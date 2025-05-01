import urllib
import numpy as np
import pyvo as vo
import astropy.units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.table import Table, QTable, MaskedColumn
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astroquery.ipac.irsa import Irsa

from euclid_tools.shared import add_magnitude


COLUMNS_TO_REMOVE = ["tileid", "x", "y", "z", "spt_ind", "htm20", "cntr"]

TABLE_MER = "euclid_q1_mer_catalogue"


def retrieve_objects(ra: float, dec: float, radius: float) -> Table:
    """
    Perform a cone search on the Euclid Q1 MER catalog and return nearby objects.

    Parameters
    ----------
    ra : float
        Right Ascension in degrees.
    dec : float
        Declination in degrees.
    radius : float
        Search radius in arcseconds.

    Returns
    -------
    Table
        Astropy Table containing catalog entries with added AB and Vega magnitudes.
        Returns `None` if no objects are found.

    Notes
    -----
    - Columns `tileid`, `x`, `y`, `z`, `spt_ind`, `htm20`, `cntr`are removed.
    - Adds VIS, Y, J, and H magnitudes using `add_magnitude`.
    - Results are sorted by angular separation in arcseconds.
    """

    position = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame="icrs")
    radius = radius * u.arcsec

    table = Irsa.query_region(coordinates=position, spatial="Cone", catalog=TABLE_MER, radius=radius)

    if len(table) == 0:
        return None

    table.remove_columns(COLUMNS_TO_REMOVE)

    add_magnitude(table, table["flux_vis_psf"], table["fluxerr_vis_psf"], "VIS")
    add_magnitude(table, table["flux_y_templfit"], table["fluxerr_y_templfit"], "Y")
    add_magnitude(table, table["flux_j_templfit"], table["fluxerr_j_templfit"], "J")
    add_magnitude(table, table["flux_h_templfit"], table["fluxerr_h_templfit"], "H")

    # Create a SkyCoord object for the target
    target_coord = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame="icrs")

    # Create a SkyCoord object for the objects in the result table
    object_coords = SkyCoord(table["ra"], table["dec"], unit=(u.deg, u.deg), frame="icrs")

    # Calculate the angular separation between the target and each object in the result table
    separations = target_coord.separation(object_coords).value

    # Convert speparations from degree to arcsec
    sparations = (separations * u.deg).to(u.arcsec)

    # Add the separations to the result table
    table["separation"] = np.round(sparations, 3)

    # Sort the result table on separation
    table.sort("separation")

    return table


def retrieve_spectrum(object_id: str, mask_bad_values: bool = False) -> QTable:
    """
    Retrieve the 1D spectrum for a given object ID from the Euclid archive.

    Parameters
    ----------
    object_id : str
        Unique identifier of the source.
    mask_bad_values : bool, default False
        Whether to mask bad flux values.

    Returns
    -------
    QTable
        A table containing wavelength (microns), flux, and flux error.
        Flux and error are scaled using the FITS header FSCALE keyword.

    Notes
    -----
    - Wavelength is converted to microns.
    - Data is masked based on the bitmask in the `MASK` column if requested.
    """

    adql = f"""
    SELECT *
      FROM euclid.objectid_spectrafile_association_q1
     WHERE objectid = {object_id}
       AND uri IS NOT NULL
    """

    table = Irsa.query_tap(adql).to_table()

    if len(table) == 0:
        return None

    file_url = urllib.parse.urljoin(Irsa.tap_url, table["uri"][0])

    with fits.open(file_url) as hdul:
        hdu = hdul[table["hdu"][0]]
        spectrum = QTable.read(hdu, format="fits")
        spec_header = hdu.header

        # Convert wavelength to microns
        wavelength = spectrum["WAVELENGTH"].to(u.micron)

        # Scale the flux and error using the scaling factor from the fits header
        flux = spectrum["SIGNAL"] * spec_header["FSCALE"]
        error = np.sqrt(spectrum["VAR"]) * spec_header["FSCALE"]

        if mask_bad_values:
            # Use the MASK column to create a boolean mask for values to ignore
            bad_mask = (spectrum["MASK"].value % 2 == 1) | (spectrum["MASK"].value >= 64)

            # Apply the mask to the spectrum data
            wavelength = MaskedColumn(wavelength, mask=bad_mask)
            flux = MaskedColumn(flux, mask=bad_mask)
            error = MaskedColumn(error, mask=bad_mask)

            # Replace masked values by NaN
            wavelength = wavelength.filled(np.nan)
            flux = flux.filled(np.nan)
            error = error.filled(np.nan)

        # Create a new QTable for the result
        result = QTable([wavelength, flux, error], names=("WAVELENGTH", "FLUX", "ERROR"))

    return result


def retrieve_cutout(ra: float, dec: float, search_radius: float, cutout_size: float, band: str) -> fits.PrimaryHDU:
    """
    Retrieve an image cutout from Euclid imaging data.

    Parameters
    ----------
    ra : float
        Right Ascension in degrees.
    dec : float
        Declination in degrees.
    search_radius : float
        Radius used to find the overlapping mosaic (in arcseconds).
    cutout_size : float
        Size of the image cutout in arcseconds.
    band : str
        Observing band to extract cutout from ('VIS', 'Y', 'J', 'H').

    Returns
    -------
    fits.HDU
        FITS HDU containing the image cutout.
        Returns `None` if no suitable image is found.

    Raises
    ------
    ValueError
        If the provided band is not one of the allowed options.
    """

    if band not in ["VIS", "Y", "J", "H"]:
        raise ValueError("Invalid band. Choose from 'VIS', 'Y', 'J', or 'H'.")

    position = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame="icrs")
    search_radius *= u.arcsec

    table = Irsa.query_sia(pos=(position, search_radius), collection="euclid_DpdMerBksMosaic")

    results = table[(table["dataproduct_subtype"] == "science")]

    if not results or len(results) == 0:
        return None

    image_url = results[results["energy_bandpassname"] == band]["access_url"][0]

    print(f"Downloading {band}-band cutout")
    with fits.open(image_url, use_fsspec=True) as hdul:
        wcs = WCS(hdul[0].header)
        cutout = Cutout2D(hdul[0].section, position=position, size=cutout_size * u.arcsec, wcs=wcs)
        hdu = fits.PrimaryHDU(data=cutout.data, header=cutout.wcs.to_header())

    return hdu


def print_catalog_info():
    """
    Print metadata information for each column in the Euclid MER catalog.

    Notes
    -----
    - Skips columns listed in `COLUMNS_TO_REMOVE`.
    - Prints column name, unit, and description.
    """

    service = vo.dal.TAPService("https://irsa.ipac.caltech.edu/TAP")
    table = service.tables[TABLE_MER]

    for col in table.columns:
        if col.name in COLUMNS_TO_REMOVE:
            continue

        print(f'{f"{col.name}":45s} {f"{col.unit}":12s} {col.description}')
