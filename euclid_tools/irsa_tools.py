import os
import tempfile
import urllib
from typing import Callable

import astropy.units as u
import numpy as np
import pyvo as vo
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.table import Table, QTable, MaskedColumn
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

    # Replace masked values with NaN
    for col in table.colnames:
        if isinstance(table[col], np.ma.MaskedArray) and np.issubdtype(table[col].dtype, np.floating):  # Only for float-type columns
            table[col] = table[col].filled(np.nan)

    return table


def retrieve_spectrum(object_id: str, save_spectrum: bool = False, output_dir: str = tempfile.gettempdir(), object_name: str = None) -> QTable:
    """
    Retrieve the 1D spectrum for a given object ID from the Euclid archive.

    Parameters
    ----------
    object_id : str
        Unique identifier of the source.
    save_spectrum : bool, optional
        If True, saves the spectrum to a FITS file in the specified output directory.
    output_dir : str, optional
        Directory where the spectrum FITS file will be saved if `save_spectrum` is True.
    object_name : str, optional
        Name of the object to use for the output file name. If not provided, defaults to 'E{object_id}'.

    Returns
    -------
    QTable
        A table containing wavelength (microns), flux, flux error, and mask.
        Flux and error are scaled using the FITS header FSCALE keyword.

    Notes
    -----
    - Wavelength is converted to microns.
    - Data is masked based on the bitmask in the `MASK` column if requested.
    - If no spectrum is found, returns `None`.
    - The output FITS file is named using the object's RA and DEC coordinates, formatted as 'E{ra}{dec}_spectrum.fits'.
    - The output directory defaults to a temporary directory if not specified.
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

        # Check if the spectrum contains valid data
        if np.isnan(spectrum["WAVELENGTH"]).all():
            return None

        # Convert wavelength to microns
        wavelength = spectrum["WAVELENGTH"].to(u.micron)

        # Scale the flux and error using the scaling factor from the fits header
        flux = spectrum["SIGNAL"] * spec_header["FSCALE"]
        error = np.sqrt(spectrum["VAR"]) * spec_header["FSCALE"]
        mask = spectrum["MASK"]

        # Create a new QTable for the result
        result = QTable([wavelength, flux, error, mask], names=("WAVELENGTH", "FLUX", "ERROR", "MASK"))

        if save_spectrum:
            if not object_name:
                # If no object name is provided, create one based on the object ID
                object_name = f"E{object_id}"

            # Create the full file path for the plot image
            output_filename = os.path.join(output_dir, object_name + "_spectrum.fits")

            # Save the spectrum to a FITS file
            result.write(output_filename, format="fits", overwrite=True)

    return result


def mask_bad_values(spectrum: QTable, mask_func: Callable[[np.ndarray], np.ndarray]) -> QTable:
    """
    Mask bad values in the spectrum data based on the MASK column.

    Parameters
    ----------
    spectrum : QTable
        The input spectrum table containing columns WAVELENGTH, FLUX, ERROR, and MASK.
    mask_func : Callable[[np.ndarray], np.ndarray]
        A function that takes the MASK column values and returns a boolean mask.
        The function should return `True` for values to ignore (bad values) and `False` for good values.
        example: `lambda mask: (mask % 2 == 1) | (mask >= 64)`.

    Returns
    -------
    QTable
        A new table with masked values replaced by NaN. The MASK column is not included in the output.

    Notes
    -----
    - The MASK column is used to identify bad values.
    - The WAVELENGTH, FLUX, and ERROR columns are masked based on the provided mask function.
    - The output table contains the same columns but with masked values replaced by NaN.
    """

    # Use the MASK column to create a boolean mask for values to ignore
    bad_mask = mask_func(spectrum["MASK"].value)

    # Apply the mask to the spectrum data
    wavelength = MaskedColumn(spectrum["WAVELENGTH"], mask=bad_mask)
    flux = MaskedColumn(spectrum["FLUX"], mask=bad_mask)
    error = MaskedColumn(spectrum["ERROR"], mask=bad_mask)

    # Replace masked values by NaN
    wavelength = wavelength.filled(np.nan)
    flux = flux.filled(np.nan)
    error = error.filled(np.nan)

    return QTable([wavelength, flux, error], names=("WAVELENGTH", "FLUX", "ERROR"))


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
