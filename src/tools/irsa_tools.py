import urllib
import numpy as np
import pyvo as vo
from enum import Enum
import astropy.units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.table import Table, QTable, MaskedColumn
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astroquery.ipac.irsa import Irsa

from tools.shared import MaskType, add_magnitude


COLUMNS_TO_REMOVE = ["tileid", "x", "y", "z", "spt_ind", "htm20", "cntr"]

TABLE_MER = "euclid_q1_mer_catalogue"


def retrieve_objects(ra: float, dec: float, radius: float) -> Table:
    """
    Perform a cone search on the Euclid archive.

    Parameters:
    ra (float): Right ascension in degrees.
    dec (float): Declination in degrees.
    radius (float): Search radius in arcseconds.

    Returns:
    Table: Astropy table containing catalog entries.
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


def retrieve_spectrum(object_id: str, maskType: Enum = MaskType.NONE) -> QTable:
    """
    Retrieve the spectrum for a given object ID.

    Parameters:
    object_id (str): The ID of the object to retrieve the spectrum for.

    Returns:
    QTable: A QTable containing the wavelength, flux, and error data.
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

        if maskType in [MaskType.FLUX, MaskType.ERROR, MaskType.BOTH]:
            # Use the MASK column to create a boolean mask for values to ignore
            bad_mask = (spectrum["MASK"].value % 2 == 1) | (spectrum["MASK"].value >= 64)

            if maskType in [MaskType.FLUX, MaskType.BOTH]:
                # Apply the mask to the flux values
                flux = MaskedColumn(flux, mask=bad_mask)

            if maskType in [MaskType.ERROR, MaskType.BOTH]:
                # Apply the mask to the error values
                error = MaskedColumn(error, mask=bad_mask)

        # Create a new QTable for the result
        result = QTable([wavelength, flux, error], names=("WAVELENGTH", "FLUX", "ERROR"))

    return result


def retrieve_cutout(ra: float, dec: float, search_radius: float, cutout_size: float, band: str) -> fits.HDUList:
    """
    Retrieve a cutout image from the Euclid archive.

    Parameters:
    ra (float): Right ascension in degrees.
    dec (float): Declination in degrees.
    search_radius (float): Search radius in arcseconds.
    cutout_size (float): Size of the cutout in arcseconds.
    band (str): Euclid band ('VIS', 'Y', 'J', 'H').

    Returns:
    fits.HDUList: HDU list containing the image cutout.
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
    service = vo.dal.TAPService("https://irsa.ipac.caltech.edu/TAP")
    table = service.tables[TABLE_MER]

    for col in table.columns:
        if col.name in COLUMNS_TO_REMOVE:
            continue

        print(f'{f"{col.name}":45s} {f"{col.unit}":12s} {col.description}')
