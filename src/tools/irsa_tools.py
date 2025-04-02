import requests
from io import BytesIO
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
import pyvo as vo
from astropy.io import fits
import astropy.units as u

import tools.shared as shr


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
    search_radius = radius * u.arcsec.to(u.deg)

    adql = f"""
    SELECT *
      FROM {shr.table_mer}
     WHERE CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra}, {dec}, {search_radius})) = 1
    """

    service = vo.dal.TAPService(shr.irsa_url + "TAP")
    results = service.search(adql)
    table = results.to_table()

    if len(table) == 0:
        return None

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


def retrieve_spectrum(object_id: str) -> fits.HDUList:
    """
    Retrieve the spectrum for a given object ID.

    Parameters:
    object_id (str): The ID of the object to retrieve the spectrum for.

    Returns:
    fits.HDUList: HDU list containing the spectrum data.
    """
    adql = f"""
    SELECT *
      FROM {shr.table_1dspectra}
     WHERE objectid = {object_id}
       AND uri IS NOT NULL
    """

    service = vo.dal.TAPService(shr.irsa_url + "TAP")
    results = service.search(adql)
    table = results.to_table()

    if len(table) == 0:
        return None

    file_url = shr.irsa_url + table["uri"][0]
    response = requests.get(file_url)

    with fits.open(BytesIO(response.content), memmap=True) as hdul:
        return hdul[table["hdu"][0]]


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

    def create_cutout(cutout_coords, cutout_size, image_url, band):
        print(f"Downloading {band} cutout")
        hdu = fits.open(image_url, use_fsspec=True)
        header = hdu[0].header
        cutout_data = Cutout2D(hdu[0].section, position=cutout_coords, size=cutout_size, wcs=WCS(hdu[0].header))
        hdu.close()
        new_hdu = fits.PrimaryHDU(data=cutout_data.data, header=header)
        new_hdu.header.update(cutout_data.wcs.to_header())
        return new_hdu

    def get_image_url(table, band):
        return table[table["energy_bandpassname"] == band]["access_url"][0]

    position = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame="icrs")
    search_radius *= u.arcsec
    cutout_size *= u.arcsec

    irsa_service = vo.dal.sia2.SIA2Service(shr.irsa_url + "SIA")

    table = irsa_service.search(pos=(position, search_radius), collection="euclid_DpdMerBksMosaic").to_table()

    results = table[(table["dataproduct_subtype"] == "science")]

    if not results or len(results) == 0:
        return None

    return create_cutout(position, cutout_size, get_image_url(results, band), band)


def print_catalog_info():
    service = vo.dal.TAPService(shr.irsa_url + "TAP")
    table = service.tables[shr.table_mer]
    for col in table.columns:
        print(f'{f"{col.name}":45s} {f"{col.unit}":12s} {col.description}')


def add_magnitude(table, flux, flux_err, band):
    table[band + "_AB_mag"], table[band + "_AB_err"] = shr.convert_flux_to_mag(
        flux, flux_err, magnitude_system=shr.MagnitudeSystem.AB
    )
    table[band + "_AB_mag"].unit = u.mag
    table[band + "_AB_err"].unit = u.mag

    table[band + "_Vega_mag"], table[band + "_Vega_err"] = shr.convert_flux_to_mag(
        flux, flux_err, magnitude_system=shr.MagnitudeSystem.Vega, band=band
    )
    table[band + "_Vega_mag"].unit = u.mag
    table[band + "_Vega_err"].unit = u.mag
