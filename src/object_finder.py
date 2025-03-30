import numpy as np
import pyvo as vo
import astropy.units as u
from astropy.coordinates import SkyCoord
from enum import Enum
import shared as shr


class MagnitudeSystem(Enum):
    AB = 'AB'
    Vega = 'Vega'


def find_objects(ra, dec, search_radius=5, nearest_object=True):
    search_radius = search_radius * u.arcsec.to(u.deg)

    adql = f"""
    SELECT *
      FROM {shr.table_mer}
     WHERE CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra}, {dec}, {search_radius})) = 1
    """

    service = vo.dal.TAPService(shr.irsa_url + 'TAP')
    result = service.search(adql)
    table = result.to_table()

    if len(table) == 0:
        return None

    def add_magnitude(flux, flux_err, band):
        table[band + '_AB_mag'], table[band + '_AB_err'] = convert_flux_to_mag(
            flux, flux_err, magnitude_system=MagnitudeSystem.AB)
        table[band + '_AB_mag'].unit = u.mag
        table[band + '_AB_err'].unit = u.mag

        table[band + '_Vega_mag'], table[band + '_Vega_err'] = convert_flux_to_mag(
            flux, flux_err, magnitude_system=MagnitudeSystem.Vega, band=band)
        table[band + '_Vega_mag'].unit = u.mag
        table[band + '_Vega_err'].unit = u.mag

    add_magnitude(result['flux_vis_psf'], result['fluxerr_vis_psf'], 'VIS')
    add_magnitude(result['flux_y_templfit'], result['fluxerr_y_templfit'], 'Y')
    add_magnitude(result['flux_j_templfit'], result['fluxerr_j_templfit'], 'J')
    add_magnitude(result['flux_h_templfit'], result['fluxerr_h_templfit'], 'H')

    # Create a SkyCoord object for the target
    target_coord = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='icrs')

    # Create a SkyCoord object for the objects in the result table
    object_coords = SkyCoord(table['ra'], table['dec'], unit=(u.deg, u.deg), frame='icrs')

    # Calculate the angular separation between the target and each object in the result table
    separations = target_coord.separation(object_coords).value

    # Convert speparations from degree to arcsec
    sparations = (separations * u.deg).to(u.arcsec)

    # Add the separations to the result table
    table['separation'] = np.round(sparations, 3)

    table.sort('separation')

    if nearest_object:
        # Find the object with the minimum separation (nearest object)
        return table[:1]
    else:
        return table


def print_catalog_info():
    service = vo.dal.TAPService(shr.irsa_url + 'TAP')
    table = service.tables[shr.table_mer]
    for col in table.columns:
        print(f'{f"{col.name}":45s} {f"{col.unit}":12s} {col.description}')


def convert_flux_to_mag(flux, flux_err, magnitude_system, band=None):
    zero_points = {
        'VIS': 2835.34,
        'Y': 1916.10,
        'J': 1370.25,
        'H': 918.35
    }

    def flux_to_mag(flux, zero_point):
        return -2.5 * np.log10(flux) + 2.5 * np.log10(zero_point)

    def flux_err_to_mag(flux, flux_err):
        return 2.5 / np.log(10) * (flux_err / flux)

    zero_point = None

    if magnitude_system == MagnitudeSystem.AB:
        zero_point = 3631

    if magnitude_system == MagnitudeSystem.Vega:
        zero_point = zero_points[band]

    zero_point = (zero_point * u.Jy).to(u.uJy).value

    mag = flux_to_mag(flux, zero_point)
    mag_err = flux_err_to_mag(flux, flux_err)

    return np.round(mag, 3), np.round(mag_err, 3)
