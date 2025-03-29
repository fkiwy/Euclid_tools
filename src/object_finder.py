import numpy as np
import pyvo as vo
import astropy.units as u
from astropy.coordinates import SkyCoord
import shared as shr
from enum import Enum


class MagnitudeSystem(Enum):
    AB = 'AB'
    Vega = 'Vega'


def find_object(ra, dec, search_radius=5):
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

    # Create a SkyCoord object for the target
    target_coord = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='icrs')

    # Create a SkyCoord object for the objects in the result table
    object_coords = SkyCoord(table['ra'], table['dec'], unit=(u.deg, u.deg), frame='icrs')

    # Calculate the angular separation between the target and each object in the result table
    separations = target_coord.separation(object_coords).value

    # Add the separations to the result table
    table['separation'] = separations

    # Find the object with the minimum separation (closest object)
    closest_object = table[separations.argmin()]

    return closest_object


def print_catalog_info():
    service = vo.dal.TAPService(shr.irsa_url + 'TAP')
    table = service.tables[shr.table_mer]
    for col in table.columns:
        print(f'{f"{col.name}":45s} {f"{col.unit}":12s} {col.description}')


def convert_flux_to_mag(flux, band=None, magnitude_system=MagnitudeSystem.AB):
    zero_points = {
        'VIS': 2835.34,
        'Y': 1916.10,
        'J': 1370.25,
        'H': 918.35
    }

    def flux_to_mag(flux, zero_point):
        return -2.5 * np.log10(flux) + 2.5 * np.log10(zero_point)

    zero_point = None

    if magnitude_system == MagnitudeSystem.AB:
        zero_point = 3631

    if magnitude_system == MagnitudeSystem.Vega:
        zero_point = zero_points[band]

    zero_point = (zero_point * u.Jy).to(u.uJy).value

    magnitude = flux_to_mag(flux, zero_point)

    return np.round(magnitude, 3)
