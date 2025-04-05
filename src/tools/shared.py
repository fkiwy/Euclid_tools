import os
import sys
import subprocess
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from enum import Enum


class MagnitudeSystem(Enum):
    AB = 1
    Vega = 2


class MaskType(Enum):
    NONE = 0
    FLUX = 1
    ERROR = 2
    BOTH = 3


def open_file(filename):
    if sys.platform == "win32":
        os.startfile(filename)
    else:
        opener = "open" if sys.platform == "darwin" else "evince"
        subprocess.call([opener, filename])


def create_object_name(ra, dec, precision=0, sep="", prefix=None, shortform=False, decimal=True):
    coords = SkyCoord(ra=ra * u.degree, dec=dec * u.degree)

    if shortform:
        coords_str = coords.to_string("hmsdms", decimal=False, sep=sep, precision=0)
        object_name = coords_str[0:4] + coords_str[7:12]
    else:
        if decimal:
            object_name = coords.to_string(decimal=True, precision=precision)
        else:
            object_name = coords.to_string("hmsdms", decimal=False, sep=sep, precision=precision).replace(" ", "")

    if prefix:
        object_name = prefix + object_name

    return str(object_name)


def add_magnitude(table, flux, flux_err, band):
    table[band + "_AB_mag"], table[band + "_AB_err"] = convert_flux_to_mag(
        flux, flux_err, magnitude_system=MagnitudeSystem.AB
    )
    table[band + "_AB_mag"].unit = u.mag
    table[band + "_AB_err"].unit = u.mag

    table[band + "_Vega_mag"], table[band + "_Vega_err"] = convert_flux_to_mag(
        flux, flux_err, magnitude_system=MagnitudeSystem.Vega, band=band
    )
    table[band + "_Vega_mag"].unit = u.mag
    table[band + "_Vega_err"].unit = u.mag


def convert_flux_to_mag(flux, flux_err, magnitude_system, band=None):
    zero_points = {"VIS": 2835.34, "Y": 1916.10, "J": 1370.25, "H": 918.35}

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


def print_results(results):
    table = results[
        "object_id",
        "ra",
        "dec",
        "VIS_AB_mag",
        "VIS_AB_err",
        "VIS_Vega_mag",
        "VIS_Vega_err",
        "Y_AB_mag",
        "Y_AB_err",
        "Y_Vega_mag",
        "Y_Vega_err",
        "J_AB_mag",
        "J_AB_err",
        "J_Vega_mag",
        "J_Vega_err",
        "H_AB_mag",
        "H_AB_err",
        "H_Vega_mag",
        "H_Vega_err",
        "separation",
    ]
    table.pprint_all()


def print_message(ra, dec, search_radius):
    print(
        f"""
No object found in Euclid Q1 MER catalog at the specified coordinates.
- Right Ascension: {ra}
- Declination: {dec}
- Search radius: {search_radius} arcsec

Please verify the coordinates or try increasing the search radius.
"""
    )
