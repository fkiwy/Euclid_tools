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
    """
    Open a file using the default system application based on the platform.

    Parameters
    ----------
    filename : str
        The path to the file to be opened.

    Notes
    -----
    - On Windows, this uses `os.startfile`.
    - On macOS, it uses the `open` command.
    - On Linux, it uses `evince` (assumes it's installed).
    """

    if sys.platform == "win32":
        os.startfile(filename)
    else:
        opener = "open" if sys.platform == "darwin" else "evince"
        subprocess.call([opener, filename])


def create_object_name(ra, dec, precision=0, sep="", prefix=None, shortform=False, decimal=True):
    """
    Generate a string-based object name from celestial coordinates.

    Parameters
    ----------
    ra : float
        Right Ascension in degrees.
    dec : float
        Declination in degrees.
    precision : int, optional
        Number of decimal places for coordinates (default is 0).
    sep : str, optional
        Separator to use in formatted output (default is "").
    prefix : str, optional
        String to prepend to the object name.
    shortform : bool, optional
        If True, returns a short HMS/DMS form like '1234+5678'.
    decimal : bool, optional
        If True, returns decimal-formatted coordinates.

    Returns
    -------
    str
        A formatted object name string.
    """

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
    """
    Add AB and Vega magnitudes and errors to a table based on flux values.

    Parameters
    ----------
    table : astropy.table.Table
        Table to which the magnitude columns will be added.
    flux : array_like
        Flux values in microJanskys (µJy).
    flux_err : array_like
        Flux uncertainties in microJanskys (µJy).
    band : str
        Photometric band name (e.g., 'VIS', 'Y', 'J', 'H').

    Notes
    -----
    This function adds four new columns to the table:
    - '<band>_AB_mag', '<band>_AB_err'
    - '<band>_Vega_mag', '<band>_Vega_err'
    """

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
    """
    Convert flux and its error to magnitude and magnitude error.

    Parameters
    ----------
    flux : array_like
        Flux values in microJanskys (µJy).
    flux_err : array_like
        Flux errors in microJanskys (µJy).
    magnitude_system : MagnitudeSystem
        Enum value specifying either AB or Vega system.
    band : str, optional
        Photometric band for Vega conversion. Required if `magnitude_system` is Vega.

    Returns
    -------
    tuple of array_like
        Tuple containing:
        - Magnitudes (rounded to 3 decimal places)
        - Magnitude errors (rounded to 3 decimal places)
    """

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
    """
    Pretty-print a summary of magnitude and coordinate data from a results table.

    Parameters
    ----------
    results : astropy.table.Table
        Table containing magnitude and positional data.

    Notes
    -----
    Only selected columns are printed, covering:
    - Object ID and coordinates
    - AB and Vega magnitudes and errors for VIS, Y, J, H bands
    - Separation (if present)
    """

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
    """
    Print a message when no object is found in the Euclid Q1 MER catalog.

    Parameters
    ----------
    ra : float
        Right Ascension in degrees.
    dec : float
        Declination in degrees.
    search_radius : float
        Search radius in arcseconds.

    Notes
    -----
    This message is intended to guide the user toward verifying or expanding the search.
    """

    print(
        f"""
No object found in Euclid Q1 MER catalog at the specified coordinates.
- Right Ascension: {ra}
- Declination: {dec}
- Search radius: {search_radius} arcsec

Please verify the coordinates or try increasing the search radius.
"""
    )
