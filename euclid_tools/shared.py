import os
import subprocess
import sys
from enum import Enum

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import QTable


class MagnitudeSystem(Enum):
    AB = 1
    Vega = 2


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


def add_magnitude(table, flux, flux_err, band, suffix=""):
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
    suffix : str, optional
        Measurement suffix appended to the band in output column names.

    Notes
    -----
    This function adds four new columns to the table:
    - '<band>_AB_mag', '<band>_AB_err'
    - '<band>_Vega_mag', '<band>_Vega_err'
    """

    column_prefix = band + suffix
    table[column_prefix + "_AB_mag"], table[column_prefix + "_AB_err"] = convert_flux_to_mag(
        flux, flux_err, magnitude_system=MagnitudeSystem.AB
    )
    table[column_prefix + "_AB_mag"].unit = u.mag
    table[column_prefix + "_AB_err"].unit = u.mag

    table[column_prefix + "_Vega_mag"], table[column_prefix + "_Vega_err"] = convert_flux_to_mag(
        flux, flux_err, magnitude_system=MagnitudeSystem.Vega, band=band
    )
    table[column_prefix + "_Vega_mag"].unit = u.mag
    table[column_prefix + "_Vega_err"].unit = u.mag


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

    if magnitude_system == MagnitudeSystem.AB:
        zero_point = 3631

    if magnitude_system == MagnitudeSystem.Vega:
        zero_points = {"VIS": 2835.34, "Y": 1916.10, "J": 1370.25, "H": 918.35}
        zero_point = zero_points[band]

    zero_point = (zero_point * u.Jy).to(u.uJy).value

    mag = -2.5 * np.log10(flux) + 2.5 * np.log10(zero_point)
    mag_err = 2.5 / np.log(10) * (flux_err / flux)

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


def _filled_column(value):
    if hasattr(value, "filled"):
        try:
            return value.filled(np.nan)
        except TypeError:
            return value.filled(0)
    return value


def normalize_spectrum_table(spectrum):
    """Return a spectrum QTable with canonical column names and NaN-filled masked values."""
    if spectrum is None:
        return None

    required = {"WAVELENGTH", "FLUX", "ERROR", "MASK"}
    missing = required - set(spectrum.colnames)
    if missing:
        raise ValueError(f"Spectrum is missing required columns: {sorted(missing)}")

    result = {
        "WAVELENGTH": _filled_column(spectrum["WAVELENGTH"]),
        "FLUX": _filled_column(spectrum["FLUX"]),
        "ERROR": _filled_column(spectrum["ERROR"]),
        "MASK": _filled_column(spectrum["MASK"]),
    }

    if "QUALITY" in spectrum.colnames:
        result["QUALITY"] = _filled_column(spectrum["QUALITY"])

    if "NDITH" in spectrum.colnames:
        result["NDITH"] = _filled_column(spectrum["NDITH"])

    return QTable(result)


def remove_bad_values(spectrum, mask_func):
    # Use the MASK column to create a boolean mask for values to ignore
    bad_mask = mask_func(spectrum["MASK"].value)

    # Keep only rows that are not flagged as bad
    good_mask = ~np.array(bad_mask, dtype=bool)

    # Also ensure wavelength and flux are finite for safety
    try:
        wl = spectrum["WAVELENGTH"].quantity if hasattr(spectrum["WAVELENGTH"], 'quantity') else spectrum["WAVELENGTH"]
    except Exception:
        wl = spectrum["WAVELENGTH"]

    # Build final mask: good and finite wavelength and flux
    finite_mask = np.isfinite(np.array(wl))
    combined_mask = good_mask & finite_mask

    # Slice the original table to keep only good rows (slice the table to drop rows)
    # Prefer explicit indexing to avoid returning MaskedColumns filled with NaN
    indices = np.nonzero(combined_mask)[0]

    # Extract wavelength preserving Quantity when available
    col_wl = spectrum["WAVELENGTH"]
    # Extract raw data and unit (if present on the column)
    wl_unit = getattr(col_wl, 'unit', None)
    wl_data = np.array(col_wl)[indices]

    col_flux = spectrum["FLUX"]
    fl_unit = getattr(col_flux, 'unit', None)
    fl_data = np.array(col_flux)[indices]

    col_err = spectrum["ERROR"]
    err_unit = getattr(col_err, 'unit', None)
    err_data = np.array(col_err)[indices]

    # Construct columns explicitly to preserve units when present
    from astropy.table import Column as _Column

    if wl_unit is not None:
        wl_col = _Column(data=wl_data, unit=wl_unit, name="WAVELENGTH")
    else:
        wl_col = _Column(data=wl_data, name="WAVELENGTH")

    if fl_unit is not None:
        fl_col = _Column(data=fl_data, unit=fl_unit, name="FLUX")
    else:
        fl_col = _Column(data=fl_data, name="FLUX")

    if err_unit is not None:
        err_col = _Column(data=err_data, unit=err_unit, name="ERROR")
    else:
        err_col = _Column(data=err_data, name="ERROR")

    return QTable([wl_col, fl_col, err_col])


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
