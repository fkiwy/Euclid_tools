"""DR1 profile photometry with optional, explicitly supplied Vega calibration.

No official release-specific Vega zero points are bundled. The legacy Q1
conversion in shared.py is deliberately separate from this opt-in contract.
"""
from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping
from numbers import Real

import astropy.units as u
import numpy as np

from euclid_tools.releases import normalize_release
from euclid_tools.shared import MagnitudeSystem, convert_flux_to_mag

BANDS = {"VIS": "psf", "Y": "templfit", "J": "templfit", "H": "templfit"}


def validate_vega_calibration(calibration, *, release="dr1"):
    """Validate and copy a user-supplied DR1 Vega convention before querying.

    ``zero_points_jy`` must give positive finite zero-magnitude flux densities
    for VIS/Y/J/H. ``source``, ``passband_reference`` and ``reference_spectrum``
    identify the adopted convention; supplying them is not official validation.
    """
    release = normalize_release(release)
    if calibration is None:
        return None
    if release != "dr1":
        raise ValueError("vega_calibration is a DR1-only opt-in; Q1 retains its legacy conversion")
    keys = {"release", "zero_points_jy", "source", "passband_reference", "reference_spectrum"}
    if not isinstance(calibration, Mapping) or set(calibration) != keys:
        raise ValueError(f"vega_calibration must specify exactly {sorted(keys)}")
    if calibration["release"] != "DR1":
        raise ValueError("Vega calibration release must be DR1")
    for key in ("source", "passband_reference", "reference_spectrum"):
        if not isinstance(calibration[key], str) or not calibration[key].strip():
            raise ValueError(f"Vega calibration requires a nonempty {key}")
    values = calibration["zero_points_jy"]
    if not isinstance(values, Mapping) or set(values) != set(BANDS):
        raise ValueError("zero_points_jy must specify VIS, Y, J and H")
    for value in values.values():
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real) or not np.isfinite(value) or value <= 0:
            raise ValueError("Vega zero_points_jy must be positive finite numbers in Jy")
    return {**calibration, "zero_points_jy": {band: float(values[band]) for band in BANDS}}


def add_dr1_magnitudes(table, *, vega_calibration=None):
    """Add profile AB columns and, on explicit calibration, profile Vega columns.

    Fluxes are in microJy. Rounding and measurement-error propagation match the
    existing AB outputs. Vega zero points are treated as fixed: their calibration
    uncertainties are not included in the reported measurement errors.
    """
    calibration = validate_vega_calibration(vega_calibration)
    for band, method in BANDS.items():
        flux = table[f"flux_{band.lower()}_{method}"]
        error = table[f"fluxerr_{band.lower()}_{method}"]
        mag, err = convert_flux_to_mag(flux, error, MagnitudeSystem.AB)
        table[f"{band}_AB_mag"], table[f"{band}_AB_err"] = mag, err
        table[f"{band}_AB_mag"].unit = table[f"{band}_AB_err"].unit = u.mag
        if calibration is not None:
            # This convention changes only the zero-magnitude reference flux.
            # It is not a recalibration of the archive's measured fluxes.
            f = np.ma.array(u.Quantity(flux, u.uJy).to_value(u.uJy),
                            mask=getattr(flux, 'mask', False))
            f0 = calibration["zero_points_jy"][band] * 1e6
            with np.errstate(divide="ignore", invalid="ignore"):
                vega = np.round(-2.5 * np.log10(f / f0), 3)
            table[f"{band}_Vega_mag"], table[f"{band}_Vega_err"] = vega, err.copy()
            for suffix in ("mag", "err"):
                column = table[f"{band}_Vega_{suffix}"]
                column.unit = u.mag
                column.description = ("User-supplied DR1 Vega convention; see VEGACAL. "
                                      "Errors include measurement uncertainty only.")
    if calibration is not None:
        encoded = json.dumps(calibration, sort_keys=True, separators=(",", ":"), allow_nan=False)
        table.meta.update(VEGACAL=encoded, VEGASHA=hashlib.sha256(encoded.encode()).hexdigest(),
                          VEGASTAT="user-supplied; DR1 applicability not independently validated")
    return table
