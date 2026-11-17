"""Retrieve one documented public DR1 spectrum; never download its parent FITS.

Run: python scripts/validate_irsa_dr1_spectrum.py --output /tmp/irsa-dr1
This checks transport, units, a reference bin and FITS roundtrip, not calibration.
"""
import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import astropy
import astropy.units as u
import astroquery
import numpy as np
from astropy.table import QTable

from euclid_tools.irsa_tools import retrieve_spectrum
from euclid_tools.irsa_spectra import FLUX_UNIT

OBJECT_ID = "2730189294682318706"


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output", type=Path, required=True)
    args = cli.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    spectrum = retrieve_spectrum(OBJECT_ID, release="dr1", grism="red", timeout=120)
    assert spectrum is not None and len(spectrum) == 531
    # Reference bin measured directly in the public SpectrumDM response on 13 Nov.
    np.testing.assert_allclose(spectrum["WAVELENGTH"][200].to_value(u.um), 1.458, rtol=1e-6)
    np.testing.assert_allclose(spectrum["FLUX"][200].to_value(FLUX_UNIT), 7.498266696929932e-16, rtol=1e-6, atol=0)
    np.testing.assert_allclose(spectrum["ERROR"][200].to_value(FLUX_UNIT), 1.662549190223217e-18, rtol=1e-6, atol=0)
    assert spectrum["MASK"][200] == 2
    destination = args.output / "spectrum.fits"
    spectrum.write(destination, format="fits", overwrite=True)
    cached = QTable.read(destination)
    np.testing.assert_array_equal(cached["MASK"], spectrum["MASK"])
    np.testing.assert_allclose(cached["FLUX"].to_value(FLUX_UNIT), spectrum["FLUX"].to_value(FLUX_UNIT))
    report = {"utc": datetime.now(timezone.utc).isoformat(), "object_id": OBJECT_ID,
              "rows": len(spectrum), "astropy": astropy.__version__, "astroquery": astroquery.__version__,
              "association_endpoint": "https://irsa.ipac.caltech.edu/TAP",
              "product": spectrum.meta, "wavelength_unit": str(spectrum["WAVELENGTH"].unit),
              "flux_unit": str(spectrum["FLUX"].unit), "all_masks_preserved": True,
              "reference_bin_200": {"wavelength_um": float(spectrum["WAVELENGTH"][200].value),
                                    "flux_cgs_per_Angstrom": float(spectrum["FLUX"][200].to_value(FLUX_UNIT)),
                                    "error_cgs_per_Angstrom": float(spectrum["ERROR"][200].to_value(FLUX_UNIT))},
              "limitations": "One public red spectrum; scientific calibration and live blue/repeat ambiguity are not validated"}
    (args.output / "validation.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"Validated {len(spectrum)} bins, physical units, preserved masks and FITS roundtrip")


if __name__ == "__main__":
    main()
