"""Exercise the masked-column operation used by older Matplotlib fill_between."""
from pathlib import Path

import astropy.units as u
import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import MaskedColumn, Table

from euclid_tools.irsa_spectra import decode_spectrum
from flux_comp.core import SED, WaveFlux


def guard_old_fill_between(monkeypatch):
    original = plt.fill_between
    calls = []

    def fill(x, lower, upper, **kwargs):
        # Reproduce the stack/slice that crashes on Astropy MaskedColumn in
        # the affected Matplotlib path, even if the installed version handles it.
        np.stack((x, lower), axis=-1)[::-1]
        np.stack((x, upper), axis=-1)[::-1]
        calls.append((x, lower, upper))
        return original(x, lower, upper, **kwargs)

    monkeypatch.setattr(plt, "fill_between", fill)
    return calls


def test_irsa_example_pipeline_renders_uncertainty_pdf(monkeypatch, tmp_path):
    path = Path(__file__).parent / "data/dr1_irsa_spectrum_excerpt.vot"
    data = decode_spectrum(path.read_bytes(), object_id="2730189294682318706",
                           path="euclid/dr1/SIR_COMBINED_SPECTRA/101839130/example.fits", hdu=6563, tileid=101839130)
    waveflux = WaveFlux(label="Spectrum", wavelength=data["WAVELENGTH"],
                        flux=data["FLUX"], uncertainty=data["ERROR"])
    waveflux.trim(1.21, 1.89)
    sed = SED("DR1 regression")
    sed.add(waveflux)
    sed.to_flux_lambda()
    before = waveflux.data.copy()
    calls = guard_old_fill_between(monkeypatch)
    sed.plot(spec_uncertainty=True, figure_size=(8, 5), show_grid=True,
             plot_format="pdf", output_dir=str(tmp_path), open_plot=False, detect_gaps=False)
    assert (tmp_path / "DR1_regression.pdf").read_bytes().startswith(b"%PDF")
    assert len(calls) == 1
    x, lower, upper = calls[0]
    assert not any(isinstance(array, MaskedColumn) for array in (x, lower, upper))
    np.testing.assert_allclose(np.asarray(x), before["Wavelength"])
    np.testing.assert_allclose(np.asarray(lower), before["Flux"] - before["Uncertainty"], rtol=1e-6)
    np.testing.assert_allclose(np.asarray(upper), before["Flux"] + before["Uncertainty"], rtol=1e-6)
    for name in before.colnames:
        np.testing.assert_array_equal(waveflux.data[name], before[name])
        assert waveflux.data[name].unit == before[name].unit


def test_masked_flux_and_errors_stay_masked_at_plot_boundary(monkeypatch, tmp_path):
    waveflux = WaveFlux()
    waveflux.photometry = False
    waveflux.data = Table([
        MaskedColumn([1., 1.1, 1.2, 1.3], mask=[False, True, False, False], unit=u.um),
        MaskedColumn([-2., 99., 0., 3.], mask=[False, True, False, False], unit=u.Jy),
        MaskedColumn([.2, .3, .4, .5], mask=[False, False, True, False], unit=u.Jy),
    ], names=["Wavelength", "Flux", "Uncertainty"])
    sed = SED("Masked regression")
    sed.add(waveflux)
    before = waveflux.data.copy()
    calls = guard_old_fill_between(monkeypatch)
    sed.plot(spec_uncertainty=True, plot_format="png", output_dir=str(tmp_path),
             open_plot=False, detect_gaps=False, offset=1)
    assert (tmp_path / "Masked_regression.png").stat().st_size > 0
    x, lower, upper = calls[0]
    np.testing.assert_array_equal(np.ma.getmaskarray(x), np.ma.getmaskarray(before["Wavelength"]))
    expected_mask = np.ma.getmaskarray(before["Flux"]) | np.ma.getmaskarray(before["Uncertainty"])
    np.testing.assert_array_equal(np.ma.getmaskarray(lower), expected_mask)
    np.testing.assert_array_equal(np.ma.getmaskarray(upper), expected_mask)
    np.testing.assert_allclose([lower[0], upper[0]], [-1.2, -.8])
    for name in before.colnames:
        np.testing.assert_array_equal(np.ma.getmaskarray(waveflux.data[name]), np.ma.getmaskarray(before[name]))
        np.testing.assert_array_equal(np.asarray(waveflux.data[name]), np.asarray(before[name]))
