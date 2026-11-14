"""Release isolation, unit scaling, all-product retrieval and failure contracts."""
import io
from unittest.mock import Mock, MagicMock, patch

import astropy.units as u
import numpy as np
import pytest
import requests
from astropy.table import QTable, Table, Column

from euclid_tools import irsa_spectra as dr1, irsa_tools
from euclid_tools.backends import get_backend

PATH = "euclid/dr1/SIR_COMBINED_SPECTRA/101839130/example.fits"
pytestmark = pytest.mark.filterwarnings("ignore:The unit .* has been deprecated:astropy.units.UnitsWarning")

ID = "-2730189294682318706"
UNIT = u.Unit(1e-16 * dr1.FLUX_UNIT)


def votable(*, blue=False, negative_variance=False):
    variance = [4., -1. if negative_variance else 9., 16.]
    table = QTable({"WAVELENGTH": ([9500, 11000, 13500] if blue else [11900, 15000, 19000]) * u.AA,
                    "SIGNAL": [-2., 0., 10.] * UNIT,
                    "VAR": variance * UNIT**2,
                    "UNCERTAINTY": [2., 3., 4.] * UNIT,
                    "MASK": [1, 64, 2], "QUALITY": [.5, .8, 1.], "NDITH": [1, 2, 3]})
    stream = io.BytesIO()
    table.write(stream, format="votable")
    return stream.getvalue()


def associations(count=1):
    return Table({"objectid": [int(ID)] * count, "tileid": [101839130] * count,
                  "path": [PATH] * count, "hdu": list(range(1, count + 1))})


def client_for(table):
    client = Mock()
    client.query_tap.return_value.to_table.return_value = table
    return client


def decode(body=None):
    return dr1.decode_spectrum(votable() if body is None else body,
                               object_id=ID, path=PATH, hdu=1, tileid=101839130)


def test_native_scaled_units_are_applied_once_and_all_bins_survive():
    spectrum = decode()
    np.testing.assert_allclose(spectrum["WAVELENGTH"].to_value(u.um), [1.19, 1.5, 1.9])
    np.testing.assert_allclose(spectrum["FLUX"].to_value(dr1.FLUX_UNIT), [-2e-16, 0, 1e-15], rtol=1e-6)
    np.testing.assert_allclose(spectrum["ERROR"].to_value(dr1.FLUX_UNIT), [2e-16, 3e-16, 4e-16], rtol=1e-6)
    assert list(spectrum["MASK"]) == [1, 64, 2]
    assert list(spectrum["NDITH"]) == [1, 2, 3]
    assert spectrum.meta["GRISM"] == "red"
    assert spectrum.meta["OBJECTID"] == ID
    assert spectrum.meta["DATASET"] == "euclid_combspec_dr1"
    assert len(spectrum.meta["RAWSHA"]) == 64


def test_invalid_variance_is_not_clipped():
    assert np.isnan(decode(votable(negative_variance=True))["ERROR"][1])


def test_all_products_and_single_product_ambiguity():
    client = client_for(associations(2))
    with patch.object(dr1, "_download", side_effect=[votable(), votable(blue=True)]):
        products = get_backend("irsa").retrieve_spectra(ID, client=client)
    assert [s.meta["GRISM"] for s in products] == ["red", "blue"]
    assert [s.meta["SRCHDU"] for s in products] == [1, 2]
    with patch.object(dr1, "_download", side_effect=[votable(), votable(blue=True)]):
        with pytest.raises(ValueError, match="Ambiguous"):
            irsa_tools.retrieve_spectrum(ID, release="dr1", client=client)
    with patch.object(dr1, "_download", side_effect=[votable(), votable(blue=True)]):
        result = get_backend("irsa").retrieve_spectrum(ID, release="DR1", grism="blue", client=client)
    assert result.meta["SRCHDU"] == 2
    query = client.query_tap.call_args.args[0]
    assert f"WHERE objectid = {ID}" in query
    assert dr1.ASSOCIATION_TABLE in query
    assert "q1" not in query


def test_repeated_same_grism_requires_exact_product_and_roundtrips(tmp_path):
    client = client_for(associations(2))
    with patch.object(dr1, "_download", return_value=votable()):
        with pytest.raises(ValueError, match="Ambiguous"):
            irsa_tools.retrieve_spectrum(ID, release="dr1", grism="red", client=client)
    with patch.object(dr1, "_download", return_value=votable()) as download:
        result = irsa_tools.retrieve_spectrum(ID, release="dr1", grism="red", product=(PATH, 2),
                                              client=client, save_spectrum=True, output_dir=str(tmp_path))
    download.assert_called_once()
    saved = QTable.read(tmp_path / f"E{ID}_spectrum.fits")
    assert saved.meta["SRCHDU"] == 2
    np.testing.assert_array_equal(saved["MASK"], result["MASK"])
    np.testing.assert_allclose(saved["FLUX"].to_value(dr1.FLUX_UNIT), result["FLUX"].to_value(dr1.FLUX_UNIT))


def test_no_products_and_unmatched_selector_are_distinct_from_errors():
    with patch.object(dr1, "_download") as download:
        assert irsa_tools.retrieve_spectrum(ID, release="dr1", client=client_for(associations(0))) is None
        assert dr1.retrieve_spectra(ID, product=(PATH, 99), client=client_for(associations())) == []
        download.assert_not_called()
    with patch.object(dr1, "_download", side_effect=requests.Timeout("test timeout")):
        with pytest.raises(requests.Timeout):
            dr1.retrieve_spectra(ID, client=client_for(associations()))
    with patch.object(dr1, "_download", side_effect=[votable(), requests.HTTPError("unavailable")]):
        with pytest.raises(requests.HTTPError):
            dr1.retrieve_spectra(ID, client=client_for(associations(2)))


@pytest.mark.parametrize("value", ["1 OR 1=1", "", "--1", "-", "1.5", "１", None, 2**63])
def test_invalid_ids_fail_before_query(value):
    client = Mock()
    with pytest.raises(ValueError):
        irsa_tools.retrieve_spectrum(value, release="dr1", client=client)
    client.query_tap.assert_not_called()


def test_bad_release_grism_metadata_and_overflow():
    client = client_for(associations())
    for kwargs in ({"release": "dr2"}, {"release": "dr1", "grism": "green"},
                   {"release": "dr1", "product": (PATH, -1)}):
        with pytest.raises(ValueError):
            irsa_tools.retrieve_spectrum(ID, client=client, **kwargs)
    client.query_tap.assert_not_called()
    with pytest.raises(ValueError, match="truncated"):
        dr1.retrieve_spectra(ID, client=client_for(associations(101)))
    rows = associations()
    rows["objectid"] = [123]
    with pytest.raises(ValueError, match="different object"):
        dr1.retrieve_spectra(ID, client=client_for(rows))
    with pytest.raises(ValueError):
        dr1.spectrum_url("euclid/dr1/../q1/file.fits", 1)


@pytest.mark.parametrize("change", ["missing_units", "wrong_units", "missing_mask", "bad_identity", "bad_error"])
def test_malformed_spectra_are_not_accepted(change):
    raw = QTable.read(io.BytesIO(votable()), format="votable")
    if change == "missing_units":
        raw.replace_column("SIGNAL", Column(np.asarray(raw["SIGNAL"].value)))
    elif change == "wrong_units":
        raw["SIGNAL"] = raw["SIGNAL"].value * u.Jy
    elif change == "missing_mask":
        raw.remove_column("MASK")
    elif change == "bad_error":
        raw["UNCERTAINTY"] *= 2
    stream = io.BytesIO()
    raw.write(stream, format="votable")
    body = stream.getvalue()
    if change == "bad_identity":
        body = body.replace(b"<TABLE>", b'<TABLE ID="Euclid_999_Combined">')
    with pytest.raises((ValueError, u.UnitConversionError)):
        decode(body)


def test_http_status_and_bounded_stream_are_enforced():
    session = MagicMock()
    response = session.get.return_value.__enter__.return_value
    response.headers = {"Content-Length": "100"}
    with pytest.raises(ValueError, match="max_bytes"):
        dr1._download("https://example.test", session=session, max_bytes=10)
    response.headers = {}
    response.iter_content.return_value = [b"123456", b"123456"]
    with pytest.raises(ValueError, match="max_bytes"):
        dr1._download("https://example.test", session=session, max_bytes=10)
    response.raise_for_status.side_effect = requests.HTTPError("503")
    with pytest.raises(requests.HTTPError):
        dr1._download("https://example.test", session=session)


def test_default_q1_retrieval_route_is_unchanged():
    client = Mock()
    client.tap_url = "https://irsa.ipac.caltech.edu/TAP"
    client.query_tap.return_value.to_table.return_value = Table({"path": ["legacy.vot"]})
    with patch.object(irsa_tools, "Irsa", client), patch.object(irsa_tools.QTable, "read", return_value=QTable.read(io.BytesIO(votable()), format="votable")), patch.object(dr1, "retrieve_spectra") as new:
        result = irsa_tools.retrieve_spectrum("123")
    assert "objectid_spectrafile_association_q1" in client.query_tap.call_args.args[0]
    assert len(result) == 3
    new.assert_not_called()


def test_public_dr1_excerpt_retains_physical_scale():
    from pathlib import Path
    path = Path(__file__).parent / "data/dr1_irsa_spectrum_excerpt.vot"
    spectrum = dr1.decode_spectrum(path.read_bytes(), object_id="2730189294682318706",
                                   path=PATH, hdu=6563, tileid=101839130)
    # Original bin 200 is row 6 in the excerpt: physical reference from the public response.
    assert float(spectrum["WAVELENGTH"][6].to_value(u.um)) == pytest.approx(1.458)
    assert float(spectrum["FLUX"][6].to_value(dr1.FLUX_UNIT)) == pytest.approx(7.498266696929932e-16, rel=1e-6, abs=0)
    assert float(spectrum["ERROR"][6].to_value(dr1.FLUX_UNIT)) == pytest.approx(1.662549190223217e-18, rel=1e-6, abs=0)
    assert spectrum["MASK"][6] == 2


@pytest.mark.parametrize("status", ["ERROR", "OVERFLOW"])
def test_service_votable_errors_are_not_valid_products(status):
    body = votable().replace(b"</RESOURCE>",
        f'<INFO name="QUERY_STATUS" value="{status}">test failure</INFO></RESOURCE>'.encode())
    with pytest.raises(ValueError, match="QUERY_STATUS"):
        decode(body)


Q1_PATH = 'euclid/q1/SIR/102158582/EUC_SIR_W-COMBSPEC_102158582_2024-11-05T16:02:47.013610Z.fits'


def test_q1_fits_association_uses_only_requested_hdu_and_preserves_photometry():
    client = Mock()
    client.query_tap.return_value.to_table.return_value = Table({'path': [Q1_PATH], 'hdu': [2338]})
    body = votable()
    raw = QTable.read(io.BytesIO(body), format='votable')
    with patch.object(irsa_tools, 'Irsa', client), patch.object(dr1, '_download', return_value=body) as download:
        result = get_backend('irsa').retrieve_spectrum('2664850113649936423', release='q1', timeout=17)
    download.assert_called_once_with(dr1.SPECTRUMDM_URL + Q1_PATH.replace(':', '%3A') + '?dataset_id=euclid_combspec_q1&hdu=2338', timeout=17)
    assert len(result) == len(raw)
    from euclid_tools.provenance import product_provenance
    provenance = product_provenance(result)
    assert provenance["release"] == "Q1"
    assert provenance["retrieval_route"] == "tap-association+spectrumdm"
    assert provenance["source_hdu"] == 2338
    assert provenance["source_file"] == Q1_PATH
    assert provenance["grism"] is None  # no bandpass metadata supplied
    np.testing.assert_allclose(result['WAVELENGTH'].to_value(u.AA), raw['WAVELENGTH'].to_value(u.AA))
    np.testing.assert_allclose(result['FLUX'].to_value(dr1.FLUX_UNIT), raw['SIGNAL'].to_value(dr1.FLUX_UNIT))
    np.testing.assert_allclose(result['ERROR'].to_value(dr1.FLUX_UNIT), raw['UNCERTAINTY'].to_value(dr1.FLUX_UNIT))
    np.testing.assert_array_equal(result['MASK'], raw['MASK'])


@pytest.mark.parametrize('paths,hdus', [([Q1_PATH], None), ([Q1_PATH], [-1]),
    ([Q1_PATH.replace('/q1/', '/dr1/')], [2338]),
    ([Q1_PATH.replace('/SIR/', '/../SIR/')], [2338]),
    ([Q1_PATH, Q1_PATH], [2338, 2339])])
def test_q1_invalid_or_ambiguous_fits_association_never_downloads(paths, hdus):
    client = Mock()
    columns = {'path': paths}
    if hdus is not None:
        columns['hdu'] = hdus
    client.query_tap.return_value.to_table.return_value = Table(columns)
    with patch.object(irsa_tools, 'Irsa', client), patch.object(dr1, '_download') as download:
        with pytest.raises(ValueError):
            irsa_tools.retrieve_spectrum('2664850113649936423')
    download.assert_not_called()
