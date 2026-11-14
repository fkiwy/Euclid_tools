"""Small in-memory mosaics test selection, range reads, WCS and calibration."""
import io
import json
from unittest.mock import Mock, patch

import astropy.units as u
import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

from euclid_tools import irsa_images as images, irsa_tools
from euclid_tools.backends import get_backend

RA, DEC = 273.167, 64.3687


def mosaic(kind="survey", band="VIS", tile=None, version="one", **changes):
    tile = tile or ("101829777" if kind == "deep" else "102158277")
    instrument = "VIS" if band == "VIS" else "NISP"
    path = f"s3://{images.BUCKET}/dr1/MER/{tile}/{instrument}/{version}.fits"
    w = WCS(naxis=2)
    w.wcs.crpix = [450, 450]
    w.wcs.crval = [RA, DEC]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    header = w.to_header()
    # A rotated CD matrix ensures the crop does not retain stale CD/PC mixtures.
    header.pop("CDELT1", None); header.pop("CDELT2", None)
    matrix = np.array([[-.9, .435889894], [.435889894, .9]]) / 36000
    for i in range(2):
        for j in range(2):
            header[f"CD{i+1}_{j+1}"] = matrix[i, j]
    pipeline = "MER_ProcessTile_EUCLID_2.0.32"
    plan = pipeline + ("-DATA_RELEASE-PLAN-test" if kind == "survey" else "-DATA_RELEASE_DEEP-PLAN-test")
    header.update(DATASETR="DR1_R1", PPLANID=plan, PIPDEFID=pipeline, PSWREL="11.3.1",
                  FILTER="VIS" if band == "VIS" else "NIR_" + band,
                  BUNIT="ADU/s" if band == "VIS" else "ELECTRON/s",
                  MAGZERO=24.5 if band == "VIS" else 29.8, EXPTIME=1234., GAIN=2.)
    header.update(changes)
    data = np.arange(900*900, dtype=np.float32).reshape(900, 900)
    hdu = fits.PrimaryHDU(data, header)
    stream = io.BytesIO(); hdu.writeto(stream, checksum=True)
    row = {"dataproduct_subtype": "science", "energy_bandpassname": band,
           "obs_id": tile + "_" + version, "cloud_access": json.dumps({"aws": {
               "bucket_name": images.BUCKET, "key": path.split(images.BUCKET + "/")[1], "region": "us-east-1"}})}
    return path, stream.getvalue(), row, hdu


class MemoryFS:
    def __init__(self, entries):
        self.entries = entries
        self.events = []

    def open(self, path, mode, **kwargs):
        assert mode == "rb" and kwargs["cache_type"] == "none"
        blob = self.entries[path]  # KeyError models an unavailable S3 object.
        events = self.events
        class File(io.BytesIO):
            def read(self, size=-1):
                events.append((path, self.tell(), size))
                assert size >= 0, "No unbounded parent-mosaic read"
                return super().read(size)
        return File(blob)


def setup(*products):
    rows = Table(rows=[p[2] for p in products])
    client = Mock()
    client.query_sia.return_value = rows
    return client, MemoryFS({p[0]: p[1] for p in products})


def test_section_read_preserves_pixels_wcs_and_calibration():
    p = mosaic()
    client, fs = setup(p)
    cutout = get_backend("irsa").retrieve_cutout(RA, DEC, search_radius=2, cutout_size=3, band="VIS",
                                                release="dr1", client=client, filesystem=fs)
    assert cutout.data.shape == (30, 30)
    x, y = cutout.header["SRCX0"], cutout.header["SRCY0"]
    np.testing.assert_array_equal(cutout.data, p[3].data[y:y+30, x:x+30])
    original_wcs, new_wcs = WCS(p[3].header), WCS(cutout.header)
    assert new_wcs.pixel_to_world(14, 17).separation(original_wcs.pixel_to_world(x+14, y+17)).arcsec < 1e-5
    assert not any(k.startswith("CD1_") for k in cutout.header)
    for key in ("BUNIT", "MAGZERO", "PPLANID", "PSWREL", "DATASETR", "EXPTIME", "GAIN", "FILTER"):
        assert cutout.header[key] == p[3].header[key]
    assert "CHECKSUM" not in cutout.header and "DATASUM" not in cutout.header
    assert cutout.header["MOSCLASS"] == "survey" and cutout.header["SRCS3"] == p[0]
    assert cutout.header["CUTBYTES"] < len(p[1]) / 10
    assert sum(size for _, _, size in fs.events) < len(p[1]) / 10
    assert client.query_sia.call_args.kwargs["collection"] == images.COLLECTIONS["wide"]


def test_wide_survey_is_distinct_from_repeat_and_deep():
    survey, repeat, deep = mosaic(), mosaic("repeat", version="repeat"), mosaic("deep")
    client, fs = setup(survey, repeat)
    found = irsa_tools.discover_mosaics(RA, DEC, client=client, filesystem=fs)
    assert list(found["kind"]) == ["survey", "repeat"]
    hdu = images.retrieve_cutout(RA, DEC, selection="repeat", cutout_size=3, client=client, filesystem=fs)
    assert hdu.header["MOSCLASS"] == "repeat" and hdu.header["SRCS3"] == repeat[0]
    client, fs = setup(deep)
    hdu = images.retrieve_cutout(RA, DEC, tiling="deep", cutout_size=3, client=client, filesystem=fs)
    assert hdu.header["MOSCLASS"] == "deep"


def test_overlap_is_ambiguous_until_tile_or_exact_product_is_selected():
    one, two = mosaic(), mosaic(tile="102158276", version="two")
    client, fs = setup(one, two)
    with pytest.raises(ValueError, match="Ambiguous"):
        images.retrieve_cutout(RA, DEC, client=client, filesystem=fs)
    chosen = images.retrieve_cutout(RA, DEC, tile_id="102158277", client=client, filesystem=fs)
    assert chosen.header["SRCS3"] == one[0]
    chosen = images.retrieve_cutout(RA, DEC, product=two[0], client=client, filesystem=fs)
    assert chosen.header["SRCS3"] == two[0]


@pytest.mark.parametrize("changes", [{"DATASETR": "Q1"}, {"FILTER": "NIR_J"},
                                      {"PPLANID": "unknown"}, {"MAGZERO": np.nan},
                                      {"CTYPE1A": "RA---TAN"}])
def test_unverified_headers_are_explicit_errors(changes):
    if changes.get("MAGZERO") is np.nan:
        # FITS rejects nonfinite header floats; a numeric-looking string is also invalid.
        changes = {"MAGZERO": "nan"}
    p = mosaic(**changes)
    client, fs = setup(p)
    found = images.discover_mosaics(RA, DEC, client=client, filesystem=fs)
    assert found["status"][0] == "error"
    with pytest.raises(ValueError, match="could not be verified"):
        images.retrieve_cutout(RA, DEC, client=client, filesystem=fs)


def test_missing_objects_are_listed_and_never_silently_replaced():
    p = mosaic()
    client, _ = setup(p)
    fs = MemoryFS({})
    found = images.discover_mosaics(RA, DEC, client=client, filesystem=fs)
    assert found["status"][0] == "error" and "KeyError" in found["error"][0]
    with pytest.raises(ValueError, match="could not be verified"):
        images.retrieve_cutout(RA, DEC, client=client, filesystem=fs)


def test_nisp_units_and_zero_point_remain_distinct():
    p = mosaic("deep", band="Y")
    client, fs = setup(p)
    hdu = images.retrieve_cutout(RA, DEC, tiling="deep", band="Y", client=client, filesystem=fs)
    assert hdu.header["BUNIT"] == "ELECTRON/s" and hdu.header["MAGZERO"] == 29.8


def test_byte_pixel_metadata_limits_and_edge_policy():
    p = mosaic()
    client, fs = setup(p)
    with pytest.raises(images.CutoutReadLimitError, match="byte budget"):
        images.retrieve_cutout(RA, DEC, max_bytes=100, client=client, filesystem=fs)
    with pytest.raises(ValueError, match="pixel budget"):
        images.retrieve_cutout(RA, DEC, cutout_size=1000, client=client, filesystem=fs)
    edge = WCS(p[3].header).pixel_to_world(2, 2)
    hdu = images.retrieve_cutout(edge.ra.deg, edge.dec.deg, cutout_size=3, client=client, filesystem=fs)
    assert hdu.header["CUTTRIM"] and min(hdu.data.shape) < 30
    with pytest.raises(ValueError):
        images.retrieve_cutout(edge.ra.deg, edge.dec.deg, mode="strict", cutout_size=3, client=client, filesystem=fs)
    client.query_sia.return_value = Table(rows=[p[2]] * images.MAX_METADATA)
    with pytest.raises(ValueError, match="metadata limit"):
        images.discover_mosaics(RA, DEC, client=client, filesystem=fs)


def test_empty_and_wrong_band_return_none_and_q1_is_unchanged():
    p = mosaic()
    client, fs = setup(p)
    assert images.retrieve_cutout(RA, DEC, band="H", client=client, filesystem=fs) is None
    with patch.object(irsa_tools.Irsa, "query_sia", return_value=Table(names=["dataproduct_subtype", "energy_bandpassname"])), patch.object(images, "retrieve_cutout") as new:
        assert irsa_tools.retrieve_cutout(RA, DEC, radius=2, cutout_size=3, band="VIS") is None
        new.assert_not_called()


def test_invalid_release_or_selection_fails_before_query():
    client = Mock()
    for kwargs in ({"release": "dr2"}, {"release": "dr1", "selection": "anything"},
                   {"release": "dr1", "tiling": "deep", "selection": "repeat"}):
        with pytest.raises(ValueError):
            irsa_tools.retrieve_cutout(RA, DEC, radius=2, cutout_size=3, band="VIS", client=client, **kwargs)
    client.query_sia.assert_not_called()


def test_live_observed_r2_dataset_is_preserved():
    p = mosaic(DATASETR="DR1_R2")
    client, fs = setup(p)
    hdu = images.retrieve_cutout(RA, DEC, client=client, filesystem=fs)
    assert hdu.header["DATASETR"] == "DR1_R2"


def test_discovery_accepts_search_radius_alias_and_empty_response():
    client = Mock()
    client.query_sia.return_value = Table()
    result = get_backend("irsa").discover_mosaics(RA, DEC, search_radius=1, client=client)
    assert len(result) == 0
    assert client.query_sia.call_args.kwargs["pos"][1] == 1*u.arcsec


def test_server_overflow_warning_cannot_become_a_partial_selection():
    import warnings
    from pyvo.dal.exceptions import DALOverflowWarning
    client = Mock()
    def overflow(**kwargs):
        warnings.warn("truncated", DALOverflowWarning)
        return Table()
    client.query_sia.side_effect = overflow
    with pytest.raises(DALOverflowWarning):
        images.discover_mosaics(RA, DEC, client=client)


@pytest.mark.parametrize("tiling,tile", [("wide", "101829777"), ("deep", "102158277")])
def test_mismatched_tiling_and_tile_id_fail_before_query(tiling, tile):
    client = Mock()
    with pytest.raises(ValueError, match="prefix disagrees"):
        images.discover_mosaics(RA, DEC, tiling=tiling, tile_id=tile, client=client)
    client.query_sia.assert_not_called()
