"""IRSA catalog release forwarding and the example's AB-only display path."""
from contextlib import redirect_stdout
from io import StringIO
from unittest.mock import Mock, patch

import numpy as np
import pytest
from astropy import units as u
from astropy.table import Table

from euclid_tools import irsa_tools
from euclid_tools.backends import get_backend
from euclid_tools.shared import print_results


def raw_dr1():
    table = Table({"object_id": [-123, -456], "ra": [80.0001, 80.], "dec": [-51., -51.],
                   "tileid": [102039700, 102039701], "flag_y": [7, 0]})
    for suffix in ("vis_psf", "y_templfit", "j_templfit", "h_templfit"):
        table[f"flux_{suffix}"] = [36.31, 3.631] * u.uJy
        table[f"fluxerr_{suffix}"] = [3.631, .3631] * u.uJy
    return table


def test_reported_backend_call_selects_dr1_and_prints_ab_only():
    client = Mock()
    client.query_region.return_value = raw_dr1()
    result = get_backend("irsa").retrieve_objects(80, -51, search_radius=5, release="dr1", client=client)
    call = client.query_region.call_args.kwargs
    assert call["catalog"] == "euclid_dr1_mer_catalogue_wide_survey"
    assert call["radius"] == 5 * u.arcsec
    assert list(result["object_id"]) == [-456, -123]
    assert result["separation"][0] == 0
    assert result.meta["release"] == "DR1"
    assert result["Y_AB_mag"].tolist() == [22.5, 20.]
    np.testing.assert_allclose(result["Y_AB_err"], [0.109, 0.109])
    assert "tileid" in result.colnames and 7 in result["flag_y"]
    assert not any("Vega" in col or "aper_AB" in col for col in result.colnames)
    output = StringIO()
    with redirect_stdout(output):
        print_results(result)
    assert "Y_AB_mag" in output.getvalue()
    assert "Vega" not in output.getvalue()


def test_id_lookup_forwarding_preserves_negative_id_and_all_rows():
    client = Mock()
    rows = raw_dr1()
    rows["object_id"] = [-123, -123]
    client.query_tap.return_value.to_table.return_value = rows
    result = get_backend("irsa").retrieve_object_by_id("-123", release="DR1", client=client)
    client.query_tap.assert_called_once_with(
        "SELECT * FROM euclid_dr1_mer_catalogue_wide_survey WHERE object_id = -123")
    assert result["object_id"].tolist() == [-123, -123]
    assert "separation" not in result.colnames
    assert result["VIS_AB_mag"][0] == 20.


def test_empty_results_return_none():
    client = Mock()
    client.query_region.return_value = Table()
    client.query_tap.return_value.to_table.return_value = Table()
    assert irsa_tools.retrieve_objects(80, -51, radius=5, release="dr1", client=client) is None
    assert irsa_tools.retrieve_object_by_id("-123", release="dr1", client=client) is None


@pytest.mark.parametrize("release", ["dr2", "", None, 1])
def test_invalid_release_fails_before_network(release):
    client = Mock()
    with pytest.raises(ValueError):
        irsa_tools.retrieve_objects(80, -51, radius=5, release=release, client=client)
    with pytest.raises(ValueError):
        irsa_tools.retrieve_object_by_id("-123", release=release, client=client)
    client.query_region.assert_not_called()
    client.query_tap.assert_not_called()


def test_q1_cone_default_and_explicit_release_keep_aperture_and_vega():
    client = Mock()
    table = raw_dr1()
    for column in irsa_tools.COLUMNS_TO_REMOVE:
        if column not in table.colnames:
            table[column] = [0, 0]
    for suffix in ("vis_4fwhm_aper", "y_2fwhm_aper", "j_2fwhm_aper", "h_2fwhm_aper"):
        table[f"flux_{suffix}"] = [36.31, 3.631]
        table[f"fluxerr_{suffix}"] = [3.631, .3631]
    for options in ({}, {"release": "Q1"}):
        client.query_region.return_value = table.copy()
        result = irsa_tools.retrieve_objects(80, -51, radius=5, client=client, **options)
        assert client.query_region.call_args.kwargs["catalog"] == irsa_tools.TABLE_MER
        assert "Y_2fwhm_aper_AB_mag" in result.colnames
        assert "Y_Vega_mag" in result.colnames


def test_catalog_metadata_uses_selected_table():
    service = Mock()
    service.tables = {irsa_tools.TABLE_MER_DR1: Mock(columns=[])}
    with patch.object(irsa_tools.vo.dal, "TAPService", return_value=service):
        irsa_tools.print_catalog_info(release="dr1")


@pytest.mark.parametrize("catalogue", ["wide_survey", "deep_survey", "wide_mode", "deep_mode"])
def test_explicit_dr1_catalogue_is_forwarded_to_cone_id_and_metadata(catalogue):
    client = Mock()
    client.query_region.return_value = raw_dr1()
    client.query_tap.return_value.to_table.return_value = raw_dr1()
    backend = get_backend("irsa")
    result = backend.retrieve_objects(80, -51, radius=5, release="dr1", catalogue=catalogue, client=client)
    expected = f"euclid_dr1_mer_catalogue_{catalogue}"
    assert client.query_region.call_args.kwargs["catalog"] == expected
    assert result.meta["catalogue_table"] == expected
    backend.retrieve_object_by_id("-123", release="dr1", catalogue=catalogue, client=client)
    assert f"FROM {expected}" in client.query_tap.call_args.args[0]
    service = Mock()
    service.tables = {expected: Mock(columns=[])}
    with patch.object(irsa_tools.vo.dal, "TAPService", return_value=service):
        irsa_tools.print_catalog_info(release="dr1", catalogue=catalogue)


@pytest.mark.parametrize("release,catalogue", [("q1", "deep_survey"), ("dr1", ""),
                                              ("dr1", "dr1.sources"), ("dr1", 1)])
def test_catalogue_rejection_precedes_network(release, catalogue):
    client = Mock()
    with pytest.raises(ValueError):
        irsa_tools.retrieve_objects(80, -51, radius=5, release=release, catalogue=catalogue, client=client)
    client.query_region.assert_not_called()


def test_no_match_message_identifies_release_and_catalogue():
    from euclid_tools.shared import print_message
    output = StringIO()
    with redirect_stdout(output):
        print_message(273, 68, 5, release="dr1", catalogue=irsa_tools.TABLE_MER_DR1)
    text = output.getvalue()
    assert "DR1" in text and "Q1" not in text
    assert irsa_tools.TABLE_MER_DR1 in text
    output = StringIO()
    with redirect_stdout(output):
        print_message(273, 68, 5)
    assert "Euclid Q1 MER catalog" in output.getvalue()
