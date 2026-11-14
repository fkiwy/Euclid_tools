"""Release isolation and arithmetic checks using a three-row public DR1 fixture."""
from pathlib import Path
from unittest.mock import Mock, patch

import numpy as np
import pytest
from astropy.table import Table

from euclid_tools import esa_tools
from euclid_tools.backends import get_backend
from euclid_tools.releases import catalog_table


@pytest.fixture
def sample():
    return Table.read(Path(__file__).parent / "data/dr1_sample.csv", format="ascii.csv")


@pytest.mark.parametrize("release", [None, "", "dr2", "dr1;drop table x", 1])
def test_invalid_release_never_queries(release):
    with patch.object(esa_tools.Euclid, "launch_job") as query:
        with pytest.raises(ValueError):
            esa_tools.retrieve_object_by_id("1", release=release)
        query.assert_not_called()


def test_release_mapping_preserves_q1_default():
    assert catalog_table() == "q1.mer_catalogue"
    assert catalog_table("Q1") == catalog_table()
    assert catalog_table("DR1") == "dr1.mer_catalogue_wide_survey"


def test_dr1_id_lookup_and_backend_forwarding(sample):
    sample = sample[:1]
    object_id = str(sample["object_id"][0])
    client = Mock()
    client.launch_job.return_value.get_results.return_value = sample
    table = get_backend("esa").retrieve_object_by_id(object_id, release="DR1", client=client)
    client.launch_job.assert_called_once_with(
        f"SELECT * FROM dr1.mer_catalogue_wide_survey WHERE object_id = {object_id}")
    assert table.meta["release"] == "DR1"
    assert "ra" in table.colnames and "dec" in table.colnames
    assert "Y_2fwhm_aper_AB_mag" not in table.colnames
    assert not any("Vega" in name for name in table.colnames)
    assert str(table["object_id"][0]) == object_id
    assert "flag_y" in table.colnames
    np.testing.assert_allclose(table["Y_AB_mag"],
                               np.round(-2.5*np.log10(np.array(table["flux_y_templfit"])*1e-6/3631), 3))
    np.testing.assert_allclose(table["Y_AB_err"],
                               np.round(2.5/np.log(10)*table["fluxerr_y_templfit"]/table["flux_y_templfit"], 3))


def test_dr1_cone_and_empty_results(sample):
    sample["dist"] = [0.001, 0.0001, 0.002]
    client = Mock()
    client.cone_search.return_value.get_results.return_value = sample
    table = esa_tools.retrieve_objects(80, -51, search_radius=5, limit=3, release="dr1", client=client)
    call = client.cone_search.call_args.kwargs
    assert call["table_name"] == catalog_table("dr1")
    assert call["ra_column_name"] == "right_ascension"
    assert client.ROW_LIMIT == 3
    np.testing.assert_allclose(table["separation"], [0.36, 3.6, 7.2])
    client.cone_search.return_value.get_results.return_value = Table()
    assert esa_tools.retrieve_objects(80, -51, radius=5, release="dr1", client=client) is None
    client.launch_job.return_value.get_results.return_value = Table()
    assert esa_tools.retrieve_object_by_id("123", release="dr1", client=client) is None


def test_dr1_product_routes_never_call_q1_endpoints():
    from euclid_tools import esa_products
    with patch.object(esa_tools.Euclid, "get_spectrum") as spectrum, patch.object(esa_tools.Euclid, "launch_job") as job:
        with patch.object(esa_products, "retrieve_spectrum", return_value=None) as dr1_spectrum:
            assert get_backend("esa").retrieve_spectrum("123", release="dr1") is None
            dr1_spectrum.assert_called_once()
        with patch.object(esa_products, "retrieve_cutout", return_value=None) as dr1_cutout:
            assert get_backend("esa").retrieve_cutout(80, -51, radius=2, cutout_size=10, band="Y", release="dr1") is None
            dr1_cutout.assert_called_once()
        spectrum.assert_not_called()
        job.assert_not_called()


def test_failed_archive_job_is_not_an_empty_catalog():
    client = Mock()
    client.launch_job.return_value = None
    with pytest.raises(RuntimeError, match="query failed"):
        esa_tools.retrieve_object_by_id("-123", release="dr1", client=client)
    client.cone_search.return_value = None
    with pytest.raises(RuntimeError, match="query failed"):
        esa_tools.retrieve_objects(80, -51, radius=2, release="dr1", client=client)


def test_dr1_catalog_metadata_selects_release():
    client = Mock()
    client.columns = []
    with patch.object(esa_tools.Euclid, "load_table", return_value=client) as load:
        esa_tools.print_catalog_info(release="dr1")
        load.assert_called_once_with("dr1.mer_catalogue_wide_survey")
