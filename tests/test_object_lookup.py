"""Offline checks for catalog ID lookups on both archives."""
import unittest
from unittest.mock import Mock, patch

import numpy as np
from astropy.table import Table

from euclid_tools import esa_tools, irsa_tools
from euclid_tools.backends import get_backend


class ObjectLookupTests(unittest.TestCase):
    def catalog_table(self, module):
        table = Table({"object_id": [123], "ra": [10.0], "dec": [20.0]})
        for column in module.COLUMNS_TO_REMOVE:
            table[column] = [0]
        for suffix in ("vis_psf", "y_templfit", "j_templfit", "h_templfit",
                       "vis_4fwhm_aper", "y_2fwhm_aper", "j_2fwhm_aper", "h_2fwhm_aper"):
            table[f"flux_{suffix}"] = [100.0]
            table[f"fluxerr_{suffix}"] = [1.0]
        if module is esa_tools:
            table.rename_column("ra", "right_ascension")
            table.rename_column("dec", "declination")
        return table

    def test_lookup_formats_results_and_queries_exact_id(self):
        for name, module, client, method, result_method, catalog in (
            ("esa", esa_tools, esa_tools.Euclid, "launch_job", "get_results", "catalogue.mer_catalogue"),
            ("irsa", irsa_tools, irsa_tools.Irsa, "query_tap", "to_table", irsa_tools.TABLE_MER),
        ):
            with self.subTest(archive=name):
                response = Mock()
                getattr(response, result_method).return_value = self.catalog_table(module)
                with patch.object(client, method, return_value=response) as query:
                    table = get_backend(name).retrieve_object_by_id("123")
                query.assert_called_once_with(f"SELECT * FROM {catalog} WHERE object_id = 123")
                self.assertEqual(table["object_id"][0], 123)
                self.assertIn("ra", table.colnames)
                self.assertIn("dec", table.colnames)
                self.assertNotIn("separation", table.colnames)
                self.assertTrue(all(col not in table.colnames for col in module.COLUMNS_TO_REMOVE))
                self.assertIn("VIS_4fwhm_aper_AB_mag", table.colnames)
                self.assertTrue(np.isfinite(table["VIS_AB_mag"][0]))
                getattr(response, result_method).return_value = Table()
                with patch.object(client, method, return_value=response):
                    self.assertIsNone(module.retrieve_object_by_id("123"))

    def test_invalid_ids_are_rejected_before_querying(self):
        for module in (esa_tools, irsa_tools):
            for object_id in ("", "1 OR 1=1", "1.0", "-1", "???", None):
                with self.subTest(module=module.__name__, object_id=object_id):
                    with self.assertRaises(ValueError):
                        module.retrieve_object_by_id(object_id)

