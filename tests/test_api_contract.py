import unittest

import numpy as np
from astropy.table import QTable

from euclid_tools.backends import get_backend, resolve_radius
from euclid_tools.shared import mask_bad_values, normalize_spectrum_table


class ApiContractTests(unittest.TestCase):
    def test_resolve_radius_accepts_aliases(self):
        self.assertEqual(resolve_radius(radius=5.0), 5.0)
        self.assertEqual(resolve_radius(search_radius=5.0), 5.0)
        self.assertEqual(resolve_radius(radius=5.0, search_radius=5.0), 5.0)

        with self.assertRaises(TypeError):
            resolve_radius()

        with self.assertRaises(ValueError):
            resolve_radius(radius=5.0, search_radius=7.0)

    def test_backend_factory(self):
        self.assertEqual(get_backend("esa").name, "esa")
        self.assertEqual(get_backend("irsa").name, "irsa")

        with self.assertRaises(ValueError):
            get_backend("unknown")

    def test_normalize_spectrum_table_uses_canonical_columns(self):
        table = QTable(
            [
                np.array([1.0, 2.0, 3.0]) * 1e-6,
                np.array([10.0, 20.0, 30.0]),
                np.array([1.0, 2.0, 3.0]),
                np.array([0, 64, 0], dtype=int),
            ],
            names=("WAVELENGTH", "FLUX", "ERROR", "MASK"),
        )

        normalized = normalize_spectrum_table(table, include_mask=True)
        self.assertEqual(list(normalized.colnames), ["WAVELENGTH", "FLUX", "ERROR", "MASK"])
        self.assertEqual(normalized["MASK"][1], 64)
        self.assertEqual(float(normalized["FLUX"][1]), 20.0)

    def test_mask_bad_values_filters_bitmask(self):
        table = QTable(
            [
                np.array([1.0, 2.0, 3.0]) * 1e-6,
                np.array([10.0, 20.0, 30.0]),
                np.array([1.0, 2.0, 3.0]),
                np.array([0, 64, 0], dtype=int),
            ],
            names=("WAVELENGTH", "FLUX", "ERROR", "MASK"),
        )

        cleaned = mask_bad_values(table, lambda mask: (mask % 2 == 1) | (mask >= 64))
        self.assertEqual(list(cleaned.colnames), ["WAVELENGTH", "FLUX", "ERROR"])
        self.assertTrue(np.isnan(cleaned["FLUX"][1]))


if __name__ == "__main__":
    unittest.main()
