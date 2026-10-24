import unittest
from contextlib import redirect_stdout
from io import StringIO

import numpy as np
from astropy import units as u
from astropy.table import QTable

from euclid_tools.backends import get_backend, resolve_radius
from euclid_tools.shared import normalize_spectrum_table, print_results, remove_bad_values


class ApiContractTests(unittest.TestCase):
    def test_print_results_only_includes_separation_when_present(self):
        columns = [
            "object_id", "ra", "dec",
            "VIS_AB_mag", "VIS_AB_err", "VIS_Vega_mag", "VIS_Vega_err",
            "VIS_4fwhm_aper_AB_mag", "VIS_4fwhm_aper_AB_err",
            "VIS_4fwhm_aper_Vega_mag", "VIS_4fwhm_aper_Vega_err",
            "Y_AB_mag", "Y_AB_err", "Y_Vega_mag", "Y_Vega_err",
            "Y_2fwhm_aper_AB_mag", "Y_2fwhm_aper_AB_err",
            "Y_2fwhm_aper_Vega_mag", "Y_2fwhm_aper_Vega_err",
            "J_AB_mag", "J_AB_err", "J_Vega_mag", "J_Vega_err",
            "J_2fwhm_aper_AB_mag", "J_2fwhm_aper_AB_err",
            "J_2fwhm_aper_Vega_mag", "J_2fwhm_aper_Vega_err",
            "H_AB_mag", "H_AB_err", "H_Vega_mag", "H_Vega_err",
            "H_2fwhm_aper_AB_mag", "H_2fwhm_aper_AB_err",
            "H_2fwhm_aper_Vega_mag", "H_2fwhm_aper_Vega_err",
        ]
        table = QTable([[0] for _ in columns], names=columns)

        output = StringIO()
        with redirect_stdout(output):
            print_results(table)
        self.assertNotIn("separation", output.getvalue())

        table["separation"] = [1.0]
        output = StringIO()
        with redirect_stdout(output):
            print_results(table)
        self.assertIn("separation", output.getvalue())

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
                np.array([1.0, 2.0, 3.0]) * u.micron,
                np.array([10.0, 20.0, 30.0]) * u.Jy,
                np.array([1.0, 2.0, 3.0]) * u.Jy,
                np.array([0, 64, 0], dtype=int),
            ],
            names=("WAVELENGTH", "FLUX", "ERROR", "MASK"),
        )

        normalized = normalize_spectrum_table(table)
        self.assertEqual(list(normalized.colnames), ["WAVELENGTH", "FLUX", "ERROR", "MASK"])
        self.assertEqual(normalized["MASK"][1], 64)
        self.assertEqual(float(normalized["FLUX"][1].value), 20.0)
        self.assertEqual(normalized["WAVELENGTH"].unit, u.micron)

    def test_remove_bad_values_filters_bitmask(self):
        table = QTable(
            [
                np.array([1.0, 2.0, 3.0]) * u.micron,
                np.array([10.0, 20.0, 30.0]) * u.Jy,
                np.array([1.0, 2.0, 3.0]) * u.Jy,
                np.array([0, 64, 0], dtype=int),
            ],
            names=("WAVELENGTH", "FLUX", "ERROR", "MASK"),
        )

        cleaned = remove_bad_values(table, lambda mask: (mask % 2 == 1) | (mask >= 64))
        self.assertEqual(list(cleaned.colnames), ["WAVELENGTH", "FLUX", "ERROR"])
        # remove_bad_values returns a sliced table containing only good rows
        self.assertEqual(len(cleaned), 2)
        self.assertEqual(float(cleaned["FLUX"][0].value), 10.0)
        self.assertEqual(float(cleaned["FLUX"][1].value), 30.0)
        self.assertEqual(cleaned["WAVELENGTH"].unit, u.micron)

    def test_normalize_spectrum_table_handles_integer_mask_columns(self):
        masked_mask = np.ma.array(np.array([0, 64, 0], dtype=int), mask=[False, True, False])
        table = QTable(
            [
                np.array([1.0, 2.0, 3.0]) * u.micron,
                np.array([10.0, 20.0, 30.0]) * u.Jy,
                np.array([1.0, 2.0, 3.0]) * u.Jy,
                masked_mask,
            ],
            names=("WAVELENGTH", "FLUX", "ERROR", "MASK"),
        )

        normalized = normalize_spectrum_table(table)
        self.assertEqual(normalized["MASK"][1], 0)
        self.assertEqual(normalized["WAVELENGTH"].unit, u.micron)


if __name__ == "__main__":
    unittest.main()
