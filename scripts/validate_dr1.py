"""Bounded public TAP validation. Run from a checkout; no bulk data or login.

Usage: python scripts/validate_dr1.py --output /tmp/dr1-validation
The sample is archive order (not representative); no scientific cuts are applied.
"""
import argparse
import csv
import io
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import requests
from astropy.table import Table

from euclid_tools.esa_tools import _prepare_catalog_table

ENDPOINT = "https://eas.esac.esa.int/tap-server/tap/sync"
COLUMNS = ["object_id", "right_ascension", "declination",
           "flux_vis_psf", "fluxerr_vis_psf", "flux_y_templfit", "fluxerr_y_templfit",
           "flux_j_templfit", "fluxerr_j_templfit", "flux_h_templfit", "fluxerr_h_templfit",
           "det_quality_flag", "flag_vis", "flag_y", "flag_j", "flag_h",
           "point_like_flag", "point_like_prob", "gaia_match_quality"]


def query(adql, maxrec=3):
    response = requests.get(ENDPOINT, params={"REQUEST": "doQuery", "LANG": "ADQL",
                            "FORMAT": "csv", "MAXREC": maxrec, "QUERY": adql}, timeout=30)
    response.raise_for_status()
    if "QUERY_STATUS" in response.text or response.text.lstrip().startswith("<"):
        raise RuntimeError(response.text[:1000])
    return response.text


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output", type=Path, required=True)
    args = cli.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    adql = f"SELECT TOP 3 {','.join(COLUMNS)} FROM dr1.mer_catalogue_wide_survey"
    raw = query(adql)
    table = Table.read(raw, format="ascii.csv")
    assert 0 < len(table) <= 3
    converted = _prepare_catalog_table(table.copy(), release="dr1")
    checks = {}
    for band, suffix in (("VIS", "psf"), ("Y", "templfit"), ("J", "templfit"), ("H", "templfit")):
        flux, err = table[f"flux_{band.lower()}_{suffix}"], table[f"fluxerr_{band.lower()}_{suffix}"]
        assert np.all(flux > 0) and np.all(err > 0)
        expected = -2.5 * np.log10(np.asarray(flux) * 1e-6 / 3631.)
        expected_err = 2.5 / np.log(10) * np.asarray(err) / flux
        np.testing.assert_allclose(converted[f"{band}_AB_mag"], expected, atol=0.000501, rtol=0)
        np.testing.assert_allclose(converted[f"{band}_AB_err"], expected_err, atol=0.000501, rtol=0)
        checks[band] = {"AB_mag": list(map(float, expected)), "AB_err": list(map(float, expected_err))}
    metadata_query = ("SELECT TOP 100 table_name,column_name,unit,description FROM TAP_SCHEMA.columns "
                      "WHERE table_name IN ('q1.mer_catalogue','dr1.mer_catalogue_wide_survey') "
                      "AND column_name IN (" + ",".join(f"'{c}'" for c in COLUMNS) + ")")
    metadata = query(metadata_query, 100)
    rows = list(csv.DictReader(io.StringIO(metadata)))
    for column in COLUMNS:
        pair = [r for r in rows if r["column_name"] == column]
        assert len(pair) == 2, column
        assert pair[0]["unit"] == pair[1]["unit"], column
        assert pair[0]["description"] == pair[1]["description"], column
    (args.output / "dr1_sample.csv").write_text(raw)
    (args.output / "catalog_metadata.csv").write_text(metadata)
    report = {"utc": datetime.now(timezone.utc).isoformat(), "endpoint": ENDPOINT,
              "query": adql, "metadata_query": metadata_query, "rows": len(table),
              "metadata_columns_equal": COLUMNS, "photometry": checks,
              "limitations": "Schema and arithmetic checks only; no calibration or selection validation"}
    (args.output / "validation.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"Validated {len(table)} DR1 rows, {len(COLUMNS)} metadata pairs and four AB bands")


if __name__ == "__main__":
    main()
