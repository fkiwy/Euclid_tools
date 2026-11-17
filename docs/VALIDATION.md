# Validation records and reproduction

The automated suite runs with this public checkout alone:

```sh
python -m pip install -e ".[dev]"
python -m pytest tests -q
```

The documentation refresh on 17 November 2026 passed **271 tests and 24 subtests**
using this repository alone.

Tests cover release selection, catalog mappings/photometry, native spectral
units, masks/variance, exact product choices, bounded transports, image WCS and
calibration preservation, provenance/FITS round trips, Q1 compatibility, import
paths and file-based plotting. They use offline fixtures and mocked transports;
a passing suite does not establish current archive access or survey calibration.

## Recorded archive checks

The records below are dated response snapshots. Their timestamps, sample sizes,
errors and historical limitations describe those checks, not guaranteed current
availability or a complete feature matrix. See [current support](DR1_FOUNDATION.md)
and the archive guides for supported functionality.

| Record | Scope |
|---|---|
| [dr1_validation.json](dr1_validation.json) | Three ESA catalog rows, 19 metadata pairs and independent AB arithmetic |
| [ESA catalog selection](esa_dr1_catalogue_selection_validation.json) | Wide/deep cone, ID and red-spectrum checks for two positions |
| [ESA deep metadata](esa_dr1_deep_catalog_metadata.json) | Native deep-survey catalog field metadata |
| [ESA products](esa_dr1_products_validation.json) | Bounded red/blue spectra, survey/repeat/deep cutouts and preserved headers |
| [ESA associations](esa_dr1_spectrum_associations_validation.json) | Three source-patch spectra and one deliberately blocked shared-patch parent |
| [ESA overlap](esa_dr1_mosaic_overlap_validation.json) | Four-band cutouts using the selected catalog tile |
| [IRSA catalog](irsa_dr1_catalog_validation.json) | Targeted cone/ID lookup and AB arithmetic |
| [IRSA coverage](irsa_dr1_catalog_coverage_validation.json) | Wide/deep coverage and backend cone checks; ID-query timeout recorded |
| [IRSA spectrum](irsa_dr1_spectrum_validation.json) | One public red spectrum and reference flux/error bin |
| [IRSA images](irsa_dr1_images_validation.json) | Small SIA/S3 cutouts, ranged GETs, WCS/calibration and byte counts |

Catalog metadata CSVs and small spectrum/header fixtures remain alongside the
records and in tests/data. Public records contain library/adapter evidence only.
No full parent spectrum bundles or survey mosaics are needed for the tests.

## Optional small live checks

These commands contact the named archives and transfer bounded samples. They
can fail or be slow when a service is unavailable; no cross-release fallback is
performed. Run only the checks relevant to your work:

```sh
python scripts/validate_dr1.py --output validation_output/esa-catalog
python scripts/validate_irsa_dr1_spectrum.py --output validation_output/irsa-spectrum
python scripts/validate_irsa_dr1_images.py --output validation_output/irsa-images
python scripts/validate_esa_dr1_products.py --output validation_output/esa-products.json
python scripts/validate_esa_catalogue_selection.py --output validation_output/esa-catalogues.json
python scripts/validate_esa_spectrum_associations.py --output validation_output/esa-associations.json
```

The catalog script records exact ADQL, endpoint, UTC metadata and independent
AB arithmetic. Product scripts preserve IDs, routes, units, hashes and transfer
limits in their reports. The guides link the official sources used to interpret
those products. Metadata agreement and transport correctness do not validate
scientific quality thresholds, Vega coefficients or release equivalence.
