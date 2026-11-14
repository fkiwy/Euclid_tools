# Explicit ESA DR1 wide/deep catalog selection

ESA helpers now accept `catalogue="wide_survey"` or `catalogue="deep_survey"` for DR1. These are different catalog/source-ID sets, not a request to improve a spectrum selected from the other catalog. The same sky object can have separate wide, deep-stack and visit IDs.

| Call | Catalog |
| --- | --- |
| No release/selector | `q1.mer_catalogue` (unchanged) |
| `release="dr1"` | `dr1.mer_catalogue_wide_survey` (unchanged) |
| `release="dr1", catalogue="wide_survey"` | `dr1.mer_catalogue_wide_survey` |
| `release="dr1", catalogue="deep_survey"` | `dr1.mer_catalogue_deep_survey` |

The selector applies to `retrieve_objects`, `retrieve_object_by_id`, and `print_catalog_info`, and is forwarded by `get_backend("esa")`. Unknown selections, visits/mode collections, and a selector supplied with Q1 fail before a query. No cross-catalog or cross-release fallback occurs.

```python
from euclid_tools import get_backend
from euclid_tools.esa_tools import print_catalog_info

backend = get_backend('esa')
release, catalogue = 'dr1', 'deep_survey'
ra, dec = 59.7913643, -47.6826163
print_catalog_info(release=release, catalogue=catalogue)
objects = backend.retrieve_objects(ra, dec, search_radius=5,
                                   release=release, catalogue=catalogue)
if objects is not None:
    # Results are sorted by separation. Inspect ambiguity before selecting a
    # counterpart for science; the library never ranks spectra by quality.
    object_id = str(objects[0]['object_id'])
    same_object = backend.retrieve_object_by_id(object_id, release=release,
                                                catalogue=catalogue)
    spectrum = backend.retrieve_spectrum(object_id, release=release, grism='red')
```

Spectra take the returned source ID, not a catalog selector: an existing wide source ID is not converted to a deep ID. The catalog table and selected survey are recorded in returned table metadata (`catalogue_table`, `catalogue`). Raw IDs/flags are preserved, with the same profile AB calculations in µJy; no aperture/Vega mapping is inferred for DR1. Q1 retains its existing aperture/Vega outputs.

For the example script, set:

```python
release = 'dr1'
catalogue = 'deep_survey'
```

Set the desired RA/Dec in the coordinates block. Changing the survey does not change the DR1 coordinates. Metadata display, cone search and ID lookup use the same selection. The example sets image `tiling="deep", selection="deep"` for the deep survey, and `tiling="wide", selection="survey"` for wide. When using image methods directly, pass these image arguments explicitly: catalog queries do not mutate backend state. Ambiguous mosaics still require a tile/exact product.

## Search integration

Use the companion `euclid_dr1/archive_dr1_esa_deep_products.json` for ESA deep catalog, spectra and deep mosaics. It explicitly sets `catalogue="deep_survey"` and the matching native table. The existing `archive_dr1_esa_products.json` remains wide by default; the Q1/IRSA configs retain their behavior. A deep table with an omitted/wide selector is rejected instead of querying the wide table silently.

Set `TEMPLATE_COLOR_SCATTER` to the same adopted value (mag) you already use; gathering still requires this explicit scientific assumption. No new value is supplied by this change.

```sh
python -m pip install -e /path/Euclid_tools
python -m euclid_dr1.gather_data --archive-config euclid_dr1/archive_dr1_esa_deep_products.json --input /tmp/targets.csv --output /tmp/deep-gathered.csv --template-color-scatter "$TEMPLATE_COLOR_SCATTER"
python -m euclid_dr1.save_spectra --archive-config euclid_dr1/archive_dr1_esa_deep_products.json --input /tmp/deep-gathered.csv --output /tmp/deep-spectra.csv --spectra-dir /tmp/deep-spectra
python -m euclid_dr1.plot_images --archive-config euclid_dr1/archive_dr1_esa_deep_products.json --input /tmp/deep-gathered.csv --output /tmp/deep-images.csv --cutouts-dir /tmp/deep-cutouts --plots-dir /tmp/deep-plots
```

Use coordinate-only inputs (or deep IDs) for deep gathering; supplying a wide ID does not authorize substituting its deep counterpart. Keep separate outputs when comparing surveys. Existing wrong-ID/backend/config cache checks remain active.

## Verification and scientific limits

On 2026-11-14, [ESA TAP](https://eas.esac.esa.int/tap-server/tap/sync) `TAP_SCHEMA.columns` showed identical names, datatypes, units and descriptions in the wide/deep catalogs for 19 required fields: object ID; RA/Dec; eight VIS/Y/J/H profile flux/error fields; detection/band flags; point-like flag/probability and Gaia match quality. [Exact metadata snapshot](esa_dr1_deep_catalog_metadata.json) records them. Fluxes are µJy and coordinates degrees. Matching metadata does not establish identical calibration, completeness or flag semantics, and none are assumed.

[Live coordinate → catalog ID → spectrum validation](esa_dr1_catalogue_selection_validation.json) checked both surveys for both targets through the public ESA helpers. An injected requests-backed native TAP cone client bypasses this environment's Astroquery HTTP/proxy limitation; the public spectrum method uses its normal native ESA transport. The checks independently verify AB photometry and source identity, with no filters or normalization:

| Position (degrees) | Wide source ID | Deep source ID |
| --- | --- | --- |
| 58.1332495, -49.1830038 | -581332416491830015 | -581332463491830110 |
| 59.7913643, -47.6826163 | -597913631476826029 | -597913580476826189 |

All four checks returned the expected IDs and 531-bin red spectra. Deep stacks have 7 zero-contribution bins each, versus 88/489 in the wide products in these samples. These counts describe delivered products; they are not new science cuts or proof of general deep/wide quality. An initial ESA spectrum request returned HTTP 503; a subsequent full run succeeded. `scripts/validate_esa_catalogue_selection.py` records catalog/product failures separately and never substitutes another product.

Other exact official sources: [ESA DR1 data access](https://www.cosmos.esa.int/web/euclid/dr1-data), [SIR processing overview](https://euclid.esac.esa.int/dr/dr1/expsup/Ch11.S1.html), and the [ESA product interface sources/limits](ESA_DR1_PRODUCTS.md). Individual visits and wide/deep mode catalog helpers are not enabled here. Production candidate-selection equivalence and detailed causes of the degraded wide spectra remain unresolved. No flux scaling, masks, fitting thresholds or scientific selection criteria changed; no bulk survey products were downloaded.


The actual companion `Archive` also retrieved the deep catalog source, saved its spectrum, and saved a 31 x 31 deep VIS cutout at `(59.7913643, -47.6826163)` using the new config. [Workflow identity record](esa_dr1_deep_workflow_validation.json) captures the native table, file/HDU, deep mosaic processing plan and hashes.

Deep queries explicitly request CSV, which Astroquery's current ESA interface also defaults to. The live full-column VOTable response has malformed `arraysize="100x*"` metadata on `calblock_id_list`; CSV avoids that parser error while retaining all native columns and exact integer IDs. The public deep ID helper was additionally checked against a full native CSV response (479 input columns), not only the 19-column validation cone. No scientific field is dropped to work around the archive issue.

Final combined validation: **224 passed, 24 subtests passed** across tools and search, including Q1 defaults, unchanged DR1-wide defaults, cone/ID/metadata selection, invalid-selector rejection, the example's stable coordinates and deep-image routing, and actual companion catalog-ID-to-spectrum integration. Both diffs pass whitespace checks.
