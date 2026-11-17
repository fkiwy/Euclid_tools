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

DR1 optionally accepts an explicitly supplied `vega_calibration` convention;
see [reference fluxes, required sources and scientific limitations](DR1_VEGA_CALIBRATION.md).
The same option works for ESA/IRSA cone searches and ID lookups. Default AB
photometry and all native measurements remain unchanged.

For the example script, keep its current defaults:

```python
release = 'dr1'
dr1_catalogue = 'deep_survey'
```

Set the desired RA/Dec in the coordinates block. Changing the survey does not change the DR1 coordinates. Metadata display, cone search and ID lookup use the same selection. The example sets image `tiling="deep", selection="deep"` for the deep survey, and `tiling="wide", selection="survey"` for wide. When using image methods directly, pass these image arguments explicitly: catalog queries do not mutate backend state. Ambiguous mosaics still require a tile/exact product.

## Verification and scientific limits

On 2026-11-14, [ESA TAP](https://eas.esac.esa.int/tap-server/tap/sync) `TAP_SCHEMA.columns` showed identical names, datatypes, units and descriptions in the wide/deep catalogs for 19 required fields: object ID; RA/Dec; eight VIS/Y/J/H profile flux/error fields; detection/band flags; point-like flag/probability and Gaia match quality. [Exact metadata snapshot](esa_dr1_deep_catalog_metadata.json) records them. Fluxes are µJy and coordinates degrees. Matching metadata does not establish identical calibration, completeness or flag semantics, and none are assumed.

[Live coordinate → catalog ID → spectrum validation](esa_dr1_catalogue_selection_validation.json) checked both surveys for both targets through the public ESA helpers. An injected requests-backed native TAP cone client bypasses this environment's Astroquery HTTP/proxy limitation; the public spectrum method uses its normal native ESA transport. The checks independently verify AB photometry and source identity, with no filters or normalization:

| Position (degrees) | Wide source ID | Deep source ID |
| --- | --- | --- |
| 58.1332495, -49.1830038 | -581332416491830015 | -581332463491830110 |
| 59.7913643, -47.6826163 | -597913631476826029 | -597913580476826189 |

All four checks returned the expected IDs and 531-bin red spectra. Deep stacks have 7 zero-contribution bins each, versus 88/489 in the wide products in these samples. These counts describe delivered products; they are not new science cuts or proof of general deep/wide quality. An initial ESA spectrum request returned HTTP 503; a subsequent full run succeeded. `scripts/validate_esa_catalogue_selection.py` records catalog/product failures separately and never substitutes another product.

Other exact official sources: [ESA DR1 data access](https://www.cosmos.esa.int/web/euclid/dr1-data), [SIR processing overview](https://euclid.esac.esa.int/dr/dr1/expsup/Ch11.S1.html), and the [ESA product interface sources/limits](ESA_DR1_PRODUCTS.md). Individual visits and wide/deep mode catalog helpers are not enabled here. Production candidate-selection equivalence and detailed causes of the degraded wide spectra remain unresolved. No flux scaling, masks, fitting thresholds or scientific selection criteria changed; no bulk survey products were downloaded.


Deep queries explicitly request CSV, which Astroquery's current ESA interface also defaults to. The live full-column VOTable response has malformed `arraysize="100x*"` metadata on `calblock_id_list`; CSV avoids that parser error while retaining all native columns and exact integer IDs. The public deep ID helper was additionally checked against a full native CSV response (479 input columns), not only the 19-column validation cone. No scientific field is dropped to work around the archive issue.


## Overlapping mosaics at 58.1332495, -49.1830038

Both deep tiles `101280483` and `101280484` cover this position in VIS/Y/J/H.
The selected deep catalog object `-581332463491830110` explicitly records
`tile_index=101280484`. The example now uses that catalog tile for DR1 images
when `image_tile_id` is None. Set `image_tile_id="101280484"` to select that tile explicitly; a user-supplied override always takes precedence. The tile
comes from the selected catalog row, not from a spectrum filename, a nearest
mosaic calculation or a quality ranking.

Direct library calls still require explicit tile/product selection when several
mosaics match. Errors now list candidate tiles, releases and exact paths; use
`backend.discover_mosaics(..., tiling="deep", selection="deep")` to inspect all
candidates. A tile with multiple product versions remains ambiguous and requires
an exact product path. Missing/masked catalog tile metadata does not authorize
arbitrary selection. Pass `tile_id` or `product` explicitly to the library methods.

On 2026-11-17 the catalog tile and the two overlapping candidates were checked
against live ESA metadata. [Overlap validation](esa_dr1_mosaic_overlap_validation.json)
records successful, small four-band cutouts from catalog tile 101280484, with
native processing/calibration/WCS checks and FITS re-serialization. Q1 routing,
scientific cuts and the library's ambiguity guard are unchanged.

See [validation](VALIDATION.md) for current local checks and dated archive evidence.
