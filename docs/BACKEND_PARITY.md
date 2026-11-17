# ESA / IRSA behavior audit

Audited 2026-11-17 using local fixtures and both real decoders. No survey data
were downloaded and archive availability/speed were not re-measured.

## Common public contract

Both backend objects and `esa_tools` / `irsa_tools` expose these methods:

| Method | Shared behavior | Archive-dependent details |
|---|---|---|
| `retrieve_objects` | Q1 default, explicit DR1; radius aliases; sorted separation in arcseconds rounded to 0.001; None for empty results | Native catalogs/raw columns, clients and default row limits |
| `retrieve_object_by_id` | Exact signed-int64 ID; matching photometry; no separation column; None for empty results | Native tables/query clients |
| `retrieve_spectrum` | Q1 default; explicit DR1; canonical QTable or None; ambiguity raises; optional FITS saving | Formats, grism evidence and enforceable selectors |
| `retrieve_spectra` | DR1 default, explicit release supported; Q1 rejected; separate products; empty list for no matches; no ranking/combining | IRSA grism filtering requires reading wavelengths |
| `discover_spectra` | DR1-only metadata Table with source_path/hdu_index, grism, route fields, retrieval_supported; no spectrum downloads | IRSA grism is unknown; grism-filtered metadata discovery is unsupported |
| `discover_mosaics` | DR1-only; 2-arcsecond default radius; both aliases; survey/deep defaults; explicit repeat selection; common product/tile/band/class/status fields | Native columns; IRSA header-verification error rows |
| `retrieve_cutout` | Q1 default; explicit radius/size/band; image HDU or None; DR1 product/tile ambiguity raises | Paths, WCS representation, rounding, formats and transport limits |

The newer DR1-only methods do not alter the singular helpers' Q1 defaults.
Release names are case-insensitive. Radius aliases must agree; radii/sizes must
be positive finite scalars. IDs never pass through floating point. Both cone
helpers accept `limit` through client ROW_LIMIT: ESA's default stays -1; IRSA's
default None preserves its client setting. Clients retain stateful row limits.

DR1 catalogs share microJy profile mapping, the existing 3631 Jy AB convention,
measurement errors and 0.001 mag rounding. Vega requires the same explicit,
source-recorded configuration. No official DR1 Vega defaults are supplied.
Q1 valid-value profile/aperture outputs are retained.

DR1 spectra provide WAVELENGTH in micrometres, FLUX, ERROR, MASK and optional
QUALITY/NDITH. Compare physical quantities with unit conversion: native scaled
units are preserved. Negative variance stays invalid; bins/flags are retained.
Missing values can be masks or NaN depending on serialization. Saving uses
`<object_name or E{object_id}>_spectrum.fits` for both archives.

## Corrections made

- IRSA plural spectra accept release. Public modules and backends expose matching
  discovery/plural methods and reject non-DR1 releases before querying.
- ESA mosaic discovery accepts search_radius. IRSA public discovery defaults to
  survey/deep products, matching cutout retrieval. Repeat visits require
  selection="repeat". Lower-level irsa_images.discover_mosaics still lists all
  native classes; unverified candidates remain visible as error rows.
- IRSA adds metadata-only spectrum discovery, deduplicates exact file/HDU rows,
  and validates all metadata before downloading. An exact product absent from a
  nonempty result raises ValueError, matching ESA, instead of looking like no data.
- Both catalog helpers reject absent/failed results and validate IDs, finite
  coordinates, radii, limits and cutout sizes. IRSA records the selected catalogue.
- Q1 spectra reject multiple ESA output files or distinct IRSA legacy VOTable
  associations instead of choosing the first. Corrupt ESA FITS errors propagate
  and opened FITS resources close.
- **ESA Q1 invalid variance:** negative variance previously became zero
  uncertainty through clipping. It now produces NaN. Valid variance, FSCALE,
  fluxes, masks, bins and thresholds are unchanged. Downstream handling of those
  previously invalid bins can change. This is an explicit validity correction,
  without recalibration or a new scientific quality threshold.

## Retained differences and limitations

- IRSA also exposes mode catalogs; coverage, IDs and collections differ.
- ESA grisms come from lambda_range. IRSA association metadata lacks grism;
  discovery leaves an empty string and rejects a grism filter with
  NotImplementedError before querying. Retrieval uses the documented wavelength
  rule. IRSA may therefore detect singular ambiguity after download while ESA
  can reject it beforehand.
- ESA needs a distinguishing native source-patch identifier; shared-patch parents
  remain blocked. IRSA enforces file/HDU with SpectrumDM. Paths/provenance differ.
- ESA server cutouts can have 31 pixels where IRSA's local 3-arcsecond cutout has
  30. Both preserve native pixels/calibration/sky coordinates, without forced
  resampling. ESA lacks parent-HDU/header evidence that IRSA can inspect; those
  fields remain unknown. Format support, size/byte limits and timeout/session/
  filesystem controls remain transport-specific.
- IRSA discovery preserves S3-header failures as error rows, including unknown
  classes; retrieval blocks until resolved. ESA checks TAP metadata at discovery
  and processing headers at retrieval.
- Legacy Q1 catalog presentation is retained: IRSA fills masked floating values
  with NaN while ESA may retain masks. DR1 raw masks are preserved by both. Q1
  image helpers retain first-mosaic behavior; strict product/tile ambiguity is DR1.
- The legacy Astroquery ESA Q1 client can log an HTTP failure and return no files.
  Its wrapper cannot distinguish that from an empty response. Direct DR1
  transports propagate failures; replacing the legacy Q1 route needs separate
  verification. Identical errors for every Q1 service failure are not guaranteed.

## Validation and references

`tests/test_backend_parity.py` compares catalog geometry/values/photometry,
defaults, inputs, empty/error/ambiguity behavior, discovery, real decoding, units,
saving and provenance. Existing tests cover bounded transport and image WCS/
calibration. Companion `euclid_dr1/tests/test_backend_integration_parity.py` runs
both actual adapters through the same manifest/FITS/cache workflow with the same
spectral data. This is not a new empirical quality/availability comparison.

Combined regression run: **368 tests and 24 subtests passed** (266 tools tests,
102 search tests). Command:

```sh
PYTHONPATH=/workspace/Euclid_tools/src:/workspace/Euclid_search \
  /workspace/.venvs/euclid-tools/bin/python -m pytest \
  /workspace/Euclid_tools/tests /workspace/Euclid_search/euclid_dr1/tests -q
```

Official formats, sources and earlier live evidence: [ESA catalogs](ESA_DR1_CATALOG.md),
[IRSA catalogs](IRSA_DR1_CATALOG.md), [ESA products](ESA_DR1_PRODUCTS.md),
[IRSA spectra](IRSA_DR1_SPECTRA.md), [IRSA images](IRSA_DR1_IMAGES.md),
[Vega assumptions](DR1_VEGA_CALIBRATION.md), [provenance](PRODUCT_PROVENANCE.md).
No undocumented mappings, photometric recalibration or new scientific selection
criteria were introduced.
