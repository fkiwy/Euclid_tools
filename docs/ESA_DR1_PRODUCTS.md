# ESA DR1 spectra and MER mosaic cutouts

ESA DR1 product retrieval is enabled through explicit `release="dr1"` calls. Q1 remains the default. Its current schema and invalid-variance/error behavior are described below. DR1 uses a separate decoder and bounded public TAP/DataLink/cutout requests, independent of Astroquery's changed `source_id`/`ids` method signatures. No Astroquery upgrade is required for this product route.

```python
from euclid_tools.backends import get_backend
backend = get_backend('esa')
ra, dec = 273.18378579, 64.37054921
object_id = '2731837857643705492'
spectrum = backend.retrieve_spectrum(object_id, release='dr1', grism='red')
image = backend.retrieve_cutout(ra, dec, radius=1, cutout_size=3, band='VIS',
                                release='dr1', tiling='wide', selection='survey',
                                tile_id='102158277')
spectra = backend.retrieve_spectra(object_id, release='dr1')
mosaics = backend.discover_mosaics(ra, dec, release='dr1', band='VIS',
                                  tiling='wide', selection='survey')
```

The runnable ESA example defaults to `release="dr1"` and
`dr1_catalogue="deep_survey"` at `(58.1332495, -49.1830038)`. Setting
`release="q1"` keeps the coordinates and automatically omits the DR1 catalog
selector. The library's DR1 catalog default remains wide survey. The example
uses the catalog's native tile for images; direct calls require an explicit tile
or product when mosaics overlap. See [catalog selection](ESA_DR1_CATALOG.md).

## Spectrum association and decoding

Query `dr1.spectra_source` for the exact integer `source_id`, retaining native `file_path`, `file_name`, `hdu_index`, and `lambda_range` (`RGS` or `BGS`). Explicit ADQL columns avoid unrelated array-valued fields. CSV metadata transport preserves identifiers as strings and avoids malformed VOTable array-size declarations encountered from this archive. HTTP/XML/query errors are not treated as an empty catalog.

Request `/sas-dd/data` with `ID="dr1 <source_id>"`, `RELEASE="dr1"`, `LINKING_PARAMETER="SOURCE_ID"`, `RETRIEVAL_TYPE="SPECTRA_RGS"` or `SPECTRA_BGS`, and `USE_ZIP_ALWAYS=true`. The response must contain exactly one bounded FITS file with a primary HDU and one extracted spectrum table. Check `EXTNAME` and any supplied source-ID/release cards. The archive's parent `hdu_index` is provenance; the extracted download's table is HDU 1.

DR1 SIGNAL and VAR already carry scaled physical units (`10**-16 erg/s/cm2/Angstrom` and its square). Per the official DPDD, **FSCALE is deprecated because its scale is embedded in the units**. The decoder preserves native flux units/values, converts wavelength to microns, computes `sqrt(VAR)` without clipping negative variance, and retains every bin and original mask/quality/NDITH value. It never multiplies FSCALE or applies Q1 masks. A nonphysical/missing unit fails explicitly.

`retrieve_spectra` returns independently identifiable products, including multiple
parents for one grism. The singular method rejects multiple choices before
downloading and lists exact selectors. `product=(full_esa_path, hdu)` selects one
parent. Use `backend.discover_spectra(object_id, release="dr1", grism="red")` to
list `source_path`, `hdu_index`, native `linking_id` and `retrieval_supported`
without downloading a spectrum. Exact duplicate file/HDU rows represent one
product; conflicting grisms for a file/HDU are rejected.

The documented `SOURCEPATCH_ID` route is now verified for native identifiers that
uniquely select one parent within the requested grism. Parse the actual
`sourcepatch_id` array from TAP, retaining exact decimal strings even beyond
int64. Do not construct identifiers by concatenating other fields. For multiple
parents, select an identifier absent from the other same-grism associations and
send `LINKING_PARAMETER=SOURCEPATCH_ID`, `ID="dr1 <native sourcepatch_id>"`.
The FITS extension still identifies the physical `source_id`, which is checked.
Unique-grism products retain the existing `SOURCE_ID` route. No ranking,
combining or archive/release fallback occurs.

Some wide parents have no distinguishing identifier: every patch identifier is
also associated with a deep stack. They remain explicitly blocked, even with an
exact file/HDU selector. Discovery marks them `retrieval_supported=False`.
The source-ID route can return multiple tables with identical extension names;
neither response order nor those names establish their parent identity. A
multi-table source-patch response is likewise rejected, never assigned to the
requested parent. All selected routes are checked before product downloads.

Provenance adds `LINKPAR`/`LINKID` (`linking_parameter`/`linking_id` in the shared
record) and uses route `tap-sourcepatch+esa-datalink` when applicable. The exact
parent path/HDU, physical source ID, returned bytes hash and retrieval URL remain
recorded. The first lexicographic unique patch identifier is only an equivalent
retrieval alias for that parent, not a scientific preference among products.

Official selector documentation:
[EuclidClass.get_spectrum](https://docs.astropy.org/projects/astroquery/en/latest/api/astroquery.esa.euclid.EuclidClass.html#astroquery.esa.euclid.EuclidClass.get_spectrum).
The bounded [live validation script](../scripts/validate_esa_spectrum_associations.py)
records results in [association validation](esa_dr1_spectrum_associations_validation.json).
Some exploratory requests timed out; timeouts are propagated, without changing
the requested parent or archive. Q1 product selection/routing is unchanged.

On 2026-11-17, the actual adapter returned both red parents of source
`-629553356493020915` (HDU 3226 in DR1_R1 and HDU 3676 in DR1_R2), and the deep
stack for source `2714057875648364582` (HDU 9906). All three have 531 bins. The
last source's wide parent (HDU 622) shares its only patch identifier with that
stack and was blocked before a product request. Total received bytes across
these checks, including metadata, were 31,019. Selected source-patch requests
took approximately 51–52 seconds at the service, so this functionality does not
guarantee the speed of the unique-source-ID route. These are small product
checks, not survey-level completeness or quality validation.

## Mosaic selection, cutouts and provenance

Discover through `dr1.mosaic_product`, including `purpose`, `processing_mode`, `data_set_release`, instrument/filter, tile, and exact file path:

| Selection | Native processing_mode | Native purpose |
| --- | --- | --- |
| wide / survey (default) | WIDE | DATA_RELEASE |
| wide / repeat | WIDE | DATA_RELEASE_DEEP |
| deep / deep | DEEP | DATA_RELEASE_DEEP |

Metadata identifies the requested class before download; overlapping products require `tile_id` or an exact `product` path returned by discovery. Repeat visits commonly require the exact path as well as a tile. ESA product paths differ from IRSA S3 URIs and are never guessed or translated.

POST the selected `FILEPATH` and `POS=CIRCLE,ra,dec,radius_degrees` to `/sas-cutout/cutout`. The cutout radius is half `cutout_size` (arcseconds). ESA controls pixel rounding: the live 3-arcsecond requests returned **31 x 31 pixels**, which are preserved without local resampling. Only a 2D, unscaled float32 primary science image with verified DR1_R1/R2 processing metadata, band, BUNIT/MAGZERO and a supported celestial WCS is accepted. Returned processing identity must agree with discovery. WCS must contain the requested position. Other planes, distortion-table/alternate-WCS formats and full parent downloads remain unsupported.

Preserve the service's WCS, native pixel values, calibration and processing headers. Remove stale FITS checksum cards after adding provenance, and remove an undefined primary `EXTNAME` emitted by ESA (otherwise Astropy cannot re-save it). No calibration/header replacement is made. The service does not identify the parent's HDU, so shared `source_hdu` remains null for images. Likewise the parent header is not downloaded: `source_header_sha256` stays null; additive `cutout_header_sha256` (`CUTHDR`) hashes the returned cutout header instead. `RAWSHA` identifies returned product content (unpacked FITS for spectra; response bytes for images); `CUTBYTES` counts image response bytes, not all discovery traffic.

Shared [product provenance](PRODUCT_PROVENANCE.md) records ESA routing, release, source/file/HDU/grism where known, mosaic processing identity and hashes. Coordinate image targets remain separate from native object identity.

Limits: 100 metadata rows, 2 MiB metadata response, 16 MiB per product/ZIP member, one million cutout pixels, `0 < cutout_size < 100` arcseconds, and search radius up to 60 arcseconds. Requests use timeouts and streaming byte caps. These are retrieval limits, not scientific selection thresholds. Transport/decode failures are explicit; only genuinely empty metadata returns no product. Authenticated use can supply an appropriate requests-compatible `session`; an injected Astroquery `client` affects metadata queries only.

## Exact official sources and validation

Sources checked 2026-11-14:
- [ESA DR1 access](https://www.cosmos.esa.int/web/euclid/dr1-data): public archive/product access and Astroquery interfaces.
- [Official Astroquery ESA guide](https://astroquery.readthedocs.io/en/latest/esa/euclid/euclid.html), sections MER Cutouts and Spectra: DataLink/download and mosaic selection workflow.
- [Official Astroquery implementation](https://github.com/astropy/astroquery/blob/main/astroquery/esa/euclid/core.py): native request parameters, changed signatures, explicit release, retrieval types and cutout radius.
- [DR1 combined-spectra DPDD](https://euclid.esac.esa.int/dr/dr1/dpdd/sirdpd/dpcards/sir_combinedspectra.html): scaled SIGNAL/VAR units, deprecated FSCALE, wavelength and diagnostic columns.
- [Live ESA TAP](https://eas.esac.esa.int/tap-server/tap/sync), `TAP_SCHEMA.columns`, `dr1.spectra_source`, `dr1.mosaic_product`, and DataLink `/sas-dd/datalink/links`: verified native fields, grism selectors, release/path/HDU and processing identities. The DataSpace user-guide URL returned a JavaScript shell; it was not used to infer missing details.

[Recorded live adapter validation](esa_dr1_products_validation.json) and `scripts/validate_esa_dr1_products.py` cover red source `2731837857643705492` (531 bins), blue source `2681886189660473837` (433 bins), wide-survey VIS/Y/J/H on tile 102158277, deep VIS on tile 101839130, and one exact wide repeat VIS. All eight passed, including FITS re-serialization, with **89,403 response bytes** total including metadata. Snapshot full-adapter times were ~0.5–1.3 s per image, ~1.3 s for blue and ~7.7 s for red (including metadata); these are not general performance guarantees.

Automated tests use a five-bin ECSV excerpt of that actual red spectrum and an actual cutout header with synthetic pixels. They exercise scaling, invalid variance, source/units, grism/ambiguity, transport/size errors, CSV integer precision, processing/WCS/calibration preservation, FITS persistence and Q1 regression. No bulk survey data were downloaded. Provenance/metadata preservation is not independent recalibration; Q1/DR1 flag equivalence and production candidate selection remain unresolved. No scientific thresholds, magnitude systems, catalog flux units, masks or fitting criteria change.

## Q1 schema compatibility after ESA's schema migration

The archive's former `sedm` alias is no longer available. Q1 helpers now explicitly pass spectrum schema `q1` and query `q1.mosaic_product`. Astroquery clients with either `source_id` or `ids` are supported by inspecting the method signature, without retrying a request under a different release. A failed mosaic job raises a clear query error; a real empty result remains no image, and an empty cutout response is handled without indexing it.

Live checks on 2026-11-14 confirmed `q1.spectra_source` and `q1.mosaic_product` through TAP_SCHEMA. Source `2664850113649936423` at `(266.4850113, 64.9936424)` returned a 531-bin Q1 red spectrum through DataLink and a 31 x 31 VIS server cutout (3 arcseconds). The returned Q1 spectrum has unscaled physical SIGNAL/VAR units and `FSCALE=1e-16`; the existing Q1 scaling calculation is retained, unlike DR1's embedded scaled units. The existing Q1 decoder successfully read that live spectrum. Endpoint tests used bounded direct HTTP requests; wrapper/API compatibility is separately regression-tested with both Astroquery signatures. No bulk data, release fallback or scientific-selection changes were introduced.

The Q1 decoder retains its FSCALE calculation for valid data. Negative variance
now stays invalid (NaN), rather than being clipped to zero uncertainty.
Ambiguous outputs and malformed FITS raise explicitly. The legacy Astroquery
client itself may still log an HTTP failure and return no files; see
[backend behavior](BACKEND_PARITY.md). Current local verification commands and
historical snapshot scope are listed in [validation](VALIDATION.md).
