# ESA DR1 spectra and MER mosaic cutouts

ESA DR1 product retrieval is enabled through explicit `release="dr1"` calls. Q1 remains the default and retains its existing decoder, query and selection behavior. DR1 uses a separate decoder and bounded public TAP/DataLink/cutout requests, independent of Astroquery's changed `source_id`/`ids` method signatures. No Astroquery upgrade is required for this product route.

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

The ESA example now has an explicit `release` variable, defaulting to Q1. Setting it to DR1 chooses this verified wide-survey target and red spectra. Image overlaps require an explicit tile or product path; no first/nearest result is selected by the DR1 methods.

## Spectrum association and decoding

Query `dr1.spectra_source` for the exact integer `source_id`, retaining native `file_path`, `file_name`, `hdu_index`, and `lambda_range` (`RGS` or `BGS`). Explicit ADQL columns avoid unrelated array-valued fields. CSV metadata transport preserves identifiers as strings and avoids malformed VOTable array-size declarations encountered from this archive. HTTP/XML/query errors are not treated as an empty catalog.

Request `/sas-dd/data` with `ID="dr1 <source_id>"`, `RELEASE="dr1"`, `LINKING_PARAMETER="SOURCE_ID"`, `RETRIEVAL_TYPE="SPECTRA_RGS"` or `SPECTRA_BGS`, and `USE_ZIP_ALWAYS=true`. The response must contain exactly one bounded FITS file with a primary HDU and one extracted spectrum table. Check `EXTNAME` and any supplied source-ID/release cards. The archive's parent `hdu_index` is provenance; the extracted download's table is HDU 1.

DR1 SIGNAL and VAR already carry scaled physical units (`10**-16 erg/s/cm2/Angstrom` and its square). Per the official DPDD, **FSCALE is deprecated because its scale is embedded in the units**. The decoder preserves native flux units/values, converts wavelength to microns, computes `sqrt(VAR)` without clipping negative variance, and retains every bin and original mask/quality/NDITH value. It never multiplies FSCALE or applies Q1 masks. A nonphysical/missing unit fails explicitly.

`retrieve_spectra` returns each uniquely associated grism; the singular method rejects multiple results. Multiple native file/HDU associations for the same grism are rejected before product download. An optional `product=(full_esa_path, hdu)` checks a unique association but cannot disambiguate several associations for one grism: the verified source-ID route cannot enforce a parent file. The source-patch route documented by newer Astroquery was **not validated** (a small test timed out); it is not implemented. No silent ranking, merging or Q1/IRSA fallback occurs.

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

## euclid_dr1 integration

Use both updated branches and install tools with `python -m pip install -e /path/Euclid_tools`. New opt-in config: `euclid_dr1/archive_dr1_esa_products.json`. Original catalog-only DR1 configurations stay disabled for products, and Q1/IRSA configurations retain their routes. The config selects ESA red spectra and official wide-survey images. Per-row grism and image tile/product selectors work as before; shared manifests and FITS caches now support ESA DR1 and validate the cached archive/dataset/selection.

```sh
python -m euclid_dr1.save_spectra --archive-config euclid_dr1/archive_dr1_esa_products.json --input /tmp/dr1-targets.csv --output /tmp/esa-spectra.csv --spectra-dir /tmp/esa-spectra
python -m euclid_dr1.plot_images --archive-config euclid_dr1/archive_dr1_esa_products.json --input /tmp/dr1-targets.csv --output /tmp/esa-images.csv --cutout-size 3 --cutouts-dir /tmp/esa-cutouts --plots-dir /tmp/esa-plots
```

Input example: `object_id,release,ra,dec,image_tile_id` followed by `2731837857643705492,DR1,273.18378579,64.37054921,102158277`. Use separate outputs or `--refresh` when switching spectrum backends; mismatched caches fail explicitly.

## Exact official sources and validation

Sources checked 2026-11-14:
- [ESA DR1 access](https://www.cosmos.esa.int/web/euclid/dr1-data): public archive/product access and Astroquery interfaces.
- [Official Astroquery ESA guide](https://astroquery.readthedocs.io/en/latest/esa/euclid/euclid.html), sections MER Cutouts and Spectra: DataLink/download and mosaic selection workflow.
- [Official Astroquery implementation](https://github.com/astropy/astroquery/blob/main/astroquery/esa/euclid/core.py): native request parameters, changed signatures, explicit release, retrieval types and cutout radius.
- [DR1 combined-spectra DPDD](https://euclid.esac.esa.int/dr/dr1/dpdd/sirdpd/dpcards/sir_combinedspectra.html): scaled SIGNAL/VAR units, deprecated FSCALE, wavelength and diagnostic columns.
- [Live ESA TAP](https://eas.esac.esa.int/tap-server/tap/sync), `TAP_SCHEMA.columns`, `dr1.spectra_source`, `dr1.mosaic_product`, and DataLink `/sas-dd/datalink/links`: verified native fields, grism selectors, release/path/HDU and processing identities. The DataSpace user-guide URL returned a JavaScript shell; it was not used to infer missing details.

[Recorded live adapter validation](esa_dr1_products_validation.json) and `scripts/validate_esa_dr1_products.py` cover red source `2731837857643705492` (531 bins), blue source `2681886189660473837` (433 bins), wide-survey VIS/Y/J/H on tile 102158277, deep VIS on tile 101839130, and one exact wide repeat VIS. All eight passed, including FITS re-serialization, with **89,403 response bytes** total including metadata. Snapshot full-adapter times were ~0.5–1.3 s per image, ~1.3 s for blue and ~7.7 s for red (including metadata); these are not general performance guarantees.

Automated tests use a five-bin ECSV excerpt of that actual red spectrum and an actual cutout header with synthetic pixels. They exercise scaling, invalid variance, source/units, grism/ambiguity, transport/size errors, CSV integer precision, processing/WCS/calibration preservation, persisted manifests/caches and Q1 regression. No bulk survey data were downloaded. Provenance/metadata preservation is not independent recalibration; Q1/DR1 flag equivalence and production candidate selection remain unresolved. No scientific thresholds, magnitude systems, catalog flux units, masks or fitting criteria change.

Final combined validation: **202 passed, 24 subtests passed** on Python 3.12; both diffs pass whitespace checks. The actual `euclid_dr1.save_spectra` and `euclid_dr1.plot_images` commands also completed against live ESA for the wide-survey target, saving the 531-bin spectrum, all four calibrated FITS cutouts and a visually inspected PDF preview. [Workflow validation record](esa_dr1_workflow_validation.json) captures the shared identities and hashes; this separate small workflow run is additional to the eight-check adapter transfer total above.
