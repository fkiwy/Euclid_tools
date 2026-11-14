# IRSA DR1 MER image cutouts

Implemented 14 November 2026. Q1 remains the default of `retrieve_cutout`.
DR1 explicitly discovers science mosaics with SIA, verifies their FITS headers,
and reads only the requested section from anonymous S3. No Q1, full-image,
nearest-tile, or repeat-visit fallback is used.

## Install and use

Update the checkout and its dependencies (new direct dependencies: fsspec/s3fs):

```sh
python -m pip install -e .
```

```python
from euclid_tools import get_backend

backend = get_backend("irsa")
# Discover each collection independently; all candidates/statuses remain visible.
for tiling in ("wide", "deep"):
    mosaics = backend.discover_mosaics(273.167, 64.3687,
        search_radius=1, band="VIS", tiling=tiling)
    print(mosaics["obs_id", "tile_id", "kind", "status", "product", "error"])

# The tutorial position overlaps two deep tiles: choose one explicitly.
hdu = backend.retrieve_cutout(273.167, 64.3687,
    search_radius=1, cutout_size=3, band="VIS", release="dr1",
    tiling="deep", tile_id="101829777")
hdu.writeto("deep_VIS_cutout.fits", overwrite=True)

# Official wide-survey processing; repeat visits require selection="repeat".
hdu = backend.retrieve_cutout(273.167, 64.3687,
    search_radius=1, cutout_size=3, band="VIS", release="dr1",
    tiling="wide", selection="survey", tile_id="102158277")
```

The same functions are available directly from `irsa_tools`.
Discovery is a new DR1-only API. It accepts `radius`/`search_radius` in arcseconds.
`tiling` selects `euclid_dr1_wide_DpdMerBksMosaic` or
`euclid_dr1_deep_DpdMerBksMosaic`; `band` is VIS/Y/J/H. Only SIA science rows in
the requested band are inspected. The source bucket/key comes from `cloud_access`,
not a guessed file name. Bucket, DR1 prefix, tile and instrument directory must
agree with the selected collection and band.

Discovery returns a Table of every matching candidate, its URI, status/error,
classification, dataset, processing plan/version, calibration and header hash.
It reads only headers. Missing files or unverified metadata remain error rows.
The singular cutout call rejects an unresolved error candidate, and raises on
multiple matches. Use `tile_id` or `product="s3://..."` (an exact URI returned by
discovery) to resolve the choice. It returns None only when no candidate matches.
No selection ranks images by depth, S/N, proximity or archive order.

## Processing identity and calibration

The official IRSA tutorial identifies the official wide-survey mosaic by
`-DATA_RELEASE-` in `PPLANID`. Wide mosaics with `DATA_RELEASE_DEEP` are repeat
visits; deep-tiling mosaics with that marker are deep stacks. Unknown plans fail
verification. Wide defaults to `selection="survey"`; deep defaults to
`selection="deep"`. `selection="repeat"` is available only with wide tiling and
usually needs an exact product, because a tile has several repeat visits.
These labels describe the processing plan, not a new scientific selection.

The code requires the header's DATASETR, PPLANID, PIPDEFID, PSWREL, FILTER, BUNIT,
MAGZERO and celestial WCS. Both `DR1_R1` and `DR1_R2` were verified in live archive
headers; other dataset labels are explicitly unsupported pending verification.
The MER pipeline/plan relationship and FILTER are checked. The implemented pixel
format is an unscaled float32 primary image with a two-dimensional celestial WCS.
Distortion tables/SIP, alternate WCS frames, scaled integer encodings and unknown
processing plans are rejected rather than interpreted without verification.

Native pixel values are copied without normalization, photometric conversion,
reprojection or resampling. The original header is retained, including BUNIT,
MAGZERO, exposure and processing metadata. **VIS and NISP calibration differ**:
the live VIS samples used ADU/s and MAGZERO 24.5; Y used ELECTRON/s and MAGZERO
29.8. Those values are preserved, not imposed on other bands or products.

The crop WCS replaces the source primary WCS keywords (removing stale CD/PC
representations), and NAXIS1/2 are updated. Source-vs-crop world coordinates are
checked at three positions to within 0.0001 arcsec. CHECKSUM/DATASUM are removed
because the original checksum no longer describes the cropped data. Source pixel
origin/dimensions, S3 URI, tile, SIA collection/observation, processing class,
requested center/size, edge policy, bytes read and source-header hash are added.
`RAWHDR` hashes the parsed source header serialization, not the entire parent file.

The default edge policy is `mode="trim"`, retaining only the overlapping image
and marking CUTTRIM. `mode="strict"` requires the complete footprint. Empty or
nonoverlapping cutouts and unsupported WCS fail explicitly. No extrapolation or
zero-filled padding is introduced.

## Transfer limits

S3 is opened anonymously with cache/read-ahead disabled. A bounded file reader
rejects unbounded reads and reads exceeding 16 MiB before fetching them. FITS
`.section` is sliced by Cutout2D; the parent HDU `.data` is never used. Requested
cutouts are limited to one million pixels. Discovery headers are limited to
64 KiB each and 100 candidates. SIA requests at most 5001 metadata rows;
server overflow warnings and a reached limit are errors, not partial selections.
A budget error is distinct from a corrupt-FITS error. There is no whole-object
GET or alternative download route on a failure.

## Validation and reproducibility

```sh
python scripts/validate_irsa_dr1_images.py --output /tmp/irsa-dr1-images
python -m pytest tests -q
```

The live validation uses the actual IRSA SIA client and anonymous S3 with an audit
hook that rejects a non-ranged GetObject request and checks HTTP 206 for every
image GET. It retrieves four 3-arcsec (30x30-pixel) samples at
RA 273.167, Dec 64.3687: deep VIS, deep Y, official wide VIS, and one explicitly
selected repeat-visit VIS. Source/crop WCS, FITS roundtrip, processing identity
and calibration headers are checked. `irsa_dr1_images_validation.json` records
exact product URIs, source headers/hashes, requested ranges, HTTP range results,
transfer totals, package versions and UTC time. The parent files are roughly
0.42 GB (deep) and 1.47 GB (wide), while each validation case transfers under
1 MiB including discovery headers.

Offline tests use small synthetic mosaics and headers sampled from the live
products. They verify pixel identity, rotated WCS, calibration preservation,
no full-array reads, official/repeat/deep separation, overlap ambiguity, missing
objects, invalid metadata, limits, edge behavior, Q1 compatibility and backend
forwarding. Companion search tests exercise the actual section reader through
archive configuration, calibrated FITS output, manifest and PDF generation.
This validates transport/metadata/geometry, not survey calibration, processing
quality, completeness or the scientific equivalence of deep and wide products.

## Search workflow and example

The updated `irsa_tools_example.py` explicitly derives wide/deep tiling from the
chosen catalog and uses its selected object's tileid. It now downloads DR1 cutouts
with the same release, rather than skipping them or querying Q1. Use discovery
and exact products if a tile/processing class remains ambiguous.

In the companion Euclid_search branch, opt in with
`euclid_dr1/archive_dr1_irsa_images.json`. It keeps the existing ESA wide catalog
mapping and IRSA spectra, and adds IRSA wide official-survey images. It is a
separate config: existing catalog-only, spectrum-only and Q1 configs keep their
behavior. Set image_tiling=deep and image_selection=deep only if you intentionally
want deep images. Per-source `image_tile_id`, per-band `image_VIS_tile_id`, and
per-band `image_VIS_product` (likewise Y/J/H) resolve overlaps/repeats explicitly.
The image workflow saves calibrated FITS files with --cutouts-dir, a PDF preview,
and per-band provenance/hashes; it does not merge images or silently mix products.

## Exact official sources

Inspected 14 November 2026:

- [IRSA DR1 MER mosaic tutorial](https://caltech-ipac.github.io/irsa-tutorials/euclid-intro-mer-images/),
  sections 1–3 and 6: SIA collections and cloud_access, PPLANID selection,
  source header example, `.section` reads, WCS and edge trimming.
- [IRSA DR1 cloud access tutorial](https://caltech-ipac.github.io/irsa-tutorials/euclid-cloud-access/).
- [IRSA DR1 overview](https://irsa.ipac.caltech.edu/data/Euclid/docs/overview_dr1.html).
- Live SIA metadata from these named collections and FITS headers from their
  `nasa-irsa-euclid-dr1/dr1/MER/...` objects; exact object keys and header hashes
  are recorded in the validation JSON. In particular the live official wide VIS
  header verifies DR1_R2, beyond the DR1_R1 deep header shown in the tutorial.

ESA DR1 product retrieval remains disabled. No scientific thresholds or catalog
selection criteria were changed by this implementation.


Final combined validation: **169 passed, 24 subtests passed** on Python 3.12.
The live search workflow saved all four VIS/Y/J/H FITS cutouts and a PDF preview,
which was visually inspected. The earlier spectrum target at
(273.018929496624, 68.231870619389) also returned a verified 30x30 VIS cutout
from deep tile 101839130 using 9,440 bytes of FITS reads. Those results are in
the validation JSON. Only science-plane pixels are retrieved; associated weight,
RMS and pixel-quality image planes are outside this implementation.
