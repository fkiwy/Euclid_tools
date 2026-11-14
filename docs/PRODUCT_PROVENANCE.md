# Shared product provenance and workflow integration

`euclid_tools.provenance.product_provenance(product, source_id=None)` returns a JSON-safe record with schema `euclid.product-provenance/v1`; `provenance_json` serializes it deterministically. Inputs can be an Astropy Table, FITS image HDU, or metadata mapping. IDs are decimal strings, retaining signed 64-bit precision in CSV/JSON. Unknown fields are JSON `null`, not guesses. `metadata_version` is null for legacy products without versioned cards.

```python
from euclid_tools.backends import get_backend
from euclid_tools.provenance import product_provenance
backend = get_backend('irsa')
spectrum = backend.retrieve_spectrum(object_id)  # Q1 remains the default
record = product_provenance(spectrum)
# DR1 remains explicit; existing grism/product selectors are unchanged.
image = backend.retrieve_cutout(ra, dec, radius=2, cutout_size=3,
                               band='VIS', release='dr1', tiling='wide',
                               selection='survey', tile_id=tile_id)
image_record = product_provenance(image, source_id=object_id)
```

The common record includes release, archive, product type, source ID, source file/HDU, grism and reference, retrieval route/URL, dataset and association table, tile and band, mosaic class, processing dataset/plan/pipeline/version, and existing raw/header hashes. It reads native metadata without modifying arrays, WCS, calibration, units, mask policy or scientific selection.

FITS storage uses additive scalar cards: `PROVVER`, `PRODTYPE`, `RETROUTE`, `PRODBAND`, plus existing `RELEASE`, `ARCHIVE`, `OBJECTID`, `SRCPATH`, `SRCHDU`, `GRISM`, etc. Long path values use FITS long-string conventions. Existing DR1 `SRCS3` remains available alongside the common `SRCPATH`. A table's metadata and an image's header survive FITS round-trips; no nested Python dictionary is written into FITS.

| Product | Recorded retrieval route | Identity evidence |
| --- | --- | --- |
| IRSA DR1 spectrum | `tap-association+spectrumdm` | Object association path/HDU; existing tutorial wavelength-based grism rule and reference |
| IRSA Q1 spectrum | `tap-association+spectrumdm` or `tap-association+votable` | Native Q1 association path/HDU and `bandpass_name`, when supplied |
| IRSA DR1 image | `sia+s3-section` | Exact S3 URI, primary HDU 0, tile, verified processing headers and survey/repeat/deep classification |
| IRSA Q1 image | `sia+fsspec-section` | Existing SIA URL, primary HDU 0, processing cards copied if present; no DR1 mosaic classification |
| ESA Q1 spectrum | `esa-get_spectrum` | Requested source ID, delivered spectrum HDU 1, archive FILENAME/GRISM only when supplied |
| ESA Q1 image | `tap-mosaic+esa-get_cutout` | Selected mosaic path/tile and returned processing cards, when supplied |

A coordinate cutout does not establish an object association: `source_id` stays null unless native product metadata supplies it. `requested_source_id` separately records the workflow target; a conflicting native source ID raises an error. An ESA server-generated cutout's parent HDU is unknown unless explicitly supplied. Local temporary filenames are not promoted to archive source paths. Requested release metadata describes routing; its presence is not an independent calibration verification.

Q1 native `RGS`/other bandpass labels are preserved, not converted to DR1 red/blue semantics. DR1's existing grism rule remains unchanged. Missing Q1 processing metadata and old caches remain unknown. Exact archive/version verification remains limited to the previously documented supported formats; provenance adds no new release/product support.

## euclid_dr1 integration

Install this tools branch and use the companion [Euclid_search PR #1](https://github.com/fkiwy/Euclid_search/pull/1). Existing commands/configuration continue to work:
- `save_spectra` adds the JSON column `spectrum_provenance` for fresh and cached products.
- `plot_images` adds `image_VIS_provenance` (likewise Y/J/H) and associates the target ID separately from the mosaic. Existing DR1 FITS cutouts and manifests remain available.
- New ESA Q1 spectrum downloads also carry the common scalar FITS metadata. Q1 output directories, retrieval choices and defaults remain unchanged.
- Failed/missing products clear their shared provenance field; a missing record is not evidence of successful retrieval. Older FITS caches remain readable without forced downloads; absent metadata remains null.
- Workflow JSON sidecars name the product provenance schema in archive settings. Existing per-product hashes and raw metadata stay available.

## Evidence and verification

This change introduces a local metadata contract, not a new undocumented archive mapping. Official sources remain [IRSA SpectrumDM OpenAPI](https://irsacloud.ipac.caltech.edu/api/spectrumdm/openapi.json), [IRSA DR1 spectrum tutorial](https://caltech-ipac.github.io/irsa-tutorials/euclid-intro-1d-spectra/), [IRSA cloud tutorial](https://caltech-ipac.github.io/irsa-tutorials/euclid-cloud-access/), and the exact official sources and live samples recorded in [spectra](IRSA_DR1_SPECTRA.md), [images](IRSA_DR1_IMAGES.md), and [catalogs](IRSA_DR1_CATALOG.md). No new scientific interpretations of Q1 flags or DR1 processing versions are introduced.

Automated validation checks real adapters/decoders/section readers with bounded fixture transports, FITS round-trips, precise IDs, units/values/masks, processing identity, Q1 defaults, missing metadata, source conflicts, workflow manifests and cache reuse. No new live survey downloads are needed for this metadata-only change; previous small public validation samples remain documented separately.
