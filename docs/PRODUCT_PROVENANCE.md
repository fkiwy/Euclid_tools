# Product provenance

`euclid_tools.provenance.product_provenance(product, source_id=None)` returns a JSON-safe record with schema `euclid.product-provenance/v1`; `provenance_json` serializes it deterministically. Inputs can be an Astropy Table, FITS image HDU, or metadata mapping. IDs are decimal strings, retaining signed 64-bit precision in CSV/JSON. Unknown fields are JSON `null`, not guesses. `metadata_version` is null for legacy products without versioned cards.

```python
from euclid_tools.backends import get_backend
from euclid_tools.provenance import product_provenance
backend = get_backend('esa')
object_id = '2731837857643705492'
ra, dec = 273.18378579, 64.37054921
spectrum = backend.retrieve_spectrum(object_id, release='dr1', grism='red')
if spectrum is not None:
    record = product_provenance(spectrum)
image = backend.retrieve_cutout(ra, dec, radius=2, cutout_size=3,
                               band='VIS', release='dr1', tiling='wide',
                               selection='survey', tile_id='102158277')
if image is not None:
    image_record = product_provenance(image, source_id=object_id)
```

The common record includes release, archive, product type, source ID, source file/HDU, grism and reference, retrieval route/URL, dataset and association table, tile and band, mosaic class, processing dataset/plan/pipeline/version, and existing raw/header hashes. It reads native metadata without modifying arrays, WCS, calibration, units, mask policy or scientific selection.

ESA DR1 spectra additionally record optional `linking_parameter`/`linking_id`
from `LINKPAR`/`LINKID`. A source-patch identifier can exceed int64 and remains an
exact decimal string, separately from the physical source ID. Earlier products
lacking these additive fields report null. `tap-sourcepatch+esa-datalink` denotes
a native patch identifier uniquely associated with the selected parent/grism.

FITS storage uses additive scalar cards: `PROVVER`, `PRODTYPE`, `RETROUTE`, `PRODBAND`, plus existing `RELEASE`, `ARCHIVE`, `OBJECTID`, `SRCPATH`, `SRCHDU`, `GRISM`, etc. Long path values use FITS long-string conventions. Existing DR1 `SRCS3` remains available alongside the common `SRCPATH`. A table's metadata and an image's header survive FITS round-trips; no nested Python dictionary is written into FITS.

| Product | Recorded retrieval route | Identity evidence |
| --- | --- | --- |
| IRSA DR1 spectrum | `tap-association+spectrumdm` | Object association path/HDU; existing tutorial wavelength-based grism rule and reference |
| IRSA Q1 spectrum | `tap-association+spectrumdm` or `tap-association+votable` | Native Q1 association path/HDU and `bandpass_name`, when supplied |
| IRSA DR1 image | `sia+s3-section` | Exact S3 URI, primary HDU 0, tile, verified processing headers and survey/repeat/deep classification |
| IRSA Q1 image | `sia+fsspec-section` | Existing SIA URL, primary HDU 0, processing cards copied if present; no DR1 mosaic classification |
| ESA DR1 spectrum | `tap-association+esa-datalink` or `tap-sourcepatch+esa-datalink` | Native source/grism association and uniquely enforceable parent file/HDU |
| ESA DR1 image | `tap-mosaic+esa-get_cutout` | Native mosaic association and returned cutout processing/calibration/WCS; parent HDU/header unknown |
| ESA Q1 spectrum | `esa-get_spectrum` | Requested source ID, delivered spectrum HDU 1, archive FILENAME/GRISM only when supplied |
| ESA Q1 image | `tap-mosaic+esa-get_cutout` | Selected mosaic path/tile and returned processing cards, when supplied |

A coordinate cutout does not establish an object association: `source_id` stays null unless native product metadata supplies it. `requested_source_id` separately records the workflow target; a conflicting native source ID raises an error. An ESA server-generated cutout's parent HDU is unknown unless explicitly supplied. Local temporary filenames are not promoted to archive source paths. Requested release metadata describes routing; its presence is not an independent calibration verification.

Q1 native `RGS`/other bandpass labels are preserved, not converted to DR1 red/blue semantics. DR1 ESA grisms use native associations; IRSA uses its documented wavelength rule. Missing Q1 processing metadata and old caches remain unknown. Exact archive/version verification remains limited to the previously documented supported formats; provenance adds no new release/product support.

## Evidence and verification

This change introduces a local metadata contract, not a new undocumented archive mapping. Official sources remain [IRSA SpectrumDM OpenAPI](https://irsacloud.ipac.caltech.edu/api/spectrumdm/openapi.json), [IRSA DR1 spectrum tutorial](https://caltech-ipac.github.io/irsa-tutorials/euclid-intro-1d-spectra/), [IRSA cloud tutorial](https://caltech-ipac.github.io/irsa-tutorials/euclid-cloud-access/), and the exact official sources and live samples recorded in [spectra](IRSA_DR1_SPECTRA.md), [images](IRSA_DR1_IMAGES.md), and [catalogs](IRSA_DR1_CATALOG.md). No new scientific interpretations of Q1 flags or DR1 processing versions are introduced.

Automated validation checks real adapters/decoders/section readers with bounded fixture transports, FITS round-trips, precise IDs, units/values/masks, processing identity, Q1 defaults, missing metadata, source conflicts, standalone serialization and legacy metadata. No new live survey downloads are needed for this metadata-only change; previous small public validation samples remain documented separately.


For ESA server cutouts, additive `cutout_header_sha256` (`CUTHDR`) identifies the
returned cutout header. The unavailable parent-header hash and parent HDU remain
unknown. These are separate evidence fields, not substitutes for an unobserved
parent product.
