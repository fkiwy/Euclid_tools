# Q1 and DR1-Foundation support

Current library support, reviewed 17 November 2026. The runnable archive examples
now default to DR1 deep survey; the singular library retrieval methods keep Q1
as their default. Public use requires only this repository and its dependencies.

## Catalog and product selection

| Feature | Q1 | DR1 |
|---|---|---|
| ESA catalogs | `q1.mer_catalogue` | Explicit wide/deep survey; wide is the library default |
| IRSA catalogs | `euclid_q1_mer_catalogue` | Wide/deep survey and mode collections; wide survey is the default |
| Profile magnitudes | AB and legacy Vega | AB; optional explicit user-supplied Vega convention |
| Aperture magnitudes | Existing VIS 4-FWHM and Y/J/H 2-FWHM outputs | Not enabled without verified calibration/units |
| ESA spectra | Current `q1` schema and Q1 FSCALE handling | Native source/grism association and bounded DataLink retrieval |
| IRSA spectra | Legacy VOTable or exact Q1 SpectrumDM file/HDU | Exact DR1 SpectrumDM file/HDU |
| ESA images | Legacy Q1 mosaic/cutout route | Verified wide-survey, repeat-visit and deep server cutouts |
| IRSA images | Legacy Q1 SIA route | SIA discovery and bounded S3 section reads |

`release="q1"` / `"dr1"` are case-insensitive. Catalog selection never changes
backend state or silently converts an ID to another release/survey. DR1 requests
do not fall back to Q1. Plural spectra and discovery methods are DR1-only and
accept explicit release selection; they default to DR1.

```python
from euclid_tools import get_backend

backend = get_backend("esa")  # Or "irsa".
objects = backend.retrieve_objects(
    58.1332495, -49.1830038, radius=5,
    release="dr1", catalogue="deep_survey",
)
if objects is not None:
    # Closest position is an example choice; inspect counterparts for science.
    source_id = str(objects[0]["object_id"])
    spectrum = backend.retrieve_spectrum(source_id, release="dr1", grism="red")
```

See [ESA catalogs](ESA_DR1_CATALOG.md), [IRSA catalogs](IRSA_DR1_CATALOG.md),
[ESA products](ESA_DR1_PRODUCTS.md), [IRSA spectra](IRSA_DR1_SPECTRA.md) and
[IRSA images](IRSA_DR1_IMAGES.md) for exact mappings, selectors and format limits.
Ambiguous spectra/images require explicit product/tile choices. Some ESA parents
cannot be enforced because their native patch identifiers are shared; they stay
blocked. IRSA grism classification needs wavelengths and is unavailable from
association metadata alone. [Backend parity](BACKEND_PARITY.md) lists differences.

## Verified profile schema and arithmetic

The original 19 shared ESA Q1/DR1-wide catalog names, units and descriptions are
recorded in [metadata](dr1_catalog_metadata.csv). Deep-survey metadata and IRSA
mappings are recorded in their catalog guides. Coordinates are degrees; VIS PSF
and Y/J/H template-fit fluxes and uncertainties are microJy. Raw flags/columns
are retained. ESA renames right_ascension/declination to ra/dec.

```
m_AB = -2.5 log10(flux_uJy * 1e-6 / 3631)
sigma_m = 2.5 / ln(10) * fluxerr / flux
```

Derived display columns round to 0.001 mag. Nonpositive flux has no ordinary
logarithmic magnitude; the existing arithmetic behavior is retained. No quality
cut, extinction correction or scientific calibration is inferred. DR1 Vega
requires an [explicit source-recorded convention](DR1_VEGA_CALIBRATION.md);
official release-specific Vega coefficients remain unverified.

## Spectral units, image metadata and unresolved assumptions

DR1 combined SIGNAL/VAR units already embed the physical scale; the DPDD marks
FSCALE deprecated. DR1 decoders preserve native scaled units and never multiply
FSCALE again. Q1 ESA spectra retain their separate FSCALE calculation. Negative
variance remains invalid, with no zero-uncertainty clipping. Retrieval retains
all bins/masks; scientific masking is an explicit downstream choice.

Image routes verify supported processing identity, native pixel units/calibration
and celestial WCS. This preserves metadata, without independently recalibrating
the survey. Unsupported planes, encodings or WCS forms fail explicitly. See the
[product provenance contract](PRODUCT_PROVENANCE.md) for evidence and unknowns.

Matching column metadata does not establish identical calibration, flag
semantics, processing, completeness or selection functions. Initial DR1 aperture
fields and fwhm lacked sufficient unit/description metadata for a scientific
migration. Aperture corrections, morphology cuts, mask policies and empirical
photometric relations need release-applicable validation. No production
candidate-selection policy is supplied by this library.

## Sources and validation

[Validation and reproduction](VALIDATION.md) distinguish current offline checks
from dated, small live archive samples. The original three-row catalog check
verified metadata and independent four-band AB arithmetic, not calibration.
Unordered TOP samples are response snapshots, not representative survey samples.
An initial installed Astroquery transport failed in the managed host while
HTTPS TAP requests succeeded; later bounded native product checks are separately
recorded. This historical environment issue is not a claim about current service
availability.

Official sources, initially inspected 12 November 2026 and extended in the
archive guides:

- [ESA DR1 release](https://www.cosmos.esa.int/web/euclid/data-release-dr1),
  [data access](https://www.cosmos.esa.int/web/euclid/dr1-data), and
  [documentation](https://www.cosmos.esa.int/web/euclid/dr1-documentation).
- [ESA TAP](https://eas.esac.esa.int/tap-server/tap/sync): TAP_SCHEMA tables/columns;
  exact catalog queries and responses are preserved in the validation records.
- [DR1 DPDD](https://euclid.esac.esa.int/dr/dr1/dpdd/) and
  [combined spectra](https://euclid.esac.esa.int/dr/dr1/dpdd/sirdpd/dpcards/sir_combinedspectra.html):
  units, FSCALE deprecation and diagnostic columns.
- [Astroquery ESA guide](https://astroquery.readthedocs.io/en/latest/esa/euclid/euclid.html):
  release and product interfaces. Tools explicitly use q1 instead of the retired
  sedm/catalogue aliases; they do not inherit those old client defaults.
- [IRSA DR1 MER tutorial](https://caltech-ipac.github.io/irsa-tutorials/euclid-intro-mer-catalog/),
  [spectra tutorial](https://caltech-ipac.github.io/irsa-tutorials/euclid-intro-1d-spectra/)
  and [image tutorial](https://caltech-ipac.github.io/irsa-tutorials/euclid-intro-mer-images/).
