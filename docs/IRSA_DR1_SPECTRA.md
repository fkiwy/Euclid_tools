# IRSA DR1 spectrum retrieval

Implemented 13 November 2026. This adds an opt-in DR1 spectrum route; Q1 remains
the default of `irsa_tools.retrieve_spectrum`. IRSA catalog helpers also support explicit DR1 selection; see
[IRSA DR1 catalogs](IRSA_DR1_CATALOG.md). IRSA DR1 images are now available via [SIA/S3 cutouts](IRSA_DR1_IMAGES.md). ESA DR1 spectra/cutouts are also available through the [dedicated ESA route](ESA_DR1_PRODUCTS.md).

## API and explicit product selection

```python
from euclid_tools import irsa_tools

# New DR1-only plural API: retains every returned product without ranking.
products = irsa_tools.retrieve_spectra("2730189294682318706")
for spectrum in products:
    print(spectrum.meta["GRISM"], spectrum.meta["SRCPATH"], spectrum.meta["SRCHDU"])

# Existing singular API, explicitly requesting DR1 and the red grism.
spectrum = irsa_tools.retrieve_spectrum("2730189294682318706",
                                        release="dr1", grism="red")

# If repeated observations/overlapping tiles remain ambiguous, select an exact
# association shown by the plural call. It must belong to this object.
path = ("euclid/dr1/SIR_COMBINED_SPECTRA/101839130/"
        "EUC_SIR_W-COMBSPEC_101839130_20250908T142948.154355Z.fits")
spectrum = irsa_tools.retrieve_spectrum("2730189294682318706", release="dr1",
                                        product=(path, 6563))
```

The IRSA backend returned by `get_backend("irsa")` also exposes these methods.
`grism` accepts `red`, `blue` or None. Following IRSA's tutorial, it is determined
from the minimum finite, unmasked wavelength: below 1.1 µm is blue, otherwise
red. The provenance records that this is a wavelength rule, not a grism header.
It does not merge overlapping red/blue spectra or select the highest-S/N result.

The plural method returns a list; the singular method returns one QTable, None
for no matching association/grism/product, or raises ValueError for multiple
matches. Any failed product fetch raises rather than silently returning a partial
set. Malformed metadata, conflicting source identity, missing/incompatible units,
VOTable service errors and inconsistent variance/uncertainty are explicit errors.

`client=` optionally supplies an Astroquery IRSA-compatible TAP client;
`session=` supplies a requests-compatible HTTP client for SpectrumDM. The default
uses the installed IRSA client and requests. SpectrumDM requests have a 120-second
read timeout. There are no automatic retries or alternate-release fallbacks.

## Verified interface and numerical behavior

The association query selects `objectid,tileid,path,hdu` from
`euclid.objectid_spectrafile_association_dr1` for an exact signed 64-bit ID.
The query requests at most 101 rows; more than 100 associations is an explicit
error, not a silently truncated product list. Each selected association is fetched
through:

```
https://irsacloud.ipac.caltech.edu/api/spectrumdm/convert/{path}?dataset_id=euclid_combspec_dr1&hdu={hdu}
```

Only the selected HDU is transferred as a SpectrumDM VOTable. No parent FITS,
SSA, S3 or Q1 fallback is used. Responses are limited to 8 MiB each and HTTP errors
propagate. DR1 SSA was still unavailable according to the official guide.

Output columns are WAVELENGTH (µm), FLUX (from SIGNAL with its original scaled
physical unit), ERROR, MASK, and available QUALITY/NDITH. **No FSCALE multiplication**
is applied. ERROR is sqrt(VAR) in the flux unit when VAR is present, with supplied
UNCERTAINTY checked for agreement where both are finite; otherwise UNCERTAINTY
is used. Negative variance stays invalid (NaN); it is never clipped to zero.
Every bin survives, including negative/zero flux and nonzero mask values. Retrieval
does not impose a scientific mask, smoothing, normalization or quality selection.
The metadata records requested release, association table, dataset, source ID,
file/HDU, tile, inferred grism, URL and SHA-256 of the raw VOTable response. It
retains returned table metadata and rejects conflicting identity when present.
These identify the retrieval route, not an independently verified processing version.

## Validation and limits

A live call through the actual installed IRSA TAP client and SpectrumDM returned
531 bins for public object `2730189294682318706`, tile `101839130`, HDU `6563`.
The response was approximately 110 KiB, rather than the large multi-object FITS.
The original SIGNAL unit was `1e-16 erg / (Angstrom s cm2)` and VAR's unit was
its square. A physical reference bin, mask preservation and FITS roundtrip passed.
The exact URL, raw-response hash, UTC timestamp, package versions and reference
values are recorded in `irsa_dr1_spectrum_validation.json`.

Reproduce this one-object check:

```sh
python scripts/validate_irsa_dr1_spectrum.py --output /tmp/irsa-dr1
python -m pytest tests -q
```

`tests/data/dr1_irsa_spectrum_excerpt.vot` contains 12 rows, at original indices
0,1,2,50,100,150,200,250,300,350,400,530, from that public response. It was
re-serialized with Astropy while preserving values, units and table identity.
Offline tests cover both grisms, repeated products, no matches, malformed input,
HTTP/size/service errors, negative variance, identity checks, Q1 compatibility,
backend forwarding and roundtrip. Cross-repository tests cover the actual decoder
through the search archive adapter, manifest and cache selection.

Only one live red spectrum was validated. Blue and repeat-product selection have
synthetic offline tests; no live blue-spectrum calibration is claimed. IRSA's
examples discuss continuum offsets and unflagged artifacts in blue spectra.
Astropy warns that Angstrom and erg are deprecated VOUnit spellings; the supplied
scaled units parse correctly and remain unchanged. Catalog selection, downstream
mask policy, spectral fitting and scientific thresholds are unchanged.

## Official sources

Read 13 November 2026:

- [IRSA DR1 SIR tutorial](https://caltech-ipac.github.io/irsa-tutorials/euclid-intro-1d-spectra/),
  sections 2–4: association fields, SpectrumDM URL/dataset, units, uncertainty,
  wavelength-based red/blue identification, multiple products and masks.
- [IRSA DR1 cloud tutorial](https://caltech-ipac.github.io/irsa-tutorials/euclid-cloud-access/),
  sections 7–8: the public validation object, source file/HDU and embedded scale.
- [IRSA data access guide](https://irsa.ipac.caltech.edu/euclid-archive-documentation/data-access/),
  spectra section: DR1 association table and current lack of DR1 SSA access.
- [IRSA DR1 overview](https://irsa.ipac.caltech.edu/data/Euclid/docs/overview_dr1.html).

## Search integration

Use the companion Euclid_search branch `codex/dr1-foundation` and
`euclid_dr1/archive_dr1_irsa_spectra.json`. It explicitly combines the existing
ESA DR1 catalog mapping with IRSA DR1 spectra and defaults to red, appropriate
for the existing 1.22–1.88 µm spectral fitting range. ESA image/spectrum methods
remain disabled; the separate [IRSA DR1 image implementation](IRSA_DR1_IMAGES.md) is opt-in.

Combined validation: **128 passed, 24 subtests passed** on Python 3.12.
The actual search `save_spectra` command also wrote a successful live manifest
and 531-bin FITS cache for the same public object. Both diffs passed whitespace checks.


## MaskedColumn plotting compatibility (14 November 2026)

`SED.plot(spec_uncertainty=True)` now passes plain NumPy masked-array copies to
Matplotlib's uncertainty-envelope call. Astropy `MaskedColumn` subclasses can
become malformed during Matplotlib's polygon stack/slice operation and raise
`AttributeError: 'MaskedColumn' object has no attribute '_mask'`. The conversion
preserves numeric values in the existing SED display units and each array's
missing-data mask; it does not change retrieval, flux scaling, mask-bit policies,
interpolation, gap detection or scientific selection.

Two regression tests reproduce the offending stack/slice with the public DR1
spectrum excerpt and with genuinely masked wavelength/flux/error columns. They
fail with the reported AttributeError before the fix and pass afterward. They
also verify uncertainty bounds, units, masks and unchanged input data for that
plotting step. The saved full 531-bin response was run through the example's
WaveFlux → trim(1.21,1.89) → SED.to_flux_lambda → uncertainty plot pipeline and
rendered successfully to PDF and PNG; the PNG was visually inspected. Tests ran
with NumPy 2.5.4, Astropy 8.0.2 and Matplotlib 3.11.2, explicitly reproducing the
polygon operation so the regression does not depend on Matplotlib's installed
internal implementation. No archive download or package downgrade was needed.

Combined validation after the plotting fix: **148 passed, 24 subtests passed**.


## Q1 FITS association compatibility

Q1 association rows can now contain a parent FITS path plus `hdu`, rather than a direct VOTable URL. These use the same bounded single-HDU converter with the explicit `euclid_combspec_q1` dataset; the original direct-VOTable route remains supported. Multiple FITS associations fail explicitly instead of choosing a spectrum arbitrarily. No parent FITS download, flux rescaling, bin filtering or Q1 scientific selection change is introduced.

Official sources checked on 2026-11-14:
- [IRSA SpectrumDM OpenAPI](https://irsacloud.ipac.caltech.edu/api/spectrumdm/openapi.json): `/convert/{file_path}` accepts `dataset_id` and `hdu`; `DataSet` lists `euclid_combspec_q1` separately from DR1.
- [Official historical Q1 tutorial](https://github.com/Caltech-IPAC/irsa-tutorials/blob/b0a451f84fee74c04a41fa9383ca839afa34f88d/tutorials/euclid/3_Euclid_intro_1D_spectra.md): Q1 VOTable spectrum interface and physical columns.
- Live IRSA TAP association `euclid.objectid_spectrafile_association_q1`, object `2664850113649936423`, coordinates `(266.48501135, 64.99364238)`, tile `102158582`: file `euclid/q1/SIR/102158582/EUC_SIR_W-COMBSPEC_102158582_2024-11-05T16:02:47.013610Z.fits`, HDU `2338`.

The exact single-HDU request returned 108,675 bytes and 531 bins. Wavelength is Angstrom before conversion to microns; SIGNAL/UNCERTAINTY carry scaled `1e-16 erg / (Angstrom s cm2)` units. The explicit Q1 dataset and legacy `euclid_combspec` alias returned identical wavelength, signal and uncertainty arrays. The actual retrieval adapter also returned all 531 bins. TAP/SSA requests were intermittently slow or timed out; the successful association lookup and bounded cloud conversion were verified separately. Automated regression tests cover the adapter, exact dataset/HDU URL, physical values, masks, invalid metadata, ambiguity and the legacy route.
