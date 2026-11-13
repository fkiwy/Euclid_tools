# DR1-Foundation support and validation (12 November 2026)

## Implemented

ESA `retrieve_objects`, `retrieve_object_by_id`, and `print_catalog_info` accept
keyword-only `release="q1"` (default) or `release="dr1"`, case insensitive.
The ESA backend forwards release arguments. Invalid releases fail before a query.
IRSA catalog/image APIs remain Q1-only. Explicit DR1 spectrum retrieval is now
available; see [IRSA DR1 spectra](IRSA_DR1_SPECTRA.md) for the 13 November extension.

Q1 retains positional arguments,
aperture/profile outputs, AB/Vega zero points and rounding. Its table is now
`q1.mer_catalogue`: a bounded live query against the retired
`catalogue.mer_catalogue` alias returned HTTP 400 (unknown table), whereas
`SELECT TOP 1 object_id FROM q1.mer_catalogue` returned a Q1 record.
This is a table-name repair, not a scientific mapping change. No global release
switch is used. An optional `client=` on catalog retrieval accepts an existing
Astroquery-compatible client, including one authenticated by the caller.

DR1 explicitly selects `dr1.mer_catalogue_wide_survey`, renames
`right_ascension`/`declination` to `ra`/`dec`, and keeps raw columns and flags.
VIS uses `flux_vis_psf`/`fluxerr_vis_psf`; Y/J/H use
`flux_<band>_templfit`/`fluxerr_<band>_templfit`. Coordinates are degrees and
flux/error are µJy. Derived DR1 magnitudes are **AB only**, using the existing
3631 Jy convention and 0.001 mag rounding:

```
m_AB = -2.5 log10(flux_uJy * 1e-6 / 3631)
sigma_m = 2.5 / ln(10) * fluxerr / flux
```

No DR1 aperture or Vega magnitudes are inferred. Nonpositive flux has no ordinary
logarithmic magnitude; existing arithmetic behavior is retained. There is no
flag filtering, threshold change, extinction correction or calibration adjustment.

```python
from euclid_tools import esa_tools
objects = esa_tools.retrieve_objects(80.0821501085, -51.2463548253,
                                     radius=2, limit=3, release="dr1")
obj = esa_tools.retrieve_object_by_id("-800821501512463548", release="dr1")
```

## What was tested

`dr1_validation.json` records the exact UTC timestamp, TAP endpoint, bounded
ADQL, metadata query and independently calculated magnitudes/errors.
`dr1_catalog_metadata.csv` records the 19 Q1/DR1 metadata pairs. The three-row
fixture is `tests/data/dr1_sample.csv`. TOP without ORDER BY is intentionally
non-deterministic: the fixture records one response, not representative sampling.

Reproduce the bounded network check (three science rows and <=100 metadata rows):

```sh
python scripts/validate_dr1.py --output /tmp/dr1-validation
python -m pytest tests -q
```

HTTPS TAP queries through `requests` succeeded; all 19 names, units and TAP
descriptions matched, and four-band AB magnitudes and errors agreed within
0.000501 mag of the rounded helper output. This verifies schema and arithmetic,
not calibration or the astrophysical validity of a selection. Installed
Astroquery's ESA transport failed with `Connection refused` in this managed
execution environment, despite HTTPS TAP working through its proxy. Direct
Astroquery end-to-end retrieval is therefore **not validated here**. Mocked
Astroquery-compatible clients test cone arguments, ID lookup, sorting, empty
responses, negative IDs, release isolation and integration; no fallback to Q1
is attempted on failure. Existing Q1 tests remain part of the regression suite.

## Product interfaces and scientific limits

ESA DR1 spectrum and cutout calls accept the release argument but raise
`NotImplementedError` **before network access**. They never use the Q1 product
paths for a DR1 request. This is the ESA catalog-only foundation; IRSA spectra are a separately enabled route.

Live TAP exposes `dr1.spectra_source`, `dr1.mosaic_product` and
`dr1.mer_cutouts_wide_survey`; their existence does not verify the DataLink
schema/retrieval selector, the release/processing identity of returned files,
cutout HDU selection, or how catalog sources link to a particular coadd.
Product metadata queried in this investigation had blank descriptions/units.
The current archive user-guide URL returned a JavaScript shell to the HTTP
client; its detailed interface contracts could not be inspected that way.

The official DR1 combined-spectra DPDD describes multiple source/dither and
combined extensions. For `m_COMBINED1D_SIGNAL`, WAVELENGTH is Angstrom, SIGNAL
has units `10**-16 erg/s/cm2/Angstrom`, VAR has the squared unit, and FSCALE is
**deprecated, with scaling embedded in the unit**. Reusing the Q1 decoder's
FSCALE multiplication could double-scale DR1 data. No ESA DR1 spectral decoder or DR1 image calibration has been enabled.
The separately documented IRSA SpectrumDM decoder preserves native units and bins.

DR1 aperture fields and `fwhm` exist in live TAP but had empty unit/description
metadata. They are not substitutes for the verified profile measurements.
Matching Q1/DR1 profile metadata does not establish equal calibration,
quality-flag meanings, processing, completeness or selection functions.
Vega offsets, aperture corrections, morphology cuts, spectral masks and
photometric type/distance calibrations still need release-specific validation.

## Exact official sources inspected

Accessed 12 November 2026:

- [ESA DR1 release page](https://www.cosmos.esa.int/web/euclid/data-release-dr1)
  identifies the release and links to the DR1 data model.
- [ESA DR1 data access](https://www.cosmos.esa.int/web/euclid/dr1-data)
  describes anonymous discovery and Astroquery/IVOA access.
- [ESA DR1 documentation](https://www.cosmos.esa.int/web/euclid/dr1-documentation)
  is the release documentation landing page.
- [Official TAP synchronous endpoint](https://eas.esac.esa.int/tap-server/tap/sync),
  including `TAP_SCHEMA.tables` and `TAP_SCHEMA.columns`; exact catalog queries
  and metadata response are checked in alongside this document.
- [DR1 DPDD index](https://euclid.esac.esa.int/dr/dr1/dpdd/) and
  [DR1 combined spectra product](https://euclid.esac.esa.int/dr/dr1/dpdd/sirdpd/dpcards/sir_combinedspectra.html),
  specifically the combined extension header and six-column table.
- [ESA archive user guide](https://euclid.dataspace.esa.int/community/help-pages/Euclid_Science_Archive_User_Guide.html)
  was attempted; only a JavaScript-required shell was returned.
- [Official Astroquery ESA API](https://astroquery.readthedocs.io/en/latest/esa/euclid/euclid.html)
  is the client reference; installed Astroquery source was inspected for
  `get_spectrum(schema=..., retrieval_type=...)` defaults (`sedm`, `ALL`).

No documentation justified altering the existing scientific cuts.

## Companion integration

Use this branch together with `fkiwy/Euclid_search` branch
`codex/dr1-foundation`; see its `euclid_dr1/DR1_FOUNDATION.md`.
Install this checkout or set the search workflow's `EUCLID_TOOLS_PATH` to it.
Run both suites in the same Python environment:

```sh
PYTHONPATH=/path/Euclid_tools/src:/path/Euclid_search python -m pytest \
  /path/Euclid_tools/tests /path/Euclid_search/euclid_dr1/tests -q
```

The search branch requires this tools change for DR1. Q1 workflows can keep
using their existing configuration. No PR should be merged automatically.

Initial catalog migration validation (12 November): **101 passed, 24 subtests passed** (Python 3.12,
12 November 2026). Both repository diffs passed `git diff --check`.
The exact search sample SQL also returned three live DR1 rows; its first row's
profile photometry passed through `gather_data.photometry` successfully.
