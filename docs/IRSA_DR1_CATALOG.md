# IRSA catalog release selection

IRSA catalog helpers now accept keyword-only `release="q1"` (default) or
`release="dr1"`, case insensitive. This fixes the backend call that previously
raised `TypeError: retrieve_objects() got an unexpected keyword argument 'release'`.

```python
from euclid_tools import get_backend
from euclid_tools.shared import print_results

backend = get_backend("irsa")
results = backend.retrieve_objects(80.2636961199, -51.2542825103,
                                  search_radius=5, release="dr1")
if results is not None:
    print_results(results)
    # ID lookup preserves every matching row: object_id is not guaranteed unique.
    rows = backend.retrieve_object_by_id(str(results["object_id"][0]), release="dr1")
```

Both direct `irsa_tools` functions and the backend support release selection,
as does `irsa_tools.print_catalog_info(release="dr1")`. `radius` and
`search_radius` retain their existing alias behavior. `client=` optionally
injects an Astroquery IRSA-compatible client for catalog retrieval.

Mappings verified against the official IRSA tutorial and live TAP metadata:

| Release | IRSA catalog | Coordinates | Photometry |
|---|---|---|---|
| Q1 (default) | `euclid_q1_mer_catalogue` | `ra`, `dec` | Existing profile/aperture AB and Vega outputs |
| DR1 | `euclid_dr1_mer_catalogue_wide_survey` | `ra`, `dec` in degrees | VIS PSF and Y/J/H template-fit AB outputs |

DR1 flux/error columns are `flux_vis_psf` / `fluxerr_vis_psf` and
`flux_<y,j,h>_templfit` / `fluxerr_<y,j,h>_templfit`, in µJy. The existing
3631 Jy AB convention and 0.001 mag rounding are used, including the matching
flux uncertainty. Raw fields, flags, exact integer IDs and `tileid` are retained.
No DR1 aperture or Vega calibration is assumed and no quality filtering occurs.
Results retain the same separation-in-arcseconds behavior as Q1 cone searches.
Empty results return None. Invalid releases fail before archive access, and an
archive error never triggers a Q1 fallback. ID lookups retain multiple matching
rows without deduplication or silently selecting one.

`print_results` now prints the available summary columns, so DR1 results do not
require the Q1-only aperture/Vega columns. Complete Q1 summaries are unchanged.
The example's `release` variable is passed to catalog metadata, cone search, ID
lookup and spectrum retrieval. For DR1 it explicitly selects red spectra and
skips images, since DR1 image retrieval is not implemented. Its default remains
Q1. Changing the release does not establish coverage at the example coordinates
or guarantee an unambiguous spectrum for the nearest object.

## Sources and validation

Read 13 November 2026:

- [Official IRSA DR1 MER catalog tutorial](https://caltech-ipac.github.io/irsa-tutorials/euclid-intro-mer-catalog/),
  sections “The MER catalog tables” and “Search the MER catalog around a position”.
- [IRSA TAP](https://irsa.ipac.caltech.edu/TAP/sync), queried against
  `TAP_SCHEMA.columns` for the selected DR1 table. The 19 verified metadata rows
  are preserved in `irsa_dr1_catalog_metadata.csv`; they cover coordinates,
  flux/error pairs and flags. Table names and units were verified directly;
  identical calibration or flag semantics across archives/releases are not assumed.

Regression tests reproduce the reported backend call, summary printing, negative
IDs and multiple matches, release rejection, empty responses, metadata selection,
and default/explicit Q1 aperture and Vega outputs. The companion Euclid_search
suite is run alongside these tests; its ESA catalog configuration is unchanged.

Live validation succeeded through the actual installed IRSA client: a 0.1-arcsec
cone around `(80.26369611992676, -51.25428251026353)` returned one DR1 record,
`-802636961512542825`. Exact-ID lookup and summary printing also succeeded;
all four AB magnitudes/errors agreed with independent µJy calculations within
0.000501 mag. Exact inputs, UTC timestamp and values are recorded in
`irsa_dr1_catalog_validation.json`. An earlier untargeted TOP 3 probe timed out
after 30 seconds; the targeted tests succeeded without a fallback or bulk data.
Combined regression validation: **137 passed, 24 subtests passed**.
