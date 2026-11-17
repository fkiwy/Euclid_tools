# Explicit DR1 Vega conventions

ESA and IRSA DR1 catalog helpers return profile AB magnitudes by default. No
official DR1-specific Vega reference-flux table was verified in the sources below
on 2026-11-17. No automatic DR1 Vega coefficients are supplied. Q1 retains its
existing profile/aperture conversions and constants.

Both `retrieve_objects` and `retrieve_object_by_id`, including backend wrappers,
accept `vega_calibration=...` for DR1. This is an explicit **user-supplied
photometric convention**, not an assertion that the archive validates it.

Copy [the JSON template](dr1_vega_calibration.example.json), replace every null
and reference placeholder with justified values, then:

```python
import json
from euclid_tools import get_backend

calibration = json.load(open("my_dr1_vega_calibration.json", encoding="utf-8"))
backend = get_backend("esa")  # The same option works with "irsa".
objects = backend.retrieve_objects(
    ra, dec, radius=5, release="dr1", catalogue="wide_survey",
    vega_calibration=calibration,
)
```

Required fields:

| Field | Meaning |
|---|---|
| `release` | Exactly `DR1`; the opt-in is rejected for Q1 |
| `zero_points_jy` | Positive finite zero-magnitude flux densities in **Jy**, for VIS/Y/J/H |
| `source` | Source/version or derivation of those reference flux densities |
| `passband_reference` | Adopted passbands and their versions |
| `reference_spectrum` | Adopted Vega spectrum/version and magnitude normalization |

The values are effective bandpass reference flux densities. They are not image
`MAGZERO` values, detector count-rate zero points, or flux corrections. A source
citation alone does not establish applicability to DR1; that remains the user's
scientific assumption. Invalid/incomplete values fail before querying either
archive, including when the catalog would be empty. No SVO lookup or release
fallback happens automatically.

For each profile band, the optional outputs are `<band>_Vega_mag` and
`<band>_Vega_err`, alongside the usual AB columns:

```
m_Vega = -2.5 log10(flux_uJy * 1e-6 / zero_points_jy[band])
m_AB - m_Vega = 2.5 log10(3631 / zero_points_jy[band])
sigma_m = 2.5 / ln(10) * fluxerr / flux
```

The error includes measurement uncertainty only, treating the supplied reference
fluxes as fixed. Zero-point uncertainty, cross-band covariance, extinction and
color-dependent recalibration are not added. Rounding stays at 0.001 mag; raw
fluxes/flags and default AB arithmetic are preserved. Masked measurements stay
masked. No DR1 aperture photometry or new quality cuts are introduced.

Table metadata `VEGACAL` contains deterministic JSON with the complete convention;
`VEGASHA` hashes that JSON, and `VEGASTAT` labels it user-supplied and not
independently validated for DR1. These fields survive FITS round trips. The
examples expose `vega_calibration=None`; load the completed JSON to opt in.

For `Euclid_search`, place this same object under `vega_calibration` in a copy of
the selected archive JSON. `Archive.objects` passes it to the shared adapter;
`gather_data` adds the optional Vega columns and calibration identity to its CSV.
The existing `vis/y/j/h_mag`, colors, photometric types and distance inputs remain
AB. The archive provenance sidecar also retains the full supplied convention.

## Official sources checked and limits

- [DR1 MER Photometry Cookbook, Magnitudes](https://euclid.esac.esa.int/dr/dr1/dpdd/merdpd/merphotometrycookbook.html):
  catalog fluxes are microJy; AB conversion uses approximately 23.9 mag. It does
  not supply a Vega reference-flux table. The existing 3631 Jy convention is kept.
- [DR1 Calibration products overview](https://euclid.esac.esa.int/dr/dr1/expsup/Ch14.S1.html)
  links the [official Schirmer et al. calibration paper](https://eceb.astro.uni-bonn.de/public/dr1/schirmer-calibration.pdf).
  Sections 7.1–7.2 discuss the survey's zero-point modeling and absolute-flux
  standards; section 8 gives AB synthetic transformations. These are not a table
  of zero-magnitude Vega flux densities for DR1 catalog conversion.
- [SVO Filter Profile Service](https://svo2.cab.inta-csic.es/theory/fps/): generic
  Euclid passband reference values do not by themselves certify a DR1-specific
  calibration or Vega convention. They were not adopted as automatic DR1 values.

Automated tests exercise both catalog entry points, archives and wide/deep survey
selections, independent arithmetic using deliberately synthetic coefficients,
unchanged AB results, provenance serialization, masked data and early rejection
of invalid/Q1 opt-ins. These verify the implementation, not official coefficients
or the empirical accuracy of a DR1 Vega calibration. Official default coefficients
remain blocked pending a documented, release-applicable convention.
