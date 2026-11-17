"""Compare a Q1 Euclid spectrum with bundled FluxComp templates.

Install from this repository with `python -m pip install -e .`, then run
`python compare_spectrum_example.py`. The script uses the library's Q1 default
through get_backend("esa") at (266.4850113, 64.9936424) degrees, with a
5-arcsecond cone. It chooses the closest catalog entry and retrieves its spectrum.
For runnable DR1 catalog/product examples, use esa_tools_example.py or
irsa_tools_example.py; their release and survey choices are explicit.

This demonstration applies a specific Q1 mask predicate and compares against
Theissen+2022 templates smoothed with smooth_window=10. An alternative
Burgasser+2017 call is commented below. SED.compare uses trim_wave=True and
metric="reduced-chi2" for one best match. These template/filter settings and the
library's fitting conventions are analysis choices, not a validated scientific
selection or an uncertainty calibration. Verify mask definitions before adapting
them to another release; retrieval itself does not apply this predicate.

The PNG comparison plot includes the spectrum, template and uncertainty display.
By default it is written to the system temporary directory and opened in an
external viewer. Set output_dir/open_plot/plot_format on sed.plot to change
these choices. For headless use, configure MPLBACKEND=Agg before launching.
Empty catalog/spectrum results print a message; transport/read errors can raise.
The legacy ESA Q1 client can log a retrieval failure and return no files.

Uses euclid_tools.shared and flux_comp.core (SED, WaveFlux, TemplateProvider).
Dependencies and the template files are installed with the package. The archive
examples and the Q1 snippets in README.md explain the separate retrieval APIs.
"""

import warnings

from astropy.utils.exceptions import AstropyWarning

from euclid_tools import get_backend
from euclid_tools.shared import create_object_name, remove_bad_values
from flux_comp.core import SED, WaveFlux, TemplateProvider

backend = get_backend('esa')

warnings.simplefilter("ignore", category=AstropyWarning)

# ------------------------------
# Compare spectrum to templates
# ------------------------------

# Define object coordinates
# ra, dec = 58.1332495, -49.1830038
# ra, dec = 59.7913643, -47.6826163
ra, dec = 266.4850113, 64.9936424

# Specify search radius in arcseconds
search_radius = 5

# Retrieve objects based on coordinates
results = backend.retrieve_objects(ra, dec, search_radius=search_radius)

if results:
    result = results[0]
    object_id = str(result["object_id"])
    ra = result["ra"]
    dec = result["dec"]

    print(f"Object found at RA: {round(ra, 7)}, Dec: {round(dec, 7)}")

    # Retrieve the spectrum for the object
    data = backend.retrieve_spectrum(object_id)

    if data and len(data) > 0:
        # Select template(s) for comparison (Theissen+2022)
        provider = TemplateProvider()
        template_name = "Theissen+2022"
        templates = provider.get_Theissen_2022_templates(smooth_window=10)

        # Alternative template set:
        # template_name = "Burgasser+2017"
        # templates = provider.get_Burgasser_2017_templates()

        # Mask bad values in the spectrum to improve comparison results
        data = remove_bad_values(data, mask_func=lambda mask: (mask % 2 == 1) | (mask >= 64))

        # Create a WaveFlux object for comparison
        spectrum = WaveFlux(
            label="Spectrum", wavelength=data["WAVELENGTH"], flux=data["FLUX"], uncertainty=data["ERROR"]
        )

        # Create object name for plotting
        object_name = create_object_name(ra, dec, precision=2, shortform=False, prefix="J", decimal=False)

        # Set up the SED object and compare the spectrum to the templates
        sed = SED(object_name + " vs. " + template_name)
        sed.compare(
            spectrum,
            templates,
            trim_wave=True,
            number_of_matches=1,
            metric="reduced-chi2",
            add_stat_to_template_label=True,
        )

        # Convert to flux lambda and plot the results
        sed.to_flux_lambda()
        sed.plot(
            reference_on_top=False,
            spec_uncertainty=True,
            colors=["tab:blue", "tab:orange"],
            figure_size=(8, 5),
            show_grid=True,
            plot_format="png",
        )
    else:
        print("No spectrum found for the given object ID")
else:
    print("No object found for the given coordinates and search radius")
