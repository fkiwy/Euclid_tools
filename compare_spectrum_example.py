"""
Example: Comparing a Euclid Spectrum to Templates Using FluxComp

This example demonstrates how to use the Euclid tools to retrieve a spectrum for a given object
and compare it to predefined templates using the `fluxcomp` tool. The script performs the following steps:

1. Retrieve an object catalog from the Euclid archive based on specified coordinates.
2. Retrieve the spectrum for the object.
3. Mask bad values in the spectrum data.
4. Select the appropriate template(s) for comparison (e.g., from "Theissen+2022").
5. Compare the spectrum to the templates using the reduced chi-squared metric.
6. Plot the results for visual comparison.

Dependencies:
- `euclid_tools.esa_tools` (for retrieving catalog objects and spectra)
- `euclid_tools.shared` (for utility functions such as object name creation)
- `flux_comp.core` (for comparing the spectrum to templates using the `SED` and `WaveFlux` classes)

Steps in the script:

1. **Retrieve Catalog Objects:**
   - The script begins by performing a search around a given RA and Dec to retrieve catalog objects from the Euclid archive. The `retrieve_objects()` function is used, and results are filtered by a specified search radius.

2. **Retrieve Spectrum:**
   - Once an object is found, the script retrieves its spectrum using the `retrieve_spectrum()` function.

3. **Mask Bad Values:**
    - The retrieved spectrum data is then masked for bad values to improve the comparison results.

3. **Template Selection:**
   - The script selects a set of templates for comparison. In this case, templates from "Theissen+2022" are used. These templates are retrieved using the `TemplateProvider()` and optionally smoothed (with a `smooth_window` parameter) before comparison.

4. **Comparison Using Reduced Chi-Squared:**
   - The comparison between the spectrum and the templates is done using the reduced chi-squared metric (`metric="reduced-chi2"`). This metric quantifies the difference between the observed spectrum and the template.

5. **Plotting the Results:**
   - The script plots the comparison results using the `SED` object's `plot()` method. The plot shows the spectrum and template(s) for visual inspection of the fit.

Notes for Users:
- **Template Selection:** This script uses templates from "Theissen+2022", but other templates (e.g., "Burgasser+2017") can be used by modifying the template name and retrieval method.
- **Comparison Metric:** The script uses the reduced chi-squared method (`metric="reduced-chi2"`) for template comparison. Users can change the comparison metric ("reduced-chi2","chi2", or "delta").
- **Plot Customization:** The plotting function visualizes the comparison between the spectrum and template. Users can adjust plotting parameters to customize the figure (e.g., size, legend, title, etc.).
  
Example Output:
- The script will generate a plot comparing the retrieved spectrum with the selected template(s). The plot will show the spectrum and template(s) overlayed with uncertainties and the reduced chi-squared fit.
"""

import warnings

from astropy.utils.exceptions import AstropyWarning

from euclid_tools.esa_tools import retrieve_objects, retrieve_spectrum, mask_bad_values
from euclid_tools.shared import create_object_name
from flux_comp.core import SED, WaveFlux, TemplateProvider

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
results = retrieve_objects(ra, dec, search_radius)

if results:
    result = results[0]
    object_id = str(result["object_id"])
    ra = result["ra"]
    dec = result["dec"]

    print(f"Object found at RA: {round(ra, 7)}, Dec: {round(dec, 7)}")

    # Retrieve the spectrum for the object
    data = retrieve_spectrum(object_id)

    if data and len(data) > 0:
        # Select template(s) for comparison (Theissen+2022)
        provider = TemplateProvider()
        template_name = "Theissen+2022"
        templates = provider.get_Theissen_2022_templates(smooth_window=10)

        # Alternative template set:
        # template_name = "Burgasser+2017"
        # templates = provider.get_Burgasser_2017_templates()

        # Mask bad values in the spectrum to improve comparison results
        data = mask_bad_values(data, mask_func=lambda mask: (mask % 2 == 1) | (mask >= 64))

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
            plot_format="png",
        )
    else:
        print("No spectrum found for the given object ID")
else:
    print("No object found for the given coordinates and search radius")
