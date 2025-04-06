# Euclid Tools

A collection of tools to work with spectral and imaging data from the **ESA Euclid mission**, specifically tailored to handle data from the **Euclid Quick Release 1 (Q1)**. This package supports retrieval, visualization, and comparison of Euclid spectra and image cutouts, as well as matching them to template libraries.

The toolkit is written in Python and is designed to be lightweight, modular, and extensible. It includes high-level convenience functions for accessing and processing data, making it suitable for both exploratory data analysis and more in-depth scientific research.

> **Note:** The tools are designed to support both **ESA** and **IRSA** data access. While functionality is similar between the two, **ESA services are significantly faster** and are preferred for most use cases.

The tools allow users to:
- Retrieve catalog objects, spectra, and image cutouts.
- Visualize spectra and images (RGB composites).
- Compare observed spectra to predefined templates.

## Functionality Overview

This toolkit provides the following core functionality:

- **Object Search and Metadata Retrieval**
  - Search for Euclid sources near specified sky coordinates.
  - Retrieve source metadata from the Euclid MER (Merged External Reference) catalog.

> **Note:** The Euclid MER catalog does **not include magnitudes for Euclid bands**. To address this, the toolkit calculates magnitudes with uncertainties from the available flux columns using appropriate zero-points. Magnitudes are computed in both **AB** and **Vega** systems, and added directly to the result table returned by `retrieve_objects()`.
> - VIS magnitude is derived from `flux_vis_psf`
> - Y magnitude from `flux_y_templfit`
> - J magnitude from `flux_j_templfit`
> - H magnitude from `flux_h_templfit`

- **Spectrum Retrieval and Visualization**
  - Retrieve NISP spectra (slitless grism) for individual sources.
  - Automatically mask unreliable flux or error values.
  - Plot flux and uncertainty versus wavelength with publication-quality visuals.

- **Image Cutout Retrieval**
  - Access small imaging cutouts in VIS, Y, J, and H bands for any Euclid-detected source.
  - Plot single- and multi-band image views (with optional RGB assignments).

- **Spectral Comparison**
  - Match retrieved Euclid spectra to template libraries (e.g., Burgasser+2017, Theissen+2022).
  - Automatically evaluate fits using reduced χ² or other metrics.
  - Generate labeled and customizable comparison plots.

## Detailed Description

### Folder Structure:
The toolset consists of several scripts and modules:
- **`compare_spectrum_example.py`**: Example script to compare a spectrum with template SEDs.
- **`esa_tools_example.py`**: Example script for retrieving and plotting data from the ESA archive.
- **`irsa_tools_example.py`**: Example script for retrieving and plotting data from the IRSA archive.
- **`euclid_tools/esa_tools.py`**: Contains functions for interacting with the ESA archive.
- **`euclid_tools/irsa_tools.py`**: Contains functions for interacting with the IRSA archive.
- **`euclid_tools/image_plotter.py`**: Functions for plotting image cutouts.
- **`euclid_tools/spectrum_plotter.py`**: Functions for plotting spectra.
- **`euclid_tools/shared.py`**: Shared utility functions used across the toolset.
- **`flux_comp/core.py`**: Core functions for spectrum comparison.

### 1. `retrieve_objects()`

Located in `euclid_tools/esa_tools.py` or `euclid_tools/irsa_tools.py` depending on which version is used.

This function allows users to retrieve catalog objects from the Euclid archive based on coordinates and a specified search radius.

#### Parameters:
- `ra`: Right Ascension of the target object.
- `dec`: Declination of the target object.
- `search_radius`: Radius (in arcseconds) around the target to search for objects.

#### Returns:
- A list of objects that are found within the specified search radius.

### 2. `retrieve_spectrum()`

Located in `euclid_tools/esa_tools.py` or `euclid_tools/irsa_tools.py`.

This function retrieves the spectrum of a given object from the Euclid archive based on the object ID. It also allows users to mask bad data points (e.g., those with erroneous error values) to improve the comparison accuracy.

#### Parameters:
- `object_id`: The unique identifier of the object.
- `maskType`: A type of mask to apply to the spectrum data (e.g., `MaskType.ERROR`).

#### Returns:
- A table containing the spectrum data, including wavelength, flux, and error.

### 3. `retrieve_cutout()`

Located in `euclid_tools/esa_tools.py` or `euclid_tools/irsa_tools.py`.

This function retrieves an image cutout for a given region around the target object and specified band (e.g., VIS, Y, J, H).

#### Parameters:
- `ra`: Right Ascension of the target object.
- `dec`: Declination of the target object.
- `search_radius`: Radius (in arcseconds) to define the region of interest.
- `cutout_size`: Size (in arcseconds) of the cutout.
- `band`: The band (e.g., "VIS", "Y", "J", "H") to retrieve.

#### Returns:
- A data object (HDUs) containing the image cutout for the specified band.

### 4. `plot_spectrum()`

Located in `euclid_tools/spectrum_plotter.py`.

This function generates a plot of the observed spectrum for a given object. The plot includes both the flux and error values, allowing users to visualize the spectrum and uncertainties.

#### Parameters:
- `data`: The spectrum data (including wavelength, flux, and error).
- `ra`: Right Ascension of the target object.
- `dec`: Declination of the target object.
- `output_dir`: Directory where the plot will be saved.
- `open_plot`: Whether to open the plot file after saving.
- `plot_format`: The format of the plot (e.g., "pdf", "png").

#### Returns:
- A saved plot of the spectrum.

### 5. `plot_images()`

Located in `euclid_tools/image_plotter.py`.

This function generates plots for multiple images (cutouts) of the object in different bands (e.g., VIS, Y, J, H). The images can be combined into RGB composites.

#### Parameters:
- `ra`: Right Ascension of the target object.
- `dec`: Declination of the target object.
- `images`: A list of dictionaries containing image data for each band.
- `img_size`: The size of the image cutout (in arcseconds).
- `image_contrast`: Contrast level to adjust the image brightness.
- `output_dir`: Directory where the plot will be saved.
- `open_plot`: Whether to open the plot file after saving.
- `plot_format`: The format of the plot (e.g., "pdf", "png").

#### Returns:
- A saved plot of the RGB images.

### 6. Spectrum Comparison Using `flux_comp`

Located in `flux_comp/core.py`.

To compare an observed spectrum to a template, the `flux_comp` package is used. The script leverages the `SED` and `WaveFlux` classes for the comparison, employing metrics such as reduced chi-squared.

#### Key steps:
1. Retrieve the catalog object based on RA, Dec.
2. Retrieve the spectrum for the object.
3. Retrieve template(s) for comparison (e.g., from "Theissen+2022").
4. Trim the spectrum to a specific wavelength range.
5. Compare the spectrum to the templates using the reduced chi-squared metric.
6. Plot the comparison.

## Example Usages

### Example 1: Retrieving and Plotting a Spectrum

This example demonstrates how to retrieve a spectrum for an object and plot it using the `plot_spectrum` function. The full example code can be found in `esa_tools_example.py`.

```python
from euclid_tools.esa_tools import retrieve_objects, retrieve_spectrum
from euclid_tools.spectrum_plotter import plot_spectrum

# Define object coordinates
ra, dec = 266.4850113, 64.9936424

# Search radius (arcsec)
search_radius = 5

# Retrieve object data from ESA
results = retrieve_objects(ra, dec, search_radius)

if results:
    result = results[0]
    object_id = str(result["object_id"])

    # Retrieve spectrum data
    spectrum_data = retrieve_spectrum(object_id)

    # Plot the spectrum
    plot_spectrum(spectrum_data, ra, dec, output_dir=".", plot_format="pdf")
else:
    print("No object found for the given coordinates and radius.")
```
![Spectrum](examples/J174556.40+645937.11_spectrum.png)

### Example 2: Retrieving and Plotting Image Cutouts

This example shows how to retrieve image cutouts for a given object using different bands, and plot them using the `plot_images` function. The full example code can be found in `esa_tools_example.py`.

```python
from euclid_tools.esa_tools import retrieve_cutout
from euclid_tools.image_plotter import plot_images

# Define object coordinates and cutout size (arcsec)
ra, dec = 266.4850113, 64.9936424
cutout_size = 20  # arcsec

# Retrieve image cutouts for different bands
search_radius = 5  # arcsec
vis_hdu = retrieve_cutout(ra, dec, search_radius, cutout_size, "VIS")
y_hdu = retrieve_cutout(ra, dec, search_radius, cutout_size, "Y")
j_hdu = retrieve_cutout(ra, dec, search_radius, cutout_size, "J")
h_hdu = retrieve_cutout(ra, dec, search_radius, cutout_size, "H")

# Prepare the images list for plotting
images = [
    {"hdu": vis_hdu, "band": "VIS", "rgb": "b"},
    {"hdu": y_hdu, "band": "Y", "rgb": "g"},
    {"hdu": j_hdu, "band": "J", "rgb": None},
    {"hdu": h_hdu, "band": "H", "rgb": "r"},
]

# Plot the images in a single figure
plot_images(ra, dec, images, cutout_size, plot_format="pdf")
```
![Images](examples/J174556.40+645937.11_images.png)

### Example 3: Compare Spectrum to Template Using `flux_comp`

This example demonstrates how to use the Euclid tools to retrieve a spectrum for a given object
and compare it to predefined templates using the `flux_comp` tool. The full example code can be found in `compare_spectrum_example.py`.

```python
from flux_comp.core import SED, WaveFlux, TemplateProvider
from euclid_tools.esa_tools import retrieve_spectrum
from euclid_tools.shared import MaskType

# Retrieve spectrum for an object
data = retrieve_spectrum(object_id, maskType=MaskType.ERROR)

# Retrieve templates for comparison
template_name = "Theissen+2022"
provider = TemplateProvider()
templates = provider.get_Theissen_2022_templates(smooth_window=10)

# Create a WaveFlux object and trim the spectrum
spectrum = WaveFlux(
    label="Spectrum", wavelength=data["WAVELENGTH"], flux=data["FLUX"], uncertainty=data["ERROR"]
)
spectrum.trim(1.22, 1.88)  # Trim to a specific wavelength range

# Compare the spectrum to the templates
sed = SED("Object vs. " + template_name)
sed.compare(spectrum, templates, trim_wave=True, metric="reduced-chi2")
sed.to_flux_lambda()
sed.plot(reference_on_top=False, spec_uncertainty=True)
```
![Comparison](examples/J174556.40+645937.11_vs._Theissen+2022.png)
![Comparison](examples/J174556.40+645937.11_vs._Burgasser+2017.png)

## Installation

To install the Euclid Tools package, clone the repository and install the dependencies:

```bash
git clone https://github.com/fkiwy/Euclid_tools.git
cd Euclid_tools
pip install .
```
