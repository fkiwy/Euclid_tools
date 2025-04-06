"""
Example: Using Euclid Tools for Object Retrieval, Spectrum Plotting, and Image Cutouts

This example demonstrates how to use tools from `tools.esa_tools` to retrieve catalog objects,
spectra, and image cutouts from the Euclid mission. It also shows how to plot the spectrum and 
the images in various bands.

The steps in the script are as follows:
1. Print catalog information.
2. Retrieve catalog objects near specified coordinates.
3. Retrieve and plot the spectrum for a selected object.
4. Retrieve image cutouts in different bands (VIS, Y, J, H).
5. Plot the retrieved images in a multi-band format.

Dependencies:
- `euclid_tools.esa_tools` (for retrieving catalog objects, spectra, and image cutouts)
- `euclid_tools.spectrum_plotter` (for spectrum plotting)
- `euclid_tools.image_plotter` (for plotting images)
- `euclid_tools.shared` (for utility functions like `MaskType`, `print_results`, and `print_message`)
- `astropy` (for units and warnings handling)

Steps in the script:
1. **Print Catalog Information:**
   - The script begins by printing information about the catalog available for querying, using the `print_catalog_info()` function.

2. **Retrieve Objects:**
   - It performs a cone search around the specified RA and Dec with a specified search radius (in arcseconds) using the `retrieve_objects()` function.
   - If objects are found, the first object is selected for further analysis. If no objects are found, the script will exit gracefully with a message.

3. **Retrieve Spectrum:**
   - The script retrieves the spectrum for the selected object using the `retrieve_spectrum()` function, with no masking applied (`MaskType.NONE`).
   - If a valid spectrum is retrieved, it will be plotted using the `plot_spectrum()` function.

4. **Retrieve Cutouts:**
   - The script retrieves image cutouts from different bands (VIS, Y, J, H) with a specified cutout size, using the `retrieve_cutout()` function.
   - Each retrieved cutout is appended to the `images` list with relevant metadata (band and RGB color mapping).

5. **Plot Images:**
   - The `plot_images()` function is called to create a multi-band image plot using the retrieved cutouts. The plot is saved as a PDF.

Notes for Users:
- **Coordinates:** The script starts with a specific set of coordinates (RA, Dec). Users can modify these to retrieve data for different objects.
- **Cutout Sizes:** The cutout size is defined in arcseconds, and users can adjust it based on their needs. For example, a `cutout_size` of 20 arcseconds is used in this script.
- **Plot Format:** The resulting plot for images is saved in PDF format. Users can change this format (e.g., PNG, JPG) by modifying the `plot_format` argument.
- **Masking Spectra:** No masking is applied to the spectrum retrieval in this example (`MaskType.NONE`). Users can apply other masking types, such as `MaskType.ERROR`, to filter poor-quality data.
- **Image Bands:** The bands (VIS, Y, J, H) are selected for cutout retrieval, but users can adjust this list to include other bands if needed.
  
Example Output:
The script will generate a plot that shows the spectrum of the selected object, as well as a multi-band image plot in the VIS, Y, J, and H bands. The plot will be saved in the specified format (PDF by default).
"""

import sys
import warnings
from astropy.utils.exceptions import AstropyWarning

from euclid_tools.esa_tools import retrieve_objects, retrieve_spectrum, retrieve_cutout, print_catalog_info
from euclid_tools.spectrum_plotter import plot_spectrum
from euclid_tools.image_plotter import plot_images
from euclid_tools.shared import MaskType, print_results, print_message

warnings.simplefilter("ignore", category=AstropyWarning)


# ==============================
# E S A  tools example
# ==============================

# ------------------------------
# Print catalog information
# ------------------------------
print_catalog_info()

# ------------------------------
# Retrieve objects
# ------------------------------

# Specify coordinates (RA, Dec) and search radius (arcsec)
ra, dec = 266.4850113, 64.9936424
search_radius = 5  # arcsec

# Perform cone search for objects
results = retrieve_objects(ra, dec, search_radius)

if not results:
    # If no results are found, print message and exit
    print_message(ra, dec, search_radius)
    sys.exit()

# Print the retrieved object information
print_results(results)

# ------------------------------
# Retrieve spectrum
# ------------------------------

# Select the first object from the results
result = results[0]
object_id = str(result["object_id"])

# Retrieve the spectrum for the selected object
table = retrieve_spectrum(object_id, maskType=MaskType.NONE)

if table and len(table) > 0:
    # Plot the spectrum if data is available
    plot_spectrum(table, ra, dec, plot_format="pdf")

# ------------------------------
# Retrieve cutouts
# ------------------------------

# Define cutout size (arcsec)
cutout_size = 20  # arcsec

# Retrieve cutouts from various bands (VIS, Y, J, H)
vis_hdu = retrieve_cutout(ra, dec, search_radius, cutout_size, "VIS")
y_hdu = retrieve_cutout(ra, dec, search_radius, cutout_size, "Y")
j_hdu = retrieve_cutout(ra, dec, search_radius, cutout_size, "J")
h_hdu = retrieve_cutout(ra, dec, search_radius, cutout_size, "H")

# Collect the images in a list with band and RGB mapping information
images = []
images.append({"hdu": vis_hdu, "band": "VIS", "rgb": "b"})
images.append({"hdu": y_hdu, "band": "Y", "rgb": "g"})
images.append({"hdu": j_hdu, "band": "J", "rgb": None})
images.append({"hdu": h_hdu, "band": "H", "rgb": "r"})

# Plot the images in a multi-band format and save the plot as a PDF
plot_images(ra, dec, images, cutout_size, plot_format="pdf")
