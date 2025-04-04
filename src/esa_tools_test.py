import warnings
from astropy.utils.exceptions import AstropyWarning

from tools.esa_tools import retrieve_objects, retrieve_spectrum, retrieve_cutout, print_catalog_info
from tools.spectrum_plotter import plot_spectrum
from tools.image_plotter import plot_images
from tools.shared import MaskType, print_results


# ==============================
# E S A  tools test
# ==============================

warnings.simplefilter("ignore", category=AstropyWarning)

# ------------------------------
# Print catalog information
# ------------------------------
print_catalog_info()


# ------------------------------
# Retrieve objects
# ------------------------------

ra, dec = 266.4850113, 64.9936424
search_radius = 5  # arcsec

results = retrieve_objects(ra, dec, search_radius)
print_results(results)


# ------------------------------
# Retrieve spectrum
# ------------------------------

result = results[0]
object_id = str(result["object_id"])

table = retrieve_spectrum(object_id, maskType=MaskType.NONE)

if table and len(table) > 0:
    plot_spectrum(table, ra, dec)


# ------------------------------
# Retrieve cutouts
# ------------------------------

search_radius = 5  # arcsec
cutout_size = 20  # arcsec

vis_hdu = retrieve_cutout(ra, dec, search_radius, cutout_size, "VIS")
y_hdu = retrieve_cutout(ra, dec, search_radius, cutout_size, "Y")
j_hdu = retrieve_cutout(ra, dec, search_radius, cutout_size, "J")
h_hdu = retrieve_cutout(ra, dec, search_radius, cutout_size, "H")

images = []

images.append({"hdu": vis_hdu, "band": "VIS", "rgb": "b"})
images.append({"hdu": y_hdu, "band": "Y", "rgb": "g"})
images.append({"hdu": j_hdu, "band": "J", "rgb": None})
images.append({"hdu": h_hdu, "band": "H", "rgb": "r"})

plot_images(ra, dec, images, cutout_size, plot_format="pdf")
