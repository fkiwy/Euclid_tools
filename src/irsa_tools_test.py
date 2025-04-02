import warnings
from astropy.table import Table
from astropy.utils.exceptions import AstropyWarning

from tools.irsa_tools import retrieve_objects, retrieve_spectrum, retrieve_cutout, print_catalog_info
from tools.spectrum_plotter import plot_spectrum
from tools.image_plotter import plot_images
import tools.shared as shr


# ==============================
# I R S A tools test
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
shr.print_results(results)


# ------------------------------
# Retrieve spectrum
# ------------------------------

result = results[0]
object_id = str(result["object_id"])

hdul = retrieve_spectrum(object_id)

if hdul and len(hdul) > 0:
    data = Table.read(hdul[0], format="fits")
    plot_spectrum(data, ra, dec)


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
