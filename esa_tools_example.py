"""Retrieve catalog photometry, a red spectrum and VIS/Y/J/H cutouts from ESA.

Setup and execution (from this repository):
    python -m pip install -e .
    python esa_tools_example.py

Configuration
-------------
The default is release="dr1", dr1_catalogue="deep_survey", at ICRS coordinates
(RA, Dec) = (58.1332495, -49.1830038) degrees. This object is present in Q1 and
DR1; (59.7913643, -47.6826163) is an alternative below. Set release="q1" to use
Q1's catalog and product routes at the same coordinates. The library API itself
still defaults to Q1. For DR1 wide products, set dr1_catalogue="wide_survey".
Changing release/survey does not imply identical source IDs, observing visits,
processing or usable spectral coverage.

The script prints the selected catalog metadata, performs a 5-arcsecond cone
search, selects its closest result and prints an ID lookup from that same
catalog. Catalog fluxes are in microJy. DR1 adds profile AB magnitudes; optional
Vega output requires an explicitly sourced vega_calibration (see
docs/DR1_VEGA_CALIBRATION.md). Q1 retains its existing AB/Vega outputs.

For DR1, a red-grism spectrum is requested for the selected catalog ID. Multiple
associations require an explicit spectrum_product=(ESA parent path, HDU), chosen
from backend.discover_spectra(object_id, release="dr1", grism="red"). Parents
without a distinguishing native source-patch identifier remain unsupported.
The spectrum display is trimmed to 1.21-1.89 micrometres with uncertainty shading;
no scientific quality-mask selection or repair of missing bins is applied.

Images use 20-arcsecond cutouts and the corresponding DR1 wide/deep collection.
The selected catalog's native tile_index resolves tile overlaps; image_tile_id
can explicitly override it. Remaining product ambiguity requires inspection via
backend.discover_mosaics and an exact product selection with the library API.
Q1 uses its legacy image route.

Outputs and dependencies
------------------------
Spectrum and image-panel plots are PNG files in the system temporary directory
and open in the external viewer by default. Set output_dir, open_plot or
plot_format on the plotting calls below to change this. Use separate output_dir
values when comparing releases to retain both sets of target-named plots. These plots
use Matplotlib Agg, avoiding a Tk GUI. No named FITS output is requested;
retrieve_spectrum(save_spectrum=True, output_dir=...) can save a spectrum.
Legacy Q1 clients may use temporary download files. An empty cone exits with
the selected catalog and coordinates. Successful retrieval depends on service
availability. DR1 failures raise; the legacy ESA Q1 client can log a failure and
return no files.

Uses euclid_tools.get_backend("esa"), euclid_tools.esa_tools for metadata,
euclid_tools.shared, euclid_tools.image_plotter, Astropy and flux_comp.
"""


import sys
import warnings

# These examples export files and use an external viewer; never create Tk GUI
# objects that an S3/async worker could later garbage-collect off the main thread.
import matplotlib
matplotlib.use("Agg")

from astropy.utils.exceptions import AstropyWarning
from euclid_tools.esa_tools import print_catalog_info
from euclid_tools.image_plotter import plot_images
from euclid_tools.shared import print_results, print_message, create_object_name
from flux_comp.core import SED, WaveFlux

from euclid_tools import get_backend
from euclid_tools.releases import catalog_table

backend = get_backend('esa')

warnings.simplefilter("ignore", category=AstropyWarning)

# ==============================
# E S A  tools example
# ==============================

# ------------------------------
# Print catalog information
# ------------------------------
release = "dr1"  # Set "q1" to compare the same target with Q1.
dr1_catalogue = "deep_survey"  # Set "wide_survey" for the DR1 wide product.
catalogue = dr1_catalogue if release == "dr1" else None
vega_calibration = None  # DR1: load a completed, source-recorded calibration JSON to opt in.
spectrum_product = None  # DR1: exact (ESA parent path, HDU) from backend.discover_spectra.
image_tile_id = None  # DR1: use selected catalog tile_index; set a string to override.
print(f"ESA catalog: {catalog_table(release, catalogue)}")
print_catalog_info(release=release, catalogue=catalogue)

# ------------------------------
# Retrieve objects
# ------------------------------

# Shared Q1/DR1 target: coordinates remain fixed when changing release/survey.
ra, dec = 58.1332495, -49.1830038
# ra, dec = 59.7913643, -47.6826163

# Specify search radius in arcseconds
search_radius = 5

# Perform cone search for objects
results = backend.retrieve_objects(ra, dec, search_radius=search_radius, release=release, catalogue=catalogue, vega_calibration=vega_calibration)

if not results:
    # If no results are found, print message and exit
    print_message(ra, dec, search_radius, release=release, catalogue=catalog_table(release, catalogue))
    sys.exit()

# Print the retrieved object information
print_results(results)

# ------------------------------
# Retrieve object by ID
# ------------------------------

# Select the first object from the cone search results
result = results[0]
object_id = str(result["object_id"])

# Retrieve the catalog entry for the selected object
object_by_id = backend.retrieve_object_by_id(object_id, release=release, catalogue=catalogue, vega_calibration=vega_calibration)
if object_by_id is not None:
    print_results(object_by_id)

# ------------------------------
# Retrieve spectrum
# ------------------------------

# Retrieve the spectrum for the selected object
spectrum_options = {"release": release}
if release == "dr1":
    spectrum_options["grism"] = "red"
    spectrum_options["product"] = spectrum_product
data = backend.retrieve_spectrum(object_id, **spectrum_options)

# Plot the spectrum if data is available
if data and len(data) > 0:
    # Create a WaveFlux object for comparison
    spectrum = WaveFlux(
        label="Spectrum", wavelength=data["WAVELENGTH"], flux=data["FLUX"], uncertainty=data["ERROR"]
    )

    # Trim spectrum to remove excessive noise at the edges
    spectrum.trim(1.21, 1.89)

    # Create object name for plotting
    object_name = create_object_name(ra, dec, precision=2, shortform=False, prefix="J", decimal=False)

    # Set up the SED object
    sed = SED(object_name)
    sed.add(spectrum)

    # Convert to flux lambda and plot the spectrum
    sed.to_flux_lambda()
    sed.plot(
        spec_uncertainty=True,
        figure_size=(8, 5),
        show_grid=True,
        plot_format="png",
    )

# ------------------------------
# Retrieve cutouts
# ------------------------------

# Define cutout size (arcsec)
cutout_size = 20  # arcsec

# Retrieve cutouts from various bands (VIS, Y, J, H)
image_options = {"release": release}
if release == "dr1":
    image_options["tiling"] = "deep" if catalogue == "deep_survey" else "wide"
    image_options["selection"] = "deep" if catalogue == "deep_survey" else "survey"
    selected_image_tile = image_tile_id
    if selected_image_tile is None and "tile_index" in result.colnames:
        import numpy as np
        catalog_tile = result["tile_index"]
        if not np.ma.is_masked(catalog_tile):
            selected_image_tile = str(catalog_tile)
    if selected_image_tile is not None:
        image_options["tile_id"] = selected_image_tile
        print(f"ESA images: using tile {selected_image_tile}")
vis_hdu = backend.retrieve_cutout(ra, dec, search_radius=search_radius, cutout_size=cutout_size, band="VIS", **image_options)
y_hdu = backend.retrieve_cutout(ra, dec, search_radius=search_radius, cutout_size=cutout_size, band="Y", **image_options)
j_hdu = backend.retrieve_cutout(ra, dec, search_radius=search_radius, cutout_size=cutout_size, band="J", **image_options)
h_hdu = backend.retrieve_cutout(ra, dec, search_radius=search_radius, cutout_size=cutout_size, band="H", **image_options)

# Collect the images in a list with band and RGB mapping information
images = []
images.append({"hdu": vis_hdu, "band": "VIS", "rgb": "b"})
images.append({"hdu": y_hdu, "band": "Y", "rgb": "g"})
images.append({"hdu": j_hdu, "band": "J", "rgb": None})
images.append({"hdu": h_hdu, "band": "H", "rgb": "r"})

# Plot the images in a multi-band format and save the plot as a PNG
plot_images(ra, dec, images, cutout_size, plot_format="png")
