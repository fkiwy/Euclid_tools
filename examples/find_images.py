import sys
import os

sys.path.append(os.path.abspath(os.path.join("..", "src")))

from image_plotter import plot_images
from image_finder import find_images


ra, dec = 273.891148, 64.526926
search_radius = 5  # arcsec
cutout_size = 20  # arcsec

vis_hdu, y_hdu, j_hdu, h_hdu = find_images(ra, dec, search_radius, cutout_size)

images = []

images.append({"hdu": vis_hdu, "band": "VIS", "rgb": "b"})
images.append({"hdu": y_hdu, "band": "Y", "rgb": "g"})
images.append({"hdu": j_hdu, "band": "J", "rgb": None})
images.append({"hdu": h_hdu, "band": "H", "rgb": "r"})

plot_images(ra, dec, images, cutout_size)
