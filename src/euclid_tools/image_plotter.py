import os
import tempfile
import warnings

import PIL.ImageOps
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.utils.exceptions import AstropyWarning
from astropy.visualization import make_lupton_rgb
from matplotlib.patches import BoxStyle, Rectangle
from reproject import reproject_interp
from reproject.mosaicking import find_optimal_celestial_wcs

import euclid_tools.shared as shr


def plot_images(
    ra, dec, images, img_size, image_contrast=10, output_dir=tempfile.gettempdir(), open_plot=True, plot_format="png"
):
    """
    Generate and save a panel of individual and RGB composite images centered on the specified coordinates.

    This function takes a list of Euclid (or other) image HDUs, reprojects them to a common WCS,
    generates grayscale band images and an RGB color composite, and annotates each with position markers
    and band labels. It also includes metadata about the target (RA/Dec, cutout size, etc.)
    in the final panel.

    Parameters:
        ra (float): Right Ascension of the target in degrees.
        dec (float): Declination of the target in degrees.
        images (list): List of dictionaries, each containing:
            - "hdu": FITS HDU with image data and WCS,
            - "band": String representing the band name (e.g., 'VIS', 'Y', 'J'),
            - "rgb": String flag indicating whether the band corresponds to R, G, or B channel (e.g., 'r', 'g', 'b').
        img_size (float): Desired cutout size in arcseconds.
        image_contrast (float, optional): Percentile range used for display scaling (default is 10).
        output_dir (str, optional): Directory to save the output image file (default: system temp directory).
        open_plot (bool, optional): Whether to open the plot after saving (default: True).
        plot_format (str, optional): Format to save the image file in (e.g., 'png', 'pdf').

    Output:
        Saves a figure in the specified format showing:
            - A grid of grayscale band images with band labels.
            - A composite RGB image based on selected bands.
            - A panel with object information (coordinates, size, orientation).
        The filename is constructed from the RA/Dec coordinates (e.g., 'JHHMMSS+DDMMSS_images.png').

    Notes:
        - The function uses astropy's WCS tools to reproject and align the images.
        - The RGB image is built using Lupton scaling via `make_lupton_rgb`.
        - Images are displayed in grayscale ('gray_r') with dynamic scaling based on median and MAD.
        - The composite image assumes all three RGB bands are provided.
    """

    def create_image(hdu, img_idx, band):
        wcs, shape = find_optimal_celestial_wcs([hdu], frame="icrs")
        data, _ = reproject_interp(hdu, wcs, shape_out=shape)
        position = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame="icrs")

        ax = fig.add_subplot(rows, cols, img_idx, projection=wcs)
        x, y = wcs.world_to_pixel(position)
        ax.plot(x, y, "ro", fillstyle="none", markersize=12, markeredgewidth=0.4)
        ax.plot(x, y, "ro", fillstyle="none", markersize=0.4, markeredgewidth=0.4)
        ax.text(
            0.03,
            0.93,
            band,
            color="black",
            fontsize=3.0,
            transform=ax.transAxes,
            bbox={"facecolor": "white", "alpha": 0.5, "linewidth": 0.1, "boxstyle": BoxStyle("Square", pad=0.3)},
        )
        ax.add_patch(Rectangle((0, 0), 1, 1, fill=False, lw=0.2, ec="black", transform=ax.transAxes))

        vmin, vmax = get_min_max(data)
        ax.imshow(data, vmin=vmin, vmax=vmax, cmap="gray_r")
        ax.axis("off")
        return data, x, y

    def create_color_image(r, g, b, img_idx, band, x, y, invert=False):
        r = Image.fromarray(get_rgb(r)).convert("L")
        g = Image.fromarray(get_rgb(g)).convert("L")
        b = Image.fromarray(get_rgb(b)).convert("L")

        if invert:
            r = PIL.ImageOps.invert(r)
            g = PIL.ImageOps.invert(g)
            b = PIL.ImageOps.invert(b)

        rgb = Image.merge("RGB", (r, g, b))

        ax = fig.add_subplot(rows, cols, img_idx)
        ax.plot(x, y, "ro", fillstyle="none", markersize=12, markeredgewidth=0.4)
        ax.plot(x, y, "ro", fillstyle="none", markersize=0.4, markeredgewidth=0.4)
        ax.text(
            0.03,
            0.93,
            band,
            color="black",
            fontsize=3.0,
            transform=ax.transAxes,
            bbox={"facecolor": "white", "alpha": 0.7, "linewidth": 0.1, "boxstyle": BoxStyle("Square", pad=0.3)},
        )
        ax.add_patch(Rectangle((0, 0), 1, 1, fill=False, lw=0.2, ec="black", transform=ax.transAxes))
        ax.imshow(rgb, origin="lower")
        ax.axis("off")

    def get_rgb(data):
        vmin, vmax = get_min_max(data)
        return make_lupton_rgb(data, data, data, minimum=vmin, stretch=vmax - vmin, Q=0)

    def get_min_max(data):
        lo = image_contrast
        hi = 100 - image_contrast
        med = np.nanmedian(data)
        mad = np.nanmedian(abs(data - med))
        dev = np.nanpercentile(data, hi) - np.nanpercentile(data, lo)
        vmin = med - 2.0 * mad
        vmax = med + 2.0 * dev
        return vmin, vmax

    warnings.simplefilter("ignore", category=AstropyWarning)

    img_idx = 1
    fontsize = 4.5
    rows, cols = 5, 6

    fig = plt.figure()
    fig.set_figheight(5)
    fig.set_figwidth(5)
    plt.subplots_adjust(wspace=0, hspace=0.1, right=1.0)
    plt.rcParams.update({"font.family": "Arial"})

    for image in images:
        hdu = image["hdu"]
        band = image["band"]
        rgb = image["rgb"]

        print(f"Creating {band}-band image")
        data, x, y = create_image(hdu, img_idx, band)

        if rgb == "r":
            r_data = data
            r_band = band
        if rgb == "g":
            g_data = data
            g_band = band
        if rgb == "b":
            b_data = data
            b_band = band
        img_idx += 1

    print("Creating color image")
    create_color_image(r_data, g_data, b_data, img_idx, f"{r_band}-{g_band}-{b_band}", x, y)
    img_idx += 1

    object_name = shr.create_object_name(ra, dec, precision=2, shortform=False, prefix="J", decimal=False)
    filename = os.path.join(output_dir, object_name + "_images." + plot_format)

    ax = fig.add_subplot(rows, cols, img_idx)
    ax.text(0.1, 0.75, object_name, fontsize=fontsize, transform=ax.transAxes)
    ax.text(0.1, 0.60, r"$\alpha$ = " + str(round(ra, 7)), fontsize=fontsize, transform=ax.transAxes)
    ax.text(0.1, 0.45, r"$\delta$ = " + str(round(dec, 7)), fontsize=fontsize, transform=ax.transAxes)
    ax.text(0.1, 0.30, "Size = " + str(int(img_size)) + " arcsec", fontsize=fontsize, transform=ax.transAxes)
    ax.text(0.1, 0.15, "North up, East left", fontsize=fontsize, transform=ax.transAxes)
    ax.axis("off")
    img_idx += 1

    plt.savefig(filename, dpi=600, bbox_inches="tight", format=plot_format)
    plt.close()

    if open_plot:
        shr.open_file(filename)
