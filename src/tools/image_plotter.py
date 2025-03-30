import os
import tempfile
import warnings
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord

# from astropy.nddata import Cutout2D
from astropy.visualization import make_lupton_rgb
from astropy.utils.exceptions import AstropyWarning
from reproject import reproject_interp
from reproject.mosaicking import find_optimal_celestial_wcs
import matplotlib.pyplot as plt
from matplotlib.patches import BoxStyle, Rectangle
from PIL import Image
import PIL.ImageOps

import tools.shared as shr


def plot_images(
    ra, dec, images, img_size, image_contrast=10, output_dir=tempfile.gettempdir(), open_plot=True, plot_format="png"
):

    def create_image(hdu, img_idx, band):
        wcs, shape = find_optimal_celestial_wcs([hdu], frame="icrs")
        data, _ = reproject_interp(hdu, wcs, shape_out=shape)

        position = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame="icrs")
        # cutout = Cutout2D(hdu.data, position, img_size*u.arcsec, wcs=wcs, mode='partial')
        # data = cutout.data
        # wcs = cutout.wcs

        ax = fig.add_subplot(rows, cols, img_idx, projection=wcs)
        x, y = wcs.world_to_pixel(position)
        ax.plot(x, y, "ro", fillstyle="none", markersize=7, markeredgewidth=0.2)
        ax.plot(x, y, "ro", fillstyle="none", markersize=0.2, markeredgewidth=0.2)
        ax.text(
            0.04,
            0.91,
            band,
            color="black",
            fontsize=1.8,
            transform=ax.transAxes,
            bbox=dict(facecolor="white", alpha=0.5, linewidth=0.1, boxstyle=BoxStyle("Square", pad=0.3)),
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
        ax.plot(x, y, "ro", fillstyle="none", markersize=7, markeredgewidth=0.2)
        ax.plot(x, y, "ro", fillstyle="none", markersize=0.2, markeredgewidth=0.2)
        ax.text(
            0.04,
            0.91,
            band,
            color="black",
            fontsize=1.8,
            transform=ax.transAxes,
            bbox=dict(facecolor="white", alpha=0.7, linewidth=0.1, boxstyle=BoxStyle("Square", pad=0.3)),
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
    fontsize = 2.4
    rows, cols = 10, 5

    fig = plt.figure()
    fig.set_figheight(5)
    fig.set_figwidth(5)
    plt.subplots_adjust(wspace=0, hspace=0.05, right=0.5)

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
    filename = os.path.join(output_dir, object_name + "." + plot_format)

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
