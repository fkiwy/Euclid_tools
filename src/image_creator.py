import os
import sys
import tempfile
import warnings
import subprocess
import numpy as np
import astropy.units as u
# from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
from astropy.visualization import make_lupton_rgb
from astropy.utils.exceptions import AstropyWarning
from reproject import reproject_interp
from reproject.mosaicking import find_optimal_celestial_wcs
import matplotlib.pyplot as plt
from matplotlib.patches import BoxStyle, Rectangle
from PIL import Image
import PIL.ImageOps


def create_images(position, b_image, g_image, r_image, b_band, g_band, r_band, img_size,
                  image_contrast=10, output_dir=tempfile.gettempdir()):

    def create_image(hdu, img_idx, band):
        wcs, shape = find_optimal_celestial_wcs([hdu], frame='icrs')
        data, _ = reproject_interp(hdu, wcs, shape_out=shape)

        # position = SkyCoord(ra*u.deg, dec*u.deg)
        # cutout = Cutout2D(hdu.data, position, img_size*u.arcsec, wcs=wcs, mode='partial')
        # data = cutout.data
        # wcs = cutout.wcs

        ax = fig.add_subplot(rows, cols, img_idx, projection=wcs)
        x, y = wcs.world_to_pixel(position)
        ax.plot(x, y, 'ro', fillstyle='none', markersize=7, markeredgewidth=0.2)
        ax.plot(x, y, 'ro', fillstyle='none', markersize=0.2, markeredgewidth=0.2)
        ax.text(0.04, 0.91, band, color='black', fontsize=1.8, transform=ax.transAxes,
                bbox=dict(facecolor='white', alpha=0.5, linewidth=0.1, boxstyle=BoxStyle('Square', pad=0.3)))
        ax.add_patch(Rectangle((0, 0), 1, 1, fill=False, lw=0.2, ec='black', transform=ax.transAxes))

        vmin, vmax = get_min_max(data)
        ax.imshow(data, vmin=vmin, vmax=vmax, cmap='gray_r')
        ax.axis('off')
        return data, x, y

    def create_color_image(r, g, b, img_idx, band, x, y, invert=False):
        r = Image.fromarray(create_lupton_rgb(r)).convert('L')
        g = Image.fromarray(create_lupton_rgb(g)).convert('L')
        b = Image.fromarray(create_lupton_rgb(b)).convert('L')

        if invert:
            r = PIL.ImageOps.invert(r)
            g = PIL.ImageOps.invert(g)
            b = PIL.ImageOps.invert(b)

        rgb = Image.merge('RGB', (r, g, b))

        ax = fig.add_subplot(rows, cols, img_idx)
        ax.plot(x, y, 'ro', fillstyle='none', markersize=7, markeredgewidth=0.2)
        ax.plot(x, y, 'ro', fillstyle='none', markersize=0.2, markeredgewidth=0.2)
        ax.text(0.04, 0.91, band, color='black', fontsize=1.8, transform=ax.transAxes,
                bbox=dict(facecolor='white', alpha=0.7, linewidth=0.1, boxstyle=BoxStyle('Square', pad=0.3)))
        ax.add_patch(Rectangle((0, 0), 1, 1, fill=False, lw=0.2, ec='black', transform=ax.transAxes))
        ax.imshow(rgb, origin='lower')
        ax.axis('off')

    def create_lupton_rgb(data):
        vmin, vmax = get_min_max(data)
        return make_lupton_rgb(data, data, data, minimum=vmin, stretch=vmax-vmin, Q=0)

    def get_min_max(data):
        lo = image_contrast
        hi = 100-image_contrast
        med = np.nanmedian(data)
        mad = np.nanmedian(abs(data - med))
        dev = np.nanpercentile(data, hi) - np.nanpercentile(data, lo)
        vmin = med - 2.0 * mad
        vmax = med + 2.0 * dev
        return vmin, vmax

    def create_obj_name(ra, dec, precision=0, sep='', prefix=None, shortform=False, decimal=True):
        coords = SkyCoord(ra=ra * u.degree, dec=dec * u.degree)

        if shortform:
            coords_str = coords.to_string('hmsdms', decimal=False, sep=sep, precision=0)
            obj_name = coords_str[0:4] + coords_str[7:12]
        else:
            if decimal:
                obj_name = coords.to_string(decimal=True, precision=precision)
            else:
                obj_name = coords.to_string('hmsdms', decimal=False, sep=sep, precision=precision).replace(' ', '')

        if prefix:
            obj_name = prefix + obj_name

        return str(obj_name)

    def start_file(filename):
        if sys.platform == 'win32':
            os.startfile(filename)
        else:
            opener = 'open' if sys.platform == 'darwin' else 'evince'
            subprocess.call([opener, filename])

    def get_vega_photometry(catalog_entry, aper_size):
        vegamag = catalog_entry[aper_size + '_vegamag']
        vegamag_err = catalog_entry[aper_size + '_vegamag_err']
        return vegamag, vegamag_err

    def get_ab_photometry(catalog_entry, aper_size):
        abmag = catalog_entry[aper_size + '_abmag']
        abmag_err = catalog_entry[aper_size + '_abmag_err']
        return abmag, abmag_err

    # --------------------------------------
    # Code for create_finder_charts function
    # --------------------------------------
    warnings.simplefilter('ignore', category=AstropyWarning)

    rows = 10
    cols = 5
    img_idx = 1
    fig = plt.figure()
    fig.set_figheight(5)
    fig.set_figwidth(5)
    plt.subplots_adjust(wspace=0, hspace=0.05, right=0.5)

    r = g = b = None

    ra = position.ra.deg
    dec = position.dec.deg

    print(f'Object {ra} {dec}')

    print(f'  Creating {b_band}-band image')
    b, x, y = create_image(b_image, img_idx, b_band)
    img_idx += 1

    print(f'  Creating {g_band}-band image')
    g, x, y = create_image(g_image, img_idx, g_band)
    img_idx += 1

    print(f'  Creating {r_band}-band image')
    r, x, y = create_image(r_image, img_idx, r_band)
    img_idx += 1

    print('  Creating color image')
    create_color_image(r, g, b, img_idx, f'{r_band}-{g_band}-{b_band}', x, y)
    img_idx += 1

    # Astrometry
    fontsize = 2.4
    ax = fig.add_subplot(rows, cols, img_idx)
    obj_name = create_obj_name(ra, dec, precision=2, shortform=False, prefix='J', decimal=False)
    ax.text(0.1, 0.75, obj_name, fontsize=fontsize, transform=ax.transAxes)
    ax.text(0.1, 0.60, r'$\alpha$ = ' + str(round(ra, 7)), fontsize=fontsize, transform=ax.transAxes)
    ax.text(0.1, 0.45, r'$\delta$ = ' + str(round(dec, 7)), fontsize=fontsize, transform=ax.transAxes)
    ax.text(0.1, 0.30, 'Size = ' + str(int(img_size)) + ' arcsec', fontsize=fontsize, transform=ax.transAxes)
    ax.text(0.1, 0.15, 'North up, East left', fontsize=fontsize, transform=ax.transAxes)
    ax.axis('off')
    img_idx += 1

    # Save and open the file
    filename = os.path.join(output_dir, obj_name + '.pdf')
    plt.savefig(filename, dpi=2000, bbox_inches='tight', format='pdf')
    plt.close()

    start_file(filename)
