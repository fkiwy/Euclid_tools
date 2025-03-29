import pyvo as vo
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astropy import units as u
from image_creator import create_images


def create_cutout(cutout_coords, cutout_size, filename):
    print(f'Opening {filename}')
    hdu = fits.open(filename, use_fsspec=True)
    header = hdu[0].header
    cutout_data = Cutout2D(
        hdu[0].section, position=cutout_coords, size=cutout_size, wcs=WCS(hdu[0].header))
    hdu.close()
    new_hdu = fits.PrimaryHDU(data=cutout_data.data, header=header)
    new_hdu.header.update(cutout_data.wcs.to_header())
    return new_hdu


def extract_file_info(df, band):
    filename = df[df['energy_bandpassname'] == band]['access_url'].to_list()[0]
    filesize = df[df['energy_bandpassname'] == band]['access_estsize'].to_list()[0]/1000000
    return filename, filesize


pd.options.mode.copy_on_write = True
pd.set_option('display.max_columns', None)
pd.set_option('display.max_colwidth', None)
# pd.reset_option('display.max_columns')
# pd.reset_option('display.max_colwidth')

ra, dec = 273.891148, 64.526926
coord = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='icrs')

search_radius = 5 * u.arcsec
cutout_size = 20 * u.arcsec

b_band, g_band, r_band = 'VIS', 'J', 'H'

irsa_service = vo.dal.sia2.SIA2Service('https://irsa.ipac.caltech.edu/SIA')

image_table = irsa_service.search(
    pos=(coord, search_radius), collection='euclid_DpdMerBksMosaic')

df_im_irsa = image_table.to_table().to_pandas()

df_im_euclid = df_im_irsa[(df_im_irsa['dataproduct_subtype'] == 'science')]

# print('Number of MER images/MER tile', len(df_im_euclid))

b_filename, b_filesize = extract_file_info(df_im_euclid, 'VIS')
# print('File name VIS:', b_filename)
# print('File size VIS:', b_filesize)

g_filename, g_filesize = extract_file_info(df_im_euclid, 'Y')
# print('File name Y:', g_filename)
# print('File size Y:', g_filesize)

r_filename, r_filesize = extract_file_info(df_im_euclid, 'H')
# print('File name H:', r_filename)
# print('File size H:', r_filesize)

b_image = create_cutout(coord, cutout_size, b_filename)
g_image = create_cutout(coord, cutout_size, g_filename)
r_image = create_cutout(coord, cutout_size, r_filename)

create_images(coord, b_image, g_image, r_image, b_band, g_band, r_band, cutout_size.value, image_contrast=5)
