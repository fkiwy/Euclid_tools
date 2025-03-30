import pyvo as vo
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astropy import units as u
import shared as shr


def find_images(ra, dec, search_radius=5, cutout_size=20):

    def create_cutout(cutout_coords, cutout_size, image_url, band):
        print(f"Downloading {band} cutout")
        hdu = fits.open(image_url, use_fsspec=True)
        header = hdu[0].header
        cutout_data = Cutout2D(hdu[0].section, position=cutout_coords, size=cutout_size, wcs=WCS(hdu[0].header))
        hdu.close()
        new_hdu = fits.PrimaryHDU(data=cutout_data.data, header=header)
        new_hdu.header.update(cutout_data.wcs.to_header())
        return new_hdu

    def get_image_url(table, band):
        return table[table["energy_bandpassname"] == band]["access_url"][0]

    position = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame="icrs")
    search_radius *= u.arcsec
    cutout_size *= u.arcsec

    irsa_service = vo.dal.sia2.SIA2Service(shr.irsa_url + "SIA")

    table = irsa_service.search(pos=(position, search_radius), collection="euclid_DpdMerBksMosaic").to_table()

    result = table[(table["dataproduct_subtype"] == "science")]

    vis_hdu = create_cutout(position, cutout_size, get_image_url(result, "VIS"), "VIS")
    y_hdu = create_cutout(position, cutout_size, get_image_url(result, "Y"), "Y")
    j_hdu = create_cutout(position, cutout_size, get_image_url(result, "J"), "J")
    h_hdu = create_cutout(position, cutout_size, get_image_url(result, "H"), "H")

    return vis_hdu, y_hdu, j_hdu, h_hdu
