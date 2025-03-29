import os
import tempfile
import requests
import warnings
import pyvo as vo
from io import BytesIO
from astropy.io import fits
from astropy.table import Table
from astropy.utils.exceptions import AstropyWarning
import matplotlib.pyplot as plt
import shared as shr


def find_spectrum(object_id, plot=False, ra=None, dec=None):
    warnings.simplefilter('ignore', category=AstropyWarning)

    adql = f"""
    SELECT *
      FROM {shr.table_1dspectra}
     WHERE objectid = {object_id}
       AND uri IS NOT NULL
    """

    service = vo.dal.TAPService(shr.irsa_url + 'TAP')
    result = service.search(adql)
    table = result.to_table()

    if len(table) == 0:
        return None

    file_url = shr.irsa_url + table['uri'][0]
    response = requests.get(file_url)

    with fits.open(BytesIO(response.content), memmap=True) as hdul:
        hdu = hdul[table['hdu'][0]]
        data = Table.read(hdu, format='fits', hdu=1)

    if plot:
        output_dir = tempfile.gettempdir()

        if ra and dec:
            object_name = shr.create_object_name(ra, dec, precision=2, shortform=False, prefix='J', decimal=False)
        else:
            object_name = object_id

        plt.rcParams.update({'font.family': 'Arial'})

        # Converting from Angstrom to microns
        plt.plot(data['WAVELENGTH']/10000., data['SIGNAL'])

        plt.xlabel('Wavelength [microns]')
        plt.ylabel('Flux [' + data['SIGNAL'].unit.to_string('latex_inline') + ']')
        plt.title(object_name)

        filename = os.path.join(output_dir, object_name + '.pdf')
        plt.savefig(filename, dpi=300, bbox_inches='tight', format='pdf')
        plt.close()

        shr.open_file(filename)

    return data
