import warnings
import requests
import pyvo as vo
from io import BytesIO
from astropy.io import fits
from astropy.table import Table
from astropy.utils.exceptions import AstropyWarning
import shared as shr


def find_spectrum(object_id):
    warnings.simplefilter("ignore", category=AstropyWarning)

    adql = f"""
    SELECT *
      FROM {shr.table_1dspectra}
     WHERE objectid = {object_id}
       AND uri IS NOT NULL
    """

    service = vo.dal.TAPService(shr.irsa_url + "TAP")
    result = service.search(adql)
    table = result.to_table()

    if len(table) == 0:
        return None

    file_url = shr.irsa_url + table["uri"][0]
    response = requests.get(file_url)

    with fits.open(BytesIO(response.content), memmap=True) as hdul:
        hdu = hdul[table["hdu"][0]]
        data = Table.read(hdu, format="fits", hdu=1)

    return data
