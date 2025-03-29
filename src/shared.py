import os
import sys
import subprocess
import astropy.units as u
from astropy.coordinates import SkyCoord

irsa_url = 'https://irsa.ipac.caltech.edu/'

table_mer = 'euclid_q1_mer_catalogue'

table_1dspectra = 'euclid.objectid_spectrafile_association_q1'


def open_file(filename):
    if sys.platform == 'win32':
        os.startfile(filename)
    else:
        opener = 'open' if sys.platform == 'darwin' else 'evince'
        subprocess.call([opener, filename])


def create_object_name(ra, dec, precision=0, sep='', prefix=None, shortform=False, decimal=True):
    coords = SkyCoord(ra=ra * u.degree, dec=dec * u.degree)

    if shortform:
        coords_str = coords.to_string('hmsdms', decimal=False, sep=sep, precision=0)
        object_name = coords_str[0:4] + coords_str[7:12]
    else:
        if decimal:
            object_name = coords.to_string(decimal=True, precision=precision)
        else:
            object_name = coords.to_string('hmsdms', decimal=False, sep=sep, precision=precision).replace(' ', '')

    if prefix:
        object_name = prefix + object_name

    return str(object_name)
