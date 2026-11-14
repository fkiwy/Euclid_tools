"""Verified ESA DR1 DataLink spectra and server-side MER cutouts.

Use public IVOA endpoints directly to avoid Astroquery's changing source_id/ids
signatures. All transfers are bounded; no Q1 or full-parent-product fallback.
"""
from __future__ import annotations

import hashlib
import csv
import io
import re
import warnings
import zipfile
from pathlib import PurePosixPath
from urllib.parse import urlencode

import astropy.units as u
import numpy as np
import requests
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import QTable, Table
from astropy.wcs import WCS

from euclid_tools.provenance import stamp
from euclid_tools.irsa_spectra import _object_id, _hdu, FLUX_UNIT
from euclid_tools.irsa_images import _metadata

BASE = 'https://eas.esac.esa.int'
TAP_URL = BASE + '/tap-server/tap/sync'
SPECTRUM_URL = BASE + '/sas-dd/data'
CUTOUT_URL = BASE + '/sas-cutout/cutout'
MAX_ROWS = 100
MAX_BYTES = 16 * 1024 * 1024
MAX_PIXELS = 1_000_000
MAX_SIZE = 100  # arcseconds, at most 1000 VIS pixels across before service rounding
SOURCE_TABLE = 'dr1.spectra_source'
MOSAIC_TABLE = 'dr1.mosaic_product'
GRISMS = {'red': 'RGS', 'blue': 'BGS'}


def _request(url, params, *, session=None, timeout=120, max_bytes=MAX_BYTES):
    if not np.isfinite(timeout) or timeout <= 0 or max_bytes <= 0:
        raise ValueError('timeout and max_bytes must be positive')
    http = requests if session is None else session
    with http.post(url, data=params, timeout=timeout, stream=True) as response:
        response.raise_for_status()
        if int(response.headers.get('Content-Length', 0)) > max_bytes:
            raise ValueError('ESA response exceeds transfer limit')
        result = bytearray()
        for chunk in response.iter_content(65536):
            if len(result) + len(chunk) > max_bytes:
                raise ValueError('ESA response exceeds transfer limit')
            result.extend(chunk)
    return bytes(result)


def _query(query, *, client=None, session=None, timeout=120):
    if client is not None:
        job = client.launch_job(query)
        if job is None:
            raise RuntimeError('ESA metadata query failed')
        result = job.get_results()
        if result is None:
            raise RuntimeError('ESA metadata query failed')
    else:
        body = _request(TAP_URL, {'REQUEST': 'doQuery', 'LANG': 'ADQL',
                                 'FORMAT': 'csv', 'QUERY': query},
                        session=session, timeout=timeout, max_bytes=2*1024*1024)
        # CSV avoids the archive's invalid VOTable arraysize=100x* on
        # sourcepatch_id. Parse strings first so large integer IDs stay exact.
        if body.lstrip().startswith(b'<'):
            raise ValueError('ESA metadata query returned an XML error instead of CSV')
        reader = csv.DictReader(io.StringIO(body.decode('utf-8')))
        if not reader.fieldnames or len(reader.fieldnames) < 2:
            raise ValueError('ESA metadata query returned no CSV header')
        records = list(reader)
        result = Table(rows=[[r[k] for k in reader.fieldnames] for r in records],
                       names=reader.fieldnames, dtype=[str]*len(reader.fieldnames))
    if len(result) > MAX_ROWS:
        raise ValueError('ESA metadata limit reached; refine the product selection')
    return result


def _path(directory, filename):
    path = str(directory).rstrip('/') + '/' + str(filename)
    if (not path.startswith('/data/euclid/repository_dr1/DR1/')
            or '..' in PurePosixPath(path).parts or not path.endswith('.fits')
            or PurePosixPath(str(filename)).name != str(filename)):
        raise ValueError('Expected an exact ESA DR1 archive FITS path')
    return path


def _single_fits(body):
    """Read one bounded FITS product in memory; never extract paths to disk."""
    if body[:2] == b'PK':
        with zipfile.ZipFile(io.BytesIO(body)) as bundle:
            members = bundle.infolist()
            if (len(members) != 1 or members[0].is_dir()
                    or not members[0].filename.endswith('.fits')
                    or members[0].file_size > MAX_BYTES):
                raise ValueError('Expected one bounded ESA FITS product; ambiguous/oversize ZIP')
            body = bundle.read(members[0])
    if not body.startswith(b'SIMPLE  '):
        raise ValueError('ESA response is not a FITS product')
    return body


def decode_spectrum(body, *, object_id, association, grism, url):
    body = _single_fits(body)
    with fits.open(io.BytesIO(body), memmap=False) as hdus:
        if len(hdus) != 2 or not isinstance(hdus[1], fits.BinTableHDU):
            raise ValueError('Expected primary plus one ESA extracted spectrum table')
        header = hdus[1].header
        if str(header.get('EXTNAME')) != str(object_id):
            raise ValueError('ESA spectrum source ID differs from requested source')
        for key in ('SOURC_ID', 'OBJECTID'):
            if key in header and str(header[key]) != str(object_id):
                raise ValueError('ESA spectrum source ID differs from requested source')
        if header.get('RELEASE', 'DR1').upper() != 'DR1':
            raise ValueError('ESA spectrum release differs from DR1 request')
        # Number is a nonstandard dimensionless diagnostic label in ESA FITS.
        number = u.def_unit('Number', u.dimensionless_unscaled)
        with u.add_enabled_units([number]), warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=u.UnitsWarning,
                                    message="'Number' did not parse as fits unit:")
            raw = QTable.read(hdus[1], format='fits')
        required = {'WAVELENGTH', 'SIGNAL', 'VAR', 'MASK'}
        if not required.issubset(raw.colnames) or not len(raw):
            raise ValueError('ESA spectrum lacks required spectral columns or bins')
        # DPDD: FSCALE is deprecated; scale is already embedded in TUNITs.
        if not all(isinstance(raw[n], u.Quantity) for n in ('WAVELENGTH', 'SIGNAL', 'VAR')):
            raise ValueError('DR1 spectrum requires physical column units')
        native_flux = raw['SIGNAL']
        native_flux.to(FLUX_UNIT)  # Validate dimensions without rescaling native values.
        if raw['VAR'].unit != native_flux.unit**2:
            raise ValueError('DR1 variance unit/scale differs from squared signal unit')
        wave = raw['WAVELENGTH'].to(u.um)
        finite = np.ma.array(wave.value, mask=getattr(wave, 'mask', False)).filled(np.nan)
        if not np.isfinite(finite).any() or np.any(finite[np.isfinite(finite)] <= 0):
            raise ValueError('Invalid DR1 wavelength coverage')
        with np.errstate(invalid='ignore'):
            error = np.sqrt(raw['VAR'])  # Negative variance remains invalid, never clipped.
        result = QTable({'WAVELENGTH': wave, 'FLUX': native_flux, 'ERROR': error,
                         'MASK': raw['MASK'].value if hasattr(raw['MASK'], 'unit') else raw['MASK']})
        for name in ('QUALITY', 'NDITH'):
            if name in raw.colnames:
                value = raw[name]
                result[name] = value.value if hasattr(value, 'unit') else value
        result.meta.update(raw.meta)
    stamp(result.meta, product_type='spectrum', release='DR1', archive='ESA',
          route='tap-association+esa-datalink', OBJECTID=str(object_id),
          SRCPATH=_path(association['file_path'], association['file_name']),
          SRCHDU=_hdu(association['hdu_index']),
          GRISM=grism, GRISMREF='ESA dr1.spectra_source.lambda_range',
          DATASET='dr1', ASSOCTBL=SOURCE_TABLE, SRCURL=url,
          RAWSHA=hashlib.sha256(body).hexdigest())
    return result


def retrieve_spectra(object_id, *, grism=None, product=None, client=None, session=None, timeout=120):
    object_id = _object_id(object_id)
    if grism not in (None, *GRISMS):
        raise ValueError('grism must be red, blue or None')
    if product is not None:
        if not isinstance(product, (tuple, list)) or len(product) != 2:
            raise ValueError('product must be an exact (ESA path, HDU) pair')
        product = (str(product[0]), _hdu(product[1]))
    columns = 'source_id,file_path,file_name,hdu_index,lambda_range'
    query = f'SELECT TOP {MAX_ROWS+1} {columns} FROM {SOURCE_TABLE} WHERE source_id={object_id}'
    if grism:
        query += f" AND lambda_range='{GRISMS[grism]}'"
    rows = _query(query, client=client, session=session, timeout=timeout)
    # DataLink's SOURCE_ID route is verified. Its source-patch route is not:
    # selecting a file locally cannot force that file on an ambiguous service.
    selected = []
    seen = set()
    for row in rows:
        if str(row['source_id']) != object_id or str(row['lambda_range']) not in GRISMS.values():
            raise ValueError('ESA association source/grism is invalid')
        actual_grism = {v:k for k,v in GRISMS.items()}[str(row['lambda_range'])]
        if grism and grism != actual_grism:
            raise ValueError('ESA association grism differs from request')
        path = _path(row['file_path'], row['file_name'])
        hdu = _hdu(row['hdu_index'])
        if actual_grism in seen:
            raise ValueError('Ambiguous ESA spectrum: multiple associations for one grism; '
                             'exact source-patch retrieval is not verified')
        seen.add(actual_grism)
        if product is None or product == (path,hdu):
            selected.append((row,actual_grism))
    results = []
    for row,actual_grism in selected:
        params = {'ID': f'dr1 {object_id}', 'RELEASE': 'dr1',
                  'LINKING_PARAMETER': 'SOURCE_ID', 'RETRIEVAL_TYPE': 'SPECTRA_'+GRISMS[actual_grism],
                  'USE_ZIP_ALWAYS': 'true', 'TAPCLIENT': 'ASTROQUERY'}
        body = _request(SPECTRUM_URL, params, session=session, timeout=timeout)
        results.append(decode_spectrum(body, object_id=object_id, association=row,
                                       grism=actual_grism, url=SPECTRUM_URL+'?'+urlencode(params)))
    return results


def retrieve_spectrum(object_id, **kwargs):
    spectra = retrieve_spectra(object_id, **kwargs)
    if not spectra:
        return None
    if len(spectra) != 1:
        raise ValueError('Ambiguous ESA DR1 spectrum; select a grism or use retrieve_spectra')
    return spectra[0]


def _image_args(ra, dec, radius, band, tiling, selection, tile_id):
    if tiling not in ('wide', 'deep') or band not in ('VIS','Y','J','H'):
        raise ValueError('Choose tiling wide/deep and band VIS/Y/J/H')
    selection = selection or ('survey' if tiling == 'wide' else 'deep')
    if selection not in (('survey','repeat') if tiling == 'wide' else ('deep',)):
        raise ValueError('Mosaic selection disagrees with wide/deep tiling')
    if not np.isfinite(radius) or radius <= 0 or radius > 60:
        raise ValueError('Mosaic search radius must be in (0,60] arcseconds')
    if tile_id is not None and not re.fullmatch('102[0-9]{6}' if tiling=='wide' else '101[0-9]{6}',str(tile_id)):
        raise ValueError('Tile ID disagrees with requested tiling')
    coord = SkyCoord(float(ra)*u.deg,float(dec)*u.deg)
    if not np.isfinite(coord.ra.deg) or not np.isfinite(coord.dec.deg):
        raise ValueError('Coordinates must be finite degrees')
    return coord, selection


def discover_mosaics(ra, dec, *, radius=2, band='VIS', tiling='wide', selection=None,
                     tile_id=None, product=None, client=None, session=None, timeout=120):
    coord, selection = _image_args(ra,dec,radius,band,tiling,selection,tile_id)
    instrument, filter_name = ('VIS','VIS') if band=='VIS' else ('NISP','NIR_'+band)
    purpose = 'DATA_RELEASE' if selection=='survey' else 'DATA_RELEASE_DEEP'
    columns = ('file_name,file_path,instrument_name,filter_name,tile_index,data_set_release,'
               'purpose,processing_mode')
    query = (f'SELECT TOP {MAX_ROWS+1} {columns} FROM {MOSAIC_TABLE} '
             f"WHERE instrument_name='{instrument}' AND filter_name='{filter_name}' "
             f"AND processing_mode='{tiling.upper()}' AND purpose='{purpose}' "
             f"AND INTERSECTS(CIRCLE('ICRS',{coord.ra.deg},{coord.dec.deg},{radius/3600}),fov)=1")
    if tile_id is not None:
        query += f' AND tile_index={int(tile_id)}'
    rows = _query(query, client=client, session=session, timeout=timeout)
    keep = []
    for index, row in enumerate(rows):
        path = _path(row['file_path'],row['file_name'])
        if (str(row['processing_mode']) != tiling.upper() or str(row['purpose']) != purpose
                or str(row['instrument_name']) != instrument or str(row['filter_name']) != filter_name
                or str(row['data_set_release']) not in ('DR1_R1','DR1_R2')
                or (tile_id is not None and str(row['tile_index']) != str(tile_id))):
            raise ValueError('ESA mosaic metadata differs from request')
        if product is None or str(product) == path:
            keep.append(index)
    rows = rows[keep]
    rows['product'] = [_path(r['file_path'],r['file_name']) for r in rows]
    return rows


def retrieve_cutout(ra, dec, *, radius=2, cutout_size=20, band='VIS', tiling='wide', selection=None,
                    tile_id=None, product=None, client=None, session=None, timeout=120):
    coord, selection = _image_args(ra,dec,radius,band,tiling,selection,tile_id)
    if cutout_size is None or not np.isfinite(cutout_size) or not 0 < cutout_size < MAX_SIZE:
        raise ValueError(f'ESA cutout_size must be positive and less than {MAX_SIZE} arcseconds')
    rows = discover_mosaics(ra,dec,radius=radius,band=band,tiling=tiling,selection=selection,
                           tile_id=tile_id,product=product,client=client,session=session,timeout=timeout)
    if not len(rows):
        return None
    if len(rows) != 1:
        raise ValueError('Ambiguous ESA mosaic; select a tile or exact product path')
    row = rows[0]
    params = {'FILEPATH':str(row['product']), 'TAPCLIENT':'ASTROQUERY',
              'POS':f'CIRCLE,{coord.ra.deg},{coord.dec.deg},{cutout_size/7200}'}
    body = _request(CUTOUT_URL,params,session=session,timeout=timeout)
    with fits.open(io.BytesIO(_single_fits(body)),memmap=False) as hdus:
        if len(hdus)!=1 or hdus[0].header.get('NAXIS')!=2:
            raise ValueError('Expected one ESA primary science image')
        header = hdus[0].header.copy()
        if header['NAXIS1']*header['NAXIS2'] > MAX_PIXELS:
            raise ValueError('ESA cutout pixel limit exceeded')
        meta = _metadata(header,band,tiling)
        if meta['kind'] != selection or meta['dataset'] != str(row['data_set_release']):
            raise ValueError('ESA image processing metadata differs from selected mosaic')
        wcs = WCS(header)
        x,y = wcs.world_to_pixel(coord)
        if not (-.5<=x<header['NAXIS1']-.5 and -.5<=y<header['NAXIS2']-.5):
            raise ValueError('ESA cutout WCS does not contain requested position')
        # ESA currently emits EXTNAME with an undefined value in primary
        # cutouts. It is not a science/WCS card and fails FITS re-serialization.
        if 'EXTNAME' in header and header['EXTNAME'] is None:
            header.remove('EXTNAME')
        result = fits.PrimaryHDU(hdus[0].data.copy(),header)
    for key in ('CHECKSUM','DATASUM'):
        result.header.remove(key,ignore_missing=True)
    stamp(result.header,product_type='image',release='DR1',archive='ESA',
          route='tap-mosaic+esa-get_cutout',SRCPATH=str(row['product']),
          SRCURL=CUTOUT_URL+'?'+urlencode(params),TILEID=str(row['tile_index']),
          PRODBAND=band,MOSCLASS=selection,CUTHDR=meta['header_sha256'],
          RAWSHA=hashlib.sha256(body).hexdigest(),CUTBYTES=len(body),
          CUTRA=coord.ra.deg,CUTDEC=coord.dec.deg,CUTSIZE=float(cutout_size))
    # A server cutout does not identify the parent's HDU; do not invent SRCHDU.
    return result
