"""Shared identity survives FITS without changing science arrays or headers."""
import json

import astropy.units as u
import numpy as np
import pytest
from astropy.io import fits
from astropy.table import QTable

from euclid_tools.provenance import product_provenance, provenance_json, stamp, SCHEMA


def test_spectrum_fits_roundtrip_preserves_exact_id_and_physical_values(tmp_path):
    table = QTable({'WAVELENGTH': [1., 2.] * u.um, 'FLUX': [-2., 3.] * u.Jy})
    stamp(table.meta, product_type='spectrum', release='q1', archive='irsa',
          route='tap-association+spectrumdm', OBJECTID='-2664850113649936423',
          SRCPATH='euclid/q1/SIR/test.fits', SRCHDU=2338, GRISM='RGS')
    before = provenance_json(table)
    destination = tmp_path / 'spectrum.fits'
    table.write(destination)
    saved = QTable.read(destination)
    assert provenance_json(saved) == before
    record = json.loads(before)
    assert record['schema'] == SCHEMA
    assert record['source_id'] == '-2664850113649936423'
    assert record['source_hdu'] == 2338 and record['grism'] == 'RGS'
    assert record['processing_plan'] is None
    np.testing.assert_array_equal(saved['FLUX'].to_value(u.Jy), [-2., 3.])


def test_image_identity_and_requested_source_are_distinct(tmp_path):
    hdu = fits.PrimaryHDU(np.arange(4,dtype=np.float32).reshape(2,2))
    hdu.header.update(CRPIX1=1., BUNIT='ADU/s', MAGZERO=24.5, DATASETR='DR1_R2',
                      PPLANID='official-plan', PIPDEFID='MER_ProcessTile', PSWREL='6.0',
                      MOSCLASS='survey')
    stamp(hdu.header, product_type='image', release='dr1', archive='irsa',
          route='sia+s3-section', SRCPATH='s3://bucket/mosaic.fits', SRCHDU=0)
    path = tmp_path / 'cutout.fits'; hdu.writeto(path)
    with fits.open(path) as hdus:
        record = product_provenance(hdus[0],source_id='1234567890123456789')
        assert record['source_id'] is None
        assert record['requested_source_id'] == '1234567890123456789'
        assert record['processing_dataset'] == 'DR1_R2'
        assert record['processing_plan'] == 'official-plan'
        assert record['processing_pipeline'] == 'MER_ProcessTile'
        assert record['processing_version'] == '6.0'
        assert record['mosaic_class'] == 'survey'
        assert record['retrieval_route'] == 'sia+s3-section'
        assert hdus[0].header['CRPIX1'] == 1 and hdus[0].header['MAGZERO'] == 24.5
        np.testing.assert_array_equal(hdus[0].data,hdu.data)


def test_missing_fields_are_unknown_and_wrong_source_is_rejected():
    record = product_provenance({},source_id='12')
    assert record['release'] is None and record['source_file'] is None
    assert record['grism'] is None and record['source_hdu'] is None
    with pytest.raises(ValueError,match='source ID'):
        product_provenance({'OBJECTID':'13'},source_id='12')


def test_q1_image_default_retains_section_selection_and_processing_identity(tmp_path, monkeypatch):
    from astropy.table import Table
    from astropy.wcs import WCS
    from unittest.mock import Mock
    from euclid_tools import irsa_tools
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [5,5]; wcs.wcs.crval = [10,20]
    wcs.wcs.cdelt = [-1/3600,1/3600]; wcs.wcs.ctype = ['RA---TAN','DEC--TAN']
    source = fits.PrimaryHDU(np.arange(100,dtype=np.float32).reshape(10,10),wcs.to_header())
    source.header.update(PPLANID='Q1-plan',PIPDEFID='Q1-pipeline',PSWREL='Q1-version',DATASETR='Q1')
    path=tmp_path/'q1.fits';source.writeto(path)
    client=Mock();client.query_sia.return_value=Table({'dataproduct_subtype':['science'],
        'energy_bandpassname':['VIS'],'access_url':[str(path)]})
    monkeypatch.setattr(irsa_tools,'Irsa',client)
    result=irsa_tools.retrieve_cutout(10,20,radius=1,cutout_size=2,band='VIS')
    record=product_provenance(result)
    assert record['release']=='Q1' and record['retrieval_route']=='sia+fsspec-section'
    assert record['source_file']==str(path) and record['source_hdu']==0
    assert record['processing_plan']=='Q1-plan'
    assert record['mosaic_class'] is None  # never assign DR1 survey/repeat meaning
    assert result.data.shape==(2,2)
    assert client.query_sia.call_args.kwargs['collection']=='euclid_DpdMerBksMosaic'
