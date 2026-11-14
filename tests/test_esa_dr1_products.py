"""Verified ESA transport, release isolation, scaled FITS and selection contracts."""
import io
from pathlib import Path
from unittest.mock import Mock, patch
import zipfile

import astropy.units as u
import numpy as np
import pytest
import requests
from astropy.io import fits
from astropy.table import QTable, Table

from euclid_tools import esa_products as esa, esa_tools
from euclid_tools.backends import get_backend
from euclid_tools.provenance import product_provenance

ID = '2731837857643705492'
DIRECTORY = '/data/euclid/repository_dr1/DR1/R2/SIR_COMBINED_SPECTRA/102158277'
FILENAME = 'EUC_SIR_W-COMBSPEC_102158277_20250902T225801.134941Z.fits'
PATH = DIRECTORY + '/' + FILENAME
DATA = Path(__file__).parent / 'data'


def association(count=1):
    return Table({'source_id':[int(ID)]*count,'file_path':[DIRECTORY]*count,
                  'file_name':[FILENAME]*count,'hdu_index':[12]*count,'lambda_range':['RGS']*count})


def spectrum_body(*, wrong_id=False, negative=False, unitless=False, zipped=True):
    t=QTable.read(DATA/'esa_dr1_spectrum_excerpt.ecsv',format='ascii.ecsv')
    if wrong_id:t.meta['EXTNAME']='123'
    if negative:t['VAR'][4]=-1*t['VAR'].unit
    if unitless:t['SIGNAL']=t['SIGNAL'].value
    stream=io.BytesIO(); t.write(stream,format='fits')
    if not zipped:return stream.getvalue()
    result=io.BytesIO()
    with zipfile.ZipFile(result,'w',compression=zipfile.ZIP_DEFLATED) as z:
        z.writestr('SPECTRA_RGS-dr1 '+ID+'.fits',stream.getvalue())
    return result.getvalue()


def image_body():
    header=fits.Header.fromtextfile(DATA/'esa_dr1_wide_cutout_header.txt')
    h=fits.PrimaryHDU(np.arange(31*31,dtype=np.float32).reshape(31,31),header)
    stream=io.BytesIO();h.writeto(stream,output_verify="ignore");return stream.getvalue()


def mosaic(count=1):
    return Table({'file_name':['mosaic.fits']*count,
        'file_path':['/data/euclid/repository_dr1/DR1/R2/MER/102158277/VIS']*count,
        'instrument_name':['VIS']*count,'filter_name':['VIS']*count,'tile_index':[102158277]*count,
        'data_set_release':['DR1_R2']*count,'purpose':['DATA_RELEASE']*count,'processing_mode':['WIDE']*count})


class Session:
    def __init__(self,*bodies):self.bodies=list(bodies);self.calls=[]
    def post(self,url,**kwargs):
        self.calls.append((url,kwargs))
        response=Mock();response.__enter__=Mock(return_value=response);response.__exit__=Mock(return_value=False)
        response.headers={};response.iter_content.return_value=[self.bodies.pop(0)]
        return response


def client(rows):
    c=Mock();c.launch_job.return_value.get_results.return_value=rows;return c


def test_actual_backend_preserves_native_scale_masks_bins_and_negative_variance(tmp_path):
    raw=QTable.read(DATA/'esa_dr1_spectrum_excerpt.ecsv',format='ascii.ecsv')
    session=Session(spectrum_body(negative=True))
    result=get_backend('esa').retrieve_spectrum(ID,release='dr1',grism='red',
        client=client(association()),session=session,save_spectrum=True,output_dir=str(tmp_path))
    np.testing.assert_allclose(result['FLUX'].to_value(esa.FLUX_UNIT),raw['SIGNAL'].to_value(esa.FLUX_UNIT))
    assert result['FLUX'][4].to_value(esa.FLUX_UNIT)==pytest.approx(-1.353587e-17,abs=1e-24)
    assert np.isnan(result['ERROR'][4]) and len(result)==5
    assert result['ERROR'][3].to_value(esa.FLUX_UNIT)**2 == pytest.approx(raw['VAR'][3].to_value(esa.FLUX_UNIT**2),rel=1e-6,abs=0)
    np.testing.assert_array_equal(result['MASK'],raw['MASK'])
    cached=QTable.read(tmp_path/f'E{ID}_spectrum.fits')
    assert product_provenance(cached)['source_file']==PATH
    assert cached.meta['SRCHDU']==12 and cached.meta['GRISM']=='red'
    assert session.calls[0][1]['data']['ID']=='dr1 '+ID
    assert session.calls[0][1]['data']['RELEASE']=='dr1'
    assert session.calls[0][1]['data']['RETRIEVAL_TYPE']=='SPECTRA_RGS'


@pytest.mark.parametrize('options,match',[({'wrong_id':True},'source ID'),({'unitless':True},'physical')])
def test_wrong_source_and_unverified_units_fail(options,match):
    with pytest.raises(ValueError,match=match):
        esa.retrieve_spectrum(ID,client=client(association()),session=Session(spectrum_body(**options)))


def test_ambiguous_same_grism_rejected_before_product_download_even_with_selector():
    session=Session()
    with pytest.raises(ValueError,match='Ambiguous ESA spectrum'):
        esa.retrieve_spectrum(ID,product=(PATH,12),client=client(association(2)),session=session)
    assert not session.calls


def test_empty_is_none_but_metadata_errors_propagate():
    assert esa.retrieve_spectrum(ID,client=client(association(0))) is None
    c=Mock();c.launch_job.return_value=None
    with pytest.raises(RuntimeError,match='query failed'):esa.retrieve_spectrum(ID,client=c)
    with pytest.raises(ValueError,match='metadata limit'):esa.retrieve_spectrum(ID,client=client(association(101)))
    with pytest.raises(ValueError):esa.retrieve_spectrum('1 OR 1=1')


def test_transport_http_failure_and_byte_budget():
    with pytest.raises(ValueError,match='transfer limit'):
        esa._request(esa.SPECTRUM_URL,{},session=Session(b'oversize'),max_bytes=2)
    s=Session(b'');orig=s.post
    def post(*a,**k):
        r=orig(*a,**k);r.raise_for_status.side_effect=requests.HTTPError('fixture 503');return r
    s.post=post
    with pytest.raises(requests.HTTPError):esa._request(esa.SPECTRUM_URL,{},session=s)


def test_ambiguous_zip_and_html_are_not_empty_spectra():
    stream=io.BytesIO()
    with zipfile.ZipFile(stream,'w') as z:
        z.writestr('one.fits',b'');z.writestr('two.fits',b'')
    with pytest.raises(ValueError,match='one bounded'):esa._single_fits(stream.getvalue())
    with pytest.raises(ValueError,match='not a FITS'):esa._single_fits(b'<html>error</html>')


def test_real_csv_transport_preserves_integer_precision_and_rejects_xml():
    columns='source_id,file_path,file_name,hdu_index,lambda_range\n'
    row=f'{ID},{DIRECTORY},{FILENAME},12,RGS\n'
    session=Session((columns+row).encode(),spectrum_body())
    result=esa.retrieve_spectrum(ID,session=session,grism='red')
    assert result.meta['OBJECTID']==ID
    assert session.calls[0][0]==esa.TAP_URL
    assert 'dr1.spectra_source' in session.calls[0][1]['data']['QUERY']
    with pytest.raises(ValueError,match='XML error'):
        esa._query('unused',session=Session(b'<VOTABLE>query error</VOTABLE>'))


def test_cutout_preserves_pixels_wcs_calibration_and_processing_identity():
    c=client(mosaic());s=Session(image_body())
    result=get_backend('esa').retrieve_cutout(273.18378579,64.37054921,radius=1,cutout_size=3,
        band='VIS',release='dr1',tile_id='102158277',client=c,session=s)
    np.testing.assert_array_equal(result.data,np.arange(961,dtype=np.float32).reshape(31,31))
    original=fits.Header.fromtextfile(DATA/'esa_dr1_wide_cutout_header.txt')
    for key in ('BUNIT','MAGZERO','PPLANID','PIPDEFID','PSWREL','CRPIX1','CRPIX2','CRVAL1','CRVAL2'):
        assert result.header[key]==original[key]
    record=product_provenance(result,source_id=ID)
    assert record['archive']=='ESA' and record['mosaic_class']=='survey'
    assert record['source_header_sha256'] is None and len(record['cutout_header_sha256'])==64
    assert record['source_id'] is None and record['requested_source_id']==ID
    assert record['source_hdu'] is None
    assert s.calls[0][0]==esa.CUTOUT_URL
    assert "purpose='DATA_RELEASE'" in c.launch_job.call_args.args[0]


def test_image_ambiguity_bad_metadata_limits_and_position():
    s=Session()
    with pytest.raises(ValueError,match='Ambiguous ESA mosaic'):
        esa.retrieve_cutout(273.18378579,64.37054921,client=client(mosaic(2)),session=s)
    assert not s.calls
    bad=mosaic();bad['purpose'][0]='WRONG'
    with pytest.raises(ValueError,match='metadata differs'):
        esa.retrieve_cutout(273.18378579,64.37054921,client=client(bad))
    with pytest.raises(ValueError,match='cutout_size'):
        esa.retrieve_cutout(1,2,cutout_size=1000)
    with pytest.raises(ValueError,match='WCS'):
        esa.retrieve_cutout(1,2,client=client(mosaic()),session=Session(image_body()))


def test_legacy_q1_route_is_preserved_and_dr1_does_not_call_it():
    with patch.object(esa_tools.Euclid,'get_spectrum',return_value=[]) as q1:
        assert esa_tools.retrieve_spectrum(ID) is None
        q1.assert_called_once()
    with patch.object(esa_tools.Euclid,'get_spectrum') as q1:
        assert esa_tools.retrieve_spectrum(ID,release='dr1',client=client(association(0))) is None
        q1.assert_not_called()


def test_distinct_grisms_plural_and_singular_contract():
    rows=association(2);rows['lambda_range'][1]='BGS'
    session=Session(spectrum_body(),spectrum_body())
    results=get_backend('esa').retrieve_spectra(ID,release='dr1',client=client(rows),session=session)
    assert [r.meta['GRISM'] for r in results]==['red','blue']
    assert [c[1]['data']['RETRIEVAL_TYPE'] for c in session.calls]==['SPECTRA_RGS','SPECTRA_BGS']
    with pytest.raises(ValueError,match='Ambiguous ESA DR1 spectrum'):
        esa.retrieve_spectrum(ID,client=client(rows),session=Session(spectrum_body(),spectrum_body()))


def test_empty_image_discovery_returns_none_without_product_request():
    session=Session()
    assert esa.retrieve_cutout(1,2,client=client(mosaic(0)),session=session) is None
    assert not session.calls
