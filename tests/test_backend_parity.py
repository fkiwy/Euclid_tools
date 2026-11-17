"""Paired archive tests of portable API behavior, independent of live services."""
import io
from unittest.mock import Mock, patch

import astropy.units as u
import numpy as np
import pytest
from astropy.io import fits
from astropy.table import QTable, Table

from euclid_tools import esa_tools, irsa_tools, irsa_spectra, esa_products
from euclid_tools.backends import get_backend
from euclid_tools.provenance import product_provenance
from test_dr1_vega_calibration import sample, calibration
from test_esa_dr1_products import association, client, ID, PATH, Session, mosaic as esa_mosaic
from test_irsa_dr1_spectra import associations, votable, client_for


@pytest.mark.parametrize('name', ['esa', 'irsa'])
@pytest.mark.parametrize('method', ['retrieve_spectra', 'discover_spectra', 'discover_mosaics'])
@pytest.mark.parametrize('release', ['q1', 'dr2', None, 1])
def test_dr1_only_helpers_reject_other_releases_before_query(name, method, release):
    c = Mock(); backend = get_backend(name)
    args = (1, 2) if method == 'discover_mosaics' else ('123',)
    with pytest.raises(ValueError):
        getattr(backend, method)(*args, release=release, client=c)
    assert not c.mock_calls
    module = esa_tools if name == 'esa' else irsa_tools
    with pytest.raises(ValueError):
        getattr(module, method)(*args, release=release, client=c)
    assert not c.mock_calls


@pytest.mark.parametrize('name', ['esa', 'irsa'])
@pytest.mark.parametrize('radius', [0, -1, np.nan, np.inf, True, '5'])
def test_invalid_radii_fail_before_any_catalog_query(name, radius):
    c = Mock()
    with pytest.raises(ValueError):
        get_backend(name).retrieve_objects(80, -51, radius=radius, client=c)
    assert not c.mock_calls


@pytest.mark.parametrize('name', ['esa', 'irsa'])
@pytest.mark.parametrize('value', ['１', '1 OR 1=1', 2**63, -(2**63)-1])
def test_catalog_and_q1_spectrum_ids_have_the_same_validation(name, value):
    c = Mock(); backend = get_backend(name)
    with pytest.raises(ValueError):
        backend.retrieve_object_by_id(value, client=c)
    assert not c.mock_calls
    module = esa_tools if name == 'esa' else irsa_tools
    archive = module.Euclid if name == 'esa' else module.Irsa
    with patch.object(archive, 'get_spectrum' if name == 'esa' else 'query_tap') as query:
        with pytest.raises(ValueError):
            backend.retrieve_spectrum(value)
        query.assert_not_called()


@pytest.mark.parametrize('name', ['esa', 'irsa'])
@pytest.mark.parametrize('release', ['q1', 'dr1'])
def test_nonfinite_coordinates_and_invalid_limits_fail_before_query(name, release):
    c = Mock(); backend = get_backend(name)
    for coords in ((np.nan, 0), (0, np.nan), (np.inf, 0)):
        with pytest.raises(ValueError):
            backend.retrieve_objects(*coords, radius=5, release=release, client=c)
    for limit in (-2, 1.2, True, np.nan):
        with pytest.raises(ValueError):
            backend.retrieve_objects(80, -51, radius=5, limit=limit, release=release, client=c)
    assert not c.mock_calls


def test_catalog_values_photometry_sorting_limit_and_metadata_agree():
    results = []
    for name in ('esa', 'irsa'):
        raw = sample(name)[[0, 0]]
        raw['object_id'] = [-123, -456]
        if name == 'esa':
            # TAP dist is angular separation, not the difference in RA.
            raw['right_ascension'] = [80.0001, 80.]; raw['dist'] = [.00006293204, 0.]
        else:
            raw['ra'] = [80.0001, 80.]
        c = Mock(); c.cone_search.return_value.get_results.return_value = raw
        c.query_region.return_value = raw
        value = get_backend(name).retrieve_objects(80, -51, search_radius=5, limit=10,
                         release='DR1', vega_calibration=calibration(), client=c)
        assert c.ROW_LIMIT == 10
        assert value.meta['catalogue'] == 'wide_survey'
        assert value.meta['release'] == 'DR1'
        assert value['object_id'].tolist() == [-456, -123]
        assert value['separation'].unit == u.arcsec
        results.append(value)
    for name in ('object_id', 'ra', 'dec', 'separation', 'flux_y_templfit', 'flag_y',
                 'VIS_AB_mag', 'Y_AB_err', 'J_Vega_mag', 'H_Vega_err'):
        np.testing.assert_array_equal(results[0][name], results[1][name])
    assert results[0].meta['VEGACAL'] == results[1].meta['VEGACAL']


@pytest.mark.parametrize('name', ['esa', 'irsa'])
def test_failed_catalog_results_are_not_empty_catalogs(name):
    c = Mock(); c.cone_search.return_value = None; c.query_region.return_value = None
    c.launch_job.return_value = None; c.query_tap.return_value = None
    b = get_backend(name)
    with pytest.raises(RuntimeError): b.retrieve_objects(80, -51, radius=5, client=c)
    with pytest.raises(RuntimeError): b.retrieve_object_by_id('123', client=c)
    c.cone_search.return_value = Mock(); c.cone_search.return_value.get_results.return_value = None
    c.launch_job.return_value = Mock(); c.launch_job.return_value.get_results.return_value = None
    c.query_tap.return_value = Mock(); c.query_tap.return_value.to_table.return_value = None
    with pytest.raises(RuntimeError): b.retrieve_objects(80, -51, radius=5, client=c)
    with pytest.raises(RuntimeError): b.retrieve_object_by_id('123', client=c)


@pytest.mark.parametrize('name', ['esa', 'irsa'])
def test_discovery_is_metadata_only_and_deduplicates_exact_parents(name):
    rows = association(2) if name == 'esa' else associations(2)
    if name == 'irsa': rows['hdu'][1] = rows['hdu'][0]
    c = client(rows) if name == 'esa' else client_for(rows)
    identity = ID if name == 'esa' else str(rows['objectid'][0])
    with patch.object(esa_products, '_retrieve_selected') as esa_download, patch.object(irsa_spectra, '_download') as irsa_download:
        result = get_backend(name).discover_spectra(identity, release='DR1', client=c)
        assert len(result) == 1 and result['retrieval_supported'][0]
        assert {'source_path', 'hdu_index', 'grism'}.issubset(result.colnames)
        esa_download.assert_not_called(); irsa_download.assert_not_called()
    if name == 'irsa':
        assert result['grism'][0] == ''
        c.reset_mock()
        with pytest.raises(NotImplementedError, match='no grism'):
            get_backend(name).discover_spectra(identity, grism='red', client=c)
        assert not c.mock_calls


@pytest.mark.parametrize('name', ['esa', 'irsa'])
def test_empty_spectrum_results_and_unassociated_selector_are_distinct(name):
    rows = association() if name == 'esa' else associations()
    c = client(rows[:0]) if name == 'esa' else client_for(rows[:0])
    identity = ID if name == 'esa' else str(rows['objectid'][0])
    b = get_backend(name)
    assert b.retrieve_spectra(identity, release='dr1', client=c) == []
    assert b.retrieve_spectrum(identity, release='dr1', client=c) is None
    c = client(rows) if name == 'esa' else client_for(rows)
    path = PATH if name == 'esa' else str(rows['path'][0])
    with pytest.raises(ValueError, match='not associated'):
        b.retrieve_spectrum(identity, release='dr1', product=(path, 99), client=c)


def test_dr1_spectrum_values_errors_masks_and_saved_files_are_portable(tmp_path):
    body = votable(negative_variance=True)
    native = QTable.read(io.BytesIO(body), format='votable')
    native.meta.update(EXTNAME=ID, SOURC_ID=int(ID))
    stream = io.BytesIO(); native.write(stream, format='fits')
    rows = associations(); rows['objectid'] = [int(ID)]
    results = []
    for name in ('esa', 'irsa'):
        directory = tmp_path/name; directory.mkdir()
        kwargs = {'client': client(association()), 'session': Session(stream.getvalue())} if name == 'esa' else {'client': client_for(rows)}
        with patch.object(irsa_spectra, '_download', return_value=body):
            result = get_backend(name).retrieve_spectrum(ID, release='DR1', grism='red',
                        save_spectrum=True, output_dir=str(directory), object_name='target', **kwargs)
        saved = QTable.read(directory/'target_spectrum.fits')
        record = product_provenance(saved)
        assert record['source_id'] == ID and record['grism'] == 'red'
        assert record['release'] == 'DR1' and len(saved) == 3
        invalid = np.ma.array(saved['ERROR'].value, mask=getattr(saved['ERROR'], 'mask', False)).filled(np.nan)
        assert np.isnan(invalid[1])
        results.append(result)
    for key in ('FLUX', 'ERROR'):
        values = [np.ma.array(r[key].to_value(irsa_spectra.FLUX_UNIT),
                              mask=getattr(r[key], 'mask', False)).filled(np.nan) for r in results]
        np.testing.assert_allclose(*values, equal_nan=True)
    for key in ('WAVELENGTH', 'MASK', 'QUALITY', 'NDITH'):
        np.testing.assert_array_equal(results[0][key], results[1][key])


@pytest.mark.parametrize('name', ['esa', 'irsa'])
def test_mosaic_radius_alias_and_release_are_forwarded_consistently(name):
    module = esa_products if name == 'esa' else __import__('euclid_tools.irsa_images', fromlist=[''])
    columns = {'kind': ['survey'], 'status': ['ready'], 'product': ['native']}
    with patch.object(module, 'discover_mosaics', return_value=Table(columns)) as discover:
        b = get_backend(name)
        b.discover_mosaics(80, -51, search_radius=3, release='DR1')
        assert discover.call_args.kwargs['radius'] == 3
        b.discover_mosaics(80, -51)
        assert discover.call_args.kwargs['radius'] == 2
        discover.reset_mock()
        with pytest.raises(ValueError): b.discover_mosaics(80, -51, radius=2, search_radius=3)
        discover.assert_not_called()


def test_mosaic_discovery_has_portable_fields_and_survey_defaults():
    from test_irsa_dr1_images import mosaic, setup, RA, DEC
    survey, repeat = mosaic(), mosaic('repeat', version='repeat')
    c, fs = setup(survey, repeat)
    irsa = get_backend('irsa').discover_mosaics(RA, DEC, release='dr1', client=c, filesystem=fs)
    esa = get_backend('esa').discover_mosaics(RA, DEC, release='dr1', client=client(esa_mosaic()))
    for table in (esa, irsa):
        assert {'product', 'tile_id', 'band', 'tiling', 'kind', 'status', 'error'}.issubset(table.colnames)
        assert len(table) == 1 and table['kind'][0] == 'survey' and table['status'][0] == 'ready'


def test_q1_esa_variance_remains_invalid_and_ambiguity_or_corruption_is_explicit(tmp_path, monkeypatch):
    table = QTable({'WAVELENGTH': [13000., 15000.] * u.AA,
                    'SIGNAL': [-2., 3.] * irsa_spectra.FLUX_UNIT,
                    'VAR': [-1., 4.] * irsa_spectra.FLUX_UNIT**2, 'MASK': [1, 0]})
    table.meta['FSCALE'] = 1e-16
    path = tmp_path/'q1.fits'; table.write(path)
    c = Mock(); c.get_spectrum.return_value = [str(path)]
    monkeypatch.setattr(esa_tools, 'Euclid', c)
    value = esa_tools.retrieve_spectrum('123')
    assert np.isnan(value['ERROR'][0]) and value['ERROR'][1].value == pytest.approx(2e-16, abs=1e-25)
    np.testing.assert_allclose(value['FLUX'].value, [-2e-16, 3e-16], atol=0)
    c.get_spectrum.return_value = [str(path), str(path)]
    assert len(esa_tools.retrieve_spectrum('123')) == 2
    c.get_spectrum.return_value = [str(path), str(tmp_path/'another.fits')]
    with pytest.raises(ValueError, match='Ambiguous'): esa_tools.retrieve_spectrum('123')
    broken = tmp_path/'broken.fits'; broken.write_bytes(b'not a FITS file')
    c.get_spectrum.return_value = [str(broken)]
    with pytest.raises(OSError): esa_tools.retrieve_spectrum('123')


def test_q1_irsa_legacy_association_ambiguity_does_not_pick_first(monkeypatch):
    c = Mock(); c.query_tap.return_value.to_table.return_value = Table({'path': ['one.vot', 'two.vot']})
    monkeypatch.setattr(irsa_tools, 'Irsa', c)
    with patch.object(irsa_tools.QTable, 'read') as download:
        with pytest.raises(ValueError, match='Ambiguous'): irsa_tools.retrieve_spectrum('123')
        download.assert_not_called()


@pytest.mark.parametrize('name', ['esa', 'irsa'])
@pytest.mark.parametrize('release', ['q1', 'dr1'])
@pytest.mark.parametrize('size,error', [(None, TypeError), (0, ValueError), (-1, ValueError), (np.nan, ValueError)])
def test_cutout_size_errors_are_consistent_before_archive_access(name, release, size, error):
    with pytest.raises(error):
        get_backend(name).retrieve_cutout(80, -51, radius=5, cutout_size=size, band='VIS', release=release)
