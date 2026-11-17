"""Explicit DR1 Vega conventions agree across archives without changing AB/Q1."""
import copy
import json
from unittest.mock import Mock

import astropy.units as u
import numpy as np
import pytest
from astropy.table import Table

from euclid_tools.backends import get_backend
from euclid_tools.photometry import validate_vega_calibration


def calibration():
    # Deliberately synthetic test values; not a scientific calibration.
    return {'release': 'DR1', 'zero_points_jy': {'VIS': 1000., 'Y': 2000., 'J': 1500., 'H': 500.},
            'source': 'Synthetic arithmetic test only', 'passband_reference': 'Test passbands',
            'reference_spectrum': 'Test reference spectrum, zero magnitude by definition'}


def sample(archive):
    table = Table({'object_id': [-123], 'flag_y': [7]})
    for name in ('vis_psf', 'y_templfit', 'j_templfit', 'h_templfit'):
        table['flux_'+name] = [36.31]*u.uJy
        table['fluxerr_'+name] = [3.631]*u.uJy
    if archive == 'esa':
        table['right_ascension'], table['declination'], table['dist'] = [80.], [-51.], [0.]
    else:
        table['ra'], table['dec'] = [80.], [-51.]
    return table


@pytest.mark.parametrize('archive', ['esa', 'irsa'])
@pytest.mark.parametrize('entry', ['cone', 'id'])
@pytest.mark.parametrize('catalogue', ['wide_survey', 'deep_survey'])
def test_opt_in_forwarding_arithmetic_provenance_and_fits(archive, entry, catalogue, tmp_path):
    raw = sample(archive); client = Mock()
    client.cone_search.return_value.get_results.side_effect = lambda: raw.copy()
    client.launch_job.return_value.get_results.side_effect = lambda: raw.copy()
    client.query_region.side_effect = lambda **kw: raw.copy()
    client.query_tap.return_value.to_table.side_effect = lambda: raw.copy()
    backend = get_backend(archive)
    def run(**options):
        kwargs = dict(release='dr1', catalogue=catalogue, client=client, **options)
        return (backend.retrieve_objects(80, -51, radius=5, **kwargs) if entry == 'cone'
                else backend.retrieve_object_by_id('-123', **kwargs))
    baseline = run(); convention = calibration(); before = copy.deepcopy(convention)
    result = run(vega_calibration=convention)
    for band, zero in convention['zero_points_jy'].items():
        assert result[f'{band}_Vega_mag'][0] == pytest.approx(round(-2.5*np.log10(36.31e-6/zero), 3))
        np.testing.assert_array_equal(result[f'{band}_AB_mag'], baseline[f'{band}_AB_mag'])
        np.testing.assert_array_equal(result[f'{band}_Vega_err'], baseline[f'{band}_AB_err'])
        assert result[f'{band}_Vega_mag'].unit == u.mag
    assert result['flag_y'][0] == 7
    np.testing.assert_array_equal(result['flux_y_templfit'], raw['flux_y_templfit'])
    assert convention == before
    assert json.loads(result.meta['VEGACAL']) == convention
    assert 'user-supplied' in result.meta['VEGASTAT']
    path = tmp_path/'photometry.fits'; result.write(path)
    saved = Table.read(path)
    assert saved.meta['VEGACAL'] == result.meta['VEGACAL']
    assert saved.meta['VEGASHA'] == result.meta['VEGASHA']


@pytest.mark.parametrize('archive', ['esa', 'irsa'])
@pytest.mark.parametrize('entry', ['cone', 'id'])
@pytest.mark.parametrize('problem', ['q1', 'missing_band', 'nan', 'negative', 'bool', 'wrong_release', 'no_source', 'no_passband', 'no_spectrum'])
def test_invalid_or_q1_calibration_fails_before_query(archive, entry, problem):
    value = calibration(); release = 'dr1'; client = Mock()
    if problem == 'q1': release = 'q1'
    elif problem == 'missing_band': del value['zero_points_jy']['H']
    elif problem in ('nan', 'negative', 'bool'):
        value['zero_points_jy']['Y'] = {'nan': float('nan'), 'negative': -1., 'bool': True}[problem]
    elif problem == 'wrong_release': value['release'] = 'Q1'
    else: value[{'no_source': 'source', 'no_passband': 'passband_reference', 'no_spectrum': 'reference_spectrum'}[problem]] = ''
    backend = get_backend(archive)
    with pytest.raises(ValueError):
        options = dict(release=release, vega_calibration=value, client=client)
        if entry == 'cone': backend.retrieve_objects(80, -51, radius=5, **options)
        else: backend.retrieve_object_by_id('-123', **options)
    assert not client.mock_calls


def test_no_calibration_and_unknown_fields():
    assert validate_vega_calibration(None, release='q1') is None
    bad = calibration(); bad['magnitude_offset'] = 1
    with pytest.raises(ValueError, match='exactly'):
        validate_vega_calibration(bad)


def test_masked_measurements_remain_masked():
    from astropy.table import MaskedColumn
    from euclid_tools.photometry import add_dr1_magnitudes
    table = sample('esa')[[0, 0]]
    for name in ('vis_psf', 'y_templfit', 'j_templfit', 'h_templfit'):
        table['flux_'+name] = MaskedColumn([36.31, 10.], mask=[False, True], unit=u.uJy)
    add_dr1_magnitudes(table, vega_calibration=calibration())
    assert table['Y_Vega_mag'].mask.tolist() == [False, True]
    assert table['Y_AB_mag'].mask.tolist() == [False, True]
