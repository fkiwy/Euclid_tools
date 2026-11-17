"""Explicit ESA survey selection must control every catalog entry point."""
from pathlib import Path
from unittest.mock import Mock, patch

import numpy as np
import pytest
from astropy.table import Table

from euclid_tools import esa_tools
from euclid_tools.backends import get_backend
from euclid_tools.releases import catalog_table


@pytest.mark.parametrize('catalogue,expected',[(None,'wide_survey'),('wide_survey','wide_survey'),('deep_survey','deep_survey')])
def test_cone_id_and_metadata_select_the_same_survey(catalogue,expected):
    raw=Table.read(Path(__file__).parent/'data/dr1_sample.csv',format='ascii.csv')[:1]
    raw['object_id']=[-597913580476826189]
    raw['dist']=[0.0001]
    c=Mock();c.cone_search.return_value.get_results.side_effect=lambda:raw.copy()
    c.launch_job.return_value.get_results.side_effect=lambda:raw.copy()
    backend=get_backend('esa')
    found=backend.retrieve_objects(59.7913643,-47.6826163,search_radius=5,
                                  release='dr1',catalogue=catalogue,client=c)
    table_name='dr1.mer_catalogue_'+expected
    assert c.cone_search.call_args.kwargs['table_name']==table_name
    if expected=='deep_survey':
        assert c.cone_search.call_args.kwargs['output_format']=='csv'
    assert found.meta['catalogue']==expected and found.meta['catalogue_table']==table_name
    assert found['separation'][0]==0.36
    identity=str(found['object_id'][0])
    by_id=backend.retrieve_object_by_id(identity,release='dr1',catalogue=catalogue,client=c)
    assert c.launch_job.call_args.args[0]==f'SELECT * FROM {table_name} WHERE object_id = {identity}'
    if expected=='deep_survey':
        assert c.launch_job.call_args.kwargs['output_format']=='csv'
    assert str(by_id['object_id'][0])=='-597913580476826189'
    np.testing.assert_allclose(found['Y_AB_mag'],by_id['Y_AB_mag'])
    assert not any('Vega' in n or 'aper_AB' in n for n in found.colnames)
    columns=Mock();columns.columns=[]
    with patch.object(esa_tools.Euclid,'load_table',return_value=columns) as load:
        esa_tools.print_catalog_info(release='dr1',catalogue=catalogue)
    load.assert_called_once_with(table_name)


@pytest.mark.parametrize('release,catalogue',[('q1','deep_survey'),('q1','wide_survey'),
    ('dr1',''),('dr1','deep_visits'),('dr1','deep_survey;drop table x'),('dr1',[]),('dr1',1)])
def test_invalid_selection_fails_before_query(release,catalogue):
    c=Mock()
    for run in [lambda:esa_tools.retrieve_object_by_id('1',release=release,catalogue=catalogue,client=c),
                lambda:esa_tools.retrieve_objects(1,2,radius=1,release=release,catalogue=catalogue,client=c)]:
        with pytest.raises(ValueError):run()
    c.launch_job.assert_not_called();c.cone_search.assert_not_called()


def test_defaults_and_q1_photometry_are_unchanged():
    assert catalog_table()=='q1.mer_catalogue'
    assert catalog_table('dr1')=='dr1.mer_catalogue_wide_survey'
    raw=Table({'object_id':[1],'right_ascension':[1.],'declination':[2.],'dist':[0.]})
    for name in ('vis_psf','y_templfit','j_templfit','h_templfit',
                 'vis_4fwhm_aper','y_2fwhm_aper','j_2fwhm_aper','h_2fwhm_aper'):
        raw['flux_'+name]=[36.31];raw['fluxerr_'+name]=[3.631]
    c=Mock();c.cone_search.return_value.get_results.return_value=raw
    result=esa_tools.retrieve_objects(1,2,radius=1,client=c)
    assert c.cone_search.call_args.kwargs['table_name']=='q1.mer_catalogue'
    assert 'Y_2fwhm_aper_AB_mag' in result.colnames and 'Y_Vega_mag' in result.colnames


@pytest.mark.parametrize('catalogue,tiling,selection',[(None,'wide','survey'),('deep_survey','deep','deep')])
@pytest.mark.parametrize('override',[False,True])
def test_example_keeps_catalog_id_and_image_collection_consistent(monkeypatch,catalogue,tiling,selection,override):
    override = ('101280483' if catalogue else '102158276') if override else None
    import euclid_tools
    from euclid_tools import shared,image_plotter
    script=Path(__file__).resolve().parents[1]/'esa_tools_example.py'
    source=script.read_text().replace('release = "q1"','release = "dr1"',1)
    source=source.replace('catalogue = None',f'catalogue = {catalogue!r}',1)
    source=source.replace('image_tile_id = None',f'image_tile_id = {override!r}',1)
    identity='-597913580476826189' if catalogue else '2731837857643705492'
    backend=Mock();backend.retrieve_objects.return_value=Table({'object_id':[int(identity)],'tile_index':[101280484 if catalogue else 102158277]})
    backend.retrieve_object_by_id.return_value=None;backend.retrieve_spectrum.return_value=None
    backend.retrieve_cutout.return_value=None
    monkeypatch.setattr(euclid_tools,'get_backend',lambda name:backend)
    info=Mock();monkeypatch.setattr(esa_tools,'print_catalog_info',info)
    monkeypatch.setattr(shared,'print_results',lambda result:None)
    monkeypatch.setattr(image_plotter,'plot_images',lambda *args,**kwargs:None)
    exec(compile(source,str(script),'exec'),{'__name__':'__main__'})
    info.assert_called_once_with(release='dr1',catalogue=catalogue)
    assert backend.retrieve_objects.call_args.kwargs['catalogue']==catalogue
    assert backend.retrieve_objects.call_args.args[:2]==(273.18378579,64.37054921)
    backend.retrieve_object_by_id.assert_called_once_with(identity,release='dr1',catalogue=catalogue,vega_calibration=None)
    backend.retrieve_spectrum.assert_called_once_with(identity,release='dr1',grism='red',product=None)
    assert backend.retrieve_cutout.call_count==4
    for call in backend.retrieve_cutout.call_args_list:
        assert call.kwargs['tiling']==tiling and call.kwargs['selection']==selection
        assert call.kwargs['tile_id']==(override or ('101280484' if catalogue else '102158277'))
