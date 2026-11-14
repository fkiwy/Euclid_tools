"""Q1 must use the current q1 schema across Astroquery API generations."""
from types import SimpleNamespace
from unittest.mock import Mock

import pytest
from astropy.table import Table

from euclid_tools import esa_tools


@pytest.mark.parametrize('current_api',[False,True])
def test_q1_spectrum_uses_explicit_schema_with_old_and_new_id_keywords(monkeypatch,current_api):
    calls=[]
    def old(*,source_id,schema,retrieval_type,output_file):
        calls.append((source_id,schema,retrieval_type));return []
    def current(*,ids,schema,retrieval_type,output_file):
        calls.append((ids,schema,retrieval_type));return []
    monkeypatch.setattr(esa_tools,'Euclid',SimpleNamespace(get_spectrum=current if current_api else old))
    assert esa_tools.retrieve_spectrum('2664850113649936423') is None
    assert calls==[('2664850113649936423','q1','SPECTRA_RGS')]


def test_q1_mosaic_query_failure_is_explicit_and_never_requests_a_cutout(monkeypatch):
    client=Mock();client.launch_job.return_value=None
    monkeypatch.setattr(esa_tools,'Euclid',client)
    with pytest.raises(RuntimeError,match='Q1 mosaic query failed'):
        esa_tools.retrieve_cutout(266.4850113,64.9936424,radius=5,cutout_size=3,band='VIS')
    query=client.launch_job.call_args.args[0]
    assert 'FROM q1.mosaic_product' in query and 'sedm.' not in query
    client.get_cutout.assert_not_called()


def test_empty_q1_mosaic_result_is_not_a_query_failure(monkeypatch):
    client=Mock();client.launch_job.return_value.get_results.return_value=Table()
    monkeypatch.setattr(esa_tools,'Euclid',client)
    assert esa_tools.retrieve_cutout(1,2,radius=5,cutout_size=3,band='VIS') is None
    client.get_cutout.assert_not_called()


def test_missing_cutout_response_does_not_raise_indexing_error(monkeypatch):
    client=Mock();client.launch_job.return_value.get_results.return_value=Table({
        'file_name':['tile.fits'],'file_path':['/native/q1'],'instrument_name':['VIS'],'tile_index':[102158582]})
    client.get_cutout.return_value=None
    monkeypatch.setattr(esa_tools,'Euclid',client)
    assert esa_tools.retrieve_cutout(1,2,radius=5,cutout_size=3,band='VIS') is None
    assert client.get_cutout.call_args.kwargs['file_path']=='/native/q1/tile.fits'
