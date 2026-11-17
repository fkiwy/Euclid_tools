"""Small ESA wide/deep cone → selected ID → spectrum validation.

A requests-backed Astroquery-compatible cone client exercises the public helper
in environments where Astroquery's http.client transport cannot use the proxy.
It uses native TAP types and the standard DISTANCE result in degrees.
"""
import argparse
import io
import json
from pathlib import Path
import time

import astropy.units as u
import numpy as np
from astropy.io.votable import parse

from euclid_tools import esa_tools
from euclid_tools.esa_products import _request, TAP_URL
from euclid_tools.provenance import product_provenance

COLUMNS=('object_id,right_ascension,declination,flux_vis_psf,fluxerr_vis_psf,'
         'flux_y_templfit,fluxerr_y_templfit,flux_j_templfit,fluxerr_j_templfit,'
         'flux_h_templfit,fluxerr_h_templfit,det_quality_flag,flag_vis,flag_y,flag_j,flag_h,'
         'point_like_flag,point_like_prob,gaia_match_quality')


class ConeClient:
    def cone_search(self,*,coordinate,radius,table_name,ra_column_name,dec_column_name,columns,async_job,output_format="csv"):
        ra,dec=coordinate.ra.deg,coordinate.dec.deg
        point=f"POINT('ICRS',{ra_column_name},{dec_column_name})"
        q=(f'SELECT TOP 10 {COLUMNS},DISTANCE({point},'
           f"POINT('ICRS',{ra},{dec})) AS dist FROM {table_name} "
           f"WHERE CONTAINS({point},CIRCLE('ICRS',{ra},{dec},{radius.to_value(u.deg)}))=1")
        body=_request(TAP_URL,{'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'votable','QUERY':q},timeout=40,max_bytes=2*1024*1024)
        result=parse(io.BytesIO(body))
        for resource in result.resources:
            for info in resource.infos:
                if info.name=='QUERY_STATUS' and info.value!='OK':raise ValueError(info.content)
        table=result.get_first_table().to_table()
        class Job:
            def get_results(self):return table
        return Job()


def main():
    cli=argparse.ArgumentParser(__doc__);cli.add_argument('--output',type=Path,required=True);args=cli.parse_args()
    targets=[(58.1332495,-49.1830038,{'wide_survey':'-581332416491830015','deep_survey':'-581332463491830110'}),
             (59.7913643,-47.6826163,{'wide_survey':'-597913631476826029','deep_survey':'-597913580476826189'})]
    checks=[]
    for ra,dec,identities in targets:
        for catalogue,expected_id in identities.items():
            start=time.monotonic()
            objects=esa_tools.retrieve_objects(ra,dec,radius=1,release='dr1',catalogue=catalogue,client=ConeClient())
            assert objects is not None and str(objects[0]['object_id'])==expected_id
            expected_y=np.round(-2.5*np.log10(float(objects['flux_y_templfit'][0])*1e-6/3631),3)
            assert float(objects['Y_AB_mag'][0])==expected_y
            check={'ra':ra,'dec':dec,'catalogue':catalogue,'catalogue_table':objects.meta['catalogue_table'],
                   'object_id':expected_id,'separation_arcsec':float(objects['separation'][0]),
                   'y_ab_mag':float(objects['Y_AB_mag'][0]),'catalogue_status':'ok'}
            try:
                spectrum=esa_tools.retrieve_spectrum(str(objects[0]['object_id']),release='dr1',grism='red',timeout=40)
                assert spectrum.meta['OBJECTID']==expected_id and len(spectrum)==531
                check.update(spectrum_bins=len(spectrum),ndith_zero_bins=int(np.sum(spectrum['NDITH']==0)),
                             maximum_ndith=int(np.max(spectrum['NDITH'])),
                             spectrum_provenance=product_provenance(spectrum),spectrum_status='ok')
            except Exception as exc:
                check.update(spectrum_status='error',spectrum_error=repr(exc))
            check['seconds']=time.monotonic()-start
            checks.append(check)
            print(catalogue,expected_id,check['spectrum_status'],flush=True)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps({'transport':'public ESA helper with requests-backed native TAP cone adapter; public spectrum method',
        'date':'2026-11-14','checks':checks},indent=2)+'\n')
    if any(c['spectrum_status']!='ok' for c in checks):raise SystemExit(1)


if __name__=='__main__':main()
