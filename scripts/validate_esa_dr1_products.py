"""Small public ESA DR1 adapter checks; no full survey products are downloaded."""
import argparse
import json
from pathlib import Path
import time

import astropy.units as u
import numpy as np
import requests

from euclid_tools.backends import get_backend
from euclid_tools.provenance import product_provenance


class AuditSession:
    def __init__(self):self.http=requests.Session();self.requests=[]
    def post(self,url,**kwargs):
        start=time.monotonic()
        response=self.http.post(url,**kwargs)
        item={'url':url,'parameters':kwargs['data'],'status':response.status_code,
              'header_seconds':time.monotonic()-start,'bytes':0}
        self.requests.append(item)
        original=response.iter_content
        def chunks(*args,**options):
            for chunk in original(*args,**options):
                item['bytes']+=len(chunk)
                yield chunk
        response.iter_content=chunks
        return response


def main():
    cli=argparse.ArgumentParser(__doc__)
    cli.add_argument('--output',type=Path,required=True)
    args=cli.parse_args()
    backend=get_backend('esa');session=AuditSession()
    report={'sources':['https://www.cosmos.esa.int/web/euclid/dr1-data',
        'https://astroquery.readthedocs.io/en/latest/esa/euclid/euclid.html',
        'https://euclid.esac.esa.int/dr/dr1/dpdd/sirdpd/dpcards/sir_combinedspectra.html'],
        'checks':[],'requests':session.requests}
    for grism,identity in [('red','2731837857643705492'),('blue','2681886189660473837')]:
        start=time.monotonic()
        try:
            result=backend.retrieve_spectrum(identity,release='dr1',grism=grism,session=session,timeout=30)
            assert result is not None and len(result)>0
            flux=result['FLUX'].to_value(u.erg/(u.s*u.cm**2*u.AA))
            report['checks'].append({'kind':'spectrum','grism':grism,'source_id':identity,'status':'ok',
                'bins':len(result),'finite_flux':int(np.isfinite(flux).sum()),
                'flux_unit':str(result['FLUX'].unit),'seconds':time.monotonic()-start,
                'provenance':product_provenance(result)})
        except Exception as exc:report['checks'].append({'kind':'spectrum','grism':grism,'status':'error','error':repr(exc)})
    cases=[{'ra':273.18378579,'dec':64.37054921,'band':b,'tiling':'wide','selection':'survey','tile_id':'102158277'} for b in ('VIS','Y','J','H')]
    cases += [{'ra':273.018929496624,'dec':68.231870619389,'band':'VIS','tiling':'deep','selection':'deep','tile_id':'101839130'},
              {'ra':273.18378579,'dec':64.37054921,'band':'VIS','tiling':'wide','selection':'repeat','tile_id':'102158277',
               'product':'/data/euclid/repository_dr1/DR1/R1/MER/102158277/VIS/EUC_MER_BGSUB-MOSAIC-VIS_TILE102158277-84DEE3_20250706T145426.302735Z_00.00.fits'}]
    for case in cases:
        start=time.monotonic()
        try:
            result=backend.retrieve_cutout(**case,radius=1,cutout_size=3,release='dr1',session=session,timeout=30)
            assert result is not None
            # Re-serialization is part of the live check (ESA can emit blank EXTNAME).
            import io
            result.writeto(io.BytesIO())
            report['checks'].append({'kind':'image',**case,'status':'ok','shape':list(result.data.shape),
                'seconds':time.monotonic()-start,'provenance':product_provenance(result),
                'bunit':result.header['BUNIT'],'magzero':result.header['MAGZERO']})
        except Exception as exc:report['checks'].append({'kind':'image',**case,'status':'error','error':repr(exc)})
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'checks':[{k:v for k,v in r.items() if k in ('kind','grism','band','selection','status','error','bins','shape')} for r in report['checks']],
                      'total_response_bytes':sum(r['bytes'] for r in session.requests)},indent=2),flush=True)
    if any(r['status']!='ok' for r in report['checks']):raise SystemExit(1)


if __name__=='__main__':main()
