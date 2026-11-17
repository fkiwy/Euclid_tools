"""Small live source-patch checks, retaining routes and bounded response sizes."""
import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
import time

from euclid_tools.backends import get_backend
from euclid_tools.provenance import product_provenance
from validate_esa_dr1_products import AuditSession


def main():
    cli = argparse.ArgumentParser(__doc__)
    cli.add_argument('--output', type=Path, required=True)
    cli.add_argument('--timeout', type=float, default=60)
    args = cli.parse_args()
    backend = get_backend('esa'); session = AuditSession()
    report = {'checked_at': datetime.now(timezone.utc).isoformat(),
              'sources': ['https://docs.astropy.org/projects/astroquery/en/latest/api/astroquery.esa.euclid.EuclidClass.html',
                          'https://eas.esac.esa.int/tap-server/tap/sync'],
              'checks': [], 'requests': session.requests}
    for identity in ('-629553356493020915', '2714057875648364582'):
        rows = backend.discover_spectra(identity, grism='red', session=session, timeout=args.timeout)
        print(identity, 'associations', len(rows), flush=True)
        for row in rows:
            item = {'source_id': identity, 'source_path': str(row['source_path']),
                    'source_hdu': int(row['hdu_index']), 'linking_id': str(row['linking_id'])}
            start = time.monotonic()
            try:
                if not row['retrieval_supported']:
                    # Exercise the blocker with the exact selector; no product request is permitted.
                    count = sum(r['url'].endswith('/data') for r in session.requests)
                    try:
                        backend.retrieve_spectrum(identity, release='dr1', grism='red',
                            product=(row['source_path'], int(row['hdu_index'])),
                            session=session, timeout=args.timeout)
                    except ValueError as exc:
                        assert 'no unique native sourcepatch_id' in str(exc)
                        assert sum(r['url'].endswith('/data') for r in session.requests) == count
                        item.update(status='blocked_as_expected', reason=str(exc))
                    else:
                        raise AssertionError('Unenforceable parent selection was not blocked')
                else:
                    result = backend.retrieve_spectrum(identity, release='dr1', grism='red',
                        product=(row['source_path'], int(row['hdu_index'])), session=session, timeout=args.timeout)
                    assert result is not None and len(result) == 531
                    assert result.meta['LINKPAR'] == 'SOURCEPATCH_ID'
                    assert result.meta['LINKID'] == str(row['linking_id'])
                    item.update(status='ok', bins=len(result), provenance=product_provenance(result))
            except Exception as exc:
                item.update(status='error', error=repr(exc))
            item['seconds'] = time.monotonic()-start
            report['checks'].append(item)
            print(identity, item['source_hdu'], item['status'], flush=True)
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(json.dumps(report, indent=2)+'\n')
    report['total_response_bytes'] = sum(r['bytes'] for r in session.requests)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    if any(c['status'] == 'error' for c in report['checks']):
        raise SystemExit(1)


if __name__ == '__main__':
    main()
