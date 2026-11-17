"""Native source-patch routing must enforce each exact parent association."""
import io

import numpy as np
import pytest
from astropy.io import fits
from astropy.table import QTable

from euclid_tools import esa_products as esa
from euclid_tools.backends import get_backend
from euclid_tools.provenance import product_provenance
from test_esa_dr1_products import association, client, Session, spectrum_body, ID, PATH, DATA


def multiple(patches=None):
    rows = association(2)
    rows['hdu_index'][1] = 13
    rows['sourcepatch_id'] = patches or [f'({ID}164)', f'({ID}171)']
    return rows


def test_discovery_retains_large_patch_identifiers_and_exact_product_pairs():
    session = Session()
    found = get_backend('esa').discover_spectra(ID, release='dr1', grism='red',
                                              client=client(multiple()), session=session)
    assert found['linking_id'].tolist() == [ID+'164', ID+'171']
    assert found['hdu_index'].tolist() == [12, 13]
    assert found['source_path'].tolist() == [PATH, PATH]
    assert found['retrieval_supported'].tolist() == [True, True]
    assert not session.calls


def test_singular_requires_selection_before_any_download():
    session = Session()
    with pytest.raises(ValueError, match='product='):
        esa.retrieve_spectrum(ID, grism='red', client=client(multiple()), session=session)
    assert not session.calls


def test_exact_selection_uses_native_patch_and_survives_fits(tmp_path):
    session = Session(spectrum_body())
    result = get_backend('esa').retrieve_spectrum(ID, release='dr1', product=(PATH, 13),
                                                client=client(multiple()), session=session)
    params = session.calls[0][1]['data']
    assert params['ID'] == 'dr1 '+ID+'171'
    assert params['LINKING_PARAMETER'] == 'SOURCEPATCH_ID'
    assert params['RETRIEVAL_TYPE'] == 'SPECTRA_RGS'
    saved = tmp_path/'selected.fits'; result.write(saved)
    record = product_provenance(QTable.read(saved))
    assert record['source_id'] == ID and record['source_hdu'] == 13
    assert record['source_file'] == PATH
    assert record['linking_id'] == ID+'171'
    assert record['linking_parameter'] == 'SOURCEPATCH_ID'
    assert record['retrieval_route'] == 'tap-sourcepatch+esa-datalink'


def test_plural_returns_each_parent_without_combining_science_arrays():
    raw = QTable.read(DATA/'esa_dr1_spectrum_excerpt.ecsv', format='ascii.ecsv')
    raw['SIGNAL'] *= 2
    stream = io.BytesIO(); raw.write(stream, format='fits')
    session = Session(spectrum_body(), stream.getvalue())
    results = esa.retrieve_spectra(ID, grism='red', client=client(multiple()), session=session)
    assert len(results) == 2 and [r.meta['SRCHDU'] for r in results] == [12, 13]
    assert all(len(r) == 5 for r in results)
    np.testing.assert_allclose(results[1]['FLUX'].value, 2*results[0]['FLUX'].value)
    assert [c[1]['data']['ID'] for c in session.calls] == ['dr1 '+ID+'164', 'dr1 '+ID+'171']


def test_shared_patch_is_blocked_but_other_parent_can_use_unique_alias():
    rows = multiple([f'({ID}164, {ID}171)', f'({ID}171)'])
    found = esa.discover_spectra(ID, client=client(rows))
    assert found['retrieval_supported'].tolist() == [True, False]
    assert found['linking_id'][0] == ID+'164'
    for product in (None, (PATH, 13)):
        session = Session()
        with pytest.raises(ValueError, match='no unique native sourcepatch_id'):
            esa.retrieve_spectra(ID, product=product, client=client(rows), session=session)
        assert not session.calls
    session = Session(spectrum_body())
    assert esa.retrieve_spectrum(ID, product=(PATH, 12), client=client(rows), session=session) is not None


def test_exact_duplicate_associations_are_one_product_on_legacy_route():
    session = Session(spectrum_body())
    result = esa.retrieve_spectrum(ID, client=client(association(2)), session=session)
    assert result.meta['LINKPAR'] == 'SOURCE_ID' and len(session.calls) == 1


def test_unknown_product_and_wrong_source_fail_explicitly():
    session = Session()
    with pytest.raises(ValueError, match='not associated'):
        esa.retrieve_spectrum(ID, product=(PATH, 99), client=client(multiple()), session=session)
    assert not session.calls
    with pytest.raises(ValueError, match='source ID'):
        esa.retrieve_spectrum(ID, product=(PATH, 12), client=client(multiple()),
                              session=Session(spectrum_body(wrong_id=True)))


def test_multi_table_patch_response_cannot_be_mislabeled_as_one_parent():
    with fits.open(io.BytesIO(spectrum_body(zipped=False))) as raw:
        hdus = fits.HDUList([raw[0].copy(), raw[1].copy(), raw[1].copy()])
        stream = io.BytesIO(); hdus.writeto(stream)
    with pytest.raises(ValueError, match='one ESA extracted spectrum'):
        esa.retrieve_spectrum(ID, product=(PATH, 12), client=client(multiple()),
                              session=Session(stream.getvalue()))


@pytest.mark.parametrize('value', ['(1.1)', '(1e20)', '(1; SELECT)', '(1,,2)'])
def test_invalid_native_patch_ids_fail_before_download(value):
    session = Session()
    with pytest.raises(ValueError, match='sourcepatch_id'):
        esa.retrieve_spectra(ID, client=client(multiple([value, '(2)'])), session=session)
    assert not session.calls


def test_negative_ids_and_native_arrays_remain_exact():
    assert esa._source_patch_ids('(-629553356493020915115, -62955335649302091553)') == (
        '-629553356493020915115', '-62955335649302091553')
    assert esa._source_patch_ids(np.array(['123456789012345678901', '2'])) == ('123456789012345678901', '2')


def test_empty_discovery_and_q1_discovery_contract():
    assert len(esa.discover_spectra(ID, client=client(association(0)))) == 0
    with pytest.raises(ValueError, match="release='dr1'"):
        get_backend('esa').discover_spectra(ID, release='q1')
