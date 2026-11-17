"""DR1 single-HDU SpectrumDM retrieval, following IRSA's official tutorial."""
from __future__ import annotations

import hashlib
import io
from pathlib import PurePosixPath
from urllib.parse import quote, urlencode

import astropy.units as u
import numpy as np
import requests
from astropy.table import QTable, Table
from astropy.io.votable import parse
from astroquery.ipac.irsa import Irsa
from euclid_tools.validation import object_id as _object_id, positive_scalar

ASSOCIATION_TABLE = "euclid.objectid_spectrafile_association_dr1"
SPECTRUMDM_URL = "https://irsacloud.ipac.caltech.edu/api/spectrumdm/convert/"
DATASET_ID = "euclid_combspec_dr1"
MAX_ASSOCIATIONS = 100
MAX_BYTES = 8 * 1024 * 1024
FLUX_UNIT = u.erg / (u.s * u.cm**2 * u.AA)


def _hdu(value):
    text = str(value)
    if not text.isascii() or not text.isdecimal():
        raise ValueError("hdu must be a nonnegative integer")
    return int(text)


def spectrum_url(path, hdu):
    """Build only a DR1 SpectrumDM URL; never fetch a parent FITS file."""
    path = str(path)
    if (not path.startswith("euclid/dr1/") or ".." in PurePosixPath(path).parts
            or not path.endswith(".fits")):
        raise ValueError("Expected an archive path under euclid/dr1/ ending in .fits")
    return SPECTRUMDM_URL + quote(path, safe="/") + "?" + urlencode(
        {"dataset_id": DATASET_ID, "hdu": _hdu(hdu)})


def _download(url, *, session=None, timeout=120, max_bytes=MAX_BYTES):
    if not np.isfinite(timeout) or timeout <= 0 or max_bytes <= 0:
        raise ValueError("timeout and max_bytes must be positive")
    http = requests if session is None else session
    with http.get(url, timeout=timeout, stream=True) as response:
        response.raise_for_status()
        size = response.headers.get("Content-Length")
        if size is not None and int(size) > max_bytes:
            raise ValueError("SpectrumDM response exceeds max_bytes")
        body = bytearray()
        for chunk in response.iter_content(chunk_size=65536):
            body.extend(chunk)
            if len(body) > max_bytes:
                raise ValueError("SpectrumDM response exceeds max_bytes")
    return bytes(body)


def decode_spectrum(body, *, object_id, path, hdu, tileid):
    """Decode supplied units once. Retain all bins and diagnostic mask values."""
    document = parse(io.BytesIO(body))
    pending = [document]
    while pending:
        node = pending.pop()
        for info in getattr(node, "infos", []):
            if str(info.name).upper() == "QUERY_STATUS" and str(info.value).upper() != "OK":
                raise ValueError(f"SpectrumDM QUERY_STATUS={info.value}: {info.content}")
        pending.extend(getattr(node, "resources", []))
    raw = QTable.read(io.BytesIO(body), format="votable")
    missing = {"WAVELENGTH", "SIGNAL", "MASK"} - set(raw.colnames)
    if missing or not ({"VAR", "UNCERTAINTY"} & set(raw.colnames)):
        raise ValueError(f"SpectrumDM response lacks spectral columns: {sorted(missing)} or error/variance")
    if len(raw) == 0:
        raise ValueError("SpectrumDM returned an empty spectrum")
    # Conflicting identity metadata must not be overwritten with requested values.
    expected_id = f"Euclid_{_object_id(object_id)}_Combined"
    for key in ("ID", "name"):
        if key in raw.meta and str(raw.meta[key]) != expected_id:
            raise ValueError(f"SpectrumDM {key} does not match the requested object")
    if "RELEASE" in raw.meta and str(raw.meta["RELEASE"]).upper() != "DR1":
        raise ValueError("SpectrumDM returned a conflicting release")
    wave, flux = raw["WAVELENGTH"], raw["SIGNAL"]
    if not isinstance(wave, u.Quantity) or not isinstance(flux, u.Quantity):
        raise ValueError("SpectrumDM wavelength and flux require physical units")
    wave = wave.to(u.um)
    if not flux.unit.is_equivalent(FLUX_UNIT):
        raise ValueError("SpectrumDM SIGNAL must have spectral flux-density units per wavelength")
    if "VAR" in raw.colnames:
        variance = raw["VAR"]
        if not isinstance(variance, u.Quantity):
            raise ValueError("SpectrumDM VAR requires physical units")
        with np.errstate(invalid="ignore"):
            error = np.sqrt(variance.to(flux.unit**2))
        if "UNCERTAINTY" in raw.colnames:
            supplied = raw["UNCERTAINTY"]
            if not isinstance(supplied, u.Quantity):
                raise ValueError("SpectrumDM UNCERTAINTY requires physical units")
            supplied = supplied.to(flux.unit)
            good = np.isfinite(error.value) & np.isfinite(supplied.value)
            if not np.allclose(error.value[good], supplied.value[good], rtol=1e-5, atol=0):
                raise ValueError("SpectrumDM variance and uncertainty disagree")
    else:
        if not isinstance(raw["UNCERTAINTY"], u.Quantity):
            raise ValueError("SpectrumDM UNCERTAINTY requires physical units")
        error = raw["UNCERTAINTY"].to(flux.unit)
    result = QTable({"WAVELENGTH": wave, "FLUX": flux, "ERROR": error, "MASK": raw["MASK"]})
    for name in ("QUALITY", "NDITH"):
        if name in raw.colnames:
            result[name] = raw[name]
    # Preserve native scaled flux units. No FSCALE multiplication or mask filtering.
    finite_wave = np.ma.array(wave.value, mask=getattr(wave, "mask", False)).filled(np.nan)
    finite_wave = finite_wave[np.isfinite(finite_wave)]
    if not len(finite_wave) or np.any(finite_wave <= 0):
        raise ValueError("SpectrumDM has no valid positive wavelength coverage")
    grism = "blue" if finite_wave.min() < 1.1 else "red"
    result.meta.update(raw.meta)
    result.meta.update(RELEASE="DR1", OBJECTID=_object_id(object_id), ARCHIVE="IRSA",
                       DATASET=DATASET_ID, ASSOCTBL=ASSOCIATION_TABLE,
                       SRCPATH=str(path), SRCHDU=_hdu(hdu), TILEID=str(tileid), GRISM=grism,
                       GRISMREF="IRSA tutorial: minimum wavelength <1.1 um identifies blue",
                       SRCURL=spectrum_url(path, hdu), RAWSHA=hashlib.sha256(body).hexdigest())
    from euclid_tools.provenance import stamp
    stamp(result.meta, product_type="spectrum", release="DR1", archive="IRSA",
          route="tap-association+spectrumdm")
    return result


def _associations(object_id, *, client=None):
    object_id = _object_id(object_id)
    client = Irsa if client is None else client
    query = (f"SELECT TOP {MAX_ASSOCIATIONS + 1} objectid,tileid,path,hdu "
             f"FROM {ASSOCIATION_TABLE} WHERE objectid = {object_id}")
    job = client.query_tap(query, maxrec=MAX_ASSOCIATIONS + 1)
    if job is None:
        raise RuntimeError("IRSA spectrum association query failed")
    rows = job.to_table()
    if rows is None:
        raise RuntimeError("IRSA spectrum association query returned no result table")
    if len(rows) > MAX_ASSOCIATIONS:
        raise ValueError("Too many spectrum associations; refusing a truncated product list")
    found = {}
    for row in rows:
        if _object_id(row["objectid"]) != object_id:
            raise ValueError("Association table returned a different object ID")
        path, hdu = str(row["path"]), _hdu(row["hdu"])
        url = spectrum_url(path, hdu)
        key = (path, hdu)
        if key in found and found[key]['tileid'] != str(row['tileid']):
            raise ValueError("Conflicting IRSA tile identities for one file/HDU")
        found[key] = {'path': path, 'hdu': hdu, 'tileid': str(row['tileid']), 'url': url}
    return object_id, list(found.values())


def discover_spectra(object_id, *, grism=None, client=None):
    """List exact DR1 associations without downloading spectral data.

    IRSA's association table supplies no grism. Discovery leaves it unknown and
    rejects a grism filter; retrieve_spectra classifies supplied wavelengths.
    """
    if grism is not None:
        if grism not in ('red', 'blue'):
            raise ValueError("grism must be red, blue or None")
        raise NotImplementedError("IRSA association metadata has no grism; omit grism for discovery. "
                                  "retrieve_spectra can filter after reading wavelengths")
    _, rows = _associations(object_id, client=client)
    return Table(rows=[('', r['path'], r['hdu'], '', '', True, r['tileid']) for r in rows],
                 names=('grism', 'source_path', 'hdu_index', 'linking_parameter', 'linking_id',
                        'retrieval_supported', 'tile_id'), dtype=(str, str, int, str, str, bool, str))


def retrieve_spectra(object_id, *, grism=None, product=None, client=None, session=None,
                     timeout=120, max_bytes=MAX_BYTES):
    """Return every matching DR1 spectrum without ranking or combining products.

    An empty association/grism result returns []; an unassociated exact product
    for a source with associations raises ValueError. Metadata is validated and
    identical file/HDU rows are deduplicated before any spectrum download.
    """
    positive_scalar(timeout, 'timeout')
    positive_scalar(max_bytes, 'max_bytes')
    if grism not in (None, "red", "blue"):
        raise ValueError("grism must be 'red', 'blue', or None")
    if product is not None:
        if not isinstance(product, (tuple, list)) or len(product) != 2:
            raise ValueError("product must be a (path, hdu) pair")
        product = (str(product[0]), _hdu(product[1]))
        spectrum_url(*product)
    object_id, rows = _associations(object_id, client=client)
    selected = [r for r in rows if product is None or product == (r['path'], r['hdu'])]
    if product is not None and rows and not selected:
        raise ValueError("Selected IRSA spectrum product is not associated with this source")
    result = []
    for row in selected:
        body = _download(row['url'], session=session, timeout=timeout, max_bytes=max_bytes)
        spectrum = decode_spectrum(body, object_id=object_id, path=row['path'], hdu=row['hdu'], tileid=row['tileid'])
        if grism is None or spectrum.meta["GRISM"] == grism:
            result.append(spectrum)
    return result
