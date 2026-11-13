"""DR1 single-HDU SpectrumDM retrieval, following IRSA's official tutorial."""
from __future__ import annotations

import hashlib
import io
from pathlib import PurePosixPath
from urllib.parse import quote, urlencode

import astropy.units as u
import numpy as np
import requests
from astropy.table import QTable
from astropy.io.votable import parse
from astroquery.ipac.irsa import Irsa

ASSOCIATION_TABLE = "euclid.objectid_spectrafile_association_dr1"
SPECTRUMDM_URL = "https://irsacloud.ipac.caltech.edu/api/spectrumdm/convert/"
DATASET_ID = "euclid_combspec_dr1"
MAX_ASSOCIATIONS = 100
MAX_BYTES = 8 * 1024 * 1024
FLUX_UNIT = u.erg / (u.s * u.cm**2 * u.AA)


def _object_id(value):
    text = str(value)
    digits = text[1:] if text.startswith("-") else text
    if not digits or not digits.isascii() or not digits.isdecimal():
        raise ValueError("object_id must be a decimal integer with an optional leading minus sign")
    number = int(text)
    if not -(2**63) <= number < 2**63:
        raise ValueError("object_id is outside the signed 64-bit range")
    return str(number)


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
    return result


def retrieve_spectra(object_id, *, grism=None, product=None, client=None, session=None,
                     timeout=120, max_bytes=MAX_BYTES):
    """Return all DR1 products, optionally filtered; fail on any download error.

    ``product=(path, hdu)`` must identify an association for this object. ``grism``
    is 'red' or 'blue', using the wavelength rule in IRSA's DR1 tutorial. Neither
    filter ranks spectra, combines products, clips variance, or removes bins.
    Returns [] only for no associations or no matching products. Raises for
    malformed metadata, transport errors, response limits, or invalid spectra.
    """
    object_id = _object_id(object_id)
    if grism not in (None, "red", "blue"):
        raise ValueError("grism must be 'red', 'blue', or None")
    if product is not None:
        if not isinstance(product, (tuple, list)) or len(product) != 2:
            raise ValueError("product must be a (path, hdu) pair")
        product = (str(product[0]), _hdu(product[1]))
        spectrum_url(*product)
    client = Irsa if client is None else client
    query = (f"SELECT TOP {MAX_ASSOCIATIONS + 1} objectid,tileid,path,hdu "
             f"FROM {ASSOCIATION_TABLE} WHERE objectid = {object_id}")
    rows = client.query_tap(query, maxrec=MAX_ASSOCIATIONS + 1).to_table()
    if len(rows) > MAX_ASSOCIATIONS:
        raise ValueError("Too many spectrum associations; refusing a truncated product list")
    result = []
    for row in rows:
        if _object_id(row["objectid"]) != object_id:
            raise ValueError("Association table returned a different object ID")
        path, hdu = str(row["path"]), _hdu(row["hdu"])
        url = spectrum_url(path, hdu)
        if product is not None and (path, hdu) != product:
            continue
        body = _download(url, session=session, timeout=timeout, max_bytes=max_bytes)
        spectrum = decode_spectrum(body, object_id=object_id, path=path, hdu=hdu, tileid=row["tileid"])
        if grism is None or spectrum.meta["GRISM"] == grism:
            result.append(spectrum)
    return result
