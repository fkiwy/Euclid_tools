"""Explicit IRSA DR1 MER mosaic discovery and bounded S3 FITS section reads."""
from __future__ import annotations

import hashlib
import io
import json
import re
import warnings
from pathlib import PurePosixPath

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from astroquery.ipac.irsa import Irsa
from pyvo.dal.exceptions import DALOverflowWarning

COLLECTIONS = {tiling: f"euclid_dr1_{tiling}_DpdMerBksMosaic" for tiling in ("wide", "deep")}
BUCKET = "nasa-irsa-euclid-dr1"
MAX_METADATA = 5001
MAX_CANDIDATES = 100
MAX_PIXELS = 1_000_000
MAX_BYTES = 16 * 1024 * 1024
_COLUMNS = ("product", "obs_id", "tile_id", "band", "tiling", "kind", "status", "error",
            "dataset", "processing_plan", "pipeline", "pipeline_version", "bunit", "magzero", "header_sha256")


class CutoutReadLimitError(RuntimeError):
    """Raised before an S3 read would exceed the cutout transfer budget."""


class _BoundedReader(io.RawIOBase):
    """Reject unbounded/oversize reads before S3 transfer; no block read-ahead."""
    def __init__(self, stream, max_bytes):
        self.stream = stream
        self.max_bytes = max_bytes
        self.bytes_read = 0
        self.ranges = []

    @property
    def name(self):
        return getattr(self.stream, "name", "mosaic.fits")

    @property
    def mode(self):
        return "rb"

    def readable(self):
        return True

    def seekable(self):
        return True

    def seek(self, offset, whence=0):
        return self.stream.seek(offset, whence)

    def tell(self):
        return self.stream.tell()

    def read(self, size=-1):
        if size < 0 or self.bytes_read + size > self.max_bytes:
            raise CutoutReadLimitError("S3 section read exceeds byte budget; no full-mosaic fallback")
        start = self.tell()
        data = self.stream.read(size)
        self.bytes_read += len(data)
        self.ranges.append((start, start + len(data)))
        return data

    def readinto(self, buffer):
        data = self.read(len(buffer))
        buffer[:len(data)] = data
        return len(data)


def _filesystem(filesystem):
    if filesystem is not None:
        return filesystem
    import s3fs
    return s3fs.S3FileSystem(anon=True, skip_instance_cache=True,
                           config_kwargs={"connect_timeout": 15, "read_timeout": 60, "retries": {"max_attempts": 1}})


def _arguments(ra, dec, radius, band, tiling, tile_id, product):
    if band not in ("VIS", "Y", "J", "H") or tiling not in COLLECTIONS:
        raise ValueError("Choose band VIS/Y/J/H and tiling wide/deep")
    if not np.isfinite(radius) or radius <= 0:
        raise ValueError("Search radius must be positive arcseconds")
    coordinate = SkyCoord(float(ra) * u.deg, float(dec) * u.deg, frame="icrs")
    if not np.isfinite(coordinate.ra.deg) or not np.isfinite(coordinate.dec.deg):
        raise ValueError("Coordinates must be finite degrees")
    if tile_id is not None and not re.fullmatch(r"10[12]\d{6}", str(tile_id)):
        raise ValueError("tile_id must be a nine-digit MER tile ID")
    if tile_id is not None and not str(tile_id).startswith("102" if tiling == "wide" else "101"):
        raise ValueError("tile_id prefix disagrees with the requested wide/deep tiling")
    if product is not None and not str(product).startswith(f"s3://{BUCKET}/dr1/MER/"):
        raise ValueError("product must be an exact DR1 S3 URI returned by discovery")
    return coordinate


def _cloud_path(row, tiling):
    aws = json.loads(str(row["cloud_access"]))["aws"]
    bucket, key = str(aws["bucket_name"]), str(aws["key"])
    parts = PurePosixPath(key).parts
    tile = str(row["obs_id"])[:9]
    if (bucket != BUCKET or len(parts) < 5 or parts[:2] != ("dr1", "MER")
            or parts[2] != tile or parts[3] != ("VIS" if str(row["energy_bandpassname"]) == "VIS" else "NISP")
            or ".." in parts or not key.endswith(".fits")
            or not re.fullmatch(r"102\d{6}" if tiling == "wide" else r"101\d{6}", tile)):
        raise ValueError("SIA cloud path/tile disagrees with the requested DR1 collection")
    if str(aws.get("region", "us-east-1")) != "us-east-1":
        raise ValueError("Unexpected DR1 bucket region")
    return f"s3://{bucket}/{key}", tile


def _metadata(header, band, tiling):
    for key in ("DATASETR", "PPLANID", "PIPDEFID", "PSWREL", "FILTER", "BUNIT", "MAGZERO"):
        if key not in header or not str(header[key]).strip():
            raise ValueError(f"Unverified MER header: missing {key}")
    if header["DATASETR"].strip() not in ("DR1_R1", "DR1_R2"):
        raise ValueError("Unverified DATASETR; expected verified DR1_R1/DR1_R2")
    # FILTER values verified from live VIS and NISP DR1 headers.
    expected_filter = "VIS" if band == "VIS" else "NIR_" + band
    if str(header["FILTER"]).strip() != expected_filter:
        raise ValueError("FITS FILTER disagrees with the selected SIA band")
    if not np.isfinite(float(header["MAGZERO"])):
        raise ValueError("Nonfinite MAGZERO")
    plan = str(header["PPLANID"])
    if not str(header["PIPDEFID"]).startswith("MER_ProcessTile") or not plan.startswith(str(header["PIPDEFID"])):
        raise ValueError("Processing plan and MER pipeline identifier disagree")
    if tiling == "wide" and "-DATA_RELEASE-" in plan:
        kind = "survey"
    elif "DATA_RELEASE_DEEP" in plan:
        kind = "repeat" if tiling == "wide" else "deep"
    else:
        raise ValueError("Unverified processing plan; cannot classify survey/repeat/deep")
    if header.get("NAXIS") != 2 or min(header.get("NAXIS1", 0), header.get("NAXIS2", 0)) <= 0:
        raise ValueError("Expected a two-dimensional primary MER image")
    if header.get("BITPIX") != -32 or header.get("BSCALE", 1) != 1 or header.get("BZERO", 0) != 0:
        raise ValueError("Unverified MER pixel encoding; require unscaled float32")
    wcs = WCS(header)
    if not wcs.has_celestial or wcs.pixel_n_dim != 2 or wcs.has_distortion:
        raise ValueError("Unverified MER WCS: require a 2D celestial WCS without distortion tables")
    if any(re.match(r"CTYPE\d[A-Z]$", key) for key in header):
        raise ValueError("Alternate WCS frames require separate cutout verification")
    scales = proj_plane_pixel_scales(wcs)
    if not np.all(np.isfinite(scales)) or np.any(scales <= 0):
        raise ValueError("Invalid WCS pixel scales")
    return {"kind": kind, "dataset": header["DATASETR"], "processing_plan": plan,
            "pipeline": header["PIPDEFID"], "pipeline_version": str(header["PSWREL"]),
            "bunit": header["BUNIT"], "magzero": float(header["MAGZERO"]),
            "header_sha256": hashlib.sha256(header.tostring().encode("ascii")).hexdigest()}


def discover_mosaics(ra, dec, *, radius=2, band="VIS", tiling="wide", tile_id=None,
                     product=None, client=None, filesystem=None):
    """List all matching science mosaics, including missing/invalid candidates.

    Reads headers only. No ranking, nearest-tile choice or collection fallback.
    Rows with status='error' retain the cause; callers must resolve them explicitly.
    """
    coordinate = _arguments(ra, dec, radius, band, tiling, tile_id, product)
    client = Irsa if client is None else client
    with warnings.catch_warnings():
        warnings.simplefilter("error", DALOverflowWarning)
        rows = client.query_sia(pos=(coordinate, radius * u.arcsec), collection=COLLECTIONS[tiling],
                                maxrec=MAX_METADATA)
    if len(rows) >= MAX_METADATA:
        raise ValueError("SIA metadata limit reached; refusing a possibly incomplete mosaic list")
    if len(rows) == 0:
        return Table(names=_COLUMNS, dtype=[str] * len(_COLUMNS))
    matches = rows[(rows["dataproduct_subtype"] == "science") & (rows["energy_bandpassname"] == band)]
    if tile_id is not None:
        matches = matches[[str(row["obs_id"]).startswith(str(tile_id) + "_") for row in matches]]
    if len(matches) > MAX_CANDIDATES:
        raise ValueError("Too many mosaic candidates; specify a tile or exact product")
    fs = _filesystem(filesystem) if len(matches) else filesystem
    records = []
    for row in matches:
        record = {key: "" for key in _COLUMNS}
        record.update(obs_id=str(row["obs_id"]), band=band, tiling=tiling, status="error", magzero=np.nan)
        try:
            uri, tile = _cloud_path(row, tiling)
            if product is not None and uri != product:
                continue
            record.update(product=uri, tile_id=tile)
            with fs.open(uri, "rb", block_size=28800, cache_type="none") as stream:
                reader = _BoundedReader(stream, 64 * 1024)
                header = fits.Header.fromfile(reader)
            record.update(_metadata(header, band, tiling), status="ready")
        except Exception as exc:
            record["error"] = f"{type(exc).__name__}: {exc}"
        records.append(record)
    return Table(rows=records, names=_COLUMNS) if records else Table(names=_COLUMNS, dtype=[str] * len(_COLUMNS))


_WCS_KEY = re.compile(r"^(WCSAXES|CRPIX\d+|CRVAL\d+|CDELT\d+|CTYPE\d+|CUNIT\d+|CROTA\d+|"
                      r"CD\d+_\d+|PC\d+_\d+|PV\d+_\d+|PS\d+_\d+|LONPOLE|LATPOLE|RADESYS|EQUINOX)$")


def retrieve_cutout(ra, dec, *, radius=2, cutout_size=20, band="VIS", tiling="wide",
                    selection=None, tile_id=None, product=None, mode="trim", client=None,
                    filesystem=None, max_bytes=MAX_BYTES):
    """Read only one unambiguous, verified DR1 MER image section from S3."""
    coordinate = _arguments(ra, dec, radius, band, tiling, tile_id, product)
    selection = ("survey" if tiling == "wide" else "deep") if selection is None else selection
    if selection not in (("survey", "repeat") if tiling == "wide" else ("deep",)):
        raise ValueError("wide selection must be survey/repeat; deep selection must be deep")
    if not np.isfinite(cutout_size) or cutout_size <= 0 or mode not in ("trim", "strict"):
        raise ValueError("Need a positive cutout_size in arcsec and mode trim/strict")
    if not isinstance(max_bytes, (int, np.integer)) or max_bytes <= 0:
        raise ValueError("max_bytes must be a positive integer")
    fs = _filesystem(filesystem)
    found = discover_mosaics(ra, dec, radius=radius, band=band, tiling=tiling, tile_id=tile_id,
                             product=product, client=client, filesystem=fs)
    if any(found["status"] == "error"):
        raise ValueError("Mosaic candidates could not be verified; inspect discover_mosaics and select an exact "
                         "available product: " + "; ".join(str(x) for x in found["error"] if str(x)))
    found = found[found["kind"] == selection]
    if len(found) == 0:
        return None
    if len(found) != 1:
        raise ValueError(f"Ambiguous DR1 mosaic: {len(found)} matches; select tile_id or exact product")
    row = found[0]
    with fs.open(str(row["product"]), "rb", block_size=28800, cache_type="none") as stream:
        reader = _BoundedReader(stream, max_bytes)
        with fits.open(reader, mode="readonly", memmap=False, lazy_load_hdus=True) as hdus:
            source_header = hdus[0].header.copy()
            metadata = _metadata(source_header, band, tiling)
            if metadata["header_sha256"] != row["header_sha256"]:
                raise ValueError("Mosaic header changed since discovery")
            wcs = WCS(source_header)
            scales = proj_plane_pixel_scales(wcs) * 3600
            pixels = np.round(cutout_size / scales).astype(int)
            if np.any(pixels < 1) or np.prod(pixels.astype(float)) > MAX_PIXELS:
                raise ValueError("Requested cutout exceeds pixel budget or is smaller than one pixel")
            # .section performs a sliced read; .data would fetch the entire mosaic.
            cutout = Cutout2D(hdus[0].section, position=coordinate, size=cutout_size * u.arcsec,
                              wcs=wcs, mode=mode, copy=True)
        bytes_read = reader.bytes_read
    header = source_header.copy()
    for key in list(header):
        if _WCS_KEY.match(key) or key in ("CHECKSUM", "DATASUM"):
            del header[key]
    header.update(cutout.wcs.to_header(relax=True))
    x0, y0 = cutout.origin_original
    header.update(RELEASE="DR1", ARCHIVE="IRSA", SRCS3=str(row["product"]), TILEID=str(row["tile_id"]),
                  SIAOBS=str(row["obs_id"]), SIACOLL=COLLECTIONS[tiling], MOSCLASS=selection,
                  RAWHDR=str(row["header_sha256"]), CUTRA=coordinate.ra.deg, CUTDEC=coordinate.dec.deg,
                  CUTSIZE=float(cutout_size), CUTMODE=mode, CUTTRIM=tuple(cutout.shape) != tuple(cutout.shape_input),
                  ORIGNX=source_header["NAXIS1"], ORIGNY=source_header["NAXIS2"],
                  SRCX0=int(x0), SRCY0=int(y0), CUTBYTES=bytes_read)
    from euclid_tools.provenance import stamp
    stamp(header, product_type="image", release="DR1", archive="IRSA",
          route="sia+s3-section", SRCPATH=str(row["product"]), SRCHDU=0,
          SRCURL=str(row["product"]), PRODBAND=band)
    result = fits.PrimaryHDU(data=cutout.data, header=header)
    new_wcs = WCS(result.header)
    x = np.array([0., (result.data.shape[1] - 1) / 2, result.data.shape[1] - 1.])
    y = np.array([0., (result.data.shape[0] - 1) / 2, result.data.shape[0] - 1.])
    separation = new_wcs.pixel_to_world(x, y).separation(wcs.pixel_to_world(x + x0, y + y0)).arcsec
    if not np.all(np.isfinite(separation)) or np.any(separation > 1e-4):
        raise ValueError("Cutout WCS failed source-pixel roundtrip verification")
    return result
