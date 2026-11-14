"""Versioned, JSON-safe product identity shared by spectra and image workflows.

Flat FITS cards remain the storage contract; missing metadata is unknown, never
inferred from a requested release, a filename, or a neighbouring product.
"""
from __future__ import annotations

import json

SCHEMA = "euclid.product-provenance/v1"
FIELDS = {
    "metadata_version": "PROVVER",
    "release": "RELEASE", "archive": "ARCHIVE", "product_type": "PRODTYPE",
    "source_id": "OBJECTID", "source_file": "SRCPATH", "source_hdu": "SRCHDU",
    "grism": "GRISM", "grism_reference": "GRISMREF", "retrieval_route": "RETROUTE",
    "retrieval_url": "SRCURL", "dataset": "DATASET", "association_table": "ASSOCTBL",
    "tile_id": "TILEID", "band": "PRODBAND", "mosaic_class": "MOSCLASS",
    "processing_dataset": "DATASETR", "processing_plan": "PPLANID",
    "processing_pipeline": "PIPDEFID", "processing_version": "PSWREL",
    "raw_sha256": "RAWSHA", "source_header_sha256": "RAWHDR",
    "cutout_header_sha256": "CUTHDR",
}


def stamp(metadata, *, product_type, release, archive, route, **cards):
    """Add identity cards without changing data, WCS, calibration or masks."""
    metadata.update(PROVVER=1, PRODTYPE=product_type, RELEASE=str(release).upper(),
                    ARCHIVE=str(archive).upper(), RETROUTE=route)
    for key, value in cards.items():
        if value is not None:
            metadata[key] = value
    return metadata


def product_provenance(product, *, source_id=None):
    """Return a stable record from a Table, HDU or metadata mapping.

    ``source_id`` is a workflow association (e.g. a cutout target), not evidence
    that the parent mosaic or spectrum metadata identifies that object.
    """
    metadata = getattr(product, "header", None)
    if metadata is None:
        metadata = getattr(product, "meta", product)
    result = {"schema": SCHEMA}
    for name, card in FIELDS.items():
        value = metadata.get(card)
        if hasattr(value, "item"):
            value = value.item()
        result[name] = value
    for key in ("source_id", "tile_id"):
        if result[key] is not None:
            result[key] = str(result[key])
    result["requested_source_id"] = None if source_id is None else str(source_id)
    if source_id is not None and result["source_id"] is not None and str(source_id) != result["source_id"]:
        raise ValueError("Product source ID differs from requested source ID")
    return result


def provenance_json(product, *, source_id=None):
    """Deterministic CSV/JSON representation; IDs retain exact decimal strings."""
    return json.dumps(product_provenance(product, source_id=source_id), sort_keys=True,
                      separators=(",", ":"), allow_nan=False)
