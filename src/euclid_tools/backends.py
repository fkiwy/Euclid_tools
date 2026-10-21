from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Optional


def resolve_radius(radius: Optional[float] = None, search_radius: Optional[float] = None) -> float:
    """Resolve the archive API radius argument standardized across backends."""
    if radius is None and search_radius is None:
        raise TypeError("Either 'radius' or 'search_radius' must be provided")
    if radius is not None and search_radius is not None and radius != search_radius:
        raise ValueError("'radius' and 'search_radius' must agree when both are given")
    return radius if radius is not None else search_radius


class EuclidArchiveBackend(ABC):
    """Small backend contract used to keep archive-specific implementations aligned."""

    name: str = "base"

    @staticmethod
    def resolve_radius(radius=None, search_radius=None):
        return resolve_radius(radius=radius, search_radius=search_radius)

    @abstractmethod
    def retrieve_objects(self, ra, dec, radius=None, search_radius=None, **kwargs):
        raise NotImplementedError

    @abstractmethod
    def retrieve_spectrum(self, object_id, **kwargs):
        raise NotImplementedError

    @abstractmethod
    def retrieve_cutout(self, ra, dec, radius=None, search_radius=None, cutout_size=None, band=None, **kwargs):
        raise NotImplementedError


class ESABackend(EuclidArchiveBackend):
    name = "esa"

    def retrieve_objects(self, ra, dec, radius=None, search_radius=None, **kwargs):
        from euclid_tools import esa_tools

        return esa_tools.retrieve_objects(ra=ra, dec=dec, radius=resolve_radius(radius, search_radius), **kwargs)

    def retrieve_spectrum(self, object_id, **kwargs):
        from euclid_tools import esa_tools

        return esa_tools.retrieve_spectrum(object_id, **kwargs)

    def retrieve_cutout(self, ra, dec, radius=None, search_radius=None, cutout_size=None, band=None, **kwargs):
        from euclid_tools import esa_tools

        return esa_tools.retrieve_cutout(
            ra=ra,
            dec=dec,
            search_radius=resolve_radius(radius, search_radius),
            cutout_size=cutout_size,
            band=band,
            **kwargs,
        )


class IRSABackend(EuclidArchiveBackend):
    name = "irsa"

    def retrieve_objects(self, ra, dec, radius=None, search_radius=None, **kwargs):
        from euclid_tools import irsa_tools

        return irsa_tools.retrieve_objects(ra=ra, dec=dec, radius=resolve_radius(radius, search_radius), **kwargs)

    def retrieve_spectrum(self, object_id, **kwargs):
        from euclid_tools import irsa_tools

        return irsa_tools.retrieve_spectrum(object_id, **kwargs)

    def retrieve_cutout(self, ra, dec, radius=None, search_radius=None, cutout_size=None, band=None, **kwargs):
        from euclid_tools import irsa_tools

        return irsa_tools.retrieve_cutout(
            ra=ra,
            dec=dec,
            search_radius=resolve_radius(radius, search_radius),
            cutout_size=cutout_size,
            band=band,
            **kwargs,
        )


def get_backend(name: str) -> EuclidArchiveBackend:
    backend_name = (name or "").lower()
    if backend_name == "esa":
        return ESABackend()
    if backend_name == "irsa":
        return IRSABackend()
    raise ValueError(f"Unsupported backend: {name!r}. Choose from 'esa' or 'irsa'.")
