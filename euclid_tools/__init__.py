"""Compatibility loader for imports from an uninstalled checkout.

Installed distributions load the package from src/ directly. This loader keeps
historical checkout imports working with the same module and class identities.
"""

from importlib.util import spec_from_file_location as _spec_from_file_location
from pathlib import Path as _Path
import sys as _sys

_package_dir = _Path(__file__).resolve().parent.parent / "src" / __name__
_spec = _spec_from_file_location(
    __name__, _package_dir / "__init__.py",
    submodule_search_locations=[str(_package_dir)],
)
_module = _sys.modules[__name__]
_module.__file__ = _spec.origin
_module.__path__ = _spec.submodule_search_locations
_module.__spec__ = _spec
_module.__loader__ = _spec.loader
_spec.loader.exec_module(_module)
del _spec_from_file_location, _Path, _sys, _package_dir, _spec, _module
