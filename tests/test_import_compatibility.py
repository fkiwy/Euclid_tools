"""Check both installed-style imports and the legacy checkout fallback."""

from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


class ImportCompatibilityTests(unittest.TestCase):
    def test_public_imports_from_checkout_and_src(self):
        # Separate interpreters prevent sys.modules from hiding loader mistakes.
        script = r'''
from pathlib import Path
import pickle
import sys
from importlib.resources import files

sys.path.insert(0, sys.argv[1])
import euclid_tools
import flux_comp
from euclid_tools import image_plotter, shared
from euclid_tools.backends import ESABackend, get_backend
from flux_comp.core import SED, WaveFlux, TemplateProvider
import flux_comp.uhs

source = Path(sys.argv[2])
assert Path(euclid_tools.__file__).resolve() == source / 'euclid_tools' / '__init__.py'
assert Path(flux_comp.__file__).resolve() == source / 'flux_comp' / '__init__.py'
assert euclid_tools.ESABackend is ESABackend
assert type(get_backend('esa')) is ESABackend
assert pickle.loads(pickle.dumps(SED)) is SED
assert pickle.loads(pickle.dumps(WaveFlux)) is WaveFlux
assert shared is sys.modules['euclid_tools.shared']
templates = Path(TemplateProvider().module_path) / 'templates'
assert Path(str(files('flux_comp'))) == source / 'flux_comp'
for family in ('Burgasser+2017', 'Kesseli+2017', 'Theissen+2022'):
    assert list((templates / family).glob('*.fits')), family
'''
        with tempfile.TemporaryDirectory() as workdir:
            for import_root in (ROOT, ROOT / "src"):
                with self.subTest(import_root=import_root):
                    result = subprocess.run(
                        [sys.executable, "-I", "-c", script, str(import_root), str(ROOT / "src")],
                        cwd=workdir,
                        capture_output=True,
                        text=True,
                        timeout=120,
                    )
                    self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
