"""Run file-export entrypoints in clean processes configured to select TkAgg."""
import os
from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("archive", ["irsa", "esa"])
def test_file_examples_render_without_tk_before_any_archive_request(tmp_path, archive):
    script = r'''
import importlib
from pathlib import Path
import runpy
import sys
import matplotlib
matplotlib.use("TkAgg")

class ReachedArchive(Exception):
    pass

def archive_boundary(*args, **kwargs):
    assert matplotlib.get_backend().lower() == "agg", matplotlib.get_backend()
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.plot([1,2], [3,4])
    fig.savefig("backend-check.pdf")
    assert type(fig.canvas).__module__ == "matplotlib.backends.backend_agg"
    plt.close(fig)
    assert Path("backend-check.pdf").read_bytes().startswith(b"%PDF")
    assert "matplotlib.backends.backend_tkagg" not in sys.modules
    raise ReachedArchive

module = importlib.import_module("euclid_tools." + sys.argv[1] + "_tools")
module.print_catalog_info = archive_boundary
try:
    runpy.run_path(sys.argv[2], run_name="__main__")
except ReachedArchive:
    pass
else:
    raise AssertionError("Did not stop at the first archive request")
'''
    env = {**os.environ, "MPLBACKEND": "TkAgg", "MPLCONFIGDIR": str(tmp_path / "mpl"),
           "PYTHONPATH": str(ROOT / "src")}
    result = subprocess.run([sys.executable, "-c", script, archive, str(ROOT / f"{archive}_tools_example.py")],
                            cwd=tmp_path, env=env, text=True, capture_output=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr


def test_library_imports_preserve_callers_backend(tmp_path):
    script = '''
import matplotlib
matplotlib.use("svg")
import euclid_tools.image_plotter
import flux_comp.core
assert matplotlib.get_backend().lower() == "svg"
'''
    result = subprocess.run([sys.executable,"-c",script], cwd=tmp_path,
                            env={**os.environ,"PYTHONPATH":str(ROOT / "src"),"MPLCONFIGDIR":str(tmp_path / "mpl")},
                            capture_output=True,text=True,timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
