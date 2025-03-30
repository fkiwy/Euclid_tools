import os
import tempfile
import astropy.units as u
import matplotlib.pyplot as plt

import tools.shared as shr


def plot_spectrum(data, ra, dec, output_dir=tempfile.gettempdir(), open_plot=True, plot_format="png"):
    object_name = shr.create_object_name(ra, dec, precision=2, shortform=False, prefix="J", decimal=False)
    filename = os.path.join(output_dir, object_name + "." + plot_format)

    plt.rcParams.update({"font.family": "Arial"})
    plt.plot(data["WAVELENGTH"].to(u.um), data["SIGNAL"])
    plt.xlabel("Wavelength [microns]")
    plt.ylabel("Flux [" + data["SIGNAL"].unit.to_string("latex_inline") + "]")
    plt.title(object_name)
    plt.savefig(filename, dpi=300, bbox_inches="tight", format=plot_format)
    plt.close()

    if open_plot:
        shr.open_file(filename)
