import os
import tempfile
import matplotlib.pyplot as plt

import tools.shared as shr


def plot_spectrum(data, ra, dec, output_dir=tempfile.gettempdir(), open_plot=True, plot_format="pdf"):
    object_name = shr.create_object_name(ra, dec, precision=2, shortform=False, prefix="J", decimal=False)
    filename = os.path.join(output_dir, object_name + "_spectrum." + plot_format)

    wavelength = data["WAVELENGTH"]
    flux = data["FLUX"]
    error = data["ERROR"]

    plt.rcParams.update({"font.family": "Arial"})
    plt.plot(wavelength, flux, color="black", label="Spectrum")
    plt.plot(wavelength, error, color="red", label="Error")
    plt.xlabel(f"Wavelength [{to_latex(wavelength.unit)}]")
    plt.ylabel(f"Flux [{to_latex(flux.unit)}]")
    plt.legend(loc="best")
    plt.title(object_name)
    plt.savefig(filename, dpi=300, bbox_inches="tight", format=plot_format)
    plt.close()

    if open_plot:
        shr.open_file(filename)


def to_latex(unit):
    return unit.to_string("latex_inline")
