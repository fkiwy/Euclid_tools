import os
import tempfile
import matplotlib.pyplot as plt

import tools.shared as shr


def plot_spectrum(data, ra, dec, output_dir=tempfile.gettempdir(), open_plot=True, plot_format="pdf"):
    """
    Plot and save the spectrum of a source at given coordinates.

    This function creates a plot of the flux and associated error as a function of wavelength,
    using a QTable with columns "WAVELENGTH", "FLUX", and "ERROR". The plot includes axis labels,
    a legend, and a title based on the target's coordinates. The result is saved to a file and
    optionally opened after saving.

    Parameters:
        data (QTable): Table containing spectral data with columns:
            - "WAVELENGTH" (Quantity): Wavelength values (with units),
            - "FLUX" (Quantity or array): Flux values,
            - "ERROR" (Quantity or array): Uncertainty in flux.
        ra (float): Right Ascension of the object (in degrees).
        dec (float): Declination of the object (in degrees).
        output_dir (str, optional): Directory where the plot will be saved (default: system temp directory).
        open_plot (bool, optional): If True, opens the plot after saving (default: True).
        plot_format (str, optional): File format to save the plot (e.g., "pdf", "png"; default: "pdf").

    Output:
        Saves a spectrum plot named after the object's coordinates, e.g., "JHHMMSS+DDMMSS_spectrum.pdf".

    Notes:
        - Uses `tools.shared.create_object_name` to generate the filename and plot title.
        - Uses LaTeX formatting for units in axis labels.
        - Assumes the "WAVELENGTH", "FLUX", and "ERROR" columns are available and properly formatted.
    """

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
