import sys
import os

sys.path.append(os.path.abspath(os.path.join("..", "src")))

from object_finder import find_objects
from object_finder import convert_flux_to_mag
from object_finder import MagnitudeSystem


def print_magnitude(flux, flux_err, band):
    mag, mag_err = convert_flux_to_mag(flux, flux_err, magnitude_system=MagnitudeSystem.AB)
    print(f"{band} AB   mag = {mag} ± {mag_err} mag")

    mag, mag_err = convert_flux_to_mag(flux, flux_err, magnitude_system=MagnitudeSystem.Vega, band=band)
    print(f"{band} Vega mag = {mag} ± {mag_err} mag")


ra, dec = 273.891148, 64.526926
search_radius = 5

results = find_objects(ra, dec, search_radius)

if results:
    result = results[0]
    print_magnitude(result["flux_vis_psf"], result["fluxerr_vis_psf"], "VIS")
    print_magnitude(result["flux_y_templfit"], result["fluxerr_y_templfit"], "Y")
    print_magnitude(result["flux_j_templfit"], result["fluxerr_j_templfit"], "J")
    print_magnitude(result["flux_h_templfit"], result["fluxerr_h_templfit"], "H")
