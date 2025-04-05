from euclid_tools.esa_tools import retrieve_objects
import euclid_tools.shared as shr


def print_magnitude(flux, flux_err, band):
    mag, mag_err = shr.convert_flux_to_mag(flux, flux_err, magnitude_system=shr.MagnitudeSystem.AB)
    print(f"{band} AB   mag = {mag} ± {mag_err} mag")

    mag, mag_err = shr.convert_flux_to_mag(flux, flux_err, magnitude_system=shr.MagnitudeSystem.Vega, band=band)
    print(f"{band} Vega mag = {mag} ± {mag_err} mag")


ra, dec = 266.4850113, 64.9936424
search_radius = 5

results = retrieve_objects(ra, dec, search_radius)

if results:
    result = results[0]
    print_magnitude(result["flux_vis_psf"], result["fluxerr_vis_psf"], "VIS")
    print_magnitude(result["flux_y_templfit"], result["fluxerr_y_templfit"], "Y")
    print_magnitude(result["flux_j_templfit"], result["fluxerr_j_templfit"], "J")
    print_magnitude(result["flux_h_templfit"], result["fluxerr_h_templfit"], "H")
