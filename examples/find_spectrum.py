import sys
import os

sys.path.append(os.path.abspath(os.path.join("..", "src")))

from object_finder import find_objects
from spectrum_finder import find_spectrum
from spectrum_plotter import plot_spectrum


ra, dec = 273.891148, 64.526926
search_radius = 5

results = find_objects(ra, dec, search_radius)

if results:
    result = results[0]
    object_id = result["object_id"]
    ra = result["ra"]
    dec = result["dec"]

    data = find_spectrum(object_id)

    if data:
        # data[:10].pprint_all()

        plot_spectrum(data, ra, dec)
    else:
        print("No spectrum found for given object ID")
else:
    print("No object found for given coordinates and search radius")
