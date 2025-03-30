import sys
import os

sys.path.append(os.path.abspath(os.path.join("..", "src")))

from object_finder import find_objects


def print_results(results):
    table = results[
        "object_id",
        "ra",
        "dec",
        "VIS_AB_mag",
        "VIS_AB_err",
        "VIS_Vega_mag",
        "VIS_Vega_err",
        # 'Y_AB_mag', 'Y_AB_err', 'Y_Vega_mag', 'Y_Vega_err',
        # 'J_AB_mag', 'J_AB_err', 'J_Vega_mag', 'J_Vega_err',
        # 'H_AB_mag', 'H_AB_err', 'H_Vega_mag', 'H_Vega_err',
        "separation",
    ]
    table.pprint_all()


ra, dec = 273.891148, 64.526926
search_radius = 5

results = find_objects(ra, dec, search_radius)
print_results(results)

results = find_objects(ra, dec, search_radius, nearest_object=False)
print_results(results)
