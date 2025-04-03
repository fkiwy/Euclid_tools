import os
import tempfile
import numpy as np
from datetime import datetime
from astroquery.esa.euclid import Euclid
from astropy.coordinates import SkyCoord
from astropy.table import Table, QTable, MaskedColumn
from astropy.io import fits
import astropy.units as u

import tools.shared as shr


COLUMNS_TO_REMOVE = ["basic_download_data_oid", "to_be_published"]


def retrieve_objects(ra: float, dec: float, radius: float) -> Table:
    """
    Perform a cone search on the Euclid archive.

    Parameters:
    ra (float): Right ascension in degrees.
    dec (float): Declination in degrees.
    radius (float): Search radius in arcseconds.

    Returns:
    Table: Astropy table containing catalog entries.
    """
    coord = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg))
    radius = u.Quantity(radius, u.arcsec)
    job = Euclid.cone_search(
        coordinate=coord,
        radius=radius,
        table_name="catalogue.mer_catalogue",
        ra_column_name="right_ascension",
        dec_column_name="declination",
        columns="*",
        async_job=False,
    )
    table = job.get_results()

    if len(table) == 0:
        return None

    table.remove_columns(COLUMNS_TO_REMOVE)

    add_magnitude(table, table["flux_vis_psf"], table["fluxerr_vis_psf"], "VIS")
    add_magnitude(table, table["flux_y_templfit"], table["fluxerr_y_templfit"], "Y")
    add_magnitude(table, table["flux_j_templfit"], table["fluxerr_j_templfit"], "J")
    add_magnitude(table, table["flux_h_templfit"], table["fluxerr_h_templfit"], "H")

    table.rename_column("right_ascension", "ra")
    table.rename_column("declination", "dec")
    table.rename_column("dist", "separation")

    separation = (table["separation"] * u.deg).to(u.arcsec)
    table["separation"] = np.round(separation, 3)

    table.sort("separation")

    return table


def retrieve_spectrum(object_id: str, ignore_bad_values: bool = False) -> fits.HDUList:
    """
    Retrieve the spectrum for a given object ID.

    Parameters:
    object_id (str): The ID of the object to retrieve the spectrum for.

    Returns:
    fits.HDUList: HDU list containing the spectrum data.
    """
    filename = generate_filename(tempfile.gettempdir(), object_id)
    results = Euclid.get_spectrum(retrieval_type="SPECTRA_RGS", source_id=object_id, output_file=filename)

    if results and len(results) > 0:
        hdu = fits.open(results[0])[1]
        spectrum = QTable.read(hdu, format="fits")
        spec_header = hdu.header

        # Convert wavelength to microns
        wavelength = spectrum["WAVELENGTH"].to(u.micron)

        # Scale the flux and error using the scaling factor from the fits header
        flux = spectrum["SIGNAL"] * spec_header["FSCALE"]
        error = np.sqrt(spectrum["VAR"]) * spec_header["FSCALE"]

        if ignore_bad_values:
            # Use the MASK column to create a boolean mask for values to ignore
            bad_mask = (spectrum["MASK"].value % 2 == 1) | (spectrum["MASK"].value >= 64)

            # Apply the mask to the flux and error values
            flux = MaskedColumn(flux, mask=bad_mask)
            error = MaskedColumn(error, mask=bad_mask)

        # Create a new QTable for the result
        result = QTable([wavelength, flux, error], names=("WAVELENGTH", "FLUX", "ERROR"))

    return result


def retrieve_cutout(ra: float, dec: float, search_radius: float, cutout_size: float, band: str) -> fits.HDUList:
    """
    Retrieve a cutout image from the Euclid archive.

    Parameters:
    ra (float): Right ascension in degrees.
    dec (float): Declination in degrees.
    search_radius (float): Search radius in arcseconds.
    cutout_size (float): Size of the cutout in arcseconds.
    band (str): Euclid band ('VIS', 'Y', 'J', 'H').

    Returns:
    fits.HDUList: HDU list containing the image cutout.
    """
    if band not in ["VIS", "Y", "J", "H"]:
        raise ValueError("Invalid band. Choose from 'VIS', 'Y', 'J', or 'H'.")

    if band == "VIS":
        instrument_name = "VIS"
        filter_name = "VIS"
    else:
        instrument_name = "NISP"
        filter_name = "NIR_" + band

    radius = search_radius * u.arcsec.to(u.deg)

    query = f"""
    SELECT file_name, 
           file_path, 
           datalabs_path, 
           mosaic_product_oid, 
           tile_index, 
           instrument_name, 
           filter_name, 
           ra, 
           dec 
      FROM sedm.mosaic_product 
     WHERE (instrument_name='{instrument_name}') 
       AND (filter_name='{filter_name}') 
       AND (((mosaic_product.fov IS NOT NULL AND INTERSECTS(CIRCLE('ICRS', {ra}, {dec},{radius}), mosaic_product.fov)=1))) 
     ORDER BY mosaic_product.tile_index ASC
    """
    job_async = Euclid.launch_job(query)
    results = job_async.get_results()

    if not results or len(results) == 0:
        return None

    result = results[0]
    file_path = result["file_path"] + "/" + result["file_name"]
    instrument = result["instrument_name"]
    obs_id = result["tile_index"]

    coord = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg))
    radius = (cutout_size / 2) * u.arcsec

    filename = generate_filename(tempfile.gettempdir(), obs_id)

    cutout_filepath = Euclid.get_cutout(
        file_path=file_path, instrument=instrument, id=obs_id, coordinate=coord, radius=radius, output_file=filename
    )

    hdul = fits.open(cutout_filepath[0])
    return hdul[0]


def print_catalog_info():
    table = Euclid.load_table("catalogue.mer_catalogue")

    for col in table.columns:
        if col.name in COLUMNS_TO_REMOVE:
            continue
        if col.name == "right_ascension":
            col.name = "ra"
        if col.name == "declination":
            col.name = "dec"

        print(f'{f"{col.name}":45s} {f"{col.unit}":12s} {col.description}')


def add_magnitude(table, flux, flux_err, band):
    table[band + "_AB_mag"], table[band + "_AB_err"] = shr.convert_flux_to_mag(
        flux, flux_err, magnitude_system=shr.MagnitudeSystem.AB
    )
    table[band + "_AB_mag"].unit = u.mag
    table[band + "_AB_err"].unit = u.mag

    table[band + "_Vega_mag"], table[band + "_Vega_err"] = shr.convert_flux_to_mag(
        flux, flux_err, magnitude_system=shr.MagnitudeSystem.Vega, band=band
    )
    table[band + "_Vega_mag"].unit = u.mag
    table[band + "_Vega_err"].unit = u.mag


def generate_filename(working_dir: str, source_id: str) -> str:
    # Get the current timestamp and format it
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    # Construct the directory path
    directory_path = os.path.join(working_dir, f"temp_{timestamp}")

    # Make sure the directory exists
    os.makedirs(directory_path, exist_ok=True)

    # Create the full file path
    file_path = os.path.join(directory_path, f"{source_id}.fits")

    return file_path
