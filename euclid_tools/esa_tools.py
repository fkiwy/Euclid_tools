import os
import tempfile
import numpy as np
from datetime import datetime
from astroquery.esa.euclid import Euclid
from astropy.coordinates import SkyCoord
from astropy.table import Table, QTable, MaskedColumn
from astropy.io import fits
import astropy.units as u

from euclid_tools.shared import add_magnitude


COLUMNS_TO_REMOVE = ["basic_download_data_oid", "to_be_published"]


def retrieve_objects(ra: float, dec: float, radius: float, limit: int = -1) -> Table:
    """
    Perform a cone search on the Euclid Q1 MER catalog and return nearby objects.

    Parameters
    ----------
    ra : float
        Right Ascension in degrees.
    dec : float
        Declination in degrees.
    radius : float
        Search radius in arcseconds.
    limit : int, default -1, meaning ESA defined limit
        The maximum number of rows to return

    Returns
    -------
    Table
        Astropy Table containing catalog entries with added AB and Vega magnitudes.
        Returns `None` if no objects are found.

    Notes
    -----
    - Columns `basic_download_data_oid` and `to_be_published` are removed.
    - Adds VIS, Y, J, and H magnitudes using `add_magnitude`.
    - Results are sorted by angular separation in arcseconds.
    """

    coord = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame="icrs")
    radius = u.Quantity(radius, u.arcsec)

    Euclid.ROW_LIMIT = limit

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


def retrieve_spectrum(object_id: str, mask_bad_values: bool = False) -> QTable:
    """
    Retrieve the 1D spectrum for a given object ID from the Euclid archive.

    Parameters
    ----------
    object_id : str
        Unique identifier of the source.
    mask_bad_values : bool, default False
        Whether to mask bad flux values.

    Returns
    -------
    QTable
        A table containing wavelength (microns), flux, and flux error.
        Flux and error are scaled using the FITS header FSCALE keyword.
        Returns `None` if no spectrum is available for the given object ID.

    Notes
    -----
    - Wavelength is converted to microns.
    - Data is masked based on the bitmask in the `MASK` column if requested.
    """

    filename = _generate_filename(tempfile.gettempdir(), object_id)
    results = Euclid.get_spectrum(retrieval_type="SPECTRA_RGS", source_id=object_id, output_file=filename)

    if results and len(results) > 0:
        try:
            hdu = fits.open(results[0])[1]
        except:
            print("No spectrum available for given object ID")
            return None

        spectrum = QTable.read(hdu, format="fits")
        spec_header = hdu.header

        # Convert wavelength to microns
        wavelength = spectrum["WAVELENGTH"].to(u.micron)

        # Scale the flux and error using the scaling factor from the fits header
        flux = spectrum["SIGNAL"] * spec_header["FSCALE"]
        error = np.sqrt(spectrum["VAR"]) * spec_header["FSCALE"]

        if mask_bad_values:
            # Use the MASK column to create a boolean mask for values to ignore
            bad_mask = (spectrum["MASK"].value % 2 == 1) | (spectrum["MASK"].value >= 64)

            # Apply the mask to the spectrum data
            wavelength = MaskedColumn(wavelength, mask=bad_mask)
            flux = MaskedColumn(flux, mask=bad_mask)
            error = MaskedColumn(error, mask=bad_mask)

            # Replace masked values by NaN
            wavelength = wavelength.filled(np.nan)
            flux = flux.filled(np.nan)
            error = error.filled(np.nan)

        # Create a new QTable for the result
        result = QTable([wavelength, flux, error], names=("WAVELENGTH", "FLUX", "ERROR"))

    return result


def retrieve_cutout(ra: float, dec: float, search_radius: float, cutout_size: float, band: str) -> fits.PrimaryHDU:
    """
    Retrieve an image cutout from Euclid imaging data.

    Parameters
    ----------
    ra : float
        Right Ascension in degrees.
    dec : float
        Declination in degrees.
    search_radius : float
        Radius used to find the overlapping mosaic (in arcseconds).
    cutout_size : float
        Size of the image cutout in arcseconds.
    band : str
        Observing band to extract cutout from ('VIS', 'Y', 'J', 'H').

    Returns
    -------
    fits.HDU
        FITS HDU containing the image cutout.
        Returns `None` if no suitable image is found.

    Raises
    ------
    ValueError
        If the provided band is not one of the allowed options.
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

    coord = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame="icrs")
    radius = (cutout_size / 2) * u.arcsec

    filename = _generate_filename(tempfile.gettempdir(), obs_id)

    print(f"Downloading {band}-band cutout")
    cutout_url = Euclid.get_cutout(
        file_path=file_path,
        instrument=instrument,
        id=obs_id,
        coordinate=coord,
        radius=radius,
        output_file=filename,
    )
    hdul = fits.open(cutout_url[0])

    return hdul[0]


def print_catalog_info():
    """
    Print metadata information for each column in the Euclid MER catalog.

    Notes
    -----
    - Skips columns listed in `COLUMNS_TO_REMOVE`.
    - Renames 'right_ascension' to 'ra' and 'declination' to 'dec' for clarity.
    - Prints column name, unit, and description.
    """

    table = Euclid.load_table("catalogue.mer_catalogue")

    for col in table.columns:
        if col.name in COLUMNS_TO_REMOVE:
            continue
        if col.name == "right_ascension":
            col.name = "ra"
        if col.name == "declination":
            col.name = "dec"

        print(f'{f"{col.name}":45s} {f"{col.unit}":12s} {col.description}')


def _generate_filename(working_dir: str, source_id: str) -> str:
    """
    Generate a unique file path in a temporary directory for storing FITS output.

    Parameters
    ----------
    working_dir : str
        Base directory in which to create the subdirectory.
    source_id : str
        Source identifier used to name the FITS file.

    Returns
    -------
    str
        Full file path to the FITS file.

    Notes
    -----
    - The filename includes a timestamp to ensure uniqueness.
    - Creates the subdirectory if it doesn't already exist.
    - Intended for internal use only.
    """

    # Get the current timestamp and format it
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S%f")

    # Construct the directory path
    directory_path = os.path.join(working_dir, f"temp_{timestamp}")

    # Make sure the directory exists
    os.makedirs(directory_path, exist_ok=True)

    # Create the full file path
    file_path = os.path.join(directory_path, f"{source_id}.fits")

    return file_path
