import os
import inspect
import tempfile
from datetime import datetime

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table, QTable
from astroquery.esa.euclid import Euclid

from euclid_tools.validation import object_id as validate_object_id, row_limit, positive_scalar
from euclid_tools.photometry import add_dr1_magnitudes, validate_vega_calibration
from euclid_tools.backends import resolve_radius
from euclid_tools.shared import add_magnitude, normalize_spectrum_table, convert_flux_to_mag, MagnitudeSystem
from euclid_tools.releases import catalog_table, normalize_release, require_product_support

def _spectrum_id_keyword(client):
    """Support Astroquery's old source_id and current ids spelling explicitly."""
    return "ids" if "ids" in inspect.signature(client.get_spectrum).parameters else "source_id"


COLUMNS_TO_REMOVE = ["basic_download_data_oid", "to_be_published"]


def retrieve_objects(ra: float, dec: float, radius: float = None, search_radius: float = None, limit: int = -1, *, release: str = "q1", catalogue=None, vega_calibration=None, client=None) -> Table:
    """
    Perform a cone search on the explicitly selected ESA MER catalog.

    ``release="q1"`` (default) preserves legacy outputs. ``release="dr1"``
    defaults to wide-survey profile photometry with AB magnitudes by default.
    ``catalogue="deep_survey"`` explicitly selects the DR1 deep-survey stack catalog.
    ``vega_calibration`` optionally supplies explicit DR1 reference fluxes and sources;
    see docs/DR1_VEGA_CALIBRATION.md. No official DR1 Vega defaults are inferred.
    ``client`` optionally supplies an authenticated Astroquery-compatible client.

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
        Catalog entries with AB magnitudes; Vega is default for Q1 and explicit for DR1.
        Returns `None` if no objects are found.

    Notes
    -----
    - Columns `basic_download_data_oid` and `to_be_published` are removed.
    - Q1 adds VIS, Y, J, and H AB/Vega magnitudes using `add_magnitude`.
    - Q1 adds aperture AB/Vega magnitudes and errors with prefixes
      `VIS_4fwhm_aper`, `Y_2fwhm_aper`, `J_2fwhm_aper`, and `H_2fwhm_aper`.
    - Results are sorted by angular separation in arcseconds.
    """

    table_name = catalog_table(release, catalogue)
    vega_calibration = validate_vega_calibration(vega_calibration, release=release)
    client = Euclid if client is None else client
    radius = u.Quantity(resolve_radius(radius=radius, search_radius=search_radius), u.arcsec)
    coord = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame="icrs")
    if not np.isfinite(coord.ra.deg) or not np.isfinite(coord.dec.deg):
        raise ValueError("Coordinates must be finite degrees")

    limit = row_limit(limit)
    if limit is not None:
        client.ROW_LIMIT = limit

    job = client.cone_search(
        coordinate=coord,
        radius=radius,
        table_name=table_name,
        ra_column_name="right_ascension",
        dec_column_name="declination",
        columns="*",
        async_job=False,
        # ESA deep-table VOTable metadata has an invalid arraysize on
        # calblock_id_list. CSV preserves all fields and exact integer IDs.
        **({"output_format": "csv"} if catalogue == "deep_survey" else {}),
    )
    if job is None:
        raise RuntimeError("ESA catalog query failed; check archive access and client diagnostics")
    table = job.get_results()

    if table is None:
        raise RuntimeError("ESA catalog query returned no result table")
    if len(table) == 0:
        return None

    table = _prepare_catalog_table(table, release=release, catalogue=catalogue, vega_calibration=vega_calibration)

    table.rename_column("dist", "separation")

    separation = (table["separation"] * u.deg).to(u.arcsec)
    table["separation"] = np.round(separation, 3)

    table.sort("separation")

    return table


def retrieve_object_by_id(object_id: str, *, release: str = "q1", catalogue=None, vega_calibration=None, client=None) -> Table:
    """Retrieve an object by numeric ID from ``release="q1"`` or ``"dr1"``.

    DR1 accepts ``catalogue="wide_survey"`` (default) or ``"deep_survey"``.
    It returns verified profile AB photometry by default; Q1 retains AB/Vega and
    aperture outputs. An optional ``client`` supplies the archive connection.

    Returns an Astropy Table with the same release-specific catalog columns and derived
    magnitudes as :func:`retrieve_objects`, without a separation column.
    Returns ``None`` if no matching object is found. Invalid IDs raise
    ``ValueError`` before querying the archive.
    """
    table_name = catalog_table(release, catalogue)
    vega_calibration = validate_vega_calibration(vega_calibration, release=release)
    client = Euclid if client is None else client
    object_id = validate_object_id(object_id)

    query = f"SELECT * FROM {table_name} WHERE object_id = {object_id}"
    job = client.launch_job(query, **({"output_format": "csv"} if catalogue == "deep_survey" else {}))
    if job is None:
        raise RuntimeError("ESA catalog query failed; check archive access and client diagnostics")
    table = job.get_results()
    if table is None:
        raise RuntimeError("ESA catalog query returned no result table")
    if len(table) == 0:
        return None
    return _prepare_catalog_table(table, release=release, catalogue=catalogue, vega_calibration=vega_calibration)


def _prepare_catalog_table(table: Table, *, release: str = "q1", catalogue=None, vega_calibration=None) -> Table:
    """Apply the catalog formatting shared by cone searches and ID lookups."""
    release = normalize_release(release)
    table.remove_columns([name for name in COLUMNS_TO_REMOVE if name in table.colnames])
    if release == "dr1":
        # Verified profile photometry; Vega requires an explicit supplied convention.
        add_dr1_magnitudes(table, vega_calibration=vega_calibration)
        table.rename_column("right_ascension", "ra")
        table.rename_column("declination", "dec")
        table.meta["release"] = "DR1"
        table.meta["catalogue_table"] = catalog_table(release, catalogue)
        table.meta["catalogue"] = "wide_survey" if catalogue is None else catalogue
        return table

    add_magnitude(table, table["flux_vis_psf"], table["fluxerr_vis_psf"], "VIS")
    add_magnitude(table, table["flux_y_templfit"], table["fluxerr_y_templfit"], "Y")
    add_magnitude(table, table["flux_j_templfit"], table["fluxerr_j_templfit"], "J")
    add_magnitude(table, table["flux_h_templfit"], table["fluxerr_h_templfit"], "H")

    for band, aperture in (("VIS", "4fwhm_aper"), ("Y", "2fwhm_aper"), ("J", "2fwhm_aper"), ("H", "2fwhm_aper")):
        flux_column = f"{band.lower()}_{aperture}"
        add_magnitude(table, table[f"flux_{flux_column}"], table[f"fluxerr_{flux_column}"], band, suffix=f"_{aperture}")

    table.rename_column("right_ascension", "ra")
    table.rename_column("declination", "dec")

    return table


def retrieve_spectrum(object_id: str, save_spectrum: bool = False, output_dir: str = tempfile.gettempdir(), object_name: str = None, *, release: str = "q1", grism=None, product=None, client=None, session=None, timeout=120) -> QTable:
    """
    Retrieve the 1D spectrum for a given object ID from the Euclid archive.

    Parameters
    ----------
    object_id : str
        Unique identifier of the source.
    save_spectrum : bool, optional
        If True, saves the spectrum to a FITS file in the specified output directory.
    output_dir : str, optional
        Directory where the spectrum FITS file will be saved if `save_spectrum` is True.
    object_name : str, optional
        Name of the object to use for the output file name. If not provided, defaults to 'E{object_id}'.

    Returns
    -------
    QTable
        A table containing wavelength (microns), flux, flux error, mask, quality flags, and number of dither positions.
        Q1 ESA uses FSCALE; DR1 preserves the scale already embedded in column units.
        Returns `None` if no spectrum is available for the given object ID.

    Notes
    -----
    - Wavelength is converted to microns.
    - Retrieval retains bins and flags; apply a scientifically justified mask separately.
    - If no spectrum is found, returns `None`.
    - Saved output is '<object_name>_spectrum.fits', defaulting to 'E<object_id>_spectrum.fits'.
    - The output directory defaults to a temporary directory if not specified.
    """

    object_id = validate_object_id(object_id)
    positive_scalar(timeout, "timeout")
    if normalize_release(release) == "dr1":
        from euclid_tools.esa_products import retrieve_spectrum as retrieve_dr1_spectrum
        result = retrieve_dr1_spectrum(object_id, grism=grism, product=product,
                                      client=client, session=session, timeout=timeout)
        if result is not None and save_spectrum:
            result.write(os.path.join(output_dir, (object_name or f"E{object_id}") + "_spectrum.fits"),
                         format="fits", overwrite=True)
        return result
    if grism is not None or product is not None or client is not None or session is not None:
        raise ValueError("ESA product selectors/injected clients require release='dr1'")
    require_product_support(release, "spectrum")
    filename = _generate_filename(tempfile.gettempdir(), object_id)
    results = Euclid.get_spectrum(retrieval_type="SPECTRA_RGS", schema="q1",
                                  output_file=filename, **{_spectrum_id_keyword(Euclid): object_id})

    result = None
    if results and len(results) > 0:
        results = list(dict.fromkeys(str(path) for path in results))
        if len(results) != 1:
            raise ValueError("Ambiguous Q1 spectrum: multiple returned products")
        with fits.open(results[0], memmap=False) as hdus:
            if len(hdus) != 2 or not isinstance(hdus[1], fits.BinTableHDU):
                raise ValueError("Expected one extracted Q1 spectrum table")
            hdu = hdus[1].copy()
            spectrum = QTable.read(hdu, format="fits")

        # Check if the spectrum contains valid data
        if np.isnan(spectrum["WAVELENGTH"]).all():
            return None

        # Extract the data and create an Astropy QTable
        wavelength = spectrum["WAVELENGTH"].to(u.micron)

        fscale = hdu.header.get("FSCALE", 1.0)

        flux = spectrum["SIGNAL"] * fscale
        # Negative variance is invalid; it must not become a zero uncertainty.
        with np.errstate(invalid="ignore"):
            error = np.sqrt(spectrum["VAR"]) * fscale

        mask = spectrum["MASK"]

        columns = [
            wavelength,
            flux,
            error,
            mask,
        ]
        names = [
            "WAVELENGTH",
            "FLUX",
            "ERROR",
            "MASK",
        ]

        if "QUALITY" in spectrum.colnames:
            columns.append(spectrum["QUALITY"])
            names.append("QUALITY")

        if "NDITH" in spectrum.colnames:
            columns.append(spectrum["NDITH"])
            names.append("NDITH")

        result = QTable(columns, names=names)
        result = normalize_spectrum_table(result)
        from euclid_tools.provenance import stamp
        stamp(result.meta, product_type="spectrum", release="Q1", archive="ESA",
              route="esa-get_spectrum", OBJECTID=str(object_id), SRCHDU=1,
              SRCPATH=hdu.header.get("FILENAME"),
              GRISM=hdu.header.get("GRISM"))

        if save_spectrum:
            if not object_name:
                # If no object name is provided, create one based on the object ID
                object_name = f"E{object_id}"

            # Create the full file path for the plot image
            output_filename = os.path.join(output_dir, object_name + "_spectrum.fits")

            # Save the spectrum to a FITS file
            result.write(output_filename, format="fits", overwrite=True)

    return result


def retrieve_cutout(ra: float, dec: float, search_radius: float = None, cutout_size: float = None, band: str = None, radius: float = None, *, release: str = "q1", tiling="wide", selection=None, tile_id=None, product=None, client=None, session=None, timeout=120) -> fits.PrimaryHDU:
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

    Notes
    -----
    Q1 is the default. DR1 requires release="dr1" and explicit wide/deep
    tiling/class choices when needed; ambiguity requires tile_id or an exact
    product. See docs/ESA_DR1_PRODUCTS.md. Native pixels, WCS and calibration
    metadata are preserved within the supported formats.

    Raises
    ------
    ValueError
        If the provided band is not one of the allowed options.
    """

    if cutout_size is None:
        raise TypeError("'cutout_size' is required")
    positive_scalar(cutout_size, "cutout_size")
    if normalize_release(release) == "dr1":
        from euclid_tools.esa_products import retrieve_cutout as retrieve_dr1_cutout
        return retrieve_dr1_cutout(ra, dec, radius=resolve_radius(radius, search_radius),
                                  cutout_size=cutout_size, band=band, tiling=tiling,
                                  selection=selection, tile_id=tile_id, product=product,
                                  client=client, session=session, timeout=timeout)
    if tiling != "wide" or selection is not None or tile_id is not None or product is not None or client is not None or session is not None:
        raise ValueError("ESA mosaic selectors/injected clients require release='dr1'")
    require_product_support(release, "cutout")
    if band not in ["VIS", "Y", "J", "H"]:
        raise ValueError("Invalid band. Choose from 'VIS', 'Y', 'J', or 'H'.")
    if cutout_size is None:
        raise TypeError("'cutout_size' is required")

    search_radius = resolve_radius(radius=radius, search_radius=search_radius)

    if band == "VIS":
        instrument_name = "VIS"
        filter_name = "VIS"
    else:
        instrument_name = "NISP"
        filter_name = "NIR_" + band

    radius_deg = search_radius * u.arcsec.to(u.deg)

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
      FROM q1.mosaic_product
     WHERE (instrument_name='{instrument_name}') 
       AND (filter_name='{filter_name}') 
       AND (((fov IS NOT NULL AND INTERSECTS(CIRCLE('ICRS', {ra}, {dec},{radius_deg}), fov)=1)))
     ORDER BY tile_index ASC
    """
    job_async = Euclid.launch_job(query)
    if job_async is None:
        raise RuntimeError("ESA Q1 mosaic query failed; no cutout requested")
    results = job_async.get_results()
    if results is None:
        raise RuntimeError("ESA Q1 mosaic query returned no result table")

    if len(results) == 0:
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
    if not cutout_url:
        return None
    hdul = fits.open(cutout_url[0])

    from euclid_tools.provenance import stamp
    stamp(hdul[0].header, product_type="image", release="Q1", archive="ESA",
          route="tap-mosaic+esa-get_cutout", SRCPATH=str(file_path), PRODBAND=band,
          TILEID=str(obs_id))
    return hdul[0]


def discover_spectra(object_id, *, release="dr1", **kwargs):
    """List ESA DR1 parent/grism associations without spectrum downloads."""
    if normalize_release(release) != 'dr1':
        raise ValueError("ESA spectrum discovery requires release='dr1'")
    from euclid_tools.esa_products import discover_spectra as discover_dr1
    return discover_dr1(object_id, **kwargs)


def retrieve_spectra(object_id, *, release="dr1", **kwargs):
    """Return independently identifiable DR1 products; Q1 retains its singular API."""
    if normalize_release(release) != 'dr1':
        raise ValueError("Plural ESA spectrum retrieval requires release='dr1'")
    from euclid_tools.esa_products import retrieve_spectra as retrieve_dr1
    return retrieve_dr1(object_id, **kwargs)


def discover_mosaics(ra, dec, radius=None, search_radius=None, *, release="dr1", **kwargs):
    """Discover ESA DR1 survey/deep mosaics, or explicitly select repeat visits."""
    if normalize_release(release) != 'dr1':
        raise ValueError("ESA mosaic discovery requires release='dr1'")
    from euclid_tools.esa_products import discover_mosaics as discover_dr1
    radius = 2 if radius is None and search_radius is None else resolve_radius(radius, search_radius)
    return discover_dr1(ra, dec, radius=radius, **kwargs)


def print_catalog_info(*, release: str = "q1", catalogue=None):
    """
    Print metadata information for each column in the Euclid MER catalog.

    Notes
    -----
    - Skips columns listed in `COLUMNS_TO_REMOVE`.
    - Renames 'right_ascension' to 'ra' and 'declination' to 'dec' for clarity.
    - Prints column name, unit, and description.
    """

    table = Euclid.load_table(catalog_table(release, catalogue))

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
