import os
import tempfile
import urllib

import astropy.units as u
import numpy as np
import pyvo as vo
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.table import Table, QTable
from astropy.wcs import WCS
from astroquery.ipac.irsa import Irsa

from euclid_tools.backends import resolve_radius
from euclid_tools.shared import add_magnitude, normalize_spectrum_table, convert_flux_to_mag, MagnitudeSystem
from euclid_tools.releases import normalize_release

COLUMNS_TO_REMOVE = ["tileid", "x", "y", "z", "spt_ind", "htm20", "cntr"]

TABLE_MER = "euclid_q1_mer_catalogue"
TABLE_MER_DR1 = "euclid_dr1_mer_catalogue_wide_survey"


DR1_CATALOGUES = {name: f"euclid_dr1_mer_catalogue_{name}" for name in
                  ("wide_survey", "deep_survey", "wide_mode", "deep_mode")}


def catalog_table(release="q1", catalogue=None):
    """Select an explicit DR1 collection; retain Q1's single-table default."""
    if normalize_release(release) == "q1":
        if catalogue is not None:
            raise ValueError("Q1 has one MER catalog; omit catalogue")
        return TABLE_MER
    name = "wide_survey" if catalogue is None else catalogue
    if not isinstance(name, str) or name not in DR1_CATALOGUES:
        raise ValueError("catalogue must be wide_survey, deep_survey, wide_mode or deep_mode")
    return DR1_CATALOGUES[name]


def retrieve_objects(ra: float, dec: float, radius: float = None, search_radius: float = None, *, release: str = "q1", catalogue=None, client=None) -> Table:
    """
    Perform a cone search on the explicitly selected IRSA MER catalog.

    Q1 (default) retains profile/aperture AB and Vega outputs. DR1 selects
    the Wide Survey by default; ``catalogue="deep_survey"`` selects the Deep Survey.
    Raw columns are preserved and derived photometry is profile AB only.
    ``client`` optionally supplies an Astroquery IRSA-compatible client.

    Parameters
    ----------
    ra : float
        Right Ascension in degrees.
    dec : float
        Declination in degrees.
    radius : float
        Search radius in arcseconds.

    Returns
    -------
    Table
        Astropy Table with profile AB magnitudes (and Q1 Vega/aperture magnitudes).
        Returns `None` if no objects are found.

    Notes
    -----
    - Q1 removes `tileid`, `x`, `y`, `z`, `spt_ind`, `htm20`, `cntr`; DR1 retains raw columns.
    - Q1 adds VIS, Y, J, and H AB/Vega magnitudes using `add_magnitude`.
    - Q1 adds aperture AB/Vega magnitudes and errors with prefixes
      `VIS_4fwhm_aper`, `Y_2fwhm_aper`, `J_2fwhm_aper`, and `H_2fwhm_aper`.
    - Results are sorted by angular separation in arcseconds.
    """

    table_name = catalog_table(release, catalogue)
    client = Irsa if client is None else client
    radius = resolve_radius(radius=radius, search_radius=search_radius) * u.arcsec
    position = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame="icrs")

    table = client.query_region(coordinates=position, spatial="Cone", catalog=table_name, radius=radius)

    if len(table) == 0:
        return None

    table = _prepare_catalog_table(table, release=release, catalogue=catalogue)

    # Create a SkyCoord object for the target
    target_coord = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame="icrs")

    # Create a SkyCoord object for the objects in the result table
    object_coords = SkyCoord(table["ra"], table["dec"], unit=(u.deg, u.deg), frame="icrs")

    # Calculate the angular separation between the target and each object in the result table
    separations = target_coord.separation(object_coords).value

    # Convert speparations from degree to arcsec
    sparations = (separations * u.deg).to(u.arcsec)

    # Add the separations to the result table
    table["separation"] = np.round(sparations, 3)

    # Sort the result table on separation
    table.sort("separation")

    return table


def retrieve_object_by_id(object_id: str, *, release: str = "q1", catalogue=None, client=None) -> Table:
    """Retrieve an IRSA Q1 (default) or explicitly selected DR1 catalog object by ID.

    Returns an Astropy Table with the same release-specific catalog columns and derived
    magnitudes as :func:`retrieve_objects`, without a separation column.
    Returns ``None`` if no matching object is found. Invalid IDs raise
    ``ValueError`` before querying the archive.
    """
    table_name = catalog_table(release, catalogue)
    client = Irsa if client is None else client
    object_id = str(object_id)
    digits = object_id[1:] if object_id.startswith("-") else object_id
    if not digits or not digits.isascii() or not digits.isdecimal():
        raise ValueError("object_id must contain decimal digits with an optional leading minus sign")

    query = f"SELECT * FROM {table_name} WHERE object_id = {object_id}"
    table = client.query_tap(query).to_table()
    if len(table) == 0:
        return None
    return _prepare_catalog_table(table, release=release, catalogue=catalogue)


def _prepare_catalog_table(table: Table, *, release: str = "q1", catalogue=None) -> Table:
    """Apply the catalog formatting shared by cone searches and ID lookups."""
    if normalize_release(release) == "dr1":
        # IRSA supplies ra/dec in degrees, and these verified flux pairs in uJy.
        # Do not infer aperture or Vega calibrations from Q1.
        for band, suffix in (("VIS", "psf"), ("Y", "templfit"), ("J", "templfit"), ("H", "templfit")):
            mag, err = convert_flux_to_mag(table[f"flux_{band.lower()}_{suffix}"],
                                           table[f"fluxerr_{band.lower()}_{suffix}"], MagnitudeSystem.AB)
            table[f"{band}_AB_mag"] = mag
            table[f"{band}_AB_err"] = err
            table[f"{band}_AB_mag"].unit = u.mag
            table[f"{band}_AB_err"].unit = u.mag
        table.meta.update(release="DR1", catalogue_table=catalog_table(release, catalogue), archive="IRSA")
        return table
    table.remove_columns(COLUMNS_TO_REMOVE)

    add_magnitude(table, table["flux_vis_psf"], table["fluxerr_vis_psf"], "VIS")
    add_magnitude(table, table["flux_y_templfit"], table["fluxerr_y_templfit"], "Y")
    add_magnitude(table, table["flux_j_templfit"], table["fluxerr_j_templfit"], "J")
    add_magnitude(table, table["flux_h_templfit"], table["fluxerr_h_templfit"], "H")

    for band, aperture in (("VIS", "4fwhm_aper"), ("Y", "2fwhm_aper"), ("J", "2fwhm_aper"), ("H", "2fwhm_aper")):
        flux_column = f"{band.lower()}_{aperture}"
        add_magnitude(table, table[f"flux_{flux_column}"], table[f"fluxerr_{flux_column}"], band, suffix=f"_{aperture}")

    # Replace masked values with NaN
    for col in table.colnames:
        if isinstance(table[col], np.ma.MaskedArray) and np.issubdtype(table[col].dtype, np.floating):  # Only for float-type columns
            table[col] = table[col].filled(np.nan)

    return table


def retrieve_spectrum(object_id: str, save_spectrum: bool = False, output_dir: str = tempfile.gettempdir(), object_name: str = None, *, release: str = "q1", grism=None, product=None,
                      client=None, session=None, timeout=120) -> QTable:
    """
    Retrieve a Q1 (default) or explicitly selected DR1 spectrum from IRSA.

    DR1 uses SpectrumDM with supplied units, returns one unambiguous product,
    and preserves all bins. Select ``grism='red'``/``'blue'`` or an exact
    ``product=(path, hdu)`` if necessary. Use ``retrieve_spectra`` for all products.
    No SSA or Q1 fallback occurs for DR1. See docs/IRSA_DR1_SPECTRA.md.

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
        Returns `None` if no spectrum is available for the given object ID.

    Notes
    -----
    - Wavelength is converted to microns.
    - Data is masked based on the bitmask in the `MASK` column if requested.
    - If no spectrum is found, returns `None`.
    - The output FITS file is named using the object's RA and DEC coordinates, formatted as 'E{ra}{dec}_spectrum.fits'.
    - The output directory defaults to a temporary directory if not specified.
    """

    if normalize_release(release) == "dr1":
        spectra = retrieve_spectra(object_id, grism=grism, product=product,
                                   client=client, session=session, timeout=timeout)
        if not spectra:
            return None
        if len(spectra) != 1:
            raise ValueError(f"Ambiguous DR1 spectrum: {len(spectra)} products; use "
                             "retrieve_spectra or select grism and exact product=(path, hdu)")
        result = spectra[0]
        if save_spectrum:
            filename = (object_name or f"E{object_id}") + "_spectrum.fits"
            result.write(os.path.join(output_dir, filename), format="fits", overwrite=True)
        return result
    if grism is not None or product is not None or client is not None or session is not None:
        raise ValueError("Product selectors and injected clients are supported only for DR1 spectra")

    adql = f"""
    SELECT *
      FROM euclid.objectid_spectrafile_association_q1
     WHERE objectid = {object_id}
       AND path IS NOT NULL
    """

    table = Irsa.query_tap(adql).to_table()

    if len(table) == 0:
        return None

    file_url = urllib.parse.urljoin(Irsa.tap_url, table["path"][0])
    spectrum = QTable.read(file_url, format="votable")

    # Check if the spectrum contains valid data
    if np.isnan(spectrum["WAVELENGTH"]).all():
        return None

    # Extract the data and create an Astropy QTable
    wavelength = spectrum["WAVELENGTH"].to(u.micron)
    flux = spectrum["SIGNAL"]
    error = spectrum["UNCERTAINTY"]
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

    if save_spectrum:
        if not object_name:
            # If no object name is provided, create one based on the object ID
            object_name = f"E{object_id}"

        # Create the full file path for the plot image
        output_filename = os.path.join(output_dir, object_name + "_spectrum.fits")

        # Save the spectrum to a FITS file
        result.write(output_filename, format="fits", overwrite=True)

    return result


def retrieve_cutout(ra: float, dec: float, search_radius: float = None, cutout_size: float = None, band: str = None, radius: float = None, *, release="q1", tiling="wide", selection=None,
                    tile_id=None, product=None, mode="trim", client=None, filesystem=None) -> fits.PrimaryHDU:
    """
    Retrieve a Q1 (default) or explicitly selected IRSA DR1 image cutout.

    DR1 uses SIA and bounded S3 section reads, preserves calibration headers and
    updates WCS. Wide defaults to official survey mosaics; deep defaults to deep
    stacks. Ambiguous results require tile_id or an exact product S3 URI.
    See docs/IRSA_DR1_IMAGES.md for discovery, selections and limits.

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

    if normalize_release(release) == "dr1":
        from euclid_tools.irsa_images import retrieve_cutout as retrieve_dr1_cutout
        return retrieve_dr1_cutout(ra, dec, radius=resolve_radius(radius, search_radius),
                                   cutout_size=cutout_size, band=band, tiling=tiling, selection=selection,
                                   tile_id=tile_id, product=product, mode=mode, client=client, filesystem=filesystem)
    if tiling != "wide" or selection is not None or tile_id is not None or product is not None or client is not None or filesystem is not None or mode != "trim":
        raise ValueError("DR1 mosaic selectors require release='dr1'")

    if band not in ["VIS", "Y", "J", "H"]:
        raise ValueError("Invalid band. Choose from 'VIS', 'Y', 'J', or 'H'.")
    if cutout_size is None:
        raise TypeError("'cutout_size' is required")

    position = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame="icrs")
    search_radius = resolve_radius(radius=radius, search_radius=search_radius) * u.arcsec

    table = Irsa.query_sia(pos=(position, search_radius), collection="euclid_DpdMerBksMosaic")

    results = table[(table["dataproduct_subtype"] == "science")]

    if not results or len(results) == 0:
        return None

    image_url = results[results["energy_bandpassname"] == band]["access_url"][0]

    print(f"Downloading {band}-band cutout")
    with fits.open(image_url, use_fsspec=True) as hdul:
        wcs = WCS(hdul[0].header)
        cutout = Cutout2D(hdul[0].section, position=position, size=cutout_size * u.arcsec, wcs=wcs)
        hdu = fits.PrimaryHDU(data=cutout.data, header=cutout.wcs.to_header())

    return hdu


def print_catalog_info(*, release: str = "q1", catalogue=None):
    """
    Print metadata information for each column in the Euclid MER catalog.

    Notes
    -----
    - Skips columns listed in `COLUMNS_TO_REMOVE`.
    - Prints column name, unit, and description.
    """

    table_name = catalog_table(release, catalogue)
    service = vo.dal.TAPService("https://irsa.ipac.caltech.edu/TAP")
    table = service.tables[table_name]

    for col in table.columns:
        if col.name in COLUMNS_TO_REMOVE:
            continue

        print(f'{f"{col.name}":45s} {f"{col.unit}":12s} {col.description}')


def retrieve_spectra(object_id, *, grism=None, product=None, client=None, session=None, timeout=120):
    """Return every matching DR1 spectrum; Q1's singular API remains unchanged."""
    from euclid_tools.irsa_spectra import retrieve_spectra as retrieve_dr1
    return retrieve_dr1(object_id, grism=grism, product=product, client=client,
                        session=session, timeout=timeout)


def discover_mosaics(ra, dec, radius=None, search_radius=None, **kwargs):
    """Discover verified DR1 wide/deep science mosaics; list all candidates."""
    from euclid_tools.irsa_images import discover_mosaics as discover_dr1
    radius = 2 if radius is None and search_radius is None else resolve_radius(radius, search_radius)
    return discover_dr1(ra, dec, radius=radius, **kwargs)
