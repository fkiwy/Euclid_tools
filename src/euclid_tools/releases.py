"""Verified ESA release mappings; endpoint environments do not select releases."""

CATALOG_TABLES = {"q1": "q1.mer_catalogue", "dr1": "dr1.mer_catalogue_wide_survey"}


def normalize_release(release):
    if not isinstance(release, str) or release.lower() not in CATALOG_TABLES:
        raise ValueError("release must be 'q1' or 'dr1'")
    return release.lower()


def catalog_table(release="q1"):
    """Select the operational Q1 or DR1-wide table without changing the API."""
    return CATALOG_TABLES[normalize_release(release)]


def require_product_support(release, product):
    if normalize_release(release) == "dr1":
        raise NotImplementedError(
            f"DR1 {product} retrieval is not verified; no fallback to Q1 products. "
            "See docs/DR1_FOUNDATION.md for official sources and limitations.")
