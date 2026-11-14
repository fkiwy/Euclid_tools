"""Verified ESA release mappings; endpoint environments do not select releases."""

CATALOG_TABLES = {"q1": "q1.mer_catalogue", "dr1": "dr1.mer_catalogue_wide_survey"}
DR1_CATALOGUES = {name: f"dr1.mer_catalogue_{name}" for name in ("wide_survey", "deep_survey")}


def normalize_release(release):
    if not isinstance(release, str) or release.lower() not in CATALOG_TABLES:
        raise ValueError("release must be 'q1' or 'dr1'")
    return release.lower()


def catalog_table(release="q1", catalogue=None):
    """Select Q1 or an explicit DR1 survey; DR1 defaults to wide_survey."""
    if normalize_release(release) == "q1":
        if catalogue is not None:
            raise ValueError("Q1 has one MER catalog; omit catalogue")
        return CATALOG_TABLES["q1"]
    name = "wide_survey" if catalogue is None else catalogue
    if not isinstance(name, str) or name not in DR1_CATALOGUES:
        raise ValueError("ESA catalogue must be wide_survey or deep_survey")
    return DR1_CATALOGUES[name]


def require_product_support(release, product):
    if normalize_release(release) == "dr1":
        raise NotImplementedError(
            f"DR1 {product} retrieval requires an explicit verified product configuration; "
            "the legacy Q1 product route is disabled for DR1. "
            "See docs/ESA_DR1_PRODUCTS.md for the ESA opt-in and limitations.")
