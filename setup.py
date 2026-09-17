from setuptools import find_packages, setup

setup(
    name="Euclid_tools",
    version="0.1.0",
    description="Package for Euclid Q1 spectra, catalog queries, and image cutouts.",
    url="https://github.com/fkiwy/Euclid_tools",
    author="Frank Kiwy",
    author_email="frank.kiwy@outlook.com",
    license="MIT",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "aiohttp",
        "astropy",
        "astroquery",
        "matplotlib",
        "numpy",
        "requests",
        "reproject",
        "scipy",
        "Pillow",
    ],
    package_data={
        "flux_comp": [
            "templates/Burgasser+2017/*.fits",
            "templates/Kesseli+2017/*.fits",
            "templates/Theissen+2022/*.fits",
        ]
    },
    zip_safe=False,
)
