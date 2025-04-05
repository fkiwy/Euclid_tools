from setuptools import setup

setup(
    name="Euclid_tools",
    version="0.1.0",
    description="A package of tools related to ESA's Euclid mission",
    url="https://github.com/fkiwy/Euclid_tools",
    author="Frank Kiwy",
    author_email="frank.kiwy@outlook.com",
    license="MIT",
    install_requires=[
        "astropy",
        "astroquery",
        "matplotlib",
        "numpy",
        "requests",
        "scipy",
        "Pillow",
    ],
    packages=["euclid_tools", "flux_comp"],
    package_dir={"euclid_tools": "tools", "flux_comp": "fluxcomp"},
    package_data={
        "flux_comp": [
            "fluxcomp/templates/Burgasser+2017/*.fits",
            "fluxcomp/templates/Kesseli+2017/*.fits",
            "fluxcomp/templates/Theissen+2022/*.fits",
        ]
    },
    zip_safe=False,
    include_package_data=True,
)
