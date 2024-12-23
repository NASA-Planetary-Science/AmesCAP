from setuptools import setup

setup(
    name="amescap",
    version="0.3",
    description="Analysis pipeline for the NASA Ames MGCM",
    url="https://github.com/NASA-Planetary-Science/AmesCAP",
    author="Mars Climate Modeling Center",
    author_email="alexandre.m.kling@nasa.gov",
    license="MIT License",
    scripts=[
        "bin/MarsPull.py",
        "bin/MarsInterp.py",
        "bin/MarsPlot.py",
        "bin/MarsVars.py",
        "bin/MarsFiles.py",
        "bin/MarsFormat.py",
        "bin/MarsCalendar.py"
    ],
    packages=["amescap"],  # Make sure this is present
    entry_points={
        'console_scripts': [
            'cap=amescap.cli:main',  # This will create the 'cap' command
        ],
    },
    install_requires=[
        "requests>=2.31.0",
        "netCDF4>=1.6.5",
        "numpy>=1.26.2",
        "matplotlib>=3.8.2",
        "scipy>=1.11.4",
        "xarray>=2023.5.0",
    ],
    data_files=[
        ("mars_data", ["mars_data/Legacy.fixed.nc"]),
        ("mars_templates", [
            "mars_templates/legacy.in",
            "mars_templates/amescap_profile"
        ])
    ],
    include_package_data=True,
    zip_safe=False,
    python_requires=">=3.7",
)
