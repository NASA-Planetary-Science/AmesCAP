Community Analysis Pipeline (CAP)
===============================

Welcome to the Mars Climate Modeling Center (MCMC) **Community Analysis Pipeline (CAP)**.

For instructions and documentation please see our `online documentation <https://amescap.readthedocs.io>`_.

About
-----
**CAP** is a set of Python3 libraries and command-line executables that streamline downloading, processing, and plotting output from the NASA Ames Mars Global Climate Models:

* NASA Ames Legacy Mars GCM
* NASA Ames Mars GCM with GFDL's FV3 dynamical core

Installation
-----------
Requirements:

* Python 3.9 through 3.11
* pip (Python package installer)

Recommended Installation
^^^^^^^^^^^^^^^^^^^^^^
For reproducible analysis, we recommend installing CAP in a dedicated virtual environment. Please reference our `installation instructions <https://amescap.readthedocs.io/en/latest/installation.html>`_ online. Briefly, a virtual environment looks like this::

    # Create a new virtual environment with pip or conda:
    python3 -m venv amescap-env # with pip
    # OR
    conda create -n amescap python=3.11 # with conda

    # Activate the environment, which varies by OS, shell, and package manager:
    source amescap-env/bin/activate     # pip + Unix/MacOS (bash) OR Windows Cygwin
    # OR
    source amescap-env/bin/activate.csh # pip + Unix/MacOS (csh/tcsh/zsh)
    # OR
    amescap-env\Scripts\activate        # pip + Windows PowerShell
    # OR
    conda activate amescap-env          # conda + Unix/MacOS (bash, csh, tcsh, zsh) OR Windows Cygwin
    # OR
    .\amescap\Scripts\Activate.ps1      # conda + Windows PowerShell

    # Install CAP and its dependencies
    pip install git+https://github.com/NASA-Planetary-Science/AmesCAP.git
    # OR install a specific branch with:
    pip install git+https://github.com/NASA-Planetary-Science/AmesCAP.git@devel

    # For spectral analysis capabilities, please follow the installation instructions.

    # Unix, macOS, or Cygwin, with pip or conda:
    cp "$(python -c 'import sys; print(sys.prefix)')/mars_templates/amescap_profile" ~/.amescap_profile
    # Windows PowerShell, with pip or conda:
    $profileDir = python -c "import sys; print(sys.prefix)"
    Copy-Item "$profileDir\mars_templates\amescap_profile" -Destination "$HOME\.amescap_profile"

This ensures consistent package versions across different systems.

For spectral analysis capabilities, please follow the `installation instructions <https://amescap.readthedocs.io/en/latest/installation.html>`_.

Available Commands
^^^^^^^^^^^^^^^
After installation, the following commands will be available:

* ``MarsInterp`` - Interpolate data to pressure or altitude coordinates
* ``MarsPull`` - Download model outputs
* ``MarsPlot`` - Create visualizations
* ``MarsVars`` - Process and analyze variables
* ``MarsFiles`` - Manage data files
* ``MarsFormat`` - Convert between model/reanalysis formats
* ``MarsCalendar`` - Handle Mars calendar calculations
* ``MarsNest`` - Map the layout of nested grids

Documentation
------------
Full documentation is available at `readthedocs.io <https://amescap.readthedocs.io>`_.

Getting Started
^^^^^^^^^^^^^
The tutorial directory contains:

* Installation instructions for Linux, MacOS, and Windows
* Documentation of CAP functions
* Practice exercises to familiarize users with CAP

  * NASA Ames MGCM Tutorial
  * Legacy GCM Tutorial

Data Sources
-----------
The tutorials use MGCM simulation outputs documented in `Haberle et al. 2019 <https://www.sciencedirect.com/science/article/pii/S0019103518305761>`_. 
Legacy files are available at `the NASA Ames Legacy GCM data directory <https://data.nas.nasa.gov/legacygcm/legacygcm/>`_, and FV3 files at `the FV3BETAOUT1 directory <https://data.nas.nasa.gov/legacygcm/fv3betaout1/fv3betaout1/>`_.

Contributing
-----------
We welcome contributions! Please see our contributing guidelines (``CONTRIBUTING.rst``) for details. Instructions for running the test suite, including its time, memory, disk, and network requirements, are in ``tests/README.rst``.

License
-------
This project is licensed under the MIT License - see the LICENSE file for details.

Citation
--------
If you use CAP in your research, please cite:
**(APA)** NASA Ames Mars Climate Modeling Center (2024). *Community Analysis Pipeline* [Computer software]. NASA Planetary Science GitHub.
