Community Analysis Pipeline (CAP)
===============================

Welcome to the Mars Climate Modeling Center (MCMC) **Community Analysis Pipeline (CAP)**.

About
-----
**CAP** is a set of Python3 libraries and command-line executables that streamline downloading, processing, and plotting output from the NASA Ames Mars Global Climate Models:

* NASA Ames Legacy Mars GCM
* NASA Ames Mars GCM with GFDL's FV3 dynamical core

Installation
-----------
Requirements:

* Python 3.7 or later
* pip (Python package installer)

Recommended Installation
^^^^^^^^^^^^^^^^^^^^^^
For reproducible analysis, we recommend installing CAP in a dedicated virtual environment::

    # Create a new virtual environment
    python3 -m venv amescap-env

    # Activate the environment
    source amescap-env/bin/activate     # On Unix/MacOS with bash
    # or
    source amescap-env/bin/activate.csh # On Unix/MacOS with csh/tcsh/zsh
    # or
    amescap-env\Scripts\activate        # On Windows

    # Install CAP and its dependencies
    pip install git+https://github.com/NASA-Planetary-Science/AmesCAP.git

This ensures consistent package versions across different systems.

Available Commands
^^^^^^^^^^^^^^^
After installation, the following commands will be available:

* ``MarsInterp`` - Interpolate data to pressure or altitude coordinates
* ``MarsPull`` - Download model outputs
* ``MarsPlot`` - Create visualizations
* ``MarsVars`` - Process and analyze variables
* ``MarsFiles`` - Manage data files
* ``MarsFormat`` - Convert between data formats
* ``MarsCalendar`` - Handle Mars calendar calculations

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
Data is available through the `MCMC Data Portal <https://data.nas.nasa.gov/mcmc/index.html>`_.

Contributing
-----------
We welcome contributions! Please see our contributing guidelines for details.

License
-------
This project is licensed under the MIT License - see the LICENSE file for details.

Citation
--------
If you use CAP in your research, please cite:
**(APA)** NASA Ames Mars Climate Modeling Center (2024). *Community Analysis Pipeline* [Computer software]. NASA Planetary Science GitHub.