Running the CAP Tests
=====================

Install CAP in a Python 3.10 or newer environment from the repository root, then
run the tests from the ``tests`` directory:

.. code-block:: bash

    pip install -e .
    cd tests
    python -m unittest -v test_*.py        # whole suite
    python -m unittest -v test_marsplot.py # one module

Each ``test_*.py`` module can also be run on its own, which is how the GitHub
Actions workflows in ``.github/workflows`` run them.

Test Data
---------

The integration tests do not download data. Each module generates synthetic
netCDF fixtures in a temporary directory, which is removed afterwards:

* ``create_ames_gcm_files.py`` writes NASA Ames MGCM-like files
  (``01336.atmos_average.nc``, ``01336.atmos_diurn.nc``, ``01336.fixed.nc``,
  etc.)
* ``create_gcm_files.py`` writes EMARS, OpenMARS, PCM, and MarsWRF-like files
  for ``test_marsformat.py``

Both scripts take a ``short`` argument that keeps every variable and dimension
but shortens the time axis. Compact fixtures are the default for the test
suite; set ``AMESCAP_FULL_FIXTURES=1`` to repeat the tests with full-size
files.

Resource Requirements
---------------------

Measured with compact fixtures on an Apple-silicon laptop (Python 3.11).
Memory is the peak resident memory of the largest process; disk is the space
used by the fixtures of a module while it runs.

====================================  =========  ===========  ============
Module                                Time       Peak memory  Disk
====================================  =========  ===========  ============
``test_fv3_utils.py``                 < 1 s      0.1 GB       none
``test_ncdf_wrapper.py``              < 1 s      0.1 GB       none
``test_marsplot_utils.py``            < 1 s      0.1 GB       none
``test_marsnest.py``                  2 s        0.1 GB       < 1 MB
``test_marspull.py``                  2 s        0.1 GB       none
``test_marscalendar.py``              6 s        0.3 GB       0.2 GB
``test_marsfiles.py``                 12 s       0.9 GB       0.2 GB
``test_marsinterp.py``                11 s       1.7 GB       0.2 GB
``test_marsvars.py``                  41 s       0.3 GB       0.2 GB
``test_marsplot.py``                  13 s       0.3 GB       0.2 GB
``test_marsformat.py``                62 s       2.4 GB       1.3 GB
**Whole suite**                       2.5 min    2.4 GB       1.3 GB
====================================  =========  ===========  ============

With ``AMESCAP_FULL_FIXTURES=1`` the fixtures need about 5.1 GB (MGCM files,
per module) and 9.7 GB (``test_marsformat.py``) of disk, and generating the
full MGCM files alone peaks at about 8.6 GB of memory.

Network Access
--------------

The default suite runs offline; ``test_marspull.py`` replaces the NAS Data
Portal with mocked responses. Tests that contact the portal are opt-in:

.. code-block:: bash

    # Directory listings and the ~350 KB FV3BETAOUT1 fixed file
    AMESCAP_LIVE_TESTS=1 python -m unittest -v test_marspull.py

    # Also download two legacy fort.11 files (~450 MB each)
    AMESCAP_LIVE_TESTS=1 AMESCAP_LARGE_DOWNLOADS=1 python -m unittest -v test_marspull.py

The small live check also runs weekly on GitHub Actions
(``.github/workflows/marspull_live_check.yml``) to detect changes to the
portal.
