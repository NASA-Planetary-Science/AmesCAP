Installation
============

*Last Updated: July 2023*

Installing CAP is done on the command line via ``git clone``. If you are installing CAP locally, please see the :ref:`Local Installation <local_install>` instructions below. If you are installing CAP on NASA's supercomputing environment, NAS, please see the :ref:`NAS Installation <local_install>` instructions below.

.. _local_install:

Installation on Local Computer
------------------------------

If you have a **pre-existing** virtual environment holding CAP, we recommend you first remove the virtual environment folder entirely:

.. code-block:: bash
    
    rm -r amescap    # change amescap to the name of your virtual environment

1. (Re)create a virtual environment in which to install CAP. Use your preferred installation of python when creating your virtual environment. Here, we create a virtual environment called ``amescap``:

.. code-block:: bash
    
    /Users/username/path/to/preferred/python3 -m venv --system-site-packages amescap

2. Activate the virtual environment using the appropriate syntax for your shell:

.. code-block:: bash

    source amescap/bin/activate         # for bash
    source amescap/bin/activate.csh     # for CSH/TCSH

3. Install ``cftime`` and ``netCDF4``, and upgrade ``pip`` **within** the virtual environment:

.. code-block:: bash
    
    pip install cftime==1.0.3.4
    pip install netCDF4==1.6.1
    python3 -m pip install --upgrade pip

4. Install CAP from the `NASA Planetary Science GitHub <https://github.com/NASA-Planetary-Science/AmesCAP>`_ using ``git clone``:

.. code-block:: bash

    pip install git+https://github.com/NASA-Planetary-Science/amescap.git

5. Test that your installation completed properly by typing within your virtual environment:

.. code-block:: bash

    MarsPlot -h

This should show the "help" documentation for MarsPlot.

6. Deactivate the virtual environment with:

.. code-block:: bash

    deactivate




.. _nas_install:

Installation on NAS
-------------------

These instructions are specific to installing CAP on Pleiades or Lou, NASA's supercomputers.

If you have a **pre-existing** virtual environment for CAP, begin by removing the virtual environment directory and all of its subdirectories:

.. code-block:: bash
    
    rm -r amescap    # or whatever the name is for you

1. (Re)create a virtual environment in which to hold CAP. Use your preferred installation of python when creating your virtual environment. Here, we create a virtual environment called ``amescap``:

.. code-block:: bash
    
    /Users/username/path/to/preferred/python3 -m venv --system-site-packages amescap

2. Load the necessary modules:

.. code-block:: bash

    module purge
    module load python3/3.9.12

3. Activate the virtual environment using the appropriate syntax for your shell:

.. code-block:: bash

    source AmesCAP/bin/activate.csh     # for CSH/TSCH
    source AmesCAP/bin/activate         # for BASH

4. Install ``cmake`` and ``setuptools``, and upgrade ``pip`` **within** the virtual environment:

.. code-block:: bash
    
    pip install cmake
    pip install setuptools --upgrade
    pip install --upgrade pip

4. Install CAP from the `NASA Planetary Science GitHub <https://github.com/NASA-Planetary-Science/AmesCAP>`_ using ``git clone``:

.. code-block:: bash

    pip install git+https://github.com/NASA-Planetary-Science/amescap.git

6. Deactivate the virtual environment with:

.. code-block:: bash

    deactivate
