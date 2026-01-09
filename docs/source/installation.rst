.. _installation:

Installation Instructions
=========================

*Last Updated: May 2025*

Installing CAP is done on the command line via ``git clone``. Here, we provide instructions for installing on Windows using Cygwin or PowerShell and pip or conda, MacOS using pip or conda, and the NASA Advanced Supercomputing (NAS) system using pip.

:ref:`MacOS Installation <mac_install>`

:ref:`Windows Installation <windows_install>`

:ref:`NAS Installation <nas_install>`

:ref:`Spectral Analysis Capabilities <spectral_analysis>`


If you are specifically seeking to use CAP's spatial filtering utilities (low-, high-, and band-pass spatial filtering, or zonal decomposition), please follow the instructions for :ref:`Spectral Analysis Capabilities <spectral_analysis>` below, as these functions require the ``pyshtools`` library for spherical harmonic transforms and other spectral analysis functions.

.. note::

   The AmesCAP package is designed to be installed in a virtual environment. This allows you to manage dependencies and avoid conflicts with other Python packages. We recommend using either `pip` or `conda` for package management.

   If you are using a conda environment, the installation process is straightforward and handles all dependencies automatically. If you are using pip, you will need to install some dependencies manually.

.. _mac_install:

Installation on MacOS
---------------------

This guide provides installation instructions for the AmesCAP package on MacOS using either ``pip`` or ``conda`` for package management.

Prerequisites
^^^^^^^^^^^^^

* A MacOS system with Python 3 installed
* Terminal access
* (Optional) `Anaconda <https://www.anaconda.com/download>`_ or `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ for conda-based installation

Installation Steps
^^^^^^^^^^^^^^^^^^

1. Remove any pre-existing CAP environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you have a **pre-existing** virtual environment holding CAP, we recommend you first remove the virtual environment folder entirely.

**NOTE:** Use the name of your virtual environment. We use ``amescap`` as an example, but you can name it whatever you like.

.. code-block:: bash

   rm -r amescap # For pip virtual environments
   # OR
   conda env remove -n amescap # For conda virtual environments

2. Create and activate a virtual environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Choose either pip or conda to create your virtual environment:

**Using pip:**

.. code-block:: bash

   /Users/username/path/to/preferred/python3 -m venv amescap
   
   source amescap/bin/activate.csh # For CSH/TCSH
   # OR
   source amescap/bin/activate # For BASH

**NOTE:** To list your Python installations, use ``which python3``. Use the path to your preferred Python installation in the pip command above.

**Using conda:**

.. code-block:: bash

   conda create -n amescap python=3.13
   conda activate amescap

3. Install CAP from GitHub
~~~~~~~~~~~~~~~~~~~~~~~~~~

Install CAP from the `NASA Planetary Science GitHub <https://github.com/NASA-Planetary-Science/AmesCAP>`_ using ``pip``:

.. code-block:: bash

   pip install git+https://github.com/NASA-Planetary-Science/AmesCAP.git

.. note::

   You can install a specific branch of the AmesCAP repository by appending ``@branch_name`` to the URL. For example, to install the ``devel`` branch, use:

   .. code-block:: bash

      pip install git+https://github.com/NASA-Planetary-Science/AmesCAP.git@devel
   
   This is useful if you want to test new features or bug fixes that are not yet in the main branch.

4. Copy the profile file to your home directory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   cp amescap/mars_templates/amescap_profile ~/.amescap_profile # For pip
   # OR
   cp /opt/anaconda3/envs/amescap/mars_templates/amescap_profile ~/.amescap_profile # For conda

5. Test your installation
~~~~~~~~~~~~~~~~~~~~~~~~~

While your virtual environment is active, run:

.. code-block:: bash

   MarsPlot -h

This should display the help documentation for MarsPlot.

6. Deactivate the virtual environment when finished
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   deactivate # For pip
   # OR
   conda deactivate # For conda

Troubleshooting Tips
^^^^^^^^^^^^^^^^^^^^

* **Python Version Issues**: Ensure you're using Python 3.6 or newer.
* **Virtual Environment Not Activating**: Verify you're using the correct activation script for your shell.
* **Package Installation Failures**: Check your internet connection and ensure you have permission to install packages.
* **Profile File Not Found**: Double-check the installation paths. The actual path may vary depending on your specific installation.
* **Shell Type**: If you're unsure which shell you're using, run ``echo $SHELL`` to determine your current shell type.

.. _windows_install:

Installation on Windows
-----------------------

This guide provides installation instructions for the AmesCAP package on Windows using either **Windows Terminal (PowerShell)** or **Cygwin**, with either ``pip`` or ``conda`` for package management.

Prerequisites
^^^^^^^^^^^^^

Choose your preferred environment:

Windows Terminal Setup
^^^^^^^^^^^^^^^^^^^^^^
* Install `Python <https://www.python.org/downloads/>`_ for Windows
* Install `Git for Windows <https://git-scm.com/download/win>`_
* Windows Terminal (pre-installed on recent Windows 10/11)
* (Optional) `Anaconda <https://www.anaconda.com/download>`_ or `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ for conda-based installation

Cygwin Setup
^^^^^^^^^^^^
* Install `Cygwin <https://www.cygwin.com/>`_ with these packages:

  * python3
  * python3-pip
  * git
  * bash
* (Optional) `Anaconda <https://www.anaconda.com/download>`_ or `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ for conda-based installation

Installation Steps
^^^^^^^^^^^^^^^^^^

1. Remove any pre-existing CAP environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you have a **pre-existing** virtual environment holding CAP, we recommend you first remove the virtual environment folder entirely.

**NOTE:** Use the name of your virtual environment. We use `amescap` as an example, but you can name it whatever you like.

Using **Windows Terminal (PowerShell):**

.. code-block:: powershell

   Remove-Item -Recurse -Force amescap # For pip virtual environments
   # OR
   conda env remove -n amescap # For conda virtual environments

Using **Cygwin:**

.. code-block:: bash

   rm -r amescap # For pip virtual environments
   # OR
   conda env remove -n amescap # For conda virtual environments

1. Create and activate a virtual environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Using **pip** with **Windows Terminal (PowerShell)**:

.. code-block:: powershell

   # Create virtual environment
   python -m venv amescap

   # Activate the environment
   .\amescap\Scripts\Activate.ps1

**NOTE:** If you get a security error about running scripts, you may need to run:

.. code-block:: powershell

   Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser

Using **pip** with **Cygwin**:

.. code-block:: bash

   # Create virtual environment (use the path to your preferred Python)
   /cygdrive/c/path/to/python3 -m venv amescap
   # Or simply use the Cygwin python:
   python3 -m venv amescap

   # Activate the environment
   source amescap/bin/activate

Using **conda** with **Windows Terminal (PowerShell)**:

.. code-block:: bash

   conda create -n amescap python=3.13
   conda activate amescap

Using **conda** with **Cygwin**:

.. code-block:: bash

   conda create -n amescap python=3.13
   conda activate amescap

1. Install CAP from GitHub
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # The same command works in both PowerShell and Cygwin
   pip install git+https://github.com/NASA-Planetary-Science/AmesCAP.git

.. note::

   You can install a specific branch of the AmesCAP repository by appending ``@branch_name`` to the URL. For example, to install the ``devel`` branch, use:

   .. code-block:: bash

      pip install git+https://github.com/NASA-Planetary-Science/AmesCAP.git@devel
   
   This is useful if you want to test new features or bug fixes that are not yet in the main branch.

4. Copy the profile file to your home directory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Using **Windows Terminal (PowerShell)**:

.. code-block:: powershell

   # For pip installation
   Copy-Item .\amescap\mars_templates\amescap_profile -Destination $HOME\.amescap_profile

   # For conda installation
   Copy-Item $env:USERPROFILE\anaconda3\envs\amescap\mars_templates\amescap_profile -Destination $HOME\.amescap_profile

Using **Cygwin**:

.. code-block:: bash

   # For pip installation
   cp amescap/mars_templates/amescap_profile ~/.amescap_profile

   # For conda installation (adjust path as needed)
   cp /cygdrive/c/Users/YourUsername/anaconda3/envs/amescap/mars_templates/amescap_profile ~/.amescap_profile

5. Test your installation
~~~~~~~~~~~~~~~~~~~~~~~~~

While your virtual environment is active, run:

.. code-block:: bash

   MarsPlot -h

This should display the help documentation for MarsPlot.

6. Deactivate the virtual environment when finished
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Using **pip** with **Windows Terminal (PowerShell)**:

.. code-block:: powershell

   deactivate

Using **pip** with **Cygwin**:

.. code-block:: bash

   deactivate

Using **conda** (both **Windows Terminal** and **Cygwin**):

.. code-block:: bash

   conda deactivate

Troubleshooting Tips
^^^^^^^^^^^^^^^^^^^^

* **Path Issues**: Windows uses backslashes (``\``) for paths, while Cygwin uses forward slashes (``/``). Make sure you're using the correct format for your environment.
* **Permission Errors**: If you encounter permission issues, try running your terminal as Administrator.
* **Virtual Environment Not Activating**: Ensure you're using the correct activation script for your shell.
* **Package Installation Failures**: Check your internet connection and ensure Git is properly installed.
* **Profile File Not Found**: Double-check the installation paths. The actual path may vary depending on your specific installation.
* **HOME**: If you encounter errors related to HOME not defined, set the variable: ``$env:HOME = $HOME`` (PowerShell) or ``export HOME="$USERPROFILE"`` (Cygwin)

.. _nas_install:

Installation in the NASA Advanced Supercomputing (NAS) Environment
------------------------------------------------------------------

This guide provides installation instructions for the AmesCAP package on NASA's Pleiades or Lou supercomputers.

Prerequisites
^^^^^^^^^^^^^

* Access to NASA's Pleiades or Lou supercomputing systems
* Familiarity with Unix command line and modules system
* Terminal access to the NAS environment

Installation Steps
^^^^^^^^^^^^^^^^^^

1. Remove any pre-existing CAP environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you have a **pre-existing** virtual environment holding CAP, we recommend you first remove the virtual environment folder entirely:

.. code-block:: bash

   rm -r amescap

**NOTE:** Use the name of your virtual environment. We use ``amescap`` as an example, but you can name it whatever you like.

2. Create and activate a virtual environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Load the required Python module, then create and activate a virtual environment:

.. code-block:: bash

   module load python3
   python3 -m venv amescap
   
   source amescap/bin/activate.csh # For CSH/TCSH
   # OR
   source amescap/bin/activate # For BASH

3. Install CAP from GitHub
~~~~~~~~~~~~~~~~~~~~~~~~~~

Install CAP from the `NASA Planetary Science GitHub <https://github.com/NASA-Planetary-Science/AmesCAP>`_ using ``pip``:

.. code-block:: bash

   pip install git+https://github.com/NASA-Planetary-Science/AmesCAP.git

.. note::

   You can install a specific branch of the AmesCAP repository by appending ``@branch_name`` to the URL. For example, to install the ``devel`` branch, use:

   .. code-block:: bash

      pip install git+https://github.com/NASA-Planetary-Science/AmesCAP.git@devel
   
   This is useful if you want to test new features or bug fixes that are not yet in the main branch.
   
4. Copy the profile file to your home directory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   cp amescap/mars_templates/amescap_profile ~/.amescap_profile

5. Test your installation
~~~~~~~~~~~~~~~~~~~~~~~~~

While your virtual environment is active, run:

.. code-block:: bash

   MarsPlot -h

This should display the help documentation for MarsPlot.

6. Deactivate the virtual environment when finished
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   deactivate

Troubleshooting Tips
^^^^^^^^^^^^^^^^^^^^

* **Module Conflicts**: If you encounter module conflicts, ensure you run ``module purge`` before loading the Python module.
* **Permission Issues**: Ensure you have the necessary permissions in your directory to create and modify virtual environments.
* **Package Installation Failures**: NAS systems may have restricted internet access. If pip installation fails, contact your system administrator.
* **Profile File Not Found**: Double-check the installation paths. The actual path may vary depending on your specific installation.
* **Python Version**: If you need a different Python version, check available modules with ``module avail python``.
* **Shell Type**: If you're unsure which shell you're using, run ``echo $SHELL`` to determine your current shell type.


.. _spectral_analysis:

Spectral Analysis Capabilities
-----------------------------

CAP includes optional spectral analysis capabilities that require additional dependencies (spatial filtering utilities). These capabilities leverage the ``pyshtools`` library for spherical harmonic transforms and other spectral analysis functions. ``pyshtools`` is a powerful library for working with spherical harmonics and it is an optional dependencies because it can be complex to install. It requires several system-level dependencies, including `libfftw3 <http://www.fftw.org/>`_ and `liblapack <http://www.netlib.org/lapack/>`_ and BLAS libraries, plus Fortran and C compilers. These dependencies are not included in the standard Python installation and may require additional setup.

If you are using a conda environment, these dependencies are automatically installed when you create the environment using the provided ``environment.yml`` file. If you are using pip, you will need to install these dependencies manually. 

Installing with Spectral Analysis Support
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are two recommended ways to install CAP with spectral analysis support:

**Method 1: Using conda (recommended)**

The conda installation method is recommended as it handles all the complex dependencies automatically. You may need to modify these instructions for your specific system, but the following should work on most systems:

.. code-block:: bash

   # Clone the repository
   git clone clone -b devel https://github.com/NASA-Planetary-Science/AmesCAP.git
   cd AmesCAP
   
   # Create conda environment with all dependencies including pyshtools
   conda env create -f environment.yml -n amescap
   
   # Activate the environment
   conda activate amescap

   # Install the package with spectral analysis support
   pip install .[spectral]

   # It is safe to remove the clone after install
   cd .. # Move out of the AmesCAP repository
   rm -rf AmesCAP # Remove the cloned repository

   # Don't forget to copy the profile file to your home directory
   cp /opt/anaconda3/envs/amescap/mars_templates/amescap_profile ~/.amescap_profile

   # To deactivate the environment, run:
   conda deactivate

**Method 1: Using pip**

The pip installation method is less recommended as it requires manual installation of the dependencies. If you choose this method, you will need to install the dependencies separately. The following command will install CAP with spectral analysis support:

.. code-block:: bash

   # Create your virtual environment with pip according to the instructions above. Make sure to follow the instructions for your operating system.

   # Activate your virtual environment

   # Install CAP with spectral analysis support
   pip install "amescap[spectral] @ git+https://github.com/NASA-Planetary-Science/AmesCAP.git"
   # OR for a specific branch, e.g., the devel branch
   pip install "amescap[spectral] @ git+https://github.com/NASA-Planetary-Science/AmesCAP.git@devel"

   # Don't forget to copy the profile file to your home directory
   cp amescap/mars_templates/amescap_profile ~/.amescap_profile
   
   # To deactivate the environment, run:
   deactivate
