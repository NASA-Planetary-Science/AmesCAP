Quick Start Guide
=================

The Community Analysis Pipeline (CAP) is a Python-based command-line tool that performs analysis and creates plots from netCDF files output by the Mars Global Climate Model (MGCM).

Model Compatibility
-------------------

CAP is currently compatible with output from the:

* MCMC Legacy MGCM
* MCMC FV3-based MGCM
* Mars Weather Research and Forecasting Model (MarsWRF)
* OpenMars

Basic Usage
------------

First, install CAP. Detailed :ref:`installation instructions are provided here <installation>`.

Once CAP is installed and you are in your virtual environment, you may type ``cap`` to view the information provided here on the command line.

You can check your CAP version and install date using ``cap version`` or ``cap info``, which returns:

.. code-block:: bash
   CAP Installation Information
   ----------------------------
   Version: 0.3
   Install Date: Fri Mar  7 11:56:48 2025
   Install Location: /Users/path/to/amescap/lib/python3.11/site-packages

Below is a list of the executables in CAP. Use this list to find the executable that performs the operation you desire.

Available Commands
^^^^^^^^^^^^^^^^^^

The CAP command-line tool includes several executables for different operations:

* **MarsCalendar** - Converts L\ :sub:`s` into day-of-year (sol) and vice versa.
* **MarsFiles** - Manipulates entire files (e.g., time-shift, regrid, filter, etc.)
* **MarsFormat** - Transforms non-MGCM model output for compatibility with CAP.
* **MarsInterp** - Interpolates files to pressure or altitude coordinates.
* **MarsPlot** - Generates plots from Custom.in template files.
* **MarsPull** - Queries data from the MGCM repository at data.nas.nasa.gov/mcmc
* **MarsVars** - Performs variable manipulations (e.g., deriving secondary variables, column-integration, etc.)

To see the arguments for each executable, use ``[-h --help]``, for example:

.. code-block:: bash
   MarsVars -h

CAP only works with netCDF files, so use CAP in a directory hosting model output in netCDF format. Then you can pass arguments to CAP executables to perform your desired operation on the file. A good place to start is to use the example command shown below your desired operation.

Example Usage
^^^^^^^^^^^^^

Let's walk through a use case for CAP. Change to a directory you want to save your data in.

1. Begin by using ``MarsPull`` to retrieve MGCM data from the `NAS Data Portal <https://data.nas.nasa.gov/mcmc>`_. 

If you check out the website at the link above, click "Reference Mars Climate Simulations" and then "FV3-based Mars GCM," you'll see a list of files. We will download the one called ``03340.atmos_average_pstd.nc`` and its associated "fixed" file, ``03340.fixed.nc``:

.. code-block:: bash
   MarsPull FV3BETAOUT1 -f 03340.atmos_average_pstd.nc 03340.fixed.nc

.. note::
   The download will take a few minutes. Actual time varies depending on your internet download speed.

While we wait for the download, let's explore how we would know to use this exact command. The :ref:`Available Commands` section above lists the executables and their functions (you can also view this in the terminal by typing ``cap``). This list tells us that we want to use ``MarsPull`` to retrieve data, and that we can use ``[-h --help]`` to view the instructions on how to use MarsPull, like so:

.. code-block:: bash
   MarsPull -h

which outputs:

.. code-block:: bash
   usage: MarsPull [-h] [-list] [-f FILENAME [FILENAME ...]] [-ls LS [LS ...]] [--debug]
                  [{FV3BETAOUT1,ACTIVECLDS,INERTCLDS,NEWBASE_ACTIVECLDS,ACTIVECLDS_NCDF}]

   Uility for downloading NASA Ames Mars Global Climate Model output files from the NAS Data Portal at:https://data.nas.nasa.gov/mcmcref/

   Requires the ``-id`` argument AND EITHER ``-f`` or ``-ls``.

   positional arguments:
   {FV3BETAOUT1,ACTIVECLDS,INERTCLDS,NEWBASE_ACTIVECLDS,ACTIVECLDS_NCDF}
                           Selects the simulation directory from the NAS data portal:
                           https://data.nas.nasa.gov/mcmcref/

                           Current options are:
                           FV3BETAOUT1
                           ACTIVECLDS
                           INERTCLDS
                           NEWBASE_ACTIVECLDS
                           ACTIVECLDS_NCDF
                           MUST be used with either ``-f`` or ``-ls``.
                           Example:
                           > MarsPull ACTIVECLDS -f fort.11_0730
                           OR
                           > MarsPull ACTIVECLDS -ls 90



   options:
   -h, --help            show this help message and exit
   -list, --list_files   Return a list of all the files available for download from:
                           https://data.nas.nasa.gov/mcmcref/

                           Example:
                           > MarsPull -list

   -f FILENAME [FILENAME ...], --filename FILENAME [FILENAME ...]
                           The name(s) of the file(s) to download.
                           Example:
                           > MarsPull ACTIVECLDS -f fort.11_0730 fort.11_0731

   -ls LS [LS ...], --ls LS [LS ...]
                           Selects the file(s) to download based on a range of solar longitudes (Ls).
                           This only works on data in the ACTIVECLDS and INERTCLDS folders.
                           Example:
                           > MarsPull ACTIVECLDS -ls 90
                           > MarsPull ACTIVECLDS -ls 180 360

   --debug               Use with any other argument to pass all Python errors and
                           status messages to the screen when running CAP.
                           Example:
                           > MarsPull ACTIVECLDS -ls 90 --debug


As we can see, MarsPull wants us to provide the simulation directory name and either one or multiple file names or an L\ :sub:`s` range. The directory name isn't very obvious, but it is listed at the end of the URL on the webpage we looked at earlier: `https://data.nas.nasa.gov/mcmcref/fv3betaout1/ <https://data.nas.nasa.gov/mcmcref/fv3betaout1/>`_.

Then, we used the ``[-f --filename]`` argument to specify which files from that page we wanted to download.

2. Inspect the file contents

Once our files are downloaded, we can look at the variables they contain using the "inspect" function in ``MarsPlot``. This is one function you'll want to remember because you'll find its always useful.

.. code-block:: bash
   MarsPlot -i 003340.atmos_average_pstd.nc

The following should be printed to your terminal:

.. image:: ./images/cli_MarsPlot_inspect.png
   :alt: Output from ``MarsPlot -i``

We can see dozens of variables in the file including surface pressure (``ps``) and atmospheric temperature (``temp``). We can use these variables to derive the CO\ :sub:`2` condensation temperature (``Tco2``). Let's derive that variable and add it to the file.

3. Derive and add ``Tco2`` to the file

Calling ``[-h --help]`` on MarsVars will return a list of variables that MarsVars can derive for you. Make sure your netCDF file has the variables required to derive your requested variable first. To add ``Tco2`` to our file, we type:

.. code-block:: bash
   MarsVars 003340.atmos_average_pstd.nc -add Tco2

When that completes, we can inspect the file to confirm that Tco2 was added:

.. code-block:: bash
   MarsPlot -i 003340.atmos_average_pstd.nc

You should see a new variable listed at the bottom of the printed output:

.. code-block:: bash
   Tco2           : ('time', 'pstd', 'lat', 'lon')= (133, 36, 90, 180), CO2 condensation temperature (derived w/CAP)  [K]

Next, let's plot the lowest-level CO\ :sub:`2` condensation temperature and the lowest level atmospheric temperature to compare them.

4. Generate some plots

CAP's plotting executable is MarsPlot, which accepts a template file called ``Custom.in`` from which it generates plots. First we need to make this template file, so we type:

.. code-block:: bash
   MarsPlot -template

This creates ``Custom.in`` in your current directory. Open ``Custom.in`` in your preferred text editor. You can set the syntax highlighting scheme to detect Python in order to make the file more readable.

The template file contains templates for several plot types. Scroll down until you see the first two templates, which are set to ``True`` by default. The default settings would create a topographical map from the ``zsurf`` variable in ``03340.fixed.nc`` and a latitude-level cross-section of the zonal wind (``ucomp``) from a ``03340.atmos_average.nc`` file, the latter of which we do not have. These are just examples that show what you might want to visualize with each plot type.

Set the second template to ``False``: ``Plot 2D lat X lev = False``

We will use the first template to generate two latitude-longitude maps showing the lowest-level temperature and CO\ :sub:`2` condensation temperature. Start by changing the ``Main Variable`` setting to call ``temp`` from the ``03340.atmos_average_pstd.nc`` file. Then, copy and paste that template right below itself, so you have two identical templates. In the second template, change ``temp`` to ``Tco2``. Your ``Custom.in`` file should look like:



Additional Information
----------------------

CAP is developed and maintained by the **Mars Climate Modeling Center (MCMC) at NASA's Ames Research Center** in Mountain View, CA. For more information, visit the `MCMC website <https://www.nasa.gov/space-science-and-astrobiology-at-ames/division-overview/planetary-systems-branch-overview-stt/mars-climate-modeling-center/>`_.
