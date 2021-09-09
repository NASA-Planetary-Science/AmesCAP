Welcome to the Mars Climate Modeling Center (MCMC) Community Analysis Pipeline (CAP). By the end of this tutorial, you will know how to download Mars Climate data from the MCMC data portal, reduce these large climate simulations to meaningful data, and make plots for Martian winds, temperature, and aerosols at specific seasons and locations.

The simulation results presented on this page are extensively documented in [Haberle et al. 2019](https://www.sciencedirect.com/science/article/pii/S0019103518305761)

# INSTALLATION {#INSTALLATION}

The analysis pipeline is entirely written in pure Python which is an intuitive and open source programming language. You may identify yourself in one the following categories:

* A. You are familiar with the Python infrastructure and would like to install the Analysis pipeline on top of your current Python installation: Check the requirements below and skip to [Installing the Pipeline](#Install_CAP). Note that you may have to manually add aliases to the `Mars***.py` executables to your search path.

* B. You have experience with Python but not with managing packages, or you are new to Python: To ensure that there is no conflict with other Python versions that may be on your system, we will install a fresh Python distribution **locally** (this does not require admin permission). Additionally, we will install the analysis pipeline in a self-contained **virtual environment**. A virtual environment is basically a stand-alone copy of your entire Python distribution minus the 'core' code that is shared with the main Python distribution.  This will allow you to use your fresh Python installation for other projects (including installing or upgrading packages) without the risk of altering CAP. It will also be safe to alter (or even completely delete) the virtual environment without breaking the main Python distribution.

##  Requirements {#Requirements}

**Python 3**
If you are already a Python user, you can install CAP on top of your current Python installation. For new users, we recommend you use the latest version of the Anaconda Python distribution because it ships with pre-compiled math and plotting packages (e.g. `numpy`, `matplotlib`) as well as pre-compiled libraries (e.g. `hdf5` headers to read `netCDF` files). You can download Anaconda using either the command-line installer or the graphical interface via the instructions [here](https://www.anaconda.com/distribution/#download-section). If you are the owner ("admin") of your system, you can choose to install Python at the system level. Otherwise, you can install Python in your home directory if you don't have permission to install at the system level.

* In MacOS and Linux, you can install Python3 **locally** from a terminal by typing the following in the command line:

```bash
$ chmod +x Anaconda3-2021.05-MacOSX-x86_64.sh   #this makes the .sh file executable
$ ./Anaconda3-2021.05MacOSX-x86_64.sh           #this runs the executable

> Welcome to Anaconda3 2021.05
>
> In order to continue the installation process, please review the license agreement.
> Please, press ENTER to continue
> >>> 
$ [ENTER]
> ...
> Do you accept the license terms? [yes|no]
> [no] >>>
$ yes
> Anaconda3 will now be installed into this location:
> /Users/username/anaconda3
> 
>  - Press ENTER to confirm the location
>  - Press CTRL-C to abort the installation
>  - Or specify a different location below
>
> [/Users/cbatters/anaconda3] >>> 
$ [ENTER]
> PREFIX=/Users/cbatters/anaconda3
> Unpacking payload ...
> Collecting package metadata (current_repodata.json): done                                                       
> Solving environment: done
> 
> ## Package Plan ##
> ...
> Preparing transaction: done
> Executing transaction: - 
> done
> installation finished.
> Do you wish the installer to initialize Anaconda3 by running conda init? [yes|no]
> [yes] >>>
$ yes
```

Read and accept the terms. Take note of the location for the installation directory. You can use the default location or change it if you would like. A good installation location is: `/Users/username/anaconda3`.

* In Windows, we recommend installing the pipeline under a Linux-type environment using [Cygwin](https://www.cygwin.com/) so that you will be able to use the pipeline as command line tools. Simply download the Windows version of Anaconda on the [Anaconda website](https://www.anaconda.com/distribution/#download-section) and follow the instructions from the installation GUI. When asked about the installation location, make sure you install Python under your emulated-Linux home directory (`/home/username`) and ***not*** in the default location (`/cygdrive/c/Users/username/anaconda3`). From the installation GUI, the path you want to select is something like `C:/Program Files/cygwin64/home/username/anaconda3` Also, make sure to check **YES** for "Add Anaconda to my `PATH` environment variable."

To make sure that your path to the Anaconda Python distribution is fully actualized, we recommend you close the current terminal, open a fresh terminal, and type 
```bash
$ python [TAB]
```

If multiple options are available (e.g. `python`, `python2`, `python 3.7`, `python.exe`), this means that you have other versions of Python sitting on your system (an old `python2` executable located in `/usr/local/bin/python`, for example).
```bash
$ python3 --version     # in bash, csh OR
$ python.exe --version  # in Cygwin/Windows
```

Do this also for the `pip command` (e.g. old  `pip`, `pip3`, `pip.exe`). Then, set your `$PATH` environment variable to point to the Anaconda Python and pip distributions we just installed and confirm this with the `which` command:
```bash
$ which python3         # in bash, csh OR
$ which python.exe      # in Cygwin/Windows
```

We are looking for a Python executable that looks like it was installed with Anaconda, such as
```bash
/username/anaconda3/bin/python3 # on MacOS/Linux, OR
/username/anaconda3/python.exe # on Cygwin/Windows
```

If `which` points to either of those locations, you are good to go and you can proceed from here using the shorthand path to your Anaconda Python distribution:
```bash
$ python3     # Linux/MacOS 
$ python.exe  # Cygwin/Windows
```

If, however, `which` points to some other location, such as `/usr/local/bin/python`, proceed from here using the FULL path to the Anaconda Python distribution like so:
```bash
$ /Users/username/anaconda3/bin/python3 # Linux/MacOS
$ /Users/username/anaconda3/python.exe  # Cygwin/Windows
```

***

### Optional Installation of `ghostscript`
If you like, you can install `ghostscript` on your local machine to allow CAP to generate multiple figures in a single PDF file. CAP supports the ability to create individual PNGs if you prefer.

To check whether you have `ghostscript` installed already, open a terminal and type
```bash
$ gs -version
``` 
If it is not installed, follow the directions on the `ghostscript` [website](https://www.ghostscript.com/download.html).

***

## Creation of the Virtual Environment {#Create_venv}

Next, we will create a virtual environment in which to download CAP. The virtual environment shares the your main Python core but branches out with its own packages. We will name it `amesGCM3` to remind ourselves that this environment shares the core Python3 structure it is derived from. To create the virtual environment, open a terminal and type:
```bash
$ python3 -m venv --system-site-packages amesGCM3 # remember to use FULL PATH to your python distribution if needed!
```

Here is the virtual environment `amesGCM3` that you just created:
```bash
anaconda3                  amesGCM3/
├── bin                    ├── bin
│   ├── pip       (copy)   │    ├── pip
│   └── python3    >>>>    │    ├──python3
└── lib                    │    ├── activate
                           │    ├── activate.csh
                           │    └── deactivate
                           └── lib             

  MAIN ENVIRONMENT           VIRTUAL ENVIRONMENT
(Leave untouched for)       (OK to mess around, will vanish
this particular project)     everytime we run 'deactivate')

```

You can activate the virtual environment by sourcing it:

```bash
$ source amesGCM3/bin/activate      # if you are using bash
$ source amesGCM3/bin/activate.csh  # if you are using csh/tcsh
```

> Note that in Cygwin/Windows, the `/bin` directory may be named `/Scripts`.

You may notice that after sourcing `amesGCM3`, your prompt changed to `(amesGCM3)>$`. This confirms that you are **inside** the virtual environment even when you navigate to different directories on your machine.

After sourcing the virtual environment, we can verify that `which python` and `which pip` unambiguously point to `amesGCM3/bin/python3` and `amesGCM3/bin/pip`, respectively. There is therefore no need to reference their full paths for the following instructions.

## Installing the Pipeline {#Install_CAP}

##### Directly from Github {#Install_CAP_GitHub}

From *inside* the virtual environment, `amesGCM3`, run:
```bash
(amesGCM3)>$ pip install git+https://github.com/alex-kling/amesgcm.git
```

##### From a .zip Archive  {#Install_CAP_zip}

If you have been provided with an archive, download and untar the `amesgcm-master.zip` archive wherever it is on your machine (e.g. in `/Downloads`). From **inside** the virtual environment, type:
```bash
(amesGCM3)>$ cd amesgcm-master
(amesGCM3)>$ pip install .
```

It is now safe to move (or remove) both the `amesgcm-master` source code and the `.zip` archive from your `/Downloads` directory since `pip` installed the pipeline inside your `amesGCM3` virtual environment.

***

To make sure the paths to the executables are correctly set in your terminal, exit the virtual environment with:
```bash
(amesGCM3)>$ deactivate
$
```

This completes the one-time installation of CAP in your virtual environment, `amesGCM3`, which now looks like:

```bash
amesGCM3/
├── bin
│   ├── MarsFiles.py
│   ├── MarsInterp.py
│   ├── MarsPlot.py
│   ├── MarsPull.py
│   ├── MarsVars.py
│   ├── activate
│   ├── activate.csh
│   ├── deactivate
│   ├── pip
│   └── python3
├── lib
│   └── python3.7
│       └── site-packages
│           ├── netCDF4
│           └── amesgcm
│               ├── FV3_utils.py
│               ├── Ncdf_wrapper.py
│               └── Script_utils.py
├── mars_data
│   └── Legacy.fixed.nc
└── mars_templates
    ├──amesgcm_profile
    └── legacy.in
```
***

**A Note:**
CAP requires the following Python packages. These were installed automatically when you installed CAP:
* `numpy` for array operations
* `matplotlib` for the MatPlotLib plotting library
* `netCDF4 Python` for handling netCDF files
* `requests` for downloading data from the MCMC Portal

**Quick Tip:**
If you prefer using the `conda` package manager for setting up your virtual environment instead of `pip`, you may use the following commands to install CAP.

First, verify (using `conda info` or `which conda`) that you are using the intented `conda` executable (two or more versions of `conda` might be present if both Python2 and Python3 are installed on your system):
```bash
$ conda create -n amesGCM3
$ conda activate amesGCM3
$ conda install pip
$ pip install git+https://github.com/alex-kling/amesgcm.git
```

The source code will be installed in `your_conda/envs/amesGCM3/`. The pipeline can then be activated and exited with
```bash
$ conda activate amesGCM3
(amesGCM3)>$ conda deactivate
$
```

***

## Routine use of the pipeline  {#Routine_Use}

Every time you want to use CAP from a new terminal session, simply source the `amesGCM3` virtual environment:
```bash
$ source amesGCM3/bin/activate      # in bash OR
$ source amesGCM3/bin/activate.csh  # in csh/tcsh
```

You can always check that the `Mars****.py` tools are installed properly by typing:
```bash
(amesGCM3)>$ Mars[TAB]
> MarsFiles.py   MarsInterp.py  MarsPlot.py    MarsPull.py    MarsVars.py
```

If no executables show up, then the paths to the executibles have not been properly set up in the virtual environment and you must use the full paths to the executables directly:

```bash
(amesGCM3)>$ ~/amesGCM3/bin/Mars[TAB]
```
or set your own aliases in `~/.bashrc` or `~/.bash_profile`:
```bash
$ vim ~/.bashrc
alias MarsPlot.py='/username/amesGCM3/bin/MarsPlot.py'  # add this line
:wq                                                     # save and quit ~/.bashrc
$ source ~/.bash_profile                                # source your ~/.bashrc
```

or in `~/.cshrc`:
```csh
$ vim ~/.cshrc
alias MarsPlot.py /username/amesGCM3/bin/MarsPlot.py    # add this line
:wq                                                     # save and quit ~/.cshrc
source ~/.cshrc                                         # source your ~/.cshrc
```

You can check the documentation for any of the executables with the `--help` option:
```bash
(amesGCM3)>$ MarsPlot.py --help  # or
(amesGCM3)>$ MarsPlot.py -h      # for short
```

Remember, you can exit the virtual environment using:
```bash
(amesGCM3)>$ deactivate
```

## Upgrade or Remove CAP {#Upgrade}
To upgrade your version of CAP to the most recent release, activate `amesGCM3` and run the `upgrade` command:
```bash
$ source amesGCM3/bin/activate      # in bash, OR
$ source amesGCM3/bin/activate.csh  # in csh/tcsh
(amesGCM3)>$ pip install git+https://github.com/alex-kling/amesgcm.git --upgrade
```

To permanently remove CAP, activate `amesGCM3` and run :
```bash
$ source amesGCM3/bin/activate      # in bash, OR
$ source amesGCM3/bin/activate.csh  # in csh/tcsh
(amesGCM3)>$ pip uninstall amesgcm
```

It is also safe to delete the entire `~/amesGCM3` directory in order to delete CAP and your virtual environment. Doing so will not alter your main Python distribution.

***

# TUTORIAL {#TUTORIAL}
## Overview of Using CAP for Analysis {#Overview}

The following steps show you how to use CAP to access Mars GCM data, reduce it, compute diagnostics, interpolate diagnostics to standard pressures levels, and visualize the results. A quick overview of what these steps are is shown in the graphic below.
![](./docs/cheat_sheet.png)

## Download Raw Legacy GCM Outputs {#Legacy_Download}
The data from the Legacy GCM are available for download on the [MCMC Data portal](https://data.nas.nasa.gov/legacygcm/data_legacygcm.php). The data are archived in 1.5 hours intervals (i.e. 16 times a day), packaged in sets of 10 sols (1 sol = 1 Martian day), and referenced by their solar longitude ("_Ls_"). Ls=0° is northern vernal equinox (beginning of Northern spring), Ls=90° is northern summer solstice, Ls=180° is northern autumnal equinox, and Ls=270° is northern winter solstice. 

Let's begin by downloading 30-sols of data starting the beginning of the Martian year (Ls=0° to Ls=15°). First, source the virtual environment. Then, navigate to the directory in which you would like to store the data and run:
```bash
(amesGCM3)>$ MarsPull.py --help      # to see what options MarsPull offers
(amesGCM3)>$ MarsPull.py --ls 0 15   # to extract the data
```

The second command will download three `LegacyGCM_Ls000***.nc` raw outputs, each ~280 MB in size. We can use the `--inspect` command from `MarsPlot.py` to see the content in the files:

```bash
(amesGCM3)>$ MarsPlot.py -i LegacyGCM_Ls000_Ls004.nc
```

Note the raw output is provided in 10-day chunks (`time`) 16 times per day (`ntod`).

## File format conversion
For analysis purposes, it is useful to reduce the data from the raw output into other formats. Format options are:
* `fixed`    static fields (e.g. surface albedo, topography)
* `average`  5-day averages
* `daily`    continuous time series
* `diurn`    5-day averages preserving each of the 16 times per day

These file formats are the new standard output format used by our newer GCM. To create the file formats listed above, we use `MarsFiles`:
```bash
(amesGCM3)>$ MarsFiles.py -h # for information about MarsFiles.py
(amesGCM3)>$ MarsFiles.py LegacyGCM_Ls* -fv3 fixed average # to create new file formats
```

Once again, use `inspect` to view the content of each of these files:
```bash
(amesGCM3)>$ MarsPlot.py -i 00000.fixed.nc
(amesGCM3)>$ MarsPlot.py -i 00000.atmos_average.nc
```

You can use cap with individual sets of files (`00000`, `00010`, and `00020` files in our example), or merge those files together into one using `MarsFiles.py` and the `--combine` call. All the utilities from the analysis pipeline (including the plotting routine) accept a list of files as input and keeping separate files can be strategic when computer memory is limited (the `daily` files are 280 MB each and there are 67 of them in one Mars year).

Since working with 5-day averages requires relatively small files, we can use `--combine` to merge several of them together along the `time` dimension, like so:
```bash
(amesGCM3)>$ MarsFiles.py *fixed.nc -c
(amesGCM3)>$ MarsFiles.py *atmos_average.nc -c
```
> Note that if you delete the `scalar_axis` variable, concatenating files will not work.

We can use `--dump` (or `--stat`) with `MarsPlot.py` to inspect the changes made to the time dimension:
```bash
(amesGCM3)>$ MarsPlot.py -i 00000.atmos_average.nc
(amesGCM3)>$ MarsPlot.py -i 00000.atmos_average.nc -dump time areo
```

## Variable Operations {#Variable_ops}
When no arguments are provided, the variable utility `MarsVars.py` has the same functionality as `MarsPlot.py -i`:
```bash
(amesGCM3)>$ MarsVars.py 00000.atmos_average.nc
```

To see what else `MarsVars.py` can do, use the `--help` option:
```bash
(amesGCM3)>$ MarsVars.py -h
```

After performing the above command, you should see that `MarsVars.py` provides an option that adds atmospheric density (rho) to the file as a derived variable calculated from the vertical grid (`pk`, `bk`), surface pressure (`ps`), and air temperature (`temp`) data. To test this utility, run:
```bash
(amesGCM3)>$ MarsVars.py 00000.atmos_average.nc -add rho
```

Check that `rho` was added to the file by running `MarsVars.py` again with no argument:
```bash
(amesGCM3)>$ MarsVars.py 00000.atmos_average.nc
```

Similarly, we can perform a column integration on the water vapor variable (`vap_mass`) using `-col`. At the same time, we will remove the dust variable (`dst_num`) and the water ice (`ice_num`) variable from the file. Since we are not planning on using these variables in our analysis, this will free up some memory.
```bash
(amesGCM3)>$ MarsVars.py 00000.atmos_average.nc -col vap_mass -rm ice_num dst_num
```

## Pressure Interpolation {#pinterp}

The Legacy GCM uses a pressure coordinate (`pfull`) in the vertical which means that a single atmospheric layer will be located at different geometric heights (and pressure levels) between the atmospheric columns. Before we do any zonal averaging, it is therefore necessary to interpolate the data to standard pressure surfaces. This operation can be done using `MarsInterp` with the `--type pstd` option:
```bash
(amesGCM3)>$ MarsInterp.py -h
(amesGCM3)>$ MarsInterp.py  00000.atmos_average.nc -t pstd
```

Inspecting the file with:
```bash
(amesGCM3)>$ MarsPlot.py -i 00000.atmos_average_pstd.nc
```
we can see that the vertical coordinate `pfull` (formerly 24 layers) has been replaced by a standard pressure coordinate `pstd`. Also, the shape of the 3D variables reflect the new shape of `pstd`.

## Plotting the results with MarsPlot {#Plotting}

While you may use the software of your choice to visualize the results (e.g. Matlab, IDL), a utility is provided in CAP that creates 2D figures and 1D line plots that are easily configured from an input template. To generate a plotting template (`Custom.in`) in the current directory, use:
```bash
(amesGCM3)>$ MarsPlot.py -h         # for instructions
(amesGCM3)>$ MarsPlot.py --template # to create the template
```

Open `Custom.in` with a text editor (you can rename the file to `something.in` if you want). You can skip the commented instructions at the top of `Custom.in` if you're following along using these instructions. Go directly to the section:
```python3
 <<<<<<<<<<<<<<<<<<<<<< Simulations >>>>>>>>>>>>>>>>>>>>>
 ref> None
 2> /u/username/path/to/GCM/data
 3>
 =======================================================
 START
```

>**Quick Tip:** `MarsPlot` uses text files with a '`.in`' extension as input files. Select "Python" as the language (in place of "Plain text") when editing the file from text editor (gedit, atom ...) to enable syntax-highlighting of keywords. If you are using `vim`, add the following lines to your `~/.vimrc` so that vim recognizes that `Custom.in` is Python code:
> ```bash
> $ vim ~/.vimrc
> # add the following lines:
> syntax on
> colorscheme default
> au BufReadPost \*.in  set syntax=python
> # write and quit the editor:
> :wq
> ```
> The next time you open `Custom.in` with `vim`, numbers and keywords like `True` or `None` will be highlighted appropriately.

To access data stored in a specific file, the `Main Variable` line accepts the syntax:
```python3
Main Variable  = XXXXX.file@N.var
```
* `XXXXX` is the sol number (e.g. "`03335`", optional)
* `file` is the file type (e.g. "`atmos_average_pstd`")
* `@N` points to the simulation listed in the section above (e.g. `@2` refers to the simulation listed under `2> /u/username/path/to/GCM/data`)
* `var` is the requested variable (e.g. `ucomp` for zonal wind)

When dimensions are omitted with `None` (e.g. `Ls 0-360 = None` or `Level [Pa/m] = None`), `MarsPlot.py` makes educated guesses for data selection. For example, if no layer is requested, CAP defaults to the surface layer. CAP will tell you exactly how the data is being processed both in the default **title of the figure** and in the **terminal output**. This behavior is detailed in the commented instructions at the top of `Custom.in`, which also includes descriptions of features such as the use of square brackets `[ ]` for variable operations and squiggly brackets `{ }` for specifying dimensions. You can use output from several simulations by specifying the directories the data are in under `<<<<< Simulations >>>>>`.

If `ghostscript` is available on your system, then you can feed the `Custom.in` template back into `MarsPlot.py` using:
```bash
(amesGCM3)>$ MarsPlot.py Custom.in
```

If you do not have `ghostscript`, use:
```bash
(amesGCM3)>$ MarsPlot.py Custom.in -o png
```

When running `MarsPlot.py`, you will see the following in the terminal:
```python3
[----------]  0 % (2D_lon_lat :fixed.zsurf)
[#####-----] 50 % (2D_lat_press :atmos_average.ucomp, Ls= (MY 1) 13.61, lon=18.0)
[##########]100 % (Done)
```

> By default, `MarsPlot.py` handles errors such as missing data by itself. It reports them after completing the plotting routine both in the **terminal** and in the relevant **figure**. To bypass this behavior, use the  `--debug` option when running `MarsPlot.py`.

Unless `-o png` is specified, `MarsPlot.py` creates a file called `Diagnostic.pdf` containing the requested plots in the current directory. `Diagnostics.pdf` can be opened with a PDF viewer such as MacOS Preview:
```bash
$ open -a Preview Diagnostic.pdf
```
The Preview App conveniently reloads PDFs automatically. The 'Skim' editor, available for download, offers a similar feature in `Preferences/Sync`. On Linux machines, view the PDF using:
```bash
$ evince Diagnostic.pdf
```
If you have used the `--output png` formatting option, the images will be located in a new directory called `/plots` created in the current directory.

You can add a new figure to the PDF by making a copy of any of the `<<<| Plot ... = True |>>>` blocks below the `HOLD ON[...]HOLD OFF` statement. The `HOLD` commands place multiple figures on the same page. To place figures on their own pages, do not wrap multiple plotting routines (`<<<| Plot ... = True |>>>`) between `HOLD ON[...]HOLD OFF`. 

For example, plot the zonal averaged and time averaged dust field for the first 10 degrees of Ls from the interpolated file `atmos_average_pstd`:

```python3
<<<<<<<<<<<<<<| Plot 2D lat X lev = True |>>>>>>>>>>>>>
Title          = This is the dust field converted to [g/kg]
Main Variable  = [atmos_average_pstd.dst_mass]*1000.  # pointing to the interpolated file, call the dust variable & convert the units to g/kg
Cmin, Cmax     = None
Ls 0-360       = 0.,10                                # plot the data from Ls 0-10 only
Lon +/-180     = all                                  # all = zonal average
2nd Variable   = None
Contours Var 2 = None
Axis Options  : Lat = [None,None] | level[Pa] = [1e3,0.2] | cmap = Wistia | scale = lin
```

Note that the square brackets `[ ]` around the variable are used to perform calculations on the variable. In this case, to change the variable units to [g/kg] from [kg/kg]. The title was also changed to reflect the variable plotted. We specify the colormap (`Wistia`) and adjust the vertical range of the axis (`level[Pa] = [1e3,0.2]`) in `Axis Options`. After saving the template, we feed it back to `MarsPlot.py`: 
```bash
(amesGCM3)>$ MarsPlot.py Custom.in
```

By default, `MarsPlot.py` runs the requested analysis on the *last* set of output files present in the directory (i.e. the highest dated file, `XXXXX.fixed.nc`). To run the analysis on a specific file or a range of files, use the `--date` option:

```bash
(amesGCM3)>$ MarsPlot.py Custom.in -d 0
```

Close and open the PDF again and you should see a new figure of the updated dust field. You can use `Custom.in` jointly with `MarsPlot.py --inspect` to add new figures and to explore the other plot types outlined in `Custom.in`. By default, plot templates are set to `= False` but can be enabled with `= True`.

# Moving Forward with Your Own Analysis {#Your_Analysis}

You can customize your own plots using the programming language of your choice. Here is a Python script to get you started. Unless you have installed python-netCDF4 and CAP on top of your main distribution, this script has to be be run from within the `amesGCM3` virtual environment. This will allow access to the `netCDF4` and `amesGCM3` packages.

Create a file called `demo.py` and copy & paste the following within it:

```python3
#======================= Import python packages ================================
import numpy as np                          # for array operations
import matplotlib.pyplot as plt             # python plotting library
from netCDF4 import Dataset                 # to read .nc files
#===============================================================================

# Open a fixed.nc file, read some variables and close it.
f_fixed=Dataset('/Users/akling/test/00000.fixed.nc','r')
lon=f_fixed.variables['lon'][:]
lat=f_fixed.variables['lat'][:]
zsurf=f_fixed.variables['zsurf'][:]  
f_fixed.close()

# Open a dataset and read the 'variables' attribute from the NETCDF FILE
f_average_pstd=Dataset('/Users/akling/test/00000.atmos_average_pstd.nc','r')
vars_list     =f_average_pstd.variables.keys()
print('The variables in the atmos files are: ',vars_list)

# Read the 'shape' and 'units' attribute from the temperature VARIABLE
Nt,Nz,Ny,Nx = f_average_pstd.variables['temp'].shape
units_txt   = f_average_pstd.variables['temp'].units
print('The data dimensions are Nt,Nz,Ny,Nx=',Nt,Nz,Ny,Nx)
# Read the pressure, time, and the temperature for an equatorial cross section
pstd       = f_average_pstd.variables['pstd'][:]   
areo       = f_average_pstd.variables['areo'][0] #solar longitude for the 1st timestep
temp       = f_average_pstd.variables['temp'][0,:,18,:] #time, press, lat, lon
f_average_pstd.close()

# Get the latitude of the cross section.
lat_cross=lat[18]

# Example of accessing  functions from the Ames Pipeline if we wanted to plot
# the data  in a different coordinate system  (0>360 instead of +/-180 )
#----
from amesgcm.FV3_utils import lon180_to_360,shiftgrid_180_to_360
lon360=lon180_to_360(lon)
temp360=shiftgrid_180_to_360(lon,temp)

# Define some contours for plotting
conts= np.linspace(150,250,32)

#Create a figure with the data
plt.close('all')
ax=plt.subplot(111)
plt.contourf(lon,pstd,temp,conts,cmap='jet',extend='both')
plt.colorbar()
# Axis labeling
ax.invert_yaxis()
ax.set_yscale("log")
plt.xlabel('Longitudes')
plt.ylabel('Pressure [Pa]')
plt.title('Temperature [%s] at Ls %03i, lat= %.2f '%(units_txt,areo,lat_cross))
plt.show()
```

Then run `demo.py` using:
```bash
(amesGCM3)>$ python3 demo.py
```

`demo.py` should produce the following image:

![](./docs/demo_figure.png)
