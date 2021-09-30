Welcome to the Mars Climate Modeling Center (MCMC) Analysis Pipeline. By the end of this tutorial, you will know how to download Mars Climate data from the MCMC's data portal, reduce these large climate simulations to meaningful data, and make plots for Martian winds, temperature, and aerosols at specific seasons and locations.

The simulation results presented on this page are extensively documented in [Haberle et al. 2019 _Documentation of the NASA/Ames Legacy Mars Global Climate Model: Simulations of the present seasonal water cycle_ ](https://www.sciencedirect.com/science/article/pii/S0019103518305761)

# INSTALLATION:
The analysis pipeline is entirely written in pure Python, which is an intuitive and open source programming language. You may identify yourself in one the following categories:

* A. You are familiar with the Python infrastructure and would like to install the Analysis pipeline on top of your current Python installation: Check the requirements below and skip to '**_Installing the pipeline_**'. Note that you may have to manually add aliases to the _Mars***.py_ executables to your search path.
* B. You have experience with Python but not with managing packages, or are new to Python: To ensure that there is no conflict with other Python versions that may be on your system, we will install a fresh Python distribution **locally** (this does not require admin permission). Additionally, we will install the analysis pipeline in a self-contained _virtual environment_  which is basically a standalone copy of your entire distribution, minus the 'core' code that is shared with the main python distribution.  This will allow you to use your fresh Python installation for other projects (including installing or upgrading packages) without the risk of altering the analysis pipeline. It will also be safe to alter (or even completely delete) that virtual environment without breaking the main distribution.
##  Requirements

**Python 3**: It you are already a Python user, you can install the Ames analysis pipeline on top of you current installation. For new users, we recommend to use the latest version of the Anaconda Python distribution as it already ships with pre-compiled math and plotting packages (e.g. _numpy_, _matplotlib_), and pre-compiled libraries (e.g. hdf5 headers to read netcdf files). You can download either the command-line installer or the graphical interface [here](https://www.anaconda.com/distribution/#download-section). If you are the owner of the system you can choose to install Python at the system level, but you can install it in your home directory if you don't have permission on your system.

* In MacOS and Linux, you can install a fresh Python3  **locally** from a terminal with:

`chmod +x Anaconda3-2020.02-MacOSX-x86_64.sh` (this make the .sh file executable)

`./Anaconda3-2020.02-MacOSX-x86_64.sh`          (this runs the executable)

Read (ENTER) and accept (yes) the terms. Take note of the location for the installation directory. You can use the default location or change it if you would like, for example _/Users/username/anaconda3_ works well.

* In Windows, we recommend installing the pipeline under a Linux-type environment using [Cygwin](https://www.cygwin.com/), so we will be able to use the pipeline as command-line tools. Simply download the _Windows_ version on the [Anaconda website](https://www.anaconda.com/distribution/#download-section) and follow the instructions from the installation GUI. When asked about the installation location, make sure you install Python under your emulated-Linux home _/home/username_ and **not** on the default location _/cygdrive/c/Users/username/anaconda3_. From the installation GUI, the path you want to select is something like: `C:/Program Files/cygwin64/home/username/anaconda3` Also make sure to check YES  for _Add Anaconda to my PATH environment variable_

***

The analysis pipeline requires the following **Python packages** which will be installed automatically in the analysis pipeline virtual environment (more on this later):
* **numpy** (array operations)
* **matplotlib** (plotting library)
* **netCDF4 Python** (handling of netcdf files)
* **requests** (for downloading data from the MCMC Portal)

Optionally, you can install:
* **ghostscript** which will allow the analysis pipeline to generate multiple figures as a single pdf file. Type
`gs -version ` to see if Ghostscript is already available on your system. If it is not, you can installing from this page: [https://www.ghostscript.com/download.html](https://www.ghostscript.com/download.html) or decide to use png images instead.

To make sure the paths are fully actualized, we recommend to close the current terminal. Then, open a fresh terminal, type `python` and hitting the TAB key. If multiple options are be available (e.g. _python_, _python2_, _python 3.7_, _python.exe_), this means that you have other versions of Python sitting on your system (e.g. an old _python2_  executable located in _/usr/local/bin/python_ ...). The same holds true for the `pip command` (e.g. old  _pip_, _pip3_, _pip.exe_). Pick the one you think may be from the Anaconda version you just installed and confirm this with the _which_ command, for example:


`python3 --version`    (`python.exe --version` in Cygwin/Windows)

`which python3`        (`which python.exe` in Cygwin/Windows)

We are looking for a python executable that looks like it was installed with Anaconda, like _/username/anaconda3/bin/python3_ (MacOS/Linux)  or _/username/anaconda3/python.exe_ (Cygwin/Windows)   If `which python3` or `which python.exe` already points to one of those locations, you are good to go. If  `which` points  to some other location (e.g. _/usr/local/bin/python_) proceed with the FULL paths to the Anaconda-Python, e.g.

`/Users/username/anaconda3/bin/python3` instead of `python3` (Linux/MacOS) or `/Users/username/anaconda3/python.exe` instead of `python.exe` (Cygwin/Windows)


## Creation of a virtual environment:
We will create a virtual environment for the Ames analysis pipeline which shares the same Python core but branches out with its own packages. We will call it  _amesGCM3_ to remind ourselves that it shares the same structure that the core python3 it is derived from. From a terminal run:

`python3 -m venv --system-site-packages amesGCM3` (remember to use FULL PATH to `python` if needed)

Here is what just happened :

```
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

Now activate the virtual environment with:

`source amesGCM3/bin/activate`      (if you are using **bash** )

`source amesGCM3/bin/activate.csh`  (if you are using **csh/tcsh** )

Note that in Cygwin/Windows, the  _bin_ directory may be named **_Scripts_**


You may notice that your prompt changed from _username>_ to _(amesGCM3)username>_ which indicates that you are INSIDE the virtual environment, even when navigating to different directories on your machine.  

After entering the virtual environment, we can verify that ```which python ``` and ```which pip ``` unambiguously point to _amesGCM3/bin/python3_ and _amesGCM3/bin/pip_ so there is no need use the full paths.

## Installing the pipeline
##### Directly from Github:

From **inside** the virtual environment, run:

`pip install git+https://github.com/alex-kling/amesgcm.git`

##### From a .zip archive:
If you have been provided with an archive, download and untar the '_amesgcm-master.zip_' archive anywhere (e.g. in our _Downloads_ directory). From **inside** the virtual environment, run:

```
cd amesgcm-master
pip install .
```
It is safe to move (or remove) both the '_amesgcm-master_' source code, and the '.zip' archive from your _Downloads_ directory since _pip_ installed the pipeline inside your _~/amesGCM3_ virtual environment.
***
To make sure the paths to the executables are correctly set in your terminal, exit the virtual environment with

`deactivate`

This complete the one-time installation of the Ames analysis pipeline:
```
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
**Quick Tip:**
If you prefer to use the _conda_ package manager to set-up your virtual environment instead of _pip_ , you may use the following commands instead to install the pipeline:
Verify with ```conda info``` or  ```which conda``` that you are using the intented _conda_ executable (two or more versions of _conda_ might be present if both python 2 and python 3 are installed on your system). Then install the pipeline with:

```
conda create -n amesGCM3
conda activate amesGCM3
conda install pip
pip install git+https://github.com/alex-kling/amesgcm.git
```

The source code will be installed in _your_conda/envs/amesGCM3/_ . The pipeline can then be activated with ```conda activate amesGCM3``` and exited with ```conda deactivate``` (more on this below).

***

## Routine use of the pipeline

Every time you want to use the analysis pipeline from a new terminal session, simply run:

`source amesGCM3/bin/activate`    (`source amesGCM3/bin/activate.csh` in **csh/tcsh**)

You can check that the tools are installed properly by typing `Mars` and hit the **TAB** key. No matter where you are on your system, you should see the following:
```
(amesGCM3) username$ Mars
MarsFiles.py   MarsInterp.py  MarsPlot.py    MarsPull.py    MarsVars.py
```
If no executable show up, the paths have not not been set-up in the virtual environment. You can use the full paths to the executable e.g. instead of using  `MarsPlot.py`, use `~/amesGCM3/bin/MarsPlot.py`. If that solution works for you, also consider setting-up your own aliases, for example:

Add `alias MarsPlot.py='/username/amesGCM3/bin/MarsPlot.py'` to your _~/.bash_profile_  and run
`source ~/.bash_profile` (in **bash**)

Add `alias MarsPlot.py /username/amesGCM3/bin/MarsPlot.py` to your _~/.cshrc_  and run
`source ~/.cshrc` (in **csh**)  

Check the documentation for any of the executables above with the `--help` option:

`MarsPlot.py --help` (`MarsPlot.py -h` for short)


After you are done with your work, you can exit the analysis pipeline with:

 `deactivate`

## Upgrade, or remove the pipeline
To upgrade the pipeline, activate the virtual environment:

`source amesGCM3/bin/activate`    (`source amesGCM3/bin/activate.csh` in **csh/tcsh**)

And run:

`pip install git+https://github.com/alex-kling/amesgcm.git --upgrade`

To permanently remove the amesgcm pipeline, activate the virtual environment and run :

`pip uninstall amesgcm`

It is also safe to delete the entire _amesGCM3_ virtual environment directory as this will not affect your main Python distribution.


## Download raw Legacy GCM outputs
The data from the Legacy GCM is archived every 1.5 hours (i.e. 16 times a day) and packaged in chunks of 10 sols (1 sol = 1 martian day). Files are available for download on the MCMC Data portal at : [https://data.nas.nasa.gov/legacygcm/data_legacygcm.php](https://data.nas.nasa.gov/legacygcm/data_legacygcm.php), and referenced by their solar longitude or "_Ls_", which  is 0° at the vernal equinox (beginning of Northern spring), 90° during the summer solstice, 180° at the autumnal equinox, and 270° during winter solstice. To download a 30-sols chunk starting at the beginning at the martian year (Ls =0 to Ls=15), navigate to a place you would like to store the data and run :
```
MarsPull.py --help
MarsPull.py --ls 0 15
```
This will download three LegacyGCM_Ls000***.nc raw outputs, each ~280MB each.

We can use the **--inspect** command of MarsPlot.py to peak into the content for one of the raw outputs:

`MarsPlot.py -i LegacyGCM_Ls000_Ls004.nc`

Note the characteristic structure for the Legacy GCM raw outputs with 10 days chunks ('_time_') , and 16 time of day ('_ntod_').

## File format conversion
For analysis purposes, it is useful to reduce the data from the raw outputs into different formats:

* **fixed**:   static fields (e.g. surface albedo, topography)
* **average**: 5 days averages
* **daily** : continuous time series
* **diurn** : 5 days average for each time of the day

New files for each of the formats listed above can be created using the _MarsFiles_ utility which handles conversions from the Legacy format to this new (FV3) format. To create _fixed_ and _average_ files for each of the 10 days output from the Legacy GCM, run:

```
MarsFiles.py -h
MarsFiles.py LegacyGCM_Ls* -fv3 fixed average
```
And check the new content for one of the files with:
```
MarsPlot.py -i 00000.fixed.nc
MarsPlot.py -i 00000.atmos_average.nc
```

Moving forward with the postprocessing pipeline, it is the user's choice to proceed with individual sets of files (00000, 00010, and 00020 files in our example), or merge those files together into one.
All the utilities from the analysis pipeline (including the plotting routine) accept a **list** of files as input, and keeping separate files can be strategic when computer memory is limited (the **daily** files remain 280MB each and there are 67 of those in one Mars year).

Since working with 5 days average involves relatively small files, we can use the **--combine** option of _MarsFiles_ to merge them together along the '_time_' dimension:

```
MarsFiles.py *fixed.nc -c
MarsFiles.py *atmos_average.nc -c
```
We can use the **--dump** (or **--stat**) option of MarsPlot to inspect the changes made to the time dimension:

```
MarsPlot.py -i 00000.atmos_average.nc
MarsPlot.py -i 00000.atmos_average.nc -dump time areo
```
