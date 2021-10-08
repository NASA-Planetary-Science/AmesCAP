![](./tutorial_images/Tutorial_Banner_Final.png)

***

# Intstalling the Community Analysis Pipeline (CAP)

## Welcome to the installation of CAP!

These are the instructions for installing the NASA Ames MCMC's Community Analysis Pipeline (CAP). **We ask that you come to the MGCM Tutorial on November 2-4 with CAP installed on your machine** so that we can jump right into using it! We'll be using the tools in CAP to analyze MGCM output on the second day of the Tutorial. Installing CAP is fairly straightforward. We will create a python virtual environment, download CAP, and then install CAP within that virtual environment. That's it!






## Creating the Virtual Environment

We'll begin by creating a virtual environment in which to install CAP. The virtual environment is an isolated Python environment cloned from an existing Python distribution. The virtual environment consists of the same directory trees as the original environment, but it includes activation and deactivation scripts that are used to move in and out of the virtual environment. Here's an illustration of how the two Python environments might look:

```
     anaconda3                    virtual_env3/
     ├── bin                      ├── bin
     │   ├── pip       (copy)     │    ├── pip
     │   └── python3    >>>>      │    ├── python3
     └── lib                      │    ├── activate
                                  │    ├── activate.csh
                                  │    └── deactivate
                                  └── lib             

  ORIGINAL ENVIRONMENT           VIRTUAL ENVIRONMENT
      (Untouched)            (Vanishes upon deactivation)
```

We can install and upgrade packages in the virtual environment without the risk of altering CAP or breaking the main Python environment. It is safe to change or even completely delete the virtual environment without breaking the main distribution. This allows us to experiment freely in the virtual environment.





### Step 1: Identify Your Preferred Python Distribution

We highly recommend using the latest version of the Anaconda Python distribution. It ships with pre-compiled math and plotting packages such as `numpy` and `matplotlib` as well as pre-compiled libraries like `hdf5 headers` for reading `netCDF` files (the preferred filetype for analysing MGCM output). 

You can install the Anaconda Python distribution via the command-line or a [graphical interface](https://www.anaconda.com/distribution/#download-section). You may install Anaconda at either the `System/` or `User/` level. For command-line installation, open a terminal and type the following:

```bash
(local)>$ chmod +x Anaconda3-2021.05-MacOSX-x86_64.sh   # creates the .sh file executable
(local)>$ ./Anaconda3-2021.05MacOSX-x86_64.sh           # runs the executable
```

Which will return the following:

```bash
> Welcome to Anaconda3 2021.05
>
> In order to continue the installation process, please review the license agreement.
> Please, press ENTER to continue
> >>> 
```

Read (ENTER) and accept (yes) the terms, choose your installation location, and initialize Anaconda3:

```bash
(local)>$ [ENTER]
> Do you accept the license terms? [yes|no]
> [no] >>>
(local)>$ yes
> Anaconda3 will now be installed into this location:
> /Users/username/anaconda3
> 
>  - Press ENTER to confirm the location
>  - Press CTRL-C to abort the installation
>  - Or specify a different location below
>
> [/Users/cbatters/anaconda3] >>> 
(local)>$ [ENTER]
> PREFIX=/Users/cbatters/anaconda3
> Unpacking payload ...
> Collecting package metadata (current_repodata.json):
>   done                                                       
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
(local)>$ yes
```

> For Windows users, we recommend installing the pipeline under a Linux-type environment using [Cygwin](https://www.cygwin.com/) so that you will be able to use the pipeline as command line tools. Simply download the Windows version of Anaconda on the [Anaconda website](https://www.anaconda.com/distribution/#download-section) and follow the instructions from the installation GUI. When asked about the installation location, make sure you install Python under your emulated-Linux home directory (`/home/username`) and ***not*** in the default location (`/cygdrive/c/Users/username/anaconda3`). From the installation GUI, the path you want to select is something like `C:/Program Files/cygwin64/home/username/anaconda3` Also, make sure to check **YES** for "Add Anaconda to my `PATH` environment variable."

Confirm that your path to the Anaconda Python distribution is fully actualized by closing the current terminal, opening a fresh terminal, and typing:

```bash
(local)>$ python[TAB]
```

If this returns multiple options (e.g. `python`, `python2`, `python 3.7`, `python.exe`), then you have other versions of Python sitting on your system (an old `python2` executable located in `/usr/local/bin/python`, for example). You can see these versions by typing:

```bash
(local)>$ python3 --version     # in bash, csh OR
(local)>$ python.exe --version  # in Cygwin/Windows
```

Do this again for the `pip` command, which could return an old  `pip`, `pip3`, or `pip.exe` in addition to the Anaconda pip distribution. Find and set your `$PATH` environment variable to point to the Anaconda Python *and* pip distributions.

```bash
# with bash:
(local)>$ echo 'export PATH=/Users/username/anaconda3/bin:$PATH' >> ~/.bash_profile
# with csh/tsch:
(local)>$ echo 'setenv PATH $PATH\:/Users/username/anaconda3/bin\:$HOME/bin\:.'  >> ~/.cshrc
```

Confirm the setting with the `which` command:

```bash
(local)>$ which python3         # in bash, csh OR
(local)>$ which python.exe      # in Cygwin/Windows
```

We are looking for a Python executable that looks like it was installed with Anaconda, such as:

```bash
/username/anaconda3/bin/python3 # on MacOS/Linux, OR
/username/anaconda3/python.exe # on Cygwin/Windows
```

If `which` points to either of those locations, you are good to go and you can proceed from here using the shorthand path to your Anaconda Python distribution:

```bash
python3     # Linux/MacOS 
python.exe  # Cygwin/Windows
```

If, however, `which` points to some other location, such as `/usr/local/bin/python`, proceed from here using the FULL path to the Anaconda Python distribution like so:

```bash
/Users/username/anaconda3/bin/python3 # Linux/MacOS
/Users/username/anaconda3/python.exe  # Cygwin/Windows
```
***






### Step 2: Set Up the Virtual Environment:

The virtual environment is created from the terminal with the following syntax:

```bash
(local)>$ python3 -m venv --system-site-packages amesGCM3` # Use FULL PATH to python if needed
```

We can now activate the virtual environment with:

```bash
(local)>$ source amesGCM3/bin/activate      # if you are using **bash**
(local)>$ source amesGCM3/bin/activate.csh  # if you are using **csh/tcsh**
```

> Note that in Cygwin/Windows, the `/bin` directory may be named `/Scripts`.

You may notice that after sourcing `amesGCM3`, your prompt changed to `(amesGCM3)>$`. This confirms that you are **inside** the virtual environment even when you navigate to different directories on your machine.

After sourcing the virtual environment, we can verify that `which python` and `which pip` unambiguously point to `amesGCM3/bin/python3` and `amesGCM3/bin/pip`, respectively. There is therefore no need to reference their full paths for the following instructions.


***




## Installing CAP

Now we can download and install cap in the `amesGCM3` virtual environment. CAP was provided to you in the tarfile `CAP_tarball.zip` that was sent along with these instructions. Download CAP_tarball.zip and leave it in `Downloads/`. Open a terminal window, activate the virtual environment, and untar the file:

```bash
(local)>$ source ~/amesGCM3/bin/activate 
(amesGCM3)>$ 
(amesGCM3)>$ cd ~/Downloads 
(amesGCM3)>$ tar -xf CAP_tarball.zip
(amesGCM3)>$ cd amesgcm-master
(amesGCM3)>$ pip install .
```

That's it! CAP is installed in `amesGCM3`, and you can see the MarsXXXX.py tools in `~/amesGCM3/bin/`:

```bash
(local)>$ ls ~/amesGCM3/bin/
> Activate.ps1     MarsPull.py      activate.csh     f2py             nc4tonc3         pip3
> MarsFiles.py     MarsVars.py      activate.fish    f2py3            ncinfo           pip3.8
> MarsInterp.py    MarsViewer.py    easy_install     f2py3.8          normalizer       python
> MarsPlot.py      activate         easy_install-3.8 nc3tonc4         pip              python3
```

It is now safe to remove both `amesgcm-master/` and the `.zip` archive from your `/Downloads` directory since `pip` installed the pipeline inside your `amesGCM3` virtual environment:

```bash
(local)>$ cd ~/Downloads
(local)>$ rm -r CAP_tarball.zip amesgcm3
```

***



Double check that the paths to the executables are correctly set in your terminal by exiting the virtual environment:

```bash
(amesGCM3)>$ deactivate
```

This completes the one-time installation of CAP in your virtual environment, `amesGCM3`, which now looks like:

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

For your reference, CAP requires the following Python packages. These were installed automatically when you installed CAP:

* matplotlib        # the MatPlotLib plotting library
* netCDF4 Python    # handling netCDF files
* requests          # downloading data from the MCMC Portal

> **Note:** If you prefer using the `conda` package manager for setting up your virtual environment instead of `pip`, you may use the following commands to install CAP.
> 
> First, verify (using `conda info` or `which conda`) that you are using the intented `conda` executable (two or more versions of `conda` might be present if both Python2 and Python3 are installed on your system):
>
>```bash
>(local)>$ conda create -n amesGCM3
>(local)>$ conda activate amesGCM3
>(amesGCM3)>$ conda install pip
>(amesGCM3)>$ pip install git+https://github.com/alex-kling/amesgcm.git
>```
>
> The source code will be installed in:
>
>```bash
>/path/to/anaconda3/envs/amesGCM3/
>```
>
> and the pipeline can then be activated and exited with conda:
>```bash
>(local)>$ conda activate amesGCM3
>(amesGCM3)>$ conda deactivate
>(local)>$ 
>```






***

### Removing CAP

To permanently remove CAP, activate the virtual environment and run the uninstall command:

```bash
(local)>$ source amesGCM3/bin/activate      # bash
(local)>$ source amesGCM3/bin/activate.csh  # csh/tcsh
(amesGCM3)>$ pip uninstall amesgcm
```

You may also delete the `amesGCM3` virtual environment directory at any time. This will uninstall CAP, remove the virtual environment from your machine, and will not affect your main Python distribution.

***




## Testing & Using CAP

Whenever you want to use CAP, simply activate the virtual environment and all of CAP's tools will be accessible from the command line:

```bash
(local)>$ source amesGCM3/bin/activate      # bash
(local)>$ source amesGCM3/bin/activate.csh` # csh/tcsh
```

You can check that the tools are installed properly by typing `Mars` adn then the **TAB** key. No matter where you are on your system, you should see the following pop up:

```bash
(amesGCM3)>$ Mars
> MarsFiles.py   MarsInterp.py  MarsPlot.py    MarsPull.py    MarsVars.py
```

If no executables show up then the paths have not been properly set in the virtual environment. You can either use the full paths to the executable:

```bash
(amesGCM3)>$ ~/amesGCM3/bin/MarsPlot.py
```

Or set up aliases in your `./bashrc` or `.cshrc`:

```bash
# with bash:
(local)>$ echo alias MarsPlot='/Users/username/amesGCM3/bin/MarsPlot.py' >> ~/.bashrc 
(local)>$ source ~/.bashrc
# with csh/tsch
(local)>$ echo alias MarsPlot /username/amesGCM3/bin/MarsPlot >> ~/.cshrc
(local)>$ source ~/.cshrc
```

***








## Tips

### Enable Syntax Highlighting for the Plot Template

The MarsPlot plotting routine needs an input template written in Python, and the template generated by MarsPlot has the `.in` file extension. If you want your text editor to colorize python syntax in the `.in` text file, select "Python" as the language in your text editor (gedit, atom ...). Enabling the correct syntax-highlighting in **vim** requires you add the following lines to `~/.vimrc`:

```bash
syntax on
colorscheme default
au BufReadPost *.in  set syntax=python
```



### Install `ghostscript` to Create Multiple-Page PDFs using MarsPlot

Installing `ghostscript` on your local machine allows CAP to generate a multiple-page PDF file when creating lots of plots. To check whether you already have `ghostscript` on your machine, open a terminal and type:

```bash
(local)>$ gs -version
> GPL Ghostscript 9.54.0 (2021-03-30)
> Copyright (C) 2021 Artifex Software, Inc.  All rights reserved.
``` 

If it is not installed, follow the directions on the `ghostscript` [website](https://www.ghostscript.com/download.html).


***
