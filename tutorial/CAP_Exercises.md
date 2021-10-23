![](./tutorial_images/Tutorial_Banner_Final.png)


<!-- TOC titleSize:2 tabSpaces:2 depthFrom:1 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 skip:0 title:1 charForUnorderedList:* -->
## Table of Contents
* [Practical: Using the Community Analysis Pipeline (CAP)](#practical-using-the-community-analysis-pipeline-cap)
  * [1. Retrieve Data](#1-retrieving-data)
    * [1.1 Download MGCM output](#11-use-marspullpy-to-download-mgcm-output)
  * [2. File Manipulations](#2-file-manipulations)
      * [2.1 `fort.11` to `netCDF` Conversion](#21-convert-the-fort11-files-into-netcdf-files-for-compatibility-with-cap)
      * [2.2 Interpolate `atmos_average` to standard pressure](#22-interpolate-atmos_average-to-standard-pressure-coordinates)
      * [2.3 Add `msf` to `atmos_average_pstd`](#23-add-mass-stream-function-msf-to-atmos_average_pstd)
      * [2.4 Add `rho` and `zfull` to `atmos_average`](#24-add-density-rho-and-mid-point-altitude-zfull-to-atmos_average)
      * [2.5 Interpolate `atmos_average` to standard altitude](#25-interpolate-atmos_average-to-standard-altitude)
      * [2.6 Time-Shift & Pressure-Interpolate the Diurn File](#26-use-marsfiles-to-time-shift-the-diurn-file-then-pressure-interpolate-the-file-with-marsinterp)
      * [2.7 Apply a Low-Pass Filter (`-lpf`) to the `atmos_daily` File](#27-Apply-a-low-pass-filter--lpf-to-the-surface-pressure-ps-and-temperature-ts-in-the-atmos_daily-file)
      * [2.8 Estimate the Magnitude of the Wind Shear](#28-estimate-the-magnitude-of-the-wind-shear-using-cap)
      * [2.9 Determine the Minimum, Mean, and Maximum Near-Surface Temperatures](#29-display-the-values-of-pfull-then-display-the-minimum-mean-and-maximum-near-surface-temperatures-temp-over-the-globe)
  * [Break!](#break)
  * [3. Plotting Routines](#3-plotting-routines)
      * [3.1 Global Map: Surface Albedo and Topography](#31-create-a-global-map-of-surface-albedo-alb-with-topography-zsurf-contoured-on-top)
      * [3.2 Zonal Mean Zonal Wind Cross-Section: RIC](#32-plot-the-zonal-mean-zonal-wind-cross-section-at-ls=270-using-altitude-as-the-vertical-coordinate)
      * [3.3 Zonal Mean Zonal Wind Cross-Section: RAC](#33-add-the-same-plot-for-the-rac-case-to-the-same-page)
      * [3.4 Overplot Temperatures](#34-overplot-temperature-in-solid-contours)
      * [3.5 Four Global Maps on One Page: `lon X lat`](#35-plot-the-following-four-global-maps-lon-x-lat-on-a-new-page)
      * [3.6 Four Global Maps on One Page: `time X lat`](#36-plot-the-following-four-time-x-lat-surface-variables-at-150-e-longitude-on-a-new-page)
      * [3.7 Four Global Maps on One Page: `time X lev`](#37-plot-the-following-four-time-x-lev-variables-at-150-e-longitude-averaged-over-all-latitudes-on-a-new-page)
      * [3.8 Two Cross-Sections on One Page](#38-plot-the-following-two-cross-sections-lat-x-lev-on-the-same-page)
      * [3.9 Zonal Mean Temperatures: RIC and RAC](#39-plot-the-zonal-mean-temperature-at-ls=270°-from-the-atmos_average-file-for-both-the-ric-and-rac-cases,-then-create-a-difference-plot)
      * [3.10 1D Temperature Profiles](#310-generate-two-1d-temperature-profiles-temp-from-the-ric-case)
      * [3.11 Tidal Analysis](#311-plot-the-filtered-and-un-filtered-surface-pressure-over-a-20-sol-period)
<!-- /TOC -->

***

# Practical: Using the Community Analysis Pipeline (CAP)

**Lecture Recap:** CAP is a Python toolkit designed to simplify post-processing and plotting MGCM output. Specifically, CAP consists of five Python executibles with the following functions:

1. `MarsPull.py` Accessing MGCM output
2. `MarsFiles.py` Reducing the files
3. `MarsVars.py` Performing variable operations
4. `MarsInterp.py` Interpolating the vertical grid
5. `MarsPlot.py` Visualizing the MGCM output

When learning to use CAP, it is useful to divide these functions into three categories and explore them in order:

1. [Retrieving Data](#1-retrieving-data) -> `MarsPull.py`
2. [File Manipulations](#2-file-manipulations) -> `MarsFiles.py`, `MarsVars.py`, & `MarsInterp.py`
3. [Plotting Routines](#3-plotting-routines) -> `MarsPlot.py`

We will practice using the functions in CAP in this order. You already have experience using `MarsPull.py` for retrieving data, which was covered at the end of [CAP_Install.md](https://github.com/alex-kling/amesgcm/blob/master/tutorial/CAP_Install.md), and we will build on that knowledge in this tutorial. You may revisit the installation instructions at any time throughout this tutorial.

## Activate CAP

As always with CAP, you must activate the `amesGCM3` virtual environment to access the Python executibles:

```bash
(local)>$ source ~/amesGCM3/bin/activate      # bash
(local)>$ source ~/amesGCM3/bin/activate.csh  # csh/tcsh
```

Remember, you can access the documentation for an executable using the `--help` argument (`-h` for short):

```bash
(amesGCM3)>$ MarsPull.py -h
```

Let's get started by reviewing the data retrieval process covered at the end of [CAP_Install.md](https://github.com/alex-kling/amesgcm/blob/master/tutorial/CAP_Install.md).




***

## 1. Retrieving Data

### 1.1 Use `MarsPull.py` to download MGCM output

`MarsPull` is a utility for accessing MGCM output files hosted on the [MCMC Data portal](https://data.nas.nasa.gov/legacygcm/data_legacygcm.php). During the installation, you were asked to create a `CAP_Tutorial/INERTCLDS/` directory, a `CAP_Tutorial/ACTIVECLDS/` directory, and to use `MarsPull` to download several `fort.11` files into each. You should have already downloaded the necessary `fort.11` files for this tutorial, so this step should be simply a review. If you haven't downloaded the files, however, you can do so now by following these instructions.

Create a `CAP_Tutorial` directory in your preferred location, and then create two subdirectories within it. Name them `INERTCLDS/` and `ACTIVECLDS/`, which represent data from the radiatively intert cloud (RIC) simulation and radiatively active cloud (RAC) simulation, respectively:

```bash
(amesGCM3)>$ cd ~/path/do/preferred/location
(amesGCM3)>$ mkdir CAP_Tutorial
(amesGCM3)>$ cd CAP_Tutorial
(amesGCM3)>$ mkdir ACTIVECLDS INERTCLDS
```

Navigate to `INERTCLDS/` and use `MarsPull` to retrieve files from the RIC simulation. Specify the simulation identifier (`INERTCLDS/`) and the range of Solar Longitudes (Ls=255° through Ls=285°) corresponding to the desired file(s) in the call to `MarsPull`:

```bash
(amesGCM3)>$ MarsPull.py -id INERTCLDS -ls 255 285
```

Then, navigate to the `ACTIVECLDS/` directory and do the same for the RAC case. After this step, you should see the following 5 `fort.11` files in each directory (`INERTCLDS/` and `ACTIVECLDS/`):

```bash
(amesGCM3)>$ ls CAP_Tutorial/INERTCLDS
> fort.11_0719 fort.11_0720 fort.11_0721 fort.11_0722 fort.11_0723
(amesGCM3)>$ ls CAP_Tutorial/ACTIVECLDS
> fort.11_0719 fort.11_0720 fort.11_0721 fort.11_0722 fort.11_0723
```

> If you have any `fort.11` files **other than the ones listed above** in **either** directory, please delete them. It will be easier to follow along during the tutorial if you work with the specific subset of files listed above.




***

## 2. File Manipulations

After retrieving the `fort.11` files from the data portal, we can process the data using CAP. CAP's post-processing capabilities include interpolating and regridding data to different vertical coordinate systems, adding derived variables to the files, and converting between filetypes, just to name a few examples.

The following exercises are designed to demonstrate how CAP can be used for post-processing MGCM output. You should follow along using the files in your `INERTCLDS/` and `ACTIVECLDS/` directories. 

**After post-processing these files, we will use them to make plots with `MarsPlot` so do *not* delete anything!**


Start with the RIC simulation (`INERTCLDS/`) and complete exercises 2.1-2.8 below. We will then provide specific instructions regarding which exercises to repeat with the RAC simulation (`ACTIVECLDS/`).




***

### 2.1 Convert the `fort.11` files into `netCDF` files for compatibility with CAP

First, go to your `INERTCLDS/` directory, and type:

```bash
(amesGCM3)>$ MarsFiles.py fort.11_* -fv3 fixed average daily diurn
```

This creates several `netCDF` files from the `fort.11` files:

```bash
(amesGCM3)>$ ls
> 07180.atmos_average.nc  07190.atmos_average.nc  07200.atmos_average.nc  07210.atmos_average.nc  07220.atmos_average.nc
> 07180.atmos_daily.nc    07190.atmos_daily.nc    07200.atmos_daily.nc    07210.atmos_daily.nc    07220.atmos_daily.nc
> 07180.atmos_diurn.nc    07190.atmos_diurn.nc    07200.atmos_diurn.nc    07210.atmos_diurn.nc    07220.atmos_diurn.nc
> 07180.fixed.nc          07190.fixed.nc          07200.fixed.nc          07210.fixed.nc          07220.fixed.nc
```

The `netCDF` filetypes and a description of their contents are listed below:

| Type                  | Description |
| --------------------- | ----------- |
| `*atmos_fixed.nc`     | static variables that **do not change over time**           |
| `*atmos_average.nc`   | **5-day averages** of MGCM output                           |
| `*atmos_diurn.nc`     | files contain **hourly** MGCM output averaged over 5 days   |
| `*atmos_daily.nc`     | **continuous time series** of the MGCM output               |

> Note that the 5-digit number at the begining of each `netCDF` file corresponds to the sol number at which that file's records begin. These files are pulled from a simulation that was warm-started from a 10 year run. `10 years x ~668 sols/year = 6680 sols`. The earliest date on our files is 07180 (the middle of the year).


For easier post-processing and plotting, we can combine like files along the `time` axis to create one of each filetype:

```bash
(amesGCM3)>$ MarsFiles.py *fixed.nc -c
(amesGCM3)>$ MarsFiles.py *average.nc -c
(amesGCM3)>$ MarsFiles.py *diurn.nc -c
(amesGCM3)>$ MarsFiles.py *daily.nc -c
```

After merging the files, our directory contains just **four** `netCDF` files:

```bash
(amesGCM3)>$ ls
> 07180.atmos_fixed.nc 07180.atmos_average.nc 07180.atmos_diurn.nc 07180.atmos_daily.nc
```





***

### 2.2 Interpolate `atmos_average` to standard pressure coordinates

This step makes use of `MarsInterp`. The documentation for `MarsInterp` can be viewed using the following command:

```bash
(amesGCM3)>$ MarsInterp.py -h
```

The documentation explains that to convert a file to standard pressure coordinates, we must use the `--tshift` (`-t`) command like so:

```bash
(amesGCM3)>$ MarsInterp.py 07180.atmos_average.nc -t pstd
```

This preserves the original file and creates a new, pressure-interpolated file named after the new vertical axis:

```bash
> 07180.atmos_average_pstd.nc
```


***

### 2.3 Add mass stream function (`msf`) to `atmos_average_pstd`

Take a look at the documentation for `MarsVars`:

```bash
(amesGCM3)>$ MarsVars.py -h # display documentation
```

Adding or removing variables from files is done using `-add` and `-rem` in the call to `MarsVars`. We can add `msf` to the pressure-interpolated average file as follows:

```bash
(amesGCM3)>$ MarsVars.py 07180.atmos_average_pstd.nc -add msf
```

`MarsVars` adds variables to the specified file (as opposed to creating a new file containing the variables). In this case, `msf` was derived from the meridional wind (`vcomp`). Note that `msf` could not be added before pressure-interpolating the file because it must be derived on vertical pressure coordinates.




***

### 2.4 Add density (`rho`) and mid-point altitude (`zfull`) to `atmos_average`


Again, add variables to files using `-add`:

```bash
(amesGCM3)>$ MarsVars.py 07180.atmos_average.nc -add rho zfull
```

Density (`rho`) was derived from the pressure and temperature variables output by the model, and mid-point altitude (`zfull`) was obtained via hydrostatic integration.

Note that we've added `rho` to the **non-interpolated** file. This is because some derived variables (such as `rho`) are computed on the native model grid in CAP, and will not be properly derived on any other vertical grid. Other variables (like `msf` from the previous step) are computed on a pressure grid, in which case we have to interpolate the file *before* adding `msf` to it.






***
### 2.5 Interpolate `atmos_average` to standard altitude

We added `rho` and `zfull` to `atmos_average` in the last step because these variables cannot be added to an interpolated file. We can now interpolate the file to standard altitude coordinates so that we have `rho` and `zfull` in an interpolated file:

```bash
(amesGCM3)>$ MarsInterp.py 07180.atmos_average.nc -t zstd   # standard altitude
```

Our directory now contains three `atmos_average` files:

```bash
> 07180.atmos_average.nc 07180.atmos_average_pstd.nc 07180.atmos_average_zstd.nc
```

The original file (`07180.atmos_average.nc`) and the altitude-interpolated file (`07180.atmos_average_zstd.nc`) contain `rho` and `zfull`. The pressure-interpolated file (`07180.atmos_average_pstd.nc`) contains `msf` but does not contain `rho` or `zfull`, because it was created before we added those variables. You can confirm that these are the variables in each file using the `--inspect` (`-i`) function from `MarsPlot`:

```bash
(amesGCM3)>$ MarsPlot.py -i 07180.atmos_average.nc          # the original file + rho, zfull
(amesGCM3)>$ MarsPlot.py -i 07180.atmos_average_zstd.nc     # the altitude interpolated file + rho, zfull
(amesGCM3)>$ MarsPlot.py -i 07180.atmos_average_pstd.nc     # the pressure interpolated file + msf
```









***

### 2.6 Use `MarsFiles` to time-shift the `diurn` file, then pressure-interpolate the file with `MarsInterp`

The variables in `07180.atmos_diurn.nc` are organized by time-of-day assuming **universal time** beginning at the Martian prime meridian. You can time-shift the fields to **uniform local time** using `MarsFiles`. This function is useful for comparing MGCM output to observations from satellites in fixed local time orbit, for example.

For this exercise, we will only time-shift the surface pressure (`ps`), surface temperature (`ts`), and atmospheric temperature (`temp`) variables in order to minimize file size and processing time. Time-shifting with `MarsFiles` and only including `ts`, `ps`, and `temp` looks like this:

```bash
(amesGCM3)>$ MarsFiles.py 07180.atmos_diurn.nc -t --include ts ps temp
```

Time-shifting can only be done on the `diurn` files since these contain hourly output. The time-shifting function creates a new file ending in `_T.nc` to differentiate it from the original file.

After time-shifting, pressure-interpolate the file using `MarsInterp` like before:

```bash
(amesGCM3)>$ MarsInterp.py 07180.atmos_diurn_T.nc -t pstd
```

This should take just over a minute. Interpolating large files (such as `daily` or `diurn` files) with CAP can take a long time because the code is written in Python. That's why we are only including three variables (`ps`, `ts`, and `temp`) in this particular demonstration. Note that we time-shifted the file before interpolating it, but that time-shifting and interpolating can be done in any order.

After time-shifting and interpolating, our directory contains three `diurn` files:

```bash
> 07180.atmos_diurn.nc 07180.atmos_diurn_T.nc 07180.atmos_diurn_T_pstd.nc
```

> **Note:** We will *not* do this here, but you can create custom vertical grids that you want CAP to interpolate to. See the documentation for `MarsInterp` for more information.




***

### 2.7 Apply a low-pass filter (`-lpf`) to the surface pressure (`ps`) and temperature (`ts`) in the `atmos_daily` file

Let's use a 10-sol cut-off frequency (`sol_max` > 10) so we can isolate synoptic-scale features. Including only `ps` and `ts` means only those two variables will be time-filtered:

```bash
(amesGCM3)>$ MarsFiles.py 07180.atmos_daily.nc -lpf 10 -include ps ts
```

This function created a new file, `atmos_daily_lpf`, which contains only `ps`, `ts`, and their dimensions.





### 2.8 Estimate the magnitude of the wind shear using CAP

Do this by adding dU/dZ and dV/dZ to `07180.atmos_average_zstd.nc`. You already know that `MarsVars` is used when adding variables, however, deriving the zonal (`ucomp`) and meridional (`vcomp`) wind shear is **not** done using the `-add` function like we've practiced. Instead, we use `-zdiff` to perform a vertical differentiation on `ucomp` and `vcomp`. `MarsVars` also has a `-col` function, which performs a column integration on specified variables. 

Apply vertical differentiation to add dU/dZ and dV/dZ to the file as follows:

```bash
(amesGCM3)>$ MarsVars.py 07180.atmos_average_zstd.nc -zdiff ucomp vcomp
```

Use `--inspect` (`-i`) to see the contents of `07180.atmos_average_zstd.nc` and find the names of the derived variables dU/dZ and dV/dZ:

```bash
(amesGCM3)>$ MarsPlot.py -i 07180.atmos_average_zstd.nc
> ===================DIMENSIONS==========================
> ['lat', 'lon', 'phalf', 'time', 'zstd']
> (etc)
> ====================CONTENT==========================
> (etc)
> d_dz_ucomp     : ('time', 'zstd', 'lat', 'lon')= (10, 45, 36, 60), vertical gradient of zonal wind  [m/s/m]]
> d_dz_vcomp     : ('time', 'zstd', 'lat', 'lon')= (10, 45, 36, 60), vertical gradient of meridional wind  [m/m]]
> (etc)
>
> Ls ranging from 255.42 to 284.19: 45.00 days
>                (MY 01)   (MY 01)
> =====================================================
```

> **Note: The `--inspect` function works on any netCDF file, not just the ones created with CAP!**








### 2.9 Display the values of `pfull`, then display the minimum, mean, and maximum near-surface temperatures `temp` over the globe

We can display values in an array using `--dump` (analog of the NCL command `ncdump`) and `MarsPlot -i`. For example, the content of the reference pressure (`pfull`) variable in `07180.atmos_average.nc` is viewed by typing:

```bash
(amesGCM3)>$ MarsPlot.py -i 07180.atmos_average.nc -dump temp
> temp=
> [8.7662227e-02 2.5499690e-01 5.4266089e-01 1.0518962e+00 1.9545468e+00
> 3.5580616e+00 6.2466631e+00 1.0509957e+01 1.7400265e+01 2.8756382e+01
> 4.7480076e+01 7.8348366e+01 1.2924281e+02 2.0770235e+02 3.0938846e+02
> 4.1609518e+02 5.1308148e+02 5.9254102e+02 6.4705731e+02 6.7754218e+02
> 6.9152936e+02 6.9731799e+02 6.9994830e+02 7.0082477e+02]
> ______________________________________________________________________
```

We can index specific values in the array using quotes and square brackets with the variable call (`'var[ ]'`). For example, the reference pressure of the first layer above the surface is indexed as follows:

```bash
(amesGCM3)>$ MarsPlot.py -i 07180.atmos_average.nc -dump 'pfull[-1]'
> pfull[-1]=
> 700.8247680664062
> ______________________________________________________________________
```

> **Note:** `-1` refers to the last array element (Python syntax)

Using `-stat` with `MarsPlot -i` calculates and displays the minimum, mean, and maximum values of a variable, which is better suited for visualizing statistics over large arrays or in data slices. Given the dimensions of the `temp` variable, `[time,pfull,lat,lon]`, we display the stats for the **near-surface air temperature over all timesteps and all locations**:

```bash
(amesGCM3)>$ MarsPlot.py -i 07180.atmos_average.nc -stat 'temp[:,-1,:,:]'
> __________________________________________________________________________
>            VAR            |      MIN      |      MEAN     |      MAX      |
> __________________________|_______________|_______________|_______________|
>             temp[:,-1,:,:]|        149.016|        202.508|         251.05|
> __________________________|_______________|_______________|_______________|
```

> **Note:** quotes `''` are necessary when browsing dimensions.

and that's it for post-processing the data in the RIC simulation! Before we move on to plotting, we need to repeat some of these steps for the RAC simulation. Feel free to repeat all of Steps 2.1-2.8 for the RAC case if you like, but **you are only required to repeat Steps 2.1, 2.2, and 2.4** for this tutorial:

* [2.1 `fort.11` to `netCDF` conversion](#21-convert-the-fort11-files-into-netcdf-files-for-compatibility-with-cap)
* [2.2 Interpolate to standard pressure](#22-interpolate-atmosaverage-to-standard-pressure-coordinates)
* [2.4 Add `rho` and `zfull` to `atmos_average`; Interpolate to standard altitude](#24-add-density-rho-and-mid-point-altitude-zfull-to-atmosaverage-then-interpolate-the-file-to-standard-altitude-zstd)


### Navigate to the `ACTIVECLDS/` directory and complete Steps 2.1, 2.2, and 2.4 before continuing







***

# Break!
Once you've completed Step 2 for **both** simulations, you are welcome to take a 15 minute break from the tutorial. You can use this time to catch up if you haven't completed Steps 1 and/or 2 already, but we highly encourage you to step away from your machine for these 15 minutes.








***

## 3. Plotting Routines

The last part of this tutorial introduces the plotting capabilities of CAP. CAP's plotting routine is `MarsPlot`, and it can create several kinds of plots:

|Type of plot      |    MarsPlot designation|
|----------------------|-------------------|
| Longitude v Latitude | Plot 2D lon X lat|
| Longitude v Time     |Plot 2D lon X time|
| Longitude v Level    |Plot 2D lon X lev|
| Latitude v Level     |Plot 2D lat X lev|
| Time v Latitude      |Plot 2D time X lat|
| Time v level         |Plot 2D time X lev |
| Any 1-dimensional line plot |Plot 1D|

`MarsPlot` is highly customizable. Plots can be displayed in PDF or image format and in landscape or portrait mode. Multiple plots can be drawn on a single page, overplotting is supported, and you can choose your axes dimensions, colormap, map projection type, contour levels, and more.

Plotting with CAP requires passing a template to `MarsPlot`. A blank template is created in the directory in which the following command is executed, so navigate to the `INERTCLDS/` directory before typing:

```bash
(amesGCM3)>$ MarsPlot.py -template
```

The blank template is called `Custom.in`. View `Custom.in` by opening it in your preferred text editor. For example:

```bash
(amesGCM3)>$ vim Custom.in
```

Scroll down until you see the first two templates shown in the image below:

![](./tutorial_images/Custom_Templates.png)

By default, `Custom.in` is set to create the two plots shown above: a global topographical map and a zonal mean wind cross-section. The plot type is indicated at the top of each template. In this case, the plot types are:

``` python
> <<<<<<<<<<<<<<| Plot 2D lon X lat = True |>>>>>>>>>>>>>
> <<<<<<<<<<<<<<| Plot 2D lat X lev = True |>>>>>>>>>>>>>
```

When the plot type is set to `True` as it is here, then `MarsPlot` will draw that plot. If `False`, `MarsPlot` will skip that plot. The variable to be plotted is `Main Variable`, which requires the variable name and the file containing it as input:

```python
 Main Variable  = fixed.zsurf # topography from the 01780.fixed.nc file
```

Without making any changes to `Custom.in`, close the file and pass it back to `MarsPlot` using the following command:

```bash
(amesGCM3)>$ MarsPlot.py Custom.in
```

This creates `Diagnostics.pdf`, a single-page PDF displaying the two plots we just discussed: global topography and zonal mean wind. Open the PDF to see the plots.

You can rename `Custom.in` and still pass it to `MarsPlot` successfully. If the template is named anything other than `Custom.in`, `MarsPlot` will produce a PDF named after the template, i.e. `myplots.in` creates `myplots.pdf`. For example:

```bash
(amesGCM3)>$ mv Custom.in myplots.in
(amesGCM3)>$ MarsPlot.py myplots.in
> Reading myplots.in
> [----------]  0 % (2D_lon_lat :fixed.zsurf)
> [#####-----] 50 % (2D_lat_lev :atmos_average.ucomp, Ls= (MY 1) 284.19, zonal avg)
> [##########]100 % (Done)
> Merging figures...
> "/username/CAP_Tutorial/INERTCLDS/myplots.pdf" was generated
```

Those are the basics of plotting with CAP. We'll create several plots in exercises 3.1-3.11 below. Begin by deleting `myplots.in` and `myplots.pdf` (if you have them), and then create a new `Custom.in` template:

```bash
(amesGCM3)>$ rm myplots.in myplots.pdf
(amesGCM3)>$ MarsPlot.py -template
```

Make sure you're in the `INERTCLDS/` directory before continuing.



***

### 3.1 Create a global map of surface albedo (`alb`) with topography (`zsurf`) contoured on top

Open `Custom.in` in your preferred text editor and make the following changes:

- Set the second default template (`Plot 2D lat X lev`) to `False` so that `MarsPlot` does not draw it (we will use it later)
- On the first template (`Plot 2D lon X lat`), set `Main Variable` to albedo (`alb`, located in the `fixed` file). This will be plotted as shaded contours
- Also set `2nd Variable` to topography (`zsurf`, located in the `fixed` file). This will be plotted as solid contours
- Set the `Title` to reflect the variable(s) being plotted

Here is what your template should look like:

```python
> <<<<<<<<<<<<<<| Plot 2D lon X lat = True |>>>>>>>>>>>>>
> Title          = 3.1: Albedo w/Topography Overplotted
> Main Variable  = fixed.alb
> Cmin, Cmax     = None
> Ls 0-360       = None
> Level [Pa/m]   = None
> 2nd Variable   = fixed.zsurf
> Contours Var 2 = None
> Axis Options  : lon = [None,None] | lat = [None,None] | cmap = binary | scale = lin | proj = cart
```

Save `Custom.in` (but don't close it!) and go back to the terminal. Pass `Custom.in` back to `MarsPlot`:

```bash
(amesGCM3)>$ MarsPlot.py Custom.in
```

Open `Diagnostics.pdf` and check to make sure it contains a global map of surface albedo with topography contoured overtop.

> Depending on the settings for your specific PDF viewer, you may have to close and reopen the file to view it.



***

### 3.2 Plot the zonal mean zonal wind cross-section at Ls=270° using altitude as the vertical coordinate

In the same template, which should still be open in your text editor:
- Set the `2D lat X lev` template to `True`
- Change `Main Variable` to point to `ucomp` stored in the `atmos_average_zstd` file
- Edit the `Title` accordingly

Save `Custom.in`, and pass it to `MarsPlot`. Again, view `Diagnostics.pdf` to see your plots!




***

### 3.3 Add the same plot for the RAC case to the same page

> **Tip:** Copy and paste the `lat x lev` plot you made in 3.2 so that you have two identical templates.

First, point `MarsPlot` to the `ACTIVECLDS/` directory. Do this by editing the `<<<<<<< Simulations <<<<<<<` section so that `2>` points to `/ACTIVECLDS` like so:

```python
> <<<<<<<<<<<<<<<<<<<<<< Simulations >>>>>>>>>>>>>>>>>>>>>
> ref> None
> 2> ../ACTIVECLDS
```

Then, edit `Main Variable` in the copy/pasted template so that the `atmos_average` file in `ACTIVECLDS/` is sourced:

```python
> Main Variable  = atmos_average@2.ucomp
```

> **Tip:** Make use of `HOLD ON` and `HOLD OFF` for these.

Save `Custom.in` and pass it to `MarsPlot`. View `Diagnostics.pdf` to see the results.




***

### 3.4 Overplot temperature in solid contours

Add `temp` as the second variable on the plots you created in 3.2 and 3.3:

```python
> <<<<<<<<<<<<<<| Plot 2D lat X lev = True |>>>>>>>>>>>>>
> > 2nd Variable     = atmos_average_zstd.temp
> (etc)
> 
> <<<<<<<<<<<<<<| Plot 2D lat X lev = True |>>>>>>>>>>>>>
> > 2nd Variable     = atmos_average_zstd@2.temp
> (etc)
```

Save `Custom.in` and pass it to `MarsPlot`. View `Diagnostics.pdf` to see the results.




***

### 3.5 Plot the following four global maps (`lon X lat`) on a new page

- Surface CO2 Ice Content (`snow`) north of 50 latitude
- Surface Temperature (`ts`). Set the colorscale (`Cmin, Cmax`) to range from 150 K to 300 K
- Surface Wind Speed (`(u^2 + v^2)/2`) (this requires the use of square brackets **and** two variables)
- Diabatic Heating Rate (`dheat`) at 50 Pa (index dimension `lev=50`)

> **Tip:** Use `HOLD ON` and `HOLD OFF` again. You can use this syntax multiple times in the same template.

Source all of the variables from the `atmos_daily` file in `INERTCLDS/`. Plot the results at Ls=270°. The general format will be:

```python
> HOLD ON
> 
> <<<<<<| Plot 2D lon X lat = True |>>>>>>
> Title    = Surface CO2 Ice (g/m2)
> (etc)
> 
> <<<<<<| Plot 2D lon X lat = True |>>>>>>
> Title    = Surface Temperature (K)
> (etc)
> 
> <<<<<<| Plot 2D lon X lat = True |>>>>>>
> Title    = Surface Wind Speed (m/s)
> (etc)
> 
> <<<<<<| Plot 2D lon X lat = True |>>>>>>
> Title    = Diabatic Heating Rate (K/sol)
> (etc)
> 
> HOLD OFF
```

> **Note: convert kg -> g using square brackets for variable operations:**
>```python
>Main Variable  = [atmos_daily.snow]*1000
>```
> **Similarly, you can add two variables together like so:**
>```python
>Main Variable  = ([atmos_daily.ucomp]**2+[atmos_daily.vcomp]**2)**0.5
>```

Name the plots accordingly. Save `Custom.in` and pass it to `MarsPlot`. You should see four plots on one page in `Diagnostics.pdf`





***

### 3.6 Plot the following four `time X lat` surface variables at 150° E longitude on a new page

- Surface temperature (`ts`). Set the colormap to `nipy_spectral`.
- Surface dust mass (`dst_mass_sfc`). Set the colormap to `jet`.
- Surface water ice mass (`ice_mass_sfc`) Set the colormap to `Spectral_r`.
- Surface water vapor mass (`vap_mass_sfc`) Set the colormap to `Spectral_r`.

Source all of the variables from the `atmos_average_pstd` file in `INERTCLDS/`. The general format will be:

```python
> HOLD ON
> 
> <<<<<<| Plot 2D time X lat = True |>>>>>>
> Title    = 150 E Surface Temperature (K)
> (etc)
> 
> <<<<<<| Plot 2D time X lat = True |>>>>>>
> Title    = 150 E Surface Dust Mass (kg/m2)
> (etc)
> 
> <<<<<<| Plot 2D time X lat = True |>>>>>>
> Title    = 150 E Surface Water Ice Mass (kg/m2)
> (etc)
> 
> <<<<<<| Plot 2D time X lat = True |>>>>>>
> Title    = 150 E Surface Water Vapor Mass (kg/m2)
> (etc)
> 
> HOLD OFF
```

> **Note: set the colormap in `Axis Options`:**
>```python
>Axis Options  : sols = [None,None] | lat = [None,None] | cmap = nipy_spectral |scale = lin
>```

Name the plots accordingly. Save `Custom.in` and pass it to `MarsPlot`.






***

### 3.7 Plot the following four `time X lev` variables at 150° E longitude, averaged over all latitudes, on a new page

- Meridional mean temperature (`temp`). Set the colormap to `nipy_spectral`.
- Meridional mean dust mass (`dst_mass`). Set the colormap to `jet`.
- Meridional mean water ice mass (`ice_mass`) Set the colormap to `Spectral_r`.
- Meridional mean water vapor mass (`vap_mass`) Set the colormap to `Spectral_r`.

Source all of the variables from the `atmos_average_pstd` file in `INERTCLDS/`. The general format will be:

```python
> HOLD ON
> 
> <<<<<<| Plot 2D time X lev = True |>>>>>>
> Title    = 150 E Temperature (K)
> (etc)
> 
> <<<<<<| Plot 2D time X lev = True |>>>>>>
> Title    = 150 E Dust Mass (kg/kg)
> (etc)
> 
> <<<<<<| Plot 2D time X lev = True |>>>>>>
> Title    = 150 E Water Ice Mass (kg/kg)
> (etc)
> 
> <<<<<<| Plot 2D time X lev = True |>>>>>>
> Title    = 150 E Water Vapor Mass (kg/kg)
> (etc)
> 
> HOLD OFF
```

Name the plots accordingly. Save `Custom.in` and pass it to `MarsPlot`.









***

### 3.8 Plot the following two cross-sections (`lat X lev`) on the same page

- Mass Streamfunction (`msf`) at Ls=270°
  - Change the colormap from `jet` to `bwr`
  - Force symmetrical contouring by setting the colorbar's minimum and maximum values to -50 and 50
  - Adjust the `y axis` limits to 1,000 Pa and 1 Pa
  - Overplot solid contours for `msf=-10` and `msf=10` *(Hint: set both `Main Variable` and `2nd Variable` to `msf`)*
- Zonal mean temperature (`temp`) at Ls=270° from the same (pressure-interpolated) file
  - Overplot the zonal wind (`ucomp`)



Don't forget to use `HOLD ON` and `HOLD OFF` and to name your plots accordingly. Save `Custom.in` and pass it to `MarsPlot`.




***

### 3.9 Plot the zonal mean temperature at Ls=270° from the `atmos_average` file for both the RIC and RAC cases, then create a difference plot

Make use of `HOLD ON` and `HOLD OFF` again here. Copy and paste a `lat x lev` template three times. For the difference plot, you'll need to use `@N` to point to the `ACTIVECLDS/` directory and use square brackets to subtract one variable from the other:

```python
> Main Variable  = [atmos_average_pstd.temp]-[atmos_average_pstd@2.temp]
```

- Set the colormap to `RdBu` for the difference plot
- For all three plots, set the vertical range to 1,000-1 Pa
- Set proper titles

Save `Custom.in` and pass it to `MarsPlot`.




***

### 3.10 Generate two 1D temperature profiles (`temp`) from the RIC case 

Both thermal profiles are at `50°N, 150°E` and Ls=270°, one at 3 AM and the other at 3 PM.

There should be two lines on one plot: the thermal profile at 3 AM and the thermal profile at 3 PM. CAP can overplot 1D data on the same graph by concatenating two 1D templates together with `ADD LINE`:

```python
> <<<<<<| Plot 1D = True |>>>>>>
> Main Variable    = var1
> (etc)
> 
> ADD LINE
> 
> <<<<<<| Plot 1D = True |>>>>>>
> Main Variable    = var2
> (etc)
```

#### *Note: You do not use `HOLD ON` and `HOLD OFF` to overplot 1D plots. `HOLD` is always for drawing separate plots on the same page*

Call `temp` from the `diurn_T_pstd` file, which is the time-shifted and pressure-interpolated version of the hourly file. 3 AM is index=`3`, 3 PM is index=`15`, so `Main Variable` will be set as:

```python
> Main Variable    = diurn_T_pstd.temp{tod=3}
```

and

```python
> Main Variable    = diurn_T_pstd.temp{tod=15}
```

You will have to specify `Level [Pa/m]` as the `y axis`:

```python
> Level [Pa/m]   = AXIS
```

Save `Custom.in` and pass it to `MarsPlot`. View `Diagnostics.pdf`.






***

### 3.11 Plot the filtered and unfiltered surface pressure over a 20 sol period

Here we're asked to compare surface pressure `ps` from the orginal file (`atmos_daily`) to the surface pressure that we time-filtered using `MarsFiles` in exercise 2.6 (`atmos_daily_lpf`). Some hints:

- Both are 1D plots so use `ADD LINE` to plot on the same axes
- Use `ps` from `07180.atmos_daily.nc` and `07180.atmos_daily_lpf.nc`
- Since we need to choose a location, set `Latitude = 50` and `Lon +/-180 = 150`
- Under `Axis Options`, set the `x axis` range (`time`) to (`Ls = [260, 280]`)
- Under `Axis Options`, set the `y axis` range (`ps`) to (`var = [850, 1000]`; Pa)

Save `Custom.in` and pass it to `MarsPlot`. View `Diagnostics.pdf`.





***



## That's a Wrap!

This concludes the practical exercise portion of the CAP tutorial. Please keep these exercises as a reference for the future!

*This document was completed in October 2021. Written by Alex Kling, Courtney Batterson, and Victoria Hartwick*

Please submit feedback to Alex Kling: alexandre.m.kling@nasa.gov




***