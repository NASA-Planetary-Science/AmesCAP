![](./tutorial_images/Tutorial_Banner_Final.png)


<!-- TOC titleSize:2 tabSpaces:2 depthFrom:1 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 skip:0 title:1 charForUnorderedList:* -->
## Table of Contents
* [Practical: Using the Community Analysis Pipeline (CAP)](#practical-using-the-community-analysis-pipeline-cap)
  * [Activate CAP](#activate-cap)
  * [1. Retrieving Data](#1-retrieving-data)
    * [Download MGCM output](#use-marspullpy-to-download-mgcm-output)
  * [2. File Manipulations](#2-file-manipulations)
      * [2.1 `fort.11` to `netCDF` Conversion](#21-convert-the-fort11-files-into-netcdf-files-for-compatibility-with-cap)
      * [2.2 Interpolate to standard pressure](#22-interpolate-atmosaverage-to-standard-pressure-coordinates)
      * [2.3 Add `rho` and `zfull` to `atmos_average`; Interpolate to Standard Altitude](#23-add-density-rho-and-mid-point-altitude-zfull-to-atmosaverage-then-interpolate-the-file-to-standard-altitude-zstd)
      * [2.4 Add `msf` to `atmos_average_pstd`](#24-add-mass-stream-function-msf-to-atmosaveragepstd)
      * [2.5 Time-Shift & Pressure-Interpolate the Diurn File](#25-use-marsfiles-to-time-shift-the-diurn-file-then-pressure-interpolate-the-file)
      * [2.6 Apply a Low-Pass Filter (`-lpf`) to the `atmos_daily` File](#26-apply-a-low-pass-filter--lpf-to-the-surface-pressure-ps-and-temperature-ts-in-the-atmosdaily-file)
      * [2.7 Estimate the Magnitude of the Wind Shear](#27-estimate-the-magnitude-of-the-wind-shear-using-cap)
      * [2.8 Determine the Minimum, Mean, and Maximum Near-Surface Temperature](#28-display-the-minimum-mean-and-maximum-near-surface-temperature)
* [Break!](#break)
  * [3. Plotting Routines](#3-plotting-routines)
      * [3.1 Global Map: Surface Albedo and Topography](#31-create-a-global-map-of-surface-albedo-alb-with-topography-zsurf-contoured-on-top)
      * [3.2 Zonal Mean Zonal Wind Cross-Section: RIC](#32-plot-the-zonal-mean-zonal-wind-cross-section-at-ls=270°-using-altitude-as-the-vertical-coordinate)
      * [3.3 Zonal Mean Zonal Wind Cross-Section: RAC](#33-create-the-same-plot-for-the-rac-case-place-these-plots-on-the-same-page)
      * [3.4 Overplot Temperatures](#34-overplot-temperature-in-solid-contours)
      * [3.5 Four Global Maps on One Page](#35-plot-the-following-four-global-maps-lon-x-lat-on-a-new-page)
      * [3.6 Two Cross-Sections on One Page](#36-plot-the-following-two-cross-sections-lat-x-lev-on-the-same-page)
      * [3.7 Zonal Mean Temperatures: RIC and RAC](#37-plot-the-zonal-mean-temperature-at-ls=270-from-the-atmosaverage-file-for-both-the-ric-and-rac-cases-also-create-a-difference-plot-for-them)
      * [3.8 1D Temperature Profiles](#38-generate-two-1d-temperature-profiles-temp-from-the-ric-case-both-at-50°n-150°e-and-ls=270-at-3-am-and-3-pm)
      * [3.9 Tidal Analysis](#39-plot-the-filtered-and-un-filtered-surface-pressure-over-a-20-sol-period)
<!-- /TOC -->

***

# Practical: Using the Community Analysis Pipeline (CAP)

**Recap:** CAP is a Python toolkit designed to simplify post-processing and plotting MGCM output. Specifically, CAP consists of five Python executibles indended to perform the following functions:

1. `MarsPull.py` Accessing MGCM output
2. `MarsFiles.py` Reducing the files
3. `MarsVars.py` Performing variable operations
4. `MarsInterp.py` Interpolating the vertical grid
5. `MarsPlot.py` Visualizing the MGCM output

When learning to use CAP, it is useful to divide its functions into three categories and explore them in order:

1. [Retrieving Data](#1-retrieving-data)
2. [File Manipulations](#2-file-manipulations)
3. [Plotting Routines](#3-plotting-routines)

We will practice using CAP for all three parts. You already have experience using CAP for Retrieving Data, which was covered at the end of the CAP installation instructions (the install asked you to use `MarsPull` to retrieve several `fort.11` files before the tutorial). Here, you will have a chance to practice using all five Python routines in CAP.

## Activate CAP

As always with CAP, you must activate the `amesGCM3` virtual environment to access the executibles (you can revisit the [installation instructions](https://github.com/alex-kling/amesgcm/blob/master/tutorial/CAP_install.md) as a refresher).

```bash
(local)>$ source ~/amesGCM3/bin/activate      # bash
(local)>$ source ~/amesGCM3/bin/activate.csh  # csh/tcsh
```

As a reminder, each Mars executable has a `--help` argument (`-h` for short) that can show you information about an executible, for example:

```bash
(amesGCM3)>$ MarsPull.py -h
```






***

## 1. Retrieving Data

### Use `MarsPull.py` to download MGCM output

`MarsPull` is a utility for accessing MGCM output files hosted on the [MCMC Data portal](https://data.nas.nasa.gov/legacygcm/data_legacygcm.php). During the installation, you were asked to use `MarsPull` to download several `fort.11` files into your `INERTCLDS/` and `ACTIVECLDS/` directories. You should have already downloaded the necessary `fort.11` files for this tutorial, so this step is likely a review for you. If you haven't downloaded the files, you can do so now by following these instructions:

Create a `CAP_Tutorial` directory in your preferred location, and then create two subdirectories, `INERTCLDS/` and `ACTIVECLDS/`:

```bash
(amesGCM3)>$ cd ~/path/do/preferred/location
(amesGCM3)>$ mkdir CAP_Tutorial
(amesGCM3)>$ cd CAP_Tutorial
(amesGCM3)>$ mkdir ACTIVECLDS INERTCLDS
```

Navigate to `INERTCLDS/` and use `MarsPull` to retrieve the files. Specify the simulation identifier (`INERTCLDS/`) and the range of Solar Longitudes (Ls=255 through Ls=285) corresponding to the desired file(s):

```bash
(amesGCM3)>$ MarsPull.py -id INERTCLDS -ls 255 285
```

Then do the same for the `ACTIVECLDS/` case.

There should now be 5 `fort.11` files in each directory (`INERTCLDS/` and `ACTIVECLDS/`):

```bash
> fort.11_0719 fort.11_0720 fort.11_0721 fort.11_0722 fort.11_0723
```

> If you have any `fort.11` files *other than the ones listed above* in *either* directory, please delete them. It will be make it easier to follow the tutorial if you work with the specific subset of files listed above.




***

## 2. File Manipulations

After retrieving output from the data portal (or recieving output from a simulation), you can process the data to create the files you need for your analysis. Post-processing includes interpolating and regridding data to different vertical coordinate systems, adding derived variables to the files, and converting between filetypes, just to name a few examples.

The following exercises are designed to demonstrate how CAP can be used for post-processing MGCM output. You should follow along in the `INERTCLDS/` and `ACTIVECLDS/` directories you created containing the `fort.11` files you downloaded. 

**After post-processing these files, we will use them to make plots with MarsPlot. Don't delete anything!**


Start with the radiatively inert cloud simulation (RIC; `INERTCLDS/`) and complete exercises 2.1-2.8 below. Then we will give you specific instructions regarding which exercises to repeat for the radiatively active cloud simulation (RAC; `ACTIVECLDS/`).




***

#### 2.1 Convert the `fort.11` files into `netCDF` files for compatibility with CAP

To do this, go to your `INERTCLDS/` directory, and type:

```bash
(amesGCM3)>$ MarsFiles.py fort.11_* -fv3 fixed average daily diurn
```

This created several `netCDF` files:

```bash
(amesGCM3)>$ ls
> 07180.atmos_average.nc  07190.atmos_average.nc  07200.atmos_average.nc  07210.atmos_average.nc  07220.atmos_average.nc
> 07180.atmos_daily.nc    07190.atmos_daily.nc    07200.atmos_daily.nc    07210.atmos_daily.nc    07220.atmos_daily.nc
> 07180.atmos_diurn.nc    07190.atmos_diurn.nc    07200.atmos_diurn.nc    07210.atmos_diurn.nc    07220.atmos_diurn.nc
> 07180.fixed.nc          07190.fixed.nc          07200.fixed.nc          07210.fixed.nc          07220.fixed.nc
```
> Note that the 5-digit number at the begining of each netCDF file corresponds to the sol number at which that file's records begin. These files are from 500 days into a simulation that was warm-started from a 10 year run as the file dates indicate: 10 years x ~668 sols/year = 6680 sols + 500 sols = 7180 sols. Hence the first date is 07180.


The `netCDF` filetypes and a description of their contents are below.

| Type                  | Description |
| --------------------- | ----------- |
| `*atmos_fixed.nc`     | static variables that **do not change over time**           |
| `*atmos_average.nc`   | **5-day averages** of MGCM output                           |
| `*atmos_diurn.nc`     | files contain **hourly** MGCM output averaged over 5 days   |
| `*atmos_daily.nc`     | **continuous time series** of the MGCM output               |

For easier post-processing and plotting, we can combine like files along the time axis to create one large file of each filetype:

```bash
(amesGCM3)>$ MarsFiles.py *fixed.nc -c
(amesGCM3)>$ MarsFiles.py *average.nc -c
(amesGCM3)>$ MarsFiles.py *diurn.nc -c
(amesGCM3)>$ MarsFiles.py *daily.nc -c
```

This merge leaves us with the following four files:

```bash
> 07180.atmos_fixed.nc 07180.atmos_average.nc 07180.atmos_diurn.nc 07180.atmos_daily.nc
```





***

#### 2.2 Interpolate `atmos_average` to standard pressure coordinates

This requires using `MarsInterp`. As a reminder, you can display documentation for MarsInterp using:

```bash
(amesGCM3)>$ MarsInterp.py -h
```

Convert to standard pressure coordinates by entering the following:
```bash
(amesGCM3)>$ MarsInterp.py 07180.atmos_average.nc -t pstd
```

which creates:

```bash
> 07180.atmos_average_pstd.nc
```



***

#### 2.3 Add density (`rho`) and mid-point altitude (`zfull`) to `atmos_average`, then interpolate the file to standard altitude (`zstd`)

Adding or removing variables from files can be done with `MarsVars`:

```bash
(amesGCM3)>$ MarsVars.py -h # display documentation
(amesGCM3)>$ MarsVars.py 07180.atmos_average.nc -add rho zfull
```

This adds the new variables to the original file. Density (`rho`) was derived from the pressure and temperature variables already present in the file, and  mid-point altitude (`zfull`) was obtained via hydrostatic integration.

> **NOTE: if you want `rho` in an interpolated file, you need to add it before performing the interpolation because. In this case, we want `rho` in an altitude-interpolated file so we've added `rho` to the original file (`atmos_average.nc`) and we perform the interpolation after.**

```bash
(amesGCM3)>$ MarsInterp.py 07180.atmos_average.nc -t zstd   # standard altitude
```

Now our directory contains three `atmos_average` files:

```bash
> 07180.atmos_average.nc 07180.atmos_average_pstd.nc 07180.atmos_average_zstd.nc
```

To see the variables in each file, use the `--inspect` function from `MarsPlot`:

```bash
(amesGCM3)>$ MarsPlot.py -i 07180.atmos_average.nc          # the original file, note that rho and zfull were added during postprocessing
(amesGCM3)>$ MarsPlot.py -i 07180.atmos_average_zstd.nc     # the pressure interpolated file
(amesGCM3)>$ MarsPlot.py -i 07180.atmos_average_pstd.nc     # the altitude interpolated file
```



***

#### 2.4 Add mass stream function (`msf`) to `atmos_average_pstd`

In this case, we add the variable after the interpolation because the mass stream function needs to be computed on a standard pressure grid.

```bash
(amesGCM3)>$ MarsVars.py 07180.atmos_average_pstd.nc -add msf
```




***

#### 2.5 Use `MarsFiles` to time-shift the diurn file, then pressure-interpolate the file

The variables in `07180.atmos_diurn.nc` are organized by time-of-day assuming **universal** time beginning at the prime martian meridian. You can time-shift the fields to **uniform** local time using `MarsFiles`. You might use this function when you want to compare MGCM output to observations from satellites in fixed local time orbit. 

For this exercise, we will only time-shift the surface pressure (`ps`), surface temperature (`ts`), and atmospheric temperature (`temp`) variables in order to minimize file size and processing time.

```bash
(amesGCM3)>$ MarsFiles.py 07180.atmos_diurn.nc -t --include ts ps temp
```

This function can only be performed on the `diurn` files since these are the only filetype containing hourly output. The `--tshift` function creates a new, time-shifted file (`07180.atmos_diurn_T.nc`) from the original file, and retains that original file. 

Now, pressure-interpolate the file using `MarsInterp` like we did for `atmos_average`:

```bash
(amesGCM3)>$ MarsInterp.py 07180.atmos_diurn_T.nc -t pstd
```

This should take just over a minute. Note that pressure-interpolating large files can take a long time, which is why we only included three variables (`ps`, `ts`, and `temp`) in this file. We now have three diurn filetypes:

```bash
> 07180.atmos_diurn.nc 07180.atmos_diurn_T.nc 07180.atmos_diurn_T_pstd.nc
```

> **Note:** We will *not* do this here, but you can create custom vertical grids that you want CAP to interpolate to. See the documentation for `MarsInterp` to learn how.




***

#### 2.6 Apply a low-pass filter (`-lpf`) to the surface pressure (`ps`) and temperature (`ts`) in the `atmos_daily` file

Use a 10-sol cut-off  frequency (specify `sol_max` > 10) to isolate synoptic-scale features. This will filter-out the pressure and save the variable in a new file:

```bash
(amesGCM3)>$ MarsFiles.py 07180.atmos_daily.nc -lpf 10 -include ps ts
```




#### 2.7 Estimate the magnitude of the wind shear using CAP

Add dU/dZ and dV/dZ to `07180.atmos_average_zstd.nc`.

In addition of adding new variables, `MarsVars` can apply certain operations, such as column integration and vertical differentiation, to existing variables. Apply vertical differentiation to add dU/dZ and dV/dZ to the file as follows:

```bash
(amesGCM3)>$ MarsVars.py 07180.atmos_average_zstd.nc -zdiff ucomp vcomp
```

You can use `--inspect` (`-i`) to find the names of the derived variables dU/dZ and dV/dZ:

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

> **The `--inspect` function works on any netCDF file, not just the ones created with CAP!**


#### 2.8 Display the minimum, mean, and maximum near-surface temperature

We can display values in an array by calling `--dump` with `MarsPlot -i` (analogue of the NCL command `ncdump`). For example, the content of the reference pressure (`pfull`) variable can be viewed by:

```bash
(amesGCM3)>$ MarsPlot.py -i 07180.atmos_average.nc -dump pfull
> pfull=
> [8.7662227e-02 2.5499690e-01 5.4266089e-01 1.0518962e+00 1.9545468e+00
> 3.5580616e+00 6.2466631e+00 1.0509957e+01 1.7400265e+01 2.8756382e+01
> 4.7480076e+01 7.8348366e+01 1.2924281e+02 2.0770235e+02 3.0938846e+02
> 4.1609518e+02 5.1308148e+02 5.9254102e+02 6.4705731e+02 6.7754218e+02
> 6.9152936e+02 6.9731799e+02 6.9994830e+02 7.0082477e+02]
> ______________________________________________________________________
```

We can index specific values in arrays by using quotes and square brackets with the variable call (`'var[ ]'`). For example, we can display the reference pressure of the first layer above the surface as follows:

```bash
(amesGCM3)>$ MarsPlot.py -i 07180.atmos_average.nc -dump 'pfull[-1]'
> pfull[-1]=
> 700.8247680664062
> ______________________________________________________________________
```

> **Note:** Use `-1` to refer to the last array element (Python syntax)

We can also use `-stat` to display the minimum, mean, and maximum values of a variable, which is better suited for visualizing statistics over large arrays or in data slices. Given the dimensions of the `temp` variable, (`[time,pfull,lat,lon]`), display the minimum, mean, and maximum near-surface air temperatures over all timesteps and all locations:

```bash
(amesGCM3)>$ MarsPlot.py -i 07180.atmos_average.nc -stat 'temp[:,-1,:,:]'
__________________________________________________________________________
           VAR            |      MIN      |      MEAN     |      MAX      |
__________________________|_______________|_______________|_______________|
            temp[:,-1,:,:]|        149.016|        202.508|         251.05|
__________________________|_______________|_______________|_______________|
```


> **Note:** quotes '' are necessary when browsing dimensions.


and that's it for post-processing the data in the `INERTCLDS/` simulation! Before we move on to plotting, we need to repeat some of these steps for the data in the `ACTIVECLDS/` simulation. Feel free to repeat all of Steps 2.1 through 2.8 for the `ACTIVECLDS/` case if you like, but **you are only required to repeat Steps 2.1, 2.2, and 2.3** for this tutorial:

* [2.1 Convert the `fort.11` files into `netCDF` files](#21-convert-the-fort11-files-into-netcdf-files-for-compatibility-with-cap)
* [2.2 Interpolate `atmos_average` to standard pressure coordinates](#22-interpolate-atmosaverage-to-standard-pressure-coordinates)
* [2.3 Add density (`rho`) and mid-point altitude (`zfull`) to the `atmos_average` file; interpolate to standard altitude](#23-add-density-rho-and-mid-point-altitude-zfull-to-atmosaverage-then-interpolate-the-file-to-standard-altitude-zstd)


### Navigate to the `ACTIVECLDS/` directory and complete Steps 2.1-2.3 before continuing





***

# Break!
Let's take a 15 minute break from the tutorial. You can use this time to catch up if you haven't completed parts 1 and 2 already, but we highly encourage you to step away from your machine for these 15 minutes.




***

## 3. Plotting Routines

The last part of this tutorial introduces the plotting capabilities of CAP. CAP can create several kinds of plots:

|Type of plot      |    MarsPlot designation|
|----------------------|-------------------|
| Longitude v Latitude | Plot 2D lon X lat|
| Longitude v Time     |Plot 2D lon X time|
| Longitude v Level    |Plot 2D lon X lev|
| Latitude v Level     |Plot 2D lat X lev|
| Time v Latitude      |Plot 2D time X lat|
| Time v level         |Plot 2D time X lev |
| Any 1-dimensional line plot |Plot 1D|

and CAP can display each plot on its own page or place multiple plots on the same page.

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

By default, `Custom.in` is set to create two plot types: a global map of topography and a zonal mean wind cross-section. The plot type is indicated at the top of each template:

``` python
<<<<<<<<<<<<<<| Plot 2D lon X lat = True |>>>>>>>>>>>>>
```

and if the plot type is set to `True` as it is here, then `MarsPlot` will draw that plot. The variable to be plotted is indicated after `Main Variable`. The variable name and the file containing it must be specified:

```python
 Main Variable  = fixed.zsurf # topography from the 01780.fixed.nc file
```

Without making any changes to `Custom.in`, close the file and pass it back to `MarsPlot` using the following command:

```bash
(amesGCM3)>$ MarsPlot.py Custom.in
```

This creates `Diagnostics.pdf`, a single-page PDF displaying the two plots we just discussed: global topography and the zonal mean wind cross-section. Open the PDF to see the plots.

You could rename `Custom.in` and still pass it to `MarsPlot` successfully. If the template is named anything other than `Custom.in`, `MarsPlot` will produce a PDF named after the template, i.e. `myplots.in` creates `myplots.pdf`. For example:

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

Those are the basics of plotting with CAP. We'll try creating several plot types in exercises 3.1--3.8 below. Begin by deleting `myplots.in` and `myplots.pdf` (if you have them), and then create a new `Custom.in` template:

```bash
(amesGCM3)>$ rm myplots.in myplots.pdf
(amesGCM3)>$ MarsPlot.py -template
```




***

#### 3.1 Create a global map of surface albedo (`alb`) with topography (`zsurf`) contoured on top

For this first plot, we'll edit `Custom.in` together. Open the template in your preferred text editor and make the following changes:

- Set the second default template (`Plot 2D lat X lev`) to `False` so that `MarsPlot` does not draw it (we will use it later)
- On the first template (`Plot 2D lon X lat`) set `Main Variable` to albedo (`alb`, located in the `fixed` file), this will be plotted as shaded contours
- Set `2nd Variable` to topography (`zsurf`, located in the `fixed` file), this will be plotted as solid contours
- Set the `Title` to reflect the variable(s) being plotted.

Here is what your template should look like:

```python
<<<<<<<<<<<<<<| Plot 2D lon X lat = True |>>>>>>>>>>>>>
Title          = 3.1: Albedo w/Topography Overplotted
Main Variable  = fixed.alb
Cmin, Cmax     = None
Ls 0-360       = None
Level [Pa/m]   = None
2nd Variable   = fixed.zsurf
Contours Var 2 = None
Axis Options  : lon = [None,None] | lat = [None,None] | cmap = binary | scale = lin | proj = cart
```

Save `Custom.in` (but don't close it!) and go back to the terminal. Pass `Custom.in` back to `MarsPlot`:

```bash
(amesGCM3)>$ MarsPlot.py Custom.in
```

Open `Diagnostics.pdf` and check to make sure it contains a global map of surface albedo with contoured topography.

> Depending on the settings for your specific PDF viewer, you may have to close and reopen the file to view it.



***

#### 3.2 Plot the zonal mean zonal wind cross-section at Ls=270° using altitude as the vertical coordinate

In the same template, which should still be open in your text editor, set the `2D lat X lev` template to `True`. Change `Main Variable` to point to `ucomp` stored in the `atmos_average_zstd` file. Edit the title accordingly, save `Custom.in`, and pass it to `MarsPlot`. 

Again, view `Diagnostics.pdf` to see your plots!




***

#### 3.3 Create the same plot for the RAC case. Place these plots on the same page

> **Tip:** Copy and paste the `lat x lev` plot so that you have two identical templates.

First, we need to point `MarsPlot` to the `ACTIVECLDS/` directory. we do this by editing the `<<<<<<< Simulations <<<<<<<` section so that `2>` points to `/ACTIVECLDS` like so:

```python
<<<<<<<<<<<<<<<<<<<<<< Simulations >>>>>>>>>>>>>>>>>>>>>
ref> None
2> ../ACTIVECLDS
```

Now we can edit `Main Variable` in one of the copy/pasted templates so that the `atmos_average` file in `ACTIVECLDS/` is sourced:

```python
Main Variable  = atmos_average@2.ucomp
```

> **Tip:** Make use of `HOLD ON` and `HOLD OFF` for these.

Save `Custom.in` and pass it to `MarsPlot`. View `Diagnostics.pdf` to see the results.




***

#### 3.4 Overplot temperature in solid contours

Add `temp` as the second variable on the plots you created in 3.2 and 3.3:

```python
<<<<<<<<<<<<<<| Plot 2D lat X lev = True |>>>>>>>>>>>>>
> 2nd Variable     = atmos_average_zstd.temp
(etc)

<<<<<<<<<<<<<<| Plot 2D lat X lev = True |>>>>>>>>>>>>>
> 2nd Variable     = atmos_average_zstd@2.temp
(etc)
```

Save `Custom.in` and pass it to `MarsPlot`. View `Diagnostics.pdf` to see the results.




***

#### 3.5 Plot the following four global maps (`lon X lat`) on a new page

> **Tip:** Use `HOLD ON` and `HOLD OFF` again. You can use this syntax multiple times in the same template.

Source all of the following variables from the `atmos_daily` file in `INERTCLDS/`. Plot the results at Ls=270.

- Surface CO2 Ice Content (`snow`) *north of 50 latitude*
- Surface Temperature (`ts`) *For this plot, set the colorscale (`Cmin, Cmax`) to range from 150 K to 300 K.*
- Surface Wind Speed (`(u^2 + v^2)/2`) (this requires the use of square brackets **and** two variables)
- Diabatic Heating Rate (`dheat`) at 50 Pa (index dimension `lev`=50).

The general format will be:

```python
HOLD ON

<<<<<<| Plot 2D lon X lat = True |>>>>>>
Title    = Surface CO2 Ice (g/m2)
(etc)

<<<<<<| Plot 2D lon X lat = True |>>>>>>
Title    = Surface Temperature (K)
(etc)

<<<<<<| Plot 2D lon X lat = True |>>>>>>
Title    = Surface Wind Speed (m/s)
(etc)

<<<<<<| Plot 2D lon X lat = True |>>>>>>
Title    = Diabatic Heating Rate (K/sol)
(etc)

HOLD OFF
```

> **Note: convert kg -> g using square brackets:**
>```python
>Main Variable  = [atmos_daily.snow]*1000
>```

> **Note: multiply two variables together like so:**
>```python
>Main Variable  = ([atmos_daily.ucomp]**2+[atmos_daily.vcomp]**2)**0.5
>```

Name the plots accordingly. Save `Custom.in` and pass it to `MarsPlot`. You should see four plots on one page in `Diagnostics.pdf`




***

#### 3.6 Plot the following two cross-sections (`lat X lev`) on the same page

- Mass Streamfunction (`msf`) at Ls=270. Change the colormap from `jet` to `bwr` and force symmetrical contouring by setting the colorbar's minimum and maximum values to -50 and 50. Adjust the y axis limits to 1,000 Pa and 1 Pa. Finally, add solid contours for `msf`=-10 and `msf`=10 on top. *Hint: set both `Main Variable` and `2nd Variable` to `msf`*
- Zonal mean temperature (`temp`) at Ls=270 from the same (pressure-interpolated) file. Overplot the zonal wind (`ucomp`).

Don't forget to use `HOLD ON` and `HOLD OFF` and to name your plots accordingly. Save `Custom.in` and pass it to `MarsPlot`.




***

#### 3.7 Plot the zonal mean temperature at Ls=270 from the `atmos_average` file for both the RIC and RAC cases. Also create a difference plot for them

Make use of `HOLD ON` and `HOLD OFF` again here. Copy and paste a `lat x lev` template three times. For the difference plot, you'll need to use `@N` to point to the `ACTIVECLDS/` directory and use square brackets to subtract one variable from the other:

```python
Main Variable  = [atmos_average_pstd.temp]-[atmos_average_pstd@2.temp]
```

Set the colormap to `RdBu` for the difference plot. For all three plots, set the vertical range to 1,000-1 Pa.

Save `Custom.in` and pass it to `MarsPlot`.




***

#### 3.8 Generate two **1D temperature profiles** (`temp`) from the RIC case, both at `50°N, 150°E` and Ls=270, at 3 AM and 3 PM 

There should be two lines on one plot. CAP can overplot 1D data on the same graph by concatenating two 1D templates together with `ADD LINE`:

```python
<<<<<<| Plot 1D = True |>>>>>>
Main Variable    = var1
(etc)

ADD LINE

<<<<<<| Plot 1D = True |>>>>>>
Main Variable    = var2
(etc)
```

#### You do not use `HOLD ON` and `HOLD OFF` to overplot 1D plots. `HOLD` is always for drawing separate plots on the same page

You'll need to call `temp` from the `diurn_T_pstd` file, which is the time-shifted and pressure-interpolated version of the hourly file. 3 AM is index=`3`, 3 PM is index=`15`, so `Main Variable` will be set as:

```python
Main Variable    = diurn_T_pstd.temp{tod=3}
```

You will have to specify `Level [Pa/m]` as the y axis:

```python
Level [Pa/m]   = AXIS
```

Save `Custom.in` and pass it to `MarsPlot`.






***

#### 3.9 Plot the filtered and un-filtered surface pressure over a 20 sol period

Some hints:
- Both are 1D plots. Use `ADD LINE` to plot on the same axes
- Use `ps` from `07180.atmos_daily.nc` and `07180.atmos_daily_lpf.nc`
- Set `Latitude = 50` and `Lon +/-180 = 150`
- Under `Axis Options`, set the x axis range (time) to 260--280 (`Ls = [260, 280]`)
- Under `Axis Options`, set the y axis range (pressure) to 850Pa--1000Pa(`var = [850, 1000]`)

Save `Custom.in` and pass it to `MarsPlot`.





***



## That's a Wrap!

This concludes the practical exercise portion of the CAP tutorial. Please keep these exercises as a reference for the future!

*This document was completed in October 2021. Written by Alex Kling, Courtney Batterson, and Victoria Hartwick*

Please submit feedback to Alex Kling: alexandre.m.kling@nasa.gov




***