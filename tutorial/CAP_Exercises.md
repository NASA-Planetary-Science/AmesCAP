![](./tutorial_images/Tutorial_Banner_Final.png)


<!-- TOC titleSize:2 tabSpaces:2 depthFrom:1 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 skip:0 title:1 charForUnorderedList:* -->
## Table of Contents
* [Practical: Using the Community Analysis Pipeline (CAP)](#practical-using-the-community-analysis-pipeline-cap)
  * [1. Retrieving Data](#1-retrieving-data)
    * [Using `MarsPull.py` to download MGCM output](#using-marspullpy-to-download-mgcm-output)
  * [2. File Manipulations](#2-file-manipulations)
      * [2.1 Convert the `fort.11` files into `netCDF` files for compatibility with CAP.](#21-convert-the-fort11-files-into-netcdf-files-for-compatibility-with-cap)
      * [2.2 Interpolate `atmos_average` to standard pressure coordinates.](#22-interpolate-atmosaverage-to-standard-pressure-coordinates)
      * [2.3 Add density (`rho`) and mid-point altitude (`zfull`) to `atmos_average`, then interpolate the file to standard altitude (`zstd`)](#23-add-density-rho-and-mid-point-altitude-zfull-to-atmosaverage-then-interpolate-the-file-to-standard-altitude-zstd)
      * [2.4 Add mass stream function (`msf`) to `atmos_average_pstd`.](#24-add-mass-stream-function-msf-to-atmosaveragepstd)
      * [2.5 Use `MarsFiles` to time-shift the diurn file, then pressure-interpolate the file.](#25-use-marsfiles-to-time-shift-the-diurn-file-then-pressure-interpolate-the-file)
      * [2.6 Apply a low-pass filter (`-lpf`) to the surface pressure (`ps`) and temperature (`ts`) in the `atmos_daily` with a 10 sols cut-off  frequency (set `sol_max` > 10) to isolate synoptic-scale feature.](#26-apply-a-low-pass-filter--lpf-to-the-surface-pressure-ps-and-temperature-ts-in-the-atmosdaily-with-a-10-sols-cut-off--frequency-set-solmax--10-to-isolate-synoptic-scale-feature)
      * [2.7 Estimate the magnitude of the wind shear using CAP. Add dU/dZ and dV/dZ to `07180.atmos_average_zstd.nc`.](#27-estimate-the-magnitude-of-the-wind-shear-using-cap-add-dudz-and-dvdz-to-07180atmosaveragezstdnc)
      * [2.8 Display the minimum, mean, and maximum near-surface temperature .](#28-display-the-minimum-mean-and-maximum-near-surface-temperature-)
    * [Remember to repeat this post-processing on the `ACTIVECLDS/` simulation as well!](#remember-to-repeat-this-post-processing-on-the-activeclds-simulation-as-well)
* [Break!](#break)
  * [3. Plotting Routines](#3-plotting-routines)
      * [3.1 Plot a global map of surface albedo (`alb`) with topography (`zsurf`) contoured on top.](#31-plot-a-global-map-of-surface-albedo-alb-with-topography-zsurf-contoured-on-top)
      * [3.2 Next, plot a cross-section of the zonal mean zonal wind at Ls=270° using altitude as the vertical coordinate.](#32-next-plot-a-cross-section-of-the-zonal-mean-zonal-wind-at-ls270-using-altitude-as-the-vertical-coordinate)
      * [3.3 Create the same plot for the radiatively active cloud case, and put both zonal mean zonal wind plots on their own page.](#33-create-the-same-plot-for-the-radiatively-active-cloud-case-and-put-both-zonal-mean-zonal-wind-plots-on-their-own-page)
      * [3.4 Add temperature as solid contours overtop of the zonal wind plot.](#34-add-temperature-as-solid-contours-overtop-of-the-zonal-wind-plot)
      * [3.5 Plot the following four global maps (`lon x lat`) on a new page:](#35-plot-the-following-four-global-maps-lon-x-lat-on-a-new-page)
      * [3.6 Plot the following two cross-sections (`lat x lev`) on the same page:](#36-plot-the-following-two-cross-sections-lat-x-lev-on-the-same-page)
      * [3.7 Plot the zonal mean temperature at Ls=270 from the average file for the inert cloud case and the active cloud case. Also create a difference plot for them.](#37-plot-the-zonal-mean-temperature-at-ls270-from-the-average-file-for-the-inert-cloud-case-and-the-active-cloud-case-also-create-a-difference-plot-for-them)
      * [3.8 Generate a **1D temperature profile** (`temp`) at `50°N, 150°E` at Ls=270 at both 3 AM and 3 PM from the radiatively inert case. Plot these on the same plot.](#38-generate-a-1d-temperature-profile-temp-at-50n-150e-at-ls270-at-both-3-am-and-3-pm-from-the-radiatively-inert-case-plot-these-on-the-same-plot)
      * [3.9 Plot the filtered and un-filtered surface pressure over a 20 sol period.](#39-plot-the-filtered-and-un-filtered-surface-pressure-over-a-20-sol-period)
  * [That's a Wrap!](#thats-a-wrap)
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

**Activate CAP**
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
### Using `MarsPull.py` to download MGCM output

`MarsPull` is a utility for accessing MGCM output files hosted on the [MCMC Data portal](https://data.nas.nasa.gov/legacygcm/data_legacygcm.php). During the installation, you were asked to use `MarsPull` to download several `fort.11` files into your `INERTCLDS/` and `ACTIVECLDS/` directories. You should have already downloaded the necessary `fort.11` files for this tutorial. If you haven't, you can do so now by following along with the instructions below.

We asked that you create a `CAP_Tutorial` directory containing two subdirectories, `INERTCLDS/` and `ACTIVECLDS/`, in `amesGCM3/`:

```bash
(amesGCM3)>$ cd ~/amesGCM3
(amesGCM3)>$ mkdir CAP_Tutorial
(amesGCM3)>$ cd CAP_Tutorial
(amesGCM3)>$ mkdir ACTIVECLDS INERTCLDS
```

Navigate to `INERTCLDS/` and use `MarsPull` to retrieve the files. Specify the simulation identifier (`INERTCLDS/`) and the range of Solar Longitudes (255 285) corresponding to the desired file(s):

```bash
(amesGCM3)>$ MarsPull.py -id INERTCLDS -ls 255 285
```

Then, do the same for the `ACTIVECLDS/` case.

There should now be 5 `fort.11` files in each directory, `INERTCLDS/` and `ACTIVECLDS/`:

```bash
> fort.11_0719 fort.11_0720 fort.11_0721 fort.11_0722 fort.11_0723
```

> If you have any `fort.11` files *other than the ones listed above* in *either* directory, please delete them. It will be make it easier to follow the tutorial if you work with the specific subset of files listed above.




***

## 2. File Manipulations

After retrieving output from the data portal or using output from a simulation you ran yourself, you will likely need to process the data to create the files you need for your analysis. Post-processing includes interpolating and regridding data to different vertical coordinate systems, adding derived variables to the files, and converting between filetypes, just to name a few examples.

The following exercises are designed to demonstrate how CAP can be used for post-processing MGCM output. You should follow along in the directories you created containing the `fort.11` files you downloaded during the installation process. After post-processing these files, **we will use them to make plots with MarsPlot**. Don't delete anything!


Start with the radiatively inert clouds simulation (RIC), `INERTCLDS/`, and complete exercises 2.1-2.8 below. Then we will give you specific instructions regarding which exercises to repeat for the radiatively active clouds (RAC) simulation `ACTIVECLDS/`. We access files from both simulations to make plots in Section 3.




***

#### 2.1 Convert the `fort.11` files into `netCDF` files for compatibility with CAP.

To do this, go to your `INERTCLDS/` directory, and type:

```bash
(amesGCM3)>$ MarsFiles.py fort.11_* -fv3 fixed average daily diurn
```

This created several `netCDF` files:

```bash
(amesGCM3)>$ ls
> 07180.atmos_average.nc  07190.atmos_average.nc  07200.atmos_average.nc  07210.atmos_average.nc  07220.atmos_average.nc  fort.11_0719            fort.11_0723
> 07180.atmos_daily.nc    07190.atmos_daily.nc    07200.atmos_daily.nc    07210.atmos_daily.nc    07220.atmos_daily.nc    fort.11_0720
> 07180.atmos_diurn.nc    07190.atmos_diurn.nc    07200.atmos_diurn.nc    07210.atmos_diurn.nc    07220.atmos_diurn.nc    fort.11_0721
> 07180.fixed.nc          07190.fixed.nc          07200.fixed.nc          07210.fixed.nc          07220.fixed.nc          fort.11_0722
```
> Note  the  five-digit  sol  numbers at the begining of each netcdf file, which corresponds to the time at the begining of each fort.11 output. Because the  simulation is issued from a 10 year run (10 x ~668 sols/year), this  particular series of  outputs start at 06690, not  00000.


The `netCDF` filetypes are:

| Type                  | Description |
| --------------------- | ----------- |
| `*atmos_fixed.nc`     | static variables that **do not change over time**           |
| `*atmos_average.nc`   | **5-day averages** of MGCM output                           |
| `*atmos_diurn.nc`     | files contain **hourly** MGCM output averaged over 5 days   |
| `*atmos_daily.nc`     | **continuous time series** of the MGCM output               |

For easier post-processing and plotting, we can combine like files along the time axis. This creates one of each filetype:

```bash
(amesGCM3)>$ MarsFiles.py *fixed.nc -c
(amesGCM3)>$ MarsFiles.py *average.nc -c
(amesGCM3)>$ MarsFiles.py *diurn.nc -c
(amesGCM3)>$ MarsFiles.py *daily.nc -c
```

This merge created the following four files:

```bash
> 07180.atmos_fixed.nc 07180.atmos_average.nc 07180.atmos_diurn.nc 07180.atmos_daily.nc
```





***

#### 2.2 Interpolate `atmos_average` to standard pressure coordinates.

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

This updates the original file to include the new variables. In this case, the density `rho` was derived from the pressure and temperature (which are already present in the file) and the mid-point altitude `zfull` was obtained through hydrostatic integration.

> **NOTE: if you want `rho` in an interpolated file, you need to add it before performing the interpolation because. In this case, we want `rho` in an altitude-interpolated file so we've added `rho` to the original file (`atmos_average.nc`) and we will perform the interpolation next .**

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

#### 2.4 Add mass stream function (`msf`) to `atmos_average_pstd`.

In this case, we add the variable after the interpolation because the mass stream function needs to be computed on a standard pressure grid.

```bash
(amesGCM3)>$ MarsVars.py 07180.atmos_average_pstd.nc -add msf
```




***

#### 2.5 Use `MarsFiles` to time-shift the diurn file, then pressure-interpolate the file.
The variables in `07180.atmos_diurn.nc` are organized by time-of-day in universal time at the prime martian meridian, but you can time-shift the fields to uniform local time using `MarsFiles`. You might use this function to allow plotting global variables at 3 AM and 3 PM, for example. We will only retain the surface pressure `ps`, surface temperature `ts` and atmospheric temperature `temp` using `--include` to minimize the size of the file and processing time.

```bash
(amesGCM3)>$ MarsFiles.py 07180.atmos_diurn.nc -t --include ts ps temp
```

This function can only be performed on `diurn` files, since only `diurn` files contain hourly output. This function creates a new, time-shifted file, `07180.atmos_diurn_T.nc`. Next, pressure interpolate the file using `MarsInterp` (like we did for `atmos_average`).

```bash
(amesGCM3)>$ MarsInterp.py 07180.atmos_diurn_T.nc -t pstd
```

This should take just over a minute. Note that pressure interpolating large files can take a long time which is why we only included `ps`, `ts`, and `temp` in this file. We now have three diurn filetypes:

```bash
> 07180.atmos_diurn.nc 07180.atmos_diurn_T.nc 07180.atmos_diurn_T_pstd.nc
```

> **Note:** We will *not* do this here, but you can specify a vertical grid to interpolate to with CAP. See the documentation for `MarsInterp.py` to learn how.




***

#### 2.6 Apply a low-pass filter (`-lpf`) to the surface pressure (`ps`) and temperature (`ts`) in the `atmos_daily` with a 10 sols cut-off  frequency (set `sol_max` > 10) to isolate synoptic-scale feature.

This will filter-out the pressure and save the variable in a new file:

```bash
(amesGCM3)>$ MarsFiles.py 07180.atmos_daily.nc -lpf 10 -include ps ts         
```




#### 2.7 Estimate the magnitude of the wind shear using CAP. Add dU/dZ and dV/dZ to `07180.atmos_average_zstd.nc`.
In addition of adding new variables, `MarsVars` can apply certain operations such as column integration or vertical differentiation to existing variables. Vertical differentiation can be done as follows:

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

> **The `--inspect` function works on any netCDF file, not just the ones created here!**


#### 2.8 Display the minimum, mean, and maximum near-surface temperature .

We can display values in an array by calling `--dump` with `MarsPlot -i` (analogue of the NCL command `ncdump`). For example, the content for the reference pressure (`pfull` variable in the file) is:

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

We can also index specific values using quotes and square brackets `'[ ]'`. For example, we can display the reference pressure in the first layer above the surface ( we use `-1` to refer to the last array element per Python convention):

```bash
(amesGCM3)>$ MarsPlot.py -i 07180.atmos_average.nc -dump 'pfull[-1]'
> pfull[-1]=
> 700.8247680664062
> ______________________________________________________________________
```

 `-stat` display the min, mean, and max values of a variable, which is better suited to display statistics over a large array or for specific data-slices. For example, to display the min, mean, and max air temperature for all timesteps, all latitudes, all longitudes, and near the surface (`[time,pfull,lat,lon]=[:,-1,:,:]`), we use:

```bash
(amesGCM3)>$ MarsPlot.py -i 07180.atmos_average.nc -stat 'temp[:,-1,:,:]'
__________________________________________________________________________
           VAR            |      MIN      |      MEAN     |      MAX      |
__________________________|_______________|_______________|_______________|
            temp[:,-1,:,:]|        149.016|        202.508|         251.05|
__________________________|_______________|_______________|_______________|
```


> **Note:** quotes '' are necessary when browsing dimensions.



### Remember to repeat this post-processing on the `ACTIVECLDS/` simulation as well!



***

# Break!
Let's take a 15 minute break from the tutorial. You can use this time to catch up if you haven't completed parts 1 and 2 already, but we highly encourage you to step away from your machine for these 15 minutes.




***

## 3. Plotting Routines

The last part of this tutorial covers the plotting capabilities within CAP. CAP can create several kinds of plots:

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

Plotting with CAP requires passing a template to `MarsPlot`. A blank template is created in the directory in which the following command is executed, so change to the `INERTCLDS/` directory and type:

```bash
(amesGCM3)>$ MarsPlot.py -template
```

The blank template is called `Custom.in`. Pass `Custom.in` back to `MarsPlot` using the following command:

```bash
(amesGCM3)>$ MarsPlot.py Custom.in
```

This will have created `Diagnostics.pdf`, a single-page PDF with a topographical plot and a cross-section of the zonal mean wind. Open the pdf to see the plots.

> You can rename `Custom.in` and still pass it to `MarsPlot` successfully:
```bash
(amesGCM3)>$ mv Custom.in myplots.in
(amesGCM3)>$ MarsPlot.py myplots.in
```
If the template is named anything other than `Custom.in`, `MarsPlot` will produce a PDF named after the renamed template, i.e. `myplots.pdf`.


Those are the basics of plotting with CAP. We'll try creating several plot types in exercises 3.8--3.8 below.




***

#### 3.1 Plot a global map of surface albedo (`alb`) with topography (`zsurf`) contoured on top.

For this first plot, we'll edit `Custom.in` together. Open the template in your preferred text editor and make the following changes:

- Change the second default template `Plot 2D lat X lev` to `False` so that `MarsPlot` does not draw it (we will use it later)
- Set the `Title` of the first default template `Plot 2D lon X lat` to reflect the variable being plotted.
- Set `Main Variable` to albedo (`alb`, located in the `fixed` file), this will be plotted as shaded contours
- Set `2nd Variable` to topography (`zsurf`, located in the `fixed` file), this will be plotted as solid contours

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

Save the template in your text editor and pass it back to `MarsPlot`:

```bash
(amesGCM3)>$ MarsPlot.py Custom.in
```

Open `Diagnostics.pdf` and check to make sure it contains a global map of surface albedo and topography.
> Depending on the settings for your specific pdf viewer, you may have to close and open the file.



***

#### 3.2 Next, plot a cross-section of the zonal mean zonal wind at Ls=270° using altitude as the vertical coordinate.

No need to create a new template, just add this plot to `Custom.in`. Use the zonal wind stored in the `atmos_average_zstd` file. Remember to set the plot template to `True`. Edit the title accordingly.

Save `Custom.in` and pass it to `MarsPlot`.




***

#### 3.3 Create the same plot for the radiatively active cloud case, and put both zonal mean zonal wind plots on their own page.

> Tip: Add to your existing template. Copy and paste the `lat x lev` plot three times. Set the plots to `True` so that `MarsPlot` recognizes them as input.

Edit the `<<<<<<< Simulations <<<<<<<` section so that
`2>` points to the `/ACTIVECLDS` directory:

```python
<<<<<<<<<<<<<<<<<<<<<< Simulations >>>>>>>>>>>>>>>>>>>>>
ref> None
2> ../ACTIVECLDS
```

Then, copy and paste the plot created in 3.2 and edit `Main Variable` to point to the correct directory using the `@N` syntax:

```python
Main Variable  = atmos_average@2.ucomp
```

> Tip: Make use of `HOLD ON` and `HOLD OFF` for these, and Copy/Paste plot types to create multiple of the same plot.

Save `Custom.in` and pass it to `MarsPlot`.




***

#### 3.4 Add temperature as solid contours overtop of the zonal wind plot.

Add `temp` as a second variable on the plots you created in 3.2 and 3.3:

```python
> 2nd Variable     = atmos_average_zstd.temp
```

Save `Custom.in` and pass it to `MarsPlot`.




***

#### 3.5 Plot the following four global maps (`lon x lat`) on a new page:

> Tip: Use `HOLD ON` and `HOLD OFF`. You can use this syntax multiple times in the same template.

All of the following variables come from `07180.atmos_daily.nc` and should be plotted at Ls=270.

- Surface CO2 ice content (`snow`) *north of 50 latitude*
- Surface temperature (`ts`) *For this plot, set the colorscale (`Cmin, Cmax`) to range from 150 K to 300 K.*
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

> *Note:* convert kg -> g using square brackets:
>```python
>Main Variable  = [atmos_daily.snow]*1000
>```
> and multiply two variables together like so:
>```python
>Main Variable  = ([atmos_daily.ucomp]**2+[atmos_daily.vcomp]**2)**0.5
>```

Name the plots accordingly. Save `Custom.in` and pass it to `MarsPlot`.




***

#### 3.6 Plot the following two cross-sections (`lat x lev`) on the same page:

- Mass Streamfunction (`msf`) at Ls=270. Change the colormap from `jet` to `bwr` and force symmetrical contouring by setting the colorbar's minimum and maximum values to -50 and 50. Adjust the y axis limits to 1,000 Pa and 1 Pa. Finally, add solid contours for `msf`=-10 and `msf`=10 on top. *Hint: set both `Main Variable` and `2nd Variable` to `msf`*
- Zonal mean temperature (`temp`) at Ls=270 from the same (pressure-interpolated) file. Overplot the zonal wind (`ucomp`).

Don't forget to use `HOLD ON` and `HOLD OFF` and to name your plots accordingly. Save `Custom.in` and pass it to `MarsPlot`.




***

#### 3.7 Plot the zonal mean temperature at Ls=270 from the average file for the inert cloud case and the active cloud case. Also create a difference plot for them.

Use `HOLD ON` and `HOLD OFF`. Copy and paste a `lat x lev` plot three times. For the difference plot, you'll need to use `@N` to point to the `ACTIVECLDS/` directory and square brackets to subtract one variable from the other:

```python
Main Variable  = [atmos_average_pstd.temp]-[atmos_average_pstd@2.temp]
```

Set the colormap to `RdBu` for the difference plot and set the vertical range to 1,000-1 Pa.

Save `Custom.in` and pass it to `MarsPlot`.




***

#### 3.8 Generate a **1D temperature profile** (`temp`) at `50°N, 150°E` at Ls=270 at both 3 AM and 3 PM from the radiatively inert case. Plot these on the same plot.

CAP can overplot 1D data on the same graph by concatenating two 1D templates together with `ADD LINE`:

```python
<<<<<<| Plot 1D = True |>>>>>>
Main Variable    = var1
(etc)

ADD LINE

<<<<<<| Plot 1D = True |>>>>>>
Main Variable    = var2
(etc)
```

> You do not need to use `HOLD ON` or `HOLD OFF` with 1D plots.

You'll need to call `temp` from the `diurn_T_pstd` file, which is the time-shifted and pressure-interpolated version of the hourly file. 3 AM is index=3, 3 PM is index=15. You will have to specify `Level [Pa/m]` as the y axis:

```python
Level [Pa/m]   = AXIS
```

Save `Custom.in` and pass it to `MarsPlot`.






***

#### 3.9 Plot the filtered and un-filtered surface pressure over a 20 sol period.

Some hints:
- Both are 1D plots. Use `ADD LINE` to plot on the same axes
- Use `ps` from the `07180.atmos_daily.nc` and `07180.atmos_daily_lpf.nc` files
- Index noon `{tod=12}`
- Set `Latitude = 50` and `Lon +/-180 = 150`
- Under `Axis Options`, set the x axis range (time) to 260--280 (`sols = [260, 280]`)
- Under `Axis Options`, set the y axis range (pressure) to 850Pa--1000Pa(`var = [850, 1000]`)

Save `Custom.in` and pass it to `MarsPlot`.





***



## That's a Wrap!

This concludes the practical exercise portion of the CAP tutorial. Please keep these exercises as a reference for the future!




***

This document was completed in October 2021. Written by Alex Kling, Courtney Batterson, and Victoria Hartwick

Please submit feedback to Alex Kling: alexandre.m.kling@nasa.gov




***
