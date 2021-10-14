![](./tutorial_images/Tutorial_Banner_Final.png)




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

**Activation of the Community Analysis Pipeline**
As always with CAP, you first start by activating the CAP virtual environment (you can revisit the [installation instructions](https://github.com/alex-kling/amesgcm/blob/master/tutorial/CAP_install.md) as a refresher).

```bash
(local)>$ source amesGCM3/bin/activate      # bash
(local)>$ source amesGCM3/bin/activate.csh  # csh/tcsh
```

For each Mars executable, we recommend you check the `--help` argument (`-h` for short) to see which documentation, for example:

```bash
(amesGCM3)>$ MarsPull.py -h
```
***

## 1. Retrieving Data
### Using `MarsPull.py` to download MGCM output

`MarsPull` is a utility for accessing MGCM output files hosted on the [MCMC Data portal](https://data.nas.nasa.gov/legacygcm/data_legacygcm.php). During the installation, you were asked to use `MarsPull` to download several `fort.11` files into your `INERTCLDS/` and `ACTIVECLDS/` directories. You can use `MarsPull` to download any file hosted on the MCMC Data portal. Specify a simulation identifier, and a Solar Longitude (L<sub>s</sub>), or a range of Solar Longitudes (L<sub>s</sub>) corresponding to the desired file(s):


```bash
(amesGCM3)>$ MarsPull.py -id INERTCLDS -ls XXX YYY
```

Where XXX and YYY are three-digit L<sub>s</sub> values. You should have already downloaded the necessary `fort.11` files for this tutorial. If you haven't, please see item five from the installation instructions, *5. Do This Before Attending the Tutorial!*, for instructions. There should be a total of ten files : fort.11_0719, fort.11_0720, fort.11_0721, fort.11_0722 and fort.11_0723 for EACH simulation.

> If you have downloaded other fort.11 files in addition of the five (ten total) listed above,  we recommend you copy the five files above in dedicated folders for the purpose of the tutorial. It  will make it easier to follow-along the tutorial if you work with a small subset of the year-long simulation.


***

## 2. File Manipulations

After retrieving output from the data portal or using output from a simulation you ran yourself, you will likely need to process the data to create the files you need for your analysis. Post-processing includes interpolating and regridding data to different vertical coordinate systems, adding derived variables to the files, and converting between filetypes, just to name a few examples.

The following exercises are designed to demonstrate how CAP can be used for post-processing MGCM output. You should follow along in the directories you created containing the `fort.11` files you downloaded during the installation process. After post-processing these files, **we will use them to make plots with MarsPlot**. Don't delete anything!


We will start with the `INERTCLDS/` simulation (radiatively inert clouds) and complete exercises 2.1-2.8 below. Repeat the exercises for the `ACTIVECLDS/` (radiatively active clouds) afterward. In section 3, we access both simulations to make plots.


***

#### 2.1 Convert the `fort.11` files into `netCDF` files for compatibility with CAP.

To do this, activate the virtual environment, go to the `INERTCLDS/` directory, and type:

```bash
(amesGCM3)>$ MarsFiles.py -h # display documentation
(amesGCM3)>$ MarsFiles.py fort.11_* -fv3 fixed average daily diurn
```

This created a bunch of `netCDF` files. Your directory should look like this:

```bash
(amesGCM3)>$ ls
> 00490.atmos_average.nc  00500.atmos_average.nc  00510.atmos_average.nc  00520.atmos_average.nc  00530.atmos_average.nc  fort.11_0719            fort.11_0723
> 00490.atmos_daily.nc    00500.atmos_daily.nc    00510.atmos_daily.nc    00520.atmos_daily.nc    00530.atmos_daily.nc    fort.11_0720
> 00490.atmos_diurn.nc    00500.atmos_diurn.nc    00510.atmos_diurn.nc    00520.atmos_diurn.nc    00530.atmos_diurn.nc    fort.11_0721
> 00490.fixed.nc          00500.fixed.nc          00510.fixed.nc          00520.fixed.nc          00530.fixed.nc          fort.11_0722
```

Several `netCDF` filetypes are located in the directory:

- `*atmos_fixed.nc` files contain static variables that **do not change over time** (e.g. albedo & topography maps)
- `*atmos_average.nc` files contain **5-day average** of MGCM output
- `*atmos_diurn.nc` files contain **hourly** MGCM output, also averaged over 5 days
- `*atmos_daily.nc` files contain **continuous time series** of the MGCM output, these are the most voluminous files

For easier post-processing and plotting, we can combine alike files along the time axis to create one of each filetype:

```bash
(amesGCM3)>$ MarsFiles.py *fixed.nc -c
(amesGCM3)>$ MarsFiles.py *average.nc -c
(amesGCM3)>$ MarsFiles.py *diurn.nc -c
(amesGCM3)>$ MarsFiles.py *daily.nc -c
```

This merged like-filetypes and created the four following files:

```bash
> 00490.atmos_fixed.nc 00490.atmos_average.nc 00490.atmos_diurn.nc 00490.atmos_daily.nc
```





***

#### 2.2 Interpolate `atmos_average` to standard pressure coordinates.

This requires using `MarsInterp`:

```bash
(amesGCM3)>$ MarsInterp.py -h # display documentation
(amesGCM3)>$ MarsInterp.py 00490.atmos_average.nc -t pstd # standard pressure
```

and creates a new file:

```bash
> 00490.atmos_average_pstd.nc
```

in the same directory as the original `00490.atmos_average.nc` file.




***

#### 2.3 Add density (`rho`) and mid-point altitude (`zfull`) to `atmos_average`, then interpolate the file to standard altitude (`zstd`)

Adding or removing variables from files can be done with `MarsVars`:

```bash
(amesGCM3)>$ MarsVars.py -h # display documentation
(amesGCM3)>$ MarsVars.py 00490.atmos_average.nc -add rho zfull
```

This updates the original file to include the new variables. In this case, the density `rho` was derived from the pressure and temperature (which are already present in the file) and the mid-point altitude `zfull` was obtained through hydrostatic integration.

> **NOTE: if you want `rho` in an interpolated file, you need to add it before performing the interpolation because. In this case, we want `rho` in an altitude-interpolated file so we've added `rho` to the original file (`atmos_average.nc`) and we will perform the interpolation next .**

```bash
(amesGCM3)>$ MarsInterp.py 00490.atmos_average.nc -t zstd   # standard altitude
```

Now our directory contains three `atmos_average` files:

```bash
> 00490.atmos_average.nc 00490.atmos_average_pstd.nc 00490.atmos_average_zstd.nc
```

To see the variables in each file, use the `--inspect` function from `MarsPlot`:

```bash
(amesGCM3)>$ MarsPlot.py -i 00490.atmos_average.nc          # the original file, note that rho and zfull were added during postprocessing
(amesGCM3)>$ MarsPlot.py -i 00490.atmos_average_zstd.nc     # the pressure interpolated file
(amesGCM3)>$ MarsPlot.py -i 00490.atmos_average_pstd.nc     # the altitude interpolated file
```



***

#### 2.4 Add mass stream function (`msf`) to `atmos_average_pstd`.

In this case, we add the variable after the interpolation because the mass stream function needs to be computed on a standard pressure grid.

```bash
(amesGCM3)>$ MarsVars.py 00490.atmos_average_pstd.nc -add msf
```




***

#### 2.5 Use `MarsFiles` to time-shift the diurn file, then pressure-interpolate the file.
The variables in `00490.atmos_diurn.nc` are organized by time-of-day in universal time at the prime martian meridian, but you can time-shift the fields to uniform local time using `MarsFiles`. You might use this function to allow plotting global variables at 3 AM and 3 PM, for example.

```bash
(amesGCM3)>$ MarsFiles.py 00490.atmos_diurn.nc -t
```

This function can only be performed on `diurn` files, since only `diurn` files contain hourly output. This function creates a new, time-shifted file, `00490.atmos_diurn_T.nc`. Next, pressure interpolate the file using `MarsInterp` (like we did for `atmos_average`):

```bash
(amesGCM3)>$ MarsInterp.py 00490.atmos_diurn_T.nc -t pstd
```

We now have three diurn filetypes:
```bash
> 00490.atmos_diurn.nc 00490.atmos_diurn_T.nc 00490.atmos_diurn_T_pstd.nc
```




***

#### 2.6 Estimate the magnitude of the wind shear using CAP. Add dU/dZ and dV/dZ to `00490.atmos_average_zstd.nc` and then display the minimum, mean, and maximum values for each variable between the surface and 10 km.

In addition of adding new variables, `MarsVars` can apply certain operations such as column integration or vertical differentiation to existing variables. We will use the later as follows:

```bash
(amesGCM3)>$ MarsVars.py 00490.atmos_average_zstd.nc -zdiff ucomp vcomp
```

You can use `--inspect` (`-i`) to find the names of the derived variables dU/dZ and dV/dZ:

```bash
(amesGCM3)>$ MarsPlot.py -i 00490.atmos_average_zstd.nc
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

Now, use `-dump` with `MarsPlot` (analogue of the NCL command `ncdump`) to print the altitude array to the terminal and determine the **index** at which `zstd` is 10 km:

```bash
(amesGCM3)>$ MarsPlot.py -i 00490.atmos_average_zstd.nc -dump zstd
> zstd=
> [ -7000.  -6000.  -5000.  -4500.  -4000.  -3500.  -3000.  -2500.  -2000.
>   -1500.  -1000.   -500.      0.    500.   1000.   1500.   2000.   2500.
>    3000.   3500.   4000.   4500.   5000.   6000.   7000.   8000.   9000.
>   10000.  12000.  14000.  16000.  18000.  20000.  25000.  30000.  35000.
>   40000.  45000.  50000.  55000.  60000.  70000.  80000.  90000. 100000.]
> ______________________________________________________________________
```

We can verify that the layer corresponding to an altitude of 10km (10000m) is the 28th element in `zstd` (i=27 since Python's indexing start at i=0)

```bash
(amesGCM3)>$ MarsPlot.py -i 00490.atmos_average_zstd.nc -dump 'zstd[27]'
>zstd[27]=
>10000.0
> ______________________________________________________________________
```

When you determine the index, use `-stat` to display the min, mean, and max values of dU/dZ, in that specific altitude slice  (`:27` following Python's convention)

```bash
(amesGCM3)>$ MarsPlot.py -i 00490.atmos_average_zstd.nc -stat 'd_dz_ucomp[:,:27,:,:]'
>__________________________________________________________________________
>           VAR            |      MIN      |      MEAN     |      MAX      |
>__________________________|_______________|_______________|_______________|
>     d_dz_ucomp[:,:27,:,:]|     -0.0161607|     0.00188501|      0.0278951|
>__________________________|_______________|_______________|_______________|
>
(amesGCM3)>$ MarsPlot.py -i 00490.atmos_average_zstd.nc -stat 'd_dz_vcomp[:,:27,:,:]'
```

Do the same for dV/dZ.

> **Note:** quotes '' are necessary when browsing dimensions.

> **Note:** We added `zfull` to the `atmos_average` file *before* interpolating to standard altitude. Then we computed the wind shear on the interpolated grid (`_zstd`).




***

#### 2.7 Compute the amplitude and phase of the semi-diurnal components of the tide (`N=2`) on the `atmos_diurn` file.
We can use `MarsFiles` with `-tidal N` (`N` denotes the tide harmonic) to create a file containing the N-diurnal components of the tide:

```bash
(amesGCM3)>$ MarsFiles.py 00490.atmos_diurn.nc -tidal 2
```

`N=1` is diurnal, `N=2` is semi diurnal, etc. By default, `MarsFiles` will perform the analysis on all available fields. To speed-up computing time, you can specify only the fields you are interested in (here the surface pressure `ps`, and atmospheric temperature `temp`)  using `--include`:

```bash
(amesGCM3)>$ MarsFiles.py 00490.atmos_diurn.nc -tidal 2 --include ps temp
```


***

#### 2.8 Apply a low-pass filter (`-lpf`) to the surface temperature (`ps`) in the `atmos_daily` with a 10 sols cut-off  frequency (set `sol_max` > 10) to isolate synoptic-scale feature.

This will filter-out the pressure and save the variable in a new file:

```bash
(amesGCM3)>$ MarsFiles.py 00490.atmos_daily.nc -lpf 10 -include ps         
```

CAP is capable of applying high-, low-, and band-pass filters to netCDF files using the syntax:

```bash
(amesGCM3)>$ MarsFiles.py file.nc -hpf --high_pass_filter sol_min          
(amesGCM3)>$ MarsFiles.py file.nc -lpf --low_pass_filter  sol_max          
(amesGCM3)>$ MarsFiles.py file.nc -bpf --band_pass_filter sol_min sol max  
```

Where `sol_min` and `sol_max` are the minimum and maximum number of days in a filtering period, respectively.


### Remember to repeat this post-processing on the `ACTIVECLDS/` simulation as well!



***

# Break!
Let's take a 15 minute break from the tutorial. You can use this time to catch up if you haven't completed parts 1 and 2 already, but we highly encourage you to step away from your machine for these 15 minutes.




***

## 3. Plotting Routines

The last part of this tutorial covers the plotting capabilities within CAP. CAP can create several kinds of plots:

- Longitude v Latitude
- Longitude v Time
- Longitude v Level
- Latitude v Level
- Time v Latitude
- Time v level
- Any 1-dimensional line plot

and CAP can display each plot on its own page or place multiple plots on the same page.

Plotting with CAP requires passing a template to `MarsPlot`. A blank template is created in the directory in which the following command is executed, so change to the `INERTCLDS/` directory and type:

```bash
(amesGCM3)>$ MarsPlot.py -template
```

The blank template is called `Custom.in`. Pass `Custom.in` back to `MarsPlot` using the following command:

```bash
(amesGCM3)>$ MarsPlot.py Custom.in
```

This will have created `Diagnostics.pdf`, a single-page PDF with a topographical plot and a cross-section of the zonal mean wind.

You can rename `Custom.in` and still pass it to `MarsPlot` successfully:

```bash
(amesGCM3)>$ mv Custom.in myplots.in
(amesGCM3)>$ MarsPlot.py myplots.in
```

If the template is named anything other than `Custom.in`, `MarsPlot` will produce a PDF named after the renamed template, i.e. `myplots.pdf`.

Those are the basics of plotting with CAP. We'll try creating several plot types in exercises 3.8--3.8 below.




***

#### 3.1 Plot a global map of surface albedo (`alb`) with topography (`zsurf`) contoured on top.

For this first plot, we'll edit `Custom.in` together. Open the template in your preferred text editor and make the following changes:

- Change the second default template `Plot 2D lat X lev` to `False` so that `MarsPlot` does not draw it.
- Set the `Title` of the first default template `Plot 2D lon X lat` to reflect the variable being plotted.
- Set `Main Variable` to thermal inertia (`thin`, located in the `fixed` file)
- Set `2nd Variable` to topography (`zsurf`, located in the `fixed` file)

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

Save the template and pass it back to `MarsPlot`:

```bash
(amesGCM3)>$ MarsPlot.py Custom.in
```

Open `Diagnostics.pdf` and check to make sure it contains a global map of surface albedo and topography.




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

All of the following variables come from `00490.atmos_daily.nc` and should be plotted at Ls=270.

- Surface CO2 ice content (`snow`) *north of 50 latitude*
- Surface temperature (`ts`) *For this plot, set the colorscale (`Cmin, Cmax`) to range from 150 K to 300 K.*
- Surface Wind Speed (`(u^2 + v^2)/2`) (this requires the use of square brackets **and** two variables)
- Diabatic Heating Rate (`dheat`) at the 14th vertical layer (index dimension `pfull`=14).

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

#### 3.4 Plot the following two cross-sections (`lat x lev`) on the same page:

- Mass Streamfunction (`msf`) at Ls=270. Change the colormap from `jet` to `bwr` and force symmetrical contouring by setting the colorbar's minimum and maximum values to -50 and 50. Adjust the y axis limits to 1,000 Pa and 1 Pa. Finally, add solid contours for `msf`=-10 and `msf`=10 on top. *Hint: set both `Main Variable` and `2nd Variable` to `msf`*
- Zonal mean temperature (`temp`) at Ls=270 from the same (pressure-interpolated) file. Overplot the zonal wind (`ucomp`).

Don't forget to use `HOLD ON` and `HOLD OFF` and to name your plots accordingly. Save `Custom.in` and pass it to `MarsPlot`.




***

#### 3.5 Plot the zonal mean temperature at Ls=270 from the daily file for the intert cloud case and the active cloud case. Also create a difference plot for them.

Use `HOLD ON` and `HOLD OFF`. Copy and paste a `lat x lev` plot three times. For the difference plot, you'll need to use `@N` to point to the `ACTIVECLDS/` directory and square brackets to subtract one variable from the other:

```python
Main Variable  = [atmos_average_pstd.temp]-[atmos_average_pstd@2.temp]
```

Set the colormap to `RdBu` for the difference plot and set the vertical range to 1,000-1 Pa.

Save `Custom.in` and pass it to `MarsPlot`.




***

#### 3.6 Generate a **1D temperature profile** (`temp`) at `50°N, 150°E` at Ls=270 at both 3 AM and 3 PM from the radiatively inert case. Plot these on the same plot.

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

#### 3.7 Plot the amplitude and phase of the thermal tide at 6 AM and Ls=270. Use orthographic projection centered over `50°N, 150°E`. Create two plots: one for each simulation.

This will be a `lat X lon` plot. Use `temp_amp` and `temp_phas` from `atmos_diurn_tidal`. Remember to use the `{tod=6}` syntax. To change the projection *and* specify the centerpoint, edit `Axis Options` with the following:

```python
proj = ortho 150, 50
```

Save `Custom.in` and pass it to `MarsPlot`.




***

#### 3.8 Plot the filtered and un-filtered noon surface temperature over a 20-Ls period.

Some hints:
- Both are 1D plots. Use `ADD LINE` to plot on the same axes
- Use `ts` from the `00490.atmos_daily.nc` and `00490.atmos_daily_lpf.nc` files
- Index noon `{tod=12}`
- Set `Latitude = 50` and `Lon +/-180 = 150`
- Under `Axis Options`, set the y axis range (temperature) to 150K--190K (`var = [150, 190]`)
- Under `Axis Options`, set the x axis range (time) to 260--280 (`sols = [260, 280]`)

Save `Custom.in` and pass it to `MarsPlot`.




***

***



## That's a Wrap!

This concludes the practical exercise portion of the CAP tutorial. Please keep these exercises as a reference for the future!




***

This document was completed in October 2021. Written by Alex Kling, Courtney Batterson, and Victoria Hartwick

Please submit feedback to Alex Kling: alexandre.m.kling@nasa.gov




***
