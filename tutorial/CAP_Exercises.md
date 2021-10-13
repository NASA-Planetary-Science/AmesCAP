![](./tutorial_images/Tutorial_Banner_Final.png)

***

# Using the Community Analysis Pipeline (CAP)

Recap: CAP is a Python toolkit designed to simplify post-processing and plotting MGCM output.

Specifically, CAP consists of five subroutines that provide tools to perform the following functions:

1. `MarsPull.py`    Access MGCM output
2. `MarsFiles.py`   Reduce the files
3. `MarsVars.py`    Perform variable operations
4. `MarsInterp.py`  Interpolate the vertical grid
5. `MarsPlot.py`    Visualize the MGCM output

When learning to use CAP, it is useful to divide the above CAP functions into three parts:

1. [Retrieving Data](#1. Retrieving Data)
2. [File Manipulations](#2.-File-Manipulations)
3. [Plotting Routines](#Plotting-Routines)

We will practice using CAP for all three parts. You already have experience using CAP for Retrieving Data, which was covered at the end of the CAP installation instructions (the install asked you to use `MarsPull` to retrieve several `fort.11` files before the tutorial). Here, you will have a chance to practice using all five Python routines in CAP.


***


## 1. Retrieving Data
### Using `MarsPull.py` to download MGCM output

`MarsPull` is a utility for accessing MGCM output files hosted on the [MCMC Data portal](https://data.nas.nasa.gov/legacygcm/data_legacygcm.php). During the installation, you were asked to use `MarsPull` to download several `fort.11` files into your `INERTCLDS/` and `ACTIVECLDS/` directories. You can use `MarsPull` to download any file hosted on the MCMC Data portal. Simply pass a specific filename, Solar Longitude (L<sub>s</sub>), or a range of Solar Longitudes (L<sub>s</sub>) corresponding to the desired file(s):

```bash
(amesGCM3)>$ MarsPull.py -f LegacyGCM_LsXXX_LsYYY.nc
(amesGCM3)>$ MarsPull.py -ls XXX 
(amesGCM3)>$ MarsPull.py -ls XXX YYY 
```

Where XXX and YYY are three-digit L<sub>s</sub> values. You should have already downloaded the necessary `fort.11` files for this tutorial. If you haven't, please see item five from the installation instructions, *5. Do This Before Attending the Tutorial!*, for instructions.


***


## 2. File Manipulations

After retrieving output from the data portal or recieving output from a simulation, you will likely need to process the data to create the files you need for your analysis. Post-processing includes interpolating and regridding data to different coordinate systems, adding derived variables to the files, and converting between filetypes, just to name a few examples.

The following exercises are designed to demonstrate how CAP can be used for post-processing MGCM output. You should follow along using the `fort.11` files you downloaded during the installation process. **We will use the files we create here to make plots with MarsPlot**.

First, create a backup copy of the `INERTCLDS/` and `ACTIVECLDS/` directories in case you make a mistake during the tutorial. We recommend saving a copy of these folders periodically throughout the tutorial. 

Begin with the `INERTCLDS/` case and complete exercises 2.1-2.8 below. Then repeat the exercises with the `ACTIVECLDS/` case.


#### 2.1 Convert the `fort.11` files into `netCDF` files for compatibility with CAP.

Activate the virtual environment, go to the `INERTCLDS/` directory, and type:

```bash
(amesGCM3)>$ MarsFiles.py fort.11_* -fv3 fixed average daily diurn
```

Several `netCDF` filetypes have been created from the `fort.11` files in the directory. These are:

- `*atmos_fixed.nc`     Variables that do not change over time (e.g. albedo & topography maps)
- `*atmos_average.nc`   5-day averages of MGCM output
- `*atmos_diurn.nc`     Hourly MGCM output
- `*atmos_daily.nc`     Daily averaged MGCM output

For easier post-processing (and plotting in the future), combine like files to create one of each filetype:

```bash
(amesGCM3)>$ MarsFiles.py *fixed.nc -c
(amesGCM3)>$ MarsFiles.py *average.nc -c
(amesGCM3)>$ MarsFiles.py *diurn.nc -c
(amesGCM3)>$ MarsFiles.py *daily.nc -c
```

This created: 

```bash
> DDDDD.atmos_fixed.nc DDDDD.atmos_average.nc DDDDD.atmos_diurn.nc DDDDD.atmos_daily.nc
```

which are compatible for use with CAP.


#### 2.2 Interpolate `atmos_average` to standard pressure coordinates.

This requires using `MarsInterp`:

```bash
(amesGCM3)>$ MarsInterp.py DDDDD.atmos_average.nc -t pstd # standard pressure
```

and creates a new file:

```bash
> DDDDD.atmos_average_pstd.nc
```

in the same directory as the original `DDDDD.atmos_average.nc` file.

#### 2.3 Add density (`rho`) and (`zfull`) to `atmos_average`, then interpolate the file to standard altitude.
Adding or removing variables from files can be done with MarsVars:

```bash
(amesGCM3)>$ MarsVars.py DDDDD.atmos_average.nc -add rho zfull
```

This updates the original file to include the new variables.
> **NOTE: if you needed `rho` in any interpolated file, you need to add it before performing the interpolation because `rho` has to be computed on the native grid. In this case, we want rho in an altitude-interpolated file, and we perform the interpolation next.**

```bash
(amesGCM3)>$ MarsInterp.py DDDDD.atmos_average.nc -t zstd # standard altitude
```

Now our directory contains three `atmos_average` files:

```bash
> DDDDD.atmos_average.nc DDDDD.atmos_average_pstd.nc DDDDD.atmos_average_zstd.nc
```

To see the variables in each file, use the inspect function from `MarsPlot`:

```bash
(amesGCM3)>$ MarsPlot.py -i DDDDD.atmos_average.nc      # the original file
(amesGCM3)>$ MarsPlot.py -i DDDDD.atmos_average_zstd.nc # the pressure interpolated file
(amesGCM3)>$ MarsPlot.py -i DDDDD.atmos_average_pstd.nc # the altitude interpolated file
```

#### 2.4 Add mass stream function (`msf`) to `atmos_average_pstd`.
In this case, we add the variable after the interpolation because we need to compute streamfunction on pressure coordinates.

```bash
(amesGCM3)>$ MarsVars.py DDDDD.atmos_average_pstd.nc -add msf
```

#### 2.5 Use `MarsFiles` to time-shift the diurn file.
The variables in `DDDDD.atmos_diurn.nc` are organized by time-of-day, but you can time-shift the field to uniform local time using `MarsFiles`. You might use this function to allow plotting global variables at 3 AM and 3 PM, for example.

```bash
(amesGCM3)>$ MarsFiles.py DDDDD.atmos_diurn.nc -t
```

This function can only be performed on `diurn` files, since only `diurn` files contain hourly output. This function creates a new, time-shifted file, `DDDDD.atmos_diurn_t.nc`.

#### 2.6 Estimate the magnitude of the wind shear using CAP. Add $\frac{dU}{dZ}$ and $\frac{dV}{dZ}$ to the `atmos_average_zstd` file, then display the minimum, mean, and maximum values of each between the surface and 10 km.

```bash
(amesGCM3)>$ MarsVars DDDDD.atmos_average_zstd.nc -zdiff ucomp vcomp
```

You can use `-dump` with `MarsPlot` (kind of like `ncdump`) to determine the **index** at which `zstd` is 10 km:

```bash
(amesGCM3)>$ MarsPlot.py -i DDDDD.atmos_average_zstd.nc -dump zstd 
```

When you determine the index, use `-stat` to display the min, mean, and max values of \frac{dU}{dZ}$ and $\frac{dV}{dZ}$:

```bash
(amesGCM3)>$ MarsPlot.py -i DDDDD.atmos_average_zstd.nc -stat 'd_dz_ucomp[:,:15,:,:]'
(amesGCM3)>$ MarsPlot.py -i DDDDD.atmos_average_zstd.nc -stat 'd_dz_vcomp[:,:15,:,:]'
```

quotes '' are necessary for browsing dimensions.

> **Note:** We added `zfull` to the `atmos_average` file *before* interpolating to standard altitude. Then we computed the wind shear on the interpolated grid (`_zstd`).

#### 2.7 Compute the amplitude and phase of the semi-diurnal components of the tide (`N=2`) on the `atmos_diurn` file.
We can use `MarsFiles` with `-tidal N` (`N` denotes the tide harmonic) to create a file containing the N-diurnal components of the tide:

```bash
(amesGCM3)>$ MarsFiles.py DDDDD.atmos_diurn.nc -tidal 2
```

`N=1` is diurnal, `N=2` is semi diurnal, etc. By default, `MarsFiles` will perform the analysis on all available fields. You can specify only the fields you want (and significantly speed-up computing time) using `--include`:

```bash
(amesGCM3)>$ MarsFiles.py DDDDD.atmos_diurn.nc -tidal 2 --include ucomp vcomp temp
```

#### 2.8 Apply a low-pass filter (`-lpf`) to the surface temperature (`ts`) in the `atmos_daily` file over a period of at least 10 sols (set `sol_max` > 10).

This will filter out noise from the surface temperature variable and save the variable in a new file:

```bash
(amesGCM3)>$ MarsFiles.py DDDDD.atmos_daily.nc -lpf 10 -include ts         
```

CAP is capable of applying high-, low-, and band-pass filters to netCDF files using the syntax:

```bash
(amesGCM3)>$ MarsFiles.py file.nc -hpf --high_pass_filter sol_min          
(amesGCM3)>$ MarsFiles.py file.nc -lpf --low_pass_filter  sol_max          
(amesGCM3)>$ MarsFiles.py file.nc -bpf --band_pass_filter sol_min sol max  
```

Where `sol_min` and `sol_max` are the minimum and maximum number of days in a filtering period, respectively.




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
(amesGCM3)>$ MarsPlot.py myplots.in
```

If the template is named anything other than `Custom.in`, `MarsPlot` will produce a PDF named after the renamed template, i.e. `myplots.pdf`.

Those are the basics of plotting with CAP. We'll try creating several plot types in exercises 3.8--3.8 below.





#### 3.1 Plot a global map of thermal inertia (`thin`) with topography (`zsurf`) contoured on top.

For this first plot, we'll edit `Custom.in` together. Open the template in your preferred text editor and make the following changes:

- Change the second default template `Plot 2D lat X lev` to `False` so that `MarsPlot` does not draw it.
- Set the `Title` to reflect the variable being plotted
- Set `Main Variable` to thermal inertia (`thin`, located in the `fixed` file)
- Set `2nd Variable` to topography (`zsurf`, located in the `fixed` file)

Here is what your template should look like:

```python
<<<<<<<<<<<<<<| Plot 2D lon X lat = True |>>>>>>>>>>>>>
Title          = 3.1: Thermal Inertia
Main Variable  = fixed.thin
Cmin, Cmax     = None
Ls 0-360       = None
Level [Pa/m]   = None
2nd Variable   = fixed.zsurf
Contours Var 2 = None
Axis Options  : lon = [None,None] | lat = [None,None] | cmap = binary | scale = lin | proj = cart
```

Save the template, then go to your terminal and run `MarsPlot`:

```bash
(amesGCM3)>$ MarsPlot.py Custom.in
```

Open `Diagnostics.pdf` and check to make sure it contains a global map of surface thermal inertia and topography.



#### 3.2 Next, plot a cross-section of the zonal mean zonal wind at Ls=270°.

No need to create a new template, just add this plot to `Custom.in`. Use the zonal wind stored in the `atmos_average` file. Edit the title accordingly.

Save `Custom.in` and pass it to `MarsPlot`.




#### 3.3 Create the same plot for the radiatively active cloud case.

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

Save `Custom.in` and pass it to `MarsPlot`.





#### 3.4 Add temperature as solid contours overtop of the zonal wind plot.

Add `temp` as a second variable on the plot you created in 3.2:

```python
> 2nd Variable     = atmos_average.temp
```

Save `Custom.in` and pass it to `MarsPlot`.





#### 3.5 Plot the following four global maps (`lon x lat`) on the same page:

> Tip: Make use of `HOLD ON` and `HOLD OFF` for these, and Copy/Paste plot types to create multiple of the same plot.

- Zonal mean surface temperature (`ts`)
- Infrared dust optical depth (`taudust_IR`). *For this plot, set the colorscale (`Cmin, Cmax`) to range from 0 to 0.1.*
- Column ice content (`cldcol`) 
- Column water vapor (`wcol`) *in units of [pr-um]* (this requires the use of square brackets for element-wise operations)

The general format will be:

```python
HOLD ON

Zonal mean surface temperature plot
IR dust optical depth plot
Column ice plot
Column water vapor plot

HOLD OFF
```

and you can convert column water vapor from its native units (kg/m$^2$) to pr-um by multiplying `wcol` by 1,000:

```python
Main Variable  = [atmos_average.wcol]*1000
```

Name the plots accordingly. Save `Custom.in` and pass it to `MarsPlot`.






#### 3.4 Plot the following two cross-sections (`lat x lev`) on the same page:

> Tip: Add to your existing template. Copy and paste the `lat x lev` plot three times. Set the plots to `True` so that MarsPlot recognizes them as input.

- Mass Streamfunction (`msf`) at Ls=270. Set the colormap to `bwr` and force symmetrical contouring by setting the colorbar's minimum and maximum values to -50 and 50. Adjust the y axis limits to 1,0000 Pa and 1 Pa. Finally, add solid contours for `msf` on top. *Hint: set both `Main Variable` and `2nd Variable` to `msf`*
- Zonal mean temperature (`temp`) at Ls=270 from the same (pressure-interpolated) file. Overplot the zonal wind (`ucomp`).

Don't forget to use `HOLD ON` and `HOLD OFF` and to name your plots accordingly. Save `Custom.in` and pass it to `MarsPlot`.





#### 3.5 Plot the zonal mean temperature at Ls=270 from the daily file for the intert cloud case and the active cloud case. Also create a difference plot for them.

Use `HOLD ON` and `HOLD OFF`, copy and paste a `lat x lev` plot three times, and set the colorbar to the range 120--250.

For the difference plot, you'll need to use `@N` to point to the `ACTIVECLDS/` directory and square brackets to subtract one variable from the other:

```python
Main Variable  = [atmos_daily@2.temp]-[atmos_daily.temp]
```

Save `Custom.in` and pass it to `MarsPlot`.





#### 3.6 Generate a **1D temperature profile** (`temp`) at `50°N, 150°E` and 3 AM and 3 PM from the radiatively active case. Plot these on the same plot.

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

You'll need to call `temp` from the `diurn_T` file, which is the time-shifted version of the hourly file. Point to the file in the `ACTIVECLDS/` directory. To index the time dimension of `temp`, use curly brackets `{}`:

```python
Main Variable  = atmos_diurn_t@2.temp{tod=15}
```

3 AM is index=3, 3 PM is index=15. Save `Custom.in` and pass it to `MarsPlot`.




#### 3.7 Plot the amplitude of the thermal tide at 6 AM and Ls=270. Use orthographic projection centered over `50°N, 150°E`. Create two plots: one for each simulation.

```python
Main Variable  = atmos_diurn_tidal.temp_amp{tod=6}
```

To change the projection, edit `Axis Options` with the following:

```python
proj = ortho -150, 50
```

Save `Custom.in` and pass it to `MarsPlot`.



#### 3.8 Plot the filtered and un-filtered noon surface temperature over a 20-Ls period.

Some hints:
- Both are 1D plots. Use `ADD LINE` to plot on the same axes
- Use `ts` from the `DDDDD.atmos_daily_lpf.nc` file
- Index noon `{tod=12}` 
- Set `Latitude = 50` and `Lon +/-180 = 150`
- Under `Axis Options`, set the y axis range (temperature) to 150K--190K (`var = [150, 190]`) 
- Under `Axis Options`, set the x axis range (time) to 260--280 (`sols = [260, 280]`)

Save `Custom.in` and pass it to `MarsPlot`.





***



## That's a Wrap!

This concludes the practical exercise portion of the CAP tutorial. Please keep these exercises as a reference for the future!


***


This document was completed in October 2021. Written by Alex Kling, Courtney Batterson, and Victoria Hartwick

Please submit feedback to Alex Kling: alexandre.m.kling@nasa.gov


