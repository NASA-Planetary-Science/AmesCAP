![](./tutorial_images/Tutorial_Banner_Final.png)



***

# Introducing the Community Analysis Pipeline (CAP)

CAP is a toolkit designed to simplify the post-processing of MGCM output. CAP is written in Python and works with existing Python libraries, allowing any Python user to install and use CAP easily and free of charge. Without CAP, plotting MGCM output requires that a user provide their own scripts for post-processing, including code for interpolating the vertical grid, computing and adding derived variables to files, converting between file types, and creating diagnostic plots. In other words, a user would be responsible for the entire post-processing effort as illustrated in Figure 1.

![Figure 1. The Typical Pipeline](./tutorial_images/Typical_Pipeline.png)

Such a process requires that users be familiar with Fortran files and be able to write (or provide) script(s) to perform file manipulations and create plots. CAP standardizes the post-processing effort by providing executables that can perform file manipulations and create diagnostic plots from the command line. This enables users of almost any skill level to post-process and plot MGCM data (Figure 2).

![Figure 2. The New Pipeline (CAP)](./tutorial_images/CAP.png)


As a foreword, we will list a few design characteristic of CAP:

* CAP is written in **Python**, an open-source programming language with extensive scientific libraries available
* CAP is installed within a Python **virtual environment**, which provides cross-platform support (MacOS, Linux and Windows), robust version control (packages updated within the main Python distribution will not affect CAP), and is not intrusive as it disappears when deactivated
* CAP is composed of a **set of libraries** (functions), callable from a user's own scripts and a collection of **five executables**, which  allows for efficient processing of model outputs from the command-line.
* CAP uses the **netCDF4 data format**, which is widely use in the climate modeling community and self-descriptive (meaning that a file contains  explicit information about its content in term of variables names, units etc...)
* CAP uses a convention for output formatting inherited from the GFDL Finite­-Volume Cubed-Sphere Dynamical Core, referred here as  "**FV3 format**": outputs may be binned and averaged in time in various ways for analysis.  
*  CAP long-term goal is to offer **multi-model support**. At the time of the writing, both the NASA Ames Legacy GCM and the NASA Ames GCM with the FV3 dynamical core are  supported. Efforts are underway to offer functionality to others Global Climate Models (e.g. eMARS, LMD, MarsWRF).


Specifically, CAP consists of five executables:

1. `MarsPull.py`    Access MGCM output
2. `MarsFiles.py`   Reduce the files
3. `MarsVars.py`    Perform variable operations
4. `MarsInterp.py`  Interpolate the vertical grid
5. `MarsPlot.py`    Visualize the MGCM output


These executables and their commonly-used functions are illustrated in the cheat sheet below in the order in which they are most often used. You should feel free to reference during and after the tutorial.

# Cheat sheet

![Figure 3. Quick Guide to Using CAP](./tutorial_images/Cheat_Sheet.png)

CAP is designed to be modular. For example, a user could post-process and plot MGCM output exclusively with CAP or a user could employ their own post-processing routine and then use CAP to plot the data. Users are free to selectively integrate CAP into their own analysis routine to the extent they see fit.



***
# The big question... How do I do this? >  <span style="color:red">Ask for help!  </span>
Use the `--help` (`-h` for short) option on any executable to display documentation and examples.
```
(amesGCM3)>$ MarsPlot.py -h
> usage: MarsPlot.py [-h] [-i INSPECT_FILE] [-d DATE [DATE ...]] [--template]
>                   [-do DO] [-sy] [-o {pdf,eps,png}] [-vert] [-dir DIRECTORY]
>                   [--debug]
>                   [custom_file]
```


***

# 1. `MarsPull.py` - Downloading Raw MGCM Output

`MarsPull` is a utility for accessing MGCM output files hosted on the [MCMC Data portal](https://data.nas.nasa.gov/legacygcm/data_legacygcm.php). MGCM data is archived in 1.5 hour intervals (16x/day) and packaged in files containing 10 sols. The files are named fort.11_XXXX in the order they were produced, but  `MarsPull` maps those files to specific solar longitudes (L<sub>s</sub>, in °). This allows users to request a file at a specific L<sub>s</sub> or for a range of L<sub>s</sub> using the `-ls` flag. Additionally the `identifier` (`-id`) flag is used to route `MarsPull` through a particular simulation. The `filename` (`-f`) flag can be used to parse specific files within a particular directory.

```bash
MarsPull.py -id INERTCLDS -ls 255 285
MarsPull.py -id ACTIVECLDS -f fort.11_0720 fort.11_0723
```
[Back to Top](#cheat-sheet)
***

# 2. `MarsFiles.py` - Reducing the Files

`MarsFiles` provides several tools for file manipulations, including code designed to create binned, averaged, and time-shifted files from MGCM output. The `-fv3` flag is used to convert fort.11 binaries to the Netcdf data format (you can select one or more of the file format listed below):

```bash
(amesGCM3)>$ MarsFiles.py fort.11* -fv3 fixed average daily diurn
```

These are the file formats that `MarsFiles` can create from the fort.11 MGCM output files.
**Primary files**

| File name | Description                                    |Timesteps for 10 sols x 16 output/sol           |Ratio to daily file (430Mb)|
|-----------|------------------------------------------------|----------------------------------------------- | ---                  |
|**atmos_daily.nc** | continuous time series                | (16 x 10)=160                                  | 1                    |
|**atmos_diurn.nc** | data binned by time of day and 5-day averaged | (16 x 2)=32                                    | x5 smaller           |
|**atmos_average.nc** | 5-day averages                              |  (1 x 2) = 2                                           | x80 smaller          |
|**fixed.nc** | statics variable such as surface albedo and topography  |  static                                        |few kB                |


**Secondary files**


| File name | description|
|-----------|------------|
|daily**_lpf**,**_hpf**,**_bpf** |low, high and band pass filtered|
|diurn**_T** |uniform local time (same time of day at all longitudes)|
|diurn**_tidal** |tidally-decomposed files into  harmonics|
|daily**_to_average**  **_to_diurn** |custom re-binning of daily files|

- `MarsFiles` can concatenate like-files together on the time dimension using the `-combine` (`-c`) flag.

```bash
> 07180.atmos_average.nc  07190.atmos_average.nc  07200.atmos_average.nc # 3 files with 10 days of output each
(amesGCM3)>$ MarsFiles.py *atmos_average.nc -c
> 07180.atmos_average.nc  # 1 file with 30 days of output
```

![Figure X. MarsFiles options](./tutorial_images/MarsFiles_diurn.png)
*3pm surface temperature before (left) and after (right) processing a diurn file with MarsFile to uniform local time*


[Back to Top](#cheat-sheet)
***

# 3. `MarsVars.py` - Performing Variable Operations

`MarsVars` provides several tools relating to variable operations such as adding and removing variables and performing column integrations. With no other arguments, passing a file to `MarsVars` displays file content much like `ncdump`:

```bash
(amesGCM3)>$ MarsVars.py 00000.atmos_average.nc
>
> ===================DIMENSIONS==========================
> ['bnds', 'time', 'lat', 'lon', 'pfull', 'scalar_axis', 'phalf']
> (etc)
> ====================CONTENT==========================
> pfull          : ('pfull',)= (30,), ref full pressure level  [Pa]
> temp           : ('time', 'pfull', 'lat', 'lon')= (4, 30, 180, 360), temperature  [K]
> (etc)
```

A typical option of `MarsVars` would be to add the atmospheric density `rho` to a file. Because the density is easily computed from the pressure and temperature fields, we do not archive in in the GCM output and instead provides a utility to add it as needed. This conservative approach to logging output allows to  minimize disk space and speed-up post processing.


```bash
(amesGCM3)>$ MarsVars.py 00000.atmos_average.nc -add rho
```

We can see that `rho` was added by calling `MarsVars` with no argument as before:

```bash
(amesGCM3)>$ MarsVars.py 00000.atmos_average.nc
>
> ===================DIMENSIONS==========================
> ['bnds', 'time', 'lat', 'lon', 'pfull', 'scalar_axis', 'phalf']
> (etc)
> ====================CONTENT==========================
> pfull          : ('pfull',)= (30,), ref full pressure level  [Pa]
> temp           : ('time', 'pfull', 'lat', 'lon')= (4, 30, 180, 360), temperature  [K]
> rho            : ('time', 'pfull', 'lat', 'lon')= (4, 30, 180, 360), density (added postprocessing)  [kg/m3]
```

The `help` (`-h`) option provides information on available variables and needed fields for each operation.

![Figure X. MarsVars](./tutorial_images/MarsVars.png)

`MarsVars` also offers the following variable operations:


| Command | flag| action|
|-----------|-----|-------|
|add | -add  | add a variable to the file|
|remove |-rm| remove a variable from a file|
|extract |-extract | extract a list of variables to a new file |
|col |-col | column integration, applicable to mixing ratios in [kg/kg] |
|zdiff |-zdiff |vertical differentiation (e.g. compute gradients)|
|zonal_detrend |-zd | zonally detrend a variable|

[Back to Top](#cheat-sheet)
***

# 4. `MarsInterp.py` - Interpolating the Vertical Grid

Native MGCM output files use pressure as the vertical coordinate (`pfull`), which means the geometric height and pressure level of an atmospheric layer varies based on location.

![Figure X. MarsInterp](./tutorial_images/MarsInterp.png)

*Pressure interpolation from the reference pressure grid to a standard pressure grid*


Climate data is usually analyzed on a standardized grid, however, and it is often necessary to interpolate the files to standard pressure coordinates. The `-type` (`-t`) argument in `MarsInterp` can interpolate files for you:

```bash
(amesGCM3)>$ MarsInterp.py  00000.atmos_average.nc -t pstd
```

An inspection of the file shows that the pressure level axis which was `pfull` (30 layers) has been replaced by a standard pressure coordinate `pstd` (36 layers), and all 3- and 4-dimensional variables reflect the new shape:

```bash
(amesGCM3)>$ MarsInterp.py  00000.atmos_average.nc -t pstd
(amesGCM3)>$ MarsVars.py 00000.atmos_average_pstd.nc
>
> ===================DIMENSIONS==========================
> ['bnds', 'time', 'lat', 'lon', 'scalar_axis', 'phalf', 'pstd']
> ====================CONTENT==========================
> pstd           : ('pstd',)= (36,), pressure  [Pa]
> temp           : ('time', 'pstd', 'lat', 'lon')= (4, 36, 180, 360), temperature  [K]
```

The following `type` (`-t` flag) of vertical interpolation are supported:

| file type | description | low-level value in a deep crater
|-----------|-----------|--------|
|_pstd | standard pressure [Pa]  |  1000Pa
|_zstd | standard altitude [m]   |  -7000m
|_zagl | standard altitude above ground level [m]   | 0 m

***

**Use of custom vertical grids**

It is also possible for the users to specify the layers for the interpolation. This is done by editing a **hidden** file `.amesgcm_profile`(note the dot '`.`) in your home directory.  

For the first use, you will need  to copy a template of `amesgcm_profile` to your /home directory:

```bash
(amesGCM3)>$ cp ~/amesGCM3/mars_templates/amesgcm_profile ~/.amesgcm_profile # Note the dot '.' !!!
```
You can open `~/.amesgcm_profile` with any text editor:

```
> <<<<<<<<<<<<<<| Pressure definitions for pstd |>>>>>>>>>>>>>

>p44=[1.0e+03, 9.5e+02, 9.0e+02, 8.5e+02, 8.0e+02, 7.5e+02, 7.0e+02,
>       6.5e+02, 6.0e+02, 5.5e+02, 5.0e+02, 4.5e+02, 4.0e+02, 3.5e+02,
>       3.0e+02, 2.5e+02, 2.0e+02, 1.5e+02, 1.0e+02, 7.0e+01, 5.0e+01,
>       3.0e+01, 2.0e+01, 1.0e+01, 7.0e+00, 5.0e+00, 3.0e+00, 2.0e+00,
>       1.0e+00, 5.0e-01, 3.0e-01, 2.0e-01, 1.0e-01, 5.0e-02, 3.0e-02,
>       1.0e-02, 5.0e-03, 3.0e-03, 5.0e-04, 3.0e-04, 1.0e-04, 5.0e-05,
>       3.0e-05, 1.0e-05]
>
>phalf_mb=[50]
```
In the example above, the user custom-defined two vertical grids, one with 44 levels (named `p44`) and one with a single layer at 50 Pa =0.5mbar(named `phalf_mb`)

You can use these by calling `MarsInterp` with the `-level` (`-l`) argument followed by the name of the new grid in `amesgcm_profile`.

```bash
(amesGCM3)>$ MarsInterp.py  00000.atmos_average.nc -t pstd -l  p44
```
[Back to Top](#cheat-sheet)
***

# 5. `MarsPlot.py` - Plotting the Results



The last component of CAP is the plotting routine, `MarsPlot`, which accepts a modifiable template (`Custom.in`) containing a list of plots to create. `MarsPlot` is useful for creating plots from MGCM output quickly, and it is designed specifically for use with the `netCDF` output files (`daily`, `diurn`, `average`, `fixed`).

The following figure shows the three components of MarsPlot:
- *MarsPlot.py*, opened in **a terminal** to inspect the netcdf files and ingest the Custom.in template
- *Custom.in* , a template opened in **a text editor**
- *Diagnostics.pdf*, refreshed in a **pdf viewer**

![Figure 4. MarsPlot workflow](./tutorial_images/MarsPlot_graphics.png)

The default template, Custom.in, can be created by passing the `-template` argument to `MarsPlot`. Custom.in is pre-populated to draw two plots on one page: a topographical plot from the fixed file and a cross-section of the zonal wind from the average file. Creating the template and passing it into `MarsPlot` creates a PDF containing the plots:

```
(amesGCM3)>$ MarsPlot.py -template
> /path/to/simulation/run_name/history/Custom.in was created
(amesGCM3)>$
(amesGCM3)>$ MarsPlot.py Custom.in
> Reading Custom.in
> [----------]  0 % (2D_lon_lat :fixed.zsurf)
> [#####-----] 50 % (2D_lat_lev :atmos_average.ucomp, Ls= (MY 2) 252.30, zonal avg)
> [##########]100 % (Done)
> Merging figures...
> /path/to/simulation/run_name/history/Diagnostics.pdf was generated
```

Specifically MarsPlot is designed to generate 2D cross - sections and 1D plots.
Let's remind ourselves that in order to create such plots from a **multi-dimensional** dataset, we first need to specify the **free** dimensions, meaning the ones that are **not** plotted.


![Figure 4. MarsPlot cross section](./tutorial_images/cross_sections.png)

*A refresher on cross-section for multi-dimensional datasets*

The data selection process to make any particular cross section is shown in the decision tree below. If an effort to make the process of generating multiple plots as **streamlined** as possible, MarsPlot selects a number of default settings for the user.

```
1.     Which simulation                                              ┌─
    (e.g. ACTIVECLDS directory)                                      │  DEFAULT    1. ref> is current directory
          │                                                          │  SETTINGS
          └── 2.   Which XXXXX epoch                                 │             2. latest XXXXX.fixed in directory
               (e.g. 00668, 07180)                                   └─
                   │                                                 ┌─
                   └── 3.   Which type of file                       │
                        (e.g. diurn, average_pstd)                   │   USER      3. provided by user
                            │                                        │ PROVIDES
                            └── 4.   Which variable                  │             4. provided by user
                                  (e.g. temp, ucomp)                 └─
                                    │                                ┌─
                                    └── 5. Which dimensions          │             5. see rule table below
                                       (e.g lat =0°,Ls =270°)        │  DEFAULT
                                           │                         │  SETTINGS
                                           └── 6. plot customization │             6. default settings
                                                  (e.g. colormap)    └─              

```

The free dimensions are set by default using day-to-day decisions from a climate modeler's perspective:


|Free dimension|Statement for default setting|Implementation|
|--|------|---------------|
|time |"*I am interested in the most recent events*"                             |time = Nt (last timestep)|
|level|"*I am more interested in the surface than any other vertical layer*"       |level = sfc|
|latitude  |"*If I have to pick a particular latitude, I would rather look at the equator*" |lat=0 (equator)|
|longitude  |"*I am more interested in a zonal average than any particular longitude*"      |lon=all (average over all values)|
|time of day| "*3pm =15hr Ok, this one is arbitrary. However if I use a diurn file, I  have a specific time of day in mind*"   |tod=15 |

*Rule table for the default setting of the free dimensions*

In practice, these cases cover 99% of the work typically done so whenever a setting is left to default (`= None` in MarsPlot's syntax) this is what is being used. This allows to considerably streamline the data selection process.

`Custom.in` can be modified using your preferred text editor (and renamed to your liking). This is an example of the code snippet in `Custom.in` used to generate a lon/lat cross-section. Note that the heading is set to `= True`, so that plot is activated for MarsPlot to process.

```python
<<<<<<<<<<<<<<| Plot 2D lon X lat = True |>>>>>>>>>>>>>
Title          = None
Main Variable  = atmos_average.temp
Cmin, Cmax     = None
Ls 0-360       = None
Level [Pa/m]   = None
2nd Variable   = None
Contours Var 2 = None
Axis Options  : lon = [None,None] | lat = [None,None] | cmap = jet | scale = lin | proj = cart

```
In the example above, we are plotting the air temperature field `temp` from the *atmos_average.nc* file as a lon/lat map. `temp` is a 4D field *(time, level, lat, lon)* but since we left the time (`Ls 0-360`) and altitude (`Level [Pa/m]`) unspecified (i.e. set to `None`) MarsPlot  will show us the *last timestep* in the file and the layer immediately adjacent to the *surface*. Similarly, MarsPlot will generate a *default title* for the figure with the variable's name (`temperature`), unit (`[K]`), selected dimensions (`last timestep, at the surface`), and makes educated choices for the range of the colormap, axis limits etc ... All those options are customizable, if desired.  Finally, note the option of adding a secondary variable as **solid contours**. For example, one may set `2nd Variable = fixed.zsurf` to plot the topography (`zsurf`) from the matching *XXXXX.fixed.nc* file.

To wrap-up (the use of `{}` to overwrite default settings is discussed later on), the following two working expressions are strictly equivalent for `Main Variable` (shaded contours) or `2nd Variable` (solid contours) fields:

```python
                     variable                                        variable
                        │                     SIMPLIFY TO               │
00668.atmos_average@1.temp{lev=1000;ls=270}     >>>      atmos_average.temp
  │         │       │              │                           │
epoch  file type simulation    free dimensions             file type
                 directory
```


|Accepted input |Meaning| Example|
|--         |-       |--|
|`None` |Use default settings from the rule table below| `Ls 0-360 = None`|
|`value`|  Return index closest to requested value in figure's unit |`Level [Pa/m]  = 50 ` (50 Pa)|
|`Val Min, Val Max`| Return the averages between two values |`Lon +/-180 = -30,30`|
|`all`| `all` is a special keyword that return the average over all values in file |`Latitude       = all`

*Accepted values for the `Ls 0-360`, `Level [Pa/m]` ,`Lon +/-180` and `Latitude` free dimensions*

> The time of day (`tod`) in diurn files is always specified using brackets`{}`, e.g. : `Main Variable = atmos_diurn.temp{tod=15,18}` for the average between 3pm and 6pm. This has allowed to streamlined all templates by not including the *time of day* free dimension, which is specific to diurn files.



# MarsPlot.py:  How to?
## Disable or add a new plot
Code blocks is set to `= True` instruct `MarsPlot` to draw those plots. Other templates in `Custom.in` are set to `= False` by default, which instructs `MarsPlot` to skip those plots. In total, `MarsPlot` is equipped to create seven plot types:

```python
<<<<<| Plot 2D lon X lat  = True |>>>>>
<<<<<| Plot 2D lon X time = True |>>>>>
<<<<<| Plot 2D lon X lev  = True |>>>>>
<<<<<| Plot 2D lat X lev  = True |>>>>>
<<<<<| Plot 2D time X lat = True |>>>>>
<<<<<| Plot 2D time X lev = True |>>>>>
<<<<<| Plot 1D            = True |>>>>> # Any 1D Plot Type (Dimension x Variable)
```

## Adjust the color range  and colormap

`Cmin, Cmax` (and `Contours Var 2`) are how the contours are set for the shaded (and solid) contours. If only two values are included, MarsPlot use 24 contours spaced between the max and min values. If more than two values are provided, MarsPlot will use those individual contours.

```python
Main Variable  = atmos_average.temp     # filename.variable *REQUIRED
Cmin, Cmax     = 240,290                # Colorbar limits (minimum, maximum)
2nd Variable   = atmos_average.ucomp    # Overplot U winds
Contours Var 2 = -200,-100,100,200      # List of contours for 2nd Variable or CMIN, CMAX
Axis Options  : Ls = [None,None] | lat = [None,None] | cmap = jet |scale = lin
```

Note the option of setting the contour spacing linearly `scale = lin` or logarithmically (`scale = log`) if the range of values spans multiple order of magnitudes.

The default colormap `cmap = jet` may be changed using any Matplotlib colormaps. A selections of those are listed below:

![Figure 4. MarsPlot workflow](./tutorial_images/all_colormaps.png)

Finally, note the use of the `_r` suffix (reverse) to reverse the order of the colormaps listed in the figure above. From example, using `cmap = jet` would have colors spanning from *blue* > *red*  and `cmap = jet_r` *red* > *blue* instead

*Supported colormap in Marsplot. The figure was generated using code from [the scipy webpage](https://scipy-lectures.org/intro/matplotlib/auto_examples/options/plot_colormaps.html) .*

## Customize Plots
`Axis Options` specify the axes limits, colormap, linestyle and color for 1D-plots, projection for certain plots :

```python
# Axis Options for 2D plots may include:
Lat         = [0,90]        # Latitude range for axes limits
Level[Pa/m] = [600,10]      # Level range for axes limits
sols        = [None,None]   # Sol range for axes limits
Lon +/-180  = [-180,180]    # Longitude range for axes limits
cmap        = jet           # Python colormap to use
scale       = lin           # Color map style ([lin]ear, [log]arithmic)
proj        = cart          # Projection ([cart]esian, [robin]son, [moll]weide, [Npole], [Spole], [ortho]graphic)
# Axis Options for 1D plots may include:

 lat,lon+/-180,[Pa/m],sols = [None,None] # range for X or Y axes limit
 var = [None,None]                       # range for displayed variables
 linestyle = -                           # Line style following matplotlib convention'-r' (solid red), '--g' (dashed green), '-ob' (solid & blue markers)
 axlabel = None                          # Change the default name for the axis
```




***
## Make a 1D-plot
The 1D plot template is different from the others in a few key ways:

- Instead of `Title`, the template requires a `Legend`. When overploting several 1D variables on top of one another, the legend option will label them insetad of changing the plot title.
- There is an additional `linestyle` axis option for the 1D plot.
- There is also a `Diurnal` option. The `Diurnal` input can only be `None` or `AXIS`, since there is syntax for selecting a specific time of day. The `AXIS` label tells `MarsPlot` which dimension serves as the X axis. `Main Variable` will dictate the Y axis.

```python
<<<<<<<<<<<<<<| Plot 1D = True |>>>>>>>>>>>>>
Legend         = None                   # Legend instead of Title
Main Variable  = atmos_average.temp
Ls 0-360       = AXIS                   #       Any of these can be selected
Latitude       = None                   #       as the X axis dimension, and
Lon +/-180     = None                   #       the free dimensions can accept
Level [Pa/m]   = None                   #       values as before. However,
Diurnal  [hr]  = None                   #   ** Diurnal can ONLY be AXIS or None **
```




***

There are several other plot customizations you can use:

* When two or more blocks are sandwiched between a `HOLD ON` and `HOLD OFF`, `MarsPlot` will draw the plots on the same page.
* Plots are created on a standard page (8.5 x 11 inches) in landscape mode, but can be drawn in portrait mode as well.
* Plots can be saved as images instead of PDFs by specifying your preferred filetype (PNG, EPS, etc.) when passing the `--output` (`-o`) argument to `MarsPlot`.
* When creating 1D plots of data spanning multiple years, you can overplot consecutive years by calling `--stack_year` (`-sy`) when submitting the template to `MarsPlot`.
* Specify which MGCM output file to use when plotting by passing the `--date` (`-d`) argument to `MarsPlot` followed by the 5-digit file prefix corresponding to the file you want to use. Alternatively, add the prefix to the filename in the template (e.g. `Main Variable = 00000.fixed.thin`).



***
## Access simulation in a different directory
The final plot-related functionality in `MarsPlot` is the simulation list, which allows you to point `MarsPlot` to different directories containing the MGCM output:

```python
<<<<<<<<<<<<<<<<<<<<<< Simulations >>>>>>>>>>>>>>>>>>>>>
ref> None
2> /path/to/another/sim                            # another simulation
3>
=======================================================
```

To access a variable from a file in another directory, just point to the correct simulation when calling `Main Variable`:

```python
Main Variable  = XXXXX.filename@N.variable`
```

Where `N` is the number in `<<< Simulations >>>` corresponding the the correct path.

## Element-wise operations
The `Main Variable` input also accepts variable operations and time-of-day selections like so:

```python
Main Variable  = [filename.variable]*1000  # multiply all values by 1000
Main Variable  = filename.variable{tod = 20}  # select the 20th hour of the day
```

At minimum, `Main Variable` requires `filename.variable` for input, but the above syntax can be combined in several ways allowing for greater plot customization. For example, to plot dust mixing ratio from the diurnal file in simulation #3 at 3 PM local time, the input is:

```python
Main Variable  = [atmos_diurn_plevs_T@2.dst_mass_micro{tod = 15}]*1.e6 # dust ppm
#                [filename@N.variable{dimension = X}]*Y
```

## Change projections


***
## Debugging
`MarsPlot` is designed to make plotting MGCM output easier and faster so it handles missing data for you. For example, when dimensions are omitted with `None`, `MarsPlot` makes educated guesses for data selection and will tell you exactly how the data is being processed both in the title for the figures (if `Title = None`), and in the terminal output. Specifics about this behavior are detailed in the instructions at the top of `Custom.in`.

> `MarsPlot` handles many errors by itself. It reports errors both in the terminal and in the generated figures. To by-pass this behavior (when debugging), use the  `--debug` option with `MarsPlot` which will raise standard Python errors.





[Back to Top](#cheat-sheet)
***
