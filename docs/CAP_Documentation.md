![](./tutorial_images/Tutorial_Banner_Final.png)

<!-- TOC titleSize:2 tabSpaces:2 depthFrom:1 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 skip:0 title:1 charForUnorderedList:* -->
## Table of Contents
* [1. `MarsPull` - Downloading Raw MGCM Output](#1-marspullpy---downloading-raw-mgcm-output)
* [2. `MarsFiles` - Reducing the Files](#2-marsfilespy---reducing-the-files)
* [3. `MarsVars` - Performing Variable Operations](#3-marsvarspy---performing-variable-operations)
* [4. `MarsInterp` - Interpolating the Vertical Grid](#4-marsinterppy---interpolating-the-vertical-grid)
* [5. `MarsPlot` - Plotting the Results](#5-marsplotpy---plotting-the-results)
<!-- /TOC -->

***


# 1. `MarsPull` - Downloading Raw MGCM Output

`MarsPull` is a utility for accessing MGCM output files hosted on the [MCMC Data portal](https://data.nas.nasa.gov/legacygcm/data_legacygcm.php). MGCM data is archived in 1.5 hour intervals (16x/day) and packaged in files containing 10 sols. The files are named fort.11_XXXX in the order they were produced, but  `MarsPull` maps those files to specific solar longitudes (L<sub>s</sub>, in °). This allows users to request a file at a specific L<sub>s</sub> or for a range of L<sub>s</sub> using the `-ls` flag. Additionally the `identifier` (`-id`) flag is used to route `MarsPull` through a particular simulation. The `filename` (`-f`) flag can be used to parse specific files within a particular directory.

```bash
MarsPull -id INERTCLDS -ls 255 285
MarsPull -id ACTIVECLDS -f fort.11_0720 fort.11_0723
```
[Back to Top](#cheat-sheet)
***

# 2. `MarsFiles` - Reducing the Files

`MarsFiles` provides several tools for file manipulations, including code designed to create binned, averaged, and time-shifted files from MGCM output. The `-fv3` flag is used to convert fort.11 binaries to the Netcdf data format (you can select one or more of the file format listed below):

```bash
(AmesCAP)>$ MarsFiles fort.11* -fv3 fixed average daily diurn
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

***

# 3. `MarsVars` - Performing Variable Operations

`MarsVars` provides several tools relating to variable operations such as adding and removing variables, and performing column integrations. With no other arguments, passing a file to `MarsVars` displays file content, much like `ncdump`:

```bash
(AmesCAP)>$ MarsVars 00000.atmos_average.nc
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
(AmesCAP)>$ MarsVars 00000.atmos_average.nc -add rho
```

We can see that `rho` was added by calling `MarsVars` with no argument as before:

```bash
(AmesCAP)>$ MarsVars 00000.atmos_average.nc
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

# 4. `MarsInterp` - Interpolating the Vertical Grid

Native MGCM output files use a terrain-following pressure coordinate as the vertical coordinate (`pfull`), which means the geometric heights and the actual mid-layer pressure of atmospheric layers vary based on the location (i.e. between adjacent grid points). In order to do any rigorous spatial averaging, it is therefore necessary to interpolate each vertical column to a same (standard) pressure grid (`_pstd` grid):

![Figure X. MarsInterp](./tutorial_images/MarsInterp.png)

*Pressure interpolation from the reference pressure grid to a standard pressure grid*

`MarsInterp` is used to perform the vertical interpolation from *reference* (`pfull`) layers to *standard* (`pstd`) layers:

```bash
(AmesCAP)>$ MarsInterp  00000.atmos_average.nc
```

An inspection of the file shows that the pressure level axis which was `pfull` (30 layers) has been replaced by a standard pressure coordinate `pstd` (36 layers), and all 3- and 4-dimensional variables reflect the new shape:

```bash
(AmesCAP)>$ MarsInterp  00000.atmos_average.nc
(AmesCAP)>$ MarsVars 00000.atmos_average_pstd.nc
>
> ===================DIMENSIONS==========================
> ['bnds', 'time', 'lat', 'lon', 'scalar_axis', 'phalf', 'pstd']
> ====================CONTENT==========================
> pstd           : ('pstd',)= (36,), pressure  [Pa]
> temp           : ('time', 'pstd', 'lat', 'lon')= (4, 36, 180, 360), temperature  [K]
```

`MarsInterp` support 3 types of vertical interpolation, which may be selected by using the `--type` (`-t` for short) flag:

| file type | description | low-level value in a deep crater
|-----------|-----------|--------|
|_pstd | standard pressure [Pa] (default) |  1000Pa
|_zstd | standard altitude [m]   |  -7000m
|_zagl | standard altitude above ground level [m]   | 0 m

***

**Use of custom vertical grids**

`MarsInterp` uses default grids for each of the interpolation listed above but it is possible for the user to specify the layers for the interpolation. This is done by editing a **hidden** file `.amescap_profile`(note the dot '`.`) in your home directory.  

For the first use, you will need to copy a template of `amescap_profile` to your /home directory:

```bash
(AmesCAP)>$ cp ~/AmesCAP/mars_templates/amescap_profile ~/.amescap_profile # Note the dot '.' !!!
```
You can open `~/.amescap_profile` with any text editor:

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

You can use these by calling `MarsInterp` with the `-level` (`-l`) argument followed by the name of the new grid defined in `.amescap_profile`.

```bash
(AmesCAP)>$ MarsInterp  00000.atmos_average.nc -t pstd -l  p44
```
[Back to Top](#cheat-sheet)
***

# 5. `MarsPlot` - Plotting the Results


[Back to Top](#cheat-sheet)
***
