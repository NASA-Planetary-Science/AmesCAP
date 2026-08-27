#!/usr/bin/env python3
"""
The MarsFormat executable is for converting non-MGCM data, such as that
from EMARS, OpenMARS, PCM, and MarsWRF, into MGCM-like netCDF data
products. The MGCM is the NASA Ames Mars Global Climate Model developed
and maintained by the Mars Climate Modeling Center (MCMC). The MGCM
data repository is available at data.nas.nasa.gov/mcmc.

The executable requires two arguments:

    * ``[input_file]``         The file to be transformed
    * ``[-gcm --gcm_name]``    The GCM from which the file originates

and optionally accepts:

    * ``[-rn --retain_names]`` Preserve original variable and dimension names
    * ``[-ba, --bin_average]`` Bin non-MGCM files like 'average' files
    * ``[-bd, --bin_diurn]``   Bin non-MGCM files like 'diurn' files

Third-party requirements:

    * ``numpy``
    * ``netCDF4``
    * ``sys``
    * ``argparse``
    * ``os``
    * ``re``
    * ``functools``
    * ``traceback``
    * ``xarray``
    * ``amescap``
"""

# Make print statements appear in color
from amescap.Script_utils import (
    Green, Yellow, Red, Blue, Cyan, Purple, Nclr,
)

# Load generic Python modules
import sys          # System commands
import argparse     # Parse arguments
import os           # Access operating system functions
import re           # Regular expressions
import numpy as np
import xarray as xr
from netCDF4 import Dataset
import functools    # For function decorators
import traceback    # For printing stack traces

# Load amesCAP modules
from amescap.Script_utils import (
    read_variable_dict_amescap_profile, reset_FV3_names
)
from amescap.FV3_utils import layers_mid_point_to_boundary

xr.set_options(keep_attrs=True)


def debug_wrapper(func):
    """
    A decorator that wraps a function with error handling
    based on the --debug flag.
    If the --debug flag is set, it prints the full traceback
    of any exception that occurs. Otherwise, it prints a
    simplified error message.

    :param func: The function to wrap.
    :type   func: function
    :return: The wrapped function.
    :rtype:  function
    :raises Exception: If an error occurs during the function call.
    :raises TypeError: If the function is not callable.
    :raises ValueError: If the function is not found.
    :raises NameError: If the function is not defined.
    :raises AttributeError: If the function does not have the
        specified attribute.
    :raises ImportError: If the function cannot be imported.
    :raises RuntimeError: If the function cannot be run.
    :raises KeyError: If the function does not have the
        specified key.
    :raises IndexError: If the function does not have the
        specified index.
    :raises IOError: If the function cannot be opened.
    :raises OSError: If the function cannot be accessed.
    :raises EOFError: If the function cannot be read.
    :raises MemoryError: If the function cannot be allocated.
    :raises OverflowError: If the function cannot be overflowed.
    :raises ZeroDivisionError: If the function cannot be divided by zero.
    :raises StopIteration: If the function cannot be stopped.
    :raises KeyboardInterrupt: If the function cannot be interrupted.
    :raises SystemExit: If the function cannot be exited.
    :raises AssertionError: If the function cannot be asserted.
    """

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        global debug
        try:
            return func(*args, **kwargs)
        except Exception as e:
            if debug:
                # In debug mode, show the full traceback
                print(f"{Red}ERROR in {func.__name__}: {str(e)}{Nclr}")
                traceback.print_exc()
            else:
                # In normal mode, show a clean error message
                print(f"{Red}ERROR in {func.__name__}: {str(e)}\nUse "
                      f"--debug for more information.{Nclr}")
            return 1  # Error exit code
    return wrapper


# ======================================================
#                  ARGUMENT PARSER
# ======================================================

parser=argparse.ArgumentParser(
    prog=('MarsFormat'),
    description=(
        f"{Yellow} Converts model output to MGCM-like format for "
        f"compatibility with CAP."
        f"{Nclr}\n\n"
    ),
    formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('input_file', nargs='+',
    type=argparse.FileType('rb'),
    help=(f"A netCDF file or list of netCDF files.\n\n"))

parser.add_argument('-gcm', '--gcm_name', type=str,
    choices=['marswrf', 'openmars', 'pcm', 'emars'],
    help=(
        f"Acceptable types include 'openmars', 'marswrf', 'emars', "
        f"and 'pcm' \n"
        f"{Green}Example:\n"
        f"> MarsFormat openmars_file.nc -gcm openmars\n"
        f"{Blue}(Creates openmars_file_daily.nc)"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-ba', '--bin_average', nargs="?", const=5,
    type=int,
    help=(
        f"Calculate 5-day averages from instantaneous data. Generates "
        f"MGCM-like 'average' files.\n"
        f"{Green}Example:\n"
        f"> MarsFormat openmars_file.nc -gcm openmars -ba\n"
        f"{Blue}(Creates openmars_file_average.nc; 5-sol bin){Green}\n"
        f"> MarsFormat openmars_file.nc -gcm openmars -ba 10\n"
        f"{Blue}(10-sol bin)"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-bd', '--bin_diurn', action='store_true',
    default=False,
    help=(
        f"Calculate 5-day averages binned by hour from instantaneous "
        f"data. Generates MGCM-like 'diurn' files.\n"
        f"Works on non-MGCM files only.\n"
        f"{Green}Example:\n"
        f"> MarsFormat openmars_file.nc -gcm openmars -bd\n"
        f"{Blue}(Creates openmars_file_diurn.nc;  5-sol bin){Green}\n"
        f"> MarsFormat openmars_file.nc -gcm openmars -bd -ba 10\n"
        f"{Blue}(10-sol bin)"
        f"{Nclr}\n\n"
    )
)

# Secondary arguments: Used with some of the arguments above

parser.add_argument('-rn', '--retain_names', action='store_true',
    default=False,
    help=(
        f"Preserves the names of the variables and dimensions in the"
        f"original file.\n"
        f"{Green}Example:\n"
        f"> MarsFormat openmars_file.nc -gcm openmars -rn\n"
        f"{Blue}(Creates openmars_file_nat_daily.nc)"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-stag', '--keep_staggered', action='store_true',
    default=False,
    help=(
        f"[marswrf only] Also keep the winds on their native Arakawa-C\n"
        f"locations: u_stag(time, pfull, lat, lon_u),\n"
        f"v_stag(time, pfull, lat_v, lon), w_half(time, phalf, lat, lon),\n"
        f"with lon_u, lat_v axes and the interface pressure phalf3D\n"
        f"(PF_PHY) and height zhalf (ZF_PHY). Mass-point ucomp, vcomp, w\n"
        f"are written regardless. u_stag, v_stag are for your own\n"
        f"scripts; CAP tools only plot the mass-point grid.\n"
        f"{Green}Example:\n"
        f"> MarsFormat wrfout_d01.nc -gcm marswrf -stag"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('--debug', action='store_true',
    help=(
        f"Use with any other argument to pass all Python errors and\n"
        f"status messages to the screen when running CAP.\n"
        f"{Green}Example:\n"
        f"> MarsFormat openmars_file.nc -gcm openmars --debug"
        f"{Nclr}\n\n"
    )
 )

args = parser.parse_args()
debug = args.debug

if args.input_file:
    for file in args.input_file:
        if not re.search(".nc", file.name):
            parser.error(f"{Red}{file.name} is not a netCDF file{Nclr}")
            exit()

# ----------------------------------------------------------------------
path2data = os.getcwd()
ref_press = 725 # TODO hard-coded reference pressure


def get_time_dimension_name(DS, model):
    """
    Find the time dimension name in the dataset.

    Updates the model object with the correct dimension name.

    :param DS: The xarray Dataset
    :type  DS: xarray.Dataset
    :param model: Model object with dimension information
    :type  model: object
    :return: The actual time dimension name found
    :rtype:  str
    :raises KeyError: If no time dimension is found
    :raises ValueError: If the model object is not defined
    :raises TypeError: If the dataset is not an xarray Dataset
    :raises AttributeError: If the model object does not have the
        specified attribute
    :raises ImportError: If the xarray module cannot be imported
    """

    # First try the expected dimension name
    if model.dim_time in DS.dims:
        return model.dim_time

    # Check alternative names
    possible_names = ['Time', 'time', 'ALSO_Time']
    for name in possible_names:
        if name in DS.dims:
            print(f"{Yellow}Notice: Using '{name}' as time dimension "
                  f"instead of '{model.dim_time}'{Nclr}")
            model.dim_time = name
            return name

    # If no time dimension is found, raise an error
    raise KeyError(f"No time dimension found in dataset. Expected one "
                   f"of: {model.dim_time}, {', '.join(possible_names)}")

# ======================================================================
#                      MarsWRF / planetWRF helpers
# ======================================================================
def _wrf_datestr_to_sol(datestr):
    """
    Convert a planetWRF date string 'YYYY-DDDDD_HH:MM:SS' to a sol
    count. DDDDD is the sol of year (1-based in the WRF convention).
    Returns float sols with year 1, sol 1, 00:00 mapping to 0.0.
    """
    datestr = str(datestr).strip()
    m = re.match(r'(\d+)-(\d+)_(\d+):(\d+):(\d+)', datestr)
    if not m:
        raise ValueError(f"Cannot parse WRF date string '{datestr}'")
    yr, sol, hh, mm, ss = (int(x) for x in m.groups())
    return (yr - 1)*669. + (sol - 1) + (hh + (mm + ss/60.)/60.)/24.


def marswrf_time_axis(DS, model):
    """
    Build the time coordinate [sols] for a planetWRF file from the
    best available clock and return (DS, source_description).

    Priority:
      1. 'Times' character array 'YYYY-DDDDD_HH:MM:SS' (absolute,
         carries the time of sol)
      2. 'XTIME' [minutes since SIMULATION_START_DATE] plus the sol of
         SIMULATION_START_DATE
      3. 'JULIAN' (fractional sol of year) plus MODEL_MARS_YEAR
      4. START_DATE global attribute plus a uniform spacing inferred
         from the L_S slope (reduced files that kept only L_S)

    The result is stored under model.time (normally 'XTIME') so that
    the rest of MarsFormat sees one time variable with units in days.
    Ls is left in model.areo ('L_S'); if it is missing it is computed
    from the sol axis with sol2ls().
    """
    tdim = model.dim_time
    ntime = DS.sizes[tdim]
    src = None
    if 'Times' in DS:
        try:
            raw = DS['Times'].values
            strs = []
            for row in raw:
                if isinstance(row, (bytes, str)):
                    strs.append(row.decode() if isinstance(row, bytes)
                                else row)
                else:
                    strs.append(b''.join(row).decode())
            sols = np.array([_wrf_datestr_to_sol(x) for x in strs])
            src = "'Times' date strings"
        except Exception as e:
            print(f"{Yellow}Could not parse 'Times' ({e}); trying XTIME")
            src = None
    if src is None and 'XTIME' in DS:
        sols = np.asarray(DS['XTIME'].values, dtype=float)/1440.
        start = DS.attrs.get('SIMULATION_START_DATE', None)
        if start is not None:
            try:
                sols = sols + _wrf_datestr_to_sol(start)
                src = "'XTIME' + SIMULATION_START_DATE"
            except ValueError:
                src = "'XTIME' (minutes since simulation start)"
        else:
            src = "'XTIME' (minutes since simulation start)"
    if src is None and 'JULIAN' in DS:
        sols = np.asarray(DS['JULIAN'].values, dtype=float)
        if 'MODEL_MARS_YEAR' in DS:
            sols = sols + 669.*(np.asarray(DS['MODEL_MARS_YEAR'].values,
                                           dtype=float) - 1)
        src = "'JULIAN' (+ MODEL_MARS_YEAR)"
    if src is None:
        # Reduced file: nothing but L_S. Assume uniform output
        # interval starting at START_DATE, spacing from the Ls slope
        # (Ls advances ~0.538 deg/sol on average).
        sol0 = 0.
        start = DS.attrs.get('START_DATE', None)
        if start is not None:
            try:
                sol0 = _wrf_datestr_to_sol(start)
            except ValueError:
                pass
        if model.areo in DS and ntime > 1:
            ls = np.asarray(DS[model.areo].values, dtype=float).ravel()
            dls = np.diff(np.unwrap(np.deg2rad(ls)))
            dsol = np.rad2deg(np.median(dls))/(360./669.)
            # snap to the nearest 1/24 sol
            dsol = np.round(dsol*24.)/24.
        else:
            dsol = 1.
        sols = sol0 + dsol*np.arange(ntime)
        src = (f"START_DATE + uniform {dsol:.4f} sol spacing inferred "
               f"from L_S (no Times/XTIME/JULIAN in file)")

    DS[model.time] = xr.DataArray(sols, dims=(tdim,))
    DS[model.time].attrs['long_name'] = 'time'
    DS[model.time].attrs['units'] = 'days since 0000-00-00 00:00:00'
    DS[model.time].attrs['description'] = (
        f'(ADDED POST-PROCESSING) sols, from {src}')
    print(f"{Cyan}Time axis [sols] built from {src}{Nclr}")

    if model.areo not in DS:
        from amescap.FV3_utils import sol2ls
        DS[model.areo] = xr.DataArray(sol2ls(sols), dims=(tdim,))
        DS[model.areo].attrs['units'] = 'degrees'
        DS[model.areo].attrs['description'] = (
            '(ADDED POST-PROCESSING) Ls from sol2ls(time)')
        print(f"{Yellow}No {model.areo} in file: Ls computed from the "
              f"sol axis with sol2ls(){Nclr}")
    return DS, src


def marswrf_fit_eta(p3d, ps, tdim_axis=0):
    """
    Recover the WRF eta coordinate from a 3D mid-level pressure field
    and surface pressure, for reduced files that did not keep
    ZNU/ZNW/P_TOP. WRF: p = p_top + eta*(ps - p_top), so per level a
    least-squares fit p = a_k + b_k*ps over all columns and times
    gives eta_k = b_k and p_top = a_k/(1 - b_k).

    p3d: array [time, lev, lat, lon] (any vertical order)
    ps : array [time, lat, lon]
    Returns (eta_mid, p_top).
    """
    nlev = p3d.shape[1]
    x = np.asarray(ps, dtype=float).ravel()
    A = np.vstack([np.ones_like(x), x]).T
    a = np.zeros(nlev)
    b = np.zeros(nlev)
    for k in range(nlev):
        y = np.asarray(p3d[:, k], dtype=float).ravel()
        coef, *_ = np.linalg.lstsq(A, y, rcond=None)
        a[k], b[k] = coef
    good = (1. - b) > 0.05
    p_top = float(np.median(a[good]/(1. - b[good]))) if good.any() else 0.
    p_top = max(p_top, 0.)
    return b, p_top


def marswrf_eta_mid_to_interfaces(eta_mid):
    """
    Interface (w-level) eta from mid (mass-level) eta, surface first:
    eta_w[0] = 1 and eta_w[k+1] = 2*eta_u[k] - eta_w[k]. The last
    value is clipped to >= 0. Input may be in either vertical order;
    output follows the input order.
    """
    eta = np.asarray(eta_mid, dtype=float)
    flip = eta[0] < eta[-1]      # top first
    if flip:
        eta = eta[::-1]
    w = np.zeros(len(eta) + 1)
    w[0] = 1.
    for k in range(len(eta)):
        w[k+1] = 2.*eta[k] - w[k]
    w[-1] = max(w[-1], 0.)
    return w[::-1] if flip else w


@debug_wrapper
def main():
    """
    Main processing function for MarsFormat.

    This function processes NetCDF files from various Mars General
    Circulation Models (GCMs)
    including MarsWRF, OpenMars, PCM, and EMARS, and reformats them for
    use in the AmesCAP
    framework.

    It performs the following operations:
        - Validates the selected GCM type and input files.
        - Loads NetCDF files and reads model-specific variable and
        dimension mappings.
        - Applies model-specific post-processing, including:
            - Unstaggering variables (for MarsWRF and EMARS).
            - Creating and orienting pressure coordinates (pfull, phalf,
            ak, bk).
            - Standardizing variable and dimension names.
            - Converting longitude ranges to 0-360 degrees east.
            - Adding scalar axes where required.
            - Handling vertical dimension orientation, especially for
            PCM files.
        - Optionally performs time binning:
            - Daily, average (over N sols), or diurnal binning.
            - Ensures correct time units and bin sizes.
            - Preserves or corrects vertical orientation after binning.
        - Writes processed datasets to new NetCDF files with appropriate
        naming conventions.

    Args:
        None. Uses global `args` for configuration and file selection.

    Raises:
        KeyError: If required dimensions or variables are missing in
        the input files.
        ValueError: If dimension swapping fails for PCM files.
        SystemExit: If no valid GCM type is specified.

    Outputs:
        Writes processed NetCDF files to disk, with suffixes indicating
        the type of processing
        (e.g., _daily, _average, _diurn, _nat).

    Note:
        This function assumes the presence of several helper functions
        and global variables,
        such as `read_variable_dict_amescap_profile`,
        `get_time_dimension_name`,  `reset_FV3_names`, and color
        constants for printing.
        """

    ext = '' # Initialize empty extension

    if args.gcm_name not in ['marswrf', 'openmars', 'pcm', 'emars']:
        print(f"{Yellow}***Notice***  No operation requested. Use "
              f"'-gcm' and specify openmars, marswrf, pcm, emars")
        exit() # Exit cleanly

    print(f"Running MarsFormat with args: {args}")
    print(f"Current working directory: {os.getcwd()}")
    print(f"Files in input_file: {[f.name for f in args.input_file]}")
    print(f"File exists check: "
          f"{all(os.path.exists(f.name) for f in args.input_file)}")

    path2data = os.getcwd()

    # Load all of the netcdf files
    file_list = [f.name for f in args.input_file]
    model_type = args.gcm_name  # e.g. 'marswrf'
    for filei in file_list:
        # Use os.path.join for platform-independent path handling
        if os.path.isabs(filei):
            fullnameIN = filei
        else:
            fullnameIN = os.path.join(path2data, filei)

        print('Processing...')
        # Load model variables, dimensions
        fNcdf = Dataset(fullnameIN, 'r')
        model = read_variable_dict_amescap_profile(fNcdf)
        fNcdf.close()

        print(f"{Cyan}Reading model attributes from ~.amescap_profile:")
        print(f"{Cyan}{vars(model)}") # Print attributes

        # Open dataset with xarray
        DS = xr.open_dataset(fullnameIN, decode_times=False)

        if model_type == 'marswrf':
            # Build the sol axis from the best clock in the file
            # (Times > XTIME > JULIAN > START_DATE+L_S) before anything
            # reads model.time; planetWRF writes XTIME with units=''.
            DS, _ = marswrf_time_axis(DS, model)

        # Store the original time values and units before any modifications
        original_time_vals = DS[model.time].values.copy()  # This will always exist
        original_time_units = DS[model.time].attrs.get('units', '')
        original_time_desc = DS[model.time].attrs.get('description', '')
        print(
            f"DEBUG: Saved original time values with units "
            f"'{original_time_units}' and description "
            f"'{original_time_desc}'"
            )

        # Find and update time dimension name
        time_dim = get_time_dimension_name(DS, model)

        # --------------------------------------------------------------
        #                      MarsWRF Processing
        # --------------------------------------------------------------
        if model_type == 'marswrf':
            # print(f"{Cyan}Current variables at top of marswrf "
            #       f"processing: \n{list(DS.variables)}{Nclr}\n")

            # First, save all variable descriptions in attrs longname
            for var_name in DS.data_vars:
                var = DS[var_name]
                if 'description' in var.attrs:
                    var.attrs['long_name'] = var.attrs['description']

            # Ensure mandatory dimensions and coordinates exist
            # Reformat Dimension Variables/Coords as Needed
            # Handle potential missing dimensions or coordinates
            if model.time not in DS:
                raise KeyError(f"Time dimension {model.time} not found")
            if model.lat not in DS:
                raise KeyError(f"Latitude dimension {model.lat} not found")
            if model.lon not in DS:
                raise KeyError(
                    f"Longitude dimension {model.lon} not found"
                    )

            # Time axis was built by marswrf_time_axis() at file open.

            # Handle latitude and longitude (2D XLAT/XLONG on a regular
            # lat-lon grid; take one row/column)
            if len(DS[model.lat].shape) > 1:
                lat = DS[model.lat][0, :, 0].values
                lon = DS[model.lon][0, 0, :].values
            else:
                lat = DS[model.lat].values
                lon = DS[model.lon].values

            # Convert longitudes to 0-360 and SORT so the coordinate
            # is monotonic (a bare modulo leaves e.g. 185..355,5..175).
            lon360 = (lon + 360) % 360
            DS = DS.assign_coords({model.dim_lon: lon360,
                                   model.dim_lat: lat})
            DS = DS.sortby(model.dim_lon)
            DS[model.lon] = xr.DataArray(DS[model.dim_lon].values,
                                         dims=(model.dim_lon,),
                                         attrs={'units': 'degrees_E',
                                                'long_name': 'longitude'})
            DS[model.lat] = xr.DataArray(lat, dims=(model.dim_lat,),
                                         attrs={'units': 'degrees_N',
                                                'long_name': 'latitude'})

            # Vertical coordinate: WRF eta (terrain-following mass)
            #   p = P_TOP + eta*(ps - P_TOP)  ->  ak = P_TOP*(1-eta),
            #   bk = eta, with eta on w-levels (ZNW) for interfaces and
            #   on mass levels (ZNU) for layer centres.
            if 'ZNW' in DS and 'P_TOP' in DS:
                p_top = float(np.asarray(DS.P_TOP.values).ravel()[0])
                eta_w = np.array(DS.ZNW.values[0, :], copy=True)
                eta_u = np.array(DS.ZNU.values[0, :], copy=True)
                print(f"{Cyan}Vertical coordinate from ZNU/ZNW, "
                      f"P_TOP = {p_top:.3e} Pa{Nclr}")
            elif 'P_PHY' in DS and 'PSFC' in DS:
                # Reduced file: recover eta by fitting P_PHY to PSFC
                pp = DS.P_PHY.transpose(model.dim_time, model.dim_pfull,
                                        ...).values
                okt = (pp > 0).reshape(pp.shape[0], -1).all(axis=1)
                pp = pp[okt]
                ps = DS.PSFC.transpose(model.dim_time, ...).values[okt]
                eta_u, p_top = marswrf_fit_eta(pp, ps)
                eta_w = marswrf_eta_mid_to_interfaces(eta_u)
                print(f"{Yellow}No ZNU/ZNW/P_TOP in file: eta fitted "
                      f"from P_PHY vs PSFC over {int(okt.sum())} "
                      f"frames, P_TOP = {p_top:.3e} Pa{Nclr}")
            else:
                raise KeyError("Cannot build the vertical coordinate: "
                               "need ZNU/ZNW/P_TOP or P_PHY/PSFC")
            P0 = float(DS.attrs.get('P0', 610.))
            phalf = p_top + eta_w*P0
            pfull = p_top + eta_u*P0

            DS = DS.assign_coords(pfull=(model.dim_pfull, pfull))
            DS = DS.assign_coords(phalf=(model.dim_phalf, phalf))

            ak = p_top*(1. - eta_w)
            bk = eta_w.copy()

            # Fill ak, bk, pfull, phalf arrays
            DS = DS.assign(ak=(model.dim_phalf, ak))
            DS = DS.assign(bk=(model.dim_phalf, bk))

            DS.phalf.attrs['description'] = (
                '(ADDED POST-PROCESSING) pressure at layer interfaces')
            DS.phalf.attrs['units'] = ('Pa')

            DS['ak'].attrs['description'] = (
               '(ADDED POST-PROCESSING)  pressure part of the hybrid coordinate')
            DS['bk'].attrs['description'] = (
            '(ADDED POST-PROCESSING) vertical coordinate sigma value')
            DS['ak'].attrs['units']='Pa'
            DS['bk'].attrs['units']='None'

            have_geopot = all(v in DS for v in ('PH', 'PHB', 'HGT'))
            if have_geopot:
                # Height of w-levels above the local surface [m]
                zagl_lvl = ((DS.PH + DS.PHB)/DS.G - DS.HGT)
            else:
                zagl_lvl = None

            # Derive full 3D pressure [Pa] and temperature [K]
            # ------------------------------------------------
            # WRF writes T as the PERTURBATION potential temperature
            # (theta = T + T0) and P as the PERTURBATION pressure
            # (p = P + PB, where PB is the base-state pressure, which
            # already includes P_TOP). Real temperature follows from
            # the Exner function evaluated with the full pressure:
            #     temp = (T + T0) * ((P + PB) / P0) ** (R_D / CP)
            # The previous form, (P_TOP + PB[0]), double-counted P_TOP
            # and dropped the perturbation pressure, giving errors of
            # several K (up to ~50 K near topography).
            #
            # Newer planetWRF files also carry the physics-grid
            # diagnostics T_PHY and P_PHY (temperature and pressure on
            # mass levels). When present and filled they are used
            # directly. They are zero on the first frame of a file
            # written straight after a restart, so fall back to the
            # reconstruction for any frame where P_PHY is not positive.
            kappa = DS.R_D / DS.CP
            have_phy = ('T_PHY' in DS) and ('P_PHY' in DS)
            if 'P' in DS and 'PB' in DS:
                pfull3D_rec = DS.P + DS.PB
            else:
                pfull3D_rec = None
            if 'T' in DS and pfull3D_rec is not None:
                temp_rec = (DS.T + DS.T0) * (pfull3D_rec / DS.P0)**kappa
            else:
                temp_rec = None

            if have_phy:
                phy_ok = (DS.P_PHY > 0).all(dim=[d for d in DS.P_PHY.dims
                                                  if d != model.dim_time])
                if temp_rec is not None:
                    temp = DS.T_PHY.where(phy_ok, temp_rec)
                    pfull3D = DS.P_PHY.where(phy_ok, pfull3D_rec)
                    n_bad = int((~phy_ok).sum())
                    if n_bad:
                        print(f"{Yellow}T_PHY/P_PHY unfilled on {n_bad} "
                              f"frame(s); reconstructed from T, P, PB "
                              f"there{Nclr}")
                else:
                    temp = DS.T_PHY
                    pfull3D = DS.P_PHY
                src = 'T_PHY, P_PHY'
            elif temp_rec is not None:
                temp = temp_rec
                pfull3D = pfull3D_rec
                src = '(T+T0)*((P+PB)/P0)**(R_D/CP)'
            else:
                raise KeyError("Cannot derive temperature: need T_PHY and "
                               "P_PHY, or T, P and PB")
            print(f"{Cyan}Temperature derived from {src}{Nclr}")

            DS = DS.assign(temp=temp)
            DS['temp'].attrs = {}
            DS['temp'].attrs['description'] = ('(ADDED POST-PROCESSING) Temperature')
            DS['temp'].attrs['long_name'] = ('(ADDED POST-PROCESSING) Temperature')
            DS['temp'].attrs['units'] = 'K'

            DS = DS.assign(pfull3D=pfull3D)
            DS['pfull3D'].attrs = {}
            DS['pfull3D'].attrs['description'] = ('(ADDED POST-PROCESSING) Pressure')
            DS['pfull3D'].attrs['long_name'] = ('(ADDED POST-PROCESSING) Pressure')
            DS['pfull3D'].attrs['units'] = 'Pa'

            # Optionally keep the native Arakawa-C winds before the
            # destaggering below overwrites U, V, W in place.
            #   u_stag: x-staggered, lon_u axis from XLONG_U
            #   v_stag: y-staggered, lat_v axis from XLAT_V
            #   w_half: on layer interfaces -> CAP's phalf dimension,
            #           with phalf3D (PF_PHY) and zhalf (ZF_PHY)
            if getattr(args, 'keep_staggered', False):
                kept = []
                if 'U' in DS and 'XLONG_U' in DS:
                    lonu = np.asarray(DS['XLONG_U'][0, 0, :].values, dtype=float)
                    DS = DS.assign_coords(west_east_stag=np.mod(lonu + 360., 360.))
                    DS = DS.sortby('west_east_stag')
                    DS['u_stag'] = DS['U'].rename({'west_east_stag': 'lon_u'})
                    DS['u_stag'].attrs['description'] = (
                        'x-wind on the native C-grid (x-staggered) location')
                    DS['u_stag'].attrs['long_name'] = DS['u_stag'].attrs['description']
                    DS = DS.assign_coords(lon_u=DS['west_east_stag'].values)
                    DS['lon_u'].attrs = {'units': 'degrees_E',
                                         'long_name': 'longitude of x-staggered points'}
                    DS = DS.drop_vars('west_east_stag')
                    kept.append('u_stag')
                if 'V' in DS and 'XLAT_V' in DS:
                    latv = np.asarray(DS['XLAT_V'][0, :, 0].values, dtype=float)
                    DS['v_stag'] = DS['V'].rename({'south_north_stag': 'lat_v'})
                    DS['v_stag'].attrs['description'] = (
                        'y-wind on the native C-grid (y-staggered) location')
                    DS['v_stag'].attrs['long_name'] = DS['v_stag'].attrs['description']
                    DS = DS.assign_coords(lat_v=latv)
                    DS['lat_v'].attrs = {'units': 'degrees_N',
                                         'long_name': 'latitude of y-staggered points'}
                    kept.append('v_stag')
                for src, dst, desc, unit in (
                        ('W', 'w_half', 'z-wind on layer interfaces', 'm s-1'),
                        ('PF_PHY', 'phalf3D', 'pressure on layer interfaces', 'Pa'),
                        ('ZF_PHY', 'zhalf', 'height above local surface on layer interfaces', 'm')):
                    if src in DS and 'bottom_top_stag' in DS[src].dims:
                        DS[dst] = DS[src].rename({'bottom_top_stag': model.dim_phalf})
                        DS[dst].attrs = {'description': f'(KEPT POST-PROCESSING) {desc}',
                                         'long_name': f'(KEPT POST-PROCESSING) {desc}',
                                         'units': unit}
                        kept.append(dst)
                if kept:
                    print(f"{Cyan}Keeping native staggered fields: "
                          f"{', '.join(kept)}{Nclr}")

            # Unstagger U, V, W, Zfull onto Regular Grid
            # ------------------------------------------
            # For variables staggered x (lon)
            #   [t,z,y,x'] -> regular mass grid [t,z,y,x]:
            # Step 1: Identify variables with the dimension
            # 'west_east_stag' and _U not in the variable name (these
            # are staggered grid identifiers)
            variables_with_west_east_stag = [var for var in DS.variables if 'west_east_stag' in DS[var].dims and '_U' not in var]

            print(
                f"{Cyan}Interpolating Staggered Variables to Standard Grid{Nclr}"
                )
            # dims_list finds the dims of the variable and replaces
            # west_east_stag with west_east
            print('From west_east_stag to west_east: ' + ', '.join(variables_with_west_east_stag))
            for var_name in variables_with_west_east_stag:
                # Inspiration: pyhton-wrf destag.py
                # https://github.com/NCAR/wrf-python/blob/57116836593b7d7833e11cf11927453c6388487b/src/wrf/destag.py#L9
                var = getattr(DS, var_name)
                dims = var.dims
                dims_list = list(dims)
                for i, dim in enumerate(dims_list):
                    if dim == 'west_east_stag':
                        dims_list[i]='west_east'
                        break # Stop the loop once the replacement is made
                new_dims = tuple(dims_list)
                # Note that XLONG_U is cyclic
                #   LON[x,0] = LON[x,-1] = 0

                transformed_var = (
                    0.5 * (var.isel(west_east_stag=slice(None, -1))
                           + var.isel(west_east_stag=slice(1, None)))
                    )
                DS[var_name] = xr.DataArray(transformed_var,
                                            dims=new_dims,
                                            coords={'XLAT':DS['XLAT']})

                print(f"\n{DS[var_name].attrs['description']}")
                DS[var_name].attrs['description'] = (
                    '(UNSTAGGERED IN POST-PROCESSING) ' + DS[var_name].attrs['description']
                    )
                DS[var_name].attrs['long_name'] = (
                    '(UNSTAGGERED IN POST-PROCESSING) ' + DS[var_name].attrs['description']
                    )
                DS[var_name].attrs['stagger'] = (
                    ('USTAGGERED IN POST-PROCESSING')
                    )

            # For variables staggered y (lat)
            #   [t,z,y',x] -> [t,z,y,x]
            variables_with_south_north_stag = [var for var in DS.variables if 'south_north_stag' in DS[var].dims and '_V' not in var]

            print('From south_north_stag to south_north: '  + ', '.join(variables_with_south_north_stag))

            for var_name in variables_with_south_north_stag:
                var = getattr(DS, var_name)
                dims = var.dims
                dims_list = list(dims)
                for i, dim in enumerate(dims_list):
                    if dim == 'south_north_stag':
                        dims_list[i]='south_north'
                        break # Stop the loop once the replacement is made
                new_dims = tuple(dims_list)

                transformed_var = (
                    0.5 * (var.isel(south_north_stag=slice(None, -1))
                           + var.isel(south_north_stag=slice(1, None)))
                    )
                DS[var_name] = xr.DataArray(transformed_var,
                                            dims=new_dims,
                                            coords={'XLONG':DS['XLONG']})

                DS[var_name].attrs['description'] = (
                    '(UNSTAGGERED IN POST-PROCESSING) ' + DS[var_name].attrs['description']
                    )
                DS[var_name].attrs['long_name'] = (
                    '(UNSTAGGERED IN POST-PROCESSING) ' + DS[var_name].attrs['description']
                    )
                DS[var_name].attrs['stagger'] = (
                    'USTAGGERED IN POST-PROCESSING'
                    )

            # For variables staggered p/z (height)
            #   [t,z',y,x] -> [t,z,y,x]
            variables_with_bottom_top_stag = [var for var in DS.variables if 'bottom_top_stag' in DS[var].dims and 'ZNW' not in var and 'phalf' not in var]

            print('From bottom_top_stag to bottom_top: '  + ', '.join(variables_with_bottom_top_stag))

            for var_name in variables_with_bottom_top_stag:
                var = getattr(DS, var_name)
                dims = var.dims
                dims_list = list(dims)
                for i, dim in enumerate(dims_list):
                    if dim == 'bottom_top_stag':
                        dims_list[i]='bottom_top'
                        break # Stop the loop once the replacement is made
                new_dims = tuple(dims_list)
                transformed_var = (
                    0.5 * (var.sel(bottom_top_stag=slice(None, -1))
                           + var.sel(bottom_top_stag=slice(1, None))))

                DS[var_name] = xr.DataArray(transformed_var, dims=new_dims)
                DS[var_name].attrs['description'] = (
                    '(UNSTAGGERED IN POST-PROCESSING) ' + DS[var_name].attrs['description']
                    )
                DS[var_name].attrs['long_name'] = (
                    '(UNSTAGGERED IN POST-PROCESSING) ' + DS[var_name].attrs['description']
                    )
                DS[var_name].attrs['stagger'] = (
                    'USTAGGERED IN POST-PROCESSING'
                    )

            # Layer-centre height above the local surface [m] (CAP zfull)
            #   1. Z_PHY (layer-centre geopotential height above the
            #      areoid, the model's own value) minus HGT, or minus a
            #      hydrostatic surface estimate from the bottom layer
            #      when HGT is absent (reduced files; ~2 m rms)
            #   2. PH, PHB, HGT: interface geopotential averaged to
            #      layer centres, minus terrain (older files)
            zfull3D = None
            zsrc = None
            if 'Z_PHY' in DS:
                zphy = DS['Z_PHY'].transpose(model.dim_time, model.dim_pfull,
                                             model.dim_lat, model.dim_lon)
                if 'HGT' in DS:
                    zs = DS['HGT'].transpose(model.dim_time, model.dim_lat,
                                             model.dim_lon).values
                    zsrc = 'Z_PHY - HGT'
                elif all(v in DS for v in ('P_PHY', 'T_PHY', 'PSFC')):
                    pp = DS['P_PHY'].transpose(model.dim_time, model.dim_pfull,
                                               model.dim_lat, model.dim_lon).values
                    tt = DS['T_PHY'].transpose(model.dim_time, model.dim_pfull,
                                               model.dim_lat, model.dim_lon).values
                    ps = DS['PSFC'].transpose(model.dim_time, model.dim_lat,
                                              model.dim_lon).values
                    kb = int(np.argmax(np.nanmean(pp, axis=(0, 2, 3))))
                    Rd = float(DS.attrs.get('R_D', 191.8))
                    g = float(DS.attrs.get('G', 3.727))
                    with np.errstate(divide='ignore', invalid='ignore'):
                        zs = (zphy.values[:, kb] - Rd*tt[:, kb]/g
                              * np.log(ps/pp[:, kb]))
                    zsrc = ('Z_PHY - hydrostatic surface estimate from '
                            'the bottom layer (no HGT)')
                    target = getattr(model, 'zsurf')
                    if target not in DS:
                        DS[target] = ((model.dim_time, model.dim_lat,
                                       model.dim_lon), zs)
                        DS[target].attrs['description'] = (
                            '(ADDED POST-PROCESSING) surface height, '
                            'hydrostatic estimate from the bottom layer')
                        DS[target].attrs['long_name'] = DS[target].attrs['description']
                        DS[target].attrs['units'] = 'm'
                if zsrc is not None:
                    zfull3D = zphy.values - zs[:, None, :, :]
            if zfull3D is None and zagl_lvl is not None:
                zagl_lvl = zagl_lvl.transpose(model.dim_time,
                                              'bottom_top_stag', ...)
                zfull3D = 0.5*(zagl_lvl.isel(bottom_top_stag=slice(None, -1)).values
                               + zagl_lvl.isel(bottom_top_stag=slice(1, None)).values)
                zsrc = 'PH, PHB, HGT'
            if zfull3D is not None and 'zhalf' in DS:
                # zhalf kept from ZF_PHY is above the areoid; make it
                # above the local surface like zfull.
                zs_da = xr.DataArray(zs if 'zs' in dir() else
                                     DS['HGT'].transpose(model.dim_time,
                                                         model.dim_lat,
                                                         model.dim_lon).values,
                                     dims=(model.dim_time, model.dim_lat,
                                           model.dim_lon))
                DS['zhalf'] = DS['zhalf'] - zs_da
            if zfull3D is not None:
                DS = DS.assign(zfull=((model.dim_time, model.dim_pfull,
                                       model.dim_lat, model.dim_lon),
                                      zfull3D))
                DS['zfull'].attrs['description'] = (
                    '(ADDED POST-PROCESSING) height above local surface')
                DS['zfull'].attrs['long_name'] = (
                    '(ADDED POST-PROCESSING) height above local surface')
                DS['zfull'].attrs['units'] = 'm'
                print(f"{Cyan}zfull from {zsrc}{Nclr}")
            else:
                print(f"{Yellow}No height information (PH/PHB/HGT or "
                      f"Z_PHY): zfull not written{Nclr}")

            if 'Times' in DS:
                print(f"{Red} Dropping 'Times' variable with non-numerical values")
                DS = DS.drop_vars("Times")

            # Winds: U, V, W sit on the Arakawa-C staggered locations
            # (W on layer interfaces). The physics-grid U_PHY, V_PHY,
            # W_PHY are the model's own values averaged to mass-point
            # layer centres, so they are preferred whenever present;
            # the 2-point destaggering above is the fallback for older
            # files that only carry U, V, W.
            # The target name is whatever the profile lookup settled on
            # ('U' when present in the file, else the CAP name 'ucomp'),
            # so the later renaming step needs no change.
            for cap, phy in (('ucomp', 'U_PHY'), ('vcomp', 'V_PHY'),
                             ('w', 'W_PHY')):
                target = getattr(model, cap)
                if phy in DS:
                    DS[target] = DS[phy]
                    DS[target].attrs['description'] = (
                        f'(FROM {phy} POST-PROCESSING) '
                        + DS[phy].attrs.get('description', ''))
                    DS[target].attrs['long_name'] = DS[target].attrs['description']
                    print(f"{Cyan}{target} taken from {phy}{Nclr}")
                elif target in DS:
                    print(f"{Cyan}{target} from destaggered {target} "
                          f"(no {phy}){Nclr}")

            # Equation of time: recent planetWRF writes LTST, the true
            # local solar time of the model sun. Its departure from the
            # mean local time UT + lon/15 is one scalar per frame
            # (it drifts through the year by tens of minutes on Mars).
            # Save it so MarsFiles -t can shift to true solar time.
            if 'LTST' in DS:
                sols = np.asarray(DS[model.time].values, dtype=float)
                ut_hr = np.mod(sols, 1.)*24.
                lon_hr = np.asarray(DS[model.dim_lon].values, dtype=float)/15.
                ltst = DS['LTST'].transpose(model.dim_time, model.dim_lat,
                                            model.dim_lon).values
                naive = ut_hr[:, None, None] + lon_hr[None, None, :]
                diff = np.mod(ltst - naive + 12., 24.) - 12.
                eot = np.nanmedian(diff.reshape(diff.shape[0], -1), axis=1)
                # Restart frames carry LTST = 0: mark and fill from
                # the nearest valid frame.
                bad = (np.nanmax(np.abs(ltst).reshape(ltst.shape[0], -1),
                                 axis=1) == 0)
                if bad.any() and (~bad).any():
                    idx = np.arange(len(eot))
                    eot[bad] = np.interp(idx[bad], idx[~bad], eot[~bad])
                DS = DS.assign(eot_offset=((model.dim_time,), eot))
                DS['eot_offset'].attrs['description'] = (
                    '(ADDED POST-PROCESSING) equation of time: true local '
                    'solar time minus mean local time (UT + lon/15), from LTST')
                DS['eot_offset'].attrs['long_name'] = DS['eot_offset'].attrs['description']
                DS['eot_offset'].attrs['units'] = 'hr'
                print(f"{Cyan}eot_offset from LTST: {eot.min()*60:+.2f} to "
                      f"{eot.max()*60:+.2f} min{Nclr}")

            # Prune variables CAP cannot carry: anything on a dimension
            # other than time/pfull/phalf/lat/lon (soil layers, dust
            # bins, radiation layers, leftover staggered grids), and
            # vertical-only metadata (ZNU, C1H, ...) that MarsInterp
            # cannot copy into an interpolated file.
            allowed = {model.dim_time, model.dim_pfull, model.dim_phalf,
                       model.dim_lat, model.dim_lon}
            keep_meta = {'ak', 'bk', 'pfull', 'phalf'}
            if getattr(args, 'keep_staggered', False):
                allowed |= {'lon_u', 'lat_v'}
                keep_meta |= {'lon_u', 'lat_v'}
            drop = []
            # Include coordinates: xarray promotes XLAT_U, XLONG_V, ...
            # to coords via the 'coordinates' attribute of U, V.
            for v in list(DS.variables):
                dims = set(DS[v].dims)
                if v in keep_meta or v in DS.dims:
                    continue
                if not dims <= allowed:
                    drop.append(v)
                elif ((model.dim_pfull in dims or model.dim_phalf in dims)
                      and not (dims & {model.dim_lat, model.dim_lon,
                                       'lat_v', 'lon_u'})):
                    # vertical-only metadata (ZNU, C1H, ...)
                    drop.append(v)
            if drop:
                print(f"{Yellow}Dropping {len(drop)} variable(s) on "
                      f"non-CAP dimensions: {', '.join(drop)}{Nclr}")
                DS = DS.drop_vars(drop)

        # --------------------------------------------------------------
        #                    OpenMars Processing
        # --------------------------------------------------------------
        elif model_type == 'openmars':
            # First save all variable FIELDNAM as longname
            var_list = list(DS.data_vars) + list(DS.coords)
            for var_name in var_list:
                var = DS[var_name]
                if 'FIELDNAM' in var.attrs:
                    var.attrs['long_name'] = var.attrs['FIELDNAM']
                if 'UNITS' in var.attrs:
                    var.attrs['units'] = var.attrs['UNITS']

            # Define Coordinates for New DataFrame
            time = DS[model.dim_time] # min since simulation start [m]
            lat = DS[model.dim_lat] # Replace DS.lat
            lon = DS[model.dim_lon]

            DS = DS.assign(pfull = DS[model.dim_pfull]*ref_press)

            DS['pfull'].attrs['long_name'] = (
                '(ADDED POST-PROCESSING) reference pressure'
                )
            DS['pfull'].attrs['units'] = ('Pa')

            # add ak,bk as variables
            # add p_half dimensions as vertical grid coordinate

            # Compute sigma values. Swap the sigma array upside down
            # twice with [::-1] because layers_mid_point_to_boundary()
            # needs (sigma[0] = 0, sigma[-1] = 1).
            # Then reorganize in the original openMars format with
            # (sigma[0] = 1, sigma[-1] = 0)
            bk = layers_mid_point_to_boundary(DS[model.dim_pfull][::-1], 1.)[::-1]
            ak = np.zeros(len(DS[model.dim_pfull]) + 1)

            DS[model.phalf] = (ak + ref_press*bk)
            DS.phalf.attrs['long_name'] = (
                '(ADDED POST-PROCESSING) pressure at layer interfaces'
                )
            DS.phalf.attrs['description'] = (
                '(ADDED POST-PROCESSING) pressure at layer interfaces'
                )
            DS.phalf.attrs['units'] = ('Pa')

            DS = DS.assign(bk=(model.dim_phalf,
                               np.array(bk)))
            DS = DS.assign(ak=(model.dim_phalf,
                               np.zeros(len(DS[model.dim_pfull]) + 1)))

            # Update Variable Description & Longname
            DS['ak'].attrs['long_name'] = (
                '(ADDED POST-PROCESSING) pressure part of the hybrid coordinate'
                )
            DS['ak'].attrs['units'] = ('Pa')
            DS['bk'].attrs['long_name'] = (
                '(ADDED POST-PROCESSING) vertical coordinate sigma value'
                )
            DS['bk'].attrs['units'] = ('None')

        # --------------------------------------------------------------
        #                     Emars Processing
        # --------------------------------------------------------------
        elif model_type == 'emars':
            # Interpolate U, V, onto Regular Mass Grid (from staggered)
            print(
                f"{Cyan}Interpolating Staggered Variables to Standard Grid"
                )

            variables_with_latu = [var for var in DS.variables if 'latu' in DS[var].dims]
            variables_with_lonv = [var for var in DS.variables if 'lonv' in DS[var].dims]

            print(f"{Cyan}Changing time units from hours to sols]")
            DS[model.time] = DS[model.time].values/24.
            DS[model.time].attrs['long_name'] = 'time'
            DS[model.time].attrs['units'] = 'days since 0000-00-00 00:00:00'

            print(f"{Cyan}Converting reference pressure to [Pa]")
            # DS[model.pfull] = DS[model.pfull].values*100
            new_pfull_vals = DS[model.pfull].values * 100
            DS = DS.assign_coords({model.pfull: new_pfull_vals})
            DS[model.pfull].attrs['units'] = 'Pa'

            # dims_list process finds dims of the variable and replaces
            # west_east_stag with west_east
            print('From latu to lat: ' + ', '.join(variables_with_latu ))
            for var_name in variables_with_latu:
                var = getattr(DS, var_name)
                dims = var.dims
                longname_txt = var.long_name
                units_txt = var.units

                # Replace latu dims with lat
                dims_list = list(dims)
                for i, dim in enumerate(dims_list):
                    if dim == 'latu':
                        dims_list[i] = 'lat'
                        break # Stop the loop when replacement made
                new_dims = tuple(dims_list)

                #TODO Using 'values' as var.isel(latu=slice(None, -1)).values and var.isel(latu=slice(None, -1)) returns array of different size
                #TODO the following reproduce the 'lat' array, but not if the keyword values is ommited
                #latu1_val=ds.latu.isel(latu=slice(None, -1)).values
                #latu2_val=ds.latu.isel(latu=slice(1,  None)).values

                #newlat_val= 0.5 * (latu1 + latu2)
                #newlat_val=np.append(newlat_val,0)
                #newlat_val ==lat
                transformed_var = (
                    0.5 * (var.isel(latu=slice(None, -1)).values
                           + var.isel(latu=slice(1, None)).values)
                    )

                # This is equal to lat[0:-1]
                # Add padding at the pole to conserve the same dimension as lat
                list_pads = []
                for i, dim in enumerate(new_dims):
                    if dim == 'lat':
                        list_pads.append((0,1))
                        # (begining, end)=(0, 1), meaning 1 padding at
                        # the end of the array
                    else:
                        list_pads.append((0,0))
                        # (begining, end)=(0, 0), no pad on that axis

                transformed_var = np.pad(transformed_var,
                                         list_pads,
                                         mode='constant',
                                         constant_values=0)

                DS[var_name] = xr.DataArray(transformed_var,
                                            dims=new_dims,
                                            coords={'lat':DS['lat']})
                DS[var_name].attrs['long_name'] = (
                    f'(UNSTAGGERED IN POST-PROCESSING) {longname_txt}'
                    )
                DS[var_name].attrs['units'] = units_txt

            for var_name in variables_with_lonv:
                var = getattr(DS, var_name)
                dims = var.dims
                longname_txt = var.long_name
                units_txt = var.units

                # Replace lonv in dimensions with lon
                dims_list = list(dims)
                for i, dim in enumerate(dims_list):
                    if dim == 'lonv':
                        dims_list[i] = 'lon'
                        break # Stop loop once the replacement is made
                new_dims = tuple(dims_list)

                transformed_var = (
                    0.5 * (var.isel(lonv=slice(None, -1)).values
                           + var.isel(lonv=slice(1, None)).values)
                    )
                # This is equal to lon[0:-1]
                # Add padding
                list_pads = []
                for i, dim in enumerate(new_dims):
                    if dim == 'lon':
                        list_pads.append((0, 1))
                        # (begining, end)=(0, 1), meaning 1 padding at
                        # the end of the array
                    else:
                        list_pads.append((0, 0))
                        # (begining, end)=(0, 0), no pad on that axis
                transformed_var = np.pad(transformed_var, list_pads,
                                         mode='wrap')
                #TODO with this method V[0] =V[-1]: can we add cyclic point before destaggering?

                DS[var_name] = xr.DataArray(transformed_var,
                                            dims=new_dims,
                                            coords={'lon':DS['lon']})

                DS[var_name].attrs['long_name'] = (
                    f'(UNSTAGGERED IN POST-PROCESSING) {longname_txt}'
                    )
                DS[var_name].attrs['units'] = units_txt
            DS.drop_vars(['latu','lonv'])

        # --------------------------------------------------------------
        #                      PCM Processing
        # --------------------------------------------------------------
        elif model_type == 'pcm':
            """
            Process PCM model output:
            1. Create pfull and phalf pressure coordinates
            2. Ensure correct vertical ordering (lowest pressure at top)
            3. Set attributes to prevent double-flipping
            """
            print(f"{Cyan}Processing pcm file")
            # Adding long_name attibutes
            for var_name in DS.data_vars:
                var = DS[var_name]
                if 'title' in var.attrs:
                    var.attrs['long_name'] = var.attrs['title']

            # Print the values for debugging
            if debug:
                print(f"ap = {DS.ap.values}")
                print(f"bp = {DS.bp.values}")

            # Create pfull variable
            pfull = (DS.aps.values + DS.bps.values*ref_press)
            DS['pfull'] = (['altitude'],  pfull)
            DS['pfull'].attrs['long_name'] = (
                '(ADDED POST-PROCESSING), reference pressure'
                )
            DS['pfull'].attrs['units'] = 'Pa'

            # Replace the PCM phalf creation section with:
            if 'ap' in DS and 'bp' in DS:
                # Calculate phalf values from ap and bp
                phalf_values = (DS.ap.values + DS.bp.values*ref_press)

                # Check the order - for vertical pressure coordinates, we want:
                # - Lowest pressure (top of atmosphere) at index 0
                # - Highest pressure (surface) at index -1
                if phalf_values[0] > phalf_values[-1]:
                    # Currently highest pressure is at index 0, so we need to flip
                    print(f"{Yellow}PCM phalf has highest pressure at index 0, flipping to standard orientation")
                    phalf_values = phalf_values[::-1]
                    DS.attrs['vertical_dimension_flipped'] = True
                else:
                    # Already in the correct orientation
                    print(f"{Green}PCM phalf values already in correct orientation (lowest at index 0)")
                    DS.attrs['vertical_dimension_flipped'] = False

                # Store phalf values in the dataset
                DS['phalf'] = (['interlayer'], phalf_values)
                DS['phalf'].attrs['long_name'] = '(ADDED POST-PROCESSING) pressure at layer interfaces'
                DS['phalf'].attrs['units'] = 'Pa'

                # Also need to fix pfull to match phalf orientation
                pfull = DS['pfull'].values
                if DS.attrs['vertical_dimension_flipped'] and pfull[0] > pfull[-1]:
                    # If we flipped phalf, also ensure pfull has lowest pressure at index 0
                    if DS['pfull'].values[0] > DS['pfull'].values[-1]:
                        DS['pfull'] = (['altitude'], DS['pfull'].values[::-1])
                        print(f"{Yellow}Also flipped pfull values to match phalf orientation")

        # --------------------------------------------------------------
        #                START PROCESSING FOR ALL MODELS
        # --------------------------------------------------------------
        if model_type == 'pcm' and 'vertical_dimension_flipped' in DS.attrs:
            print(f"{Cyan}Using PCM-specific vertical orientation handling")
            # Skip automatic flipping - we've already handled it in PCM processing
        else:
            # Standard vertical processing for other models
            if DS[model.pfull].values[0] != DS[model.pfull].values.min():
                # Flip vertical dimensions using explicit index arrays
                # This approach is more robust across xarray versions 
                # than slice(None, None, -1) which can cause dimension 
                # tracking issues in xarray >= 2025.12.0
                n_pfull = DS.sizes[model.dim_pfull]
                DS = DS.isel(**{model.dim_pfull: list(range(n_pfull - 1, -1, -1))})
                # Flip phalf:
                n_phalf = DS.sizes[model.dim_phalf]
                DS = DS.isel(**{model.dim_phalf: list(range(n_phalf - 1, -1, -1))})
                print(f"{Red}NOTE: all variables flipped along vertical dimension. "
                      f"Top of the atmosphere is now index = 0")
                 
        # Reorder dimensions
        print(f"{Cyan} Transposing variable dimensions to match order "
              f"expected in CAP")
        DS = DS.transpose(model.dim_time, model.dim_pfull, model.dim_lat,
                          model.dim_lon, ...)
        # Native staggered fields kept by MarsFormat -stag (marswrf) use
        # their own horizontal/vertical axes; order them explicitly
        for v, order in (
                ('v_stag', (model.dim_time, model.dim_pfull, 'lat_v', model.dim_lon)),
                ('w_half', (model.dim_time, model.dim_phalf, model.dim_lat, model.dim_lon)),
                ('phalf3D', (model.dim_time, model.dim_phalf, model.dim_lat, model.dim_lon)),
                ('zhalf', (model.dim_time, model.dim_phalf, model.dim_lat, model.dim_lon))):
            if v in DS:
                DS[v] = DS[v].transpose(*order)

        # Change longitude from -180-179 to 0-360
        if min(DS[model.dim_lon]) < 0:
            tmp = np.array(DS[model.dim_lon])
            tmp = np.where(tmp<0, tmp + 360, tmp)
            # Use assign_coords for robust coordinate update across 
            # xarray versions. Direct assignment 
            # (DS[model.dim_lon] = tmp) can fail in xarray >= 2025.12.0
            DS = DS.assign_coords({model.dim_lon: tmp})
            DS = DS.sortby(model.dim_lon)
            DS[model.lon].attrs['long_name'] = (
                '(MODIFIED POST-PROCESSING) longitude'
                )
            DS[model.lon].attrs['units'] = ('degrees_E')
            print(f"{Red} NOTE: Longitude changed to 0-360E and all variables "
                  f"appropriately reindexed")

        # Add scalar axis to areo [time, scalar_axis])
        inpt_dimlist = DS.dims
        # First check if dims are correct - don't need to be modified
        if 'scalar_axis' not in inpt_dimlist:
            # If scalar axis is a dimension
            scalar_axis = DS.assign_coords(scalar_axis=1)
        if DS[model.areo].dims != (model.time,scalar_axis):
            DS[model.areo] = DS[model.areo].expand_dims('scalar_axis', axis=1)
            DS[model.areo].attrs['long_name'] = (
                '(SCALAR AXIS ADDED POST-PROCESSING)'
                )

            print(f"{Red}NOTE: scalar axis added to aerocentric longitude")

        # Standardize variables names if requested
        if args.retain_names:
            print(f"{Purple}Preserving original names for variable and "
                  f"dimensions")
            ext = f'{ext}_nat'
        else:
            print(f"{Purple}Using standard FV3 names for variables and "
                  f"dimensions")

            # Model has {'ucomp':'U','temp':'T', 'dim_lon'='XLON'...}
            # Create a reversed dictionary, e.g. {'U':'ucomp','T':'temp'...}
            # to revert to the original variable names before archiving
            # Note that model.() is constructed only from
            # .amescap_profile and may include names not present in file
            model_dims = dict()
            model_vars = dict()

            # Loop over optential dims and vars in model
            for key_i in model.__dict__.keys():
                val = getattr(model,key_i)
                if key_i[0:4] != 'dim_':
                    # Potential variables
                    if val in (list(DS.keys()) + list(DS.coords)):
                        # Check if the key is a variable (e.g. temp)
                        # or coordinate (e.g. lat)
                        model_vars[val] = key_i
                else:
                    # Potential dimensions
                    if val in list(DS.dims):
                        model_dims[val] = key_i[4:]

            # Sort the dictionaries to remove duplicates
            # Remove key/val duplicates: e.g if {'lat':'lat'}, there is
            # no need to rename that variable
            model_dims = {key: val for key, val in model_dims.items() if key != val}
            model_vars = {key: val for key, val in model_vars.items() if key != val}

            # Avoiding conflict with derived 'temp' variable in MarsWRF
            # T is perturb T in MarsWRF, and temp was derived earlier
            if (model_type == 'marswrf' and
                'T' in model_vars and
                model_vars['T'] == 'temp'):
                print(f"{Yellow}Note: Removing 'T' from variable mapping for "
                      f"MarSWRF to avoid conflict with derived 'temp'")
                del model_vars['T']  # Remove the T -> temp mapping

            print(f"DEBUG: Model dimensions: {model_dims}")
            print(f"DEBUG: Model variables: {model_vars}")

            # Special handling for PCM to avoid dimension swap errors
            dimension_swap_failed = False
            if model_type == 'pcm':
                try:
                    DS = DS.swap_dims(dims_dict = model_dims)
                except ValueError as e:
                    if "replacement dimension" in str(e):
                        print(f"{Yellow}Warning: PCM dimension swap failed. "
                            f"Automatically using retain_names approach for dimensions.{Nclr}")
                        # Skip the dimension swap but continue with variable renaming
                        dimension_swap_failed = True
                    else:
                        # Re-raise other errors
                        raise
            else:
                # Normal processing for other GCM types
                DS = DS.swap_dims(dims_dict = model_dims)

            # Continue with variable renaming regardless of dimension swap status
            if not dimension_swap_failed:
                DS = DS.rename_vars(name_dict = model_vars)
                # print(f"{Cyan}Renamed variables:\n{list(DS.variables)}{Nclr}\n")
                # Update CAP's internal variables dictionary
                model = reset_FV3_names(model)
            else:
                # If dimension swap failed, still rename variables but handle as if using retain_names
                DS = DS.rename_vars(name_dict = model_vars)
                # print(f"{Cyan}Renamed variables (with original dimensions):\n{list(DS.variables)}{Nclr}\n")
                # Add the _nat suffix as if -rn was used, but we still renamed variables
                if '_nat' not in ext:
                    ext = f'{ext}_nat'

        # --------------------------------------------------------------
        # CREATE ATMOS_DAILY, ATMOS_AVERAGE, & ATMOS_DIURN FILES
        # --------------------------------------------------------------
        if args.bin_average and not args.bin_diurn:
            ext = f'{ext}_average'
            nday = args.bin_average

            # Calculate time step from original unmodified values
            dt_in = float(original_time_vals[1] - original_time_vals[0])
            print(f"DEBUG: Using original time values with dt_in = {dt_in}")

            # Convert time step to days based on original units
            dt_days = dt_in
            if 'minute' in original_time_units.lower() or 'minute' in original_time_desc.lower():
                dt_days = dt_in / 1440.0  # Convert minutes to days
                print(f"DEBUG: Converting {dt_in} minutes to {dt_days} days")
            elif 'hour' in original_time_units.lower() or 'hour' in original_time_desc.lower():
                dt_days = dt_in / 24.0  # Convert hours to days
                print(f"DEBUG: Converting {dt_in} hours to {dt_days} days")
            else:
                print(f"DEBUG: No time unit found in original attributes, assuming 'days'")

            # Check if bin size is appropriate
            if dt_days >= nday:
                print(f"{Red}***Error***: Requested bin size ({nday} days) is smaller than or equal to "
                    f"the time step in the data ({dt_days:.2f} days)")
                continue  # Skip to next file

            # Calculate samples per day and samples per bin
            samples_per_day = 1.0 / dt_days
            samples_per_bin = nday * samples_per_day

            # Need at least one sample per bin
            if samples_per_bin < 1:
                print(f"{Red}***Error***: Time sampling in file ({1.0/samples_per_day:.2f} days "
                    f"between samples) is too coarse for {nday}-day bins")
                continue  # Skip to next file

            # Round to nearest integer for coarsen function
            combinedN = max(1, int(round(samples_per_bin)))
            print(f"DEBUG: Using {combinedN} time steps per {nday}-day bin")

            # Coarsen and average
            DS_average = DS.coarsen(**{model.dim_time:combinedN}, boundary='trim').mean()

            # Update the time coordinate attribute
            DS_average[model.dim_time].attrs['long_name'] = (
                f'time averaged over {nday} sols'
                )

            # For PCM files, ensure vertical orientation is preserved after averaging
            if model_type == 'pcm':
                # Check phalf values and orientation
                phalf_vals = DS_average['phalf'].values
                if len(phalf_vals) > 1:  # Only check if we have more than one value
                    # Correct orientation: lowest pressure at index 0, highest at index -1
                    if phalf_vals[0] > phalf_vals[-1]:
                        print(f"{Yellow}Warning: phalf orientation incorrect after binning, fixing...")
                        DS_average['phalf'] = (['interlayer'], phalf_vals[::-1])

            # Create New File, set time dimension as unlimitted
            base_name = os.path.splitext(fullnameIN)[0]
            fullnameOUT = f"{base_name}{ext}.nc"
            DS_average.to_netcdf(fullnameOUT, unlimited_dims=model.dim_time,
                                 format='NETCDF4_CLASSIC')

        elif args.bin_diurn:
            ext = f'{ext}_diurn'
            print(f"Doing diurnal binning")

            # Custom number of sols
            if args.bin_average:
                nday = args.bin_average
            else:
                nday = 5

            print(f"Using {nday}-day bins then diurnal binning")

            # Calculate time step from original unmodified values
            dt_in = float(original_time_vals[1] - original_time_vals[0])
            print(f"DEBUG: Using original time values with dt_in = {dt_in}")

            # Convert time step to days based on original units
            dt_days = dt_in
            if 'minute' in original_time_units.lower() or 'minute' in original_time_desc.lower():
                dt_days = dt_in / 1440.0  # Convert minutes to days
                print(f"DEBUG: Converting {dt_in} minutes to {dt_days} days")
            elif 'hour' in original_time_units.lower() or 'hour' in original_time_desc.lower():
                dt_days = dt_in / 24.0  # Convert hours to days
                print(f"DEBUG: Converting {dt_in} hours to {dt_days} days")
            else:
                print(f"DEBUG: No time unit found in original attributes, assuming 'days'")

            # Calculate samples per day and check if valid
            samples_per_day = 1.0 / dt_days
            if samples_per_day < 1:
                print(f"{Red}***Error***: Operation not permitted because "
                    f"time sampling in file < one time step per day")
                continue  # Skip to next file

            # Calculate number of steps to bin
            iperday = int(round(samples_per_day))
            combinedN = int(iperday * nday)
            print(f"DEBUG: Using {combinedN} time steps per {nday}-day bin with {iperday} samples per day")

            # Output Binned Data to New  **atmos_diurn.nc file
            # create a new time of day dimension
            tod_name = f"time_of_day_{iperday:02d}"
            days = len(DS[model.dim_time]) / iperday

            # Initialize the new dataset
            DS_diurn = None

            for i in range(0, int(days/nday)):
                # Slice original dataset in 5 sol increments
                downselect = (
                    DS.isel(**{model.dim_time:slice(i*combinedN,
                                                    i*combinedN + combinedN)})
                    )

                # Rename time dimension to time of day and find the
                # local time equivalent
                downselect = downselect.rename({model.dim_time: tod_name})
                downselect[tod_name] = (
                    np.mod(downselect[tod_name]*24, 24).values
                    )

                # Average all instances of the same local time
                idx = downselect.groupby(tod_name).mean()

                # Add back in the time dimensionn
                idx = idx.expand_dims({model.dim_time: [i]})

                # Concatenate into new diurn array with a local time
                # and time dimension (stack along time)
                if DS_diurn is None:
                    DS_diurn = idx
                else:
                    DS_diurn = xr.concat([DS_diurn, idx], dim=model.dim_time)
            #TODO
            # ==== Overwrite the ak, bk arrays=== [AK]
            #For some reason I can't track down, the ak('phalf') and bk('phalf')
            # turn into ak ('time', 'time_of_day_12','phalf'), in PCM which messes
            # the pressure interpolation.
            # Safe approach to fix the dimensions for ak/bk arrays
            # First check if these variables exist in the diurn dataset
            ak_dims = None
            bk_dims = None

            if model.ak in DS_diurn and model.bk in DS_diurn:
                # Get the dimensions and print for debugging
                ak_dims = DS_diurn[model.ak].dims
                bk_dims = DS_diurn[model.bk].dims
                print(f"DEBUG: {ak_dims} and {bk_dims}")

                # Use explicit dimension checks (safer than len(dims) > 1)
                if any(dim != model.dim_phalf for dim in ak_dims):
                    print(f"DEBUG: Fixing dimensions for {model.ak} and {model.bk} in diurn file")
                    # Ensure we're assigning the correct structure
                    DS_diurn[model.ak] = DS[model.ak].copy()
                    DS_diurn[model.bk] = DS[model.bk].copy()

            # replace the time dimension with the time dimension from DS_average
            time_DS = DS[model.dim_time]
            time_avg_DS = time_DS.coarsen(**{model.dim_time:combinedN},boundary='trim').mean()
            DS_diurn[model.dim_time] = time_avg_DS[model.dim_time]

            # Update the time coordinate attribute
            DS_diurn[model.dim_time].attrs['long_name'] = (
                f'time averaged over {nday} sols'
                )
            
            # Safe phalf check for PCM files
            # During diurn processing, phalf can gain extra dimensions from 
            # groupby().mean() operations (similar issue to ak/bk above).
            # We need to restore phalf from the original DS if it has wrong dimensions.
            if model_type == 'pcm' and DS_diurn is not None and 'phalf' in DS_diurn:
                try:
                    phalf_dims = DS_diurn['phalf'].dims
                    # Check if phalf has acquired extra dimensions (should only have 1)
                    if len(phalf_dims) > 1:
                        print(f"DEBUG: phalf has extra dimensions {phalf_dims}, restoring from original DS")
                        # Restore phalf from original dataset
                        if 'phalf' in DS:
                            DS_diurn['phalf'] = DS['phalf'].copy()
                    
                    # Now check orientation using the corrected phalf
                    phalf_vals = DS_diurn['phalf'].values
                    # Ensure we have a 1D array before checking orientation
                    if phalf_vals.ndim == 1 and len(phalf_vals) > 1:
                        # Extract actual values and convert to regular Python floats
                        first_val = float(phalf_vals[0])
                        last_val = float(phalf_vals[-1])
                        if first_val > last_val:
                            print(f"{Yellow}Warning: phalf orientation incorrect in diurn file, fixing...")
                            DS_diurn['phalf'] = (DS_diurn['phalf'].dims, phalf_vals[::-1])
                except Exception as e:
                    print(f"{Yellow}Note: Could not check phalf orientation: {str(e)}")

            # Create New File, set time dimension as unlimitted
            fullnameOUT = f'{fullnameIN[:-3]}{ext}.nc'
            DS_diurn.to_netcdf(fullnameOUT, unlimited_dims=model.dim_time,
                               format='NETCDF4_CLASSIC')

        else:
            ext = f'{ext}_daily'
            fullnameOUT = f'{fullnameIN[:-3]}{ext}.nc'
            DS.to_netcdf(fullnameOUT, unlimited_dims=model.dim_time,
                         format='NETCDF4_CLASSIC')
        print(f"{Cyan}{fullnameOUT} was created")


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
