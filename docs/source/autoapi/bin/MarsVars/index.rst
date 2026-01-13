:py:mod:`bin.MarsVars`
======================

.. py:module:: bin.MarsVars

.. autoapi-nested-parse::

   The MarsVars executable is for performing variable manipulations in
   existing files. Most often, it is used to derive and add variables to
   existing files, but it also differentiates variables with respect to
   (w.r.t) the Z axis, column-integrates variables, converts aerosol
   opacities from opacity per Pascal to opacity per meter, removes and
   extracts variables from files, and enables scaling variables or editing
   variable names, units, etc.

   The executable requires:

       * ``[input_file]``           The file to be transformed

   and optionally accepts:

       * ``[-add --add_variable]``          Derive and add variable to file
       * ``[-zdiff --differentiate_wrt_z]`` Differentiate variable w.r.t. Z axis
       * ``[-col --column_integrate]``      Column-integrate variable
       * ``[-zd --zonal_detrend]``          Subtract zonal mean from variable
       * ``[-to_dz --dp_to_dz]``            Convert aerosol opacity op/Pa -> op/m
       * ``[-to_dp --dz_to_dp]``            Convert aerosol opacity op/m -> op/Pa
       * ``[-rm --remove_variable]``        Remove variable from file
       * ``[-extract --extract_copy]``      Copy variable to new file
       * ``[-edit --edit_variable]``        Edit variable attributes or scale it

   Third-party Requirements:

       * ``sys``
       * ``argparse``
       * ``os``
       * ``warnings``
       * ``re``
       * ``numpy``
       * ``netCDF4``
       * ``shutil``
       * ``functools``
       * ``traceback``
       * ``matplotlib``
       * ``time``
       * ``io``
       * ``locale``
       * ``amescap``



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   bin.MarsVars.add_help
   bin.MarsVars.check_dependencies
   bin.MarsVars.check_variable_exists
   bin.MarsVars.compute_DP_3D
   bin.MarsVars.compute_DZ_3D
   bin.MarsVars.compute_DZ_full_pstd
   bin.MarsVars.compute_Ek
   bin.MarsVars.compute_Ep
   bin.MarsVars.compute_MF
   bin.MarsVars.compute_N
   bin.MarsVars.compute_Tco2
   bin.MarsVars.compute_Vg_sed
   bin.MarsVars.compute_WMFF
   bin.MarsVars.compute_dustref_per_pa
   bin.MarsVars.compute_dustref_per_z
   bin.MarsVars.compute_mmr
   bin.MarsVars.compute_p_3D
   bin.MarsVars.compute_rho
   bin.MarsVars.compute_scorer
   bin.MarsVars.compute_theta
   bin.MarsVars.compute_w
   bin.MarsVars.compute_w_net
   bin.MarsVars.compute_xzTau
   bin.MarsVars.compute_zfull
   bin.MarsVars.compute_zhalf
   bin.MarsVars.debug_wrapper
   bin.MarsVars.ensure_file_closed
   bin.MarsVars.force_close_netcdf_files
   bin.MarsVars.get_existing_var_name
   bin.MarsVars.main
   bin.MarsVars.patched_print_message
   bin.MarsVars.process_add_variables
   bin.MarsVars.safe_copy_replace
   bin.MarsVars.safe_move_file
   bin.MarsVars.safe_print
   bin.MarsVars.safe_remove_file



Attributes
~~~~~~~~~~

.. autoapisummary::

   bin.MarsVars.C_dst
   bin.MarsVars.C_ice
   bin.MarsVars.Cp
   bin.MarsVars.Kb
   bin.MarsVars.M_co2
   bin.MarsVars.N
   bin.MarsVars.Na
   bin.MarsVars.Qext_dst
   bin.MarsVars.Qext_ice
   bin.MarsVars.R
   bin.MarsVars.Rd
   bin.MarsVars.Reff_dst
   bin.MarsVars.Reff_ice
   bin.MarsVars.S0
   bin.MarsVars.T0
   bin.MarsVars.Tpole
   bin.MarsVars.amu
   bin.MarsVars.amu_co2
   bin.MarsVars.args
   bin.MarsVars.cap_str
   bin.MarsVars.debug
   bin.MarsVars.exit_code
   bin.MarsVars.filepath
   bin.MarsVars.fill_value
   bin.MarsVars.g
   bin.MarsVars.mass_co2
   bin.MarsVars.master_list
   bin.MarsVars.n0
   bin.MarsVars.original_print_message
   bin.MarsVars.parser
   bin.MarsVars.psrf
   bin.MarsVars.rgas
   bin.MarsVars.rho_air
   bin.MarsVars.rho_dst
   bin.MarsVars.rho_ice
   bin.MarsVars.sigma


.. py:function:: add_help(var_list)

   Create a help string for the add_variable argument.

   :param var_list: Dictionary of variables and their attributes
   :type  var_list: dict
   :return: Formatted help string
   :rtype:  str


.. py:function:: check_dependencies(f, var, master_list, add_missing=True, dependency_chain=None)

   Check for variable dependencies in a file, add missing dependencies.

   :param f: NetCDF file object
   :param var: Variable to check deps. for
   :param master_list: Dict of supported vars and their deps.
   :param add_missing: Whether to try adding missing deps. (default: True)
   :param dependency_chain: List of vars in the current dep. chain (for detecting cycles)
   :return: True if all deps. are present or successfully added, False otherwise
   :raises RuntimeError: If the variable is not in the master list
   :raises Exception: If any other error occurs


.. py:function:: check_variable_exists(var_name, file_vars)

   Check if a variable exists in a file.

   Considers alternative naming conventions.

   :param var_name: Variable name to check
   :type  var_name: str
   :param file_vars: Set of variable names in the file
   :type  file_vars: set
   :return: True if the variable exists, False otherwise
   :rtype:  bool
   :raises ValueError: If the input dimensions are not compatible
   :raises TypeError: If the input types are not compatible
   :raises Exception: If any other error occurs


.. py:function:: compute_DP_3D(ps, ak, bk, shape_out)

   Calculate the thickness of a layer in pressure units.

   :param ps: Surface pressure (Pa)
   :type  ps: array [time, lat, lon]
   :param ak: Vertical coordinate pressure value (Pa)
   :type  ak: array [phalf]
   :param bk: Vertical coordinate sigma value (None)
   :type  bk: array [phalf]
   :param shape_out: Determines how to handle the dimensions of DP_3D.
       If len(time) = 1 (one timestep), DP_3D is returned as
       [1, lev, lat, lon] as opposed to [lev, lat, lon]
   :type  shape_out: float
   :return: ``DP`` Layer thickness in pressure units (Pa)
   :rtype:  array [time, lev, lat, lon]
   :raises ValueError: If the input dimensions are not compatible
   :raises TypeError: If the input types are not compatible
   :raises Exception: If any other error occurs


.. py:function:: compute_DZ_3D(ps, ak, bk, temp, shape_out)

   Calculate the thickness of a layer in altitude units.

   :param ps: Surface pressure (Pa)
   :type  ps: array [time, lat, lon]
   :param ak: Vertical coordinate pressure value (Pa)
   :type  ak: array [phalf]
   :param bk: Vertical coordinate sigma value (None)
   :type  bk: array [phalf]
   :param shape_out: Determines how to handle the dimensions of DZ_3D.
       If len(time) = 1 (one timestep), DZ_3D is returned as
       [1, lev, lat, lon] as opposed to [lev, lat, lon]
   :type  shape_out: float
   :return: ``DZ`` Layer thickness in altitude units (m)
   :rtype:  array [time, lev, lat, lon]
   :raises ValueError: If the input dimensions are not compatible
   :raises TypeError: If the input types are not compatible
   :raises Exception: If any other error occurs


.. py:function:: compute_DZ_full_pstd(pstd, T, ftype='average')

   Calculate layer thickness.

   Computes from the midpoint of the standard pressure levels (``pstd``).

   In this context, ``pfull=pstd`` with the layer interfaces
   defined somewhere in between successive layers::

       --- Nk --- TOP       ========  phalf
       --- Nk-1 ---
                            --------  pfull = pstd    ^
                                                      | DZ_full_pstd
                            ========  phalf           |
       --- 1 ---            --------  pfull = pstd    v
       --- 0 --- SFC        ========  phalf
                             / / / /

   :param pstd: Vertical coordinate (pstd; Pa)
   :type  pstd: array [lev]
   :param T: Temperature (K)
   :type  T: array [time, lev, lat, lon]
   :param f_type: The FV3 file type: diurn, daily, or average
   :type  f_stype: str
   :return: DZ_full_pstd, Layer thicknesses (Pa)
   :rtype:  array [time, lev, lat, lon]
   :raises ValueError: If the input dimensions are not compatible
   :raises TypeError: If the input types are not compatible
   :raises Exception: If any other error occurs
   :raises RuntimeError: If the layer thickness calculation fails
   :raises ZeroDivisionError: If the temperature is zero


.. py:function:: compute_Ek(ucomp, vcomp)

   Calculate wave kinetic energy

   Calculation::

       Ek = 1/2 (u'**2+v'**2)

   :param ucomp: Zonal wind (m/s)
   :type  ucomp: array [time, lev, lat, lon]
   :param vcomp: Meridional wind (m/s)
   :type  vcomp: array [time, lev, lat, lon]
   :return: ``Ek`` Wave kinetic energy (J/kg)
   :rtype:  array [time, lev, lat, lon]
   :raises ValueError: If the input dimensions are not compatible
   :raises TypeError: If the input types are not compatible
   :raises Exception: If any other error occurs


.. py:function:: compute_Ep(temp)

   Calculate wave potential energy.

   Calculation::

       Ep = 1/2 (g/N)^2 (temp'/temp)^2

   :param temp: Temperature (K)
   :type  temp: array [time, lev, lat, lon]
   :return: ``Ep`` Wave potential energy (J/kg)
   :rtype:  array [time, lev, lat, lon]
   :raises ValueError: If the input dimensions are not compatible
   :raises TypeError: If the input types are not compatible
   :raises Exception: If any other error occurs


.. py:function:: compute_MF(UVcomp, w)

   Calculate zonal or meridional momentum fluxes.

   :param UVcomp: Zonal or meridional wind (ucomp or vcomp)(m/s)
   :type  UVcomp: array
   :param w: Vertical wind (m/s)
   :type  w: array [time, lev, lat, lon]
   :return: ``u'w'`` or ``v'w'``, Zonal/meridional momentum flux (J/kg)
   :rtype:  array [time, lev, lat, lon]
   :raises ValueError: If the input dimensions are not compatible
   :raises TypeError: If the input types are not compatible
   :raises Exception: If any other error occurs


.. py:function:: compute_N(theta, zfull)

   Calculate the Brunt Vaisala freqency.

   :param theta: Potential temperature (K)
   :type  theta: array [time, lev, lat, lon]
   :param zfull: Altitude above ground level at the layer midpoint (m)
   :type  zfull: array [time, lev, lat, lon]
   :return: ``N``, Brunt Vaisala freqency [rad/s]
   :rtype:  array [time, lev, lat, lon]
   :raises ValueError: If the input dimensions are not compatible
   :raises TypeError: If the input types are not compatible
   :raises Exception: If any other error occurs


.. py:function:: compute_Tco2(P_3D)

   Calculate the frost point of CO2.

   Adapted from Fannale (1982) - Mars: The regolith-atmosphere cap
   system and climate change. Icarus.

   :param P_3D: The full 3D pressure array (Pa)
   :type  p_3D: array [time, lev, lat, lon]
   :return: CO2 frost point [K]
   :rtype:  array [time, lev, lat, lon]
   :raises ValueError: If the input dimensions are not compatible
   :raises TypeError: If the input types are not compatible
   :raises Exception: If any other error occurs


.. py:function:: compute_Vg_sed(xTau, nTau, T)

   Calculate the sedimentation rate of the dust.
   [Courtney Batterson, 2023]

   :param xTau: Dust or ice MASS mixing ratio (ppm)
   :type  xTau: array [time, lev, lat, lon]
   :param nTau: Dust or ice NUMBER mixing ratio (None)
   :type  nTau: array [time, lev, lat, lon]
   :param T: Temperature (K)
   :type  T: array [time, lev, lat, lon]
   :return: ``Vg`` Dust sedimentation rate (m/s)
   :rtype:  array [time, lev, lat, lon]
   :raises ValueError: If the input dimensions are not compatible
   :raises TypeError: If the input types are not compatible
   :raises Exception: If any other error occurs
   :raises RuntimeError: If the sedimentation rate calculation fails


.. py:function:: compute_WMFF(MF, rho, lev, interp_type)

   Calculate the zonal or meridional wave-mean flow forcing.

   Calculation::

       ax = -1/rho d(rho u'w')/dz
       ay = -1/rho d(rho v'w')/dz

   If interp_type == ``pstd``, then::

       [du/dz = (du/dp).(dp/dz)] > [du/dz = -rho*g * (du/dp)]

   where::

       dp/dz = -rho*g
       [du/dz = (du/dp).(-rho*g)] > [du/dz = -rho*g * (du/dp)]

   :param MF: Zonal/meridional momentum flux (J/kg)
   :type  MF: array [time, lev, lat, lon]
   :param rho: Atmospheric density (kg/m^3)
   :type  rho: array [time, lev, lat, lon]
   :param lev: Array for the vertical grid (zagl, zstd, pstd, or pfull)
   :type  lev: array [lev]
   :param interp_type: The vertical grid type (``zagl``, ``zstd``,
       ``pstd``, or ``pfull``)
   :type  interp_type: str
   :return: The zonal or meridional wave-mean flow forcing (m/s2)
   :rtype:  array [time, lev, lat, lon]
   :raises ValueError: If the input dimensions are not compatible
   :raises TypeError: If the input types are not compatible
   :raises Exception: If any other error occurs
   :raises RuntimeError: If the wave-mean flow forcing calculation fails
   :raises ZeroDivisionError: If rho is zero


.. py:function:: compute_dustref_per_pa(dustref, delp)

   Computes visible dust opacity per Pascal from dustref and delp. 
   dustref is visible dust opacity per level (model layer).

   opacity per Pa = opacity per layer / layer thickness in Pa::

       dustref_per_pa = dustref/delp

   [Courtney Batterson, 2025]

   :param dustref: Visible dust opacity [op/model layer]
   :type  dustref: array [time, lev, lat, lon]
   :param delp: Layer thickness [Pa]
   :type  delp: array [time, lev, lat, lon]
   :return: `dustref_per_pa` Visible dust opacity [op/Pa]
   :rtype:  array [time, lev, lat, lon]
   :raises ValueError: If the input dimensions are not compatible
   :raises TypeError: If the input types are not compatible
   :raises Exception: If any other error occurs


.. py:function:: compute_dustref_per_z(dustref, delz)

   Computes visible dust opacity per kilometer from dustref and delz. 
   dustref is visible dust opacity per level (model layer).

   opacity per km = opacity per layer / layer thickness in m * 1000::

       dustref_per_z = dustref/delz*1000

   [Courtney Batterson, 2025]

   :param dustref: Visible dust opacity [op/model layer]
   :type  dustref: array [time, lev, lat, lon]
   :param delz: Layer thickness [m]
   :type  delz: array [time, lev, lat, lon]
   :return: `dustref_per_z` Visible dust opacity [op/km]
   :rtype:  array [time, lev, lat, lon]
   :raises ValueError: If the input dimensions are not compatible
   :raises TypeError: If the input types are not compatible
   :raises Exception: If any other error occurs


.. py:function:: compute_mmr(xTau, temp, lev, const, f_type)

   Compute the dust or ice mixing ratio.

   Adapted from Heavens et al. (2011) observations from MCS (JGR).
   [Courtney Batterson, 2023]

   :param xTau: Dust or ice extinction rate (km-1)
   :type  xTau: array [time, lev, lat, lon]
   :param temp: Temperature (K)
   :type  temp: array [time, lev, lat, lon]
   :param lev: Vertical coordinate (e.g., pstd) (e.g., Pa)
   :type  lev: array [lev]
   :param const: Dust or ice constant
   :type  const: array
   :param f_type: The FV3 file type: diurn, daily, or average
   :type  f_stype: str
   :return: ``q``, Dust or ice mass mixing ratio (ppm)
   :rtype:  array [time, lev, lat, lon]
   :raises ValueError: If the input dimensions are not compatible
   :raises TypeError: If the input types are not compatible
   :raises Exception: If any other error occurs
   :raises RuntimeError: If the mixing ratio calculation fails


.. py:function:: compute_p_3D(ps, ak, bk, shape_out)

   Compute the 3D pressure at layer midpoints.

   :param ps: Surface pressure (Pa)
   :type  ps: array [time, lat, lon]
   :param ak: Vertical coordinate pressure value (Pa)
   :type  ak: array [phalf]
   :param bk: Vertical coordinate sigma value (None)
   :type  bk: array [phalf]
   :param shape_out: Determines how to handle the dimensions of p_3D.
       If ``len(time) = 1`` (one timestep), ``p_3D`` is returned as
       [1, lev, lat, lon] as opposed to [lev, lat, lon]
   :type  shape_out: float
   :return: ``p_3D`` The full 3D pressure array (Pa)
   :rtype:  array [time, lev, lat, lon]
   :raises ValueError: If the input dimensions are not compatible
   :raises TypeError: If the input types are not compatible
   :raises Exception: If any other error occurs
   :raises RuntimeError: If the pressure calculation fails


.. py:function:: compute_rho(p_3D, temp)

   Compute density.

   :param p_3D: Pressure (Pa)
   :type  p_3D: array [time, lev, lat, lon]
   :param temp: Temperature (K)
   :type  temp: array [time, lev, lat, lon]
   :return: Density (kg/m^3)
   :rtype:  array [time, lev, lat, lon]
   :raises ValueError: If the input dimensions are not compatible
   :raises TypeError: If the input types are not compatible
   :raises Exception: If any other error occurs


.. py:function:: compute_scorer(N, ucomp, zfull)

   Calculate the Scorer wavelength.

   :param N: Brunt Vaisala freqency (rad/s)
   :type  N: float [time, lev, lat, lon]
   :param ucomp: Zonal wind (m/s)
   :type  ucomp: array [time, lev, lat, lon]
   :param zfull: Altitude above ground level at the layer midpoint (m)
   :type  zfull: array [time, lev, lat, lon]
   :return: ``scorer_wl`` Scorer horizontal wavelength (m)
   :rtype:  array [time, lev, lat, lon]
   :raises ValueError: If the input dimensions are not compatible
   :raises TypeError: If the input types are not compatible
   :raises Exception: If any other error occurs


.. py:function:: compute_theta(p_3D, ps, T, f_type)

   Compute the potential temperature.

   :param p_3D: The full 3D pressure array (Pa)
   :type  p_3D: array [time, lev, lat, lon]
   :param ps: Surface pressure (Pa)
   :type  ps: array [time, lat, lon]
   :param T: Temperature (K)
   :type  T: array [time, lev, lat, lon]
   :param f_type: The FV3 file type: diurn, daily, or average
   :type  f_type: str
   :return: Potential temperature (K)
   :rtype:  array [time, lev, lat, lon]
   :raises ValueError: If the input dimensions are not compatible
   :raises TypeError: If the input types are not compatible
   :raises Exception: If any other error occurs


.. py:function:: compute_w(rho, omega)

   Compute the vertical wind using the omega equation.

   Under hydrostatic balance, omega is proportional to the vertical
   wind velocity (``w``)::

       omega = dp/dt = (dp/dz)(dz/dt) = (dp/dz) * w

   Under hydrostatic equilibrium::

       dp/dz = -rho * g

   So ``omega`` can be calculated as::

       omega = -rho * g * w

   :param rho: Atmospheric density (kg/m^3)
   :type  rho: array [time, lev, lat, lon]
   :param omega: Rate of change in pressure at layer midpoint (Pa/s)
   :type  omega: array [time, lev, lat, lon]
   :return: vertical wind (m/s)
   :rtype:  array [time, lev, lat, lon]
   :raises ValueError: If the input dimensions are not compatible
   :raises TypeError: If the input types are not compatible
   :raises Exception: If any other error occurs
   :raises RuntimeError: If the vertical wind calculation fails
   :raises ZeroDivisionError: If rho or omega is zero
   :raises OverflowError: If the calculation results in an overflow
   :raises Exception: If any other error occurs


.. py:function:: compute_w_net(Vg, wvar)

   Computes the net vertical wind.

   w_net = vertical wind (w) - sedimentation rate (``Vg_sed``)::

       w_net = w - Vg_sed

   [Courtney Batterson, 2023]

   :param Vg: Dust sedimentation rate (m/s)
   :type  Vg: array [time, lev, lat, lon]
   :param wvar: Vertical wind (m/s)
   :type  wvar: array [time, lev, lat, lon]
   :return: `w_net` Net vertical wind speed (m/s)
   :rtype:  array [time, lev, lat, lon]
   :raises ValueError: If the input dimensions are not compatible
   :raises TypeError: If the input types are not compatible
   :raises Exception: If any other error occurs


.. py:function:: compute_xzTau(q, temp, lev, const, f_type)

   Compute the dust or ice extinction rate.

   Adapted from Heavens et al. (2011) observations from MCS (JGR).
   [Courtney Batterson, 2023]

   :param q: Dust or ice mass mixing ratio (ppm)
   :type  q: array [time, lev, lat, lon]
   :param temp: Temperature (K)
   :type  temp: array [time, lev, lat, lon]
   :param lev: Vertical coordinate (e.g., pstd) (e.g., Pa)
   :type  lev: array [lev]
   :param const: Dust or ice constant
   :type  const: array
   :param f_type: The FV3 file type: diurn, daily, or average
   :type  f_stype: str
   :return: ``xzTau`` Dust or ice extinction rate (km-1)
   :rtype:  array [time, lev, lat, lon]
   :raises ValueError: If the input dimensions are not compatible
   :raises TypeError: If the input types are not compatible
   :raises Exception: If any other error occurs
   :raises RuntimeError: If the extinction rate calculation fails


.. py:function:: compute_zfull(ps, ak, bk, T)

   Calculate the altitude of the layer midpoints above ground level.

   :param ps: Surface pressure (Pa)
   :type  ps: array [time, lat, lon]
   :param ak: Vertical coordinate pressure value (Pa)
   :type  ak: array [phalf]
   :param bk: Vertical coordinate sigma value (None)
   :type  bk: array [phalf]
   :param T: Temperature (K)
   :type  T: array [time, lev, lat, lon]
   :return: ``zfull`` (m)
   :rtype:  array [time, lev, lat, lon]
   :raises ValueError: If the input dimensions are not compatible
   :raises TypeError: If the input types are not compatible
   :raises Exception: If any other error occurs


.. py:function:: compute_zhalf(ps, ak, bk, T)

   Calculate the altitude of the layer interfaces above ground level.

   :param ps: Surface pressure (Pa)
   :type  ps: array [time, lat, lon]
   :param ak: Vertical coordinate pressure value (Pa)
   :type  ak: array [phalf]
   :param bk: Vertical coordinate sigma value (None)
   :type  bk: array [phalf]
   :param T: Temperature (K)
   :type  T: array [time, lev, lat, lon]
   :return: ``zhalf`` (m)
   :rtype:  array [time, lev, lat, lon]
   :raises ValueError: If the input dimensions are not compatible
   :raises TypeError: If the input types are not compatible
   :raises Exception: If any other error occurs


.. py:function:: debug_wrapper(func)

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


.. py:function:: ensure_file_closed(filepath, delay=0.5)

   Try to ensure a file is not being accessed by the system.

   This is especially helpful for Windows environments.

   :param filepath: Path to the file
   :param delay: Delay in seconds to wait for handles to release
   :return: None
   :rtype:  None
   :raises FileNotFoundError: If the file does not exist
   :raises OSError: If the file is locked or cannot be accessed
   :raises Exception: If any other error occurs
   :raises TypeError: If the filepath is not a string
   :raises ValueError: If the filepath is empty
   :raises RuntimeError: If the file cannot be closed


.. py:function:: force_close_netcdf_files(file_or_dir, delay=1.0)

   Aggressively try to ensure netCDF files are closed on Windows systems.

   :param file_or_dir: Path to the file or directory to process
   :param delay: Delay in seconds after forcing closure
   :return: None
   :rtype:  None
   :raises FileNotFoundError: If the file or directory does not exist
   :raises OSError: If the file is locked or cannot be accessed
   :raises Exception: If any other error occurs
   :raises TypeError: If the file_or_dir is not a string
   :raises ValueError: If the file_or_dir is empty
   :raises RuntimeError: If the file or directory cannot be processed
   :raises ImportError: If the netCDF4 module is not available
   :raises AttributeError: If the netCDF4 module does not have the required attributes
   :raises Exception: If any other error occurs


.. py:function:: get_existing_var_name(var_name, file_vars)

   Get the actual variable name that exists in the file.

   Considers alternative naming conventions.

   :param var_name: Variable name to check
   :type  var_name: str
   :param file_vars: Set of variable names in the file
   :type  file_vars: set
   :return: Actual variable name in the file
   :rtype:  str
   :raises ValueError: If the input dimensions are not compatible
   :raises TypeError: If the input types are not compatible
   :raises Exception: If any other error occurs


.. py:function:: main()

   Main function for variable manipulations in NetCDF files.

   This function performs a sequence of operations on one or more
   NetCDF files, as specified by command-line arguments. The operations
   include removing variables, extracting variables, adding variables,
   vertical differentiation, zonal detrending, opacity conversions,
   column integration, and editing variable metadata or values.

   Workflow:
       - Iterates over all input NetCDF files.
       - For each file, performs the following operations as requested
         by arguments:
           * Remove specified variables and update the file.
           * Extract specified variables into a new file.
           * Add new variables using provided methods.
           * Compute vertical derivatives of variables with respect to
             height or pressure.
           * Remove zonal mean (detrend) from specified variables.
           * Convert variables between dp/dz and dz/dp representations.
           * Perform column integration of variables.
           * Edit variable metadata (name, long_name, units) or scale
             values.

   Arguments:
       args: Namespace
           Parsed command-line arguments specifying which operations
           to perform and their parameters.
       master_list: list
           List of available variables and their properties (used for
           adding variables).
       debug: bool
           If True, prints detailed error messages and stack traces.
   Notes:
       - Handles both Unix and Windows file operations for safe file
       replacement.
       - Uses helper functions for NetCDF file manipulation, variable
       existence checks, and error handling.
       - Assumes global constants and utility functions (e.g., Dataset,
       Ncdf, check_file_tape, etc.) are defined elsewhere.
       - Uses global variables lev_T and lev_T_out for axis
       manipulation in vertical operations.

   Raises:
       Exceptions are caught and logged for each operation; files are
       cleaned up on error.


.. py:function:: patched_print_message(self, message, file=None)

   Patched version of _print_message that handles Unicode encoding errors.

   :param self: The ArgumentParser instance
   :param message: The message to print
   :param file: The file to print to (default is sys.stdout)
   :type  file: file-like object
   :return: None
   :rtype:  None
   :raises UnicodeEncodeError: If the message cannot be encoded
   :raises TypeError: If the message is not a string
   :raises ValueError: If the message is empty
   :raises Exception: If any other error occurs
   :raises RuntimeError: If the message cannot be printed


.. py:function:: process_add_variables(file_name, add_list, master_list, debug=False)

   Process a list of variables to add.

   Dependent variables are added in the correct order.
   If a variable is already in the file, it is skipped.
   If a variable cannot be added, an error message is printed.

   :param file_name: Input file path
   :param add_list: List of variables to add
   :param master_list: Dictionary of supported variables and their dependencies
   :param debug: Whether to show debug information
   :type  debug: bool
   :raises ValueError: If the input dimensions are not compatible
   :raises TypeError: If the input types are not compatible
   :raises Exception: If any other error occurs
   :raises RuntimeError: If the variable cannot be added


.. py:function:: safe_copy_replace(src_file, dst_file, max_attempts=5, delay=1.0)

   Windows-specific approach to copy file contents and replace destination.

   This avoids move operations which are more likely to fail with locking

   :param src_file: Source file path
   :param dst_file: Destination file path
   :param max_attempts: Maximum number of retry attempts
   :param delay: Base delay between attempts (increases with retries)
   :return: True if successful, False otherwise
   :rtype:  bool
   :raises FileNotFoundError: If the source file does not exist
   :raises OSError: If the file is locked or cannot be accessed
   :raises Exception: If any other error occurs
   :raises TypeError: If the src_file or dst_file is not a string
   :raises ValueError: If the src_file or dst_file is empty
   :raises RuntimeError: If the file cannot be copied or replaced


.. py:function:: safe_move_file(src_file, dst_file, max_attempts=5, delay=1)

   Safely move a file with retries for Windows file locking issues.

   :param src_file: Source file path
   :param dst_file: Destination file path
   :param max_attempts: Number of attempts to make
   :param delay: Delay between attempts in seconds
   :return: True if successful, False otherwise
   :rtype:  bool
   :raises FileNotFoundError: If the source file does not exist
   :raises OSError: If the file is locked or cannot be accessed
   :raises Exception: If any other error occurs
   :raises TypeError: If the src_file or dst_file is not a string
   :raises ValueError: If the src_file or dst_file is empty
   :raises RuntimeError: If the file cannot be moved


.. py:function:: safe_print(text)

   Print text safely, handling encoding issues on Windows.

   :param text: Text to print
   :type  text: str
   :return: None
   :rtype:  None
   :raises UnicodeEncodeError: If the text cannot be encoded
   :raises TypeError: If the text is not a string
   :raises ValueError: If the text is empty
   :raises Exception: If any other error occurs
   :raises RuntimeError: If the text cannot be printed


.. py:function:: safe_remove_file(filepath, max_attempts=5, delay=1)

   Safely remove a file with retries for Windows file locking issues

   :param filepath: Path to the file to remove
   :param max_attempts: Number of attempts to make
   :param delay: Delay between attempts in seconds
   :return: True if successful, False otherwise
   :rtype:  bool
   :raises FileNotFoundError: If the file does not exist
   :raises OSError: If the file is locked or cannot be accessed
   :raises Exception: If any other error occurs
   :raises TypeError: If the filepath is not a string
   :raises ValueError: If the filepath is empty
   :raises RuntimeError: If the file cannot be removed


.. py:data:: C_dst

   

.. py:data:: C_ice

   

.. py:data:: Cp
   :value: 735.0

   

.. py:data:: Kb

   

.. py:data:: M_co2
   :value: 0.044

   

.. py:data:: N
   :value: 0.01

   

.. py:data:: Na

   

.. py:data:: Qext_dst
   :value: 0.35

   

.. py:data:: Qext_ice
   :value: 0.773

   

.. py:data:: R
   :value: 8.314

   

.. py:data:: Rd
   :value: 192.0

   

.. py:data:: Reff_dst
   :value: 1.06

   

.. py:data:: Reff_ice
   :value: 1.41

   

.. py:data:: S0
   :value: 222

   

.. py:data:: T0
   :value: 273.15

   

.. py:data:: Tpole
   :value: 150.0

   

.. py:data:: amu

   

.. py:data:: amu_co2
   :value: 44.0

   

.. py:data:: args

   

.. py:data:: cap_str
   :value: ' (derived w/CAP)'

   

.. py:data:: debug

   

.. py:data:: exit_code

   

.. py:data:: filepath

   

.. py:data:: fill_value
   :value: 0.0

   

.. py:data:: g
   :value: 3.72

   

.. py:data:: mass_co2

   

.. py:data:: master_list

   

.. py:data:: n0

   

.. py:data:: original_print_message

   

.. py:data:: parser

   

.. py:data:: psrf
   :value: 610.0

   

.. py:data:: rgas
   :value: 189.0

   

.. py:data:: rho_air

   

.. py:data:: rho_dst
   :value: 2500.0

   

.. py:data:: rho_ice
   :value: 900

   

.. py:data:: sigma
   :value: 0.63676

   

