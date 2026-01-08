:py:mod:`amescap.FV3_utils`
===========================

.. py:module:: amescap.FV3_utils

.. autoapi-nested-parse::

   FV3_utils contains internal Functions for processing data in MGCM
   output files such as vertical interpolation.

   These functions can be used on their own outside of CAP if they are
   imported as a module::

       from /u/path/FV3_utils import fms_press_calc

   Third-party Requirements:

       * ``numpy``
       * ``warnings``
       * ``scipy``



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   amescap.FV3_utils.MGStau_ls_lat
   amescap.FV3_utils.MGSzmax_ls_lat
   amescap.FV3_utils.UT_LTtxt
   amescap.FV3_utils.add_cyclic
   amescap.FV3_utils.alt_KM
   amescap.FV3_utils.area_weights_deg
   amescap.FV3_utils.areo_avg
   amescap.FV3_utils.axis_interp
   amescap.FV3_utils.azimuth2cart
   amescap.FV3_utils.broadcast
   amescap.FV3_utils.cart_to_azimut_TR
   amescap.FV3_utils.compute_uneven_sigma
   amescap.FV3_utils.daily_to_average
   amescap.FV3_utils.daily_to_diurn
   amescap.FV3_utils.dvar_dh
   amescap.FV3_utils.expand_index
   amescap.FV3_utils.find_n
   amescap.FV3_utils.find_n0
   amescap.FV3_utils.fms_Z_calc
   amescap.FV3_utils.fms_press_calc
   amescap.FV3_utils.frontogenesis
   amescap.FV3_utils.gauss_profile
   amescap.FV3_utils.get_trend_2D
   amescap.FV3_utils.interp_KDTree
   amescap.FV3_utils.layers_mid_point_to_boundary
   amescap.FV3_utils.lin_interp
   amescap.FV3_utils.lon180_to_360
   amescap.FV3_utils.lon360_to_180
   amescap.FV3_utils.ls2sol
   amescap.FV3_utils.mass_stream
   amescap.FV3_utils.mollweide2cart
   amescap.FV3_utils.ortho2cart
   amescap.FV3_utils.polar2XYZ
   amescap.FV3_utils.polar_warming
   amescap.FV3_utils.press_pa
   amescap.FV3_utils.press_to_alt_atmosphere_Mars
   amescap.FV3_utils.ref_atmosphere_Mars_PTD
   amescap.FV3_utils.regression_2D
   amescap.FV3_utils.robin2cart
   amescap.FV3_utils.second_hhmmss
   amescap.FV3_utils.shiftgrid_180_to_360
   amescap.FV3_utils.shiftgrid_360_to_180
   amescap.FV3_utils.sol2ls
   amescap.FV3_utils.sol_hhmmss
   amescap.FV3_utils.spherical_curl
   amescap.FV3_utils.spherical_div
   amescap.FV3_utils.swinbank
   amescap.FV3_utils.time_shift_calc
   amescap.FV3_utils.transition
   amescap.FV3_utils.vinterp
   amescap.FV3_utils.vw_from_MSF
   amescap.FV3_utils.zonal_detrend



.. py:function:: MGStau_ls_lat(ls, lat)

   Return the max altitude for the dust from "MGS scenario" from
   Montmessin et al. (2004), Origin and role of water ice clouds in
   the Martian water cycle as inferred from a general circulation model

   :param ls: solar longitude [°]
   :type  ls: array
   :param lat : latitude [°]
   :type  lat: array
   :return: top altitude for the dust [km]


.. py:function:: MGSzmax_ls_lat(ls, lat)

   Return the max altitude for the dust from "MGS scenario" from
   Montmessin et al. (2004), Origin and role of water ice clouds in
   the Martian water cycle as inferred from a general circulation model

   :param ls: solar longitude [°]
   :type  ls: array
   :param lat : latitude [°]
   :type  lat: array
   :return: top altitude for the dust [km]


.. py:function:: UT_LTtxt(UT_sol, lon_180=0.0, roundmin=None)

   Returns the time in HH:MM:SS at a certain longitude.

   :param time_sol: the time in sols
   :type  time_sol: float
   :param lon_180: the center longitude in -180/180 coordinates.
       Increments by 1hr every 15°
   :type  lon_180: float
   :param roundmin: round to the nearest X minute. Typical values are
       ``roundmin = 1, 15, 60``
   :type  roundmin: int

   .. note::
       If ``roundmin`` is requested, seconds are not shown


.. py:function:: add_cyclic(data, lon)

   Add a cyclic (overlapping) point to a 2D array. Useful for azimuth
   and orthographic projections.

   :param data: variable of size ``[nlat, nlon]``
   :type  data: array
   :param lon: longitudes
   :type  lon: array
   :return: a 2D array of size ``[nlat, nlon+1]`` with last column
       identical to the 1st; and a 1D array of longitudes
       size [nlon+1] where the last element is ``lon[-1] + dlon``


.. py:function:: alt_KM(press, scale_height_KM=8.0, reference_press=610.0)

   Gives the approximate altitude [km] for a given pressure

   :param press: the pressure [Pa]
   :type  press: 1D array
   :param scale_height_KM: scale height [km] (default is 8 km, an
       isothermal at 155K)
   :type  scale_height_KM: float
   :param reference_press: reference surface pressure [Pa] (default is
       610 Pa)
   :type  reference_press: float
   :return: ``z_KM`` the equivalent altitude for that pressure [km]

   .. note::
       Scale height is ``H = rT/g``


.. py:function:: area_weights_deg(var_shape, lat_c, axis=-2)

   Returns weights scaled so that np.mean(var*W) gives an area-weighted 
   average. This works because grid cells near the poles have smaller 
   areas than those at the equator, so they should contribute less to 
   a global average.

   Expected dimensions are:

   [lat] ``axis`` not needed
   [lat, lon] ``axis = -2`` or ``axis = 0``
   [time, lat, lon] ``axis = -2`` or ``axis = 1``
   [time, lev, lat, lon] ``axis = -2`` or ``axis = 2``
   [time, time_of_day_24, lat, lon] ``axis = -2`` or ``axis = 2``
   [time, time_of_day_24, lev, lat, lon] ``axis = -2`` or ``axis = 3``

   Because ``dlat`` is computed as ``lat_c[1]-lat_c[0]``, ``lat_c``
   may be truncated on either end (e.g., ``lat = [-20 ..., 0... 50]``)
   but must be continuous.

   :param var_shape: the shape/dimensions of your data array
   :type  var_shape: tuple
   :param lat_c: latitude cell centers in degrees [°]
   :type  lat_c: float
   :param axis: which dimension contains latitude, default: 2nd-to-last
   :type  axis: int
   :return: ``W`` weights for the variable ready for standard
       averaging as ``np.mean(var*W)`` [condensed form] or
       ``np.average(var, weights=W)`` [expanded form]

   .. note::
       Given a variable var:

       ``var = [v1, v2, ...vn]``

       The regular average is

       ``AVG = (v1 + v2 + ... vn) / N``

       and the weighted average is

       ``AVG_W = (v1*w1 + v2*w2 + ... vn*wn) / (w1 + w2 + ...wn)``

       This function returns

       ``W = [w1, w2, ... , wn]*N / (w1 + w2 + ...wn)``

       Therfore taking a regular average of (``var*W``) with
       ``np.mean(var*W)`` or ``np.average(var, weights=W)``

       returns the weighted average of the variable. Use
       ``np.average(var, weights=W, axis = X)`` to average over a
       specific axis.


.. py:function:: areo_avg(VAR, areo, Ls_target, Ls_angle, symmetric=True)

   Return a value average over a central solar longitude

   EX::

       ``Ls_target = 90.``
       ``Ls_angle = 10.``

   Nominally, the time average is done over solar longitudes
   ``85 < Ls_target < 95`` (10°).

   If ``symmetric = True`` and the input data range is Ls = 88-100°
   then ``88 < Ls_target < 92`` (4°, symmetric)

   If ``symmetric = False`` and the input data range is Ls = 88-100°
   then ``88 < Ls_target < 95`` (7°, assymetric)

   :param VAR: a variable with ``time`` in the 1st dimension
   :type  VAR: ND array
   :param areo: solar longitude of the input variable (0-720)
   :type  areo: 1D array
   :param Ls_target: central solar longitude of interest
   :type  Ls_target: float
   :param Ls_angle: requested window angle centered at ``Ls_target``
   :type  Ls_angle: float
   :param symmetric: If ``True`` and the requested window is out of range,
       ``Ls_angle`` is reduced. If False, the time average is performed
       on the data available
   :type  symmetric: bool (defaults to True)
   :return: the variable averaged over solar longitudes
       ``Ls_target-Ls_angle/2`` to ``Ls_target+Ls_angle/2``

   .. note::
       The routine can bin data from muliples Mars years


.. py:function:: axis_interp(var_IN, x, xi, axis, reverse_input=False, type_int='lin', modulo=None)

   One dimensional linear/logarithmic interpolation along one axis.

   :param var_IN: Variable on a REGULAR grid (e.g.,
       ``[lev, lat, lon]`` or ``[time, lev, lat, lon]``)
   :type  var_IN: ND array
   :param x: Original position array (e.g., ``time``)
   :type  x: 1D array
   :param xi: Target array to interpolate the array on
   :type  xi: 1D array
   :param axis: Position of the interpolation axis (e.g., ``0`` for a
       temporal interpolation on ``[time, lev, lat, lon]``)
   :type  axis: int
   :param reverse_input: Reverse input arrays (e.g., if
       ``zfull(0)``= 120 km, ``zfull(N)``= 0 km, which is typical)
   :type  reverse_input: bool
   :param type_int: "log" for logarithmic (typically pressure),
       "lin" for linear
   :type  type_int: str
   :param modulo: For "lin" interpolation only, use cyclic input
       (e.g., when using ``modulo = 24`` for time of day, 23.5 and
       00 am are considered 30 min apart, not 23.5 hr apart)
   :type  modulo: float
   :return: ``VAR_OUT`` interpolated data on the requested axis

   .. note::
       This routine is similar but simpler than the vertical
       interpolation ``vinterp()`` as the interpolation axis is
       assumed to be fully defined by a 1D array such as ``time``,
       ``pstd`` or ``zstd`` rather than 3D arrays like ``pfull`` or
       ``zfull``.
       For lon/lat interpolation, consider using ``interp_KDTree()``.

   Calculation::

       X_OUT = Xn*A + (1-A)*Xn + 1
       with ``A = log(xi/xn + 1) / log(xn/xn + 1)`` in "log" mode
       or ``A = (xi-xn + 1)/(xn-xn + 1)`` in "lin" mode


.. py:function:: azimuth2cart(LAT, LON, lat0, lon0=0)

   Azimuthal equidistant projection. Converts from latitude-longitude
   to cartesian coordinates.

   :param LAT: latitudes[°] size [nlat]
   :type  LAT: 1D or 2D array
   :param LON: longitudes [°] size [nlon]
   :type  LON: 1D or 2D array
   :param lat0: latitude coordinate of the pole
   :type  lat0: float
   :param lon0: longitude coordinate of the pole
   :type  lon0: float
   :return: the cartesian coordinates for the latitudes and longitudes


.. py:function:: broadcast(var_1D, shape_out, axis)

   Broadcast a 1D array based on a variable's dimensions

   :param var_1D: variable (e.g., ``lat`` size = 36, or ``time``
       size = 133)
   :type  var_1D: 1D array
   :param shape_out: broadcasting shape (e.g.,
       ``temp.shape = [133, lev, 36, lon]``)
   :type  shape_out: list
   :return: (ND array) broadcasted variables (e.g., size =
       ``[1,36,1,1]`` for ``lat`` or ``[133,1,1,1]`` for ``time``)


.. py:function:: cart_to_azimut_TR(u, v, mode='from')

   Convert cartesian coordinates or wind vectors to radians using azimuth angle.

   :param x: the cartesian coordinate
   :type  x: 1D array
   :param y: the cartesian coordinate
   :type  y: 1D array
   :param mode: "to" for the direction that the vector is pointing,
       "from" for the direction from which the vector is coming
   :type  mode: str
   :return: ``Theta`` [°] and ``R`` the polar coordinates


.. py:function:: compute_uneven_sigma(num_levels, N_scale_heights, surf_res, exponent, zero_top)

   Construct an initial array of sigma based on the number of levels
   and an exponent

   :param num_levels: the number of levels
   :type  num_levels: float
   :param N_scale_heights: the number of scale heights to the top of
       the model (e.g., ``N_scale_heights`` = 12.5 ~102 km assuming an
       8 km scale height)
   :type  N_scale_heights: float
   :param surf_res: the resolution at the surface
   :type  surf_res: float
   :param exponent: an exponent to increase the thickness of the levels
   :type  exponent: float
   :param zero_top: if True, force the top pressure boundary
       (in N = 0) to 0 Pa
   :type  zero_top: bool
   :return: an array of sigma layers


.. py:function:: daily_to_average(varIN, dt_in, nday=5, trim=True)

   Bin a variable from an ``atmos_daily`` file format to the
   ``atmos_average`` file format.

   :param varIN: variable with ``time`` dimension first (e.g.,
       ``ts[time, lat, lon]``)
   :type  varIN: ND array
   :param dt_in: delta of time betwen timesteps in sols (e.g.,
       ``dt_in = time[1] - time[0]``)
   :type  dt_in: float
   :param nday: bining period in sols, default is 5 sols
   :type  nday: int
   :param trim: whether to discard any leftover data at the end of file
       before binning
   :type  trim: bool
   :return: the variable bin over ``nday``

   .. note::
       If ``varIN[time, lat, lon]`` from ``atmos_daily`` is
       ``[160, 48, 96]`` and has 4 timesteps per day (every 6 hours),
       then the resulting variable for ``nday = 5`` is
       ``varOUT(160/(4*5), 48, 96) = varOUT(8, 48, 96)``

   .. note::
       If the daily file has 668 sols, then there are
       ``133 x 5 + 3`` sols leftover. If ``trim = True``, then the
       time is 133 and last 3 sols the are discarded. If
       ``trim = False``, the time is 134 and last bin contains only
       3 sols of data.


.. py:function:: daily_to_diurn(varIN, time_in)

   Bin a variable from an ``atmos_daily`` file into the
   ``atmos_diurn`` format.

   :param varIN: variable with time dimension first (e.g.,
       ``[time, lat, lon]``)
   :type  varIN: ND array
   :param time_in: time array in sols. Only the first N elements
       are actually required if saving memory is important
   :type  time_in: ND array
   :return: the variable binned in the ``atmos_diurn`` format
       (``[time, time_of_day, lat, lon]``) and the time of day array
       [hr]

   .. note::
       If ``varIN[time, lat, lon]`` from ``atmos_daily`` is
       ``[40, 48, 96]`` and has 4 timestep per day (every 6 hours),
       then the resulting variable is
       ``varOUT[10, 4, 48, 96] = [time, time_of_day, lat, lon]`` and
       ``tod = [0., 6., 12., 18.]``.

   .. note::
       Since the time dimension is first, the output variables
       may be passed to the ``daily_to_average()`` function for
       further binning.


.. py:function:: dvar_dh(arr, h=None)

   Differentiate an array ``A[dim1, dim2, dim3...]`` w.r.t ``h``. The
   differentiated dimension must be the first dimension.

   EX: Compute ``dT/dz`` where ``T[time, lev, lat, lon]`` is the
   temperature and ``Zkm`` is the array of  level heights [km].

   First, transpose ``T`` so the vertical dimension comes first:
   ``T[lev, time, lat, lon]``.

   Then transpose back to get ``dTdz[time, lev, lat, lon]``::

       dTdz = dvar_dh(t.transpose([1, 0, 2, 3]),
                      Zkm).transpose([1, 0, 2, 3])

   If ``h`` is 1D, then ``h``and ``dim1`` must have the same length

   If ``h`` is 2D, 3D or 4D, then ``arr`` and ``h`` must have the
   same shape

   :param arr: variable
   :type  arr: ND array
   :param h: the dimension (``Z``, ``P``, ``lat``, ``lon``)
   :type  h: str
   :return: d_arr: the array differentiated w.r.t ``h``, e.g., d(array)/dh


.. py:function:: expand_index(Nindex, VAR_shape_axis_FIRST, axis_list)

   Repeat interpolation indices along an axis.

   :param Nindex: Interpolation indices, size is (``n_axis``,
       ``Nfull = [time, lat, lon]``)
   :type  Nindex: idx
   :param VAR_shape_axis_FIRST: Shape for the variable to interpolate
       with interpolation axis first (e.g., ``[tod, time, lev, lat, lon]``)
   :type  VAR_shape_axis_FIRST: tuple
   :param axis_list: Position or list of positions for axis to insert
       (e.g., ``2`` for ``lev`` in ``[tod, time, lev, lat, lon]``, ``[2, 4]``
       for ``lev`` and ``lon``). The axis positions are those for the final
       shape (``VAR_shape_axis_FIRST``) and must be INCREASING
   :type  axis_list: int or list
   :return: ``LFULL`` a 2D array (size ``n_axis``,
       ``NfFULL = [time, lev, lat, lon]``) with the indices expanded
       along the ``lev`` dimension and flattened

   .. note::
       Example of application:
       Observational time of day may be the same at all vertical levels
       so the interpolation of a 5D variable ``[tod, time, lev, lat, lon]``
       only requires the interpolation indices for ``[tod, time, lat, lon]``.
       This routine expands the indices from ``[tod, time, lat, lon]`` to
       ``[tod, time, lev, lat, lon]`` with ``Nfull = [time, lev, lat, lon]``
       for use in interpolation.


.. py:function:: find_n(X_IN, X_OUT, reverse_input=False, modulo=None)

   Maps the closest index from a 1D input array to a ND output array
   just below the input values.

   :param X_IN: Source level [Pa] or [m]
   :type  X_IN: float or 1D array
   :param X_OUT: Desired pressure [Pa] or altitude [m] at layer
       midpoints. Level dimension is FIRST
   :type  X_OUT: array
   :param reverse_input: If input array is decreasing (e.g., if z(0)
       = 120 km, z(N) = 0 km, which is typical, or if data is
       p(0) = 1000 Pa, p(N) = 0 Pa, which is uncommon)
   :type  reverse_input: bool
   :return: The index for the level(s) where the pressure < ``plev``


.. py:function:: find_n0(Lfull_IN, Llev_OUT, reverse_input=False)

   Return the index for the level(s) just below ``Llev_OUT``.
   This assumes ``Lfull_IN`` is increasing in the array
   (e.g., ``p(0) = 0``, ``p(N) = 1000`` [Pa]).

   :param Lfull_IN: Input pressure [Pa] or altitude [m] at layer
       midpoints. ``Level`` dimension is FIRST
   :type  Lfull_IN: array
   :param Llev_OUT: Desired level type for interpolation [Pa] or [m]
   :type  Llev_OUT: float or 1D array
   :param reverse_input: Reverse array (e.g., if ``z(0) = 120 km``,
       ``z(N) = 0km`` -- which is typical -- or if input data is
       ``p(0) = 1000Pa``, ``p(N) = 0Pa``)
   :type  reverse_input: bool
   :return: ``n`` index for the level(s) where the pressure is just
       below ``plev``

   .. note::
       If ``Lfull_IN`` is a 1D array and ``Llev_OUT`` is a float
       then ``n`` is a float.

   .. note::
       If ``Lfull_IN`` is ND ``[lev, time, lat, lon]`` and
       ``Llev_OUT`` is a 1D array of size ``klev`` then ``n`` is an
       array of size ``[klev, Ndim]`` with ``Ndim = [time, lat, lon]``


.. py:function:: fms_Z_calc(psfc, ak, bk, T, topo=0.0, lev_type='full')

   Returns the 3D altitude field [m] AGL (or above aeroid).

   :param psfc: The surface pressure [Pa] or array of surface
       pressures (1D, 2D, or 3D)
   :type  psfc: array
   :param ak: 1st vertical coordinate parameter
   :type  ak: array
   :param bk: 2nd vertical coordinate parameter
   :type  bk: array
   :param T: The air temperature profile. 1D array (for a single grid
       point), ND array with VERTICAL AXIS FIRST
   :type  T: 1D array or ND array
   :param topo: The surface elevation. Same dimension as ``psfc``.
       If None is provided, AGL is returned
   :type  topo: array
   :param lev_type: "full" (layer midpoint) or "half" (layer
       interfaces). Defaults to "full"
   :type  lev_type: str
   :return: The layer altitude at the full level ``Z_f(:, :, Nk-1)``
       or half-level ``Z_h(:, :, Nk)`` [m]. ``Z_f`` and ``Z_h`` are
       AGL if ``topo = None``. ``Z_f`` and ``Z_h`` are above aeroid
       if topography is not None.

   Calculation::

       --- 0 --- TOP        ========  z_half
       --- 1 ---
                           --------  z_full

                           ========  z_half
       ---Nk-1---          --------  z_full
       --- Nk --- SFC      ========  z_half
                           / / / / /

   .. note::
       Expands to the time dimension using::

           topo = np.repeat(zsurf[np.newaxis, :], ps.shape[0], axis = 0)

   Calculation is derived from
       ``./atmos_cubed_sphere_mars/Mars_phys.F90``::

           # (dp/dz = -rho g) => (dz = dp/(-rho g)) and
           # (rho = p/(r T)) => (dz = rT/g * (-dp/p))

           # Define log-pressure (``u``) as:
           u = ln(p)

           # Then:
           du = {du/dp}*dp = {1/p)*dp} = dp/p

           # Finally, ``dz`` for the half-layers:
           (dz = rT/g * -(du)) => (dz = rT/g * (+dp/p))
           # with ``N`` layers defined from top to bottom.

   Z_half calculation::

       # Hydrostatic relation within the layer > (P(k+1)/P(k) =
       # exp(-DZ(k)/H))

       # layer thickness:
       DZ(k) = rT/g * -(du)

       # previous layer altitude + thickness of layer:
       Z_h k) = Z_h(k+1)  +DZ_h(h)

   Z_full calculation::

       # previous altitude + half the thickness of previous layer and
       # half of current layer
       Z_f(k) = Z_f(k+1) + (0.5 DZ(k) + 0.5 DZ(k+1))

       # Add ``+0.5 DZ(k)-0.5 DZ(k)=0`` and re-organiz the equation
       Z_f(k) = Z_f(k+1) + DZ(k) + 0.5 (DZ(k+1) - DZ(k))
       Z_f(k) = Z_h(k+1) + 0.5 (DZ(k+1) - DZ(k))

   The specific heat ratio:
   ``γ = cp/cv (cv = cp-R)`` => ``γ = cp/(cp-R)`` Also ``(γ-1)/γ=R/cp``

   The dry adiabatic lapse rate:
   ``Γ = g/cp`` => ``Γ = (gγ)/R``

   The isentropic relation:
   ``T2 = T1(p2/p1)**(R/cp)``

   Therefore::

       line 1) =====Thalf=====zhalf[k]          line 2)                                   line 3)                                    line 4) -----Tfull-----zfull[k]     \ T(z)= To-Γ (z-zo)
       line 5)                                      line 6)                                       line 7) =====Thalf=====zhalf[k+1]      
   Line 1: T_half[k+1]/Tfull[k] = (p_half[k+1]/p_full[k])**(R/Cp)

   Line 4: From the lapse rate, assume T decreases linearly within the
   layer so ``T_half[k+1] = T_full[k] + Γ(Z_full[k]-Z_half[k+1])``
   and (``Tfull < Thalf`` and ``Γ > 0``)

   Line 7: ``Z_full[k] = Z_half[k] + (T_half[k+1]-T_full[k])/Γ``
   Pulling out ``Tfull`` from above equation and using ``Γ = (gγ)/R``::

       Z_full[k] = (Z_half[k+1] + (R Tfull[k]) / (gγ)(T_half[k+1]
       / T_full[k] - 1))

   Using the isentropic relation above::

       Z_full = (Z_half[k+1] + (R Tfull[k]) / (gγ)(p_half[k+1]
       / p_full[k])**(R/Cp)-1))


.. py:function:: fms_press_calc(psfc, ak, bk, lev_type='full')

   Returns the 3D pressure field from the surface pressure and the
   ak/bk coefficients.

   :param psfc: the surface pressure [Pa] or an array of surface
       pressures (1D, 2D, or 3D if time dimension)
   :type  psfc: array
   :param ak: 1st vertical coordinate parameter
   :type  ak: array
   :param bk: 2nd vertical coordinate parameter
   :type  bk: array:
   :param lev_type: "full" (layer midpoints) or "half"
       (layer interfaces). Defaults to "full."
   :type  lev_type: str
   :return: the 3D pressure field at the full levels
       ``PRESS_f(Nk-1:,:,:)`` or half-levels ``PRESS_h(Nk,:,:,)`` [Pa]

   Calculation::

       --- 0 --- TOP        ========  p_half
       --- 1 ---
                            --------  p_full

                            ========  p_half
       ---Nk-1---           --------  p_full
       --- Nk --- SFC       ========  p_half
                           / / / / /

   .. note::
       Some literature uses pk (pressure) instead of ak with
       ``p3d = ps * bk + P_ref * ak`` instead of ``p3d = ps * bk + ak``


.. py:function:: frontogenesis(U, V, theta, lon_deg, lat_deg, R=3400 * 1000.0, spacing='varying')

   Compute the frontogenesis (local change in potential temperature
   gradient near a front) following Richter et al. 2010: Toward a
   Physically Based Gravity Wave Source Parameterization in a General
   Circulation Model, JAS 67.

   We have ``Fn = 1/2 D(Del Theta)^2/Dt`` [K/m/s]

   :param U: wind field with ``lat`` SECOND TO LAST and ``lon`` as last
       dimensions (e.g., ``[lat, lon]`` or ``[time, lev, lat, lon``]
       etc.)
   :type  U: array
   :param V: wind field with ``lat`` SECOND TO LAST and ``lon`` as last
       dimensions (e.g., ``[lat, lon]`` or ``[time, lev, lat, lon``]
       etc.)
   :type  V: array
   :param theta: potential temperature [K]
   :type  theta: array
   :param lon_deg: longitude [°] (2D if irregularly-spaced)
   :type  lon_deg: 1D array
   :param lat_deg: latitude [°] (2D if irregularly-spaced)
   :type  lat_deg: 1D array
   :param R: planetary radius [m]
   :type  R: float
   :param spacing: when ``lon`` and ``lat`` are 1D arrays, using
       spacing = "varying" differentiates latitude and longitude. When
       spacing = "regular", ``dx = lon[1]-lon[0]``,
       `` dy=lat[1]-lat[0]`` and the ``numpy.gradient()`` method are
       used
   :type  spacing: str (defaults to "varying")
   :return: the frontogenesis field [m-1]


.. py:function:: gauss_profile(x, alpha, x0=0.0)

   Return Gaussian line shape at x. This can be used to generate a
   bell-shaped mountain.


.. py:function:: get_trend_2D(VAR, LON, LAT, type_trend='wmean')

   Extract spatial trends from the data. The output can be directly
   subtracted from the original field.

   :param VAR: Variable for decomposition. ``lat`` is SECOND TO LAST
       and ``lon`` is LAST  (e.g., ``[time, lat, lon]`` or
       ``[time, lev, lat, lon]``)
   :type  VAR: ND array
   :param LON: longitude coordinates
   :type  LON: 2D array
   :param LAT: latitude coordinates
   :type  LAT: 2D array
   :param type_trend: type of averaging to perform:
       "mean" - use a constant average over all lat/lon
       "wmean" - use a area-weighted average over all lat/lon
       "zonal" - detrend over the zonal axis only
       "2D" - use a 2D planar regression (not area-weighted)
   :type  type_trend: str
   :return: the trend, same size as ``VAR``


.. py:function:: interp_KDTree(var_IN, lat_IN, lon_IN, lat_OUT, lon_OUT, N_nearest=10)

   Inverse distance-weighted interpolation using nearest neighboor for
   ND variables. Alex Kling, May 2021

   :param var_IN: ND variable to regrid (e.g., ``[lev, lat, lon]``,
       ``[time, lev, lat, lon]`` with ``[lat, lon]`` dimensions LAST
       [°])
   :type  var_IN: ND array
   :param lat_IN: latitude [°] (``LAT[y, x]`` array for
       irregular grids)
   :type  lat_IN: 1D or 2D array
   :param lon_IN: latitude [°] (``LAT[y, x]`` array for
       irregular grids)
   :type  lon_IN: 1D or 2D array
   :param lat_OUT: latitude [°] for the TARGET grid structure
       (or ``LAT1[y,x]`` for irregular grids)
   :type  lat_OUT: 1D or 2D array
   :param lon_OUT: longitude [°] for the TARGET grid structure
       (or ``LON1[y,x]`` for irregular grids)
   :type  lon_OUT: 1D or 2D array
   :param N_nearest: number of nearest neighbours for the search
   :type  N_nearest: int
   :return: ``VAR_OUT`` interpolated data on the target grid

   .. note::
       This implementation is much FASTER than ``griddata`` and
       it supports unstructured grids like an MGCM tile.
       The nearest neighbour interpolation is only done on the lon/lat
       axis (not level). Although this interpolation works well on the
       3D field [x, y, z], this is typically not what is expected. In
       a 4°x4° run, the closest points in all directions (N, E, S, W)
       on the target grid are 100's of km away while the closest
       points in the vertical are a few 10's -100's meter in the PBL.
       This would result in excessive weighting in the vertical.


.. py:function:: layers_mid_point_to_boundary(pfull, sfc_val)

   A general description for the layer boundaries is::

       p_half = ps*bk + pk

   This routine converts the coordinate of the layer MIDPOINTS,
   ``p_full`` or ``bk``, into the coordinate of the layer BOUNDARIES
   ``p_half``. The surface value must be provided.

   :param p_full: Pressure/sigma values for the layer MIDPOINTS,
       INCREASING with ``N`` (e.g., [0.01 -> 720] or [0.001 -> 1])
   :type  p_full: 1D array
   :param sfc_val: The surface value for the lowest layer's boundary
       ``p_half[N]`` (e.g., ``sfc_val`` = 720 Pa or ``sfc_val`` = 1 in
       sigma coordinates)
   :type  sfc_val: float
   :return: ``p_half`` the pressure at the layer boundaries
       (size = ``N+1``)

   Structure::

       --- 0 --- TOP   ========  p_half
       --- 1 ---
                       --------  p_full

                       ========  p_half
       ---Nk-1---      --------  p_full
       --- Nk --- SFC  ========  p_half
                       / / / / /

   We have::

       pfull[N] = ((phalf[N]-phalf[N-1]) / np.log(phalf[N]/phalf[N-1]))
       => phalf[N-1] - pfull[N] log(phalf[N-1])
       = phalf[N] - pfull[N] log(phalf[N])

   We want to solve for ``phalf[N-1] = X``::

       v                v                             v
       X      - pfull[N]       log(X)   =             B

   ``=> X= -pfull[N] W{-exp(-B/pfull[N])/pfull[N]}``

   with ``B = phalf[N] - pfull[N] log(phalf[N])`` (known at N) and

   ``W`` is the product-log (Lambert) function.

   This was tested on an L30 simulation: The values of ``phalf`` are
   reconstructed from ``pfull`` with a max error of:

   ``100*(phalf - phalf_reconstruct)/phalf < 0.4%`` at the top.


.. py:function:: lin_interp(X_in, X_ref, Y_ref)

   Simple linear interpolation with no dependance on scipy

   :param X_in: input values
   :type  X_in: float or array
   :param X_ref x values
   :type  X_ref: array
   :param Y_ref y values
   :type  Y_ref: array
   :return: y value linearly interpolated at ``X_in``


.. py:function:: lon180_to_360(lon)

   Transform a float or an array from the -180/180 coordinate system
   to 0-360

   :param lon: longitudes in the -180/180 coordinate system
   :type  lon: float, 1D array, or 2D array
   :return: the equivalent longitudes in the 0-360 coordinate system


.. py:function:: lon360_to_180(lon)

   Transform a float or an array from the 0-360 coordinate system to
       -180/180.

   :param lon: longitudes in the 0-360 coordinate system
   :type  lon: float, 1D array, or 2D array
   :return: the equivalent longitudes in the -180/180 coordinate system


.. py:function:: ls2sol(Ls_in)

   Ls to sol converter.

   :param Ls_in: solar longitudes (0-360...720)
   :type  Ls_in: float or 1D array
   :return: the corresponding sol number

   .. note::
       This function simply uses a numerical solver on the
       ``sol2ls()`` function.


.. py:function:: mass_stream(v_avg, lat, level, type='pstd', psfc=700, H=8000.0, factor=1e-08)

   Compute the mass stream function::

                               P
                               ⌠
       Ph i= (2 pi a) cos(lat)/g ⎮vz_tavg dp
                               ⌡
                               p_top

   :param v_avg: zonal wind [m/s] with ``lev`` dimension FIRST and
       ``lat`` dimension SECOND (e.g., ``[pstd, lat]``,
       ``[pstd, lat, lon]`` or ``[pstd, lat, lon, time]``)
   :type  v_avg: ND array
   :param lat: latitudes [°]
   :type  lat: 1D array
   :param level: interpolated layers [Pa] or [m]
   :type  level: 1D array
   :param type: interpolation type (``pstd``, ``zstd`` or ``zagl``)
   :type  type: str
   :param psfc: reference surface pressure [Pa]
   :type  psfc: float
   :param H: reference scale height [m] when pressures are used
   :type  H: float
   :param factor: normalize the mass stream function by a factor, use
       ``factor = 1`` for [kg/s]
   :type  factor: int
   :return: ``MSF`` the meridional mass stream function (in
       ``factor * [kg/s]``)

   .. note::
       This routine allows the time and zonal averages to be
       computed before OR after the MSF calculation.

   .. note::
       The expressions for MSF use log(pressure) Z coordinates,
       which integrate better numerically.

       With ``p = p_sfc exp(-Z/H)`` and ``Z = H log(p_sfc/p)``
       then ``dp = -p_sfc/H exp(-Z/H) dZ`` and we have::

                                           Z_top
                                           ⌠
           Phi = +(2pi a)cos(lat)psfc/(gH) ⎮v_rmv exp(-Z/H)dZ
                                           ⌡
                                           Z
       With ``p = p_sfc exp(-Z/H)``

       The integral is calculated using trapezoidal rule::

               n
               ⌠
           .g. ⌡ f(z)dz = (Zn-Zn-1){f(Zn) + f(Zn-1)}/2
             n-1


.. py:function:: mollweide2cart(LAT, LON)

   Mollweide projection. Converts from latitude-longitude to
   cartesian coordinates.

   :param LAT: latitudes[°] size [nlat]
   :type  LAT: 1D or 2D array
   :param LON: longitudes [°] size [nlon]
   :type  LON: 1D or 2D array
   :param lat0: latitude coordinate of the pole
   :type  lat0: float
   :param lon0: longitude coordinate of the pole
   :type  lon0: float
   :return: the cartesian coordinates for the latitudes and longitudes


.. py:function:: ortho2cart(LAT, LON, lat0, lon0=0)

   Orthographic projection. Converts from latitude-longitude to
   cartesian coordinates.

   :param LAT: latitudes[°] size [nlat]
   :type  LAT: 1D or 2D array
   :param LON: longitudes [°] size [nlon]
   :type  LON: 1D or 2D array
   :param lat0: latitude coordinate of the pole
   :type  lat0: float
   :param lon0: longitude coordinate of the pole
   :type  lon0: float
   :return: the cartesian coordinates for the latitudes and longitudes;
       and a mask (NaN array) that hides the back side of the planet


.. py:function:: polar2XYZ(lon, lat, alt, Re=3400 * 10**3)

   Spherical to cartesian coordinate transformation.

   :param lon: Longitude in radians
   :type  lon: ND array
   :param lat: Latitude in radians
   :type  lat: ND array
   :param alt: Altitude [m]
   :type  alt: ND array
   :param Re: Planetary radius [m], defaults to 3400*10^3
   :type  Re: float
   :return: ``X``, ``Y``, ``Z`` in cartesian coordinates [m]

   .. note::
       This is a classic polar coordinate system with
       ``colatitude = pi/2 - lat`` where ``cos(colat) = sin(lat)``


.. py:function:: polar_warming(T, lat, outside_range=np.nan)

   Return the polar warming, following McDunn et al. 2013:
   Characterization of middle-atmosphere polar warming at Mars, JGR
   Alex Kling

   :param T: temperature with the lat dimension FIRST (transpose as
       needed)
   :type  T: ND array
   :param lat: latitude array
   :type  lat: 1D array
   :param outside_range: values to set the polar warming to when
       outside pf the range. Default = NaN but 0 may be desirable.
   :type  outside_range: float
   :return: The polar warming [K]

   .. note::
       ``polar_warming()`` concatenates the results from both
       hemispheres obtained from the nested function
       ``PW_half_hemisphere()``


.. py:function:: press_pa(alt_KM, scale_height_KM=8.0, reference_press=610.0)

   Gives the approximate altitude [km] for a given pressure

   :param alt_KM: the altitude [km]
   :type  alt_KM: 1D array
   :param scale_height_KM: scale height [km] (default is 8 km, an
       isothermal at 155K)
   :type  scale_height_KM: float
   :param reference_press: reference surface pressure [Pa] (default is
       610 Pa)
   :type  reference_press: float
   :return: ``press_pa`` the equivalent pressure at that altitude [Pa]

   .. note::
       Scale height is ``H = rT/g``


.. py:function:: press_to_alt_atmosphere_Mars(Pi)

   Return the altitude [m] as a function of pressure from the
   analytical calculation above.

   :param Pi: input pressure [Pa] (<= 610 Pa)
   :type  Pi: float or 1D array
   :return: the corresponding altitude [m] (float or 1D array)


.. py:function:: ref_atmosphere_Mars_PTD(Zi)

   Analytical atmospheric model for Martian pressure, temperature, and
   density. Alex Kling, June 2021

   :param Zi: input altitude [m] (must be >= 0)
   :type  Zi: float or 1D array
   :return: tuple of corresponding pressure [Pa], temperature [K],
   and density [kg/m3] floats or arrays

   .. note::
       This model was obtained by fitting globally and annually
       averaged reference temperature profiles derived from the Legacy
       GCM, MCS observations, and Mars Climate Database.

       The temperature fit was constructed using quadratic temperature
       ``T(z) = T0 + gam(z-z0) + a*(z-z0)^2`` over 4 segments (0>57 km,
       57>110 km, 110>120 km and 120>300 km).

       From the ground to 120 km, the pressure is obtained by
       integrating (analytically) the hydrostatic equation:

       ``dp/dz=-g. p/(rT)`` with ``T(z) = T0 + gam(z-z0) + a*(z-z0)^2``

       Above ~120 km, ``P = P0 exp(-(z-z0)g/rT)`` is not a good
       approximation as the fluid is in molecula regime. For those
       altitudes, we provide a fit in the form of
       ``P = P0 exp(-az-bz^2)`` based on diurnal average of the MCD
       database at lat = 0, Ls = 150.


.. py:function:: regression_2D(X, Y, VAR, order=1)

   Linear and quadratic regression on the plane.

   :param X: first coordinate
   :type  X: 2D array
   :param Y: second coordinate
   :type  Y: 2D array
   :param VAR: variable of the same size as X
   :type  VAR: 2D array
   :param order: 1 (linear) or 2 (quadratic)
   :type  order: int

   .. note::
       When ``order = 1``, the equation is: ``aX + bY + C = Z``.
       When ``order = 2``, the equation is:
       ``aX^2 + 2bX*Y + cY^2 + 2dX + 2eY + f = Z``

   For the linear case::, ``ax + by + c = z`` is re-written as
   ``A X = b`` with::

               |x0   y0   1|        |a      |z0
           A = |x1   y1   1|    X = |b   b= |
               |      ...  |        |c      |...
               |xn   yn   1|                |zn

                   [n,3]           [3]       [n]

   The least-squares regression provides the solution that that
   minimizes ``||b – A x||^2``


.. py:function:: robin2cart(LAT, LON)

   Robinson projection. Converts from latitude-longitude to cartesian
   coordinates.

   :param LAT: latitudes[°] size [nlat]
   :type  LAT: 1D or 2D array
   :param LON: longitudes [°] size [nlon]
   :type  LON: 1D or 2D array
   :param lat0: latitude coordinate of the pole
   :type  lat0: float
   :param lon0: longitude coordinate of the pole
   :type  lon0: float
   :return: the cartesian coordinates for the latitudes and longitudes


.. py:function:: second_hhmmss(seconds, lon_180=0.0)

   Given the time [sec], return local true solar time at a
   certain longitude.

   :param seconds: the time [sec]
   :type  seconds: float
   :param lon_180: the longitude in -180/180 coordinate
   :type  lon_180: float
   :return: the local time [float] or a tuple (hours, minutes, seconds)


.. py:function:: shiftgrid_180_to_360(lon, data)

   This function shifts ND data from a -180/180 to a 0-360 grid.

   :param lon: longitudes in the 0-360 coordinate system
   :type  lon: 1D array
   :param data: variable with ``lon`` in the last dimension
   :type  data: ND array
   :return: shifted data


.. py:function:: shiftgrid_360_to_180(lon, data)

   This function shifts ND data from a 0-360 to a -180/180 grid.

   :param lon: longitudes in the 0-360 coordinate system
   :type  lon: 1D array
   :param data: variable with ``lon`` in the last dimension
   :type  data: ND array
   :return: shifted data

   .. note::
       Use ``np.ma.hstack`` instead of ``np.hstack`` to keep the
       masked array properties


.. py:function:: sol2ls(jld, cumulative=False)

   Return the solar longitude (Ls) as a function of the sol number.
   Sol=0 is the spring equinox.

   :param jld: sol number after perihelion
   :type  jld: float or 1D array
   :param cumulative: if True, result is cumulative
       (Ls=0-360, 360-720 etc..)
   :type  cumulative: bool
   :return: the corresponding solar longitude


.. py:function:: sol_hhmmss(time_sol, lon_180=0.0)

   Given the time in days, return return local true solar time at a
   certain longitude.

   :param time_sol: the time in sols
   :type  seconds: float
   :param lon_180: the longitude in -180/180 coordinate
   :type  lon_180: float
   :return: the local time [float] or a tuple (hours, minutes, seconds)


.. py:function:: spherical_curl(U, V, lon_deg, lat_deg, R=3400 * 1000.0, spacing='varying')

   Compute the vertical component of the relative vorticity using
   finite difference::

       curl = dv/dx -du/dy
            = 1/(r cos lat)[d(v)/dlon + d(u(cos lat)/dlat]

   :param U: wind field with ``lat`` SECOND TO LAST and ``lon`` as last
       dimensions (e.g., ``[lat, lon]`` or ``[time, lev, lat, lon``]
       etc.)
   :type  U: array
   :param V: wind field with ``lat`` SECOND TO LAST and ``lon`` as last
       dimensions (e.g., ``[lat, lon]`` or ``[time, lev, lat, lon``]
       etc.)
   :type  V: array
   :param lon_deg: longitude [°] (2D if irregularly-spaced)
   :type  lon_deg: 1D array
   :param lat_deg: latitude [°] (2D if irregularly-spaced)
   :type  lat_deg: 1D array
   :param R: planetary radius [m]
   :type  R: float
   :param spacing: when ``lon`` and ``lat`` are 1D arrays, using
       spacing = "varying" differentiates latitude and longitude. When
       spacing = "regular", ``dx = lon[1]-lon[0]``,
       `` dy=lat[1]-lat[0]`` and the ``numpy.gradient()`` method are
       used
   :type  spacing: str (defaults to "varying")
   :return: the vorticity of the wind field [m-1]


.. py:function:: spherical_div(U, V, lon_deg, lat_deg, R=3400 * 1000.0, spacing='varying')

   Compute the divergence of the wind fields using finite difference::

       div = du/dx + dv/dy
       -> = 1/(r cos lat)[d(u)/dlon + d(v cos lat)/dlat]

   :param U: wind field with ``lat`` SECOND TO LAST and ``lon`` as last
       dimensions (e.g., ``[lat, lon]`` or ``[time, lev, lat, lon``]
       etc.)
   :type  U: array
   :param V: wind field with ``lat`` SECOND TO LAST and ``lon`` as last
       dimensions (e.g., ``[lat, lon]`` or ``[time, lev, lat, lon``]
       etc.)
   :type  V: array
   :param lon_deg: longitude [°] (2D if irregularly-spaced)
   :type  lon_deg: 1D array
   :param lat_deg: latitude [°] (2D if irregularly-spaced)
   :type  lat_deg: 1D array
   :param R: planetary radius [m]
   :type  R: float
   :param spacing: when ``lon`` and ``lat`` are 1D arrays, using
       spacing = "varying" differentiates latitude and longitude. When
       spacing = "regular", ``dx = lon[1]-lon[0]``,
       `` dy=lat[1]-lat[0]`` and the ``numpy.gradient()`` method are
       used
   :type  spacing: str (defaults to "varying")
   :return: the horizonal divergence of the wind field [m-1]


.. py:function:: swinbank(plev, psfc, ptrans=1.0)

   Compute ``ak`` and ``bk`` values with a transition based on Swinbank

   :param plev: the pressure levels [Pa]
   :type  plev: 1D array
   :param psfc: the surface pressure [Pa]
   :type  psfc: 1D array
   :param ptrans: the transition pressure [Pa]
   :type  ptrans: 1D array
   :return: the coefficients for the new layers


.. py:function:: time_shift_calc(var_in, lon, tod, target_times=None)

   Conversion to uniform local time.

   Mars rotates approx. 14.6° lon per Mars-hour (360° ÷ 24.6 hr)
   Each 14.6° shift in lon represents a 1-hour shift in local time
   This code uses the more precise calculation: lon_shift = 24.0 * lon / 360.

   :param var_in: variable to be shifted. Assume ``lon`` is the first
       dimension and ``tod`` is the last dimension
   :type  var_in: ND array
   :param lon: longitude
   :type  lon: 1D array
   :param tod: ``time_of_day`` index from the input file
   :type  tod: 1D array
   :param target_times: local time(s) [hr] to shift to (e.g., ``"3. 15."``)
   :type  target_times: float (optional)
   :return: the array shifted to uniform local time

   .. note::
       If ``target_times`` is not specified, the file is interpolated
       on the same ``tod`` as the input


.. py:function:: transition(pfull, p_sigma=0.1, p_press=0.05)

   Return the transition factor to construct ``ak`` and ``bk``

   In the MGCM code, the full pressures are computed from::

                      del(phalf)
        pfull = -----------------------------
                log(phalf(k+1/2)/phalf(k-1/2))

   :param pfull: the pressure [Pa]
   :type  pfull: 1D array
   :param p_sigma: the pressure level where the vertical grid starts
       transitioning from sigma to pressure
   :type  p_sigma: float
   :param p_press: the pressure level above which the vertical grid is
       pure (constant) pressure
   :type  p_press: float
   :return: the transition factor. = 1 for pure sigma, = 0 for pure
       pressure and =0-1 for the transition


.. py:function:: vinterp(varIN, Lfull, Llev, type_int='log', reverse_input=False, masktop=True, index=None)

   Vertical linear or logarithmic interpolation for pressure or altitude.

   :param varIN: Variable to interpolate (VERTICAL AXIS FIRST)
   :type  varIN: ND array
   :param Lfull: Pressure [Pa] or altitude [m] at full layers, same
       dimensions as ``varIN``
   :type  Lfull: array
   :param Llev: Desired level for interpolation [Pa] or [m]. May be
       increasing or decreasing as the output levels are processed one
       at a time
   :type  Llev: 1D array
   :param type_int: "log" for logarithmic (typically pressure),
       "lin" for linear (typically altitude)
   :type  type_int: str
   :param reverse_input: Reverse input arrays. e.g., if
       ``zfull[0]`` = 120 km then ``zfull[N]`` = 0km (typical) or if
       input data is ``pfull[0]``=1000 Pa, ``pfull[N]``=0 Pa
   :type  reverse_input: bool
   :param masktop: Set to NaN values if above the model top
   :type  masktop: bool
   :param index: Indices for the interpolation, already processed as
       ``[klev, Ndim]``. Indices calculated if not provided
   :type  index: None or array
   :return: ``varOUT`` variable interpolated on the ``Llev`` pressure
       or altitude levels

   .. note::
       This interpolation assumes pressure decreases with height::

           --  0  -- TOP  [0 Pa]   : [120 km]| X_OUT = Xn*A + (1-A)*Xn + 1
           --  1  --               :         |
                                   :         |
           --  n  -- pn   [30 Pa]  : [800 m] | Xn
                                   :         |
           --  k  -- Llev [100 Pa] : [500 m] | X_OUT
           -- n+1 -- pn+1 [200 Pa] : [200 m] | Xn+1

           -- SFC --
           / / / / / /

       with ``A = log(Llev/pn + 1) / log(pn/pn + 1)`` in "log" mode
       or ``A = (zlev-zn + 1) / (zn-zn + 1)`` in "lin" mode


.. py:function:: vw_from_MSF(msf, lat, lev, ztype='pstd', norm=True, psfc=700, H=8000.0)

   Return the V and W components of the circulation from the mass
   stream function.

   :param msf: the mass stream function with ``lev`` SECOND TO
       LAST and the ``lat`` dimension LAST (e.g., ``[lev, lat]``,
       ``[time, lev, lat]``, ``[time, lon, lev, lat]``)
   :type  msf: ND array
   :param lat: latitude [°]
   :type  lat: 1D array
   :param lev: level [Pa] or [m] (``pstd``, ``zagl``, ``zstd``)
   :type  lev: 1D array
   :param ztype: Use ``pstd`` for pressure so vertical
       differentation is done in log space.
   :type  ztype: str
   :param norm: if True, normalize ``lat`` and ``lev`` before
       differentiating to avoid having to rescale manually the
       vectors in quiver plots
   :type  norm: bool
   :param psfc: surface pressure for pseudo-height when
       ``ztype = pstd``
   :type  psfc: float
   :param H: scale height for pseudo-height when ``ztype = pstd``
   :type  H: float
   :return: the meditional and altitude components of the mass stream
       function for plotting as a quiver or streamlines.

   .. note::
       The components are:
       ``[v]=  g/(2 pi cos(lat)) dphi/dz``
       ``[w]= -g/(2 pi cos(lat)) dphi/dlat``


.. py:function:: zonal_detrend(VAR)

   Substract the zonal average mean value from a field.

   :param VAR: variable with detrending dimension last
   :type  VAR: ND array
   :return: detrented field (same size as input)

   .. note::
       ``RuntimeWarnings`` are expected if the slice contains
       only NaNs which is the case below the surface and above the
       model top in the interpolated files. This routine disables such
       warnings temporarily.


