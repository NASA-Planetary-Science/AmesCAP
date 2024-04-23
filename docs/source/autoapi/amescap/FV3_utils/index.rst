:orphan:

:py:mod:`amescap.FV3_utils`
===========================

.. py:module:: amescap.FV3_utils


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   amescap.FV3_utils.fms_press_calc
   amescap.FV3_utils.fms_Z_calc
   amescap.FV3_utils.find_n0
   amescap.FV3_utils.find_n
   amescap.FV3_utils.expand_index
   amescap.FV3_utils.vinterp
   amescap.FV3_utils.axis_interp
   amescap.FV3_utils.layers_mid_point_to_boundary
   amescap.FV3_utils.polar2XYZ
   amescap.FV3_utils.interp_KDTree
   amescap.FV3_utils.cart_to_azimut_TR
   amescap.FV3_utils.sfc_area_deg
   amescap.FV3_utils.area_meridional_cells_deg
   amescap.FV3_utils.area_weights_deg
   amescap.FV3_utils.areo_avg
   amescap.FV3_utils.mass_stream
   amescap.FV3_utils.vw_from_MSF
   amescap.FV3_utils.alt_KM
   amescap.FV3_utils.press_pa
   amescap.FV3_utils.shiftgrid_360_to_180
   amescap.FV3_utils.shiftgrid_180_to_360
   amescap.FV3_utils.second_hhmmss
   amescap.FV3_utils.sol_hhmmss
   amescap.FV3_utils.UT_LTtxt
   amescap.FV3_utils.dvar_dh
   amescap.FV3_utils.zonal_detrend
   amescap.FV3_utils.get_trend_2D
   amescap.FV3_utils.regression_2D
   amescap.FV3_utils.daily_to_average
   amescap.FV3_utils.daily_to_diurn
   amescap.FV3_utils.gauss_profile
   amescap.FV3_utils.compute_uneven_sigma
   amescap.FV3_utils.transition
   amescap.FV3_utils.swinbank
   amescap.FV3_utils.polar_warming
   amescap.FV3_utils.tshift
   amescap.FV3_utils.lin_interp
   amescap.FV3_utils.add_cyclic
   amescap.FV3_utils.spherical_div
   amescap.FV3_utils.spherical_curl
   amescap.FV3_utils.frontogenesis
   amescap.FV3_utils.MGSzmax_ls_lat
   amescap.FV3_utils.MGStau_ls_lat
   amescap.FV3_utils.broadcast
   amescap.FV3_utils.ref_atmosphere_Mars_PTD
   amescap.FV3_utils.press_to_alt_atmosphere_Mars
   amescap.FV3_utils.azimuth2cart
   amescap.FV3_utils.ortho2cart
   amescap.FV3_utils.mollweide2cart
   amescap.FV3_utils.robin2cart
   amescap.FV3_utils.sol2ls
   amescap.FV3_utils.ls2sol



.. py:function:: fms_press_calc(psfc, ak, bk, lev_type='full')

   Returns the 3D pressure field from the surface pressure and the ak/bk coefficients.

   Args:
       psfc:       the surface pressure in [Pa] or
                   an array of surface pressures (1D, 2D, or 3D if time dimension)
       ak:         1st vertical coordinate parameter
       bk:         2nd vertical coordinate parameter
       lev_type:   "full" (layer midpoints) or  "half" (layer interfaces).
                   Defaults to "full."
   Returns:
       The 3D pressure field at the full PRESS_f(Nk-1:,:,:) or half-levels PRESS_h(Nk,:,:,) in [Pa]
   --- 0 --- TOP        ========  p_half
   --- 1 ---
                        --------  p_full

                        ========  p_half
   ---Nk-1---           --------  p_full
   --- Nk --- SFC       ========  p_half
                       / / / / /

   *NOTE*
       Some literature uses pk (pressure) instead of ak.
       With (p3d = ps * bk + P_ref * ak) vs the current (p3d = ps * bk + ak)


.. py:function:: fms_Z_calc(psfc, ak, bk, T, topo=0.0, lev_type='full')

   Returns the 3D altitude field in [m] AGL or above aeroid.

   Args:
       psfc:       The surface pressure [Pa] or array of surface pressures (1D, 2D, or 3D).
       ak:         1st vertical coordinate parameter.
       bk:         2nd vertical coordinate parameter.
       T:          The air temperature profile. 1D array (for a single grid point),
                   N-dimensional array with VERTICAL AXIS FIRST.
       topo:       The surface elevation. Same dimension as 'psfc'. If None is provided,
                   AGL is returned.
       lev_type:   "full" (layer midpoint) or "half" (layer interfaces). Defaults to "full".
   Returns:
       The layer altitude at the full level Z_f(:, :, Nk-1) or half-level Z_h(:, :, Nk) in [m].
       Z_f and Z_h are AGL if topo = None.
       Z_f and Z_h are above aeroid if topo is provided.

   --- 0 --- TOP        ========  z_half
   --- 1 ---
                        --------  z_full

                        ========  z_half
   ---Nk-1---           --------  z_full
   --- Nk --- SFC       ========  z_half
                       / / / / /

   *NOTE*
       Expands to the time dimension using:
           topo = np.repeat(zsurf[np.newaxis, :], ps.shape[0], axis = 0)

   *NOTE*
       Expands topo to the time dimension using:
           topo = np.repeat(zsurf[np.newaxis, :], ps.shape[0], axis = 0)

       Calculation is derived from ./atmos_cubed_sphere_mars/Mars_phys.F90:
           (dp/dz = -rho g) => (dz = dp/(-rho g)) and
           (rho= p/(r T)) => (dz=rT/g * (-dp/p))
       Define log-pressure (u) as:
           u = ln(p)
       Then:
           du = {du/dp}*dp = {1/p)*dp} = dp/p
       Finally, dz for the half-layers:
           (dz = rT/g * -(du)) => (dz = rT/g *(+dp/p))
           with N layers defined from top to bottom.

   Z_half calculation:
   ------------------
   Hydrostatic relation within the layer > (P(k+1)/P(k) = exp(-DZ(k)/H))
   > DZ(k) = rT/g * -(du)         (layer thickness)
   > Z_h(k) = Z_h(k+1)  +DZ_h(h)  (previous layer altitude + thickness of layer)

   Z_full calculation:
   ------------------
   Z_f(k) = Z_f(k+1) + (0.5 DZ(k) + 0.5 DZ(k+1)) (previous altitude + half the thickness
                                                  of previous layer and half of current layer)
          = Z_f(k+1) + DZ(k) + 0.5 (DZ(k+1) - DZ(k)) (Added '+0.5 DZ(k)-0.5 DZ(k)=0' and
                                                      re-organized the equation)
          = Z_h(k+1) + 0.5 (DZ(k+1) - DZ(k))

   The  specific heat ratio       γ  = cp/cv (cv = cp-R) => γ = cp/(cp-R) Also (γ-1)/γ=R/cp
   The dry adiabatic lapse rate   Γ  = g/cp => Γ = (gγ)/R
   The isentropic relation        T2 = T1(p2/p1)**(R/cp)

   therefore, T_half[k+1]/Tfull[k] = (p_half[k+1]/p_full[k])**(R/Cp)            =====Thalf=====zhalf[k]                                                                                                                                                                                                                             From the lapse rate, assume T decreases linearly within the layer:           -----Tfull-----zfull[k]     \ T(z)= To-Γ (z-zo)
       T_half[k+1] = T_full[k] + Γ(Z_full[k]-Z_half[k+1])                                                            (Tfull < Thalf and Γ > 0)                                                                                      Z_full[k] = Z_half[k] + (T_half[k+1]-T_full[k])/Γ                        =====Thalf=====zhalf[k+1]          Pulling out Tfull from above equation and using  Γ = (gγ)/R:
       Z_full[k] = Z_half[k+1] + (R Tfull[k])/(gγ)(T_half[k+1]/T_full[k] - 1)
   Using the isentropic relation above:
       Z_full = Z_half[k+1] + (R Tfull[k])/(gγ)(p_half[k+1]/p_full[k])**(R/Cp)-1)


.. py:function:: find_n0(Lfull_IN, Llev_OUT, reverse_input=False)

   Return the index for the level(s) just below Llev_OUT.
   This assumes Lfull_IN is increasing in the array (e.g p(0) = 0Pa, p(N) = 1000Pa)

   Args:
       Lfull_IN (array):               input pressure [pa] or altitude [m] at layer midpoints. 'Level' dimension is FIRST.
       Llev_OUT (float or 1D array):   Desired level type for interpolation [Pa] or [m].
       reverse_input (boolean):        Reverse array (e.g if z(0) = 120 km, z(N) = 0km -- which is typical -- or if
                                       input data is p(0) = 1000Pa, p(N) = 0Pa).
   Returns:
       n:    index for the level(s) where the pressure is just below 'plev'.
   ***NOTE***
       - if Lfull_IN is 1D array and Llev_OUT is a float then n is a float.
       - if Lfull_IN is ND [lev, time, lat, lon] and Llev_OUT is a 1D array of size 'klev' then
         'n' is an array of size [klev, Ndim] with 'Ndim' = (time x lat x lon).


.. py:function:: find_n(X_IN, X_OUT, reverse_input=False, modulo=None)

   Map the closest index from a 1D input array to a ND output array just below the input values.
   Args:
       X_IN (float or 1D array):   source level [Pa] or [m].
       X_OUT (ND array):           desired pressure [pa] or altitude [m] at layer midpoints. 'Level' dimension is FIRST.
       reverse_input (boolean):    if input array is decreasing (e.g if z(0) = 120 km, z(N)=0km -- which is typical -- or
                                   data is p(0) = 1000Pa, p(N) = 0Pa -- which is uncommon in MGCM output
   Returns:
       n:    index for the level(s) where the pressure is just below 'plev'.

       Case 1:      Case 2:      Case 3:        Case 4:
       (ND)   (1D)  (1D)  (1D)   (1D)   (ND)    (ND)        (ND)
      |x|x|         |x|            |x|            |x|x|
      |x|x| > |x|   |x| > |x|      |x| > |x|x|    |x|x|  >  |x|x|
      |x|x|   |x|   |x|   |x|      |x|   |x|x|    |x|x|     |x|x|  (case 4, must have same number of
      |x|x|   |x|   |x|   |x|      |x|   |x|x|    |x|x|     |x|x|  elements along the other dimensions)

      *** Note on cyclic values ***

      *** Note ***
      Cyclic array are handled naturally (e.g. time of day 0.5 ..23.5 > 0.5) or longitudes 0 >... 359 >0
      >>> if first (0) array element is above requested value, (e.g 0.2 is requested from [0.5 1.5... 23.5], n is set to 0-1=-1 which refers to the last element, 23.5 here)
      >>> last element in array is always inferior or equal to selected value: (e.g 23.8 is requested from [0.5 1.5... 23.5], 23.5 will be selected

      Therefore, the cyclic values must therefore be handled during the interpolation but not at this stage.

      


.. py:function:: expand_index(Nindex, VAR_shape_axis_FIRST, axis_list)

   Repeat interpolation indices along an axis.
   Args:
       Nindex: inteprolation indices, size is (n_axis, Nfull= time x lat x lon)
       VAR_shape_axis_FIRST: shape for the variable to interpolate with interpolation axis first e.g. (tod,time,lev,lat,lon)
       axis_list (int or list): position or list of positions for axis to insert, e.g. '2' for LEV in (tod,time,LEV,lat,lon), '[2,4]' for LEV and LON
                                The axis position are those for the final shape (VAR_shape_axis_FIRST) and must be INCREASING
   Returns:
       LFULL: a 2D array size(n_axis, NfFULL= time x LEV x lat x lon) with the indices expended  along the LEV dimensions and flattened
   ***NOTE***
   Example of application:
    Observational time of day  may the same at all vertical levels so the interpolation of a 5D variable
   (tod,time,LEV,lat,lon) only requires the interpolation indices for (tod,time,lat,lon).
   This routines expands the indices from (tod,time,lat,lon) to (tod,time,LEV,lat,lon) with Nfull=time x lev x lat x lon for use in interpolation



.. py:function:: vinterp(varIN, Lfull, Llev, type_int='log', reverse_input=False, masktop=True, index=None)

   Vertical linear or logarithmic interpolation for pressure or altitude.   Alex Kling 5-27-20
   Args:
       varIN: variable to interpolate (N-dimensional array with VERTICAL AXIS FIRST)
       Lfull: pressure [Pa] or altitude [m] at full layers same dimensions as varIN
       Llev : desired level for interpolation as a 1D array in [Pa] or [m] May be either increasing or decreasing as the output levels are processed one at the time.
       reverse_input (boolean) : reverse input arrays, e.g if zfull(0)=120 km, zfull(N)=0km (which is typical) or if your input data is pfull(0)=1000Pa, pfull(N)=0Pa
       type_int : 'log' for logarithmic (typically pressure), 'lin' for linear (typically altitude)
       masktop: set to NaN values if above the model top
       index: indices for the interpolation, already processed as [klev,Ndim]
              Indices will be recalculated if not provided.
   Returns:
       varOUT: variable interpolated on the Llev pressure or altitude levels

   *** IMPORTANT NOTE***
   This interpolation assumes pressure are increasing downward, i.e:

       ---  0  --- TOP   [0 Pa]   : [120 km]|    X_OUT= Xn*A + (1-A)*Xn+1
       ---  1  ---                :         |
                                  :         |
       ---  n  ---  pn   [30 Pa]  : [800 m] | Xn
                                  :         |
   >>> ---  k  ---  Llev [100 Pa] : [500 m] | X_OUT
       --- n+1 ---  pn+1 [200 Pa] : [200 m] | Xn+1

       --- SFC ---
       / / / / / /

   with A = log(Llev/pn+1)/log(pn/pn+1) in 'log' mode
        A =    (zlev-zn+1)/(zn-zn+1)    in 'lin' mode




.. py:function:: axis_interp(var_IN, x, xi, axis, reverse_input=False, type_int='lin', modulo=None)

   One dimensional linear /log interpolation along one axis. [Alex Kling, May 2021]
   Args:
       var_IN (N-D array): N-Dimensional variable, e.g.  [lev,lat,lon],[time,lev,lat,lon] on a REGULAR grid.
       x (1D array)      : original position array (e.g. time)
       xi (1D array)     : target array to interpolate the array on
       axis (int)        : position of  the interpolation axis (e.g. 0 if time interpolation for [time,lev,lat,lon])
       reverse_input (boolean) : reverse input arrays, e.g if zfull(0)=120 km, zfull(N)=0km (which is typical)
       type_int : 'log' for logarithmic (typically pressure), 'lin' for linear
       modulo (float)    : for 'lin' interpolation only, use cyclic input (e.g when using modulo = 24 for time of day, 23.5 and 00am are considered 30 min appart, not 23.5hr)
   Returns:
       VAR_OUT: interpolated data on the requested axis

   ***NOTE***
   > This routine is similar, but simpler, than the vertical interpolation vinterp()  as the  interpolation axis is assumed to be fully defined by a 1D
    array (e.g. 'time', 'pstd' or 'zstd) unlike  pfull or zfull which are 3D arrays.

   > For lon/lat interpolation, you may consider using  interp_KDTree() instead

   We have:

   X_OUT= Xn*A + (1-A)*Xn+1
   with A = log(xi/xn+1)/log(xn/xn+1) in 'log' mode
        A =    (xi-xn+1)/(xn-xn+1)    in 'lin' mode


.. py:function:: layers_mid_point_to_boundary(pfull, sfc_val)

   A general description for the layer boundaries is p_half= ps*bk +pk
   This routine convert  p_full or bk, the coordinate  of the layer MIDPOINTS into the coordinate of the layers
   BOUNDARIES, p_half. The surface value must be provided    [A. Kling, 2022]

   Args:
       p_full : 1D array of presure/sigma values for the layers's MIDPOINTS, INCREASING with N (e.g. [0.01... 720] or [0.001.. 1])
       sfc_val : the surface value for the lowest layer's boundary p_half[N], e.g. sfc_val=720Pa or sfc_val=1. for sigma coordinates
   Returns:
       p_half: the pressure at the layers boundaries, the size is N+1

   ***NOTE***

   --- 0 --- TOP        ========  p_half
   --- 1 ---
                        --------  p_full

                        ========  p_half
   ---Nk-1---           --------  p_full
   --- Nk --- SFC       ========  p_half
                       / / / / /
   We have pfull[N]= (phalf[N]-phalf[N-1])/np.log(phalf[N]/phalf[N-1])

   => phalf[N-1]- pfull[N] log(phalf[N-1])= phalf[N]-pfull[N] log(phalf[N]) . We want to solve for phalf[N-1]=X
         v                v                             v
         X      - pfull[N]       log(X)   =             B

   ==> X= - pfull[N] W{-exp(-B/pfull[N])/pfull[N]}  with B = phalf[N] - pfull[N] log(phalf[N]) (known at N)
   and W the product-log (Lambert)   function

   Though the product-log function is available in python, we use an approximation for portability
   (e.g. Appendix in Kling et al. 2020, Icarus)

   This was tested on a L30 simulation:
   The value of phalf are reconstruted from pfull with a max error of 100*(phalh-phalf_reconstruct)/phalf < 0.4% at the top.


.. py:function:: polar2XYZ(lon, lat, alt, Re=3400 * 10**3)

   Spherical to cartesian coordinates transformation
   Args:
       lon,lat (ND array): longitude and latitude, in [rad]
       alt (ND array): altitude in [m]
   Return:
       X,Y,Z in cartesian coordinates [m]
   ***NOTE***
   This is a classic polar coordinate system with colat = pi/2 -lat,  the colatitude and cos(colat) = sin(lat)


.. py:function:: interp_KDTree(var_IN, lat_IN, lon_IN, lat_OUT, lon_OUT, N_nearest=10)

   Inverse-distance-weighted interpolation using nearest neighboor for ND variables.  [Alex Kling , May 2021]
   Args:
       var_IN: N-Dimensional variable to regrid, e.g.  [lev,lat,lon],[time,lev,lat,lon]... with [lat, lon] dimensions LAST in [deg]
       lat_IN,lon_IN        (1D or 2D):   lat, lon 1D arrays or LAT[y,x] LON[y,x] for irregular grids in [deg]
       lat_OUT,lon_OUT(1D or 2D):lat,lon for the TARGET grid structure , e.g. lat1,lon1 or LAT1[y,x], LON1[y,x] for irregular grids in [deg]
       N_nearest: integer, number of nearest neighbours for the search.
   Returns:
       VAR_OUT: interpolated data on the target grid

   ***NOTE***
   > This implementation is much FASTER than griddata and supports unstructured grids (e.g. FV3 tile)
   > The nearest neighbour interpolation is only done on the lon/lat axis, (not level).  Although this interpolation work well on the 3D field (x,y,z),
   this is typically not what is expected: In a 4°x4° run, the closest points East, West, North and South, on the target grid  are 100's of km away
   while the closest points in the vertical are a few 10's -100's meter in the PBL, which would results in excessive weighting in the vertical.


.. py:function:: cart_to_azimut_TR(u, v, mode='from')

   Convert cartesian coordinates or wind vectors to radian,using azimut angle.

   Args:
       x,y: 1D arrays for the cartesian coordinate
       mode='to' direction towards the vector is pointing, 'from': direction from the vector is coming
   Returns:
       Theta [deg], R the polar coordinates


.. py:function:: sfc_area_deg(lon1, lon2, lat1, lat2, R=3390000.0)

   Return the surface between two set of latitudes/longitudes
   S= Int[R**2 dlon cos(lat) dlat]     _____lat2
   Args:                               \            lon1,lon2: in [degree]           \____\lat1
       lat1,lat2: in [degree]        lon1    lon2
       R: planetary radius in [m]
   *** NOTE***
   Lon and Lat define the corners of the area, not the grid cells' centers



.. py:function:: area_meridional_cells_deg(lat_c, dlon, dlat, normalize=False, R=3390000.0)

   Return area of invidual cells for a meridional band of thickness dlon
   S= Int[R**2 dlon cos(lat) dlat]
   with  sin(a)-sin(b)=2 cos((a+b)/2)sin((a+b)/2)
   >>> S= 2 R**2 dlon 2 cos(lat)sin(dlat/2)         _________lat+dlat/2
   Args:                                            \    lat \             ^
       lat_c: latitude of cell center in [degree]    \lon +   \            | dlat
       dlon : cell angular width  in [degree]         \________\lat-dlat/2 v
       dlat : cell angular height in [degree]   lon-dlon/2      lon+dlon/2
       R: planetary radius in [m]                       <------>
       normalize: if True, sum of output elements is 1.   dlon
   Returns:
       S: areas of the cells, same size as lat_c in [m2] or normalized by the total area


.. py:function:: area_weights_deg(var_shape, lat_c, axis=-2)

   Return weights for averaging of the variable var.
   Args:
       var_shape: variable's shape, e.g. [133,36,48,46] typically obtained with 'var.shape'
       Expected dimensions are:                      (lat) [axis not needed]
                                                (lat, lon) [axis=-2 or axis=0]
                                          (time, lat, lon) [axis=-2 or axis=1]
                                     (time, lev, lat, lon) [axis=-2 or axis=2]
                          (time, time_of_day_24, lat, lon) [axis=-2 or axis=2]
                     (time, time_of_day_24, lev, lat, lon) [axis=-2 or axis=3]

       lat_c: latitude of cell centers in [degree]
       axis: Position of the latitude axis for 2D and higher-dimensional arrays. The default is the SECOND TO LAST dimension, e.g: axis=-2
          >>> Because dlat is computed as lat_c[1]-lat_c[0] lat_c may be truncated on either end (e.g. lat= [-20 ...,0... +50]) but must be contineous.
   Returns:
       W: weights for var, ready for standard averaging as np.mean(var*W) [condensed form] or np.average(var,weights=W) [expended form]

   ***NOTE***
   Given a variable var:
       var= [v1,v2,...vn]
   Regular average is:
       AVG = (v1+v2+... vn)/N
   Weighted average is:
       AVG_W= (v1*w1+v2*w2+... vn*wn)/(w1+w2+...wn)

   This function returns:
       W= [w1,w2,... ,wn]*N/(w1+w2+...wn)

   >>> Therfore taking a regular average of (var*W) with np.mean(var*W) or np.average(var,weights=W) returns the weighted-average of var
   Use np.average(var,weights=W,axis=X) to average over a specific axis


.. py:function:: areo_avg(VAR, areo, Ls_target, Ls_angle, symmetric=True)

   Return a value average over a central solar longitude

   Args:
       VAR: a ND variable variable with the 1st dimensions being the time, e.g (time,lev,lat,lon)
       areo: 1D array of solar longitude of the input variable in degree (0->720)
       Ls_target: central solar longitude of interest.
       Ls_angle:  requested window angle centered around    Ls_target
       symmetric: a boolean (default =True) If True, and if the requested window is out of range, Ls_angle is reduced
                                            If False, the time average is done on the data available
   Returns:
       The variable VAR averaged over solar longitudes  Ls_target-Ls_angle/2 to Ls_target+Ls_angle/2
        E.g in our example the size would (lev,lat,lon)

   Expl:  Ls_target= 90.
          Ls_angle = 10.

          ---> Nominally, the time average is done over solar longitudes      85 <Ls_target < 95 (10 degree)

          ---> If  symmetric =True and the input data ranges from Ls 88 to 100     88 <Ls_target < 92 (4  degree, symmetric)
               If  symmetric =False and the input data ranges from Ls 88 to 100    88 <Ls_target < 95 (7  degree, assymetric)
   *NOTE*

   [Alex] The routine will bin data from muliples Mars years if provided



.. py:function:: mass_stream(v_avg, lat, level, type='pstd', psfc=700, H=8000.0, factor=1e-08)

   Compute the mass stream function.
                           P
                           ⌠
   Phi=(2 pi a) cos(lat)/g ⎮vz_tavg dp
                           ⌡
                           p_top
   Args:

       v_avg:  zonal winds  [m/s] with 'level' dimensions FIRST and 'lat' dimension SECOND e.g (pstd,lat), (pstd,lat,lon) or (pstd,lat,lon,time)
                >> This routine is set-up so the time and zonal averages may be done either ahead or after the MSF calculation.
       lat  :1D array of latitudes in [degree]
       level: interpolated layers in [Pa] or [m]
       type : interpolation type, i.e. 'pstd', 'zstd' or 'zagl'
       psfc : reference surface pressure in [Pa]
       H    : reference scale height in [m] when pressure are used.
       factor: normalize the mass stream function by a factor, use factor =1. to obtain [kg/s]
   Returns:
       MSF: The meridional mass stream function in factor*[kg/s]
   ***NOTE***
   [Alex. K] : The expressions for the MSF I have seen uses log(pressure) Z coordinate, which I assume integrates better numerically.

   With p=p_sfc exp(-Z/H)  i.e. Z= H log(p_sfc/p) ==> dp= -p_sfc/H exp(-Z/H) dZ, we have:

                                     Z_top
                                    ⌠
   Phi=+(2 pi a) cos(lat)psfc/(g H) ⎮v_rmv exp(-Z/H) dZ  With p=p_sfc exp(-Z/H)
                                    ⌡
                                    Z
                                                            n
                                                           ⌠
   The integral is calculated using trapezoidal rule, e.g. ⌡ f(z)dz  = (Zn-Zn-1){f(Zn)+f(Zn-1)}/2
                                                           n-1


.. py:function:: vw_from_MSF(msf, lat, lev, ztype='pstd', norm=True, psfc=700, H=8000.0)

   Return the [v] and [w] component of the circulation from the mass stream function.

   Args:
       msf  : the mass stream function with 'level' SECOND to LAST and the 'latitude' dimension LAST, e.g. (lev,lat), (time,lev,lat), (time,lon,lev,lat)...
       lat  : 1D latitude array in [degree]
       lev  : 1D level array  in [Pa] or [m]  e.g. pstd, zagl, zstd
       ztype: Use 'pstd' for pressure so vertical differentation is done in log space.
       norm : if  True, normalize  the lat and lev before differentiation avoid having to rescale manually  the vectors in quiver plots
       psfc : surface  pressure for pseudo-height when ztype ='pstd'
       H    : scale height for pseudo-height when ztype ='pstd'
   Return:
       V,W the meditional and altitude component of the mass stream function, to be plotted as quiver or streamlines.

   ***NOTE***
   The components are:
       [v]=  g/(2 pi cos(lat)) dphi/dz
       [w]= -g/(2 pi cos(lat)) dphi/dlat


.. py:function:: alt_KM(press, scale_height_KM=8.0, reference_press=610.0)

   Gives the approximate altitude in km for a given pressure
   Args:
       press: the pressure in [Pa]
       scale_height_KM: a scale height in [km], (default is 8 km, an isothermal at 155K)
       reference_press: reference surface pressure in [Pa], (default is 610 Pa)
   Returns:
       z_KM: the equivalent altitude for that pressure level in [km]

   ***NOTE***
   Scale height is H=rT/g


.. py:function:: press_pa(alt_KM, scale_height_KM=8.0, reference_press=610.0)

   Gives the approximate altitude in km for a given pressure
   Args:
       alt_KM: the altitude in  [km]
       scale_height_KM: a scale height in [km], (default is 8 km, an isothermal at 155K)
       reference_press: reference surface pressure in [Pa], (default is 610 Pa)
   Returns:
        press_pa: the equivalent pressure at that altitude in [Pa]
   ***NOTE***
   Scale height is H=rT/g


.. py:function:: shiftgrid_360_to_180(lon, data)

   This function shift N dimensional data a 0->360 to a -180/+180 grid.
   Args:
       lon: 1D array of longitude 0->360
       data: ND array with last dimension being the longitude (transpose first if necessary)
   Returns:
       data: shifted data
   Note: Use np.ma.hstack instead of np.hstack to keep the masked array properties


.. py:function:: shiftgrid_180_to_360(lon, data)

   This function shift N dimensional data a -180/+180 grid to a 0->360
   Args:
       lon: 1D array of longitude -180/+180
       data: ND array with last dimension being the longitude (transpose first if necessary)
   Returns:
       data: shifted data


.. py:function:: second_hhmmss(seconds, lon_180=0.0)

   Given the time in seconds return Local true Solar Time at a certain longitude
   Args:
       seconds: a float, the time in seconds
       lon_180: a float, the longitude in -/+180 coordinate
   Returns:
       hours: float, the local time or  (hours,minutes, seconds)



.. py:function:: sol_hhmmss(time_sol, lon_180=0.0)

   Given the time in days, return the Local true Solar Time at a certain longitude
   Args:
       time_sol: a float, the time, eg. sols 2350.24
       lon_180: a float, the longitude in a -/+180 coordinate
   Returns:
       hours: float, the local time or  (hours,minutes, seconds)


.. py:function:: UT_LTtxt(UT_sol, lon_180=0.0, roundmin=None)

   Returns the time in HH:MM:SS format at a certain longitude.
   Args:
       time_sol: a float, the time, eg. sols 2350.24
       lon_180: a float, the center longitude in  -/+180 coordinate. Increment by 1hr every 15 deg
       roundmin: round to the nearest X minute  Typical values are roundmin=1,15,60
   ***Note***
   If roundmin is requested, seconds are not shown


.. py:function:: dvar_dh(arr, h=None)

   Differentiate an array A(dim1,dim2,dim3...) with respect to h. The differentiated dimension must be the first dimension.
   > If h is 1D: h and dim1 must have the same length
   > If h is 2D, 3D or 4D, arr and h must have the same shape
   Args:
       arr:   an array of dimension n
       h:     the dimension, eg Z, P, lat, lon

   Returns:
       d_arr: the array differentiated with respect to h, e.g d(array)/dh

   *Example*
    #Compute dT/dz where T[time,LEV,lat,lon] is the temperature and Zkm is the array of  level heights in Km:
    #First we transpose t so the vertical dimension comes first as T[LEV,time,lat,lon] and then we transpose back to get dTdz[time,LEV,lat,lon].
    dTdz=dvar_dh(t.transpose([1,0,2,3]),Zkm).transpose([1,0,2,3])



.. py:function:: zonal_detrend(VAR)

   Substract zonnally averaged mean value from a field
   Args:
       VAR: ND-array with detrending dimension last (e.g time,lev,lat,lon)
   Returns:
       OUT: detrented field (same size as input)

   ***NOTE***
   RuntimeWarnings are expected if the slice contains only NaN, which is the case below the surface
   and above the model's top in the interpolated files. We will disable those warnings temporarily


.. py:function:: get_trend_2D(VAR, LON, LAT, type_trend='wmean')

   Extract spatial trend from data. The output can be directly substracted from the original field.
   Args:
       VAR:  Variable for decomposition, latitude is SECOND to LAST and longitude is LAST  e.g. (time,lat,lon) or (time,lev,lat,lon)
       LON,LAT: 2D arrays of coordinates
       type_trend:  'mean' > use a constant average over all latitude/longitude
                    'wmean'> use a area-weighted average over all latitude/longitude
                    'zonal'> detrend over the zonal axis only
                    '2D'   > use a 2D planar regression (not area-weighted)
   Returns:
       TREND      : trend, same size as VAR e.g. (time,lev,lat,lon)


.. py:function:: regression_2D(X, Y, VAR, order=1)

   Linear and quadratic regression on the plane.
   Args:
       X: 2D array of first coordinate
       Y: 2D array of decond coordinate
       VAR: 2D array, same size as X
       order : 1 (linear) or 2 (quadratic)


   ***NOTE***
   With order =1, the equation is: aX + bY + C = Z
   With order =2, the equation is:  a X**2 + 2b X*Y +c Y**2 +2dX +2eY+f = Z

   For the linear case:
   > ax + by + c = z is re-writtent as A X =b with:
       |x0   y0   1|        |a      |z0
   A = |x1   y1   1|    X = |b   b= |
       |      ...  |        |c      |...
       |xn   yn   1|                |zn

           [n,3]           [3]       [n]

   The least square regression provides the solution that that minimizes  ||b – A x||**2


.. py:function:: daily_to_average(varIN, dt_in, nday=5, trim=True)

   Bin a variable from an atmos_daily file to the atmos_average format.
   Args:
       varIN: ND-array with time dimension first (e.g ts(time,lat,lon))
       dt_in: Delta of time betwen timesteps in sols, e.g. dt_in=time[1]-time[0]
       nday : bining period in sols, default is 5 sols
       trim : discard any leftover data at the end of file before binning.

   Returns:
       varOUT: the variable bin over nday

   ***NOTE***

   1) If varIN(time,lat,lon) from atmos_daily = (160,48,96) and has 4 timestep per day (every 6 hours), the resulting variable  for nday=5 is
      varOUT(160/(4x5),48,96)=varOUT(8,48,96)

   2) If daily file is 668 sols, that is =133x5 +3 leftover sols.
      >If trim=True,  the time is 133 and last 3 sols the are discarded
      >If trim=False, the time is 134 and last bin contains only 3 sols of data


.. py:function:: daily_to_diurn(varIN, time_in)

   Bin a variable from an atmos_daily file into the atmos_diurn format.
   Args:
       varIN: ND-array with time dimension first (e.g ts(time,lat,lon))
       time_in: Time array in sols. Only the first N elements (e.g. time[0:N]) are actually needed (if saving memory is important).
   Returns:
       varOUT: the variable bined in the atmos_diurn format, e.g. ts(time,time_of_day,lat,lon)
       tod   : time of day in [hours]

   ***NOTE***
   1) If varIN(time,lat,lon) from atmos_daily = (40,48,96) and has 4 timestep per day (every 6 hours):
   > The resulting variable is varOUT(10,4,48,96)=(time,time_of_day,lat,lon)
   > tod=[0.,6.,12.,18.] (for example)

   2) Since the time dimension remains first, the output variables may be passed to the daily_to_average() function for further binning.


.. py:function:: gauss_profile(x, alpha, x0=0.0)

   Return Gaussian line shape at x This can be used to generate a bell-shaped mountain


.. py:function:: compute_uneven_sigma(num_levels, N_scale_heights, surf_res, exponent, zero_top)

   Construct an initial array of sigma based on the number of levels, an exponent
   Args:
       num_levels: the number of levels
       N_scale_heights: the number of scale heights to the top of the model (e.g scale_heights =12.5 ~102km assuming 8km scale height)
       surf_res: the resolution at the surface
       exponent: an exponent to increase th thickness of the levels
       zero_top: if True, force the top pressure boundary (in N=0) to 0 Pa
   Returns:
       b: an array of sigma layers



.. py:function:: transition(pfull, p_sigma=0.1, p_press=0.05)

   Return the transition factor to construct the ak and bk
   Args:
       pfull: the pressure in Pa
       p_sigma: the pressure level where the vertical grid starts transitioning from sigma to pressure
       p_press: the pressure level above those  the vertical grid is pure (constant) pressure
   Returns:
       t: the transition factor =1 for pure sigma, 0 for pure pressure and 0<t<1 for the transition

   NOTE:
   In the FV code full pressure are computed from:
                      del(phalf)
        pfull = -----------------------------
                log(phalf(k+1/2)/phalf(k-1/2))


.. py:function:: swinbank(plev, psfc, ptrans=1.0)

   Compute ak and bk values with a transition based on Swinbank
   Args:
       plev: the pressure levels in Pa
       psfc: the surface pressure in Pa
       ptrans:the transition pressure in Pa
   Returns:
        aknew, bknew,ks: the coefficients for the new layers


.. py:function:: polar_warming(T, lat, outside_range=np.NaN)

   Return the polar warming, following  [McDunn et al. 2013]: Characterization of middle-atmosphere polar warming at Mars, JGR
   A. Kling
   Args:
       T:   temperature array, 1D, 2D or ND, with the latitude dimension FIRST (transpose as needed)
       lat: latitude array
       outside_range: values to set the polar warming to outside the range. Default is Nan but 'zero' may be desirable.
   Returns:
       DT_PW:   The polar warming in [K]


   *NOTE*  polar_warming() concatenates the results from both hemispheres obtained from the nested function PW_half_hemisphere()


.. py:function:: tshift(array, lon, timeo, timex=None)

   Conversion to uniform local time.
   Args:
       array: variable to be shifted. Assume longitude is the first dimension and time_of_day is the last dimension
       lon: longitude
       timeo : time_of_day index from input file
       timex (optional) : local time (hr) to shift to, e.g. '3. 15.'
   Returns:
       tshift: array shifted to uniform local time.

   ***Note***
   If timex is not specified, the file is interpolated on the same time_of_day as the input


.. py:function:: lin_interp(X_in, X_ref, Y_ref)

   Simple linear interpolation with no dependance on scipy
   Args:
       X_in (float or array): input values
       X_ref (array): x values
       Y_ref (array): y values
   Returns:
       Y_out: y value linearly interpolated at X_in


.. py:function:: add_cyclic(data, lon)

   Add an additional cyclic (overlapping) point to a 2D array, useful for azimuth and orthographic projections
   Args:
       data: 2D array of size (nlat,nlon)
       lon: 1D array of longitudes
   Returns:
       data_c: 2D array of size (nlat,nlon+1), with last column identical to the 1st
       lon_c: 1D array of longitudes size nlon+1 where the last element is lon[-1]+dlon



.. py:function:: spherical_div(U, V, lon_deg, lat_deg, R=3400 * 1000.0, spacing='varying')

   Compute the divergence of the wind fields using finite difference.
   div = du/dx + dv/dy = 1/(r cos lat)[d(u)/dlon +d(v cos lat)/dlat]
   Args:
       U,V    : wind field with latitude second to last and longitude as last dimensions  e.g. (lat,lon) or (time,lev,lat,lon)...
       lon_deg: 1D array of longitude in [degree] or 2D (lat,lon) if irregularly-spaced
       lat_deg: 1D array of latitude  in [degree] or 2D (lat,lon) if irregularly-spaced
       R      : planetary radius in [m]
       spacing : When lon, lat are  1D arrays, using spacing ='varying' differentiate lat and lon (default)
                 If spacing='regular', only uses uses dx=lon[1]-lon[0], dy=lat[1]-lat[0] and the numpy.gradient() method
   Return:
       div: the horizonal divergence of the wind field   in [m-1]



.. py:function:: spherical_curl(U, V, lon_deg, lat_deg, R=3400 * 1000.0, spacing='varying')

   Compute the vertical component of the relative vorticy using finite difference.
   curl = dv/dx -du/dy  = 1/(r cos lat)[d(v)/dlon +d(u(cos lat)/dlat]
   Args:
       U,V    : wind fields with latitude second to last and longitude as last dimensions  e.g. (lat,lon) or (time,lev,lat,lon)...
       lon_deg: 1D array of longitude in [degree] or 2D (lat,lon) if irregularly-spaced
       lat_deg: 1D array of latitude  in [degree] or 2D (lat,lon) if irregularly-spaced
       R      : planetary radius in [m]
       spacing : When lon, lat are  1D arrays, using spacing ='varying' differentiate lat and lon (default)
                 If spacing='regular', only uses uses dx=lon[1]-lon[0], dy=lat[1]-lat[0] and the numpy.gradient() method
   Return:
       curl: the vorticity of the wind field in [m-1]



.. py:function:: frontogenesis(U, V, theta, lon_deg, lat_deg, R=3400 * 1000.0, spacing='varying')

   Compute the frontogenesis,i.e. local change in potential temperature gradient near a front.
   Following Richter et al. 2010 Toward a Physically Based Gravity Wave Source Parameterization in
    a General Circulation Model, JAS 67 we have Fn= 1/2 D(Del Theta)**2/Dt in [K/m/s]

   Args:
       U,V    : wind fields with latitude second to last and longitude as last dimensions  e.g. (lat,lon) or (time,lev,lat,lon)...
       theta  : potential temperature [K]
       lon_deg: 1D array of longitude in [degree] or 2D (lat,lon) if irregularly-spaced
       lat_deg: 1D array of latitude  in [degree] or 2D (lat,lon) if irregularly-spaced
       R      : planetary radius in [m]
       spacing : When lon, lat are  1D arrays, using spacing ='varying' differentiate lat and lon (default)
                 If spacing='regular', only uses uses dx=lon[1]-lon[0], dy=lat[1]-lat[0] and the numpy.gradient() method
   Return:
       Fn: the frontogenesis field in [m-1]



.. py:function:: MGSzmax_ls_lat(ls, lat)

   Return the max altitude for the dust from "MGS scenario"
   from Montmessin et al. (2004), Origin and role of water ice clouds in the Martian
                                  water cycle as inferred from a general circulation model

   Args:
       ls  : solar longitude in degree
       lat : latitude in degree
   Returns:
       zmax : top altitude for the dusk in [km]


.. py:function:: MGStau_ls_lat(ls, lat)

   Return the max altitude for the dust from "MGS scenario"
   from Montmessin et al. (2004), Origin and role of water ice clouds in the Martian
                                  water cycle as inferred from a general circulation model

   Args:
       ls  : solar longitude in degree
       lat : latitude in degree
   Returns:
       zmax : top altitude for the dusk in [km]


.. py:function:: broadcast(var_1D, shape_out, axis)

   Broadcast a 1D array based on a variable's dimensions
   Args:
       var_1D     (1D array), e.g. lat size (36), or time size (133)
       shape_out (ND list) : braodcasting shape e.g temp.shape= [133,(lev),36,(lon)]
   Return:
       var_b (ND array): broadcasted variables, e.g. size [1,36,1,1] for lat or [133,1,1,1] for time


.. py:function:: ref_atmosphere_Mars_PTD(Zi)

   Analytical atmospheric model for Martian pressure, temperature and density,  [Alex Kling, June 2021]
   Args:
       Zi (float or 1D array): input altitude in m (must be >= 0)
   Return:
       P,T,D (floats ot 1D arrays): tuple of corresponding pressure [Pa], temperature [K] and density  [kg/m3]

   ***NOTE***

   This model was obtained by  fitting globally and annually averaged reference temperature profiles derived from the Legacy GCM, MCS observations,
   and Mars Climate  Database.

   The temperature fit was constructed using quadratic temperature T(z)= T0+ gam(z-z0)+a*(z-z0)**2
   over 4 segments (0>57 km, 57>110km, 110>120 km and 120>300km)

   From the ground to 120km, he pressure is obtained be integrating analytically  the hydrostatic equation:
    dp/dz=-g. p/(rT) with T(z)= T0+ gam(z-z0)+a*(z-z0)**2 .

   Above ~120km P=P0 exp(-(z-z0)g/rT) is not a good approximation as the fluid is in molecula regime. For those altitude, we  provide
   fit in the form P=P0 exp(-az-bz**2),based on diurnal average of the MCD database at lat 0, Ls 150.


.. py:function:: press_to_alt_atmosphere_Mars(Pi)

   Return the altitude in m as a function of pressure from the analytical calculations derived above.
   Args:
       Pi (float or 1D array): input pressure in Pa (must be <=610 Pa)
   Return:
       Z (float ot 1D array): corresponding altitude in m


.. py:function:: azimuth2cart(LAT, LON, lat0, lon0=0)

   Azimuthal equidistant projection, convert from lat/lon to cartesian coordinates
   Args:
       LAT,LON: 1D or 2D array of latitudes, longitudes in degree, size [nlat,nlon]
       lat0,lon0:(floats) coordinates of the pole
   Returns:
       X,Y: cartesian coordinates for the latitudes and longitudes


.. py:function:: ortho2cart(LAT, LON, lat0, lon0=0)

   Orthographic projection, convert from lat/lon to cartesian coordinates
   Args:
       LAT,LON: 1D or 2D array of latitudes, longitudes in degree, size [nlat,nlon]
       lat0,lon0:(floats) coordinates of the pole
   Returns:
       X,Y: cartesian coordinates for the latitudes and longitudes
       MASK: NaN array that is used to hide the back side of the planet


.. py:function:: mollweide2cart(LAT, LON)

   Mollweide projection, convert from lat/lon to cartesian coordinates
   Args:
       LAT,LON: 1D or 2D array of latitudes, longitudes in degree, size [nlat,nlon]
   Returns:
       X,Y: cartesian coordinates for the latitudes and longitudes


.. py:function:: robin2cart(LAT, LON)

   Robinson projection, convert from lat/lon to cartesian coordinates
   Args:
       LAT,LON: floats, 1D or 2D array (nalt,nlon) of latitudes, longitudes in degree
   Returns:
       X,Y: cartesian coordinates for the latitudes and longitudes


.. py:function:: sol2ls(jld, cummulative=False)

   Return the solar longitude Ls as a function of the sol number. Sol 0 is spring equinox.
   Args:
       jld [float or 1D array]: sol number after perihelion
       cummulative [bool]     : if True, result is cummulative Ls 0>360>720 etc..

   Returns:
       Ls: The corresponding solar longitude Ls


.. py:function:: ls2sol(Ls_in)

   Ls to sol converter.
   Args:
       Ls_in (float or 1D array) : Solar longitudes 0-360...720
   Return:
       sol: the corresponding sol number
   ***NOTE***
   This function simply uses a numerical solver on the sol2ls() function.


