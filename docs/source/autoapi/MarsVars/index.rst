:py:mod:`MarsVars`
==================

.. py:module:: MarsVars

.. autoapi-nested-parse::

   The MarsVars executable is for performing variable manipulations in
   existing files. Most often, it is used to derive and add variables to
   existing files, but it also differentiates variables with respect to
   (w.r.t) the Z axis, column-integrates variables, converts aerosol
   opacities from opacity per Pascal to opacity per meter, removes and
   extracts variables from files, and enables scaling variables or editing
   variable names, units, etc.

   The executable requires:
       * ``[input_file]``          the file to be transformed

   and optionally accepts:
       * ``[-add --add]``           derive and add variable to file
       * ``[-zdiff --zdiff]``       differentiate variable w.r.t. Z axis
       * ``[-col --col]``           column-integrate variable
       * ``[-zd --zonal_detrend]``  subtract zonal mean from variable
       * ``[-dp_to_dz --dp_to_dz]`` convert aerosol opacity op/Pa -> op/m
       * ``[-dz_to_dp --dz_to_dp]`` convert aerosol opacity op/m -> op/Pa
       * ``[-rm --remove]``         remove variable from file
       * ``[-extract --extract]``   copy variable to new file
       * ``[-edit --edit]``         edit variable attributes or scale it

   Third-party Requirements:
       * ``numpy``
       * ``netCDF4``
       * ``argparse``
       * ``os``
       * ``subprocess``
       * ``matplotlib``



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   MarsVars.compute_p_3D
   MarsVars.compute_rho
   MarsVars.compute_xzTau
   MarsVars.compute_mmr
   MarsVars.compute_Vg_sed
   MarsVars.compute_w_net
   MarsVars.compute_theta
   MarsVars.compute_w
   MarsVars.compute_zfull
   MarsVars.compute_zhalf
   MarsVars.compute_DZ_full_pstd
   MarsVars.compute_N
   MarsVars.compute_Tco2
   MarsVars.compute_scorer
   MarsVars.compute_DP_3D
   MarsVars.compute_DZ_3D
   MarsVars.compute_Ep
   MarsVars.compute_Ek
   MarsVars.compute_MF
   MarsVars.compute_WMFF



.. py:function:: compute_p_3D(ps, ak, bk, shape_out)

   Compute the 3D pressure at layer midpoints.

   :param ps: Surface pressure (Pa)
   :type ps: array [time, lat, lon]
   :param ak: Vertical coordinate pressure value (Pa)
   :type ak: array [phalf]
   :param bk: Vertical coordinate sigma value (None)
   :type bk: array [phalf]
   :param shape_out: Determines how to handle the dimensions of p_3D.
       If ``len(time) = 1`` (one timestep), ``p_3D`` is returned as
       [1, lev, lat, lon] as opposed to [lev, lat, lon]
   :type shape_out: float

   :raises:

   :return: ``p_3D`` The full 3D pressure array (Pa)
   :rtype: array [time, lev, lat, lon]


.. py:function:: compute_rho(p_3D, temp)

   Compute density.

   :param p_3D: Pressure (Pa)
   :type p_3D: array [time, lev, lat, lon]
   :param temp: Temperature (K)
   :type temp: array [time, lev, lat, lon]

   :raises:

   :return: Density (kg/m^3)
   :rtype: array [time, lev, lat, lon]


.. py:function:: compute_xzTau(q, temp, lev, const, f_type)

   Compute the dust or ice extinction rate.
   Adapted from Heavens et al. (2011) observations from MCS (JGR).

   :param q: Dust or ice mass mixing ratio (ppm)
   :type q: array [time, lev, lat, lon]
   :param temp: Temperature (K)
   :type temp: array [time, lev, lat, lon]
   :param lev: Vertical coordinate (e.g., pstd) (e.g., Pa)
   :type lev: array [lev]
   :param const: Dust or ice constant
   :type const: array
   :param f_type: The FV3 file type: diurn, daily, or average
   :type f_stype: str

   :raises:

   :return: ``xzTau`` Dust or ice extinction rate (km-1)
   :rtype: array [time, lev, lat, lon]


.. py:function:: compute_mmr(xTau, temp, lev, const, f_type)

   Compute the dust or ice mixing ratio.
   Adapted from Heavens et al. (2011) observations from MCS (JGR).

   :param xTau: Dust or ice extinction rate (km-1)
   :type xTau: array [time, lev, lat, lon]
   :param temp: Temperature (K)
   :type temp: array [time, lev, lat, lon]
   :param lev: Vertical coordinate (e.g., pstd) (e.g., Pa)
   :type lev: array [lev]
   :param const: Dust or ice constant
   :type const: array
   :param f_type: The FV3 file type: diurn, daily, or average
   :type f_stype: str

   :raises:

   :return: ``q``, Dust or ice mass mixing ratio (ppm)
   :rtype: array [time, lev, lat, lon]


.. py:function:: compute_Vg_sed(xTau, nTau, temp)

   Calculate the sedimentation rate of the dust.

   :param xTau: Dust or ice MASS mixing ratio (ppm)
   :type xTau: array [time, lev, lat, lon]
   :param nTau: Dust or ice NUMBER mixing ratio (None)
   :type nTau: array [time, lev, lat, lon]
   :param temp: Temperature (K)
   :type temp: array [time, lev, lat, lon]

   :raises:

   :return: ``Vg`` Dust sedimentation rate (m/s)
   :rtype: array [time, lev, lat, lon]


.. py:function:: compute_w_net(Vg, wvar)

   Computes the net vertical wind, which is the vertical wind (w)
   minus the sedimentation rate (``Vg_sed``)::

       w_net = w - Vg_sed

   :param Vg: Dust sedimentation rate (m/s)
   :type Vg: array [time, lev, lat, lon]
   :param wvar: Vertical wind (m/s)
   :type wvar: array [time, lev, lat, lon]

   :raises:

   :return: `w_net` Net vertical wind speed (m/s)
   :rtype: array [time, lev, lat, lon]


.. py:function:: compute_theta(p_3D, ps, temp, f_type)

   Compute the potential temperature. 

   :param p_3D: The full 3D pressure array (Pa)
   :type p_3D: array [time, lev, lat, lon]
   :param ps: Surface pressure (Pa)
   :type ps: array [time, lat, lon]
   :param temp: Temperature (K)
   :type temp: array [time, lev, lat, lon]
   :param f_type: The FV3 file type: diurn, daily, or average
   :type f_type: str

   :raises:

   :return: Potential temperature (K)
   :rtype: array [time, lev, lat, lon]


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
   :type rho: array [time, lev, lat, lon]
   :param omega: Rate of change in pressure at layer midpoint (Pa/s)
   :type omega: array [time, lev, lat, lon]

   :raises:

   :return: vertical wind (m/s)
   :rtype: array [time, lev, lat, lon]


.. py:function:: compute_zfull(ps, ak, bk, temp)

   Calculate the altitude of the layer midpoints above ground level.

   :param ps: Surface pressure (Pa)
   :type ps: array [time, lat, lon]
   :param ak: Vertical coordinate pressure value (Pa)
   :type ak: array [phalf]
   :param bk: Vertical coordinate sigma value (None)
   :type bk: array [phalf]
   :param temp: Temperature (K)
   :type temp: array [time, lev, lat, lon]

   :raises:

   :return: ``zfull`` (m)
   :rtype: array [time, lev, lat, lon]


.. py:function:: compute_zhalf(ps, ak, bk, temp)

   Calculate the altitude of the layer interfaces above ground level.

   :param ps: Surface pressure (Pa)
   :type ps: array [time, lat, lon]
   :param ak: Vertical coordinate pressure value (Pa)
   :type ak: array [phalf]
   :param bk: Vertical coordinate sigma value (None)
   :type bk: array [phalf]
   :param temp: Temperature (K)
   :type temp: array [time, lev, lat, lon]

   :raises:

   :return: ``zhalf`` (m)
   :rtype: array [time, lev, lat, lon]


.. py:function:: compute_DZ_full_pstd(pstd, temp, ftype='average')

   Calculate the thickness of a layer from the midpoint of the
   standard pressure levels (``pstd``).

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
   :type pstd: array [lev]
   :param temp: Temperature (K)
   :type temp: array [time, lev, lat, lon]
   :param f_type: The FV3 file type: diurn, daily, or average
   :type f_stype: str

   :raises:

   :return: DZ_full_pstd, Layer thicknesses (Pa)
   :rtype: array [time, lev, lat, lon]


.. py:function:: compute_N(theta, zfull)

   Calculate the Brunt Vaisala freqency.

   :param theta: Potential temperature (K)
   :type theta: array [time, lev, lat, lon]
   :param zfull: Altitude above ground level at the layer midpoint (m)
   :type zfull: array [time, lev, lat, lon]

   :raises:

   :return: ``N``, Brunt Vaisala freqency [rad/s]
   :rtype: array [time, lev, lat, lon]


.. py:function:: compute_Tco2(P_3D)

   Calculate the frost point of CO2.
   Adapted from Fannale (1982) - Mars: The regolith-atmosphere cap
   system and climate change. Icarus.

   :param P_3D: The full 3D pressure array (Pa)
   :type p_3D: array [time, lev, lat, lon]

   :raises:

   :return: CO2 frost point [K]
   :rtype: array [time, lev, lat, lon]


.. py:function:: compute_scorer(N, ucomp, zfull)

   Calculate the Scorer wavelength.

   :param N: Brunt Vaisala freqency (rad/s)
   :type N: float [time, lev, lat, lon]
   :param ucomp: Zonal wind (m/s)
   :type ucomp: array [time, lev, lat, lon]
   :param zfull: Altitude above ground level at the layer midpoint (m)
   :type zfull: array [time, lev, lat, lon]

   :raises:

   :return: ``scorer_wl`` Scorer horizontal wavelength (m)
   :rtype: array [time, lev, lat, lon]


.. py:function:: compute_DP_3D(ps, ak, bk, shape_out)

   Calculate the thickness of a layer in pressure units.

   :param ps: Surface pressure (Pa)
   :type ps: array [time, lat, lon]
   :param ak: Vertical coordinate pressure value (Pa)
   :type ak: array [phalf]
   :param bk: Vertical coordinate sigma value (None)
   :type bk: array [phalf]
   :param shape_out: Determines how to handle the dimensions of DP_3D.
       If len(time) = 1 (one timestep), DP_3D is returned as
       [1, lev, lat, lon] as opposed to [lev, lat, lon]
   :type shape_out: float

   :raises:

   :return: ``DP`` Layer thickness in pressure units (Pa)
   :rtype: array [time, lev, lat, lon]


.. py:function:: compute_DZ_3D(ps, ak, bk, temp, shape_out)

   Calculate the thickness of a layer in altitude units.

   :param ps: Surface pressure (Pa)
   :type ps: array [time, lat, lon]
   :param ak: Vertical coordinate pressure value (Pa)
   :type ak: array [phalf]
   :param bk: Vertical coordinate sigma value (None)
   :type bk: array [phalf]
   :param shape_out: Determines how to handle the dimensions of DZ_3D.
       If len(time) = 1 (one timestep), DZ_3D is returned as
       [1, lev, lat, lon] as opposed to [lev, lat, lon]
   :type shape_out: float

   :raises:

   :return: ``DZ`` Layer thickness in altitude units (m)
   :rtype: array [time, lev, lat, lon]


.. py:function:: compute_Ep(temp)

   Calculate wave potential energy::

       Ep = 1/2 (g/N)^2 (temp'/temp)^2

   :param temp: Temperature (K)
   :type temp: array [time, lev, lat, lon]

   :raises:

   :return: ``Ep`` Wave potential energy (J/kg)
   :rtype: array [time, lev, lat, lon]


.. py:function:: compute_Ek(ucomp, vcomp)

   Calculate wave kinetic energ::

       Ek = 1/2 (u'**2+v'**2)

   :param ucomp: Zonal wind (m/s)
   :type ucomp: array [time, lev, lat, lon]
   :param vcomp: Meridional wind (m/s)
   :type vcomp: array [time, lev, lat, lon]

   :raises:

   :return: ``Ek`` Wave kinetic energy (J/kg)
   :rtype: array [time, lev, lat, lon]


.. py:function:: compute_MF(UVcomp, w)

   Calculate zonal or meridional momentum fluxes.

   :param UVcomp: Zonal or meridional wind (ucomp or vcomp)(m/s)
   :type UVcomp: array
   :param w: Vertical wind (m/s)
   :type w: array [time, lev, lat, lon]

   :raises:

   :return: ``u'w'`` or ``v'w'``, Zonal/meridional momentum flux (J/kg)
   :rtype: array [time, lev, lat, lon]


.. py:function:: compute_WMFF(MF, rho, lev, interp_type)

   Calculate the zonal or meridional wave-mean flow forcing::

       ax = -1/rho d(rho u'w')/dz
       ay = -1/rho d(rho v'w')/dz

   If interp_type == ``pstd``, then::

       [du/dz = (du/dp).(dp/dz)] > [du/dz = -rho*g * (du/dp)]

   where::

       dp/dz = -rho*g
       [du/dz = (du/dp).(-rho*g)] > [du/dz = -rho*g * (du/dp)]

   :param MF: Zonal/meridional momentum flux (J/kg)
   :type MF: array [time, lev, lat, lon]
   :param rho: Atmospheric density (kg/m^3)
   :type rho: array [time, lev, lat, lon]
   :param lev: Array for the vertical grid (zagl, zstd, pstd, or pfull)
   :type lev: array [lev]
   :param interp_type: The vertical grid type (``zagl``, ``zstd``,
       ``pstd``, or ``pfull``)
   :type interp_type: str

   :raises:

   :return: The zonal or meridional wave-mean flow forcing (m/s2)
   :rtype: array [time, lev, lat, lon]


