:orphan:

:py:mod:`MarsVars`
==================

.. py:module:: MarsVars


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

   Return the 3D pressure field at the layer midpoint.
   *** NOTE***
   The shape_out argument ensures that, when time=1 (one timestep) results are returned as (1,lev,lat,lon), not (lev,lat,lon)


.. py:function:: compute_rho(p_3D, temp)

   Return the density in [kg/m3]


.. py:function:: compute_xzTau(q, temp, lev, const, f_type)

   Return dust or ice extinction in [km-1]
   Adapted from Heavens et al. 2011, observations by MCS (JGR)


.. py:function:: compute_mmr(xTau, temp, lev, const, f_type)

   Return dust or ice mixing ratio [kg/kg]
   Adapted from Heavens et al. 2011, observations by MCS (JGR)


.. py:function:: compute_Vg_sed(xTau, nTau, temp)

   Return sedimentation rate for dust


.. py:function:: compute_w_net(Vg, wvar)

   Return net vertical wind: w - sedimentation rate (Vg_sed)


.. py:function:: compute_theta(p_3D, ps, temp, f_type)

   Return the potential temperature in [K]


.. py:function:: compute_zfull(ps, ak, bk, temp)

   Compute the altitude AGL in [m]


.. py:function:: compute_zhalf(ps, ak, bk, temp)

   Compute the altitude AGL in [m]


.. py:function:: compute_DZ_full_pstd(pstd, temp, ftype='average')

   Return the distance between two layers  mid-point from the standard pressure levels

   Args:
       pstd: 1D  array of standard pressure in [Pa]
       temp : 3D array of temperature
       ftype: 'daily', 'aveage' or 'diurn'
   Returns:
       DZ_full_pstd: 3D array of distancez  between adjacent layers

   *** NOTE***
   In this context p_full = p_std, with the half layers boundaries defined somewhere in between successive layers

   --- Nk --- TOP        ========  p_half
   --- Nk-1 ---
                        --------  p_full = p_std   ^
                                                   | DZ_full_pstd
                        ========  p_half           |
   --- 1 ---            --------  p_full = p_std   v
   --- 0 --- SFC        ========  p_half
                       / / / /


.. py:function:: compute_N(theta, zfull)

   Compute the Brunt Vaisala freqency in [rad/s]


.. py:function:: compute_Tco2(P_3D, temp)

   Compute the frost point of CO2 in [K]
   From [Fannale 1982] Mars: The regolit-atmosphere cap system and climate change. Icarus


.. py:function:: compute_scorer(N, ucomp, zfull)

   Compute the Scorer wavelenght in [m]


.. py:function:: compute_DP_3D(ps, ak, bk, shape_out)

   Compute the thickness of a layer in [Pa]


.. py:function:: compute_DZ_3D(ps, ak, bk, temp, shape_out)

   Compute the thickness of a layer in [Pa]


.. py:function:: compute_Ep(temp)

   Return the wave potential energy: Ep= 1/2 (g/N)**2 (T'/T)**2 in [J/kg]


.. py:function:: compute_Ek(ucomp, vcomp)

   Return the wave kinetic energy: Ek= 1/2 (u'**2+v'**2) in[J/kg]


.. py:function:: compute_MF(UVcomp, w)

   Return the zonal or meridional momentum fluxes u'w' or v'w'


.. py:function:: compute_WMFF(MF, rho, lev, interp_type)

   Return the zonal or meridional wave-mean flow forcing ax= -1/rho d(rho u'w')/dz in [m/s/s]
   ***NOTE***                                            ay= -1/rho d(rho v'w')/dz in [m/s/s]
   For pstd, we have:
       du/dz= (du/dp).(dp/dz) > du/dz=-rho g (du/dp) with dp/dz = -rho g


