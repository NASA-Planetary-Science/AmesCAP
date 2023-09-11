:py:mod:`MarsVars`
==================

.. py:module:: MarsVars

.. autoapi-nested-parse::

   The MarsVars executable is for ...

   The executable requires x arguments:
       * [-x --x]      define

   Third-party Requirements:
       * numpy
       * argparse
       * requests

   List of Functions:
       * x



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
   The shape_out argument ensures that when time = 1 (one timestep), results are returned
   as (1, lev, lat, lon) not (lev, lat, lon)


.. py:function:: compute_rho(p_3D, temp)

   Returns density in [kg/m3].


.. py:function:: compute_xzTau(q, temp, lev, const, f_type)

   Returns dust or ice extinction in [km-1].
   Adapted from Heavens et al. 2011, observations by MCS (JGR).


.. py:function:: compute_mmr(xTau, temp, lev, const, f_type)

   Return dust or ice mixing ratio [kg/kg]
   Adapted from Heavens et al. 2011. observations by MCS (JGR)


.. py:function:: compute_Vg_sed(xTau, nTau, temp)

   Returns the dust sedimentation rate.


.. py:function:: compute_w_net(Vg, wvar)

   Returns the net vertical wind (subtracts the sedimentation rate (Vg_sed)
   from the vertical wind (w))
   w_net = w - Vg_sed


.. py:function:: compute_theta(p_3D, ps, temp, f_type)

   Returns the potential temperature in [K].


.. py:function:: compute_zfull(ps, ak, bk, temp)

   Returns the altitude of the layer midpoints AGL in [m].


.. py:function:: compute_zhalf(ps, ak, bk, temp)

   Returns the altitude of the layer interfaces AGL in [m]


.. py:function:: compute_DZ_full_pstd(pstd, temp, ftype='average')

   Returns the thickness of a layer (distance between two layers) from the
   midpoint of the standard pressure levels ('pstd').

   Args:
       pstd:   1D array of standard pressure in [Pa]
       temp:   3D array of temperature
       ftype: 'daily', 'aveage', or 'diurn'
   Returns:
       DZ_full_pstd: 3D array of thicknesses

   *** NOTE***
   In this context, 'pfull' = 'pstd' with the layer interfaces defined somewhere
   in between successive layers.

   --- Nk --- TOP       ========  phalf
   --- Nk-1 ---
                        --------  pfull = pstd    ^
                                                  | DZ_full_pstd
                        ========  phalf           |
   --- 1 ---            --------  pfull = pstd    v
   --- 0 --- SFC        ========  phalf
                       / / / /


.. py:function:: compute_N(theta, zfull)

   Returns the Brunt Vaisala freqency in [rad/s].


.. py:function:: compute_Tco2(P_3D, temp)

   Returns the frost point of CO2 in [K].
   Adapted from Fannale, 1982. Mars: The regolith-atmosphere cap system and climate change. Icarus.


.. py:function:: compute_scorer(N, ucomp, zfull)

   Returns the Scorer wavelength in [m].


.. py:function:: compute_DP_3D(ps, ak, bk, shape_out)

   Returns the thickness of a layer in [Pa].


.. py:function:: compute_DZ_3D(ps, ak, bk, temp, shape_out)

   Returns the layer thickness in [Pa].


.. py:function:: compute_Ep(temp)

   Returns the wave potential energy (Ep) in [J/kg].
   Ep = 1/2 (g/N)**2 (T'/T)**2


.. py:function:: compute_Ek(ucomp, vcomp)

   Returns the wave kinetic energy (Ek) in [J/kg].
   Ek= 1/2 (u'**2+v'**2)


.. py:function:: compute_MF(UVcomp, w)

   Returns the zonal or meridional momentum fluxes (u'w' or v'w').


.. py:function:: compute_WMFF(MF, rho, lev, interp_type)

   Returns the zonal or meridional wave-mean flow forcing
   ax= -1/rho d(rho u'w')/dz in [m/s/s]
   ay= -1/rho d(rho v'w')/dz in [m/s/s]

   For 'pstd':
       [du/dz = (du/dp).(dp/dz)] > [du/dz = -rho g (du/dp)] with dp/dz = -rho g


