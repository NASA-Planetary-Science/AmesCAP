Derived Variables
=================

``MarsVars -add`` derives the variables below from fields already in a
file. Run ``MarsVars -h`` for the list in your installed version. Each
entry gives the definition, units, the fields it needs, the vertical
grids it supports, and the reference for the method. The
``compute_*`` function named in each entry is documented in the
:doc:`MarsVars API <autoapi/bin/MarsVars/index>`.

Grids are ``pfull`` (native model levels), ``pstd`` (standard pressure),
``zstd`` (standard altitude), and ``zagl`` (altitude above ground level);
files on the last three are made with ``MarsInterp``. Planetary constants
(gravity :math:`g`, gas constant :math:`R`, :math:`c_p`, radius
:math:`a`) are read from the file when available (e.g., planetWRF
output) and otherwise default to Mars values.

Notation: a prime denotes the deviation from the zonal mean,
:math:`X' = X - \overline{X}`, where the overbar is the zonal mean.

Thermodynamics and vertical structure
-------------------------------------

``pfull3D`` -- Mid-layer pressure [Pa]
   Pressure at the layer midpoints. The interface pressures come from
   the hybrid vertical coordinate, :math:`p_{k+1/2} = a_k + b_k p_s`,
   and the midpoint pressure is the FV3 log-mean of the two bounding
   interfaces.
   Needs ``ps``, ``temp`` and the ``ak``, ``bk`` coefficients. Grid:
   ``pfull``.

``DP`` -- Layer thickness in pressure [Pa]
   :math:`\Delta p_k = p_{k+1/2} - p_{k-1/2}` between the layer
   interfaces. Needs ``ps``, ``temp``. Grid: ``pfull``.
   Function: ``compute_DP_3D``.

``zfull`` -- Mid-layer altitude above ground [m]
   Integrated upward from the surface with the hypsometric equation,
   :math:`\Delta z = (R T/g)\ln(p_{k-1/2}/p_{k+1/2})`.
   Needs ``ps``, ``temp``. Grid: ``pfull``. Function: ``compute_zfull``.

``DZ`` -- Layer thickness in altitude [m]
   Thickness of each layer from the hypsometric equation above.
   Needs ``ps``, ``temp``. Grid: ``pfull``. Function: ``compute_DZ_3D``.

``rho`` -- Density [kg/m\ :sup:`3`]
   Ideal gas law, :math:`\rho = p/(R T)`.
   Needs ``ps``, ``temp``. Grid: ``pfull``. Function: ``compute_rho``.

``theta`` -- Potential temperature [K]
   :math:`\theta = T (p_s/p)^{\kappa}` with :math:`\kappa = R/c_p`.
   The reference pressure is the local surface pressure :math:`p_s`.
   Needs ``ps``, ``temp``. Grid: ``pfull``. Function: ``compute_theta``.

``N`` -- Brunt-Väisälä frequency [rad/s]
   :math:`N = \sqrt{(g/\theta)\, \partial\theta/\partial z}`.
   Needs ``ps``, ``temp``. Grid: ``pfull``. Function: ``compute_N``.
   Reference: Holton and Hakim (2013).

``Ri`` -- Richardson number [dimensionless]
   :math:`Ri = N^2 / \left[(\partial u/\partial z)^2 + (\partial v/\partial z)^2\right]`.
   Needs ``ps``, ``temp``, ``ucomp``, ``vcomp``. Grid: ``pfull``.
   Reference: Holton and Hakim (2013).

``Tco2`` -- CO\ :sub:`2` condensation temperature [K]
   Below the CO\ :sub:`2` triple-point pressure (518 kPa),
   :math:`T = -3167.8 / [\ln(0.01\,p) - 23.23]`, with :math:`p` in Pa.
   Needs ``ps``, ``temp``. Grids: ``pfull``, ``pstd``.
   Function: ``compute_Tco2``. Reference: Fanale et al. (1982).

Winds and circulation
---------------------

``w`` -- Vertical wind [m/s]
   From the pressure velocity under hydrostatic balance,
   :math:`w = -\omega/(\rho g)`.
   Needs ``ps``, ``temp``, ``omega``. Grid: ``pfull``.
   Function: ``compute_w``.

``wspeed`` -- Horizontal wind speed [m/s]
   :math:`\sqrt{u^2 + v^2}`. Needs ``ucomp``, ``vcomp``. All grids.

``wdir`` -- Wind direction [degrees]
   Direction the wind blows from, clockwise from north (meteorological
   convention). Needs ``ucomp``, ``vcomp``. All grids.

``div`` -- Horizontal divergence [s\ :sup:`-1`]
   :math:`\nabla\cdot\mathbf{v} = \frac{1}{a\cos\phi}\left[\frac{\partial u}{\partial\lambda} + \frac{\partial (v\cos\phi)}{\partial\phi}\right]`,
   with longitude :math:`\lambda` and latitude :math:`\phi`. Staggered
   C-grid winds are used when present in the file. Needs ``ucomp``,
   ``vcomp``. All grids. Reference: Holton and Hakim (2013).

``curl`` -- Relative vorticity [s\ :sup:`-1`]
   :math:`\zeta = \frac{1}{a\cos\phi}\left[\frac{\partial v}{\partial\lambda} - \frac{\partial (u\cos\phi)}{\partial\phi}\right]`.
   Staggered C-grid winds are used when present. Needs ``ucomp``,
   ``vcomp``. All grids. Reference: Holton and Hakim (2013).

``msf`` -- Mass streamfunction [10\ :sup:`8` kg/s]
   :math:`\Psi(p) = \frac{2\pi a\cos\phi}{g}\int_0^{p} v\,dp'`,
   integrated with the trapezoidal rule in log-pressure height
   :math:`Z = H\ln(p_s/p)` (:math:`H` = 8 km, :math:`p_s` = 700 Pa).
   On ``pstd`` the integral runs downward from the model top. On
   ``zstd`` and ``zagl`` it runs upward from the lowest level, with
   :math:`e^{-Z/H}` as the density weighting. Usually applied to
   time- and zonal-mean meridional wind. Needs ``vcomp``, ``temp``.
   Grids: ``pstd``, ``zstd``, ``zagl``. Function: ``mass_stream``.
   Reference: Holton and Hakim (2013).

``fn`` -- Frontogenesis [K m\ :sup:`-1` s\ :sup:`-1`]
   :math:`F_n = \tfrac{1}{2}\, D(|\nabla\theta|^2)/Dt`.
   Needs ``ucomp``, ``vcomp``, ``theta``. Grids: ``pstd``, ``zstd``,
   ``zagl``. Function: ``frontogenesis``. Reference: Richter et al.
   (2010).

Waves
-----

``tp_t`` -- Normalized temperature perturbation [dimensionless]
   :math:`T'/T`. Needs ``temp``. Grids: ``pstd``, ``zstd``, ``zagl``.

``ek`` -- Wave kinetic energy [J/kg]
   :math:`E_k = \tfrac{1}{2}\left(u'^2 + v'^2\right)`.
   Needs ``ucomp``, ``vcomp``. Grids: ``pstd``, ``zstd``, ``zagl``.
   Function: ``compute_Ek``. Reference: Andrews et al. (1987).

``ep`` -- Wave potential energy [J/kg]
   :math:`E_p = \tfrac{1}{2}\,(g/N)^2\,(T'/T)^2`, with a constant
   :math:`N` = 0.01 rad/s. Needs ``temp``. Grids: ``pstd``, ``zstd``,
   ``zagl``. Function: ``compute_Ep``. Reference: Andrews et al. (1987).

``mx``, ``my`` -- Vertical flux of zonal and meridional momentum [m\ :sup:`2`/s\ :sup:`2`]
   :math:`u'w'` and :math:`v'w'`. Needs ``ucomp`` or ``vcomp`` and
   ``w``. Grids: ``pstd``, ``zstd``, ``zagl``. Function: ``compute_MF``.
   Reference: Andrews et al. (1987).

``ax``, ``ay`` -- Zonal and meridional wave-mean flow forcing [m/s\ :sup:`2`]
   :math:`a_x = -\frac{1}{\rho}\frac{\partial(\rho\,u'w')}{\partial z}`
   and likewise for :math:`a_y` with :math:`v'w'`. On ``pstd``, the
   vertical derivative uses :math:`\partial/\partial z = -\rho g\,\partial/\partial p`.
   Needs ``ucomp`` or ``vcomp``, ``w``, ``rho``. Grids: ``pstd``,
   ``zstd``, ``zagl``. Function: ``compute_WMFF``.
   Reference: Andrews et al. (1987).

``scorer_wl`` -- Scorer horizontal wavelength [m]
   :math:`\lambda = 2\pi/l`, with the Scorer parameter
   :math:`l^2 = N^2/u^2 - (\partial^2 u/\partial z^2)/u`.
   Needs ``ps``, ``temp``, ``ucomp``. Grid: ``pfull``.
   Function: ``compute_scorer``. Reference: Scorer (1949).

Aerosols
--------

``dzTau``, ``izTau`` -- Dust and water ice extinction [km\ :sup:`-1`]
   :math:`\beta = \rho\, q / C` with
   :math:`C = \tfrac{4}{3}\,\rho_p\, r_\mathrm{eff} / Q_\mathrm{ext}`,
   where :math:`q` is the mass mixing ratio and :math:`\rho = p/(R T)`.
   Dust uses :math:`\rho_p` = 2500 kg/m\ :sup:`3`,
   :math:`r_\mathrm{eff}` = 1.06 µm, :math:`Q_\mathrm{ext}` = 0.35;
   water ice uses :math:`\rho_p` = 900 kg/m\ :sup:`3`,
   :math:`r_\mathrm{eff}` = 1.41 µm, :math:`Q_\mathrm{ext}` = 0.773.
   Needs ``dst_mass_mom`` or ``ice_mass_mom`` and ``temp``. Grid:
   ``pfull``. Function: ``compute_xzTau``. References: Heavens et al.
   (2010, 2011), Kleinböhl et al. (2011).

``dst_mass_mom``, ``ice_mass_mom`` -- Dust and water ice mass mixing ratio [kg/kg]
   The inverse of the extinction relation above,
   :math:`q = C\,\beta/\rho`. Needs ``dzTau`` or ``izTau`` and
   ``temp``. Grid: ``pfull``. Function: ``compute_mmr``.

``Vg_sed`` -- Dust sedimentation velocity [m/s]
   Stokes settling with the Cunningham slip correction,
   :math:`V_g = \frac{2\rho_p g r^2}{9\eta}\left[1 + \mathrm{Kn}\left(A + B e^{-C/\mathrm{Kn}}\right)\right]`,
   with :math:`A` = 1.246, :math:`B` = 0.42, :math:`C` = 0.87, Knudsen
   number :math:`\mathrm{Kn} = \ell/r` and gas mean free path
   :math:`\ell = 2\eta/(\rho\bar{c})`, where
   :math:`\bar{c} = \sqrt{3 k_B T/m_{CO_2}}` and :math:`\rho` is a
   constant reference density (610 Pa, 150 K). The particle radius :math:`r`
   comes from the dust mass and number mixing ratios assuming a
   log-normal size distribution. The viscosity :math:`\eta` follows
   Sutherland's law,
   :math:`\eta = \eta_0 (T/T_0)^{3/2} (T_0 + S)/(T + S)`, with the
   CO\ :sub:`2` values :math:`\eta_0` = 1.37×10\ :sup:`-5` N s/m\ :sup:`2`,
   :math:`T_0` = 273.15 K, :math:`S` = 222 K.
   Needs ``dst_mass_mom``, ``dst_num_mom`` (or the ``_micro``
   versions), ``temp``. All grids. Function: ``compute_Vg_sed``.
   References: Kasten (1968), White (1991).

``w_net`` -- Net vertical velocity of dust [m/s]
   :math:`w - V_g`. Needs ``Vg_sed``, ``w``. All grids.

``dustref_per_pa`` -- Visible dust opacity per pascal [Pa\ :sup:`-1`]
   ``dustref/delp``, the layer opacity divided by the layer pressure
   thickness. Needs ``dustref``, ``delp``. Grid: ``pfull``.

``dustref_per_km`` -- Visible dust opacity per kilometer [km\ :sup:`-1`]
   ``dustref/delz``, the layer opacity divided by the layer thickness.
   Needs ``dustref``, ``delz``. Grid: ``pfull``.

References
----------

Andrews, D. G., Holton, J. R., and Leovy, C. B. (1987), *Middle
Atmosphere Dynamics*, International Geophysics Series, vol. 40,
Academic Press, San Diego.

Fanale, F. P., Salvail, J. R., Banerdt, W. B., and Saunders, R. S.
(1982), Mars: The regolith-atmosphere-cap system and climate change,
*Icarus*, 50, 381-407, https://doi.org/10.1016/0019-1035(82)90131-2

Heavens, N. G., et al. (2010), Water ice clouds over the Martian tropics
during northern summer, *Geophys. Res. Lett.*, 37, L18202,
https://doi.org/10.1029/2010GL044610

Heavens, N. G., et al. (2011), The vertical distribution of dust in the
Martian atmosphere during northern spring and summer: Observations by
the Mars Climate Sounder and analysis of zonal average vertical dust
profiles, *J. Geophys. Res.*, 116, E04003,
https://doi.org/10.1029/2010JE003691

Holton, J. R., and Hakim, G. J. (2013), *An Introduction to Dynamic
Meteorology*, 5th ed., Academic Press,
https://doi.org/10.1016/C2009-0-63394-8

Kasten, F. (1968), Falling speed of aerosol particles, *J. Appl.
Meteor.*, 7, 944-947,
https://doi.org/10.1175/1520-0450(1968)007<0944:FSOAP>2.0.CO;2

Kleinböhl, A., Schofield, J. T., Abdou, W. A., Irwin, P. G. J., and de
Kok, R. J. (2011), A single-scattering approximation for infrared
radiative transfer in limb geometry in the Martian atmosphere, *J.
Quant. Spectrosc. Radiat. Transfer*, 112, 1568-1580,
https://doi.org/10.1016/j.jqsrt.2011.03.006

Richter, J. H., Sassi, F., and Garcia, R. R. (2010), Toward a physically
based gravity wave source parameterization in a general circulation
model, *J. Atmos. Sci.*, 67, 136-156,
https://doi.org/10.1175/2009JAS3112.1

Scorer, R. S. (1949), Theory of waves in the lee of mountains, *Q. J. R.
Meteorol. Soc.*, 75, 41-56, https://doi.org/10.1002/qj.49707532308

White, F. M. (1991), *Viscous Fluid Flow*, 2nd ed., McGraw-Hill,
Table 1-2.
