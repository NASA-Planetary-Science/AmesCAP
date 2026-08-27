planetWRF and MarsWRF Output
============================

``MarsFormat -gcm marswrf`` converts output from planetWRF (and its Mars
configuration, MarsWRF) into the CAP file layout, after which every other
CAP executable works on it unchanged: ``MarsInterp``, ``MarsFiles``,
``MarsVars`` and ``MarsPlot``.

Two kinds of input are supported:

* **Full** ``wrfout`` **files** as written by the model (hundreds of
  variables, staggered winds, eta-coordinate metadata, ``Times`` strings).
* **Reduced files**, such as the subsets published alongside papers, that
  keep only a handful of physics-grid fields (``T_PHY``, ``P_PHY``,
  ``Z_PHY``, ``U_PHY``, ``V_PHY``, ``PSFC``, ``TSK``, ``L_S``, ...) and no
  grid metadata.

Quick start
-----------

.. code-block:: bash

    (amesCAP)$ MarsFormat wrfout_d01_0001-00144_00:00:00 -gcm marswrf
    (amesCAP)$ MarsInterp wrfout_d01_0001-00144_00:00:00_daily.nc -t pstd
    (amesCAP)$ MarsFiles  wrfout_d01_0001-00144_00:00:00_daily.nc -bd -ba 1
    (amesCAP)$ MarsFiles  wrfout_d01_0001-00144_00:00:00_daily_to_average_to_diurn.nc -t "3 15"
    (amesCAP)$ MarsPlot -template

Notes on the flags:

* ``MarsFiles -bd`` averages the diurnal bins over 5 sols by default; use
  ``-ba 1`` (or another period) for short files.
* Files from a restart have physics diagnostics equal to zero on their
  first frame; ``MarsFormat`` detects this and reconstructs that frame
  from the dynamics fields.
* ``MarsFormat`` streams large files in time chunks when ``dask`` is
  installed (``pip install dask``, or the ``[large]`` extra); without it
  the whole file is held in memory. ``MarsInterp`` interpolates in time
  chunks (``-chunk N`` to set the size). A one-year, 7 GB planetWRF file
  converts in about 20 s with a 2 GB peak and interpolates to ``pstd``
  in about 80 s with a 5 GB peak.

Where each quantity comes from
------------------------------

``MarsFormat`` takes each quantity from the best source present in the
file and falls back to older or reconstructed forms when it is missing.
The source used is printed for every file.

.. list-table::
   :header-rows: 1
   :widths: 22 32 46

   * - Quantity
     - First choice
     - Fallbacks
   * - ``temp``
     - ``T_PHY``
     - ``(T + T0) * ((P + PB) / P0) ** (R_D / CP)``
   * - ``pfull3D``
     - ``P_PHY``
     - ``P + PB``
   * - ``time`` [sols]
     - ``Times`` strings (``YYYY-DDDDD_HH:MM:SS``)
     - ``XTIME`` + ``SIMULATION_START_DATE``; ``JULIAN`` + ``MODEL_MARS_YEAR``;
       ``START_DATE`` + uniform spacing inferred from ``L_S``
   * - ``areo`` [Ls]
     - ``L_S``
     - ``sol2ls(time)``
   * - ``ak``, ``bk``, ``pfull``, ``phalf``
     - ``ZNW``, ``ZNU``, ``P_TOP`` (``ak = P_TOP * (1 - eta)``, ``bk = eta``)
     - per-level least-squares fit of ``P_PHY`` against ``PSFC``
   * - ``ucomp``, ``vcomp``, ``w``
     - ``U_PHY``, ``V_PHY``, ``W_PHY`` (model's own mass-point values)
     - two-point destaggering of ``U``, ``V``, ``W``
   * - ``zfull``, ``zsurf``
     - ``Z_PHY - HGT``
     - ``Z_PHY`` minus a hydrostatic surface estimate from the bottom layer
       (no ``HGT``); ``(PH + PHB) / g - HGT``
   * - local time
     - ``LTST`` (true local solar time; saved as ``eot_offset``)
     - mean local time ``UT + lon / 15``
   * - ``rgas3D``, ``cp3D``
     - ``8314.46 / MW_AIR_3D``, ``rgas3D / RCP_3D``
     - global attributes ``R_D``, ``CP``
   * - ``g3D``
     - ``G * (1 - Z_PHY / RADIUS) ** 2`` when ``DO_VARIABLE_GRAVITY = 1``
     - global attribute ``G``

Variables on dimensions CAP cannot carry (soil layers, dust bins,
radiation layers, vertical-only metadata such as ``ZNU`` or ``C1H``) are
dropped with a message.

Planet constants downstream
---------------------------

``MarsVars`` reads the gas constant, heat capacity, gravity and planetary
radius from the file when it carries the WRF attribute set ``R_D``,
``CP``, ``G`` and ``RADIUS`` (all four are required, so files from other
models keep CAP's defaults). The 3D fields ``rgas3D``, ``cp3D`` and
``g3D`` take precedence when present and are used in ``zfull``, ``rho``,
``theta`` and the column integrals.

Native C-grid winds
-------------------

``MarsFormat ... -gcm marswrf -stag`` additionally keeps the winds on their
Arakawa-C locations:

* ``u_stag(time, pfull, lat, lon_u)`` and ``v_stag(time, pfull, lat_v, lon)``,
  with the ``lon_u`` and ``lat_v`` axes from ``XLONG_U`` and ``XLAT_V``;
* ``w_half(time, phalf, lat, lon)`` on the layer interfaces, together with
  the interface pressure ``phalf3D`` (``PF_PHY``) and height ``zhalf``
  (``ZF_PHY``).

``w_half``, ``phalf3D`` and ``zhalf`` go through ``MarsInterp`` like any
interface field. ``u_stag`` and ``v_stag`` are intended for your own
scripts; CAP plots only the mass-point grid, and ``MarsInterp`` skips them.
When they are present, ``MarsVars -add div curl`` forms the divergence and
relative vorticity exactly as the model does on the C grid instead of
using centred differences of the mass-point winds.

True local solar time
---------------------

Recent planetWRF runs with ``ra_use_true_solar_time`` write ``LTST``, the
true local solar time of the model sun. Its departure from mean local
time (``UT + lon / 15``) is the equation of time, a single value per
frame. ``MarsFormat`` saves it as ``eot_offset(time)`` [hr] and
``MarsFiles -t`` applies it sol by sol, so the shifted output is on true
solar time. Files without ``LTST`` are shifted to mean local time, which
for older planetWRF versions is the model's solar time by construction.
