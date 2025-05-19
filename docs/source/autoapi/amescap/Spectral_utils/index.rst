:py:mod:`amescap.Spectral_utils`
================================

.. py:module:: amescap.Spectral_utils

.. autoapi-nested-parse::

   Spectral_utils contains wave analysis routines. Note the dependencies on
   scipy.signal.

   Third-party Requirements:

       * ``numpy``
       * ``amescap.Script_utils``
       * ``scipy.signal``
       * ``ImportError``
       * ``Exception``
       



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   amescap.Spectral_utils.diurn_extract
   amescap.Spectral_utils.reconstruct_diurn
   amescap.Spectral_utils.space_time
   amescap.Spectral_utils.zeroPhi_filter


.. py:function:: diurn_extract(VAR, N, tod, lon)

   Extract the diurnal component of a field. Original code by John
   Wilson. Adapted by Alex Kling April, 2021

   :param VAR: field to process. Time of day dimension must be first
       (e.g., ``[tod, time, lat, lon]`` or ``[tod]``
   :type VAR: 1D or ND array

   :param N: number of harmonics to extract (``N=1`` for diurnal,
       ``N=2`` for diurnal AND semi diurnal, etc.)
   :type N: int

   :param tod: universal time of day in sols (``0-1.``) If provided in
       hours (``0-24``), it will be normalized.
   :type tod: 1D array

   :param lon: longitudes ``0-360``
   :type lon: 1D array or float

   :return: the amplitudes & phases for the Nth first harmonics,
       (e.g., size ``[Nh, time, lat, lon]``)
   :rtype: ND arrays



.. py:function:: reconstruct_diurn(amp, phas, tod, lon, sumList=[])

   Reconstructs a field wave based on its diurnal harmonics

   :param amp: amplitude of the signal. Harmonics dimension FIRST
       (e.g., ``[N, time, lat, lon]``)
   :type amp: array

   :param phas: phase of the signal [hr]. Harmonics dimension FIRST
   :type phas: array

   :param tod: time of day in universal time [hr]
   :type tod: 1D array

   :param lon: longitude for converting universal -> local time
   :type lon: 1D array or float

   :param sumList: the harmonics to include when reconstructing the
       wave (e.g., ``sumN=[1, 2, 4]``), defaults to ``[]``
   :type sumList: list, optional

   :return: a variable with reconstructed harmonics with N dimension
       FIRST and time of day SECOND (``[N, tod, time, lat, lon]``). If
       sumList is provided, the wave output harmonics will be
       aggregated (i.e., size = ``[tod, time, lat, lon]``)
   :rtype: _type_



.. py:function:: space_time(lon, timex, varIN, kmx, tmx)

   Obtain west and east propagating waves. This is a Python
       implementation of John Wilson's ``space_time`` routine.
       Alex Kling 2019.

   :param lon: longitude [°] (0-360)
   :type lon: 1D array

   :param timex: time [sol] (e.g., 1.5 days sampled every hour is
       ``[0/24, 1/24, 2/24,.. 1,.. 1.5]``)
   :type timex: 1D array

   :param varIN: variable for the Fourier analysis. First axis must be
       ``lon`` and last axis must be ``time`` (e.g.,
       ``varIN[lon, time]``, ``varIN[lon, lat, time]``, or
       ``varIN[lon, lev, lat, time]``)
   :type varIN: array

   :param kmx: the number of longitudinal wavenumbers to extract
       (max = ``nlon/2``)
   :type kmx: int

   :param tmx: the number of tidal harmonics to extract
       (max = ``nsamples/2``)
   :type tmx: int

   :return: (ampe) East propagating wave amplitude [same unit as
       varIN]; (ampw) West propagating wave amplitude [same unit as
       varIN]; (phasee) East propagating phase [°]; (phasew) West
       propagating phase [°]

   .. NOTE::   1. ``ampe``, ``ampw``, ``phasee``, and ``phasew`` have
       dimensions ``[kmx, tmx]`` or ``[kmx, tmx, lat]`` or
       ``[kmx, tmx, lev, lat]`` etc.

       2. The x and y axes may be constructed as follows, which will
       display the eastern and western modes::

           klon = np.arange(0, kmx) # [wavenumber] [cycle/sol]
           ktime = np.append(-np.arange(tmx, 0, -1), np.arange(0, tmx))
           KTIME, KLON = np.meshgrid(ktime, klon)
           amplitude = np.concatenate((ampw[:, ::-1], ampe), axis = 1)
           phase = np.concatenate((phasew[:, ::-1], phasee), axis = 1)
           


.. py:function:: zeroPhi_filter(VAR, btype, low_highcut, fs, axis=0, order=4, add_trend=False)

   A temporal filter that uses a forward and backward pass to prevent
   phase shift. Alex Kling 2020.

   :param VAR: values for filtering 1D or ND array. Filtered dimension
       must be FIRST. Adjusts axis as necessary.
   :type VAR: array

   :param btype: filter type (i.e., "low", "high" or "band")
   :type btype: str

   :param low_high_cut: low, high, or [low, high] cutoff frequency
       depending on the filter [Hz or m-1]
   :type low_high_cut: int

   :param fs: sampling frequency [Hz or m-1]
   :type fs: int

   :param axis: if data is an ND array, this identifies the filtering
       dimension
   :type axis: int

   :param order: order for the filter
   :type order: int

   :param add_trend: if True, return the filtered output. If false,
       return the trend and filtered output.
   :type add_trend: bool

   :return: the filtered data

   .. NOTE:: ``Wn=[low, high]`` are expressed as a function of the
       Nyquist frequency
       
