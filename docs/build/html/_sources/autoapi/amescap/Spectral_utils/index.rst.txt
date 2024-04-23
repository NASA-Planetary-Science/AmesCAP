:orphan:

:py:mod:`amescap.Spectral_utils`
================================

.. py:module:: amescap.Spectral_utils


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   amescap.Spectral_utils.init_shtools
   amescap.Spectral_utils.diurn_extract
   amescap.Spectral_utils.reconstruct_diurn
   amescap.Spectral_utils.space_time
   amescap.Spectral_utils.zeroPhi_filter
   amescap.Spectral_utils.zonal_decomposition
   amescap.Spectral_utils.zonal_construct



.. py:function:: init_shtools()

   The following code simply loads the pyshtools module and provides adequate referencing.
   Since dependencies may need to be solved by the user, the module import is wrapped in a function
   that may be called when needed.


.. py:function:: diurn_extract(VAR, N, tod, lon)

   Extract the diurnal component of a field. Original code by J.Wilson adapted by A. Kling. April, 2021
   Args:
       VAR (1D or ND array)   : field to process with time of day dimension FIRST, e.g (tod,time,lat,lon) or (tod)
       N   (int)              : number of harmonics to extract (N=1 for diurnal,N=2  for diurnal + semi diurnal etc...)
       tod (1D array)         : universal time of day in sols (0>1.) If provided in hours (0>24), it  will  be  normalized.
       lon (1D array or float): longitudes 0>360
   Return:
       amp (ND array)  :  the amplitudes for the Nth first harmonics, e.g. size (Nh,time,lat,lon)
       phase (ND array):  the phases for the Nth first harmonics, e.g. size (Nh,time,lat,lon)


.. py:function:: reconstruct_diurn(amp, phas, tod, lon, sumList=[])

   Reconstruct a field wave based on its diurnal harmonics
   Args:
       amp   : amplitude of the signal, with  harmonics dimension FIRST, e.g. (N,time,lat,lon)
       phas : phase of the signal, in [hr], with harmonics dimension FIRST
       tod   : 1D array for the time of day, in UT [hr]
       lon   : 1D array  or float for the longitudes, used to convert UT to LT
       sumList : (optional) list containing the harmonics to include when reconstructing the wave, e.g. sumN=[1,2,4]
   Return:
       VAR   : a variable with reconstructed harmonics with N dimension FIRST and time of day SECOND, e.g. (N,tod,time,lat,lon)
               if  sumList is provided, the wave output has the harmonics already agregated, e.g. size is    (tod,time,lat,lon)


.. py:function:: space_time(lon, timex, varIN, kmx, tmx)

   Obtain west and east propagating waves. This is a Python implementation of John Wilson's  space_time routine by [A. Kling, 2019]
   Args:
       lon:   longitude array in [degrees]   0->360
       timex: 1D time array in units of [day]. Expl 1.5 days sampled every hour is  [0/24,1/24, 2/24,.. 1,.. 1.5]
       varIN: input array for the Fourier analysis.
              First axis must be longitude and last axis must be time.  Expl: varIN[lon,time] varIN[lon,lat,time],varIN[lon,lev,lat,time]
       kmx: an integer for the number of longitudinal wavenumber to extract   (max allowable number of wavenumbers is nlon/2)
       tmx: an integer for the number of tidal harmonics to extract           (max allowable number of harmonics  is nsamples/2)

   Returns:
       ampe:   East propagating wave amplitude [same unit as varIN]
       ampw:   West propagating wave amplitude [same unit as varIN]
       phasee: East propagating phase [degree]
       phasew: West propagating phase [degree]



   *NOTE*  1. ampe,ampw,phasee,phasew have dimensions [kmx,tmx] or [kmx,tmx,lat] [kmx,tmx,lev,lat] etc...
           2. The x and y axis may be constructed as follow to display the easter and western modes:

               klon=np.arange(0,kmx)  [wavenumber]  [cycle/sol]
               ktime=np.append(-np.arange(tmx,0,-1),np.arange(0,tmx))
               KTIME,KLON=np.meshgrid(ktime,klon)

               amplitude=np.concatenate((ampw[:,::-1], ampe), axis=1)
               phase=    np.concatenate((phasew[:,::-1], phasee), axis=1)



.. py:function:: zeroPhi_filter(VAR, btype, low_highcut, fs, axis=0, order=4, no_trend=False)

   Temporal filter: use a forward pass and a backward pass to prevent phase shift. [A. Kling, 2020]
   Args:
       VAR:  values to filter 1D or ND array. Filtered dimension is FIRST, otherwise, adjust axis
       btype: filter type: 'low', 'high' or 'band'
       low_high_cut: low , high or [low,high] cutoff frequency depending on the filter [Hz or m-1]
       fs:     sampling frequency [Hz or m-1]
       axis:  if data is N-dimensional array, the filtering dimension
       order: order for the filter
       no_trend: if True, only return the filtered-output, not TREND+ FILTER

   Returns:
       out: the filtered data

   ***NOTE***
   Wn=[low, high] are expressed as a function of the Nyquist frequency


.. py:function:: zonal_decomposition(VAR)

   Decomposition into spherical harmonics. [A. Kling, 2020]
   Args:
       VAR:  Detrend variable for decomposition, latitude is SECOND to LAST and longitude is LAST  e.g. (time,lat,lon) or (time,lev,lat,lon)
   Returns:
       COEFFS      : coefficient for harmonic decomposion, shape is flatten e.g. (time,2,lat/2, lat/2) (time x lev,2,lat/2, lat/2)
       power_per_l : power spectral density, shape is re-organized, e.g. (time, lat/2) or  (time,lev,lat/2)
   ***NOTE***
   Output size is (...,lat/2, lat/2) as latitude is the smallest dimension and to match the Nyquist frequency


.. py:function:: zonal_construct(COEFFS_flat, VAR_shape, btype=None, low_highcut=None)

   Recomposition into spherical harmonics
   Args:
       COEFFS_flat:  Spherical harmonics coefficients as flatten array, e.g. (time,lat,lon) or (time x lev,2, lat,lon)
       VAR_shape:    Shape of the original variable e.g. VAR_shape=temp.shape
       btype: filter type: 'low', 'high' or 'band'. If None, returns array as reconstructed using all zonal wavenumber
       low_high_cut: low , high or [low,high] cutoff zonal wavenumber(s)
   Returns:
       VAR      : reconstructed output, size is same as original detrened variable e.g. (time,lev,lat,lon)

   ***NOTE***
   # The minimum and maximum wavelenghts in [km] may be computed as follows:
   dx=2*np.pi*3400
   L_min=(1./kmax)*dx
   L_max=1./max(kmin,1.e-20)*dx ; if L_max>1.e20:L_max=np.inf
   print('(kmin,kmax)=(#g,#g)>>dx min= #g km,dx max= #g km'#(kmin,kmax,L_min,L_max))


