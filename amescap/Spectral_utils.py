# !/usr/bin/env python3
"""
Spectral_utils contains wave analysis routines. Note the dependencies on
scipy.signal.

Third-party Requirements:

    * ``numpy``
    * ``amescap.Script_utils``
    * ``scipy.signal``
    * ``ImportError``
    * ``Exception``
    
"""

# Load generic Python modules
import numpy as np

from amescap.Script_utils import Yellow, Cyan, Nclr, progress

try:
    from scipy.signal import butter, filtfilt, detrend
except ImportError as error_msg:
    print(f"{Yellow}Error while importing modules from scipy.signal{Nclr}")
    exit()
except Exception as exception:
    # Output unexpected Exceptions.
    print(exception.__class__.__name__ + ": ", exception)
    exit()

# ======================================================================
#                           DEFINITIONS
# ======================================================================
# Try to import pyshtools with proper error handling
try:
    import pyshtools
    PYSHTOOLS_AVAILABLE = True
except ImportError:
    PYSHTOOLS_AVAILABLE = False
    print(
        f"{Yellow}__________________\n"
        f"Zonal decomposition relies on the pyshtools library, "
        f"referenced at:\n\n"
        f"Mark A. Wieczorek and Matthias Meschede (2018). "
        f"SHTools - Tools for working with spherical harmonics,"
        f"Geochemistry, Geophysics, Geosystems, 2574-2592, "
        f"doi:10.1029/2018GC007529\n\nPlease consult pyshtools  "
        f"documentation at:\n"
        f" {Cyan}https://pypi.org/project/pyshtools\n"
        f"{Yellow}And installation instructions for CAP with pyshtools:\n"
        f" {Cyan}https://amescap.readthedocs.io/en/latest/installation."
        f"html#_spectral_analysis{Yellow}\n"
        f"__________________{Nclr}\n\n"
    )

def diurn_extract(VAR, N, tod, lon):
    """
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
    
    """
    dimsIN = VAR.shape
    nsteps = len(tod)
    period = 24
    rnorm = 1/nsteps
    delta = period/nsteps

    # Be sure that the local time grid units = sols
    if max(tod) > 1:
        tod = tod/24.

    freq = (delta/period) * 2.0*np.pi
    arg = tod * 2*np.pi
    # Reshape array for matrix operations
    arg = arg.reshape([len(tod), 1])

    # Dimensions for the output
    if len(dimsIN) == 1:
        # if size(VAR) = [tod, time, lat, lon] then
        # size(dimsOUT) = [N, time, lat, lon]
        dimsOUT=[N, 1]
    else:
        dimsOUT=np.append([N], dimsIN[1:])

    # Reshape VAR as a 2D array (tod, Nelements) for generalization
    # Ndim is the product of all dimensions except tod axis
    # (e.g., ``time x lat x lon``)
    Ndim = int(np.prod(dimsIN[1:]))
    # Set the shape of the flattened arrays
    dimsFLAT = np.append([nsteps], [Ndim])
    dimsOUT_flat = np.append([N], [Ndim])
    # Flatten array to (tod, Nelements)
    VAR = VAR.reshape(dimsFLAT)

    # Initialize output arrays
    amp = np.zeros(dimsOUT_flat)
    phas = np.zeros(dimsOUT_flat)

    if len(dimsIN) == 1:
        # If ``nargin > 3`` (Initial code, assume lon is provided)
        corr = np.array([lon])
    else:
        # Repeat lon array to match the size of the input variables
        # minus the first (tod) axis.
        # Create an axis to expand the lon array
        tilenm = np.append(dimsIN[1:-1], 1)
        # If VAR = [tod, time, lat, lon] then lonN is the lon repeated
        # as [time,lat,lon]
        lonND = np.tile(lon, tilenm)
        corr = lonND.flatten()

    for nn in range(1,N+1):
        # Python does zero-indexing, start at nn-1
        cosser = np.dot(VAR[:, ...].T, np.cos(nn * arg)).squeeze()
        sinser = np.dot(VAR[:, ...].T, np.sin(nn * arg)).squeeze()

        amp[nn-1, :] = 2 * rnorm * np.sqrt(cosser**2 + sinser**2)
        phas[nn-1, :] = (180/np.pi) * np.arctan2(sinser, cosser)

        # Apply local time correction to the phase
        phas[nn-1, :] = phas[nn-1, :] + 360 + nn*corr[:]
        phas[nn-1, :] = (24/(nn)/360) * np.mod(phas[nn-1, :], 360)

    # Return the phase and amplitude
    return amp.reshape(dimsOUT), phas.reshape(dimsOUT)

def reconstruct_diurn(amp, phas, tod, lon, sumList=[]):
    """
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
    
    """
    dimsIN = amp.shape
    N = dimsIN[0]
    dimsSUM = np.append([len(tod)], dimsIN[1:])
    # Create dimensional array for 4D variables (e.g., [0, 1, 2, 3])
    dimAXIS = np.arange(len(dimsSUM))
    dimsOUT = np.append([dimsIN[0], len(tod)], dimsIN[1:])
    varOUT = np.zeros(dimsOUT)

    # Special case for station data (lon = float)
    if len(np.atleast_1d(lon)) == 1:
        lon = np.array([lon])

    #Reshape  lon array for broadcasting, e.g. lon[96] -> [1,1,1,96]
    dimAXIS = np.arange(len(dimsSUM))
    dimAXIS[:] = 1
    dimAXIS[-1] = len(lon)
    lon = lon.reshape(dimAXIS)
    
    # Reshape tod array
    dimAXIS = np.arange(len(dimsSUM))
    dimAXIS[:] = 1
    dimAXIS[0] = len(tod)
    tod = tod.reshape(dimAXIS)

    # Shift in phase due to local time
    DT = lon/360 * 24

    varSUM  = np.zeros(dimsSUM)
    for nn in range(1, N+1): # Compute each harmonic
        varOUT[nn-1, ...] = (
            amp[nn-1, ...] 
            * np.cos(nn * (tod-phas[nn-1, ...] + DT) / 24*2*np.pi)
            )
        # If a sum of harmonics is requested, sum it
        if nn in sumList:
            varSUM += varOUT[nn-1, ...]

    if sumList:
        # Return the aggregated harmonic
        return varSUM
    else:
        # Return harmonics separately
        return varOUT

def space_time(lon, timex, varIN, kmx, tmx):
    """
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
        ``[kmx, tmx, lev, lat]`` etc.\n
        2. The x and y axes may be constructed as follows, which will
        display the eastern and western modes::

            klon = np.arange(0, kmx) # [wavenumber] [cycle/sol]
            ktime = np.append(-np.arange(tmx, 0, -1), np.arange(0, tmx))
            KTIME, KLON = np.meshgrid(ktime, klon)
            amplitude = np.concatenate((ampw[:, ::-1], ampe), axis = 1)
            phase = np.concatenate((phasew[:, ::-1], phasee), axis = 1)
            
    """
    # Get input variable dimensions
    dims = varIN.shape
    lon_id = dims[0]
    time_id = dims[-1]

    # Additional dimensions stacked in the middle
    dim_sup_id = dims[1:-1]

    # jd = total number of dimensions in the middle (``varIN > 3D``)
    jd = int(np.prod(dim_sup_id))

    # Flatten the middle dimensions, if any
    varIN = np.reshape(varIN, (lon_id, jd, time_id))

    # Initialize 4 empty arrays
    ampw, ampe, phasew, phasee = ([np.zeros((kmx, tmx, jd)) for _x
                                   in range(0, 4)])

    #TODO not implemented yet:
    # zamp, zphas = [np.zeros((jd, tmx)) for _x in range(0, 2)]

    tpi = 2*np.pi
    # Normalize longitude array
    argx = lon * 2*np.pi / 360
    rnorm = 2. / len(argx)
    # If timex = [0/24, 1/24, 2/24,.. 1] arg cycles for m [0, 2Pi]
    arg =  timex * 2*np.pi
    # Nyquist cut off
    rnormt =  2. / len(arg)

    for kk in range(0, kmx):
        progress(kk, kmx)
        cosx = np.cos(kk * argx) * rnorm
        sinx = np.sin(kk * argx) * rnorm

        # Inner product calculates the Fourier coefficients of the
        # cosine and sine contributions of the spatial variation
        acoef = np.dot(varIN.T, cosx)
        bcoef = np.dot(varIN.T, sinx)

        for nn in range(0, tmx):
            # Get the cosine and sine series expansions of the temporal
            # variations of the acoef and bcoef spatial terms
            cosray = rnormt * np.cos(nn * arg)
            sinray = rnormt * np.sin(nn * arg)

            cosA = np.dot(acoef.T, cosray)
            sinA = np.dot(acoef.T, sinray)
            cosB = np.dot(bcoef.T, cosray)
            sinB = np.dot(bcoef.T, sinray)

            wr = 0.5*(cosA - sinB)
            wi = 0.5*(-sinA - cosB)
            er = 0.5*(cosA + sinB)
            ei = 0.5*(sinA - cosB)

            aw = np.sqrt(wr**2 + wi**2)
            ae = np.sqrt(er**2 + ei**2)
            pe = np.arctan2(ei, er) * 180/np.pi
            pw = np.arctan2(wi, wr) * 180/np.pi

            pe = np.mod(-np.arctan2(ei, er) + tpi, tpi) * 180/np.pi
            pw = np.mod(-np.arctan2(wi, wr) + tpi, tpi) * 180/np.pi

            ampw[kk, nn, :] = aw.T
            ampe[kk, nn, :] = ae.T
            phasew[kk, nn, :] = pw.T
            phasee[kk, nn, :] = pe.T
    #End loop

    ampw = np.reshape(ampw, (kmx, tmx) + dim_sup_id)
    ampe = np.reshape(ampe, (kmx, tmx) + dim_sup_id)
    phasew = np.reshape(phasew, (kmx, tmx) + dim_sup_id)
    phasee = np.reshape(phasee, (kmx, tmx) + dim_sup_id)

    # TODO implement zonal mean: zamp, zphas (standing wave k = 0,
    # zonally  averaged) stamp, stphs (stationary component ktime = 0)

    # # varIN = reshape(varIN, dims)
    # # if nargout < 5:
    # #     # only ampe, ampw, phasee, phasew are requested
    # #     return

    # #   Now calculate the axisymmetric tides  zamp,zphas

    # zvarIN = np.mean(varIN, axis=0)
    # zvarIN = np.reshape(zvarIN, (jd, time_id))

    # arg = timex * 2*np.pi
    # arg = np.reshape(arg, (len(arg), 1))
    # rnorm = 2/time_id

    # for nn in range(0, tmx):
    #     cosray = rnorm * np.cos(nn*arg)
    #     sinray = rnorm * np.sin(nn*arg)
    #     cosser =  np.dot(zvarIN, cosray)
    #     sinser =  np.dot(zvarIN, sinray)

    #     zamp[:, nn] = np.sqrt(cosser[:]**2 + sinser[:]**2).T
    #     zphas[:, nn] = np.mod(-np.arctan2(sinser, cosser)
    #                           + tpi, tpi).T * 180/np.pi

    # zamp = zamp.T
    # zphas = zphas.T

    # if len(dims) > 2:
    #     zamp = np.reshape(zamp, (tmx,) + dim_sup_id)
    #     zamp = np.reshape(zphas, (tmx,) + dim_sup_id)

    # # if nargout < 7:
    # #   return

    # sxx = np.mean(varIN, ndims(varIN))
    # [stamp, stphs] = amp_phase(sxx, lon, kmx)

    # if len(dims) > 2:
    #     stamp = reshape(stamp, [kmx dims(2:end-1)])
    #     stphs = reshape(stphs, [kmx dims(2:end-1)])

    return ampe, ampw, phasee, phasew

def zeroPhi_filter(VAR, btype, low_highcut, fs, axis=0, order=4,
                   add_trend=False):
    """
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
        
    """
    # Create the filter
    low_highcut = np.array(low_highcut)
    nyq = 0.5*fs
    b, a = butter(order, low_highcut/nyq, btype = btype)

    # Detrend the data, this is the equivalent of doing linear
    # regressions across the time axis at each grid point
    VAR_detrend = detrend(VAR, axis = axis, type = "linear")

    # Trend = variable - detrend array
    VAR_trend = VAR - VAR_detrend

    VAR_f = filtfilt(b, a, VAR_detrend, axis = axis)

    if add_trend:
        return VAR_trend + VAR_f
    else:
        return VAR_f

def zonal_decomposition(VAR):
    """
    Decomposition into spherical harmonics. [A. Kling, 2020]

    :param VAR: Detrend variable for decomposition. Lat is SECOND to
        LAST dimension and lon is LAST (e.g., ``[time,lat,lon]`` or
        ``[time,lev,lat,lon]``)

    :return: (COEFFS) coefficient for harmonic decomposion, shape is
        flattened (e.g., ``[time, 2, lat/2, lat/2]``
        ``[time x lev, 2, lat/2, lat/2]``);
        (power_per_l) power spectral density, shape is re-organized
        (e.g., [time, lat/2] or [time, lev, lat/2])

    .. NOTE:: Output size is (``[...,lat/2, lat/2]``) as lat is the
        smallest dimension. This matches the Nyquist frequency.
        
    """
    if not PYSHTOOLS_AVAILABLE:
        raise ImportError(
            "This function requires pyshtools. Install with:\n"
            "conda install -c conda-forge pyshtools\n"
            "or\n"
            "pip install amescap[spectral]"
        )

    var_shape = np.array(VAR.shape)

    # Flatten array (e.g., [10, 36, lat, lon] -> [360, lat, lon])
    nflatten = int(np.prod(var_shape[:-2]))
    reshape_flat = np.append(nflatten, var_shape[-2:])
    VAR = VAR.reshape(reshape_flat)

    coeff_out_shape = np.append(var_shape[0:-2],
                                [2, var_shape[-2]//2, var_shape[-2]//2])
    psd_out_shape = np.append(var_shape[0:-2], var_shape[-2]//2)
    coeff_flat_shape = np.append(nflatten, coeff_out_shape[-3:])
    COEFFS = np.zeros(coeff_flat_shape)

    psd = np.zeros((nflatten,var_shape[-2]//2))

    for ii in range(0,nflatten):
        COEFFS[ii,...] = pyshtools.expand.SHExpandDH(VAR[ii,...], sampling = 2)
        psd[ii,:] = pyshtools.spectralanalysis.spectrum(COEFFS[ii,...])

    return  COEFFS, psd.reshape(psd_out_shape)

def zonal_construct(COEFFS_flat, VAR_shape, btype=None, low_highcut=None):
    """
    Recomposition into spherical harmonics

    :param COEFFS_flat: Spherical harmonic coefficients as a flattened
        array, (e.g., ``[time, lat, lon]`` or
        ``[time x lev, 2, lat, lon]``)
    :type COEFFS_flat: array
    
    :param VAR_shape: shape of the original variable
    :type VAR_shape: tuple
    
    :param btype: filter type: "low", "high", or "band". If None,
        returns reconstructed array using all zonal wavenumbers
    :type btype: str or None
    
    :param low_high_cut: low, high or [low, high] zonal wavenumber
    :type low_high_cut: int

    :return: detrended variable reconstructed to original size
        (e.g., [time, lev, lat, lon])

    .. NOTE:: The minimum and maximum wavelenghts in [km] are computed::
        dx = 2*np.pi * 3400
        L_min = (1./kmax) * dx
        L_max = 1./max(kmin, 1.e-20) * dx
        if L_max > 1.e20:
            L_max = np.inf
        print("(kmin,kmax) = ({kmin}, {kmax})
              -> dx min = {L_min} km, dx max = {L_max} km")
              
    """
    if not PYSHTOOLS_AVAILABLE:
        raise ImportError(
            "This function requires pyshtools. Install with:\n"
            "conda install -c conda-forge pyshtools\n"
            "or\n"
            "pip install amescap[spectral]"
        )

    # Initialization
    nflatten = COEFFS_flat.shape[0]
    kmin = 0
    kmax = COEFFS_flat.shape[-1]

    VAR = np.zeros((nflatten, VAR_shape[-2], VAR_shape[-1]))

    if btype == "low":
        kmax= int(low_highcut)
    if btype == "high":
        kmin= int(low_highcut)
    if btype == "band":
        kmin, kmax= int(low_highcut[0]), int(low_highcut[1])

    for ii in range(0, nflatten):
        # Filtering
        COEFFS_flat[ii, :, :kmin, :] = 0.
        COEFFS_flat[ii, :, kmax:, :] = 0.
        VAR[ii, :, :] = pyshtools.expand.MakeGridDH(COEFFS_flat[ii, :, :],
                                                  sampling = 2)
    return  VAR.reshape(VAR_shape)


def init_shtools():
    """
    Check if pyshtools is available and return its availability status
    
    :return: True if pyshtools is available, False otherwise
    :rtype: bool
    """
    return PYSHTOOLS_AVAILABLE