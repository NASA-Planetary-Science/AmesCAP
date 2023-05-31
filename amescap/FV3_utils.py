import numpy as np
from netCDF4 import Dataset
import os           # access operating systems function
import warnings     # suppress certain errors when dealing with NaN arrays

# NOTE p_half = half-level = layer interfaces
# NOTE p_full = full-level = layer midpoints

def fms_press_calc(psfc, ak, bk, lev_type='full'):
    """
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
    """

    Nk = len(ak)
    # If 'psfc' is a float (e.g. psfc = 700.), make it a 1-element array (e.g. psfc = [700])
    if len(np.atleast_1d(psfc)) == 1:
        psfc = np.array([np.squeeze(psfc)])

    # Flatten psfc array to generalize it to N dimensions
    psfc_flat = psfc.flatten()

    # Expand the dimensions. Vectorized calculations:
    psfc_v  = np.repeat(psfc_flat[:, np.newaxis], Nk, axis=1)   # (Np) -> (Np, Nk)
    ak_v    = np.repeat(ak[np.newaxis, :], len(
                        psfc_flat), axis=0)                     # (Nk) -> (Np, Nk)
    bk_v    = np.repeat(bk[np.newaxis, :], 1, axis=0)           # (Nk) -> (1,  Nk)

    # Pressure at layer interfaces. The size of Z axis is 'Nk'
    PRESS_h = psfc_v*bk_v+ak_v

    # Pressure at layer midpoints. The size of Z axis is 'Nk-1'
    PRESS_f = np.zeros((len(psfc_flat), Nk-1))
    # Top layer (1st element is i = 0 in Python)
    if ak[0] == 0 and bk[0] == 0:
        PRESS_f[:, 0] = 0.5*(PRESS_h[:, 0]+PRESS_h[:, 1])
    else:
        PRESS_f[:, 0] = (PRESS_h[:, 1]-PRESS_h[:, 0]) / \
            np.log(PRESS_h[:, 1]/PRESS_h[:, 0])

    # The rest of the column (i = 1 ... Nk).
    # [2:] goes from the 3rd element to 'Nk' and
    # [1:-1] goes from the 2nd element to 'Nk-1'
    PRESS_f[:, 1:] = (PRESS_h[:, 2:]-PRESS_h[:, 1:-1]) / \
        np.log(PRESS_h[:, 2:]/PRESS_h[:, 1:-1])

    # First, transpose PRESS(:, Nk) to PRESS(Nk, :),
    # then reshape PRESS(Nk, :) to the original pressure shape PRESS(Nk, :, :, :) (resp. 'Nk-1')

    if lev_type == "full":
        new_dim_f = np.append(Nk-1, psfc.shape)
        return np.squeeze(PRESS_f.T.reshape(new_dim_f))
    elif lev_type == "half":
        new_dim_h = np.append(Nk, psfc.shape)
        return np.squeeze(PRESS_h.T.reshape(new_dim_h))
    else:
        raise Exception(
            """Pressure level type not recognized by press_lev(): use 'full' or 'half' """)


def fms_Z_calc(psfc, ak, bk, T, topo=0., lev_type='full'):
    """
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

    therefore, T_half[k+1]/Tfull[k] = (p_half[k+1]/p_full[k])**(R/Cp)            =====Thalf=====zhalf[k]  \
                                                                                                           \
                                                                                                            \
    From the lapse rate, assume T decreases linearly within the layer:           -----Tfull-----zfull[k]     \ T(z)= To-Γ (z-zo)
        T_half[k+1] = T_full[k] + Γ(Z_full[k]-Z_half[k+1])                                                    \
        (Tfull < Thalf and Γ > 0)                                                                              \
        Z_full[k] = Z_half[k] + (T_half[k+1]-T_full[k])/Γ                        =====Thalf=====zhalf[k+1]      \
    Pulling out Tfull from above equation and using  Γ = (gγ)/R:
        Z_full[k] = Z_half[k+1] + (R Tfull[k])/(gγ)(T_half[k+1]/T_full[k] - 1)
    Using the isentropic relation above:
        Z_full = Z_half[k+1] + (R Tfull[k])/(gγ)(p_half[k+1]/p_full[k])**(R/Cp)-1)
    """
    g       = 3.72      # acc. m/s2
    r_co2   = 191.00    # kg/mol
    Nk      = len(ak)

    # Get the half and full pressure levels from fms_press_calc
    PRESS_f = fms_press_calc(psfc, ak, bk, 'full')  # Z axis is first
    PRESS_h = fms_press_calc(psfc, ak, bk, 'half')  # Z axis is first

    # If 'psfc' is a float, turn it into a 1-element array:
    if len(np.atleast_1d(psfc)) == 1:
        psfc = np.array([np.squeeze(psfc)])
    if len(np.atleast_1d(topo)) == 1:
        topo = np.array([np.squeeze(topo)])

    psfc_flat = psfc.flatten()
    topo_flat = topo.flatten()

    # Reshape arrays for vector calculations and compute the log pressure
    PRESS_h = PRESS_h.reshape((Nk, len(psfc_flat)))
    PRESS_f = PRESS_f.reshape((Nk-1, len(psfc_flat)))
    T = T.reshape((Nk-1, len(psfc_flat)))

    logPPRESS_h = np.log(PRESS_h)

    # Initialize the output arrays
    Z_f = np.zeros((Nk-1, len(psfc_flat)))
    Z_h = np.zeros((Nk, len(psfc_flat)))

    # First half-layer is equal to the surface elevation
    Z_h[-1, :] = topo_flat

    # Other layers from the bottom-up:
    # Isothermal within the layer, we have Z = Z0 + r*T0/g*ln(P0/P)
    for k in range(Nk-2, -1, -1):
        Z_h[k, :] = Z_h[k+1, :]+(r_co2*T[k, :]/g) * \
            (logPPRESS_h[k+1, :]-logPPRESS_h[k, :])
        Z_f[k, :] = Z_h[k+1, :]+(r_co2*T[k, :]/g) * \
            (1-PRESS_h[k, :]/PRESS_f[k, :])

    # Return the arrays
    if lev_type == "full":
        new_dim_f = np.append(Nk-1, psfc.shape)
        return Z_f.reshape(new_dim_f)
    elif lev_type == "half":
        new_dim_h = np.append(Nk, psfc.shape)
        return Z_h.reshape(new_dim_h)
    # Return the levels in Z coordinates [m]
    else:
        raise Exception(
            """Altitude level type not recognized: use 'full' or 'half' """)

# TODO: delete: Former version of find_n() : only provides 1D >1D and ND > 1D mapping


def find_n0(Lfull_IN, Llev_OUT, reverse_input=False):
    '''
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
    '''
    # Number of original layers
    Lfull_IN = np.array(Lfull_IN)
    Nlev = len(np.atleast_1d(Llev_OUT))
    if Nlev == 1:
        Llev_OUT = np.array([Llev_OUT])
    dimsIN = Lfull_IN.shape  # Get input variable dimensions
    Nfull = dimsIN[0]
    dimsOUT = tuple(np.append(Nlev, dimsIN[1:]))
    # 'Ndim' is the product of all the dimensions except for the vertical axis
    Ndim = np.int(np.prod(dimsIN[1:]))
    Lfull_IN = np.reshape(Lfull_IN, (Nfull, Ndim))

    if reverse_input:
        Lfull_IN = Lfull_IN[::-1, :]

    ncol = Lfull_IN.shape[-1]
    n = np.zeros((Nlev, ncol), dtype=int)

    for i in range(0, Nlev):
        for j in range(0, ncol):
            n[i, j] = np.argmin(np.abs(Lfull_IN[:, j]-Llev_OUT[i]))
            if Lfull_IN[n[i, j], j] > Llev_OUT[i]:
                n[i, j] = n[i, j]-1
    return n
# =========================================================================================


def find_n(X_IN, X_OUT, reverse_input=False, modulo=None):
    '''
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

       '''
    # number of original layers
    if type(X_IN) != np.ndarray:
        X_IN = np.array(X_IN)
    # If one value is requested, convert float to array
    if len(np.atleast_1d(X_OUT)) == 1:
        X_OUT = np.array([X_OUT])
    elif type(X_OUT) != np.ndarray:  # Convert list to numpy array as needed
        X_OUT = np.array(X_OUT)

    dimsIN = X_IN.shape  # get input variable dimensions
    dimsOUT = X_OUT.shape  # get output variable dimensions
    N_IN = dimsIN[0]  # Size of interpolation axis
    N_OUT = dimsOUT[0]

    # Number of element in arrays other than interpolation axis
    NdimsIN = np.int(np.prod(dimsIN[1:]))
    NdimsOUT = np.int(np.prod(dimsOUT[1:]))

    if (NdimsIN > 1 and NdimsOUT > 1) and (NdimsIN != NdimsOUT):
        print('*** Error in find_n(): dimensions of arrays other than the interpolated (first) axis must be 1 or identical***')

    # Ndims_IN and Ndims_OUT are either '1' or equal
    Ndim = max(NdimsIN, NdimsOUT)

    X_IN = np.reshape(X_IN, (N_IN, NdimsIN))
    X_OUT = np.reshape(X_OUT, (N_OUT, NdimsOUT))

    # Reverse input array if monotically decreasing
    if reverse_input:
        X_IN = X_IN[::-1, :]

    n = np.zeros((N_OUT, Ndim), dtype=int)

    # Some repetition below but that allows to keep the 'if' statement out of the big loop over all array elements
    if len(dimsIN) == 1:
        for i in range(N_OUT):
            for j in range(Ndim):
                n[i, j] = np.argmin(np.abs(X_OUT[i, j]-X_IN[:]))
                if X_IN[n[i, j]] > X_OUT[i, j]:
                    n[i, j] = n[i, j]-1
    elif len(dimsOUT) == 1:
        for i in range(N_OUT):
            for j in range(Ndim):
                n[i, j] = np.argmin(np.abs(X_OUT[i]-X_IN[:, j]))
                if X_IN[n[i, j], j] > X_OUT[i]:
                    n[i, j] = n[i, j]-1
    else:
        for i in range(N_OUT):
            for j in range(Ndim):
                n[i, j] = np.argmin(np.abs(X_OUT[i, j]-X_IN[:, j]))
                if X_IN[n[i, j], j] > X_OUT[i, j]:
                    n[i, j] = n[i, j]-1

    if len(dimsOUT) == 1:
        n = np.squeeze(n)
    return n


def expand_index(Nindex, VAR_shape_axis_FIRST, axis_list):
    '''
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

    '''
    # If one element, turn axis to list
    if len(np.atleast_1d(axis_list)) == 1:
        axis_list = [axis_list]

    # size for the interpolation, e.g. (tod, time x lev x lat x lon)
    Nfull = Nindex.shape[0]
    # Desired output size with LEV axis repeated(tod, time x lev x lat x lon)
    dimsOUT_flat = tuple(np.append(Nfull, np.prod(VAR_shape_axis_FIRST[1:])))

    dimsIN = []  # Reconstruct  the initial (un-flattened) size of Nindex
    # using the dimenions from VAR_shape_axis_FIRST
    for ii, len_axis in enumerate(VAR_shape_axis_FIRST):
        if ii > 0 and not ii in axis_list:
            # skip the first (interpolated) axis and add to the list of initial dims unless the axis are the ones we are expending.
            dimsIN.append(len_axis)
    dimsIN = np.insert(dimsIN, 0, Nfull)  # Initial shape for :Nindex
    # Reshape Nindex from its iniatial flattened sahpe (Nfull,time x lat x lon) to a  ND array (tod, time , lat , lon)
    Nindex = np.reshape(Nindex, dimsIN)
    # Repeat the interpolation indices on the requested axis
    for ii, len_axis in enumerate(VAR_shape_axis_FIRST):
        if ii in axis_list:
            Nindex = np.repeat(np.expand_dims(Nindex, axis=ii), len_axis, ii)
    #e.g. Nindex is now (tod, time ,LEV, lat , lon)
    # Return the new, flattened version of Nindex
    return Nindex.reshape(dimsOUT_flat)


def vinterp(varIN, Lfull, Llev, type_int='log', reverse_input=False, masktop=True, index=None):
    '''
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


    '''
    # Special case where only 1 layer is requested
    Nlev = len(np.atleast_1d(Llev))
    if Nlev == 1:
        Llev = np.array([Llev])

    dimsIN = varIN.shape  # get input variable dimensions
    Nfull = dimsIN[0]

    # Special case where varIN and Lfull are a single profile
    if len(varIN.shape) == 1:
        varIN = varIN.reshape([Nfull, 1])
    if len(Lfull.shape) == 1:
        Lfull = Lfull.reshape([Nfull, 1])

    dimsIN = varIN.shape  # repeat in case varIN and Lfull were reshaped

    dimsOUT = tuple(np.append(Nlev, dimsIN[1:]))
    # Ndim is the product  of all dimensions but the vertical axis
    Ndim = np.int(np.prod(dimsIN[1:]))
    # flatten the other dimensions to (Nfull, Ndim)
    varIN = np.reshape(varIN, (Nfull, Ndim))
    # flatten the other dimensions to (Nfull, Ndim)
    Lfull = np.reshape(Lfull, (Nfull, Ndim))
    varOUT = np.zeros((Nlev, Ndim))
    Ndimall = np.arange(0, Ndim)  # all indices (does not change)

    #
    if reverse_input:
        Lfull = Lfull[::-1, :]
        varIN = varIN[::-1, :]

    for k in range(0, Nlev):
        # Find nearest layer to Llev[k]
        if np.any(index):
            # index have been pre-computed:
            n = index[k, :]
        else:
            # Compute index on the fly for that layer.
            # Note that reversed_input is always set to False as if desired, Lfull was reversed earlier
            n = np.squeeze(find_n(Lfull, Llev[k], False))
        # ==Slower method (but explains what is done below): loop over Ndim======
        # for ii in range(Ndim):
        #     if n[ii]<Nfull-1:
        #         alpha=np.log(Llev[k]/Lfull[n[ii]+1,ii])/np.log(Lfull[n[ii],ii]/Lfull[n[ii]+1,ii])
        #         varOUT[k,ii]=varIN[n[ii],ii]*alpha+(1-alpha)*varIN[n[ii]+1,ii]

        # =================    Fast method  no loop  =======================
        # Convert the layers n to indexes, for a 2D matrix using nindex=i*ncol+j
        nindex = n*Ndim+Ndimall  # n
        nindexp1 = (n+1)*Ndim+Ndimall  # n+1

        # initialize alpha, size is [Ndim]
        alpha = np.NaN*Ndimall
        # Only calculate alpha  where the indices are <Nfull
        Ndo = Ndimall[nindexp1 < Nfull*Ndim]
        if type_int == 'log':
            alpha[Ndo] = np.log(Llev[k]/Lfull.flatten()[nindexp1[Ndo]])/np.log(
                Lfull.flatten()[nindex[Ndo]]/Lfull.flatten()[nindexp1[Ndo]])
        elif type_int == 'lin':
            alpha[Ndo] = (Llev[k]-Lfull.flatten()[nindexp1[Ndo]]) / \
                (Lfull.flatten()[nindex[Ndo]] - Lfull.flatten()[nindexp1[Ndo]])

        # Mask if Llev[k]<model top for the pressure interpolation
        if masktop:
            alpha[Llev[k] < Lfull.flatten()[nindex]] = np.NaN

        # Here, we need to make sure n+1 is never> Nfull by setting n+1=Nfull, if it is the case.
        # This does not affect the calculation as alpha is set to NaN for those values.
        nindexp1[nindexp1 >= Nfull*Ndim] = nindex[nindexp1 >= Nfull*Ndim]

        varOUT[k, :] = varIN.flatten()[nindex]*alpha+(1-alpha) * \
            varIN.flatten()[nindexp1]

    return np.reshape(varOUT, dimsOUT)


def axis_interp(var_IN, x, xi, axis, reverse_input=False, type_int='lin', modulo=None):
    '''
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
    '''
    # Convert list to numpy array as needed
    if type(var_IN) != np.ndarray:
        var_IN = np.array(var_IN)
    # Move interpolated axis to 1st axis:
    var_IN = np.moveaxis(var_IN, axis, 0)
    if reverse_input:
        var_IN = var_IN[::-1, ...]
        x = x[::-1]

    # This is called everytime as it is fast on a 1D array
    index = find_n(x, xi, False)

    dimsIN = var_IN.shape
    dimsOUT = tuple(np.append(len(xi), dimsIN[1:]))
    print(dimsOUT)
    var_OUT = np.zeros(dimsOUT)

    for k in range(0, len(index)):
        n = index[k]
        np1 = n+1
        # Treatment of edge cases where the interpolated value is outside the domain, i.e. n is the last element  and n+1 does not exist
        if np1 >= len(x):
            # If looping around (e.g. longitude, time of day...)replace n+1 by the first element
            if modulo is not None:
                np1 = 0
            else:
                # This will set the interpolated value to NaN in xi as last value  as x[n] - x[np1] =0
                np1 -= 1
        # Also set n=n+1 (which results in NaN) if n =-1 (requested value is samller than first element array) and the values are NOT cyclic
        if n == -1 and modulo is None:
            n = 0
        if type_int == 'log':
            alpha = np.log(xi[k]/x[np1])/np.log(x[n]/x[np1])
        elif type_int == 'lin':
            if modulo is None:
                alpha = (xi[k]-x[np1])/(x[n] - x[np1])
            else:
                alpha = np.mod(xi[k]-x[np1]+modulo, modulo) / \
                    np.mod(x[n] - x[np1]+modulo, modulo)
        var_OUT[k, :] = var_IN[n, ...]*alpha+(1-alpha)*var_IN[np1, ...]

    return np.moveaxis(var_OUT, 0, axis)


def layers_mid_point_to_boundary(pfull, sfc_val):
    '''
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
    '''

    def lambertW_approx(x):
        # Internal function. Uniform approximation for the product log-function
        A = 2.344
        B = 0.8842
        C = 0.9294
        D = 0.5106
        E = -1.213
        y = np.sqrt(2*np.e*x+2)
        return (2*np.log(1+B*y)-np.log(1+C*np.log(1+D*y))+E)/(1+1./(2*np.log(1+B*y)+2*A))

    N = len(pfull)
    phalf = np.zeros(N+1)
    phalf[N] = sfc_val

    for i in range(N, 0, -1):
        B = phalf[i]-pfull[i-1]*np.log(phalf[i])
        phalf[i-1] = -pfull[i-1] * \
            lambertW_approx(-np.exp(-B/pfull[i-1])/pfull[i-1])

    return phalf


def polar2XYZ(lon, lat, alt, Re=3400*10**3):  # radian
    '''
    Spherical to cartesian coordinates transformation
    Args:
        lon,lat (ND array): longitude and latitude, in [rad]
        alt (ND array): altitude in [m]
    Return:
        X,Y,Z in cartesian coordinates [m]
    ***NOTE***
    This is a classic polar coordinate system with colat = pi/2 -lat,  the colatitude and cos(colat) = sin(lat)
    '''
    lon = np.array(lon)
    lat = np.array(lat)
    alt = np.array(alt)
    R = Re+alt
    X = R*np.cos(lon)*np.cos(lat)
    Y = R*np.sin(lon)*np.cos(lat)
    # added in case broadcasted variables are used (e.g. [1,lat,1], [1,1,lon])
    Z = R*np.sin(lat)*np.ones_like(lon)
    return X, Y, Z


def interp_KDTree(var_IN, lat_IN, lon_IN, lat_OUT, lon_OUT, N_nearest=10):
    '''
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
    '''
    from scipy.spatial import cKDTree  # TODO Import called each time. May be moved out of the routine is scipy is a requirement for the pipeline

    dimsIN = var_IN.shape
    nlon_IN = dimsIN[-1]
    nlat_IN = dimsIN[-2]

    # If var is 2D, extend the dimensions for generality
    if len(dimsIN) == 2:
        var_IN = var_IN.reshape(1, nlat_IN, nlon_IN)

    # If input/output latitudes/longitudes are 1D, extend the dimensions for generality:
    if len(lat_IN.shape) == 1:
        lon_IN, lat_IN = np.meshgrid(lon_IN, lat_IN)  # TODO broadcast instead?
    if len(lat_OUT.shape) == 1:
        lon_OUT, lat_OUT = np.meshgrid(lon_OUT, lat_OUT)

    nlat_OUT = lat_OUT.shape[0]
    nlon_OUT = lon_OUT.shape[1]

    # If lat, lon are 1D, broadcast dimensions:

    # Ndim is the product of all input dimensions but lat & lon
    Ndim = np.int(np.prod(dimsIN[0:-2]))
    dims_IN_reshape = tuple(np.append(Ndim, nlon_IN*nlat_IN))
    dims_OUT_reshape = tuple(np.append(Ndim, nlat_OUT*nlon_OUT))
    # Needed if var is (lat,lon)
    dims_OUT = np.append(dimsIN[0:-2], [nlat_OUT, nlon_OUT]).astype(int)
    var_OUT = np.zeros(dims_OUT_reshape)  # Initialization
    Ndimall = np.arange(0, Ndim)  # all indices (does not change)

    # Compute cartesian coordinate for source and target files  polar2XYZ(lon,lat,lev)
    xs, ys, zs = polar2XYZ(lon_IN*np.pi/180, lat_IN*np.pi/180, 0., Re=1.)
    xt, yt, zt = polar2XYZ(lon_OUT*np.pi/180, lat_OUT*np.pi/180, 0., Re=1.)

    tree = cKDTree(list(zip(xs.flatten(), ys.flatten(), zs.flatten())))
    d, inds = tree.query(
        list(zip(xt.flatten(), yt.flatten(), zt.flatten())), k=N_nearest)
    # Inverse distance
    w = 1.0 / d**2
    # sum the weights  and normalize
    var_OUT = np.sum(w*var_IN.reshape(dims_IN_reshape)
                     [:, inds], axis=2)/np.sum(w, axis=1)
    return var_OUT.reshape(dims_OUT)


def cart_to_azimut_TR(u, v, mode='from'):
    '''
    Convert cartesian coordinates or wind vectors to radian,using azimut angle.

    Args:
        x,y: 1D arrays for the cartesian coordinate
        mode='to' direction towards the vector is pointing, 'from': direction from the vector is coming
    Returns:
        Theta [deg], R the polar coordinates
    '''
    if mode == 'from':
        cst = 180
    if mode == 'to':
        cst = 0.
    return np.mod(np.arctan2(u, v)*180/np.pi+cst, 360), np.sqrt(u**2+v**2)


def sfc_area_deg(lon1, lon2, lat1, lat2, R=3390000.):
    '''
    Return the surface between two set of latitudes/longitudes
    S= Int[R**2 dlon cos(lat) dlat]     _____lat2
    Args:                               \    \
        lon1,lon2: in [degree]           \____\lat1
        lat1,lat2: in [degree]        lon1    lon2
        R: planetary radius in [m]
    *** NOTE***
    Lon and Lat define the corners of the area, not the grid cells' centers

    '''
    lat1 *= np.pi/180
    lat2 *= np.pi/180
    lon1 *= np.pi/180
    lon2 *= np.pi/180
    return (R**2)*np.abs(lon1-lon2)*np.abs(np.sin(lat1)-np.sin(lat2))


def area_meridional_cells_deg(lat_c, dlon, dlat, normalize=False, R=3390000.):
    '''
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
    '''
    # Initialize
    area_tot = 1.
    # Compute total area in a longitude band extending from lat[0]-dlat/2 to lat_c[-1]+dlat/2
    if normalize:
        area_tot = sfc_area_deg(-dlon/2, dlon/2,
                                lat_c[0]-dlat/2, lat_c[-1]+dlat/2, R)
    # Now convert to radians
    lat_c = lat_c*np.pi/180
    dlon *= np.pi/180
    dlat *= np.pi/180
    return 2.*R**2*dlon*np.cos(lat_c)*np.sin(dlat/2.)/area_tot


def area_weights_deg(var_shape, lat_c, axis=-2):
    '''
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
    '''

    # var or lat is a scalar, do nothing
    if len(np.atleast_1d(lat_c)) == 1 or len(np.atleast_1d(var_shape)) == 1:
        return np.ones(var_shape)
    else:
        # Then, lat has at least 2 elements
        dlat = lat_c[1]-lat_c[0]
        # Calculate cell areas. Since it is normalized, we can use dlon= 1 and R=1 without changing the result
        # Note that sum(A)=(A1+A2+...An)=1
        A = area_meridional_cells_deg(lat_c, 1, dlat, normalize=True, R=1)
        # var is a 1D array. of size (lat). Easiest case since (w1+w2+...wn)=sum(A)=1 and N=len(lat)
        if len(var_shape) == 1:
            W = A*len(lat_c)
        else:
            # Generate the appropriate shape for the area A, e.g  (time, lev, lat, lon) > (1, 1, lat, 1)
            # In this case, N=time*lev*lat*lon and  (w1+w2+...wn) =time*lev*lon*sum(A) , therefore N/(w1+w2+...wn)=lat
            reshape_shape = [1 for i in range(0, len(var_shape))]
            reshape_shape[axis] = len(lat_c)
            W = A.reshape(reshape_shape)*len(lat_c)
        return W*np.ones(var_shape)


def areo_avg(VAR, areo, Ls_target, Ls_angle, symmetric=True):
    """
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

    """
    # Take the modulo of solar longitude
    areo = np.mod(areo, 360)

    shape_out = VAR.shape[1:]  # All dimensions but time
    VAR = VAR.reshape((len(areo), np.prod(shape_out)))  # flatten array

    # compute bounds from Ls_target and Ls_angle
    Ls_min = Ls_target-Ls_angle/2.
    Ls_max = Ls_target+Ls_angle/2.

    if (Ls_min < 0.):
        Ls_min += 360.
    if (Ls_max > 360.):
        Ls_max -= 360.

    # Initialize output array
    VAR_avg = np.zeros(np.prod(shape_out))

    # This was removed, for exemple, if 10 degree of data are requested around Ls 0:
    #                       Ls 355 <-- (0.00)--> 5
    #   and the file is       Ls 1 <-- (180)--> 357   the  data selected should be 1>5 and 355 > 357
    '''
    #check is the Ls of interest is within the data provided, raise execption otherwise
    if Ls_target <= areo.min() or Ls_target >=areo.max() :
        raise Exception("Error \nNo data found, requested  data :       Ls %.2f <-- (%.2f)--> %.2f\n However, data in file only ranges      Ls %.2f <-- (%.2f)--> %.2f"%(Ls_min,Ls_target,Ls_max,areo.min(),(areo.min()+areo.max())/2.,areo.max()))
    '''

    if Ls_min < areo.min() or Ls_max > areo.max():
        print("In areo_avg() Warning: \nRequested  data ranging    Ls %.2f <-- (%.2f)--> %.2f" %
              (Ls_min, Ls_target, Ls_max))
        if symmetric:  # Case 1: reduce the window
            if Ls_min < areo.min():
                Ls_min = areo.min()
                Ls_angle = 2*(Ls_target-Ls_min)
                Ls_max = Ls_target+Ls_angle/2.

            if Ls_max > areo.max():
                Ls_max = areo.max()
                Ls_angle = 2*(Ls_max-Ls_target)
                Ls_min = Ls_target-Ls_angle/2.

            print("Reshaping data ranging     Ls %.2f <-- (%.2f)--> %.2f" %
                  (Ls_min, Ls_target, Ls_max))
        else:  # Case 2: Use all data available
            print("I am only using            Ls %.2f <-- (%.2f)--> %.2f \n" %
                  (max(areo.min(), Ls_min), Ls_target, min(areo.max(), Ls_max)))
    count = 0

    for t in range(len(areo)):
        # special case Ls around Ls =0 (wrap around)
        if (Ls_min <= areo[t] <= Ls_max):
            VAR_avg += VAR[t, ...]
            count += 1

    if count > 0:
        VAR_avg /= count
    return VAR_avg.reshape(shape_out)


def mass_stream(v_avg, lat, level, type='pstd', psfc=700, H=8000., factor=1.e-8):
    '''
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
    '''
    g = 3.72  # m/s2
    a = 3400*1000  # m
    nlev = len(level)
    shape_out = v_avg.shape

    # If size is (pstd,lat), turns to (pstd,lat,1) for generality
    if len(shape_out) == 2:
        v_avg = v_avg.reshape(nlev, len(lat), 1)

    # Flatten array
    v_avg = v_avg.reshape((nlev, len(lat), np.prod(v_avg.shape[2:])))
    MSF = np.zeros_like(v_avg)

    # Sum variable, same dimensions as v_avg but for the first dimension
    I = np.zeros(v_avg.shape[2:])

    # Make note of NaN positions and replace by zero for downward integration
    isNan = False
    if np.isnan(v_avg).any():
        isNan = True
        mask = np.isnan(v_avg)
        v_avg[mask] = 0.

    # Missing data may also be masked instead of set to NaN:
    isMasked = False
    if np.ma.is_masked(v_avg):
        isMasked = True
        mask0 = np.ma.getmaskarray(v_avg)
        mask = mask0.copy()  # Make a standalone copy of the mask array
        # Set masked elements to 0. Note that this effectively unmask the array as 0. is a valid entry.
        v_avg[mask0] = 0.

    if type == 'pstd':
        Z = H*np.log(psfc/level)
    else:  # Copy zagl or zstd instead of using a pseudo height
        Z = level.copy()

    for k0 in range(nlev-2, 0, -1):
        I[:] = 0.
        for k in range(nlev-2, k0, -1):
            zn = Z[k]
            znp1 = Z[k+1]
            fn = v_avg[k, :, ...] * np.exp(-zn/H)
            fnp1 = v_avg[k+1, :, ...]*np.exp(-znp1/H)
            I = I+0.5*(znp1-zn)*(fnp1+fn)
        MSF[k0, :, ...] = 2*np.pi*a*psfc / \
            (g*H)*np.cos(np.pi/180*lat).reshape([len(lat), 1])*I*factor

    # Replace NaN where they initially were:
    if isNan:
        MSF[mask] = np.NaN
    if isMasked:
        MSF = np.ma.array(MSF, mask=mask)

    return MSF.reshape(shape_out)


def vw_from_MSF(msf, lat, lev, ztype='pstd', norm=True, psfc=700, H=8000.):
    '''
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
    '''
    g = 3.72  # m/s2

    lat = lat*np.pi/180
    var_shape = msf.shape

    xx = lat.copy()
    zz = lev.copy()

    if ztype == 'pstd':
        zz = H*np.log(psfc/lev)

    if norm:
        xx = (xx-xx.min())/(xx.max()-xx.min())
        zz = (zz-zz.min())/(zz.max()-zz.min())

    # Extend broadcasting dimensions for the latitude, e.g  [1,1,lat] if msf is size (time,lev,lat)
    reshape_shape = [1 for i in range(0, len(var_shape))]
    reshape_shape[-1] = lat.shape[0]
    lat1d = lat.reshape(reshape_shape)

    # Transpose shapes:
    T_array = np.arange(len(msf.shape))

    # one permutation only: lat is passed to the 1st dimension
    T_latIN = np.append(T_array[-1], T_array[0:-1])
    # one permutation only: lat is passed to the 1st dimension
    T_latOUT = np.append(T_array[1:], T_array[0])

    T_levIN = np.append(np.append(T_array[-2], T_array[0:-2]), T_array[-1])
    T_levOUT = np.append(np.append(T_array[1:-1], T_array[0]), T_array[-1])

    V = g/(2*np.pi*np.cos(lat1d)) * \
        dvar_dh(msf.transpose(T_levIN), zz).transpose(T_levOUT)
    W = -g/(2*np.pi*np.cos(lat1d)) * \
        dvar_dh(msf.transpose(T_latIN), xx).transpose(T_latOUT)

    return V, W


def alt_KM(press, scale_height_KM=8., reference_press=610.):
    """
    Gives the approximate altitude in km for a given pressure
    Args:
        press: the pressure in [Pa]
        scale_height_KM: a scale height in [km], (default is 8 km, an isothermal at 155K)
        reference_press: reference surface pressure in [Pa], (default is 610 Pa)
    Returns:
        z_KM: the equivalent altitude for that pressure level in [km]

    ***NOTE***
    Scale height is H=rT/g
    """
    return -scale_height_KM*np.log(press/reference_press)  # p to altitude in km


def press_pa(alt_KM, scale_height_KM=8., reference_press=610.):
    """
    Gives the approximate altitude in km for a given pressure
    Args:
        alt_KM: the altitude in  [km]
        scale_height_KM: a scale height in [km], (default is 8 km, an isothermal at 155K)
        reference_press: reference surface pressure in [Pa], (default is 610 Pa)
    Returns:
         press_pa: the equivalent pressure at that altitude in [Pa]
    ***NOTE***
    Scale height is H=rT/g
    """
    return reference_press*np.exp(-alt_KM/scale_height_KM)  # p to altitude in km


def lon180_to_360(lon):
    lon = np.array(lon)
    """
    Transform a float or an array from the -180/+180 coordinate system to 0-360
    Args:
        lon: a float, 1D or 2D array of longitudes in the 180/+180 coordinate system
    Returns:
        lon: the equivalent longitudes in the 0-360 coordinate system

    """
    if len(np.atleast_1d(lon)) == 1:  # lon180 is a float
        if lon < 0:
            lon += 360
    else:  # lon180 is an array
        lon[lon < 0] += 360
        # reogranize lon by increasing values
        lon = np.append(lon[lon < 180], lon[lon >= 180])
    return lon


def lon360_to_180(lon):
    lon = np.array(lon)
    """
    Transform a float or an array from the 0-360 coordinate system to -180/+180
    Args:
        lon: a float, 1D or 2D array of longitudes in the 0-360 coordinate system
    Returns:
        lon: the equivalent longitudes in the -180/+180 coordinate system

    """
    if len(np.atleast_1d(lon)) == 1:  # lon is a float
        if lon > 180:
            lon -= 360
    else:  # lon is an array
        lon[lon > 180] -= 360
        # reogranize lon by increasing values
        lon = np.append(lon[lon < 0], lon[lon >= 0])
    return lon


def shiftgrid_360_to_180(lon, data):  # longitude is LAST
    '''
    This function shift N dimensional data a 0->360 to a -180/+180 grid.
    Args:
        lon: 1D array of longitude 0->360
        data: ND array with last dimension being the longitude (transpose first if necessary)
    Returns:
        data: shifted data
    Note: Use np.ma.hstack instead of np.hstack to keep the masked array properties
    '''
    lon = np.array(lon)
    lon[lon > 180] -= 360.  # convert to +/- 180
    data = np.concatenate(
        (data[..., lon < 0], data[..., lon >= 0]), axis=-1)  # stack data
    return data


def shiftgrid_180_to_360(lon, data):  # longitude is LAST
    '''
    This function shift N dimensional data a -180/+180 grid to a 0->360
    Args:
        lon: 1D array of longitude -180/+180
        data: ND array with last dimension being the longitude (transpose first if necessary)
    Returns:
        data: shifted data
    '''
    lon = np.array(lon)
    lon[lon < 0] += 360.  # convert to 0-360
    data = np.concatenate(
        (data[..., lon < 180], data[..., lon >= 180]), axis=-1)  # stack data
    return data


def second_hhmmss(seconds, lon_180=0.):
    """
    Given the time in seconds return Local true Solar Time at a certain longitude
    Args:
        seconds: a float, the time in seconds
        lon_180: a float, the longitude in -/+180 coordinate
    Returns:
        hours: float, the local time or  (hours,minutes, seconds)

    """
    hours = seconds // (60*60)
    seconds %= (60*60)
    minutes = seconds // 60
    seconds %= 60
    # Add timezone offset (1hr/15 degree)
    hours = np.mod(hours+lon_180/15., 24)

    return np.int(hours), np.int(minutes), np.int(seconds)


def sol_hhmmss(time_sol, lon_180=0.):
    """
    Given the time in days, return the Local true Solar Time at a certain longitude
    Args:
        time_sol: a float, the time, eg. sols 2350.24
        lon_180: a float, the longitude in a -/+180 coordinate
    Returns:
        hours: float, the local time or  (hours,minutes, seconds)
    """
    return second_hhmmss(time_sol*86400., lon_180)


def UT_LTtxt(UT_sol, lon_180=0., roundmin=None):
    '''
    Returns the time in HH:MM:SS format at a certain longitude.
    Args:
        time_sol: a float, the time, eg. sols 2350.24
        lon_180: a float, the center longitude in  -/+180 coordinate. Increment by 1hr every 15 deg
        roundmin: round to the nearest X minute  Typical values are roundmin=1,15,60
    ***Note***
    If roundmin is requested, seconds are not shown
    '''
    def round2num(number, interval):
        # Internal function to round a number to the closest  range.
        # e.g. round2num(26,5)=25 ,round2num(28,5)=30
        return round(number / interval) * interval

    hh, mm, ss = sol_hhmmss(UT_sol, lon_180)

    if roundmin:
        sec = hh*3600+mm*60+ss
        # Round to the nearest increment (in seconds) and run a second pass
        rounded_sec = round2num(sec, roundmin*60)
        hh, mm, ss = second_hhmmss(rounded_sec, lon_180)
        return '%02d:%02d' % (hh, mm)
    else:
        return '%02d:%02d:%02d' % (hh, mm, ss)


def dvar_dh(arr, h=None):
    '''
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

    '''
    h = np.array(h)
    d_arr = np.zeros_like(arr)
    if h.any():
        # h is provided as a 1D array
        if len(h.shape) == 1:
            reshape_shape = np.append(
                [arr.shape[0]-2], [1 for i in range(0, arr.ndim - 1)])
            d_arr[0, ...] = (arr[1, ...]-arr[0, ...])/(h[1]-h[0])
            d_arr[-1, ...] = (arr[-1, ...]-arr[-2, ...])/(h[-1]-h[-2])
            d_arr[1:-1, ...] = (arr[2:, ...]-arr[0:-2, ...]) / \
                (np.reshape(h[2:]-h[0:-2], reshape_shape))
        # h has the same dimension as var
        elif h.shape == arr.shape:
            d_arr[0, ...] = (arr[1, ...]-arr[0, ...])/(h[1, ...]-h[0, ...])
            d_arr[-1, ...] = (arr[-1, ...]-arr[-2, ...]) / \
                (h[-1, ...]-h[-2, ...])
            d_arr[1:-1, ...] = (arr[2:, ...]-arr[0:-2, ...]
                                )/(h[2:, ...]-h[0:-2, ...])
        else:
            print('Error,h.shape=', h.shape, 'arr.shape=', arr.shape)
    # h is not defined, we return only d_var, not d_var/dh
    else:
        d_arr[0, ...] = arr[1, ...]-arr[0, ...]
        d_arr[-1, ...] = arr[-1, ...]-arr[-2, ...]
        # > Note the 0.5 factor since differentiation uses a central scheme
        d_arr[1:-1, ...] = 0.5*(arr[2:, ...]-arr[0:-2, ...])

    return d_arr


def zonal_detrend(VAR):
    '''
    Substract zonnally averaged mean value from a field
    Args:
        VAR: ND-array with detrending dimension last (e.g time,lev,lat,lon)
    Returns:
        OUT: detrented field (same size as input)

    ***NOTE***
    RuntimeWarnings are expected if the slice contains only NaN, which is the case below the surface
    and above the model's top in the interpolated files. We will disable those warnings temporarily
    '''
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return VAR-np.nanmean(VAR, axis=-1)[..., np.newaxis]


def get_trend_2D(VAR, LON, LAT, type_trend='wmean'):
    '''
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
    '''
    var_shape = np.array(VAR.shape)

    # Type 'zonal' is the easiest as averaging is performed over 1 dimension only.
    if type_trend == 'zonal':
        return np.repeat(np.nanmean(VAR, axis=-1)[..., np.newaxis], var_shape[-1], axis=-1)

    # The following options involve avering over both latitude and longitude dimensions:

    # Flatten array e.g. turn (10,36,lat,lon) to (360,lat,lon)
    nflatten = int(np.prod(var_shape[:-2]))
    reshape_flat = np.append(nflatten, var_shape[-2:])
    VAR = VAR.reshape(reshape_flat)

    TREND = np.zeros(reshape_flat)
    for ii in range(nflatten):
        if type_trend == 'mean':
            TREND[ii, ...] = np.mean(VAR[ii, ...].flatten())
        elif type_trend == 'wmean':
            W = area_weights_deg(var_shape[-2:], LAT[:, 0])
            TREND[ii, ...] = np.mean((VAR[ii, ...]*W).flatten())
        elif type_trend == '2D':
            TREND[ii, ...] = regression_2D(LON, LAT, VAR[ii, :, :], order=1)
        else:
            print("Error, in area_trend, type '%s' not recognized" % (type_trend))
            return None
    return TREND.reshape(var_shape)


def regression_2D(X, Y, VAR, order=1):
    '''
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
    '''
    if order == 1:
        A = np.array([X.flatten(), Y.flatten(), np.ones_like(X.flatten())]).T

        # An Equivalent notation is:
        # A=np.c_[X.flatten(),Y.flatten(),np.ones_like(X.flatten())]

        b = VAR.flatten()

        # P is the solution of  A X =b, ==> P[0] x + P[1]y + P[2] = z
        P, residuals, rank, s = np.linalg.lstsq(A, b, rcond=None)

        Z = P[0]*X + P[1]*Y+np.ones_like(X)*P[2]

    elif order == 2:
        # best-fit quadratic curve: a X**2 + 2b X*Y +c Y**2 +2dX +2eY+f
        XX = X.flatten()
        YY = Y.flatten()
        ZZ = VAR.flatten()
        data = np.zeros((len(XX), 3))
        data[:, 0] = XX
        data[:, 1] = YY
        data[:, 2] = ZZ

        A = np.c_[np.ones(data.shape[0]), data[:, :2], np.prod(
            data[:, :2], axis=1), data[:, :2]**2]
        P, _, _, _ = np.linalg.lstsq(A, data[:, 2], rcond=None)

        # evaluate it on a grid (using vector product)
        Z = np.dot(np.c_[np.ones(XX.shape), XX, YY, XX *
                   YY, XX**2, YY**2], P).reshape(X.shape)
    return Z


def daily_to_average(varIN, dt_in, nday=5, trim=True):
    '''
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
    '''
    vshape_in = varIN.shape
    Nin = vshape_in[0]  # time dimension

    iperday = int(np.round(1/dt_in))
    combinedN = int(iperday*nday)
    N_even = Nin//combinedN
    N_left = Nin % combinedN

    # Nin/(ndayxiperday) is not a round number
    if N_left != 0 and not trim:
        # Do the average on the even part
        vreshape = np.append([-1, combinedN], vshape_in[1:]).astype(int)
        var_even = np.mean(
            varIN[0:N_even*combinedN, ...].reshape(vreshape), axis=1)

        # Left over time steps
        var_left = np.mean(
            varIN[N_even*combinedN:, ...], axis=0, keepdims=True)
        # Combine both
        varOUT = np.concatenate((var_even, var_left), axis=0)

    # Nin/(ndayxiperday) is a round number or we request to trim the array
    else:
        vreshape = np.append([-1, combinedN], vshape_in[1:]).astype(int)
        varOUT = np.mean(
            varIN[0:N_even*combinedN, ...].reshape(vreshape), axis=1)
    return varOUT


def daily_to_diurn(varIN, time_in):
    '''
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
    '''
    dt_in = time_in[1]-time_in[0]
    iperday = int(np.round(1/dt_in))
    vshape_in = varIN.shape
    vreshape = np.append([-1, iperday], vshape_in[1:]).astype(int)
    varOUT = varIN.reshape(vreshape)

    # Get time of day in hours
    tod = np.mod(time_in[0:iperday]*24, 24)

    # Sort by time of day, e.g. if  tod=[6.,12.,18.,0.] re-arange  into [0.,6.,12.,18.]
    # every element is array  is greater than the one on its left.
    if not np.all(tod[1:] >= tod[:-1], axis=0):

        # This returns the permutation, e.g. if tod=[6.,12.,18.,0.], i_sort = [3, 0, 1, 2]
        i_sort = sorted(range(len(tod)), key=lambda k: tod[k])
        tod = tod[i_sort]
        varOUT = varOUT[:, i_sort, ...]
    return varOUT
# =========================================================================
# =======================vertical grid utilities===========================
# =========================================================================


def gauss_profile(x, alpha, x0=0.):
    """ Return Gaussian line shape at x This can be used to generate a bell-shaped mountain"""
    return np.sqrt(np.log(2) / np.pi) / alpha\
        * np.exp(-((x-x0) / alpha)**2 * np.log(2))


def compute_uneven_sigma(num_levels, N_scale_heights, surf_res, exponent, zero_top):
    """
    Construct an initial array of sigma based on the number of levels, an exponent
    Args:
        num_levels: the number of levels
        N_scale_heights: the number of scale heights to the top of the model (e.g scale_heights =12.5 ~102km assuming 8km scale height)
        surf_res: the resolution at the surface
        exponent: an exponent to increase th thickness of the levels
        zero_top: if True, force the top pressure boundary (in N=0) to 0 Pa
    Returns:
        b: an array of sigma layers

    """
    b = np.zeros(int(num_levels)+1)
    for k in range(0, num_levels):
        zeta = 1.-k/np.float(num_levels)  # zeta decreases with k
        z = surf_res*zeta + (1.0 - surf_res)*(zeta**exponent)
        b[k] = np.exp(-z*N_scale_heights)  # z goes from 1 to 0
    b[-1] = 1.0
    if(zero_top):
        b[0] = 0.0
    return b


def transition(pfull, p_sigma=0.1, p_press=0.05):
    """
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
    """
    t = np.zeros_like(pfull)
    for k in range(0, len(pfull)):
        if(pfull[k] <= p_press):
            t[k] = 0.0
        elif (pfull[k] >= p_sigma):
            t[k] = 1.0
        else:
            x = pfull[k] - p_press
            xx = p_sigma - p_press
            t[k] = (np.sin(0.5*np.pi*x/xx))**2

    return t


def swinbank(plev, psfc, ptrans=1.):
    """
    Compute ak and bk values with a transition based on Swinbank
    Args:
        plev: the pressure levels in Pa
        psfc: the surface pressure in Pa
        ptrans:the transition pressure in Pa
    Returns:
         aknew, bknew,ks: the coefficients for the new layers
    """

    # ks= number of pure pressure levels
    ktrans = np.argmin(np.abs(plev - ptrans))
    km = len(plev)-1

    aknew = np.zeros(len(plev))
    bknew = np.zeros(len(plev))

    #   pnorm= 1.e5;
    pnorm = psfc
    eta = plev / pnorm

    ep = eta[ktrans+1]  # ks= number of pure pressure levels
    es = eta[-1]
    rnorm = 1. / (es-ep)**2

    #   Compute alpha, beta, and gamma using Swinbank's formula
    alpha = (ep**2 - 2.*ep*es) / (es-ep)**2
    beta = 2.*ep*es**2 / (es-ep)**2
    gamma = -(ep*es)**2 / (es-ep)**2

    #   Pure Pressure levels
    aknew = eta * pnorm

    #  Hybrid pressure-sigma levels
    kdex = range(ktrans+1, km)
    aknew[kdex] = alpha*eta[kdex] + beta + gamma/eta[kdex]
    aknew[kdex] = aknew[kdex] * pnorm
    aknew[-1] = 0.0

    bknew[kdex] = (plev[kdex] - aknew[kdex])/psfc
    bknew[-1] = 1.0

    # find the transition level ks where (bk[ks]>0)
    ks = 0
    while bknew[ks] == 0.:
        ks += 1
    # ks is the one that would be use in fortran indexing in fv_eta.f90
    return aknew, bknew, ks


def polar_warming(T, lat, outside_range=np.NaN):
    """
    Return the polar warming, following  [McDunn et al. 2013]: Characterization of middle-atmosphere polar warming at Mars, JGR
    A. Kling
    Args:
        T:   temperature array, 1D, 2D or ND, with the latitude dimension FIRST (transpose as needed)
        lat: latitude array
        outside_range: values to set the polar warming to outside the range. Default is Nan but 'zero' may be desirable.
    Returns:
        DT_PW:   The polar warming in [K]


    *NOTE*  polar_warming() concatenates the results from both hemispheres obtained from the nested function PW_half_hemisphere()
    """

    def PW_half_hemisphere(T_half, lat_half, outside_range=np.NaN):

        # Easy case, T is a 1D on the latitude direction only
        if len(T_half.shape) == 1:
            imin = np.argmin(T_half)
            imax = np.argmax(T_half)

            # Note that we compute polar warming at ALL latitudes and then set NaN the latitudes outside the desired range.
            # We test on the absolute values (np.abs) of the latitudes, therefore the function is usable on both hemispheres
            DT_PW_half = T_half-T_half[imin]
            exclude = np.append(np.where(np.abs(
                lat_half)-np.abs(lat_half[imin]) < 0), np.where(np.abs(lat_half[imax])-np.abs(lat_half) < 0))
            DT_PW_half[exclude] = outside_range  # set to NaN
            return DT_PW_half

        # General case for N dimensions
        else:

            # Flatten the diemnsions but the first dimension (the latitudes)
            arr_flat = T_half.reshape(
                [T_half.shape[0], np.prod(T_half.shape[1:])])
            LAT_HALF = np.repeat(
                lat_half[:, np.newaxis], arr_flat.shape[1], axis=1)

            imin = np.argmin(arr_flat, axis=0)
            imax = np.argmax(arr_flat, axis=0)

            # Initialize four working arrays
            tmin0, tmax0, latmin0, latmax0 = [
                np.zeros_like(arr_flat) for _ in range(4)]

            # get the min/max temperature and latitudes
            for i in range(0, arr_flat.shape[1]):
                tmax0[:, i] = arr_flat[imax[i], i]
                tmin0[:, i] = arr_flat[imin[i], i]
                latmin0[:, i] = lat_half[imin[i]]
                latmax0[:, i] = lat_half[imax[i]]

            # Compute polar warming for that hemisphere
            DT_PW_half = arr_flat-tmin0

            # Set to NaN values outside the range.
            tuple_lower_than_latmin = np.where(
                np.abs(LAT_HALF)-np.abs(latmin0) < 0)
            tuple_larger_than_latmax = np.where(
                np.abs(latmax0)-np.abs(LAT_HALF) < 0)

            DT_PW_half[tuple_lower_than_latmin] = outside_range
            DT_PW_half[tuple_larger_than_latmax] = outside_range

            return DT_PW_half.reshape(T_half.shape)

    # ======================================================
    # ======Actual calculations for both hemispheres========
    # ======================================================
    T_SH = T[0:len(lat)//2]
    lat_SH = lat[0:len(lat)//2]

    T_NH = T[len(lat)//2:]
    lat_NH = lat[len(lat)//2:]

    return np.concatenate((PW_half_hemisphere(T_SH, lat_SH, outside_range), PW_half_hemisphere(T_NH, lat_NH, outside_range)), axis=0)


def tshift(array, lon, timeo, timex=None):
    '''
    Conversion to uniform local time.
    Args:
        array: variable to be shifted. Assume longitude is the first dimension and time in the last dimension
        lon: longitude
        timeo : time_of_day index from input file
        timex (optional) : local time (hr) to shift to, e.g. '3. 15.'
    Returns:
        tshift: array shifted to uniform local time.

    ***Note***
    If timex is not specified, the file is interpolated on the same time_of_day as the input
    '''
    if np.shape(array) == len(array):
        print('Need longitude and time dimensions')
        return

    dims = np.shape(array)  # get dimensions of array
    end = len(dims)-1
    id = dims[0]  # number of longitudes in file
    nsteps = len(timeo)   # number of timesteps per day in input

    nsf = np.float_(nsteps)   # number of timesteps per day in input

    timeo = np.squeeze(timeo)

    # array dimensions for output
    if timex is None:       # time shift all local times
        nsteps_out = nsteps
    else:
        nsteps_out = len(timex)

    # Assuming time is last dimension, check if it is local time timex
    # If not, reshape the array into (stuff, days, local time)
    if dims[end] != nsteps:
        ndays = dims[end] / nsteps
        if ndays*nsteps != dims[end]:
            print('Time dimensions do not conform')
            return
        array = np.reshape(array, (dims[0, end-1], nsteps, ndays))
        newdims = np.linspace(len(dims+1), dtype=np.int32)
        newdims[len(dims)-1] = len(dims)
        newdims[len(dims)] = len(dims)-1
        array = np.transpose(array, newdims)

    dims = np.shape(array)  # get new dims of array if reshaped

    if len(dims) > 2:
        recl = np.prod(dims[1:len(dims)-1])
    else:
        recl = 1

    array = np.reshape(array, (id, recl, nsteps))

    # create output array
    narray = np.zeros((id, recl, nsteps_out))

    dt_samp = 24.0/nsteps  # Time increment of input data (in hours)

    # time increment of output
    if timex is None:       # match dimensions of output file to input
        dt_save = dt_samp  # Time increment of output data (in hours)
    else:
        dt_save = 1.    # assume output time increament in 1 hour

    # calculate interpolation indeces
    # convert east longitude to equivalent hours
    xshif = 24.0*lon/360.
    kk = np.where(xshif < 0)
    xshif[kk] = xshif[kk]+24.

    fraction = np.zeros((id, nsteps_out))
    imm = np.zeros((id, nsteps_out))
    ipp = np.zeros((id, nsteps_out))

    for nd in range(nsteps_out):
        #dtt = nd*dt_save - xshif - timex[0] + dt_samp
        if timex is None:
            dtt = nd*dt_save-xshif - timeo[0] + dt_samp
        else:
            # time_out - xfshif - tod[0] + hrs/stpe in input
            dtt = timex[nd] - xshif

        #      insure that data local time is bounded by [0,24] hours
        kk = np.where(dtt < 0.)
        dtt[kk] = dtt[kk] + 24.

        im = np.floor(dtt/dt_samp)  # this is index into the data aray
        fraction[:, nd] = dtt-im*dt_samp
        kk = np.where(im < 0.)
        im[kk] = im[kk] + nsf

        ipa = im + 1.
        kk = np.where(ipa >= nsf)
        ipa[kk] = ipa[kk] - nsf

        imm[:, nd] = im[:]
        ipp[:, nd] = ipa[:]

    fraction = fraction / dt_samp  # assume uniform tinc between input data samples

    #           Now carry out the interpolation
    for nd in range(nsteps_out):  # Number of output time levels
        for i in range(id):  # Number of longitudes
            im = np.int(imm[i, nd]) % 24
            ipa = np.int(ipp[i, nd])
            frac = fraction[i, nd]
            narray[i, :, nd] = (1.-frac)*array[i, :, im] + \
                frac*array[i, :, ipa]

    narray = np.squeeze(narray)
    ndimsfinal = np.zeros(len(dims), dtype=int)
    for nd in range(end):
        ndimsfinal[nd] = dims[nd]
    ndimsfinal[end] = nsteps_out
    narray = np.reshape(narray, ndimsfinal)

    return narray


def lin_interp(X_in, X_ref, Y_ref):
    '''
    Simple linear interpolation with no dependance on scipy
    Args:
        X_in (float or array): input values
        X_ref (array): x values
        Y_ref (array): y values
    Returns:
        Y_out: y value linearly interpolated at X_in
    '''
    X_ref = np.array(X_ref)
    Y_ref = np.array(Y_ref)
    # ===Definition of the interpolating function=====

    def lin_oneElement(x, X_ref, Y_ref):
        if x < X_ref.min() or x > X_ref.max():
            return np.NaN
        # Find closest left-hand size index
        n = np.argmin(np.abs(x-X_ref))
        if X_ref[n] > x or n == len(X_ref):
            n -= 1
        a = (Y_ref[n+1]-Y_ref[n])/(X_ref[n+1]-X_ref[n])
        b = Y_ref[n]-a*X_ref[n]
        return a*x+b

    # ======Wrapper to the function above=====
    if len(np.atleast_1d(X_in)) == 1:
        Y_out = lin_oneElement(X_in, X_ref, Y_ref)
    else:
        X_in = np.array(X_in)
        Y_out = np.zeros_like(X_in)
        for i, x_in in enumerate(X_in):
            Y_out[i] = lin_oneElement(x_in, X_ref, Y_ref)
    return Y_out


def add_cyclic(data, lon):
    """
    Add an additional cyclic (overlapping) point to a 2D array, useful for azimuth and orthographic projections
    Args:
        data: 2D array of size (nlat,nlon)
        lon: 1D array of longitudes
    Returns:
        data_c: 2D array of size (nlat,nlon+1), with last column identical to the 1st
        lon_c: 1D array of longitudes size nlon+1 where the last element is lon[-1]+dlon

    """
    # Compute increment
    dlon = lon[1]-lon[0]
    # Create new array, size [nlon+1]
    data_c = np.zeros((data.shape[0], data.shape[1]+1), np.float)
    data_c[:, 0:-1] = data[:, :]
    data_c[:, -1] = data[:, 0]
    return data_c, np.append(lon, lon[-1]+dlon)


def spherical_div(U, V, lon_deg, lat_deg, R=3400*1000., spacing='varying'):
    '''
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

    '''
    lon = lon_deg*np.pi/180
    lat = lat_deg*np.pi/180
    var_shape = U.shape

    # Transpose shapes:
    T_array = np.arange(len(U.shape))
    # one permutation only: lon is passsed to the 1st dimension
    T_lonIN = np.append(T_array[-1], T_array[0:-1])
    # one permutation only: lon is passsed to the 1st dimension
    T_lonOUT = np.append(T_array[1:], T_array[0])
    T_latIN = np.append(np.append(T_array[-2], T_array[0:-2]), T_array[-1])
    T_latOUT = np.append(np.append(T_array[1:-1], T_array[0]), T_array[-1])

    # ----lon, lat are 1D arrays---
    if len(lon.shape) == 1:
        # Extend broadcasting dimensions for the latitude, e.g  [1,1,lat,1] if U is size (time,lev,lat,lon)
        reshape_shape = [1 for i in range(0, len(var_shape))]
        reshape_shape[-2] = lat.shape[0]
        lat_b = lat.reshape(reshape_shape)
        if spacing == 'regular':
            out = 1/(R*np.cos(lat_b))*(np.gradient(U, axis=-1) /
                                       (lon[1]-lon[0])+np.gradient(V*np.cos(lat_b), axis=-2)/(lat[1]-lat[0]))
        else:
            out = 1/(R*np.cos(lat_b))*(dvar_dh(U.transpose(T_lonIN), lon).transpose(T_lonOUT) +
                                       dvar_dh((V*np.cos(lat_b)).transpose(T_latIN), lat).transpose(T_latOUT))
    # ----lon, lat are 2D array---
    else:
        # if U is (time,lev,lat,lon), also reshape lat ,lon to (time,lev,lat,lon)
        if var_shape != lon.shape:
            # (time,lev,lat,lon)> (time,lev) and reverse, so first lev, then time
            for ni in var_shape[:-2][::-1]:
                lat = np.repeat(lat[np.newaxis, ...], ni, axis=0)
                lon = np.repeat(lon[np.newaxis, ...], ni, axis=0)

        out = 1/(R*np.cos(lat))*(dvar_dh(U.transpose(T_lonIN), lon.transpose(T_lonIN)).transpose(T_lonOUT) +
                                 dvar_dh((V*np.cos(lat)).transpose(T_latIN), lat.transpose(T_latIN)).transpose(T_latOUT))
    return out


def spherical_curl(U, V, lon_deg, lat_deg, R=3400*1000., spacing='varying'):
    '''
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

    '''
    lon = lon_deg*np.pi/180
    lat = lat_deg*np.pi/180

    var_shape = U.shape

    # Transpose shapes:
    T_array = np.arange(len(U.shape))
    # one permutation only: lon is passsed to the 1st dimension
    T_lonIN = np.append(T_array[-1], T_array[0:-1])
    # one permutation only: lon is passsed to the 1st dimension
    T_lonOUT = np.append(T_array[1:], T_array[0])
    T_latIN = np.append(np.append(T_array[-2], T_array[0:-2]), T_array[-1])
    T_latOUT = np.append(np.append(T_array[1:-1], T_array[0]), T_array[-1])

    # ----lon, lat are 1D arrays---
    if len(lon.shape) == 1:
        # Extend broadcasting dimensions for the latitude, e.g  [1,1,lat,1] if U is size (time,lev,lat,lon)
        reshape_shape = [1 for i in range(0, len(var_shape))]
        reshape_shape[-2] = lat.shape[0]
        lat_b = lat.reshape(reshape_shape)
        if spacing == 'regular':
            out = 1/(R*np.cos(lat_b))*(np.gradient(V, axis=-1) /
                                       (lon[1]-lon[0])-np.gradient(U*np.cos(lat_b), axis=-2)/(lat[1]-lat[0]))
        else:
            out = 1/(R*np.cos(lat_b))*(dvar_dh(V.transpose(T_lonIN), lon).transpose(T_lonOUT) -
                                       dvar_dh((U*np.cos(lat_b)).transpose(T_latIN), lat).transpose(T_latOUT))

    # ----lon, lat are 2D array---
    else:
        # if U is (time,lev,lat,lon), also reshape lat ,lon to (time,lev,lat,lon)
        if var_shape != lon.shape:
            # (time,lev,lat,lon)> (time,lev) and reverse, so first lev, then time
            for ni in var_shape[:-2][::-1]:
                lat = np.repeat(lat[np.newaxis, ...], ni, axis=0)
                lon = np.repeat(lon[np.newaxis, ...], ni, axis=0)

        out = 1/(R*np.cos(lat))*(dvar_dh(V.transpose(T_lonIN), lon.transpose(T_lonIN)).transpose(T_lonOUT) -
                                 dvar_dh((U*np.cos(lat)).transpose(T_latIN), lat.transpose(T_latIN)).transpose(T_latOUT))
    return out


def frontogenesis(U, V, theta, lon_deg, lat_deg, R=3400*1000., spacing='varying'):
    '''
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

    '''
    lon = lon_deg*np.pi/180
    lat = lat_deg*np.pi/180

    var_shape = U.shape

    # Transpose shapes:
    T_array = np.arange(len(U.shape))
    # one permutation only: lon is passsed to the 1st dimension
    T_lonIN = np.append(T_array[-1], T_array[0:-1])
    # one permutation only: lon is passsed to the 1st dimension
    T_lonOUT = np.append(T_array[1:], T_array[0])
    T_latIN = np.append(np.append(T_array[-2], T_array[0:-2]), T_array[-1])
    T_latOUT = np.append(np.append(T_array[1:-1], T_array[0]), T_array[-1])

    # ----lon, lat are 1D arrays---
    if len(lon.shape) == 1:
        # Extend broadcasting dimensions for the colatitude, e.g  [1,1,lat,1] if U is size (time,lev,lat,lon)
        reshape_shape = [1 for i in range(0, len(var_shape))]
        reshape_shape[-2] = lat.shape[0]
        lat_b = lat.reshape(reshape_shape)
        if spacing == 'regular':

            du_dlon = np.gradient(U, axis=-1)/(lon[1]-lon[0])
            dv_dlon = np.gradient(V, axis=-1)/(lon[1]-lon[0])
            dtheta_dlon = np.gradient(theta, axis=-1)/(lon[1]-lon[0])

            du_dlat = np.gradient(U, axis=-2)/(lat[1]-lat[0])
            dv_dlat = np.gradient(V, axis=-2)/(lat[1]-lat[0])
            dtheta_dlat = np.gradient(theta, axis=-2)/(lat[1]-lat[0])

        else:

            du_dlon = dvar_dh(U.transpose(T_lonIN), lon).transpose(T_lonOUT)
            dv_dlon = dvar_dh(V.transpose(T_lonIN), lon).transpose(T_lonOUT)
            dtheta_dlon = dvar_dh(theta.transpose(
                T_lonIN), lon).transpose(T_lonOUT)

            du_dlat = dvar_dh(U.transpose(T_latIN), lat).transpose(T_latOUT)
            dv_dlat = dvar_dh(V.transpose(T_latIN), lat).transpose(T_latOUT)
            dtheta_dlat = dvar_dh(theta.transpose(
                T_latIN), lat).transpose(T_latOUT)

    # ----lon, lat are 2D array---
    else:
        # if U is (time,lev,lat,lon), also reshape lat ,lon to (time,lev,lat,lon)
        if var_shape != lon.shape:
            lat_b = lat.copy()
            # (time,lev,lat,lon)> (time,lev) and reverse, so first lev, then time
            for ni in var_shape[:-2][::-1]:
                lat = np.repeat(lat[np.newaxis, ...], ni, axis=0)
                lon = np.repeat(lon[np.newaxis, ...], ni, axis=0)

            du_dlon = dvar_dh(U.transpose(T_lonIN), lon.transpose(
                T_lonIN)).transpose(T_lonOUT)
            dv_dlon = dvar_dh(V.transpose(T_lonIN), lon.transpose(
                T_lonIN)).transpose(T_lonOUT)
            dtheta_dlon = dvar_dh(theta.transpose(
                T_lonIN), lon.transpose(T_lonIN)).transpose(T_lonOUT)

            du_dlat = dvar_dh(U.transpose(T_latIN), lat.transpose(
                T_latIN)).transpose(T_latOUT)
            dv_dlat = dvar_dh(V.transpose(T_latIN), lat.transpose(
                T_latIN)).transpose(T_latOUT)
            dtheta_dlat = dvar_dh(theta.transpose(
                T_latIN), lat.transpose(T_latIN)).transpose(T_latOUT)

    out = -(1/(R*np.cos(lat_b))*dtheta_dlon)**2 *\
        (1/(R*np.cos(lat_b))*du_dlon - V*np.tan(lat_b)/R) -\
        (1/R*dtheta_dlat)**2*(1/R*dv_dlat) -\
        (1/(R*np.cos(lat_b))*dtheta_dlon)*(1/R*dtheta_dlat) *\
        (1/(R*np.cos(lat_b))*dv_dlon+1/R*du_dlat+U*np.tan(lat_b)/R)

    return out


def MGSzmax_ls_lat(ls, lat):
    '''
    Return the max altitude for the dust from "MGS scenario"
    from Montmessin et al. (2004), Origin and role of water ice clouds in the Martian
                                   water cycle as inferred from a general circulation model

    Args:
        ls  : solar longitude in degree
        lat : latitude in degree
    Returns:
        zmax : top altitude for the dusk in [km]
    '''
    lat = np.array(lat)*np.pi/180
    ls_p = (np.array(ls)-158)*np.pi/180

    return 60+18*np.sin(ls_p)-(32+18*np.sin(ls_p))*np.sin(lat)**4-8*np.sin(ls_p)*np.sin(lat)**5


def MGStau_ls_lat(ls, lat):
    '''
    Return the max altitude for the dust from "MGS scenario"
    from Montmessin et al. (2004), Origin and role of water ice clouds in the Martian
                                   water cycle as inferred from a general circulation model

    Args:
        ls  : solar longitude in degree
        lat : latitude in degree
    Returns:
        zmax : top altitude for the dusk in [km]
    '''
    lat = np.array(lat)
    ls_p = (np.array(ls)-250)*np.pi/180

    tn = 0.1
    teq = 0.2+0.3*np.cos(0.5*ls_p)**14
    ts = 0.1+0.4*np.cos(0.5*ls_p)**14

    # We have tanh(-x)=-tanh(x)
    t_north = tn+0.5*(teq-tn)*(1+np.tanh(4.5-lat/10))
    t_south = ts+0.5*(teq-ts)*(1+np.tanh(4.5+lat/10))

    # One latitude
    if len(np.atleast_1d(lat)) == 1:
        tau = t_north if lat >= 0 else t_south
    else:
        tau = np.zeros_like(lat)
        tau[lat <= 0] = t_south[lat <= 0]
        tau[lat > 0] = t_north[lat > 0]

    return tau


def broadcast(var_1D, shape_out, axis):
    '''
    Broadcast a 1D array based on a variable's dimensions
    Args:
        var_1D     (1D array), e.g. lat size (36), or time size (133)
        shape_out (ND list) : braodcasting shape e.g temp.shape= [133,(lev),36,(lon)]
    Return:
        var_b (ND array): broadcasted variables, e.g. size [1,36,1,1] for lat or [133,1,1,1] for time
    '''
    var_1D = np.atleast_1d(
        var_1D)  # Special case where var_1D has only one element
    reshape_shape = [1 for i in range(0, len(shape_out))]
    reshape_shape[axis] = len(var_1D)  # e.g [28,1,1,1]
    return var_1D.reshape(reshape_shape)


def ref_atmosphere_Mars_PTD(Zi):
    '''
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
    '''
    # =================================
    # =======Internal functions========
    # =================================

    def alt_to_temp_quad(Z, Z0, T0, gam, a):
        '''
        Return the a representative globally and annually average temperature in the form T(z)= a(z-z0)**2 +b(z-z0) +T0
        Args:
            Z:  Altitude in [m]
            Z0: Starting altitude [m]
            T0,gam,a: quadratic coefficients.
        Returns:
            T: temperature at altitude Z [K]
        '''
        return T0+gam*(Z-Z0)+a*(Z-Z0)**2

    def alt_to_press_quad(Z, Z0, P0, T0, gam, a, rgas, g):
        '''
        Return the pressure in Pa in the troposphere as a function of height for a quadratic temperature profile.
        T(z)= T0+ gam(z-z0)+a*(z-z0)**2
        Args:
            Z:  Altitude in [m]
            Z0: Referecence altitude in [m]
            P0: reference pressure  at Z0[Pa]
            T0: reference temperature  at Z0[Pa]
            gam,a: linear and quadratic coefficients for the temperature
            rgas,g:specific gas constant [J/kg/K] and gravity [m/s2]
        Returns:
            P: pressure at alitude Z [Pa]
        '''
        delta = 4*a*T0-gam**2

        if delta >= 0:
            sq_delta = np.sqrt(4*a*T0-gam**2)
            return P0*np.exp(-2*g/(rgas*sq_delta)*(np.arctan((2*a*(Z-Z0)+gam)/sq_delta)-np.arctan(gam/sq_delta)))
        else:
            delta = -delta
            sq_delta = np.sqrt(delta)
            return P0*(((2*a*(Z-Z0)+gam)-sq_delta)*(gam+sq_delta)/(((2*a*(Z-Z0)+gam)+sq_delta)*(gam-sq_delta)))**(-g/(rgas*sq_delta))

    def P_mars_120_300(Zi, Z0=120000., P0=0.00012043158397922564, p1=1.09019694e-04, p2=-3.37385416e-10):
        '''
        Martian pressure from 120 to 300km . Above ~120km P=P0 exp(-(z-z0)g/rT) is not a good approximation as the fluid is in molecular
        regime
        Fit from a diurnal average of the MCD database at lat 0, Ls 150. from Alex K.
        Consistently to what is done for the Earth, we use P=P0 exp(-az-bz**2-cz**c-d**4 ...)
        Args:
            Z: altitude in [m]
        Returns:
            P: pressure in [Pa]
        '''
        return P0*np.exp(-p1*(Zi-Z0)-p2*(Zi-Z0)**2)
    # The following is a best fit of globally averaged temperature profile from various sources: Legacy GCM, MCS, MCD

    def T_analytic_scalar(Zi):
        if Zi <= 57000:
            return alt_to_temp_quad(Zi, Z0=0, T0=225.9, gam=-0.00213479, a=1.44823e-08)
        elif 57000 < Zi <= 110000:
            return alt_to_temp_quad(Zi, Z0=57000, T0=151.2, gam=-0.000367444, a=-6.8256e-09)
        elif 110000 < Zi <= 170000:
            return alt_to_temp_quad(Zi, Z0=110000, T0=112.6, gam=0.00212537, a=-1.81922e-08)
        elif 170000 <= Zi:
            return 174.6
    # Analytical solution for the pressure derived from the temperature profile

    def P_analytic_scalar(Zi):
        if Zi <= 57000:
            return alt_to_press_quad(Zi, Z0=0, P0=610, T0=225.9, gam=-0.00213479, a=1.44823e-08, rgas=192, g=3.72)
        elif 57000 < Zi <= 110000:
            return alt_to_press_quad(Zi, Z0=57000, P0=1.2415639872674782, T0=151.2, gam=-0.000367444, a=-6.8256e-09, rgas=192, g=3.72)
        # The following must be discarded above 120 km when we enter the molecular regime
        elif 110000 < Zi <= 120000:
            return alt_to_press_quad(Zi, Z0=110000, P0=0.0005866878792825923, T0=112.6, gam=0.00212537, a=-1.81922e-08, rgas=192, g=3.72)
        elif 120000 <= Zi:
            return P_mars_120_300(Zi)

    # Vectorize array as needed:
    if len(np.atleast_1d(Zi)) > 1:
        P_analytic_scalar = np.vectorize(P_analytic_scalar)
        T_analytic_scalar = np.vectorize(T_analytic_scalar)
    return P_analytic_scalar(Zi), T_analytic_scalar(Zi),  P_analytic_scalar(Zi)/(192*T_analytic_scalar(Zi))


def press_to_alt_atmosphere_Mars(Pi):
    '''
    Return the altitude in m as a function of pressure from the analytical calculations derived above.
    Args:
        Pi (float or 1D array): input pressure in Pa (must be <=610 Pa)
    Return:
        Z (float ot 1D array): corresponding altitude in m
    '''
    # =================================
    # =======Internal functions========
    # =================================

    def press_to_alt_quad(P, Z0, P0, T0, gam, a, rgas, g):
        '''
        Return the altitude in m as a function of pressure for a quadratic temperature profile.
        T(z)= T0+ gam(z-z0)+a*(z-z0)**2
        Args:
            P:  Pressure in [Pa]
            Z0: Referecence altitude in [m]
            P0: reference pressure  at Z0[Pa]
            T0: reference temperature  at Z0[Pa]
            gam,a: linear and quadratic coefficients for the temperature
            rgas,g:specific gas constant [J/kg/K] and gravity [m/s2]
        Returns:
            Z: altitude in m
        '''
        delta = 4*a*T0-gam**2

        if delta >= 0:
            sq_delta = np.sqrt(4*a*T0-gam**2)
            return Z0+sq_delta/(2*a)*np.tan(np.arctan(gam/sq_delta)+np.log(P0/P)*rgas*sq_delta/(2*g))-gam/(2*a)
        else:
            delta = -delta
            sq_delta = np.sqrt(delta)
            return Z0+(gam-sq_delta)/(2*a)*(1-(P/P0)**(-rgas*sq_delta/g))/((gam-sq_delta)/(gam+sq_delta)*(P/P0)**(-rgas*sq_delta/g)-1)

    def press_to_alt_mars_120_300(P, Z0=120000., P0=0.00012043158397922564, p1=1.09019694e-04, p2=-3.37385416e-10):
        '''
        Martian altitude as a function of pressure from 120 to 300km .
        Args:
            P: pressure in [m]
        Returns:
            Z: altitude in [m]
        '''
        delta = p1**2-4*p2*np.log(P/P0)  # delta >0 on this pressure interval
        return (-p1+np.sqrt(delta))/(2*p2) + Z0

    def alt_analytic_scalar(Pi):
        '''
        Analytical solution for altitude as a function of pressure.
        '''
        if Pi >= 610:
            return 0.
        elif 610 > Pi >= 1.2415639872674782:  # This is the pressure from alt_to_press_quad at 57000. m
            return press_to_alt_quad(Pi, Z0=0, P0=610, T0=225.9, gam=-0.00213479, a=1.44823e-08, rgas=192, g=3.72)
        elif 1.2415639872674782 > Pi >= 0.0005866878792825923:  # 57000 to 110000m
            return press_to_alt_quad(Pi, Z0=57000, P0=1.2415639872674782, T0=151.2, gam=-0.000367444, a=-6.8256e-09, rgas=192, g=3.72)
        elif 0.0005866878792825923 > Pi >= 0.00012043158397922564:  # 110000m to 120000 m
            return press_to_alt_quad(Pi, Z0=110000, P0=0.0005866878792825923, T0=112.6, gam=0.00212537, a=-1.81922e-08, rgas=192, g=3.72)
        elif 0.00012043158397922564 > Pi:  # 120000m to 300000
            return press_to_alt_mars_120_300(Pi, Z0=120000., P0=0.00012043158397922564, p1=1.09019694e-04, p2=-3.37385416e-10)

    if len(np.atleast_1d(Pi)) > 1:
        alt_analytic_scalar = np.vectorize(alt_analytic_scalar)
    return alt_analytic_scalar(Pi)


# ==================================Projections==================================
'''
The projections below were implemented by Alex Kling, following:
An Album of Map Projections,
USGS  Professional Paper 1453, (1994)
>> https://pubs.usgs.gov/pp/1453/report.pdf
'''
# ===============================================================================


def azimuth2cart(LAT, LON, lat0, lon0=0):
    '''
    Azimuthal equidistant projection, convert from lat/lon to cartesian coordinates
    Args:
        LAT,LON: 1D or 2D array of latitudes, longitudes in degree, size [nlat,nlon]
        lat0,lon0:(floats) coordinates of the pole
    Returns:
        X,Y: cartesian coordinates for the latitudes and longitudes
    '''

    LAT = LAT*np.pi/180
    lat0 = lat0*np.pi/180
    LON = LON*np.pi/180
    lon0 = lon0*np.pi/180

    c = np.arccos(np.sin(lat0) * np.sin(LAT) + np.cos(lat0)
                  * np.cos(LAT) * np.cos(LON-lon0))
    k = c / np.sin(c)
    X = k * np.cos(LAT) * np.sin(LON-lon0)
    Y = k * (np.cos(lat0)*np.sin(LAT) - np.sin(lat0)
             * np.cos(LAT)*np.cos(LON-lon0))
    return X, Y


def ortho2cart(LAT, LON, lat0, lon0=0):
    '''
    Orthographic projection, convert from lat/lon to cartesian coordinates
    Args:
        LAT,LON: 1D or 2D array of latitudes, longitudes in degree, size [nlat,nlon]
        lat0,lon0:(floats) coordinates of the pole
    Returns:
        X,Y: cartesian coordinates for the latitudes and longitudes
        MASK: NaN array that is used to hide the back side of the planet
    '''

    LAT = LAT*np.pi/180
    lat0 = lat0*np.pi/180
    LON = LON*np.pi/180
    lon0 = lon0*np.pi/180
    MASK = np.ones_like(LON)

    X = np.cos(LAT) * np.sin(LON-lon0)
    Y = np.cos(lat0)*np.sin(LAT) - np.sin(lat0)*np.cos(LAT)*np.cos(LON-lon0)

    # Filter values on the other side of the globe, i.e cos(c)<0
    cosc = np.sin(lat0) * np.sin(LAT) + np.cos(lat0) * \
        np.cos(LAT) * np.cos(LON-lon0)
    MASK[cosc < 0] = np.NaN
    return X, Y, MASK


def mollweide2cart(LAT, LON):
    '''
    Mollweide projection, convert from lat/lon to cartesian coordinates
    Args:
        LAT,LON: 1D or 2D array of latitudes, longitudes in degree, size [nlat,nlon]
    Returns:
        X,Y: cartesian coordinates for the latitudes and longitudes
    '''

    LAT = np.array(LAT)*np.pi/180
    LON = np.array(LON)*np.pi/180
    lon0 = 0

    def compute_theta(lat):
        '''
        Internal function to compute theta, lat is in radians here
        '''
        theta0 = lat
        sum = 0
        running = True
        # Solve for theta using Newton–Raphson
        while running and sum <= 100:
            theta1 = theta0-(2*theta0+np.sin(2*theta0)-np.pi *
                             np.sin(lat))/(2+2*np.cos(2*theta0))
            sum += 1
            if np.abs((theta1-theta0)) < 10**-3:
                running = False
            theta0 = theta1
        if sum == 100:
            print("Warning,in mollweide2cart():  Reached Max iterations")
        return theta1

    # Float or 1D array
    if len(np.atleast_1d(LAT).shape) == 1:
        nlat = len(np.atleast_1d(LAT))
        LAT = LAT.reshape((nlat))
        LON = LON.reshape((nlat))
        THETA = np.zeros((nlat))
        for i in range(0, nlat):
            THETA[i] = compute_theta(LAT[i])

    else:  # 2D array
        nlat = LAT.shape[0]
        nlon = LAT.shape[1]
        theta = np.zeros((nlat))
        for i in range(0, nlat):
            theta[i] = compute_theta(LAT[i, 0])
        THETA = np.repeat(theta[:, np.newaxis], nlon, axis=1)

    X = 2*np.sqrt(2)/np.pi*(LON-lon0)*np.cos(THETA)
    Y = np.sqrt(2)*np.sin(THETA)
    return np.squeeze(X), np.squeeze(Y)


def robin2cart(LAT, LON):
    '''
    Robinson projection, convert from lat/lon to cartesian coordinates
    Args:
        LAT,LON: floats, 1D or 2D array (nalt,nlon) of latitudes, longitudes in degree
    Returns:
        X,Y: cartesian coordinates for the latitudes and longitudes
    '''
    lon0 = 0.

    LAT = np.array(LAT)*np.pi/180
    LON = np.array(LON)*np.pi/180

    lat_ref = np.array([0, 5, 10, 15, 20, 25, 30, 35, 40,
                       45, 50, 55, 60, 65, 70, 75, 80, 85, 90.])*np.pi/180
    x_ref = np.array([1.0000, 0.9986, 0.9954, 0.9900, 0.9822, 0.9730, 0.9600, 0.9427, 0.9216,
                     0.8962, 0.8679, 0.8350, 0.7986, 0.7597, 0.7186, 0.6732, 0.6213, 0.5722, 0.5322])
    y_ref = np.array([0.0000, 0.0620, 0.1240, 0.1860, 0.2480, 0.3100, 0.3720, 0.4340, 0.4958,
                     0.5571, 0.6176, 0.6769, 0.7346, 0.7903, 0.8435, 0.8936, 0.9394, 0.9761, 1.0000])

    # Float or 1D array
    if len(np.atleast_1d(LAT).shape) == 1:
        X1 = lin_interp(np.abs(LAT), lat_ref, x_ref)
        Y1 = np.sign(LAT)*lin_interp(np.abs(LAT), lat_ref, y_ref)
    else:
        # 2D array
        nlat = LAT.shape[0]
        nlon = LAT.shape[1]
        lat = LAT[:, 0]
        x1 = lin_interp(np.abs(lat), lat_ref, x_ref)
        y1 = np.sign(lat)*lin_interp(np.abs(lat), lat_ref, y_ref)

        X1 = np.repeat(x1[:, np.newaxis], nlon, axis=1)
        Y1 = np.repeat(y1[:, np.newaxis], nlon, axis=1)

    X = 0.8487*X1*(LON-lon0)
    Y = 1.3523*Y1

    return X, Y

# ===================== (End projections section) ================================


def sol2ls(jld, cummulative=False):
    """
    Return the solar longitude Ls as a function of the sol number. Sol 0 is spring equinox.
    Args:
        jld [float or 1D array]: sol number after perihelion
        cummulative [bool]     : if True, result is cummulative Ls 0>360>720 etc..

    Returns:
        Ls: The corresponding solar longitude Ls
    """

    # constants
    year = 668.0  # This is in martian days  ie sols
    zero_date = 488.  # time of perihelion passage
    equinox = 180  # day of northern equinox:  (for 668 sol year)
    small_value = 1.0e-7
    pi = np.math.pi
    degrad = pi/180.0

    # if jld is a scalar, reshape as a 1-element arra
    jld = np.array(jld).astype(float).reshape(len(np.atleast_1d(jld)))

    # =============================================================================
    # ===Internal function 1 : calculate Ls with 0>360 using a numerical solver====
    # =============================================================================
    def sol2ls_mod(jld):
        '''
        Based on Tanguys' aerols.py
        Also useful link: http://www.jgiesen.de/kepler/kepler.html
        '''
        # Specify orbit eccentricity and angle of planets's inclination

        ec = .093  # orbit eccentricity
        er = ((1.0+ec)/(1.0-ec))**0.5

        # Initialize working arrays
        w, rad, als, areols, MY = [np.zeros_like(jld) for _ in range(0, 5)]

        # date= days since last perihelion passage
        date = jld - zero_date

        # determine true anomaly at equinox:  eq

        qq = 2.0 * pi * equinox / year   # qq is the mean anomaly
        e = 1.0
        diff = 1.0
        while (diff > small_value):
            ep = e - (e-ec*np.math.sin(e)-qq) / (1.0-ec*np.cos(e))
            diff = abs(ep-e)
            e = ep

        eq1 = 2.0 * np.math.atan(er * np.math.tan(0.5*e))

        # determine true anomaly at current date:  w

        for i in range(0, len(jld)):
            e = 1.0
            diff = 1.0
            em = 2. * pi * date[i] / year
            while (diff > small_value):
                ep = e - (e - ec * np.sin(e) - em) / (1.0 - ec * np.cos(e))
                diff = abs(ep-e)
                e = ep
            w[i] = 2.0 * np.math.atan(er * np.math.tan(0.5*e))

        als = w - eq1  # Aerocentric Longitude
        areols = als/degrad
        areols[areols < 0.] += 360.
        return areols

    # =============================================================================
    # ===Internal function 2 : calculate cumulative Ls with 0>720 cumulative=======
    # =============================================================================

    def sol2ls_cum(jld):
        '''
        Calculate cumulative Ls.
        Contineous solar longitude Ls_c =0-359- 361.. 720  are obtained by adding +360 to the Ls at every equinox based on the sol number.
        Since sol2ls function uses a numerical solver, the equinox may return either 359.9999 or 0.0001. +360 should only be added in the later case to mitigate outlier points.

        For those edges cases where Ls is close to 359.9, the routine calculate again the Ls at a later time (say 1 sols) to check for outlier points.
        '''
        # Calculate cummulative Ls using sol2ls function() and adding +360 for every mars year
        areols, MY = [np.zeros_like(jld) for _ in range(0, 2)]
        date = jld - zero_date
        MY = (date-equinox)//(year)+1  # MY=(date-equinox)//(year)
        Ls_mod = sol2ls_mod(jld)
        Ls_cum = Ls_mod+MY*360.

        # Check indexes where the returned Ls is close to 360.
        # The [0] turns tuple from np.where into a list
        index = np.where(Ls_mod >= 359.9)[0]

        for ii in index:
            jld_plus1 = jld[ii] + \
                1.  # compute Ls one day after (arbitrary length)
            Ls_plus1 = sol2ls(jld_plus1)
            date_plus1 = jld_plus1 - zero_date
            MY_plus1 = (date_plus1-equinox)//(668.)+1  # Compute MY
            Ls_cum_plus1 = Ls_plus1+MY_plus1 * \
                360.  # Cummulative Ls 1 day after.
            # If things are smooth the Ls should go from [359>361]. If it reads [359>721], we need to update the MY for those indices
            # difference between two consecutive Ls, should be small unless Ls_cum was too big at the first place
            diff = Ls_cum_plus1-Ls_cum[ii]
            if diff < 0:
                MY[ii] -= 1
        # Recompute one more more time with updated MY
        Ls_cum = Ls_mod+MY*360.
        return Ls_cum

    if cummulative:
        return sol2ls_cum(jld)
    else:
        return sol2ls_mod(jld)
# ====End sol 2 Ls function====
