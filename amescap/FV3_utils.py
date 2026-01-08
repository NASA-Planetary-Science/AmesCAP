# !/usr/bin/env python3
"""
FV3_utils contains internal Functions for processing data in MGCM
output files such as vertical interpolation.

These functions can be used on their own outside of CAP if they are
imported as a module::

    from /u/path/FV3_utils import fms_press_calc

Third-party Requirements:

    * ``numpy``
    * ``warnings``
    * ``scipy``
"""

# Load generic Python modules
import warnings     # suppress errors triggered by NaNs
import numpy as np
from scipy import optimize
from scipy.spatial import cKDTree

# p_half = half-level = layer interfaces
# p_full = full-level = layer midpoints


def fms_press_calc(psfc, ak, bk, lev_type='full'):
    """
    Returns the 3D pressure field from the surface pressure and the
    ak/bk coefficients.

    :param psfc: the surface pressure [Pa] or an array of surface
        pressures (1D, 2D, or 3D if time dimension)
    :type  psfc: array
    :param ak: 1st vertical coordinate parameter
    :type  ak: array
    :param bk: 2nd vertical coordinate parameter
    :type  bk: array:
    :param lev_type: "full" (layer midpoints) or "half"
        (layer interfaces). Defaults to "full."
    :type  lev_type: str
    :return: the 3D pressure field at the full levels
        ``PRESS_f(Nk-1:,:,:)`` or half-levels ``PRESS_h(Nk,:,:,)`` [Pa]

    Calculation::

        --- 0 --- TOP        ========  p_half
        --- 1 ---
                             --------  p_full

                             ========  p_half
        ---Nk-1---           --------  p_full
        --- Nk --- SFC       ========  p_half
                            / / / / /

    .. note::
        Some literature uses pk (pressure) instead of ak with
        ``p3d = ps * bk + P_ref * ak`` instead of ``p3d = ps * bk + ak``
    """

    Nk = len(ak)
    # If ``psfc`` is a float (e.g., ``psfc = 700.``), make it a
    # 1-element array (e.g., ``psfc = [700]``)
    if len(np.atleast_1d(psfc)) == 1:
        psfc = np.array([np.squeeze(psfc)])

    # Flatten ``psfc`` array to generalize it to N dimensions
    psfc_flat = psfc.flatten()

    # Expand the dimensions. Vectorized calculations:
    # (Np) -> (Np, Nk)
    psfc_v = np.repeat(psfc_flat[:, np.newaxis], Nk, axis = 1)
    # (Nk) -> (Np, Nk)
    ak_v = np.repeat(ak[np.newaxis, :], len(psfc_flat), axis = 0)
    # (Nk) -> (1,  Nk)
    bk_v = np.repeat(bk[np.newaxis, :], 1, axis = 0)
    # (Nk) -> (1,  Nk)

    # Pressure at layer interfaces. The size of Z axis is ``Nk``
    PRESS_h = psfc_v*bk_v + ak_v

    # Pressure at layer midpoints. The size of Z axis is ``Nk-1``
    PRESS_f = np.zeros((len(psfc_flat), Nk-1))

    if ak[0] == 0 and bk[0] == 0:
        # Top layer (1st element is ``i = 0`` in Python)
        PRESS_f[:, 0] = 0.5 * (PRESS_h[:, 0]+PRESS_h[:, 1])
    else:
        PRESS_f[:, 0] = (
            (PRESS_h[:, 1] - PRESS_h[:, 0])
            / np.log(PRESS_h[:, 1] / PRESS_h[:, 0])
            )

    # The rest of the column (``i = 1 ... Nk``). [2:] goes from the 3rd
    # element to ``Nk`` and [1:-1] goes from the 2nd element to ``Nk-1``
    PRESS_f[:, 1:] = ((PRESS_h[:, 2:] - PRESS_h[:, 1:-1])
                      / np.log(PRESS_h[:, 2:] / PRESS_h[:, 1:-1]))

    # First, transpose ``PRESS(:, Nk)`` to ``PRESS(Nk, :)``. Then
    # reshape ``PRESS(Nk, :)`` to the original pressure shape
    # ``PRESS(Nk, :, :, :)`` (resp. ``Nk-1``)
    #TODO
    if lev_type == "full":
        new_dim_f = np.append(Nk-1, psfc.shape)
        return PRESS_f.T.reshape(new_dim_f)
    elif lev_type == "half":
        new_dim_h = np.append(Nk, psfc.shape)
        return PRESS_h.T.reshape(new_dim_h)
    else:
        raise Exception("Pressure level type not recognized by "
                        "``press_lev()``: use 'full' or 'half' ")


def fms_Z_calc(psfc, ak, bk, T, topo=0., lev_type="full"):
    """
    Returns the 3D altitude field [m] AGL (or above aeroid).

    :param psfc: The surface pressure [Pa] or array of surface
        pressures (1D, 2D, or 3D)
    :type  psfc: array
    :param ak: 1st vertical coordinate parameter
    :type  ak: array
    :param bk: 2nd vertical coordinate parameter
    :type  bk: array
    :param T: The air temperature profile. 1D array (for a single grid
        point), ND array with VERTICAL AXIS FIRST
    :type  T: 1D array or ND array
    :param topo: The surface elevation. Same dimension as ``psfc``.
        If None is provided, AGL is returned
    :type  topo: array
    :param lev_type: "full" (layer midpoint) or "half" (layer
        interfaces). Defaults to "full"
    :type  lev_type: str
    :return: The layer altitude at the full level ``Z_f(:, :, Nk-1)``
        or half-level ``Z_h(:, :, Nk)`` [m]. ``Z_f`` and ``Z_h`` are
        AGL if ``topo = None``. ``Z_f`` and ``Z_h`` are above aeroid
        if topography is not None.

    Calculation::

        --- 0 --- TOP        ========  z_half
        --- 1 ---
                            --------  z_full

                            ========  z_half
        ---Nk-1---          --------  z_full
        --- Nk --- SFC      ========  z_half
                            / / / / /

    .. note::
        Expands to the time dimension using::

            topo = np.repeat(zsurf[np.newaxis, :], ps.shape[0], axis = 0)

    Calculation is derived from
        ``./atmos_cubed_sphere_mars/Mars_phys.F90``::

            # (dp/dz = -rho g) => (dz = dp/(-rho g)) and
            # (rho = p/(r T)) => (dz = rT/g * (-dp/p))

            # Define log-pressure (``u``) as:
            u = ln(p)

            # Then:
            du = {du/dp}*dp = {1/p)*dp} = dp/p

            # Finally, ``dz`` for the half-layers:
            (dz = rT/g * -(du)) => (dz = rT/g * (+dp/p))
            # with ``N`` layers defined from top to bottom.

    Z_half calculation::

        # Hydrostatic relation within the layer > (P(k+1)/P(k) =
        # exp(-DZ(k)/H))

        # layer thickness:
        DZ(k) = rT/g * -(du)

        # previous layer altitude + thickness of layer:
        Z_h k) = Z_h(k+1)  +DZ_h(h)

    Z_full calculation::

        # previous altitude + half the thickness of previous layer and
        # half of current layer
        Z_f(k) = Z_f(k+1) + (0.5 DZ(k) + 0.5 DZ(k+1))

        # Add ``+0.5 DZ(k)-0.5 DZ(k)=0`` and re-organiz the equation
        Z_f(k) = Z_f(k+1) + DZ(k) + 0.5 (DZ(k+1) - DZ(k))
        Z_f(k) = Z_h(k+1) + 0.5 (DZ(k+1) - DZ(k))

    The specific heat ratio:
    ``γ = cp/cv (cv = cp-R)`` => ``γ = cp/(cp-R)`` Also ``(γ-1)/γ=R/cp``

    The dry adiabatic lapse rate:
    ``Γ = g/cp`` => ``Γ = (gγ)/R``

    The isentropic relation:
    ``T2 = T1(p2/p1)**(R/cp)``

    Therefore::

        line 1) =====Thalf=====zhalf[k]  \
        line 2)                           \
        line 3)                            \
        line 4) -----Tfull-----zfull[k]     \ T(z)= To-Γ (z-zo)
        line 5)                              \
        line 6)                               \
        line 7) =====Thalf=====zhalf[k+1]      \

    Line 1: T_half[k+1]/Tfull[k] = (p_half[k+1]/p_full[k])**(R/Cp)

    Line 4: From the lapse rate, assume T decreases linearly within the
    layer so ``T_half[k+1] = T_full[k] + Γ(Z_full[k]-Z_half[k+1])``
    and (``Tfull < Thalf`` and ``Γ > 0``)

    Line 7: ``Z_full[k] = Z_half[k] + (T_half[k+1]-T_full[k])/Γ``
    Pulling out ``Tfull`` from above equation and using ``Γ = (gγ)/R``::

        Z_full[k] = (Z_half[k+1] + (R Tfull[k]) / (gγ)(T_half[k+1]
        / T_full[k] - 1))

    Using the isentropic relation above::

        Z_full = (Z_half[k+1] + (R Tfull[k]) / (gγ)(p_half[k+1]
        / p_full[k])**(R/Cp)-1))
    """

    g = 3.72 # acc. m/s2
    r_co2 = 191.00 # kg/mol
    Nk = len(ak)

    # Get the half and full pressure levels from ``fms_press_calc``
    # Z axis is first
    PRESS_f = fms_press_calc(psfc, ak, bk, 'full')
    PRESS_h = fms_press_calc(psfc, ak, bk, 'half')

    # If ``psfc`` is a float, turn it into a 1-element array:
    if len(np.atleast_1d(psfc)) == 1:
        psfc = np.array([np.squeeze(psfc)])
    if len(np.atleast_1d(topo)) == 1:
        topo = np.array([np.squeeze(topo)])

    psfc_flat = psfc.flatten()
    topo_flat = topo.flatten()

    # Reshape arrays for vector calculations and compute log pressure
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
    # Isothermal within the layer, we have ``Z = Z0 + r*T0/g*ln(P0/P)``
    for k in range(Nk-2, -1, -1):
        Z_h[k, :] = (Z_h[k+1, :]
                     + (r_co2*T[k, :]/g)
                     * (logPPRESS_h[k+1, :]-logPPRESS_h[k, :]))
        Z_f[k, :] = (Z_h[k+1, :]
                     + (r_co2*T[k, :]/g)
                     * (1-PRESS_h[k, :]/PRESS_f[k, :]))

    # Return the arrays
    if lev_type == "full":
        new_dim_f = np.append(Nk-1, psfc.shape)
        return Z_f.reshape(new_dim_f)
    elif lev_type == "half":
        new_dim_h = np.append(Nk, psfc.shape)
        return Z_h.reshape(new_dim_h)
    else: # Return the levels in Z coordinates [m]
        raise Exception("Altitude level type not recognized: use 'full' "
                        "or 'half'")


# TODO: delete: Former version of find_n(): only provides 1D > 1D and
# ND > 1D mapping

def find_n0(Lfull_IN, Llev_OUT, reverse_input=False):
    """
    Return the index for the level(s) just below ``Llev_OUT``.
    This assumes ``Lfull_IN`` is increasing in the array
    (e.g., ``p(0) = 0``, ``p(N) = 1000`` [Pa]).

    :param Lfull_IN: Input pressure [Pa] or altitude [m] at layer
        midpoints. ``Level`` dimension is FIRST
    :type  Lfull_IN: array
    :param Llev_OUT: Desired level type for interpolation [Pa] or [m]
    :type  Llev_OUT: float or 1D array
    :param reverse_input: Reverse array (e.g., if ``z(0) = 120 km``,
        ``z(N) = 0km`` -- which is typical -- or if input data is
        ``p(0) = 1000Pa``, ``p(N) = 0Pa``)
    :type  reverse_input: bool
    :return: ``n`` index for the level(s) where the pressure is just
        below ``plev``

    .. note::
        If ``Lfull_IN`` is a 1D array and ``Llev_OUT`` is a float
        then ``n`` is a float.

    .. note::
        If ``Lfull_IN`` is ND ``[lev, time, lat, lon]`` and
        ``Llev_OUT`` is a 1D array of size ``klev`` then ``n`` is an
        array of size ``[klev, Ndim]`` with ``Ndim = [time, lat, lon]``
    """

    # Number of original layers
    Lfull_IN = np.array(Lfull_IN)
    Nlev = len(np.atleast_1d(Llev_OUT))

    if Nlev == 1:
        Llev_OUT = np.array([Llev_OUT])

    # Get input variable dimensions
    dimsIN = Lfull_IN.shape
    Nfull = dimsIN[0]
    dimsOUT = tuple(np.append(Nlev, dimsIN[1:]))
    # ``Ndim`` is the product of all the dimensions except for the
    # vertical axis
    Ndim = int(np.prod(dimsIN[1:]))
    Lfull_IN = np.reshape(Lfull_IN, (Nfull, Ndim))

    if reverse_input:
        Lfull_IN = Lfull_IN[::-1, :]

    ncol = Lfull_IN.shape[-1]
    n = np.zeros((Nlev, ncol), dtype=int)

    for i in range(0, Nlev):
        for j in range(0, ncol):
            n[i, j] = np.argmin(np.abs(Lfull_IN[:, j] - Llev_OUT[i]))
            if Lfull_IN[n[i, j], j] > Llev_OUT[i]:
                n[i, j] = n[i, j] - 1
    return n


def find_n(X_IN, X_OUT, reverse_input=False, modulo=None):
    """
    Maps the closest index from a 1D input array to a ND output array
    just below the input values.

    :param X_IN: Source level [Pa] or [m]
    :type  X_IN: float or 1D array
    :param X_OUT: Desired pressure [Pa] or altitude [m] at layer
        midpoints. Level dimension is FIRST
    :type  X_OUT: array
    :param reverse_input: If input array is decreasing (e.g., if z(0)
        = 120 km, z(N) = 0 km, which is typical, or if data is
        p(0) = 1000 Pa, p(N) = 0 Pa, which is uncommon)
    :type  reverse_input: bool
    :return: The index for the level(s) where the pressure < ``plev``
    """

    if type(X_IN) != np.ndarray:
        # Number of original layers
        X_IN = np.array(X_IN)

    if len(np.atleast_1d(X_OUT)) == 1:
        # If one value is requested, convert float to array
        X_OUT = np.array([X_OUT])
    elif type(X_OUT) != np.ndarray:
        # Convert list to numpy array as needed
        X_OUT = np.array(X_OUT)

    # Get input variable dimensions
    dimsIN = X_IN.shape
    # Get output variable dimensions
    dimsOUT = X_OUT.shape
    # Get size of interpolation axis
    N_IN = dimsIN[0]
    N_OUT = dimsOUT[0]

    # Get number of elements in arrays other than interpolation axis
    NdimsIN = int(np.prod(dimsIN[1:]))
    NdimsOUT = int(np.prod(dimsOUT[1:]))

    if ((NdimsIN > 1 and NdimsOUT > 1) and
        (NdimsIN != NdimsOUT)):
        print("*** Error in ``find_n()``: dimensions of arrays other "
              "than the interpolated (first) axis must be 1 or identical ***")

    # Ndims_IN and Ndims_OUT are either 1 or identical
    Ndim = max(NdimsIN, NdimsOUT)
    X_IN = np.reshape(X_IN, (N_IN, NdimsIN))
    X_OUT = np.reshape(X_OUT, (N_OUT, NdimsOUT))

    # Reverse input array if monotically decreasing
    if reverse_input:
        X_IN = X_IN[::-1, :]

    n = np.zeros((N_OUT, Ndim), dtype = int)

    # Some redundancy below but this allows keeping the "if" statement
    # out of the larger loop over all of the array elements
    if len(dimsIN) == 1:
        for i in range(N_OUT):
            for j in range(Ndim):
                # Handle the case where j might be out of bounds
                if j < NdimsOUT:
                    n[i, j] = np.argmin(np.abs(X_OUT[i, j] - X_IN[:]))
                    if X_IN[n[i, j]] > X_OUT[i, j]:
                        n[i, j] = n[i, j] - 1
                else:
                    # For indices beyond the available dimensions, use index 0
                    n[i, j] = 0
    elif len(dimsOUT) == 1:
        for i in range(N_OUT):
            for j in range(Ndim):
                # Handle the case where j might be out of bounds for X_IN
                if j < NdimsIN:
                    n[i, j] = np.argmin(np.abs(X_OUT[i] - X_IN[:, j]))
                    if X_IN[n[i, j], j] > X_OUT[i]:
                        n[i, j] = n[i, j] - 1
                else:
                    # For indices beyond the available dimensions, use index 0
                    n[i, j] = 0
    else:
        for i in range(N_OUT):
            for j in range(Ndim):
                # Handle the case where j might be out of bounds for either array
                if j < NdimsIN and j < NdimsOUT:
                    n[i, j] = np.argmin(np.abs(X_OUT[i, j] - X_IN[:, j]))
                    if X_IN[n[i, j], j] > X_OUT[i, j]:
                        n[i, j] = n[i, j] - 1
                else:
                    # For indices beyond the available dimensions, use index 0
                    n[i, j] = 0

    if len(dimsOUT) == 1:
        n = np.squeeze(n)
    return n


def expand_index(Nindex, VAR_shape_axis_FIRST, axis_list):
    """
    Repeat interpolation indices along an axis.

    :param Nindex: Interpolation indices, size is (``n_axis``,
        ``Nfull = [time, lat, lon]``)
    :type  Nindex: idx
    :param VAR_shape_axis_FIRST: Shape for the variable to interpolate
        with interpolation axis first (e.g., ``[tod, time, lev, lat, lon]``)
    :type  VAR_shape_axis_FIRST: tuple
    :param axis_list: Position or list of positions for axis to insert
        (e.g., ``2`` for ``lev`` in ``[tod, time, lev, lat, lon]``, ``[2, 4]``
        for ``lev`` and ``lon``). The axis positions are those for the final
        shape (``VAR_shape_axis_FIRST``) and must be INCREASING
    :type  axis_list: int or list
    :return: ``LFULL`` a 2D array (size ``n_axis``,
        ``NfFULL = [time, lev, lat, lon]``) with the indices expanded
        along the ``lev`` dimension and flattened

    .. note::
        Example of application:
        Observational time of day may be the same at all vertical levels
        so the interpolation of a 5D variable ``[tod, time, lev, lat, lon]``
        only requires the interpolation indices for ``[tod, time, lat, lon]``.
        This routine expands the indices from ``[tod, time, lat, lon]`` to
        ``[tod, time, lev, lat, lon]`` with ``Nfull = [time, lev, lat, lon]``
        for use in interpolation.
    """

    # If one element, turn axis to list
    if len(np.atleast_1d(axis_list)) == 1:
        axis_list = [axis_list]

    # Size for the interpolation (e.g., [tod, time, lev, lat, lon]
    Nfull = Nindex.shape[0]
    # Desired output size with ``lev`` axis repeated
    # [tod, time, lev, lat, lon]
    dimsOUT_flat = tuple(np.append(Nfull, np.prod(VAR_shape_axis_FIRST[1:])))

    # Reconstruct the initial (un-flattened) size of ``Nindex``
    dimsIN = []
    for ii, len_axis in enumerate(VAR_shape_axis_FIRST):
        # Use the dimenions from VAR_shape_axis_FIRST
        if (ii > 0 and not
            ii in axis_list):
            # Skip the first (interpolated) axis and add to the list of
            # initial dims unless the axis is the one we are expending.
            dimsIN.append(len_axis)

    # Initial shape for ``Nindex``
    dimsIN = np.insert(dimsIN, 0, Nfull)
    # Reshape ``Nindex`` from its iniatial flattened sahpe
    # ``[Nfull, time, lat, lon]`` -> ND array ``[tod, time, lat, lon]``
    Nindex = np.reshape(Nindex, dimsIN)
    # Repeat the interpolation indices on the requested axis
    for ii, len_axis in enumerate(VAR_shape_axis_FIRST):
        if ii in axis_list:
            # e.g., ``Nindex`` is now ``[tod, time, lev, lat, lon]``
            Nindex = np.repeat(np.expand_dims(Nindex, axis = ii), len_axis, ii)
    # Return the new, flattened version of ``Nindex``
    return Nindex.reshape(dimsOUT_flat)


def vinterp(varIN, Lfull, Llev, type_int="log", reverse_input=False,
            masktop=True, index=None):
    """
    Vertical linear or logarithmic interpolation for pressure or altitude.

    :param varIN: Variable to interpolate (VERTICAL AXIS FIRST)
    :type  varIN: ND array
    :param Lfull: Pressure [Pa] or altitude [m] at full layers, same
        dimensions as ``varIN``
    :type  Lfull: array
    :param Llev: Desired level for interpolation [Pa] or [m]. May be
        increasing or decreasing as the output levels are processed one
        at a time
    :type  Llev: 1D array
    :param type_int: "log" for logarithmic (typically pressure),
        "lin" for linear (typically altitude)
    :type  type_int: str
    :param reverse_input: Reverse input arrays. e.g., if
        ``zfull[0]`` = 120 km then ``zfull[N]`` = 0km (typical) or if
        input data is ``pfull[0]``=1000 Pa, ``pfull[N]``=0 Pa
    :type  reverse_input: bool
    :param masktop: Set to NaN values if above the model top
    :type  masktop: bool
    :param index: Indices for the interpolation, already processed as
        ``[klev, Ndim]``. Indices calculated if not provided
    :type  index: None or array
    :return: ``varOUT`` variable interpolated on the ``Llev`` pressure
        or altitude levels

    .. note::
        This interpolation assumes pressure decreases with height::

            --  0  -- TOP  [0 Pa]   : [120 km]| X_OUT = Xn*A + (1-A)*Xn + 1
            --  1  --               :         |
                                    :         |
            --  n  -- pn   [30 Pa]  : [800 m] | Xn
                                    :         |
            --  k  -- Llev [100 Pa] : [500 m] | X_OUT
            -- n+1 -- pn+1 [200 Pa] : [200 m] | Xn+1

            -- SFC --
            / / / / / /

        with ``A = log(Llev/pn + 1) / log(pn/pn + 1)`` in "log" mode
        or ``A = (zlev-zn + 1) / (zn-zn + 1)`` in "lin" mode
    """

    Nlev = len(np.atleast_1d(Llev))
    if Nlev == 1:
        # Special case where only 1 layer is requested
        Llev = np.array([Llev])

    # Get input variable dimensions
    dimsIN = varIN.shape
    Nfull = dimsIN[0]

    if len(varIN.shape) == 1:
        # Special case where ``varIN`` is a single profile
        varIN = varIN.reshape([Nfull, 1])
    if len(Lfull.shape) == 1:
        # Special case where ``Lfull`` is a single profile
        Lfull = Lfull.reshape([Nfull, 1])

    # Repeat in case ``varIN`` and ``Lfull`` were reshaped above
    dimsIN = varIN.shape

    dimsOUT = tuple(np.append(Nlev, dimsIN[1:]))
    # ``Ndim`` is the product of all dimensions except the vertical axis
    Ndim = int(np.prod(dimsIN[1:]))
    # Flatten the other dimensions to ``[Nfull, Ndim]``
    varIN = np.reshape(varIN, (Nfull, Ndim))
    # Flatten the other dimensions to ``[Nfull, Ndim]``
    Lfull = np.reshape(Lfull, (Nfull, Ndim))
    varOUT = np.zeros((Nlev, Ndim))
    # All indices (does not change)
    Ndimall = np.arange(0, Ndim)

    if reverse_input:
        Lfull = Lfull[::-1, :]
        varIN = varIN[::-1, :]

    for k in range(0, Nlev):
        # Find nearest layer to Llev[k]
        if np.any(index):
            # Indices have been pre-computed:
            n = index[k, :]
        else:
            # Compute index on the fly for that layer. Note that
            # ``reversed_input`` is always set to ``False``. ``Lfull``
            # was reversed earlier
            n = np.squeeze(find_n(Lfull, Llev[k], False))

        # Fast method (no loop)
        # Convert the layers (``nindex``) to indices for a 2D matrix
        # using nindex = i*ncol + j
        nindex = n*Ndim + Ndimall
        nindexp1 = (n + 1)*Ndim + Ndimall

        # Initialize alpha (size = ``[Ndim]``)
        alpha = np.nan * Ndimall
        # Only calculate alpha where ``nindex < Nfull``
        Ndo = Ndimall[nindexp1 < Nfull*Ndim]
        if type_int == 'log':
            alpha[Ndo] = (np.log(Llev[k] / Lfull.flatten()[nindexp1[Ndo]])
                          / np.log(Lfull.flatten()[nindex[Ndo]]
                                   / Lfull.flatten()[nindexp1[Ndo]]))
        elif type_int == 'lin':
            alpha[Ndo] = ((Llev[k] - Lfull.flatten()[nindexp1[Ndo]])
                          / (Lfull.flatten()[nindex[Ndo]]
                             - Lfull.flatten()[nindexp1[Ndo]]))

        # Mask if ``Llev[k]`` < model top for the pressure interpolation
        if masktop:
            alpha[Llev[k] < Lfull.flatten()[nindex]] = np.nan

        # Ensure ``n+1`` is never > ``Nfull`` by setting ``n+1 = Nfull``
        # if ever ``n+1 > Nfull``. This does not affect the calculation
        # as alpha is set to NaN for those values.
        nindexp1[nindexp1 >= Nfull*Ndim] = nindex[nindexp1 >= Nfull*Ndim]

        varOUT[k, :] = (varIN.flatten()[nindex] * alpha
                        + (1-alpha) * varIN.flatten()[nindexp1])
    return np.reshape(varOUT, dimsOUT)


def axis_interp(var_IN, x, xi, axis, reverse_input=False, type_int="lin",
                modulo=None):
    """
    One dimensional linear/logarithmic interpolation along one axis.

    :param var_IN: Variable on a REGULAR grid (e.g.,
        ``[lev, lat, lon]`` or ``[time, lev, lat, lon]``)
    :type  var_IN: ND array
    :param x: Original position array (e.g., ``time``)
    :type  x: 1D array
    :param xi: Target array to interpolate the array on
    :type  xi: 1D array
    :param axis: Position of the interpolation axis (e.g., ``0`` for a
        temporal interpolation on ``[time, lev, lat, lon]``)
    :type  axis: int
    :param reverse_input: Reverse input arrays (e.g., if
        ``zfull(0)``= 120 km, ``zfull(N)``= 0 km, which is typical)
    :type  reverse_input: bool
    :param type_int: "log" for logarithmic (typically pressure),
        "lin" for linear
    :type  type_int: str
    :param modulo: For "lin" interpolation only, use cyclic input
        (e.g., when using ``modulo = 24`` for time of day, 23.5 and
        00 am are considered 30 min apart, not 23.5 hr apart)
    :type  modulo: float
    :return: ``VAR_OUT`` interpolated data on the requested axis

    .. note::
        This routine is similar but simpler than the vertical
        interpolation ``vinterp()`` as the interpolation axis is
        assumed to be fully defined by a 1D array such as ``time``,
        ``pstd`` or ``zstd`` rather than 3D arrays like ``pfull`` or
        ``zfull``.
        For lon/lat interpolation, consider using ``interp_KDTree()``.

    Calculation::

        X_OUT = Xn*A + (1-A)*Xn + 1
        with ``A = log(xi/xn + 1) / log(xn/xn + 1)`` in "log" mode
        or ``A = (xi-xn + 1)/(xn-xn + 1)`` in "lin" mode
    """

    # Convert list to numpy array as needed
    if type(var_IN) != np.ndarray:
        var_IN = np.array(var_IN)
    # Move interpolated axis to 1st axis:
    var_IN = np.moveaxis(var_IN, axis, 0)
    if reverse_input:
        var_IN = var_IN[::-1, ...]
        x = x[::-1]

    # This is called every time as it is fast on a 1D array
    index = find_n(x, xi, False)

    dimsIN = var_IN.shape
    dimsOUT = tuple(np.append(len(xi), dimsIN[1:]))
    var_OUT = np.zeros(dimsOUT)

    for k in range(0, len(index)):
        n = index[k]
        np1 = n+1
        # Treatment of edge cases where the interpolated value is
        # outside the domain (i.e. ``n`` is the last element and ``n+1``
        # does not exist)
        if np1 >= len(x):
            if modulo is not None:
                # If looping around (e.g., longitude, time of day...)
                # replace ``n+1`` by the first element
                np1 = 0
            else:
                # Sets the interpolated value to NaN in ``xi`` as last
                # value as ``x[n] - x[np1] = 0``
                np1 -= 1

        if n == -1 and modulo is None:
            n = 0
            # Also set ``n = n+1`` (which results in NaN) if ``n = -1``
            # (i.e., if the requested value is samller than first
            # element array) and the values are NOT cyclic.

        if type_int == "log":
            # Add error handling to avoid logarithm and division issues
            if x[np1] <= 0 or xi[k] <= 0 or x[n] <= 0 or x[n] == x[np1]:
                alpha = 0  # Default to 0 if we can't compute logarithm
                var_OUT[k, :] = var_IN[np1, ...]  # Use nearest value
            else:
                try:
                    alpha = (np.log(xi[k]/x[np1]) / np.log(x[n]/x[np1]))
                    var_OUT[k, :] = (var_IN[n, ...]*alpha + (1-alpha)*var_IN[np1, ...])
                except:
                    # Handle any other errors by using nearest value
                    alpha = 0
                    var_OUT[k, :] = var_IN[np1, ...]
        elif type_int == "lin":
            if modulo is None:
                if x[n] == x[np1]:
                    # Avoid division by zero
                    alpha = 0
                    var_OUT[k, :] = var_IN[np1, ...]
                else:
                    alpha = (xi[k] - x[np1]) / (x[n] - x[np1])
                    var_OUT[k, :] = (var_IN[n, ...]*alpha + (1-alpha)*var_IN[np1, ...])
            else:
                # Handle modulo case with similar error checking
                denom = np.mod(x[n]-x[np1] + modulo, modulo)
                if denom == 0:
                    # Avoid division by zero
                    alpha = 0
                    var_OUT[k, :] = var_IN[np1, ...]
                else:
                    alpha = np.mod(xi[k]-x[np1] + modulo, modulo) / denom
                    var_OUT[k, :] = (var_IN[n, ...]*alpha + (1-alpha)*var_IN[np1, ...])

    return np.moveaxis(var_OUT, 0, axis)


def layers_mid_point_to_boundary(pfull, sfc_val):
    """
    A general description for the layer boundaries is::

        p_half = ps*bk + pk

    This routine converts the coordinate of the layer MIDPOINTS,
    ``p_full`` or ``bk``, into the coordinate of the layer BOUNDARIES
    ``p_half``. The surface value must be provided.

    :param p_full: Pressure/sigma values for the layer MIDPOINTS,
        INCREASING with ``N`` (e.g., [0.01 -> 720] or [0.001 -> 1])
    :type  p_full: 1D array
    :param sfc_val: The surface value for the lowest layer's boundary
        ``p_half[N]`` (e.g., ``sfc_val`` = 720 Pa or ``sfc_val`` = 1 in
        sigma coordinates)
    :type  sfc_val: float
    :return: ``p_half`` the pressure at the layer boundaries
        (size = ``N+1``)

    Structure::

        --- 0 --- TOP   ========  p_half
        --- 1 ---
                        --------  p_full

                        ========  p_half
        ---Nk-1---      --------  p_full
        --- Nk --- SFC  ========  p_half
                        / / / / /

    We have::

        pfull[N] = ((phalf[N]-phalf[N-1]) / np.log(phalf[N]/phalf[N-1]))
        => phalf[N-1] - pfull[N] log(phalf[N-1])
        = phalf[N] - pfull[N] log(phalf[N])

    We want to solve for ``phalf[N-1] = X``::

        v                v                             v
        X      - pfull[N]       log(X)   =             B

    ``=> X= -pfull[N] W{-exp(-B/pfull[N])/pfull[N]}``

    with ``B = phalf[N] - pfull[N] log(phalf[N])`` (known at N) and

    ``W`` is the product-log (Lambert) function.

    This was tested on an L30 simulation: The values of ``phalf`` are
    reconstructed from ``pfull`` with a max error of:

    ``100*(phalf - phalf_reconstruct)/phalf < 0.4%`` at the top.
    """

    def lambertW_approx(x):
        # Internal Function. Uniform approximation for the product-log
        # function
        A = 2.344
        B = 0.8842
        C = 0.9294
        D = 0.5106
        E = -1.213
        y = np.sqrt(2*np.e*x + 2)
        return ((2*np.log(1+B*y) - np.log(1 + C*np.log(1+D*y)) + E)
                / (1 + 1./(2*np.log(1+B*y) + 2*A)))

    N = len(pfull)
    phalf = np.zeros(N+1)
    phalf[N] = sfc_val

    for i in range(N, 0, -1):
        B = phalf[i] - pfull[i-1]*np.log(phalf[i])
        phalf[i-1] = (-pfull[i-1] * lambertW_approx(-np.exp(-B / pfull[i-1])
                                                    / pfull[i-1]))
    return phalf


def polar2XYZ(lon, lat, alt, Re=3400*10**3):
    """
    Spherical to cartesian coordinate transformation.

    :param lon: Longitude in radians
    :type  lon: ND array
    :param lat: Latitude in radians
    :type  lat: ND array
    :param alt: Altitude [m]
    :type  alt: ND array
    :param Re: Planetary radius [m], defaults to 3400*10^3
    :type  Re: float
    :return: ``X``, ``Y``, ``Z`` in cartesian coordinates [m]

    .. note::
        This is a classic polar coordinate system with
        ``colatitude = pi/2 - lat`` where ``cos(colat) = sin(lat)``
    """
    lon = np.array(lon)
    lat = np.array(lat)
    alt = np.array(alt)
    R = Re + alt
    X = R * np.cos(lon)*np.cos(lat)
    Y = R * np.sin(lon)*np.cos(lat)
    # Added in case broadcasted variables are used (e.g.,
    # ``[1, lat, 1]`` or ``[1, 1, lon]``)
    Z = R * np.sin(lat)*np.ones_like(lon)
    return X, Y, Z


def interp_KDTree(var_IN, lat_IN, lon_IN, lat_OUT, lon_OUT, N_nearest=10):
    """
    Inverse distance-weighted interpolation using nearest neighboor for
    ND variables. Alex Kling, May 2021

    :param var_IN: ND variable to regrid (e.g., ``[lev, lat, lon]``,
        ``[time, lev, lat, lon]`` with ``[lat, lon]`` dimensions LAST
        [°])
    :type  var_IN: ND array
    :param lat_IN: latitude [°] (``LAT[y, x]`` array for
        irregular grids)
    :type  lat_IN: 1D or 2D array
    :param lon_IN: latitude [°] (``LAT[y, x]`` array for
        irregular grids)
    :type  lon_IN: 1D or 2D array
    :param lat_OUT: latitude [°] for the TARGET grid structure
        (or ``LAT1[y,x]`` for irregular grids)
    :type  lat_OUT: 1D or 2D array
    :param lon_OUT: longitude [°] for the TARGET grid structure
        (or ``LON1[y,x]`` for irregular grids)
    :type  lon_OUT: 1D or 2D array
    :param N_nearest: number of nearest neighbours for the search
    :type  N_nearest: int
    :return: ``VAR_OUT`` interpolated data on the target grid

    .. note::
        This implementation is much FASTER than ``griddata`` and
        it supports unstructured grids like an MGCM tile.
        The nearest neighbour interpolation is only done on the lon/lat
        axis (not level). Although this interpolation works well on the
        3D field [x, y, z], this is typically not what is expected. In
        a 4°x4° run, the closest points in all directions (N, E, S, W)
        on the target grid are 100's of km away while the closest
        points in the vertical are a few 10's -100's meter in the PBL.
        This would result in excessive weighting in the vertical.
    """

    dimsIN = var_IN.shape
    nlon_IN = dimsIN[-1]
    nlat_IN = dimsIN[-2]

    # If var is 2D, extend the dimensions for generality
    if len(dimsIN) == 2:
        var_IN = var_IN.reshape(1, nlat_IN, nlon_IN)

    # If input/output latitudes/longitudes are 1D, extend the
    # dimensions for generality:
    if len(lat_IN.shape) == 1:
        # TODO broadcast instead?
        lon_IN, lat_IN = np.meshgrid(lon_IN, lat_IN)
    if len(lat_OUT.shape) == 1:
        lon_OUT, lat_OUT = np.meshgrid(lon_OUT, lat_OUT)

    nlat_OUT = lat_OUT.shape[0]
    nlon_OUT = lon_OUT.shape[1]

    # If ``lat``, ``lon`` are 1D, broadcast dimensions:
    # ``Ndim`` = product of all input dimensions but ``lat`` & ``lon``
    Ndim = int(np.prod(dimsIN[0:-2]))
    dims_IN_reshape = tuple(np.append(Ndim, nlon_IN*nlat_IN))
    dims_OUT_reshape = tuple(np.append(Ndim, nlat_OUT*nlon_OUT))
    # Needed if var is ``lat``, ``lon``
    dims_OUT = np.append(dimsIN[0:-2], [nlat_OUT, nlon_OUT]).astype(int)
    # Initialization
    var_OUT = np.zeros(dims_OUT_reshape)
    # All indices (does not change)
    Ndimall = np.arange(0, Ndim)

    # Compute cartesian coordinate for source and target files
    # ``polar2XYZ(lon, lat, lev)``
    xs, ys, zs = polar2XYZ(lon_IN*np.pi/180, lat_IN*np.pi/180, 0., Re = 1.)
    xt, yt, zt = polar2XYZ(lon_OUT*np.pi/180, lat_OUT*np.pi/180, 0., Re = 1.)

    tree = cKDTree(list(zip(xs.flatten(), ys.flatten(), zs.flatten())))
    d, inds = tree.query(list(zip(xt.flatten(), yt.flatten(), zt.flatten())),
                         k = N_nearest)
    # Inverse distance
    w = 1.0 / d**2
    # Sum the weights and normalize
    var_OUT = (np.sum(w*var_IN.reshape(dims_IN_reshape)[:, inds], axis = 2)
               / np.sum(w, axis = 1))
    return var_OUT.reshape(dims_OUT)


def cart_to_azimut_TR(u, v, mode="from"):
    """
    Convert cartesian coordinates or wind vectors to radians using azimuth angle.

    :param x: the cartesian coordinate
    :type  x: 1D array
    :param y: the cartesian coordinate
    :type  y: 1D array
    :param mode: "to" for the direction that the vector is pointing,
        "from" for the direction from which the vector is coming
    :type  mode: str
    :return: ``Theta`` [°] and ``R`` the polar coordinates
    """

    if mode == "from":
        cst = 180
    if mode == "to":
        cst = 0.
    return (np.mod(np.arctan2(u, v)*180/np.pi+  cst, 360),
            np.sqrt(u**2 + v**2))


def area_weights_deg(var_shape, lat_c, axis = -2):
    """
    Returns weights scaled so that np.mean(var*W) gives an area-weighted 
    average. This works because grid cells near the poles have smaller 
    areas than those at the equator, so they should contribute less to 
    a global average.

    Expected dimensions are:

    [lat] ``axis`` not needed
    [lat, lon] ``axis = -2`` or ``axis = 0``
    [time, lat, lon] ``axis = -2`` or ``axis = 1``
    [time, lev, lat, lon] ``axis = -2`` or ``axis = 2``
    [time, time_of_day_24, lat, lon] ``axis = -2`` or ``axis = 2``
    [time, time_of_day_24, lev, lat, lon] ``axis = -2`` or ``axis = 3``

    Because ``dlat`` is computed as ``lat_c[1]-lat_c[0]``, ``lat_c``
    may be truncated on either end (e.g., ``lat = [-20 ..., 0... 50]``)
    but must be continuous.

    :param var_shape: the shape/dimensions of your data array
    :type  var_shape: tuple
    :param lat_c: latitude cell centers in degrees [°]
    :type  lat_c: float
    :param axis: which dimension contains latitude, default: 2nd-to-last
    :type  axis: int
    :return: ``W`` weights for the variable ready for standard
        averaging as ``np.mean(var*W)`` [condensed form] or
        ``np.average(var, weights=W)`` [expanded form]

    .. note::
        Given a variable var:

        ``var = [v1, v2, ...vn]``

        The regular average is

        ``AVG = (v1 + v2 + ... vn) / N``

        and the weighted average is

        ``AVG_W = (v1*w1 + v2*w2 + ... vn*wn) / (w1 + w2 + ...wn)``

        This function returns

        ``W = [w1, w2, ... , wn]*N / (w1 + w2 + ...wn)``

        Therfore taking a regular average of (``var*W``) with
        ``np.mean(var*W)`` or ``np.average(var, weights=W)``

        returns the weighted average of the variable. Use
        ``np.average(var, weights=W, axis = X)`` to average over a
        specific axis.
    """
    # Check if either the lat array or var shape is essentially
    # scalar (single value)
    # np.atleast_1d() ensures the input is treated as at least a 1D 
    # array for checking its length
    if (len(np.atleast_1d(lat_c)) == 1 or
        len(np.atleast_1d(var_shape)) == 1):
        # If either is scalar, returns an array of 1s matching the var 
        # shape (no weighting needed since there's nothing to weight across)
        return np.ones(var_shape)
    else:
        # Calculates lat spacing by taking the difference between the 
        # first two latitude points. Assumes uniform grid spacing
        dlat = lat_c[1]-lat_c[0]
        
        # Calculate cell areas 
        # Use dlon = 1 (arbitrary lon spacing and planet 
        # radius because normalization makes the absolute values 
        # irrelevant), then normalize
        R = 1.
        dlon = 1.
        
        # Compute total surface area for normalization
        lon1_deg = -dlon/2
        lon2_deg = dlon/2
        lat1_deg = lat_c[0] - dlat/2
        lat2_deg = lat_c[-1] + dlat/2
        
        # Convert to radians for area calculation
        lat1_rad = lat1_deg * np.pi/180
        lat2_rad = lat2_deg * np.pi/180
        lon1_rad = lon1_deg * np.pi/180
        lon2_rad = lon2_deg * np.pi/180
        
        area_tot = ((R**2) 
                    * np.abs(lon1_rad - lon2_rad) 
                    * np.abs(np.sin(lat1_rad) - np.sin(lat2_rad)))
        
        # Convert to radians
        lat_c_rad = lat_c * np.pi/180
        dlon_rad = dlon * np.pi/180
        dlat_rad = dlat * np.pi/180
        
        # Calculate normalized areas. Areas sum to 1
        A = (2. * R**2 * dlon_rad * np.cos(lat_c_rad) * 
             np.sin(dlat_rad/2.) / area_tot)

        # Check if var is 1D (just lat values)
        if len(var_shape) == 1:
            # For 1D case: multiply normalized areas by the number of 
            # lat points. Creates weights where sum(W) = len(lat_c), 
            # allowing np.mean(var*W) to give the area-weighted average
            W = A * len(lat_c)
        else:
            # For multidimensional data: create a list of 1s with 
            # length matching the number of dims in the var 
            # (e.g., [1, 1, 1, 1] for 4D data).
            reshape_shape = [1 for i in range(0, len(var_shape))]
            # Set the lat dim to the actual number of lat points 
            # (e.g., [1, 1, lat, 1] for 4D data where lat = third dim)
            reshape_shape[axis] = len(lat_c)
            # Reshape the 1D area array to match the var's 
            # dimensionality (broadcasting-ready shape), then multiply 
            # by the number of lat points. Allows the weights to 
            # broadcast correctly across all other dims.
            W = A.reshape(reshape_shape)*len(lat_c)
        
        # Expand the weights to full var shape by multiplying with an 
        # array of 1s. Creates the final weight array that can be 
        # directly multiplied with other data for area-weighted avg.
        return W*np.ones(var_shape)


def areo_avg(VAR, areo, Ls_target, Ls_angle, symmetric=True):
    """
    Return a value average over a central solar longitude

    EX::

        ``Ls_target = 90.``
        ``Ls_angle = 10.``

    Nominally, the time average is done over solar longitudes
    ``85 < Ls_target < 95`` (10°).

    If ``symmetric = True`` and the input data range is Ls = 88-100°
    then ``88 < Ls_target < 92`` (4°, symmetric)

    If ``symmetric = False`` and the input data range is Ls = 88-100°
    then ``88 < Ls_target < 95`` (7°, assymetric)

    :param VAR: a variable with ``time`` in the 1st dimension
    :type  VAR: ND array
    :param areo: solar longitude of the input variable (0-720)
    :type  areo: 1D array
    :param Ls_target: central solar longitude of interest
    :type  Ls_target: float
    :param Ls_angle: requested window angle centered at ``Ls_target``
    :type  Ls_angle: float
    :param symmetric: If ``True`` and the requested window is out of range,
        ``Ls_angle`` is reduced. If False, the time average is performed
        on the data available
    :type  symmetric: bool (defaults to True)
    :return: the variable averaged over solar longitudes
        ``Ls_target-Ls_angle/2`` to ``Ls_target+Ls_angle/2``

    .. note::
        The routine can bin data from muliples Mars years
    """

    # Take the modulo of solar longitude
    areo = np.mod(areo, 360)
    # All dimensions but time
    shape_out = VAR.shape[1:]
    # Flatten array
    VAR = VAR.reshape((len(areo), np.prod(shape_out)))

    # Compute bounds from ``Ls_target`` and ``Ls_angle``
    Ls_min = Ls_target - Ls_angle/2.
    Ls_max = Ls_target + Ls_angle/2.

    if (Ls_min < 0.):
        Ls_min += 360.
    if (Ls_max > 360.):
        Ls_max -= 360.

    # Initialize output array
    VAR_avg = np.zeros(np.prod(shape_out))

    # EX: This is removed if 10° of data are requested around Ls 0:
    #                    Ls 355 <-- (0.00) --> 5
    # and the file is    Ls 1   <-- (180)  --> 357
    # the data selected should be 1 > 5 and 355 > 357

    # Check if the Ls of interest is within the range of data, raise
    # execption otherwise
    if Ls_target <= areo.min() or Ls_target >= areo.max():
        raise Exception(
            f"Error\nNo data found, requested data range:\n"
            f"Ls {Ls_min:.3} <-- ({Ls_target:.3})--> {Ls_max:.3}\n"
            f"However, the available data range is:\n"
            f"Ls {areo.min():.3} <-- ({(areo.min()+areo.max())/2.:.3}) --> "
            f"{areo.max():.3}")

    if Ls_min < areo.min() or Ls_max > areo.max():
        print(f"In ``areo_avg()`` Warning:\nRequested data range:\n"
              f"Ls {Ls_min:.3} <-- ({Ls_target:.3})--> {Ls_max:.3}")
        if symmetric:
            # Case 1: reduce the window
            if Ls_min < areo.min():
                Ls_min = areo.min()
                Ls_angle = 2*(Ls_target - Ls_min)
                Ls_max = Ls_target + Ls_angle/2.

            if Ls_max > areo.max():
                Ls_max = areo.max()
                Ls_angle = 2*(Ls_max - Ls_target)
                Ls_min = Ls_target - Ls_angle/2.

            print(f"Reshaping data ranging Ls "
                  f"{Ls_min:.3} <-- ({Ls_target:.3})--> {Ls_max:.3}")
        else:
            # Case 2: use all data available
            print(f"Only using data ranging Ls "
                  f"{max(areo.min(), Ls_min):.3} <-- ({Ls_target:.3})--> "
                  f"{min(areo.max(), Ls_max):.3} \n")
    count = 0

    for t in range(len(areo)):
        if (Ls_min <= areo[t] <= Ls_max):
            # Special case: Ls around Ls = 0 (wrap around)
            VAR_avg += VAR[t, ...]
            count += 1

    if count > 0:
        VAR_avg /= count
    return VAR_avg.reshape(shape_out)


def mass_stream(v_avg, lat, level, type="pstd", psfc=700, H=8000.,
                factor=1.e-8):
    """
    Compute the mass stream function::

                                P
                                ⌠
        Ph i= (2 pi a) cos(lat)/g ⎮vz_tavg dp
                                ⌡
                                p_top

    :param v_avg: zonal wind [m/s] with ``lev`` dimension FIRST and
        ``lat`` dimension SECOND (e.g., ``[pstd, lat]``,
        ``[pstd, lat, lon]`` or ``[pstd, lat, lon, time]``)
    :type  v_avg: ND array
    :param lat: latitudes [°]
    :type  lat: 1D array
    :param level: interpolated layers [Pa] or [m]
    :type  level: 1D array
    :param type: interpolation type (``pstd``, ``zstd`` or ``zagl``)
    :type  type: str
    :param psfc: reference surface pressure [Pa]
    :type  psfc: float
    :param H: reference scale height [m] when pressures are used
    :type  H: float
    :param factor: normalize the mass stream function by a factor, use
        ``factor = 1`` for [kg/s]
    :type  factor: int
    :return: ``MSF`` the meridional mass stream function (in
        ``factor * [kg/s]``)

    .. note::
        This routine allows the time and zonal averages to be
        computed before OR after the MSF calculation.

    .. note::
        The expressions for MSF use log(pressure) Z coordinates,
        which integrate better numerically.

        With ``p = p_sfc exp(-Z/H)`` and ``Z = H log(p_sfc/p)``
        then ``dp = -p_sfc/H exp(-Z/H) dZ`` and we have::

                                            Z_top
                                            ⌠
            Phi = +(2pi a)cos(lat)psfc/(gH) ⎮v_rmv exp(-Z/H)dZ
                                            ⌡
                                            Z
        With ``p = p_sfc exp(-Z/H)``

        The integral is calculated using trapezoidal rule::

                n
                ⌠
            .g. ⌡ f(z)dz = (Zn-Zn-1){f(Zn) + f(Zn-1)}/2
              n-1
    """
    reverting = False
    if level[0] < level[-1]:
        reverting = True
        print("Reversing pstd array for mass stream function calculation...")
        level = level[::-1]
        v_avg = v_avg[::-1, ...]
                    
    g = 3.72 # m/s2
    a = 3400*1000 # m
    nlev = len(level)
    shape_out = v_avg.shape

    # If size is ``[pstd, lat]``, convert to ``[pstd, lat, 1]`` for
    # generality
    if len(shape_out) == 2:
        v_avg = v_avg.reshape(nlev, len(lat), 1)

    # Flatten array
    v_avg = v_avg.reshape((nlev, len(lat), np.prod(v_avg.shape[2:])))
    MSF = np.zeros_like(v_avg)

    # Sum variable, same dimensions as ``v_avg`` but for first dimension
    I = np.zeros(v_avg.shape[2:])

    # Replace NaN with 0 for downward integration
    isNan = False
    if np.isnan(v_avg).any():
        isNan = True
        mask = np.isnan(v_avg)
        v_avg[mask] = 0.

    isMasked = False
    if np.ma.is_masked(v_avg):
        # Missing data may be masked instead of set to NaN
        isMasked = True
        mask0 = np.ma.getmaskarray(v_avg)
        # Make a standalone copy of the mask array
        mask = mask0.copy()
        # Set masked elements to ``0.`` Note that this effectively
        # unmasks the array as ``0.`` is a valid entry.
        v_avg[mask0] = 0.

    if type == "pstd":
        Z = H * np.log(psfc/level)
    else:
        # Copy ``zagl`` or ``zstd`` instead of using a pseudo height
        Z = level.copy()

    for k0 in range(nlev-2, 0, -1):
        I[:] = 0.
        for k in range(nlev-2, k0, -1):
            zn = Z[k]
            znp1 = Z[k+1]
            fn = v_avg[k, :, ...] * np.exp(-zn/H)
            fnp1 = v_avg[k+1, :, ...] * np.exp(-znp1/H)
            I = I + 0.5 * (znp1-zn) * (fnp1+fn)
        MSF[k0, :, ...] = (2 * np.pi * a * psfc
                           / (g*H)
                           * np.cos(np.pi/180*lat).reshape([len(lat), 1])
                           * I * factor)

    # Put NaNs back to where they initially were
    if isNan:
        MSF[mask] = np.nan
    if isMasked:
        MSF = np.ma.array(MSF, mask = mask)
    
    if reverting:
        print("Reversing pstd dimension of MSF array for compatibility...")
        MSF = MSF[::-1, ...]
        
    return MSF.reshape(shape_out)


def vw_from_MSF(msf, lat, lev, ztype="pstd", norm=True, psfc=700, H=8000.):
    """
    Return the V and W components of the circulation from the mass
    stream function.

    :param msf: the mass stream function with ``lev`` SECOND TO
        LAST and the ``lat`` dimension LAST (e.g., ``[lev, lat]``,
        ``[time, lev, lat]``, ``[time, lon, lev, lat]``)
    :type  msf: ND array
    :param lat: latitude [°]
    :type  lat: 1D array
    :param lev: level [Pa] or [m] (``pstd``, ``zagl``, ``zstd``)
    :type  lev: 1D array
    :param ztype: Use ``pstd`` for pressure so vertical
        differentation is done in log space.
    :type  ztype: str
    :param norm: if True, normalize ``lat`` and ``lev`` before
        differentiating to avoid having to rescale manually the
        vectors in quiver plots
    :type  norm: bool
    :param psfc: surface pressure for pseudo-height when
        ``ztype = pstd``
    :type  psfc: float
    :param H: scale height for pseudo-height when ``ztype = pstd``
    :type  H: float
    :return: the meditional and altitude components of the mass stream
        function for plotting as a quiver or streamlines.

    .. note::
        The components are:
        ``[v]=  g/(2 pi cos(lat)) dphi/dz``
        ``[w]= -g/(2 pi cos(lat)) dphi/dlat``
    """

    g = 3.72 # m/s2

    lat = lat * np.pi/180
    var_shape = msf.shape

    xx = lat.copy()
    zz = lev.copy()

    if ztype == "pstd":
        zz = H * np.log(psfc/lev)

    if norm:
        xx = (xx-xx.min()) / (xx.max()-xx.min())
        zz = (zz-zz.min()) / (zz.max()-zz.min())

    # Extend broadcasting dimensions for the latitude (``[1, 1, lat]``
    # if ``msf`` is size ``[time, lev, lat]``)
    reshape_shape = [1 for i in range(0, len(var_shape))]
    reshape_shape[-1] = lat.shape[0]
    lat1d = lat.reshape(reshape_shape)

    # Transpose shapes:
    T_array = np.arange(len(msf.shape))

    # One permutation only: ``lat`` is passed to the 1st dimension
    T_latIN = np.append(T_array[-1], T_array[0:-1])
    T_latOUT = np.append(T_array[1:], T_array[0])

    T_levIN = np.append(np.append(T_array[-2], T_array[0:-2]), T_array[-1])
    T_levOUT = np.append(np.append(T_array[1:-1], T_array[0]), T_array[-1])

    V = (g / (2*np.pi*np.cos(lat1d))
         * dvar_dh(msf.transpose(T_levIN), zz).transpose(T_levOUT))
    W = (-g / (2*np.pi*np.cos(lat1d))
         * dvar_dh(msf.transpose(T_latIN), xx).transpose(T_latOUT))
    return V, W


def alt_KM(press, scale_height_KM=8., reference_press=610.):
    """
    Gives the approximate altitude [km] for a given pressure

    :param press: the pressure [Pa]
    :type  press: 1D array
    :param scale_height_KM: scale height [km] (default is 8 km, an
        isothermal at 155K)
    :type  scale_height_KM: float
    :param reference_press: reference surface pressure [Pa] (default is
        610 Pa)
    :type  reference_press: float
    :return: ``z_KM`` the equivalent altitude for that pressure [km]

    .. note::
        Scale height is ``H = rT/g``
    """

    # Pressure -> altitude [km]
    return (-scale_height_KM * np.log(press/reference_press))


def press_pa(alt_KM, scale_height_KM=8., reference_press=610.):
    """
    Gives the approximate altitude [km] for a given pressure

    :param alt_KM: the altitude [km]
    :type  alt_KM: 1D array
    :param scale_height_KM: scale height [km] (default is 8 km, an
        isothermal at 155K)
    :type  scale_height_KM: float
    :param reference_press: reference surface pressure [Pa] (default is
        610 Pa)
    :type  reference_press: float
    :return: ``press_pa`` the equivalent pressure at that altitude [Pa]

    .. note::
        Scale height is ``H = rT/g``
    """

    return (reference_press * np.exp(-alt_KM/scale_height_KM))


def lon180_to_360(lon):
    """
    Transform a float or an array from the -180/180 coordinate system
    to 0-360

    :param lon: longitudes in the -180/180 coordinate system
    :type  lon: float, 1D array, or 2D array
    :return: the equivalent longitudes in the 0-360 coordinate system
    """

    lon = np.array(lon)

    if len(np.atleast_1d(lon)) == 1:
        # ``lon`` is a float
        if lon < 0:
            lon += 360
    else:
        # ``lon`` is an array
        lon[lon < 0] += 360
        # Reogranize lon by increasing values
        lon = np.append(lon[lon <= 180], lon[lon > 180])
    return lon


def lon360_to_180(lon):
    """
    Transform a float or an array from the 0-360 coordinate system to
        -180/180.

    :param lon: longitudes in the 0-360 coordinate system
    :type  lon: float, 1D array, or 2D array
    :return: the equivalent longitudes in the -180/180 coordinate system
    """

    lon = np.array(lon)
    if len(np.atleast_1d(lon)) == 1:
        # ``lon`` is a float
        if lon > 180:
            lon -= 360
    else:
        # ``lon`` is an array
        lon[lon > 180] -= 360
        # Reogranize lon by increasing values
        lon = np.append(lon[lon < 0], lon[lon >= 0])
    return lon


def shiftgrid_360_to_180(lon, data):
    """
    This function shifts ND data from a 0-360 to a -180/180 grid.

    :param lon: longitudes in the 0-360 coordinate system
    :type  lon: 1D array
    :param data: variable with ``lon`` in the last dimension
    :type  data: ND array
    :return: shifted data

    .. note::
        Use ``np.ma.hstack`` instead of ``np.hstack`` to keep the
        masked array properties
    """

    lon = np.array(lon)
    # convert to +/- 180
    lon[lon > 180] -= 360.
    # stack data
    data = np.concatenate((data[..., lon < 0], data[..., lon >= 0]), axis = -1)
    return data


def shiftgrid_180_to_360(lon, data):
    """
    This function shifts ND data from a -180/180 to a 0-360 grid.

    :param lon: longitudes in the 0-360 coordinate system
    :type  lon: 1D array
    :param data: variable with ``lon`` in the last dimension
    :type  data: ND array
    :return: shifted data
    """

    lon = np.array(lon)
    # convert to 0-360
    lon[lon < 0] += 360.
    # stack data
    data = np.concatenate((data[..., lon <= 180], data[..., lon > 180]),
                          axis = -1)
    return data


def second_hhmmss(seconds, lon_180=0.):
    """
    Given the time [sec], return local true solar time at a
    certain longitude.

    :param seconds: the time [sec]
    :type  seconds: float
    :param lon_180: the longitude in -180/180 coordinate
    :type  lon_180: float
    :return: the local time [float] or a tuple (hours, minutes, seconds)
    """

    hours = seconds // (60*60)
    seconds %= (60 * 60)
    minutes = seconds // 60
    seconds %= 60
    # Add timezone offset (1hr/15°)
    hours = np.mod(hours + lon_180/15., 24)
    return (np.int32(hours), np.int32(minutes), np.int32(seconds))


def sol_hhmmss(time_sol, lon_180=0.):
    """
    Given the time in days, return return local true solar time at a
    certain longitude.

    :param time_sol: the time in sols
    :type  seconds: float
    :param lon_180: the longitude in -180/180 coordinate
    :type  lon_180: float
    :return: the local time [float] or a tuple (hours, minutes, seconds)
    """

    return second_hhmmss(time_sol * 86400., lon_180)


def UT_LTtxt(UT_sol, lon_180=0., roundmin=None):
    """
    Returns the time in HH:MM:SS at a certain longitude.

    :param time_sol: the time in sols
    :type  time_sol: float
    :param lon_180: the center longitude in -180/180 coordinates.
        Increments by 1hr every 15°
    :type  lon_180: float
    :param roundmin: round to the nearest X minute. Typical values are
        ``roundmin = 1, 15, 60``
    :type  roundmin: int

    .. note::
        If ``roundmin`` is requested, seconds are not shown
    """

    def round2num(number, interval):
        # Internal Function to round a number to the closest range.
        # e.g., ``round2num(26, 5) = 25``, ``round2num(28, 5) = 30``
        return round(number / interval) * interval

    hh, mm, ss = sol_hhmmss(UT_sol, lon_180)

    if roundmin:
        sec = hh*3600 + mm*60 + ss
        # Round to the nearest increment [sec] and run a second pass
        rounded_sec = round2num(sec, roundmin*60)
        hh, mm, ss = second_hhmmss(rounded_sec, lon_180)
        return (f"{hh:02}:{mm:02}")
    else:
        return (f"{hh:02}:{mm:02}:{ss:02}")


def dvar_dh(arr, h=None):
    """
    Differentiate an array ``A[dim1, dim2, dim3...]`` w.r.t ``h``. The
    differentiated dimension must be the first dimension.

    EX: Compute ``dT/dz`` where ``T[time, lev, lat, lon]`` is the
    temperature and ``Zkm`` is the array of  level heights [km].

    First, transpose ``T`` so the vertical dimension comes first:
    ``T[lev, time, lat, lon]``.

    Then transpose back to get ``dTdz[time, lev, lat, lon]``::

        dTdz = dvar_dh(t.transpose([1, 0, 2, 3]),
                       Zkm).transpose([1, 0, 2, 3])

    If ``h`` is 1D, then ``h``and ``dim1`` must have the same length

    If ``h`` is 2D, 3D or 4D, then ``arr`` and ``h`` must have the
    same shape

    :param arr: variable
    :type  arr: ND array
    :param h: the dimension (``Z``, ``P``, ``lat``, ``lon``)
    :type  h: str
    :return: d_arr: the array differentiated w.r.t ``h``, e.g., d(array)/dh
    """

    h = np.array(h)
    d_arr = np.zeros_like(arr)
    if h.any():
        if len(h.shape) == 1:
            # ``h`` is a 1D array
            reshape_shape = np.append(
                [arr.shape[0]-2], [1 for i in range(0, arr.ndim - 1)])
            d_arr[0, ...] = ((arr[1, ...] - arr[0, ...])
                             / (h[1] - h[0]))
            d_arr[-1, ...] = ((arr[-1, ...] - arr[-2, ...])
                              / (h[-1] - h[-2]))
            d_arr[1:-1, ...] = ((arr[2:, ...] - arr[0:-2, ...])
                                / (np.reshape(h[2:] - h[0:-2], reshape_shape)))
        elif h.shape == arr.shape:
            # ``h`` has the same shape as the variable
            d_arr[0, ...] = ((arr[1, ...] - arr[0, ...])
                             / (h[1, ...] - h[0, ...]))
            d_arr[-1, ...] = ((arr[-1, ...] - arr[-2, ...])
                              / (h[-1, ...] - h[-2, ...]))
            d_arr[1:-1, ...] = ((arr[2:, ...] - arr[0:-2, ...])
                                / (h[2:, ...] - h[0:-2, ...]))
        else:
            print(f"Error, ``h.shape=``{h.shape}, ``arr.shape=``{arr.shape}")
    else:
        # ``h`` is not defined, return ``d_var``, not ``d_var/dh``
        d_arr[0, ...] = arr[1, ...] - arr[0, ...]
        d_arr[-1, ...] = arr[-1, ...] - arr[-2, ...]
        # 0.5 factor since differentiation uses a centered scheme
        d_arr[1:-1, ...] = 0.5*(arr[2:, ...] - arr[0:-2, ...])
    return d_arr


def zonal_detrend(VAR):
    """
    Substract the zonal average mean value from a field.

    :param VAR: variable with detrending dimension last
    :type  VAR: ND array
    :return: detrented field (same size as input)

    .. note::
        ``RuntimeWarnings`` are expected if the slice contains
        only NaNs which is the case below the surface and above the
        model top in the interpolated files. This routine disables such
        warnings temporarily.
    """

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category = RuntimeWarning)
        return (VAR - np.nanmean(VAR, axis = -1)[..., np.newaxis])


def get_trend_2D(VAR, LON, LAT, type_trend="wmean"):
    """
    Extract spatial trends from the data. The output can be directly
    subtracted from the original field.

    :param VAR: Variable for decomposition. ``lat`` is SECOND TO LAST
        and ``lon`` is LAST  (e.g., ``[time, lat, lon]`` or
        ``[time, lev, lat, lon]``)
    :type  VAR: ND array
    :param LON: longitude coordinates
    :type  LON: 2D array
    :param LAT: latitude coordinates
    :type  LAT: 2D array
    :param type_trend: type of averaging to perform:
        "mean" - use a constant average over all lat/lon
        "wmean" - use a area-weighted average over all lat/lon
        "zonal" - detrend over the zonal axis only
        "2D" - use a 2D planar regression (not area-weighted)
    :type  type_trend: str
    :return: the trend, same size as ``VAR``
    """

    var_shape = np.array(VAR.shape)

    # Type "zonal" is the easiest as averaging is performed over 1
    # dimension only.
    if type_trend == "zonal":
        return (np.repeat(np.nanmean(VAR, axis = -1)[..., np.newaxis],
                          var_shape[-1], axis = -1))

    # The following options involve averaging over both lat and lon
    # dimensions

    # Flatten array (``[10, 36, lat, lon]`` -> ``[360, lat, lon]``)
    nflatten = int(np.prod(var_shape[:-2]))
    reshape_flat = np.append(nflatten, var_shape[-2:])
    VAR = VAR.reshape(reshape_flat)

    TREND = np.zeros(reshape_flat)
    for ii in range(nflatten):
        if type_trend == "mean":
            TREND[ii, ...] = np.mean(VAR[ii, ...].flatten())
        elif type_trend == "wmean":
            W = area_weights_deg(var_shape[-2:], LAT[:, 0])
            TREND[ii, ...] = np.mean((VAR[ii, ...] * W).flatten())
        elif type_trend == "2D":
            TREND[ii, ...] = regression_2D(LON, LAT, VAR[ii, :, :], order = 1)
        else:
            print(f"Error, in ``area_trend``, type '{type_trend}' not "
                  f"recognized")
            return None
    return TREND.reshape(var_shape)


def regression_2D(X, Y, VAR, order=1):
    """
    Linear and quadratic regression on the plane.

    :param X: first coordinate
    :type  X: 2D array
    :param Y: second coordinate
    :type  Y: 2D array
    :param VAR: variable of the same size as X
    :type  VAR: 2D array
    :param order: 1 (linear) or 2 (quadratic)
    :type  order: int

    .. note::
        When ``order = 1``, the equation is: ``aX + bY + C = Z``.
        When ``order = 2``, the equation is:
        ``aX^2 + 2bX*Y + cY^2 + 2dX + 2eY + f = Z``

    For the linear case::, ``ax + by + c = z`` is re-written as
    ``A X = b`` with::

                |x0   y0   1|        |a      |z0
            A = |x1   y1   1|    X = |b   b= |
                |      ...  |        |c      |...
                |xn   yn   1|                |zn

                    [n,3]           [3]       [n]

    The least-squares regression provides the solution that that
    minimizes ``||b – A x||^2``
    """

    if order == 1:
        A = np.array([X.flatten(), Y.flatten(), np.ones_like(X.flatten())]).T
        # An Equivalent notation is:
        # A = np.c_[X.flatten(), Y.flatten(), np.ones_like(X.flatten())]

        b = VAR.flatten()

        # ``P`` is the solution of  ``A X =b``
        # ==> ``P[0] x + P[1]y + P[2] = z``
        P, residuals, rank, s = np.linalg.lstsq(A, b, rcond = None)
        Z = P[0]*X + P[1]*Y + np.ones_like(X)*P[2]

    elif order == 2:
        # Best-fit quadratic curve:
        # ``aX^2 + 2bX*Y + cY^2 + 2dX + 2eY + f``
        XX = X.flatten()
        YY = Y.flatten()
        ZZ = VAR.flatten()
        data = np.zeros((len(XX), 3))
        data[:, 0] = XX
        data[:, 1] = YY
        data[:, 2] = ZZ

        A = np.c_[np.ones(data.shape[0]), data[:, :2], np.prod(
            data[:, :2], axis = 1), data[:, :2]**2]
        P, _, _, _ = np.linalg.lstsq(A, data[:, 2], rcond = None)

        # Evaluate it on a grid (using vector product)
        Z = np.dot(np.c_[np.ones(XX.shape),
                         XX, YY, XX * YY, XX**2, YY**2], P).reshape(X.shape)
    return Z


def daily_to_average(varIN, dt_in, nday=5, trim=True):
    """
    Bin a variable from an ``atmos_daily`` file format to the
    ``atmos_average`` file format.

    :param varIN: variable with ``time`` dimension first (e.g.,
        ``ts[time, lat, lon]``)
    :type  varIN: ND array
    :param dt_in: delta of time betwen timesteps in sols (e.g.,
        ``dt_in = time[1] - time[0]``)
    :type  dt_in: float
    :param nday: bining period in sols, default is 5 sols
    :type  nday: int
    :param trim: whether to discard any leftover data at the end of file
        before binning
    :type  trim: bool
    :return: the variable bin over ``nday``

    .. note::
        If ``varIN[time, lat, lon]`` from ``atmos_daily`` is
        ``[160, 48, 96]`` and has 4 timesteps per day (every 6 hours),
        then the resulting variable for ``nday = 5`` is
        ``varOUT(160/(4*5), 48, 96) = varOUT(8, 48, 96)``

    .. note::
        If the daily file has 668 sols, then there are
        ``133 x 5 + 3`` sols leftover. If ``trim = True``, then the
        time is 133 and last 3 sols the are discarded. If
        ``trim = False``, the time is 134 and last bin contains only
        3 sols of data.
    """

    vshape_in = varIN.shape
    # 0 is the time dimension
    Nin = vshape_in[0]

    # Add safety check for dt_in
    if np.isclose(dt_in, 0.0):
        print("Error: Time difference dt_in is zero or very close to zero.")
        return None

    iperday = int(np.round(1 / dt_in))
    combinedN = int(iperday * nday)
    N_even = Nin // combinedN
    N_left = Nin % combinedN


    if N_left != 0 and not trim:
        # If ``Nin/(nday * iperday)`` is not a round number
        # Do the average on the even part
        vreshape = np.append([-1, combinedN], vshape_in[1:]).astype(int)
        var_even = np.mean(varIN[0:N_even*combinedN, ...].reshape(vreshape),
                           axis = 1)
        # Left over time steps
        var_left = np.mean(varIN[N_even*combinedN:, ...], axis = 0,
                           keepdims = True)
        # Combine both
        varOUT = np.concatenate((var_even, var_left), axis = 0)
    else:
        # If ``Nin/(nday * iperday)`` is a round number, otherwise trim
        # the array
        vreshape = np.append([-1, combinedN], vshape_in[1:]).astype(int)
        varOUT = np.mean(varIN[0:N_even*combinedN, ...].reshape(vreshape),
                         axis = 1)
    return varOUT


def daily_to_diurn(varIN, time_in):
    """
    Bin a variable from an ``atmos_daily`` file into the
    ``atmos_diurn`` format.

    :param varIN: variable with time dimension first (e.g.,
        ``[time, lat, lon]``)
    :type  varIN: ND array
    :param time_in: time array in sols. Only the first N elements
        are actually required if saving memory is important
    :type  time_in: ND array
    :return: the variable binned in the ``atmos_diurn`` format
        (``[time, time_of_day, lat, lon]``) and the time of day array
        [hr]

    .. note::
        If ``varIN[time, lat, lon]`` from ``atmos_daily`` is
        ``[40, 48, 96]`` and has 4 timestep per day (every 6 hours),
        then the resulting variable is
        ``varOUT[10, 4, 48, 96] = [time, time_of_day, lat, lon]`` and
        ``tod = [0., 6., 12., 18.]``.

    .. note::
        Since the time dimension is first, the output variables
        may be passed to the ``daily_to_average()`` function for
        further binning.
    """

    dt_in = time_in[1] - time_in[0]

    # Add safety check for dt_in
    if np.isclose(dt_in, 0.0):
        print("Error: Time difference dt_in is zero or very close to zero.")
        return None

    iperday = int(np.round(1/dt_in))
    vshape_in = varIN.shape

    # Add safety check for integer sols
    if not(np.mod(vshape_in[0],iperday) == 0):
        print("Error: File cannot be split evenly into sols")
        return None

    vreshape = np.append([-1, iperday], vshape_in[1:]).astype(int)
    varOUT = varIN.reshape(vreshape)

    # Get time of day in hours
    tod = np.mod(time_in[0:iperday]*24, 24)

    # Sort by time of day (e.g., if ``tod = [6., 12., 18., 0.]``,
    # re-arange into ``[0., 6., 12., 18.]``. Every element in array is
    # greater than the one to its left
    if not np.all(tod[1:] >= tod[:-1], axis = 0):
        # This returns the permutation (if ``tod = [6., 12., 18., 0.]``,
        # ``i_sort = [3, 0, 1, 2]``)
        i_sort = sorted(range(len(tod)), key=lambda k: tod[k])
        tod = tod[i_sort]
        varOUT = varOUT[:, i_sort, ...]
    return varOUT


# ======================================================================
#                       Vertical Grid Utilities
# ======================================================================


def gauss_profile(x, alpha, x0=0.):
    """
    Return Gaussian line shape at x. This can be used to generate a
    bell-shaped mountain.
    """

    return (np.sqrt(np.log(2) / np.pi) / alpha
            * np.exp(-((x-x0) / alpha)**2 * np.log(2)))


def compute_uneven_sigma(num_levels, N_scale_heights, surf_res,
                         exponent, zero_top):
    """
    Construct an initial array of sigma based on the number of levels
    and an exponent

    :param num_levels: the number of levels
    :type  num_levels: float
    :param N_scale_heights: the number of scale heights to the top of
        the model (e.g., ``N_scale_heights`` = 12.5 ~102 km assuming an
        8 km scale height)
    :type  N_scale_heights: float
    :param surf_res: the resolution at the surface
    :type  surf_res: float
    :param exponent: an exponent to increase the thickness of the levels
    :type  exponent: float
    :param zero_top: if True, force the top pressure boundary
        (in N = 0) to 0 Pa
    :type  zero_top: bool
    :return: an array of sigma layers
    """

    b = np.zeros(int(num_levels)+1)
    for k in range(0, num_levels):
        # zeta decreases with k
        zeta = 1. - k/float(num_levels)
        z = surf_res*zeta + (1.0 - surf_res)*(zeta**exponent)
        # z goes from 1 to 0
        b[k] = np.exp(-z*N_scale_heights)
    b[-1] = 1.0
    if(zero_top):
        b[0] = 0.0
    return b


def transition(pfull, p_sigma=0.1, p_press=0.05):
    """
    Return the transition factor to construct ``ak`` and ``bk``

    In the MGCM code, the full pressures are computed from::

                       del(phalf)
         pfull = -----------------------------
                 log(phalf(k+1/2)/phalf(k-1/2))

    :param pfull: the pressure [Pa]
    :type  pfull: 1D array
    :param p_sigma: the pressure level where the vertical grid starts
        transitioning from sigma to pressure
    :type  p_sigma: float
    :param p_press: the pressure level above which the vertical grid is
        pure (constant) pressure
    :type  p_press: float
    :return: the transition factor. = 1 for pure sigma, = 0 for pure
        pressure and =0-1 for the transition
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
            t[k] = (np.sin(0.5 * np.pi * x/xx))**2
    return t


def swinbank(plev, psfc, ptrans=1.):
    """
    Compute ``ak`` and ``bk`` values with a transition based on Swinbank

    :param plev: the pressure levels [Pa]
    :type  plev: 1D array
    :param psfc: the surface pressure [Pa]
    :type  psfc: 1D array
    :param ptrans: the transition pressure [Pa]
    :type  ptrans: 1D array
    :return: the coefficients for the new layers
    """

    # ``ks`` = number of pure pressure levels
    ktrans = np.argmin(np.abs(plev - ptrans))
    km = len(plev)-1

    aknew = np.zeros(len(plev))
    bknew = np.zeros(len(plev))

    # ``pnorm`` = 1.e5;
    pnorm = psfc
    eta = plev / pnorm

    # ``ks`` = number of pure pressure levels
    ep = eta[ktrans+1]
    es = eta[-1]
    rnorm = 1. / (es-ep)**2

    # Compute ``alpha``, ``beta``, and ``gamma`` using Swinbank formula
    alpha = (ep**2 - 2.*ep*es) / (es-ep)**2
    beta = 2. * ep*es**2 / (es-ep)**2
    gamma = -(ep*es)**2 / (es-ep)**2

    # Pure Pressure levels
    aknew = eta * pnorm

    # Hybrid pressure-sigma levels
    kdex = range(ktrans+1, km)
    aknew[kdex] = alpha*eta[kdex] + beta + gamma/eta[kdex]
    aknew[kdex] = aknew[kdex] * pnorm
    aknew[-1] = 0.0

    bknew[kdex] = (plev[kdex] - aknew[kdex])/psfc
    bknew[-1] = 1.0

    # Find the transition level ``ks`` where ``bk[ks]>0``
    ks = 0
    while bknew[ks] == 0.:
        ks += 1
    # ``ks`` would be used for fortran indexing in ``fv_eta.f90``
    return aknew, bknew, ks


def polar_warming(T, lat, outside_range=np.nan):
    """
    Return the polar warming, following McDunn et al. 2013:
    Characterization of middle-atmosphere polar warming at Mars, JGR
    Alex Kling

    :param T: temperature with the lat dimension FIRST (transpose as
        needed)
    :type  T: ND array
    :param lat: latitude array
    :type  lat: 1D array
    :param outside_range: values to set the polar warming to when
        outside pf the range. Default = NaN but 0 may be desirable.
    :type  outside_range: float
    :return: The polar warming [K]

    .. note::
        ``polar_warming()`` concatenates the results from both
        hemispheres obtained from the nested function
        ``PW_half_hemisphere()``
    """

    def PW_half_hemisphere(T_half, lat_half, outside_range=np.nan):

        # Easy case, T is a 1D on the latitude direction only
        if len(T_half.shape) == 1:
            imin = np.argmin(T_half)
            imax = np.argmax(T_half)

            # Note that we compute polar warming at ALL latitudes and
            # then set NaN the latitudes outside the desired range.
            # We test on the absolute values (``np.abs``) of the
            # latitudes, therefore the function is usable on both
            # hemispheres
            DT_PW_half = T_half - T_half[imin]
            exclude = np.append(np.where(np.abs(lat_half)
                                         - np.abs(lat_half[imin]) < 0),
                                np.where(np.abs(lat_half[imax])
                                         - np.abs(lat_half) < 0))
            # set to NaN
            DT_PW_half[exclude] = outside_range
            return DT_PW_half

        else:
            # General case for N dimensions
            # Flatten all dimension but the first (lat)
            arr_flat = T_half.reshape([T_half.shape[0],
                                       np.prod(T_half.shape[1:])])
            LAT_HALF = np.repeat(lat_half[:, np.newaxis],
                                 arr_flat.shape[1],
                                 axis = 1)

            imin = np.argmin(arr_flat, axis = 0)
            imax = np.argmax(arr_flat, axis = 0)

            # Initialize four working arrays
            tmin0, tmax0, latmin0, latmax0 = [
                np.zeros_like(arr_flat) for _ in range(4)
            ]

            # Get the min/max temperature and latitudes
            for i in range(0, arr_flat.shape[1]):
                tmax0[:, i] = arr_flat[imax[i], i]
                tmin0[:, i] = arr_flat[imin[i], i]
                latmin0[:, i] = lat_half[imin[i]]
                latmax0[:, i] = lat_half[imax[i]]

            # Compute polar warming for that hemisphere
            DT_PW_half = arr_flat - tmin0

            # Set to NaN values outside the range
            tuple_lower_than_latmin = np.where(np.abs(LAT_HALF)
                                               - np.abs(latmin0) < 0)
            tuple_larger_than_latmax = np.where(np.abs(latmax0)
                                                - np.abs(LAT_HALF) < 0)

            DT_PW_half[tuple_lower_than_latmin] = outside_range
            DT_PW_half[tuple_larger_than_latmax] = outside_range
            return DT_PW_half.reshape(T_half.shape)

    # Actual calculations for both hemispheres
    T_SH = T[0:len(lat)//2]
    lat_SH = lat[0:len(lat)//2]

    T_NH = T[len(lat)//2:]
    lat_NH = lat[len(lat)//2:]

    return (np.concatenate((PW_half_hemisphere(T_SH, lat_SH, outside_range),
                            PW_half_hemisphere(T_NH, lat_NH, outside_range)),
                           axis = 0))


def time_shift_calc(var_in, lon, tod, target_times=None):
    """
    Conversion to uniform local time.

    Mars rotates approx. 14.6° lon per Mars-hour (360° ÷ 24.6 hr)
    Each 14.6° shift in lon represents a 1-hour shift in local time
    This code uses the more precise calculation: lon_shift = 24.0 * lon / 360.

    :param var_in: variable to be shifted. Assume ``lon`` is the first
        dimension and ``tod`` is the last dimension
    :type  var_in: ND array
    :param lon: longitude
    :type  lon: 1D array
    :param tod: ``time_of_day`` index from the input file
    :type  tod: 1D array
    :param target_times: local time(s) [hr] to shift to (e.g., ``"3. 15."``)
    :type  target_times: float (optional)
    :return: the array shifted to uniform local time

    .. note::
        If ``target_times`` is not specified, the file is interpolated
        on the same ``tod`` as the input
    """

    if np.shape(var_in) == len(var_in):
        print('Need longitude and time dimensions')
        return

    # Get dimensions of var_in
    dims_in = np.shape(var_in)
    n_dims_in = len(dims_in) - 1

    # Number of longitudes in file
    n_lon = dims_in[0]

    # Number of timesteps per day in input
    n_tod_in = len(tod)
    if n_tod_in == 0:
        print('No time steps in input (time_shift_calc in FV3_utils.py)')
        exit()

    # Store as float for calculations but keep integer version for reshaping
    n_tod_in_float = float(n_tod_in)

    tod = np.squeeze(tod)

    # Array dimensions for output
    if target_times is None:
        # Time shift all local times
        n_tod_out = n_tod_in
    else:
        n_tod_out = len(target_times)

    # Assuming ``time`` is the last dimension, check if it is a local
    # time ``target_times``. If not, reshape the array into
    # ``[..., days, local time]``
    if dims_in[n_dims_in] != n_tod_in:
        n_days = dims_in[n_dims_in] // n_tod_in  # Integer division
        if (n_days * n_tod_in) != dims_in[n_dims_in]:
            print("Time dimensions do not conform")
            return

        # Fix the incorrect indexing
        var_in = np.reshape(var_in, (dims_in[0], dims_in[n_dims_in - 1], n_tod_in, n_days))
        dims_out = np.linspace(len(dims_in) + 1, dtype=np.int32)
        dims_out[len(dims_in)-1] = len(dims_in)
        dims_out[len(dims_in)] = len(dims_in)-1
        var_in = np.transpose(var_in, dims_out)

    # Get new dims_in of var_in if reshaped
    dims_in = np.shape(var_in)

    if len(dims_in) > 2:
        # Assuming lon is the first dimension and time is the last
        # dimension, we need to reshape the array
        # into ``[lon, time, ...]``. The ``...`` dimension is
        # assumed to be the same size as the lon dimension.
        # If there are more than 2 dimensions, we need to
        # reshape the array into ``[lon, combined_dims, time]``
        # where ``combined_dims`` is the product of all dimensions
        # except first (lon) and last (time) = total # data points for
        # each longitude-time combination
        combined_dims = int(np.prod(dims_in[1:len(dims_in)-1]))  # Ensure integer
    else:
        combined_dims = 1

    # Use integer n_tod_in for reshaping
    var_in = np.reshape(var_in, (n_lon, combined_dims, n_tod_in))

    # Create output array
    var_out = np.zeros((n_lon, combined_dims, n_tod_out))

    # Time increment of input data (in hours)
    dt_in = 24.0 / n_tod_in_float  # Use float version for calculations

    # Time increment of output
    if target_times is None:
        # Preserve original time sampling pattern (in hours) but shift
        # it for each lon so # timesteps in output = # timesteps in input
        dt_out = dt_in
    else:
        # Interpolate to all local times
        dt_out = 1.

    # Calculate interpolation indices
    # Convert east longitude to equivalent hours
    lon_shift = 24.0 * lon / 360.
    kk = np.where(lon_shift < 0)
    lon_shift[kk] = lon_shift[kk] + 24.

    fraction = np.zeros((n_lon, n_tod_out))
    lower_indices = np.zeros((n_lon, n_tod_out))
    upper_indices = np.zeros((n_lon, n_tod_out))

    # Core calculation
    # target_times[n]: The target local Mars time we want (e.g., 15:00 local Mars time)
    # lon_shift: The offset in Mars-hours due to Martian longitude
    # The result dtt (delta time transform) tells us which time indices
    # in the original Mars data to interpolate between
    for n in range(n_tod_out):
        # dtt = n*dt_out - lon_shift - target_times[0] + dt_in
        if target_times is None:
            dtt = n*dt_out - lon_shift - tod[0] + dt_in
        else:
            # For specifying target local times
            # ``time_out - xfshif - tod[0] + hrs/stpe`` in input
            dtt = target_times[n] - lon_shift

        # Ensure that local time is bounded by [0, 24] hours
        kk = np.where(dtt < 0.)
        # dtt: time in OG data corresponding to time we want
        dtt[kk] = dtt[kk] + 24.

        # This is the index into the data aray
        lower_idx = np.floor(dtt/dt_in) # time step before target local time
        fraction[:, n] = dtt - lower_idx*dt_in
        kk = np.where(lower_idx < 0.)
        lower_idx[kk] = lower_idx[kk] + n_tod_in_float  # Use float version

        upper_idx = lower_idx + 1. # time step after target local time
        kk = np.where(upper_idx >= n_tod_in_float)  # Use float version
        upper_idx[kk] = upper_idx[kk] - n_tod_in_float  # Use float version

        # Store lower_idx and upper_idx for each lon point and output time
        lower_indices[:, n] = lower_idx[:]
        upper_indices[:, n] = upper_idx[:]

    # Assume uniform time between input data samples
    fraction = fraction / dt_in

    # Now carry out the interpolation
    for n in range(n_tod_out):
        # Number of output time levels
        for i in range(n_lon):
            # Number of longitudes
            lower_idx = np.int32(lower_indices[i, n]) % n_tod_in  # Use modulo with integer n_tod_in
            upper_idx = np.int32(upper_indices[i, n])
            frac = fraction[i, n]
            # Interpolate between the two time levels
            var_out[i, :, n] = (
                (1.-frac) * var_in[i, :, lower_idx]
                + frac * var_in[i, :, upper_idx]
                )

    var_out = np.squeeze(var_out)
    dims_out = np.zeros(len(dims_in), dtype=int)
    for d in range(n_dims_in):
        dims_out[d] = dims_in[d]
    dims_out[n_dims_in] = n_tod_out
    var_out = np.reshape(var_out, dims_out)
    return var_out


def lin_interp(X_in, X_ref, Y_ref):
    """
    Simple linear interpolation with no dependance on scipy

    :param X_in: input values
    :type  X_in: float or array
    :param X_ref x values
    :type  X_ref: array
    :param Y_ref y values
    :type  Y_ref: array
    :return: y value linearly interpolated at ``X_in``
    """

    X_ref = np.array(X_ref)
    Y_ref = np.array(Y_ref)

    # Definition of the interpolating function
    def lin_oneElement(x, X_ref, Y_ref):
        if x < X_ref.min() or x > X_ref.max():
            return np.nan

        # Find closest left-hand size index
        n = np.argmin(np.abs(x-X_ref))
        if X_ref[n] > x or n == len(X_ref):
            n -= 1
        a = (Y_ref[n+1] - Y_ref[n])/(X_ref[n+1] - X_ref[n])
        b = Y_ref[n] - a*X_ref[n]
        return a*x + b

    # Wrapper to the function above
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
    Add a cyclic (overlapping) point to a 2D array. Useful for azimuth
    and orthographic projections.

    :param data: variable of size ``[nlat, nlon]``
    :type  data: array
    :param lon: longitudes
    :type  lon: array
    :return: a 2D array of size ``[nlat, nlon+1]`` with last column
        identical to the 1st; and a 1D array of longitudes
        size [nlon+1] where the last element is ``lon[-1] + dlon``
    """

    # Compute increment
    dlon = lon[1]-lon[0]
    # Create new array, size ``[nlon + 1]``
    data_c = np.zeros((data.shape[0], data.shape[1]+1), float)
    data_c[:, 0:-1] = data[:, :]
    data_c[:, -1] = data[:, 0]
    return data_c, np.append(lon, lon[-1] + dlon)


def spherical_div(U, V, lon_deg, lat_deg, R=3400*1000., spacing="varying"):
    """
    Compute the divergence of the wind fields using finite difference::

        div = du/dx + dv/dy
        -> = 1/(r cos lat)[d(u)/dlon + d(v cos lat)/dlat]

    :param U: wind field with ``lat`` SECOND TO LAST and ``lon`` as last
        dimensions (e.g., ``[lat, lon]`` or ``[time, lev, lat, lon``]
        etc.)
    :type  U: array
    :param V: wind field with ``lat`` SECOND TO LAST and ``lon`` as last
        dimensions (e.g., ``[lat, lon]`` or ``[time, lev, lat, lon``]
        etc.)
    :type  V: array
    :param lon_deg: longitude [°] (2D if irregularly-spaced)
    :type  lon_deg: 1D array
    :param lat_deg: latitude [°] (2D if irregularly-spaced)
    :type  lat_deg: 1D array
    :param R: planetary radius [m]
    :type  R: float
    :param spacing: when ``lon`` and ``lat`` are 1D arrays, using
        spacing = "varying" differentiates latitude and longitude. When
        spacing = "regular", ``dx = lon[1]-lon[0]``,
        `` dy=lat[1]-lat[0]`` and the ``numpy.gradient()`` method are
        used
    :type  spacing: str (defaults to "varying")
    :return: the horizonal divergence of the wind field [m-1]
    """

    lon = lon_deg * np.pi/180
    lat = lat_deg * np.pi/180
    var_shape = U.shape

    # Transpose shapes:
    T_array = np.arange(len(U.shape))
    # One permutation only: ``lon`` is passsed to the 1st dimension
    T_lonIN = np.append(T_array[-1], T_array[0:-1])
    # One permutation only: ``lon`` is passsed to the 1st dimension
    T_lonOUT = np.append(T_array[1:], T_array[0])
    T_latIN = np.append(np.append(T_array[-2], T_array[0:-2]), T_array[-1])
    T_latOUT = np.append(np.append(T_array[1:-1], T_array[0]), T_array[-1])

    if len(lon.shape) == 1:
        # ``lon`` and ``lat`` are 1D arrays
        # Extend broadcasting dimensions for the latitude (e.g.,
        # ``[1, 1, lat, 1]`` if ``U`` is size ``[time, lev, lat, lon]``)
        reshape_shape = [1 for i in range(0, len(var_shape))]
        reshape_shape[-2] = lat.shape[0]
        lat_b = lat.reshape(reshape_shape)
        if spacing == 'regular':
            out = (1 / (R*np.cos(lat_b))
                   * (np.gradient(U, axis = -1) / (lon[1]-lon[0])
                      + np.gradient(V*np.cos(lat_b), axis = -2)
                      / (lat[1]-lat[0])))
        else:
            out = (1 / (R*np.cos(lat_b))
                   * (dvar_dh(U.transpose(T_lonIN),
                              lon).transpose(T_lonOUT)
                      + dvar_dh((V * np.cos(lat_b)).transpose(T_latIN),
                                lat).transpose(T_latOUT)))

    else:
        # ``lon`` and ``lat`` are 2D arrays
        # If ``U`` is ``[time, lev, lat, lon]``, reshape ``lat`` and
        # ``lon`` to ``[time, lev, lat, lon]``
        if var_shape != lon.shape:
            # ``[time, lev, lat, lon] > [time, lev]`` and reverse, so
            # first ``lev``, then ``time``
            for ni in var_shape[:-2][::-1]:
                lat = np.repeat(lat[np.newaxis, ...], ni, axis = 0)
                lon = np.repeat(lon[np.newaxis, ...], ni, axis = 0)

        out = (1 / (R*np.cos(lat))
               * (dvar_dh(U.transpose(T_lonIN),
                          lon.transpose(T_lonIN)).transpose(T_lonOUT)
                  + dvar_dh((V*np.cos(lat)).transpose(T_latIN),
                            lat.transpose(T_latIN)).transpose(T_latOUT)))
    return out


def spherical_curl(U, V, lon_deg, lat_deg, R=3400*1000., spacing="varying"):
    """
    Compute the vertical component of the relative vorticity using
    finite difference::

        curl = dv/dx -du/dy
             = 1/(r cos lat)[d(v)/dlon + d(u(cos lat)/dlat]

    :param U: wind field with ``lat`` SECOND TO LAST and ``lon`` as last
        dimensions (e.g., ``[lat, lon]`` or ``[time, lev, lat, lon``]
        etc.)
    :type  U: array
    :param V: wind field with ``lat`` SECOND TO LAST and ``lon`` as last
        dimensions (e.g., ``[lat, lon]`` or ``[time, lev, lat, lon``]
        etc.)
    :type  V: array
    :param lon_deg: longitude [°] (2D if irregularly-spaced)
    :type  lon_deg: 1D array
    :param lat_deg: latitude [°] (2D if irregularly-spaced)
    :type  lat_deg: 1D array
    :param R: planetary radius [m]
    :type  R: float
    :param spacing: when ``lon`` and ``lat`` are 1D arrays, using
        spacing = "varying" differentiates latitude and longitude. When
        spacing = "regular", ``dx = lon[1]-lon[0]``,
        `` dy=lat[1]-lat[0]`` and the ``numpy.gradient()`` method are
        used
    :type  spacing: str (defaults to "varying")
    :return: the vorticity of the wind field [m-1]
    """

    lon = lon_deg * np.pi/180
    lat = lat_deg * np.pi/180

    var_shape = U.shape

    # Transpose shapes:
    T_array = np.arange(len(U.shape))
    # One permutation only: ``lon`` is passsed to the 1st dimension
    T_lonIN = np.append(T_array[-1], T_array[0:-1])
    # One permutation only: ``lon`` is passsed to the 1st dimension
    T_lonOUT = np.append(T_array[1:], T_array[0])
    T_latIN = np.append(np.append(T_array[-2], T_array[0:-2]), T_array[-1])
    T_latOUT = np.append(np.append(T_array[1:-1], T_array[0]), T_array[-1])

    if len(lon.shape) == 1:
        # lon, lat are 1D arrays
        # Extend broadcasting dimensions for the latitude (e.g.,
        # ``[1 ,1, lat, 1]`` if ``U`` is size ``[time, lev, lat, lon``)
        reshape_shape = [1 for i in range(0, len(var_shape))]
        reshape_shape[-2] = lat.shape[0]
        lat_b = lat.reshape(reshape_shape)
        if spacing == "regular":
            out = (1 / (R*np.cos(lat_b))
            * (np.gradient(V, axis = -1) / (lon[1]-lon[0])
               - np.gradient(U*np.cos(lat_b), axis = -2)/(lat[1] - lat[0])))
        else:
            out = (1 / (R*np.cos(lat_b))
            * (dvar_dh(V.transpose(T_lonIN),
                       lon).transpose(T_lonOUT)
               - dvar_dh((U*np.cos(lat_b)).transpose(T_latIN),
                         lat).transpose(T_latOUT)))

    else:
        # ``lon``, ``lat`` are 2D arrays
        # If ``U`` is ``[time, lev, lat, lon]``, reshape ``lat`` and
        # ``lon`` to ``[time, lev, lat, lon]``)
        if var_shape != lon.shape:
            # ``[time, lev, lat, lon]`` > ``[time, lev]`` and reverse,
            # so first ``lev`` then ``time``
            for ni in var_shape[:-2][::-1]:
                lat = np.repeat(lat[np.newaxis, ...], ni, axis = 0)
                lon = np.repeat(lon[np.newaxis, ...], ni, axis = 0)

        out = (1 / (R*np.cos(lat))
               * (dvar_dh(V.transpose(T_lonIN),
                          lon.transpose(T_lonIN)).transpose(T_lonOUT)
                  - dvar_dh((U*np.cos(lat)).transpose(T_latIN),
                            lat.transpose(T_latIN)).transpose(T_latOUT)))
    return out


def frontogenesis(U, V, theta, lon_deg, lat_deg, R=3400*1000.,
                  spacing="varying"):
    """
    Compute the frontogenesis (local change in potential temperature
    gradient near a front) following Richter et al. 2010: Toward a
    Physically Based Gravity Wave Source Parameterization in a General
    Circulation Model, JAS 67.

    We have ``Fn = 1/2 D(Del Theta)^2/Dt`` [K/m/s]

    :param U: wind field with ``lat`` SECOND TO LAST and ``lon`` as last
        dimensions (e.g., ``[lat, lon]`` or ``[time, lev, lat, lon``]
        etc.)
    :type  U: array
    :param V: wind field with ``lat`` SECOND TO LAST and ``lon`` as last
        dimensions (e.g., ``[lat, lon]`` or ``[time, lev, lat, lon``]
        etc.)
    :type  V: array
    :param theta: potential temperature [K]
    :type  theta: array
    :param lon_deg: longitude [°] (2D if irregularly-spaced)
    :type  lon_deg: 1D array
    :param lat_deg: latitude [°] (2D if irregularly-spaced)
    :type  lat_deg: 1D array
    :param R: planetary radius [m]
    :type  R: float
    :param spacing: when ``lon`` and ``lat`` are 1D arrays, using
        spacing = "varying" differentiates latitude and longitude. When
        spacing = "regular", ``dx = lon[1]-lon[0]``,
        `` dy=lat[1]-lat[0]`` and the ``numpy.gradient()`` method are
        used
    :type  spacing: str (defaults to "varying")
    :return: the frontogenesis field [m-1]
    """

    lon = lon_deg*np.pi/180
    lat = lat_deg*np.pi/180

    var_shape = U.shape

    # Transpose shapes:
    T_array = np.arange(len(U.shape))
    # One permutation only: ``lon`` is passsed to the 1st dimension
    T_lonIN = np.append(T_array[-1], T_array[0:-1])
    # One permutation only: ``lon`` is passsed to the 1st dimension
    T_lonOUT = np.append(T_array[1:], T_array[0])
    T_latIN = np.append(np.append(T_array[-2], T_array[0:-2]), T_array[-1])
    T_latOUT = np.append(np.append(T_array[1:-1], T_array[0]), T_array[-1])

    if len(lon.shape) == 1:
        # lon, lat are 1D arrays
        # Extend broadcasting dimensions for the colatitude
        # (e.g., ``[1 ,1, lat, 1]`` if ``U`` is size
        # ``[time, lev, lat, lon``)
        reshape_shape = [1 for i in range(0, len(var_shape))]
        reshape_shape[-2] = lat.shape[0]
        lat_b = lat.reshape(reshape_shape)
        if spacing == "regular":

            du_dlon = np.gradient(U, axis = -1)/(lon[1] - lon[0])
            dv_dlon = np.gradient(V, axis = -1)/(lon[1] - lon[0])
            dtheta_dlon = np.gradient(theta, axis = -1)/(lon[1] - lon[0])

            du_dlat = np.gradient(U, axis = -2)/(lat[1] - lat[0])
            dv_dlat = np.gradient(V, axis = -2)/(lat[1] - lat[0])
            dtheta_dlat = np.gradient(theta, axis = -2)/(lat[1] - lat[0])
        else:
            du_dlon = dvar_dh(U.transpose(T_lonIN), lon).transpose(T_lonOUT)
            dv_dlon = dvar_dh(V.transpose(T_lonIN), lon).transpose(T_lonOUT)
            dtheta_dlon = dvar_dh(theta.transpose(T_lonIN),
                                  lon).transpose(T_lonOUT)

            du_dlat = dvar_dh(U.transpose(T_latIN), lat).transpose(T_latOUT)
            dv_dlat = dvar_dh(V.transpose(T_latIN), lat).transpose(T_latOUT)
            dtheta_dlat = dvar_dh(theta.transpose(T_latIN),
                                  lat).transpose(T_latOUT)
    else:
        # lon, lat are 2D array
        # If ``U`` is ``[time, lev, lat, lon]``, reshape ``lat`` and
        # ``lon`` to ``[time, lev, lat, lon]``)
        if var_shape != lon.shape:
            lat_b = lat.copy()
            # ``[time, lev, lat, lon]`` > ``[time, lev]`` and reverse,
            # so first ``lev`` then ``time``
            for ni in var_shape[:-2][::-1]:
                lat = np.repeat(lat[np.newaxis, ...], ni, axis = 0)
                lon = np.repeat(lon[np.newaxis, ...], ni, axis = 0)

            du_dlon = (dvar_dh(U.transpose(T_lonIN),
                               lon.transpose(T_lonIN)).transpose(T_lonOUT))
            dv_dlon = (dvar_dh(V.transpose(T_lonIN),
                               lon.transpose(T_lonIN)).transpose(T_lonOUT))
            dtheta_dlon = (dvar_dh(theta.transpose(T_lonIN),
                                   lon.transpose(T_lonIN)).transpose(T_lonOUT))

            du_dlat = (dvar_dh(U.transpose(T_latIN),
                               lat.transpose(T_latIN)).transpose(T_latOUT))
            dv_dlat = (dvar_dh(V.transpose(T_latIN),
                               lat.transpose(T_latIN)).transpose(T_latOUT))
            dtheta_dlat = (dvar_dh(theta.transpose(T_latIN),
                                   lat.transpose(T_latIN)).transpose(T_latOUT))

    out = (-(1 / (R*np.cos(lat_b)) * dtheta_dlon)**2
           * (1 / (R*np.cos(lat_b)) * du_dlon - V*np.tan(lat_b)/R)
           - (1 / R*dtheta_dlat)**2 * (1/R*dv_dlat)
           - (1 / (R*np.cos(lat_b)) * dtheta_dlon) * (1/R*dtheta_dlat)
           * (1 / (R*np.cos(lat_b)) * dv_dlon + 1/R*du_dlat
              + U*np.tan(lat_b)/R))
    return out


def MGSzmax_ls_lat(ls, lat):
    """
    Return the max altitude for the dust from "MGS scenario" from
    Montmessin et al. (2004), Origin and role of water ice clouds in
    the Martian water cycle as inferred from a general circulation model

    :param ls: solar longitude [°]
    :type  ls: array
    :param lat : latitude [°]
    :type  lat: array
    :return: top altitude for the dust [km]
    """

    lat = np.array(lat) * np.pi/180
    ls_p = (np.array(ls)-158) * np.pi/180

    return (60
            + 18*np.sin(ls_p)
            - (32 + 18*np.sin(ls_p)) * np.sin(lat)**4
            - 8*np.sin(ls_p) * np.sin(lat)**5)


def MGStau_ls_lat(ls, lat):
    """
    Return the max altitude for the dust from "MGS scenario" from
    Montmessin et al. (2004), Origin and role of water ice clouds in
    the Martian water cycle as inferred from a general circulation model

    :param ls: solar longitude [°]
    :type  ls: array
    :param lat : latitude [°]
    :type  lat: array
    :return: top altitude for the dust [km]
    """

    lat = np.array(lat)
    ls_p = (np.array(ls)-250) * np.pi/180

    tn = 0.1
    teq = 0.2 + 0.3*np.cos(0.5*ls_p)**14
    ts = 0.1 + 0.4*np.cos(0.5*ls_p)**14

    # We have ``tanh(-x)=-tanh(x)``
    t_north = tn + 0.5*(teq-tn) * (1 + np.tanh(4.5 - lat/10))
    t_south = ts + 0.5*(teq-ts) * (1 + np.tanh(4.5 + lat/10))

    if len(np.atleast_1d(lat)) == 1:
        # One latitude
        if lat >= 0:
            tau = t_north
        else:
            t_south
    else:
        tau = np.zeros_like(lat)
        tau[lat <= 0] = t_south[lat <= 0]
        tau[lat > 0] = t_north[lat > 0]
    return tau


def broadcast(var_1D, shape_out, axis):
    """
    Broadcast a 1D array based on a variable's dimensions

    :param var_1D: variable (e.g., ``lat`` size = 36, or ``time``
        size = 133)
    :type  var_1D: 1D array
    :param shape_out: broadcasting shape (e.g.,
        ``temp.shape = [133, lev, 36, lon]``)
    :type  shape_out: list
    :return: (ND array) broadcasted variables (e.g., size =
        ``[1,36,1,1]`` for ``lat`` or ``[133,1,1,1]`` for ``time``)
    """
    # Special case where var_1D has only one element
    var_1D = np.atleast_1d(var_1D)
    reshape_shape = [1 for i in range(0, len(shape_out))]
    # Reshape, e.g., [28, 1, 1, 1]
    reshape_shape[axis] = len(var_1D)
    return var_1D.reshape(reshape_shape)


def ref_atmosphere_Mars_PTD(Zi):
    """
    Analytical atmospheric model for Martian pressure, temperature, and
    density. Alex Kling, June 2021

    :param Zi: input altitude [m] (must be >= 0)
    :type  Zi: float or 1D array
    :return: tuple of corresponding pressure [Pa], temperature [K],
    and density [kg/m3] floats or arrays

    .. note::
        This model was obtained by fitting globally and annually
        averaged reference temperature profiles derived from the Legacy
        GCM, MCS observations, and Mars Climate Database.

        The temperature fit was constructed using quadratic temperature
        ``T(z) = T0 + gam(z-z0) + a*(z-z0)^2`` over 4 segments (0>57 km,
        57>110 km, 110>120 km and 120>300 km).

        From the ground to 120 km, the pressure is obtained by
        integrating (analytically) the hydrostatic equation:

        ``dp/dz=-g. p/(rT)`` with ``T(z) = T0 + gam(z-z0) + a*(z-z0)^2``

        Above ~120 km, ``P = P0 exp(-(z-z0)g/rT)`` is not a good
        approximation as the fluid is in molecula regime. For those
        altitudes, we provide a fit in the form of
        ``P = P0 exp(-az-bz^2)`` based on diurnal average of the MCD
        database at lat = 0, Ls = 150.
    """


    # Internal Functions
    def alt_to_temp_quad(Z, Z0, T0, gam, a):
        """
        Return the a representative globally and annually averaged
        temperature in the form ``T(z) = a(z-z0)^2 + b(z-z0) + T0``

        :param Z: altitude [m]
        :type  Z: float
        :param Z0: starting altitude [m]
        :type  Z0: float
        :param T0: quadratic coefficient
        :type  T0: float
        :param gam: quadratic coefficient
        :type  gam: float
        :param a: quadratic coefficient
        :type  a: float
        :return: temperature at altitude Z [K]
        """

        return (T0 + gam*(Z-Z0) + a*(Z-Z0)**2)


    def alt_to_press_quad(Z, Z0, P0, T0, gam, a, rgas, g):
        """
        Return the pressure [Pa] in the troposphere as a function of
        height for a quadratic temperature profile.

        ``T(z) = T0 + gam(z-z0) + a*(z-z0)^2``

        :param Z: altitude [m]
        :type  Z: float
        :param Z0: starting altitude [m]
        :type  Z0: float
        :param P0: reference pressure at Z0[Pa]
        :type  P0: float
        :param T0: reference temperature at Z0[Pa]
        :type  T0: float
        :param gam: linear and quadratic coeff for the temperature
        :type  gam: float
        :param a: linear and quadratic coeff for the temperature
        :type  a: float
        :param rgas: specific gas constant [J/kg/K]
        :type  rgas: float
        :param rg: gravity [m/s2]
        :type  rg: float
        :return: pressure at alitude Z [Pa]
        """

        delta = (4*a*T0 - gam**2)

        if delta >= 0:
            sq_delta = np.sqrt(4*a*T0 - gam**2)
            return (P0
                    *np.exp(-2*g / (rgas*sq_delta)
                            * (np.arctan((2*a*(Z-Z0)+gam) / sq_delta)
                               - np.arctan(gam/sq_delta))))
        else:
            delta = -delta
            sq_delta = np.sqrt(delta)
            return (P0 * (((2*a*(Z-Z0) + gam) - sq_delta)
                       * (gam + sq_delta)
                       / (((2*a*(Z-Z0)+gam) + sq_delta) * (gam-sq_delta)))
                    **(-g / (rgas*sq_delta)))


    def P_mars_120_300(Zi, Z0=120000., P0=0.00012043158397922564,
                       p1=1.09019694e-04, p2=-3.37385416e-10):
        """
        Martian pressures from 120-300 km. Above ~120 km,
        ``P = P0 exp(-(z-z0)g/rT)`` is not a good approximation as the
        fluid is in a molecular regime. Alex Kling

        Fit from a diurnal average of the MCD database at ``lat = 0``,
        ``Ls = 150``.

        To be consistent with Earth physics, we use
        ``P = P0 exp(-az - bz^2 - cz^c-d^4 ...)``

        :param Z: altitude [m]
        :type  Z: float
        :return: the equivalent pressure at altitude [Pa]
        """

        return (P0 * np.exp(-p1*(Zi-Z0) - p2*(Zi-Z0)**2))


    def T_analytic_scalar(Zi):
        """
        A best fit of globally averaged temperature profiles from
        various sources including: Legacy MGCM, MCS, & MCD
        """

        if Zi <= 57000:
            return alt_to_temp_quad(Zi, Z0 = 0, T0 = 225.9,
                                    gam = -0.00213479, a = 1.44823e-08)
        elif 57000 < Zi <= 110000:
            return alt_to_temp_quad(Zi, Z0 = 57000, T0 = 151.2,
                                    gam = -0.000367444, a = -6.8256e-09)
        elif 110000 < Zi <= 170000:
            return alt_to_temp_quad(Zi, Z0 = 110000, T0 = 112.6,
                                    gam = 0.00212537, a = -1.81922e-08)
        elif 170000 <= Zi:
            return 174.6


    def P_analytic_scalar(Zi):
        """
        Analytic solution for a pressure derived from a temperature
        profile
        """

        if Zi <= 57000:
            return alt_to_press_quad(Zi, Z0 = 0, P0 = 610, T0 = 225.9,
                                     gam = -0.00213479, a = 1.44823e-08,
                                     rgas = 192, g = 3.72)
        elif 57000 < Zi <= 110000:
            return alt_to_press_quad(Zi, Z0 = 57000, P0 = 1.2415639872674782,
                                     T0 = 151.2, gam = -0.000367444,
                                     a = -6.8256e-09, rgas = 192, g = 3.72)
        elif 110000 < Zi <= 120000:
            # Discarded above 120 km when we enter the molecular regime
            return alt_to_press_quad(Zi, Z0 = 110000,
                                     P0 = 0.0005866878792825923, T0 = 112.6,
                                     gam = 0.00212537, a = -1.81922e-08,
                                     rgas = 192, g = 3.72)
        elif 120000 <= Zi:
            return P_mars_120_300(Zi)


    if len(np.atleast_1d(Zi)) > 1:
        # Vectorize array
        P_analytic_scalar = np.vectorize(P_analytic_scalar)
        T_analytic_scalar = np.vectorize(T_analytic_scalar)
    return (P_analytic_scalar(Zi), T_analytic_scalar(Zi),
            P_analytic_scalar(Zi) / (192*T_analytic_scalar(Zi)))


def press_to_alt_atmosphere_Mars(Pi):
    """
    Return the altitude [m] as a function of pressure from the
    analytical calculation above.

    :param Pi: input pressure [Pa] (<= 610 Pa)
    :type  Pi: float or 1D array
    :return: the corresponding altitude [m] (float or 1D array)
    """

    # Internal Functions
    def press_to_alt_quad(P, Z0, P0, T0, gam, a, rgas, g):
        """
        Return the altitude [m] as a function of pressure for a
        quadratic temperature profile::

            T(z) = T0 + gam(z-z0) + a*(z-z0)^2

        :param P: pressure [Pa]
        :type  P: float

        :param Z0: referecence altitude [m]
        :type  Z0: float
        :param P0: reference pressure at Z0[Pa]
        :type  P0: float
        :param T0: reference temperature at Z0[Pa]
        :type  T0: float
        :param gam: linear and quadratic coefficients for temperature
        :type  gam: float
        :param a: linear and quadratic coefficients for temperature
        :type  a: float
        :param rgas: specific gas constant [J/kg/K]
        :type  rgas: float
        :param g: gravity [m/s2]
        :type  g:
        :return: the corresponding altitude [m] (float or 1D array)
        """

        delta = (4*a*T0 - gam**2)

        if delta >= 0:
            sq_delta = np.sqrt(4*a*T0 - gam**2)
            return (Z0
                    + sq_delta/(2*a)
                    * np.tan(np.arctan(gam/sq_delta)
                             + np.log(P0/P) * rgas * sq_delta / (2*g))
                    - gam / (2*a))
        else:
            delta = -delta
            sq_delta = np.sqrt(delta)
            return (Z0
                    + (gam-sq_delta)/(2*a)
                    * (1-(P/P0)**(-rgas * sq_delta/g))
                    / ((gam - sq_delta)
                       / (gam + sq_delta)
                       * (P/P0)**(-rgas * sq_delta/g) - 1))


    def press_to_alt_mars_120_300(P, Z0=120000., P0=0.00012043158397922564,
                                  p1=1.09019694e-04, p2=-3.37385416e-10):
        """
        Martian altitude as a function of pressure from 120-300 km.

        :param P: pressure [m]
        :type  P: float
        :return: altitude [m]
        """

        # ``delta > 0`` on this pressure interval
        delta = (p1**2 - 4*p2*np.log(P/P0))
        return ((-p1 + np.sqrt(delta)) / (2*p2) + Z0)


    def alt_analytic_scalar(Pi):
        """
        Analytic solution for the altitude as a function of pressure.
        """

        if Pi >= 610:
            return 0.
        elif 610 > Pi >= 1.2415639872674782:
            # The pressure from ``alt_to_press_quad`` at 57,000 m
            return press_to_alt_quad(Pi, Z0 = 0,
                                     P0 = 610,
                                     T0 = 225.9,
                                     gam = -0.00213479,
                                     a = 1.44823e-08,
                                     rgas = 192,
                                     g = 3.72)
        elif 1.2415639872674782 > Pi >=  0.0005866878792825923:
            # 57,000-110,000 m
            return press_to_alt_quad(Pi, Z0 = 57000,
                                     P0 = 1.2415639872674782,
                                     T0 = 151.2,
                                     gam = -0.000367444,
                                     a = -6.8256e-09,
                                     rgas = 192,
                                     g = 3.72)
        elif 0.0005866878792825923 > Pi >=  0.00012043158397922564:
            # 110,000-120,000 m
            return press_to_alt_quad(Pi,
                                     Z0 = 110000,
                                     P0 = 0.0005866878792825923,
                                     T0 = 112.6,
                                     gam = 0.00212537,
                                     a = -1.81922e-08,
                                     rgas = 192,
                                     g = 3.72)
        elif 0.00012043158397922564 > Pi:
            # 120,000-300,000 m
            return press_to_alt_mars_120_300(Pi,
                                             Z0 = 120000.,
                                             P0 = 0.00012043158397922564,
                                             p1 = 1.09019694e-04,
                                             p2 = -3.37385416e-10)
    if len(np.atleast_1d(Pi)) > 1:
        alt_analytic_scalar = np.vectorize(alt_analytic_scalar)
    return alt_analytic_scalar(Pi)


# ============================ Projections =============================
# The projections below were implemented by Alex Kling following "An
# Album of Map Projections," USGS Professional Paper 1453, (1994)
# https://pubs.usgs.gov/pp/1453/report.pdf


def azimuth2cart(LAT, LON, lat0, lon0=0):
    """
    Azimuthal equidistant projection. Converts from latitude-longitude
    to cartesian coordinates.

    :param LAT: latitudes[°] size [nlat]
    :type  LAT: 1D or 2D array
    :param LON: longitudes [°] size [nlon]
    :type  LON: 1D or 2D array
    :param lat0: latitude coordinate of the pole
    :type  lat0: float
    :param lon0: longitude coordinate of the pole
    :type  lon0: float
    :return: the cartesian coordinates for the latitudes and longitudes
    """

    # Convert to radians
    LAT = LAT * np.pi/180
    lat0 = lat0 * np.pi/180
    LON = LON * np.pi/180
    lon0 = lon0 * np.pi/180

    c = np.arccos(np.sin(lat0) * np.sin(LAT)
                  + np.cos(lat0) * np.cos(LAT) * np.cos(LON - lon0))
    k = c / np.sin(c)
    X = k * np.cos(LAT) * np.sin(LON - lon0)
    Y = k * (np.cos(lat0) * np.sin(LAT)
             - np.sin(lat0) * np.cos(LAT) * np.cos(LON - lon0))
    return X, Y


def ortho2cart(LAT, LON, lat0, lon0=0):
    """
    Orthographic projection. Converts from latitude-longitude to
    cartesian coordinates.

    :param LAT: latitudes[°] size [nlat]
    :type  LAT: 1D or 2D array
    :param LON: longitudes [°] size [nlon]
    :type  LON: 1D or 2D array
    :param lat0: latitude coordinate of the pole
    :type  lat0: float
    :param lon0: longitude coordinate of the pole
    :type  lon0: float
    :return: the cartesian coordinates for the latitudes and longitudes;
        and a mask (NaN array) that hides the back side of the planet
    """

    # Convert to radians
    LAT = LAT  * np.pi/180
    lat0 = lat0  * np.pi/180
    LON = LON  * np.pi/180
    lon0 = lon0  * np.pi/180
    MASK = np.ones_like(LON)

    X = np.cos(LAT) * np.sin(LON - lon0)
    Y = (np.cos(lat0) * np.sin(LAT)
         - np.sin(lat0) * np.cos(LAT) * np.cos(LON - lon0))

    # Filter values on the opposite side of Mars (i.e., ``cos(c) < 0``)
    cosc = (np.sin(lat0) * np.sin(LAT)
            + np.cos(lat0) * np.cos(LAT) * np.cos(LON - lon0))
    MASK[cosc < 0] = np.nan
    return X, Y, MASK


def mollweide2cart(LAT, LON):
    """
    Mollweide projection. Converts from latitude-longitude to
    cartesian coordinates.

    :param LAT: latitudes[°] size [nlat]
    :type  LAT: 1D or 2D array
    :param LON: longitudes [°] size [nlon]
    :type  LON: 1D or 2D array
    :param lat0: latitude coordinate of the pole
    :type  lat0: float
    :param lon0: longitude coordinate of the pole
    :type  lon0: float
    :return: the cartesian coordinates for the latitudes and longitudes
    """

    # Convert to radians
    LAT = np.array(LAT) * np.pi/180
    LON = np.array(LON) * np.pi/180
    lon0 = 0

    def compute_theta(lat):
        """
        Internal Function to compute theta. Latitude is in radians here.
        """

        theta0 = lat
        sum = 0
        running = True

        while running and sum <= 100:
            # Solve for theta using Newton–Raphson
            theta1 = (theta0
                      - (2*theta0 + np.sin(2*theta0) - np.pi*np.sin(lat))
                      / (2 + 2*np.cos(2*theta0)))
            sum += 1
            if np.abs((theta1-theta0)) < 10**(-3):
                running = False
            theta0 = theta1
        if sum == 100:
            print("Warning in ``mollweide2cart()``: Reached maximum "
                  "number of iterations")
        return theta1


    if len(np.atleast_1d(LAT).shape) == 1:
        # Float or 1D array
        nlat = len(np.atleast_1d(LAT))
        LAT = LAT.reshape((nlat))
        LON = LON.reshape((nlat))
        THETA = np.zeros((nlat))
        for i in range(0, nlat):
            THETA[i] = compute_theta(LAT[i])
    else:
        # 2D array
        nlat = LAT.shape[0]
        nlon = LAT.shape[1]
        theta = np.zeros((nlat))
        for i in range(0, nlat):
            theta[i] = compute_theta(LAT[i, 0])
        THETA = np.repeat(theta[:, np.newaxis], nlon, axis = 1)

    X = 2 * np.sqrt(2) / np.pi * (LON - lon0) * np.cos(THETA)
    Y = np.sqrt(2) * np.sin(THETA)
    return np.squeeze(X), np.squeeze(Y)


def robin2cart(LAT, LON):
    """
    Robinson projection. Converts from latitude-longitude to cartesian
    coordinates.

    :param LAT: latitudes[°] size [nlat]
    :type  LAT: 1D or 2D array
    :param LON: longitudes [°] size [nlon]
    :type  LON: 1D or 2D array
    :param lat0: latitude coordinate of the pole
    :type  lat0: float
    :param lon0: longitude coordinate of the pole
    :type  lon0: float
    :return: the cartesian coordinates for the latitudes and longitudes
    """

    # Convert to radians
    lon0 = 0.
    LAT = np.array(LAT) * np.pi/180
    LON = np.array(LON) * np.pi/180

    lat_ref = np.array([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65,
                        70, 75, 80, 85, 90.]) * np.pi/180
    x_ref = np.array([1.0000, 0.9986, 0.9954, 0.9900, 0.9822, 0.9730, 0.9600,
                      0.9427, 0.9216, 0.8962, 0.8679, 0.8350, 0.7986, 0.7597,
                      0.7186, 0.6732, 0.6213, 0.5722, 0.5322])
    y_ref = np.array([0.0000, 0.0620, 0.1240, 0.1860, 0.2480, 0.3100, 0.3720,
                      0.4340, 0.4958, 0.5571, 0.6176, 0.6769, 0.7346, 0.7903,
                      0.8435, 0.8936, 0.9394, 0.9761, 1.0000])

    if len(np.atleast_1d(LAT).shape) == 1:
        # Float or 1D array
        X1 = lin_interp(np.abs(LAT), lat_ref, x_ref)
        Y1 = np.sign(LAT) * lin_interp(np.abs(LAT), lat_ref, y_ref)
    else:
        # 2D array
        nlat = LAT.shape[0]
        nlon = LAT.shape[1]
        lat = LAT[:, 0]
        x1 = lin_interp(np.abs(lat), lat_ref, x_ref)
        y1 = np.sign(lat)*lin_interp(np.abs(lat), lat_ref, y_ref)

        X1 = np.repeat(x1[:, np.newaxis], nlon, axis = 1)
        Y1 = np.repeat(y1[:, np.newaxis], nlon, axis = 1)

    X = 0.8487 * X1 * (LON - lon0)
    Y = 1.3523 * Y1
    return X, Y


# ======================= End Projections Section ======================


def sol2ls(jld, cumulative=False):
    """
    Return the solar longitude (Ls) as a function of the sol number.
    Sol=0 is the spring equinox.

    :param jld: sol number after perihelion
    :type  jld: float or 1D array
    :param cumulative: if True, result is cumulative
        (Ls=0-360, 360-720 etc..)
    :type  cumulative: bool
    :return: the corresponding solar longitude
    """

    # Constants
    # Year in sols
    year = 668.0
    # Date of perihelion
    zero_date = 488.
    # Date of northern equinox for a 668-sol year
    equinox = 180
    small_value = 1.0e-7
    pi = np.pi
    degrad = pi/180.0

    # If ``jld`` is a scalar, reshape to a 1-element array
    jld = np.array(jld).astype(float).reshape(len(np.atleast_1d(jld)))


    # ==================================================================
    # Internal Function 1: Calculate Ls 0-360 using a numerical solver
    # ==================================================================
    def sol2ls_mod(jld):
        """
        Based on Tanguy's ``aerols.py``. Useful link:
        http://www.jgiesen.de/kepler/kepler.html
        """

        # Specify orbit eccentricity
        ec = .093
        # Specify angle of planet inclination
        er = ((1.0 + ec) / (1.0 - ec))**0.5

        # Initialize working arrays
        w, rad, als, areols, MY = [np.zeros_like(jld) for _ in range(0, 5)]

        # Days since last perihelion passage
        date = jld - zero_date

        # Determine true anomaly at equinox (``eq1``)
         # ``qq`` is the mean anomaly
        qq = 2.0 * pi * equinox / year
        e = 1.0
        diff = 1.0
        while (diff > small_value):
            ep = e - (e-ec*np.sin(e)-qq) / (1.0-ec*np.cos(e))
            diff = abs(ep - e)
            e = ep

        eq1 = 2.0 * np.arctan(er * np.tan(0.5 * e))

        # Determine true anomaly at current date (``w``)
        for i in range(0, len(jld)):
            e = 1.0
            diff = 1.0
            em = 2. * pi * date[i] / year
            while (diff > small_value):
                ep = e - (e-ec*np.sin(e)-em) / (1.0-ec*np.cos(e))
                diff = abs(ep - e)
                e = ep
            w[i] = 2.0 * np.arctan(er * np.tan(0.5*e))

        # Aerocentric Longitude (``als``)
        als = w - eq1
        areols = als / degrad
        areols[areols < 0.] += 360.
        return areols


    # ==================================================================
    # Internal Function 2: Calculate cumulative Ls 0-720
    # ==================================================================


    def sol2ls_cumu(jld):
        """
        Calculate cumulative Ls. Continuous solar longitudes
        (``Ls_c``=0-359, 360-720...) are obtained by adding 360 to the
        Ls at every equinox based on the sol number. Since ``sol2ls``
        uses a numerical solver, the equinox may return either
        359.9999 or 0.0001. Adding 360 should only be done in the latter
        case to mitigate outlier points.

        For those edge cases where Ls is close to 359.9, the routine
        recalculates the Ls at a later time (say 1 sols) to check for
        outlier points.
        """

        # Calculate cumulative Ls using ``sol2ls`` and adding 360 for
        # every Mars year
        areols, MY = [np.zeros_like(jld) for _ in range(0, 2)]
        date = jld - zero_date
        MY = (date - equinox)//(year) + 1
        Ls_mod = sol2ls_mod(jld)
        Ls_cumu = Ls_mod + MY*360.

        # Check indices where the returned Ls is close to 360. The [0]
        # turns a tuple from ``np.where`` into a list
        index = np.where(Ls_mod >= 359.9)[0]
        for ii in index:
            # Compute Ls one day after (arbitrary length)
            jld_plus1 = jld[ii] + 1.
            Ls_plus1 = sol2ls(jld_plus1)
            date_plus1 = jld_plus1 - zero_date
            MY_plus1 = (date_plus1 - equinox)//(668.) + 1
            # cumulative Ls 1 day after
            Ls_cumu_plus1 = Ls_plus1+MY_plus1 * 360.
            # If things are smooth, the Ls should go from 359-361. If
            # it is 359-721, we need to update the MY for those indices.
            # Difference between two consecutive Ls should be small
            # unless ``Ls_cumu`` was too big in the first place
            diff = Ls_cumu_plus1 - Ls_cumu[ii]
            if diff < 0:
                MY[ii] -= 1

        # Recompute one more time with updated MY
        Ls_cumu = Ls_mod + MY*360.
        return Ls_cumu

    if cumulative:
        return sol2ls_cumu(jld)
    else:
        return sol2ls_mod(jld)


def ls2sol(Ls_in):
    """
    Ls to sol converter.

    :param Ls_in: solar longitudes (0-360...720)
    :type  Ls_in: float or 1D array
    :return: the corresponding sol number

    .. note::
        This function simply uses a numerical solver on the
        ``sol2ls()`` function.
    """

    def internal_func(Ls_in):
        func_int = lambda x: sol2ls(x, cumulative = True)
        # Distance to the target function
        errfunc = lambda x, y: func_int(x) - y
        # Initial guess for the parameters
        p0 = [0.]
        p1, success = optimize.leastsq(errfunc, p0[:], args = (Ls_in))
        return  p1
    Ls_in = np.array(Ls_in)

    if len(np.atleast_1d(Ls_in)) == 1:
        return internal_func(Ls_in)[0]
    else:
        sol_all = []
        for ii in Ls_in:
            sol_all.append(internal_func(ii)[0])
        return   sol_all
