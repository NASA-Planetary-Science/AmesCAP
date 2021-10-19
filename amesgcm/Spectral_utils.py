#===================================================================================
#   This files contains wave analysis routine. Note the dependencies on scipy.signal
#===================================================================================
import numpy as np
from amesgcm.Script_utils import prYellow,prCyan,prRed,prGreen
from amesgcm.FV3_utils import area_weights_deg

try:
    from scipy.signal import butter,filtfilt,detrend
except ImportError as error_msg:
    prYellow("Error while importing modules from scipy.signal")
    exit()
except Exception as exception:
    # Output unexpected Exceptions.
    print(exception.__class__.__name__ + ": ", exception)
    exit()
    
    
# def init_shtools():
#     '''
#     The following code simply loads the pyshtools module and provides adequate referencing.
#     Since dependencies may need to be solved by the user, the module import is wrapped in a function
#     that may be called when needed.
#     '''
#     try:
#         import  pyshtools
#     except ImportError as error_msg:
#         prYellow("__________________")
#         prYellow("Zonal decomposition relies on the pyshtools library, referenced at:")
#         prYellow('')
#         prYellow("Mark A. Wieczorek and Matthias Meschede (2018). SHTools - Tools for working with spherical harmonics, Geochemistry, Geophysics, Geosystems, 19, 2574-2592, doi:10.1029/2018GC007529")
#         prYellow("Please consult installation instructions at:")
#         prCyan("https://pypi.org/project/pyshtools")
#         prYellow("And install with:")
#         prCyan("pip install pyshtools")
#         exit()
#     except Exception as exception:
#         # Output unexpected Exceptions.
#         print(exception.__class__.__name__ + ": ", exception)
#  


def diurn_extract(VAR,N,tod,lon):
    '''
    Extract the diurnal component of a field. Original code by J.Wilson adapted by A. Kling. April, 2021
    Args:
        VAR (1D or ND array)   : field to process with time of day dimension FIRST, e.g (tod,time,lat,lon) or (tod)
        N   (int)              : number of harmonics to extract (N=1 for diurnal,N=2  for diurnal + semi diurnal etc...)
        tod (1D array)         : universal time of day in sols (0>1.) If provided in hours (0>24), it  will  be  normalized.
        lon (1D array or float): longitudes 0>360
    Return:
        amp (ND array)  :  the amplitudes for the Nth first harmonics, e.g. size (Nh,time,lat,lon)    
        phase (ND array):  the phases for the Nth first harmonics, e.g. size (Nh,time,lat,lon) 
    '''
    dimsIN=VAR.shape
    
    nsteps= len(tod)    
    period= 24;   rnorm= 1/nsteps;   delta= period/nsteps;
    
    #   Be sure that the local time grid is in units of days
    if max(tod) > 1:   tod= tod/24.
    
    freq= (delta/period)*2.0*np.pi 
    arg= tod * 2* np.pi 
    arg= arg.reshape([len(tod), 1] ) #reshape array for matrix operations
        
    #Dimensions for the output
    if len(dimsIN) == 1:
        dimsOUT=[N,1] # if VAR is size (tod,time,lat,lon) dimsOUT is size (N,time,lat,lon)
    else:
        dimsOUT=np.append([N],dimsIN[1:])
        
    #Reshape input variable VAR as a 2D array (tod,Nelements) for generalization 
    
    Ndim= np.int(np.prod(dimsIN[1:]))   #Ndim is the product of all dimensions but the time of day axis, e.g. time x lat x lon
    dimsFLAT=np.append([nsteps],[Ndim])  # Shape of flattened array        
    dimsOUT_flat=np.append([N],[Ndim])   # Shape of flattened array    
    VAR= VAR.reshape(dimsFLAT)     #Flatten array to  (tod,Nelements)
    
    #Initialize output arrays
    amp,phas=np.zeros(dimsOUT_flat),np.zeros(dimsOUT_flat) 
    
    #if nargin > 3 (Initial code, we will assume lon is always provided)
    
    if len(dimsIN) == 1:
        corr= np.array([lon])
    else : #Repeat longitude array to match the size of the input VARIABLES, minus the first (tod) axis.
        tilenm= np.append(dimsIN[1:-1],1)  # Create axis to expend the longitude array
        lonND=np.tile(lon,tilenm)      # if VAR is (tod,time,lat,lon) lonN is longitude repeated as (time,lat,lon)
        corr=lonND.flatten()
         
    # Python indexing starts at 0 so we index at nn-1
    for nn in range(1,N+1):
        cosser=   np.dot(VAR[:,...].T ,np.cos(nn*arg )).squeeze()
        sinser=   np.dot(VAR[:,...].T ,np.sin(nn*arg )).squeeze()
    
        amp[nn-1,:]= 2*rnorm*np.sqrt( cosser**2 + sinser**2)
        phas[nn-1,:]= (180/np.pi) * np.arctan2( sinser, cosser)
    
        #Apply local time correction to the phase
        phas[nn-1,:]= phas[nn-1,:] + 360 + nn*corr[:]
        phas[nn-1,:]= (24/(nn)/360) * np.mod( phas[nn-1,:],360 )
         
    
    # Return the phase and amplitude
    return  amp.reshape( dimsOUT), phas.reshape( dimsOUT )
    
    
def reconstruct_diurn(amp,phas,tod,lon,sumList=[]):
    '''
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
    '''
    
    #Initialization
    dimsIN=amp.shape
    N=dimsIN[0]
    dimsSUM=np.append([len(tod)],dimsIN[1:])
    dimAXIS=np.arange(len(dimsSUM)) #Create dimensions array, e.g. [0,1,2,3] for a 4D variables 
  
    dimsOUT=np.append([dimsIN[0],len(tod)],dimsIN[1:])
    varOUT=np.zeros(dimsOUT)
    
    #Special case for station data (lon is a float)
    if len(np.atleast_1d(lon))==1:lon=np.array([lon])
    
    #Reshape  lon array for broadcasting, e.g. lon[96] to  [1,1,1,96]
    dimAXIS=np.arange(len(dimsSUM));dimAXIS[:]=1;dimAXIS[-1]=len(lon)#
    lon=lon.reshape(dimAXIS)
    #Reshape tod array
    dimAXIS=np.arange(len(dimsSUM));dimAXIS[:]=1;dimAXIS[0]=len(tod) 
    tod=tod.reshape(dimAXIS)    

    # Shift in phase due to local time
    DT=lon/360*24
    
    varSUM =np.zeros(dimsSUM)
    for nn in range(1,N+1):
        #Compute each harmonic
        varOUT[nn-1,...]=amp[nn-1,...]*np.cos(nn*(tod-phas[nn-1,...]+DT)/24*2*np.pi)
        # If a sum of harmonics is requested, sum it
        if nn in sumList:varSUM+=varOUT[nn-1,...]
        
    if sumList:
        #Return the agregated harmonics 
        return varSUM
    else:
        #Return all harmonics individually
        return varOUT     
        
def space_time(lon,timex, varIN,kmx,tmx):
    """
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

    """           
    
    dims= varIN.shape             #get input variable dimensions
    
    lon_id= dims[0]    # lon          
    time_id= dims[-1]  # time     
    dim_sup_id=dims[1:-1] #additional dimensions stacked in the middle
    jd= np.int(np.prod( dim_sup_id))     #jd is the total number of dimensions in the middle is varIN>3D
    
    varIN= np.reshape(varIN, (lon_id, jd, time_id) )   #flatten the middle dimensions if any
    
    #Initialize 4 empty arrays
    ampw, ampe,phasew,phasee =[np.zeros((kmx,tmx,jd)) for _x in range(0,4)]
    
    #TODO not implemented yet: zamp,zphas=[np.zeros((jd,tmx)) for _x in range(0,2)] 
    
    tpi= 2*np.pi
    argx= lon * 2*np.pi/360  #nomalize longitude array
    rnorm= 2./len(argx)
    
    arg= timex * 2* np.pi
    rnormt= 2./len(arg)
    
    #
    for kk in range(0,kmx): 
        progress(kk,kmx) 
        cosx= np.cos( kk*argx )*rnorm
        sinx= np.sin( kk*argx )*rnorm
        
    #   Inner product to calculate the Fourier coefficients of the cosine
    #   and sine contributions of the spatial variation
        acoef = np.dot(varIN.T,cosx) 
        bcoef = np.dot(varIN.T,sinx)

    # Now get the cos/sine series expansions of the temporal
    #variations of the acoef and bcoef spatial terms.
        for nn in range(0,tmx):
            cosray= rnormt*np.cos(nn*arg )
            sinray= rnormt*np.sin(nn*arg )
        
            cosA=  np.dot(acoef.T,cosray)
            sinA=  np.dot(acoef.T,sinray)
            cosB=  np.dot(bcoef.T,cosray)
            sinB=  np.dot(bcoef.T,sinray)
        
        
            wr= 0.5*(  cosA - sinB )
            wi= 0.5*( -sinA - cosB )
            er= 0.5*(  cosA + sinB )
            ei= 0.5*(  sinA - cosB )
        
            aw= np.sqrt( wr**2 + wi**2 )
            ae= np.sqrt( er**2 + ei**2)
            pe= np.arctan2(ei,er) * 180/np.pi
            pw= np.arctan2(wi,wr) * 180/np.pi
        
            pe= np.mod( -np.arctan2(ei,er) + tpi, tpi ) * 180/np.pi
            pw= np.mod( -np.arctan2(wi,wr) + tpi, tpi ) * 180/np.pi
        
            ampw[kk,nn,:]= aw.T
            ampe[kk,nn,:]= ae.T
            phasew[kk,nn,:]= pw.T
            phasee[kk,nn,:]= pe.T
    #End loop
    
    
    ampw=   np.reshape( ampw,    (kmx,tmx)+dim_sup_id )
    ampe=   np.reshape( ampe,    (kmx,tmx)+dim_sup_id )
    phasew= np.reshape( phasew,  (kmx,tmx)+dim_sup_id )
    phasee= np.reshape( phasee,  (kmx,tmx)+dim_sup_id )

    #TODO implement zonal mean: zamp,zphas,stamp,stphs
    '''
    #  varIN= reshape( varIN, dims );
    
    #if nargout < 5;  return;  end ---> only  ampe,ampw,phasee,phasew are requested
    
    
    #   Now calculate the axisymmetric tides  zamp,zphas
    
    zvarIN= np.mean(varIN,axis=0)
    zvarIN= np.reshape( zvarIN, (jd, time_id) )
    
    arg= timex * 2* np.pi
    arg= np.reshape( arg, (len(arg), 1 ))
    rnorm= 2/time_id
    
    for nn in range(0,tmx):
        cosray= rnorm*np.cos( nn*arg )
        sinray= rnorm*np.sin( nn*arg )
    
        cosser=  np.dot(zvarIN,cosray)
        sinser=  np.dot(zvarIN,sinray)
    
        zamp[:,nn]= np.sqrt( cosser[:]**2 + sinser[:]**2 ).T
        zphas[:,nn]= np.mod( -np.arctan2( sinser, cosser )+tpi, tpi ).T * 180/np.pi
    
    
    zamp=  zamp.T #np.permute( zamp,  (2 1) )
    zphas= zphas.T #np.permute( zphas, (2,1) )
    
    if len(dims)> 2:
        zamp=  np.reshape( zamp,  (tmx,)+dim_sup_id )
        zamp=  np.reshape( zphas, (tmx,)+dim_sup_id ) 
    
    
    
    #if nargout < 7;  return;  end
    
    sxx= np.mean(varIN,ndims(varIN));
    [stamp,stphs]= amp_phase( sxx, lon, kmx );
    
    if len(dims)> 2;
        stamp= reshape( stamp, [kmx dims(2:end-1)] );
        stphs= reshape( stphs, [kmx dims(2:end-1)] );
    end

    '''
    return ampe,ampw,phasee,phasew
    
    
def zeroPhi_filter(VAR, btype, low_highcut, fs,axis=0,order=4,no_trend=False):
    '''
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
    '''  
    
    #Create the filter
    low_highcut=np.array(low_highcut)
    nyq = 0.5 * fs
    b, a = butter(order, low_highcut/nyq, btype=btype)
    
    #Detrend the data, this is the equivalent of doing linear regressions across the time axis at each grid point
    VAR_detrend=detrend(VAR, axis=axis, type='linear')
    VAR_trend=VAR-VAR_detrend #By substracting the detrend array from the variable, we get the trend
    
    VAR_f= filtfilt(b, a, VAR_detrend,axis=axis)
    
    if no_trend:
        return VAR_f
    else:
        return VAR_trend +VAR_f 

    
def zonal_decomposition(VAR):
    '''
    Decomposition into spherical harmonics. [A. Kling, 2020]
    Args:
        VAR:  Detrend variable for decomposition, latitude is SECOND to LAST and longitude is LAST  e.g. (time,lat,lon) or (time,lev,lat,lon)
    Returns:
        COEFFS      : coefficient for harmonic decomposion, shape is flatten e.g. (time,2,lat/2, lat/2) (time x lev,2,lat/2, lat/2)
        power_per_l : power spectral density, shape is re-organized, e.g. (time, lat/2) or  (time,lev,lat/2)
    ***NOTE***    
    Output size is (...,lat/2, lat/2) as latitude is the smallest dimension and to match the Nyquist frequency
    '''  
    init_shtools() #Not optimal but prevent issues when library is not installed
    var_shape=np.array(VAR.shape)
    
    #Flatten array e.g. turn (10,36,lat,lon) to (360,lat,lon)
    nflatten=int(np.prod(var_shape[:-2]))
    reshape_flat=np.append(nflatten,var_shape[-2:])
    VAR=VAR.reshape(reshape_flat)
    
    coeff_out_shape=np.append(var_shape[0:-2],[2,var_shape[-2]//2,var_shape[-2]//2])
    psd_out_shape=np.append(var_shape[0:-2],var_shape[-2]//2)
    coeff_flat_shape=np.append(nflatten,coeff_out_shape[-3:])
    COEFFS=np.zeros(coeff_flat_shape)
    
    psd=np.zeros((nflatten,var_shape[-2]//2))
    
    for ii in range(0,nflatten):
        COEFFS[ii,...]=pyshtools.expand.SHExpandDH(VAR[ii,...], sampling=2)
        psd[ii,:]=pyshtools.spectralanalysis.spectrum(COEFFS[ii,...])
        
    return  COEFFS, psd.reshape(psd_out_shape)  

def zonal_construct(COEFFS_flat,VAR_shape,btype=None,low_highcut=None):
    '''
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
    '''  
    init_shtools() #Not optimal but prevent issues when library is not installed
    #--Initialization---
    nflatten=COEFFS_flat.shape[0]
    kmin=0
    kmax=COEFFS_flat.shape[-1]
    
    VAR=np.zeros((nflatten,VAR_shape[-2],VAR_shape[-1]))
    
    if btype=='low':  kmax= int(low_highcut)
    if btype=='high': kmin= int(low_highcut)
    if btype=='band': kmin,kmax= int(low_highcut[0]), int(low_highcut[1])

    for ii in range(0,nflatten):
        #=========Filtering===========
        COEFFS_flat[ii,:, :kmin, :] = 0.
        COEFFS_flat[ii,:, kmax:, :] = 0.
        VAR[ii,:,:] = pyshtools.expand.MakeGridDH(COEFFS_flat[ii,:,:], sampling=2)
    return  VAR.reshape(VAR_shape)