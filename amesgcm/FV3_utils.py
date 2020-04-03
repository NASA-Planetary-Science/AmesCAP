import numpy as np

def fms_press_calc(psfc,ak,bk,lev_type='full'):
    """
    Return the 3d pressure field from the surface pressure and the ak/bk coefficients.

    Args:
        psfc: the surface pressure in [Pa] or array of surface pressures 1D or 2D, or 3D (if time dimension)
        ak: 1st vertical coordinate parameter
        bk: 2nd vertical coordinate parameter
        lev_type: "full" (centers of the levels) or "half" (layer interfaces)
                  Default is "full"
    Returns:
        The 3D pressure field at the full PRESS_f(:,:,Nk-1) or half levels PRESS_h(:,:,Nk) in [Pa]
    --- 0 --- TOP        ========  p_half
    --- 1 ---                    
                         --------  p_full
    
                         ========  p_half
    ---Nk-1---           --------  p_full
    --- Nk --- SFC       ========  p_half 
                        / / / / /
    
    *NOTE* 
        Some litterature uses pk (pressure) instead of ak. 
        With p3d=  ps*bk +pref*ak  vs the current  p3d= ps*bk +ak      
        
         
    """
    
    Nk=len(ak)
    # If psfc is a float (e.g. psfc=7.) make it a one element array (e.g. psfc=[7.])
    if len(np.atleast_1d(psfc))==1: psfc=np.array([np.squeeze(psfc)])
        
    #Flatten the pressure array to generalize to N dimensions
    psfc_flat=psfc.flatten()
    
    # Expands the dimensions vectorized calculations: 
    psfc_v=np.repeat(psfc_flat[:,np.newaxis],Nk, axis=1)    #(Np) ->(Np,Nk)
    ak_v=np.repeat(ak[np.newaxis,:],len(psfc_flat), axis=0) #(Nk) ->(Np,Nk)
    bk_v=np.repeat(bk[np.newaxis,:],1, axis=0)              #(Nk) ->(1, Nk)
    
    #Pressure at half level = layers interfaces. The size of z axis is Nk
    PRESS_h=psfc_v*bk_v+ak_v
    
    #Pressure at full levels = centers of the levels. The size of z axis is Nk-1
    PRESS_f=np.zeros((len(psfc_flat),Nk-1))
    
    #Top layer (1st element is i=0 in Python)
    if ak[0]==0 and bk[0]==0: 
        PRESS_f[:,0]= 0.5*(PRESS_h[:,0]+PRESS_h[:,1])
    else:
        PRESS_f[:,0] = (PRESS_h[:,1]-PRESS_h[:,0])/np.log(PRESS_h[:,1]/PRESS_h[:,0])
    
    #Rest of the column (i=1..Nk).
    #[2:] goes from the 3rd element to Nk and [1:-1] goes from the 2nd element to Nk-1
    PRESS_f[:,1:]= (PRESS_h[:,2:]-PRESS_h[:,1:-1])/np.log(PRESS_h[:,2:]/PRESS_h[:,1:-1])
    
    # Reshape PRESS(:,Nk) to the original pressure shape PRESS(:,:,Nk) (resp. Nk-1)
                       
    if lev_type=="full":
        new_dim_f=np.append(psfc.shape,Nk-1)
        return np.squeeze(PRESS_f.reshape(new_dim_f))
    elif lev_type=="half" :  
        new_dim_h=np.append(psfc.shape,Nk)
        return np.squeeze(PRESS_h.reshape(new_dim_h))
    else: 
        raise Exception("""Pressure levels type not recognized in press_lev(): use 'full' or 'half' """)

def fms_Z_calc(psfc,ak,bk,T,topo=0.,lev_type='full'):
    """
    Return the 3d altitude field in [m]

    Args:
        psfc: the surface pressure in [Pa] or array of surface pressures 1D or 2D, or 3D (if time dimension)
        ak: 1st vertical coordinate parameter
        bk: 2nd vertical coordinate parameter
        T : the air temperature profile, 1D array (for a single grid point) or 2D, 3D 4D 
        topo: the surface elevation, same dimension as psfc
        lev_type: "full" (centers of the levels) or "half" (layer interfaces)
                  Default is "full"
    Returns:
        The layers' altitude  at the full Z_f(:,:,Nk-1) or half levels Z_h(:,:,Nk) in [m]
    
    --- 0 --- TOP        ========  z_half
    --- 1 ---                    
                         --------  z_full
    
                         ========  z_half
    ---Nk-1---           --------  z_full
    --- Nk --- SFC       ========  z_half 
                        / / / / /
    
    
    *NOTE* 
        Calculation is derived from ./atmos_cubed_sphere_mars/Mars_phys.F90:
        We have dp/dz = -rho g => dz= dp/(-rho g) and rho= p/(r T)  => dz=rT/g *(-dp/p) 
        Let's definethe log-pressure u as u = ln(p). We have du = du/dp *dp = (1/p)*dp =dp/p
        
        Finally , we have dz for the half layers:  dz=rT/g *-(du) => dz=rT/g *(+dp/p)   with N the layers defined from top to bottom.
    """
    g=3.72 #acc. m/s2
    r_co2= 191.00 # kg/mol
    Nk=len(ak)
    #===get the half and full pressure levels from fms_press_calc==
    
    PRESS_f=fms_press_calc(psfc,ak,bk,'full') 
    PRESS_h=fms_press_calc(psfc,ak,bk,'half')
    
    # If psfc is a float, turn it into a one-element array:
    if len(np.atleast_1d(psfc))==1: 
        psfc=np.array([np.squeeze(psfc)])
    if len(np.atleast_1d(topo))==1:    
        topo=np.array([np.squeeze(topo)])
        
    psfc_flat=psfc.flatten()
    topo_flat=topo.flatten()
    
    #  reshape arrays for vector calculations and compute the log pressure====
    
    PRESS_h=PRESS_h.reshape((len(psfc_flat),Nk))
    PRESS_f=PRESS_f.reshape((len(psfc_flat),Nk-1))
    T=T.reshape((len(psfc_flat),Nk-1))
    
    logPPRESS_h=np.log(PRESS_h)
    
    #===Initialize the output arrays===
    Z_f=np.zeros((len(psfc_flat),Nk-1))
    Z_h=np.zeros((len(psfc_flat),Nk))

    #First half layer is equal to the surface elevation
    
    Z_h[:,-1] = topo_flat
    
    # Other layes, from the bottom-ip:
    for k in range(Nk-2,-1,-1):
        Z_h[:,k] = Z_h[:,k+1]+(r_co2*T[:,k]/g)*(logPPRESS_h[:,k+1]-logPPRESS_h[:,k])
        Z_f[:,k] = Z_h[:,k+1]+(r_co2*T[:,k]/g)*(1-PRESS_h[:,k]/PRESS_f[:,k])
        
    #return the arrays
    if lev_type=="full":
        new_dim_f=np.append(psfc.shape,Nk-1)
        return Z_f.reshape(new_dim_f)
    elif lev_type=="half" : 
        new_dim_h=np.append(psfc.shape,Nk)
        return  Z_h.reshape(new_dim_h)
    #=====return the levels in Z coordinates [m]====
    else: 
        raise Exception("""Altitudes levels type not recognized: use 'full' or 'half' """)
        

def find_n(Lfull,Llev):
    '''
    Return the index for the level(s) just below Llev.
    Args:
        Lfull (array)         : input pressure [pa] or altitude [m] at full levels, level dimension is FIRST 
        Llev (float or 1D array) : desired level for interpolation [Pa] or [m]
    Returns:
        n:    index for the level(s) where the pressure is just below plev.
        alpha: alpha coefficients for the interpolation
    ***NOTE***
        - if Lfull is 1D array and  Llev is a float        > n is a float 
        - if Lfull is ND [lev,time,lat,lon] and Llev is a 1D array of size klev > n is an array of size[klev,Ndim]
          with Ndim =time x lat x lon
    '''
                   #number of original layers
                   
    Nlev=len(np.atleast_1d(Llev))
    if Nlev==1:Llev=np.array([Llev])
    dimsIN=Lfull.shape                         #get input variable dimensions
    Nfull=dimsIN[0]  
    dimsOUT=tuple(np.append(Nlev,dimsIN[1:]))
    Ndim= np.int(np.prod(dimsIN[1:]))           #Ndim is the product  of all dimensions but the vertical axis
    Lfull= np.reshape(Lfull, (Nfull, Ndim) )               
    ncol=Lfull.shape[-1]
    n=np.zeros((Nlev,ncol),dtype=int)
    
    for i in range(0,Nlev):
        for j in range(0,ncol) :
            n[i,j]=np.argmin(np.abs(Lfull[:,j]-Llev[i]))
            if Lfull[n[i,j],j]>Llev[i]:n[i,j]=n[i,j]-1
    return n
    
       
def pinterp(varIN,pfull,plev,masktop=True,index=None):
    '''
    Logarithmic interpolation pressure interpolation.   Alex Kling 3-26-20
    Args:
        varIN: variable to interpolate (N-dimensional array with vertical axis first)
        pfull: pressure at full layers same dimensions as varIN
        plev : desired level for interpolation as a 1D array 
        masktop: set to NaN values if above the model top
        index: indices for the interpolation, already procesed as [klev,Ndim] 
               Indices will be recalculated in not provided.
    Returns:
        varOUT: variable interpolated on the plev pressure levels
    
        ---  0  --- TOP    [e.g]   |    X_OUT= Xn*A + (1-A)*Xn+1
        ---  1  ---                |          
                                   |    with A = log(plev/pn)/log(pn+1/pn)
        ---  n  ---  pn   [30 Pa]  |Xn
                                   |
    >>> ---  k  ---  plev [100 Pa] |X_OUT       
        --- n+1 ---  pn+1 [200 Pa] |Xn+1
    
        --- SFC ---      
        / / / / / /
        
    '''
    #Special case where only 1 layer is requested
    Nlev=len(np.atleast_1d(plev))
    if Nlev==1:Llev=np.array([plev])
 
    dimsIN=varIN.shape               #get input variable dimensions
    Nfull=dimsIN[0]
    
    #Special case where varIN and pfull are a single profile            
    if len(varIN.shape )==1:varIN=varIN.reshape([Nfull,1])
    if len(pfull.shape )==1:pfull=pfull.reshape([Nfull,1])
       
    dimsIN=varIN.shape       #repeat in case varIN and pfull were reshaped
               
    dimsOUT=tuple(np.append(Nlev,dimsIN[1:]))
    Ndim= np.int(np.prod(dimsIN[1:]))          #Ndim is the product  of all dimensions but the vertical axis
    varIN= np.reshape(varIN, (Nfull, Ndim))    #flatten the other dimensions to (Nfull, Ndim)
    pfull= np.reshape(pfull, (Nfull, Ndim) )   #flatten the other dimensions to (Nfull, Ndim)
    varOUT=np.zeros((Nlev, Ndim))

    Ndimall=np.arange(0,Ndim)                   #all indices (does not change)
    
    for k in range(0,Nlev):
        #progress(k,Nlev)   #Display progress bar
        #Find nearest layer to plev[k]
        if np.any(index):
            #index have been pre-computed:  
            n= index[k,:]
        else:
            #Compute index on the fly for that layer
            n= np.squeeze(find_n(pfull,plev[k]))
                
        #==Slower method (but explains what is done below): loop over Ndim======
        # for ii in range(Ndim):
        #     if n[ii]<Nfull-1:
        #         alpha=np.log(plev[k]/pfull[n[ii]+1,ii])/np.log(pfull[n[ii],ii]/pfull[n[ii]+1,ii])
        #         varOUT[k,ii]=varIN[n[ii],ii]*alpha+(1-alpha)*varIN[n[ii]+1,ii]
        
        #=================    Fast method  no loop  =======================
        #Convert the layers n to indexes, for a 2D matrix using nindex=i*ncol+j
        nindex  =    n*Ndim+Ndimall  # n
        nindexp1=(n+1)*Ndim+Ndimall  # n+1
        
        #initialize alpha, size is [Ndim]
        alpha=np.NaN*Ndimall
        #Only calculate alpha  where the indices are <Nfull
        Ndo=Ndimall[nindexp1<Nfull*Ndim]
        alpha[Ndo]=np.log(plev[k]/pfull.flatten()[nindexp1[Ndo]])/np.log(pfull.flatten()[nindex[Ndo]]/pfull.flatten()[nindexp1[Ndo]])
        
        #Mask if plev[k]<model top
        if masktop: alpha[plev[k]<pfull.flatten()[nindex]]=np.NaN
        #Here, we need to make sure n+1 is never> Nfull by setting n+1=Nfull, if it it the case.
        #This does not affect the calculation as alpha is set to NaN for those values. 
        nindexp1[nindexp1>=Nfull*Ndim]=nindex[nindexp1>=Nfull*Ndim]

        varOUT[k,:]=varIN.flatten()[nindex]*alpha+(1-alpha)*varIN.flatten()[nindexp1]
        
    return np.reshape(varOUT,dimsOUT)


        
def akbk_loader(NLAY,data_dir='/u/mkahre/MCMC/data_files'):
    """
    Return the ak and bk values given a number of layers for standards resolutions 
    Default directory is /lou/s2n/mkahre/MCMC/data_files/ 
    Args:
        NLAY: the number of layers (float or integer)
    Returns:
        ak: 1st vertical coordinate parameter [Pa]
        bk: 2nd vertical coordinate parameter [none]
    
    *NOTE*    ak,bk have a size NLAY+1 since they define the position of the layer interfaces (half layers):
              p_half = ak + bk*p_sfc 
    """  
        
    from netCDF4 import Dataset
    NLAY=int(NLAY)
    file=Dataset(data_dir+'/akbk_L%i.nc'%(NLAY), 'r', format='NETCDF4_CLASSIC')
    ak=file.variables['pk'][:]
    bk=file.variables['bk'][:]
    file.close()
    return ak,bk 
    
    
def zonal_avg_P_lat(Ls,var,Ls_target,Ls_angle,symmetric=True):
    """
    Return the zonally averaged mean value of a pressure interpolated 4D variable.

    Args:
        Ls: 1D array of solar longitude of the input variable in degree (0->360)
        var: a 4D variable var [time,levels,lat,lon] interpolated on the pressure levels (f_average_plevs file)
        Ls_target: central solar longitude of interest.     
        Ls_angle:  requested window angle centered around   Expl:  Ls_angle = 10.  (Window will go from Ls 85  
        symmetric: a boolean (default =True) If True, and if the requested window is out of range, Ls_angle is reduced
                                             If False, the time average is done on the data available
    Returns:
        The zonnally and latitudinally-averaged field zpvar[level,lat]
    
    Expl:  Ls_target= 90.
           Ls_angle = 10.  
           
           ---> Nominally, the time average is done over solar longitudes      85 <Ls_target < 95 (10 degree)
           
           ---> If  symmetric =True and the input data ranges from Ls 88 to 100     88 <Ls_target < 92 (4  degree, symmetric)
                If  symmetric =False and the input data ranges from Ls 88 to 100    88 <Ls_target < 95 (7  degree, assymetric)
    *NOTE* 
    
    [Alex] as of 6/8/18, the routine will bin data from muliples Mars years if provided
         
    """
    #compute bounds from Ls_target and Ls_angle
    Ls_min= Ls_target-Ls_angle/2.
    Ls_max= Ls_target+Ls_angle/2.
    
    if (Ls_min<0.):Ls_min+=360.
    if (Ls_max>360.):Ls_max-=360. 
    
    #Initialize output array
    zpvar=np.zeros((var.shape[1],var.shape[2])) #nlev, nlat
    
    #check is the Ls of interest is within the data provided, raise execption otherwise
    if Ls_target <= Ls.min() or Ls_target >=Ls.max() :
        raise Exception("Error \nNo data found, requested  data :       Ls %.2f <-- (%.2f)--> %.2f\nHowever, data in file only ranges      Ls %.2f <-- (%.2f)--> %.2f"%(Ls_min,Ls_target,Ls_max,Ls.min(),(Ls.min()+Ls.max())/2.,Ls.max()))

    
    else : #If only some of the requested data is outside the ranges, process this data
        if Ls_min <Ls.min() or Ls_max >Ls.max():
            print("In zonal_avg_P_lat() Warning: \nRequested  data ranging    Ls %.2f <-- (%.2f)--> %.2f"%(Ls_min,Ls_target,Ls_max))
            if symmetric: #Case 1: reduce the window
                if Ls_min <Ls.min():
                    Ls_min =Ls.min()
                    Ls_angle=2*(Ls_target-Ls_min)
                    Ls_max= Ls_target+Ls_angle/2.
                    
                if Ls_max >Ls.max():
                    Ls_max =Ls.max()
                    Ls_angle=2*(Ls_max-Ls_target)
                    Ls_min= Ls_target-Ls_angle/2.
                    
                print("Reshaping data ranging     Ls %.2f <-- (%.2f)--> %.2f"%(Ls_min,Ls_target,Ls_max))        
            else: #Case 2: Use all data available
                print("I am only using            Ls %.2f <-- (%.2f)--> %.2f \n"%(max(Ls.min(),Ls_min),Ls_target,min(Ls.max(),Ls_max)))
    count=0
    #perform longitude average on the field
    zvar= np.mean(var,axis=3)
    
    for t in xrange(len(Ls)):
    #special case Ls around Ls =0 (wrap around)
        if (Ls_min<=Ls[t] <= Ls_max):
            zpvar[:,:]=zpvar[:,:]+zvar[t,:,:]
            count+=1
            
    if  count>0:
        zpvar/=count
    return zpvar
    

    
def alt_KM(press,scale_height_KM=8.,reference_press=610.):
    """
    Gives the approximate altitude in km for a given pressure
    Args:
        press: the pressure in [Pa]
        scale_height_KM: a scale height in [km], (default is 10 km)
        reference_press: reference surface pressure in [Pa], (default is 610 Pa)
    Returns:
        z_KM: the equivalent altitude for that pressure level in [km]
   
    """      
    return -scale_height_KM*np.log(press/reference_press) # p to altitude in km      
    
def press_pa(alt_KM,scale_height_KM=8.,reference_press=610.):
    """
    Gives the approximate altitude in km for a given pressure
    Args:
        alt_KM: the altitude in  [km]
        scale_height_KM: a scale height in [km], (default is 8 km)
        reference_press: reference surface pressure in [Pa], (default is 610 Pa)
    Returns:
         press_pa: the equivalent pressure at that altitude in [Pa]
   
    """      
    return reference_press*np.exp(-alt_KM/scale_height_KM) # p to altitude in km 
     
def lon180_to_360(lon):
    lon=np.array(lon)
    """
    Transform a float or an array from the -180/+180 coordinate system to 0-360
    Args:
        lon: a float, 1D or 2D array of longitudes in the 180/+180 coordinate system
    Returns:
        lon: the equivalent longitudes in the 0-360 coordinate system
   
    """ 
    if len(np.atleast_1d(lon))==1: #lon180 is a float
        if lon<0:lon+=360 
    else:                            #lon180 is an array
        lon[lon<0]+=360
        lon=np.append(lon[lon<180],lon[lon>=180]) #reogranize lon by increasing values
    return lon
    
def lon360_to_180(lon):
    lon=np.array(lon)
    """
    Transform a float or an array from the 0-360 coordinate system to -180/+180
    Args:
        lon: a float, 1D or 2D array of longitudes in the 0-360 coordinate system
    Returns:
        lon: the equivalent longitudes in the -180/+180 coordinate system
   
    """ 
    if len(np.atleast_1d(lon))==1:   #lon is a float
        if lon>180:lon-=360 
    else:                            #lon is an array
        lon[lon>180]-=360
        lon=np.append(lon[lon<0],lon[lon>=0]) #reogranize lon by increasing values
    return lon    
     
def shiftgrid_360_to_180(lon,data): #longitude is LAST
    '''
    This function shift N dimensional data a 0->360 to a -180/+180 grid.
    Args:
        lon: 1D array of longitude 0->360
        data: ND array with last dimension being the longitude (transpose first if necessary)
    Returns:
        data: shifted data
    Note: Use np.ma.hstack instead of np.hstack to keep the masked array properties
    '''
    lon=np.array(lon)
    lon[lon>180]-=360. #convert to +/- 180
    data=np.concatenate((data[...,lon<0],data[...,lon>=0]),axis=-1) #stack data
    return data


def shiftgrid_180_to_360(lon,data): #longitude is LAST
    '''
    This function shift N dimensional data a -180/+180 grid to a 0->360
    Args:
        lon: 1D array of longitude -180/+180
        data: ND array with last dimension being the longitude (transpose first if necessary)
    Returns:
        data: shifted data
    Note: Use np.ma.hstack instead of np.hstack to keep the masked array properties
    '''
    lon=np.array(lon)
    lon[lon<0]+=360. #convert to 0-360
    data=np.concatenate((data[...,lon<180],data[...,lon>=180]),axis=-1) #stack data
    return data
        
def second_hhmmss(seconds,lon_180=0.,show_mmss=True):
    """
    Given the time seconds return Local true Solar Time at a certain longitude
    Args:
        seconds: a float, the time in seconds
        lon_180: a float, the longitude in a -/+180 coordinate
        show_mmss: returns min and second if true
    Returns:
        hours: float, the local time or  (hours,minutes, seconds)
   
    """ 
    hours = seconds // (60*60)
    seconds %= (60*60)
    minutes = seconds // 60
    seconds %= 60
    #Add timezone offset (1hr/15 degree)
    hours=np.mod(hours+lon_180/15.,24)
    
    if show_mmss:
        return np.int(hours), np.int(minutes), np.int(seconds)
    else:
        return np.int(hours)   


def sol2LTST(time_sol,lon_180=0.,show_minute=False):
    """
    Given the time in days, return the Local true Solar Time at a certain longitude
    Args:
        time_sol: a float, the time, eg. sols 2350.24
        lon_180: a float, the longitude in a -/+180 coordinate
        show_minute: show minutes if true, otherwise show whole hours
    Returns:
        hours: float, the local time or  (hours,minutes, seconds)
   
    """ 
    return second_hhmmss(time_sol*86400.,lon_180,show_minute)

def space_time(lon,timex, varIN,kmx,tmx):
    """
    Obtain west and east propagating waves. This is a Python implementation of John Wilson's  space_time routine by Alex
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
    h=np.array(h)
    if h.any():
        # h is provided as a 1D array
        if len(h.shape)==1:
            d_arr = np.copy(arr)
            reshape_shape=np.append([arr.shape[0]-2],[1 for i in range(0,arr.ndim -1)]) 
            d_arr[0,...] = (arr[1,...]-arr[0,...])/(h[1]-h[0])
            d_arr[-1,...] = (arr[-1,...]-arr[-2,...])/(h[-1]-h[-2])
            d_arr[1:-1,...] = (arr[2:,...]-arr[0:-2,...])/(np.reshape(h[2:]-h[0:-2],reshape_shape))
         
        #h has the same dimension as var    
        elif h.shape==arr.shape:
            d_arr = np.copy(arr)
            d_arr[0,...] = (arr[1,...]-arr[0,...])/(h[1,...]-h[0,...])
            d_arr[-1,...] = (arr[-1,...]-arr[-2,...])/(h[-1,...]-h[-2,...])
            d_arr[1:-1,...] = (arr[2:,...]-arr[0:-2,...])/(h[2:,...]-h[0:-2,...])
        else:     
            print('Error,h.shape=', h.shape,'arr.shape=',arr.shape)
        
    # h is not defined, we return only d_var, not d_var/dh
    else:
        d_arr = np.copy(arr)
        reshape_shape=np.append([arr.shape[0]-2],[1 for i in range(0,arr.ndim -1)]) 
        d_arr[0,...] = arr[1,...]-arr[0,...]
        d_arr[-1,...] = arr[-1,...]-arr[-2,...]
        d_arr[1:-1,...] = arr[2:,...]-arr[0:-2,...]
        
    
    return d_arr

#========================================================================= 
#=============Wrapper for creation of netcdf files========================
#=========================================================================

class Ncdf(object):
    '''
    Alex K.
    NetCdf wrapper for quick archiving of data into netcdf format 
    
    USAGE: 
    
    from netcdf_wrapper import Ncdf

    Fgeo= 0.03 #W/m2, a constant
    TG=np.ones((24,8)) #ground temperature

    #---create file---
    filename="/lou/s2n/mkahre/MCMC/analysis/working/myfile.nc"
    description="results from new simulation, Alex 01-01-19"
    Log=Ncdf(filename,description)
    
    #---Save the constant to the file---
    Log.add_constant('Fgeo',Fgeo,"geothermal flux","W/m2")
    
    #---Save the TG array to the file---
    Log.add_dimension('Nx',8)
    Log.add_dimension('time',24)
    
    Log.log_variable('TG',TG,('time','Nx'),'soil temperature','K')
    
    Log.close()
    
    
    '''
    def __init__(self,filename=None,description_txt="",action='w'):
        if filename: 
            if filename[-3:]!=".nc":
            #assume that only path is provided so make a name for the file
                import datetime;now = datetime.datetime.now()
                filename=filename+\
                '/run_%02d-%02d-%04d_%i-%i-%i.nc'%(now.day,now.month,now.year,now.hour,now.minute,now.second)
        else:   #create a default file name  if path and filename are not provided
            import os #use a default path if not provided
            pathname=os.getcwd()+'/'
            import datetime;now = datetime.datetime.now()
            filename=pathname+\
            'run_%02d-%02d-%04d_%i-%i-%i.nc'%(now.day,now.month,now.year,now.hour,now.minute,now.second)
        self.filename=filename
        from netCDF4 import Dataset 
        if action=='w':
            self.f_Ncdf = Dataset(filename, 'w', format='NETCDF4_CLASSIC')
            self.f_Ncdf.description = description_txt
        elif action=='a': #append to file
            self.f_Ncdf = Dataset(filename, 'a', format='NETCDF4_CLASSIC')
        #create dictionaries to hold dimensions and variables
        self.dim_dict=dict()
        self.var_dict=dict()
        #print(filename+ " was created")
        
    def close(self):
        self.f_Ncdf.close()
        print(self.filename+" was closed")
        
    def add_dimension(self,dimension_name,length):
        self.dim_dict[dimension_name]= self.f_Ncdf.createDimension(dimension_name,length)
        
    def print_dimension(self):
        print(self.dim_dict.items())
    def print_variable(self):
        print(self.var_dict.keys())    
        
    def add_constant(self,variable_name,value,longname_txt="",unit_txt=""):
        if not any('constant' in s for s in self.dim_dict.keys()):
            self.add_dimension('constant',1)
        longname_txt =longname_txt+' (%g)'%(value)   #add the value to the longname
        self._def_variable(variable_name,('constant'),longname_txt,unit_txt)
        self.var_dict[variable_name][:]=value
    #=====Private definitions=====   
    def _def_variable(self,variable_name,dim_array,longname_txt="",unit_txt=""):
        self.var_dict[variable_name]= self.f_Ncdf.createVariable(variable_name,'f4',dim_array)    
        self.var_dict[variable_name].units=unit_txt
        self.var_dict[variable_name].long_name=longname_txt
        self.var_dict[variable_name].dim_name=str(dim_array)  

    def _def_axis1D(self,variable_name,dim_array,longname_txt="",unit_txt="",cart_txt=""):
        self.var_dict[variable_name]= self.f_Ncdf.createVariable(variable_name,'f8',dim_array)
        self.var_dict[variable_name].units=unit_txt
        self.var_dict[variable_name].long_name=longname_txt
        self.var_dict[variable_name].cartesian_axis=cart_txt
    #================================    
    def log_variable(self,variable_name,DATAin,dim_array,longname_txt="",unit_txt=""):
        if not any(variable_name in s for s in self.var_dict.keys()):
            self._def_variable(variable_name,dim_array,longname_txt,unit_txt)
        self.var_dict[variable_name].long_name=longname_txt
        self.var_dict[variable_name].dim_name=str(dim_array)  
        self.var_dict[variable_name].units=unit_txt
        self.var_dict[variable_name][:]=DATAin 
        
    def log_axis1D(self,variable_name,DATAin,dim_array,longname_txt="",unit_txt="",cart_txt="",):
        if not any(variable_name == s for s in self.var_dict.keys()):
            self._def_axis1D(variable_name,dim_array,longname_txt,unit_txt,cart_txt)
        self.var_dict[variable_name].long_name=longname_txt
        self.var_dict[variable_name].units=unit_txt
        self.var_dict[variable_name].cartesian_axis=cart_txt
        self.var_dict[variable_name][:]=DATAin
        
    #Function to define a dimension and add a variable with at the same time
    #lon_array=np.linspace(0,360)
    #Log.add_dim_with_content('lon',lon_array,'longitudes','degree')
    def add_dim_with_content(self,dimension_name,DATAin,longname_txt="",unit_txt=""):
        self.add_dimension(dimension_name,len(DATAin))
        #---If no longname is provided, use dimension_name as default---
        if longname_txt=="":longname_txt=dimension_name
        self.log_axis1D(dimension_name,DATAin,(dimension_name),longname_txt,unit_txt) 
         
    #Copy a netcdf DIMENSION variable e.g Ncdim is:  f.variables['lon']
    # if the dimension for that variable does not exist, it will be created
    def copy_Ncdim_with_content(self,Ncdim_var):
        longname_txt=getattr(Ncdim_var,'long_name',Ncdim_var.name)
        unit_txt=    getattr(Ncdim_var,'units','')
        self.add_dim_with_content(Ncdim_var.name,Ncdim_var[:],longname_txt,unit_txt)
               
    #Copy a netcdf variable from another file, e.g Ncvar is: f.variables['ucomp']
    def copy_Ncvar(self,Ncvar):
        dim_array=Ncvar.dimensions
        longname_txt=getattr(Ncvar,'long_name',Ncvar.name)
        unit_txt=    getattr(Ncvar,'units','')
        self._def_variable(Ncvar.name,Ncvar.dimensions,longname_txt,unit_txt)
        self.log_variable(Ncvar.name,Ncvar[:],Ncvar.dimensions,longname_txt,unit_txt)
    
#=====TEST ONLY=======
'''       
fname='/Users/akling/test/00030.atmos_average.nc'

f=Dataset(fname,'r')
varnames=f.variables.keys()
lat=f.variables['lat']
lon=f.variables['lon'][:]  
ucomp=f.variables['ucomp']  


Log=Ncdf('/Users/akling/mytest2.nc')
Fgeo= 0.03 #W/m2, a constant

#Log.add_dimension('Nx',8)
Log.add_dim_with_content('lon',lon,'longitudes','degree')
Log.copy_Ncdim_with_content(f.variables['lat'])
Log.copy_Ncdim_with_content(f.variables['pfull'])
Log.copy_Ncdim_with_content(f.variables['time'])
Log.copy_Ncvar(f.variables['ucomp']  )
Log.close()  
'''
  
  
    
#========================================================================= 
#=======================vertical grid utilities===========================
#=========================================================================

def gauss_profile(x, alpha,x0=0.):
    """ Return Gaussian line shape at x This can be used to generate a bell-shaped mountain"""
    return np.sqrt(np.log(2) / np.pi) / alpha\
                             * np.exp(-((x-x0) / alpha)**2 * np.log(2))

def compute_uneven_sigma(num_levels, N_scale_heights, surf_res, exponent, zero_top ):
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
    b=np.zeros(int(num_levels)+1)
    for k in range(0,num_levels):
        zeta = 1.-k/np.float(num_levels) #zeta decreases with k
        z  = surf_res*zeta + (1.0 - surf_res)*(zeta**exponent)
        b[k] = np.exp(-z*N_scale_heights)
    b[-1] = 1.0
    if(zero_top):  b[0] = 0.0
    return b


def transition( pfull, p_sigma=0.1, p_press=0.05):
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
    t=np.zeros_like(pfull)
    for k in range(0,len(pfull)):
        if( pfull[k] <= p_press): 
            t[k] = 0.0
        elif ( pfull[k] >= p_sigma) :
            t[k] = 1.0
        else:
            x  = pfull[k]    - p_press
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

    ktrans= np.argmin(np.abs( plev- ptrans) ) # ks= number of pure pressure levels
    km= len(plev)-1
    
    aknew=np.zeros(len(plev))
    bknew=np.zeros(len(plev))
    
    #   pnorm= 1.e5; 
    pnorm= psfc 
    eta= plev / pnorm
    
    ep= eta[ktrans+1]       #  ks= number of pure pressure levels
    es= eta[-1]
    rnorm= 1. / (es-ep)**2
    
    #   Compute alpha, beta, and gamma using Swinbank's formula
    alpha = (ep**2 - 2.*ep*es) / (es-ep)**2
    beta  =        2.*ep*es**2 / (es-ep)**2
    gamma =        -(ep*es)**2 / (es-ep)**2
    
    #   Pure Pressure levels 
    aknew= eta * pnorm
    
    #  Hybrid pressure-sigma levels
    kdex= range(ktrans+1,km) 
    aknew[kdex] = alpha*eta[kdex] + beta + gamma/eta[kdex] 
    aknew[kdex]= aknew[kdex] * pnorm
    aknew[-1]= 0.0
    
    bknew[kdex] = (plev[kdex] - aknew[kdex])/psfc 
    bknew[-1] = 1.0 
    
    #find the transition level ks where (bk[ks]>0)
    ks=0
    while bknew[ks]==0. :
        ks+=1
    #ks is the one that would be use in fortran indexing in fv_eta.f90    
    return  aknew, bknew,ks
   

def polar_warming(T,lat,outside_range=np.NaN):
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
    
    def PW_half_hemisphere(T_half,lat_half,outside_range=np.NaN):

        #Easy case, T is a 1D on the latitude direction only
        if len(T_half.shape)==1:
            imin=np.argmin(T_half)
            imax=np.argmax(T_half)
            
            #Note that we compute polar warming at ALL latitudes and then set NaN the latitudes outside the desired range. 
            #We test on the absolute values (np.abs) of the latitudes, therefore the function is usable on both hemispheres
            DT_PW_half=T_half-T_half[imin] 
            exclude=np.append(np.where(np.abs(lat_half)-np.abs(lat_half[imin])<0),np.where(np.abs(lat_half[imax])-np.abs(lat_half)<0))
            DT_PW_half[exclude]=outside_range #set to NaN
            return DT_PW_half
        
        #General case for N dimensions        
        else:
            
            #Flatten the diemnsions but the first dimension (the latitudes)
            arr_flat=T_half.reshape([T_half.shape[0], np.prod(T_half.shape[1:])])
            LAT_HALF=np.repeat(lat_half[:,np.newaxis],arr_flat.shape[1],axis=1)
            
            imin=np.argmin(arr_flat,axis=0)
            imax=np.argmax(arr_flat,axis=0)
            
            #Initialize four working arrays
            tmin0,tmax0,latmin0,latmax0=[np.zeros_like(arr_flat) for _ in range(4)]
            
            #get the min/max temperature and latitudes
            for i in range(0,arr_flat.shape[1]):
                tmax0[:,i]=arr_flat[imax[i],i]
                tmin0[:,i]=arr_flat[imin[i],i]
                latmin0[:,i]=lat_half[imin[i]]
                latmax0[:,i]=lat_half[imax[i]]
            
            #Compute polar warming for that hemisphere
            DT_PW_half=arr_flat-tmin0
            
            # Set to NaN values outside the range. 
            tuple_lower_than_latmin=np.where(np.abs(LAT_HALF)-np.abs(latmin0)<0)
            tuple_larger_than_latmax=np.where(np.abs(latmax0)-np.abs(LAT_HALF)<0)
            
            DT_PW_half[tuple_lower_than_latmin]=outside_range
            DT_PW_half[tuple_larger_than_latmax]=outside_range
        
            return DT_PW_half.reshape(T_half.shape)

    #======================================================
    #======Actual calculations for both hemispheres========
    #======================================================
    T_SH=T[0:len(lat)//2]
    lat_SH=lat[0:len(lat)//2]
    
    T_NH=T[len(lat)//2:]
    lat_NH=lat[len(lat)//2:]

    return np.concatenate((PW_half_hemisphere(T_SH,lat_SH,outside_range),PW_half_hemisphere(T_NH,lat_NH,outside_range)),axis=0) 
