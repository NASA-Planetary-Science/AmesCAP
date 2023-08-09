#!/usr/bin/env python3
"""
The MarsFormat executable is for ...

The executable requires x arguments:
    * [-x --x]      define

Third-party Requirements:
    * numpy
    * argparse
    * requests

List of Functions:
    * x
"""

# make print statements appear in color
from amescap.Script_utils import prCyan

# load generic Python modules
import argparse   # parse arguments
import numpy as np
import xarray as xr
import os

# load amesCAP modules
from amescap.FV3_utils import layers_mid_point_to_boundary

xr.set_options(keep_attrs=True)

#---
# fit2FV3.py
# Routine to Transform Model Input (variable names, dimension names, array order)
# to expected configuration CAP

# ======================================================
#                  ARGUMENT PARSER
# ======================================================

parser = argparse.ArgumentParser(
   description="""\033[93m fit2FV3.py  Used to convert model output to FV3 format  \n \033[00m""",
   formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('input_file', nargs='+',  # sys.stdin
                    help='***.nc file or list of ***.nc files ')


parser.add_argument('-openmars', '--openmars', nargs='+',
                    help="""Produce a FV3-like daily file \n"""
                    """> Usage: MarsFormat.py fileIN*.nc --model \n"""
                    """> Available options are:                     \n"""
                    """(-openmars)  --openmars  daily \n"""
                    """(-marswrf)   --marswrf   daily \n"""
                    """\n""")

parser.add_argument('-marswrf', '--marswrf', nargs='+',
                    help=argparse.SUPPRESS)

parser.add_argument('-legacy', '--legacy', nargs='+',
                    help=argparse.SUPPRESS)

# ======================================================
#                  MAIN PROGRAM
# ======================================================

def main():
   path2data = os.getcwd()
   # Open a single File
   file_list=parser.parse_args().input_file
   #path_inpt
   for filei in file_list:
      #Add path unless full path is provided
      if not ('/' in filei):
         fullnameIN = path2data + '/' + filei
      else:
         fullnameIN=filei
      fullnameOUT = fullnameIN[:-3]+'_atmos_daily.nc'

      print('Processing...')
      #dataDIR = path+filename+'.nc'
      DS = xr.open_dataset(fullnameIN, decode_times=False)

      #=================================================================
      # ===================OpenMars Specific Processing==================
      #=================================================================
      if parser.parse_args().marswrf:
         #TODO longname is 'description' for MarsWRF
         '''
         print('Input File content (description) and (description) attibutes:')
         print('------')
         for ivar in  DS.keys():
            print(ivar,DS[ivar].attrs['description'],DS[ivar].attrs['units'])
         print('------')
         '''
         #==================================================================
         # Find Shape of Coordinates
         #==================================================================
         # [t,z,y,x] = 100,43,90,180
         ppt_dims = np.shape(DS.T)
         lmax = ppt_dims[3] # x = 180
         jmax = ppt_dims[2] # y = 90
         tmax = ppt_dims[0] # t = 100
         pmax = ppt_dims[1] # z = 43 (layer)

         #==================================================================
         # Define Coordinates for New DataFrame
         #==================================================================
         time        = DS.XTIME/ 60/ 24         # minutes since simulation start [m]
         lat = DS.XLAT[0,:,0]
         lon2D = DS.XLONG[0,:]
         lon = np.squeeze(lon2D[0,:])

         # Derive half and full reference pressure levels (Pa)
         pfull = DS.P_TOP[0]+ DS.ZNU[0,:]* DS.P0
         phalf = DS.P_TOP[0]+ DS.ZNW[0,:]* DS.P0

         #==================================================================
         # Calculate *Level* Heights above the Surface (i.e. above topo)
         #==================================================================
         zagl_lvl = (DS.PH[:,:pmax,:,:] + DS.PHB[0,:pmax,:,:]) / DS.G - DS.HGT[0,:,:]

         #==================================================================
         # Find Layer Pressures [Pa]
         #==================================================================
         try:
            pfull3D = DS.P_TOP + DS.PB[0,:] # perturb. pressure + base state pressure, time-invariant
         except NameError:
            pfull3D = DS.PSFC[:,:jmax,:lmax] * DS.ZNU[:,:pmax]

         #==================================================================
         # Interpolate U, V, W, Zfull onto Regular Mass Grid (from staggered)
         #==================================================================
         # For variables staggered x (lon) [t,z,y,x'] -> regular mass grid [t,z,y,x]:
         ucomp = 0.5 * (DS.U[..., :-1] + DS.U[..., 1:])

         # For variables staggered y (lat) [t,z,y',x] -> regular mass grid [t,z,y,x]:
         vcomp = 0.5 * (DS.V[:,:,:-1,:] + DS.V[:,:,1:,:])

         # For variables staggered p/z (height) [t,z',y,x] -> regular mass grid [t,z,y,x]:
         w = 0.5 * (DS.W[:,:-1,:,:] + DS.W[:,1:,:,:])

         # ALSO INTERPOLATE TO FIND *LAYER* HEIGHTS ABOVE THE SURFACE (i.e., above topography; m)
         zfull3D = 0.5 * (zagl_lvl[:,:-1,:,:] + zagl_lvl[:,1:,:,:])

         #==================================================================
         # Derive attmospherice temperature [K]
         #==================================================================
         gamma = DS.CP / (DS.CP - DS.R_D)
         temp = (DS.T + DS.T0) * (pfull3D / DS.P0)**((gamma-1.) / gamma)

         #==================================================================
         # Derive ak, bk
         #==================================================================
         ak  = np.zeros(len(phalf))
         bk = np.zeros(len(phalf))
         ak[-1]=DS.P_TOP[0]  #MarsWRF comes with pressure increasing with N
         bk[:]=DS.ZNW[0,:]

         #==================================================================
         # Create New DataFrame
         #==================================================================
         coords = {'time': np.array(time), 'phalf':np.array(phalf),'pfull': np.array(pfull), 'lat': np.array(lat), 'lon': np.array(lon)}  # Coordinates dictionary

         archive_vars = {
            'ak' :         [ak, ['phalf'],'pressure part of the hybrid coordinate','Pa'],
            'bk' :         [bk, ['phalf'],'vertical coordinate sigma value','none'],
            'areo' :       [DS.L_S, ['time'],'solar longitude','degree'],
            'ps' :         [DS.PSFC, ['time','lat','lon'],'surface pressure','Pa'],
            'zsurf':       [DS.HGT[0,:],['lat','lon'],'surface height','m'],
            'ucomp':       [ucomp, ['time', 'pfull', 'lat', 'lon'],'zonal winds','m/sec'],
            'vcomp':       [vcomp, ['time', 'pfull', 'lat', 'lon'],'meridional winds','m/sec'],
            'w':           [w, ['time', 'pfull', 'lat', 'lon'],'vertical winds','m/s'],
            'pfull3D':     [pfull3D, ['time', 'pfull', 'lat', 'lon'],'pressure','Pa'],
            'temp':        [temp, ['time', 'pfull', 'lat', 'lon'],'temperature','K'],
            'h2o_ice_sfc': [DS.H2OICE, ['time','lat','lon'],'Surface H2O Ice','kg/m2'],
            'co2_ice_sfc': [DS.CO2ICE, ['time','lat','lon'],'Surface CO2 Ice','kg/m2'],
            'ts':          [DS.TSK, ['time','lat','lon'],'surface temperature','K'],
         }


      #=================================================================
      # ===================OpenMars Specific Processing==================
      #=================================================================
      elif parser.parse_args().openmars:
         '''
         print('Input File content (FIELDNAM) and (UNITS) attibutes:')
         print('------')
         for ivar in  DS.keys():
            print(ivar,DS[ivar].attrs['FIELDNAM'],DS[ivar].attrs['UNITS'])
         print('------')
         '''
         #==================================================================
         # Define Coordinates for New DataFrame
         #==================================================================
         ref_press=720 #TODO this is added on to create ak/bk
         time        = DS.time         # minutes since simulation start [m]
         lat = DS.lat
         lon = DS.lon
         # Derive half and full reference pressure levels (Pa)
         pfull = DS.lev*ref_press


         #==================================================================
         # add p_half dimensions and ak, bk vertical grid coordinates
         #==================================================================

         #DS.expand_dims({'p_half':len(pfull)+1})
         #Compute sigma values. Swap the sigma array upside down twice  with [::-1] since the layers_mid_point_to_boundary() needs sigma[0]=0, sigma[-1]=1) and then to reorganize the array in the original openMars format with sigma[0]=1, sigma[-1]=0
         DS['bk'] = layers_mid_point_to_boundary(DS.lev[::-1],1.)[::-1]
         DS['ak'] = np.zeros(len(pfull)+1) #Pure sigma model, set bk to zero
         phalf= np.array(DS['ak']) + ref_press*np.array(DS['bk'])  #compute phalf


         #==================================================================
         # Make New DataFrame
         #==================================================================

         coords = {'time': np.array(time), 'phalf': np.array(phalf), 'pfull': np.array(pfull), 'lat':np.array(lat), 'lon': np.array(lon)}

         #Variable to archive [name, values, dimensions, longname,units]
         archive_vars = {
         'bk' :          [DS.bk, ['phalf'],'vertical coordinate sigma value','none'],
         'ak' :          [DS.ak, ['phalf'], 'pressure part of the hybrid coordinate','Pa'],
         'areo' :        [DS.Ls, ['time'],'solar longitude','degree'],
         'ps' :          [DS.ps, ['time','lat','lon'],'surface pressure','Pa'],
         'ucomp':        [DS.u, ['time', 'pfull', 'lat', 'lon'],'zonal winds','m/sec'],
         'vcomp':        [DS.v, ['time', 'pfull', 'lat', 'lon'],'meridional wind','m/sec'],
         'temp':         [DS.temp, ['time', 'pfull', 'lat', 'lon'],'temperature','K'],
         'dust_mass_col':[DS.dustcol, ['time','lat','lon'],'column integration of dust','kg/m2'],
         'co2_ice_sfc':  [DS.co2ice, ['time','lat','lon'],'surace CO2 ice','kg/m2'],
         'ts':           [DS.tsurf, ['time','lat','lon'],'Surface Temperature','K']
         }



      #==================================================================
      #========Create output file (common to all models)=================
      #==================================================================
      archive_coords = {
      'time':  ['time','days'],
      'pfull': ['ref full pressure level','Pa'],
      'lat':   ['latitudes' ,'degrees_N'],
      'lon':   ['longitudes','degrees_E'],
      'phalf': ['ref pressure at layer boundaries','Pa']
      }

      # Empty xarray dictionary
      data_vars = {}
      # Assign description and units attributes to the xarray dictionary
      for ivar in archive_vars.keys():
         data_vars[ivar] = xr.DataArray(np.array(archive_vars[ivar][0]), dims=archive_vars[ivar][1])
         data_vars[ivar].attrs['long_name'] = archive_vars[ivar][2]
         data_vars[ivar].attrs['units'] =  archive_vars[ivar][3]

      # Create the dataset with the data variables and assigned attributes
      DF = xr.Dataset(data_vars, coords=coords)

      #Add longname and units attibutes to the coordiate variables
      for ivar in archive_coords.keys():
         DF[ivar].attrs['long_name']=archive_coords[ivar][0]
         DF[ivar].attrs['units']=archive_coords[ivar][1]



      #==================================================================
      # check that vertical grid starts at toa with highest level at surface
      #==================================================================
      if DF.pfull[0] != DF.pfull.min(): # if toa, lev = 0 is surface then flip
          DF=DF.isel(pfull=slice(None, None, -1)) # regrids DS based on pfull
          DF=DF.isel(phalf=slice(None, None, -1)) #Also flip phalf,ak, bk


      #==================================================================
      # reorder dimensions
      #==================================================================
      DF = DF.transpose("phalf","time", "pfull" ,"lat","lon")


      #==================================================================
      # change longitude from -180-179 to 0-360
      #==================================================================
      if min(DF.lon)<0:
            tmp = np.array(DF.lon)
            tmp = np.where(tmp<0,tmp+360,tmp)
            DF=DF.assign_coords({'lon':('lon',tmp,DF.lon.attrs)})
            DF = DF.sortby("lon")

      #==================================================================
      # add scalar axis to areo [time, scalar_axis])
      #==================================================================
      inpt_dimlist = DF.dims
      # first check if dimensions are correct and don't need to be modified
      if 'scalar_axis' not in inpt_dimlist:           # first see if scalar axis is a dimension
            scalar_axis = DF.assign_coords(scalar_axis=1)
      if DF.areo.dims != ('time',scalar_axis):
            DF['areo'] = DF.areo.expand_dims('scalar_axis', axis=1)




      #==================================================================
      # Output Processed Data to New NC File
      #==================================================================
      DF.to_netcdf(fullnameOUT)
      prCyan(fullnameOUT +' was created')

# ======================================================
#                  END OF PROGRAM
# ======================================================

if __name__ == '__main__':
    main()
