#!/usr/bin/env python3

#Load generic Python Modules

import argparse   # parse arguments
import numpy as np
import xarray as xr
import os
#---Use in-script function for now---
from amesgcm.FV3_utils import layers_mid_point_to_boundary
from amesgcm.Script_utils import prCyan
#---

# Routine to Transform Model Input (variable names, dimension names, array order)
# to expected configuration CAP

parser = argparse.ArgumentParser(description="""\033[93m openMars2FV3.py  Used to convert openMars output to FV3 format  \n \033[00m""",
                                formatter_class=argparse.RawTextHelpFormatter)


parser.add_argument('input_file', nargs='+', #sys.stdin
                             help='***.nc file or list of ***.nc files ')



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


        # keep attributes of variables
        xr.set_options(keep_attrs=True)

        DS = xr.open_dataset(fullnameIN, decode_times=False)

        # Make a list of variables and dimensions of the input
        inpt_varlist = list(DS.keys())
        inpt_dimlist = list(DS.coords)

        #===================================================================
        # change dimension names in input file to the expected names for FV3
        #===================================================================
        # list of possible dimensionn name variations, always set FV3 name as first element in list
        dimvar_list =  [['pfull', 'lev', 'level', 'pstd','dimvert'],['lat','lats','latitude','dimlat'],['lon','lons','longitude','dimlon'],['time','sols']]

        # Make a dictionary of the input dimension name and the expected dimension name for CAP
        dims_dict = {}
        for e in dimvar_list:
            for f in inpt_dimlist:
                if f in e[1:]:   # skip first element so we don't try to rename dimensions that are already GTG
                    dims_dict[f]=e[0]

        # Rename Dimensions/Coordinates
        DS = DS.rename_dims(dims_dict=dims_dict) # keeps all the dimension attributes
        DS = DS.rename(dims_dict)

        #==================================================================
        # check that vertical grid starts at toa with highest level at surface
        #==================================================================
        if DS.pfull[0] != DS.pfull.min(): # if toa, lev = 0 is surface then flip
            DS=DS.isel(pfull=slice(None, None, -1)) # regrids DS based on pfull

        #==================================================================
        # change variable names in input file to the expected names for FV3
        #==================================================================

        # list of possible variable name variations, always set FV3 name as first element in list
        varvar_list = [['pfull','lev'],['lat','latitude','latitudes'],['lon','lons','longitude','longitudes'],['areo',         'Ls'],['ucomp','u'],['vcomp','v'],['ts','tsurf'],['dust_mass_col','dustcol'],['frost','cos2ice']]

        # make dictionary key=input variable name, value = preferred variable name
        var_dict = {}
        for e in varvar_list:
            for f in inpt_varlist:
                if f in e[1:]:   # skip first element so we don't try to rename dimensions that are already GTG
                    var_dict[f]=e[0]

        # rename variables
        DS = DS. rename_vars(name_dict = var_dict)

        #==================================================================
        # reorder dimensions
        #==================================================================
        DS = DS.transpose("time", "pfull", "lat","lon")

        #==================================================================
        # change longitude from -180-179 to 0-360
        #==================================================================
        if min(DS.lon)<0:
        	tmp = np.array(DS.lon)
        	tmp = np.where(tmp<0,tmp+360,tmp)
        	DS=DS.assign_coords({'lon':('lon',tmp,DS.lon.attrs)})
        	DS = DS.sortby("lon")

        #==================================================================
        # add scalar axis to areo [time, scalar_axis])
        #==================================================================
        # first check if dimensions are correct and don't need to be modified
        if 'scalar_axis' not in inpt_dimlist:		# first see if scalar axis is a dimension
        	scalar_axis = DS.assign_coords(scalar_axis=1)
        if DS.areo.dims != ('time',scalar_axis):
        	DS['areo'] = DS.areo.expand_dims('scalar_axis', axis=1)


        #==================================================================
        # OUTDATED (no longer needed after update)
        # Check for All Necessary Attributes
        #==================================================================
        new_varlist = list(DS.keys())
        new_dimlist = list(DS.coords)
        attrs_list = list(DS.attrs)
        if 'long_name' not in attrs_list:
        	for i in new_varlist:
        		DS[i].attrs['long_name'] = DS[i].attrs['FIELDNAM']
        	for i in new_dimlist:
        		DS[i].attrs['long_name'] = DS[i].attrs['FIELDNAM']
        if 'units' not in attrs_list:
        	for i in new_varlist:
        		DS[i].attrs['units'] = DS[i].attrs['UNITS']
        	for i in new_dimlist:
        		DS[i].attrs['units'] = DS[i].attrs['UNITS']


        #==================================================================
        # PlACEHOOLDER: check units of vertical grid and time
        #==================================================================
        #if DS['pfull'].attrs['units'] != 'Pa':

        #==================================================================
        # check if vertical dimension is interpolated pressure, native
        # grid or sigma levels
        #==================================================================
        if DS['pfull'].attrs['FIELDNAM'][:5]=='sigma' or DS['pfull'].attrs['FIELDNAM'][:5]=='Sigma':
                DS['bk'] = layers_mid_point_to_boundary(DS.pfull,1.) #Compute the bk_half (p_half) from sigma (pfull)
                DS['pk'] = np.zeros(len(DS.pfull)+1) #Pure sigma model, this is zero
                DS['pfull']=DS.pfull*610 #TODO set to 610Pa, for consistency should be set to the model

        #==================================================================
        # Output Processed Data to New NC File
        #==================================================================
        DS.to_netcdf(fullnameOUT)
        prCyan(fullnameOUT + ' was created')
        #==================================================================
        # Add Dummy Fixed File if Necessary
        #==================================================================
        #if os.path.exists(str(path_inpt+'/'+filename[:-3]+'.fixed.nc')):
        #	exit
        #else:
        #	tmp = DS.pk.values
        #	var1 = xr.DataArray(tmp, coords={'pfull': tmp},dims=['pfull'])
        #	var1.name = 'pk'
        #
        #	tmp = DS.bk.values
        #	var2 = xr.DataArray(tmp, coords={'pfull': tmp},dims=['pfull'])
        #	var2.name = 'bk'
        #	NEW_DF = xr.merge([var1,var2])
        #	NEW_DF.to_netcdf(str(path_inpt+'/'+filename[:-3]+'.fixed.nc'))
        #f = open(path+'/'+filename[:-3]+'.fixed.nc', "a")
        #f.close()


if __name__ == '__main__':
    main()
