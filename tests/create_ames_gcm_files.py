#!/usr/bin/env python3
"""
Script to create test NetCDF files for Mars Global Climate Model (MGCM) data.
This script generates files with variables that exactly match the specifications
in the mgcm_contents.txt file.
"""

import numpy as np
from netCDF4 import Dataset
import sys
import os

def create_mgcm_fixed():
    """Create fixed.nc with the exact variables and structure as specified."""
    nc_file = Dataset('01336.fixed.nc', 'w', format='NETCDF4')
    
    # Define dimensions based on the provided data
    lat_dim = nc_file.createDimension('lat', 48)
    lon_dim = nc_file.createDimension('lon', 96)
    phalf_dim = nc_file.createDimension('phalf', 31)
    bnds_dim = nc_file.createDimension('bnds', 2)
    
    # Create and populate lat variable
    lat_var = nc_file.createVariable('lat', 'f4', ('lat',))
    lat_var.long_name = 'latitude'
    lat_var.units = 'degrees_N'
    lat_values = np.array([-88.125, -84.375, -80.625, -76.875, -73.125, -69.375, -65.625, -61.875, -58.125,
                         -54.375, -50.625, -46.875, -43.125, -39.375, -35.625, -31.875, -28.125, -24.375,
                         -20.625, -16.875, -13.125, -9.375, -5.625, -1.875, 1.875, 5.625, 9.375,
                         13.125, 16.875, 20.625, 24.375, 28.125, 31.875, 35.625, 39.375, 43.125,
                         46.875, 50.625, 54.375, 58.125, 61.875, 65.625, 69.375, 73.125, 76.875,
                         80.625, 84.375, 88.125])
    lat_var[:] = lat_values
    
    # Create and populate grid_yt_bnds variable
    grid_yt_bnds_var = nc_file.createVariable('grid_yt_bnds', 'f4', ('lat', 'bnds'))
    grid_yt_bnds_values = np.array([
        [-90.0, -86.25], [-86.25, -82.5], [-82.5, -78.75], [-78.75, -75.0],
        [-75.0, -71.25], [-71.25, -67.5], [-67.5, -63.75], [-63.75, -60.0],
        [-60.0, -56.25], [-56.25, -52.5], [-52.5, -48.75], [-48.75, -45.0],
        [-45.0, -41.25], [-41.25, -37.5], [-37.5, -33.75], [-33.75, -30.0],
        [-30.0, -26.25], [-26.25, -22.5], [-22.5, -18.75], [-18.75, -15.0],
        [-15.0, -11.25], [-11.25, -7.5], [-7.5, -3.75], [-3.75, 0.0],
        [0.0, 3.75], [3.75, 7.5], [7.5, 11.25], [11.25, 15.0],
        [15.0, 18.75], [18.75, 22.5], [22.5, 26.25], [26.25, 30.0],
        [30.0, 33.75], [33.75, 37.5], [37.5, 41.25], [41.25, 45.0],
        [45.0, 48.75], [48.75, 52.5], [52.5, 56.25], [56.25, 60.0],
        [60.0, 63.75], [63.75, 67.5], [67.5, 71.25], [71.25, 75.0],
        [75.0, 78.75], [78.75, 82.5], [82.5, 86.25], [86.25, 90.0]
    ])
    grid_yt_bnds_var[:] = grid_yt_bnds_values
    grid_yt_bnds_var.long_name = 'T-cell latitude'
    grid_yt_bnds_var.units = 'degrees_N'
    
    # Create and populate lon variable
    lon_var = nc_file.createVariable('lon', 'f4', ('lon',))
    lon_var.long_name = 'longitude'
    lon_var.units = 'degrees_E'
    lon_values = np.array([1.875, 5.625, 9.375, 13.125, 16.875, 20.625, 24.375, 28.125, 31.875,
                         35.625, 39.375, 43.125, 46.875, 50.625, 54.375, 58.125, 61.875, 65.625,
                         69.375, 73.125, 76.875, 80.625, 84.375, 88.125, 91.875, 95.625, 99.375,
                         103.125, 106.875, 110.625, 114.375, 118.125, 121.875, 125.625, 129.375, 133.125,
                         136.875, 140.625, 144.375, 148.125, 151.875, 155.625, 159.375, 163.125, 166.875,
                         170.625, 174.375, 178.125, 181.875, 185.625, 189.375, 193.125, 196.875, 200.625,
                         204.375, 208.125, 211.875, 215.625, 219.375, 223.125, 226.875, 230.625, 234.375,
                         238.125, 241.875, 245.625, 249.375, 253.125, 256.875, 260.625, 264.375, 268.125,
                         271.875, 275.625, 279.375, 283.125, 286.875, 290.625, 294.375, 298.125, 301.875,
                         305.625, 309.375, 313.125, 316.875, 320.625, 324.375, 328.125, 331.875, 335.625,
                         339.375, 343.125, 346.875, 350.625, 354.375, 358.125])
    lon_var[:] = lon_values
    
    # Create and populate grid_xt_bnds variable
    grid_xt_bnds_var = nc_file.createVariable('grid_xt_bnds', 'f4', ('lon', 'bnds'))
    grid_xt_bnds_values = np.array([
        [0.0, 3.75], [3.75, 7.5], [7.5, 11.25], [11.25, 15.0],
        [15.0, 18.75], [18.75, 22.5], [22.5, 26.25], [26.25, 30.0],
        [30.0, 33.75], [33.75, 37.5], [37.5, 41.25], [41.25, 45.0],
        [45.0, 48.75], [48.75, 52.5], [52.5, 56.25], [56.25, 60.0],
        [60.0, 63.75], [63.75, 67.5], [67.5, 71.25], [71.25, 75.0],
        [75.0, 78.75], [78.75, 82.5], [82.5, 86.25], [86.25, 90.0],
        [90.0, 93.75], [93.75, 97.5], [97.5, 101.25], [101.25, 105.0],
        [105.0, 108.75], [108.75, 112.5], [112.5, 116.25], [116.25, 120.0],
        [120.0, 123.75], [123.75, 127.5], [127.5, 131.25], [131.25, 135.0],
        [135.0, 138.75], [138.75, 142.5], [142.5, 146.25], [146.25, 150.0],
        [150.0, 153.75], [153.75, 157.5], [157.5, 161.25], [161.25, 165.0],
        [165.0, 168.75], [168.75, 172.5], [172.5, 176.25], [176.25, 180.0],
        [180.0, 183.75], [183.75, 187.5], [187.5, 191.25], [191.25, 195.0],
        [195.0, 198.75], [198.75, 202.5], [202.5, 206.25], [206.25, 210.0],
        [210.0, 213.75], [213.75, 217.5], [217.5, 221.25], [221.25, 225.0],
        [225.0, 228.75], [228.75, 232.5], [232.5, 236.25], [236.25, 240.0],
        [240.0, 243.75], [243.75, 247.5], [247.5, 251.25], [251.25, 255.0],
        [255.0, 258.75], [258.75, 262.5], [262.5, 266.25], [266.25, 270.0],
        [270.0, 273.75], [273.75, 277.5], [277.5, 281.25], [281.25, 285.0],
        [285.0, 288.75], [288.75, 292.5], [292.5, 296.25], [296.25, 300.0],
        [300.0, 303.75], [303.75, 307.5], [307.5, 311.25], [311.25, 315.0],
        [315.0, 318.75], [318.75, 322.5], [322.5, 326.25], [326.25, 330.0],
        [330.0, 333.75], [333.75, 337.5], [337.5, 341.25], [341.25, 345.0],
        [345.0, 348.75], [348.75, 352.5], [352.5, 356.25], [356.25, 360.0]
    ])
    grid_xt_bnds_var[:] = grid_xt_bnds_values
    grid_xt_bnds_var.long_name = 'T-cell longitude'
    grid_xt_bnds_var.units = 'degrees_E'
    
    # Create and populate other variables
    zsurf_var = nc_file.createVariable('zsurf', 'f4', ('lat', 'lon'))
    zsurf_var.long_name = 'surface height'
    zsurf_var.units = 'm'
    zsurf_var[:] = np.random.uniform(-7.1e+03, 1.1e+04, size=(48, 96))
    
    thin_var = nc_file.createVariable('thin', 'f4', ('lat', 'lon'))
    thin_var.long_name = 'Surface Thermal Inertia'
    thin_var.units = 'mks'
    thin_var[:] = np.random.uniform(40.6, 1037.6, size=(48, 96))
    
    alb_var = nc_file.createVariable('alb', 'f4', ('lat', 'lon'))
    alb_var.long_name = 'Surface Albedo'
    alb_var.units = 'mks'
    alb_var[:] = np.random.uniform(0.1, 0.3, size=(48, 96))
    
    emis_var = nc_file.createVariable('emis', 'f4', ('lat', 'lon'))
    emis_var.long_name = 'Surface Emissivity'
    emis_var[:] = np.random.uniform(0.9, 1.0, size=(48, 96))
    
    gice_var = nc_file.createVariable('gice', 'f4', ('lat', 'lon'))
    gice_var.long_name = 'GRS Ice'
    gice_var[:] = np.random.uniform(-57.9, 58.6, size=(48, 96))
    
    bk_var = nc_file.createVariable('bk', 'f4', ('phalf',))
    bk_var.long_name = 'vertical coordinate sigma value'
    bk_values = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.00193664, 0.00744191, 0.01622727, 
                          0.02707519, 0.043641, 0.0681068, 0.1028024, 0.14971954, 
                          0.20987134, 0.28270233, 0.3658161, 0.4552023, 0.545936, 
                          0.6331097, 0.7126763, 0.7819615, 0.8397753, 0.88620347, 
                          0.9222317, 0.94934535, 0.9691962, 0.98337257, 0.9932694, 
                          0.996, 0.999, 1.0])
    bk_var[:] = bk_values
    
    phalf_var = nc_file.createVariable('phalf', 'f4', ('phalf',))
    phalf_var.long_name = 'ref half pressure level'
    phalf_var.units = 'mb'
    phalf_values = np.array([1.94482759e-04, 5.57983413e-04, 1.90437332e-03, 5.75956606e-03,
                            1.52282217e-02, 3.74336530e-02, 7.93855540e-02, 1.42458016e-01,
                            2.19247580e-01, 3.35953688e-01, 5.07962971e-01, 7.51680775e-01,
                            1.08112355e+00, 1.50342598e+00, 2.01470398e+00, 2.59814516e+00,
                            3.22560498e+00, 3.86251679e+00, 4.47443525e+00, 5.03295321e+00,
                            5.51929991e+00, 5.92512245e+00, 6.25102343e+00, 6.50392232e+00,
                            6.69424554e+00, 6.83358777e+00, 6.93309849e+00, 7.00256879e+00,
                            7.02180194e+00, 7.04295002e+00, 7.05000000e+00])
    phalf_var[:] = phalf_values
    
    pk_var = nc_file.createVariable('pk', 'f4', ('phalf',))
    pk_var.long_name = 'pressure part of the hybrid coordinate'
    pk_var.units = 'Pa'
    pk_values = np.array([1.9448277e-02, 5.5798341e-02, 1.9043733e-01, 5.7595658e-01, 1.5228221e+00,
                          2.3780346e+00, 2.6920066e+00, 2.8055782e+00, 2.8367476e+00, 2.8284638e+00,
                          2.7810004e+00, 2.6923854e+00, 2.5600791e+00, 2.3833103e+00, 2.1652553e+00,
                          1.9141655e+00, 1.6428767e+00, 1.3668065e+00, 1.1011865e+00, 8.5853612e-01,
                          6.4712650e-01, 4.7065818e-01, 3.2891041e-01, 2.1889719e-01, 1.3609630e-01,
                          7.5470544e-02, 3.2172799e-02, 1.9448276e-03, 1.9448275e-04, 1.9448275e-06,
                          0.0000000e+00])
    pk_var[:] = pk_values
    
    nc_file.close()
    print("Created 01336.fixed.nc")

def create_mgcm_atmos_average(short=False):
    """Create atmos_average.nc with the exact variables and structure as specified."""
    nc_file = Dataset('01336.atmos_average.nc', 'w', format='NETCDF4')

    # Shorten file length if wanted
    if short:
        len_time = 2
    else:
        len_time = 133

    # Define dimensions
    time_dim = nc_file.createDimension('time', len_time)
    lat_dim = nc_file.createDimension('lat', 48)
    lon_dim = nc_file.createDimension('lon', 96)
    pfull_dim = nc_file.createDimension('pfull', 30)
    phalf_dim = nc_file.createDimension('phalf', 31)
    scalar_axis_dim = nc_file.createDimension('scalar_axis', 1)
    
    # Create and populate variables
    
    # Time, lat, lon, scalar_axis variables
    time_var = nc_file.createVariable('time', 'f4', ('time',))
    time_var.long_name = 'time'
    time_var.units = 'days'
    time_var[:] = np.linspace(1338.5, 1998.5, len_time)
    
    lat_var = nc_file.createVariable('lat', 'f4', ('lat',))
    lat_var.long_name = 'latitude'
    lat_var.units = 'degrees_N'
    lat_var[:] = np.array([-88.125, -84.375, -80.625, -76.875, -73.125, -69.375, -65.625, -61.875, -58.125,
                         -54.375, -50.625, -46.875, -43.125, -39.375, -35.625, -31.875, -28.125, -24.375,
                         -20.625, -16.875, -13.125, -9.375, -5.625, -1.875, 1.875, 5.625, 9.375,
                         13.125, 16.875, 20.625, 24.375, 28.125, 31.875, 35.625, 39.375, 43.125,
                         46.875, 50.625, 54.375, 58.125, 61.875, 65.625, 69.375, 73.125, 76.875,
                         80.625, 84.375, 88.125])
    
    lon_var = nc_file.createVariable('lon', 'f4', ('lon',))
    lon_var.long_name = 'longitude'
    lon_var.units = 'degrees_E'
    lon_var[:] = np.array([1.875, 5.625, 9.375, 13.125, 16.875, 20.625, 24.375, 28.125, 31.875,
                         35.625, 39.375, 43.125, 46.875, 50.625, 54.375, 58.125, 61.875, 65.625,
                         69.375, 73.125, 76.875, 80.625, 84.375, 88.125, 91.875, 95.625, 99.375,
                         103.125, 106.875, 110.625, 114.375, 118.125, 121.875, 125.625, 129.375, 133.125,
                         136.875, 140.625, 144.375, 148.125, 151.875, 155.625, 159.375, 163.125, 166.875,
                         170.625, 174.375, 178.125, 181.875, 185.625, 189.375, 193.125, 196.875, 200.625,
                         204.375, 208.125, 211.875, 215.625, 219.375, 223.125, 226.875, 230.625, 234.375,
                         238.125, 241.875, 245.625, 249.375, 253.125, 256.875, 260.625, 264.375, 268.125,
                         271.875, 275.625, 279.375, 283.125, 286.875, 290.625, 294.375, 298.125, 301.875,
                         305.625, 309.375, 313.125, 316.875, 320.625, 324.375, 328.125, 331.875, 335.625,
                         339.375, 343.125, 346.875, 350.625, 354.375, 358.125])
    
    scalar_axis_var = nc_file.createVariable('scalar_axis', 'f4', ('scalar_axis',))
    scalar_axis_var.long_name = 'none'
    scalar_axis_var[:] = np.array([0.0])
    
    # Pfull, phalf, bk, pk variables
    pfull_var = nc_file.createVariable('pfull', 'f4', ('pfull',))
    pfull_var.long_name = 'ref full pressure level'
    pfull_var.units = 'mb'
    pfull_values = np.array([3.44881953e-04, 1.09678471e-03, 3.48347419e-03, 9.73852715e-03,
                           2.46886197e-02, 5.58059295e-02, 1.07865789e-01, 1.78102296e-01,
                           2.73462607e-01, 4.16048919e-01, 6.21882691e-01, 9.06446215e-01,
                           1.28069137e+00, 1.74661073e+00, 2.29407255e+00, 2.90057273e+00,
                           3.53450185e+00, 4.16097960e+00, 4.74822077e+00, 5.27238854e+00,
                           5.71981194e+00, 6.08661884e+00, 6.37663706e+00, 6.59862648e+00,
                           6.76367744e+00, 6.88322325e+00, 6.96777592e+00, 7.01218097e+00,
                           7.03237068e+00, 7.04647442e+00])
    pfull_var[:] = pfull_values
    
    phalf_var = nc_file.createVariable('phalf', 'f4', ('phalf',))
    phalf_var.long_name = 'ref half pressure level'
    phalf_var.units = 'mb'
    phalf_values = np.array([1.94482759e-04, 5.57983413e-04, 1.90437332e-03, 5.75956606e-03,
                           1.52282217e-02, 3.74336530e-02, 7.93855540e-02, 1.42458016e-01,
                           2.19247580e-01, 3.35953688e-01, 5.07962971e-01, 7.51680775e-01,
                           1.08112355e+00, 1.50342598e+00, 2.01470398e+00, 2.59814516e+00,
                           3.22560498e+00, 3.86251679e+00, 4.47443525e+00, 5.03295321e+00,
                           5.51929991e+00, 5.92512245e+00, 6.25102343e+00, 6.50392232e+00,
                           6.69424554e+00, 6.83358777e+00, 6.93309849e+00, 7.00256879e+00,
                           7.02180194e+00, 7.04295002e+00, 7.05000000e+00])
    phalf_var[:] = phalf_values
    
    bk_var = nc_file.createVariable('bk', 'f4', ('phalf',))
    bk_var.long_name = 'vertical coordinate sigma value'
    bk_values = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.00193664, 0.00744191, 0.01622727, 
                         0.02707519, 0.043641, 0.0681068, 0.1028024, 0.14971954, 
                         0.20987134, 0.28270233, 0.3658161, 0.4552023, 0.545936, 
                         0.6331097, 0.7126763, 0.7819615, 0.8397753, 0.88620347, 
                         0.9222317, 0.94934535, 0.9691962, 0.98337257, 0.9932694, 
                         0.996, 0.999, 1.0])
    bk_var[:] = bk_values
    
    pk_var = nc_file.createVariable('pk', 'f4', ('phalf',))
    pk_var.long_name = 'pressure part of the hybrid coordinate'
    pk_var.units = 'Pa'
    pk_values = np.array([1.9448277e-02, 5.5798341e-02, 1.9043733e-01, 5.7595658e-01, 1.5228221e+00,
                         2.3780346e+00, 2.6920066e+00, 2.8055782e+00, 2.8367476e+00, 2.8284638e+00,
                         2.7810004e+00, 2.6923854e+00, 2.5600791e+00, 2.3833103e+00, 2.1652553e+00,
                         1.9141655e+00, 1.6428767e+00, 1.3668065e+00, 1.1011865e+00, 8.5853612e-01,
                         6.4712650e-01, 4.7065818e-01, 3.2891041e-01, 2.1889719e-01, 1.3609630e-01,
                         7.5470544e-02, 3.2172799e-02, 1.9448276e-03, 1.9448275e-04, 1.9448275e-06,
                         0.0000000e+00])
    pk_var[:] = pk_values
    
    # Other variables specific to atmos_average.nc
    areo_var = nc_file.createVariable('areo', 'f4', ('time', 'scalar_axis'))
    areo_var.long_name = 'areo'
    areo_var.units = 'degrees'
    areo_vals = np.linspace(722.2, 1076.5, len_time)
    areo_data = np.zeros((len_time, 1))  # Create a 2D array with shape (len_time, 1)
    for i in range(len_time):
        areo_data[i, 0] = areo_vals[i]
    areo_var[:] = areo_data
    
    cldcol_var = nc_file.createVariable('cldcol', 'f4', ('time', 'lat', 'lon'))
    cldcol_var.long_name = 'ice column'
    cldcol_var[:] = np.random.uniform(1.2e-11, 4.1e-02, size=(len_time, 48, 96))
    
    dst_mass_micro_var = nc_file.createVariable('dst_mass_micro', 'f4', ('time', 'pfull', 'lat', 'lon'))
    dst_mass_micro_var.long_name = 'dust_mass'
    dst_mass_micro_var[:] = np.random.uniform(1.5e-17, 2.5e-04, size=(len_time, 30, 48, 96))
    
    dst_num_micro_var = nc_file.createVariable('dst_num_micro', 'f4', ('time', 'pfull', 'lat', 'lon'))
    dst_num_micro_var.long_name = 'dust_number'
    dst_num_micro_var[:] = np.random.uniform(-3.8e-15, 6.3e+10, size=(133, 30, 48, 96))
    
    ice_mass_micro_var = nc_file.createVariable('ice_mass_micro', 'f4', ('time', 'pfull', 'lat', 'lon'))
    ice_mass_micro_var.long_name = 'ice_mass'
    ice_mass_micro_var[:] = np.random.uniform(-5.8e-34, 3.1e-03, size=(133, 30, 48, 96))
    
    omega_var = nc_file.createVariable('omega', 'f4', ('time', 'pfull', 'lat', 'lon'))
    omega_var.long_name = 'vertical wind'
    omega_var.units = 'Pa/s'
    omega_var[:] = np.random.uniform(-0.045597, 0.0806756, size=(133, 30, 48, 96))
    
    ps_var = nc_file.createVariable('ps', 'f4', ('time', 'lat', 'lon'))
    ps_var.long_name = 'surface pressure'
    ps_var.units = 'Pa'
    ps_var[:] = np.random.uniform(176.8, 1318.8, size=(len_time, 48, 96))
    
    r_var = nc_file.createVariable('r', 'f4', ('time', 'pfull', 'lat', 'lon'))
    r_var.long_name = 'specific humidity'
    r_var.units = 'kg/kg'
    r_var[:] = np.random.uniform(8.6e-13, 3.4e-03, size=(len_time, 30, 48, 96))
    
    taudust_IR_var = nc_file.createVariable('taudust_IR', 'f4', ('time', 'lat', 'lon'))
    taudust_IR_var.long_name = 'Dust opacity IR'
    taudust_IR_var.units = 'op'
    taudust_IR_var[:] = np.random.uniform(0.0, 0.5, size=(len_time, 48, 96))
    
    temp_var = nc_file.createVariable('temp', 'f4', ('time', 'pfull', 'lat', 'lon'))
    temp_var.long_name = 'temperature'
    temp_var.units = 'K'
    temp_var[:] = np.random.uniform(104.1, 258.8, size=(len_time, 30, 48, 96))
    
    ts_var = nc_file.createVariable('ts', 'f4', ('time', 'lat', 'lon'))
    ts_var.long_name = 'Surface Temperature'
    ts_var.units = 'K'
    ts_var[:] = np.random.uniform(143.4, 258.7, size=(len_time, 48, 96))
    
    ucomp_var = nc_file.createVariable('ucomp', 'f4', ('time', 'pfull', 'lat', 'lon'))
    ucomp_var.long_name = 'zonal wind'
    ucomp_var.units = 'm/sec'
    ucomp_var[:] = np.random.uniform(-268.7, 212.7, size=(len_time, 30, 48, 96))
    
    vcomp_var = nc_file.createVariable('vcomp', 'f4', ('time', 'pfull', 'lat', 'lon'))
    vcomp_var.long_name = 'meridional wind'
    vcomp_var.units = 'm/sec'
    vcomp_var[:] = np.random.uniform(-97.5, 109.6, size=(len_time, 30, 48, 96))
    
    nc_file.close()
    print("Created 01336.atmos_average.nc")

def create_mgcm_atmos_daily(short=False):
    """Create atmos_daily.nc with the exact variables and structure as specified."""
    nc_file = Dataset('01336.atmos_daily.nc', 'w', format='NETCDF4')

    # Shorten file length if wanted
    if short:
        len_time = 2
    else:
        len_time = 2672

    # Define dimensions
    time_dim = nc_file.createDimension('time', len_time)
    lat_dim = nc_file.createDimension('lat', 48)
    lon_dim = nc_file.createDimension('lon', 96)
    pfull_dim = nc_file.createDimension('pfull', 30)
    scalar_axis_dim = nc_file.createDimension('scalar_axis', 1)
    
    # Create variables
    time_var = nc_file.createVariable('time', 'f4', ('time',))
    time_var.long_name = 'time'
    time_var.units = 'days'
    time_var[:] = np.linspace(1336.2, 2004.0, len_time)
    
    lat_var = nc_file.createVariable('lat', 'f4', ('lat',))
    lat_var.long_name = 'latitude'
    lat_var.units = 'degrees_N'
    lat_var[:] = np.array([-88.125, -84.375, -80.625, -76.875, -73.125, -69.375, -65.625, -61.875, -58.125,
                          -54.375, -50.625, -46.875, -43.125, -39.375, -35.625, -31.875, -28.125, -24.375,
                          -20.625, -16.875, -13.125, -9.375, -5.625, -1.875, 1.875, 5.625, 9.375,
                          13.125, 16.875, 20.625, 24.375, 28.125, 31.875, 35.625, 39.375, 43.125,
                          46.875, 50.625, 54.375, 58.125, 61.875, 65.625, 69.375, 73.125, 76.875,
                          80.625, 84.375, 88.125])
    
    lon_var = nc_file.createVariable('lon', 'f4', ('lon',))
    lon_var.long_name = 'longitude'
    lon_var.units = 'degrees_E'
    lon_var[:] = np.array([1.875, 5.625, 9.375, 13.125, 16.875, 20.625, 24.375, 28.125, 31.875,
                          35.625, 39.375, 43.125, 46.875, 50.625, 54.375, 58.125, 61.875, 65.625,
                          69.375, 73.125, 76.875, 80.625, 84.375, 88.125, 91.875, 95.625, 99.375,
                          103.125, 106.875, 110.625, 114.375, 118.125, 121.875, 125.625, 129.375, 133.125,
                          136.875, 140.625, 144.375, 148.125, 151.875, 155.625, 159.375, 163.125, 166.875,
                          170.625, 174.375, 178.125, 181.875, 185.625, 189.375, 193.125, 196.875, 200.625,
                          204.375, 208.125, 211.875, 215.625, 219.375, 223.125, 226.875, 230.625, 234.375,
                          238.125, 241.875, 245.625, 249.375, 253.125, 256.875, 260.625, 264.375, 268.125,
                          271.875, 275.625, 279.375, 283.125, 286.875, 290.625, 294.375, 298.125, 301.875,
                          305.625, 309.375, 313.125, 316.875, 320.625, 324.375, 328.125, 331.875, 335.625,
                          339.375, 343.125, 346.875, 350.625, 354.375, 358.125])
    
    scalar_axis_var = nc_file.createVariable('scalar_axis', 'f4', ('scalar_axis',))
    scalar_axis_var.long_name = 'none'
    scalar_axis_var[:] = np.array([0.0])
    
    pfull_var = nc_file.createVariable('pfull', 'f4', ('pfull',))
    pfull_var.long_name = 'ref full pressure level'
    pfull_var.units = 'mb'
    pfull_values = np.array([3.44881953e-04, 1.09678471e-03, 3.48347419e-03, 9.73852715e-03,
                            2.46886197e-02, 5.58059295e-02, 1.07865789e-01, 1.78102296e-01,
                            2.73462607e-01, 4.16048919e-01, 6.21882691e-01, 9.06446215e-01,
                            1.28069137e+00, 1.74661073e+00, 2.29407255e+00, 2.90057273e+00,
                            3.53450185e+00, 4.16097960e+00, 4.74822077e+00, 5.27238854e+00,
                            5.71981194e+00, 6.08661884e+00, 6.37663706e+00, 6.59862648e+00,
                            6.76367744e+00, 6.88322325e+00, 6.96777592e+00, 7.01218097e+00,
                            7.03237068e+00, 7.04647442e+00])
    pfull_var[:] = pfull_values
    
    # Add specific variables for atmos_daily.nc
    areo_var = nc_file.createVariable('areo', 'f4', ('time', 'scalar_axis'))
    areo_var.long_name = 'areo'
    areo_var.units = 'degrees'
    areo_vals = np.linspace(720.3, 1079.8, len_time)
    areo_data = np.zeros((len_time, 1))  # Create a 2D array with shape (133, 1)
    for i in range(len_time):
        areo_data[i, 0] = areo_vals[i]
    areo_var[:] = areo_data
    
    ps_var = nc_file.createVariable('ps', 'f4', ('time', 'lat', 'lon'))
    ps_var.long_name = 'surface pressure'
    ps_var.units = 'Pa'
    ps_var[:] = np.random.uniform(170.3, 1340.2, size=(len_time, 48, 96))
    
    temp_var = nc_file.createVariable('temp', 'f4', ('time', 'pfull', 'lat', 'lon'))
    temp_var.long_name = 'temperature'
    temp_var.units = 'K'
    temp_var[:] = np.random.uniform(101.6, 283.9, size=(len_time, 30, 48, 96))
    
    nc_file.close()
    print("Created 01336.atmos_daily.nc")

def create_mgcm_atmos_average_pstd(short=False):
    """Create atmos_average_pstd.nc with the exact variables and structure as specified."""
    nc_file = Dataset('01336.atmos_average_pstd.nc', 'w', format='NETCDF4')

    # Shorten file length if wanted
    if short:
        len_time = 2
    else:
        len_time = 133
    
    # Define dimensions
    time_dim = nc_file.createDimension('time', len_time)
    pstd_dim = nc_file.createDimension('pstd', 44)
    lat_dim = nc_file.createDimension('lat', 48)
    lon_dim = nc_file.createDimension('lon', 96)
    scalar_axis_dim = nc_file.createDimension('scalar_axis', 1)
    
    # Create and populate key variables
    pstd_var = nc_file.createVariable('pstd', 'f4', ('pstd',))
    pstd_var.long_name = 'standard pressure'
    pstd_var.units = 'Pa'
    # Using exactly the pstd values from the file
    pstd_values = np.array([1.0e-05, 3.0e-05, 5.0e-05, 1.0e-04, 3.0e-04, 5.0e-04, 3.0e-03, 5.0e-03, 1.0e-02,
                           3.0e-02, 5.0e-02, 1.0e-01, 2.0e-01, 3.0e-01, 5.0e-01, 1.0e+00, 2.0e+00, 3.0e+00,
                           5.0e+00, 7.0e+00, 1.0e+01, 2.0e+01, 3.0e+01, 5.0e+01, 7.0e+01, 1.0e+02, 1.5e+02,
                           2.0e+02, 2.5e+02, 3.0e+02, 3.5e+02, 4.0e+02, 4.5e+02, 5.0e+02, 5.3e+02, 5.5e+02,
                           5.9e+02, 6.0e+02, 6.3e+02, 6.5e+02, 6.9e+02, 7.0e+02, 7.5e+02, 8.0e+02])
    pstd_var[:] = pstd_values
    
    lon_var = nc_file.createVariable('lon', 'f4', ('lon',))
    lon_var.long_name = 'longitude'
    lon_var.units = 'degrees_E'
    lon_var[:] = np.array([1.875, 5.625, 9.375, 13.125, 16.875, 20.625, 24.375, 28.125, 31.875,
                          35.625, 39.375, 43.125, 46.875, 50.625, 54.375, 58.125, 61.875, 65.625,
                          69.375, 73.125, 76.875, 80.625, 84.375, 88.125, 91.875, 95.625, 99.375,
                          103.125, 106.875, 110.625, 114.375, 118.125, 121.875, 125.625, 129.375, 133.125,
                          136.875, 140.625, 144.375, 148.125, 151.875, 155.625, 159.375, 163.125, 166.875,
                          170.625, 174.375, 178.125, 181.875, 185.625, 189.375, 193.125, 196.875, 200.625,
                          204.375, 208.125, 211.875, 215.625, 219.375, 223.125, 226.875, 230.625, 234.375,
                          238.125, 241.875, 245.625, 249.375, 253.125, 256.875, 260.625, 264.375, 268.125,
                          271.875, 275.625, 279.375, 283.125, 286.875, 290.625, 294.375, 298.125, 301.875,
                          305.625, 309.375, 313.125, 316.875, 320.625, 324.375, 328.125, 331.875, 335.625,
                          339.375, 343.125, 346.875, 350.625, 354.375, 358.125])
    
    lat_var = nc_file.createVariable('lat', 'f4', ('lat',))
    lat_var.long_name = 'latitude'
    lat_var.units = 'degrees_N'
    lat_var[:] = np.array([-88.125, -84.375, -80.625, -76.875, -73.125, -69.375, -65.625, -61.875, -58.125,
                          -54.375, -50.625, -46.875, -43.125, -39.375, -35.625, -31.875, -28.125, -24.375,
                          -20.625, -16.875, -13.125, -9.375, -5.625, -1.875, 1.875, 5.625, 9.375,
                          13.125, 16.875, 20.625, 24.375, 28.125, 31.875, 35.625, 39.375, 43.125,
                          46.875, 50.625, 54.375, 58.125, 61.875, 65.625, 69.375, 73.125, 76.875,
                          80.625, 84.375, 88.125])
    
    time_var = nc_file.createVariable('time', 'f4', ('time',))
    time_var.long_name = 'time'
    time_var.units = 'days'
    time_var[:] = np.linspace(1338.5, 1998.5, len_time)
    
    scalar_axis_var = nc_file.createVariable('scalar_axis', 'f4', ('scalar_axis',))
    scalar_axis_var.long_name = 'none'
    scalar_axis_var[:] = np.array([0.0])
    
    # Create other variables specific to atmos_average_pstd.nc
    areo_var = nc_file.createVariable('areo', 'f4', ('time', 'scalar_axis'))
    areo_var.long_name = 'areo'
    areo_var.units = 'degrees'
    areo_vals = np.linspace(723.7, 1076.9, len_time)
    areo_data = np.zeros((len_time, 1))  # Create a 2D array with shape (len_time, 1)
    for i in range(len_time):
        areo_data[i, 0] = areo_vals[i]
    areo_var[:] = areo_data
    
    cldcol_var = nc_file.createVariable('cldcol', 'f4', ('time', 'lat', 'lon'))
    cldcol_var.long_name = 'ice column'
    cldcol_var[:] = np.random.uniform(1.2e-11, 4.1e-02, size=(len_time, 48, 96))
    
    dst_mass_micro_var = nc_file.createVariable('dst_mass_micro', 'f4', ('time', 'pstd', 'lat', 'lon'))
    dst_mass_micro_var.long_name = 'dust_mass'
    dst_mass_micro_var[:] = np.random.uniform(2.5e-16, 2.0e-04, size=(len_time, 44, 48, 96))
    
    omega_var = nc_file.createVariable('omega', 'f4', ('time', 'pstd', 'lat', 'lon'))
    omega_var.long_name = 'vertical wind'
    omega_var.units = 'Pa/s'
    omega_var[:] = np.random.uniform(-0.045597, 0.0806756, size=(133, 44, 48, 96))
    
    w_var = nc_file.createVariable('w', 'f4', ('time', 'pstd', 'lat', 'lon'))
    w_var.long_name = 'w'
    w_var.units = 'm/s'
    w_var[:] = np.random.uniform(-2.02603, 1.58804, size=(133, 44, 48, 96))
     
    ps_var = nc_file.createVariable('ps', 'f4', ('time', 'lat', 'lon'))
    ps_var.long_name = 'surface pressure'
    ps_var.units = 'Pa'
    ps_var[:] = np.random.uniform(176.8, 1318.8, size=(len_time, 48, 96))
    
    r_var = nc_file.createVariable('r', 'f4', ('time', 'pstd', 'lat', 'lon'))
    r_var.long_name = 'specific humidity'
    r_var.units = 'kg/kg'
    r_var[:] = np.random.uniform(9.2e-13, 3.4e-03, size=(len_time, 44, 48, 96))
    
    taudust_IR_var = nc_file.createVariable('taudust_IR', 'f4', ('time', 'lat', 'lon'))
    taudust_IR_var.long_name = 'Dust opacity IR'
    taudust_IR_var.units = 'op'
    taudust_IR_var[:] = np.random.uniform(0.0, 0.5, size=(len_time, 48, 96))
    
    temp_var = nc_file.createVariable('temp', 'f4', ('time', 'pstd', 'lat', 'lon'))
    temp_var.long_name = 'temperature'
    temp_var.units = 'K'
    temp_var[:] = np.random.uniform(104.8, 258.5, size=(len_time, 44, 48, 96))
    
    ts_var = nc_file.createVariable('ts', 'f4', ('time', 'lat', 'lon'))
    ts_var.long_name = 'Surface Temperature'
    ts_var.units = 'K'
    ts_var[:] = np.random.uniform(143.4, 258.7, size=(len_time, 48, 96))
    
    ucomp_var = nc_file.createVariable('ucomp', 'f4', ('time', 'pstd', 'lat', 'lon'))
    ucomp_var.long_name = 'zonal wind'
    ucomp_var.units = 'm/sec'
    ucomp_var[:] = np.random.uniform(-258.5, 209.6, size=(len_time, 44, 48, 96))
    
    vcomp_var = nc_file.createVariable('vcomp', 'f4', ('time', 'pstd', 'lat', 'lon'))
    vcomp_var.long_name = 'meridional wind'
    vcomp_var.units = 'm/sec'
    vcomp_var[:] = np.random.uniform(-94.7, 108.6, size=(len_time, 44, 48, 96))
    
    nc_file.close()
    print("Created 01336.atmos_average_pstd.nc")

def create_mgcm_atmos_diurn_pstd(short=False):
    """Create atmos_diurn_pstd.nc with the exact variables and structure as specified."""
    nc_file = Dataset('01336.atmos_diurn_pstd.nc', 'w', format='NETCDF4')

    # Shorten file length if wanted
    if short:
        len_time = 2
    else:
        len_time = 133
    
    # Define dimensions
    time_dim = nc_file.createDimension('time', len_time)
    time_of_day_24_dim = nc_file.createDimension('time_of_day_24', 24)
    lat_dim = nc_file.createDimension('lat', 48)
    lon_dim = nc_file.createDimension('lon', 96)
    scalar_axis_dim = nc_file.createDimension('scalar_axis', 1)
    
    # Create and populate key variables
    time_var = nc_file.createVariable('time', 'f4', ('time',))
    time_var.long_name = 'time'
    time_var.units = 'days'
    time_var[:] = np.linspace(1338.5, 1998.5, len_time)
    
    time_of_day_24_var = nc_file.createVariable('time_of_day_24', 'f4', ('time_of_day_24',))
    time_of_day_24_var.long_name = 'time of day'
    time_of_day_24_var.units = 'hours since 0000-00-00 00:00:00'
    time_of_day_24_values = np.array([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5,
                                     14.5, 15.5, 16.5, 17.5, 18.5, 19.5, 20.5, 21.5, 22.5, 23.5])
    time_of_day_24_var[:] = time_of_day_24_values
    
    lat_var = nc_file.createVariable('lat', 'f4', ('lat',))
    lat_var.long_name = 'latitude'
    lat_var.units = 'degrees_N'
    lat_var[:] = np.array([-88.125, -84.375, -80.625, -76.875, -73.125, -69.375, -65.625, -61.875, -58.125,
                          -54.375, -50.625, -46.875, -43.125, -39.375, -35.625, -31.875, -28.125, -24.375,
                          -20.625, -16.875, -13.125, -9.375, -5.625, -1.875, 1.875, 5.625, 9.375,
                          13.125, 16.875, 20.625, 24.375, 28.125, 31.875, 35.625, 39.375, 43.125,
                          46.875, 50.625, 54.375, 58.125, 61.875, 65.625, 69.375, 73.125, 76.875,
                          80.625, 84.375, 88.125])
    
    lon_var = nc_file.createVariable('lon', 'f4', ('lon',))
    lon_var.long_name = 'longitude'
    lon_var.units = 'degrees_E'
    lon_var[:] = np.array([1.875, 5.625, 9.375, 13.125, 16.875, 20.625, 24.375, 28.125, 31.875,
                          35.625, 39.375, 43.125, 46.875, 50.625, 54.375, 58.125, 61.875, 65.625,
                          69.375, 73.125, 76.875, 80.625, 84.375, 88.125, 91.875, 95.625, 99.375,
                          103.125, 106.875, 110.625, 114.375, 118.125, 121.875, 125.625, 129.375, 133.125,
                          136.875, 140.625, 144.375, 148.125, 151.875, 155.625, 159.375, 163.125, 166.875,
                          170.625, 174.375, 178.125, 181.875, 185.625, 189.375, 193.125, 196.875, 200.625,
                          204.375, 208.125, 211.875, 215.625, 219.375, 223.125, 226.875, 230.625, 234.375,
                          238.125, 241.875, 245.625, 249.375, 253.125, 256.875, 260.625, 264.375, 268.125,
                          271.875, 275.625, 279.375, 283.125, 286.875, 290.625, 294.375, 298.125, 301.875,
                          305.625, 309.375, 313.125, 316.875, 320.625, 324.375, 328.125, 331.875, 335.625,
                          339.375, 343.125, 346.875, 350.625, 354.375, 358.125])
    
    scalar_axis_var = nc_file.createVariable('scalar_axis', 'f4', ('scalar_axis',))
    scalar_axis_var.long_name = 'none'
    scalar_axis_var[:] = np.array([0.0])
    
    # Create other variables specific to atmos_diurn_pstd.nc
    areo_var = nc_file.createVariable('areo', 'f4', ('time', 'time_of_day_24', 'scalar_axis'))
    areo_var.long_name = 'areo'
    areo_var.units = 'degrees'
    
    # Create base values for areo dimension
    areo_base = np.linspace(721.2, 1077.3, len_time)

    # Create 3D array with shape (len_time, 24, 1)
    areo_data = np.zeros((len_time, 24, 1))

    # Fill array with increasing values
    for t in range(len_time):
        # Base value for this areo
        base_val = areo_base[t]
        
        # Daily oscillation (values increase slightly throughout the day)
        # Starting with a small offset and incrementing by a small amount
        for tod in range(24):
            # Small daily oscillation of ~0.4 degrees
            daily_increment = (tod / 24.0) * 0.4
            areo_data[t, tod, 0] = base_val - 0.2 + daily_increment

    # Assign the data to the variable
    areo_var[:] = areo_data
    
    ps_var = nc_file.createVariable('ps', 'f4', ('time', 'time_of_day_24', 'lat', 'lon'))
    ps_var.long_name = 'surface pressure'
    ps_var.units = 'Pa'
    ps_var[:] = np.random.uniform(167.9, 1338.7, size=(len_time, 24, 48, 96))
    
    nc_file.close()
    print("Created 01336.atmos_diurn_pstd.nc")

def create_mgcm_atmos_diurn(short=False):
    """Create atmos_diurn.nc with the exact variables and structure as specified."""
    nc_file = Dataset('01336.atmos_diurn.nc', 'w', format='NETCDF4')

    # Shorten file length if wanted
    if short:
        len_time = 2
    else:
        len_time = 133
    
    # Define dimensions
    time_dim = nc_file.createDimension('time', len_time)
    time_of_day_24_dim = nc_file.createDimension('time_of_day_24', 24)
    pfull_dim = nc_file.createDimension('pfull', 30)
    lat_dim = nc_file.createDimension('lat', 48)
    lon_dim = nc_file.createDimension('lon', 96)
    scalar_axis_dim = nc_file.createDimension('scalar_axis', 1)
    
    # Create key variables
    time_var = nc_file.createVariable('time', 'f4', ('time',))
    time_var.long_name = 'time'
    time_var.units = 'days'
    time_var[:] = np.linspace(1338.5, 1998.5, len_time)
    
    time_of_day_24_var = nc_file.createVariable('time_of_day_24', 'f4', ('time_of_day_24',))
    time_of_day_24_var.long_name = 'time of day'
    time_of_day_24_var.units = 'hours since 0000-00-00 00:00:00'
    time_of_day_24_values = np.array([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5,
                                     14.5, 15.5, 16.5, 17.5, 18.5, 19.5, 20.5, 21.5, 22.5, 23.5])
    time_of_day_24_var[:] = time_of_day_24_values
    
    lat_var = nc_file.createVariable('lat', 'f4', ('lat',))
    lat_var.long_name = 'latitude'
    lat_var.units = 'degrees_N'
    lat_var[:] = np.array([-88.125, -84.375, -80.625, -76.875, -73.125, -69.375, -65.625, -61.875, -58.125,
                          -54.375, -50.625, -46.875, -43.125, -39.375, -35.625, -31.875, -28.125, -24.375,
                          -20.625, -16.875, -13.125, -9.375, -5.625, -1.875, 1.875, 5.625, 9.375,
                          13.125, 16.875, 20.625, 24.375, 28.125, 31.875, 35.625, 39.375, 43.125,
                          46.875, 50.625, 54.375, 58.125, 61.875, 65.625, 69.375, 73.125, 76.875,
                          80.625, 84.375, 88.125])
    
    lon_var = nc_file.createVariable('lon', 'f4', ('lon',))
    lon_var.long_name = 'longitude'
    lon_var.units = 'degrees_E'
    lon_var[:] = np.array([1.875, 5.625, 9.375, 13.125, 16.875, 20.625, 24.375, 28.125, 31.875,
                          35.625, 39.375, 43.125, 46.875, 50.625, 54.375, 58.125, 61.875, 65.625,
                          69.375, 73.125, 76.875, 80.625, 84.375, 88.125, 91.875, 95.625, 99.375,
                          103.125, 106.875, 110.625, 114.375, 118.125, 121.875, 125.625, 129.375, 133.125,
                          136.875, 140.625, 144.375, 148.125, 151.875, 155.625, 159.375, 163.125, 166.875,
                          170.625, 174.375, 178.125, 181.875, 185.625, 189.375, 193.125, 196.875, 200.625,
                          204.375, 208.125, 211.875, 215.625, 219.375, 223.125, 226.875, 230.625, 234.375,
                          238.125, 241.875, 245.625, 249.375, 253.125, 256.875, 260.625, 264.375, 268.125,
                          271.875, 275.625, 279.375, 283.125, 286.875, 290.625, 294.375, 298.125, 301.875,
                          305.625, 309.375, 313.125, 316.875, 320.625, 324.375, 328.125, 331.875, 335.625,
                          339.375, 343.125, 346.875, 350.625, 354.375, 358.125])
    
    scalar_axis_var = nc_file.createVariable('scalar_axis', 'f4', ('scalar_axis',))
    scalar_axis_var.long_name = 'none'
    scalar_axis_var[:] = np.array([0.0])
    
    pfull_var = nc_file.createVariable('pfull', 'f4', ('pfull',))
    pfull_var.long_name = 'ref full pressure level'
    pfull_var.units = 'mb'
    pfull_values = np.array([3.44881953e-04, 1.09678471e-03, 3.48347419e-03, 9.73852715e-03,
                            2.46886197e-02, 5.58059295e-02, 1.07865789e-01, 1.78102296e-01,
                            2.73462607e-01, 4.16048919e-01, 6.21882691e-01, 9.06446215e-01,
                            1.28069137e+00, 1.74661073e+00, 2.29407255e+00, 2.90057273e+00,
                            3.53450185e+00, 4.16097960e+00, 4.74822077e+00, 5.27238854e+00,
                            5.71981194e+00, 6.08661884e+00, 6.37663706e+00, 6.59862648e+00,
                            6.76367744e+00, 6.88322325e+00, 6.96777592e+00, 7.01218097e+00,
                            7.03237068e+00, 7.04647442e+00])
    pfull_var[:] = pfull_values
    
    # Create specific variables for atmos_diurn.nc
    areo_var = nc_file.createVariable('areo', 'f4', ('time', 'time_of_day_24', 'scalar_axis'))
    areo_var.long_name = 'areo'
    areo_var.units = 'degrees'
    # Create base values for areo dimension
    areo_base = np.linspace(721.2, 1077.3, len_time)

    # Create 3D array with shape (len_time, 24, 1)
    areo_data = np.zeros((len_time, 24, 1))

    # Fill array with increasing values
    for t in range(len_time):
        # Base value for this areo
        base_val = areo_base[t]
        
        # Daily oscillation (values increase slightly throughout the day)
        # Starting with a small offset and incrementing by a small amount
        for tod in range(24):
            # Small daily oscillation of ~0.4 degrees
            daily_increment = (tod / 24.0) * 0.4
            areo_data[t, tod, 0] = base_val - 0.2 + daily_increment
    
    # Assign the data to the variable
    areo_var[:] = areo_data
    
    ps_var = nc_file.createVariable('ps', 'f4', ('time', 'time_of_day_24', 'lat', 'lon'))
    ps_var.long_name = 'surface pressure'
    ps_var.units = 'Pa'
    ps_var[:] = np.random.uniform(167.9, 1338.7, size=(len_time, 24, 48, 96))
    
    temp_var = nc_file.createVariable('temp', 'f4', ('time', 'time_of_day_24', 'pfull', 'lat', 'lon'))
    temp_var.long_name = 'temperature'
    temp_var.units = 'K'
    temp_var[:] = np.random.uniform(101.6, 286.5, size=(len_time, 24, 30, 48, 96))
    
    nc_file.close()
    print("Created 01336.atmos_diurn.nc")

def create_mgcm_atmos_average_pstd_c48(short=False):
    """Create atmos_average_pstd_c48.nc with the exact variables and structure as specified."""
    nc_file = Dataset('01336.atmos_average_pstd_c48.nc', 'w', format='NETCDF4')

    # Shorten file length if wanted
    if short:
        len_time = 2
    else:
        len_time = 133
    
    # Define dimensions - note this file has different lat/lon dimensions
    time_dim = nc_file.createDimension('time', len_time)
    pstd_dim = nc_file.createDimension('pstd', 48)
    lat_dim = nc_file.createDimension('lat', 90)
    lon_dim = nc_file.createDimension('lon', 180)
    scalar_axis_dim = nc_file.createDimension('scalar_axis', 1)
    
    # Create key variables
    pstd_var = nc_file.createVariable('pstd', 'f4', ('pstd',))
    pstd_var.long_name = 'pressure'
    pstd_var.units = 'Pa'
    # Using exactly the pstd values from the file but extending to 48 elements
    pstd_values = np.array([1.0e-05, 3.0e-05, 5.0e-05, 1.0e-04, 3.0e-04, 5.0e-04, 3.0e-03, 5.0e-03, 1.0e-02,
                           3.0e-02, 5.0e-02, 1.0e-01, 2.0e-01, 3.0e-01, 5.0e-01, 1.0e+00, 2.0e+00, 3.0e+00,
                           5.0e+00, 7.0e+00, 1.0e+01, 2.0e+01, 3.0e+01, 5.0e+01, 7.0e+01, 1.0e+02, 1.5e+02,
                           2.0e+02, 2.5e+02, 3.0e+02, 3.5e+02, 4.0e+02, 4.5e+02, 5.0e+02, 5.3e+02, 5.5e+02,
                           5.9e+02, 6.0e+02, 6.3e+02, 6.5e+02, 6.9e+02, 7.0e+02, 7.5e+02, 8.0e+02, 8.5e+02,
                           9.0e+02, 9.5e+02, 1.0e+03])
    pstd_var[:] = pstd_values
    
    # NOTE: This file uses different lat and lon values than the other files
    lat_var = nc_file.createVariable('lat', 'f4', ('lat',))
    lat_var.long_name = 'latitude'
    lat_var.units = 'degrees_N'
    lat_values = np.array([-89.0, -87.0, -85.0, -83.0, -81.0, -79.0, -77.0, -75.0, -73.0, -71.0, -69.0, -67.0, -65.0, -63.0,
                          -61.0, -59.0, -57.0, -55.0, -53.0, -51.0, -49.0, -47.0, -45.0, -43.0, -41.0, -39.0, -37.0, -35.0,
                          -33.0, -31.0, -29.0, -27.0, -25.0, -23.0, -21.0, -19.0, -17.0, -15.0, -13.0, -11.0, -9.0, -7.0,
                          -5.0, -3.0, -1.0, 1.0, 3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0, 17.0, 19.0, 21.0,
                          23.0, 25.0, 27.0, 29.0, 31.0, 33.0, 35.0, 37.0, 39.0, 41.0, 43.0, 45.0, 47.0, 49.0,
                          51.0, 53.0, 55.0, 57.0, 59.0, 61.0, 63.0, 65.0, 67.0, 69.0, 71.0, 73.0, 75.0, 77.0,
                          79.0, 81.0, 83.0, 85.0, 87.0, 89.0])
    lat_var[:] = lat_values
    
    lon_var = nc_file.createVariable('lon', 'f4', ('lon',))
    lon_var.long_name = 'longitude'
    lon_var.units = 'degrees_E'
    lon_values = np.array([1.0, 3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0, 17.0, 19.0, 21.0, 23.0, 25.0, 27.0,
                          29.0, 31.0, 33.0, 35.0, 37.0, 39.0, 41.0, 43.0, 45.0, 47.0, 49.0, 51.0, 53.0, 55.0,
                          57.0, 59.0, 61.0, 63.0, 65.0, 67.0, 69.0, 71.0, 73.0, 75.0, 77.0, 79.0, 81.0, 83.0,
                          85.0, 87.0, 89.0, 91.0, 93.0, 95.0, 97.0, 99.0, 101.0, 103.0, 105.0, 107.0, 109.0, 111.0,
                          113.0, 115.0, 117.0, 119.0, 121.0, 123.0, 125.0, 127.0, 129.0, 131.0, 133.0, 135.0, 137.0, 139.0,
                          141.0, 143.0, 145.0, 147.0, 149.0, 151.0, 153.0, 155.0, 157.0, 159.0, 161.0, 163.0, 165.0, 167.0,
                          169.0, 171.0, 173.0, 175.0, 177.0, 179.0, 181.0, 183.0, 185.0, 187.0, 189.0, 191.0, 193.0, 195.0,
                          197.0, 199.0, 201.0, 203.0, 205.0, 207.0, 209.0, 211.0, 213.0, 215.0, 217.0, 219.0, 221.0, 223.0,
                          225.0, 227.0, 229.0, 231.0, 233.0, 235.0, 237.0, 239.0, 241.0, 243.0, 245.0, 247.0, 249.0, 251.0,
                          253.0, 255.0, 257.0, 259.0, 261.0, 263.0, 265.0, 267.0, 269.0, 271.0, 273.0, 275.0, 277.0, 279.0,
                          281.0, 283.0, 285.0, 287.0, 289.0, 291.0, 293.0, 295.0, 297.0, 299.0, 301.0, 303.0, 305.0, 307.0,
                          309.0, 311.0, 313.0, 315.0, 317.0, 319.0, 321.0, 323.0, 325.0, 327.0, 329.0, 331.0, 333.0, 335.0,
                          337.0, 339.0, 341.0, 343.0, 345.0, 347.0, 349.0, 351.0, 353.0, 355.0, 357.0, 359.0])
    lon_var[:] = lon_values
    
    time_var = nc_file.createVariable('time', 'f4', ('time',))
    time_var.long_name = 'time'
    time_var.units = 'days'
    time_var[:] = np.linspace(670.5, 1330.5, len_time)
    
    scalar_axis_var = nc_file.createVariable('scalar_axis', 'f4', ('scalar_axis',))
    scalar_axis_var.long_name = 'none'
    scalar_axis_var[:] = np.array([0.0])
    
    # Create specific variables for atmos_average_pstd_c48.nc
    areo_var = nc_file.createVariable('areo', 'f4', ('time', 'scalar_axis'))
    areo_var.long_name = 'areo'
    areo_var.units = 'degrees'
    areo_vals = np.linspace(362.1, 716.6, len_time)
    areo_data = np.zeros((len_time, 1))  # Create a 2D array with shape (len_time, 1)
    for i in range(len_time):
        areo_data[i, 0] = areo_vals[i]
    areo_var[:] = areo_data
    
    temp_var = nc_file.createVariable('temp', 'f4', ('time', 'pstd', 'lat', 'lon'))
    temp_var.long_name = 'temperature'
    temp_var.units = 'K'
    temp_var[:] = np.random.uniform(106.9, 260.6, size=(len_time, 48, 90, 180))
    
    nc_file.close()
    print("Created 01336.atmos_average_pstd_c48.nc")

def main(short=False):
    """Main function to create all MGCM test files."""
    if short:
        print("Making short GCM files")
    create_mgcm_fixed()
    create_mgcm_atmos_average(short)
    create_mgcm_atmos_daily(short)
    create_mgcm_atmos_average_pstd(short)
    create_mgcm_atmos_diurn_pstd(short)
    create_mgcm_atmos_diurn(short)
    create_mgcm_atmos_average_pstd_c48(short)
    
    print("All MGCM test NetCDF files created successfully.")

if __name__ == "__main__":
    short_flag = False
    if len(sys.argv) > 1:
        for arg in sys.argv:
            if arg.lower() == "short":
                short_flag = True
    main(short=short_flag)