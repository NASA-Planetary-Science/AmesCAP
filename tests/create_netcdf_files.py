import numpy as np
from netCDF4 import Dataset
import os

# This script creates NetCDF files for Mars climate model data
# based on the specifications provided
import argparse     # Parse arguments

parser = argparse.ArgumentParser(
    prog=('create_netcdf_files'),
    description=(
        f"Create NetCDF files for Mars climate model data.\n"
        f"Files are created in the test directory for the executable.\n"
    ),
    formatter_class=argparse.RawTextHelpFormatter
)

parser.add_argument('directory', default=os.getcwd(),
    type=str,
    help=(
        f"Path to save the NetCDF files.\n"
    )
)

args = parser.parse_args()

# Create directory for the output files if it doesn't exist
output_dir = args.directory
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Common latitude and longitude arrays
lat_array = np.array([-88.125, -84.375, -80.625, -76.875, -73.125, -69.375, -65.625, -61.875, -58.125, 
                     -54.375, -50.625, -46.875, -43.125, -39.375, -35.625, -31.875, -28.125, -24.375,
                     -20.625, -16.875, -13.125, -9.375, -5.625, -1.875, 1.875, 5.625, 9.375, 13.125, 
                     16.875, 20.625, 24.375, 28.125, 31.875, 35.625, 39.375, 43.125, 46.875, 50.625, 
                     54.375, 58.125, 61.875, 65.625, 69.375, 73.125, 76.875, 80.625, 84.375, 88.125])

lon_array = np.array([1.875, 5.625, 9.375, 13.125, 16.875, 20.625, 24.375, 28.125, 31.875, 35.625, 
                      39.375, 43.125, 46.875, 50.625, 54.375, 58.125, 61.875, 65.625, 69.375, 73.125, 
                      76.875, 80.625, 84.375, 88.125, 91.875, 95.625, 99.375, 103.125, 106.875, 110.625, 
                      114.375, 118.125, 121.875, 125.625, 129.375, 133.125, 136.875, 140.625, 144.375, 
                      148.125, 151.875, 155.625, 159.375, 163.125, 166.875, 170.625, 174.375, 178.125, 
                      181.875, 185.625, 189.375, 193.125, 196.875, 200.625, 204.375, 208.125, 211.875, 
                      215.625, 219.375, 223.125, 226.875, 230.625, 234.375, 238.125, 241.875, 245.625, 
                      249.375, 253.125, 256.875, 260.625, 264.375, 268.125, 271.875, 275.625, 279.375, 
                      283.125, 286.875, 290.625, 294.375, 298.125, 301.875, 305.625, 309.375, 313.125, 
                      316.875, 320.625, 324.375, 328.125, 331.875, 335.625, 339.375, 343.125, 346.875, 
                      350.625, 354.375, 358.125])

# c48 versions
lat_c48 = np.array([-89., -87., -85., -83., -81., -79., -77., -75., -73., -71., -69., -67., -65., -63., 
                  -61., -59., -57., -55., -53., -51., -49., -47., -45., -43., -41., -39., -37., -35., 
                  -33., -31., -29., -27., -25., -23., -21., -19., -17., -15., -13., -11., -9., -7., 
                  -5., -3., -1., 1., 3., 5., 7., 9., 11., 13., 15., 17., 19., 21., 23., 25., 27., 
                  29., 31., 33., 35., 37., 39., 41., 43., 45., 47., 49., 51., 53., 55., 57., 59., 
                  61., 63., 65., 67., 69., 71., 73., 75., 77., 79., 81., 83., 85., 87., 89.])

lon_c48 = np.array([1., 3., 5., 7., 9., 11., 13., 15., 17., 19., 21., 23., 25., 27., 29., 31., 33., 
                  35., 37., 39., 41., 43., 45., 47., 49., 51., 53., 55., 57., 59., 61., 63., 65., 
                  67., 69., 71., 73., 75., 77., 79., 81., 83., 85., 87., 89., 91., 93., 95., 97., 
                  99., 101., 103., 105., 107., 109., 111., 113., 115., 117., 119., 121., 123., 125., 
                  127., 129., 131., 133., 135., 137., 139., 141., 143., 145., 147., 149., 151., 153., 
                  155., 157., 159., 161., 163., 165., 167., 169., 171., 173., 175., 177., 179., 181., 
                  183., 185., 187., 189., 191., 193., 195., 197., 199., 201., 203., 205., 207., 209., 
                  211., 213., 215., 217., 219., 221., 223., 225., 227., 229., 231., 233., 235., 237., 
                  239., 241., 243., 245., 247., 249., 251., 253., 255., 257., 259., 261., 263., 265., 
                  267., 269., 271., 273., 275., 277., 279., 281., 283., 285., 287., 289., 291., 293., 
                  295., 297., 299., 301., 303., 305., 307., 309., 311., 313., 315., 317., 319., 321., 
                  323., 325., 327., 329., 331., 333., 335., 337., 339., 341., 343., 345., 347., 349., 
                  351., 353., 355., 357., 359.])

# Pressure standard level array
pstd = np.array([1.0e+03, 9.5e+02, 9.0e+02, 8.5e+02, 8.0e+02, 7.5e+02, 7.0e+02, 6.5e+02, 6.0e+02, 
                5.5e+02, 5.0e+02, 4.5e+02, 4.0e+02, 3.5e+02, 3.0e+02, 2.5e+02, 2.0e+02, 1.5e+02, 
                1.0e+02, 7.0e+01, 5.0e+01, 3.0e+01, 2.0e+01, 1.0e+01, 7.0e+00, 5.0e+00, 3.0e+00, 
                2.0e+00, 1.0e+00, 5.0e-01, 3.0e-01, 2.0e-01, 1.0e-01, 5.0e-02, 3.0e-02, 1.0e-02, 
                5.0e-03, 3.0e-03, 5.0e-04, 3.0e-04, 1.0e-04, 5.0e-05, 3.0e-05, 1.0e-05])

# Other common arrays
pfull = np.linspace(0.0, 7.0, 13)  # 13 pressure levels
phalf = np.linspace(0.0, 7.0, 14)  # 14 half pressure levels
bk = np.linspace(0.0, 1.0, 14)     # vertical coordinate sigma values
pk = np.linspace(0.0, 2.8, 14)     # pressure part of hybrid coordinate
scalar_axis = np.array([0.0])      # scalar axis
bnds = np.array([0, 1])            # bounds dimension

# Time arrays
time_avg = np.linspace(1338.5, 1998.5, 22)     # For average files
time_daily = np.linspace(1336.2, 2004.0, 23)   # For daily files
time_pstd = np.linspace(670.5, 1330.5, 22)     # For pstd files
time_of_day_24 = np.linspace(0.5, 23.5, 24)    # Time of day

# ----------------------------------------------------------------------
# 1. Create atmos_average.nc
# ----------------------------------------------------------------------
def create_atmos_average():
    filename = os.path.join(output_dir, "01336.atmos_average.nc")
    with Dataset(filename, 'w', format='NETCDF4') as nc:
        # Create dimensions
        time_dim = nc.createDimension('time', len(time_avg))
        lat_dim = nc.createDimension('lat', len(lat_array))
        lon_dim = nc.createDimension('lon', len(lon_array))
        pfull_dim = nc.createDimension('pfull', len(pfull))
        phalf_dim = nc.createDimension('phalf', len(phalf))
        scalar_dim = nc.createDimension('scalar_axis', 1)
        
        # Create variables
        time_var = nc.createVariable('time', 'f4', ('time',))
        lat_var = nc.createVariable('lat', 'f4', ('lat',))
        lon_var = nc.createVariable('lon', 'f4', ('lon',))
        pfull_var = nc.createVariable('pfull', 'f4', ('pfull',))
        phalf_var = nc.createVariable('phalf', 'f4', ('phalf',))
        bk_var = nc.createVariable('bk', 'f4', ('phalf',))
        pk_var = nc.createVariable('pk', 'f4', ('phalf',))
        scalar_var = nc.createVariable('scalar_axis', 'f4', ('scalar_axis',))
        
        # Create data variables
        areo_var = nc.createVariable('areo', 'f4', ('time', 'scalar_axis'))
        ps_var = nc.createVariable('ps', 'f4', ('time', 'lat', 'lon'))
        temp_var = nc.createVariable('temp', 'f4', ('time', 'pfull', 'lat', 'lon'))
        ts_var = nc.createVariable('ts', 'f4', ('time', 'lat', 'lon'))
        ucomp_var = nc.createVariable('ucomp', 'f4', ('time', 'pfull', 'lat', 'lon'))
        vcomp_var = nc.createVariable('vcomp', 'f4', ('time', 'pfull', 'lat', 'lon'))
        r_var = nc.createVariable('r', 'f4', ('time', 'pfull', 'lat', 'lon'))
        cldcol_var = nc.createVariable('cldcol', 'f4', ('time', 'lat', 'lon'))
        dst_mass_micro_var = nc.createVariable('dst_mass_micro', 'f4', ('time', 'pfull', 'lat', 'lon'))
        taudust_IR_var = nc.createVariable('taudust_IR', 'f4', ('time', 'lat', 'lon'))
        
        # Fill coordinate variables
        time_var[:] = time_avg
        lat_var[:] = lat_array
        lon_var[:] = lon_array
        pfull_var[:] = pfull
        phalf_var[:] = phalf
        bk_var[:] = bk
        pk_var[:] = pk
        scalar_var[:] = scalar_axis
        
        # Set attributes
        time_var.units = 'days'
        time_var.long_name = 'time'
        lat_var.units = 'degrees_N'
        lat_var.long_name = 'latitude'
        lon_var.units = 'degrees_E'
        lon_var.long_name = 'longitude'
        pfull_var.units = 'mb'
        pfull_var.long_name = 'ref full pressure level'
        phalf_var.units = 'mb'
        phalf_var.long_name = 'ref half pressure level'
        bk_var.long_name = 'vertical coordinate sigma value'
        pk_var.units = 'Pa'
        pk_var.long_name = 'pressure part of the hybrid coordinate'
        scalar_var.long_name = 'none'
        
        # Fill data variables with random values within specified ranges
        areo_var[:, :] = np.random.uniform(721.3, 1077.2, size=(len(time_avg), 1))
        areo_var.units = 'degrees'
        areo_var.long_name = 'areo'
        
        ps_var[:] = np.random.uniform(176.8, 1318.8, size=(len(time_avg), len(lat_array), len(lon_array)))
        ps_var.units = 'Pa'
        ps_var.long_name = 'surface pressure'
        
        temp_var[:] = np.random.uniform(104.1, 258.8, size=(len(time_avg), len(pfull), len(lat_array), len(lon_array)))
        temp_var.units = 'K'
        temp_var.long_name = 'temperature'
        
        ts_var[:] = np.random.uniform(143.4, 258.7, size=(len(time_avg), len(lat_array), len(lon_array)))
        ts_var.units = 'K'
        ts_var.long_name = 'Surface Temperature'
        
        ucomp_var[:] = np.random.uniform(-268.7, 212.7, size=(len(time_avg), len(pfull), len(lat_array), len(lon_array)))
        ucomp_var.units = 'm/sec'
        ucomp_var.long_name = 'zonal wind'
        
        vcomp_var[:] = np.random.uniform(-97.5, 109.6, size=(len(time_avg), len(pfull), len(lat_array), len(lon_array)))
        vcomp_var.units = 'm/sec'
        vcomp_var.long_name = 'meridional wind'
        
        r_var[:] = np.random.uniform(8.6e-13, 3.4e-03, size=(len(time_avg), len(pfull), len(lat_array), len(lon_array)))
        r_var.units = 'kg/kg'
        r_var.long_name = 'specific humidity'
        
        cldcol_var[:] = np.random.uniform(1.2e-11, 4.1e-02, size=(len(time_avg), len(lat_array), len(lon_array)))
        cldcol_var.long_name = 'ice column'
        
        dst_mass_micro_var[:] = np.random.uniform(1.5e-17, 2.5e-04, size=(len(time_avg), len(pfull), len(lat_array), len(lon_array)))
        dst_mass_micro_var.long_name = 'dust_mass'
        
        taudust_IR_var[:] = np.random.uniform(0.0, 0.5, size=(len(time_avg), len(lat_array), len(lon_array)))
        taudust_IR_var.units = 'op'
        taudust_IR_var.long_name = 'Dust opacity IR'

# ----------------------------------------------------------------------
# 2. Create atmos_daily.nc
# ----------------------------------------------------------------------
def create_atmos_daily():
    filename = os.path.join(output_dir, "01336.atmos_daily.nc")
    with Dataset(filename, 'w', format='NETCDF4') as nc:
        # Create dimensions
        time_dim = nc.createDimension('time', len(time_daily))
        lat_dim = nc.createDimension('lat', len(lat_array))
        lon_dim = nc.createDimension('lon', len(lon_array))
        pfull_dim = nc.createDimension('pfull', len(pfull))
        scalar_dim = nc.createDimension('scalar_axis', 1)
        
        # Create variables
        time_var = nc.createVariable('time', 'f4', ('time',))
        lat_var = nc.createVariable('lat', 'f4', ('lat',))
        lon_var = nc.createVariable('lon', 'f4', ('lon',))
        pfull_var = nc.createVariable('pfull', 'f4', ('pfull',))
        scalar_var = nc.createVariable('scalar_axis', 'f4', ('scalar_axis',))
        
        # Create data variables
        areo_var = nc.createVariable('areo', 'f4', ('time', 'scalar_axis'))
        ps_var = nc.createVariable('ps', 'f4', ('time', 'lat', 'lon'))
        temp_var = nc.createVariable('temp', 'f4', ('time', 'pfull', 'lat', 'lon'))
        
        # Fill coordinate variables
        time_var[:] = time_daily
        lat_var[:] = lat_array
        lon_var[:] = lon_array
        pfull_var[:] = pfull
        scalar_var[:] = scalar_axis
        
        # Set attributes
        time_var.units = 'days'
        time_var.long_name = 'time'
        lat_var.units = 'degrees_N'
        lat_var.long_name = 'latitude'
        lon_var.units = 'degrees_E'
        lon_var.long_name = 'longitude'
        pfull_var.units = 'mb'
        pfull_var.long_name = 'ref full pressure level'
        scalar_var.long_name = 'none'
        
        # Fill data variables with random values within specified ranges
        areo_var[:, :] = np.random.uniform(720.1, 1080.0, size=(len(time_daily), 1))
        areo_var.units = 'degrees'
        areo_var.long_name = 'areo'
        
        ps_var[:] = np.random.uniform(170.3, 1340.2, size=(len(time_daily), len(lat_array), len(lon_array)))
        ps_var.units = 'Pa'
        ps_var.long_name = 'surface pressure'
        
        temp_var[:] = np.random.uniform(101.6, 283.9, size=(len(time_daily), len(pfull), len(lat_array), len(lon_array)))
        temp_var.units = 'K'
        temp_var.long_name = 'temperature'

# ----------------------------------------------------------------------
# 3. Create atmos_diurn.nc
# ----------------------------------------------------------------------
def create_atmos_diurn():
    filename = os.path.join(output_dir, "01336.atmos_diurn.nc")
    with Dataset(filename, 'w', format='NETCDF4') as nc:
        # Create dimensions
        time_dim = nc.createDimension('time', len(time_avg))
        time_of_day_dim = nc.createDimension('time_of_day_24', len(time_of_day_24))
        lat_dim = nc.createDimension('lat', len(lat_array))
        lon_dim = nc.createDimension('lon', len(lon_array))
        pfull_dim = nc.createDimension('pfull', len(pfull))
        scalar_dim = nc.createDimension('scalar_axis', 1)
        
        # Create variables
        time_var = nc.createVariable('time', 'f4', ('time',))
        time_of_day_var = nc.createVariable('time_of_day_24', 'f4', ('time_of_day_24',))
        lat_var = nc.createVariable('lat', 'f4', ('lat',))
        lon_var = nc.createVariable('lon', 'f4', ('lon',))
        pfull_var = nc.createVariable('pfull', 'f4', ('pfull',))
        scalar_var = nc.createVariable('scalar_axis', 'f4', ('scalar_axis',))
        
        # Create data variables
        areo_var = nc.createVariable('areo', 'f4', ('time', 'time_of_day_24', 'scalar_axis'))
        ps_var = nc.createVariable('ps', 'f4', ('time', 'time_of_day_24', 'lat', 'lon'))
        temp_var = nc.createVariable('temp', 'f4', ('time', 'time_of_day_24', 'pfull', 'lat', 'lon'))
        
        # Fill coordinate variables
        time_var[:] = time_avg
        time_of_day_var[:] = time_of_day_24
        lat_var[:] = lat_array
        lon_var[:] = lon_array
        pfull_var[:] = pfull
        scalar_var[:] = scalar_axis
        
        # Set attributes
        time_var.units = 'days'
        time_var.long_name = 'time'
        time_of_day_var.units = 'hours since 0000-00-00 00:00:00'
        time_of_day_var.long_name = 'time of day'
        lat_var.units = 'degrees_N'
        lat_var.long_name = 'latitude'
        lon_var.units = 'degrees_E'
        lon_var.long_name = 'longitude'
        pfull_var.units = 'mb'
        pfull_var.long_name = 'ref full pressure level'
        scalar_var.long_name = 'none'
        
        # Fill data variables with random values within specified ranges
        areo_var[:, :, :] = np.random.uniform(721.1, 1077.4, size=(len(time_avg), len(time_of_day_24), 1))
        areo_var.units = 'degrees'
        areo_var.long_name = 'areo'
        
        ps_var[:] = np.random.uniform(167.9, 1338.7, size=(len(time_avg), len(time_of_day_24), len(lat_array), len(lon_array)))
        ps_var.units = 'Pa'
        ps_var.long_name = 'surface pressure'
        
        temp_var[:] = np.random.uniform(101.6, 286.5, size=(len(time_avg), len(time_of_day_24), len(pfull), len(lat_array), len(lon_array)))
        temp_var.units = 'K'
        temp_var.long_name = 'temperature'

# ----------------------------------------------------------------------
# 4. Create fixed.nc
# ----------------------------------------------------------------------
def create_fixed():
    filename = os.path.join(output_dir, "01336.fixed.nc")
    with Dataset(filename, 'w', format='NETCDF4') as nc:
        # Create dimensions
        lat_dim = nc.createDimension('lat', len(lat_array))
        lon_dim = nc.createDimension('lon', len(lon_array))
        phalf_dim = nc.createDimension('phalf', len(phalf))
        bnds_dim = nc.createDimension('bnds', len(bnds))
        
        # Create variables
        lat_var = nc.createVariable('lat', 'f4', ('lat',))
        lon_var = nc.createVariable('lon', 'f4', ('lon',))
        phalf_var = nc.createVariable('phalf', 'f4', ('phalf',))
        bk_var = nc.createVariable('bk', 'f4', ('phalf',))
        pk_var = nc.createVariable('pk', 'f4', ('phalf',))
        
        # Bounds variables
        grid_yt_bnds_var = nc.createVariable('grid_yt_bnds', 'f4', ('lat', 'bnds'))
        grid_xt_bnds_var = nc.createVariable('grid_xt_bnds', 'f4', ('lon', 'bnds'))
        
        # Create data variables
        zsurf_var = nc.createVariable('zsurf', 'f4', ('lat', 'lon'))
        thin_var = nc.createVariable('thin', 'f4', ('lat', 'lon'))
        alb_var = nc.createVariable('alb', 'f4', ('lat', 'lon'))
        emis_var = nc.createVariable('emis', 'f4', ('lat', 'lon'))
        gice_var = nc.createVariable('gice', 'f4', ('lat', 'lon'))
        
        # Fill coordinate variables
        lat_var[:] = lat_array
        lon_var[:] = lon_array
        phalf_var[:] = phalf
        bk_var[:] = bk
        pk_var[:] = pk
        
        # Create bounds arrays
        for i in range(len(lat_array)):
            grid_yt_bnds_var[i, 0] = max(-90.0, lat_array[i] - 1.875)
            grid_yt_bnds_var[i, 1] = min(90.0, lat_array[i] + 1.875)
            
        for i in range(len(lon_array)):
            grid_xt_bnds_var[i, 0] = max(0.0, lon_array[i] - 1.875)
            grid_xt_bnds_var[i, 1] = min(360.0, lon_array[i] + 1.875)
        
        # Set attributes
        lat_var.units = 'degrees_N'
        lat_var.long_name = 'latitude'
        lon_var.units = 'degrees_E'
        lon_var.long_name = 'longitude'
        phalf_var.units = 'mb'
        phalf_var.long_name = 'ref half pressure level'
        bk_var.long_name = 'vertical coordinate sigma value'
        pk_var.units = 'Pa'
        pk_var.long_name = 'pressure part of the hybrid coordinate'
        grid_yt_bnds_var.units = 'degrees_N'
        grid_yt_bnds_var.long_name = 'T-cell latitude'
        grid_xt_bnds_var.units = 'degrees_E'
        grid_xt_bnds_var.long_name = 'T-cell longitude'
        
        # Fill data variables with random values within specified ranges
        zsurf_var[:] = np.random.uniform(-7.1e+03, 1.1e+04, size=(len(lat_array), len(lon_array)))
        zsurf_var.units = 'm'
        zsurf_var.long_name = 'surface height'
        
        thin_var[:] = np.random.uniform(40.6, 1037.6, size=(len(lat_array), len(lon_array)))
        thin_var.units = 'mks'
        thin_var.long_name = 'Surface Thermal Inertia'
        
        alb_var[:] = np.random.uniform(0.1, 0.3, size=(len(lat_array), len(lon_array)))
        alb_var.units = 'mks'
        alb_var.long_name = 'Surface Albedo'
        
        emis_var[:] = np.random.uniform(0.9, 1.0, size=(len(lat_array), len(lon_array)))
        emis_var.long_name = 'Surface Emissivity'
        
        gice_var[:] = np.random.uniform(-57.9, 58.6, size=(len(lat_array), len(lon_array)))
        gice_var.long_name = 'GRS Ice'

# ----------------------------------------------------------------------
# 5. Create atmos_average_pstd_c48.nc
# ----------------------------------------------------------------------
def create_atmos_average_pstd_c48():
    filename = os.path.join(output_dir, "01336.atmos_average_pstd_c48.nc")
    with Dataset(filename, 'w', format='NETCDF4') as nc:
        # Create dimensions
        time_dim = nc.createDimension('time', len(time_pstd))
        lat_dim = nc.createDimension('lat', len(lat_c48))
        lon_dim = nc.createDimension('lon', len(lon_c48))
        pstd_dim = nc.createDimension('pstd', len(pstd))
        scalar_dim = nc.createDimension('scalar_axis', 1)
        
        # Create variables
        time_var = nc.createVariable('time', 'f4', ('time',))
        lat_var = nc.createVariable('lat', 'f4', ('lat',))
        lon_var = nc.createVariable('lon', 'f4', ('lon',))
        pstd_var = nc.createVariable('pstd', 'f4', ('pstd',))
        scalar_var = nc.createVariable('scalar_axis', 'f4', ('scalar_axis',))
        
        # Create data variables
        areo_var = nc.createVariable('areo', 'f4', ('time', 'scalar_axis'))
        temp_var = nc.createVariable('temp', 'f4', ('time', 'pstd', 'lat', 'lon'))
        
        # Fill coordinate variables
        time_var[:] = time_pstd
        lat_var[:] = lat_c48
        lon_var[:] = lon_c48
        pstd_var[:] = pstd
        scalar_var[:] = scalar_axis
        
        # Set attributes
        time_var.units = 'days'
        time_var.long_name = 'time'
        lat_var.units = 'degrees_N'
        lat_var.long_name = 'latitude'
        lon_var.units = 'degrees_E'
        lon_var.long_name = 'longitude'
        pstd_var.units = 'Pa'
        pstd_var.long_name = 'pressure'
        scalar_var.long_name = 'none'
        
        # Fill data variables with random values within specified ranges
        areo_var[:, :] = np.random.uniform(361.3, 717.2, size=(len(time_pstd), 1))
        areo_var.units = 'degrees'
        areo_var.long_name = 'areo'
        
        temp_var[:] = np.random.uniform(106.9, 260.6, size=(len(time_pstd), len(pstd), len(lat_c48), len(lon_c48)))
        temp_var.units = 'K'
        temp_var.long_name = 'temperature'

# ----------------------------------------------------------------------
# 6. Create atmos_average_pstd.nc
# ----------------------------------------------------------------------
def create_atmos_average_pstd():
    filename = os.path.join(output_dir, "01336.atmos_average_pstd.nc")
    with Dataset(filename, 'w', format='NETCDF4') as nc:
        # Create dimensions
        time_dim = nc.createDimension('time', len(time_pstd))
        lat_dim = nc.createDimension('lat', len(lat_array))
        lon_dim = nc.createDimension('lon', len(lon_array))
        pstd_dim = nc.createDimension('pstd', len(pstd))
        scalar_dim = nc.createDimension('scalar_axis', 1)
        
        # Create variables
        time_var = nc.createVariable('time', 'f4', ('time',))
        lat_var = nc.createVariable('lat', 'f4', ('lat',))
        lon_var = nc.createVariable('lon', 'f4', ('lon',))
        pstd_var = nc.createVariable('pstd', 'f4', ('pstd',))
        scalar_var = nc.createVariable('scalar_axis', 'f4', ('scalar_axis',))
        
        # Create data variables
        areo_var = nc.createVariable('areo', 'f4', ('time', 'scalar_axis'))
        temp_var = nc.createVariable('temp', 'f4', ('time', 'pstd', 'lat', 'lon'))
        ucomp_var = nc.createVariable('ucomp', 'f4', ('time', 'pstd', 'lat', 'lon'))
        vcomp_var = nc.createVariable('vcomp', 'f4', ('time', 'pstd', 'lat', 'lon'))
        r_var = nc.createVariable('r', 'f4', ('time', 'pstd', 'lat', 'lon'))
        dst_mass_micro_var = nc.createVariable('dst_mass_micro', 'f4', ('time', 'pstd', 'lat', 'lon'))
        ps_var = nc.createVariable('ps', 'f4', ('time', 'lat', 'lon'))
        ts_var = nc.createVariable('ts', 'f4', ('time', 'lat', 'lon'))
        cldcol_var = nc.createVariable('cldcol', 'f4', ('time', 'lat', 'lon'))
        taudust_IR_var = nc.createVariable('taudust_IR', 'f4', ('time', 'lat', 'lon'))
        
        # Fill coordinate variables
        time_var[:] = time_pstd
        lat_var[:] = lat_array
        lon_var[:] = lon_array
        pstd_var[:] = pstd
        scalar_var[:] = scalar_axis
        
        # Set attributes
        time_var.units = 'days'
        time_var.long_name = 'time'
        lat_var.units = 'degrees_N'
        lat_var.long_name = 'latitude'
        lon_var.units = 'degrees_E'
        lon_var.long_name = 'longitude'
        pstd_var.units = 'Pa'
        pstd_var.long_name = 'pressure'
        scalar_var.long_name = 'none'
        
        # Fill data variables with random values within specified ranges
        areo_var[:, :] = np.random.uniform(361.3, 717.2, size=(len(time_pstd), 1))
        areo_var.units = 'degrees'
        areo_var.long_name = 'areo'
        
        temp_var[:] = np.random.uniform(106.9, 260.6, size=(len(time_pstd), len(pstd), len(lat_array), len(lon_array)))
        temp_var.units = 'K'
        temp_var.long_name = 'temperature'
        
        dst_mass_micro_var[:] = np.random.uniform(2.5e-16, 2.0e-04, size=(len(time_pstd), len(pstd), len(lat_array), len(lon_array)))
        dst_mass_micro_var.units = ''
        dst_mass_micro_var.long_name = 'dust_mass'
        
        r_var[:] = np.random.uniform(9.2e-13, 3.4e-03, size=(len(time_pstd), len(pstd), len(lat_array), len(lon_array)))
        r_var.units = 'kg/kg'
        r_var.long_name = 'specific humidity'
        
        ucomp_var[:] = np.random.uniform(-258.5, 209.6, size=(len(time_pstd), len(pstd), len(lat_array), len(lon_array)))
        ucomp_var.units = 'm/sec'
        ucomp_var.long_name = 'zonal wind'
        
        vcomp_var[:] = np.random.uniform(-94.7, 108.6, size=(len(time_pstd), len(pstd), len(lat_array), len(lon_array)))
        vcomp_var.units = 'm/sec'
        vcomp_var.long_name = 'meridional wind'
        
        ps_var[:] = np.random.uniform(167.9, 1338.7, size=(len(time_avg), len(lat_array), len(lon_array)))
        ps_var.units = 'Pa'
        ps_var.long_name = 'surface pressure'
        
        ts_var[:] = np.random.uniform(143.4, 258.7, size=(len(time_avg), len(lat_array), len(lon_array)))
        ts_var.units = 'K'
        ts_var.long_name = 'surface temperature'
        
        taudust_IR_var[:] = np.random.uniform(0.0, 0.5, size=(len(time_avg), len(lat_array), len(lon_array)))
        taudust_IR_var.units = 'op'
        taudust_IR_var.long_name = 'Dust opacity IR'
        
        cldcol_var[:] = np.random.uniform(1.2e-11, 4.1e-02, size=(len(time_avg), len(lat_array), len(lon_array)))
        cldcol_var.units = ''
        cldcol_var.long_name = 'ice column'
          
# ----------------------------------------------------------------------
# 7. Create atmos_diurn_pstd.nc
# ----------------------------------------------------------------------
def create_atmos_diurn_pstd():
    filename = os.path.join(output_dir, "01336.atmos_diurn_pstd.nc")
    with Dataset(filename, 'w', format='NETCDF4') as nc:
        # Create dimensions
        time_dim = nc.createDimension('time', len(time_avg))
        time_of_day_dim = nc.createDimension('time_of_day_24', len(time_of_day_24))
        lat_dim = nc.createDimension('lat', len(lat_array))
        lon_dim = nc.createDimension('lon', len(lon_array))
        scalar_dim = nc.createDimension('scalar_axis', 1)
        
        # Create variables
        time_var = nc.createVariable('time', 'f4', ('time',))
        time_of_day_var = nc.createVariable('time_of_day_24', 'f4', ('time_of_day_24',))
        lat_var = nc.createVariable('lat', 'f4', ('lat',))
        lon_var = nc.createVariable('lon', 'f4', ('lon',))
        scalar_var = nc.createVariable('scalar_axis', 'f4', ('scalar_axis',))
        
        # Create data variables
        areo_var = nc.createVariable('areo', 'f4', ('time', 'time_of_day_24', 'scalar_axis'))
        ps_var = nc.createVariable('ps', 'f4', ('time', 'time_of_day_24', 'lat', 'lon'))
        
        # Fill coordinate variables
        time_var[:] = time_avg
        time_of_day_var[:] = time_of_day_24
        lat_var[:] = lat_array
        lon_var[:] = lon_array
        scalar_var[:] = scalar_axis
        
        # Set attributes
        time_var.units = 'days'
        time_var.long_name = 'time'
        time_of_day_var.units = 'hours since 0000-00-00 00:00:00'
        time_of_day_var.long_name = 'time of day'
        lat_var.units = 'degrees_N'
        lat_var.long_name = 'latitude'
        lon_var.units = 'degrees_E'
        lon_var.long_name = 'longitude'
        scalar_var.long_name = 'none'
        
        # Fill data variables with random values within specified ranges
        areo_var[:, :, :] = np.random.uniform(721.1, 1077.4, size=(len(time_avg), len(time_of_day_24), 1))
        areo_var.units = 'degrees'
        areo_var.long_name = 'areo'
        
        ps_var[:] = np.random.uniform(167.9, 1338.7, size=(len(time_avg), len(time_of_day_24), len(lat_array), len(lon_array)))
        ps_var.units = 'Pa'
        ps_var.long_name = 'surface pressure'

# ----------------------------------------------------------------------
# Main execution
# ----------------------------------------------------------------------
def main():
    print("Creating NetCDF files for Mars climate model...")
    
    # Create all 6 files (removed the duplicate)
    create_atmos_average()
    print("Created 01336.atmos_average.nc")
    
    create_atmos_daily()
    print("Created 01336.atmos_daily.nc")
    
    create_atmos_diurn()
    print("Created 01336.atmos_diurn.nc")
    
    create_fixed()
    print("Created 01336.fixed.nc")
    
    create_atmos_average_pstd_c48()
    print("Created 01336.atmos_average_pstd_c48.nc")
    
    create_atmos_average_pstd()
    print("Created 01336.atmos_average_pstd.nc")
    
    create_atmos_diurn_pstd()
    print("Created 01336.atmos_diurn_pstd.nc")
    
    print(f"All files created successfully in the {output_dir} directory.")

if __name__ == "__main__":
    main()