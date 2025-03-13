#!/usr/bin/env python3
"""
Integration tests for MarsFormat.py

These tests verify the functionality of MarsFormat.py for converting
different Mars climate model file formats.
"""

import os
import sys
import unittest
import shutil
import tempfile
import subprocess
import netCDF4 as nc
import numpy as np

class TestMarsFormat(unittest.TestCase):
    """Integration test suite for MarsFormat"""
    
    @classmethod
    def setUpClass(cls):
        """Set up the test environment"""
        # Create a temporary directory for test files
        cls.test_dir = os.path.join(os.path.expanduser('~'), 'MarsFormat_test')
        os.makedirs(cls.test_dir, exist_ok=True)
        
        # Project root directory
        cls.project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    
    def setUp(self):
        """Create dummy netCDF files for each test"""
        # Change to test directory
        os.chdir(self.test_dir)
        
        # Create dummy netCDF files for EMARS, OpenMARS, PCM, MarsWRF
        self.emars_file = os.path.join(self.test_dir, "emars_test.nc")
        self.openmars_file = os.path.join(self.test_dir, "openmars_test.nc")
        self.pcm_file = os.path.join(self.test_dir, "pcm_test.nc")
        self.marswrf_file = os.path.join(self.test_dir, "marswrf_test.nc")

        self.create_dummy_emars()
        self.create_dummy_openmars()
        self.create_dummy_pcm()
        self.create_dummy_marswrf()
    
    def tearDown(self):
        """Clean up any generated files after each test"""
        # Clean up any generated output files after each test
        output_files = [
            os.path.join(self.test_dir, "emars_test_daily.nc"), 
            os.path.join(self.test_dir, "emars_test_average.nc"), 
            os.path.join(self.test_dir, "emars_test_diurn.nc"),
            os.path.join(self.test_dir, "openmars_test_daily.nc"), 
            os.path.join(self.test_dir, "openmars_test_average.nc"), 
            os.path.join(self.test_dir, "openmars_test_diurn.nc"),
            os.path.join(self.test_dir, "pcm_test_daily.nc"), 
            os.path.join(self.test_dir, "pcm_test_average.nc"), 
            os.path.join(self.test_dir, "pcm_test_diurn.nc"),
            os.path.join(self.test_dir, "marswrf_test_daily.nc"), 
            os.path.join(self.test_dir, "marswrf_test_average.nc"), 
            os.path.join(self.test_dir, "marswrf_test_diurn.nc")
        ]
        for file in output_files:
            if os.path.exists(file):
                os.remove(file)
        
        # Clean up dummy input files
        # dummy_files = [
        #     self.emars_file, self.openmars_file, self.pcm_file, self.marswrf_file
        # ]
        # for file in dummy_files:
        #     if os.path.exists(file):
        #         os.remove(file)
    
    @classmethod
    def tearDownClass(cls):
        """Clean up the test environment"""
        try:
            shutil.rmtree(cls.test_dir, ignore_errors=True)
        except Exception:
            print(f"Warning: Could not remove test directory {cls.test_dir}")
    
    def run_mars_format(self, args):
        """
        Run MarsFormat using subprocess
        
        :param args: List of arguments to pass to MarsFormat
        :return: The subprocess result object
        """
        # Construct the full command to run MarsFormat
        cmd = [sys.executable, os.path.join(self.project_root, "bin", "MarsFormat.py")] + args
        
        # Run the command
        result = subprocess.run(
            cmd, 
            capture_output=True, 
            text=True, 
            cwd=self.test_dir,  # Run in the test directory
            env=dict(os.environ, PWD=self.test_dir)  # Ensure current working directory is set
        )
        
        return result
        
    def create_dummy_emars(self):
        """Create a dummy EMARS netCDF file with the necessary attributes."""
        dataset = nc.Dataset(self.emars_file, 'w', format='NETCDF4')
        
        # Create dimensions
        dataset.createDimension('time', 840)
        dataset.createDimension('lat', 36)
        dataset.createDimension('lon', 60)
        dataset.createDimension('pfull', 28)
        dataset.createDimension('latu', 36)
        dataset.createDimension('lonv', 60)
        dataset.createDimension('phalf', 29)

        # Create variables and fill with realistic data
        lat = dataset.createVariable('lat', 'f4', ('lat',))
        lat[:] = np.linspace(-89.0, 89.0, 36)
        lat.units = 'degree_N'
        # Note the change from longname to long_name - this is critical
        lat.long_name = 'latitude'

        lon = dataset.createVariable('lon', 'f4', ('lon',))
        lon[:] = np.linspace(3.0, 357.0, 60)
        lon.units = 'degree_E'
        lon.long_name = 'longitude'

        latu = dataset.createVariable('latu', 'f4', ('latu',))
        latu[:] = np.linspace(-87.4, 87.4, 36)
        latu.units = 'degree_N'
        latu.long_name = 'latitude at u points'

        lonv = dataset.createVariable('lonv', 'f4', ('lonv',))
        lonv[:] = np.linspace(0.0, 354.0, 60)
        lonv.units = 'degree_E'
        lonv.long_name = 'longitude at v points'

        pfull = dataset.createVariable('pfull', 'f4', ('pfull',))
        pfull[:] = np.linspace(0.0, 7.7, 28)
        pfull.units = 'mb'
        pfull.long_name = 'ref full pressure level'

        phalf = dataset.createVariable('phalf', 'f4', ('phalf',))
        phalf[:] = np.linspace(0.0, 7.7, 29)
        phalf.units = 'mb'
        phalf.long_name = 'ref half pressure level'

        ak = dataset.createVariable('ak', 'f4', ('phalf',))
        ak[:] = np.linspace(0.0, 2.9, 29)
        ak.units = 'pascal'
        ak.long_name = 'pressure part of the hybrid coordinate'

        bk = dataset.createVariable('bk', 'f4', ('phalf',))
        bk[:] = np.linspace(0.0, 1.0, 29)
        bk.units = 'none'
        bk.long_name = 'vertical coordinate sigma value'

        time = dataset.createVariable('time', 'f4', ('time',))
        time[:] = np.linspace(0.0, 839.0, 840)
        time.units = 'Martian hour'
        time.long_name = 'number of hours since start of file'

        Ls = dataset.createVariable('Ls', 'f4', ('time',))
        Ls[:] = np.linspace(103.5, 120.0, 840)
        Ls.units = 'deg'
        Ls.long_name = 'areocentric longitude'

        MY = dataset.createVariable('MY', 'f4', ('time',))
        MY[:] = np.ones(840) * 24.0
        MY.units = 'Martian year'
        MY.long_name = 'Mars Year'

        earth_day = dataset.createVariable('earth_day', 'f4', ('time',))
        earth_day[:] = np.random.randint(1, 31, 840)
        earth_day.units = 'Earth day'
        earth_day.long_name = 'Earth day of the month'

        earth_hour = dataset.createVariable('earth_hour', 'f4', ('time',))
        earth_hour[:] = np.random.randint(0, 23, 840)
        earth_hour.units = 'Earth hour'
        earth_hour.long_name = 'Earth hour of the day'

        earth_minute = dataset.createVariable('earth_minute', 'f4', ('time',))
        earth_minute[:] = np.random.randint(0, 59, 840)
        earth_minute.units = 'Earth minute'
        earth_minute.long_name = 'Earth minute of the hour'

        earth_month = dataset.createVariable('earth_month', 'f4', ('time',))
        earth_month[:] = np.random.randint(2, 4, 840)
        earth_month.units = 'Earth month'
        earth_month.long_name = 'Earth month of the year'

        earth_second = dataset.createVariable('earth_second', 'f4', ('time',))
        earth_second[:] = np.random.uniform(0.1, 60.0, 840)
        earth_second.units = 'Earth second'
        earth_second.long_name = 'Earth second and fractional second of the minute'

        earth_year = dataset.createVariable('earth_year', 'f4', ('time',))
        earth_year[:] = np.ones(840) * 1999.0
        earth_year.units = 'Earth year'
        earth_year.long_name = 'Earth year AD'

        emars_sol = dataset.createVariable('emars_sol', 'f4', ('time',))
        emars_sol[:] = np.linspace(1075.0, 1109.0, 840)
        emars_sol.units = 'Martian sol'
        emars_sol.long_name = 'sols after MY 22 perihelion'

        macda_sol = dataset.createVariable('macda_sol', 'f4', ('time',))
        macda_sol[:] = np.linspace(223.0, 257.0, 840)
        macda_sol.units = 'Martian year'
        macda_sol.long_name = 'sols after the start of MY 24'

        mars_hour = dataset.createVariable('mars_hour', 'f4', ('time',))
        mars_hour[:] = np.random.randint(0, 23, 840)
        mars_hour.units = 'Martian hour'
        mars_hour.long_name = 'hour of the Martian day'

        mars_soy = dataset.createVariable('mars_soy', 'f4', ('time',))
        mars_soy[:] = np.linspace(223.0, 257.0, 840)
        mars_soy.units = 'Martian sol'
        mars_soy.long_name = 'sols after the last Martian vernal equinox'

        Surface_geopotential = dataset.createVariable('Surface_geopotential', 'f4', ('lat', 'lon',))
        Surface_geopotential[:] = np.random.uniform(-24000, 26000, (36, 60))
        Surface_geopotential.units = 'm^2/s/s'
        Surface_geopotential.long_name = 'surface geopotential height'

        ps = dataset.createVariable('ps', 'f4', ('time', 'lat', 'lon',))
        ps[:] = np.random.uniform(285.5, 1295.8, (840, 36, 60))
        ps.units = 'pascal'
        ps.long_name = 'surface pressure'

        T = dataset.createVariable('T', 'f4', ('time', 'pfull', 'lat', 'lon',))
        T[:] = np.random.uniform(105.7, 253.7, (840, 28, 36, 60))
        T.units = 'K'
        T.long_name = 'Temperature'

        U = dataset.createVariable('U', 'f4', ('time', 'pfull', 'latu', 'lon',))
        U[:] = np.random.uniform(-159.6, 198.1, (840, 28, 36, 60))
        U.units = 'm/s'
        U.long_name = 'zonal wind'

        V = dataset.createVariable('V', 'f4', ('time', 'pfull', 'lat', 'lonv',))
        V[:] = np.random.uniform(-197.3, 158.3, (840, 28, 36, 60))
        V.units = 'm/s'
        V.long_name = 'meridional wind'

        dataset.close()

    def create_dummy_openmars(self):
        dataset = nc.Dataset(self.openmars_file, 'w', format='NETCDF4')

        # Create dimensions based on GCM_netcdf_file_structures.txt
        dataset.createDimension('lon', 60)
        dataset.createDimension('lat', 36)
        dataset.createDimension('lev', 30)
        dataset.createDimension('time', 30)

        # Create variables and fill with realistic data
        lon = dataset.createVariable('lon', 'f4', ('lon',))
        lon[:] = np.linspace(-180.0, 175.0, 60)
        lon.FIELDNAM = 'longitude'
        lon.UNITS = 'degrees_east'
        
        lat = dataset.createVariable('lat', 'f4', ('lat',))
        lat[:] = np.linspace(-87.0, 87.0, 36)
        lat.FIELDNAM = 'latitude'
        lat.UNITS = 'degrees_north'
        
        lev = dataset.createVariable('lev', 'f4', ('lev',))
        lev[:] = np.linspace(5.1e-05, 1.0e+00, 30)
        lev.FIELDNAM = 'sigma level'
        lev.UNITS = 'none'
        
        time = dataset.createVariable('time', 'f4', ('time',))
        time[:] = np.linspace(2911.1, 2941.0, 30)
        time.FIELDNAM = 'time'
        time.UNITS = 'hours'
        
        Ls = dataset.createVariable('Ls', 'f4', ('time',))
        Ls[:] = np.linspace(110.0, 124.3, 30)
        Ls.FIELDNAM = 'solar longitude'
        Ls.UNITS = 'degrees'
        
        MY = dataset.createVariable('MY', 'f4', ('time',))
        MY[:] = np.ones(30) * 28.0
        MY.FIELDNAM = 'Mars Year'
        MY.UNITS = 'year'
        
        ps = dataset.createVariable('ps', 'f4', ('time', 'lat', 'lon',))
        ps[:] = np.random.uniform(178.2, 1347.1, (30, 36, 60))
        ps.FIELDNAM = 'surface pressure'
        ps.UNITS = 'Pa'
        
        tsurf = dataset.createVariable('tsurf', 'f4', ('time', 'lat', 'lon',))
        tsurf[:] = np.random.uniform(142.2, 288.2, (30, 36, 60))
        tsurf.FIELDNAM = 'surface temperature'
        tsurf.UNITS = 'K'
        
        co2ice = dataset.createVariable('co2ice', 'f4', ('time', 'lat', 'lon',))
        co2ice[:] = np.random.uniform(0.0, 6790.3, (30, 36, 60))
        co2ice.FIELDNAM = 'CO2 ice'
        co2ice.UNITS = 'kg/m2'
        
        dustcol = dataset.createVariable('dustcol', 'f4', ('time', 'lat', 'lon',))
        dustcol[:] = np.random.uniform(0.0, 1.1, (30, 36, 60))
        dustcol.FIELDNAM = 'dust column'
        dustcol.UNITS = 'none'
        
        u = dataset.createVariable('u', 'f4', ('time', 'lev', 'lat', 'lon',))
        u[:] = np.random.uniform(-194.2, 211.8, (30, 30, 36, 60))
        u.FIELDNAM = 'zonal wind'
        u.UNITS = 'm/s'
        
        v = dataset.createVariable('v', 'f4', ('time', 'lev', 'lat', 'lon',))
        v[:] = np.random.uniform(-203.3, 181.0, (30, 30, 36, 60))
        v.FIELDNAM = 'meridional wind'
        v.UNITS = 'm/s'
        
        temp = dataset.createVariable('temp', 'f4', ('time', 'lev', 'lat', 'lon',))
        temp[:] = np.random.uniform(99.1, 258.7, (30, 30, 36, 60))
        temp.FIELDNAM = 'temperature'
        temp.UNITS = 'K'
        
        dataset.close()

    def create_dummy_pcm(self):
        dataset = nc.Dataset(self.pcm_file, 'w', format='NETCDF4')

        # Create dimensions based on PCM_diagfi_small.nc
        dataset.createDimension('time', 10)
        dataset.createDimension('altitude', 24)
        dataset.createDimension('latitude', 32)
        dataset.createDimension('longitude', 64)
        dataset.createDimension('interlayer', 24)
        dataset.createDimension('subsurface_layers', 10)
        dataset.createDimension('index', 10)

        # Create variables and fill with realistic data
        time = dataset.createVariable('time', 'f4', ('time',))
        time[:] = np.linspace(0.5, 3.0, 10)
        time.units = 'hours since start of simulation'
        time.title = 'Time'

        Ls = dataset.createVariable('Ls', 'f4', ('time',))
        Ls[:] = np.linspace(9.3, 10.6, 10)
        Ls.units = 'deg'
        Ls.title = 'Solar longitude'

        Sols = dataset.createVariable('Sols', 'f4', ('time',))
        Sols[:] = np.linspace(687.5, 690.0, 10)
        Sols.units = 'sols'
        Sols.title = 'Mars sols'

        altitude = dataset.createVariable('altitude', 'f4', ('altitude',))
        altitude[:] = np.linspace(0.0, 80.5, 24)
        altitude.units = 'km'
        altitude.title = 'Altitude'

        latitude = dataset.createVariable('latitude', 'f4', ('latitude',))
        latitude[:] = np.linspace(-90.0, 90.0, 32)
        latitude.units = 'degrees_north'
        latitude.title = 'North latitude'

        longitude = dataset.createVariable('longitude', 'f4', ('longitude',))
        longitude[:] = np.linspace(-180.0, 180.0, 64)
        longitude.units = 'degrees_east'
        longitude.title = 'East longitude'

        aps = dataset.createVariable('aps', 'f4', ('altitude',))
        aps[:] = np.linspace(0.0, 9.0, 24)
        aps.units = 'Pa'
        aps.title = 'Pressure levels a coefficient'

        bps = dataset.createVariable('bps', 'f4', ('altitude',))
        bps[:] = np.linspace(0.0, 1.0, 24)
        bps.units = ''
        bps.title = 'Pressure levels b coefficient'

        temp = dataset.createVariable('temp', 'f4', ('time', 'altitude', 'latitude', 'longitude',))
        temp[:] = np.random.uniform(108.7, 250.7, (10, 24, 32, 64))
        temp.units = 'K'
        temp.title = 'Temperature'

        ps = dataset.createVariable('ps', 'f4', ('time', 'latitude', 'longitude',))
        ps[:] = np.random.uniform(149.3, 1198.9, (10, 32, 64))
        ps.units = 'Pa'
        ps.title = 'Surface pressure'

        u = dataset.createVariable('u', 'f4', ('time', 'altitude', 'latitude', 'longitude',))
        u[:] = np.random.uniform(-124.5, 191.5, (10, 24, 32, 64))
        u.units = 'm.s-1'
        u.title = 'Zonal wind'

        v = dataset.createVariable('v', 'f4', ('time', 'altitude', 'latitude', 'longitude',))
        v[:] = np.random.uniform(-96.9, 105.6, (10, 24, 32, 64))
        v.units = 'm.s-1'
        v.title = 'Meridional wind'

        h2o_ice = dataset.createVariable('h2o_ice', 'f4', ('time', 'altitude', 'latitude', 'longitude',))
        h2o_ice[:] = np.random.uniform(-3.2e-21, 3.2e-04, (10, 24, 32, 64))
        h2o_ice.units = 'kg/kg'
        h2o_ice.title = 'Water ice mass mixing ratio'

        h2o_vap = dataset.createVariable('h2o_vap', 'f4', ('time', 'altitude', 'latitude', 'longitude',))
        h2o_vap[:] = np.random.uniform(-2.3e-23, 3.2e-04, (10, 24, 32, 64))
        h2o_vap.units = 'kg/kg'
        h2o_vap.title = 'Water vapor mass mixing ratio'

        dataset.close()

    def create_dummy_marswrf(self):
        """Create a MarsWRF file with explicitly defined coordinate arrays based on real MarsWRF data."""
        dataset = nc.Dataset(self.marswrf_file, 'w', format='NETCDF4')
        
        # Create dimensions according to actual MarsWRF format
        dataset.createDimension('Time', 5)
        dataset.createDimension('DateStrLen', 19)  
        dataset.createDimension('south_north', 86)
        dataset.createDimension('west_east', 180)
        dataset.createDimension('south_north_stag', 87)
        dataset.createDimension('west_east_stag', 181)
        dataset.createDimension('bottom_top', 10)
        dataset.createDimension('bottom_top_stag', 11)
        dataset.createDimension('soil_layers_stag', 5)
        
        # Create Times variable with string dates
        Times = dataset.createVariable('Times', 'S1', ('Time', 'DateStrLen'))
        time_strings = []
        for i in range(5):
            time_str = f"2020-01-{(i+1):02d}_00:00:00"
            time_strings.append([c for c in time_str])
        Times[:] = time_strings
        
        # Create required time variables
        XTIME = dataset.createVariable('XTIME', 'f4', ('Time',))
        XTIME[:] = np.arange(5) * 24.0 * 60.0  # minutes
        XTIME.units = 'minutes'
        XTIME.description = 'minutes since simulation start'
        
        # Create explicit coordinate arrays exactly as seen in real MarsWRF files
        
        # Create latitude arrays (XLAT and staggered versions)
        XLAT = dataset.createVariable('XLAT', 'f4', ('Time', 'south_north', 'west_east'))
        XLAT.units = 'degree_north'
        XLAT.description = 'latitude'
        
        # Fill XLAT properly - constant along west_east dimension
        for t in range(5):
            lat_values = np.linspace(-89, 89, 86)
            for i, lat in enumerate(lat_values):
                XLAT[t, i, :] = lat
        
        # Create XLAT_U - same as XLAT but one more point in west_east dimension
        XLAT_U = dataset.createVariable('XLAT_U', 'f4', ('Time', 'south_north', 'west_east_stag'))
        XLAT_U.units = 'degree_north'
        XLAT_U.description = 'latitude at u-points'
        
        # XLAT_U matches XLAT in first 180 points, repeats first value at end
        for t in range(5):
            for i in range(86):
                XLAT_U[t, i, 0:180] = XLAT[t, i, :]
                XLAT_U[t, i, 180] = XLAT[t, i, 0]
        
        # Create XLAT_V - staggered in south_north
        XLAT_V = dataset.createVariable('XLAT_V', 'f4', ('Time', 'south_north_stag', 'west_east'))
        XLAT_V.units = 'degree_north'
        XLAT_V.description = 'latitude at v-points'
        
        # First fill south pole values
        for t in range(5):
            XLAT_V[t, 0, :] = -90.0
            
            # Fill interior points as averages
            for i in range(1, 86):
                XLAT_V[t, i, :] = (XLAT[t, i-1, :] + XLAT[t, i, :]) / 2.0
                
            # Fill north pole values
            XLAT_V[t, 86, :] = 90.0
        
        # Create longitude arrays (XLONG and staggered versions)
        XLONG = dataset.createVariable('XLONG', 'f4', ('Time', 'south_north', 'west_east'))
        XLONG.units = 'degree_east'
        XLONG.description = 'longitude'
        
        # Fill XLONG properly
        # Create array from 1 to 179, then -179 to -1
        lon_values = np.concatenate([np.arange(1, 180, 2), np.arange(-179, 0, 2)])
        for t in range(5):
            for i in range(86):
                XLONG[t, i, :] = lon_values
        
        # Create XLONG_U - staggered in west_east
        XLONG_U = dataset.createVariable('XLONG_U', 'f4', ('Time', 'south_north', 'west_east_stag'))
        XLONG_U.units = 'degree_east'
        XLONG_U.description = 'longitude at u-points'
        
        # XLONG_U contains midpoints between XLONG values, starting at 0
        lon_u_values = np.concatenate([np.arange(0, 180, 2), np.arange(-178, 2, 2), [0]])
        for t in range(5):
            for i in range(86):
                XLONG_U[t, i, :] = lon_u_values
        
        # Create XLONG_V - same as XLONG in west_east dimension
        XLONG_V = dataset.createVariable('XLONG_V', 'f4', ('Time', 'south_north_stag', 'west_east'))
        XLONG_V.units = 'degree_east'
        XLONG_V.description = 'longitude at v-points'
        
        # Copy XLONG to XLONG_V for each north_south_stag point
        for t in range(5):
            for i in range(87):
                if i < 86:
                    XLONG_V[t, i, :] = XLONG[t, i, :]
                else:
                    XLONG_V[t, i, :] = XLONG[t, 0, :]  # Use first latitude for extra staggered point
        
        # Required thermodynamic constants
        P_TOP = dataset.createVariable('P_TOP', 'f4', ('Time',))
        P_TOP[:] = np.ones(5) * 5.0
        P_TOP.units = 'Pa'
        P_TOP.description = 'pressure at model top'
        
        P0 = dataset.createVariable('P0', 'f4', ())
        P0[:] = 610.0
        P0.units = 'Pa'
        P0.description = 'reference pressure'
        
        G = dataset.createVariable('G', 'f4', ())
        G[:] = 3.72
        G.units = 'm s-2'
        G.description = 'gravitational acceleration'
        
        CP = dataset.createVariable('CP', 'f4', ())
        CP[:] = 735.0
        CP.units = 'J kg-1 K-1'
        CP.description = 'specific heat capacity'
        
        R_D = dataset.createVariable('R_D', 'f4', ())
        R_D[:] = 192.0
        R_D.units = 'J kg-1 K-1'
        R_D.description = 'gas constant for dry air'
        
        T0 = dataset.createVariable('T0', 'f4', ('Time',))
        T0[:] = np.ones(5) * 170.0
        T0.units = 'K'
        T0.description = 'reference temperature'
        
        # Create vertical coordinate variables
        ZNU = dataset.createVariable('ZNU', 'f4', ('Time', 'bottom_top'))
        ZNU[:] = np.tile(np.linspace(0, 1, 10).reshape(1, 10), (5, 1))
        ZNU.description = 'eta values on full (mass) levels'
        
        ZNW = dataset.createVariable('ZNW', 'f4', ('Time', 'bottom_top_stag'))
        ZNW[:] = np.tile(np.linspace(0, 1, 11).reshape(1, 11), (5, 1))
        ZNW.description = 'eta values on half (w) levels'
        
        # Core meteorological variables - minimal set for testing
        T = dataset.createVariable('T', 'f4', ('Time', 'bottom_top', 'south_north', 'west_east'))
        T[:] = np.random.uniform(-170.9, 349.1, (5, 10, 86, 180))
        T.units = 'K'
        T.description = 'perturbation potential temperature'
        
        PSFC = dataset.createVariable('PSFC', 'f4', ('Time', 'south_north', 'west_east'))
        PSFC[:] = np.random.uniform(85.0, 1330.0, (5, 86, 180))
        PSFC.units = 'Pa'
        PSFC.description = 'surface pressure'
        
        # Wind variables (staggered) - minimal set
        U = dataset.createVariable('U', 'f4', ('Time', 'bottom_top', 'south_north', 'west_east_stag'))
        U[:] = np.random.uniform(-514.7, 519.5, (5, 10, 86, 181))
        U.units = 'm s-1'
        U.description = 'x-wind component'
        
        V = dataset.createVariable('V', 'f4', ('Time', 'bottom_top', 'south_north_stag', 'west_east'))
        V[:] = np.random.uniform(-205.6, 238.3, (5, 10, 87, 180))
        V.units = 'm s-1'
        V.description = 'y-wind component'
        
        W = dataset.createVariable('W', 'f4', ('Time', 'bottom_top_stag', 'south_north', 'west_east'))
        W[:] = np.random.uniform(-40.7, 39.0, (5, 11, 86, 180))
        W.units = 'm s-1'
        W.description = 'z-wind component'
        
        # Surface variables - minimal set
        HGT = dataset.createVariable('HGT', 'f4', ('Time', 'south_north', 'west_east'))
        HGT[:] = np.random.uniform(-7400.0, 19000.0, (5, 86, 180))
        HGT.units = 'm'
        HGT.description = 'terrain height'
        
        TSK = dataset.createVariable('TSK', 'f4', ('Time', 'south_north', 'west_east'))
        TSK[:] = np.random.uniform(142.0, 308.7, (5, 86, 180))
        TSK.units = 'K'
        TSK.description = 'surface skin temperature'
        
        # Mars-specific variables - minimal set
        L_S = dataset.createVariable('L_S', 'f4', ('Time',))
        L_S[:] = np.linspace(0, 90, 5)
        L_S.units = 'degrees'
        L_S.description = 'Mars solar longitude'
        
        CO2ICE = dataset.createVariable('CO2ICE', 'f4', ('Time', 'south_north', 'west_east'))
        CO2ICE[:] = np.random.uniform(0.0, 1697.5, (5, 86, 180))
        CO2ICE.units = 'kg/m^2'
        CO2ICE.description = 'CO2 ice on surface'
        
        TAU_OD2D = dataset.createVariable('TAU_OD2D', 'f4', ('Time', 'south_north', 'west_east'))
        TAU_OD2D[:] = np.random.uniform(0.0, 2.1, (5, 86, 180))
        TAU_OD2D.units = 'Unitless'
        TAU_OD2D.description = 'Dust optical depth'
        
        # Add geopotential variables - THESE WERE MISSING
        PH = dataset.createVariable('PH', 'f4', ('Time', 'bottom_top_stag', 'south_north', 'west_east'))
        PH[:] = np.random.uniform(-7.3e+04, 1.6e+04, (5, 11, 86, 180))
        PH.units = 'm2 s-2'
        PH.description = 'perturbation geopotential'
        
        PHB = dataset.createVariable('PHB', 'f4', ('Time', 'bottom_top_stag', 'south_north', 'west_east'))
        PHB[:] = np.random.uniform(-2.8e+04, 2.7e+05, (5, 11, 86, 180))
        PHB.units = 'm2 s-2'
        PHB.description = 'base-state geopotential'
        
        # Add pressure-related variables we may have missed
        MU = dataset.createVariable('MU', 'f4', ('Time', 'south_north', 'west_east'))
        MU[:] = np.random.uniform(-267.5, 87.1, (5, 86, 180))
        MU.units = 'Pa'
        MU.description = 'perturbation dry air mass in column'
        
        MUB = dataset.createVariable('MUB', 'f4', ('Time', 'south_north', 'west_east'))
        MUB[:] = np.random.uniform(128.4, 1244.7, (5, 86, 180))
        MUB.units = 'Pa'
        MUB.description = 'base state dry air mass in column'
        
        P = dataset.createVariable('P', 'f4', ('Time', 'bottom_top', 'south_north', 'west_east'))
        P[:] = np.random.uniform(-267.2, 87.0, (5, 10, 86, 180))
        P.units = 'Pa'
        P.description = 'perturbation pressure'
        
        PB = dataset.createVariable('PB', 'f4', ('Time', 'bottom_top', 'south_north', 'west_east'))
        PB[:] = np.random.uniform(1.6, 1243.7, (5, 10, 86, 180))
        PB.units = 'Pa'
        PB.description = 'base state pressure'
        
        dataset.close()
    
    def test_daily_average(self):
        """Test the daily average functionality of MarsFormat.py for EMARS file."""
        # Run MarsFormat with the correct arguments
        result = self.run_mars_format([os.path.basename(self.emars_file), "-gcm", "emars"])
        
        # Check that the command executed successfully
        self.assertEqual(result.returncode, 0, f"MarsFormat.py failed: {result.stderr}")
        
        # Check that the output file was created
        output_file = os.path.join(self.test_dir, "emars_test_daily.nc")
        self.assertTrue(os.path.exists(output_file), f"Output file {output_file} was not created.")
        
        # Open the output file and check that it contains the expected variables
        dataset = nc.Dataset(output_file, 'r')
        
        # Debug - print all variable names to examine what's actually in the file
        print(f"Variables in {output_file}: {list(dataset.variables.keys())}")
        
        # Check that key variables are present - note that 'T' may have been renamed to 'temp'
        self.assertIn('temp', dataset.variables, "Temperature variable not found in output file.")
        self.assertIn('ps', dataset.variables, "Surface pressure variable not found in output file.")
        
        dataset.close()

    def test_diurn_average(self):
        """Test the diurnal average functionality of MarsFormat.py for EMARS file."""
        # Run MarsFormat with the correct arguments
        result = self.run_mars_format([os.path.basename(self.emars_file), "-gcm", "emars", "-bd"])
        
        # Check that the command executed successfully
        self.assertEqual(result.returncode, 0, f"MarsFormat.py failed: {result.stderr}")
        
        # Check that the output file was created
        output_file = os.path.join(self.test_dir, "emars_test_diurn.nc")
        self.assertTrue(os.path.exists(output_file), f"Output file {output_file} was not created.")
        
        # Open the output file and check appropriate variables
        dataset = nc.Dataset(output_file, 'r')
        
        # Debug - print all variable names 
        print(f"Variables in {output_file}: {list(dataset.variables.keys())}")
        
        # Check for key variables we expect to see in diurnal average - note variable names
        self.assertIn('temp', dataset.variables, "Temperature variable not found in output file.")
        self.assertIn('ps', dataset.variables, "Surface pressure variable not found in output file.")
        
        dataset.close()

    def test_openmars_format(self):
        """Test MarsFormat.py with OpenMARS file."""
        # Run MarsFormat with the correct arguments
        result = self.run_mars_format([os.path.basename(self.openmars_file), "-gcm", "openmars"])
        
        # Check that the command executed successfully
        self.assertEqual(result.returncode, 0, f"MarsFormat.py failed with OpenMARS: {result.stderr}")
        
        # Check that the output file was created
        output_file = os.path.join(self.test_dir, "openmars_test_daily.nc")
        self.assertTrue(os.path.exists(output_file), f"Output file {output_file} was not created.")
        
        # Check output contains expected variables
        dataset = nc.Dataset(output_file, 'r')
        self.assertIn('temp', dataset.variables, "Temperature variable not found in output file.")
        self.assertIn('ps', dataset.variables, "Surface pressure variable not found in output file.")
        dataset.close()

    def test_pcm_format(self):
        """Test MarsFormat.py with PCM file."""
        # Run MarsFormat with the correct arguments
        result = self.run_mars_format([os.path.basename(self.pcm_file), "-gcm", "pcm"])
        
        # Check that the command executed successfully
        self.assertEqual(result.returncode, 0, f"MarsFormat.py failed with PCM: {result.stderr}")
        
        # Check that the output file was created
        output_file = os.path.join(self.test_dir, "pcm_test_daily.nc")
        self.assertTrue(os.path.exists(output_file), f"Output file {output_file} was not created.")
        
        # Check output contains expected variables
        dataset = nc.Dataset(output_file, 'r')
        self.assertIn('temp', dataset.variables, "Temperature variable not found in output file.")
        self.assertIn('ps', dataset.variables, "Surface pressure variable not found in output file.")
        dataset.close()

    def test_marswrf_format(self):
        """Test MarsFormat.py with MarsWRF file."""
        # Create a simplified MarsWRF test that's less likely to cause coordinate issues
        self.create_dummy_marswrf()
        
        # Run MarsFormat with the correct arguments
        result = self.run_mars_format([os.path.basename(self.marswrf_file), "-gcm", "marswrf"])
        
        # Check that the command executed successfully - if it fails, print more detailed error
        if result.returncode != 0:
            print(f"STDERR: {result.stderr}")
            print(f"STDOUT: {result.stdout}")
        
        self.assertEqual(result.returncode, 0, f"MarsFormat.py failed with MarsWRF: {result.stderr}")
        
        # Check that the output file was created
        output_file = os.path.join(self.test_dir, "marswrf_test_daily.nc")
        self.assertTrue(os.path.exists(output_file), f"Output file {output_file} was not created.")
        
        # Check output contains expected variables - note that variable names may have changed
        dataset = nc.Dataset(output_file, 'r')
        print(f"Variables in {output_file}: {list(dataset.variables.keys())}")
        
        # Look for appropriate variables, flexibly handling variable name changes
        temp_var_found = False
        for var_name in dataset.variables:
            if 'temp' in var_name.lower() or 't' == var_name:
                temp_var_found = True
                break
        
        self.assertTrue(temp_var_found, "No temperature variable found in output file")
        
        # Check for pressure variable similarly
        pressure_var_found = False
        for var_name in dataset.variables:
            if 'ps' in var_name.lower() or 'psfc' in var_name.lower():
                pressure_var_found = True
                break
        
        self.assertTrue(pressure_var_found, "No pressure variable found in output file")
        
        dataset.close()


if __name__ == '__main__':
    unittest.main()