#!/usr/bin/env python3
"""
Integration tests for MarsFiles.py

These tests verify the functionality of MarsFiles for manipulating netCDF files.
"""

import os
import sys
import unittest
import shutil
import subprocess
import tempfile
import glob
import numpy as np
from netCDF4 import Dataset

class TestMarsFiles(unittest.TestCase):
    """Integration test suite for MarsFiles"""
    
    @classmethod
    def setUpClass(cls):
        """Set up the test environment"""
        # Create a temporary directory for the tests
        cls.test_dir = tempfile.mkdtemp(prefix='MarsFiles_test_')
        
        # Project root directory
        cls.project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        
        # Run the script to create test netCDF files
        cls.create_test_files()
    
    @classmethod
    def create_test_files(cls):
        """Create test netCDF files using create_ames_gcm_files.py"""
        # Get path to create_ames_gcm_files.py script
        create_files_script = os.path.join(cls.project_root, "tests", "create_ames_gcm_files.py")
        
        # Execute the script to create test files
        cmd = [sys.executable, create_files_script, cls.test_dir]
        
        # Print the command being executed
        print(f"Creating test files with command: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                cwd=cls.project_root  # Run in the project root directory
            )
            
            # Print output for debugging
            print(f"File creation output: {result.stdout}")
            
            if result.returncode != 0:
                raise Exception(f"Failed to create test files: {result.stderr}")
            
        except Exception as e:
            raise Exception(f"Error running create_ames_gcm_files.py: {e}")
            
        # Verify files were created
        expected_files = [
            '01336.atmos_average.nc',
            '01336.atmos_average_pstd_c48.nc',
            '01336.atmos_daily.nc',
            '01336.atmos_diurn.nc',
            '01336.atmos_diurn_pstd.nc',
            '01336.fixed.nc'
        ]
        
        for filename in expected_files:
            filepath = os.path.join(cls.test_dir, filename)
            if not os.path.exists(filepath):
                raise Exception(f"Test file {filename} was not created")
    
    def setUp(self):
        """Change to temporary directory before each test"""
        os.chdir(self.test_dir)
    
    def tearDown(self):
        """Clean up after each test"""
        # Clean up any generated output files after each test but keep input files
        output_patterns = [
            '*_T.nc',
            '*_to_average.nc',
            '*_to_diurn.nc',
            '*_tide_decomp*.nc',
            '*_hpt*.nc',
            '*_lpt*.nc',
            '*_bpt*.nc',
            '*_regrid.nc',
            '*_zavg*.nc',
            '*_Ls*_*.nc',
            '*_lat_*_*.nc'
        ]
        
        for pattern in output_patterns:
            for file_path in glob.glob(os.path.join(self.test_dir, pattern)):
                try:
                    os.remove(file_path)
                    print(f"Removed file: {file_path}")
                except Exception as e:
                    print(f"Warning: Could not remove file {file_path}: {e}")
        
        # Return to test_dir
        os.chdir(self.test_dir)
    
    @classmethod
    def tearDownClass(cls):
        """Clean up the test environment"""
        try:
            shutil.rmtree(cls.test_dir, ignore_errors=True)
        except Exception:
            print(f"Warning: Could not remove test directory {cls.test_dir}")
    
    def run_mars_files(self, args):
        """
        Run MarsFiles using subprocess
        
        :param args: List of arguments to pass to MarsFiles
        :return: subprocess result object
        """
        # Convert any relative file paths to absolute paths
        abs_args = []
        for arg in args:
            if isinstance(arg, str) and arg.endswith('.nc'):
                abs_args.append(os.path.join(self.test_dir, arg))
            else:
                abs_args.append(arg)
        
        # Construct the full command to run MarsFiles
        cmd = [sys.executable, '-m', 'bin.MarsFiles'] + abs_args
        
        # Print the command being executed
        print(f"Running command: {' '.join(cmd)}")
        print(f"Working directory: {self.test_dir}")
        
        # Run the command
        try:
            result = subprocess.run(
                cmd, 
                capture_output=True, 
                text=True, 
                cwd=self.project_root,  # Run in the project root directory
                env=dict(os.environ, PWD=self.test_dir)  # Ensure current working directory is set
            )
            
            # Print both stdout and stderr to help debug
            print(f"STDOUT: {result.stdout}")
            print(f"STDERR: {result.stderr}")
            
            return result
        except Exception as e:
            self.fail(f"Failed to run MarsFiles: {e}")
    
    def check_file_exists(self, filename):
        """
        Check if a file exists and is not empty
        
        :param filename: Filename to check
        """
        filepath = os.path.join(self.test_dir, filename)
        self.assertTrue(os.path.exists(filepath), f"File {filename} does not exist")
        self.assertGreater(os.path.getsize(filepath), 0, f"File {filename} is empty")
        return filepath
    
    def verify_netcdf_has_variable(self, filename, variable):
        """
        Verify that a netCDF file has a specific variable
        
        :param filename: Path to the netCDF file
        :param variable: Variable name to check for
        """
        nc = Dataset(filename, 'r')
        try:
            self.assertIn(variable, nc.variables, f"Variable {variable} not found in {filename}")
        finally:
            nc.close()
    
    def test_help_message(self):
        """Test that help message can be displayed"""
        result = self.run_mars_files(['-h'])
        
        # Check for successful execution
        self.assertEqual(result.returncode, 0, "Help command failed")
        
        # Check for typical help message components
        help_checks = [
            'usage:',
            '--bin_files',
            '--concatenate',
            '--split',
            '--time_shift',
            '--bin_average',
            '--bin_diurn',
            '--tide_decomp',
            '--regrid',
            '--zonal_average'
        ]
        
        for check in help_checks:
            self.assertIn(check, result.stdout, f"Help message missing '{check}'")
    
    def test_time_shift(self):
        """Test time_shift operation on diurn file"""
        result = self.run_mars_files(['01336.atmos_diurn_pstd.nc', '-t'])
        
        # Check for successful execution
        self.assertEqual(result.returncode, 0, "Time shift command failed")
        
        # Check that output file was created
        output_file = self.check_file_exists('01336.atmos_diurn_pstd_T.nc')
        
        # Verify that the output file has expected structure
        self.verify_netcdf_has_variable(output_file, 'time')
        self.verify_netcdf_has_variable(output_file, 'time_of_day_24')  # Default is 24 time of day bins
    
    def test_time_shift_specific_times(self):
        """Test time_shift operation with specific local times"""
        result = self.run_mars_files(['01336.atmos_diurn_pstd.nc', '-t', '3 15'])
        
        # Check for successful execution
        self.assertEqual(result.returncode, 0, "Time shift with specific times command failed")
        
        # Check that output file was created
        output_file = self.check_file_exists('01336.atmos_diurn_pstd_T.nc')
        
        # Verify that the output file has expected structure with only 2 time of day values
        nc = Dataset(output_file, 'r')
        try:
            # Should have 'time_of_day_02' dimension
            self.assertIn('time_of_day_02', nc.dimensions, "No time_of_day_02 dimension found")
            # Should have 2 times of day approximately at 3 and 15 hours
            tod_var = None
            for var_name in nc.variables:
                if 'time_of_day' in var_name:
                    tod_var = nc.variables[var_name]
                    break
            
            self.assertIsNotNone(tod_var, "No time_of_day variable found")
            self.assertEqual(len(tod_var), 2, "Expected 2 time of day values")
            
            # Check that values are close to 3 and 15
            values = sorted(tod_var[:])
            self.assertAlmostEqual(values[0], 3.0, delta=0.5)
            self.assertAlmostEqual(values[1], 15.0, delta=0.5)
        finally:
            nc.close()
    
    def test_bin_average(self):
        """Test bin_average operation on daily file"""
        result = self.run_mars_files(['01336.atmos_daily.nc', '-ba'])
        
        # Check for successful execution
        self.assertEqual(result.returncode, 0, "Bin average command failed")
        
        # Check that output file was created
        output_file = self.check_file_exists('01336.atmos_daily_to_average.nc')
        
        # Verify that the output file has expected structure
        self.verify_netcdf_has_variable(output_file, 'time')
        
        # Verify that time dimension is smaller in output (binned) file
        nc_in = Dataset(os.path.join(self.test_dir, '01336.atmos_daily.nc'), 'r')
        nc_out = Dataset(output_file, 'r')
        try:
            self.assertLess(
                len(nc_out.dimensions['time']), 
                len(nc_in.dimensions['time']), 
                "Output time dimension should be smaller than input"
            )
        finally:
            nc_in.close()
            nc_out.close()
    
    def test_bin_diurn(self):
        """Test bin_diurn operation on daily file"""
        result = self.run_mars_files(['01336.atmos_daily.nc', '-bd'])
        
        # Check for successful execution
        self.assertEqual(result.returncode, 0, "Bin diurn command failed")
        
        # Check that output file was created
        output_file = self.check_file_exists('01336.atmos_daily_to_diurn.nc')
        
        # Verify that the output file has expected structure
        nc = Dataset(output_file, 'r')
        try:
            # Should have a time of day dimension
            found_tod = False
            for dim_name in nc.dimensions:
                if 'time_of_day' in dim_name:
                    found_tod = True
                    break
            
            self.assertTrue(found_tod, "No time_of_day dimension found in output file")
        finally:
            nc.close()
    
    def test_split_file_by_areo(self):
        """Test split operation on average file by Ls (areo)"""
        # First check what Ls values are available in the file
        nc = Dataset(os.path.join(self.test_dir, '01336.atmos_average.nc'), 'r')
        try:
            areo_values = nc.variables['areo'][:]
            ls_min = np.min(areo_values) % 360
            ls_max = np.max(areo_values) % 360
            
            # Choose a range within the available values
            ls_range = [ls_min, ls_max]
        finally:
            nc.close()
        
        # Run split command
        result = self.run_mars_files([
            '01336.atmos_average.nc', 
            '--split', 
            str(ls_range[0]), 
            str(ls_range[1]), 
            '-dim', 
            'areo'
        ])
        
        # Check for successful execution
        self.assertEqual(result.returncode, 0, "Split by areo command failed")
        
        # Check that output file was created - name will depend on Ls values
        ls_files = glob.glob(os.path.join(self.test_dir, '*_Ls*_*.nc'))
        self.assertTrue(len(ls_files) > 0, "No output files from split operation found")
    
    def test_split_file_by_lat(self):
        """Test split operation on average file by latitude"""
        result = self.run_mars_files([
            '01336.atmos_average.nc', 
            '--split', 
            '-45', 
            '45', 
            '-dim', 
            'lat'
        ])
        
        # Check for successful execution
        self.assertEqual(result.returncode, 0, "Split by latitude command failed")
        
        # Check that output file was created
        lat_files = glob.glob(os.path.join(self.test_dir, '*_lat_*_*.nc'))
        self.assertTrue(len(lat_files) > 0, "No output files from split by latitude operation found")
        
        # Verify the latitude range in the output file
        if lat_files:
            nc = Dataset(lat_files[0], 'r')
            try:
                lat_values = nc.variables['lat'][:]
                self.assertGreaterEqual(min(lat_values), -45)
                self.assertLessEqual(max(lat_values), 45)
            finally:
                nc.close()
    
    def test_temporal_high_pass_filter(self):
        """Test high-pass temporal filtering"""
        result = self.run_mars_files(['01336.atmos_daily.nc', '-hpt', '10'])
        
        # Check for successful execution
        self.assertEqual(result.returncode, 0, "High-pass temporal filter command failed")
        
        # Check that output file was created
        output_file = self.check_file_exists('01336.atmos_daily_hpt.nc')
        
        # Verify that the output file has expected structure
        nc = Dataset(output_file, 'r')
        try:
            # Should have a constant for filter cutoff
            found_constant = False
            for attr_name in nc.ncattrs():
                if 'sol_min' in attr_name:
                    found_constant = True
                    break
            
            self.assertTrue(found_constant, "No sol_min constant found in output file")
        finally:
            nc.close()
    
    def test_temporal_low_pass_filter(self):
        """Test low-pass temporal filtering"""
        result = self.run_mars_files(['01336.atmos_daily.nc', '-lpt', '0.75'])
        
        # Check for successful execution
        self.assertEqual(result.returncode, 0, "Low-pass temporal filter command failed")
        
        # Check that output file was created
        output_file = self.check_file_exists('01336.atmos_daily_lpt.nc')
    
    def test_temporal_band_pass_filter(self):
        """Test band-pass temporal filtering"""
        result = self.run_mars_files(['01336.atmos_daily.nc', '-bpt', '0.75', '10'])
        
        # Check for successful execution
        self.assertEqual(result.returncode, 0, "Band-pass temporal filter command failed")
        
        # Check that output file was created
        output_file = self.check_file_exists('01336.atmos_daily_bpt.nc')
    
    def test_temporal_filter_with_trend(self):
        """Test temporal filtering with trend added back"""
        result = self.run_mars_files(['01336.atmos_daily.nc', '-hpt', '10', '-add_trend'])
        
        # Check for successful execution
        self.assertEqual(result.returncode, 0, "Temporal filter with trend command failed")
        
        # Check that output file was created
        output_file = self.check_file_exists('01336.atmos_daily_hpt_trended.nc')
        
    def test_tide_decomposition(self):
        """Test tidal decomposition on diurn file"""
        result = self.run_mars_files(['01336.atmos_diurn.nc', '-tide', '2', '-incl', 'ps', 'temp'])
        
        # Check for successful execution
        self.assertEqual(result.returncode, 0, "Tide decomposition command failed")
        
        # Check that output file was created
        output_file = self.check_file_exists('01336.atmos_diurn_tide_decomp.nc')
        
        # Verify that the output file has amplitude and phase variables
        self.verify_netcdf_has_variable(output_file, 'ps_amp')
        self.verify_netcdf_has_variable(output_file, 'ps_phas')
        self.verify_netcdf_has_variable(output_file, 'temp_amp')
        self.verify_netcdf_has_variable(output_file, 'temp_phas')
    
    def test_tide_decomposition_with_normalize(self):
        """Test tidal decomposition with normalization"""
        result = self.run_mars_files(['01336.atmos_diurn.nc', '-tide', '2', '-incl', 'ps', '-norm'])
        
        # Check for successful execution
        self.assertEqual(result.returncode, 0, "Tide decomposition with normalization command failed")
        
        # Check that output file was created
        output_file = self.check_file_exists('01336.atmos_diurn_tide_decomp_norm.nc')
        
        # Verify output has normalized amplitude
        nc = Dataset(output_file, 'r')
        try:
            # Check if amplitude variable uses percentage units
            self.assertIn('ps_amp', nc.variables)
            amp_var = nc.variables['ps_amp']
            # Modified check: Either the units contain '%' or the variable has a 'normalized' attribute
            has_percent = False
            if 'units' in amp_var.ncattrs():
                units = getattr(amp_var, 'units', '')
                has_percent = '%' in units.lower()
            
            has_normalized_attr = False
            if 'normalized' in amp_var.ncattrs():
                has_normalized_attr = getattr(amp_var, 'normalized') in [True, 1, 'true', 'yes']
            
            self.assertTrue(has_percent or has_normalized_attr, 
                        "Normalized amplitude should have '%' in units or a 'normalized' attribute")
        finally:
            nc.close()
    
    def test_tide_decomposition_with_reconstruct(self):
        """Test tidal decomposition with reconstruction"""
        result = self.run_mars_files(['01336.atmos_diurn.nc', '-tide', '2', '-incl', 'ps', '-recon'])
        
        # Check for successful execution
        self.assertEqual(result.returncode, 0, "Tide decomposition with reconstruction command failed")
        
        # Check that output file was created
        output_file = self.check_file_exists('01336.atmos_diurn_tide_decomp_reconstruct.nc')
        
        # Verify output has reconstructed harmonics
        self.verify_netcdf_has_variable(output_file, 'ps_N1')
        self.verify_netcdf_has_variable(output_file, 'ps_N2')
    
    def test_regrid(self):
        """Test regridding operation"""
        result = self.run_mars_files([
            '01336.atmos_average_pstd.nc', 
            '-regrid', 
            '01336.atmos_average_pstd_c48.nc'
        ])
        
        # Check for successful execution
        self.assertEqual(result.returncode, 0, "Regrid command failed")
        
        # Check that output file was created
        output_file = self.check_file_exists('01336.atmos_average_regrid.nc')
        
        # Verify that the grid dimensions match the target file
        nc_target = Dataset(os.path.join(self.test_dir, '01336.atmos_average_pstd_c48.nc'), 'r')
        nc_output = Dataset(output_file, 'r')
        try:
            self.assertEqual(
                len(nc_output.dimensions['lat']), 
                len(nc_target.dimensions['lat']), 
                "Latitude dimension doesn't match target file"
            )
            self.assertEqual(
                len(nc_output.dimensions['lon']), 
                len(nc_target.dimensions['lon']), 
                "Longitude dimension doesn't match target file"
            )
        finally:
            nc_target.close()
            nc_output.close()
    
    def test_zonal_average(self):
        """Test zonal averaging operation"""
        result = self.run_mars_files(['01336.atmos_average.nc', '-zavg'])
        
        # Check for successful execution
        self.assertEqual(result.returncode, 0, "Zonal average command failed")
        
        # Check that output file was created
        output_file = self.check_file_exists('01336.atmos_average_zavg.nc')
        
        # Verify that the longitude dimension is 1 in the output file
        nc = Dataset(output_file, 'r')
        try:
            self.assertEqual(len(nc.dimensions['lon']), 1, "Longitude dimension should have size 1")
        finally:
            nc.close()
    
    def test_custom_extension(self):
        """Test using custom extension"""
        result = self.run_mars_files(['01336.atmos_average.nc', '-zavg', '-ext', '_custom'])
        
        # Check for successful execution
        self.assertEqual(result.returncode, 0, "Zonal average with custom extension command failed")
        
        # Check that output file was created with custom extension
        output_file = self.check_file_exists('01336.atmos_average_zavg_custom.nc')
    
    def test_include_vars(self):
        """Test including only specific variables"""
        result = self.run_mars_files(['01336.atmos_average.nc', '-zavg', '-incl', 'ps', 'temp'])
        
        # Check for successful execution
        self.assertEqual(result.returncode, 0, "Zonal average with include command failed")
        
        # Check that output file was created
        output_file = self.check_file_exists('01336.atmos_average_zavg.nc')
        
        # Verify that only included variables (plus dimensions) are in the output
        nc = Dataset(output_file, 'r')
        try:
            # Should have ps and temp
            self.assertIn('ps', nc.variables, "Variable ps not found in output")
            self.assertIn('temp', nc.variables, "Variable temp not found in output")
            
            # Should not have other variables that might be in the original file
            # This check depends on what's in your test files, adjust as needed
            all_vars = set(nc.variables.keys())
            expected_vars = {'ps', 'temp', 'time', 'lat', 'lon'}
            # Add any dimension variables
            for dim in nc.dimensions:
                expected_vars.add(dim)
            
            # Check if there are unexpected variables
            # Skip this test if we don't know what all variables should be
            # This is just an example of a possible check
            unexpected_vars = all_vars - expected_vars
            for var in unexpected_vars:
                # Skip dimension variables and coordinate variables
                if var in nc.dimensions or var in ['areo', 'scalar_axis'] or var.startswith('time_of_day'):
                    continue
                # Skip typical dimension bounds variables
                if var.endswith('_bnds') or var.endswith('_bounds'):
                    continue
                # Skip typical dimension coordinate variables
                if var in ['pfull', 'phalf', 'pstd', 'zstd', 'zagl', 'pk', 'bk']:
                    continue
                
                self.fail(f"Unexpected variable {var} found in output")
        finally:
            nc.close()


if __name__ == '__main__':
    unittest.main()