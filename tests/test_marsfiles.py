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
from base_test import BaseTestCase
import time

# Check if pyshtools is available
try:
    import pyshtools
    HAVE_PYSHTOOLS = True
except ImportError:
    HAVE_PYSHTOOLS = False
    
class TestMarsFiles(BaseTestCase):
    """Integration test suite for MarsFiles"""

    # Class attribute for storing modified files
    modified_files = {}
    
    @classmethod
    def setUpClass(cls):
        """Set up the test environment once for all tests"""
        # Create a temporary directory for the tests
        cls.test_dir = tempfile.mkdtemp(prefix='MarsFiles_test_')
        print(f"Created temporary test directory: {cls.test_dir}")

        # Project root directory
        cls.project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        print(f"Project root directory: {cls.project_root}")

        # Start timing for test file creation
        start_time = time.time()
        
        # Run the script to create test netCDF files
        cls.create_test_files()
        
        # Report how long file creation took
        elapsed = time.time() - start_time
        print(f"Test file creation completed in {elapsed:.2f} seconds")
        
        # Dictionary to keep track of modified files
        cls.modified_files = {}
        
        # Dictionary to track files created in each test
        cls.test_created_files = {}
        
        # Initialize modified_files dictionary with original files
        expected_files = [
            '01336.atmos_average.nc',
            '01336.atmos_average_pstd_c48.nc',
            '01336.atmos_daily.nc',
            '01336.atmos_diurn.nc',
            '01336.atmos_diurn_pstd.nc',
            '01336.fixed.nc'
        ]
        
        for filename in expected_files:
            cls.modified_files[filename] = os.path.join(cls.test_dir, filename)

    @classmethod
    def create_test_files(cls):
        """Create test netCDF files using create_ames_gcm_files.py"""
        # Get path to create_ames_gcm_files.py script
        create_files_script = os.path.join(cls.project_root, "tests", "create_ames_gcm_files.py")

        # Execute the script to create test files - Important: pass the test_dir as argument
        cmd = [sys.executable, create_files_script, cls.test_dir]

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                cwd=cls.test_dir,  # Run in the test directory to ensure files are created there
                timeout=300  # Add a timeout to prevent hanging
            )

            if result.returncode != 0:
                raise Exception(f"Failed to create test files: {result.stderr}")

        except subprocess.TimeoutExpired:
            raise Exception("Timeout creating test files - process took too long and was terminated")
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
                raise Exception(f"Test file {filename} was not created in {cls.test_dir}")
            # Print file size to help with debugging
            file_size = os.path.getsize(filepath) / (1024 * 1024)  # Size in MB
            print(f"Created {filename}: {file_size:.2f} MB")


    def setUp(self):
        """Change to temporary directory before each test"""
        os.chdir(self.test_dir)
        
        # Store the current test method name
        self.current_test = self.id().split('.')[-1]
        
        # Print test start time for debugging
        print(f"\nStarting test: {self.current_test} at {time.strftime('%H:%M:%S')}")
        self.start_time = time.time()
        
        # Initialize an empty list to track files created by this test
        self.__class__.test_created_files[self.current_test] = []
        
        # Get a snapshot of files in the directory before the test runs
        self.files_before_test = set(os.listdir(self.test_dir))

    def tearDown(self):
        """Clean up after each test"""
        # Calculate and print test duration
        elapsed = time.time() - self.start_time
        print(f"Test {self.current_test} completed in {elapsed:.2f} seconds")
        
        # Get files that exist after the test
        files_after_test = set(os.listdir(self.test_dir))
        
        # Find new files created during this test
        new_files = files_after_test - self.files_before_test
        
        # Store these new files in our tracking dictionary
        for filename in new_files:
            file_path = os.path.join(self.test_dir, filename)
            self.__class__.test_created_files[self.current_test].append(file_path)
            
            # Also track in modified_files if it's a netCDF file we want to keep
            if filename.endswith('.nc') and filename not in self.modified_files:
                self.modified_files[filename] = file_path
        
        # Get the list of files to clean up (files created by this test that aren't in modified_files)
        files_to_clean = []
        for file_path in self.__class__.test_created_files[self.current_test]:
            filename = os.path.basename(file_path)
            # If this is a permanent file we want to keep, skip it
            if file_path in self.modified_files.values():
                continue
            # Clean up temporary files
            files_to_clean.append(file_path)
        
        # Remove temporary files created by this test
        for file_path in files_to_clean:
            try:
                if os.path.exists(file_path):
                    os.remove(file_path)
                    print(f"Cleaned up: {os.path.basename(file_path)}")
            except Exception as e:
                print(f"Warning: Could not remove file {file_path}: {e}")

        # Return to test_dir
        os.chdir(self.test_dir)

    @classmethod
    def tearDownClass(cls):
        """Clean up the test environment"""
        try:
            shutil.rmtree(cls.test_dir, ignore_errors=True)
        except Exception as e:
            print(f"Warning: Could not remove test directory {cls.test_dir}: {e}")

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
                # Check if we have a modified version of this file
                base_filename = os.path.basename(arg)
                if base_filename in self.modified_files:
                    abs_args.append(self.modified_files[base_filename])
                else:
                    abs_args.append(os.path.join(self.test_dir, arg))
            else:
                abs_args.append(arg)

        # Construct the full command to run MarsFiles
        cmd = [sys.executable, os.path.join(self.project_root, "bin", "MarsFiles.py")] + abs_args

        # Get a snapshot of files before running the command
        files_before = set(os.listdir(self.test_dir))

        # Run the command with a timeout
        start_time = time.time()
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                cwd=self.test_dir,  # Run in the test directory
                env=dict(os.environ, PWD=self.test_dir),  # Ensure current working directory is set
                timeout=300  # Set a reasonable timeout (5 minutes) per subprocess
            )
            elapsed = time.time() - start_time
            print(f"Subprocess completed in {elapsed:.2f} seconds")

            # If command succeeded, find any new files that were created
            if result.returncode == 0:
                files_after = set(os.listdir(self.test_dir))
                new_files = files_after - files_before
                
                # Track these new files
                for filename in new_files:
                    file_path = os.path.join(self.test_dir, filename)
                    # Add to test tracking
                    self.__class__.test_created_files[self.current_test].append(file_path)
                    
                    # Also update the modified_files dictionary for output files we need to track
                    if filename.endswith('.nc'):
                        # Track the file in our modified_files dictionary
                        self.modified_files[filename] = file_path

            return result
        except subprocess.TimeoutExpired:
            print(f"ERROR: Subprocess timed out after {time.time() - start_time:.2f} seconds")
            self.fail("Subprocess timed out")
        except Exception as e:
            self.fail(f"Failed to run MarsFiles: {e}")

    def check_file_exists(self, filename):
        """
        Check if a file exists and is not empty

        :param filename: Filename to check
        """
        # First check if we have this file in our modified_files dictionary
        if filename in self.modified_files:
            filepath = self.modified_files[filename]
        else:
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
        
        print("Help message displayed successfully")

    def test_time_shift(self):
        """Test time_shift operation on diurn file"""
        result = self.run_mars_files(['01336.atmos_diurn_pstd.nc', '-t', '--debug'])

        # Check for successful execution
        self.assertEqual(result.returncode, 0, "Time shift command failed")

        # Check that output file was created
        output_file = self.check_file_exists('01336.atmos_diurn_pstd_T.nc')

        # Verify that the output file has expected structure
        self.verify_netcdf_has_variable(output_file, 'time')
        self.verify_netcdf_has_variable(output_file, 'time_of_day_24')  # Default is 24 time of day bins
        
        print("Time shift operation succeeded")

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
        
        print("Time shift with specific times succeeded")

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
        
        print("Bin average operation succeeded")

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
        
        print("Bin diurn operation succeeded")

    def test_split_file_by_areo(self):
        """Test split operation on average file by Ls (areo)"""
        # First check what Ls values are available in the file
        nc = Dataset(os.path.join(self.test_dir, '01336.atmos_average.nc'), 'r')
        try:
            areo_values = nc.variables['areo'][:]
            ls_min = np.min(areo_values) % 360
            ls_max = np.max(areo_values) % 360

            # Choose a range within the available values
            ls_range = [ls_min + 50 , ls_max - 50]
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
        
        print(f"Split file by areo succeeded (Ls range: {ls_range[0]}-{ls_range[1]})")

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
        
        print("Split file by latitude succeeded")
        
    def test_split_file_by_lon(self):
        """Test split operation on average file by longitude"""
        result = self.run_mars_files([
            '01336.atmos_average.nc',
            '--split',
            '10',
            '20',
            '-dim',
            'lon'
        ])

        # Check for successful execution
        self.assertEqual(result.returncode, 0, "Split by longitude command failed")

        # Check that output file was created
        lon_files = glob.glob(os.path.join(self.test_dir, '*_lon_*_*.nc'))
        self.assertTrue(len(lon_files) > 0, "No output files from split by longitude operation found")

        # Verify the longitude range in the output file
        if lon_files:
            nc = Dataset(lon_files[0], 'r')
            try:
                lon_values = nc.variables['lon'][:]
                self.assertGreaterEqual(min(lon_values), 10)
                self.assertLessEqual(max(lon_values), 20)
            finally:
                nc.close()
        
        print("Split file by longitude succeeded")

    def test_temporal_filters(self):
        """Test all temporal filtering operations"""
        # High-pass filter
        result = self.run_mars_files(['01336.atmos_daily.nc', '-hpt', '10', '-incl', 'temp'])
        self.assertEqual(result.returncode, 0, "High-pass temporal filter command failed")
        high_pass_file = self.check_file_exists('01336.atmos_daily_hpt.nc')
        self.verify_netcdf_has_variable(high_pass_file, 'temp')
        print("High-pass temporal filter succeeded")
        
        # Low-pass filter
        result = self.run_mars_files(['01336.atmos_daily.nc', '-lpt', '0.75', '-incl', 'temp'])
        self.assertEqual(result.returncode, 0, "Low-pass temporal filter command failed")
        low_pass_file = self.check_file_exists('01336.atmos_daily_lpt.nc')
        self.verify_netcdf_has_variable(low_pass_file, 'temp')
        print("Low-pass temporal filter succeeded")
        
        # Band-pass filter
        result = self.run_mars_files(['01336.atmos_daily.nc', '-bpt', '0.75', '10', '-add_trend', '-incl', 'temp'])
        self.assertEqual(result.returncode, 0, "Band-pass temporal filter with trend command failed")
        band_pass_file = self.check_file_exists('01336.atmos_daily_bpt_trended.nc')
        self.verify_netcdf_has_variable(band_pass_file, 'temp')
        print("Band-pass temporal filter succeeded")
    
    def test_spatial_filters(self):
        """Test all spatial filtering operations"""
        if not HAVE_PYSHTOOLS:
            self.skipTest("pyshtools is not available in this version of CAP."
                          "Install it to run this test.")
            
        # High-pass filter
        result = self.run_mars_files(['01336.atmos_daily.nc', '-hps', '10', '-incl', 'temp'])
        self.assertEqual(result.returncode, 0, "High-pass spatial filter command failed")
        high_pass_file = self.check_file_exists('01336.atmos_daily_hps.nc')
        self.verify_netcdf_has_variable(high_pass_file, 'temp')
        print("High-pass spatial filter succeeded")
        
        # Low-pass filter
        result = self.run_mars_files(['01336.atmos_daily.nc', '-lps', '20', '-incl', 'temp'])
        self.assertEqual(result.returncode, 0, "Low-pass spatial filter command failed")
        low_pass_file = self.check_file_exists('01336.atmos_daily_lps.nc')
        self.verify_netcdf_has_variable(low_pass_file, 'temp')
        print("Low-pass spatial filter succeeded")
        
        # Band-pass filter
        result = self.run_mars_files(['01336.atmos_daily.nc', '-bps', '10', '20', '-incl', 'temp'])
        self.assertEqual(result.returncode, 0, "Band-pass spatial filter command failed")
        band_pass_file = self.check_file_exists('01336.atmos_daily_bps.nc')
        self.verify_netcdf_has_variable(band_pass_file, 'temp')
        print("Band-pass spatial filter succeeded")

    def test_tide_decomposition(self):
        """Test tidal decomposition on diurn file"""
        if not HAVE_PYSHTOOLS:
            self.skipTest("pyshtools is not available in this version of CAP."
                          "Install it to run this test.")
        
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
        
        print("Tide decomposition succeeded")

    def test_tide_decomposition_with_normalize(self):
        """Test tidal decomposition with normalization"""
        if not HAVE_PYSHTOOLS:
            self.skipTest("pyshtools is not available in this version of CAP."
                          "Install it to run this test.")
            
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
            
        print("Tide decomposition with normalization succeeded")

    def test_tide_decomposition_with_reconstruct(self):
        """Test tidal decomposition with reconstruction"""
        if not HAVE_PYSHTOOLS:
            self.skipTest("pyshtools is not available in this version of CAP."
                          "Install it to run this test.")
            
        result = self.run_mars_files(['01336.atmos_diurn.nc', '-tide', '2', '-incl', 'ps', '-recon'])

        # Check for successful execution
        self.assertEqual(result.returncode, 0, "Tide decomposition with reconstruction and include commands failed")

        # Check that output file was created
        output_file = self.check_file_exists('01336.atmos_diurn_tide_decomp_reconstruct.nc')

        # Verify output has reconstructed harmonics
        self.verify_netcdf_has_variable(output_file, 'ps_N1')
        self.verify_netcdf_has_variable(output_file, 'ps_N2')
        
        print("Tide decomposition with reconstruction succeeded")

        # Verify that only included variables (plus dimensions) are in the output
        nc = Dataset(output_file, 'r')
        try:
            # Should have ps and temp
            self.assertIn('ps_N1', nc.variables, "Variable ps_N1 not found in output")

            # Should not have other variables that might be in the original file
            # This check depends on what's in the test files, adjust as needed
            all_vars = set(nc.variables.keys())
            expected_vars = {'ps_N1', 'ps_N2', 'time', 'lat', 'lon'}
            # Add any dimension variables
            for dim in nc.dimensions:
                expected_vars.add(dim)

            # Check if there are unexpected variables
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
        
        print("Include argument succeeded")
            
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
        output_file = self.check_file_exists('01336.atmos_average_pstd_regrid.nc')

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
        
        print("Regrid operation succeeded")

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
        
        print("Zonal average operation succeeded")

    def test_custom_extension(self):
        """Test using custom extension"""
        result = self.run_mars_files(['01336.atmos_average.nc', '-zavg', '-ext', 'custom'])

        # Check for successful execution
        self.assertEqual(result.returncode, 0, "Zonal average with custom extension command failed")

        # Check that output file was created with custom extension
        output_file = self.check_file_exists('01336.atmos_average_zavg_custom.nc')
        
        print("Custom extension operation succeeded")

    def test_zzz_output_cleanup_stats(self):
        """
        This test runs last (due to alphabetical sorting of 'zzz') and outputs statistics
        about files created and cleaned during testing.
        """
        if not hasattr(self.__class__, 'test_created_files'):
            self.skipTest("No file tracking information available")
            
        # Calculate total files created
        total_files = sum(len(files) for files in self.__class__.test_created_files.values())
        
        # Calculate files per test
        files_per_test = {test: len(files) for test, files in self.__class__.test_created_files.items()}
        
        # Find the test that created the most files
        max_files_test = max(files_per_test.items(), key=lambda x: x[1]) if files_per_test else (None, 0)
        
        # Output statistics
        print("\n===== File Cleanup Statistics =====")
        print(f"Total files created during testing: {total_files}")
        print(f"Average files per test: {total_files / len(self.__class__.test_created_files) if self.__class__.test_created_files else 0:.1f}")
        print(f"Test creating most files: {max_files_test[0]} ({max_files_test[1]} files)")
        print("==================================\n")
        
        # Test passes if we reach this point
        print("Selective cleanup system is working")
        
        # This isn't really a test, but we'll assert True to make it pass
        self.assertTrue(True)

if __name__ == '__main__':
    unittest.main()