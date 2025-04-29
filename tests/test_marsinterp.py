#!/usr/bin/env python3
"""
Integration tests for MarsInterp.py

These tests verify the functionality of MarsInterp for interpolating netCDF files
to various pressure and altitude coordinates.
"""

import os
import sys
import unittest
import shutil
import platform
import subprocess
import tempfile
import glob
import re
import numpy as np
from netCDF4 import Dataset

class TestMarsInterp(unittest.TestCase):
    """Integration test suite for MarsInterp"""
    
    @classmethod
    def setUpClass(cls):
        """Set up the test environment"""
        # Create a temporary directory for the tests
        cls.test_dir = tempfile.mkdtemp(prefix='MarsInterp_test_')
        print(f"Created temporary test directory: {cls.test_dir}")
        
        # Project root directory
        cls.project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        print(f"Project root directory: {cls.project_root}")
        
        # Run script to create test netCDF files (assumed to exist)
        cls.create_test_files()
    
    @classmethod
    def create_test_files(cls):
        """Create test netCDF files using create_ames_gcm_files.py"""
        # Get path to create_ames_gcm_files.py script
        create_files_script = os.path.join(cls.project_root, "tests", "create_ames_gcm_files.py")
        
        # Execute the script to create test files
        cmd = [sys.executable, create_files_script, cls.test_dir, 'short']
        
        print(f"Creating test files with command: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                cwd=cls.test_dir
            )
            
            print(f"File creation STDOUT: {result.stdout}")
            print(f"File creation STDERR: {result.stderr}")
            
            if result.returncode != 0:
                raise Exception(f"Failed to create test files: {result.stderr}")
            
        except Exception as e:
            raise Exception(f"Error running create_ames_gcm_files.py: {e}")
        
        # List files in the temp directory to debug
        print(f"Files in test directory after creation: {os.listdir(cls.test_dir)}")
            
        # Verify files were created
        expected_files = [
            '01336.atmos_average.nc',
            '01336.atmos_daily.nc',
            '01336.atmos_diurn.nc',
            '01336.fixed.nc'
        ]
        
        for filename in expected_files:
            filepath = os.path.join(cls.test_dir, filename)
            if not os.path.exists(filepath):
                raise Exception(f"Test file {filename} was not created in {cls.test_dir}")
            else:
                print(f"Confirmed test file exists: {filepath}")
    
    def setUp(self):
        """Change to temporary directory before each test"""
        os.chdir(self.test_dir)
        print(f"Changed to test directory: {os.getcwd()}")
    
    def tearDown(self):
        """Clean up after each test"""
        # Clean up any generated output files after each test but keep input files
        output_patterns = [
            '*_pstd*.nc',
            '*_zstd*.nc',
            '*_zagl*.nc',
            '*.txt'
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
            # List files in temp directory before deleting to debug
            print(f"Files in test directory before cleanup: {os.listdir(cls.test_dir)}")
            shutil.rmtree(cls.test_dir, ignore_errors=True)
            print(f"Removed test directory: {cls.test_dir}")
        except Exception as e:
            print(f"Warning: Could not remove test directory {cls.test_dir}: {e}")
    
    def run_mars_interp(self, args, expected_success=True):
        """
        Run MarsInterp using subprocess
        
        :param args: List of arguments to pass to MarsInterp
        :param expected_success: Whether the command is expected to succeed
        :return: subprocess result object
        """
        # Construct the full command to run MarsInterp
        cmd = [sys.executable, os.path.join(self.project_root, "bin", "MarsInterp.py")] + args
        
        # Print debugging info
        print(f"Running command: {' '.join(cmd)}")
        print(f"Working directory: {self.test_dir}")
        print(f"File exists check: {os.path.exists(os.path.join(self.project_root, 'bin', 'MarsInterp.py'))}")
        
        # Run the command
        try:
            result = subprocess.run(
                cmd, 
                capture_output=True, 
                text=True, 
                cwd=self.test_dir,
                env=dict(os.environ, PWD=self.test_dir)
            )
            
            # Print both stdout and stderr to help debug
            print(f"STDOUT: {result.stdout}")
            print(f"STDERR: {result.stderr}")
            
            if expected_success:
                self.assertEqual(result.returncode, 0, f"MarsInterp command failed with code {result.returncode}")
            else:
                self.assertNotEqual(result.returncode, 0, "MarsInterp command succeeded but was expected to fail")
            
            return result
        except Exception as e:
            self.fail(f"Failed to run MarsInterp: {e}")
    
    def check_file_exists(self, filename):
        """
        Check if a file exists and is not empty
        
        :param filename: Filename to check
        """
        filepath = os.path.join(self.test_dir, filename)
        self.assertTrue(os.path.exists(filepath), f"File {filename} does not exist")
        self.assertGreater(os.path.getsize(filepath), 0, f"File {filename} is empty")
        return filepath
    
    def check_netcdf_structure(self, filepath, expected_dimension):
        """
        Check if a netCDF file has the expected structure
        
        :param filepath: Path to the netCDF file
        :param expected_dimension: Expected vertical dimension name
        """
        try:
            nc = Dataset(filepath, 'r')
            
            # Check if the expected dimension exists
            self.assertIn(expected_dimension, nc.dimensions, 
                          f"Expected dimension {expected_dimension} not found in {filepath}")
            
            # Check if some typical 3D/4D variables are present
            # This depends on your specific setup; adjust as needed
            has_typical_vars = False
            for var_name in nc.variables:
                var = nc.variables[var_name]
                dims = var.dimensions
                if expected_dimension in dims and len(dims) >= 3:
                    has_typical_vars = True
                    break
            
            self.assertTrue(has_typical_vars, 
                           f"No variables with dimension {expected_dimension} found in {filepath}")
            
            nc.close()
        except Exception as e:
            self.fail(f"Error checking netCDF structure: {e}")
    
    def test_help_message(self):
        """Test that help message can be displayed"""
        result = self.run_mars_interp(['-h'])
        
        # Check for typical help message components
        help_checks = [
            'usage:',
            'input_file',
            '--interp_type',
            '--vertical_grid',
            '--include',
            '--extension',
            '--print_grid'
        ]
        
        for check in help_checks:
            self.assertTrue(any(check in line for line in result.stdout.split('\n')), 
                          f"Help message missing '{check}'")
    
    def test_print_grid(self):
        """Test printing a vertical grid without interpolation"""
        result = self.run_mars_interp(['01336.atmos_average.nc', '-t', 'pstd', '-print'])
        
        # Check that numeric values were printed
        # The output should contain floating point numbers
        has_numbers = bool(re.search(r'\d+(\.\d+)?', result.stdout))
        self.assertTrue(has_numbers, "No grid values found in output")
    
    def test_interpolate_to_pstd(self):
        """Test basic pressure interpolation (pstd)"""
        result = self.run_mars_interp(['01336.atmos_average.nc'])
        
        # Check for successful execution based on typical output
        self.assertIn("Completed in", result.stdout, "Missing completion message")
        
        # Check that the output file was created
        output_file = self.check_file_exists("01336.atmos_average_pstd.nc")
        
        # Check that the file has the expected structure
        self.check_netcdf_structure(output_file, "pstd")
    
    def test_interpolate_to_zstd(self):
        """Test interpolation to standard altitude (zstd)"""
        result = self.run_mars_interp(['01336.atmos_average.nc', '-t', 'zstd'])
        
        # Check for successful execution
        self.assertIn("Completed in", result.stdout, "Missing completion message")
        
        # Check that the output file was created
        output_file = self.check_file_exists("01336.atmos_average_zstd.nc")
        
        # Check that the file has the expected structure
        self.check_netcdf_structure(output_file, "zstd")
    
    def test_interpolate_to_zagl(self):
        """Test interpolation to altitude above ground level (zagl)"""
        result = self.run_mars_interp(['01336.atmos_average.nc', '-t', 'zagl'])
        
        # Check for successful execution
        self.assertIn("Completed in", result.stdout, "Missing completion message")
        
        # Check that the output file was created
        output_file = self.check_file_exists("01336.atmos_average_zagl.nc")
        
        # Check that the file has the expected structure
        self.check_netcdf_structure(output_file, "zagl")
    
    def test_custom_vertical_grid(self):
        """Test interpolation to a custom vertical grid"""
        # Note: This assumes a custom grid named 'pstd_default' exists in the profile
        result = self.run_mars_interp(['01336.atmos_average.nc', '-t', 'pstd', '-v', 'pstd_default'])
        
        # Check for successful execution
        self.assertIn("Completed in", result.stdout, "Missing completion message")
        
        # Check that the output file was created
        output_file = self.check_file_exists("01336.atmos_average_pstd.nc")
        
        # Check that the file has the expected structure
        self.check_netcdf_structure(output_file, "pstd")
    
    def test_include_specific_variables(self):
        """Test including only specific variables in interpolation"""
        # First we need to know valid variable names from the file
        nc = Dataset(os.path.join(self.test_dir, '01336.atmos_average.nc'), 'r')
        try:
            # Get variables that are likely 3D/4D (have 'pfull' dimension)
            var_names = []
            for var in nc.variables:
                if 'pfull' in nc.variables[var].dimensions:
                    var_names.append(var)
                    if len(var_names) >= 2:  # Get at least 2 variables
                        break
            
            if len(var_names) < 2:
                self.skipTest("Not enough variables with pfull dimension found")
                
        finally:
            nc.close()
        
        # Run interpolation with only these variables
        result = self.run_mars_interp(['01336.atmos_average.nc', '-incl'] + var_names)
        
        # Check for successful execution
        self.assertIn("Completed in", result.stdout, "Missing completion message")
        
        # Check that the output file was created
        output_file = self.check_file_exists("01336.atmos_average_pstd.nc")
        
        # Open the output file and confirm only specified variables are there
        nc = Dataset(output_file, 'r')

        # Remove 'pfull' from the list to check
        if 'pfull' in var_names:
            var_names.remove('pfull')

        try:
            # Check that each specified variable is present
            for var_name in var_names:
                self.assertIn(var_name, nc.variables, f"Variable {var_name} missing from output")
            
            # Count how many variables with pstd dimension are there
            pstd_var_count = 0
            for var in nc.variables:
                if 'pstd' in nc.variables[var].dimensions and var not in ['pstd', 'time', 'lon', 'lat']:
                    pstd_var_count += 1
            
            # Should only have our specified variables (plus dimensions/1D variables)
            self.assertEqual(pstd_var_count, len(var_names), 
                           f"Expected {len(var_names)} variables with pstd dimension, found {pstd_var_count}")
        finally:
            nc.close()
    
    def test_custom_extension(self):
        """Test creating output with a custom extension"""
        extension = "custom_test"
        result = self.run_mars_interp(['01336.atmos_average.nc', '-ext', extension])
        
        # Check for successful execution
        self.assertIn("Completed in", result.stdout, "Missing completion message")
        
        # Check that the output file was created with the custom extension
        output_file = self.check_file_exists(f"01336.atmos_average_pstd_{extension}.nc")
        
        # Check that the file has the expected structure
        self.check_netcdf_structure(output_file, "pstd")
    
    def test_multiple_files(self):
        """Test interpolating multiple files at once"""
        input_files = ['01336.atmos_average.nc', '01336.atmos_daily.nc']
        result = self.run_mars_interp(input_files)
        
        # Check for successful execution
        self.assertIn("Completed in", result.stdout, "Missing completion message")
        
        # Check that both output files were created
        for input_file in input_files:
            output_file = input_file.replace('.nc', '_pstd.nc')
            self.check_file_exists(output_file)
            self.check_netcdf_structure(os.path.join(self.test_dir, output_file), "pstd")
    
    def test_debug_mode(self):
        """Test running in debug mode with an error"""
        # Create an invalid netCDF file
        with open(os.path.join(self.test_dir, "Invalid.nc"), "w") as f:
            f.write("This is not a netCDF file")
        
        # This should fail with detailed error in debug mode
        result = self.run_mars_interp(['Invalid.nc', '--debug'], expected_success=False)
        
        # Look for traceback in output
        self.assertTrue('Traceback' in result.stdout or 'Traceback' in result.stderr,
                      "No traceback found in debug output")
    
    def test_invalid_interpolation_type(self):
        """Test error handling with invalid interpolation type"""
        result = self.run_mars_interp(['01336.atmos_average.nc', '-t', 'invalid_type'], expected_success=True)
        
        # Check for error message about unsupported interpolation type
        error_indicators = [
            'not supported',
            'use `pstd`, `zstd` or `zagl`'
        ]
        
        # At least one of these should appear in stderr or stdout
        all_output = result.stdout + result.stderr
        self.assertTrue(any(indicator in all_output for indicator in error_indicators),
                      "No error message about invalid interpolation type")
    
    def test_invalid_netcdf_file(self):
        """Test error handling with invalid netCDF file"""
        # Create an invalid netCDF file
        with open(os.path.join(self.test_dir, "Invalid.nc"), "w") as f:
            f.write("This is not a netCDF file")
        
        # This should fail
        result = self.run_mars_interp(['Invalid.nc'], expected_success=False)
        
        # Check for error message
        self.assertTrue(
            any(word in result.stdout for word in ['ERROR', 'error', 'failed', 'Failed', 'Unknown']) or
            any(word in result.stderr for word in ['ERROR', 'error', 'failed', 'Failed', 'Unknown']),
            "No error message for invalid netCDF file"
        )
    
    def test_diurnal_file_interpolation(self):
        """Test interpolation of a diurnal file"""
        result = self.run_mars_interp(['01336.atmos_diurn.nc'])
        
        # Check for successful execution
        self.assertIn("Completed in", result.stdout, "Missing completion message")
        
        # Check that the output file was created
        output_file = self.check_file_exists("01336.atmos_diurn_pstd.nc")
        
        # Check that the file has the expected structure
        self.check_netcdf_structure(output_file, "pstd")
        
        # For diurnal files, we should also check that the time_of_day dimension is preserved
        nc = Dataset(output_file, 'r')
        try:
            # Find the time_of_day dimension (name might vary)
            tod_dim = None
            for dim in nc.dimensions:
                if 'time_of_day' in dim:
                    tod_dim = dim
                    break
            
            self.assertIsNotNone(tod_dim, "No time_of_day dimension found in diurnal output file")
        finally:
            nc.close()


if __name__ == '__main__':
    unittest.main()
