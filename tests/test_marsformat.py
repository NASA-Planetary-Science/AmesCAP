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
        # Create a temporary directory for test files instead of a fixed path
        cls.test_dir = tempfile.mkdtemp(prefix='MarsFormat_test_')
        
        # Project root directory
        cls.project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    
    def setUp(self):
        """Create test netCDF files using create_gcm_files.py"""
        # Change to test directory
        os.chdir(self.test_dir)
        
        # Define file paths
        self.emars_file = os.path.join(self.test_dir, "emars_test.nc")
        self.openmars_file = os.path.join(self.test_dir, "openmars_test.nc")
        self.pcm_file = os.path.join(self.test_dir, "pcm_test.nc")
        self.marswrf_file = os.path.join(self.test_dir, "marswrf_test.nc")

        # Get path to create_gcm_files.py script
        create_files_script = os.path.join(self.project_root, "tests", "create_gcm_files.py")
        
        # Execute the script to create test files
        result = subprocess.run(
            [sys.executable, create_files_script],
            capture_output=True,
            text=True,
            cwd=self.test_dir
        )
        
        # Print output for debugging
        print(f"File creation output: {result.stdout}")
        
        # Check files were created
        for test_file in [self.emars_file, self.openmars_file, self.pcm_file, self.marswrf_file]:
            if not os.path.exists(test_file):
                print(f"Warning: Test file {test_file} was not created!")
                if result.stderr:
                    print(f"Error output: {result.stderr}")
    
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
        dummy_files = [
            self.emars_file, self.openmars_file, self.pcm_file, self.marswrf_file
        ]
        for file in dummy_files:
            if os.path.exists(file):
                os.remove(file)
    
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
        # Convert any relative file paths to absolute paths
        abs_args = []
        for arg in args:
            if arg.endswith('.nc'):
                abs_args.append(os.path.join(self.test_dir, arg))
            else:
                abs_args.append(arg)
        
        # Construct the full command to run MarsFormat
        cmd = [sys.executable, os.path.join(self.project_root, "bin", "MarsFormat.py")] + abs_args
        
        # Print debugging info
        print(f"Running command: {' '.join(cmd)}")
        print(f"Working directory: {self.test_dir}")
        print(f"File exists check: {os.path.exists(os.path.join(self.project_root, 'bin', 'MarsFormat.py'))}")
        
        # Run the command
        result = subprocess.run(
            cmd, 
            capture_output=True, 
            text=True, 
            cwd=self.test_dir,  # Run in the test directory
            env=dict(os.environ, PWD=self.test_dir)  # Ensure current working directory is set
        )
        
        # Print both stdout and stderr to help debug
        print(f"STDOUT: {result.stdout}")
        print(f"STDERR: {result.stderr}")
        
        return result
    
    # Test methods remain unchanged
    def test_daily_average(self):
        """Test the daily average functionality of MarsFormat.py for EMARS file."""
        # Run MarsFormat with the correct arguments
        result = self.run_mars_format([os.path.basename(self.emars_file), "-gcm", "emars"])
        
        # Check that the command executed successfully
        self.assertEqual(result.returncode, 0, f"MarsFormat.py failed: {result.stderr}")
        
        # Check that the output file was created
        output_file = os.path.join(self.test_dir, "emars_test_daily.nc")
        
        print(f"Checking for file at: {output_file}")
        print(f"Current directory: {os.getcwd()}")
        print(f"File exists: {os.path.exists(output_file)}")
        
        self.assertTrue(os.path.exists(output_file), f"Output file {output_file} was not created.")
            
        # Open the output file and check that it contains the expected variables
        dataset = nc.Dataset(output_file, 'r')
        
        # Debug - print all variable names to examine what's actually in the file
        print(f"Variables in {output_file}: {list(dataset.variables.keys())}")
        
        # Check that key variables are present
        self.assertIn('temp', dataset.variables, "Temperature variable not found in output file.")
        self.assertIn('ps', dataset.variables, "Surface pressure variable not found in output file.")
        
        # Check that coordinate variables are present
        self.assertIn('lat', dataset.variables, "Latitude variable not found in output file.")
        self.assertIn('lon', dataset.variables, "Longitude variable not found in output file.")
        self.assertIn('time', dataset.variables, "Time variable not found in output file.")
        self.assertIn('pfull', dataset.variables, "Pressure levels variable not found in output file.")
        
        # Check that variables have correct units/attributes
        self.assertEqual(dataset.variables['temp'].units, 'K', "Temperature units incorrect.")
        self.assertEqual(dataset.variables['ps'].units, 'pascal', "Surface pressure units incorrect.")
        
        dataset.close()

    # The rest of the test methods remain unchanged
    # ...
    
    # Include all the remaining test methods from the original file
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
        
        # Check for key variables we expect to see in diurnal average
        self.assertIn('temp', dataset.variables, "Temperature variable not found in output file.")
        self.assertIn('ps', dataset.variables, "Surface pressure variable not found in output file.")
        
        # Check that coordinate variables are present
        self.assertIn('lat', dataset.variables, "Latitude variable not found in output file.")
        self.assertIn('lon', dataset.variables, "Longitude variable not found in output file.")
        self.assertIn('time', dataset.variables, "Time variable not found in output file.")
        
        # Check that time_of_day dimension exists
        time_of_day_found = False
        for dim_name in dataset.dimensions:
            if 'time_of_day' in dim_name:
                time_of_day_found = True
                break
        self.assertTrue(time_of_day_found, "No time_of_day dimension found in output file.")
        
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
        
        # Debug - print all variable names 
        print(f"Variables in {output_file}: {list(dataset.variables.keys())}")
        
        # Check for core variables and coordinates
        self.assertIn('temp', dataset.variables, "Temperature variable not found in output file.")
        self.assertIn('ps', dataset.variables, "Surface pressure variable not found in output file.")
        self.assertIn('lat', dataset.variables, "Latitude variable not found in output file.")
        self.assertIn('lon', dataset.variables, "Longitude variable not found in output file.")
        self.assertIn('time', dataset.variables, "Time variable not found in output file.")
        self.assertIn('pfull', dataset.variables, "Pressure levels variable not found in output file.")
        
        # Check that hybrid coordinate variables exist
        self.assertIn('ak', dataset.variables, "Hybrid coordinate 'ak' not found in output file.")
        self.assertIn('bk', dataset.variables, "Hybrid coordinate 'bk' not found in output file.")
        
        # Check that coordinate values have been properly converted
        # OpenMARS typically has longitudes from -180 to 180, but should be 0 to 360 after conversion
        self.assertGreaterEqual(dataset.variables['lon'][0], 0, "Longitude values not converted to 0-360 range")
        
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
        
        # Print variables for debugging
        print(f"Variables in {output_file}: {list(dataset.variables.keys())}")
        
        # Check for key variables
        self.assertIn('temp', dataset.variables, "Temperature variable not found in output file.")
        self.assertIn('ps', dataset.variables, "Surface pressure variable not found in output file.")
        
        # Check for coordinate variables
        self.assertIn('latitude', dataset.variables, "Latitude variable not found in output file.")
        self.assertIn('longitude', dataset.variables, "Longitude variable not found in output file.")
        self.assertIn('time', dataset.variables, "Time variable not found in output file.")
        self.assertIn('pfull', dataset.variables, "Pfull variable not found in output file.")
        
        # Check that longitude values have been properly converted (if needed)
        if min(dataset.variables['longitude'][:]) < 0:
            self.fail("Longitude values not converted to 0-360 range")
            
        dataset.close()

    def test_marswrf_format(self):
        """Test MarsFormat.py with MarsWRF file."""
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
        
        # Check output contains expected variables
        dataset = nc.Dataset(output_file, 'r')
        print(f"Variables in {output_file}: {list(dataset.variables.keys())}")
        
        # Check for key variables
        # Temperature could be named 'temp' or 'T' depending on settings
        temp_var_found = False
        for var_name in dataset.variables:
            if var_name.lower() in ['temp', 't']:
                temp_var_found = True
                break
        self.assertTrue(temp_var_found, "No temperature variable found in output file")
        
        # Surface pressure could be 'ps' or 'PSFC'
        pressure_var_found = False
        for var_name in dataset.variables:
            if var_name.lower() in ['ps', 'psfc']:
                pressure_var_found = True
                break
        self.assertTrue(pressure_var_found, "No surface pressure variable found in output file")
        
        # Check for required hybrid coordinate parameters
        self.assertIn('ak', dataset.variables, "Hybrid coordinate 'ak' not found in output file.")
        self.assertIn('bk', dataset.variables, "Hybrid coordinate 'bk' not found in output file.")
        self.assertIn('pfull', dataset.variables, "Pressure levels (pfull) not found in output file.")
        self.assertIn('phalf', dataset.variables, "Pressure interfaces (phalf) not found in output file.")
        
        # Check for core coordinates
        self.assertIn('lat', dataset.variables, "Latitude variable not found in output file.")
        self.assertIn('lon', dataset.variables, "Longitude variable not found in output file.")
        self.assertIn('time', dataset.variables, "Time variable not found in output file.")
        
        dataset.close()
        
    def test_bin_average(self):
        """Test the binning and averaging functionality for EMARS."""
        # Run MarsFormat with the correct arguments including -ba for bin_average
        result = self.run_mars_format([os.path.basename(self.emars_file), "-gcm", "emars", "-ba", "2"])
        
        # Check that the command executed successfully
        self.assertEqual(result.returncode, 0, f"MarsFormat.py failed during bin average: {result.stderr}")
        
        # Check that the output file was created
        output_file = os.path.join(self.test_dir, "emars_test_average.nc")
        self.assertTrue(os.path.exists(output_file), f"Output file {output_file} was not created.")
        
        # Open the output file and check that it contains expected variables
        dataset = nc.Dataset(output_file, 'r')
        print(f"Variables in {output_file}: {list(dataset.variables.keys())}")
        
        # Check key variables
        self.assertIn('temp', dataset.variables, "Temperature variable not found in average output file.")
        self.assertIn('ps', dataset.variables, "Surface pressure variable not found in average output file.")
        
        # Check that time dimension has fewer steps due to averaging
        # Original had 840 time steps, with 2-day bins we should have fewer
        self.assertLess(len(dataset.dimensions['time']), 840, 
                        "Time dimension not properly reduced by binning.")
        
        dataset.close()
        
    def test_retain_names(self):
        """Test the -rn (retain_names) option with EMARS."""
        # Run MarsFormat with the -rn option
        result = self.run_mars_format([os.path.basename(self.emars_file), "-gcm", "emars", "-rn"])
        
        # Check that the command executed successfully
        self.assertEqual(result.returncode, 0, f"MarsFormat.py failed with retain_names: {result.stderr}")
        
        # Check that the output file was created with the expected name
        output_file = os.path.join(self.test_dir, "emars_test_nat_daily.nc")
        self.assertTrue(os.path.exists(output_file), f"Output file {output_file} was not created.")
        
        # Check that the original variable names were retained
        dataset = nc.Dataset(output_file, 'r')
        print(f"Variables in {output_file}: {list(dataset.variables.keys())}")
        
        # Should find original names like 'T' instead of 'temp'
        self.assertIn('T', dataset.variables, "Original temperature variable 'T' not found in output file.")
        
        dataset.close()
        
    def test_debug_flag(self):
        """Test the --debug flag functionality."""
        # Run MarsFormat with the --debug flag
        result = self.run_mars_format([os.path.basename(self.emars_file), "-gcm", "emars", "--debug"])
        
        # Check that the command executed successfully
        self.assertEqual(result.returncode, 0, f"MarsFormat.py failed with debug flag: {result.stderr}")
        
        # With debug flag, the output should contain more detailed information
        # Check for certain debug messages that would only appear with --debug
        self.assertIn("Current variables at top", result.stdout, 
                     "Debug output not found with --debug flag.")
                     
        # Check that the output file was created properly despite the debug flag
        output_file = os.path.join(self.test_dir, "emars_test_daily.nc")
        self.assertTrue(os.path.exists(output_file), 
                       f"Output file not created with --debug flag.")


if __name__ == '__main__':
    unittest.main()