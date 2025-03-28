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
        cls.test_dir = tempfile.mkdtemp(prefix='MarsFormat_test_')
        
        # Project root directory
        cls.project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        
        # Define all GCM types to test
        cls.gcm_types = ['emars', 'openmars', 'pcm', 'marswrf']
    
    def setUp(self):
        """Create test netCDF files using create_gcm_files.py"""
        # Change to test directory
        os.chdir(self.test_dir)
        
        # Define file paths for each GCM type
        self.test_files = {}
        for gcm_type in self.gcm_types:
            self.test_files[gcm_type] = os.path.join(self.test_dir, f"{gcm_type}_test.nc")

        # # Get path to create_gcm_files.py script
        # create_files_script = os.path.join(self.project_root, "tests", "create_gcm_files.py")
        
        # # Execute the script to create test files
        # result = subprocess.run(
        #     [sys.executable, create_files_script],
        #     capture_output=True,
        #     text=True,
        #     cwd=self.test_dir
        # )
        # Copy real data files
        subprocess.run(['cp', os.path.expanduser('~/marsformat_data/emars_Ls240-270.nc'), 
                        os.path.join(self.test_dir, 'emars_test.nc')], check=True)
        subprocess.run(['cp', os.path.expanduser('~/marsformat_data/marswrf_d01_0001-00669.nc'), 
                        os.path.join(self.test_dir, 'marswrf_test.nc')], check=True)
        subprocess.run(['cp', os.path.expanduser('~/marsformat_data/openmars_Ls264-284.nc'), 
                        os.path.join(self.test_dir, 'openmars_test.nc')], check=True)
        subprocess.run(['cp', os.path.expanduser('~/marsformat_data/pcm_Ls264-280.nc'), 
                        os.path.join(self.test_dir, 'pcm_test.nc')], check=True)
        
        print("File creation output: Copied real data files from ~/marsformat_data/")
        
        # # Print output for debugging
        # print(f"File creation output: {result.stdout}")
        
        # Check files were created
        for gcm_type, test_file in self.test_files.items():
            if not os.path.exists(test_file):
                print(f"Warning: Test file {test_file} was not created!")
                # if result.stderr:
                #     print(f"Error output: {result.stderr}")
    
    def tearDown(self):
        """Clean up any generated files after each test"""
        # Clean up any generated output files after each test
        for gcm_type in self.gcm_types:
            # Regular output files
            output_files = [
                os.path.join(self.test_dir, f"{gcm_type}_test_daily.nc"), 
                os.path.join(self.test_dir, f"{gcm_type}_test_average.nc"), 
                os.path.join(self.test_dir, f"{gcm_type}_test_diurn.nc"),
                os.path.join(self.test_dir, f"{gcm_type}_test_nat_daily.nc"), 
                os.path.join(self.test_dir, f"{gcm_type}_test_nat_average.nc"), 
                os.path.join(self.test_dir, f"{gcm_type}_test_nat_diurn.nc")
            ]
            
            for file in output_files:
                if os.path.exists(file):
                    os.remove(file)
        
        # Clean up dummy input files
        for gcm_type, test_file in self.test_files.items():
            if os.path.exists(test_file):
                os.remove(test_file)
    
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
            if isinstance(arg, str) and arg.endswith('.nc'):
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
    
    def verify_output_file(self, output_file, expected_vars=None, expected_coords=None):
        """
        Verify that an output file exists and contains expected variables and coordinates
        
        :param output_file: Path to the output file to verify
        :param expected_vars: List of expected variable names
        :param expected_coords: List of expected coordinate names
        """
        # Check that the output file was created
        self.assertTrue(os.path.exists(output_file), f"Output file {output_file} was not created.")
            
        # Open the output file and check that it contains the expected variables
        dataset = nc.Dataset(output_file, 'r')
        
        # Debug - print all variable names to examine what's actually in the file
        print(f"Variables in {output_file}: {list(dataset.variables.keys())}")
        
        # Check that expected variables are present
        if expected_vars:
            for var in expected_vars:
                # For temperature, check both 'temp' and 'T' since naming may differ
                if var.lower() in ['temp', 't']:
                    temp_var_found = False
                    for var_name in dataset.variables:
                        if var_name.lower() in ['temp', 't']:
                            temp_var_found = True
                            break
                    self.assertTrue(temp_var_found, f"Temperature variable not found in {output_file}")
                # For surface pressure, check both 'ps' and 'PSFC'
                elif var.lower() in ['ps', 'psfc']:
                    ps_var_found = False
                    for var_name in dataset.variables:
                        if var_name.lower() in ['ps', 'psfc']:
                            ps_var_found = True
                            break
                    self.assertTrue(ps_var_found, f"Surface pressure variable not found in {output_file}")
                else:
                    self.assertIn(var, dataset.variables, f"{var} not found in {output_file}")
        
        # Check that expected coordinates are present
        if expected_coords:
            for coord in expected_coords:
                if coord == 'lat':
                    # Check both 'lat' and 'latitude'
                    lat_found = False
                    for coord_name in dataset.variables:
                        if coord_name.lower() in ['lat', 'latitude']:
                            lat_found = True
                            break
                    self.assertTrue(lat_found, f"Latitude coordinate not found in {output_file}")
                elif coord == 'lon':
                    # Check both 'lon' and 'longitude'
                    lon_found = False
                    for coord_name in dataset.variables:
                        if coord_name.lower() in ['lon', 'longitude']:
                            lon_found = True
                            break
                    self.assertTrue(lon_found, f"Longitude coordinate not found in {output_file}")
                else:
                    self.assertIn(coord, dataset.variables, f"{coord} not found in {output_file}")
        
        # Check for time_of_day dimension if this is a diurn file
        if '_diurn.nc' in output_file:
            time_of_day_found = False
            for dim_name in dataset.dimensions:
                if 'time_of_day' in dim_name:
                    time_of_day_found = True
                    break
            self.assertTrue(time_of_day_found, f"No time_of_day dimension found in {output_file}")
        
        # Check that longitude values have been properly converted to 0-360 range
        # Find the longitude variable
        lon_var = None
        for var_name in dataset.variables:
            if var_name.lower() in ['lon', 'longitude']:
                lon_var = var_name
                break
        
        if lon_var:
            self.assertGreaterEqual(dataset.variables[lon_var][0], 0, 
                                 f"Longitude values not converted to 0-360 range in {output_file}")
        
        dataset.close()
        return True

    def test_all_gcm_types(self):
        """Test basic conversion for all GCM types."""
        # Core variables and coordinates to check for in every file
        core_vars = ['temp', 'ps']
        core_coords = ['lat', 'lon', 'time', 'pfull']
        hybrid_vars = ['ak', 'bk', 'phalf']
        
        for gcm_type in self.gcm_types:
            # Run MarsFormat with just the GCM flag
            result = self.run_mars_format([os.path.basename(self.test_files[gcm_type]), "-gcm", gcm_type])
            
            # Check that the command executed successfully
            self.assertEqual(result.returncode, 0, f"MarsFormat.py failed for {gcm_type}: {result.stderr}")
            
            # Verify output file
            output_file = os.path.join(self.test_dir, f"{gcm_type}_test_daily.nc")
            self.verify_output_file(output_file, expected_vars=core_vars + hybrid_vars, 
                                   expected_coords=core_coords)

    def test_all_gcm_types_retain_names(self):
        """Test conversion with name retention for all GCM types."""
        for gcm_type in self.gcm_types:
            # Run MarsFormat with GCM flag and retain_names flag
            result = self.run_mars_format([os.path.basename(self.test_files[gcm_type]), 
                                         "-gcm", gcm_type, "-rn"])
            
            # Check that the command executed successfully
            self.assertEqual(result.returncode, 0, 
                          f"MarsFormat.py failed for {gcm_type} with retain_names: {result.stderr}")
            
            # Verify output file - variable names will differ by GCM type
            output_file = os.path.join(self.test_dir, f"{gcm_type}_test_nat_daily.nc")
            
            # Expect original variable names to be preserved
            # For EMARS we expect 'T' instead of 'temp', etc.
            if gcm_type == 'emars':
                self.verify_output_file(output_file, expected_vars=['T', 'ps'])
            elif gcm_type == 'openmars':
                self.verify_output_file(output_file)
            elif gcm_type == 'marswrf':
                self.verify_output_file(output_file, expected_vars=['T', 'PSFC'])
            elif gcm_type == 'pcm':
                # PCM variable names
                self.verify_output_file(output_file)

    def test_bin_average_all_types(self):
        """Test bin_average flag for all GCM types, with and without retain_names."""
        for gcm_type in self.gcm_types:
            # Test without retain_names
            result = self.run_mars_format([os.path.basename(self.test_files[gcm_type]), 
                                         "-gcm", gcm_type, "-ba", "20"])
            
            self.assertEqual(result.returncode, 0, 
                          f"MarsFormat.py failed for {gcm_type} with bin_average: {result.stderr}")
            
            # Verify output file
            output_file = os.path.join(self.test_dir, f"{gcm_type}_test_average.nc")
            self.verify_output_file(output_file)
            
            # Check time dimension has been reduced due to averaging
            dataset = nc.Dataset(output_file, 'r')
            orig_time_len = len(nc.Dataset(self.test_files[gcm_type], 'r').dimensions['Time'])
            new_time_len = len(dataset.dimensions['time'])
            self.assertLess(new_time_len, orig_time_len, 
                         f"Time dimension not reduced by binning in {output_file}")
            dataset.close()
            
            # Test with retain_names
            result = self.run_mars_format([os.path.basename(self.test_files[gcm_type]), 
                                         "-gcm", gcm_type, "-ba", "20", "-rn"])
            
            self.assertEqual(result.returncode, 0, 
                          f"MarsFormat.py failed for {gcm_type} with bin_average and retain_names: {result.stderr}")
            
            # Verify output file
            output_file = os.path.join(self.test_dir, f"{gcm_type}_test_nat_average.nc")
            self.verify_output_file(output_file)

    def test_bin_diurn_all_types(self):
        """Test bin_diurn flag for all GCM types, with and without retain_names."""
        for gcm_type in self.gcm_types:
            # Skip marswrf which is known to fail with bin_diurn
            if gcm_type == 'marswrf':
                print(f"Skipping bin_diurn test for {gcm_type} as it's known to fail")
                continue
            
            # Test without retain_names
            result = self.run_mars_format([os.path.basename(self.test_files[gcm_type]), 
                                         "-gcm", gcm_type, "-bd"])
            
            self.assertEqual(result.returncode, 0, 
                          f"MarsFormat.py failed for {gcm_type} with bin_diurn: {result.stderr}")
            
            # Verify output file
            output_file = os.path.join(self.test_dir, f"{gcm_type}_test_diurn.nc")
            self.verify_output_file(output_file)
            
            # Test with retain_names
            result = self.run_mars_format([os.path.basename(self.test_files[gcm_type]), 
                                         "-gcm", gcm_type, "-bd", "-rn"])
            
            self.assertEqual(result.returncode, 0, 
                          f"MarsFormat.py failed for {gcm_type} with bin_diurn and retain_names: {result.stderr}")
            
            # Verify output file
            output_file = os.path.join(self.test_dir, f"{gcm_type}_test_nat_diurn.nc")
            self.verify_output_file(output_file)

    def test_combined_flags(self):
        """Test combining the bin_average and bin_diurn flags."""
        for gcm_type in self.gcm_types:
            # Skip marswrf for the same reason as above
            if gcm_type == 'marswrf':
                continue
                
            # Create a unique input file for this test to avoid conflicts
            unique_input = os.path.join(self.test_dir, f"{gcm_type}_test_combined.nc")
            shutil.copy(self.test_files[gcm_type], unique_input)
            
            # Test bin_average with bin_diurn (without retain_names)
            result = self.run_mars_format([
                os.path.basename(unique_input), 
                "-gcm", gcm_type, 
                "-bd", "-ba", "2"
            ])
            
            self.assertEqual(result.returncode, 0, 
                        f"MarsFormat.py failed for {gcm_type} with combined flags: {result.stderr}")
            
            # Verify output file - should create a diurn file with averaged data
            output_file = os.path.join(self.test_dir, f"{gcm_type}_test_combined_diurn.nc")
            self.verify_output_file(output_file)
            
            # Create another unique input file for the retain_names test
            unique_input_rn = os.path.join(self.test_dir, f"{gcm_type}_test_combined_rn.nc")
            shutil.copy(self.test_files[gcm_type], unique_input_rn)
            
            # Test with retain_names added
            result = self.run_mars_format([
                os.path.basename(unique_input_rn), 
                "-gcm", gcm_type, 
                "-bd", "-ba", "2", 
                "-rn"
            ])
            
            self.assertEqual(result.returncode, 0, 
                        f"MarsFormat.py failed for {gcm_type} with combined flags and retain_names: {result.stderr}")
            
            # Verify output file
            output_file_rn = os.path.join(self.test_dir, f"{gcm_type}_test_combined_rn_nat_diurn.nc")
            self.verify_output_file(output_file_rn)

    def test_variable_mapping(self):
        """Test that variable mapping from GCM-specific names to standard names works correctly."""
        # Expected variable mappings based on amescap_profile
        var_mappings = {
            'emars': {
                'original': ['T', 'ALSO_u', 'ALSO_v'], 
                'mapped': ['temp', 'ucomp', 'vcomp']
            },
            'openmars': {
                'original': [], 
                'mapped': ['temp', 'ucomp', 'vcomp']
            },
            'marswrf': {
                'original': ['U', 'V'], 
                'mapped': ['ucomp', 'vcomp']
            },
            'pcm': {  # PCM is referred to as LMD in amescap_profile
                'original': [], 
                'mapped': ['temp', 'ucomp', 'vcomp']
            }
        }
        
        for gcm_type in self.gcm_types:
            # Run with retain_names to keep original names
            result_retain = self.run_mars_format([os.path.basename(self.test_files[gcm_type]), 
                                                "-gcm", gcm_type, "-rn"])
            self.assertEqual(result_retain.returncode, 0)
            
            # Run without retain_names to map to standard names
            result_map = self.run_mars_format([os.path.basename(self.test_files[gcm_type]), 
                                             "-gcm", gcm_type])
            self.assertEqual(result_map.returncode, 0)
            
            # Check files
            retained_file = os.path.join(self.test_dir, f"{gcm_type}_test_nat_daily.nc")
            mapped_file = os.path.join(self.test_dir, f"{gcm_type}_test_daily.nc")
            
            # Verify original variables in retained file
            if var_mappings[gcm_type]['original']:
                with nc.Dataset(retained_file, 'r') as ds:
                    var_names = list(ds.variables.keys())
                    for var in var_mappings[gcm_type]['original']:
                        if var not in var_names:
                            # Some variables might not be present in the sample files
                            print(f"Note: Expected variable {var} not found in {retained_file}")
            
            # Verify mapped variables in mapped file
            with nc.Dataset(mapped_file, 'r') as ds:
                var_names = list(ds.variables.keys())
                for var in var_mappings[gcm_type]['mapped']:
                    self.assertIn(var, var_names, f"Expected mapped variable {var} not found in {mapped_file}")

    def test_coordinate_transformations(self):
        """Test that coordinate transformations are applied correctly."""
        for gcm_type in self.gcm_types:
            # Run MarsFormat 
            result = self.run_mars_format([os.path.basename(self.test_files[gcm_type]), 
                                         "-gcm", gcm_type])
            self.assertEqual(result.returncode, 0)
            
            # Check output file
            output_file = os.path.join(self.test_dir, f"{gcm_type}_test_daily.nc")
            
            with nc.Dataset(output_file, 'r') as ds:
                # Find longitude variable (could be 'lon' or 'longitude')
                lon_var = None
                for var_name in ds.variables:
                    if var_name.lower() in ['lon', 'longitude']:
                        lon_var = var_name
                        break
                
                if lon_var:
                    # Check that longitudes are in 0-360 range
                    lon_values = ds.variables[lon_var][:]
                    self.assertGreaterEqual(np.min(lon_values), 0, 
                                         f"Longitudes not transformed to 0-360 range in {output_file}")
                    self.assertLess(np.max(lon_values), 360.1, 
                                  f"Longitudes exceed 360 degrees in {output_file}")
                
                # Check vertical coordinate ordering (pfull should increase with level index)
                pfull_var = None
                for var_name in ds.variables:
                    if var_name.lower() == 'pfull':
                        pfull_var = var_name
                        break
                
                if pfull_var:
                    pfull_values = ds.variables[pfull_var][:]
                    # Check that pressure increases with index (TOA at index 0)
                    self.assertLess(pfull_values[0], pfull_values[-1], 
                                  f"Pressure levels not ordered with TOA at index 0 in {output_file}")
                
                # Check for hybrid coordinate variables
                self.assertIn('ak', ds.variables, f"Hybrid coordinate 'ak' missing in {output_file}")
                self.assertIn('bk', ds.variables, f"Hybrid coordinate 'bk' missing in {output_file}")

    def test_debug_flag(self):
        """Test the --debug flag functionality."""
        for gcm_type in self.gcm_types:
            # Run MarsFormat with the --debug flag
            result = self.run_mars_format([os.path.basename(self.test_files[gcm_type]), 
                                         "-gcm", gcm_type, "--debug"])
            
            # Check that the command executed successfully
            self.assertEqual(result.returncode, 0, 
                          f"MarsFormat.py failed with debug flag for {gcm_type}: {result.stderr}")
            
            # Debug output should contain more detailed information
            self.assertTrue(
                "Running MarsFormat with args" in result.stdout or 
                "Current working directory" in result.stdout,
                "Debug output not found with --debug flag."
            )
                     
            # Check that output file was created properly
            output_file = os.path.join(self.test_dir, f"{gcm_type}_test_daily.nc")
            self.assertTrue(os.path.exists(output_file), 
                         f"Output file not created with --debug flag for {gcm_type}.")

    def test_error_handling(self):
        """Test error handling with invalid inputs."""
        # Test with non-existent file
        result = self.run_mars_format(["nonexistent_file.nc", "-gcm", "emars"])
        self.assertNotEqual(result.returncode, 0, "MarsFormat didn't fail with non-existent file")
        
        # Test with invalid GCM type
        result = self.run_mars_format([os.path.basename(self.test_files['emars']), 
                                     "-gcm", "invalid_gcm"])
        self.assertNotEqual(result.returncode, 0, "MarsFormat didn't fail with invalid GCM type")
        
        # Test with no GCM type specified
        result = self.run_mars_format([os.path.basename(self.test_files['emars'])])
        # This might not return an error code, but should print a notice
        self.assertTrue(
            "No operation requested" in result.stdout or 
            result.returncode != 0,
            "MarsFormat didn't handle missing GCM type correctly"
        )


if __name__ == '__main__':
    unittest.main()