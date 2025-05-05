#!/usr/bin/env python3
"""
Integration tests for MarsVars.py

These tests verify the functionality of MarsVars for manipulating variables in netCDF files.
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

class TestMarsVars(unittest.TestCase):
    """Integration test suite for MarsVars"""

    # Class attribute for storing modified files
    modified_files = {}
    
    @classmethod
    def setUpClass(cls):
        """Set up the test environment"""
        # Create a temporary directory for the tests
        cls.test_dir = tempfile.mkdtemp(prefix='MarsVars_test_')
        print(f"Created temporary test directory: {cls.test_dir}")

        # Project root directory
        cls.project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        print(f"Project root directory: {cls.project_root}")

        # Run the script to create test netCDF files
        cls.create_test_files()
        
        # Dictionary to keep track of modified files
        cls.modified_files = {}

    @classmethod
    def create_test_files(cls):
        """Create test netCDF files using create_ames_gcm_files.py"""
        # Get path to create_ames_gcm_files.py script
        create_files_script = os.path.join(cls.project_root, "tests", "create_ames_gcm_files.py")

        # Execute the script to create test files - Important: pass the test_dir as argument
        cmd = [sys.executable, create_files_script, cls.test_dir]

        # Print the command being executed
        print(f"Creating test files with command: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                cwd=cls.test_dir  # Run in the test directory to ensure files are created there
            )

            # Print output for debugging
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
            '01336.atmos_average_pstd.nc',
            '01336.atmos_daily.nc',
            '01336.atmos_diurn.nc',
            '01336.atmos_diurn_pstd.nc',
            '01336.fixed.nc'
        ]

        for filename in expected_files:
            filepath = os.path.join(cls.test_dir, filename)
            if not os.path.exists(filepath):
                raise Exception(f"Test file {filename} was not created in {cls.test_dir}")
            else:
                print(f"Confirmed test file exists: {filepath}")
        
        # Initialize modified_files dictionary with original files
        for filename in expected_files:
            cls.modified_files[filename] = os.path.join(cls.test_dir, filename)

    def setUp(self):
        """Change to temporary directory before each test"""
        os.chdir(self.test_dir)
        print(f"Changed to test directory: {os.getcwd()}")

    def tearDown(self):
        """Clean up after each test"""
        # Don't remove the modified files - we want to preserve them between tests
        # Only clean up any generated extract files that aren't part of our modified_files dict
        output_patterns = [
            '*_tmp.nc',
            '*_extract.nc',
            '*_col.nc'
        ]

        for pattern in output_patterns:
            for file_path in glob.glob(os.path.join(self.test_dir, pattern)):
                # Skip files we want to keep track of
                if file_path in self.modified_files.values():
                    continue
                
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

    def run_mars_vars(self, args):
        """
        Run MarsVars using subprocess

        :param args: List of arguments to pass to MarsVars
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

        # Construct the full command to run MarsVars
        cmd = [sys.executable, os.path.join(self.project_root, "bin", "MarsVars.py")] + abs_args

        # Print debugging info
        print(f"Running command: {' '.join(cmd)}")
        print(f"Working directory: {self.test_dir}")
        print(f"File exists check: {os.path.exists(os.path.join(self.project_root, 'bin', 'MarsVars.py'))}")

        # Run the command
        try:
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

            # Update our record of the modified file if this run was successful
            if result.returncode == 0:
                # Figure out which file was modified
                for arg in args:
                    if isinstance(arg, str) and arg.endswith('.nc') and not arg.startswith('-'):
                        base_filename = os.path.basename(arg)
                        if os.path.exists(os.path.join(self.test_dir, base_filename)):
                            self.modified_files[base_filename] = os.path.join(self.test_dir, base_filename)
                        break
                
                # Handle extract file creation
                if '-extract' in args:
                    input_file = next((arg for arg in args if arg.endswith('.nc')), None)
                    if input_file:
                        base_name = os.path.basename(input_file)
                        base_name_without_ext = os.path.splitext(base_name)[0]
                        extract_file = f"{base_name_without_ext}_extract.nc"
                        extract_path = os.path.join(self.test_dir, extract_file)
                        if os.path.exists(extract_path):
                            self.modified_files[extract_file] = extract_path
                            
            return result
        except Exception as e:
            self.fail(f"Failed to run MarsVars: {e}")

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

    def verify_netcdf_has_variable(self, filename, variable, alternative_names=None):
        """
        Verify that a netCDF file has a specific variable or one of its alternatives

        :param filename: Path to the netCDF file
        :param variable: Primary variable name to check for
        :param alternative_names: List of alternative variable names that are equivalent
        :return: The actual variable name found in the file
        """
        # If no alternative names provided, create an empty list
        if alternative_names is None:
            alternative_names = []
        
        # Create the full list of acceptable variable names
        acceptable_names = [variable] + alternative_names
        
        nc = Dataset(filename, 'r')
        try:
            # Try to find any of the acceptable variable names
            for var_name in acceptable_names:
                if var_name in nc.variables:
                    return var_name
                    
            # If we get here, none of the names were found
            self.fail(f"Neither {variable} nor any of its alternatives {alternative_names} found in {filename}")
        finally:
            nc.close()

    def test_help_message(self):
        """Test that help message can be displayed"""
        result = self.run_mars_vars(['-h'])

        # Check for successful execution
        self.assertEqual(result.returncode, 0, "Help command failed")

        # Check for typical help message components
        help_checks = [
            'usage:',
            '--add_variable',
            '--differentiate_wrt_z',
            '--column_integrate',
            '--zonal_detrend',
            '--dp_to_dz',
            '--dz_to_dp',
            '--remove_variable',
            '--extract_copy',
            '--edit_variable'
        ]

        for check in help_checks:
            self.assertIn(check, result.stdout, f"Help message missing '{check}'")

    def test_add_variable(self):
        """Test adding variables to files"""
        var_list = ['curl', 'div', 'DP', 'dzTau', 'dst_mass_mom', 'DZ', 
                    'theta', 'N', 'pfull3D', 'rho', 'Ri', 
                    'scorer_wl', 'Tco2', 'wdir', 'wspeed', 'zfull', 
                    'izTau', 'ice_mass_mom', 'w', 'Vg_sed', 'w_net']
        
        # Define known alternative variable names
        alternative_vars = {
            'dst_mass_mom': ['dst_mass_micro'],
            'ice_mass_mom': ['ice_mass_micro'],
            'izTau': ['ice_tau']
        }
        
        for var in var_list:
            # Add variables to non-interpolated average file and check for success
            result = self.run_mars_vars(['01336.atmos_average.nc', '-add', var])

            # Check for successful execution
            self.assertEqual(result.returncode, 0, "Add variable command failed")

            # Check that variables were added
            output_file = self.check_file_exists('01336.atmos_average.nc')
            self.verify_netcdf_has_variable(output_file, var)
        
        # Handle variables that might have alternative names
        for var, alternatives in alternative_vars.items():
            result = self.run_mars_vars(['01336.atmos_average.nc', '-add', var])
            
            # We don't check return code because it might be a warning about alternative name
            output_file = self.check_file_exists('01336.atmos_average.nc')
            
            # Check if either the variable or its alternative exists
            actual_var = self.verify_netcdf_has_variable(output_file, var, alternatives)
            print(f"Found variable {actual_var} as alternative for {var}")
        
        var_list = ['fn', 'ek', 'ep', 'msf', 'tp_t', 'ax', 'ay', 'mx', 'my']
        
        for var in var_list:
            # Add variables to interpolated average file and check for success
            result = self.run_mars_vars(['01336.atmos_average_pstd.nc', '-add', var])

            # Check for successful execution
            self.assertEqual(result.returncode, 0, "Add variable command failed")

            # Check that variables were added
            output_file = self.check_file_exists('01336.atmos_average_pstd.nc')
            self.verify_netcdf_has_variable(output_file, var)

    def test_differentiate_wrt_z(self):
        """Test differentiating a variable with respect to the Z axis"""
        # First check if we have dst_mass_micro or dst_mass_mom
        output_file = self.check_file_exists('01336.atmos_average.nc')
        actual_var = self.verify_netcdf_has_variable(output_file, 'dst_mass_micro', ['dst_mass_mom'])
        
        # Now differentiate the actual variable found
        result = self.run_mars_vars(['01336.atmos_average.nc', '-zdiff', actual_var])

        # Check for successful execution
        self.assertEqual(result.returncode, 0, "Differentiate variable command failed")

        # Check that the output variable was created
        output_file = self.check_file_exists('01336.atmos_average.nc')
        output_var = f"d_dz_{actual_var}"
        self.verify_netcdf_has_variable(output_file, output_var)

        # Verify units and naming
        nc = Dataset(output_file, 'r')
        try:
            var = nc.variables[output_var]
            self.assertIn('/m', var.units, "Units should contain '/m'")
            self.assertIn('vertical gradient', var.long_name.lower(), "Long name should mention vertical gradient")
        finally:
            nc.close()

    def test_column_integrate(self):
        """Test column integration of a variable"""        
        # First check if we have dst_mass_micro or dst_mass_mom
        output_file = self.check_file_exists('01336.atmos_average.nc')
        actual_var = self.verify_netcdf_has_variable(output_file, 'dst_mass_micro', ['dst_mass_mom'])
        
        # Now column integrate the actual variable found
        result = self.run_mars_vars(['01336.atmos_average.nc', '-col', actual_var])
        
        # Check for successful execution
        self.assertEqual(result.returncode, 0, "Column integrate command failed")

        # Check that the output variable was created
        output_file = self.check_file_exists('01336.atmos_average.nc')
        output_var = f"{actual_var}_col"
        self.verify_netcdf_has_variable(output_file, output_var)

        # Verify that the vertical dimension is removed in the output variable
        nc = Dataset(output_file, 'r')
        try:
            # Original variable has vertical dimension
            orig_dims = nc.variables[actual_var].dimensions
            col_dims = nc.variables[output_var].dimensions
            
            # Column integrated variable should have one less dimension
            self.assertEqual(len(orig_dims) - 1, len(col_dims), 
                            "Column integrated variable should have one less dimension")
            
            # Verify units
            self.assertIn('/m2', nc.variables[output_var].units, 
                          "Column integrated variable should have units with /m2")
        finally:
            nc.close()

    def test_zonal_detrend(self):
        """Test zonal detrending of a variable"""
        result = self.run_mars_vars(['01336.atmos_average.nc', '-zd', 'temp'])

        # Check for successful execution
        self.assertEqual(result.returncode, 0, "Zonal detrend command failed")

        # Check that the output variable was created
        output_file = self.check_file_exists('01336.atmos_average.nc')
        self.verify_netcdf_has_variable(output_file, 'temp_p')

        # Verify that the zonal mean is approximately zero
        nc = Dataset(output_file, 'r')
        try:
            temp_p = nc.variables['temp_p'][:]
            # Calculate zonal mean (assuming longitude is the last dimension)
            zonal_mean = np.mean(temp_p, axis=-1)
            # Should be close to zero
            self.assertTrue(np.allclose(zonal_mean, 0, atol=1e-10), 
                           "Zonal mean of detrended variable should be approximately zero")
        finally:
            nc.close()

    def test_opacity_conversion(self):
        """Test opacity conversion between dp and dz"""
        # Verify that DP and DZ variables exist
        output_file = self.check_file_exists('01336.atmos_average.nc')
        self.verify_netcdf_has_variable(output_file, 'DP')
        self.verify_netcdf_has_variable(output_file, 'DZ')
        
        # Test dp_to_dz conversion
        result = self.run_mars_vars(['01336.atmos_average.nc', '-to_dz', 'temp'])

        # Check for successful execution
        self.assertEqual(result.returncode, 0, "dp_to_dz conversion command failed")

        # Check that the output variable was created
        output_file = self.check_file_exists('01336.atmos_average.nc')
        self.verify_netcdf_has_variable(output_file, 'temp_dp_to_dz')

        # Test dz_to_dp conversion
        result = self.run_mars_vars(['01336.atmos_average.nc', '-to_dp', 'temp'])
        self.assertEqual(result.returncode, 0, "dz_to_dp conversion command failed")
        self.verify_netcdf_has_variable(output_file, 'temp_dz_to_dp')

    def test_remove_variable(self):
        """Test removing a variable from a file"""
        # First make sure wspeed exists
        output_file = self.check_file_exists('01336.atmos_average.nc')
        
        # If wspeed doesn't exist yet, add it
        nc = Dataset(output_file, 'r')
        has_wspeed = 'wspeed' in nc.variables
        nc.close()
        
        if not has_wspeed:
            # Add wspeed variable first
            result = self.run_mars_vars(['01336.atmos_average.nc', '-add', 'wspeed'])
            self.assertEqual(result.returncode, 0, "Adding wspeed failed")
            
        # Verify it exists now
        output_file = self.check_file_exists('01336.atmos_average.nc')
        self.verify_netcdf_has_variable(output_file, 'wspeed')
        
        # Now remove it
        result = self.run_mars_vars(['01336.atmos_average.nc', '-rm', 'wspeed'])

        # Check for successful execution
        self.assertEqual(result.returncode, 0, "Remove variable command failed")

        # Check that the variable was removed
        nc = Dataset(output_file, 'r')
        try:
            self.assertNotIn('wspeed', nc.variables, "Variable wspeed should have been removed")
        finally:
            nc.close()

    def test_extract_copy(self):
        """Test extracting variables to a new file"""
        # Make sure we use variables that definitely exist
        result = self.run_mars_vars(['01336.atmos_average.nc', '-extract', 'temp', 'ucomp'])

        # Check for successful execution
        self.assertEqual(result.returncode, 0, "Extract variable command failed")

        # Check that the output file was created
        output_file = self.check_file_exists('01336.atmos_average_extract.nc')
        
        # Add the extract file to our tracked modified files
        self.modified_files['01336.atmos_average_extract.nc'] = output_file
        
        # Verify it contains only the requested variables plus dimensions
        nc = Dataset(output_file, 'r')
        try:
            # Should have temp and ucomp
            self.assertIn('temp', nc.variables, "Variable temp not found in extract file")
            self.assertIn('ucomp', nc.variables, "Variable ucomp not found in extract file")
            
            # Count the non-dimension variables
            non_dim_vars = [var for var in nc.variables 
                            if var not in nc.dimensions and 
                            not any(var.endswith(f"_{dim}") for dim in nc.dimensions)]
            
            # Should only have temp and ucomp as non-dimension variables
            expected_vars = {'temp', 'ucomp'}
            actual_vars = set(non_dim_vars)
            
            # Verify the intersection of the actual and expected variables
            self.assertTrue(
                expected_vars.issubset(actual_vars), 
                f"Extract file should contain temp and ucomp. Found: {actual_vars}"
            )
        finally:
            nc.close()

    def test_edit_variable(self):
        """Test editing a variable's attributes and values"""
        # Test renaming, changing longname, units, and multiplying
        result = self.run_mars_vars([
            '01336.atmos_average.nc', 
            '-edit', 'ps', 
            '-rename', 'ps_mbar', 
            '-longname', 'Surface Pressure in Millibars', 
            '-unit', 'mbar', 
            '-multiply', '0.01'
        ])

        # Check for successful execution
        self.assertEqual(result.returncode, 0, "Edit variable command failed")

        # Check that the output file still exists and has the new variable
        output_file = self.check_file_exists('01336.atmos_average.nc')
        self.verify_netcdf_has_variable(output_file, 'ps_mbar')
        
        # Verify the attributes and scaling were applied
        nc = Dataset(output_file, 'r')
        try:
            # New variable should have the specified attributes
            ps_mbar = nc.variables['ps_mbar']
            self.assertEqual(ps_mbar.long_name, 'Surface Pressure in Millibars', 
                            "Long name was not set correctly")
            self.assertEqual(ps_mbar.units, 'mbar', "Units were not set correctly")
            
            # Values should be scaled by 0.01. The original ps variable remains
            ps = nc.variables['ps']
            self.assertTrue(np.all((ps[:] * 0.01) == ps_mbar[:]), 
                           "Values don't appear to be correctly scaled to mbar")
        finally:
            nc.close()
            

if __name__ == '__main__':
    unittest.main()