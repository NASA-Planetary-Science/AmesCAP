#!/usr/bin/env python3
"""
Integration tests for MarsPlot.py

These tests verify the functionality of MarsPlot for visualizing netCDF files.
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

class TestMarsPlot(unittest.TestCase):
    """Integration test suite for MarsPlot"""
    
    @classmethod
    def setUpClass(cls):
        """Set up the test environment"""
        # Create a temporary directory for the tests
        cls.test_dir = tempfile.mkdtemp(prefix='MarsPlot_test_')
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
        cmd = [sys.executable, create_files_script, cls.test_dir]
        
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
        
        # Create a test template file
        cls.create_test_template()
    
    @classmethod
    def create_test_template(cls):
        """Create a simple Custom.in test file"""
        template_content = """# Test Custom.in template for MarsPlot
# Simple template with one 2D plot
<<<<<<<<<<<<<<<<<<<<<< Simulations >>>>>>>>>>>>>>>>>>>>>
ref> None
2> 
3>
=======================================================
START
HOLD ON
<<<<<<<<<<<<<<| Plot 2D lon X lat = True |>>>>>>>>>>>>>
Title          = None
Main Variable  = fixed.zsurf
Cmin, Cmax     = None
Ls 0-360       = None
Level Pa/m     = None
2nd Variable   = None
Contours Var 2 = None
Axis Options  : Lon = [None,None] | Lat = [None,None] | cmap = jet | scale = lin | proj = cart 
HOLD OFF
"""
        template_path = os.path.join(cls.test_dir, "Custom.in")
        with open(template_path, "w") as f:
            f.write(template_content)
        
        print(f"Created test template file: {template_path}")
    
    def setUp(self):
        """Change to temporary directory before each test"""
        os.chdir(self.test_dir)
        print(f"Changed to test directory: {os.getcwd()}")
    
    def tearDown(self):
        """Clean up after each test"""
        # Clean up any generated output files after each test but keep input files
        output_patterns = [
            '*_figure*.pdf',
            '*_figure*.png',
            '*_figure*.eps',
            '*.pdf',
            '*.png',
            '*.eps'
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
    
    def run_mars_plot(self, args, expected_success=True):
        """
        Run MarsPlot using subprocess
        
        :param args: List of arguments to pass to MarsPlot
        :param expected_success: Whether the command is expected to succeed
        :return: subprocess result object
        """
        # Construct the full command to run MarsPlot
        cmd = [sys.executable, os.path.join(self.project_root, "bin", "MarsPlot.py")] + args
        
        # Print debugging info
        print(f"Running command: {' '.join(cmd)}")
        print(f"Working directory: {self.test_dir}")
        print(f"File exists check: {os.path.exists(os.path.join(self.project_root, 'bin', 'MarsPlot.py'))}")
        
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
                self.assertEqual(result.returncode, 0, f"MarsPlot command failed with code {result.returncode}")
            else:
                self.assertNotEqual(result.returncode, 0, "MarsPlot command succeeded but was expected to fail")
            
            return result
        except Exception as e:
            self.fail(f"Failed to run MarsPlot: {e}")
    
    def check_file_exists(self, filename):
        """
        Check if a file exists and is not empty
        
        :param filename: Filename to check
        """
        filepath = os.path.join(self.test_dir, filename)
        self.assertTrue(os.path.exists(filepath), f"File {filename} does not exist")
        self.assertGreater(os.path.getsize(filepath), 0, f"File {filename} is empty")
        return filepath
    
    def test_help_message(self):
        """Test that help message can be displayed"""
        result = self.run_mars_plot(['-h'])
        
        # Check for typical help message components
        help_checks = [
            'usage:',
            'template_file',
            '--inspect_file',
            '--generate_template',
            '--date',
            '--stack_years',
            '--figure_filetype',
            '--portrait_mode',
            '--pixel_width',
            '--directory'
        ]
        
        for check in help_checks:
            self.assertIn(check, result.stdout, f"Help message missing '{check}'")
    
    def test_generate_template(self):
        """Test generating a template file"""
        # Delete any existing Custom.in to ensure we're testing the generation
        if os.path.exists(os.path.join(self.test_dir, "Custom.in")):
            os.remove(os.path.join(self.test_dir, "Custom.in"))
        
        result = self.run_mars_plot(['-template'])
        
        # Check that the template file was created
        template_file = self.check_file_exists("Custom.in")
        
        # Verify that the template has expected content
        with open(template_file, 'r') as f:
            content = f.read()
        
        expected_sections = [
            "Simulations",
            "START",
            "HOLD",
            "Title",
            "Plot",
            "Variable"
        ]
        
        for section in expected_sections:
            self.assertIn(section, content, f"Template missing expected section: {section}")
    
    def test_generate_trimmed_template(self):
        """Test generating a trimmed template file"""
        # Delete any existing Custom.in to ensure we're testing the generation
        if os.path.exists(os.path.join(self.test_dir, "Custom.in")):
            os.remove(os.path.join(self.test_dir, "Custom.in"))
        
        result = self.run_mars_plot(['-template', '-trim'])
        
        # Check that the template file was created
        template_file = self.check_file_exists("Custom.in")
        
        # Verify that the template has expected content but is shorter
        with open(template_file, 'r') as f:
            content = f.read()
        
        # The trimmed version should still have these sections
        expected_sections = [
            "Simulations",
            "START",
            "HOLD",
            "Title",
            "Plot",
            "Variable"
        ]
        
        for section in expected_sections:
            self.assertIn(section, content, f"Template missing expected section: {section}")
    
    def test_inspect_file(self):
        """Test inspecting a netCDF file"""
        result = self.run_mars_plot(['-i', '01336.atmos_average.nc'])
        
        # Check for typical inspect output
        inspect_checks = [
            'DIMENSIONS',
            'CONTENT'
        ]
        
        for check in inspect_checks:
            self.assertIn(check, result.stdout, f"Inspect output missing '{check}'")
    
    def test_inspect_with_values(self):
        """Test inspecting a netCDF file with values for a variable"""
        # First we need to know a valid variable name from the file
        nc = Dataset(os.path.join(self.test_dir, '01336.atmos_average.nc'), 'r')
        try:
            # Get the first variable that's not a dimension
            var_name = None
            for var in nc.variables:
                if var not in nc.dimensions:
                    var_name = var
                    break
            
            if var_name is None:
                self.skipTest("No suitable variable found for test_inspect_with_values")
        finally:
            nc.close()
        
        result = self.run_mars_plot(['-i', '01336.atmos_average.nc', '-values', var_name])
        
        # Check for the variable name in output
        self.assertIn(var_name, result.stdout, f"Variable {var_name} not found in inspect output")
    
    def test_inspect_with_statistics(self):
        """Test inspecting a netCDF file with statistics for a variable"""
        # First we need to know a valid variable name from the file
        nc = Dataset(os.path.join(self.test_dir, '01336.atmos_average.nc'), 'r')
        try:
            # Get the first variable that's not a dimension
            var_name = None
            for var in nc.variables:
                if var not in nc.dimensions:
                    var_name = var
                    break
            
            if var_name is None:
                self.skipTest("No suitable variable found for test_inspect_with_statistics")
        finally:
            nc.close()
        
        result = self.run_mars_plot(['-i', '01336.atmos_average.nc', '-stats', var_name])
        
        # Check for statistics in output
        stats_checks = [
            'min',
            'max',
            'mean'
        ]
        
        for check in stats_checks:
            self.assertIn(check, result.stdout.lower(), f"Statistics output missing '{check}'")
    
    def test_run_template(self):
        """Test running MarsPlot with a template file"""
        result = self.run_mars_plot(['Custom.in'])
        
        # Check for successful execution based on typical output
        success_checks = [
            'Reading Custom.in',
            'generated'  # This might need adjustment based on actual output
        ]
        
        for check in success_checks:
            self.assertIn(check, result.stdout, f"Output missing expected message: '{check}'")
        
        # Check that at least one output file was created
        # The exact name depends on implementation, so we check for any PDF
        pdf_files = glob.glob(os.path.join(self.test_dir, "*.pdf"))
        self.assertTrue(len(pdf_files) > 0, "No PDF output file generated")
    
    def test_run_with_specific_date(self):
        """Test running MarsPlot with a specific date"""
        result = self.run_mars_plot(['Custom.in', '-d', '01336'])
        
        # Check for successful execution
        success_checks = [
            'Reading Custom.in',
            'Done' 
        ]
        
        for check in success_checks:
            self.assertIn(check, result.stdout, f"Output missing expected message: '{check}'")
    
    def test_run_with_png_output(self):
        """Test running MarsPlot with PNG output format"""
        result = self.run_mars_plot(['Custom.in', '-ftype', 'png'])
        
        # Check for successful execution
        self.assertIn('Reading Custom.in', result.stdout)
        
        # Check that a PNG file was created
        png_files = glob.glob(os.path.join(self.test_dir, "plots/*.png"))
        self.assertTrue(len(png_files) > 0, "No PNG output file generated")
    
    def test_run_in_portrait_mode(self):
        """Test running MarsPlot in portrait mode"""
        result = self.run_mars_plot(['Custom.in', '-portrait', '--debug'])
        
        # Check for successful execution
        self.assertIn('Reading Custom.in', result.stdout)
        
        # Since we can't easily verify the aspect ratio of the output,
        # we just check that a file was created
        output_files = glob.glob(os.path.join(self.test_dir, "*.pdf"))
        self.assertTrue(len(output_files) > 0, "No output file generated in portrait mode")
    
    def test_run_with_custom_pixel_width(self):
        """Test running MarsPlot with custom pixel width"""
        result = self.run_mars_plot(['Custom.in', '-pw', '1000'])
        
        # Check for successful execution
        self.assertIn('Reading Custom.in', result.stdout)
        
        # Since we can't easily verify the dimensions of the output,
        # we just check that a file was created
        output_files = glob.glob(os.path.join(self.test_dir, "*.pdf"))
        self.assertTrue(len(output_files) > 0, "No output file generated with custom pixel width")
    
    def test_stack_years_option(self):
        """Test running MarsPlot with stack years option"""
        result = self.run_mars_plot(['Custom.in', '-sy'])
        
        # Check for successful execution
        self.assertIn('Reading Custom.in', result.stdout)     
    
    def test_debug_mode(self):
        """Test running MarsPlot in debug mode"""
        result = self.run_mars_plot(['Invalid.txt', '--debug'],expected_success=False)
        
        # Debug mode releases the errors to standard output
        debug_indicators = [
            'error',
            'ERROR',
            'Failed',
            'failed'
        ]
        
        # At least one of these should appear in the output
        self.assertTrue(any(indicator in result.stderr for indicator in debug_indicators),
                       "No indication that debug mode was active")
    
    def test_invalid_template_extension(self):
        """Test error handling with invalid template extension"""
        # Create a test file with wrong extension
        with open(os.path.join(self.test_dir, "Invalid.txt"), "w") as f:
            f.write("Some content")
        
        # This should fail because the template must be a .in file
        result = self.run_mars_plot(['Invalid.txt'], expected_success=False)
        
        # Check for error message
        self.assertIn('not a \'.in\'', result.stderr, "Missing error message about invalid file extension")
    
    def test_invalid_netcdf_file(self):
        """Test error handling with invalid netCDF file"""
        # Create a test file that's not a valid netCDF
        with open(os.path.join(self.test_dir, "Invalid.nc"), "w") as f:
            f.write("This is not a netCDF file")
        
        # This should fail because the file is not a valid netCDF
        result = self.run_mars_plot(['-i', 'Invalid.nc'], expected_success=False)
        
        # Check for error message (may vary depending on implementation)
        error_indicators = [
            'error',
            'failed',
            'invalid'
        ]
        
        # At least one of these should appear in stderr
        error_output = result.stdout.lower()
        print(error_output)
        self.assertTrue(any(indicator in error_output for indicator in error_indicators),
                       "No indication of error with invalid netCDF file")


if __name__ == '__main__':
    unittest.main()
