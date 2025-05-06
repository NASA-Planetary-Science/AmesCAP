#!/usr/bin/env python3
"""
Shared test class methods and functions
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

class BaseTestCase(unittest.TestCase):
    """Base class for integration tests with common setup methods"""
    
    PREFIX = "Default_test_"
    FILESCRIPT = "create_ames_gcm_files.py"
    SHORTFILE = "short"

    @classmethod
    def setUpClass(cls):
        """Set up the test environment"""
        # Create a temporary directory for the tests
        cls.test_dir = tempfile.mkdtemp(prefix=cls.PREFIX)
        print(f"Created temporary test directory: {cls.test_dir}")
        # Project root directory
        cls.project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        print(f"Project root directory: {cls.project_root}")
        # Create test files
        cls.create_test_files()

    @classmethod
    def create_test_files(cls):
        """Create test netCDF files using create_ames_gcm_files.py"""
        # Get path to create_ames_gcm_files.py script
        create_files_script = os.path.join(cls.project_root, "tests", cls.FILESCRIPT)

        # Execute the script to create test files - Important: pass the test_dir as argument
        cmd = [sys.executable, create_files_script, cls.test_dir, cls.SHORTFILE]

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
            # List files in temp directory before deleting to debug
            print(f"Files in test directory before cleanup: {os.listdir(cls.test_dir)}")
            shutil.rmtree(cls.test_dir, ignore_errors=True)
            print(f"Removed test directory: {cls.test_dir}")
        except Exception as e:
            print(f"Warning: Could not remove test directory {cls.test_dir}: {e}")

