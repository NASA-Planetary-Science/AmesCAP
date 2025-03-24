#!/usr/bin/env python3
"""
Integration tests for MarsCalendar.py

These tests verify the functionality of MarsCalendar for converting between
Martian solar longitude (Ls) and sol values.
"""

import os
import sys
import unittest
import shutil
import tempfile
import argparse
import subprocess
import re

class TestMarsCalendar(unittest.TestCase):
    """Integration test suite for MarsCalendar"""
    
    @classmethod
    def setUpClass(cls):
        """Set up the test environment"""
        # Create a temporary directory in the user's home directory
        cls.test_dir = os.path.join(os.path.expanduser('~'), 'MarsCalendar_test')
        os.makedirs(cls.test_dir, exist_ok=True)
        
        # Project root directory
        cls.project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    
    def setUp(self):
        """Change to temporary directory before each test"""
        os.chdir(self.test_dir)
    
    @classmethod
    def tearDownClass(cls):
        """Clean up the test environment"""
        try:
            shutil.rmtree(cls.test_dir, ignore_errors=True)
        except Exception:
            print(f"Warning: Could not remove test directory {cls.test_dir}")
    
    def run_mars_calendar(self, args):
        """
        Run MarsCalendar using subprocess to avoid import-time argument parsing
        
        :param args: List of arguments to pass to MarsCalendar
        """
        # Construct the full command to run MarsCalendar
        cmd = [sys.executable, '-m', 'bin.MarsCalendar'] + args
        
        # Run the command
        try:
            result = subprocess.run(
                cmd, 
                capture_output=True, 
                text=True, 
                cwd=self.test_dir,
                env=dict(os.environ, PWD=self.test_dir)
            )
            
            # Check if the command was successful
            if result.returncode != 0:
                self.fail(f"MarsCalendar failed with error: {result.stderr}")
            
            return result
        except Exception as e:
            self.fail(f"Failed to run MarsCalendar: {e}")
    
    def extract_values_from_output(self, output):
        """
        Extract the values from the MarsCalendar output table
        
        Returns a list of tuples (input, output) from the rows in the table
        """
        # Split the output by lines and find the lines after the header
        lines = output.strip().split('\n')
        table_start = 0
        for i, line in enumerate(lines):
            if '-----' in line:
                table_start = i + 1
                break
        
        values = []
        for i in range(table_start, len(lines)):
            if lines[i].strip():  # Skip empty lines
                # Parse the line with format " 350.00  |  626.17"
                match = re.search(r'(\d+\.\d+)\s+\|\s+(\d+\.\d+)', lines[i])
                if match:
                    input_val = float(match.group(1))
                    output_val = float(match.group(2))
                    values.append((input_val, output_val))
        
        return values
    
    def test_ls_to_sol_single(self):
        """Test converting a single Ls value to sol"""
        result = self.run_mars_calendar(['-ls', '90'])
        
        # Extract the results
        values = self.extract_values_from_output(result.stdout)
        
        # Check that we got a result
        self.assertTrue(len(values) > 0, "No values found in the output")
        
        # Check that the input Ls is 90
        self.assertAlmostEqual(values[0][0], 90.0, places=1)
        
        # Check that the output sol is approximately 167
        # (This value should be verified against expected results for Ls=90)
        self.assertGreater(values[0][1], 160)
        self.assertLess(values[0][1], 175)
    
    def test_ls_to_sol_range(self):
        """Test converting a range of Ls values to sols"""
        result = self.run_mars_calendar(['-ls', '0', '90', '30'])
        
        # Extract the results
        values = self.extract_values_from_output(result.stdout)
        
        # Check that we got the expected number of results (0, 30, 60, 90)
        self.assertEqual(len(values), 4, "Expected 4 values in output")
        
        # Check that the Ls values are as expected
        expected_ls = [0.0, 30.0, 60.0, 90.0]
        for i, (ls_val, _) in enumerate(values):
            self.assertAlmostEqual(ls_val, expected_ls[i], places=1)
    
    def test_sol_to_ls_single(self):
        """Test converting a single sol value to Ls"""
        result = self.run_mars_calendar(['-sol', '167'])
        
        # Extract the results
        values = self.extract_values_from_output(result.stdout)
        
        # Check that we got a result
        self.assertTrue(len(values) > 0, "No values found in the output")
        
        # Check that the input sol is 167
        self.assertAlmostEqual(values[0][0], 167.0, places=1)
        
        # Check that the output Ls is approximately 90
        # (This value should be verified against expected results for sol=167)
        self.assertGreater(values[0][1], 85)
        self.assertLess(values[0][1], 95)
    
    def test_sol_to_ls_range(self):
        """Test converting a range of sol values to Ls"""
        result = self.run_mars_calendar(['-sol', '0', '300', '100'])
        
        # Extract the results
        values = self.extract_values_from_output(result.stdout)
        
        # Check that we got the expected number of results (0, 100, 200, 300)
        self.assertEqual(len(values), 4, "Expected 4 values in output")
        
        # Check that the sol values are as expected
        expected_sols = [0.0, 100.0, 200.0, 300.0]
        for i, (sol_val, _) in enumerate(values):
            self.assertAlmostEqual(sol_val, expected_sols[i], places=1)
    
    def test_mars_year_option(self):
        """Test using the Mars Year option"""
        # Test with MY=1
        result = self.run_mars_calendar(['-ls', '90', '-my', '1'])
        
        # Extract the results
        values = self.extract_values_from_output(result.stdout)
        
        # Check that the output sol is shifted by 668 (one Mars year)
        base_result = self.run_mars_calendar(['-ls', '90'])
        base_values = self.extract_values_from_output(base_result.stdout)
        
        # The sol for MY=1 should be approximately 668 more than for MY=0
        self.assertAlmostEqual(values[0][1], base_values[0][1] + 668, delta=1)
    
    def test_continuous_option(self):
        """Test using the continuous Ls option"""
        # Test with a sol value that should give Ls > 360 with continuous option
        result = self.run_mars_calendar(['-sol', '700', '-c'])
        
        # Extract the results
        values = self.extract_values_from_output(result.stdout)
        
        # With continuous flag, Ls should be > 360
        self.assertGreater(values[0][1], 360)
        
        # Compare with non-continuous option
        regular_result = self.run_mars_calendar(['-sol', '700'])
        regular_values = self.extract_values_from_output(regular_result.stdout)
        
        # Without continuous flag, Ls should be < 360
        self.assertLess(regular_values[0][1], 360)
        
        # The difference should be approximately 360
        self.assertAlmostEqual(values[0][1], regular_values[0][1] + 360, delta=1)
    
    def test_help_message(self):
        """Test that help message can be displayed"""
        result = self.run_mars_calendar(['-h'])
        
        # Check that something was printed
        self.assertTrue(len(result.stdout) > 0, "No help message generated")
        
        # Check for typical help message components
        help_checks = [
            'usage:',
            '-ls',
            '-sol',
            '-my',
            '--continuous',
            '--debug'
        ]
        
        for check in help_checks:
            self.assertIn(check, result.stdout.lower(), f"Help message missing '{check}'")


if __name__ == '__main__':
    unittest.main()