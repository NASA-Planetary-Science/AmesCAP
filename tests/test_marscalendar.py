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
            if lines[i].strip() and not lines[i].strip().startswith('\n'):  # Skip empty lines
                # Updated regex pattern to match actual output format
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
        
        # Verify exact value based on your example output (sol=192.57 for Ls=90)
        self.assertAlmostEqual(values[0][1], 192.57, places=1)
        
    def test_ls_to_sol_multiple_values(self):
        """Test converting multiple individual Ls values to sols"""
        # Test specific known Ls values instead of using a range
        test_cases = [
            (0.0, 1.0, 1.0),     # Ls=0 should give sol ~0
            (30.0, 55.0, 65.0),  # Ls=30 should give sol ~60.81 (allow range)
            (60.0, 120.0, 130.0), # Ls=60 should give sol ~125.95 (allow range)
            (90.0, 185.0, 200.0)  # Ls=90 should give sol ~192.57 (allow range)
        ]
        
        for ls_val, min_sol, max_sol in test_cases:
            result = self.run_mars_calendar(['-ls', str(ls_val)])
            values = self.extract_values_from_output(result.stdout)
            
            # Check that we got a result
            self.assertTrue(len(values) > 0, f"No values found for Ls={ls_val}")
            
            # Check that the input Ls is correct
            self.assertAlmostEqual(values[0][0], ls_val, places=1)
            
            # Check that the output sol is in the expected range
            self.assertGreaterEqual(values[0][1], min_sol, 
                                f"Sol value {values[0][1]} for Ls={ls_val} is less than expected minimum {min_sol}")
            self.assertLessEqual(values[0][1], max_sol, 
                            f"Sol value {values[0][1]} for Ls={ls_val} is greater than expected maximum {max_sol}")
    
    def test_sol_to_ls_single(self):
        """Test converting a single sol value to Ls"""
        result = self.run_mars_calendar(['-sol', '167'])
        
        # Extract the results
        values = self.extract_values_from_output(result.stdout)
        
        # Check that we got a result
        self.assertTrue(len(values) > 0, "No values found in the output")
        
        # Check that the input sol is 167
        self.assertAlmostEqual(values[0][0], 167.0, places=1)
        
        # Verify exact value based on your example output (Ls=78.46 for sol=167)
        self.assertAlmostEqual(values[0][1], 78.46, places=1)
    
    def test_sol_to_ls_multiple_values(self):
        """Test converting multiple individual sol values to Ls"""
        # Test specific known sol values instead of using a range
        test_cases = [
            (0.0, 0.0, 1.0),      # sol=0 should give Ls ~0
            (100.0, 45.0, 55.0),  # sol=100 should give Ls ~50 (allow range)
            (200.0, 85.0, 95.0),  # sol=200 should give Ls ~90 (allow range)
            (300.0, 140.0, 150.0) # sol=300 should give Ls ~145 (allow range)
        ]
        
        for sol_val, min_ls, max_ls in test_cases:
            result = self.run_mars_calendar(['-sol', str(sol_val)])
            values = self.extract_values_from_output(result.stdout)
            
            # Check that we got a result
            self.assertTrue(len(values) > 0, f"No values found for sol={sol_val}")
            
            # Check that the input sol is correct
            self.assertAlmostEqual(values[0][0], sol_val, places=1)
            
            # Check that the output Ls is in the expected range
            self.assertGreaterEqual(values[0][1], min_ls, 
                                f"Ls value {values[0][1]} for sol={sol_val} is less than expected minimum {min_ls}")
            self.assertLessEqual(values[0][1], max_ls, 
                            f"Ls value {values[0][1]} for sol={sol_val} is greater than expected maximum {max_ls}")
    
    def test_mars_year_option(self):
        """Test using the Mars Year option"""
        # Test with base case MY=0
        base_result = self.run_mars_calendar(['-ls', '90'])
        base_values = self.extract_values_from_output(base_result.stdout)
        
        # Test with MY=1
        result1 = self.run_mars_calendar(['-ls', '90', '-my', '1'])
        values1 = self.extract_values_from_output(result1.stdout)
        
        # Test with MY=2
        result2 = self.run_mars_calendar(['-ls', '90', '-my', '2'])
        values2 = self.extract_values_from_output(result2.stdout)
        
        # Check that the output matches the expected value for MY=2
        self.assertAlmostEqual(values2[0][1], 1528.57, places=1)
        
        # Test that MY=0 and MY=1 produce the same results
        # This test verifies the behavior that MY=0 is treated the same as MY=1
        # It's a behavior we want to document but discourage
        my0_result = self.run_mars_calendar(['-ls', '90', '-my', '0'])
        my0_values = self.extract_values_from_output(my0_result.stdout)
        
        # MY=0 should produce the same result as the default (no MY specified)
        self.assertAlmostEqual(my0_values[0][1], base_values[0][1], places=1)
        
        # MY=1 should produce a result approximately 668 sols greater than default
        self.assertAlmostEqual(values1[0][1], base_values[0][1] + 668, delta=2)
        
        # MY=2 should produce a result approximately 2*668 sols greater than default
        self.assertAlmostEqual(values2[0][1], base_values[0][1] + 2*668, delta=2)
        
        # Additional check to verify MY=0 and default (no MY) produce the same result
        # This test documents the potentially confusing behavior
        self.assertAlmostEqual(my0_values[0][1], base_values[0][1], places=1)
    
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