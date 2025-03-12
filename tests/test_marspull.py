#!/usr/bin/env python3
"""
Tests for MarsPull.py

This test module mocks HTTP requests to test the functionality of MarsPull
without actually downloading files from the NAS Data Portal.
"""

import os
import sys
import unittest
from unittest.mock import patch, MagicMock
import tempfile
import shutil
import subprocess
import re

# For finding the MarsPull executable
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


class TestMarsPull(unittest.TestCase):
    """Test suite for MarsPull"""
    
    @classmethod
    def setUpClass(cls):
        """Set up the test environment"""
        # Create a temporary directory for test outputs
        cls.test_dir = tempfile.mkdtemp()
        os.chdir(cls.test_dir)  # Change to test directory
        
        # Find MarsPull executable
        cls.mars_pull_exec = cls.find_mars_pull_executable()
    
    @classmethod
    def tearDownClass(cls):
        """Clean up the test environment"""
        # Change back to original directory and remove test directory
        os.chdir(SCRIPT_DIR)
        shutil.rmtree(cls.test_dir)
    
    @classmethod
    def find_mars_pull_executable(cls):
        """Find the MarsPull executable"""
        # Check common locations
        possible_paths = [
            os.path.join(SCRIPT_DIR, 'MarsPull'),
            os.path.join(SCRIPT_DIR, '..', 'bin', 'MarsPull'),
            'MarsPull',  # Check if in PATH
        ]
        
        for path in possible_paths:
            if os.path.exists(path) or shutil.which(path):
                return path
        
        # If we can't find it, use the most likely path and let the tests fail if needed
        return 'MarsPull'
    
    def run_mars_pull(self, args, expect_success=True):
        """Run MarsPull command with given arguments"""
        cmd = [self.mars_pull_exec] + args
        
        # Run the command
        try:
            result = subprocess.run(cmd, check=False, capture_output=True, text=True)
            if expect_success:
                self.assertEqual(result.returncode, 0, f"MarsPull failed: {result.stderr}")
            else:
                self.assertNotEqual(result.returncode, 0, "MarsPull should have failed but didn't")
            return result.stdout, result.stderr
        except Exception as e:
            self.fail(f"Failed to run MarsPull: {e}")
    
    @patch('requests.get')
    def test_list_files_option(self, mock_get):
        """Test the -list option"""
        # Mock the responses for directory listings
        legacy_response = MagicMock()
        legacy_response.text = """
        <html>
        <a href="https://data.nas.nasa.gov/legacygcm/legacygcmdata/ACTIVECLDS/">ACTIVECLDS</a>
        <a href="https://data.nas.nasa.gov/legacygcm/legacygcmdata/INERTCLDS/">INERTCLDS</a>
        </html>
        """
        
        activeclds_response = MagicMock()
        activeclds_response.text = """
        <html>
        <a download="fort.11_0730">fort.11_0730</a>
        <a download="fort.11_0731">fort.11_0731</a>
        </html>
        """
        
        fv3_response = MagicMock()
        fv3_response.text = """
        <html>
        <a href="https://example.com/file1.nc">file1.nc</a>
        <a href="https://example.com/file2.nc">file2.nc</a>
        </html>
        """
        
        # Configure mock to return different responses for different URLs
        def side_effect(url):
            if 'legacygcm/' in url and not 'legacygcmdata/' in url:
                return legacy_response
            elif 'legacygcmdata/' in url:
                return activeclds_response
            elif 'fv3betaout1/' in url:
                return fv3_response
            return MagicMock()
        
        mock_get.side_effect = side_effect
        
        # Run MarsPull with -list option
        stdout, stderr = self.run_mars_pull(['-list'])
        
        # Check that the command tried to fetch directory listings
        self.assertGreaterEqual(mock_get.call_count, 2)
        
        # Check that the output includes file listings
        self.assertIn("Available directories", stdout)
        self.assertIn("Available files", stdout)
    
    @patch('requests.get')
    def test_download_file_by_name(self, mock_get):
        """Test downloading a specific file by name"""
        # Mock the response for file download
        response = MagicMock()
        response.status_code = 200
        response.headers.get.return_value = '1024'  # Content length
        response.iter_content.return_value = [b'test data']
        mock_get.return_value = response
        
        # Run MarsPull to download a specific file
        stdout, stderr = self.run_mars_pull(['ACTIVECLDS', '-f', 'fort.11_0730'])
        
        # Check that the command tried to download the file
        mock_get.assert_called_with('https://data.nas.nasa.gov/legacygcm/legacygcmdata/ACTIVECLDS/fort.11_0730', stream=True)
        
        # Check that the file was created
        self.assertTrue(os.path.exists('fort.11_0730'))
        
        # Check file content
        with open('fort.11_0730', 'rb') as f:
            content = f.read()
            self.assertEqual(content, b'test data')
    
    @patch('requests.get')
    def test_download_file_by_ls(self, mock_get):
        """Test downloading files by Ls value"""
        # Mock the response for file download
        response = MagicMock()
        response.status_code = 200
        response.headers.get.return_value = '1024'  # Content length
        response.iter_content.return_value = [b'test data']
        mock_get.return_value = response
        
        # Run MarsPull to download a file by Ls
        stdout, stderr = self.run_mars_pull(['ACTIVECLDS', '-ls', '90'])
        
        # Check that the command tried to download at least one file
        self.assertGreaterEqual(mock_get.call_count, 1)
        
        # Check that at least one file was created
        files = os.listdir('.')
        self.assertGreater(len(files), 0)
        
        # Look for fort.11_* files which would be created by the Ls download
        fort_files = [f for f in files if re.match(r'fort\.11_\d+', f)]
        self.assertGreater(len(fort_files), 0)
    
    def test_missing_required_args(self):
        """Test error handling when required arguments are missing"""
        # Run MarsPull without required arguments
        stdout, stderr = self.run_mars_pull([], expect_success=False)
        
        # Check for usage/error message
        self.assertTrue('usage:' in stdout.lower() or 'error:' in stdout.lower())
        
        # Run MarsPull with directory but no -ls or -f
        stdout, stderr = self.run_mars_pull(['ACTIVECLDS'], expect_success=False)
        
        # Check for specific error message
        self.assertIn('ERROR No file requested', stdout)
    
    def test_help_message(self):
        """Test that help message displays correctly"""
        stdout, stderr = self.run_mars_pull(['-h'])
        
        # Check for typical help message components
        self.assertIn('Usage:', stdout)
        self.assertIn('-ls', stdout)
        self.assertIn('-f', stdout)
        self.assertIn('--list', stdout)


if __name__ == "__main__":
    unittest.main()