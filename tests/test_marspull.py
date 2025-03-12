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

# Add project root to Python path
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

class TestMarsPull(unittest.TestCase):
    """Test suite for MarsPull"""
    
    @classmethod
    def setUpClass(cls):
        """Set up the test environment"""
        # Create a temporary directory for test outputs
        cls.test_dir = tempfile.mkdtemp()
        os.chdir(cls.test_dir)  # Change to test directory
        
        # Find MarsPull script
        cls.mars_pull_script = os.path.join(PROJECT_ROOT, 'bin', 'MarsPull.py')
    
    @classmethod
    def tearDownClass(cls):
        """Clean up the test environment"""
        # Change back to original directory and remove test directory
        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        shutil.rmtree(cls.test_dir)
    
    def run_mars_pull(self, args, expect_success=True):
        """Run MarsPull command with given arguments"""
        cmd = [sys.executable, self.mars_pull_script] + args
        
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
        
        # Configure mock to return different responses for different URLs
        def side_effect(url, **kwargs):
            response = MagicMock()
            if url == 'https://data.nas.nasa.gov/mcmcref/legacygcm/':
                response.text = """<html>URLs for legacygcm</html>"""
            elif url == 'https://data.nas.nasa.gov/mcmcref/fv3betaout1/':
                response.text = """<html>URLs for fv3betaout1</html>"""
            else:
                response.text = """<html>Some directory content</html>"""
            return response
        
        mock_get.side_effect = side_effect
        
        # Run MarsPull with -list option
        stdout, stderr = self.run_mars_pull(['-list'])
        
        # Check that the command tried to fetch directory listings
        self.assertGreaterEqual(mock_get.call_count, 1)
        
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
        
        def side_effect(url, **kwargs):
            # Simulate finding multiple files for Ls 90
            if 'fort.11_0090' in url or 'fort.11_0091' in url:
                return response
            return MagicMock()
        
        mock_get.side_effect = side_effect
        
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
        # Test without any arguments
        stdout, stderr = self.run_mars_pull([], expect_success=False)
        self.assertIn("Error: You must specify either -list or a directory", stdout)
        
        # Test with directory but no file specification
        stdout, stderr = self.run_mars_pull(['ACTIVECLDS'], expect_success=False)
        self.assertIn("ERROR No file requested", stdout)
    
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