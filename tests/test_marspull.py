import unittest
import subprocess
import os
import tempfile
from unittest.mock import patch, MagicMock
import sys

# Add the path containing MarsPull script
# Adjust this path based on your project structure
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

try:
    # Try to import the script for mocking tests
    from bin.MarsPull import download, main
    CAN_IMPORT = True
except ImportError:
    CAN_IMPORT = False
    print("Note: Could not import MarsPull directly. Some tests will be skipped.")


class TestMarsPullExecutable(unittest.TestCase):
    """Test executing MarsPull as a subprocess"""
    
    def test_help_command(self):
        """Test if the help command works correctly"""
        result = subprocess.run(["python", "bin/MarsPull", "-h"], 
                               capture_output=True, text=True)
        self.assertEqual(result.returncode, 0)
        self.assertIn("Uility for downloading NASA Ames Mars Global Climate", result.stdout)
    
    def test_list_command(self):
        """Test if the list command returns available files"""
        # This test will make actual web requests, so it might be slow
        result = subprocess.run(["python", "bin/MarsPull", "-list"], 
                               capture_output=True, text=True)
        self.assertEqual(result.returncode, 0)
        self.assertIn("Available directories", result.stdout)
        self.assertIn("Available files", result.stdout)
    
    def test_missing_required_args(self):
        """Test if error is raised when required args are missing"""
        result = subprocess.run(["python", "bin/MarsPull", "ACTIVECLDS"], 
                               capture_output=True, text=True)
        # Should fail because -ls or -f is required
        self.assertNotEqual(result.returncode, 0)


@unittest.skipIf(not CAN_IMPORT, "Cannot import MarsPull directly")
class TestMarsPullFunctions(unittest.TestCase):
    """Test MarsPull's internal functions with mocks"""
    
    @patch('bin.MarsPull.requests.get')
    def test_download_function(self, mock_get):
        """Test the download function with mocked requests"""
        # Setup mock response
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.headers.get.return_value = '1024'  # Content length
        mock_response.iter_content.return_value = [b'test data']
        mock_get.return_value = mock_response
        
        # Create a temp file for testing
        with tempfile.NamedTemporaryFile() as tmp:
            # Call the function
            download('https://example.com/test.nc', tmp.name)
            
            # Verify the function called requests.get
            mock_get.assert_called_once_with('https://example.com/test.nc', stream=True)
            
            # Verify the file was written
            with open(tmp.name, 'rb') as f:
                content = f.read()
                self.assertEqual(content, b'test data')
    
    @patch('bin.MarsPull.download')
    @patch('bin.MarsPull.argparse.ArgumentParser.parse_args')
    def test_main_with_ls_param(self, mock_args, mock_download):
        """Test main function with -ls parameter"""
        # Setup mock args
        args = MagicMock()
        args.list_files = False
        args.directory_name = 'ACTIVECLDS'
        args.ls = [90]
        args.filename = None
        mock_args.return_value = args
        
        # Call main function
        main()
        
        # Verify download was called
        self.assertTrue(mock_download.called)
        # We'd need to verify the exact URL in a more complete test


if __name__ == '__main__':
    unittest.main()
