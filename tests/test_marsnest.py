#!/usr/bin/env python3
"""
Integration tests for MarsNest using small synthetic grid_spec files
"""
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np
from netCDF4 import Dataset


def write_grid_spec(path, lon_range, lat_range, n=9):
    """Write a grid_spec file with a regular (n+1) x (n+1) corner grid"""
    lon, lat = np.meshgrid(np.linspace(*lon_range, n + 1),
                           np.linspace(*lat_range, n + 1))
    with Dataset(path, 'w') as f:
        f.createDimension('nyp', n + 1)
        f.createDimension('nxp', n + 1)
        f.createVariable('grid_lon', 'f8', ('nyp', 'nxp'))[:] = lon
        f.createVariable('grid_lat', 'f8', ('nyp', 'nxp'))[:] = lat


class TestMarsNest(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp(prefix='MarsNest_test_')
        self.project_root = os.path.dirname(
            os.path.dirname(os.path.abspath(__file__)))

    def tearDown(self):
        shutil.rmtree(self.test_dir, ignore_errors=True)

    def run_mars_nest(self, args):
        return subprocess.run(
            [sys.executable, os.path.join(self.project_root, 'bin',
                                          'MarsNest.py')] + args,
            capture_output=True, text=True, cwd=self.test_dir,
            stdin=subprocess.DEVNULL, timeout=300,
        )

    def test_layout_pdf_for_parent_and_nest(self):
        history = os.path.join(self.test_dir, 'history')
        os.mkdir(history)
        write_grid_spec(os.path.join(history, 'grid_spec.tile1.nc'),
                        (0., 90.), (-45., 45.))
        write_grid_spec(os.path.join(history, 'grid_spec.nest02.tile7.nc'),
                        (30., 60.), (-15., 15.))

        result = self.run_mars_nest([history])

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        output = os.path.join(history, 'nest_layout.pdf')
        self.assertTrue(os.path.exists(output))
        self.assertGreater(os.path.getsize(output), 0)

    def test_missing_grid_spec_files_is_an_error(self):
        result = self.run_mars_nest([self.test_dir])

        self.assertNotEqual(result.returncode, 0)
        self.assertIn('No grid_spec files found', result.stdout)

    def test_help_documents_usage(self):
        result = self.run_mars_nest(['-h'])

        self.assertEqual(result.returncode, 0)
        self.assertIn('input_path', result.stdout)
        self.assertIn('-topo', result.stdout)


if __name__ == '__main__':
    unittest.main()
