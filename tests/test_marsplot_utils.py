import os
import sys
import unittest

import numpy as np

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)
_argv = sys.argv
try:
    # MarsPlot parses the command line at import time
    sys.argv = ['MarsPlot']
    from bin import MarsPlot
finally:
    sys.argv = _argv


class FakeVariable:
    def __init__(self, data):
        self.data = np.ma.asarray(data)

    def __getitem__(self, key):
        return self.data[key]


class TestAreoByTime(unittest.TestCase):
    def test_average_file_areo(self):
        areo = FakeVariable(np.array([[10.0], [20.0], [30.0]]))
        np.testing.assert_array_equal(MarsPlot.areo_by_time(areo),
                                      [10.0, 20.0, 30.0])

    def test_diurn_file_uses_midnight_value(self):
        areo = FakeVariable(np.array([[[10.0], [10.5]], [[20.0], [20.5]]]))
        np.testing.assert_array_equal(MarsPlot.areo_by_time(areo),
                                      [10.0, 20.0])

    def test_single_time_step_diurn_file_keeps_time_axis(self):
        # Regression: np.squeeze dropped the time axis, which broke
        # plots of single-record diurn files (e.g., after MarsFiles -split)
        areo = FakeVariable(np.array([[[268.1], [268.4]]]))
        ls = MarsPlot.areo_by_time(areo)
        self.assertEqual(ls.shape, (1,))
        self.assertEqual(float(ls[0]), 268.1)


class TestGetTimeIndex(unittest.TestCase):
    def test_single_time_step(self):
        ti, txt = MarsPlot.get_time_index(270, np.array([268.1]))
        self.assertEqual(int(ti), 0)
        self.assertIn('268.10', txt)

    def test_scalar_time_step(self):
        ti, txt = MarsPlot.get_time_index(270, np.float64(268.1))
        self.assertEqual(int(ti), 0)
        self.assertIn('268.10', txt)


class TestTopographyForOverlay(unittest.TestCase):
    lon = np.array([45.0, 135.0, 225.0, 315.0])
    zsurf = np.array([[1.0, 2.0, 3.0, 4.0],
                      [5.0, 6.0, 7.0, 8.0]])

    def test_matching_grid_is_shifted_like_the_data(self):
        # Regression: topography was never drawn on lon X lat plots
        _, var = MarsPlot.shift_data(self.lon, self.zsurf * 10.)
        topo = MarsPlot.topography_for_overlay(self.lon, self.zsurf, var)
        _, expected = MarsPlot.shift_data(self.lon, self.zsurf)
        np.testing.assert_array_equal(topo, expected)

    def test_no_fixed_file(self):
        var = np.zeros((2, 4))
        self.assertIsNone(
            MarsPlot.topography_for_overlay(self.lon, None, var))

    def test_different_grid_is_skipped(self):
        # e.g., data regridded to a different resolution
        var = np.zeros((3, 4))
        self.assertIsNone(
            MarsPlot.topography_for_overlay(self.lon, self.zsurf, var))


class TestLongitudeRange(unittest.TestCase):
    lons = np.arange(1.875, 360., 3.75)    # FV3-based MGCM grid

    def test_range_across_prime_meridian(self):
        # Regression: -45,45 averaged the 270 deg on the far side
        loni, txt = MarsPlot.get_lon_index(np.array([-45., 45.]), self.lons)
        self.assertEqual(len(loni), 25)
        self.assertAlmostEqual(float(self.lons[loni[0]]), 313.125)
        self.assertAlmostEqual(float(self.lons[loni[-1]]), 43.125)
        self.assertIn('-46.9<->43.1', txt)

    def test_range_across_dateline_keeps_direction(self):
        loni, txt = MarsPlot.get_lon_index(np.array([160., -40.]), self.lons)
        self.assertEqual(len(loni), 44)
        self.assertIn('159.4<->-39.4', txt)

    def test_range_within_eastern_hemisphere(self):
        loni, _ = MarsPlot.get_lon_index(np.array([10., 50.]), self.lons)
        self.assertEqual(len(loni), 12)


if __name__ == '__main__':
    unittest.main()
