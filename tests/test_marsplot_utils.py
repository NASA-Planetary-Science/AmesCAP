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


if __name__ == '__main__':
    unittest.main()
