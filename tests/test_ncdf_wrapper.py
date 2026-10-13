import unittest
from unittest.mock import patch

import numpy as np

from amescap.Ncdf_wrapper import Fort


class FakeVariable:
    def __init__(self, data, dimensions, long_name='', units=''):
        self.data = np.asarray(data)
        self.dimensions = dimensions
        self.long_name = long_name
        self.units = units

    def __getitem__(self, key):
        return self.data[key]


class TestFortAverage(unittest.TestCase):
    @patch('amescap.Ncdf_wrapper.Ncdf')
    @patch('amescap.Ncdf_wrapper.daily_to_average', return_value=np.array([2.5]))
    def test_average_time_axis_is_binned_from_time_coordinate(
        self, daily_average, ncdf
    ):
        fort = Fort.__new__(Fort)
        fort.path = '/tmp'
        fort.fdate = '00000'
        fort.variables = {
            'lat': FakeVariable([0.0], ('lat',)),
            'lon': FakeVariable([0.0], ('lon',)),
            'pfull': FakeVariable([1.0], ('pfull',)),
            'phalf': FakeVariable([0.0, 2.0], ('phalf',)),
            'zgrid': FakeVariable([0.0], ('zgrid',)),
            'pk': FakeVariable([0.0], ('phalf',)),
            'bk': FakeVariable([0.0], ('phalf',)),
            'time': FakeVariable(
                np.arange(10.0), ('time',), 'time', 'sols since 00000'
            ),
            'temperature': FakeVariable(
                np.ones((10, 1)), ('time', 'lat'), 'temperature', 'K'
            ),
        }

        fort.write_to_average()

        self.assertIs(daily_average.call_args_list[0].kwargs['varIN'],
                      fort.variables['time'])
        ncdf.return_value.close.assert_called_once_with()


if __name__ == '__main__':
    unittest.main()