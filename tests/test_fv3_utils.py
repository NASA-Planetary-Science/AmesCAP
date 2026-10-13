import unittest

import numpy as np

from amescap.FV3_utils import (
    R_CO2, R_REF_ATMOS, fms_Z_calc, interp_KDTree, lon180_to_360,
    lon360_to_180, mass_stream, press_to_alt_atmosphere_Mars,
    ref_atmosphere_Mars_PTD,
    shiftgrid_180_to_360, shiftgrid_360_to_180
)


class TestCO2GasConstant(unittest.TestCase):
    psfc = np.array([600.0])
    ak = np.array([0.0, 0.0])
    bk = np.array([0.5, 1.0])
    temperature = np.array([[200.0]])

    def test_shared_co2_constant_value(self):
        self.assertEqual(R_CO2, 189.0)

    def test_half_level_altitude_matches_hypsometric_equation(self):
        # Isothermal layer between 300 Pa and 600 Pa:
        # dz = R T / g * ln(p_bottom / p_top)
        altitude = fms_Z_calc(self.psfc, self.ak, self.bk, self.temperature,
                              lev_type='half')
        expected_top = R_CO2 * 200.0 / 3.72 * np.log(600.0 / 300.0)

        self.assertAlmostEqual(float(altitude[0, 0]), expected_top, places=8)
        self.assertAlmostEqual(float(altitude[1, 0]), 0.0, places=8)

    def test_default_gas_constant_is_shared_constant(self):
        default = fms_Z_calc(self.psfc, self.ak, self.bk, self.temperature)
        explicit = fms_Z_calc(self.psfc, self.ak, self.bk, self.temperature,
                              rgas=R_CO2)

        np.testing.assert_allclose(default, explicit)


class TestReferenceAtmosphere(unittest.TestCase):
    def test_reference_profile_is_continuous_at_segment_boundaries(self):
        # The segment boundary pressures were derived with R_REF_ATMOS;
        # using a different gas constant opens ~10 % jumps here
        for boundary in (57000.0, 110000.0, 120000.0):
            below, _, _ = ref_atmosphere_Mars_PTD(boundary - 1.0)
            above, _, _ = ref_atmosphere_Mars_PTD(boundary + 1.0)
            self.assertLess(abs(float(below) / float(above) - 1.0), 1e-3,
                            f'discontinuity at {boundary} m')

    def test_density_uses_reference_gas_constant(self):
        pressure, temperature, density = ref_atmosphere_Mars_PTD(0.0)
        self.assertAlmostEqual(float(pressure), 610.0, places=6)
        self.assertAlmostEqual(
            float(density), 610.0 / (R_REF_ATMOS * float(temperature)),
            places=12)

    def test_pressure_to_altitude_inverts_profile(self):
        altitudes = np.array([1000.0, 30000.0, 80000.0, 115000.0])
        pressures, _, _ = ref_atmosphere_Mars_PTD(altitudes)
        np.testing.assert_allclose(press_to_alt_atmosphere_Mars(pressures),
                                   altitudes, rtol=1e-6)


class TestLongitudeConversion(unittest.TestCase):
    legacy_lon = np.arange(-180., 180., 6.)    # Legacy MGCM grid
    fv3_lon = np.arange(1.875, 360., 3.75)    # FV3-based MGCM grid

    def test_legacy_grid_to_360_is_increasing(self):
        # Regression: -180 became 180 but stayed first, which garbled
        # MarsPlot maps of Legacy data
        lon = lon180_to_360(self.legacy_lon)
        self.assertTrue(np.all(np.diff(lon) > 0))
        np.testing.assert_array_equal(lon[29:32], [174., 180., 186.])

    def test_legacy_data_follows_its_longitudes(self):
        # Store each column's own longitude as the data value
        data = np.vstack([self.legacy_lon, self.legacy_lon])
        shifted = shiftgrid_180_to_360(self.legacy_lon, data)
        np.testing.assert_array_equal(np.mod(shifted[0], 360.),
                                      lon180_to_360(self.legacy_lon))

    def test_round_trip_to_180(self):
        lon = lon360_to_180(self.fv3_lon)
        self.assertTrue(np.all(np.diff(lon) > 0))
        data = shiftgrid_360_to_180(self.fv3_lon, self.fv3_lon[None, :])
        np.testing.assert_array_equal(np.where(data[0] > 180,
                                               data[0] - 360, data[0]), lon)

    def test_grids_already_in_target_range_are_unchanged(self):
        np.testing.assert_array_equal(lon180_to_360(self.fv3_lon),
                                      self.fv3_lon)
        np.testing.assert_array_equal(lon360_to_180(self.legacy_lon),
                                      self.legacy_lon)

    def test_masked_data_stays_masked(self):
        data = np.ma.masked_array(self.legacy_lon[None, :],
                                  mask=(self.legacy_lon == 0)[None, :])
        shifted = shiftgrid_180_to_360(self.legacy_lon, data)
        self.assertTrue(np.ma.isMaskedArray(shifted))
        self.assertTrue(shifted.mask[0, 0])    # 0 deg is now first


class TestInterpKDTree(unittest.TestCase):
    lat = np.array([-45., 0., 45.])
    lon = np.array([0., 90., 180., 270.])
    var = np.arange(12.).reshape(3, 4)

    def test_same_grid_returns_input(self):
        out = interp_KDTree(self.var, self.lat, self.lon, self.lat, self.lon)
        self.assertFalse(np.isnan(out).any())
        np.testing.assert_allclose(out, self.var)

    def test_coincident_points_use_source_value(self):
        # Mix coincident and in-between target points, with extra
        # leading dimensions
        var = np.stack([self.var, 2*self.var])
        lat_out = np.array([0., 22.5])
        lon_out = np.array([90., 135.])
        out = interp_KDTree(var, self.lat, self.lon, lat_out, lon_out)
        self.assertEqual(out.shape, (2, 2, 2))
        self.assertFalse(np.isnan(out).any())
        np.testing.assert_allclose(out[:, 0, 0], [5., 10.])


class TestMassStream(unittest.TestCase):
    # Constant northward wind on 51 levels
    lat = np.array([-30., 0., 30.])
    H, psfc, g, a = 8000., 700., 3.72, 3400.e3

    def msf(self, level, ztype):
        v_avg = np.ones((len(level), len(self.lat)))
        return mass_stream(v_avg, self.lat, level, type=ztype,
                           psfc=self.psfc, H=self.H, factor=1.,
                           g=self.g, a=self.a)

    def segments(self, Z):
        # Trapezoidal integral of exp(-Z/H) over each layer
        f = np.exp(-Z/self.H)
        return 0.5 * np.diff(Z) * (f[1:] + f[:-1])

    def scale(self):
        return (2*np.pi*self.a*self.psfc/(self.g*self.H)
                * np.cos(np.deg2rad(self.lat))[None, :])

    def test_pstd_integrates_every_layer(self):
        # Pressure levels from the surface up; MSF is integrated
        # downward from the top
        pstd = self.psfc * np.exp(-np.linspace(0., 50000., 51)/self.H)
        seg = self.segments(self.H * np.log(self.psfc/pstd))
        I = np.append(np.cumsum(seg[::-1])[::-1], 0.)
        msf = self.msf(pstd, "pstd")
        np.testing.assert_allclose(msf, I[:, None]*self.scale(),
                                   rtol=1e-12)
        self.assertTrue((msf[0, :] > 0).all())
        np.testing.assert_array_equal(msf[-1, :], 0.)

    def test_zagl_integrates_every_layer(self):
        # Altitude levels from the surface up; MSF is integrated
        # upward from the surface
        zagl = np.linspace(0., 50000., 51)
        I = -np.append(0., np.cumsum(self.segments(zagl)))
        msf = self.msf(zagl, "zagl")
        np.testing.assert_allclose(msf, I[:, None]*self.scale(),
                                   rtol=1e-12)
        np.testing.assert_array_equal(msf[0, :], 0.)
        self.assertTrue((msf[-1, :] < 0).all())

if __name__ == '__main__':
    unittest.main()
