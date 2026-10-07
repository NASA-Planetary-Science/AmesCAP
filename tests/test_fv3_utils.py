import unittest

import numpy as np

from amescap.FV3_utils import (
    R_CO2, R_REF_ATMOS, fms_Z_calc, press_to_alt_atmosphere_Mars,
    ref_atmosphere_Mars_PTD
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


if __name__ == '__main__':
    unittest.main()
