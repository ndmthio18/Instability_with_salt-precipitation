import math
import unittest

from porous_salt_growth import PorousSaltModel


class PorousSaltGrowthTests(unittest.TestCase):
    def setUp(self):
        self.model = PorousSaltModel(0.2634, 2165.0, 0.43, 2.4729683670962094e-6, math.sqrt(156.0))

    def test_transport_factory_reproduces_recomputed_parameters(self):
        eps = 0.43
        model = PorousSaltModel.from_transport_parameters(
            equilibrium_mass_fraction=0.2634,
            crystal_density_kg_m3=2165.0,
            top_porosity=eps,
            effective_diffusivity_m2_s=math.sqrt(eps) * 1.3e-9,
            reaction_coefficient_m_s=0.0023,
            specific_area_m_inv=36.0 * (1.0 - eps) / 0.00058,
            beta=2.0,
        )
        self.assertAlmostEqual(model.reaction_length_m, 2.4729683670962094e-6, places=17)
        self.assertAlmostEqual(model.sqrt_da_star**2, 169.77504515830265, places=10)

    def test_prat_finite_condition_is_recovered(self):
        for h_over_ell in (0.05, 0.5, 2.0, 20.0):
            response = self.model.prat_finite_response(h_over_ell * self.model.reaction_length_m, 1.0e-4)
            x = self.model.sqrt_da_star
            expected = x * math.tanh(h_over_ell) / (1.0 + x * math.tanh(h_over_ell))
            self.assertAlmostEqual(response.top_partition, expected, places=12)

    def test_finite_layer_signed_reaction_and_mass_budget(self):
        response = self.model.finite_layer_response(
            8.0 * self.model.reaction_length_m, 8.0e-5, supply_ratio=0.60
        )
        self.assertGreater(response.internal_precipitation_kg_m2_s, 0.0)
        self.assertGreater(response.internal_dissolution_kg_m2_s, 0.0)
        closed = (response.top_solid_rate_kg_m2_s
                  + response.internal_precipitation_kg_m2_s
                  - response.internal_dissolution_kg_m2_s)
        self.assertAlmostEqual(closed, response.net_solid_rate_kg_m2_s, places=17)

    def test_thick_layer_partition_and_kinematics(self):
        response = self.model.thick_layer_response(8.0e-5, supply_ratio=0.60)
        closed = (response.top_solid_rate_kg_m2_s
                  + response.internal_precipitation_kg_m2_s
                  - response.internal_dissolution_kg_m2_s)
        self.assertAlmostEqual(closed, response.net_solid_rate_kg_m2_s, places=17)
        self.assertAlmostEqual(response.thickness_velocity_m_s,
                               0.60 * response.top_velocity_m_s, places=18)

    def test_constant_forcing_time_integration(self):
        response = self.model.thick_layer_response(8.0e-5, supply_ratio=0.60)
        final = self.model.integrate_thick_layer(
            [0.0, 50.0, 100.0], [8.0e-5] * 3, [0.60] * 3
        )[-1]
        self.assertAlmostEqual(final.top_position_m, 100.0 * response.top_velocity_m_s, places=18)
        self.assertAlmostEqual(final.base_position_m, 100.0 * response.bottom_retreat_velocity_m_s, places=18)
        self.assertAlmostEqual(final.thickness_m, 100.0 * response.thickness_velocity_m_s, places=18)
        self.assertAlmostEqual(final.cumulative_net_solid_kg_m2,
                               100.0 * response.net_solid_rate_kg_m2_s, places=17)

    def test_disconnected_thick_closure_is_rejected(self):
        disconnected = PorousSaltModel(0.2634, 2165.0, 0.43, 1.0e-6, 1.0)
        with self.assertRaises(ValueError):
            disconnected.thick_layer_response(1.0e-4, supply_ratio=0.8)


if __name__ == "__main__":
    unittest.main()
