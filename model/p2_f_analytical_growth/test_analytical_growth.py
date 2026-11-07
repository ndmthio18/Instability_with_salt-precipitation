import math
import pathlib
import unittest

import analytical_growth as ag


class AnalyticalGrowthTests(unittest.TestCase):
    def setUp(self):
        self.p = ag.Parameters()
        self.h = 1.0e-5
        self.j = 1.8e-5

    def test_solubility_conversion_is_explicit(self):
        self.assertAlmostEqual(self.p.equilibrium_mass_fraction,
                               0.2649408697676898, delta=1e-15)

    def test_zero_evaporation_means_zero_growth(self):
        self.assertEqual(ag.growth_velocity(self.p, self.h, 0.0), 0.0)

    def test_growth_is_nonnegative(self):
        for j in (0.0, 1e-8, 1e-6, self.j):
            self.assertGreaterEqual(ag.growth_velocity(self.p, self.h, j), 0.0)

    def test_zero_internal_reaction_limit(self):
        # In Prat's finite-Da factor, sqrt(Da* Da) == Da_top.  Letting a_v -> 0
        # gives tanh(sqrt(Da))/sqrt(Da) -> 1, hence Da_top/(1+Da_top).
        target = self.p.damkohler_top(self.h)
        for area in (1e3, 1e1, 1e-1, 1e-3):
            da = (self.h**2 * (1-self.p.equilibrium_mass_fraction) * area
                  * self.p.reaction_coefficient
                  / (self.p.porosity*self.p.effective_diffusivity))
            da_star = (self.p.beta**2 * (1-self.p.porosity)**2
                       * (1-self.p.equilibrium_mass_fraction)
                       * self.p.reaction_coefficient
                       / (self.p.porosity*self.p.effective_diffusivity*area))
            s = math.sqrt(da_star) * math.tanh(math.sqrt(da))
            self.assertAlmostEqual(math.sqrt(da*da_star), target, delta=1e-12)
            if area == 1e-3:
                self.assertAlmostEqual(s/(1+s), target/(1+target), delta=2e-8)

    def test_independent_bvp_reduces_at_small_peclet(self):
        j = 1e-9
        simple = ag.growth_velocity(self.p, self.h, j)
        independent = ag.independent_advective_velocity(self.p, self.h, j)
        self.assertLess(abs(simple-independent)/simple, 1e-4)

    def test_advection_is_not_silently_ignored_when_large(self):
        p = ag.Parameters(maximum_solution_peclet=1e-8)
        state = ag.activation_trial(p, 1e-3)
        trial = ag.evaluate_trial(state, p, 1e-3, 1.0)
        flags = ag.validity(p, state, trial, p.equilibrium_mass_fraction, 1.0)
        self.assertFalse(flags["solution_advection_negligible"])
        self.assertNotAlmostEqual(
            ag.independent_advective_velocity(p, state.height, 1e-3),
            trial.velocity, delta=1e-14)

    def test_large_da_top_reference(self):
        h = 1e3 / self.p.height_coefficient
        v = ag.growth_velocity(self.p, h, self.j)
        vinf = ag.large_damkohler_top_velocity(self.p, self.j)
        self.assertAlmostEqual(v/vinf, 1000/1001, delta=1e-15)

    def test_constant_porosity_and_no_internal_precipitation(self):
        source = pathlib.Path(ag.__file__).read_text(encoding="utf-8")
        self.assertNotIn("solveProfile", source)
        self.assertNotIn("specific_area", source)
        self.assertEqual(self.p.porosity, ag.Parameters().porosity)

    def test_inventory_fluxes_match_state_change(self):
        state = ag.activation_trial(self.p, self.h)
        before_salt = state.total_salt
        before_water = state.water
        trial = ag.evaluate_trial(state, self.p, self.j, 2.0)
        self.assertAlmostEqual(
            trial.salt_flux*2.0,
            trial.proposed.total_salt-before_salt, delta=1e-18)
        self.assertAlmostEqual(
            (trial.water_flux-self.j)*2.0,
            trial.proposed.water-before_water, delta=1e-18)

    def test_transactional_lifecycle(self):
        state = ag.activation_trial(self.p, self.h)
        a = ag.evaluate_trial(state, self.p, self.j, 1.0)
        b = ag.evaluate_trial(state, self.p, self.j, 1.0)
        self.assertEqual(a.proposed, b.proposed)
        ag.rollback(state, a)
        with self.assertRaises(ag.ModelDomainExit):
            ag.commit(state, a)
        ag.commit(state, b)
        with self.assertRaises(ag.ModelDomainExit):
            ag.commit(state, b)

    def test_euler_converges_to_independent_height_reference(self):
        t_end = 2000.0
        exact = ag.analytical_height(self.p, self.h, self.j, t_end)
        errors = []
        for dt in (20.0, 10.0, 5.0):
            state = ag.activation_trial(self.p, self.h)
            for _ in range(round(t_end/dt)):
                trial = ag.evaluate_trial(state, self.p, self.j, dt)
                ag.commit(state, trial)
            errors.append(abs(state.height-exact))
        self.assertGreater(errors[0]/errors[1], 1.9)
        self.assertGreater(errors[1]/errors[2], 1.9)

    def test_supply_flags_do_not_clip_growth(self):
        state = ag.activation_trial(self.p, self.h)
        trial = ag.evaluate_trial(state, self.p, self.j, 1.0)
        invalid = ag.validity(self.p, state, trial,
                              0.95*self.p.equilibrium_mass_fraction, 0.5)
        self.assertFalse(invalid["valid"])
        self.assertGreater(trial.velocity, 0.0)


if __name__ == "__main__":
    unittest.main()
