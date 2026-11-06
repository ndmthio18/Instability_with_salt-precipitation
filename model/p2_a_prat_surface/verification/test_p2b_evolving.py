#!/usr/bin/env python3
"""Verification tests for the P2-B evolving-porosity crust.

Coverage:

* reduction to Prat equation (28) in the frozen-uniform limit
                                              -> TestReductionToClosedForm
* variable-coefficient transport solve        -> TestTransportSolve
* porosity evolution from equation (2)        -> TestPorosityEvolution
* water and salt mass balances, solid split   -> TestMassBalance
* sealing and the growth-before-sealing law   -> TestSealing
* declared a_v(eps) and D(eps) sensitivities  -> TestDeclaredClosures

Tolerances come from `config/acceptance_criteria_p2b.json`, registered before
the runs, together with three predictions registered before measurement.
"""

from __future__ import annotations

import math
import unittest

from reference_config import PARAMS, reference_film
from reference_config_p2b import (
    P2B,
    TOL_P2B,
    evolving_parameters,
    run_to_sealing,
    seeded_state,
    snapshot,
)
from run_p2b import predicted_growth
from surface import evolving_crust as ec
from surface import prat_crust as pc

C_EQ = PARAMS["equilibrium_mass_fraction"]
H_SEED = PARAMS["seed_height"]


def _rel(a: float, b: float) -> float:
    return abs(a - b) / max(abs(a), abs(b), 1.0e-300)


def _flux() -> float:
    return reference_film().flux_from_mass_fraction(C_EQ)


class TestReductionToClosedForm(unittest.TestCase):
    """The generalised solver must reproduce the audited closed form."""

    def setUp(self):
        self.j = _flux()
        # rho_v = 0 so the comparison is against equation (28) exactly as the
        # audit states it; the vapour sweep is checked separately.
        self.frozen = evolving_parameters(vapour_density=0.0).frozen()

    def _velocity(self, cells: int, linearised: bool) -> float:
        params = evolving_parameters(
            vapour_density=0.0, linearised_top_closure=linearised
        )
        state = seeded_state(params, self.j, cells=cells)
        return ec.solve_profile(state, params, self.j).growth_velocity

    def test_linearised_closure_converges_to_equation_28(self):
        reference = pc.growth_velocity(self.frozen, H_SEED, self.j)
        errors = []
        for cells in (16, 32, 64, 128, 256, 512):
            errors.append(_rel(self._velocity(cells, True), reference))
        for coarse, fine in zip(errors, errors[1:]):
            order = math.log2(coarse / fine)
            self.assertGreater(order, TOL_P2B["eq28_reduction_observed_order_min"])
        self.assertLess(errors[-1], TOL_P2B["eq28_reduction_relative_at_n512"])

    def test_exact_closure_offset_equals_the_top_linearisation(self):
        """Without the C ~ C_eq step the answer shifts by C*_top / C_eq."""
        reference = pc.growth_velocity(self.frozen, H_SEED, self.j)
        measured = _rel(self._velocity(512, False), reference)
        predicted = pc.top_supersaturation(self.frozen, H_SEED, self.j) / C_EQ
        self.assertLess(
            _rel(measured, predicted),
            TOL_P2B["exact_closure_offset_predicted_relative"],
        )

    def test_frozen_parameters_round_trip(self):
        params = evolving_parameters()
        frozen = params.frozen()
        self.assertEqual(frozen.crust_porosity, params.deposition_porosity)
        self.assertEqual(frozen.specific_area, params.specific_area)
        self.assertEqual(frozen.equilibrium_mass_fraction, C_EQ)


class TestTransportSolve(unittest.TestCase):
    """The variable-coefficient two-point solve itself."""

    def setUp(self):
        self.params = evolving_parameters()
        self.j = _flux()
        self.state = seeded_state(self.params, self.j)

    def test_profile_is_monotone_and_bounded(self):
        solution = ec.solve_profile(self.state, self.params, self.j)
        self.assertGreater(solution.profile[0], 0.0)
        for lower, upper in zip(solution.profile, solution.profile[1:]):
            self.assertGreater(upper, lower)
        self.assertLess(solution.profile[-1], solution.top_supersaturation)

    def test_solution_mass_balance_across_the_crust(self):
        solution = ec.solve_profile(self.state, self.params, self.j)
        self.assertLess(
            _rel(
                solution.solution_mass_flux_base,
                solution.solution_mass_flux_top + solution.internal_reaction,
            ),
            1.0e-14,
        )

    def test_inactive_or_zero_height_is_refused(self):
        empty = ec.EvolvingCrustState(cell_size=H_SEED / 32)
        with self.assertRaises(pc.CrustModelDomainExit):
            ec.solve_profile(empty, self.params, self.j)
        with self.assertRaises(pc.CrustModelDomainExit):
            ec.step(empty, self.params, self.j, 100.0)

    def test_zero_evaporation_gives_no_supersaturation(self):
        solution = ec.solve_profile(self.state, self.params, 0.0)
        self.assertEqual(solution.top_supersaturation, 0.0)
        self.assertEqual(solution.growth_velocity, 0.0)
        self.assertEqual(solution.internal_reaction, 0.0)

    def test_lower_porosity_raises_top_supersaturation(self):
        """A tighter crust impedes the diffusive back-flux."""
        loose = ec.solve_profile(self.state, self.params, self.j)
        tight = seeded_state(self.params, self.j)
        for i in range(len(tight.solid)):
            tight.solid[i] = (1.0 - 0.1) * tight.thickness(i)
        tightened = ec.solve_profile(tight, self.params, self.j)
        self.assertGreater(tightened.top_supersaturation, loose.top_supersaturation)
        self.assertGreater(tightened.growth_velocity, loose.growth_velocity)


class TestPorosityEvolution(unittest.TestCase):
    """Equation (2) applied pointwise on the material coordinate z."""

    def setUp(self):
        self.params = evolving_parameters()
        self.j = _flux()

    def test_porosity_decreases_everywhere_under_evaporation(self):
        state = seeded_state(self.params, self.j)
        previous = state.porosity_profile()
        for _ in range(20):
            ec.step(state, self.params, self.j, 200.0)
            current = state.porosity_profile()
            for index in range(len(previous)):
                self.assertLessEqual(current[index], previous[index] + 1.0e-18)
            previous = current

    def test_new_material_enters_at_the_deposition_porosity(self):
        state = seeded_state(self.params, self.j)
        for _ in range(50):
            ec.step(state, self.params, self.j, 200.0)
        profile = state.porosity_profile()
        self.assertLess(profile[-1], self.params.deposition_porosity)
        self.assertGreater(profile[-1], profile[len(profile) // 2])

    def test_porosity_rate_matches_equation_2(self):
        state = seeded_state(self.params, self.j)
        solution = ec.solve_profile(state, self.params, self.j)
        before = state.porosity_profile()
        dt = 10.0
        ec.step(state, self.params, self.j, dt)
        after = state.porosity_profile()
        for index in range(len(before) - 1):
            expected = (
                self.params.area(solution.porosity[index])
                * self.params.rho
                * self.params.reaction_coefficient
                * solution.profile[index]
                / self.params.crystal_density
            )
            measured = (before[index] - after[index]) / dt
            self.assertLess(_rel(measured, expected), 1.0e-6)

    def test_the_porosity_minimum_sits_at_the_original_seed_top(self):
        """Registered prediction: the oldest ever-top material clogs first."""
        state = seeded_state(self.params, self.j, cells=128)
        for _ in range(2000):
            ec.step(state, self.params, self.j, 100.0)
            if state.sealed:
                break
        profile = state.porosity_profile()
        index = min(range(len(profile)), key=lambda i: profile[i])
        self.assertLess(
            abs(state.centre(index) - H_SEED) / H_SEED,
            TOL_P2B["clogging_depth_relative"],
        )


class TestMassBalance(unittest.TestCase):
    """Budgets must still close once porosity evolves."""

    def test_salt_and_water_budgets_close(self):
        params = evolving_parameters()
        state, initial, _ = run_to_sealing(params, _flux(), dt=200.0)
        budget = ec.budget(state, initial)
        self.assertLess(
            budget["salt_residual_relative"], TOL_P2B["salt_budget_relative"]
        )
        self.assertLess(
            budget["water_residual_relative"], TOL_P2B["water_budget_relative"]
        )

    def test_solid_splits_exactly_into_top_and_internal(self):
        """Now legitimate: z carries the deposition-age information."""
        params = evolving_parameters()
        state, initial, _ = run_to_sealing(params, _flux(), dt=200.0)
        budget = ec.budget(state, initial)
        self.assertGreater(budget["internal_reaction_kg_m2"], 0.0)
        self.assertGreater(budget["top_skeleton_kg_m2"], 0.0)
        self.assertLess(
            budget["solid_split_residual_relative"],
            TOL_P2B["solid_split_residual_relative"],
        )

    def test_budget_closes_under_time_step_refinement(self):
        params = evolving_parameters()
        for dt in (400.0, 200.0, 100.0):
            state, initial, _ = run_to_sealing(params, _flux(), dt=dt)
            budget = ec.budget(state, initial)
            self.assertLess(
                budget["salt_residual_relative"], TOL_P2B["salt_budget_relative"]
            )

    def test_salt_flux_is_positive_while_solid_accumulates(self):
        params = evolving_parameters()
        state = seeded_state(params, _flux())
        for _ in range(30):
            record = ec.step(state, params, _flux(), 200.0)
            self.assertGreater(record["salt_flux_soil_kg_m2_s"], 0.0)
            self.assertGreater(record["top_skeleton_kg_m2_s"], 0.0)

    def test_seed_inventory_is_charged_to_the_soil(self):
        params = evolving_parameters()
        state = seeded_state(params, _flux())
        self.assertGreater(state.seed_salt, 0.0)
        self.assertGreater(state.seed_water, 0.0)
        self.assertLess(
            _rel(state.cumulative_salt_from_soil, state.seed_salt), 1.0e-15
        )


class TestSealing(unittest.TestCase):
    """Sealing behaviour and the registered growth-before-sealing law."""

    def test_the_crust_seals_at_the_declared_reference(self):
        params = evolving_parameters()
        state, _, _ = run_to_sealing(params, _flux(), dt=200.0)
        self.assertTrue(state.sealed)
        self.assertIsNotNone(state.sealed_time)
        self.assertLessEqual(state.minimum_porosity, params.clogging_porosity)

    def test_stepping_a_sealed_crust_is_refused(self):
        params = evolving_parameters()
        state, _, _ = run_to_sealing(params, _flux(), dt=200.0)
        with self.assertRaises(pc.CrustModelDomainExit):
            ec.step(state, params, _flux(), 100.0)

    def test_growth_before_sealing_matches_beta_eps0_over_av(self):
        """Registered prediction, checked in the thin-crust regime."""
        low, high = TOL_P2B["sealing_law_ratio_band_thin_crust"]
        for a_v in (1.0e5, 1.0e6):
            params = evolving_parameters(specific_area=a_v)
            state, _, _ = run_to_sealing(params, _flux())
            growth = state.height - H_SEED
            predicted = params.beta * params.deposition_porosity / a_v
            ratio = growth / predicted
            self.assertGreater(ratio, low, f"a_v={a_v} ratio={ratio}")
            self.assertLess(ratio, high, f"a_v={a_v} ratio={ratio}")

    def test_growth_before_sealing_is_independent_of_evaporation(self):
        params = evolving_parameters()
        heights = []
        for humidity in (0.2, 0.4, 0.6):
            flux = reference_film(relative_humidity=humidity).flux_from_mass_fraction(
                C_EQ
            )
            state, _, _ = run_to_sealing(params, flux)
            heights.append(state.height - H_SEED)
        for value in heights[1:]:
            self.assertLess(
                _rel(value, heights[0]),
                TOL_P2B["sealing_law_evaporation_independence_relative"],
            )

    def test_growth_before_sealing_is_nearly_independent_of_kinetics(self):
        heights = []
        for k_r in (1.0e-8, 1.0e-7, 1.0e-6):
            params = evolving_parameters(reaction_coefficient=k_r)
            state, _, _ = run_to_sealing(params, _flux())
            heights.append(state.height - H_SEED)
        for value in heights[1:]:
            self.assertLess(
                _rel(value, heights[0]),
                TOL_P2B["sealing_law_kinetics_independence_relative"],
            )

    def test_sealing_time_converges_under_grid_refinement(self):
        params = evolving_parameters()
        times = []
        for cells in (32, 64, 128):
            state, _, _ = run_to_sealing(params, _flux(), dt=200.0, cells=cells)
            times.append(state.sealed_time)
        self.assertLess(
            _rel(times[-1], times[-2]), TOL_P2B["sealing_time_grid_relative"]
        )

    def test_sealing_time_is_insensitive_to_the_time_step(self):
        params = evolving_parameters()
        times = []
        for dt in (400.0, 200.0, 100.0):
            state, _, _ = run_to_sealing(params, _flux(), dt=dt)
            times.append(state.sealed_time)
        self.assertLess(
            _rel(times[-1], times[0]), TOL_P2B["sealing_time_timestep_relative"]
        )


class TestGeneralSealingLaw(unittest.TestCase):
    """The a_v(eps) power-law generalisation, registered before measurement."""

    def test_prediction_holds_for_every_declared_exponent(self):
        slack = TOL_P2B["sealing_law_lower_bound_slack"]
        for exponent in (0.0, 1.0, 2.0):
            params = evolving_parameters(specific_area_exponent=exponent)
            state, _, _ = run_to_sealing(
                params, _flux(), dt=200.0, maximum_steps=20000
            )
            self.assertTrue(state.sealed, f"m={exponent} did not seal")
            growth = state.height - H_SEED
            predicted = predicted_growth(params, params.clogging_porosity)
            ratio = growth / predicted
            # Registered as a lower bound: deeper material sees a smaller C*
            # and seals later, so the thin-crust estimate under-predicts.
            self.assertGreater(ratio, 1.0 - slack, f"m={exponent} ratio={ratio}")

    def test_prediction_is_tight_in_the_thin_crust_limit(self):
        params = evolving_parameters()
        state, _, _ = run_to_sealing(params, _flux(), dt=200.0)
        growth = state.height - H_SEED
        predicted = predicted_growth(params, params.clogging_porosity)
        self.assertLess(
            abs(growth / predicted - 1.0),
            TOL_P2B["sealing_law_thin_crust_relative"],
        )

    def test_prediction_is_independent_of_evaporation_for_every_exponent(self):
        for exponent in (0.0, 1.0):
            params = evolving_parameters(specific_area_exponent=exponent)
            heights = []
            for humidity in (0.2, 0.6):
                state, _, _ = run_to_sealing(
                    params,
                    reference_film(
                        relative_humidity=humidity
                    ).flux_from_mass_fraction(C_EQ),
                    maximum_steps=20000,
                )
                heights.append(state.height - H_SEED)
            self.assertLess(
                _rel(heights[0], heights[1]),
                TOL_P2B["sealing_law_evaporation_independence_relative"],
            )

    def test_under_resolved_cases_are_identified_and_improve_on_refinement(self):
        """The grid must resolve h/sqrt(Da), not the growth increment."""
        params = evolving_parameters(
            specific_area=1.0e7, reaction_coefficient=1.0e-6
        )
        cells_default = int(P2B["cells_per_seed_height"])
        layer = H_SEED / math.sqrt(params.frozen().damkohler(H_SEED))
        coarse_cells = layer / (H_SEED / cells_default)
        self.assertLess(coarse_cells, TOL_P2B["boundary_layer_cells_minimum"])

        predicted = predicted_growth(params, params.clogging_porosity)
        errors = []
        for cells in (cells_default, 8 * cells_default):
            state, _, _ = run_to_sealing(
                params, _flux(), cells=cells, maximum_steps=20000
            )
            errors.append(abs((state.height - H_SEED) / predicted - 1.0))
        self.assertLess(errors[1], errors[0])


class TestDeclaredClosures(unittest.TestCase):
    """The optional a_v(eps) and D(eps) closures, default off."""

    def test_defaults_are_constant(self):
        params = evolving_parameters()
        self.assertEqual(params.area(0.01), params.specific_area)
        self.assertEqual(
            params.diffusivity(0.01), params.effective_salt_diffusivity
        )

    def test_shrinking_specific_area_delays_sealing(self):
        baseline, _, _ = run_to_sealing(evolving_parameters(), _flux(), dt=200.0)
        delayed, _, _ = run_to_sealing(
            evolving_parameters(specific_area_exponent=2.0),
            _flux(),
            dt=200.0,
            maximum_steps=4000,
        )
        self.assertGreater(delayed.height - H_SEED, baseline.height - H_SEED)

    def test_shrinking_diffusivity_raises_the_growth_velocity(self):
        """A tighter crust traps supersaturation at the top."""
        flux = _flux()
        plain = evolving_parameters()
        tortuous = evolving_parameters(diffusivity_exponent=1.5)
        state_a = seeded_state(plain, flux)
        state_b = seeded_state(tortuous, flux)
        for _ in range(100):
            ec.step(state_a, plain, flux, 200.0)
            ec.step(state_b, tortuous, flux, 200.0)
        self.assertGreater(
            state_b.history[-1]["growth_velocity_m_s"],
            state_a.history[-1]["growth_velocity_m_s"],
        )

    def test_exponents_do_not_change_the_initial_state(self):
        flux = _flux()
        plain = ec.solve_profile(
            seeded_state(evolving_parameters(), flux), evolving_parameters(), flux
        )
        tuned = evolving_parameters(
            diffusivity_exponent=1.5, specific_area_exponent=2.0
        )
        shaped = ec.solve_profile(seeded_state(tuned, flux), tuned, flux)
        self.assertLess(
            _rel(plain.growth_velocity, shaped.growth_velocity), 1.0e-12
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
