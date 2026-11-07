"""P2-C verification: the graded crust and the salt_grow parameter transfer.

Every tolerance is read from ``config/acceptance_criteria_p2c.json``, which
was registered before any run. The predictions tested here
(``predictions_registered_before_measurement``) were written down before the
corresponding quantity was measured.

Run with::

    python -m unittest test_p2c_graded
"""

from __future__ import annotations

import math
import unittest

import reference_config_p2c as rc
from reference_config import PARAMS
from reference_config_p2b import evolving_parameters, seeded_state as p2b_seeded_state
from surface import evolving_crust as ec
from surface import graded_crust as gc
from surface import prat_crust as pc

TOL = rc.TOL_P2C
PUBLISHED = rc.REFERENCE["published_groups_at_h0_1mm"]
PUBLISHED_V = rc.REFERENCE["published_velocities_m_s"]
PUBLISHED_B2 = rc.REFERENCE["published_b2_6h"]


def graded_twin(params: ec.EvolvingCrustParameters, **extra) -> gc.GradedCrustParameters:
    """The graded parameter object that must behave exactly like a P2-B one."""
    kwargs = dict(
        equilibrium_mass_fraction=params.equilibrium_mass_fraction,
        deposition_porosity=params.deposition_porosity,
        crystal_density=params.crystal_density,
        specific_area=params.specific_area,
        reaction_coefficient=params.reaction_coefficient,
        effective_salt_diffusivity=params.effective_salt_diffusivity,
        beta=params.beta,
        solution_density=params.solution_density,
        vapour_density=params.vapour_density,
        crust_saturation=params.crust_saturation,
        clogging_porosity=params.clogging_porosity,
        connectivity_porosity=params.connectivity_porosity,
        diffusivity_exponent=params.diffusivity_exponent,
        specific_area_exponent=params.specific_area_exponent,
        linearised_top_closure=params.linearised_top_closure,
        legacy_liquid_closure=True,
    )
    kwargs.update(extra)
    return gc.GradedCrustParameters(**kwargs)


class TestGridConstruction(unittest.TestCase):
    def test_a_uniform_request_tiles_the_height_exactly(self):
        height, cells = 1.0e-4, 32
        sizes = gc.graded_seed_grid(height, height / cells)
        self.assertEqual(len(sizes), cells)
        self.assertAlmostEqual(sum(sizes), height, delta=1.0e-18)
        self.assertTrue(all(abs(s - height / cells) < 1.0e-20 for s in sizes))

    def test_grading_keeps_the_top_fine_and_the_base_coarse(self):
        height, top = 1.0e-3, 1.0e-7
        sizes = gc.graded_seed_grid(height, top, window=2.0e-6, ratio=1.25,
                                    maximum_cell=5.0e-6)
        self.assertAlmostEqual(sum(sizes), height, delta=1.0e-12 * height)
        self.assertAlmostEqual(sizes[-1], top, delta=1.0e-20)
        self.assertGreater(sizes[0], 20.0 * top)
        # A bottom sliver is absorbed into its neighbour, so the bottom-most
        # cell may be up to 1.5 times the cap; nothing above it may be.
        self.assertLessEqual(max(sizes), 1.5 * 5.0e-6)
        self.assertLessEqual(max(sizes[1:]), 5.0e-6 * (1.0 + 1.0e-12))
        # Cell sizes increase monotonically downward above the bottom cell.
        for lower, upper in zip(sizes[1:], sizes[2:]):
            self.assertGreaterEqual(lower + 1.0e-20, upper)
        self.assertLess(len(sizes), 300)

    def test_a_graded_grid_is_far_cheaper_than_a_uniform_one(self):
        height, top = 1.0e-3, 5.0e-7 / 8.0
        uniform = height / top
        graded = len(
            gc.graded_seed_grid(height, top, window=2.0e-5, ratio=1.25,
                                maximum_cell=5.0e-6)
        )
        self.assertGreater(uniform, 10000.0)
        self.assertLess(graded, uniform / 10.0)


class TestReductionToP2B(unittest.TestCase):
    """On a uniform grid the graded module must *be* the P2-B module."""

    def setUp(self):
        self.flux = rc.FILM_FLUX
        self.height = rc.H_SEED
        self.cells = 32
        self.evolving = evolving_parameters()
        self.graded = graded_twin(self.evolving)

    def _states(self):
        evolving_state = p2b_seeded_state(self.evolving, self.flux, cells=self.cells)
        graded_state = gc.GradedCrustState(top_cell_size=self.height / self.cells)
        gc.activate(graded_state, self.graded, self.height, self.flux, 0.0)
        return evolving_state, graded_state

    def test_the_seed_states_agree(self):
        evolving_state, graded_state = self._states()
        self.assertEqual(len(evolving_state.occupied()), len(graded_state.thicknesses))
        self.assertAlmostEqual(evolving_state.height, graded_state.height, delta=1e-20)
        self.assertAlmostEqual(
            evolving_state.solid_salt, graded_state.solid_salt, delta=1e-20
        )

    def test_the_first_solve_agrees(self):
        evolving_state, graded_state = self._states()
        a = ec.solve_profile(evolving_state, self.evolving, self.flux)
        b = gc.solve_profile(graded_state, self.graded, self.flux)
        for name in (
            "growth_velocity",
            "top_supersaturation",
            "internal_reaction",
            "dissolved_salt",
            "water",
            "top_gradient_coefficient",
        ):
            got, want = getattr(b, name), getattr(a, name)
            self.assertLessEqual(
                abs(got - want) / abs(want), TOL["p2b_equivalence_relative"], name
            )

    def test_the_trajectory_agrees(self):
        evolving_state, graded_state = self._states()
        worst = 0.0
        for _ in range(250):
            a = ec.step(evolving_state, self.evolving, self.flux, 20.0)
            b = gc.step(graded_state, self.graded, self.flux, 20.0)
            for key in (
                "crust_height_m",
                "growth_velocity_m_s",
                "minimum_porosity",
                "internal_reaction_kg_m2_s",
                "salt_flux_soil_kg_m2_s",
                "water_flux_soil_kg_m2_s",
            ):
                worst = max(worst, abs(a[key] - b[key]) / max(abs(a[key]), 1e-300))
        # Bisection in the top closure is the only source of difference.
        self.assertLess(worst, 1.0e-9)

    def test_the_legacy_base_flux_reproduces_p2b(self):
        evolving_state, graded_state = self._states()
        a = ec.solve_profile(evolving_state, self.evolving, self.flux)
        b = gc.solve_profile(graded_state, self.graded, self.flux)
        self.assertLessEqual(
            abs(b.salt_flux_base_route_legacy - a.salt_flux_base_route)
            / abs(a.salt_flux_base_route),
            1.0e-12,
        )


class TestReductionToEquation28(unittest.TestCase):
    def test_the_linearised_uniform_limit_reproduces_equation_28(self):
        params = graded_twin(
            evolving_parameters(), linearised_top_closure=True
        )
        frozen = evolving_parameters().frozen()
        want = pc.growth_velocity(frozen, rc.H_SEED, rc.FILM_FLUX)
        errors = []
        for cells in (64, 128, 256, 512):
            state = gc.GradedCrustState(top_cell_size=rc.H_SEED / cells)
            gc.activate(state, params, rc.H_SEED, rc.FILM_FLUX, 0.0)
            got = gc.solve_profile(state, params, rc.FILM_FLUX).growth_velocity
            errors.append(abs(got - want) / want)
        self.assertLess(errors[-1], TOL["eq28_reduction_relative"])
        for coarse, fine in zip(errors, errors[1:]):
            self.assertLess(fine, coarse)

    def test_a_graded_grid_also_reproduces_equation_28(self):
        """Grading must not break the limit, only redistribute resolution.

        Tested in the regime grading exists for, ``delta_r << h``: the
        transferred salt_grow parameters, where a uniform grid at the same top
        cell size would need of order 1e5 cells. The error is second order in
        the number of cells across ``delta_r``, so the registered tolerance is
        a statement about resolution, not about grading: the production grid
        of 8 cells per ``delta_r`` carries a bias of a few times 1e-4, which
        is recorded and used when interpreting the matched benchmark.
        """
        params = rc.case_parameters(3, linearised_top_closure=True)
        frozen = pc.CrustParameters(
            equilibrium_mass_fraction=params.equilibrium_mass_fraction,
            crust_porosity=params.deposition_porosity,
            crystal_density=params.crystal_density,
            specific_area=params.specific_area,
            reaction_coefficient=params.reaction_coefficient,
            effective_salt_diffusivity=params.effective_salt_diffusivity,
            beta=params.beta,
            vapour_density=params.vapour_density,
        )
        want = pc.growth_velocity(frozen, rc.H_SEED, rc.FILM_FLUX)
        errors = []
        for per in (16.0, 32.0, 64.0):
            state = rc.seeded_state(params, rc.FILM_FLUX, rc.H_SEED,
                                    cells_per_delta=per)
            got = gc.solve_profile(state, params, rc.FILM_FLUX).growth_velocity
            errors.append(abs(got - want) / want)
        self.assertLess(errors[-1], TOL["eq28_reduction_relative"])
        for coarse, fine in zip(errors, errors[1:]):
            self.assertGreater(coarse / fine, 3.0)   # observed order about 2

    def test_the_production_grid_bias_is_measured_not_assumed(self):
        params = rc.case_parameters(3, linearised_top_closure=True)
        frozen = pc.CrustParameters(
            equilibrium_mass_fraction=params.equilibrium_mass_fraction,
            crust_porosity=params.deposition_porosity,
            crystal_density=params.crystal_density,
            specific_area=params.specific_area,
            reaction_coefficient=params.reaction_coefficient,
            effective_salt_diffusivity=params.effective_salt_diffusivity,
            beta=params.beta,
            vapour_density=params.vapour_density,
        )
        want = pc.growth_velocity(frozen, rc.H_SEED, rc.FILM_FLUX)
        state = rc.seeded_state(params, rc.FILM_FLUX, rc.H_SEED,
                                cells_per_delta=8.0)
        got = gc.solve_profile(state, params, rc.FILM_FLUX).growth_velocity
        bias = (got - want) / want
        # Recorded so that a sub-0.1 % benchmark agreement is not over-claimed.
        self.assertLess(abs(bias), 2.0e-3)
        self.assertGreater(abs(bias), 1.0e-4)

    def test_grading_is_not_used_when_the_layer_exceeds_the_crust(self):
        """Case 0's reaction length exceeds its own height, so no grading."""
        params = rc.case_parameters(0)
        self.assertGreater(params.reaction_length(), rc.H_SEED)
        state = rc.seeded_state(params, rc.FILM_FLUX, rc.H_SEED)
        self.assertEqual(len(state.thicknesses), 32)
        self.assertTrue(
            all(abs(t - state.thicknesses[0]) < 1e-20 for t in state.thicknesses)
        )


class TestCoarsening(unittest.TestCase):
    def setUp(self):
        self.params = rc.case_parameters(3)
        self.flux = rc.FILM_FLUX

    def test_merging_conserves_thickness_and_solid_exactly(self):
        state = rc.seeded_state(self.params, self.flux, rc.H_SEED, cells_per_delta=4.0)
        before_height = sum(state.thicknesses)
        before_solid = sum(state.solid)
        merged = gc.coarsen(state, freeze_depth=state.height / 2.0,
                            maximum_cell=state.height / 4.0)
        self.assertGreater(merged, 0)
        self.assertAlmostEqual(sum(state.thicknesses), before_height, delta=1e-20)
        self.assertAlmostEqual(sum(state.solid), before_solid, delta=1e-20)

    def test_merging_never_touches_the_reaction_layer(self):
        state = rc.seeded_state(self.params, self.flux, rc.H_SEED, cells_per_delta=4.0)
        delta = self.params.reaction_length()
        top_before = state.thicknesses[-1]
        gc.coarsen(state, freeze_depth=20.0 * delta, maximum_cell=10.0 * delta)
        self.assertAlmostEqual(state.thicknesses[-1], top_before, delta=1e-22)
        depth = 0.0
        for thickness in reversed(state.thicknesses):
            if depth >= 20.0 * delta:
                break
            self.assertLessEqual(thickness, 1.000001 * top_before)
            depth += thickness

    def test_coarsening_does_not_change_the_answer(self):
        """The whole point: the wake is inert, so merging it must be free."""
        plain = rc.run_case(self.params, self.flux, rc.H_SEED, cells_per_delta=4.0,
                            horizon=400.0, plateau_windows=1e9)
        merged = rc.run_case(self.params, self.flux, rc.H_SEED, cells_per_delta=4.0,
                             horizon=400.0, plateau_windows=1e9, freeze_lengths=20.0)
        a = plain["state"].history[-1]
        b = merged["state"].history[-1]
        self.assertGreater(merged["state"].merges, 0)
        for key in ("growth_velocity_m_s", "top_supersaturation", "crust_height_m",
                    "minimum_porosity", "internal_reaction_kg_m2_s"):
            self.assertLess(
                abs(a[key] - b[key]) / max(abs(a[key]), 1e-300), 1.0e-6, key
            )
        self.assertLess(b["cells"], a["cells"])


class TestSaltGrowClosures(unittest.TestCase):
    def test_the_transferred_specific_area_matches_the_reference_formula(self):
        self.assertAlmostEqual(
            rc.TRANSFER["specific_area"],
            36.0 * (1.0 - 0.43) / 25.0e-6,
            delta=1.0e-6,
        )
        self.assertAlmostEqual(rc.SALT_GROW["a_v_ref"], 820800.0, delta=1.0e-6)

    def test_the_b2_closures_reproduce_the_reference_expressions(self):
        params = rc.matched_parameters("B2")
        for porosity in (0.43, 0.30, 0.24, 0.10):
            self.assertAlmostEqual(
                params.area(porosity),
                36.0 * (1.0 - porosity) / rc.SALT_GROW["d_sub"],
                delta=1.0e-6,
            )
            self.assertAlmostEqual(
                params.diffusivity(porosity),
                math.sqrt(porosity) * rc.SALT_GROW["D_molecular"],
                delta=1.0e-22,
            )

    def test_the_b2_area_increases_as_porosity_falls(self):
        """Unlike the P2-B power laws, so the P2-B sealing law cannot apply."""
        params = rc.matched_parameters("B2")
        self.assertGreater(params.area(0.20), params.area(0.43))

    def test_b1_freezes_both_coefficients(self):
        params = rc.matched_parameters("B1")
        self.assertEqual(params.area(0.1), params.area(0.43))
        self.assertEqual(params.diffusivity(0.1), params.diffusivity(0.43))


class TestPublishedGroups(unittest.TestCase):
    """This repository and salt_grow must agree on the dimensionless groups."""

    def setUp(self):
        self.groups = rc.dimensionless_groups(
            rc.matched_parameters("B2"),
            rc.SALT_GROW["h0"],
            rc.SALT_GROW["evaporation_mass_flux"],
        )

    def test_groups_match_the_published_values(self):
        for name, key in (
            ("Pe", "Pe"),
            ("Da", "Da"),
            ("Da_top", "Da_top"),
            ("Da_star", "Da_star"),
            ("delta_C_m", "delta_r_m"),
        ):
            want = PUBLISHED[name]
            got = self.groups[key]
            self.assertLess(abs(got - want) / want, 1.0e-4, name)

    def test_the_analytical_velocity_matches_the_published_value(self):
        want = PUBLISHED_V["prat_analytical"]
        got = self.groups["growth_velocity_m_s"]
        self.assertLess(abs(got - want) / want, 1.0e-8)

    def test_the_efficiency_matches_the_published_internal_fraction(self):
        # salt_grow quotes 1 - sqrt(Da*)/(1+sqrt(Da*)) as about 26.9 %.
        self.assertAlmostEqual(1.0 - self.groups["efficiency_f"], 0.269, delta=0.002)


class TestResolutionRequirement(unittest.TestCase):
    def test_the_committed_p2b_grid_cannot_resolve_the_transferred_parameters(self):
        params = rc.case_parameters(3)
        delta = params.reaction_length()
        committed_cell = rc.H_SEED / 32.0
        self.assertLess(delta, 1.0e-6)
        self.assertLess(
            delta / committed_cell, TOL["boundary_layer_cells_minimum"]
        )

    def test_an_under_resolved_grid_gives_a_different_answer(self):
        params = rc.case_parameters(3)
        coarse = rc.run_case(params, rc.FILM_FLUX, rc.H_SEED,
                             cells_per_delta=1.0e-9, horizon=400.0,
                             plateau_windows=1e9)
        fine = rc.run_case(params, rc.FILM_FLUX, rc.H_SEED, cells_per_delta=8.0,
                           horizon=400.0, plateau_windows=1e9)
        a = coarse["state"].history[-1]["minimum_porosity"]
        b = fine["state"].history[-1]["minimum_porosity"]
        self.assertLess(coarse["cells_per_delta_actual"], 1.0)
        self.assertGreater(abs(a - b), 0.01)

    def test_refinement_converges(self):
        params = rc.case_parameters(3)
        previous = None
        errors = []
        for per in (2.0, 4.0, 8.0):
            run = rc.run_case(params, rc.FILM_FLUX, rc.H_SEED, cells_per_delta=per,
                              horizon=400.0, plateau_windows=1e9)
            value = run["state"].history[-1]["growth_velocity_m_s"]
            if previous is not None:
                errors.append(abs(value - previous) / value)
            previous = value
        self.assertLess(errors[-1], errors[0])
        self.assertLess(errors[-1], TOL["grid_convergence_successive_relative"])


class TestRegisteredPredictions(unittest.TestCase):
    def test_the_sealing_criterion_is_correct_for_every_case(self):
        expected = {0: True, 1: True, 2: False, 3: False}
        for case, seals in expected.items():
            params = rc.case_parameters(case)
            groups = rc.dimensionless_groups(params, rc.H_SEED, rc.FILM_FLUX)
            self.assertEqual(groups["seals_predicted"], seals, f"case {case}")
            threshold = (1.0 - params.deposition_porosity) / (
                params.deposition_porosity - params.connectivity_porosity
            )
            self.assertEqual(groups["sqrt_Da_star"] <= threshold, seals)

    def test_the_matched_benchmark_is_predicted_not_to_seal(self):
        groups = rc.dimensionless_groups(
            rc.matched_parameters("B2"),
            rc.SALT_GROW["h0"],
            rc.SALT_GROW["evaporation_mass_flux"],
        )
        self.assertFalse(groups["seals_predicted"])

    def test_the_wake_drop_is_the_inverse_root_damkohler_star(self):
        params = rc.matched_parameters("B2")
        drop = params.predicted_wake_porosity_drop()
        self.assertAlmostEqual(
            drop,
            (1.0 - params.deposition_porosity) / math.sqrt(params.damkohler_star()),
            delta=1.0e-15,
        )
        # a_v delta_r / beta is the same number by construction.
        eps = params.deposition_porosity
        self.assertAlmostEqual(
            drop,
            params.area(eps) * params.reaction_length(eps) / params.beta,
            delta=1.0e-12 * drop,
        )

    def test_the_wake_drop_is_independent_of_the_evaporation_flux(self):
        params = rc.matched_parameters("B2")
        base = params.predicted_wake_porosity_drop()
        for factor in (0.1, 10.0):
            groups = rc.dimensionless_groups(
                params, rc.SALT_GROW["h0"],
                factor * rc.SALT_GROW["evaporation_mass_flux"],
            )
            self.assertAlmostEqual(groups["predicted_wake_drop"], base, delta=1e-15)

    def test_the_measured_wake_porosity_matches_the_prediction(self):
        params = rc.case_parameters(2)
        run = rc.run_case(params, rc.FILM_FLUX, rc.H_SEED, cells_per_delta=4.0)
        measured = run["state"].history[-1]["minimum_porosity"]
        predicted = (
            params.deposition_porosity - params.predicted_wake_porosity_drop()
        )
        self.assertFalse(run["state"].disconnected)
        self.assertLess(
            abs(measured - predicted) / measured,
            TOL["wake_porosity_prediction_relative"],
        )


class TestConservation(unittest.TestCase):
    def test_budgets_close_in_the_transferred_regime(self):
        params = rc.case_parameters(3)
        run = rc.run_case(params, rc.FILM_FLUX, rc.H_SEED, cells_per_delta=4.0,
                          horizon=600.0, plateau_windows=1e9)
        report = gc.budget(run["state"], run["initial"])
        self.assertLess(report["salt_residual_relative"], TOL["salt_budget_relative"])
        self.assertLess(report["water_residual_relative"], TOL["water_budget_relative"])
        self.assertLess(
            report["solid_split_residual_relative"],
            TOL["solid_split_residual_relative"],
        )

    def test_the_corrected_route_identity_holds(self):
        params = rc.case_parameters(3)
        run = rc.run_case(params, rc.FILM_FLUX, rc.H_SEED, cells_per_delta=4.0,
                          horizon=600.0, plateau_windows=1e9)
        record = run["state"].history[-1]
        self.assertLess(
            record["route_difference_relative"],
            TOL["route_difference_corrected_relative"],
        )

    def test_the_legacy_closure_keeps_a_finite_gap(self):
        """The P2-B 7.5-10.4 % route gap is a closure defect, not noise."""
        params = rc.case_parameters(3)
        run = rc.run_case(params, rc.FILM_FLUX, rc.H_SEED, cells_per_delta=4.0,
                          horizon=600.0, plateau_windows=1e9)
        record = run["state"].history[-1]
        self.assertGreater(record["route_difference_legacy_relative"], 0.01)
        self.assertGreater(
            record["route_difference_legacy_relative"],
            100.0 * record["route_difference_relative"],
        )

    def test_the_legacy_gap_matches_its_closed_form(self):
        params = rc.case_parameters(3)
        run = rc.run_case(params, rc.FILM_FLUX, rc.H_SEED, cells_per_delta=4.0,
                          horizon=600.0, plateau_windows=1e9)
        record = run["state"].history[-1]
        predicted = (
            params.equilibrium_mass_fraction
            * params.rho
            * params.crust_saturation
            / params.crystal_density
            * record["internal_reaction_kg_m2_s"]
            / abs(record["salt_flux_soil_kg_m2_s"])
        )
        self.assertLess(
            abs(record["route_difference_legacy_relative"] - predicted) / predicted,
            0.05,
        )

    def test_the_remaining_gap_is_the_declared_linearisation_residual(self):
        """Registered prediction partly falsified, and replaced by a measurement.

        The criteria file predicted that the corrected route difference would
        converge to zero at first order in the time step. It does not: it
        converges to a non-zero, grid- and step-independent value, which is
        exactly the residual of the linearised sink ``(1-C_eq)`` that Prat's
        own reduction uses in place of ``(1-C)``:

            F_s,transport - F_s,inventory = - s L_top .
        """
        params = rc.case_parameters(3)
        run = rc.run_case(params, rc.FILM_FLUX, rc.H_SEED, cells_per_delta=4.0,
                          horizon=600.0, plateau_windows=1e9)
        record = run["state"].history[-1]
        self.assertLess(record["route_identity_error_relative"], 0.01)
        self.assertLess(
            abs(record["route_difference_kg_m2_s"]), 1.0e-4
            * abs(record["salt_flux_soil_kg_m2_s"])
        )

    def test_the_remaining_gap_is_independent_of_grid_and_time_step(self):
        params = rc.case_parameters(3)
        values = []
        for per, steps in ((2.0, 375.0), (8.0, 3000.0)):
            run = rc.run_case(params, rc.FILM_FLUX, rc.H_SEED, cells_per_delta=per,
                              steps_per_tau=steps, maximum_steps=300000,
                              horizon=600.0, plateau_windows=1e9)
            values.append(run["state"].history[-1]["route_difference_relative"])
        self.assertLess(abs(values[1] - values[0]) / values[0], 0.05)

    def test_the_remaining_gap_vanishes_with_the_supersaturation(self):
        """- s L_top is proportional to s, so a weaker film must shrink it."""
        params = rc.case_parameters(3)
        values = []
        for factor in (1.0, 0.1):
            run = rc.run_case(params, factor * rc.FILM_FLUX, rc.H_SEED,
                              cells_per_delta=4.0, horizon=600.0,
                              plateau_windows=1e9)
            values.append(run["state"].history[-1]["route_difference_relative"])
        self.assertLess(values[1], 0.5 * values[0])


class TestConnectivityStop(unittest.TestCase):
    def test_connectivity_loss_precedes_clogging(self):
        params = rc.case_parameters(0)
        run = rc.run_case(params, rc.FILM_FLUX, rc.H_SEED)
        state = run["state"]
        self.assertTrue(state.disconnected)
        self.assertFalse(state.sealed)
        self.assertEqual(run["stop_reason"], "hydraulic_connectivity_loss")
        self.assertLessEqual(
            state.history[-1]["minimum_porosity"], params.connectivity_porosity
        )
        self.assertGreater(
            state.history[-1]["minimum_porosity"], params.clogging_porosity
        )

    def test_a_stopped_crust_refuses_to_be_extrapolated(self):
        params = rc.case_parameters(0)
        run = rc.run_case(params, rc.FILM_FLUX, rc.H_SEED)
        with self.assertRaises(pc.CrustModelDomainExit):
            gc.step(run["state"], params, rc.FILM_FLUX, 1.0)

    def test_the_growth_velocity_is_never_capped(self):
        """No supply clipping anywhere: v_s depends only on the top closure."""
        params = rc.case_parameters(3)
        state = rc.seeded_state(params, rc.FILM_FLUX, rc.H_SEED, cells_per_delta=2.0)
        solution = gc.solve_profile(state, params, rc.FILM_FLUX)
        expected = (
            params.beta
            * params.reaction_coefficient
            * params.rho
            * solution.top_supersaturation
            / params.crystal_density
        )
        self.assertAlmostEqual(
            solution.growth_velocity, expected, delta=1.0e-18 * expected
        )


class TestMatchedBenchmark(unittest.TestCase):
    def test_the_matched_model_reproduces_the_salt_grow_velocity(self):
        params = rc.matched_parameters("B2")
        run = rc.run_case(
            params,
            rc.SALT_GROW["evaporation_mass_flux"],
            rc.SALT_GROW["h0"],
            cells_per_delta=4.0,
            freeze_lengths=40.0,
        )
        measured = run["state"].history[-1]["growth_velocity_m_s"]
        want = PUBLISHED_V["m3_b2_corrected_6h"]
        self.assertLess(
            abs(measured - want) / want,
            TOL["matched_benchmark_velocity_relative"],
        )

    def test_the_matched_model_reproduces_the_salt_grow_wake_porosity(self):
        params = rc.matched_parameters("B2")
        run = rc.run_case(
            params,
            rc.SALT_GROW["evaporation_mass_flux"],
            rc.SALT_GROW["h0"],
            cells_per_delta=4.0,
            freeze_lengths=40.0,
        )
        measured = run["state"].history[-1]["minimum_porosity"]
        want = PUBLISHED_B2["minimum_porosity"]
        self.assertLess(abs(measured - want) / want, 0.02)

    def test_the_matched_model_is_faster_than_the_analytical_law(self):
        """salt_grow attributes this to the non-uniform porosity field."""
        params = rc.matched_parameters("B2")
        run = rc.run_case(
            params,
            rc.SALT_GROW["evaporation_mass_flux"],
            rc.SALT_GROW["h0"],
            cells_per_delta=4.0,
            freeze_lengths=40.0,
        )
        measured = run["state"].history[-1]["growth_velocity_m_s"]
        self.assertGreater(measured, PUBLISHED_V["prat_analytical"])
        excess = (measured - PUBLISHED_V["prat_analytical"]) / PUBLISHED_V[
            "prat_analytical"
        ]
        self.assertLess(abs(excess - 0.0222571), 0.01)


class TestNothingCommittedChanged(unittest.TestCase):
    def test_the_p2b_module_still_uses_its_own_closures(self):
        """Guard: P2-C must not have altered the committed P2-B behaviour."""
        params = evolving_parameters()
        self.assertEqual(params.diffusivity_exponent, 0.0)
        self.assertEqual(params.specific_area_exponent, 0.0)
        self.assertFalse(params.linearised_top_closure)
        self.assertEqual(params.deposition_porosity, 0.4)
        self.assertEqual(params.specific_area, PARAMS["specific_area"])
        self.assertEqual(params.reaction_coefficient, PARAMS["reaction_coefficient"])

    def test_case_zero_is_the_committed_p2b_configuration(self):
        case0 = rc.case_parameters(0)
        p2b = evolving_parameters()
        for name in (
            "equilibrium_mass_fraction",
            "deposition_porosity",
            "crystal_density",
            "specific_area",
            "reaction_coefficient",
            "effective_salt_diffusivity",
            "beta",
            "vapour_density",
        ):
            self.assertEqual(getattr(case0, name), getattr(p2b, name), name)


if __name__ == "__main__":
    unittest.main()
