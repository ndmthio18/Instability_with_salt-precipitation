#!/usr/bin/env python3
"""Verification tests for the P2-A surface subgrid.

Coverage maps one-to-one onto P2-A item 10:

* water activity versus NaCl concentration     -> TestWaterActivity
* vapour driving force and evaporation flux    -> TestFilmEvaporation
* Prat finite-Da and large-Da growth           -> TestPratGrowthLaw
* crust height evolution                       -> TestHeightEvolution
* water and salt mass balances                 -> TestMassBalance
* response to insufficient liquid/salt supply  -> TestSupplyLimitation

plus the structural requirements of items 6, 8 and 9:

* explicit activation gate                     -> TestActivationGate
* conservative soil/crust exchange             -> TestInterfaceReconstruction
* evaporation-interface options A and B        -> TestInterfaceComparison

All tolerances come from `config/acceptance_criteria.json`, registered before
the runs. No test validates physics; they bound implementation and
discretisation error.
"""

from __future__ import annotations

import math
import unittest

from reference_config import (  # noqa: E402
    PARAMS,
    TOL,
    reference_crust,
    reference_film,
    reference_patch,
    reference_seed,
)
from surface import coupling as cp  # noqa: E402
from surface import prat_crust as pc  # noqa: E402
from surface import thermo  # noqa: E402
from surface import water_activity as wa  # noqa: E402

C_EQ = PARAMS["equilibrium_mass_fraction"]


def _rel(a: float, b: float) -> float:
    scale = max(abs(a), abs(b), 1.0e-300)
    return abs(a - b) / scale


class TestWaterActivity(unittest.TestCase):
    """Item 10: water activity versus NaCl concentration."""

    def test_pure_water_activity_is_one(self):
        self.assertEqual(wa.water_activity(0.0), 1.0)
        self.assertEqual(wa.osmotic_coefficient(0.0), 1.0)

    def test_strictly_decreasing_in_molality(self):
        values = [wa.water_activity(m) for m in [i * 0.05 for i in range(0, 161)]]
        for lower, upper in zip(values, values[1:]):
            self.assertLess(upper, lower)

    def test_saturated_value_in_declared_band(self):
        low, high = TOL["saturated_water_activity_band"]
        a_w = wa.water_activity_from_mass_fraction(C_EQ)
        self.assertGreater(a_w, low)
        self.assertLess(a_w, high)

    def test_mole_and_mass_fraction_routes_agree(self):
        for x in (1.0e-6, 0.0036, 0.02, 0.0999):
            w = thermo.mass_fraction_from_mole_fraction(x)
            self.assertLess(
                _rel(
                    wa.water_activity_from_mole_fraction(x),
                    wa.water_activity_from_mass_fraction(w),
                ),
                1.0e-14,
            )

    def test_mass_mole_conversion_round_trip(self):
        for x in (1.0e-8, 0.0036, 0.05, 0.3):
            w = thermo.mass_fraction_from_mole_fraction(x)
            self.assertLess(_rel(thermo.mole_fraction_from_mass_fraction(w), x), 1e-14)

    def test_molality_limit_is_enforced(self):
        with self.assertRaises(ValueError):
            wa.water_activity(wa.MOLALITY_LIMIT + 1.0e-9)
        with self.assertRaises(ValueError):
            wa.water_activity(-1.0e-12)

    def test_frozen_conversion_reference_values(self):
        """The P1-R2 registered conversions must still reproduce.

        The quoted task-3 references carry 8-9 significant digits, so the
        comparison is made at that precision, not bit-for-bit.
        """
        w = thermo.mass_fraction_from_mole_fraction(0.0036)
        self.assertLess(_rel(w, 0.011584506), 1.0e-8)
        rho = thermo.BrineDensityLaw().density(w)
        self.assertLess(_rel(rho, 1006.024104), 1.0e-8)
        mixture = (1.0 - 0.0036) * thermo.M_H2O + 0.0036 * thermo.M_NACL
        self.assertLess(_rel(mixture, 0.018160808992), 1.0e-9)


class TestFilmEvaporation(unittest.TestCase):
    """Item 10: vapour driving force and evaporation flux."""

    def setUp(self):
        self.film = reference_film()

    def test_hand_computed_flux(self):
        a_w = wa.water_activity_from_mass_fraction(C_EQ)
        y_s = a_w * self.film.saturation_pressure / self.film.atmospheric_pressure
        y_ref = (
            self.film.relative_humidity
            * self.film.saturation_pressure
            / self.film.atmospheric_pressure
        )
        rho_ref = thermo.ideal_gas_density(
            self.film.temperature, self.film.atmospheric_pressure, y_ref
        )
        d_v = thermo.vapour_diffusion_coefficient(
            self.film.temperature, self.film.atmospheric_pressure
        )
        expected = (
            rho_ref
            * d_v
            / self.film.thickness
            * (
                thermo.gas_water_mass_fraction(y_s)
                - thermo.gas_water_mass_fraction(y_ref)
            )
        )
        self.assertLess(
            _rel(self.film.flux_from_mass_fraction(C_EQ), expected),
            TOL["film_flux_hand_check_relative"],
        )

    def test_zero_driving_force_gives_exactly_zero_flux(self):
        film = reference_film(relative_humidity=0.5)
        self.assertLessEqual(
            abs(film.flux_from_activity(0.5)), TOL["zero_driving_absolute_kg_m2_s"]
        )

    def test_condensation_is_signed_not_clipped(self):
        film = reference_film(relative_humidity=0.9)
        flux = film.flux_from_mass_fraction(C_EQ)
        self.assertLess(flux, 0.0)

    def test_flux_decreases_with_salinity(self):
        fluxes = [
            self.film.flux_from_mass_fraction(w)
            for w in [0.0, 0.05, 0.1, 0.2, C_EQ]
        ]
        for lower, upper in zip(fluxes, fluxes[1:]):
            self.assertLess(upper, lower)

    def test_total_rate_is_not_normalised_to_a_fixed_value(self):
        """Item 4: the total evaporation rate must change with salinity."""
        fresh = self.film.flux_from_mass_fraction(0.0)
        salty = self.film.flux_from_mass_fraction(C_EQ)
        self.assertGreater(abs(fresh - salty) / fresh, 0.1)

    def test_flux_is_inversely_proportional_to_film_thickness(self):
        thin = reference_film(thickness=PARAMS["film_thickness"] / 3.0)
        self.assertLess(
            _rel(
                thin.flux_from_mass_fraction(C_EQ),
                3.0 * self.film.flux_from_mass_fraction(C_EQ),
            ),
            1.0e-13,
        )

    def test_no_normalisation_helper_was_imported(self):
        """Item 4: the upstream weight normaliser must stay excluded."""
        self.assertFalse(hasattr(wa, "normalizedEvaporationWeights"))
        self.assertFalse(hasattr(wa, "normalized_evaporation_weights"))

    def test_invalid_physical_inputs_are_rejected(self):
        for kwargs in (
            {"temperature": -1.0},
            {"relative_humidity": 1.0},
            {"thickness": 0.0},
            {"atmospheric_pressure": 0.0},
        ):
            with self.assertRaises(ValueError):
                reference_film(**kwargs)


class TestPratGrowthLaw(unittest.TestCase):
    """Item 10: Prat finite-Da and large-Da growth."""

    def setUp(self):
        self.crust = reference_crust()
        self.j = reference_film().flux_from_mass_fraction(C_EQ)

    def test_damkohler_identity(self):
        for h in (1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2):
            da = self.crust.damkohler(h)
            da_top = self.crust.damkohler_top(h)
            self.assertLess(
                _rel(da_top * da_top / da, self.crust.damkohler_star),
                TOL["damkohler_identity_relative"],
            )

    def test_damkohler_star_is_height_independent(self):
        values = {self.crust.damkohler_star for _ in range(3)}
        self.assertEqual(len(values), 1)

    def test_three_routes_to_growth_velocity_agree(self):
        """eq (28) == eq (8) on the analytic profile == the eq (20) kinetic law."""
        for h in (1.0e-5, 1.0e-4, 1.0e-3, 5.0e-3, 2.0e-2):
            v_28 = pc.growth_velocity(self.crust, h, self.j)
            c_star = pc.top_supersaturation(self.crust, h, self.j)
            u = math.sqrt(self.crust.damkohler(h))
            grad_top = c_star * u / math.tanh(u) / h
            v_eq8 = (
                C_EQ * self.j
                - self.crust.rho
                * self.crust.crust_porosity
                * self.crust.effective_salt_diffusivity
                * grad_top
            ) / (
                self.crust.crystal_density
                * (1.0 - self.crust.crust_porosity)
                * (1.0 - C_EQ)
            )
            v_kin = (
                self.crust.beta
                * self.crust.reaction_coefficient
                * self.crust.rho
                * c_star
                / self.crust.crystal_density
            )
            self.assertLess(_rel(v_28, v_eq8), TOL["three_route_growth_relative"])
            self.assertLess(_rel(v_28, v_kin), TOL["three_route_growth_relative"])

    def test_large_da_limit_is_approached(self):
        da_star = self.crust.damkohler_star
        f_inf = pc.efficiency_large_damkohler(da_star)
        self.assertLess(
            _rel(pc.efficiency(10.0, da_star), f_inf),
            TOL["large_da_limit_relative_at_da_10"],
        )
        self.assertLess(
            _rel(pc.efficiency(1000.0, da_star), f_inf),
            TOL["large_da_limit_relative_at_da_1000"],
        )

    def test_efficiency_is_bounded_and_monotone(self):
        da_star = self.crust.damkohler_star
        previous = -1.0
        for da in [0.0, 0.1, 1.0, 10.0, 100.0, 1.0e4]:
            f = pc.efficiency(da, da_star)
            self.assertGreaterEqual(f, 0.0)
            self.assertLess(f, 1.0)
            self.assertGreater(f, previous)
            previous = f
        self.assertLess(pc.efficiency_large_damkohler(da_star), 1.0)

    def test_growth_velocity_is_proportional_to_evaporation_flux(self):
        """Large-Da form: v_s is strictly linear in j (eq 29/32)."""
        h = 1.0e-3
        v1 = pc.growth_velocity(self.crust, h, self.j, pc.GrowthLaw.LARGE_DA)
        v2 = pc.growth_velocity(self.crust, h, 2.0 * self.j, pc.GrowthLaw.LARGE_DA)
        self.assertLess(_rel(v2, 2.0 * v1), 1.0e-14)

    def test_profile_satisfies_boundary_conditions(self):
        h = 1.0e-3
        base = pc.concentration_profile(self.crust, h, self.j, 0.0)
        top = pc.concentration_profile(self.crust, h, self.j, 1.0)
        self.assertLess(abs(base - C_EQ), TOL["profile_boundary_absolute"])
        self.assertLess(
            abs(top - pc.top_mass_fraction(self.crust, h, self.j)),
            TOL["profile_boundary_absolute"],
        )

    def test_profile_satisfies_the_reduced_ode(self):
        """d2 C*/dz*2 = Da C* checked by a central difference."""
        h = 1.0e-3
        da = self.crust.damkohler(h)
        dz = 1.0e-4
        for z in (0.25, 0.5, 0.75):
            values = [
                pc.concentration_profile(self.crust, h, self.j, z + k * dz) - C_EQ
                for k in (-1, 0, 1)
            ]
            second = (values[0] - 2.0 * values[1] + values[2]) / (dz * dz)
            self.assertLess(_rel(second, da * values[1]), TOL["profile_ode_relative"])

    def test_internal_over_top_approaches_inverse_sqrt_damkohler_star(self):
        """Derived identity: internal/top -> 1/sqrt(Da*) for Da >> 1.

        A direct consequence is ``f_inf = 1/(1 + internal/top)``, i.e. the Prat
        efficiency is the top-growth share of total precipitation. Prat (2026)
        does not state this; it is derived here from the analytic profile and
        is checked numerically.
        """
        da_star = self.crust.damkohler_star
        target = 1.0 / math.sqrt(da_star)
        previous = math.inf
        for h in (1.0e-3, 1.0e-2, 1.0e-1):
            flags = pc.validity_flags(self.crust, h, self.j)
            error = abs(flags["internal_over_top_precipitation"] - target) / target
            self.assertLessEqual(error, previous)
            previous = error
        self.assertLess(previous, 1.0e-6)
        self.assertLess(
            _rel(
                pc.efficiency_large_damkohler(da_star),
                1.0 / (1.0 + target),
            ),
            1.0e-14,
        )

    def test_internal_reaction_matches_quadrature(self):
        h = 1.0e-3
        n = 20001
        total = 0.0
        for i in range(n):
            z = i / (n - 1)
            weight = 0.5 if i in (0, n - 1) else 1.0
            c_star = pc.concentration_profile(self.crust, h, self.j, z) - C_EQ
            total += weight * c_star
        total *= h / (n - 1)
        expected = (
            self.crust.specific_area
            * self.crust.rho
            * self.crust.reaction_coefficient
            * total
        )
        self.assertLess(
            _rel(pc.internal_reaction_rate(self.crust, h, self.j), expected),
            TOL["quadrature_relative"],
        )

    def test_inventories_match_quadrature(self):
        h = 1.0e-3
        n = 20001
        total = 0.0
        for i in range(n):
            z = i / (n - 1)
            weight = 0.5 if i in (0, n - 1) else 1.0
            total += weight * pc.concentration_profile(self.crust, h, self.j, z)
        total *= h / (n - 1)
        expected = self.crust.crust_porosity * self.crust.rho * total
        self.assertLess(
            _rel(pc.dissolved_salt_inventory(self.crust, h, self.j), expected),
            TOL["quadrature_relative"],
        )
        pore = self.crust.crust_porosity * self.crust.rho * h
        self.assertLess(
            _rel(
                pc.water_inventory(self.crust, h, self.j),
                pore - expected,
            ),
            TOL["quadrature_relative"],
        )

    def test_parameter_sweep_covers_both_regimes(self):
        """a_v and k_r are declared, so both Da regimes must be reachable."""
        small = reference_crust(specific_area=1.0e3, reaction_coefficient=1.0e-9)
        large = reference_crust(specific_area=1.0e7, reaction_coefficient=1.0e-5)
        self.assertLess(small.damkohler(1.0e-3), 10.0)
        self.assertGreater(large.damkohler(1.0e-3), 10.0)


class TestActivationGate(unittest.TestCase):
    """Item 6: growth only after an explicitly declared seed event."""

    def test_growth_law_refuses_zero_height(self):
        crust = reference_crust()
        with self.assertRaises(pc.CrustModelDomainExit):
            pc.growth_velocity(crust, 0.0, 1.0e-5)
        with self.assertRaises(pc.CrustModelDomainExit):
            crust.damkohler(0.0)

    def test_patch_is_inert_below_the_declared_trigger(self):
        patch = reference_patch()
        started = patch.consider_activation(0.5 * C_EQ, 0.0)
        self.assertFalse(started)
        self.assertFalse(patch.state.active)
        chain = patch.diagnose(0.5 * C_EQ)
        self.assertEqual(chain.growth_velocity, 0.0)
        self.assertEqual(chain.height, 0.0)
        self.assertGreater(chain.evaporation_flux, 0.0)

    def test_inactive_step_transfers_no_salt(self):
        patch = reference_patch()
        record = patch.step(0.5 * C_EQ, 100.0)
        self.assertFalse(record["active"])
        self.assertEqual(record["salt_flux_soil_kg_m2_s"], 0.0)
        self.assertEqual(record["crust_height_m"], 0.0)
        self.assertAlmostEqual(
            record["water_flux_soil_kg_m2_s"],
            record["evaporation_flux_kg_m2_s"],
            places=20,
        )

    def test_activation_requires_the_declared_trigger(self):
        patch = reference_patch()
        self.assertTrue(patch.consider_activation(C_EQ, 1234.0))
        self.assertTrue(patch.state.active)
        self.assertEqual(patch.state.activation_time, 1234.0)
        self.assertEqual(patch.state.eligibility_time, 1234.0)
        self.assertEqual(patch.state.height, PARAMS["seed_height"])
        self.assertFalse(patch.consider_activation(C_EQ, 2000.0))

    def test_seed_rejects_zero_height(self):
        with self.assertRaises(ValueError):
            reference_seed(seed_height=0.0)

    def test_seed_porosity_must_match_the_crust(self):
        patch = reference_patch()
        patch.seed = reference_seed(seed_porosity=0.3)
        with self.assertRaises(ValueError):
            patch.consider_activation(C_EQ, 0.0)


class TestInterfaceReconstruction(unittest.TestCase):
    """Item 8: base fluxes from the moving-interface balances."""

    def setUp(self):
        self.crust = reference_crust()
        self.j = reference_film().flux_from_mass_fraction(C_EQ)

    def test_route_difference_matches_the_closed_form_identity(self):
        """The gap is a cancellation of two large fluxes.

        It is therefore compared against the flux scale, not against itself;
        measuring a catastrophic cancellation relative to its own magnitude
        would only report floating-point noise amplification.
        """
        for h in (1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2):
            f = pc.interface_fluxes(self.crust, h, self.j)
            error = abs(f.route_difference - f.route_difference_identity)
            self.assertLess(
                error / abs(f.salt_flux),
                TOL["route_difference_identity_relative"],
            )

    def test_salt_flux_is_positive_while_solid_accumulates(self):
        """Item 8: zero soil salt flux with accumulating solid is forbidden."""
        for h in (1.0e-5, 1.0e-4, 1.0e-3):
            f = pc.interface_fluxes(self.crust, h, self.j)
            self.assertGreater(f.top_skeleton_rate, 0.0)
            self.assertGreater(f.salt_flux, 0.0)
            self.assertGreaterEqual(f.salt_flux, f.top_skeleton_rate)

    def test_water_flux_identity(self):
        """F_w = rho U(0) - F_s and F_w = j - rho_v v_s + dM_w/dt."""
        for h in (1.0e-4, 1.0e-3):
            f = pc.interface_fluxes(self.crust, h, self.j)
            self.assertLess(
                _rel(f.water_flux, f.solution_mass_flux_base - f.salt_flux), 1.0e-13
            )
            self.assertLess(
                _rel(
                    f.water_flux,
                    f.evaporation_flux + f.vapour_sweep_rate + f.new_pore_water_rate,
                ),
                1.0e-13,
            )

    def test_crust_top_does_not_equal_the_base_concentration(self):
        """Item 8: C_top == C_b must not be assumed."""
        f = pc.interface_fluxes(self.crust, 1.0e-3, self.j)
        self.assertGreater(f.top_supersaturation, 0.0)
        self.assertGreater(f.top_mass_fraction, C_EQ)

    def test_solution_mass_balance_across_the_crust(self):
        f = pc.interface_fluxes(self.crust, 1.0e-3, self.j)
        self.assertLess(
            _rel(
                f.solution_mass_flux_base,
                f.solution_mass_flux_top + f.internal_reaction,
            ),
            1.0e-14,
        )

    def test_disconnected_crust_exits_instead_of_growing(self):
        crust = reference_crust(crust_saturation=0.5)
        with self.assertRaises(pc.CrustModelDomainExit):
            pc.interface_fluxes(crust, 1.0e-3, self.j)


class TestHeightEvolution(unittest.TestCase):
    """Item 10: crust height evolution."""

    def test_large_da_height_is_exactly_linear_in_time(self):
        patch = reference_patch(law=pc.GrowthLaw.LARGE_DA, activated=True)
        h0 = patch.state.height
        v_s = pc.growth_velocity(
            patch.crust, h0, patch.film.flux_from_mass_fraction(C_EQ),
            pc.GrowthLaw.LARGE_DA,
        )
        dt, steps = 500.0, 400
        for _ in range(steps):
            patch.step(C_EQ, dt)
        self.assertLess(
            _rel(patch.state.height, h0 + v_s * dt * steps),
            TOL["analytic_linear_height_relative"],
        )

    def test_finite_da_height_grows_and_saturates_towards_the_large_da_rate(self):
        crust = reference_crust()
        j = reference_film().flux_from_mass_fraction(C_EQ)
        v_inf = pc.growth_velocity(crust, 1.0, j, pc.GrowthLaw.LARGE_DA)
        previous = 0.0
        for h in (1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2):
            v = pc.growth_velocity(crust, h, j)
            self.assertGreater(v, previous)
            self.assertLess(v, v_inf * (1.0 + 1.0e-12))
            previous = v

    def test_time_step_convergence_order(self):
        def final_height(dt: float, scheme: str) -> float:
            patch = reference_patch(activated=True)
            total, elapsed = 2.0e5, 0.0
            while elapsed < total - 1.0e-9:
                step = min(dt, total - elapsed)
                patch.step(C_EQ, step, scheme=scheme)
                elapsed += step
            return patch.state.height

        reference = final_height(12.5, "heun")
        for scheme, key in (
            ("heun", "time_step_observed_order_min_heun"),
            ("euler", "time_step_observed_order_min_euler"),
        ):
            coarse = abs(final_height(1600.0, scheme) - reference)
            fine = abs(final_height(800.0, scheme) - reference)
            order = math.log2(coarse / fine)
            self.assertGreater(order, TOL[key], f"{scheme} order {order}")

    def test_height_is_monotone_under_evaporation(self):
        patch = reference_patch(activated=True)
        heights = [patch.state.height]
        for _ in range(50):
            patch.step(C_EQ, 1000.0)
            heights.append(patch.state.height)
        for lower, upper in zip(heights, heights[1:]):
            self.assertGreater(upper, lower)


class TestMassBalance(unittest.TestCase):
    """Item 10: water and salt mass balances."""

    def _run(self, interface, law, steps=500, dt=400.0, scheme="heun"):
        patch = reference_patch(interface=interface, law=law)
        patch.consider_activation(C_EQ, 0.0)
        initial = patch.state.copy()
        for _ in range(steps):
            patch.step(C_EQ, dt, scheme=scheme)
        return patch, patch.budget(initial)

    def test_salt_and_water_budgets_close(self):
        for interface in cp.EvaporationInterface:
            for law in pc.GrowthLaw:
                with self.subTest(interface=interface, law=law):
                    _, budget = self._run(interface, law)
                    self.assertLess(
                        budget["salt_residual_relative"], TOL["salt_budget_relative"]
                    )
                    self.assertLess(
                        budget["water_residual_relative"], TOL["water_budget_relative"]
                    )

    def test_salt_entering_the_crust_equals_the_inventory_increase(self):
        patch, budget = self._run(
            cp.EvaporationInterface.SOIL_TOP, pc.GrowthLaw.FINITE_DA
        )
        self.assertGreater(budget["delta_solid_salt_kg_m2"], 0.0)
        self.assertGreater(budget["integrated_salt_from_soil_kg_m2"], 0.0)
        self.assertLess(
            _rel(
                budget["integrated_salt_from_soil_kg_m2"],
                budget["delta_solid_salt_kg_m2"]
                + budget["delta_dissolved_salt_kg_m2"],
            ),
            TOL["salt_budget_relative"],
        )

    def test_solid_inventory_exceeds_the_geometric_skeleton(self):
        """Total solid = new skeleton + internal reaction, never less.

        At the declared reference ``a_v`` and ``k_r`` the internal reaction is
        in fact the dominant term. That is a reportable property of the
        declared parameters, not a defect: the reduced solution integrates the
        volumetric reaction of equations (1)-(2) explicitly, which is option 1
        of PRAT_PHYSICS_AUDIT section 6. Its consequence for the
        constant-porosity assumption is checked separately.
        """
        patch, _ = self._run(cp.EvaporationInterface.SOIL_TOP, pc.GrowthLaw.LARGE_DA)
        geometric = (
            (1.0 - patch.crust.crust_porosity)
            * patch.crust.crystal_density
            * patch.state.height
        )
        self.assertGreater(patch.state.solid_salt, geometric)
        record = patch.state.history[-1]
        self.assertGreater(record["validity_internal_over_top_precipitation"], 0.0)

    def test_constant_porosity_assumption_is_audited_not_assumed(self):
        """Equation (2) forces porosity to change wherever R_int > 0."""
        crust = reference_crust()
        j = reference_film().flux_from_mass_fraction(C_EQ)
        h = PARAMS["seed_height"]
        t_eps = pc.porosity_consumption_timescale(crust, h, j)
        t_growth = pc.growth_timescale(crust, h, j)
        self.assertGreater(t_eps, 0.0)
        self.assertTrue(math.isfinite(t_eps))
        flags = pc.validity_flags(crust, h, j)
        self.assertIn("constant_porosity_defensible", flags)
        self.assertEqual(
            flags["constant_porosity_defensible"], t_eps / t_growth >= 10.0
        )

    def test_slow_kinetics_recover_a_skeleton_dominated_crust(self):
        """With a small a_v k_r product the top skeleton dominates again."""
        crust = reference_crust(specific_area=1.0e2, reaction_coefficient=1.0e-9)
        j = reference_film().flux_from_mass_fraction(C_EQ)
        flags = pc.validity_flags(crust, 1.0e-4, j)
        self.assertLess(flags["internal_over_top_precipitation"], 1.0)

    def test_quasi_steady_storage_gap_is_reported_and_small(self):
        patch, _ = self._run(cp.EvaporationInterface.SOIL_TOP, pc.GrowthLaw.FINITE_DA)
        record = patch.state.history[-1]
        gap = abs(record["dissolved_storage_gap_kg_m2"])
        self.assertGreater(gap, 0.0)
        self.assertLess(gap / record["dissolved_salt_kg_m2"], 1.0e-2)

    def test_budget_closes_across_the_activation_event(self):
        """The declared seed inventory must be charged to the soil.

        Item 8 forbids accumulating solid salt in the crust while the soil salt
        flux stays zero. The seed is a declared input, but its mass still has
        to come from somewhere, so it is booked as a single instantaneous
        soil-to-crust transfer and the budget must close across it.
        """
        patch = reference_patch()
        initial = patch.state.copy()
        w_start = thermo.mass_fraction_from_mole_fraction(0.0036)
        ramp, dt = 5.0e4, 250.0
        for i in range(int(2.0e5 / dt)):
            t = i * dt
            w = min(C_EQ, w_start + (C_EQ - w_start) * t / ramp)
            patch.consider_activation(w, t)
            patch.step(w, dt)
        self.assertTrue(patch.state.active)
        budget = patch.budget(initial)
        self.assertGreater(budget["seed_salt_kg_m2"], 0.0)
        self.assertGreater(budget["seed_water_kg_m2"], 0.0)
        self.assertLess(budget["salt_residual_relative"], TOL["salt_budget_relative"])
        self.assertLess(
            budget["water_residual_relative"], TOL["water_budget_relative"]
        )

    def test_budget_closes_under_time_step_refinement(self):
        for dt in (1600.0, 400.0, 100.0):
            patch = reference_patch()
            patch.consider_activation(C_EQ, 0.0)
            initial = patch.state.copy()
            for _ in range(int(2.0e5 / dt)):
                patch.step(C_EQ, dt)
            budget = patch.budget(initial)
            self.assertLess(
                budget["salt_residual_relative"], TOL["salt_budget_relative"]
            )


class TestSupplyLimitation(unittest.TestCase):
    """Item 10: response to insufficient liquid/salt supply."""

    def setUp(self):
        self.crust = reference_crust()
        self.j = reference_film().flux_from_mass_fraction(C_EQ)
        self.fluxes = pc.interface_fluxes(self.crust, 1.0e-3, self.j)

    def test_sufficient_supply_is_not_flagged(self):
        diag = pc.supply_limitation_diagnostic(
            self.fluxes,
            2.0 * self.fluxes.salt_flux,
            2.0 * self.fluxes.water_flux,
        )
        self.assertFalse(diag.salt_limited)
        self.assertFalse(diag.water_limited)

    def test_insufficient_salt_supply_is_flagged_not_clipped(self):
        reduced = 0.25 * self.fluxes.salt_flux
        diag = pc.supply_limitation_diagnostic(
            self.fluxes, reduced, 2.0 * self.fluxes.water_flux
        )
        self.assertTrue(diag.salt_limited)
        self.assertLess(_rel(diag.salt_ratio, 4.0), 1.0e-12)
        # Nothing was modified: the growth law still returns its own value.
        self.assertEqual(diag.salt_demand, self.fluxes.salt_flux)
        self.assertEqual(
            pc.growth_velocity(self.crust, 1.0e-3, self.j),
            self.fluxes.growth_velocity,
        )

    def test_shortfall_can_raise_a_declared_model_domain_exit(self):
        with self.assertRaises(pc.CrustModelDomainExit):
            pc.supply_limitation_diagnostic(
                self.fluxes,
                0.5 * self.fluxes.salt_flux,
                2.0 * self.fluxes.water_flux,
                raise_on_shortfall=True,
            )

    def test_insufficient_water_supply_is_flagged(self):
        diag = pc.supply_limitation_diagnostic(
            self.fluxes, 2.0 * self.fluxes.salt_flux, 0.5 * self.fluxes.water_flux
        )
        self.assertTrue(diag.water_limited)

    def test_loss_of_connectivity_stops_the_closure(self):
        crust = reference_crust(crust_saturation=0.8, connectivity_threshold=0.95)
        with self.assertRaises(pc.CrustModelDomainExit):
            pc.interface_fluxes(crust, 1.0e-3, self.j)


class TestInterfaceComparison(unittest.TestCase):
    """Item 9: evaporation-interface options A and B."""

    def test_option_b_reduces_the_flux(self):
        film, crust = reference_film(), reference_crust()
        comparison = cp.compare_interfaces(film, crust, 1.0e-3, C_EQ)
        self.assertLess(comparison.flux_b, comparison.flux_a)
        self.assertGreater(comparison.top_supersaturation_b, 0.0)

    def test_difference_grows_with_crust_height(self):
        film, crust = reference_film(), reference_crust()
        previous = -1.0
        for h in (1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2):
            rel = cp.compare_interfaces(
                film, crust, h, C_EQ
            ).relative_flux_difference
            self.assertGreater(rel, previous)
            previous = rel

    def test_options_coincide_as_supersaturation_vanishes(self):
        film = reference_film()
        fast = reference_crust(specific_area=1.0e9, reaction_coefficient=1.0e-3)
        comparison = cp.compare_interfaces(film, fast, 1.0e-5, C_EQ)
        self.assertLess(comparison.relative_flux_difference, 1.0e-4)

    def test_option_b_fixed_point_converges(self):
        patch = reference_patch(
            interface=cp.EvaporationInterface.CRUST_TOP, activated=True, height=1.0e-3
        )
        chain = patch.diagnose(C_EQ)
        self.assertGreater(chain.fixed_point_iterations, 0)
        self.assertLess(abs(chain.fixed_point_residual), 1.0e-16)
        reconstructed = patch.film.flux_from_mass_fraction(
            pc.top_mass_fraction(patch.crust, 1.0e-3, chain.evaporation_flux)
        )
        self.assertLess(_rel(reconstructed, chain.evaporation_flux), 1.0e-9)

    def test_options_coincide_before_activation(self):
        a = reference_patch(interface=cp.EvaporationInterface.SOIL_TOP)
        b = reference_patch(interface=cp.EvaporationInterface.CRUST_TOP)
        w = 0.5 * C_EQ
        self.assertEqual(
            a.diagnose(w).evaporation_flux, b.diagnose(w).evaporation_flux
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
