"""P2-D verification: the 1D Richards column, the film boundary and the
soil/crust interface.

Tolerances come from `config/acceptance_criteria_p2d.json`, registered before
the runs. The predictions tested here were written down before the
corresponding quantity was measured.

Run with::

    python -m unittest test_p2d_coupling
"""

from __future__ import annotations

import math
import re
import shutil
import unittest
from pathlib import Path

import reference_config_p2d as rc
from column import interface, soil
from surface import graded_crust as gc

TOL = rc.TOL
CPP = rc.P2D_ROOT / "cpp"


_STAGE1_CACHE: dict = {}


def evaporate_to_solubility(cells: int = 30, dt: float = 800.0):
    """Stage 1: run the activity-controlled film until the surface saturates.

    Cached per (cells, dt) because several test classes need the same state
    and the run costs about thirty seconds. The declared pair below is the
    coarsest that stays inside the Newton solve's step limit; at dt = 2000 s
    on 20 cells the surface concentration moves far enough in one step for the
    block Jacobian to become singular, which is a declared limitation of the
    scheme rather than something to tune around silently.
    """
    key = (cells, dt)
    if key not in _STAGE1_CACHE:
        state = rc.hydrostatic_state(cells=cells)
        boundary = soil.BoundaryCondition(
            top_flux=rc.film_flux(rc.film()), bottom="no_flux"
        )
        _STAGE1_CACHE[key] = soil.run(
            rc.SOIL, state, boundary, dt, 900, record_stride=50,
            stop=rc.solubility_stop,
        )
    report = _STAGE1_CACHE[key]
    return {
        "state": report["state"].copy(),
        "initial": report["initial"].copy(),
        "transferred": dict(report["transferred"]),
        "history": report["history"],
        "stop_reason": report["stop_reason"],
    }


class TestSoilClosuresMirrorP1R2(unittest.TestCase):
    """Term-by-term transcription checks against main.cc."""

    def test_saturation_and_its_inverse_agree(self):
        for sw in (0.15, 0.25, 0.4, 0.7, 0.95):
            self.assertAlmostEqual(rc.SOIL.sw(rc.SOIL.pw(sw)), sw, delta=1e-12)

    def test_saturation_is_the_van_genuchten_expression(self):
        p = rc.SOIL
        pw = 9.0e4
        pc = soil.P_NON_WETTING_REFERENCE - pw
        m = 1.0 - 1.0 / p.pc_n
        want = p.residual_saturation + (1.0 - p.residual_saturation) / (
            (1.0 + (p.pc_alpha * pc) ** p.pc_n) ** m
        )
        self.assertAlmostEqual(p.sw(pw), want, delta=1e-15)

    def test_saturation_is_one_above_the_reference_pressure(self):
        self.assertEqual(rc.SOIL.sw(soil.P_NON_WETTING_REFERENCE + 1.0), 1.0)

    def test_relative_permeability_uses_the_exponent_form(self):
        """Case 36877 sets KrwMualemExponent = true, so krw(Se=1) = 1."""
        p = rc.SOIL
        self.assertTrue(p.krw_mualem_exponent)
        self.assertAlmostEqual(p.krw(1.0), 1.0, delta=1e-12)
        self.assertEqual(p.krw(p.residual_saturation), 0.0)
        self.assertLess(p.krw(0.3), p.krw(0.6))

    def test_effective_diffusion_is_millington_quirk(self):
        p = rc.SOIL
        for sw in (0.2, 0.5, 0.9):
            want = p.porosity**1.5 * sw**3.5 * p.salt_diffusivity
            self.assertAlmostEqual(
                p.effective_diffusion_coefficient(sw), want, delta=1e-30
            )

    def test_mass_and_mole_fraction_invert(self):
        for x in (1e-6, 0.0036, 0.05, 0.1, 0.3):
            w = rc.SOIL.mass_fraction(x)
            self.assertAlmostEqual(rc.SOIL.mole_fraction(w), x, delta=1e-14)

    def test_the_declared_equilibrium_matches_the_frozen_solubility(self):
        """C_eq is the frozen Transport.SolubilityMoleFraction converted."""
        self.assertAlmostEqual(
            rc.C_EQ, rc.SOIL.mass_fraction(rc.SOLUBILITY_MOLE_FRACTION), delta=1e-18
        )
        self.assertAlmostEqual(rc.C_EQ, 0.2649408697676898, delta=1e-15)

    def test_density_law_is_the_frozen_one(self):
        p = rc.SOIL
        w = p.mass_fraction(0.0036)
        want = 998.203 + 1000.0 * w * (0.670804 + 0.373854 * w)
        self.assertAlmostEqual(p.rho(0.0036), want, delta=1e-12)

    def test_storage_is_the_two_p1r2_equations(self):
        p = rc.SOIL
        total, salt = soil.storage(p, 9.0e4, 0.0036)
        theta_rho = p.porosity * p.sw(9.0e4) * p.rho(0.0036)
        self.assertAlmostEqual(total, theta_rho, delta=1e-15)
        self.assertAlmostEqual(salt, theta_rho * p.mass_fraction(0.0036), delta=1e-15)


class TestHydrostaticSteadyState(unittest.TestCase):
    def test_the_p1r2_initial_profile_is_an_exact_steady_state(self):
        state = rc.hydrostatic_state(cells=20)
        boundary = soil.BoundaryCondition(top_flux=rc.no_flux, bottom="no_flux")
        residual = soil.residual(rc.SOIL, state, state.copy(), boundary, 1.0)
        scale = soil.inventory(rc.SOIL, state)["liquid_kg_m2"]
        self.assertLess(max(abs(v) for v in residual), 1e-12 * scale)

    def test_it_does_not_drift(self):
        state = rc.hydrostatic_state(cells=20)
        before = soil.inventory(rc.SOIL, state)
        boundary = soil.BoundaryCondition(top_flux=rc.no_flux, bottom="no_flux")
        report = soil.run(rc.SOIL, state, boundary, 1000.0, 20, record_stride=20)
        after = soil.inventory(rc.SOIL, state)
        for key in before:
            self.assertLess(
                abs(after[key] - before[key]) / abs(before[key]),
                TOL["steady_state_drift_relative"],
                key,
            )
        self.assertEqual(report["stop_reason"], "steps")


class TestPrescribedFluxBudget(unittest.TestCase):
    def test_component_budgets_close(self):
        state = rc.hydrostatic_state(cells=20)
        boundary = soil.BoundaryCondition(
            top_flux=rc.prescribed_flux(1.0e-5, 1.0e-7), bottom="no_flux"
        )
        report = soil.run(rc.SOIL, state, boundary, 50.0, 40, record_stride=40)
        report_audit = soil.audit(
            rc.SOIL, report["state"], report["initial"], report["transferred"]
        )
        for name in ("liquid", "salt", "water"):
            self.assertLess(
                report_audit[f"{name}_residual_relative"],
                TOL["budget_relative"],
                name,
            )

    def test_the_transfer_is_the_flux_the_scheme_applied(self):
        """Backward Euler uses the END-of-step flux; the audit must too."""
        state = rc.hydrostatic_state(cells=20)
        boundary = soil.BoundaryCondition(
            top_flux=rc.prescribed_flux(1.0e-5, 1.0e-7), bottom="no_flux"
        )
        info = soil.step(rc.SOIL, state, boundary, 100.0)
        self.assertAlmostEqual(info["top_total_kg_m2_s"], 1.0e-5, delta=1e-18)
        self.assertAlmostEqual(info["top_salt_kg_m2_s"], 1.0e-7, delta=1e-20)
        self.assertEqual(info["bottom_total_kg_m2_s"], 0.0)


class TestFilmEvaporation(unittest.TestCase):
    """Task item 4 and 5."""

    @classmethod
    def setUpClass(cls):
        cls.report = evaporate_to_solubility()

    def test_the_surface_reaches_solubility_in_finite_time(self):
        self.assertEqual(self.report["stop_reason"], "soil_top_reached_solubility")
        self.assertGreater(self.report["state"].time, 0.0)
        self.assertLess(self.report["state"].time, 1.0e7)

    def test_a_too_large_step_fails_loudly_rather_than_silently(self):
        """The declared step limit is enforced by the solver, not assumed."""
        state = rc.hydrostatic_state(cells=20)
        boundary = soil.BoundaryCondition(
            top_flux=rc.film_flux(rc.film()), bottom="no_flux"
        )
        report = soil.run(
            rc.SOIL, state, boundary, 2000.0, 900, record_stride=100,
            stop=rc.solubility_stop,
        )
        self.assertNotEqual(report["stop_reason"], "steps")
        self.assertIn("Error", report["stop_reason"])

    def test_evaporation_falls_monotonically(self):
        fluxes = [row["top_total_kg_m2_s"] for row in self.report["history"]]
        self.assertGreater(len(fluxes), 3)
        for earlier, later in zip(fluxes, fluxes[1:]):
            self.assertLess(later, earlier)

    def test_the_flux_is_the_activity_controlled_film(self):
        state = self.report["state"]
        w = rc.SOIL.mass_fraction(state.mole_fraction[-1])
        want = rc.film().flux_from_mass_fraction(w)
        got = rc.film_flux(rc.film())(rc.SOIL, state)[0]
        self.assertAlmostEqual(got, want, delta=1e-18)

    def test_nothing_is_prescribed_or_normalised(self):
        """The initial flux is a physical prediction, not an input."""
        fresh = rc.hydrostatic_state(cells=20)
        flux = rc.film_flux(rc.film())(rc.SOIL, fresh)[0]
        self.assertGreater(flux, 0.0)
        self.assertLess(flux, rc.film().flux_from_mass_fraction(0.0))

    def test_salt_is_conserved_exactly(self):
        report_audit = soil.audit(
            rc.SOIL, self.report["state"], self.report["initial"],
            self.report["transferred"],
        )
        self.assertEqual(report_audit["transferred_salt_kg_m2"], 0.0)
        self.assertLess(
            report_audit["salt_residual_relative"], TOL["budget_relative"]
        )

    def test_water_and_liquid_budgets_close(self):
        report_audit = soil.audit(
            rc.SOIL, self.report["state"], self.report["initial"],
            self.report["transferred"],
        )
        for name in ("liquid", "water"):
            self.assertLess(
                report_audit[f"{name}_residual_relative"],
                TOL["budget_relative"],
                name,
            )

    def test_the_molality_limit_is_refused_not_extrapolated(self):
        state = rc.hydrostatic_state(cells=4)
        state.mole_fraction = [0.4] * 4
        with self.assertRaises(soil.ColumnDomainExit):
            rc.film_flux(rc.film())(rc.SOIL, state)


class TestCrustBaseGeneralisation(unittest.TestCase):
    """Task item 2: the verified graded_crust model must be preserved."""

    def setUp(self):
        self.params = rc.crust_parameters()
        self.film = rc.film()
        delta = self.params.reaction_length()
        self.state = gc.GradedCrustState(top_cell_size=delta / 8.0)
        gc.activate(
            self.state, self.params, rc.DECLARED["seed_height_m"],
            self.film.flux_from_mass_fraction(rc.C_EQ), 0.0,
            window=40.0 * delta, ratio=1.25, maximum_cell=10.0 * delta,
        )

    def test_default_and_explicit_c_eq_are_identical(self):
        flux = self.film.flux_from_mass_fraction(rc.C_EQ)
        a = gc.solve_profile(self.state, self.params, flux)
        b = gc.solve_profile(self.state, self.params, flux, rc.C_EQ)
        for name in ("top_supersaturation", "growth_velocity",
                     "internal_reaction", "dissolved_salt", "water",
                     "solution_mass_flux_base", "salt_flux_base_route"):
            self.assertEqual(getattr(a, name), getattr(b, name), name)

    def test_a_richer_base_raises_the_salt_flux_into_the_crust(self):
        flux = self.film.flux_from_mass_fraction(rc.C_EQ)
        lean = gc.solve_profile(self.state, self.params, flux, 0.99 * rc.C_EQ)
        rich = gc.solve_profile(self.state, self.params, flux, 1.01 * rc.C_EQ)
        self.assertLess(lean.salt_flux_base_route, rich.salt_flux_base_route)

    def test_an_undersaturated_base_dissolves(self):
        flux = self.film.flux_from_mass_fraction(rc.C_EQ)
        lean = gc.solve_profile(self.state, self.params, flux, 0.99 * rc.C_EQ)
        self.assertLess(min(lean.profile), 0.0)

    def test_prat_equation_17_is_the_area_closure(self):
        for eps in (0.43, 0.30, 0.20):
            self.assertAlmostEqual(
                self.params.area(eps), 36.0 * (1.0 - eps) / 25.0e-6, delta=1e-6
            )

    def test_the_tortuosity_closure_is_sqrt_porosity(self):
        for eps in (0.43, 0.30, 0.20):
            self.assertAlmostEqual(
                self.params.diffusivity(eps), math.sqrt(eps) * 1.3e-9, delta=1e-22
            )


class TestSeed(unittest.TestCase):
    """Task item 6."""

    def setUp(self):
        self.report = evaporate_to_solubility()
        self.column = self.report["state"]
        self.before = soil.inventory(rc.SOIL, self.column)["salt_kg_m2"]
        self.crust, self.salt, self.water = interface.seed_from_soil(
            rc.SOIL, self.column, rc.crust_parameters(), rc.film(),
            rc.DECLARED["seed_height_m"],
        )

    def test_the_seed_is_budget_neutral(self):
        after = soil.inventory(rc.SOIL, self.column)["salt_kg_m2"]
        self.assertLess(
            abs((after + self.salt) - self.before) / self.before,
            TOL["budget_relative"],
        )

    def test_the_seed_comes_out_of_the_soil(self):
        after = soil.inventory(rc.SOIL, self.column)["salt_kg_m2"]
        self.assertLess(after, self.before)
        self.assertGreater(self.salt, 0.0)

    def test_the_seed_is_a_connected_saturated_crust(self):
        params = rc.crust_parameters()
        self.assertEqual(params.crust_saturation, 1.0)
        self.assertGreater(self.crust.minimum_porosity, params.connectivity_porosity)
        self.assertAlmostEqual(
            self.crust.height, rc.DECLARED["seed_height_m"], delta=1e-12
        )
        self.assertGreater(
            self.crust.height / params.reaction_length(), 4.0
        )

    def test_an_oversized_seed_is_refused(self):
        with self.assertRaises(interface.InterfaceDomainExit):
            interface.seed_from_soil(
                rc.SOIL, self.column.copy(), rc.crust_parameters(), rc.film(),
                1.0e-2,
            )


class TestInterfaceCoupling(unittest.TestCase):
    """Task items 7, 8 and 9."""

    @classmethod
    def setUpClass(cls):
        report = evaporate_to_solubility()
        cls.column = report["state"]
        cls.initial_column = cls.column.copy()
        cls.params = rc.crust_parameters()
        cls.film = rc.film()
        cls.crust, cls.seed_salt, cls.seed_water = interface.seed_from_soil(
            rc.SOIL, cls.column, cls.params, cls.film,
            rc.DECLARED["seed_height_m"],
        )

    def fresh_state(self):
        return interface.CoupledState(
            column=self.column.copy(),
            crust=self.crust.copy(),
            seed_salt_kg_m2=self.seed_salt,
            seed_water_kg_m2=self.seed_water,
        )

    def test_the_residual_changes_sign_across_equilibrium(self):
        """One unknown, one equation: the root exists and is bracketed."""
        state = self.fresh_state()
        record = interface.step(rc.SOIL, state, self.params, self.film, 5.0)
        self.assertLess(
            record["interface_residual_relative"],
            TOL["interface_closure_relative"],
        )

    def test_the_crust_pins_the_interface_to_halite_saturation(self):
        state = self.fresh_state()
        record = interface.step(rc.SOIL, state, self.params, self.film, 5.0)
        self.assertAlmostEqual(record["interface_over_c_eq"], 1.0, delta=1e-4)

    def test_the_interface_residual_is_monotone(self):
        """Monotone residual means the compatible concentration is unique."""
        values = []
        for factor in (0.9995, 0.9998, 1.0, 1.0002, 1.0005):
            response = interface.crust_response(
                self.crust, self.params, self.film, factor * rc.C_EQ
            )
            soil_flux = interface.soil_interface_salt_flux(
                rc.SOIL, self.column, response["total_mass_flux"],
                factor * rc.C_EQ,
            )
            values.append(soil_flux - response["salt_mass_flux"])
        for earlier, later in zip(values, values[1:]):
            self.assertLess(later, earlier)
        self.assertGreater(values[0], 0.0)
        self.assertLess(values[-1], 0.0)

    def test_neither_concentration_nor_flux_is_prescribed(self):
        """Item 7: the solved interface is not an input anywhere."""
        state = self.fresh_state()
        record = interface.step(rc.SOIL, state, self.params, self.film, 5.0)
        self.assertNotEqual(record["interface_mass_fraction"], rc.C_EQ)
        self.assertNotEqual(
            record["interface_mass_fraction"],
            rc.SOIL.mass_fraction(self.column.mole_fraction[-1]),
        )

    def test_the_soil_is_given_the_inventory_route(self):
        state = self.fresh_state()
        record = interface.step(rc.SOIL, state, self.params, self.film, 5.0)
        self.assertLess(
            record["route_2_minus_3_relative"], TOL["interface_closure_relative"]
        )

    def test_the_transport_route_is_an_independent_audit(self):
        """It differs by the declared linearisation residual, not by chance."""
        state = self.fresh_state()
        record = interface.step(rc.SOIL, state, self.params, self.film, 5.0)
        measured = record["route_1_minus_2_relative"]
        predicted = record["route_1_minus_2_predicted_relative"]
        self.assertGreater(measured, 0.0)
        self.assertLess(measured, 1.0e-4)
        self.assertLess(abs(measured - predicted) / measured, 1.0)

    def test_the_reduction_stays_inside_pe_much_less_than_da(self):
        state = self.fresh_state()
        record = interface.step(rc.SOIL, state, self.params, self.film, 5.0)
        self.assertLess(record["peclet_over_damkohler"], 1.0e-3)

    def test_the_coupled_budget_closes(self):
        state = self.fresh_state()
        # The combined budget starts from the PRE-SEED soil, because the seed
        # was formed out of that soil rather than added to the system.
        initial = self.initial_column.copy()
        for _ in range(4):
            interface.step(rc.SOIL, state, self.params, self.film, 5.0)
        report = interface.audit(
            rc.SOIL, state, initial, self.seed_salt, self.seed_water
        )
        self.assertLess(report["salt_residual_relative"], 1.0e-8)
        self.assertLess(report["water_residual_relative"], 1.0e-8)

    def test_an_undersaturated_soil_is_refused(self):
        """Item 9: outside the halite-equilibrated envelope the model stops."""
        state = self.fresh_state()
        for index in range(state.column.cells):
            state.column.mole_fraction[index] = 1.0e-4
        with self.assertRaises(interface.InterfaceDomainExit):
            interface.step(rc.SOIL, state, self.params, self.film, 5.0)
        self.assertTrue(state.stopped)
        self.assertIsNotNone(state.stop_reason)

    def test_a_stopped_model_refuses_to_continue(self):
        state = self.fresh_state()
        state.stopped = True
        state.stop_reason = "declared test stop"
        with self.assertRaises(interface.InterfaceDomainExit):
            interface.step(rc.SOIL, state, self.params, self.film, 5.0)

    def test_the_growth_velocity_is_never_capped(self):
        response = interface.crust_response(
            self.crust, self.params, self.film, rc.C_EQ
        )
        solution = response["solution"]
        expected = (
            self.params.beta
            * self.params.reaction_coefficient
            * self.params.rho
            * solution.top_supersaturation
            / self.params.crystal_density
        )
        self.assertLess(
            abs(solution.growth_velocity - expected) / expected, 1.0e-14
        )


class TestCppPort(unittest.TestCase):
    """Task item 3. The module cannot be built here; what can be checked is."""

    def test_no_cpp_toolchain_is_available(self):
        """Recorded as a fact of this environment, not assumed."""
        for tool in ("g++", "clang++", "cl", "cmake", "make", "ninja"):
            self.assertIsNone(shutil.which(tool), tool)

    def test_the_header_and_its_checker_exist(self):
        for name in ("pratcrust.hh", "verify_vectors.cc", "CMakeLists.txt",
                     "reference_vectors.txt", "README.md"):
            self.assertTrue((CPP / name).is_file(), name)

    def test_the_header_states_that_it_is_unverified(self):
        text = (CPP / "pratcrust.hh").read_text(encoding="utf-8")
        self.assertIn("NEVER BEEN COMPILED", text)

    def test_every_reference_vector_is_computed_by_the_checker(self):
        """A silent drift between the two sides would void the whole check."""
        reference = (CPP / "reference_vectors.txt").read_text(encoding="utf-8")
        names = {
            line.split()[0]
            for line in reference.splitlines()
            if line and not line.startswith("#")
        }
        source = (CPP / "verify_vectors.cc").read_text(encoding="utf-8")
        plain = set(re.findall(r'add\(v,\s*"([A-Za-z0-9_]+)"', source))
        suffixes = set(re.findall(r'add\(v,\s*prefix\s*\+\s*"([A-Za-z0-9_]+)"', source))
        bases = re.findall(r'\{"([A-Za-z0-9_]+)",\s*[-0-9.*a-zA-Z ]+\}', source)
        computed = set(plain)
        for base in bases:
            for suffix in suffixes:
                computed.add(f"profile_{base}_{suffix}")
        self.assertEqual(names, computed)
        self.assertGreaterEqual(len(names), 60)

    def test_the_declared_vector_configuration_matches_on_both_sides(self):
        source = (CPP / "verify_vectors.cc").read_text(encoding="utf-8")
        import generate_cpp_vectors as gen
        for name, value in (
            ("SEED_HEIGHT", gen.SEED_HEIGHT),
            ("CELLS_PER_DELTA", gen.CELLS_PER_DELTA),
            ("WINDOW_LENGTHS", gen.WINDOW_LENGTHS),
            ("RATIO", gen.RATIO),
            ("MAXIMUM_CELL_LENGTHS", gen.MAXIMUM_CELL_LENGTHS),
            ("DT", gen.DT),
        ):
            match = re.search(rf"{name}\s*=\s*([0-9.eE+-]+);", source)
            self.assertIsNotNone(match, name)
            self.assertAlmostEqual(float(match.group(1)), value, delta=1e-15, msg=name)
        match = re.search(r"STEPS\s*=\s*(\d+);", source)
        self.assertEqual(int(match.group(1)), gen.STEPS)

    def test_the_vectors_are_current(self):
        """Regenerating must not change the committed file."""
        import generate_cpp_vectors as gen
        rows = gen.vectors()
        reference = (CPP / "reference_vectors.txt").read_text(encoding="utf-8")
        stored = {}
        for line in reference.splitlines():
            if line and not line.startswith("#"):
                name, value = line.split()
                stored[name] = float(value)
        self.assertEqual(set(stored), {name for name, _ in rows})
        for name, value in rows:
            denominator = max(abs(value), 1e-300)
            self.assertLess(abs(stored[name] - value) / denominator, 1e-13, name)


if __name__ == "__main__":
    unittest.main()
