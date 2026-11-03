"""Numerical tests for the command-15 acceptance and comparison logic.

Pure post-processing; no solver, no network, no launch path is exercised.
"""
import copy,json,math,unittest
from pathlib import Path
import numpy as np
import wp8_full4t_campaign as f
from analyze_wp8 import rows,pairs
from single_wavelength_post import widths,params
from wp6_paired_stopping import same_time,validate_leg

OUT=f.OUT;CACHE=f.CACHE;STAGE=f.STAGE;ROOT=f.ROOT
RUN=OUT/STAGE
HIST=ROOT/'analysis/wp8_strict_background_20261103'


class TestTimeAxis(unittest.TestCase):
    def test_accepted_times_are_the_exact_double_accumulation(self):
        plan=json.loads((RUN/'plan.json').read_text());ev=rows(RUN/'surface_events.csv')
        t=0.
        for i,r in enumerate(ev):
            if i:t+=plan['dt_s']
            self.assertEqual(float(r['time_s']),t)

    def test_endpoint_cannot_equal_decimal_TEnd_but_stays_inside_the_budget(self):
        plan=json.loads((RUN/'plan.json').read_text());ev=rows(RUN/'surface_events.csv')
        last=float(ev[-1]['time_s']);H=plan['horizon_s']
        self.assertNotEqual(last,H)
        self.assertFalse(same_time(last,H))
        self.assertLessEqual(abs(H-last),plan['steps_cap']*math.ulp(H))
        self.assertLessEqual(last,H+1e-7)

    def test_clipped_endpoint_is_excluded_rather_than_interpolated(self):
        a=[0.,100.,501695.8170227767];b=[0.,100.,501695.81702277]
        self.assertFalse(same_time(a[2],b[2]))
        self.assertEqual(pairs(a,b),[(0,0),(1,1)])


class TestAcceptance(unittest.TestCase):
    def setUp(self):
        self.v=json.loads((OUT/'validation.json').read_text())

    def test_every_gate_passed_and_none_was_dropped(self):
        self.assertTrue(self.v['passed'])
        self.assertEqual(self.v['failures'],[])
        for key in ['exact_horizontal_equality','salt_closure','water_closure','independent_budget',
                    'one_Newton_success_per_step','no_failed_Newton_attempt','accepted_scaled_contract',
                    'reached_requested_horizon','frozen_physics','no_replay_or_restart','hashes']:
            self.assertIn(key,self.v['checks']);self.assertTrue(self.v['checks'][key],key)

    def test_horizon_reached_without_a_threshold_event(self):
        stop=json.loads((RUN/'stop_record.json').read_text())
        self.assertEqual(stop['reason'],'requested_horizon')
        self.assertEqual(stop['accepted_steps'],4800)
        self.assertIsNone(stop['surface_threshold_high_s'])
        self.assertFalse(self.v['censored_by_threshold'])
        self.assertLess(self.v['metrics']['max_domain_x'],0.1)

    def test_a_fabricated_post_threshold_accepted_state_is_rejected(self):
        ev=rows(RUN/'surface_events.csv');stop=json.loads((RUN/'stop_record.json').read_text())
        validate_leg(ev,stop,0)
        bad=copy.deepcopy(ev);bad[-2]['global_threshold_reached']='1';bad[-2]['x_domain_max']='0.2'
        with self.assertRaises(ValueError):validate_leg(bad,stop,0)

    def test_scaled_criterion_subcriteria_all_bind_at_the_old_failure_step(self):
        plan=json.loads((RUN/'plan.json').read_text());c=plan['criterion']
        di=rows(RUN/'convergence_v2.csv')
        acc=[r for r in di if int(r['accepted']) and int(r['step'])==4414]
        self.assertEqual(len(acc),1)
        r=acc[0]
        self.assertLessEqual(float(r['relative_shift']),c['relative_shift'])
        self.assertLessEqual(float(r['scaled_nonlinear']),c['scaled_residual'])
        self.assertLessEqual(float(r['linear_eta_component']),c['componentwise_linear_backward_error'])
        self.assertLessEqual(float(r['regional_salt_error_m']),c['regional_step_salt_m'])
        self.assertLessEqual(float(r['regional_water_error_m']),c['regional_step_water_m'])
        # the step the old fixed absolute limit could not accept is accepted here
        self.assertGreater(float(r['salt_max']),2e-20)


class TestIndependentBudget(unittest.TestCase):
    def test_storage_reconstruction_matches_the_recorded_inventory(self):
        a=dict(np.load(CACHE/STAGE/'post/profiles.npz'));bud=rows(RUN/'transport_budget.csv')
        par=params(RUN/'params.input');phi=float(par['Problem.Porosity'])
        z=a['z_m'];w=widths(z);g='upper_0.950000'
        rr=[r for r in bud if r['region']==g];mask=z>=float(rr[0]['actual_z_over_H'])*.2
        for i in [0,1200,2400,3600,4413,4800]:
            self.assertLess(abs(phi*a['sx_mean'][i,mask]@w[mask]-float(rr[i]['inventory_salt_m'])),1e-14)

    def test_budget_identity_is_exact_in_the_recorded_columns(self):
        bud=rows(RUN/'transport_budget.csv');g='free_domain'
        rr=[r for r in bud if r['region']==g];I0=float(rr[0]['inventory_salt_m']);C=0.
        for i,r in enumerate(rr):
            if i:C+=float(r['dt_s'])*float(r['cut_salt_total_down_m_s'])
            self.assertEqual(C,float(r['cumulative_cut_salt_m']))
            self.assertAlmostEqual(float(r['inventory_salt_m'])-I0+C,float(r['salt_residual_m']),delta=1e-18)


class TestHistoricalComparison(unittest.TestCase):
    def setUp(self):
        self.c=json.loads((OUT/'historical_comparison.json').read_text())

    def test_byte_identity_is_measured_not_required(self):
        r=self.c['references']['historical_4T_prefix']
        self.assertFalse(self.c['byte_identity_required'])
        self.assertLess(r['identical_field_fraction'],1.0)
        self.assertGreater(r['identical_field_sha256'],0)
        self.assertLess(r['max_abs_dx'],1e-14)

    def test_comparison_uses_matched_times_only(self):
        r=self.c['references']['historical_4T_prefix']
        self.assertFalse(self.c['interpolation_used'])
        self.assertEqual(r['shared_accepted_states'],r['reference_accepted_steps']+1)
        self.assertAlmostEqual(r['last_shared_over_T'],3.6775,places=6)

    def test_corrected_solver_passes_the_old_failure_window(self):
        w=self.c['old_failure_window']
        self.assertTrue(w['passed'])
        self.assertGreaterEqual(w['historical_failed_min_salt_residual'],w['historical_fixed_absolute_salt_threshold'])
        self.assertLess(w['corrected_iterations'],w['historical_failed_iterations'])

    def test_newton_cost_strictly_improved(self):
        n=self.c['newton_iteration_counts']
        self.assertLess(n['corrected_max'],n['historical_max'])
        self.assertLess(n['corrected_mean'],n['historical_mean'])


class TestScientificGuards(unittest.TestCase):
    def setUp(self):
        self.s=json.loads((OUT/'scientific_analysis.json').read_text())

    def test_no_effective_time_or_surrogate_label_is_assigned(self):
        for k in ['t_eff','crystallization_delay','surrogate_label']:
            self.assertIsNone(self.s[k])

    def test_no_solubility_event_in_any_compared_run(self):
        for e in self.s['solubility_threshold_events']:
            self.assertFalse(e['threshold_observed'])
            self.assertLess(e['max_x_domain_max'],e['threshold'])

    def test_grid_sensitivity_is_reported_separately_from_validation(self):
        self.assertIn('NOT a convergence proof',self.s['grid_sensitivity']['scope'])
        self.assertNotIn('grid',json.dumps(self.s['background']))

    def test_response_below_the_solver_criterion_floor_is_not_claimed(self):
        floor=self.s['solver_criterion_caveat']['bound_from_background_to_background_difference_max_abs_dx']
        early=self.s['comparisons']['R24_B']['checkpoints']['T']['mean_response_Linf']
        self.assertLess(abs(early),floor)
        late=self.s['comparisons']['R24_B']['checkpoints']['4T']['mean_response_Linf']
        self.assertGreater(abs(late),1e3*floor)


if __name__=='__main__':unittest.main()
