"""WP9 numerical tests. Pure checks on saved products and on the mode/time algebra."""
import json,math,unittest
from pathlib import Path
import numpy as np
import wp9_common as c
import wp9_r96_protocol as P


def load(name):
    p=c.OUT/name
    return json.loads(p.read_text(encoding='utf-8')) if p.exists() else None


class TestModeAlgebra(unittest.TestCase):
    def test_box_mode_to_lsa_index_is_exact(self):
        self.assertEqual(c.NX_PER_M,14.0)
        self.assertLess(c.WIDTH_CONSISTENCY,1e-15)
        for m,nx in [(1,14),(2,28),(6,84),(10,140)]:
            self.assertEqual(c.lsa_nx(m),nx)

    def test_wavelength_and_wavenumber_agree_between_conventions(self):
        for m in (1,2,3,6,10,24):
            lsa_lambda=c.FULL_WIDTH_M/c.lsa_nx(m)
            self.assertAlmostEqual(lsa_lambda,c.wavelength(m),places=15)
            self.assertLess(abs(lsa_lambda/c.wavelength(m)-1),1e-15)
            self.assertAlmostEqual(2*math.pi/c.wavelength(m),c.k_phys(m),places=9)

    def test_canonical_mode_spans_one_box_width(self):
        self.assertAlmostEqual(c.wavelength(c.CANONICAL_M),c.WIDTH_M,places=18)
        self.assertEqual(c.wavelength(c.CANONICAL_M),c.WIDTH_M)

    def test_cells_per_wavelength(self):
        self.assertAlmostEqual(c.cells_per_wavelength(10,24),4.8)
        self.assertAlmostEqual(c.cells_per_wavelength(10,48),9.6)
        self.assertAlmostEqual(c.cells_per_wavelength(10,96),19.2)

    def test_global_critical_mode_does_not_fit_the_box(self):
        self.assertEqual(c.lsa_nx(1),14.0)            # the smallest box mode is nx=14
        self.assertLess(5,c.lsa_nx(1))                # the stored critical nx=5 is below it
        self.assertGreater(c.FULL_WIDTH_M/5,c.WIDTH_M)  # its wavelength exceeds the box width


class TestTimeAxis(unittest.TestCase):
    def test_half_step_is_an_exact_binary_halving(self):
        self.assertEqual(2*(c.DT_S/2),c.DT_S)

    def test_paired_accumulation_stays_inside_the_registered_budget(self):
        ev=P.time_axis_evidence()
        self.assertTrue(ev['exact_halving'])
        self.assertLessEqual(ev['max_paired_representation_error_s'],ev['accumulation_budget_s'])
        self.assertGreater(ev['max_paired_representation_error_s'],32*math.ulp(c.HORIZON_S))

    def test_frame_index_identity(self):
        a=0.0;coarse=[]
        for _ in range(11):coarse.append(a);a+=c.DT_S
        b=0.0;fine=[]
        for _ in range(21):fine.append(b);b+=c.DT_S/2
        for j in range(11):
            self.assertLess(abs(coarse[j]-fine[2*j]),1e-9)

    def test_step_count_reaches_the_horizon(self):
        self.assertAlmostEqual(9600*(c.DT_S/2),c.HORIZON_S,places=6)


class TestAcceptedLegs(unittest.TestCase):
    """Gates that can only be checked once a leg has run."""

    def legs(self):
        out={}
        for name,rel in [('M1','r96_validation.json'),('M0','r96_m0/validation.json')]:
            p=c.OUT/rel
            if p.exists():out[name]=json.loads(p.read_text(encoding='utf-8'))
        return out

    def test_every_gate_passed(self):
        legs=self.legs()
        if not legs:self.skipTest('no leg validated yet')
        for k,v in legs.items():
            self.assertTrue(v['passed'],f'{k}: {v["failures"]}')

    def test_only_the_last_step_is_truncated(self):
        for k,v in self.legs().items():
            m=v['metrics']
            self.assertTrue(v['checks']['fixed_dt_on_every_step_but_the_last'],k)
            self.assertGreaterEqual(m['final_step_truncation_s'],0.0,k)
            self.assertLessEqual(m['final_step_truncation_s'],m['time_axis_budget_s'],k)
            self.assertEqual(m['fixed_dt_accumulation_drift_s'],0.0,k)

    def test_endpoint_lands_on_the_requested_horizon(self):
        for k,v in self.legs().items():
            if v['censored_by_threshold']:continue
            self.assertLessEqual(v['metrics']['endpoint_error_s'],
                                 v['metrics']['time_axis_budget_s'],k)

    def test_shared_time_axis_with_R48(self):
        for k,v in self.legs().items():
            m=v['metrics']
            self.assertLessEqual(m['max_paired_time_difference_s'],m['paired_time_budget_s'],k)

    def test_both_legs_use_the_same_executable(self):
        legs=self.legs()
        if len(legs)<2:self.skipTest('both legs not validated yet')
        for k,v in legs.items():
            self.assertTrue(v['checks']['elf_is_the_R24_R48_executable'],k)

    def test_the_pair_runs_on_two_hosts_without_a_protocol_difference(self):
        legs=self.legs()
        if len(legs)<2:self.skipTest('both legs not validated yet')
        self.assertNotEqual(legs['M1']['host'],legs['M0']['host'])
        port=load('hpcc_portability/portability.json')
        ts=load('hpcc_portability/thread_scaling.json')
        self.assertTrue(port and port['numerically_identical'])
        self.assertTrue(ts and ts['all_bitwise_identical'])


class TestProtocol(unittest.TestCase):
    def setUp(self):
        self.proto=load('r96_protocol.json');self.adm=load('r96_admission.json')
        if self.proto is None:self.skipTest('protocol not built')

    def test_only_the_registered_keys_changed(self):
        self.assertEqual(sorted(self.adm['input_key_changes']),self.adm['allowed_key_changes'])

    def test_executable_is_the_R24_R48_one(self):
        self.assertEqual(self.proto['numerics']['executable']['sha256'],
                         '112076759fdff0942000307f729c2bd3c7e96963229b43c1b67c1f12bbf9df81')

    def test_eigenprofile_is_the_original_file(self):
        self.assertEqual(self.adm['eigenprofile_sha256'],
                         'd721e605d126e555f92b43200aa73b999ed4e2b61397edd1e98c0d9fddb36575')

    def test_forbidden_list_covers_the_command(self):
        for item in ['R96-M0','R192','multiwavelength','full-domain','8T']:
            self.assertIn(item,self.proto['forbidden'])

    def test_hard_wall_is_conservative(self):
        r=self.proto['resources']
        self.assertGreaterEqual(r['hard_wall_s'],1.5*r['projected_worker_wall_s'])

    def test_the_timestep_confound_is_registered(self):
        names=[k['name'] for k in self.proto['known_confounds']]
        self.assertIn('simultaneous timestep change',names)


class TestCorrections(unittest.TestCase):
    def setUp(self):
        self.d=load('corrections_v2.json')
        if self.d is None:self.skipTest('corrections not built')

    def test_whole_domain_increase_is_bottom_input(self):
        w=self.d['correction_1_inventory']['whole_domain']
        self.assertLess(w['attribution_relative_error'],1e-9)
        self.assertGreater(w['delta_salt_m'],0)

    def test_water_budget_closes(self):
        w=self.d['correction_1_inventory']['whole_domain']
        self.assertLess(abs(w['water_balance_residual_m']),1e-12)

    def test_evaporation_removes_no_salt(self):
        w=self.d['correction_1_inventory']['whole_domain']
        self.assertGreater(w['evaporated_water_m'],0)
        self.assertAlmostEqual(w['delta_salt_m'],w['cumulative_bottom_salt_input_m'],places=15)

    def test_every_plane_balance_closes(self):
        for r in self.d['correction_1_inventory']['per_plane']:
            self.assertLess(abs(r['balance_residual_m']),1e-14)

    def test_censoring_is_recorded_as_censoring(self):
        cz=self.d['correction_2_censoring']
        self.assertEqual(cz['status'],'right_censored')
        for r in cz['trajectories']:
            self.assertFalse(r['crossed'])
            self.assertLess(r['max_local_x'],c.SOLUBILITY)

    def test_v2_executable_guard_is_present_in_source(self):
        for e in self.d['correction_3_v2_executable']['evidence']:
            self.assertTrue(e['present'])

    def test_executable_gap_is_round_off(self):
        g=self.d['elf_gap']
        self.assertLess(g['max_abs_x_top_area_mean'],1e-12)
        self.assertTrue(g['times_within_budget'])


class TestLsaAudit(unittest.TestCase):
    def setUp(self):
        self.d=load('lsa_audit.json')
        if self.d is None:self.skipTest('audit not built')

    def test_stored_anchors_reproduced(self):
        self.assertTrue(self.d['anchors_reproduced'])
        a=self.d['reproduced_anchors']
        self.assertLess(a['sigma_c_tc_at2p0']['relative_difference'],1e-9)
        self.assertEqual(a['critical_rayleigh_at_tc']['recomputed_critical_nx'],
                         a['critical_rayleigh_at_tc']['stored_critical_nx'])

    def test_wide_scan_peaks_are_interior(self):
        for k,v in self.d['per_checkpoint'].items():
            self.assertFalse(v['wide_scan_peak_at_scan_edge'],k)

    def test_fastest_box_mode_does_not_decrease_with_time(self):
        order=['t_c','2t_c','T','2T','3T','3.6775T','4T']
        seq=[self.d['per_checkpoint'][k]['fastest_box_mode_m'] for k in order if k in self.d['per_checkpoint']]
        self.assertEqual(seq,sorted(seq))

    def test_m10_is_stable_at_tc_and_unstable_at_4T(self):
        self.assertLess(self.d['per_checkpoint']['t_c']['sigma_per_s_m10'],0)
        self.assertGreater(self.d['per_checkpoint']['4T']['sigma_per_s_m10'],0)

    def test_assumptions_are_stated(self):
        self.assertGreaterEqual(len(self.d['quasi_steady_assumptions']),4)
        self.assertFalse(self.d['box_mode_algebra']['global_critical_nx_representable_in_box'])


class TestLsaVersusDns(unittest.TestCase):
    def setUp(self):
        self.d=load('lsa_versus_dns.json')
        if self.d is None:self.skipTest('comparison not built')

    def test_only_admitted_points_are_compared(self):
        for r in self.d['rows']:
            if not r['admitted_small_perturbation']:
                self.assertNotEqual(r['comparison'],'quantitative')

    def test_admitted_points_agree_within_ten_percent(self):
        self.assertLess(self.d['agreement']['max_absolute_relative_difference'],0.10)

    def test_window_ends_before_four_T(self):
        for leg,w in self.d['windows'].items():
            self.assertLess(w['last_time_over_T'],4.0,leg)


class TestComparison(unittest.TestCase):
    def setUp(self):
        self.d=load('resolution_comparison.json')
        if self.d is None:self.skipTest('comparison not built')

    def test_time_axes_within_budget(self):
        for k,v in self.d['time_axis'].items():
            self.assertTrue(v['within_budget'],k)

    def test_uncertainty_sources_are_named(self):
        u=self.d['uncertainty']
        for k in ['x_resolvability_floor','modal_energy_floor','timestep_prior_on_growth',
                  'transport_timestep_prior']:
            self.assertIn(k,u)

    def test_R48_is_compared_to_R96_directly(self):
        self.assertIn('R48 and R96 are compared directly',self.d['independence_note'])

    def test_discriminator_has_a_threshold(self):
        self.assertIsNotNone(self.d['discriminator']['R48_over_R24_m10'])


class TestScientificGuards(unittest.TestCase):
    def test_no_effective_time_or_crystallization_label_is_asserted(self):
        banned={'t_eff','crystallization_delay','crystallization_delay_label','surrogate_label',
                'training_labels','t_eff_s'}
        def walk(node,where):
            if isinstance(node,dict):
                for k,v in node.items():
                    if k in banned:self.assertIsNone(v,f'{where}/{k}')
                    walk(v,f'{where}/{k}')
            elif isinstance(node,list):
                for i,v in enumerate(node):walk(v,f'{where}[{i}]')
        for p in sorted(c.OUT.glob('*.json')):
            walk(json.loads(p.read_text(encoding='utf-8')),p.name)

    def test_grid_sensitivity_is_not_called_convergence(self):
        d=load('resolution_comparison.json')
        if d is None:self.skipTest('comparison not built')
        self.assertIn('does not establish a convergence order',
                      json.loads((c.OUT/'r96_protocol.json').read_text(encoding='utf-8'))['prediction']['caveat'])


if __name__=='__main__':unittest.main(verbosity=2)
