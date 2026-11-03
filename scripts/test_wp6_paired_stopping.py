"""Synthetic stopping semantics tests, independent of solver execution."""
import contextlib,csv,io,json,shutil,tempfile
from pathlib import Path
import unittest
from wp6_paired_stopping import validate_leg, paired_times, budget_rows_at_events


def fixture(times=(0,1,2,3), event=None):
    rows=[]
    for i,t in enumerate(times):
        crossing=i==len(times)-1 and event is not None
        top=.101 if crossing and event=='surface' else .03
        domain=.102 if crossing else top
        rows.append(dict(step=i,time_s=t,dt_s=t-times[i-1] if i else 0,x_top_max=top,x_domain_max=domain,surface_threshold_reached=int(top>=.1),global_threshold_reached=int(domain>=.1)))
    stop=dict(initial_diagnostics_written=True,accepted_steps=len(times)-1,last_accepted_time_s=times[-1],requested_end_s=10 if event else times[-1],reason=event+'_solubility_threshold' if event else 'requested_horizon')
    for kind in ['surface','global']:
        observed=event is not None and (kind=='global' or event=='surface')
        stop[kind+'_threshold_low_s']=times[max(0,len(times)-2)] if observed else None
        stop[kind+'_threshold_high_s']=times[-1] if observed else None
    return rows,stop


class StoppingTests(unittest.TestCase):
    def leg(self,times=(0,1,2,3),event=None): return validate_leg(*fixture(times,event))
    def test_m0_first(self):
        p=paired_times(self.leg((0,1,2),'surface'),self.leg((0,1,2,3),'surface'))
        self.assertEqual(p['shared_times_s'],[0,1]); self.assertEqual(p['surface_time_difference']['lower_s'],0)
        self.assertEqual(p['surface_time_difference']['upper_s'],2)
    def test_m1_first(self):
        p=paired_times(self.leg((0,1,2,3),'surface'),self.leg((0,1,2),'surface'))
        self.assertEqual(p['shared_times_s'],[0,1]); self.assertEqual(p['surface_time_difference']['upper_s'],0)
    def test_both_censored_different_horizons(self):
        p=paired_times(self.leg(),self.leg((0,1,2)))
        self.assertEqual(p['shared_times_s'],[0,1,2]); self.assertEqual(p['surface_time_difference']['status'],'both_right_censored')
    def test_event_types_separate(self):
        p=paired_times(self.leg((0,1,2),'interior'),self.leg((0,1,2,3),'surface'))
        self.assertTrue(p['different_first_event_types']); self.assertEqual(p['surface_time_difference']['status'],'upper_bound_only')
        self.assertEqual(p['surface_time_difference']['upper_s'],1)
        self.assertEqual(p['any_local_time_difference']['status'],'interval_observed')
    def test_adaptive_no_interpolation(self):
        p=paired_times(self.leg((0,.5,1,2),'surface'),self.leg((0,.25,1,1.5,3),'surface'))
        self.assertEqual(p['shared_times_s'],[0,1])
    def test_censor_bounds(self):
        p=paired_times(self.leg((0,1,2),'surface'),self.leg((0,1,2,3)))
        self.assertEqual(p['surface_time_difference']['lower_s'],1)
        self.assertIsNone(p['surface_time_difference']['upper_s'])
    def test_initial_crossing_empty_comparison(self):
        p=paired_times(self.leg((0,),'surface'),self.leg())
        self.assertEqual(p['shared_times_s'],[])
    def test_crossing_overshoot_excluded(self):
        p=paired_times(self.leg((0,1,2),'surface'),self.leg((0,1,2),'surface'))
        self.assertEqual(p['shared_times_s'],[0,1])
    def test_reject_invalid_stop(self):
        a,b=fixture((0,1,2),'interior');b['reason']='surface_solubility_threshold'
        with self.assertRaises(ValueError):validate_leg(a,b)
    def test_reject_missing_bracket(self):
        a,b=fixture((0,1,2),'surface');b['surface_threshold_low_s']=None
        with self.assertRaises(ValueError):validate_leg(a,b)
    def test_reject_failed_exit(self):
        with self.assertRaises(ValueError):validate_leg(*fixture(),returncode=1)
    def test_reject_bad_time_dt(self):
        for field,value in [('dt_s',-1),('time_s',0),('dt_s',.1)]:
            a,b=fixture();a[1][field]=value
            with self.assertRaises(ValueError):validate_leg(a,b)
    def test_reject_continuation_after_crossing(self):
        a,b=fixture((0,1,2),'surface');a[1].update(x_top_max=.101,x_domain_max=.101,surface_threshold_reached=1,global_threshold_reached=1)
        with self.assertRaises(ValueError):validate_leg(a,b)
    def test_near_but_not_shared_times(self):
        p=paired_times(self.leg((0,1,2)),self.leg((0,1.00000001,2)))
        self.assertEqual(p['shared_times_s'],[0,2])
    def test_terminal_tolerance_does_not_snap_accepted_time(self):
        a,b=fixture();b['requested_end_s']+=1e-8
        v=validate_leg(a,b);self.assertEqual(v['times'][-1],3)
    def test_reject_missing_budget_state(self):
        with self.assertRaises(ValueError):budget_rows_at_events([],self.leg(),'upper')
    def test_initial_surface_event_interval(self):
        p=paired_times(self.leg((0,),'surface'),self.leg((0,),'surface'))
        self.assertEqual(p['surface_time_difference']['lower_s'],0)
        self.assertEqual(p['surface_time_difference']['upper_s'],0)


class PostprocessorIntegrationTests(unittest.TestCase):
    """Synthetic CSV bundles exercise the public full postprocessor, not a solver."""
    def setUp(self):
        import analyze_wp6_review as analysis
        self.analysis=analysis;self.old_out=analysis.OUT
        self.root=Path(tempfile.mkdtemp(prefix='_test_wp6_',dir=analysis.ROOT/'analysis'))
    def tearDown(self):
        self.analysis.OUT=self.old_out
        target=self.root.resolve();assert target.is_relative_to((self.analysis.ROOT/'analysis').resolve()) and target.name.startswith('_test_wp6_')
        shutil.rmtree(target)
    def bundle(self,role,times,event):
        d=self.root/role;(d/'post').mkdir(parents=True)
        ev,stop=fixture(times,event)
        (d/'finished.json').write_text(json.dumps(dict(solver=dict(stop=stop,returncode=0))))
        (d/'post/independent_checks.json').write_text('{"passed":true}')
        cuts=[dict(name='upper_'+str(i),actual_z_over_H=.5+.05*i) for i in range(7)]
        (d/'transport_metadata.json').write_text(json.dumps(dict(cuts=cuts)))
        budget=[]
        for t in times:
            for cut in cuts:
                r=dict(region=cut['name'],time_s=t,inventory_salt_m=1e-5,inventory_initial_salt_m=1e-5)
                for k in ['salt_residual_m','water_residual_m','whole_salt_residual_m','whole_water_residual_m','cumulative_cut_salt_m','cumulative_cut_adv_m','cumulative_cut_diff_m','cut_salt_total_down_m_s','cut_salt_adv_down_m_s','cut_salt_diff_down_m_s']:r[k]=0
                budget.append(r)
        for name,rr in [('surface_events.csv',ev),('transport_budget.csv',budget)]:
            with (d/name).open('w',newline='') as f:
                w=csv.DictWriter(f,fieldnames=list(rr[0]));w.writeheader();w.writerows(rr)
        return d
    def exercise(self,t0,e0,t1,e1,expected):
        dirs=dict(M0=self.bundle('M0',t0,e0),M1=self.bundle('M1',t1,e1))
        self.analysis.OUT=self.root/'output'
        with contextlib.redirect_stdout(io.StringIO()):self.analysis.main(dirs)
        result=json.loads((self.analysis.OUT/'async_stopping_checks.json').read_text())
        self.assertEqual(result['pair']['shared_times_s'],expected)
        self.assertTrue(all(v['valid'] for v in result['legs'].values()))
    def test_m0_first_full_pipeline(self):self.exercise((0,1,2),'surface',(0,1,2,3),'surface',[0,1])
    def test_m1_first_full_pipeline(self):self.exercise((0,1,2,3),'surface',(0,1,2),'surface',[0,1])
    def test_both_censored_full_pipeline(self):self.exercise((0,1,2,3),None,(0,1,2),None,[0,1,2])
    def test_different_events_full_pipeline(self):self.exercise((0,1,2),'interior',(0,1,2,3),'surface',[0,1])
    def test_empty_prethreshold_full_pipeline(self):self.exercise((0,),'surface',(0,1,2),None,[])


if __name__=='__main__':unittest.main()
