import unittest
import json,copy
from analyze_wp8_strict import persistent_crossing
from wp6_paired_stopping import same_time
from analyze_wp8 import pairs
from analyze_wp8 import rows
from wp8_strict_campaign import OUT
from validate_wp8_strict import accepted_prefix_audit

class TestStrictDiagnostics(unittest.TestCase):
    def test_isolated_noise_not_persistent_departure(self):
        r=persistent_crossing(range(8),[0,2,0,2,0,2,0,2],1,n=3)
        self.assertEqual(r['status'],'not_observed')
    def test_departure_bracket_and_confirmation(self):
        r=persistent_crossing([0,10,20,30,40],[0,.5,2,3,4],1,n=3)
        self.assertEqual((r['low_s'],r['high_s'],r['confirmed_s']),(10,20,40))
    def test_exact_shared_window_excludes_clipped_endpoint(self):
        t=136423.0637723484;other=136423.06377235046
        self.assertFalse(same_time(t,other));self.assertEqual(pairs([0,100,t],[0,100,other]),[(0,0),(1,1)])
    def test_saved_prefix_keeps_rejected_step_out(self):
        d=OUT/'4T';ev=rows(d/'surface_events.csv');stop=json.loads((d/'stop_record.json').read_text());log=(d/'solver.log').read_text();self.assertTrue(accepted_prefix_audit(ev,stop,log))
        bad=copy.deepcopy(ev);bad.append(copy.deepcopy(bad[-1]));self.assertFalse(accepted_prefix_audit(bad,stop,log))
        bad=copy.deepcopy(ev);bad[-1]['global_threshold_reached']='1';self.assertFalse(accepted_prefix_audit(bad,stop,log))
if __name__=='__main__':unittest.main()
