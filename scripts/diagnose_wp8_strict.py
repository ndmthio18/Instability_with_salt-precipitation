"""Read-only failure localization. No solver tolerance edits or resubmissions."""
import re,json
from wp8_strict_campaign import OUT,T,save
from analyze_wp8 import rows
from single_wavelength_post import savecsv
def main():
    log=(OUT/'4T/solver.log').read_text();ev=rows(OUT/'4T/surface_events.csv');trace=[];step=0;block=[];hist={}
    for line in log.splitlines():
        m=re.search(r'WP8_ITER (\d+) water_residual=([^ ]+) salt_residual=([^ ]+) shift=([^ ]+)',line)
        if m:
            it=int(m[1])
            if it==0:step+=1;block=[]
            row=dict(attempted_step=step,iteration=it,water_residual=float(m[2]),salt_residual=float(m[3]),relative_shift=float(m[4]),accepted=step<len(ev));trace.append(row);block.append(row)
        m=re.search(r'WP8_CONSTRAINT iterations=(\d+)',line)
        if m:hist[m[1]]=hist.get(m[1],0)+1
    failed=[r for r in trace if not r['accepted']];late=[r for r in failed if r['iteration']>=2]
    savecsv(OUT/'newton_trace.csv',trace);savecsv(OUT/'failed_step_iterations.csv',failed)
    result=dict(last_accepted_s=float(ev[-1]['time_s']),last_accepted_over_T=float(ev[-1]['time_s'])/T,accepted_steps=len(ev)-1,failed_attempt_step=step,proposed_unaccepted_time_s=float(ev[-1]['time_s'])+float(ev[-1]['dt_s']),newton_success_histogram=hist,failed_attempt_iterations=len(failed),late_failed_salt_min=min(r['salt_residual'] for r in late),late_failed_salt_max=max(r['salt_residual'] for r in late),late_failed_water_max=max(r['water_residual'] for r in late),late_failed_shift_max=max(r['relative_shift'] for r in late),fixed_internal_limits=dict(water=2e-18,salt=2e-20,shift=1e-10),interpretation='salt algebraic residual stagnation just above fixed cutoff; compatible with rounding-floor limitation, not proof of physical divergence or permission to relax gates',automatic_retry=False,source_changed=False)
    save(OUT/'failure_diagnosis.json',result);print(json.dumps(result,indent=2))
if __name__=='__main__':main()
