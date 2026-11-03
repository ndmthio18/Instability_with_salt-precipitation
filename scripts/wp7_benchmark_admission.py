"""Admit the R48 pair only from measured benchmark cost and frozen limits."""
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'analysis/wp7_shortwave_validation_20261103'
def main():
    d=OUT/'runs/WP7-R48-M1-bench10';f=json.loads((d/'finished.json').read_text());s=f['solver'];steps=s['stop']['accepted_steps'];assert steps==10 and s['returncode']==0
    seconds_per_step=s['wall_s']/steps;projected_solver=seconds_per_step*4800
    # R24 post/archive was ~689 s; scale by R48/R24 field size with 25% margin.
    projected_post=689.193*2*1.25;projected_worker=projected_solver+projected_post
    pair_core_h=2*projected_worker/3600;admitted=projected_worker<=7200 and pair_core_h<=4.373333333333333
    result=dict(admitted=admitted,benchmark_steps=steps,benchmark_solver_wall_s=s['wall_s'],seconds_per_step=seconds_per_step,
      projected_solver_wall_s_per_leg=projected_solver,projected_post_archive_wall_s_per_leg=projected_post,
      projected_worker_wall_s_per_leg=projected_worker,projected_pair_core_h=pair_core_h,
      gates=dict(worker_limit_s=7200,pair_planning_cap_core_h=4.373333333333333),method='linear 10-step solver extrapolation plus 2x R24 post/archive with 25% margin; admission estimate, not confidence bound')
    (OUT/'benchmark_admission.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2));assert admitted
if __name__=='__main__':main()
