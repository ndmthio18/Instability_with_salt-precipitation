"""Convert the triggered R96 benchmark into a bounded planning estimate."""
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'analysis/wp7_shortwave_validation_20261103';NAME='WP7-R96-M1-bench10'
def main():
    d=json.loads((OUT/'runs'/NAME/'finished.json').read_text());s=d['solver'];assert s['stop']['accepted_steps']==10 and s['returncode']==0
    r48=[]
    for role in ['M0','M1']:r48.append(json.loads((OUT/'runs'/f'WP7-R48-{role}-4T/finished.json').read_text()))
    seconds_per_step=s['wall_s']/10;solver=seconds_per_step*4800;post=max(x['post_archive_wall_s'] for x in r48)*2*1.25;worker=solver+post;result=dict(benchmark_solver_wall_s=s['wall_s'],seconds_per_step=seconds_per_step,projected_solver_wall_s_per_leg=solver,projected_post_archive_wall_s_per_leg=post,projected_worker_wall_s_per_leg=worker,projected_pair_core_h=2*worker/3600,method='R96 10-step linear solver extrapolation plus twice measured worst R48 post/archive with 25 percent margin',authorization='cost evidence only; no full R96 launch')
    (OUT/'r96_cost_estimate.json').write_text(json.dumps(result,indent=2)+'\n')
    comparison=json.loads((OUT/'r48_comparison.json').read_text());decision=comparison['decision']
    decision.update(R96_cost_benchmark_completed=True,R96_projected_pair_core_h=result['projected_pair_core_h'],R96_projected_worker_wall_s_per_leg=result['projected_worker_wall_s_per_leg'],R96='PROPOSE_ONLY_NO_GO; full run exceeds current 2 h per-leg cap')
    comparison['decision']=decision
    (OUT/'r48_comparison.json').write_text(json.dumps(comparison,indent=2)+'\n')
    (OUT/'decision.json').write_text(json.dumps(decision,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
