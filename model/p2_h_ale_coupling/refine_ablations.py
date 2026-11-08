"""Temporal replay check; same initial state and exact recorded flux intervals."""
import json,os,pathlib,subprocess
assert os.environ.get('SLURM_JOB_ID')
root=pathlib.Path.cwd()/'ablations'
for name in ('eq0_n320_native','eq0_n320_replay_total'):
    source=root/name;target=root/(name+'_half');target.mkdir(exist_ok=False)
    args=json.loads((source/'command.json').read_text());args[-1]='.0125'
    (target/'command.json').write_text(json.dumps(args,indent=2))
    with (target/'solver.out').open('w') as f:
        result=subprocess.run(args,cwd=target,stdout=f,stderr=subprocess.STDOUT,timeout=600)
    (target/'exit_status.txt').write_text(str(result.returncode))
    print(name,result.returncode,flush=True)
