"""Slurm entry point for one WP7 run; solver remains one rank and one core."""
import os,sys
from wp7_worker import main
if __name__=='__main__':
    allowed=sorted(os.sched_getaffinity(0))
    if 'SLURM_JOB_ID' not in os.environ or not allowed:raise RuntimeError('WP7 sal worker requires a Slurm allocation')
    os.environ['WP6_MANAGED']='1';os.environ['WP6_CPU']=str(allowed[0])
    raise SystemExit(main(sys.argv[1],sys.argv[2]))
