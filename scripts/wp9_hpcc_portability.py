"""Portability check: does the frozen sal bundle reproduce itself on hpcc?

Scope: one ten-step R96 benchmark, identical input and identical executable, run as a
single-core Slurm job on the 探索1000 cluster. It is an infrastructure check, not a
scientific result: no new physics, no new campaign, and nothing here authorizes a
production leg. The success criterion is that the diagnostics match the sal benchmark
that was already recorded, so a later leg could be placed on either host without a
host-dependent protocol difference.

Actions: stage | submit | status | collect
"""
import argparse,json,subprocess
from pathlib import Path
import wp9_common as c

OUT=c.OUT/'hpcc_portability'
REMOTE='/WORK/hulm_work/zhixin/paper2_wp9_portability_20261104'
BENCH=c.ROOT/'analysis/wp7_shortwave_validation_20261103/runs/WP7-R96-M1-bench10'
BUNDLE=Path('C:/Users/12608/AppData/Local/Temp/claude/8ac5c3d6dacf0260/'
            '9089512e-4ea2-4e2c-8139-aa64039ad641/scratchpad/bundle_wp9.tar.gz')
BUNDLE_SHA='af885bc127bce044af72989347209aed35a4bd0fbba2f6722dc836c681967b01'
ELF_SHA='112076759fdff0942000307f729c2bd3c7e96963229b43c1b67c1f12bbf9df81'
SLURM='export PATH=/rmprog/slurm/v22.05.7/bin:$PATH'
COMPARE=['surface_events.csv','transport_budget.csv','weighted_online.csv',
         'WP7-R96-M1-bench10_sigmatop.csv']


def ssh(cmd,timeout=1800):
    q=subprocess.run(['ssh','-o','ConnectTimeout=25','-o','BatchMode=yes','hpcc',
                      f"bash -lc '{SLURM}; {cmd}'"],capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=timeout)
    if q.returncode:raise RuntimeError((q.stderr or '')+(q.stdout or ''))
    return q.stdout


def scp(src,dst,timeout=1800):
    subprocess.run(['scp','-o','ConnectTimeout=25','-o','BatchMode=yes',str(src),dst],
                   check=True,timeout=timeout)


def stage():
    OUT.mkdir(parents=True,exist_ok=True)
    assert c.sha(BUNDLE)==BUNDLE_SHA,'the staged bundle is not the one taken from sal'
    ssh(f'mkdir -p {REMOTE}')
    scp(BUNDLE,f'hpcc:{REMOTE}/bundle.tar.gz')
    for n in ['params.input','eigenprofile.dat']:scp(BENCH/n,f'hpcc:{REMOTE}/{n}')
    scp(c.ROOT/'scripts/wp9_hpcc_job.sbatch',f'hpcc:{REMOTE}/job.sbatch')
    probe=ssh(f'cd {REMOTE} && tar -xzf bundle.tar.gz && rm -f bundle.tar.gz && '
              f'bash -n job.sbatch && sha256sum bundle/bin/paper2_transport params.input '
              f'eigenprofile.dat job.sbatch && find bundle -type f | wc -l && quota -s 2>&1 | tail -2')
    lines=[l for l in probe.strip().splitlines() if l.strip()]
    sums={}
    count=None
    for l in lines:
        parts=l.split()
        if len(parts)==2 and len(parts[0])==64:sums[parts[1].lstrip('*')]=parts[0]
        elif l.strip().isdigit():count=int(l.strip())
    rec=dict(host='hpcc',login_node='ibcln01',remote=REMOTE,
             bundle_sha256=BUNDLE_SHA,
             bundle_source='the live sal WP9 run root, tarred after the orted relocation',
             checksums=sums,bundle_files=count,
             elf_matches_sal=sums.get('bundle/bin/paper2_transport')==ELF_SHA,
             input_matches_bench=sums.get('params.input')==c.sha(BENCH/'params.input'),
             profile_matches_bench=sums.get('eigenprofile.dat')==c.sha(BENCH/'eigenprofile.dat'),
             quota_tail=[l for l in lines if '/WORK' in l or 'balance' in l],
             scope='infrastructure portability only; authorizes no production leg')
    assert rec['elf_matches_sal'] and rec['input_matches_bench'] and rec['profile_matches_bench']
    c.save(OUT/'stage.json',rec)
    print(json.dumps({k:rec[k] for k in ['elf_matches_sal','input_matches_bench',
                                         'profile_matches_bench','bundle_files']},indent=2))


def submit():
    assert json.loads((OUT/'stage.json').read_text(encoding='utf-8'))['elf_matches_sal']
    assert not (OUT/'job.json').exists(),'already submitted'
    job=ssh(f'cd {REMOTE} && sbatch --parsable job.sbatch').strip().splitlines()[-1]
    c.save(OUT/'job.json',dict(job_id=job,remote=REMOTE,partition='cnmix',cpus=1,
                               scope='ten-step R96 benchmark, portability only'))
    print(json.dumps(dict(job_id=job),indent=2))


def status():
    job=json.loads((OUT/'job.json').read_text(encoding='utf-8'))['job_id']
    q=ssh(f'squeue -j {job} -h -o "%T %M %R" || true')
    s=ssh(f'sacct -j {job} -X -n -P -o JobID,State,ExitCode,Elapsed,NodeList,MaxRSS || true')
    tail=ssh(f'cd {REMOTE} && tail -3 solver.log 2>/dev/null || true')
    print(json.dumps(dict(squeue=q.strip(),sacct=s.strip(),
                          tail=tail.strip().splitlines()[-2:]),indent=2))


def collect():
    for n in ['environment.txt','modules.txt','solver.log','stop_record.json',
              'transport_metadata.json']+COMPARE:
        q=subprocess.run(['scp','-o','ConnectTimeout=25','-o','BatchMode=yes',
                          f'hpcc:{REMOTE}/{n}',str(OUT/n)],capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=1800)
        if q.returncode:print('missing',n)
    job=json.loads((OUT/'job.json').read_text(encoding='utf-8'))['job_id']
    sac=ssh(f'sacct -j {job} -X -n -P -o JobID,State,ExitCode,Elapsed,AllocCPUS,NodeList,MaxRSS')
    (OUT/'sacct.txt').write_text(sac,encoding='utf-8')
    env=dict(l.split('=',1) for l in (OUT/'environment.txt').read_text().splitlines() if '=' in l)
    stop=json.loads((OUT/'stop_record.json').read_text())
    comparison={}
    for n in COMPARE:
        a=OUT/n;b=BENCH/n
        if not a.exists() or not b.exists():comparison[n]=dict(present=False);continue
        ta=a.read_text().splitlines();tb=b.read_text().splitlines()
        worst=0.0;cells=0;mismatched=0
        for la,lb in zip(ta,tb):
            fa=la.lstrip('#').replace(',',' ').split();fb=lb.lstrip('#').replace(',',' ').split()
            for x,y in zip(fa,fb):
                try:fx,fy=float(x),float(y)
                except ValueError:continue
                cells+=1
                if fx!=fy:
                    mismatched+=1
                    worst=max(worst,abs(fx-fy)/max(abs(fy),1e-300))
        comparison[n]=dict(present=True,bitwise_identical=c.sha(a)==c.sha(b),
                           lines=len(ta),line_count_matches=len(ta)==len(tb),
                           numeric_cells=cells,differing_cells=mismatched,
                           max_relative_difference=worst)
    rec=dict(host='hpcc',node=env.get('hostname'),kernel=env.get('kernel'),
             host_glibc=env.get('host_glibc'),python=env.get('python'),
             returncode=int(env.get('returncode','1')),
             solver_wall_s=float(env.get('solver_wall_s','nan')),
             vtu_files=int(env.get('vtu_files','0')),
             elf_sha256=env.get('elf_sha256'),elf_matches_sal=env.get('elf_sha256')==ELF_SHA,
             loader_sha256=env.get('loader_sha256'),libc_sha256=env.get('libc_sha256'),
             stop=stop,sacct=sac.strip(),comparison=comparison,
             sal_benchmark=dict(run=str(BENCH.relative_to(c.ROOT)),
                 solver_wall_s=json.loads((BENCH/'finished.json').read_text())['solver']['wall_s']))
    rec['speedup_versus_sal']=(rec['sal_benchmark']['solver_wall_s']/rec['solver_wall_s']
                               if rec['solver_wall_s']>0 else None)
    rec['portable']=bool(rec['returncode']==0 and rec['elf_matches_sal']
                         and stop['reason']=='requested_horizon' and stop['accepted_steps']==10)
    rec['numerically_identical']=all(v.get('bitwise_identical') for v in comparison.values()
                                     if v.get('present'))
    rec['conclusion']=('the frozen bundle runs unchanged on CentOS 7 and reproduces the sal '
                       'benchmark bit for bit' if rec['portable'] and rec['numerically_identical']
                       else 'a host-dependent difference exists and must be registered before any '
                            'leg is placed on this cluster; see comparison')
    rec['scope']='infrastructure only; no scientific result and no production authorization'
    c.save(OUT/'portability.json',rec)
    print(json.dumps({k:rec[k] for k in ['portable','numerically_identical','node','returncode',
                                         'solver_wall_s','speedup_versus_sal','vtu_files',
                                         'conclusion']},indent=2,default=str))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('action',
        choices=['stage','submit','status','collect'])
    globals()[ap.parse_args().action]()
