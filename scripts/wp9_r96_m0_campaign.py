"""The paired R96-M0 leg: freeze, stage, submit, watch and collect it on hpcc.

M0 is the unperturbed companion of the R96-M1 leg already running on sal. Every other
resolution in the campaign has such a pair (R24_M0/M1, R48_M0/M1); R96 did not, which is
why the registered short-wave discriminator had to say that the R96 numerical floor was
unmeasured. The input differs from the frozen M1 input in exactly the five keys that
define the perturbation, the same five that separate M0 from M1 at R24 and R48.

Host and thread count were measured not to be protocol differences
(hpcc_portability/portability.json and thread_scaling.json), so this leg stays
comparable with R24, R48 and the R96-M1 leg.

Actions: freeze | stage | submit | status | collect | reduce
"""
import argparse,json,shlex,subprocess
from pathlib import Path
import wp9_common as c

OUT=c.OUT/'r96_m0'
PENDING=OUT/'pending'
NAME='WP9-R96-M0-4T'
M1=c.OUT/'pending_r96'
REMOTE='/WORK/hulm_work/zhixin/paper2_wp9_r96_m0_20261104'
PORTREMOTE='/WORK/hulm_work/zhixin/paper2_wp9_portability_20261104'
VIEW=f'{NAME}-frames'
ELF_SHA='112076759fdff0942000307f729c2bd3c7e96963229b43c1b67c1f12bbf9df81'
SLURM='export PATH=/rmprog/slurm/v22.05.7/bin:$PATH'
THREADS=4
STAGED=['wp9_hpcc_worker.py','wp9_hpcc_m0_job.sbatch','wp9_post.py','wp9_pair_frames.py',
        'single_wavelength_post.py','wp6_4t_metrics.py','wp7_spatial.py','wp8_fields.py']
SMALL=['params.input','plan.json','started.json','solver.log','solver_finished.json',
       'stop_record.json','finished.json','surface_events.csv','transport_budget.csv',
       'transport_metadata.json','weighted_online.csv','post.log','eigenprofile.dat',
       f'{NAME}_sigmatop.csv',f'{NAME}.pvd',
       'post/field_audit.json','post/field_hashes.json','post/independent_checks.json',
       'post/independent_inventories.csv','post/weighted_diagnostics.csv',
       'post/profiles.npz','post/representative_fields.npz','failure.json']
# the five keys that define the unperturbed leg, exactly as at R24 and R48
PERTURBATION={'Problem.Name':NAME,'Problem.Amplitude':'0',
              'Problem.StandardDeviationSaltMoleFraction':'0',
              'Problem.PerturbationType':'Constant'}
DROP=['Problem.EigenprofilePath']


def ssh(cmd,timeout=3600):
    # the inner command carries its own quotes, so it is shell-quoted as one argument
    q=subprocess.run(['ssh','-o','ConnectTimeout=25','-o','BatchMode=yes','hpcc',
                      'bash -lc '+shlex.quote(f'{SLURM}; {cmd}')],
                     capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=timeout)
    if q.returncode:raise RuntimeError((q.stderr or '')+(q.stdout or ''))
    return q.stdout


def ssh_json(cmd,timeout=3600):
    """The login shell prints a quota banner before anything else, so take the last
    line rather than the whole stream."""
    out=[l for l in ssh(cmd,timeout).splitlines() if l.strip()]
    return json.loads(out[-1])


def scp(src,dst,timeout=3600):
    subprocess.run(['scp','-o','ConnectTimeout=25','-o','BatchMode=yes',str(src),dst],
                   check=True,timeout=timeout)


def freeze():
    import wp9_r96_protocol as P
    PENDING.mkdir(parents=True,exist_ok=True)
    text=[];changed={}
    for e in P.read_input(M1/'params.input'):
        if e[0]=='raw':text.append(e[3] if len(e)>3 else e[1]);continue
        key,value=e[1],e[2]
        if key in DROP:changed[key]=(value,None);continue
        if key in PERTURBATION:
            changed[key]=(value,PERTURBATION[key]);value=PERTURBATION[key]
        text.append(f'{key.split(".",1)[1]} = {value}')
    (PENDING/'params.input').write_text('\n'.join(text)+'\n',encoding='utf-8')
    assert set(changed)==set(PERTURBATION)|set(DROP),f'unexpected changes {sorted(changed)}'
    import shutil;shutil.copyfile(M1/'eigenprofile.dat',PENDING/'eigenprofile.dat')
    # the reference pair: the same five differences separate M0 from M1 at R48
    r48=lambda role:P.flat(c.ROOT/f'analysis/wp7_shortwave_validation_20261103/runs/WP7-R48-{role}-4T/params.input')
    a,b=r48('M0'),r48('M1')
    reference={k:(b.get(k),a.get(k)) for k in set(a)|set(b) if a.get(k)!=b.get(k)}
    assert set(reference)==set(PERTURBATION)|set(DROP),'the M0/M1 difference is not the R48 one'

    m1=json.loads((c.OUT/'r96_protocol.json').read_text(encoding='utf-8'))
    ts=json.loads((c.OUT/'hpcc_portability/thread_scaling.json').read_text(encoding='utf-8'))
    port=json.loads((c.OUT/'hpcc_portability/portability.json').read_text(encoding='utf-8'))
    proto=dict(schema_version=1,id='WP9-R96-M0',created_date='2026-11-04',
        authority='command/16.md for the comparison, user instruction of 2026-11-04 for this leg',
        authority_sha256=c.sha(c.ROOT/'command/16.md'),
        decision='FROZEN_AND_PENDING_LAUNCH',launched=False,
        scope=['exactly one unperturbed R96 leg, strictly paired with the R96-M1 leg'],
        purpose=('measure the R96 numerical floor, which the registered short-wave discriminator '
                 'records as unmeasured, and complete the A/B/C comparison at R96'),
        forbidden=['R192','multiwavelength','full-domain','8T','restart','retry',
                   'any change of the solver formulation'],
        pairing=dict(partner='WP9-R96-M1-4T',
            partner_input_sha256=m1['numerics'] and c.sha(M1/'params.input'),
            differences={k:dict(m1=v[0],m0=v[1]) for k,v in changed.items()},
            identical_to_the_R48_pair_difference=True,
            same_grid=True,same_dt=True,same_time_axis=True,same_output_interval=True),
        physics=dict(m1['physics'],perturbation=dict(type='Constant',amplitude=0.0,
            note='no imposed perturbation; the leg measures what the discretization alone produces')),
        numerics=m1['numerics'],
        execution=dict(host='hpcc',partition='cnmix',ranks=1,cpus=THREADS,threads=THREADS,
            thread_backend='TBB, compiled into the frozen executable',
            hard_wall_s=21600,projected_solver_wall_s=9600*0.842,
            projection_source='the four-thread ten-step measurement, 0.842 s per step',
            minimum_free_disk_bytes=100*1024**3,minimum_available_memory_bytes=16*1024**3,
            scratch=REMOTE,filesystem='Lustre /WORK',
            io_caveat=('the sal legs wrote to a RAM-backed tmpfs; 4801 files on Lustre may be '
                       'slower, so the projection is a floor, not a promise')),
        host_and_threads_are_not_protocol_differences=dict(
            executable_identical=port['elf_matches_sal'],
            diagnostics_bitwise_identical_to_sal=port['numerically_identical'],
            diagnostics_bitwise_identical_across_thread_counts=ts['all_bitwise_identical'],
            evidence=['hpcc_portability/portability.json','hpcc_portability/thread_scaling.json'],
            rule='neither the host nor the thread count enters any protocol hash'),
        acceptance=m1['physical_acceptance'],
        postprocessing=m1['postprocessing'],
        known_confounds=m1['known_confounds'],
        failure_policy=m1['failure_policy'])
    c.save(OUT/'protocol.json',proto)
    adm=dict(schema_version=1,protocol_sha256=c.sha(OUT/'protocol.json'),
        input_sha256=c.sha(PENDING/'params.input'),
        eigenprofile_sha256=c.sha(PENDING/'eigenprofile.dat'),
        reference_input='analysis/wp9_r96_convergence_20261104/pending_r96/params.input',
        reference_input_sha256=c.sha(M1/'params.input'),
        input_key_changes={k:dict(frm=v[0],to=v[1]) for k,v in changed.items()},
        r48_pair_key_changes={k:dict(frm=v[0],to=v[1]) for k,v in reference.items()},
        status='PENDING')
    c.save(OUT/'admission.json',adm)
    auth=dict(schema_version=1,id='WP9-A02',created_date='2026-11-04',
        question='多线程跑哪一条？',answer='R96-M0（推荐）',
        user_words='嗯嗯，那就在hpcc上跑一个多线程',
        granted=dict(runs=1,role='M0',nx_cells=96,host='hpcc',threads=THREADS,
                     dt_s=proto['numerics']['dt_s'],steps=proto['numerics']['steps']),
        not_granted=['R192','multiwavelength','full-domain','8T','restart','retry'],
        protocol_sha256=adm['protocol_sha256'],input_sha256=adm['input_sha256'])
    c.save(OUT/'user_authorization.json',auth)
    print(json.dumps(dict(input_sha256=adm['input_sha256'],
                          changes=adm['input_key_changes'],
                          matches_R48_pair=set(adm['input_key_changes'])==set(adm['r48_pair_key_changes'])),
                     indent=2))


def plan():
    proto=json.loads((OUT/'protocol.json').read_text(encoding='utf-8'))
    adm=json.loads((OUT/'admission.json').read_text(encoding='utf-8'))
    return dict(id=NAME,role='M0',stage='r96_m0',case_id=int(c.CASE),
        nx_cells=96,nz_cells=200,T_s=c.HORIZON_S,dt_s=proto['numerics']['dt_s'],
        steps=proto['numerics']['steps'],output_interval=proto['numerics']['output_interval'],
        input_sha256=adm['input_sha256'],profile_sha256=adm['eigenprofile_sha256'],
        elf_sha256=ELF_SHA,preregistration_sha256=adm['protocol_sha256'],
        user_authorization_sha256=c.sha(OUT/'user_authorization.json'),
        hard_wall_limit_s=proto['execution']['hard_wall_s'],
        minimum_free_disk_bytes=proto['execution']['minimum_free_disk_bytes'],
        minimum_available_memory_bytes=proto['execution']['minimum_available_memory_bytes'],
        threads=THREADS,host='hpcc')


PROBE=r'''
import hashlib,json,os
from pathlib import Path
root=Path("%s");b=root/"bundle"
man=json.loads((b/"runtime_manifest.json").read_text())
bad=[n for n,v in man.items() if hashlib.sha256((b/n).read_bytes()).hexdigest()!=v["sha256"]]
st=os.statvfs(root)
avail=int(next(s.split()[1] for s in Path("/proc/meminfo").read_text().splitlines() if s.startswith("MemAvailable:")))*1024
print(json.dumps(dict(root=str(root),
  elf_sha256=hashlib.sha256((b/"bin/paper2_transport").read_bytes()).hexdigest(),
  runtime_manifest_mismatches=bad,
  input_sha256=hashlib.sha256((root/"runs/%s/params.input").read_bytes()).hexdigest(),
  profile_sha256=hashlib.sha256((root/"runs/%s/eigenprofile.dat").read_bytes()).hexdigest(),
  preregistration_sha256=hashlib.sha256((root/"preregistration.json").read_bytes()).hexdigest(),
  authorization_sha256=hashlib.sha256((root/"user_authorization.json").read_bytes()).hexdigest(),
  free_disk_bytes=st.f_bavail*st.f_frsize,available_memory_bytes=avail)))
'''


def stage():
    assert not (OUT/'stage.json').exists(),'one authorized run only'
    p=plan();c.save(PENDING/'plan.json',p)
    ssh(f'test ! -e {REMOTE} && mkdir -p {REMOTE}/runs/{NAME}')
    ssh(f'cp -a {PORTREMOTE}/bundle {REMOTE}/bundle')
    for n in STAGED:scp(c.ROOT/'scripts'/n,f'hpcc:{REMOTE}/{n}')
    scp(OUT/'protocol.json',f'hpcc:{REMOTE}/preregistration.json')
    scp(OUT/'user_authorization.json',f'hpcc:{REMOTE}/user_authorization.json')
    for n in ['params.input','eigenprofile.dat','plan.json']:
        scp(PENDING/n,f'hpcc:{REMOTE}/runs/{NAME}/{n}')
    probe=ssh_json(f"python3 - <<'PY'\n{PROBE%(REMOTE,NAME,NAME)}\nPY")
    checks=dict(
        elf_matches_the_campaign=probe['elf_sha256']==ELF_SHA,
        runtime_manifest_clean=not probe['runtime_manifest_mismatches'],
        input_matches=probe['input_sha256']==p['input_sha256'],
        profile_matches=probe['profile_sha256']==p['profile_sha256'],
        protocol_matches=probe['preregistration_sha256']==p['preregistration_sha256'],
        authorization_matches=probe['authorization_sha256']==p['user_authorization_sha256'],
        free_disk_above_minimum=probe['free_disk_bytes']>p['minimum_free_disk_bytes'],
        memory_above_minimum=probe['available_memory_bytes']>p['minimum_available_memory_bytes'])
    probe.update(checks=checks,admitted=all(checks.values()),plan=p)
    c.save(OUT/'stage.json',probe)
    print(json.dumps(dict(admitted=probe['admitted'],checks=checks,
                          free_disk_gib=probe['free_disk_bytes']/1024**3,
                          available_memory_gib=probe['available_memory_bytes']/1024**3),indent=2))
    assert probe['admitted'],'admission failed; nothing submitted'


def submit():
    assert json.loads((OUT/'stage.json').read_text(encoding='utf-8'))['admitted']
    assert not (OUT/'job.json').exists(),'already submitted'
    job=ssh(f'cd {REMOTE} && sbatch --parsable wp9_hpcc_m0_job.sbatch {REMOTE} {NAME}'
            ).strip().splitlines()[-1]
    c.save(OUT/'job.json',dict(job_id=job,run=NAME,remote=REMOTE,host='hpcc',threads=THREADS))
    print(json.dumps(dict(job_id=job),indent=2))


STATUS=r'''
import json,os
from pathlib import Path
r=Path("%s/runs/%s")
rows=0
if (r/"surface_events.csv").exists():
    rows=sum(1 for _ in (r/"surface_events.csv").open())-1
tail=(r/"solver.log").read_text().splitlines()[-3:] if (r/"solver.log").exists() else None
print(json.dumps(dict(steps_recorded=rows,vtu=len(list(r.glob("*.vtu"))),
  started=(r/"started.json").exists(),solver_finished=(r/"solver_finished.json").exists(),
  finished=(r/"finished.json").exists(),failure=(r/"failure.json").exists(),
  stop=json.loads((r/"stop_record.json").read_text()) if (r/"stop_record.json").exists() else None,
  tail=tail)))
'''


def status():
    job=json.loads((OUT/'job.json').read_text(encoding='utf-8'))['job_id']
    q=ssh(f'squeue -j {job} -h -o "%T %M %L %R" || true')
    st=ssh_json(f"python3 - <<'PY'\n{STATUS%(REMOTE,NAME)}\nPY")
    st['squeue']=q.strip().splitlines()[-1] if q.strip() else ''
    st['progress']=st['steps_recorded']/9600
    print(json.dumps(st,indent=2,default=str))
    return st


def collect():
    run=OUT/'run';run.mkdir(parents=True,exist_ok=True)
    known=json.loads((run/'finished.json').read_text())['files'] if (run/'finished.json').exists() else {}
    for n in SMALL:
        dst=run/n
        if n in known and dst.exists() and c.sha(dst)==known[n]['sha256']:continue
        dst.parent.mkdir(parents=True,exist_ok=True)
        q=subprocess.run(['scp','-o','ConnectTimeout=25','-o','BatchMode=yes',
                          f'hpcc:{REMOTE}/runs/{NAME}/{n}',str(dst)],
                         capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=3600)
        if q.returncode and n!='failure.json':print('missing',n)
    fin=json.loads((run/'finished.json').read_text())
    pulled=[n for n in fin['files'] if (run/n).exists()]
    bad=[n for n in pulled if c.sha(run/n)!=fin['files'][n]['sha256']]
    assert not bad,f'transfer hash mismatch: {bad}'
    rec=ssh_json(f"python3 - <<'PY'\nimport hashlib,json\nfrom pathlib import Path\n"
        f"a=Path('{REMOTE}/{NAME}_fields.tar.gz')\n"
        f"m=json.loads(Path('{REMOTE}/runs/{NAME}/post/field_hashes.json').read_text())\n"
        f"print(json.dumps(dict(archive=str(a),bytes=a.stat().st_size,"
        f"sha256=hashlib.sha256(a.read_bytes()).hexdigest(),frames=len(m))))\nPY",timeout=3600)
    rec.update(archive_matches_finished=rec['sha256']==fin['archive']['sha256'],
               field_manifest_sha256=c.sha(run/'post/field_hashes.json'),
               independent_checks=json.loads((run/'post/independent_checks.json').read_text()),
               retained_remotely_because='the workstation has little free disk; only derived products are pulled',
               local_files={str(n.relative_to(run)).replace('\\','/'):dict(sha256=c.sha(n),bytes=n.stat().st_size)
                            for n in sorted(run.rglob('*')) if n.is_file()})
    c.save(OUT/'archive_receipt.json',rec)
    print(json.dumps(dict(frames=rec['frames'],bytes=rec['bytes'],
                          archive_matches_finished=rec['archive_matches_finished'],
                          independent_passed=rec['independent_checks']['passed']),indent=2))


def reduce():
    interval=json.loads((OUT/'protocol.json').read_text(encoding='utf-8'))['numerics']['output_interval']
    pairing=ssh_json(f'cd {REMOTE} && python3 wp9_pair_frames.py runs/{NAME} runs/{VIEW} '
                     f'{interval}',timeout=3600)
    assert pairing['sha256_values_unchanged']
    c.save(OUT/'frame_pairing_receipt.json',pairing)
    outs={}
    for script,sub in [('wp8_fields.py','reduced'),('wp7_spatial.py','spatial')]:
        log=ssh(f'cd {REMOTE} && env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 '
                f'python3 {script} runs/{VIEW} {sub}_{NAME} 2>&1 | tail -3',timeout=7200)
        outs[script]=log.strip().splitlines()[-1]
    dest=OUT/'fields';dest.mkdir(parents=True,exist_ok=True)
    for remote,local in [(f'reduced_{NAME}/timeseries.csv','reduced_timeseries.csv'),
                         (f'reduced_{NAME}/audit.json','reduced_audit.json'),
                         (f'spatial_{NAME}/spatial_metrics.csv','spatial_metrics.csv'),
                         (f'spatial_{NAME}/spatial_audit.json','spatial_spatial_audit.json'),
                         (f'spatial_{NAME}/selected_vertical_profiles.csv','spatial_selected_vertical_profiles.csv')]:
        scp(f'hpcc:{REMOTE}/{remote}',str(dest/local))
    c.CACHE.mkdir(parents=True,exist_ok=True)
    scp(f'hpcc:{REMOTE}/reduced_{NAME}/profiles_spectra.npz',str(c.CACHE/'R96_M0_profiles_spectra.npz'))
    audit=json.loads((dest/'reduced_audit.json').read_text())
    assert audit['all_field_hashes_verified']
    assert audit['files']['profiles_spectra.npz']==c.sha(c.CACHE/'R96_M0_profiles_spectra.npz')
    c.save(OUT/'reduction_receipt.json',dict(remote_logs=outs,reduced_audit=audit,
        spatial_audit=json.loads((dest/'spatial_spatial_audit.json').read_text()),
        local={p.name:dict(sha256=c.sha(p),bytes=p.stat().st_size)
               for p in sorted(dest.iterdir()) if p.is_file()},
        spectra_npz=dict(path=str(c.CACHE/'R96_M0_profiles_spectra.npz'),
                         sha256=c.sha(c.CACHE/'R96_M0_profiles_spectra.npz'))))
    print(json.dumps(dict(frames=audit['frames'],verified=audit['all_field_hashes_verified']),indent=2))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('action',
        choices=['freeze','stage','submit','status','collect','reduce'])
    globals()[ap.parse_args().action]()
