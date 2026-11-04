"""Stage, submit, watch and collect the single authorized WP9 R96-M1 leg on sal.

No solver code, no input and no reviewed worker is modified here. The run uses the
bundle that already produced R24 and R48, copied byte-for-byte on the remote host.
"""
import argparse,json,subprocess,sys
from pathlib import Path
import wp9_common as c

OUT=c.OUT
PENDING=OUT/'pending_r96'
NAME='WP9-R96-M1-4T'
REMOTE='/tmp/zhixin_paper2_wp9_r96_20261104'
WP7_REMOTE='/tmp/zhixin_paper2_wp7_shortwave_20261103_r01'
ARCHIVE_HOME='paper2_wp9_r96_20261104'
ELF_SHA='112076759fdff0942000307f729c2bd3c7e96963229b43c1b67c1f12bbf9df81'
STAGED=['wp7_worker.py','wp7_sal_worker.py','wp7_sal_job.sh','wp6_post.py',
        'single_wavelength_post.py','wp6_4t_metrics.py','wp7_spatial.py','wp8_fields.py',
        'wp9_pair_frames.py']
VIEW=f'{NAME}-frames'
SMALL=['params.input','plan.json','started.json','solver.log','solver_finished.json',
       'stop_record.json','finished.json','surface_events.csv','transport_budget.csv',
       'transport_metadata.json','weighted_online.csv','post.log','eigenprofile.dat',
       f'{NAME}_sigmatop.csv',f'{NAME}.pvd',
       'post/field_audit.json','post/field_hashes.json','post/independent_checks.json',
       'post/independent_inventories.csv','post/weighted_diagnostics.csv',
       'post/profiles.npz','post/representative_fields.npz','failure.json']


def ssh(cmd,capture=True,timeout=1800):
    q=subprocess.run(['ssh','-o','ConnectTimeout=20','-o','BatchMode=yes','sal',cmd],
                     capture_output=capture,text=True,timeout=timeout)
    if q.returncode:raise RuntimeError((q.stderr or '')+(q.stdout or ''))
    return q.stdout if capture else ''


def scp(src,dst,timeout=1800):
    subprocess.run(['scp','-o','ConnectTimeout=20','-o','BatchMode=yes',str(src),dst],
                   check=True,timeout=timeout)


def plan():
    auth=OUT/'user_authorization.json'
    assert auth.exists(),'no authorization record; command/16 forbids launching'
    a=json.loads(auth.read_text(encoding='utf-8'))
    proto=json.loads((OUT/'r96_protocol.json').read_text(encoding='utf-8'))
    adm=json.loads((OUT/'r96_admission.json').read_text(encoding='utf-8'))
    assert a['protocol_sha256']==c.sha(OUT/'r96_protocol.json')==adm['protocol_sha256']
    assert a['input_sha256']==c.sha(PENDING/'params.input')==adm['input_sha256']
    assert a['granted']['runs']==1 and a['granted']['role']=='M1'
    return dict(id=NAME,role='M1',stage='r96',case_id=int(c.CASE),
        nx_cells=96,nz_cells=200,T_s=c.HORIZON_S,dt_s=proto['numerics']['dt_s'],
        steps=proto['numerics']['steps'],output_interval=proto['numerics']['output_interval'],
        scalar_unknowns=2*(96+1)*(200+1),
        input_sha256=adm['input_sha256'],profile_sha256=adm['eigenprofile_sha256'],
        elf_sha256=ELF_SHA,preregistration_sha256=adm['protocol_sha256'],
        user_authorization_sha256=c.sha(auth),
        authority_sha256=proto['authority_sha256'],
        hard_wall_limit_s=proto['resources']['hard_wall_s'],
        minimum_free_disk_bytes=proto['resources']['minimum_free_scratch_bytes'],
        runtime_relocation_sha256=c.sha(c.ROOT/'analysis/wp7_shortwave_validation_20261103/runtime_relocation_receipt.json'),
        reference_run='analysis/wp7_shortwave_validation_20261103/runs/WP7-R48-M1-4T')


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
  runtime_manifest_entries=len(man),runtime_manifest_mismatches=bad,
  input_sha256=hashlib.sha256((root/"runs/%s/params.input").read_bytes()).hexdigest(),
  profile_sha256=hashlib.sha256((root/"runs/%s/eigenprofile.dat").read_bytes()).hexdigest(),
  preregistration_sha256=hashlib.sha256((root/"preregistration.json").read_bytes()).hexdigest(),
  free_disk_bytes=st.f_bavail*st.f_frsize,available_memory_bytes=avail,
  load=os.getloadavg())))
'''


def stage():
    assert not (OUT/'stage.json').exists(),'one authorized run only'
    p=plan();c.save(PENDING/'plan.json',p)
    ssh(f'test ! -e {REMOTE} && mkdir -p {REMOTE}/runs/{NAME}')
    ssh(f'cp -a {WP7_REMOTE}/bundle {REMOTE}/bundle')
    for n in STAGED:scp(c.ROOT/'scripts'/n,f'sal:{REMOTE}/{n}')
    scp(OUT/'r96_protocol.json',f'sal:{REMOTE}/preregistration.json')
    scp(OUT/'user_authorization.json',f'sal:{REMOTE}/user_authorization.json')
    for n in ['params.input','eigenprofile.dat','plan.json']:
        scp(PENDING/n,f'sal:{REMOTE}/runs/{NAME}/{n}')
    probe=json.loads(ssh(f"python3 - <<'PY'\n{PROBE%(REMOTE,NAME,NAME,)}\nPY"))
    checks=dict(
        elf_matches_R24_R48=probe['elf_sha256']==ELF_SHA,
        runtime_manifest_clean=not probe['runtime_manifest_mismatches'],
        input_matches=probe['input_sha256']==p['input_sha256'],
        profile_matches=probe['profile_sha256']==p['profile_sha256'],
        protocol_matches=probe['preregistration_sha256']==p['preregistration_sha256'],
        free_scratch_above_minimum=probe['free_disk_bytes']>p['minimum_free_disk_bytes'],
        memory_above_worker_gate=probe['available_memory_bytes']>16*1024**3,
        reduced_spectra_present=all((c.ROOT/'tmp/wp8_zero_mode_20261103'/k/'profiles_spectra.npz').exists()
                                    for k in ['R24_M0','R24_M1','R48_M0','R48_M1']))
    probe.update(checks=checks,admitted=all(checks.values()),plan=p)
    c.save(OUT/'stage.json',probe)
    print(json.dumps(dict(admitted=probe['admitted'],checks=checks,
                          free_disk_gib=probe['free_disk_bytes']/1024**3,
                          available_memory_gib=probe['available_memory_bytes']/1024**3),indent=2))
    assert probe['admitted'],'admission failed; nothing submitted'


def submit():
    assert json.loads((OUT/'stage.json').read_text())['admitted']
    assert not (OUT/'job.json').exists(),'already submitted'
    job=ssh(f'cd {REMOTE} && sbatch --parsable --job-name=WP9-R96 --time=14:00:00 '
            f'wp7_sal_job.sh {REMOTE} {NAME}').strip()
    c.save(OUT/'job.json',dict(job_id=job,run=NAME,remote=REMOTE,
                               submitted_from='wp9_r96_campaign.py'))
    print(json.dumps(dict(job_id=job),indent=2))


STATUS=r'''
import json,os
from pathlib import Path
r=Path("%s/runs/%s")
def head(n,k=3):
    p=r/n
    return p.read_text().splitlines()[-k:] if p.exists() else None
rows=0
if (r/"surface_events.csv").exists():
    rows=sum(1 for _ in (r/"surface_events.csv").open())-1
print(json.dumps(dict(steps_recorded=rows,
  vtu=len(list(r.glob("*.vtu"))),
  started=(r/"started.json").exists(),
  solver_finished=(r/"solver_finished.json").exists(),
  finished=(r/"finished.json").exists(),
  failure=(r/"failure.json").exists(),
  stop=json.loads((r/"stop_record.json").read_text()) if (r/"stop_record.json").exists() else None,
  tail=head("solver.log",4))))
'''


def status():
    q=json.loads((OUT/'job.json').read_text())
    squeue=ssh(f'squeue -j {q["job_id"]} -h -o "%T %M %L" || true')
    st=json.loads(ssh(f"python3 - <<'PY'\n{STATUS%(REMOTE,NAME)}\nPY"))
    st['squeue']=squeue.strip()
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
        q=subprocess.run(['scp','-o','ConnectTimeout=20','-o','BatchMode=yes',
                          f'sal:{REMOTE}/runs/{NAME}/{n}',str(dst)],
                         capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=3600)
        if q.returncode and n!='failure.json':print('missing',n)
    fin=json.loads((run/'finished.json').read_text())
    pulled=[n for n in fin['files'] if (run/n).exists()]
    bad=[n for n in pulled if c.sha(run/n)!=fin['files'][n]['sha256']]
    assert not bad,f'transfer hash mismatch: {bad}'
    assert len(pulled)>=len(SMALL)-2,f'only {len(pulled)} of the recorded files arrived'
    # The field archive stays on sal: the workstation has no room for it. The registered
    # plan was to move it to the NFS home directory so it survives a reboot of the
    # RAM-backed /tmp, but that account's home quota is exhausted, so the move is
    # attempted and its failure recorded rather than hidden.
    moved=True;move_error=None
    try:
        ssh(f'mkdir -p ~/{ARCHIVE_HOME} && mv {REMOTE}/{NAME}_fields.tar.gz ~/{ARCHIVE_HOME}/ '
            f'&& cp {REMOTE}/runs/{NAME}/post/field_hashes.json ~/{ARCHIVE_HOME}/')
    except RuntimeError as exc:
        moved=False;move_error=str(exc).strip() or 'the remote move failed with no message'
        quota=ssh(f'quota -s 2>&1 | tail -3 || true').strip()
        move_error+=f' | quota: {quota}'
    where=(f'~/{ARCHIVE_HOME}/{NAME}_fields.tar.gz' if moved
           else f'{REMOTE}/{NAME}_fields.tar.gz')
    rec=json.loads(ssh(f"python3 - <<'PY'\nimport hashlib,json,os\nfrom pathlib import Path\n"
        f"a=Path(os.path.expanduser('{where}'))\n"
        f"m=json.loads(Path('{REMOTE}/runs/{NAME}/post/field_hashes.json').read_text())\n"
        f"print(json.dumps(dict(archive=str(a),bytes=a.stat().st_size,"
        f"sha256=hashlib.sha256(a.read_bytes()).hexdigest(),frames=len(m))))\nPY",timeout=3600))
    rec.update(archive_matches_finished=rec['sha256']==fin['archive']['sha256'],
               field_manifest_sha256=c.sha(run/'post/field_hashes.json'),
               independent_checks=json.loads((run/'post/independent_checks.json').read_text()),
               moved_to_durable_storage=moved,move_error=move_error,
               volatility_warning=(None if moved else
                   'sal /tmp is a RAM-backed tmpfs, so this archive does not survive a reboot; '
                   'the per-frame hash manifest and every derived product are held locally'),
               retained_remotely_because='the workstation has under 20 GiB free; only derived products are pulled',
               local_files={str(n.relative_to(run)).replace(chr(92),'/'):dict(sha256=c.sha(n),bytes=n.stat().st_size)
                            for n in sorted(run.rglob('*')) if n.is_file()})
    c.save(OUT/'archive_receipt.json',rec)
    print(json.dumps(dict(frames=rec['frames'],bytes=rec['bytes'],
                          archive_matches_finished=rec['archive_matches_finished']),indent=2))


def reduce():
    """Run the two reviewed field passes remotely and pull their small products.

    Both passes assume one archived frame per recorded accepted step. This leg writes a
    frame every second step, so they run against a symlink view built by
    wp9_pair_frames.py, which carries every frame hash over unchanged.
    """
    scp(c.ROOT/'scripts/wp9_pair_frames.py',f'sal:{REMOTE}/wp9_pair_frames.py')
    pairing=json.loads(ssh(f'cd {REMOTE} && python3 wp9_pair_frames.py runs/{NAME} runs/{VIEW} '
                           f'{json.loads((OUT/"r96_protocol.json").read_text(encoding="utf-8"))["numerics"]["output_interval"]}',
                           timeout=1800).strip().splitlines()[-1])
    assert pairing['sha256_values_unchanged'] and pairing['frames']==pairing['kept_steps']
    c.save(OUT/'frame_pairing_receipt.json',pairing)
    outs={}
    for script,sub in [('wp8_fields.py','reduced'),('wp7_spatial.py','spatial')]:
        log=ssh(f'cd {REMOTE} && env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 '
                f'python3 {script} runs/{VIEW} {sub}_{NAME} 2>&1 | tail -3',timeout=5400)
        outs[script]=log.strip().splitlines()[-1]
    dest=OUT/'fields';dest.mkdir(parents=True,exist_ok=True)
    for remote,local in [(f'reduced_{NAME}/timeseries.csv','reduced_timeseries.csv'),
                         (f'reduced_{NAME}/audit.json','reduced_audit.json'),
                         (f'spatial_{NAME}/spatial_metrics.csv','spatial_metrics.csv'),
                         (f'spatial_{NAME}/spatial_audit.json','spatial_spatial_audit.json'),
                         (f'spatial_{NAME}/selected_vertical_profiles.csv','spatial_selected_vertical_profiles.csv')]:
        scp(f'sal:{REMOTE}/{remote}',str(dest/local))
    c.CACHE.mkdir(parents=True,exist_ok=True)
    scp(f'sal:{REMOTE}/reduced_{NAME}/profiles_spectra.npz',str(c.CACHE/'R96_M1_profiles_spectra.npz'))
    scp(f'sal:{REMOTE}/spatial_{NAME}/selected_fields.npz',str(dest/'spatial_selected_fields.npz'))
    audit=json.loads((dest/'reduced_audit.json').read_text())
    assert audit['all_field_hashes_verified']
    assert audit['files']['profiles_spectra.npz']==c.sha(c.CACHE/'R96_M1_profiles_spectra.npz')
    c.save(OUT/'reduction_receipt.json',dict(remote_logs=outs,reduced_audit=audit,
        spatial_audit=json.loads((dest/'spatial_spatial_audit.json').read_text()),
        local={p.name:dict(sha256=c.sha(p),bytes=p.stat().st_size) for p in sorted(dest.iterdir()) if p.is_file()},
        spectra_npz=dict(path=str(c.CACHE/'R96_M1_profiles_spectra.npz'),
                         sha256=c.sha(c.CACHE/'R96_M1_profiles_spectra.npz'))))
    print(json.dumps(dict(frames=audit['frames'],verified=audit['all_field_hashes_verified']),indent=2))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('action',
        choices=['stage','submit','status','collect','reduce'])
    a=ap.parse_args();globals()[a.action]()
