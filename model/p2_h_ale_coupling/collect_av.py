"""Build the lightweight Git delivery while retaining raw evidence on HPCC."""
import hashlib,json,pathlib,shutil

root=pathlib.Path.cwd()
target=root/'delivery'/'analysis'/'p2_h_av_sensitivity_20261108'
if target.exists():
    raise RuntimeError(f'refusing to overwrite {target}')
shutil.copytree(root/'analysis',target)
for folder in ('unit_exact','stage_a_exact','mechanism_exact'):
    shutil.copytree(root/folder,target/folder)
logs=target/'logs'; logs.mkdir()
for pattern in ('build-*.out','smoke-*.out','matrix-*.out','build-exact-*.out','exact-resume-*.out','audit-*.out'):
    for path in root.glob(pattern): shutil.copy2(path,logs/path.name)
case_root=target/'cases'
for stage in ('smoke_exact','primary_exact','dt_half_exact','ale_refine_exact','ale256_exact'):
    for case in json.loads((root/stage/'cases.json').read_text()):
        source=root/stage/case['name']; destination=case_root/stage/case['name']; destination.mkdir(parents=True)
        for name in ('command.json','params.input','exit_status.txt','activation_config.csv','coupling_events.csv',
                     'donor_attempts.csv','donor_profiles.csv','stop_record.json','transport_metadata.json'):
            if (source/name).exists(): shutil.copy2(source/name,destination/name)
        if (source/'solver.out').exists():
            lines=(source/'solver.out').read_text(errors='replace').splitlines()
            (destination/'solver.tail.txt').write_text('\n'.join(lines[-120:])+'\n')
source_target=target/'source'; source_target.mkdir()
for name in ('ale.hh','bdf.hh','soil_methods.inc','grid_diagnostics.inc','donor.hh','prepare_dumux.py',
             'av_parameter_test.cc','av_mechanism.cc','run_av.py','analyze_av.py','verify_av.py','plot_av.py',
             'AV_SENSITIVITY.md','build_av_exact.slurm','resume_av_exact.slurm','matrix_av.slurm','audit_av.slurm'):
    shutil.copy2(root/'repo/model/p2_h_ale_coupling'/name,source_target/name)
archives={}
for path in sorted(root.glob('*r*.sha256')):
    archives[path.name]=path.read_text().strip()
provenance=dict(base_commit='301097e3d2de3541bf63741cd729bd5cda22329a',
    branch='research/p2-h-av-sensitivity',hpcc_root=str(root),source_archive_hashes=archives,
    exact_binary_hashes=(root/'analysis/executable_sha256.txt').read_text().splitlines(),
    historical_default='ReactionCoefficient=2.3e-3, AvMultiplier=1, splitReactionLedgers=false',
    sensitivity_configuration='ReactionCoefficient=1e-3, AvMultiplier in {0.1,1,10}, splitReactionLedgers=true',
    preliminary_endpoint_quadrature='retained in raw primary/dt_half/ale_refine; excluded from final quantitative tables',
    raw_evidence='full_av_evidence.tar.gz retained on HPCC; see archive SHA-256')
(target/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
manifest=[]
for path in sorted(target.rglob('*')):
    if path.is_file(): manifest.append(f'{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(target).as_posix()}')
(target/'delivery.sha256').write_text('\n'.join(manifest)+'\n')
print(target)
