"""Lightweight Git deliverables, full evidence retained separately on HPCC."""
import csv,gzip,hashlib,json,os,pathlib,shutil,subprocess,tarfile
import pandas as pd
assert os.environ.get('SLURM_JOB_ID')
root=pathlib.Path.cwd();tag=os.environ.get('P2H_DELIVERY_TAG','');suffix='_'+tag if tag else ''
out=root/('delivery'+suffix);out.mkdir(exist_ok=False)
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()
baseline=root.parent/'p2_h_crystal_boundary_20261108'
extra={}
actual=root/'ale_diagnostic/historical320_ale/coupled_history.csv'
extra['original_320_coupled_entire_history_equal']=actual.read_bytes()==(baseline/'matrix/crystal_n320/coupled_history.csv').read_bytes()
current=(root/'correction/repo/model/p2_h_ale_coupling/bdf.hh').read_text()
current=current.replace('    std::function<void(const Vec&,double)> rejectedState;\n','')
current=current.replace('        try { model.validate(y); } catch(...) { if(rejectedState) rejectedState(y,reached); throw; }\n        return y;','        model.validate(y); return y;')
extra['BDF_identical_after_removing_observer']=current==(baseline/'repo/model/p2_h_ale_coupling/bdf.hh').read_text()
assert all(extra.values()),extra
(root/'audit/additional_regressions.json').write_text(json.dumps(extra,indent=2))
shutil.copytree(root/'audit',out,dirs_exist_ok=True)
shutil.copytree(root/'correction/audit',out/'corrected_audit')
stages=['smoke','diagnostic','ale_diagnostic','ablations','correction/fixed_smoke','correction/fixed_tests','correction/fixed_pre','correction/fixed_coupled']
manifest=[];raw=[]
for stage in stages:
    folder=root/stage
    for d in folder.iterdir():
        if not d.is_dir():
            if d.suffix in ('.json','.csv'):raw.append(d)
            continue
        if not (d/'exit_status.txt').exists():raise RuntimeError('unfinished case '+str(d))
        target=out/'cases'/stage/d.name;target.mkdir(parents=True)
        for p in d.iterdir():
            if p.suffix not in ('.csv','.json','.txt','.out','.input'):continue
            raw.append(p);manifest.append(dict(path=str(p.relative_to(root)),bytes=p.stat().st_size,sha256=sha(p)))
            if p.name in ('coupled_history.csv','history.csv'):
                df=pd.read_csv(p);keep=df.index%100==0;keep|=df.index==len(df)-1
                if 'active' in df:
                    keep|=df.time_s==df.activation_s;keep|=df.rejected_trials.diff().fillna(0)!=0
                if 'status' in df:keep|=df.status!='accepted'
                df[keep].to_csv(target/p.name,index=False)
            elif p.name=='grid_profiles.csv':
                df=pd.read_csv(p);times=sorted(df.time_s.unique());wanted=[0,1200,1900,2500,times[-1]]
                h=pd.read_csv(d/'coupled_history.csv')
                if (h.active==1).any():wanted.append(h[h.active==1].time_s.iloc[0])
                df[df.time_s.isin(wanted)].to_csv(target/p.name,index=False)
            elif p.name=='ale_diagnostics.csv':
                df=pd.read_csv(p);accepted=df[df.status=='accepted'];birth=accepted.time_s.min()
                keep=(df.status!='accepted')|(df.time_s-birth<=10)|(df.time_s%100<1e-8)|(df.time_s==accepted.time_s.max())
                df[keep].to_csv(target/p.name,index=False)
            elif p.name in ('grid_budgets.csv','fields.csv','accepted_steps.csv','transport_budget.csv','component_budget.csv','cut_fluxes.csv','surface_events.csv','soil_fields.csv','crust_fields.csv') or p.name.endswith('_sigmatop.csv'):
                # Entire files remain in the checksum-pinned full archive.
                continue
            elif p.name=='solver.out':
                lines=p.read_text(errors='replace').splitlines()
                (target/'solver.tail.txt').write_text('\n'.join(lines[-45:])+'\n')
            else:shutil.copyfile(p,target/p.name)
    for name in ('cases.json','outcomes.json'):
        if (folder/name).exists():
            dest=out/'cases'/stage/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(folder/name,dest)
pd.DataFrame(manifest).to_csv(out/'raw_output_sha256.csv',index=False)
sources=[]
for prefix in ('repo','correction/repo','ablation_build'):
    for p in sorted((root/prefix).rglob('*')):
        if p.is_file() and '__pycache__' not in p.parts and p.name!='driver':sources.append(dict(path=str(p.relative_to(root)),bytes=p.stat().st_size,sha256=sha(p)))
pd.DataFrame(sources).to_csv(out/'runtime_source_sha256.csv',index=False)
for prefix in ('build','correction/build','ablation_build'):
    for p in (root/prefix).glob('*.log'):
        dest=out/'logs'/(prefix.replace('/','_')+'_'+p.name+'.gz');dest.parent.mkdir(exist_ok=True)
        with gzip.open(dest,'wb') as f:f.write(p.read_bytes())
for folder in (root,root/'correction'):
    for pattern in ('*.out','*.sha256'):
        for p in folder.glob(pattern):
            dest=out/'logs'/((('correction_' if folder!=root else '')+p.name)+('.gz' if p.suffix=='.out' else ''))
            if p.suffix=='.out':
                with gzip.open(dest,'wb') as f:f.write(p.read_bytes())
            else:shutil.copyfile(p,dest)
jobs='16380077,16380079,16380102,16380103,16380104,16380110,16380113,16380114,16380116,16380117,16380118,16380119,16380120,16380122,16380200,16380204,16380219,16380222,16380225,16380240,16380253,'+os.environ['SLURM_JOB_ID']
r=subprocess.run(['sacct','-j',jobs,'--format=JobID,JobName,State,ExitCode,Elapsed,AllocCPUS,NodeList','-P'],capture_output=True,text=True,check=True)
(out/'accounting.txt').write_text(r.stdout)
protected=['model/p2_h_ale_coupling/ale.hh','model/p2_h_ale_coupling/donor.hh','model/p2_d_coupling/cpp/pratcrust.hh']
protected += [str(p.relative_to(root/'correction/repo')) for p in (root/'correction/repo/model/p1_r2_conservative_richards/source/appl').rglob('*') if p.is_file()]
baseline=root.parent/'p2_h_crystal_boundary_20261108/repo'
for name in protected:assert sha(root/'correction/repo'/name)==sha(baseline/name),name
provenance=dict(base_commit='301097e3d2de3541bf63741cd729bd5cda22329a',branch='research/p2-h-grid-precipitation-fix',
    hpcc_root=str(root),diagnostic_binary_sha256=sha(root/'build/soil_solver'),corrected_binary_sha256=sha(root/'correction/build/soil_solver'),
    ablation_binary_sha256=sha(root/'ablation_build/driver'),protected_sources_unchanged=protected,
    source_archives={p.name:sha(p) for p in root.glob('*.tar.gz') if p.name.startswith(('source','overlay'))},
    stages=stages,raw_archive='full_grid_evidence.tar.gz',compact_archive='grid_delivery'+suffix+'.tar.gz',
    accepted_vs_rejected='Explicit status column; failed histories never extrapolated',
    compatibility_default='projection; corrected matrix explicitly selects implicit_surface',
    threads='DUMUX_NUM_THREADS=1; independent cases parallel in Slurm',scientific_status='See REPORT_ZH.md and gate CSVs')
provenance['compiled_main_sha256']={prefix:sha(root/prefix/'compiled_main.cc') for prefix in ('build','correction/build')}
provenance['historical_reference_sha256']={f'crystal_n{n}/{name}':sha(root.parent/f'p2_h_crystal_boundary_20261108/matrix/crystal_n{n}'/name)
    for n in (40,80) for name in ('coupled_history.csv','donor_attempts.csv')}
(out/'provenance.json').write_text(json.dumps(provenance,indent=2))
if suffix:
    old=pd.read_csv(root/'delivery/raw_output_sha256.csv')
    assert dict(zip(old.path,old.sha256))=={r['path']:r['sha256'] for r in manifest},'raw case files changed during visual QA'
    assert (root/'full_grid_evidence.tar.gz').exists()
else:
    with tarfile.open(root/'full_grid_evidence.tar.gz','w:gz') as tar:
        for p in sorted(set(raw)):tar.add(p,arcname=str(p.relative_to(root)))
compact='grid_delivery'+suffix+'.tar.gz'
with tarfile.open(root/compact,'w:gz') as tar:tar.add(out,arcname='delivery')
manifest_file=root/('delivery'+suffix+'.sha256')
manifest_file.write_text(''.join(sha(root/name)+'  '+name+'\n' for name in ('full_grid_evidence.tar.gz',compact)))
print(manifest_file.read_text())
