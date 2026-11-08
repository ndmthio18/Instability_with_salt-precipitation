"""Keep compact auditable evidence in Git and full accepted histories separately."""
import csv,gzip,hashlib,json,pathlib,shutil,tarfile
root=pathlib.Path.cwd(); out=root/'deliverables'; out.mkdir(exist_ok=True)
shutil.copytree(root/'audit',out/'audit',dirs_exist_ok=True)
manifest={}
for dirname in ('repo','controlled','smoke','matrix','refusal','build','failures'):
    for p in sorted((root/dirname).rglob('*')):
        if p.is_file() and '__pycache__' not in p.parts:
            manifest[str(p.relative_to(root))]=hashlib.sha256(p.read_bytes()).hexdigest()
for p in root.glob('*.tar.gz'):manifest[p.name]=hashlib.sha256(p.read_bytes()).hexdigest()
(out/'remote_manifest.json').write_text(json.dumps(manifest,indent=2))
for stage in ('controlled','smoke','matrix','refusal'):
    for p in (root/stage).glob('*'):
        target=out/stage/p.name
        if p.is_file(): target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,target);continue
        target.mkdir(parents=True,exist_ok=True)
        for name in ('command.json','params.input','exit_status.txt','donor_profiles.csv'):
            if (p/name).exists():shutil.copyfile(p/name,target/name)
        for name in ('coupled_history.csv','donor_attempts.csv'):
            if not (p/name).exists():continue
            with (p/name).open() as f: rows=list(csv.DictReader(f))
            selected=[]
            for i,r in enumerate(rows):
                keep=i%100==0 or i==len(rows)-1
                if 'active' in r:
                    keep=keep or r['time_s']==r['activation_s'] or (i>0 and (r['active']!=rows[i-1]['active'] or r['rejected_trials']!=rows[i-1]['rejected_trials']))
                else:keep=keep or float(r['sufficient'])==1
                if keep:selected.append(r)
            if rows:
                with (target/name).open('w',newline='') as f:
                    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(selected)
        for name in ('coupling_events.csv','solver.out'):
            if (p/name).exists():
                lines=(p/name).read_text(errors='replace').splitlines()
                (target/(name+'.tail.txt')).write_text('\n'.join(lines[-35:])+'\n')
for p in (root/'build').glob('*.log'):
    dest=out/'compiler_logs'/(p.name+'.gz');dest.parent.mkdir(exist_ok=True)
    with gzip.open(dest,'wb') as f:f.write(p.read_bytes())
for pattern in ('*.out','*.txt','*.sha256'):
    for p in root.glob(pattern):shutil.copyfile(p,out/p.name)
if (root/'resume_metadata.json').exists():shutil.copyfile(root/'resume_metadata.json',out/'resume_metadata.json')
with tarfile.open(root/'full_accepted_histories.tar.gz','w:gz') as tar:
    for stage in ('controlled','smoke','matrix','refusal','failures'):
        for p in sorted((root/stage).rglob('*')):
            if p.is_file() and p.suffix in ('.csv','.json','.txt'):tar.add(p,arcname=str(p.relative_to(root)))
with tarfile.open(root/'compact_deliverables.tar.gz','w:gz') as tar:tar.add(out,arcname='deliverables')
(root/'transfer.sha256').write_text(''.join(hashlib.sha256((root/name).read_bytes()).hexdigest()+'  '+name+'\n' for name in ('full_accepted_histories.tar.gz','compact_deliverables.tar.gz')))
print((root/'transfer.sha256').read_text())
