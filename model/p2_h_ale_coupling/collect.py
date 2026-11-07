"""Create lightweight evidence plus separately archived full accepted histories."""
import csv
import gzip
import hashlib
import json
import pathlib
import shutil
import tarfile

root=pathlib.Path.cwd()
out=root/'deliverables'
out.mkdir(exist_ok=True)
manifest={}
for dirname in ('source','repo','reference','stage_a','stage_b','stage_c','failures','build'):
    for path in sorted((root/dirname).rglob('*')):
        if path.is_file() and '__pycache__' not in path.parts:
            manifest[str(path.relative_to(root))]=hashlib.sha256(path.read_bytes()).hexdigest()
for path in root.glob('*.slurm'):
    manifest[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
for path in root.glob('*.tar*'):
    if path.name not in ('full_accepted_histories.tar.gz','compact_deliverables.tar.gz'):
        manifest[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
(out/'remote_manifest.json').write_text(json.dumps(manifest,indent=2))

for stage in ('stage_a','stage_b','stage_c'):
    for name in ('verification.json','node.txt','job.txt','binary.sha256','binaries.sha256','outcomes.json'):
        path=root/stage/name
        if path.exists():
            dest=out/stage/name; dest.parent.mkdir(parents=True,exist_ok=True); shutil.copyfile(path,dest)
    for path in (root/stage).iterdir():
        if not path.is_dir(): continue
        dest=out/stage/path.name; dest.mkdir(parents=True,exist_ok=True)
        for name in ('params.input','command.json','exit_status.txt'):
            if (path/name).exists(): shutil.copyfile(path/name,dest/name)
        for name in ('history.csv','coupled_history.csv'):
            if not (path/name).exists(): continue
            with (path/name).open() as f: rows=list(csv.DictReader(f))
            chosen=[]
            for i,r in enumerate(rows):
                keep=(i%40==0 or i==len(rows)-1 or float(r.get('activation_s',-1))==float(r['time_s']))
                if i and 'active' in r:
                    keep=keep or r['active']!=rows[i-1]['active'] or r['rejected_trials']!=rows[i-1]['rejected_trials']
                if keep: chosen.append(r)
            with (dest/name).open('w',newline='') as f:
                writer=csv.DictWriter(f,fieldnames=rows[0]); writer.writeheader(); writer.writerows(chosen)
        for name in ('coupling_events.csv','solver.out'):
            if (path/name).exists():
                # Full logs retained remotely; light log tail is sufficient to see exit/rejection.
                lines=(path/name).read_text(errors='replace').splitlines()
                (dest/(name+'.tail.txt')).write_text('\n'.join(lines[-60:])+'\n')
for name in ('regression.log','operator.log','soil_build.log','sundials.log'):
    path=root/'build'/name
    if path.exists():
        dest=out/'compiler_logs'/(name+'.gz'); dest.parent.mkdir(exist_ok=True)
        with gzip.open(dest,'wb') as f: f.write(path.read_bytes())
for path in root.glob('stage-*.out'):
    dest=out/'job_logs'/path.name; dest.parent.mkdir(exist_ok=True); shutil.copyfile(path,dest)
for path in (root/'failures').rglob('verification.log'):
    dest=out/'failures'/str(path.relative_to(root/'failures'))
    dest.parent.mkdir(parents=True,exist_ok=True); shutil.copyfile(path,dest)
for name in ('accounting.txt','source_checksums.txt','postprocess_job.txt','postprocess_node.txt'):
    if (root/name).exists(): shutil.copyfile(root/name,out/name)

with tarfile.open(root/'full_accepted_histories.tar.gz','w:gz') as tar:
    for stage in ('stage_a','stage_b','stage_c'):
        for path in sorted((root/stage).rglob('*')):
            if path.is_file() and path.suffix in ('.csv','.json','.txt'):
                tar.add(path,arcname=str(path.relative_to(root)))
with tarfile.open(root/'compact_deliverables.tar.gz','w:gz') as tar:
    tar.add(out,arcname='deliverables')
(root/'transfer.sha256').write_text(''.join(hashlib.sha256((root/name).read_bytes()).hexdigest()+'  '+name+'\n'
    for name in ('full_accepted_histories.tar.gz','compact_deliverables.tar.gz')))
print((root/'transfer.sha256').read_text())
