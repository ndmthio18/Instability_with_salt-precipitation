"""Read-only copy of the registered frodo ELF, its runtime and reference inputs.

Run on frodo. Writes only a new export directory under the Pilot-0 run root.
"""
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import tarfile

ROOT=Path('/temp/zhixin/paper2_transport_20261102/pilot0_integrated_20261102_r01')
OUT=ROOT/'sal_export'; OUT.mkdir(exist_ok=False)
bundle=OUT/'bundle';bundle.mkdir()
plan=json.loads((ROOT/'M1/plan.json').read_text())
exe=Path(plan['executable']['path'])
assert hashlib.sha256(exe.read_bytes()).hexdigest()==plan['executable']['sha256']
manifest={};seen=set()


def copy(source,target):
    target=bundle/target;target.parent.mkdir(parents=True,exist_ok=True)
    h=hashlib.sha256(source.read_bytes()).hexdigest()
    if target.exists():assert hashlib.sha256(target.read_bytes()).hexdigest()==h
    else:shutil.copy2(source,target)
    manifest[str(target.relative_to(bundle))]=dict(source=str(source),sha256=h,bytes=source.stat().st_size)


def dependencies(path):
    if str(path) in seen:return
    seen.add(str(path))
    output=subprocess.check_output(['ldd',str(path)],text=True)
    if 'not found' in output:raise ValueError(output)
    for line in output.splitlines():
        found=re.search(r'(?:=>\s+)?(/[^\s]+)',line)
        if found:
            p=Path(found.group(1));copy(p,'lib/'+p.name)


copy(exe,'bin/paper2_transport');dependencies(exe)
orted=Path('/usr/bin/orted')
if orted.exists():copy(orted,'bin/orted.real');dependencies(orted)
plugins=Path('/usr/lib/x86_64-linux-gnu/openmpi/lib/openmpi3')
for p in sorted(plugins.glob('*.so')):
    copy(p,'openmpi/'+p.name);dependencies(p)
for relative in ('M0','M1'):
    for p in (ROOT/relative).iterdir():
        if p.is_file():copy(p,relative+'/'+p.name)
m0=Path(plan['M0_remote'])
copy(next(m0.glob('*-00000.vtu')),'M0/'+next(m0.glob('*-00000.vtu')).name)
for name,h in plan['source_sha256'].items():
    src=Path(plan['remote_source'])/name
    assert hashlib.sha256(src.read_bytes()).hexdigest()==h
    copy(src,'source/'+name)
copy(ROOT/'pilot0_field_review.py','pilot0_field_review.py')
(bundle/'runtime_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
with tarfile.open(OUT/'bundle.tar.gz','w:gz',compresslevel=1) as tar:
    tar.add(bundle,arcname='bundle')
record=dict(archive=str(OUT/'bundle.tar.gz'),sha256=hashlib.sha256((OUT/'bundle.tar.gz').read_bytes()).hexdigest(),
    bytes=(OUT/'bundle.tar.gz').stat().st_size,files=len(manifest),elf_sha256=plan['executable']['sha256'])
(OUT/'export.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record))
