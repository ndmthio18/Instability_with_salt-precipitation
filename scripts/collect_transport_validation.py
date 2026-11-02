"""Retrieve and verify completed small validation evidence; never launch simulations."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tarfile

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'analysis/transport_validation_20261102'
REMOTE = '/temp/zhixin/paper2_transport_20261102/validation_20261102_r01'
BRIDGE = '/mnt/d/3OutputFile/1Research/Article/22Surrogate_model/paper/2/scripts/frodo_bridge.py'
ARCHIVE = ROOT/'data/generated/validation_20261102_r01.tar.gz'
RAW = ROOT/'data/generated/transport_validation_20261102/r01'


def remote(code):
    p = subprocess.run(['wsl', 'python3', BRIDGE, 'run', 'python3 -'], input=code.encode(), capture_output=True, timeout=60)
    if p.returncode:
        raise RuntimeError(p.stderr.decode(errors='replace'))
    return json.loads(p.stdout)


plan = json.loads((OUT/'plan.json').read_text())
upstream = json.loads((ROOT/'model/transport_v1/upstream_manifest.json').read_text())
provenance = remote('''
import hashlib,json,datetime,tarfile
from pathlib import Path
root=Path(REMOTE)
finished=json.loads((root/'finished.json').read_text())
assert finished['completed']==finished['planned']==8
with tarfile.open(str(root)+'.tar.gz','w:gz') as t:
 for name in ['plan.json','summary.json','finished.json','worker.py','worker.log','runs']:t.add(root/name,arcname=name)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
plan=PLAN
upstream=UPSTREAM
shared={n:sha(Path(upstream['remote_path'])/n) for n in upstream['files']}
isolated={n:sha(Path('/temp/zhixin/paper2_transport_20261102/source')/n) for n in plan['build_source_sha256']}
print(json.dumps({'checked_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
 'shared_source_unchanged':all(shared[n]==h['sha256'] for n,h in upstream['files'].items()),
 'isolated_source_matches_build':isolated==plan['build_source_sha256'],
 'executable_matches_build':sha(Path(plan['executable']['path']))==plan['executable']['sha256'],
 'shared_sha256':shared,'isolated_sha256':isolated,'finished':finished,
 'archive_sha256':sha(Path(str(root)+'.tar.gz'))}))
'''.replace('REMOTE', repr(REMOTE)).replace('PLAN', repr(plan)).replace('UPSTREAM', repr(upstream)))
for key in ('shared_source_unchanged', 'isolated_source_matches_build', 'executable_matches_build'):
    assert provenance[key], key
ARCHIVE.parent.mkdir(parents=True, exist_ok=True)
wsl_target = '/mnt/'+ARCHIVE.drive[0].lower()+ARCHIVE.as_posix()[2:]
subprocess.run(['wsl', 'python3', BRIDGE, 'pull', REMOTE+'.tar.gz', wsl_target], check=True, timeout=60)
assert hashlib.sha256(ARCHIVE.read_bytes()).hexdigest() == provenance['archive_sha256']
RAW.mkdir(parents=True, exist_ok=True)
with tarfile.open(ARCHIVE) as t:
    for member in t.getmembers():
        if member.isfile():
            existing = RAW/member.name
            if existing.exists():
                assert existing.read_bytes() == t.extractfile(member).read(), 'Refuse to overwrite changed evidence: '+str(existing)
    t.extractall(RAW, filter='data')
summary = json.loads((RAW/'summary.json').read_text())
count = 0
for s in summary:
    for name, record in s['files'].items():
        p = RAW/'runs'/s['id']/name
        assert hashlib.sha256(p.read_bytes()).hexdigest() == record['sha256'], str(p)
        count += 1
for p in RAW.rglob('*'):
    if p.is_file() and p.suffix not in ('.vtu', '.pvd'):
        dest = OUT/'evidence'/p.relative_to(RAW)
        dest.parent.mkdir(parents=True, exist_ok=True)
        if dest.exists():
            assert dest.read_bytes() == p.read_bytes(), str(dest)
        else:
            shutil.copyfile(p, dest)
(OUT/'source_after_validation.json').write_text(json.dumps(provenance, indent=2), encoding='utf-8')
(OUT/'transfer_check.json').write_text(json.dumps(dict(archive_sha256=provenance['archive_sha256'], verified_files=count, hash_mismatches=[]), indent=2), encoding='utf-8')
print(json.dumps(dict(verified_files=count, shared_source_unchanged=provenance['shared_source_unchanged'], completed=8)))
