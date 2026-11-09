"""Verify transferred x009 evidence and the exact comparison source files."""
import hashlib,json,pathlib
repo=pathlib.Path(__file__).resolve().parents[2]
out=repo/'analysis/prat_analytical_boundary_20261109/x009_comparison'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
manifest=json.loads((out/'checksums.json').read_text())
for name,value in manifest.items(): assert sha(out/name)==value,name
sources=0
for name,value in json.loads((out/'provenance.json').read_text())['sha256'].items():
    if '/repo/' in name:
        local=repo/name.split('/repo/',1)[1]; assert sha(local)==value,name; sources+=1
print('PASS {} comparison artifacts and {} source files'.format(len(manifest),sources))
