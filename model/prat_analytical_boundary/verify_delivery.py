"""Metadata-only local verification; numerical validation belongs on HPCC."""
import hashlib,json,pathlib
repo=pathlib.Path(__file__).resolve().parents[2]
out=repo/'analysis/prat_analytical_boundary_20261109'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
count=0
for name,value in json.loads((out/'checksums.json').read_text()).items():
    assert sha(out/name)==value,('delivery checksum mismatch',name)
    count+=1
sources=0
for name,value in json.loads((out/'provenance.json').read_text())['sha256'].items():
    if name.startswith('repo/'):
        assert sha(repo/name[5:])==value,('tested source mismatch',name)
        sources+=1
print('PASS {} artifact hashes and {} tested source hashes'.format(count,sources))
# Includes the final human-authored report and post-job scheduler accounting.
manifest={str(p.relative_to(out)).replace('\\','/'):sha(p) for p in out.rglob('*') if p.is_file() and p.name!='DELIVERY_SHA256.json'}
(out/'DELIVERY_SHA256.json').write_text(json.dumps(manifest,indent=2)+'\n')
