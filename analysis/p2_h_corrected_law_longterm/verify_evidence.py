"""Verify retrieved archives and the production source/version crosswalk."""
from pathlib import Path
import hashlib
import json
R=Path(__file__).resolve().parent
record=json.loads((R/'execution_record.json').read_text())
def digest(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
checks={}
for name,v in record['versions'].items():
    path=R/v['driver_copy'];checks[name]=digest(path)==v['driver_sha256']
for name,expected in record['shared_operator_sha256'].items():
    checks[name]=digest(R.parent.parent/'model/p2_h_ale_coupling'/name)==expected
checks['forcing.csv']=digest(R.parent/'p2_h_corrected_vs_ale_x008_20261112/results/x008_forcing.csv')==record['forcing_sha256']
for line in (R/'final_sha256.txt').read_text().splitlines():
    expected,filename=line.split(maxsplit=1)
    if filename.startswith('cases') and filename.endswith('.csv'):
        path=R/filename if (R/filename).exists() else R/'provenance'/filename
        checks[filename]=digest(path)==expected
for manifest in ['compact.sha256','compact_p1.sha256','compact_refine.sha256']:
    expected,filename=(R/manifest).read_text().split()
    checks[filename]=digest(R/filename)==expected
out={'all_pass':all(checks.values()),'checks':checks,
     'full_archive':'Remote retained; full archive digest is recorded, not falsely claimed locally verified.'}
(R/'results/evidence_verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
assert out['all_pass'],checks
