"""Collect read-only provenance/initial-field audit, not a simulation."""
import base64
import hashlib
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'analysis/matched_pair_36877_20261102'
OUT.mkdir(exist_ok=True)
script = ROOT/'scripts/audit_matched_pair_remote.py'
bridge = '/mnt/d/3OutputFile/1Research/Article/22Surrogate_model/paper/2/scripts/frodo_bridge.py'
p = subprocess.run(['wsl', 'python3', bridge, 'run',
                    'env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -'],
                   input=script.read_bytes(), capture_output=True, timeout=90)
if p.returncode:
    raise RuntimeError(p.stderr.decode(errors='replace'))
j = json.loads(p.stdout)
m0_input = base64.b64decode(j.pop('M0_input_base64'))
(OUT/'M0_actual_params.input').write_bytes(m0_input)
small = j.pop('legacy_files_base64')
folder = OUT/'legacy_nx28'
folder.mkdir(exist_ok=True)
for name, value in small.items():
    assert Path(name).name == name
    data = base64.b64decode(value)
    expected = next(v['sha256'] for k, v in j['files'].items() if k.endswith('/'+name))
    assert hashlib.sha256(data).hexdigest() == expected
    (folder/name).write_bytes(data)
j['audit_script_sha256'] = hashlib.sha256(script.read_bytes()).hexdigest()
(OUT/'remote_audit.json').write_text(json.dumps(j, indent=2), encoding='utf-8')
print(json.dumps({k:j[k] for k in ('executable','initial_checks','initial_field','profile')}, indent=2))
