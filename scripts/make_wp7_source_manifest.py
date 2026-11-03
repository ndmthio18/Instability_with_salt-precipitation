"""Mechanically freeze hashes of the isolated WP7 source tree."""
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'model/transport_wp7_localflux_v1'
SOURCE=BASE/'source'
records={str(path.relative_to(SOURCE)).replace('\\','/'):{
    'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
    'bytes':path.stat().st_size,
} for path in sorted(SOURCE.rglob('*')) if path.is_file() and 'build-cmake' not in path.parts and '.git' not in path.parts}
(BASE/'source_manifest.json').write_text(json.dumps(records,indent=2)+'\n')
print(len(records),BASE/'source_manifest.json')
