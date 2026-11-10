"""Restore original conversation sources from a verified, lossless payload.

Uses only Python's standard library. Writes only this archive directory.
The script does not access credentials, run network requests, or modify model code.
"""
from __future__ import annotations
import base64
import hashlib
import io
import json
from pathlib import Path
import tarfile
import zipfile

ROOT = Path(__file__).resolve().parent
EXPECTED_ARCHIVE = '6bf0643f5f1ec32122748eee1c711d60eb5cfe8f84504fc52ba73f13731bb921'
EXPECTED = {
    'Prat_crust_salt_balance_editable.svg': (48244, 'afee2cefa9309442d23932e187c5c4094c9c7809f073b537158901ccd5107ba5'),
    'prat_salt_balance.svg': (48244, 'afee2cefa9309442d23932e187c5c4094c9c7809f073b537158901ccd5107ba5'),
    'prat_editable_direct.svg': (19877, 'ebe97e62b0ad0a43eff0a7bf56baa6f671c49ebe349f6d5d0234871468c2f353'),
    'prat_editable_direct.txt': (19877, 'ebe97e62b0ad0a43eff0a7bf56baa6f671c49ebe349f6d5d0234871468c2f353'),
    'prat_interface_editable.svg': (53026, '4871a9f7918d3aaddf5717e1927fbc3fbe3998c8837d329a05e09f7cc8aa7c75'),
    'prat_interface_source.txt': (53026, '4871a9f7918d3aaddf5717e1927fbc3fbe3998c8837d329a05e09f7cc8aa7c75'),
}

def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()

def checked_write(path: Path, data: bytes) -> None:
    if path.exists() and path.read_bytes() != data:
        raise RuntimeError(f'Refusing to overwrite changed file: {path}')
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)

def main() -> None:
    parts = [ROOT / '_payload' / f'{i:02d}.b64' for i in range(17)]
    payload = base64.b64decode(''.join(p.read_text(encoding='ascii').strip() for p in parts), validate=True)
    if len(payload) != 12120 or sha256(payload) != EXPECTED_ARCHIVE:
        raise RuntimeError('Compressed payload checksum mismatch')
    sources: dict[str, bytes] = {}
    with tarfile.open(fileobj=io.BytesIO(payload), mode='r:xz') as archive:
        for member in archive.getmembers():
            if not member.isfile() or member.name not in EXPECTED or member.name in sources:
                raise RuntimeError(f'Unexpected archive member: {member.name}')
            size, digest = EXPECTED[member.name]
            if member.size != size:
                raise RuntimeError(f'Wrong size: {member.name}')
            stream = archive.extractfile(member)
            if stream is None:
                raise RuntimeError(f'Unreadable member: {member.name}')
            data = stream.read(size + 1)
            if len(data) != size or sha256(data) != digest:
                raise RuntimeError(f'File checksum mismatch: {member.name}')
            sources[member.name] = data
    if set(sources) != set(EXPECTED):
        raise RuntimeError('Missing source files')
    report = []
    for name, data in sources.items():
        path = ROOT / 'editable' / name
        checked_write(path, data)
        report.append({'path': str(path.relative_to(ROOT)), 'bytes': len(data), 'sha256': sha256(data), 'status': 'VERIFIED_ORIGINAL'})
    # Reproduce the small original ZIP byte-for-byte, including its ZIP metadata.
    buff = io.BytesIO()
    info = zipfile.ZipInfo('prat_salt_balance.svg', (2026, 11, 10, 9, 13, 34))
    info.create_system = 3
    info.create_version = 20
    info.extract_version = 20
    info.external_attr = 2175008768
    with zipfile.ZipFile(buff, 'w') as archive:
        archive.writestr(info, sources['prat_salt_balance.svg'], compress_type=zipfile.ZIP_DEFLATED, compresslevel=6)
    small_zip = buff.getvalue()
    if sha256(small_zip) != 'bafaea275ba7046187558d8d7317a3965aa0d19d7f0a1fe157f046f31661b396':
        raise RuntimeError('Original ZIP reproduction mismatch')
    checked_write(ROOT / 'original_archives/prat_salt_balance.zip', small_zip)
    report.append({'path': 'original_archives/prat_salt_balance.zip', 'bytes': len(small_zip), 'sha256': sha256(small_zip), 'status': 'VERIFIED_ORIGINAL'})
    picture = ROOT / 'images/figure9_physical_columns.png'
    if picture.exists():
        data = picture.read_bytes()
        if sha256(data) != '9d7a58764f56bf545de6e2f5612a968268f8422e211079c850d25847f69b2c3b':
            raise RuntimeError('Original figure9 checksum mismatch')
        report.append({'path': 'images/figure9_physical_columns.png', 'bytes': len(data), 'sha256': sha256(data), 'status': 'VERIFIED_ORIGINAL'})
    # New download bundle; not the old conversation package containing PNG previews.
    buff = io.BytesIO()
    with zipfile.ZipFile(buff, 'w') as archive:
        for row in report:
            info = zipfile.ZipInfo(row['path'], (2026, 11, 10, 0, 0, 0))
            archive.writestr(info, (ROOT / row['path']).read_bytes(), compress_type=zipfile.ZIP_DEFLATED, compresslevel=6)
    checked_write(ROOT / 'verified_download_bundle.zip', buff.getvalue())
    record = {'status': 'PARTIAL_CONVERSATION_ARCHIVE', 'verified_original_count': len(report), 'payload_sha256': EXPECTED_ARCHIVE, 'files': report, 'note': 'Original AI-generated PNGs and the large original preview package are not included. See UPLOAD_STATUS.md.'}
    checked_write(ROOT / 'VERIFIED_FILES.json', (json.dumps(record, ensure_ascii=False, indent=2) + '\n').encode('utf-8'))
    print(json.dumps(record, ensure_ascii=False, indent=2))

if __name__ == '__main__':
    main()
