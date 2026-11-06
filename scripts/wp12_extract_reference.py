"""Extract the dense 0--6T R96 field prefix from the immutable WP11 archive."""
from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import socket
import sys
import tarfile
from datetime import datetime, timezone
from pathlib import Path


EXPECTED_SOURCE = "a8bafcd3f1963813e9fc2a373914790218bd401084dc7dab08d8066c1f3716cb"
FRAME = re.compile(r"-(\d{5})\.vtu$")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def save_in_place(path: Path, value: dict) -> None:
    with path.open("w") as handle:
        json.dump(value, handle, indent=2)
        handle.flush()
        os.fsync(handle.fileno())


def main(source: Path, destination: Path, manifest: Path) -> None:
    scratch = Path(os.environ.get("TMPDIR", "/tmp")) / f"wp12_reference_{os.environ['SLURM_JOB_ID']}"
    scratch.mkdir()
    local_source = scratch / "fields.tar.gz"
    local_output = scratch / "R96_reference_0_6T_fields.tar.gz"
    shutil.copyfile(source, local_source)
    source_hash = sha256(local_source)
    if source_hash != EXPECTED_SOURCE:
        raise RuntimeError(f"source hash drift: {source_hash}")

    names = []
    with tarfile.open(local_source, "r|gz") as incoming, tarfile.open(local_output, "w:gz", compresslevel=1) as outgoing:
        pvd_names = []
        for member in incoming:
            match = FRAME.search(member.name)
            is_pvd = member.name.endswith(".pvd")
            if not is_pvd and (not match or int(match.group(1)) > 144):
                continue
            stream = incoming.extractfile(member)
            if stream is None:
                raise RuntimeError(f"could not read {member.name}")
            outgoing.addfile(member, stream)
            names.append(member.name)
            if is_pvd:
                pvd_names.append(member.name)
    field_names = [name for name in names if FRAME.search(name)]
    selected_indices = sorted(int(FRAME.search(name).group(1)) for name in field_names)
    if len(field_names) != 145 or selected_indices != list(range(145)) or len(pvd_names) != 1:
        raise RuntimeError(f"unexpected field selection: {len(names)} {names[:1]} {names[-1:]}")

    shutil.copyfile(local_output, destination)
    result = {
        "schema_version": 1,
        "utc": datetime.now(timezone.utc).isoformat(),
        "job_id": os.environ["SLURM_JOB_ID"],
        "node": socket.gethostname(),
        "source": str(source),
        "source_sha256": source_hash,
        "selected_frames": len(names),
        "selected_vtu_frames": len(field_names),
        "pvd_member": pvd_names[0],
        "first_member": next(name for name in field_names if name.endswith("-00000.vtu")),
        "last_member": next(name for name in field_names if name.endswith("-00144.vtu")),
        "destination_bytes": destination.stat().st_size,
        "destination_sha256": sha256(destination),
    }
    save_in_place(manifest, result)
    shutil.rmtree(scratch)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main(*(Path(value) for value in sys.argv[1:4]))
