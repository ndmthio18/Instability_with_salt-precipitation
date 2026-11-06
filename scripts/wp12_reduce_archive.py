"""Reduce trusted WP12-compatible archives in node-local scratch."""
from __future__ import annotations

import hashlib
import json
import os
import shutil
import socket
import subprocess
import sys
import tarfile
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def safe_extract(archive: Path, target: Path) -> None:
    target_resolved = target.resolve()
    with tarfile.open(archive, "r:gz") as handle:
        members = handle.getmembers()
        for member in members:
            resolved = (target / member.name).resolve()
            if target_resolved not in resolved.parents and resolved != target_resolved:
                raise RuntimeError(f"unsafe archive member {member.name}")
        handle.extractall(target, members=members)


def truncate_pvd(run: Path) -> dict:
    pvd_files = list(run.glob("*.pvd"))
    if len(pvd_files) != 1:
        raise RuntimeError(f"expected one PVD, got {len(pvd_files)}")
    available = {path.name for path in run.glob("*.vtu")}
    tree = ET.parse(pvd_files[0])
    collection = next(tree.getroot().iter("Collection"))
    before = len(list(collection))
    for dataset in list(collection):
        if Path(dataset.attrib["file"]).name not in available:
            collection.remove(dataset)
    after = len(list(collection))
    if after != len(available):
        raise RuntimeError(f"PVD/VTU mismatch after truncation: {after}/{len(available)}")
    tree.write(pvd_files[0], encoding="utf-8", xml_declaration=True)
    return {"datasets_before": before, "datasets_after": after, "vtu_files": len(available)}


def save_in_place(path: Path, value: dict) -> None:
    with path.open("w") as handle:
        json.dump(value, handle, indent=2)
        handle.flush()
        os.fsync(handle.fileno())


def run_checked(command: list[str], log: Path) -> None:
    result = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    with log.open("a") as handle:
        handle.write("$ " + " ".join(command) + "\n")
        handle.write(result.stdout)
    if result.returncode:
        raise RuntimeError(f"command failed ({result.returncode}): {' '.join(command)}")


def main(script_root: Path, fields: Path, products: Path, output: Path, manifest: Path, label: str, amplitude: str) -> None:
    job = os.environ["SLURM_JOB_ID"]
    scratch = Path(os.environ.get("TMPDIR", "/tmp")) / f"wp12_reduce_{label}_{job}"
    run, metrics = scratch / "run", scratch / "metrics"
    run.mkdir(parents=True)
    metrics.mkdir()
    local_fields, local_products = scratch / "fields.tar.gz", scratch / "products.tar.gz"
    shutil.copyfile(fields, local_fields)
    shutil.copyfile(products, local_products)
    source_hashes = {"fields": sha256(local_fields), "products": sha256(local_products)}
    safe_extract(local_products, run)
    safe_extract(local_fields, run)
    mapping = truncate_pvd(run)
    post = run / "post_reduced_0_6T"
    log = metrics / "reduction.log"
    run_checked([sys.executable, str(script_root / "wp9_post.py"), str(run), str(post)], log)
    run_checked([sys.executable, str(script_root / "wp12_scalar_metrics.py"), str(run), str(metrics), label,
                 amplitude, "--max-T", "6", "--post-dir", str(post)], log)
    run_checked([sys.executable, str(script_root / "wp12_field_metrics.py"), str(run), str(metrics), label,
                 "--max-T", "6", "--post-dir", str(post)], log)
    evidence = fields.parent
    acceptance = None
    if (evidence / "production_started.json").exists() and (evidence / "plan.json").exists():
        plan = json.loads((evidence / "plan.json").read_text())
        expected_job = (evidence / "submission.txt").read_text().strip()
        acceptance_dir = metrics / "acceptance"
        run_checked([sys.executable, str(script_root / "wp11_acceptance.py"), str(evidence), str(run),
                     str(acceptance_dir), label, plan["input_sha256"], expected_job], log)
        acceptance = json.loads((acceptance_dir / "acceptance.json").read_text())["passed"]
    with tarfile.open(scratch / "metrics.tar.gz", "w:gz", compresslevel=6) as handle:
        for path in sorted(metrics.rglob("*")):
            if path.is_file():
                handle.add(path, arcname=str(path.relative_to(metrics)))
    shutil.copyfile(scratch / "metrics.tar.gz", output)
    result = {
        "schema_version": 1,
        "utc": datetime.now(timezone.utc).isoformat(),
        "job_id": job,
        "node": socket.gethostname(),
        "label": label,
        "amplitude": float(amplitude),
        "source_archives": {"fields": str(fields), "products": str(products)},
        "source_hashes": source_hashes,
        "pvd_mapping": mapping,
        "numerical_acceptance_passed": acceptance,
        "output_bytes": output.stat().st_size,
        "output_sha256": sha256(output),
    }
    save_in_place(manifest, result)
    shutil.rmtree(scratch)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4]), Path(sys.argv[5]), sys.argv[6], sys.argv[7])
