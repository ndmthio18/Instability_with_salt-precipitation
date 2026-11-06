"""Reduce the WP14 dt/2 archives and run the WP14-specific acceptance gate."""
from __future__ import annotations

import json
import os
import shutil
import socket
import sys
import tarfile
from datetime import datetime, timezone
from pathlib import Path

from wp12_reduce_archive import run_checked, safe_extract, save_in_place, sha256, truncate_pvd


def main(script_root: Path, evidence: Path, output: Path, manifest: Path) -> None:
    job = os.environ["SLURM_JOB_ID"]
    scratch = Path(os.environ.get("TMPDIR", "/tmp")) / f"wp14_reduce_R48_{job}"
    run, metrics = scratch / "run", scratch / "metrics"
    run.mkdir(parents=True)
    metrics.mkdir()
    local_fields, local_products = scratch / "fields.tar.gz", scratch / "products.tar.gz"
    shutil.copyfile(evidence / "fields.tar.gz", local_fields)
    shutil.copyfile(evidence / "products.tar.gz", local_products)
    source_hashes = {"fields": sha256(local_fields), "products": sha256(local_products)}
    safe_extract(local_products, run)
    safe_extract(local_fields, run)
    mapping = truncate_pvd(run)
    post = run / "post_reduced_0_6T"
    log = metrics / "reduction.log"
    run_checked([sys.executable, str(script_root / "wp9_post.py"), str(run), str(post)], log)
    run_checked([sys.executable, str(script_root / "wp12_scalar_metrics.py"), str(run), str(metrics),
                 "R48", "0.0001", "--max-T", "6", "--post-dir", str(post)], log)
    run_checked([sys.executable, str(script_root / "wp12_field_metrics.py"), str(run), str(metrics),
                 "R48", "--max-T", "6", "--post-dir", str(post)], log)
    shutil.copyfile(run / "transport_budget.csv", metrics / "transport_budget.csv")
    shutil.copyfile(run / "surface_events.csv", metrics / "surface_events.csv")
    shutil.copyfile(post / "profiles.npz", metrics / "profiles.npz")
    shutil.copyfile(post / "representative_fields.npz", metrics / "representative_fields.npz")
    acceptance_dir = metrics / "acceptance"
    acceptance_dir.mkdir()
    run_checked([sys.executable, str(script_root / "wp14_archive_acceptance.py"), str(evidence),
                 str(acceptance_dir / "acceptance.json")], log)
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
    main(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4]))
