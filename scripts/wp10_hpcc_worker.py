"""Run the single authorized WP10 trajectory in node-local scratch and stage archives."""
import csv, hashlib, json, os, resource, shutil, signal, socket, subprocess, sys, tarfile, time, traceback
from datetime import datetime, timezone
from pathlib import Path

ELF = "112076759fdff0942000307f729c2bd3c7e96963229b43c1b67c1f12bbf9df81"
INPUT = "c1441be8b69a11a23e9c8ce36839c127866d1d5567207fb0c63b6645c33d1bd1"
PROFILE = "d721e605d126e555f92b43200aa73b999ed4e2b61397edd1e98c0d9fddb36575"
RUN = "WP10-R24-M1-long50T"
SOLVER_TIMEOUT_S = 4 * 3600

def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024*1024), b""): h.update(block)
    return h.hexdigest()

def save(path, value):
    tmp = Path(str(path)+".tmp")
    tmp.write_text(json.dumps(value, indent=2))
    tmp.replace(path)

def progress(run, durable, process, started):
    path = run/"surface_events.csv"; rows = 0; tail = None
    if path.exists():
        with path.open(newline="") as f:
            data = list(csv.DictReader(f))
        rows = len(data); tail = data[-1] if data else None
    save(durable/"progress.json", {
        "utc": datetime.now(timezone.utc).isoformat(), "pid": process.pid,
        "elapsed_s": time.monotonic()-started, "accepted_rows_including_initial": rows,
        "last_surface_event": tail,
    })

def archive(run, path, fields):
    with tarfile.open(path, "w:gz", compresslevel=1) as tar:
        for f in run.rglob("*"):
            if not f.is_file(): continue
            isfield = f.suffix in (".vtu", ".pvd")
            if isfield == fields: tar.add(f, arcname=str(f.relative_to(run)))

def main(root):
    durable = Path(root); threads = int(os.environ["WP10_THREADS"])
    assert os.environ.get("SLURM_JOB_ID") and threads == 4
    assert not (durable/"production_started.json").exists(), "single-run guard"
    scratch_root = Path(os.environ.get("TMPDIR", "/tmp"))
    usage = shutil.disk_usage(scratch_root)
    mem = int(next(x.split()[1] for x in Path("/proc/meminfo").read_text().splitlines()
                   if x.startswith("MemAvailable:"))) * 1024
    assert usage.free >= 20*1024**3 and mem >= 16*1024**3
    run = scratch_root/f"wp10_{os.environ['SLURM_JOB_ID']}"/RUN
    run.mkdir(parents=True)
    shutil.copytree(durable/"bundle", run/"bundle", symlinks=False)
    for name in ("params.input", "eigenprofile.dat"):
        shutil.copy2(durable/name, run/name)
    assert sha(run/"bundle/bin/paper2_transport") == ELF
    assert sha(run/"params.input") == INPUT and sha(run/"eigenprofile.dat") == PROFILE
    started_record = {
        "utc": datetime.now(timezone.utc).isoformat(), "job_id": os.environ["SLURM_JOB_ID"],
        "node": socket.gethostname(), "affinity": sorted(os.sched_getaffinity(0)),
        "threads": threads, "scratch": str(run), "scratch_free_bytes": usage.free,
        "available_memory_bytes": mem, "elf_sha256": ELF, "input_sha256": INPUT,
        "profile_sha256": PROFILE, "solver_timeout_s": SOLVER_TIMEOUT_S,
    }
    save(durable/"production_started.json", started_record)
    env = dict(os.environ, DUMUX_NUM_THREADS=str(threads), OPENBLAS_NUM_THREADS="1",
               OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    def limit(): resource.setrlimit(resource.RLIMIT_AS, (8*1024**3, 8*1024**3))
    started = time.monotonic()
    with (run/"solver.log").open("w") as log:
        proc = subprocess.Popen([sys.executable, str(run/"bundle/launch_solver.py"), "params.input"],
            cwd=run, stdout=log, stderr=subprocess.STDOUT, env=env, preexec_fn=limit)
        while proc.poll() is None:
            progress(run, durable, proc, started)
            if time.monotonic()-started > SOLVER_TIMEOUT_S:
                proc.send_signal(signal.SIGTERM); time.sleep(20)
                if proc.poll() is None: proc.kill()
                break
            time.sleep(60)
        rc = proc.wait()
    progress(run, durable, proc, started)
    stop = json.loads((run/"stop_record.json").read_text()) if (run/"stop_record.json").exists() else None
    solver = {"returncode": rc, "wall_s": time.monotonic()-started, "stop": stop,
              "peak_rss_kib": resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss}
    save(run/"solver_finished.json", solver)
    if rc == 0 and stop and stop["reason"] in ("requested_horizon", "surface_solubility_threshold",
                                                "interior_solubility_threshold"):
        q = subprocess.run([sys.executable, str(durable/"wp9_post.py"), str(run), str(run/"post")],
                           capture_output=True, text=True, timeout=7200)
        (run/"post.log").write_text(q.stdout+q.stderr)
        if q.returncode: solver["post_returncode"] = q.returncode
    fields = durable/"fields.tar.gz"; products = durable/"products.tar.gz"
    archive(run, fields, True); archive(run, products, False)
    result = {"solver": solver,
              "fields": {"bytes": fields.stat().st_size, "sha256": sha(fields)},
              "products": {"bytes": products.stat().st_size, "sha256": sha(products)},
              "complete": rc == 0 and (run/"post/independent_checks.json").exists(),
              "scratch_removed": False}
    save(durable/"finished.json", result)
    shutil.rmtree(run.parent)
    result["scratch_removed"] = True; save(durable/"finished.json", result)
    return 0 if result["complete"] else 1

if __name__ == "__main__":
    root = Path(sys.argv[1])
    try: raise SystemExit(main(root))
    except Exception:
        save(root/"failure.json", {"utc": datetime.now(timezone.utc).isoformat(),
                                   "traceback": traceback.format_exc()})
        raise

