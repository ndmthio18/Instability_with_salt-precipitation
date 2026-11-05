"""Run one immutable WP11 grid leg in node-local scratch and stage two archives."""
import csv,hashlib,json,os,resource,shutil,signal,socket,subprocess,sys,tarfile,time,traceback
from collections import deque
from datetime import datetime,timezone
from pathlib import Path

ELF="112076759fdff0942000307f729c2bd3c7e96963229b43c1b67c1f12bbf9df81"
PROFILE="d721e605d126e555f92b43200aa73b999ed4e2b61397edd1e98c0d9fddb36575"

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""):h.update(block)
    return h.hexdigest()
def save(path,value):
    # WP11 starts while the account is above its soft file-count quota and the grace
    # period may expire during R96. All durable targets are pre-created before submit;
    # overwrite the same inode instead of creating a temporary inode for every update.
    path=Path(path)
    with path.open("w") as f:
        json.dump(value,f,indent=2);f.flush();os.fsync(f.fileno())
def last_csv(path):
    if not path.exists():return 0,None
    with path.open(newline="") as f:
        reader=csv.DictReader(f);tail=deque(reader,maxlen=1);rows=reader.line_num-1
    return rows,(tail[0] if tail else None)
def archive(run,path,fields):
    with tarfile.open(path,"w:gz",compresslevel=1) as tar:
        for f in run.rglob("*"):
            if not f.is_file():continue
            isfield=f.parent==run and f.suffix in (".vtu",".pvd")
            if isfield==fields:tar.add(f,arcname=str(f.relative_to(run)))

def main(root,grid):
    root=Path(root);durable=root/grid;plan=json.loads((durable/"plan.json").read_text())
    threads=int(os.environ["WP11_THREADS"]);job=os.environ.get("SLURM_JOB_ID")
    assert job and threads==plan["threads"] and not (durable/"production_started.json").exists()
    scratch_root=Path(os.environ.get("TMPDIR","/tmp"));disk=shutil.disk_usage(scratch_root)
    mem=int(next(x.split()[1] for x in Path("/proc/meminfo").read_text().splitlines() if x.startswith("MemAvailable:")))*1024
    assert disk.free>=plan["minimum_scratch_free_bytes"] and mem>=plan["minimum_available_memory_bytes"]
    scratch=scratch_root/f"wp11_{grid}_{job}";run=scratch/plan["run_name"]
    run.mkdir(parents=True);shutil.copytree(root/"bundle",scratch/"bundle",symlinks=False)
    for n in ("params.input","eigenprofile.dat"):shutil.copy2(durable/n,run/n)
    assert sha(scratch/"bundle/bin/paper2_transport")==ELF
    assert sha(run/"params.input")==plan["input_sha256"] and sha(run/"eigenprofile.dat")==PROFILE
    started={"utc":datetime.now(timezone.utc).isoformat(),"job_id":job,"grid":grid,"node":socket.gethostname(),
             "affinity":sorted(os.sched_getaffinity(0)),"threads":threads,"scratch":str(run),
             "scratch_free_bytes":disk.free,"available_memory_bytes":mem,"elf_sha256":ELF,
             "input_sha256":plan["input_sha256"],"profile_sha256":PROFILE,
             "solver_timeout_s":plan["solver_timeout_s"]}
    save(durable/"production_started.json",started)
    env=dict(os.environ,DUMUX_NUM_THREADS=str(threads),OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",MKL_NUM_THREADS="1")
    def limit():resource.setrlimit(resource.RLIMIT_AS,(12*1024**3,12*1024**3))
    began=time.monotonic();timed_out=False
    with (run/"solver.log").open("w") as log:
        proc=subprocess.Popen([sys.executable,str(scratch/"bundle/launch_solver.py"),"params.input"],cwd=run,
                              stdout=log,stderr=subprocess.STDOUT,env=env,preexec_fn=limit)
        while proc.poll() is None:
            rows,tail=last_csv(run/"surface_events.csv")
            save(durable/"progress.json",{"utc":datetime.now(timezone.utc).isoformat(),"elapsed_s":time.monotonic()-began,
                                          "accepted_rows_including_initial":rows,"last_surface_event":tail})
            if time.monotonic()-began>plan["solver_timeout_s"]:
                timed_out=True;proc.send_signal(signal.SIGTERM);time.sleep(20)
                if proc.poll() is None:proc.kill()
                break
            time.sleep(60)
        rc=proc.wait()
    rows,tail=last_csv(run/"surface_events.csv")
    save(durable/"progress.json",{"utc":datetime.now(timezone.utc).isoformat(),"elapsed_s":time.monotonic()-began,
                                  "accepted_rows_including_initial":rows,"last_surface_event":tail})
    stop=json.loads((run/"stop_record.json").read_text()) if (run/"stop_record.json").exists() else None
    solver={"returncode":rc,"timed_out":timed_out,"wall_s":time.monotonic()-began,"stop":stop,
            "peak_rss_kib":resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss}
    save(run/"solver_finished.json",solver)
    accepted=rc==0 and stop and stop["reason"] in ("requested_horizon","surface_solubility_threshold","interior_solubility_threshold")
    if accepted:
        q=subprocess.run([sys.executable,str(root/"wp9_post.py"),str(run),str(run/"post")],capture_output=True,text=True,timeout=10800)
        (run/"post.log").write_text(q.stdout+q.stderr);solver["post_returncode"]=q.returncode
    fields=durable/"fields.tar.gz";products=durable/"products.tar.gz"
    archive(run,fields,True);archive(run,products,False)
    complete=accepted and (run/"post/independent_checks.json").exists() and json.loads((run/"post/independent_checks.json").read_text())["passed"]
    result={"solver":solver,"fields":{"bytes":fields.stat().st_size,"sha256":sha(fields)},
            "products":{"bytes":products.stat().st_size,"sha256":sha(products)},"complete":complete,"scratch_removed":False}
    save(durable/"finished.json",result);shutil.rmtree(scratch);result["scratch_removed"]=True;save(durable/"finished.json",result)
    return 0 if complete else 1

if __name__=="__main__":
    root=Path(sys.argv[1]);grid=sys.argv[2]
    try:raise SystemExit(main(root,grid))
    except Exception:
        save(root/grid/"failure.json",{"utc":datetime.now(timezone.utc).isoformat(),"traceback":traceback.format_exc()});raise
