"""Reduce one WP11 same-node, same-grid one-versus-four-thread benchmark."""
import hashlib,json,os,socket,sys
from pathlib import Path

scratch,out,grid,job=Path(sys.argv[1]),Path(sys.argv[2]),sys.argv[3],sys.argv[4]
files=["surface_events.csv","transport_budget.csv","weighted_online.csv"]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
runs={}
for n in (1,4):
    r=scratch/f"t{n}";start,end=map(float,(r/"times.txt").read_text().splitlines())
    stop=json.loads((r/"stop_record.json").read_text())
    runs[str(n)]={"threads":n,"wall_s":end-start,"s_per_step":(end-start)/stop["accepted_steps"],
                  "stop":stop,"files":{f:sha(r/f) for f in files}}
identical={f:runs["1"]["files"][f]==runs["4"]["files"][f] for f in files}
result={"scope":"infrastructure benchmark only; no scientific result","grid":grid,"steps":100,
        "job_id":job,"node":socket.gethostname(),"affinity":sorted(os.sched_getaffinity(0)),
        "elf_sha256":sha(scratch/"bundle/bin/paper2_transport"),"runs":runs,
        "bitwise_identical":identical,"all_bitwise_identical":all(identical.values()),
        "speedup_1_to_4":runs["1"]["wall_s"]/runs["4"]["wall_s"]}
(out/"benchmark_result.json").write_text(json.dumps(result,indent=2))
for n in (1,4):(out/f"benchmark_t{n}.log").write_bytes((scratch/f"t{n}/solver.log").read_bytes())
print(json.dumps(result))
