"""Run/retrieve bounded field analysis on frodo; never executes a solver."""
import base64,hashlib,json,sys
from pathlib import Path
from collect_post_onset import remote
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'analysis/wp6_review_20261103';CACHE=ROOT/'tmp/wp6_review_20261103'
REMOTE='/temp/zhixin/paper2_wp6_review_20261103_r01'


def main(role):
    assert role in ['baseline','half']
    run=('/temp/zhixin/paper2_wp6_execution_20261103_r01/runs/WP6-M0-extend' if role=='baseline' else REMOTE+'/runs/WP6-M0-dt_half')
    data=(ROOT/'scripts/wp6_background_structure.py').read_bytes();digest=hashlib.sha256(data).hexdigest()
    code="""import json,base64,hashlib,subprocess,sys,os
from pathlib import Path
r=Path(REMOTE);script=r/'wp6_background_structure.py';data=base64.b64decode(DATA)
assert not script.exists() or script.read_bytes()==data
if not script.exists():script.write_bytes(data)
assert hashlib.sha256(script.read_bytes()).hexdigest()==DIGEST
assert (Path(RUN)/'finished.json').exists(), 'solver/archive postprocessing not finished'
out=r/('background_'+ROLE);assert os.getloadavg()[0]<8
os.sched_setaffinity(0,{2})
p=subprocess.run([sys.executable,str(script),RUN,str(out)],capture_output=True,text=True,timeout=600)
assert p.returncode==0,p.stderr
files={f.name:dict(data=base64.b64encode(f.read_bytes()).decode(),sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in out.iterdir() if f.is_file()}
print(json.dumps(dict(role=ROLE,files=files,stdout=p.stdout)))
""".replace('REMOTE',repr(REMOTE)).replace('DATA',repr(base64.b64encode(data).decode())).replace('DIGEST',repr(digest)).replace('RUN',repr(run)).replace('ROLE',repr(role))
    result=remote(code,timeout=650)
    records={}
    for name,v in result['files'].items():
        b=base64.b64decode(v['data']);assert hashlib.sha256(b).hexdigest()==v['sha256']
        d=(CACHE if name.endswith('.npz') else OUT)/('background_'+role);d.mkdir(exist_ok=True,parents=True);(d/name).write_bytes(b)
        records[name]=dict(sha256=v['sha256'],bytes=len(b))
    (OUT/('background_'+role)/'transfer.json').write_text(json.dumps(records,indent=2))
    print(result['stdout'])


if __name__=='__main__':main(sys.argv[1])
