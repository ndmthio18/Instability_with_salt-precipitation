"""Execute bounded postprocessing and collect small products, no new physics."""
import base64,hashlib,json,sys
from pathlib import Path
from collect_post_onset import remote
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'analysis/wp6_4T_exploratory_20261103';CACHE=ROOT/'tmp/wp6_4T_exploratory_20261103';REMOTE='/temp/zhixin/paper2_wp6_4T_exploratory_20261103_r01'
def main(role='M1'):
    assert role in ['M0','M1'];folder='spatial' if role=='M1' else 'spatial_M0'
    data={name:base64.b64encode((ROOT/'scripts'/name).read_bytes()).decode() for name in ['wp6_4t_spatial.py','wp6_4t_metrics.py']}
    code='''import json,base64,hashlib,os,subprocess,sys
from pathlib import Path
r=Path(REMOTE);script=r/'wp6_4t_spatial.py'
for name,encoded in DATA.items():
 p=r/name;content=base64.b64decode(encoded);assert not p.exists() or p.read_bytes()==content
 if not p.exists():p.write_bytes(content)
run=r/'runs'/('WP6-'+ROLE+'-explore4T');assert (run/'finished.json').exists()
cpu=json.loads((r/'explore4T_manager_resources.json').read_text())['cpus'][0 if ROLE=='M0' else 1];os.sched_setaffinity(0,{cpu})
out=r/FOLDER;p=subprocess.run([sys.executable,str(script),str(run),str(out)],capture_output=True,text=True,timeout=900);assert p.returncode==0,p.stderr
files={f.name:dict(data=base64.b64encode(f.read_bytes()).decode(),sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in out.iterdir() if f.is_file()}
print(json.dumps(dict(files=files,summary=p.stdout)))
'''.replace('REMOTE',repr(REMOTE)).replace('DATA',repr(data)).replace('ROLE',repr(role)).replace('FOLDER',repr(folder))
    r=remote(code,timeout=950);receipt={}
    for name,v in r['files'].items():
        b=base64.b64decode(v['data']);assert hashlib.sha256(b).hexdigest()==v['sha256'];d=(CACHE if name.endswith('.npz') else OUT)/folder;d.mkdir(parents=True,exist_ok=True);(d/name).write_bytes(b);receipt[name]=dict(sha256=v['sha256'],bytes=len(b))
    (OUT/folder/'transfer.json').write_text(json.dumps(receipt,indent=2));print(r['summary'])
if __name__=='__main__':main(sys.argv[1] if len(sys.argv)>1 else 'M1')
