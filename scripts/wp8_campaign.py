"""WP8 postprocessing only, using existing archived or resident fields."""
import subprocess,json,hashlib,tarfile,tempfile,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'analysis/wp8_zero_mode_20261103';CACHE=ROOT/'tmp/wp8_zero_mode_20261103';REMOTE='/tmp/zhixin_paper2_wp8_zero_mode_20261103';OLD='/tmp/zhixin_paper2_wp7_shortwave_20261103_r01'
RUNS={'R24_M0':'analysis/wp6_4T_exploratory_20261103/runs/WP6-M0-explore4T','R24_M1':'analysis/wp6_4T_exploratory_20261103/runs/WP6-M1-explore4T','R48_M0':'analysis/wp7_shortwave_validation_20261103/runs/WP7-R48-M0-4T','R48_M1':'analysis/wp7_shortwave_validation_20261103/runs/WP7-R48-M1-4T','B':'analysis/wp7_shortwave_validation_20261103/runs/WP7-1D-background-T'}
def run(a):return subprocess.check_output(a,text=True)
def ssh(s):return run(['ssh','-o','ConnectTimeout=15','sal',s])
def scp(a,b):subprocess.run(['scp',str(a),str(b)],check=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
def launch():
    OUT.mkdir(exist_ok=True,parents=True);CACHE.mkdir(exist_ok=True,parents=True)
    assert not (OUT/'jobs.json').exists();ssh(f'mkdir -p {REMOTE}')
    for n in ['wp8_fields.py','single_wavelength_post.py','wp6_4t_metrics.py','wp8_post.sh']:scp(ROOT/'scripts'/n,f'sal:{REMOTE}/{n}')
    jobs={}
    for key,rel in RUNS.items():
        local=ROOT/rel;name=local.name;remote=f'{OLD}/runs/{name}'
        if key.startswith('R24'):
            remote=f'{REMOTE}/inputs/{name}';ssh(f'mkdir -p {remote}')
            archive=ROOT/'tmp/wp6_4T_exploratory_20261103'/f'{name}_fields.tar.gz'
            assert sha(archive)==json.loads((local/'finished.json').read_text())['archive']['sha256']
            scp(archive,f'sal:{REMOTE}/{archive.name}');ssh(f'tar -xzf {REMOTE}/{archive.name} -C {remote}')
            with tempfile.TemporaryDirectory() as td:
                arc=Path(td)/'meta.tar.gz'
                with tarfile.open(arc,'w:gz') as tar:
                    for n in ['params.input','surface_events.csv','post/field_hashes.json']:tar.add(local/n,arcname=n)
                scp(arc,f'sal:{REMOTE}/{key}_meta.tar.gz');ssh(f'tar -xzf {REMOTE}/{key}_meta.tar.gz -C {remote}')
        jobs[key]=ssh(f'cd {REMOTE} && sbatch --parsable --job-name=WP8-{key} wp8_post.sh {REMOTE} {remote} {key}').strip()
        save(OUT/'jobs.json',jobs)
    print(jobs)
def status():
    jobs=json.loads((OUT/'jobs.json').read_text());print(ssh('squeue -h -j '+','.join(jobs.values())+" -o '%i|%T|%M|%R'"))
def collect():
    jobs=json.loads((OUT/'jobs.json').read_text());assert not ssh('squeue -h -j '+','.join(jobs.values())).strip()
    for key in jobs:
        dest=OUT/key;dest.mkdir(exist_ok=True);cache=CACHE/key;cache.mkdir(exist_ok=True)
        for n in ['audit.json','timeseries.csv']:scp(f'sal:{REMOTE}/{key}/{n}',dest/n)
        scp(f'sal:{REMOTE}/{key}/profiles_spectra.npz',cache/'profiles_spectra.npz')
        audit=json.loads((dest/'audit.json').read_text());assert sha(cache/'profiles_spectra.npz')==audit['files']['profiles_spectra.npz'];assert sha(dest/'timeseries.csv')==audit['files']['timeseries.csv']
    print('all five products transferred and SHA verified')
if __name__=='__main__':globals()[sys.argv[1]]()
