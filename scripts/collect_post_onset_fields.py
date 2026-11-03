import sys;sys.path.insert(0,'scripts')
from analyze_post_onset import *
from collect_post_onset import remote
rows=readcsv(OUT/'run_metrics.csv');inv=json.load(open(ROOT/'analysis/phase0_20261102/remote_frodo_inventory.json'))['runs'];bykey={r['run_key']:r for r in inv}
selected=[]
for cid in ['25748','27710','31932','30858','47647','45164','29422','35733','19674']:
 rr=[r for r in rows if r['case_id']==cid and r['source_group']=='main' and r['series_status']=='hash_verified']
 r=max(rr,key=lambda r:num(r.get('G_end'),0));selected.append(bykey[r['run_key']])
# Deliberate quarantined comparator, never admitted to physical claims.
selected.append(next(r for r in inv if r['run_key']=='frodo:hf1/runs/shard-019/0079--case-38852--native'))
savejson(OUT/'sparse_field_selection.json',dict(reason='growth leaders, relative-amplitude leaders, near-threshold negative and known provenance-conflict comparator; dense cases reused separately',runs=selected))
cache=ROOT/'tmp/post_onset_amplification_20261103/fields';cache.mkdir(parents=True,exist_ok=True)
code=(ROOT/'scripts/post_onset_field_worker.py').read_text(encoding='utf-8-sig')
for i,s in enumerate(selected):
 target=cache/(str(i)+'.json')
 if target.exists() and json.loads(target.read_text())['run_key']==s['run_key']:continue
 result=remote(code+'\nprint(json.dumps(inspect('+repr(s)+')))\n',timeout=240)
 import base64
 data=base64.b64decode(result.pop('fields_base64'));(cache/(str(i)+'.npz')).write_bytes(data);result['npz_sha256']=hashlib.sha256(data).hexdigest();savejson(target,result)
 print(i,s['run_key'],[(r['consistent'],r['x_max']) for r in result['records']],flush=True)
