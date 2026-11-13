import csv,json
from pathlib import Path
names=['fast_grid_x008_n512','parallel_x008_n512']
h=[]
for name in names:
 with open('raw/'+name+'_history.csv') as f:
  r=list(csv.DictReader(f));h.append({float(x['t']):x for x in r})
t=sorted(h[0].keys()&h[1].keys());equal=all(h[0][x]==h[1][x] for x in t)
result=dict(cases=names,common_end_s=t[-1],samples=len(t),bitwise_equal_saved_rows=equal)
Path('parallel_prefix_gate.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
assert equal and t[-1]>=1000
