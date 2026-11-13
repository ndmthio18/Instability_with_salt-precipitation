import csv,json
from pathlib import Path
ceq=.264940869768;cr=2165.;phi=.8898899218857486
def read(n):return {float(r['t']):{k:float(v) for k,v in r.items()} for r in csv.DictReader(open('raw/'+n+'_history.csv'))}
report=[]
for lab in ['u00','u03','u06','u09','x008']:
 a=f'fast_production_{lab}_n256' if lab in ['u00','u03'] else f'fast_p2_{lab}_n256_k1' if lab!='x008' else 'fast_grid_x008_n256'
 b=f'fast_production_{lab}_n512' if lab in ['u00','u03'] else f'fast_grid_{lab}_n512'
 if not all(Path('raw/'+c+'_status.txt').exists() for c in [a,b]):continue
 x,y=map(read,[a,b]);tt=sorted(x.keys()&y.keys());row={'case':lab,'end_s':tt[-1]}
 for k in ['vt','vb','vnet','zt','zb','hc']:
  row[k+'_max_scaled']=max(abs(x[t][k]-y[t][k])/(ceq*y[t]['j']/((1-ceq)*cr*phi) if k.startswith('v') else max(abs(y[t][k]),1e-9)) for t in tt)
 report.append(row)
Path('partial_grid.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
