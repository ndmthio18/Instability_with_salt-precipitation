import csv,json,os
from pathlib import Path

def rows(name):
 with open('raw/'+name+'_history.csv') as f:return {float(r['t']):{k:float(v) for k,v in r.items()} for r in csv.DictReader(f)}
def complete(name,end):return dict(x.split('=',1) for x in Path('raw/'+name+'_status.txt').read_text().splitlines())['status']=='end' and max(rows(name))>=end
prefix=os.getenv('CASE_PREFIX','setup_')
ceq=.264940869768;cr=2165.;phi=.8898899218857486
vel=[]
for lab,end in [('u06',180000),('u09',180000),('x008',160951.5)]:
 names=[f'{prefix}p0_ld_{lab}_r1e-09',f'{prefix}p0_ld_{lab}_r1e-10'];a,b=map(rows,names);t=sorted(a.keys()&b.keys());worst=0
 for ti in t:
  ref=ceq*b[ti]['j']/((1-ceq)*cr*phi)
  worst=max(worst,*[abs(a[ti][k]-b[ti][k])/ref for k in ['vt','vb','vnet']])
 vel.append(dict(case=lab,complete=all(complete(n,end) for n in names),worst_velocity_normalized=worst))
flux=[]
for lab in ['u06','x008']:
 names=[f'{prefix}early_{lab}_n{n}_b1e-09' for n in [128,256,512]]+[f'{prefix}early_{lab}_n{n}_tight' for n in [128,256,512]]
 hs=list(map(rows,names));tt=sorted(set.intersection(*[set(h) for h in hs]));err=[]
 for t in tt:
  e=0
  for i,j in [(0,1),(1,2),(0,3),(1,4),(2,5)]:
   ref=ceq*hs[j][t]['j']/(1-ceq)
   e=max(e,*[abs(hs[i][t][k]-hs[j][t][k])/ref for k in ['Fplus','Fminus']])
  err.append(e)
 earliest=next((t for i,t in enumerate(tt) if max(err[i:])<.01),None)
 flux.append(dict(case=lab,complete=all(complete(n,100) for n in names),earliest_s=earliest,common_end=tt[-1]))
all_early=list(csv.DictReader(open('fast_p01_cases.csv' if prefix=='fast_' else 'retry_cases.csv')))
failures=[r['id'] for r in all_early if not complete(r['id'],float(r['end']))]
result=dict(velocity=vel,flux=flux,failures=failures,P0_pass=all(v['complete'] and v['worst_velocity_normalized']<.01 for v in vel),P1_pass=all(v['complete'] and v['earliest_s'] is not None for v in flux) and not failures)
Path('fast_p01_gate.json' if prefix=='fast_' else 'p01_gate.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
