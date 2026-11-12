"""Freeze P0 cases before launching; P1 is gated on P0 findings."""
import csv
from pathlib import Path

ROOT=Path(__file__).resolve().parent
rows=[]
for u in [0,.3,.6,.9,-1]:
    label='x008' if u<0 else f'u{u:g}'
    settings=[(32,1.5,1e-8,30),(64,1.5,1e-8,30),(128,1.5,1e-8,30),
              (64,1.5,1e-9,15),(64,1,1e-8,30),(64,2,1e-8,30)]
    for n,p,r,dt in settings:
        rows.append([f'{label}_n{n}_p{p:g}_r{r:g}_dt{dt}',n,p,r,dt,u,2,1,1])
# Follow-up resolving the observed 50 h top-layer under-resolution.
for n in [128,256,512]:
    rows.append([f'u0_n{n}_p2_r1e-08_dt30',n,2,1e-8,30,0,2,1,1])
with (ROOT/'cases.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.writer(f,lineterminator='\n');w.writerow(['id','n','power','rtol','maxstep','u','beta','av','jfactor']);w.writerows(rows)
print(f'{len(rows)} frozen P0 cases')
