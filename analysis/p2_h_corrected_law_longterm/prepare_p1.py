import csv
from pathlib import Path
R=Path(__file__).resolve().parent
rows=[]
settings=[('beta05',.5,1,1),('beta_unit_da',2/(3.1750355738459297**.5),1,1),
          ('av_unit_da',2,3.1750355738459297,1),('av16',2,16,1),
          ('j1000',2,1,1000),('j10000',2,1,10000)]
for label,beta,av,jf in settings:
    for u in [0,.6]:
        for n in [128,256]:
            rows.append([f'p1_{label}_u{u:g}_n{n}',n,2,1e-8,30/jf,u,beta,av,jf])
with (R/'cases_p1.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.writer(f,lineterminator='\n');w.writerow(['id','n','power','rtol','maxstep','u','beta','av','jfactor']);w.writerows(rows)
print(len(rows))
