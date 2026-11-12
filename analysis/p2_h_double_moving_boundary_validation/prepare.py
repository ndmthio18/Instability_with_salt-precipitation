import csv
from pathlib import Path
R=Path(__file__).resolve().parent
fields=['id','nc','nl','kappa','gap','rtol','maxstep','u','end']
pilot=[]
for n in [32,64,128]:pilot.append([f'pilot_u06_n{n}',n,n//4,1,1e-9,1e-8,10,.6,3600])
for gap in [1e-10,1e-8]:pilot.append([f'pilot_seed{gap:g}',64,16,1,gap,1e-8,10,.6,3600])
pilot.append(['pilot_u09_n64',64,16,1,1e-9,1e-8,10,.9,3600])
for row in pilot:row[0]+='_qjac'
cases=[]
for u in [0,.3,.6,.9,-1]:
    for n in [64,128,256]:
        name='x008' if u<0 else f'u{u:g}'
        cases.append([f'{name}_n{n}',n,n//4,1,1e-9,1e-8,30,u,180000])
for kb in [.5,2]:cases.append([f'u0.6_k{kb:g}',128,32,kb,1e-9,1e-8,30,.6,180000])
for gap in [1e-10,1e-8]:cases.append([f'u0.6_seed{gap:g}',128,32,1,gap,1e-8,30,.6,180000])
for u in [0,.6,-1]:
    name='x008' if u<0 else f'u{u:g}'
    cases.append([f'{name}_strict',128,32,1,1e-9,1e-9,15,u,180000])
cases.append(['u0.6_liquid64',128,64,1,1e-9,1e-8,30,.6,180000])
for name,rows in [('pilot_cases.csv',pilot),('cases.csv',cases)]:
    with (R/name).open('w',newline='') as f:
        w=csv.writer(f,lineterminator='\n');w.writerow(fields);w.writerows(rows)
print(len(pilot),'pilots;',len(cases),'production cases registered')
