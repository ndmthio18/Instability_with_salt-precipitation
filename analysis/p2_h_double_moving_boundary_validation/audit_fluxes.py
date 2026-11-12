"""Do not infer instantaneous flux convergence from inventory conservation."""
from pathlib import Path
import numpy as np
import pandas as pd
from analyze import CEQ
R=Path(__file__).resolve().parent
h=pd.read_csv(R/'results/history.csv.gz')
cv=pd.read_csv(R/'results/convergence.csv')
pairs=[(a.id,a.reference_id,'grid') for a in cv.itertuples()]
pairs += [(a,b,'time') for a,b in [('u0_n128','u0_strict'),('x008_n128','x008_strict'),('u0.6_n128','u0.6_strict'),('u0.6_n128','u0.6_strict_tighter')]]
rows=[]
for first,second,kind in pairs:
    a=h[h.id==first];b=h[h.id==second]
    for start in [1,60,3600]:
        q=a[(a.t>=start)&(a.t<=b.t.max())]
        if q.empty:continue
        row={'id':first,'reference_id':second,'kind':kind,'start_s':start,'common_end_s':q.t.max()}
        scale=CEQ*q.j/(1-CEQ)
        for field in ['Fplus','Fminus','Splus','Sminus','Jsoil','Cb']:
            delta=abs(q[field]-np.interp(q.t,b.t,b[field]))
            row[field+'_max_abs']=float(delta.max())
            row[field+'_scaled_max']=float((delta/(CEQ if field=='Cb' else scale)).max())
        rows.append(row)
pd.DataFrame(rows).to_csv(R/'results/flux_resolution.csv',index=False)
print('Wrote',len(rows),'windowed flux comparisons, including unresolved startup transients.')
