import csv,json
from pathlib import Path
name='grid_x008_n512'
with open('raw/'+name+'_history.csv') as f:rows=list(csv.DictReader(f))
last=float(rows[-1]['t']);prev=float(rows[-2]['t'])
Path('raw/'+name+'_status.txt').write_text(f'id={name}\nstatus=censored_superseded_numerical_version\nlast_s={last}\nprevious_s={prev}\ntarget_s=160951.5\ntransfers=0\nreason=v6 duplicate replaced by uniformly matched v7 primary; full v6 P0 bridge is retained; this accepted prefix is not a numerical failure\n')
print(last)
