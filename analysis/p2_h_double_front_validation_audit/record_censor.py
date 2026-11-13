import csv,json
from pathlib import Path
name='parallel_x008_n512'
with open('raw/'+name+'_history.csv') as f:rows=list(csv.DictReader(f))
last=float(rows[-1]['t']);prev=float(rows[-2]['t'])
Path('raw/'+name+'_status.txt').write_text(f'id={name}\nstatus=censored_duplicate_performance_trial\nlast_s={last}\nprevious_s={prev}\ntarget_s=160951.5\ntransfers=0\nreason=8-thread trial had no observed throughput gain on this allocation; identical prefix retained; primary serial run continues\n')
print(last)
