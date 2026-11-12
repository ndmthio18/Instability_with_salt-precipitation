"""Adaptive P1 spatial checks after the measured 128/256 differences."""
from pathlib import Path
import csv
R=Path(__file__).resolve().parent
with (R/'cases_events.csv').open() as f:
    reader=csv.DictReader(f);fields=reader.fieldnames
    rows=[r for r in reader if r['n']=='256']
with (R/'cases_p1.csv').open() as f:
    rows += [r for r in csv.DictReader(f) if r['id']=='p1_j1000_u0_n256']
for r in rows:r['id']=r['id'].replace('n256','n512');r['n']='512'
with (R/'cases_refine.csv').open('w',newline='') as f:
    writer=csv.DictWriter(f,fields,lineterminator='\n');writer.writeheader();writer.writerows(rows)
