"""Coordinate-only nested-grid cut review; no mesh construction or solver run."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'analysis/matched_pair_36877_20261102/remote_audit.json'
OUT=ROOT/'analysis/pilot0_integrated_20261102'
z=np.array(json.loads(SOURCE.read_text())['vertical_vertices_m']);mid=(z[:-1]+z[1:])/2
assert len(z)==201 and np.all(np.diff(z)>0)
records=[]
for subdivisions in (1,2,3):
    fine=np.r_[np.concatenate([np.linspace(a,b,subdivisions+1)[:-1] for a,b in zip(z[:-1],z[1:])]),z[-1]]
    cuts=(fine[:-1]+fine[1:])/2
    for q in (.9,.75,.5):
        old=mid[np.argmin(abs(mid-q*.2))];near=cuts[np.argmin(abs(cuts-q*.2))]
        records.append(dict(subdivisions=subdivisions,Nz=len(fine)-1,target_z_over_H=q,baseline_actual_z_over_H=old/.2,
            fine_nearest_target_z_over_H=near/.2,cut_shift_m=near-old,distance_to_nearest_original_cut_m=float(min(abs(cuts-old))),
            old_cut_is_representable=bool(min(abs(cuts-old))<1e-14),bottom_cell_m=fine[1]-fine[0],top_cell_m=fine[-1]-fine[-2],
            status='coordinate_precheck_only_no_solver'))
with (OUT/'vertical_cut_precheck.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
(OUT/'vertical_cut_provenance.json').write_text(json.dumps(dict(source=str(SOURCE),sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
    interpretation='Hypothetical equal subdivision of actual baseline intervals; no claim that setting Grid.Cells1 alone implements these coordinates',
    solver_executed=False,transport_error_estimated=False),indent=2)+'\n',encoding='utf-8')
print('PASS coordinate precheck; Nz400 shifts cuts, Nz600 represents original cuts but default nearest-target selection can differ')
