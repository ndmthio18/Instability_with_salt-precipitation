"""Independent review of actual M0/M1 initial VTK fields; never runs a solver."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np


def field(path):
    arrays = {a.get('Name'): np.fromstring(a.text or '', sep=' ') for a in ET.parse(path).iter('DataArray')
              if a.get('Name') in ('Coordinates','moleFraction','sw','pw') and a.get('format')=='ascii'}
    xyz=arrays.pop('Coordinates').reshape(-1,3)
    order=np.lexsort((xyz[:,0],xyz[:,1]))
    return xyz[order], {k:v[order] for k,v in arrays.items()}


def review(m0, m1):
    plan=json.loads((m1/'plan.json').read_text()); p=plan['parameters']
    f0=next(m0.glob('*-00000.vtu')); f1=next(m1.glob('*-00000.vtu'))
    xyz,a=field(f0); xyz1,b=field(f1)
    assert np.array_equal(xyz,xyz1), 'Actual grids differ'
    xs,xi=np.unique(xyz[:,0],return_inverse=True); zs,zi=np.unique(xyz[:,1],return_inverse=True)
    assert len(xs)==601 and len(zs)==201 and len(xyz)==601*201
    width=xs[-1]-xs[0];height=zs[-1]-zs[0]
    wx=np.diff(np.r_[xs[0],(xs[1:]+xs[:-1])/2,xs[-1]])
    wz=np.diff(np.r_[zs[0],(zs[1:]+zs[:-1])/2,zs[-1]])
    tab=np.loadtxt(m1/'eigenprofile.dat')
    expected=float(p['Problem.SaltMoleFraction'])+1e-6*np.interp(xyz[:,1],tab[:,0],tab[:,1])*np.cos(2*np.pi/(width/28)*xyz[:,0])
    checks=dict(actual_grid_600x200=True,
        actual_pressure_identical=bool(np.array_equal(a['pw'],b['pw'])),
        actual_saturation_identical=bool(np.array_equal(a['sw'],b['sw'])),
        actual_M1_formula=bool(np.max(abs(b['moleFraction']-expected))<1e-14),
        M0_uniform=bool(np.ptp(a['moleFraction'])<1e-12))
    weights=float(p['Problem.Porosity'])*wx[xi]*wz[zi]/width
    initial={}
    for label,run in [('M0',m0),('M1',m1)]:
        with (run/'transport_budget.csv').open() as f:
            initial[label]={r['region']:r for r in csv.DictReader(f) if r['step']=='0'}
    meta=json.loads((m1/'transport_metadata.json').read_text())
    inventories=[]
    for cut in meta['cuts']:
        mask=xyz[:,1]>cut['actual_z_over_H']*height
        vals=[math.fsum((weights[mask]*v['sw'][mask]*v['moleFraction'][mask]).tolist()) for v in (a,b)]
        errors=[vals[i]-float(initial[label][cut['name']]['inventory_salt_m']) for i,label in enumerate(('M0','M1'))]
        inventories.append(dict(region=cut['name'],actual_z_over_H=cut['actual_z_over_H'],
            M0_actual_initial_inventory_m=vals[0],M1_actual_initial_inventory_m=vals[1],
            measured_initial_inventory_difference_m=vals[1]-vals[0],
            M0_reconstruction_error_m=errors[0],M1_reconstruction_error_m=errors[1]))
    checks['initial_inventory_reconstruction']=all(abs(r[k])<1e-14 for r in inventories for k in ('M0_reconstruction_error_m','M1_reconstruction_error_m'))
    return dict(passed=all(checks.values()),checks=checks,inventories=inventories,
        M1_max_formula_error=float(np.max(abs(b['moleFraction']-expected))),
        M1_initial_min=float(b['moleFraction'].min()),M1_initial_max=float(b['moleFraction'].max()),
        sources={str(q):hashlib.sha256(q.read_bytes()).hexdigest() for q in (f0,f1,m1/'eigenprofile.dat')},
        scope='Measured new M1 and existing M0 fields, not substituted historical nx28 inventories')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--m0',type=Path,required=True);ap.add_argument('--m1',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args(); result=review(args.m0,args.m1)
    args.out.write_text(json.dumps(result,indent=2)+'\n')
    raise SystemExit(0 if result['passed'] else 1)
