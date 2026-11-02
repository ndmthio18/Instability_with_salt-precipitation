"""Read-only frodo audit; no solver execution and no remote file writes."""
from pathlib import Path
import base64
import csv
import hashlib
import json
import math
import xml.etree.ElementTree as ET
import numpy as np

BASE = Path('/temp/zhixin/paper2_transport_20261102')
M0 = BASE/'unperturbed_36877_20261102_r01'
OLD = Path('/temp/zhixin/c1b/runs_600/shard-001/0017--case-36877--native')


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def field(p):
    names = ('Coordinates', 'moleFraction', 'sw', 'pw')
    a = {v.attrib.get('Name'): np.fromstring(v.text or '', sep=' ')
         for v in ET.parse(p).iter('DataArray')
         if v.attrib.get('Name') in names and v.attrib.get('format') == 'ascii'}
    xyz = a.pop('Coordinates').reshape(-1, 3)
    ix = np.lexsort((xyz[:, 0], xyz[:, 1]))
    return xyz[ix], {k: v[ix] for k, v in a.items()}


plan = json.loads((M0/'plan.json').read_text())
exe = Path(plan['executable']['path'])
profile = next(OLD.glob('*_eigenprofile.dat'))
table = np.loadtxt(profile)
assert len(table) >= 2 and np.all(np.diff(table[:, 0]) > 0)
f0 = M0/'unperturbed_36877-00000.vtu'
f1 = next(OLD.glob('*-00000.vtu'))
xyz, a0 = field(f0)
xyz1, a1 = field(f1)
assert np.array_equal(xyz, xyz1)
xs, xi = np.unique(xyz[:, 0], return_inverse=True)
zs, zi = np.unique(xyz[:, 1], return_inverse=True)
assert len(xs) == 601 and len(zs) == 201
assert len(xyz) == len(xs)*len(zs)
wx = np.diff(np.r_[xs[0], (xs[1:]+xs[:-1])/2, xs[-1]])
wz = np.diff(np.r_[zs[0], (zs[1:]+zs[:-1])/2, zs[-1]])
v = plan['parameters']
width = float(v['Grid.Positions0'].split()[-1])
phi = float(v['Problem.Porosity'])
phase = np.cos(2*np.pi/(width/28)*xyz[:, 0])
delta = 1e-6*np.interp(xyz[:, 1], table[:, 0], table[:, 1])*phase
predicted_x = float(v['Problem.SaltMoleFraction']) + delta
weight = phi*a0['sw']*wx[xi]*wz[zi]/width
meta = json.loads((M0/'transport_metadata.json').read_text())
budget = list(csv.DictReader((M0/'transport_budget.csv').open()))
initial = {r['region']: r for r in budget if r['step'] == '0'}
initial_checks = []
for cut in meta['cuts']:
    mask = xyz[:, 1] > cut['actual_z_over_H']*.2
    i0 = math.fsum((weight[mask]*a0['moleFraction'][mask]).tolist())
    ip = math.fsum((weight[mask]*predicted_x[mask]).tolist())
    direct = math.fsum((weight[mask]*delta[mask]).tolist())
    old = math.fsum((phi*a1['sw'][mask]*wx[xi[mask]]*wz[zi[mask]]/width*a1['moleFraction'][mask]).tolist())
    initial_checks.append(dict(region=cut['name'], actual_z_over_H=cut['actual_z_over_H'],
        M0_inventory_m=i0, planned_M1_inventory_m=ip,
        planned_M1_minus_M0_m=ip-i0, perturbation_integral_m=direct,
        legacy_M1_inventory_m=old, legacy_minus_M0_m=old-i0,
        M0_csv_reconstruction_error_m=i0-float(initial[cut['name']]['inventory_salt_m'])))

source_hashes = {name: sha(BASE/'source'/name) for name in plan['source_sha256']}
files = {}
for p in [M0/'params.input', M0/'plan.json', M0/'finished.json', M0/'transport_budget.csv',
          M0/'surface_events.csv', M0/'transport_metadata.json', f0, f1, profile,
          *[OLD/n for n in ['case_manifest.json', 'run_meta.json']], *OLD.glob('*.input')]:
    files[str(p)] = dict(bytes=p.stat().st_size, sha256=sha(p))
small = {p.name: base64.b64encode(p.read_bytes()).decode()
         for p in [profile, *OLD.glob('*.input'), OLD/'case_manifest.json', OLD/'run_meta.json']}
print(json.dumps(dict(
    scope='read_only_no_solver', M0_input_base64=base64.b64encode((M0/'params.input').read_bytes()).decode(),
    executable=dict(path=str(exe),sha256=sha(exe)),
    source_hashes=source_hashes, files=files, legacy_files_base64=small,
    initial_checks=initial_checks, vertical_vertices_m=zs.tolist(),
    horizontal_box_weighted_cos_mean=float(np.sum(wx*np.cos(2*np.pi/(width/28)*xs))/width),
    profile=dict(rows=len(table),z_min=float(table[0,0]),z_max=float(table[-1,0]),
                 A_min=float(table[:,1].min()),A_max=float(table[:,1].max()),
                 A_bottom=float(table[0,1]), A_top=float(table[-1,1])),
    initial_field=dict(legacy_vs_formula_max_x_error=float(np.max(np.abs(a1['moleFraction']-predicted_x))),
        legacy_vs_M0_max_pressure_error_Pa=float(np.max(np.abs(a1['pw']-a0['pw']))),
        legacy_vs_M0_max_Sw_error=float(np.max(np.abs(a1['sw']-a0['sw']))),
        planned_x_min=float(predicted_x.min()), planned_x_max=float(predicted_x.max()),
        planned_top_node_mean=float(predicted_x[zi==len(zs)-1].mean()),
        planned_top_weighted_mean=float(np.sum(predicted_x[zi==len(zs)-1]*wx)/width)),
    legacy_manifest=json.loads((OLD/'case_manifest.json').read_text()),
    legacy_log_tail=next(OLD.glob('*.log')).read_text()[-1800:]
), indent=2))
