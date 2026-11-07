"""Independent pinned-Python operator and saved M3 trajectory comparisons."""
import csv
import hashlib
import json
import pathlib
import subprocess
import sys
import numpy as np

root=pathlib.Path.cwd()
ref=root/'reference/prat_dynamic_growth'
sys.path.insert(0,str(ref/'src'))
from prat_growth.dynamic import DynamicALEModel
from prat_growth.parameters import NumericalParameters, PhysicalParameters

checks=[]
for label,j,excess in [('growth',1200*.01/86400,2e-7),('reverse_mesh',0.,-2e-7)]:
    p=PhysicalParameters(U_evap=j/1200)
    m=DynamicALEModel(p,NumericalParameters(n_cells=12,grid_power=3),'B2')
    y=m.initial_state(); y[:12]=.43-.1*m.centers
    y[12:24]+=excess*np.exp(-10*(1-m.centers))
    path=root/'stage_a'/('operator_'+label+'.txt')
    path.write_text('12 '+repr(j)+'\n'+' '.join(map(repr,y.tolist())))
    run=subprocess.run([str(root/'build/operator_check'),str(path)],capture_output=True,text=True,check=True)
    (root/'stage_a'/('operator_'+label+'.out')).write_text(run.stdout)
    actual=np.fromstring(run.stdout,sep=' ')
    f=m.flux_state(y[:12],y[12:24],y[24])
    expected=np.concatenate(([f['c_top'],f['velocity']],f['liquid_flux'],f['salt_flux'],f['epsilon_mesh_flux'],m.rhs(0,y)))
    # Near-zero residuals are compared with a mixed absolute/relative bound.
    scaled=np.abs(actual-expected)/(1e-12+np.abs(expected))
    assert np.max(scaled)<1e-5,(label,float(np.max(scaled)))
    checks.append(dict(case=label,max_scaled_difference=float(np.max(scaled)),
                       max_absolute_difference=float(np.max(abs(actual-expected)))))

with (root/'stage_a/n640/history.csv').open() as f:
    rows=list(csv.DictReader(f))
archive=np.load(ref/'outputs/m3/runs/grid_640.npz')
meta=json.loads((ref/'outputs/m3/runs/grid_640.json').read_text())
assert hashlib.sha256((ref/'outputs/m3/runs/grid_640.npz').read_bytes()).hexdigest()==meta['data_sha256']
print('Frozen NPZ keys:',archive.files)
last={k:float(v) for k,v in rows[-1].items()}
targets=dict(height_m=.00164919695,epsilon_min=.23612005,internal_fraction=.25378548,velocity_m_s=3.00578141e-8)
differences={k:abs(last[k]-v)/abs(v) for k,v in targets.items()}
assert all(v<5e-5 for v in differences.values()),differences
assert last['time_s']==21600
first={k:float(v) for k,v in rows[0].items()}
salt_scale=first['dissolved_salt_kg_m2']+2170*.57*.001
water_scale=first['pore_water_kg_m2']
budgets=dict(salt=max(abs(float(r['salt_residual_kg_m2'])) for r in rows)/salt_scale,
             water=max(abs(float(r['water_residual_kg_m2'])) for r in rows)/water_scale)
assert max(budgets.values())<1e-8,budgets
times=np.array([float(r['time_s']) for r in rows])
assert np.array_equal(times,archive['time'])
height_error=float(np.max(np.abs(np.array([float(r['height_m']) for r in rows])/archive['height']-1)))
minimum_error=float(np.max(np.abs(np.array([float(r['epsilon_min']) for r in rows])-archive['epsilon'].min(axis=0))))
assert height_error<5e-5 and minimum_error<5e-5,(height_error,minimum_error)
with (root/'stage_a/n640/fields.csv').open() as f:
    fields=list(csv.DictReader(f))
epsilon_error=concentration_error=0.
for t in sorted({float(r['time_s']) for r in fields}):
    index=int(np.flatnonzero(archive['time']==t)[0])
    snap=[r for r in fields if float(r['time_s'])==t]
    epsilon_error=max(epsilon_error,float(np.max(np.abs(np.array([float(r['epsilon']) for r in snap])-archive['epsilon'][:,index]))))
    concentration_error=max(concentration_error,float(np.max(np.abs(np.array([float(r['concentration']) for r in snap])-archive['concentration'][:,index]))))
assert epsilon_error<5e-5 and concentration_error<1e-7,(epsilon_error,concentration_error)
result=dict(operator_checks=checks,final=last,relative_target_differences=differences,
            maximum_relative_budgets=budgets,reference_commit='8c0df6eddeea95a8a76c61f245cef2f891255c86',
            history_height_relative_error=height_error,history_epsilon_min_absolute_error=minimum_error,
            snapshot_epsilon_absolute_error=epsilon_error,snapshot_concentration_absolute_error=concentration_error,
            integrator='SUNDIALS 6.7.0 CVODE BDF/Newton in E,Q inventories; mapped tolerance coefficients, distinct error norm from primitive SciPy BDF')
(root/'stage_a/verification.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
