"""Analytic-field check of the postprocessor, not a simulation or scientific result."""
import hashlib,json,subprocess,sys,tempfile
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory(prefix='paper2_fieldcheck_') as temp:
    base=Path(temp);source=base/'source';source.mkdir();out=base/'analysis';out.mkdir()
    params={'Grid.Positions0':'0 2','Grid.Positions1':'0 1','Grid.Cells0':'2','Grid.Cells1':'4',
        'Problem.Porosity':'.3','Problem.PcN':'2','Problem.PcAlpha':'.0001','Problem.Swr':'.1',
        'Problem.MolarMassH2O':'.01801528','Problem.MolarMassNaCl':'.05844','Problem.PureWaterDensity':'998.203',
        'Problem.DensityLinearCoefficient':'.670804','Problem.DensityQuadraticCoefficient':'.373854'}
    lines=[]
    for section in ['Grid','Problem']:
        lines.append('['+section+']');lines.extend(k.split('.')[1]+' = '+v for k,v in params.items() if k.startswith(section+'.'))
    (source/'test.input').write_text('\n'.join(lines)+'\n')
    xs=np.linspace(0,2,3);zs=np.linspace(0,1,5);xx,zz=np.meshgrid(xs,zs)
    sw0=.1+.9/np.sqrt(2);sig=[];pvd=[]
    for frame in range(101):
        x=.01+frame*1e-4*zz+1e-4*np.cos(np.pi*xx/2)
        w=x*.05844/((1-x)*.01801528+x*.05844);rho=998.203+1000*w*(.670804+.373854*w)
        arrays={'Coordinates':np.column_stack([xx.ravel(),zz.ravel(),np.zeros(xx.size)]),
                'moleFraction':x,'sw':np.full(xx.shape,sw0),'pw':np.full(xx.shape,90000),'rhoW':rho}
        xml='<VTKFile><UnstructuredGrid><Piece NumberOfCells="8">'+''.join('<DataArray Name="'+k+'" format="ascii">'+' '.join(format(v,'.17g') for v in a.ravel())+'</DataArray>' for k,a in arrays.items())+'</Piece></UnstructuredGrid></VTKFile>'
        name=f'test-{frame:05d}.vtu';(source/name).write_text(xml)
        pvd.append(f'<DataSet timestep="{frame}" file="{name}"/>')
        sig.append([frame,np.mean(x[-1]),np.std(x[-1]),3])
    (source/'test.pvd').write_text('<VTKFile><Collection>'+''.join(pvd)+'</Collection></VTKFile>')
    np.savetxt(source/'test_sigmatop.csv',sig)
    hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [source/'test.input',source/'test.pvd']}
    plan={'runs':[dict(id='analytic_field',case_id=-1,run_key='analytic_not_simulation',path=str(source),parameters=params,
            small_file_hashes=hashes,sigmatop_sha256=hashlib.sha256((source/'test_sigmatop.csv').read_bytes()).hexdigest())]}
    (out/'plan.json').write_text(json.dumps(plan))
    result=subprocess.run([sys.executable,str(ROOT/'scripts/dense_history_worker.py'),str(out)],capture_output=True,text=True)
    assert result.returncode==0,result.stderr
    audit=json.loads((out/'summary.json').read_text())[0]
    assert audit['status']=='field_input_consistent_binary_unverified',audit
    import csv
    frames=list(csv.DictReader((out/'runs/analytic_field/frame_metrics.csv').open()))
    expected_initial=.3*sw0*.01;expected_final=.3*sw0*.015
    assert abs(float(frames[0]['whole_salt_m'])-expected_initial)<1e-16
    assert abs(float(frames[-1]['whole_salt_m'])-expected_final)<1e-16
    assert abs(float(frames[-1]['top_x_mean_area'])-.02)<1e-16
    regional=list(csv.DictReader((out/'runs/analytic_field/inventory_timeseries.csv').open()))
    for row in regional:
        if int(row['frame'])==0:
            expected=.3*sw0*.01*(1-float(row['actual_z_over_H']))
            assert abs(float(row['inventory_salt_m'])-expected)<1e-16
    report=dict(passed=True,kind='analytic sampled field, no time integration or PDE solver',
                frames=101,known_full_inventory_initial=expected_initial,known_full_inventory_final=expected_final,
                known_final_surface_mean=.02,independent_quadrature_relative_error=audit['max_quadrature_relative_error'])
    dest=ROOT/'analysis/dense_history_20261102/postprocessing_selfcheck.json'
    dest.write_text(json.dumps(report,indent=2));print(json.dumps(report))
