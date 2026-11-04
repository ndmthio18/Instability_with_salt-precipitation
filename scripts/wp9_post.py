"""Stride-aware independent VTU quadrature and accepted-step budget checks.

This is `wp6_post.py` with one change: it no longer assumes that archived frame i is
accepted step i. `single_wavelength_post.post` already records each frame's true step,
matched through the sigma-top clock, and that step is used here to address the budget
and online rows. The R96 leg advances with dt/2 and archives every second step, so the
original indexing would have compared frame i against step i and failed.

Two consequences, both strictly stronger than the original:

  * the cumulative flux re-integration now walks every accepted step, not only the
    archived ones, so it reproduces the recorded cumulative column exactly;
  * the independent inventory closure is evaluated at the archived steps, against the
    cumulative sum carried through the steps in between.

For a run that archives every step this reduces to the original computation; the
equivalence is checked against the committed R48 output in `wp9_post_equivalence.json`.
Thresholds are unchanged.
"""
import csv,json,sys,hashlib,xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
from single_wavelength_post import post,params,widths,savecsv


def main(run,out):
    run=Path(run);out=Path(out);audit=post(run,out)
    rows=list(csv.DictReader((run/'transport_budget.csv').open()))
    on=list(csv.DictReader((run/'weighted_online.csv').open()))
    weighted=list(csv.DictReader((out/'weighted_diagnostics.csv').open()))
    meta=json.loads((run/'transport_metadata.json').read_text());cuts=meta['cuts'];assert len(cuts)==8
    p=params(run/'params.input');phi=float(p['Problem.Porosity'])
    by={}
    for r in rows:by.setdefault(int(r['step']),{})[r['region']]=r
    online_by={int(r['step']):r for r in on}
    frame_step=[int(r['step']) for r in weighted]
    assert frame_step==sorted(frame_step) and len(set(frame_step))==len(frame_step)
    stride=(frame_step[1]-frame_step[0]) if len(frame_step)>1 else 1
    assert all(b-a==stride for a,b in zip(frame_step,frame_step[1:])),'frames are not evenly spaced'

    inv=[];online_error={k:0 for k in ['weighted_top_mean','weighted_top_sigma',
                                       'weighted_fourier','weighted_projection']};grids=None
    for i,ds in enumerate(ET.parse(next(run.glob('*.pvd'))).getroot().iter('DataSet')):
        step=frame_step[i]
        a={v.attrib['Name']:np.fromstring(v.text or '',sep=' ')
           for v in ET.parse(run/ds.attrib['file']).getroot().iter('DataArray')
           if v.attrib.get('Name') in ['Coordinates','moleFraction','sw']}
        xyz=a['Coordinates'].reshape(-1,3)
        if grids is None:
            xs=np.unique(xyz[:,0]);zs=np.unique(xyz[:,1])
            order=np.argsort(np.searchsorted(zs,xyz[:,1])*len(xs)+np.searchsorted(xs,xyz[:,0]))
            vol=widths(zs)[:,None]*widths(xs)[None,:]/(xs[-1]-xs[0]);grids=True
        x=a['moleFraction'][order].reshape(len(zs),len(xs));sw=a['sw'][order].reshape(x.shape)
        for cdef in cuts:
            mask=zs>zs[0]+cdef['actual_z_over_H']*(zs[-1]-zs[0]);r=by[step][cdef['name']]
            salt=float(np.sum((phi*sw*x*vol)[mask]));water=float(np.sum((phi*sw*vol)[mask]))
            inv.append(dict(step=step,frame=i,time_s=float(r['time_s']),region=cdef['name'],
                independent_salt_m=salt,independent_water_m=water,
                salt_difference_m=salt-float(r['inventory_salt_m']),
                water_difference_m=water-float(r['inventory_water_m'])))
        o=online_by[step]
        assert abs(float(o['time_s'])-float(weighted[i]['time_s']))<1e-7
        for k in online_error:
            if k=='weighted_projection' and o['projection_defined']=='0':continue
            online_error[k]=max(online_error[k],abs(float(o[k])-float(weighted[i][k])))

    budget_checks={}
    for cdef in cuts:
        rr=[r for r in rows if r['region']==cdef['name']]
        ii={r['step']:r for r in inv if r['region']==cdef['name']}
        first=ii[frame_step[0]]['independent_salt_m'];firstw=ii[frame_step[0]]['independent_water_m']
        cumS=cumW=0.;maxS=maxW=maxCS=maxCW=0.
        for r in rr:
            step=int(r['step']);dt=float(r['dt_s'])
            cumS+=dt*(float(r['cut_salt_total_down_m_s'])+float(r['other_salt_out_m_s']))
            cumW+=dt*(float(r['cut_water_down_m_s'])+float(r['other_water_out_m_s']))
            maxCS=max(maxCS,abs(cumS-float(r['cumulative_cut_salt_m'])))
            maxCW=max(maxCW,abs(cumW-float(r['cumulative_net_water_out_m'])))
            if step in ii:
                maxS=max(maxS,abs(ii[step]['independent_salt_m']-first+cumS))
                maxW=max(maxW,abs(ii[step]['independent_water_m']-firstw+cumW))
        budget_checks[cdef['name']]=dict(salt_closure_m=maxS,water_closure_m=maxW,
            cumulative_salt_reintegration_error=maxCS,cumulative_water_reintegration_error=maxCW)

    result=dict(frames=len(weighted),recorded_steps=len(on),frame_stride=stride,
        max_inventory_salt_difference_m=max(abs(r['salt_difference_m']) for r in inv),
        max_inventory_water_difference_m=max(abs(r['water_difference_m']) for r in inv),
        online_absolute_errors=online_error,independent_budgets=budget_checks,
        scope='independent VTU storage + independent integration of recorded conservative flux; not a second flux-kernel implementation')
    result['passed']=(result['max_inventory_salt_difference_m']<1e-15
        and result['max_inventory_water_difference_m']<1e-12
        and max(online_error.values())<1e-12
        and all(v['salt_closure_m']<1e-14 and v['water_closure_m']<1e-11
                and v['cumulative_salt_reintegration_error']<1e-18
                and v['cumulative_water_reintegration_error']<1e-16
                for v in budget_checks.values()))
    savecsv(out/'independent_inventories.csv',inv)
    (out/'independent_checks.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result));assert result['passed']


if __name__=='__main__':main(sys.argv[1],sys.argv[2])
