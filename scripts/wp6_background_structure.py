"""Read every archived accepted M0 VTU; hash-check and quantify horizontal noise.

Runs on the compute host, writes small summaries and an ignored plotting cache.
This diagnoses numerical trajectories, not a physical instability growth rate.
"""
import csv,hashlib,json,sys,time
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def widths(x): return np.diff(np.r_[x[0],(x[1:]+x[:-1])/2,x[-1]])
def main(run,out):
    run=Path(run);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    hashes=json.loads((run/'post/field_hashes.json').read_text())
    weighted=list(csv.DictReader((run/'post/weighted_diagnostics.csv').open()))
    assert len(hashes)==len(weighted)
    rows=[];top=[];mode=[];energies=[];depth=[];snaps={};grid=None
    for k,(rec,w) in enumerate(zip(hashes,weighted)):
        p=Path(rec['path']);assert p.resolve().is_relative_to(run.resolve())
        raw=p.read_bytes();assert hashlib.sha256(raw).hexdigest()==rec['sha256']
        a={v.attrib['Name']:np.fromstring(v.text or '',sep=' ') for v in ET.fromstring(raw).iter('DataArray') if v.attrib.get('Name') in ['Coordinates','moleFraction','pw','sw','rhoW']}
        xyz=a['Coordinates'].reshape(-1,3)
        if grid is None:
            grid=xyz;xs=np.unique(xyz[:,0]);zs=np.unique(xyz[:,1]);nx=len(xs);nz=len(zs);wx=widths(xs);wz=widths(zs);W=wx.sum();H=wz.sum()
            order=np.argsort(np.searchsorted(zs,xyz[:,1])*nx+np.searchsorted(xs,xyz[:,0]));assert len(xyz)==nx*nz
            cm=np.cos(np.pi*np.arange(1,nx)[:,None]*(xs-xs[0])/W)
        else:assert np.array_equal(xyz,grid)
        row=dict(frame=k,time_s=float(w['time_s']))
        for name in ['moleFraction','pw','sw','rhoW']:
            q=a[name][order].reshape(nz,nx);mean=q@wx/W;d=q-mean[:,None];sd=np.sqrt(d*d@wx/W)
            row.update({name+'_top_sigma':sd[-1],name+'_volume_rms':np.sqrt(np.sum(sd*sd*wz)/H),name+'_max_endpoint':np.max(abs(q[:,0]-q[:,-1])),name+'_max_anomaly':np.max(abs(d)),name+'_max_mirror':np.max(abs(q-q[:,::-1]))})
            if name=='moleFraction':
                coef=2*(d*wx)@cm.T/W;en=.5*coef**2;en[:,-1]*=.5;en=np.sum(en*wz[:,None],axis=0)/H
                top.append(d[-1]);mode.append(coef[-1]);energies.append(en);depth.append(sd)
                row['x_top_signed_endpoint']=q[-1,0]-q[-1,-1]
                row['x_peak_endpoint_depth_m']=zs[-1]-zs[np.argmax(abs(q[:,0]-q[:,-1]))]
                row['x_peak_sigma_depth_m']=zs[-1]-zs[np.argmax(sd)]
                row['x_top_range_ulps']=(q[-1].max()-q[-1].min())/np.spacing(mean[-1])
                assert abs(sd[-1]-float(w['weighted_top_sigma']))<1e-18
                if k in [0,1,len(hashes)//2,len(hashes)-1]:snaps['x_anomaly_'+str(k)]=d
        rows.append(row)
    with (out/'noise_timeseries.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    top=np.array(top);mode=np.array(mode);energies=np.array(energies);depth=np.array(depth)
    np.savez_compressed(out/'noise_fields.npz',time_s=[r['time_s'] for r in rows],xs=xs,zs=zs,top_anomaly=top,top_neumann_coeff=mode,volume_mode_energy=energies,depth_sigma=depth,**snaps)
    dominant=np.argmax(energies[-1])+1; total=energies.sum(axis=1)
    summary=dict(run=str(run),frames=len(rows),all_vtu_hashes_verified=True,input_sha256=sha(run/'params.input'),field_manifest_sha256=sha(run/'post/field_hashes.json'),script_sha256=sha(Path(__file__)),final=rows[-1],initial=rows[0],first_step=rows[1],maxima={key:max(r[key] for r in rows) for key in rows[0] if key not in ['frame','time_s']},dominant_final_volume_neumann_mode=int(dominant),final_volume_mode_energy_fractions=(energies[-1]/total[-1]).tolist(),max_top_noise_sigma=float(max(r['moleFraction_top_sigma'] for r in rows)),frozen_endpoint_screen_passed=bool(max(r['moleFraction_max_endpoint'] for r in rows)<=1e-12),npz_sha256=sha(out/'noise_fields.npz'))
    (out/'noise_summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(dict(frames=len(rows),dominant_mode=int(dominant),max_endpoint=summary['maxima']['moleFraction_max_endpoint'])))


if __name__=='__main__':main(*sys.argv[1:])
