"""Hash-verified accepted M1 fields: modal structure and initial-shape control.

Run remotely after the solver/postprocessor finishes. No solver is invoked.
"""
import csv,hashlib,json,sys
from pathlib import Path
import numpy as np
import xml.etree.ElementTree as ET
from single_wavelength_post import widths
from wp6_4t_metrics import energy_depth,neumann_spectrum

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(run,out):
    run=Path(run);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    hashes=json.loads((run/'post/field_hashes.json').read_text())
    with (run/'surface_events.csv').open() as f:ev=list(csv.DictReader(f))
    times=np.array([float(r['time_s']) for r in ev]);assert len(times)==len(hashes)
    targets=[0,136423.0637723484,272846.1275446968,409269.1913170452,545692.2550893936]
    selected={int(np.argmin(abs(times-t))):t for t in targets if t<=times[-1]+1e-7};selected.setdefault(len(times)-1,float(times[-1]))
    profile=np.loadtxt(run/'eigenprofile.dat');rows=[];grid=None;initial=None;snap={};snapmeta=[];depthprofiles=[]
    for i,record in enumerate(hashes):
        path=Path(record['path']);assert path.resolve().is_relative_to(run.resolve());raw=path.read_bytes();assert hashlib.sha256(raw).hexdigest()==record['sha256']
        a={v.attrib['Name']:np.fromstring(v.text or '',sep=' ') for v in ET.fromstring(raw).iter('DataArray') if v.attrib.get('Name') in ['Coordinates','moleFraction']}
        xyz=a['Coordinates'].reshape(-1,3)
        if grid is None:
            grid=xyz;xs=np.unique(xyz[:,0]);zs=np.unique(xyz[:,1]);nx=len(xs);nz=len(zs);wx=widths(xs);wz=widths(zs);W=wx.sum();H=wz.sum();vol=wz[:,None]*wx[None,:]/W
            order=np.argsort(np.searchsorted(zs,xyz[:,1])*nx+np.searchsorted(xs,xyz[:,0]));assert len(xyz)==nx*nz
            cm=np.cos(np.pi*np.arange(1,nx)[:,None]*(xs-xs[0])/W);basis=np.interp(zs,profile[:,0],profile[:,1])[:,None]*cm[1][None,:];bnorm=np.sum(basis*basis*vol)
        else:assert np.array_equal(grid,xyz)
        x=a['moleFraction'][order].reshape(nz,nx);mean=x@wx/W;d=x-mean[:,None];sd=np.sqrt(d*d@wx/W);coeff,en=neumann_spectrum(xs,zs,d);total=np.sum(sd**2*wz)
        if initial is None:initial=d.copy();initial_energy=total
        proj=np.sum(d*basis*vol)/bnorm
        initial_scale=np.sum(d*initial*vol)/initial_energy if initial_energy>0 else None
        shape_res=float(np.sqrt(np.sum((d-initial_scale*initial)**2*vol)/total)) if initial_scale is not None and total>0 else None
        fraction=lambda value:float(value/total) if total>0 else None
        row=dict(frame=i,time_s=float(times[i]),top_sigma=float(sd[-1]),top_relative_sigma=float(sd[-1]/mean[-1]),volume_anomaly_rms=float(np.sqrt(total/H)),energy_relative_initial=float(total/initial_energy) if initial_energy>0 else None,
            fundamental_volume_energy_fraction=fraction(en[1]),second_harmonic_volume_energy_fraction=fraction(en[3]),third_harmonic_volume_energy_fraction=fraction(en[5]),
            even_higher_harmonic_fraction=fraction(sum(en[(np.arange(1,nx)%2==0)&(np.arange(1,nx)>=4)])),odd_neumann_energy_fraction=fraction(sum(en[np.arange(1,nx)%2==1])),
            dominant_neumann_mode=int(np.argmax(en)+1) if total>0 else None,dominant_cells_per_wavelength=float(2*(nx-1)/(np.argmax(en)+1)) if total>0 else None,energy_fraction_at_four_or_fewer_cells_per_wavelength=fraction(sum(en[np.arange(1,nx)>=12])),
            top_fundamental_amplitude=float(abs(coeff[-1,1])),eigenmode_projection=float(proj),initial_shape_relative_residual=shape_res,max_relative_anomaly=float(np.max(abs(d)/mean[:,None])),
            depth90_mm=energy_depth(zs,sd,.9),depth99_mm=energy_depth(zs,sd,.99),
            max_mirror_difference=float(np.max(abs(x-x[:,::-1]))),max_endpoint_difference=float(np.max(abs(x[:,0]-x[:,-1]))))
        rows.append(row)
        if i in selected:
            snap['x_'+str(i)]=x;snap['anomaly_'+str(i)]=d;snap['mean_'+str(i)]=mean;snap['sigma_'+str(i)]=sd;snap['mode_energy_'+str(i)]=en
            snapmeta.append(dict(frame=i,actual_time_s=float(times[i]),nominal_target_s=selected[i],source_path=str(path),sha256=record['sha256']))
            for z,mu,sigma in zip(zs,mean,sd):depthprofiles.append(dict(frame=i,time_s=float(times[i]),z_m=float(z),depth_mm=float((zs[-1]-z)*1000),mean=float(mu),sigma=float(sigma)))
    for name,rr in [('spatial_metrics.csv',rows),('selected_vertical_profiles.csv',depthprofiles)]:
        with (out/name).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rr[0]));w.writeheader();w.writerows(rr)
    np.savez_compressed(out/'selected_fields.npz',xs=xs,zs=zs,**snap)
    result=dict(frames=len(rows),selected=snapmeta,all_field_hashes_verified=True,script_sha256=sha(Path(__file__)),metrics_helper_sha256=sha(Path(__file__).with_name('wp6_4t_metrics.py')),field_manifest_sha256=sha(run/'post/field_hashes.json'),npz_sha256=sha(out/'selected_fields.npz'),initial=rows[0],final=rows[-1],convention='Neumann m=2 fundamental; initial shape residual also includes linear profile evolution, not proof of nonlinearity')
    (out/'spatial_audit.json').write_text(json.dumps(result,indent=2));print(json.dumps(dict(frames=len(rows),final=rows[-1])))
if __name__=='__main__':main(*sys.argv[1:])
