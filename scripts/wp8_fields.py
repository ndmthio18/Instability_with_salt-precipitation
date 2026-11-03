"""Read-only hash-verified field reduction; no solver calls."""
import csv, hashlib, json, sys, time
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
from single_wavelength_post import widths, params
from wp6_4t_metrics import neumann_spectrum

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def reduce_fields(run, out):
    run,out=Path(run),Path(out);out.mkdir(parents=True,exist_ok=True);start=time.time()
    records=json.loads((run/'post/field_hashes.json').read_text())
    with (run/'surface_events.csv').open() as f:times=[float(r['time_s']) for r in csv.DictReader(f)]
    assert len(times)==len(records)
    phi=float(params(run/'params.input')['Problem.Porosity'])
    arrays={k:[] for k in ['x_mean','sw_mean','salt_mean','covariance','x_sigma','energy','top_coeff']};summary=[];snap={};grid=None
    for i,rec in enumerate(records):
        p=run/Path(rec['path']).name;raw=p.read_bytes();assert hashlib.sha256(raw).hexdigest()==rec['sha256']
        a={v.attrib['Name']:np.fromstring(v.text or '',sep=' ') for v in ET.fromstring(raw).iter('DataArray') if v.attrib.get('Name') in ['Coordinates','moleFraction','sw']}
        xyz=a['Coordinates'].reshape(-1,3)
        if grid is None:
            grid=xyz;xs=np.unique(xyz[:,0]);zs=np.unique(xyz[:,1]);wx=widths(xs);wz=widths(zs);W=wx.sum();order=np.argsort(np.searchsorted(zs,xyz[:,1])*len(xs)+np.searchsorted(xs,xyz[:,0]));k=np.arange(1,len(xs))*np.pi/W
        else:assert np.array_equal(xyz,grid)
        x=a['moleFraction'][order].reshape(len(zs),len(xs));sw=a['sw'][order].reshape(x.shape)
        mean=x@wx/W;sm=sw@wx/W;d=x-mean[:,None];cov=((sw-sm[:,None])*d)@wx/W;salt=phi*(sw*x)@wx/W
        coeff,en=neumann_spectrum(xs,zs,d);total=en.sum();sigma=np.sqrt(d*d@wx/W)
        assert np.max(abs(salt-phi*(sm*mean+cov)))<1e-16
        assert abs(total-np.sum(sigma*sigma*wz))<max(1e-25,total*1e-10)
        for key,val in zip(arrays,[mean,sm,salt,cov,sigma,en,coeff[-1]]):arrays[key].append(val)
        summary.append(dict(frame=i,time_s=times[i],nonzero_energy=total,k_peak_rad_m=k[np.argmax(en)] if total>1e-26 else None,k_E_rad_m=float(k@en/total) if total>1e-26 else None,canonical_rms=float(np.sqrt(en[1]/wz.sum())) if len(en)>1 else None,max_covariance=float(np.max(abs(cov))),integrated_covariance_salt_m=float(phi*cov@wz),mean_top=mean[-1],max_top=x[-1].max()))
        if i%1200==0 or i==len(records)-1:snap['anomaly_'+str(i)]=d;snap['x_'+str(i)]=x
    np.savez_compressed(out/'profiles_spectra.npz',time_s=times,z_m=zs,xs=xs,k_rad_m=k,phi=phi,**arrays,**snap)
    with (out/'timeseries.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(summary[0]));w.writeheader();w.writerows(summary)
    audit=dict(run=str(run),frames=len(records),all_field_hashes_verified=True,field_manifest_sha256=sha(run/'post/field_hashes.json'),input_sha256=sha(run/'params.input'),script_sha256=sha(__file__),wall_s=time.time()-start,files={n:sha(out/n) for n in ['profiles_spectra.npz','timeseries.csv']})
    (out/'audit.json').write_text(json.dumps(audit,indent=2));print(json.dumps(audit))
if __name__=='__main__':reduce_fields(*sys.argv[1:])
