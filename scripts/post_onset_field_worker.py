"""Remote stdin-only sparse VTU inspection; all original files are read-only."""
import base64,hashlib,io,json,os
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np

def sha(b):return hashlib.sha256(b).hexdigest()
def inspect(spec):
 source=Path(spec['path']);p=spec['parameters'];checks={};hashes={}
 for name,expected in spec['small_file_hashes'].items():
  b=(source/name).read_bytes();hashes[name]=sha(b);checks[name]=hashes[name]==expected
 assert all(checks.values()),'source changed'
 sigpath=next(source.glob('*_sigmatop.csv'));b=sigpath.read_bytes();assert sha(b)==spec['sigmatop_sha256'];sig=np.loadtxt(io.BytesIO(b),ndmin=2)
 datasets=list(ET.parse(next(source.glob('*.pvd'))).getroot().iter('DataSet'))
 profile=np.loadtxt(next(source.glob('*eigenprofile.dat')))
 scalar=lambda k:float(p['Problem.'+k]);mw,ms,r0,lin,quad,n,alpha,swr=[scalar(k) for k in ['MolarMassH2O','MolarMassNaCl','PureWaterDensity','DensityLinearCoefficient','DensityQuadraticCoefficient','PcN','PcAlpha','Swr']]
 fields={};records=[]
 for fi in sorted(set([0,len(datasets)-1])):
  ds=datasets[fi];f=(source/ds.attrib['file']).resolve();assert f.is_relative_to(source.resolve());st=f.stat();b=f.read_bytes();assert (st.st_size,st.st_mtime_ns)==(f.stat().st_size,f.stat().st_mtime_ns)
  tree=ET.fromstring(b);arrays={a.attrib.get('Name'):np.fromstring(a.text or '',sep=' ') for a in tree.iter('DataArray') if a.attrib.get('Name') in ['Coordinates','pw','sw','moleFraction','rhoW'] and a.attrib.get('format')=='ascii'}
  xyz=arrays['Coordinates'].reshape(-1,3);xs=np.unique(xyz[:,0]);zs=np.unique(xyz[:,1]);nx=len(xs);nz=len(zs);ix=np.searchsorted(xs,xyz[:,0]);iz=np.searchsorted(zs,xyz[:,1]);order=np.argsort(iz*nx+ix);assert np.array_equal((iz*nx+ix)[order],np.arange(nx*nz))
  x,sw,pw,rho=[arrays[k][order].reshape(nz,nx) for k in ['moleFraction','sw','pw','rhoW']];assert all(np.isfinite(v).all() for v in [x,sw,pw,rho])
  w=x*ms/((1-x)*mw+x*ms);rhoerr=float(np.max(abs(rho-r0-1000*w*(lin+quad*w))));swcalc=swr+(1-swr)/(1+(alpha*np.maximum(1e5-pw,0))**n)**(1-1/n);swerr=float(np.max(abs(sw-swcalc)))
  wx=np.diff(np.r_[xs[0],.5*(xs[1:]+xs[:-1]),xs[-1]]);mean=x@wx/(xs[-1]-xs[0]);sd=np.sqrt(((x-mean[:,None])**2)@wx/(xs[-1]-xs[0]))
  expected=scalar('SaltMoleFraction')+scalar('Amplitude')*np.interp(zs,profile[:,0],profile[:,1])[:,None]*np.cos(2*np.pi*scalar('WaveNumber')*xs[None,:]/(xs[-1]-xs[0]))
  j=0 if fi==0 else len(sig)-1;terr=abs(float(ds.attrib['timestep'])-sig[j,0]);same_round=float(format(sig[j,0],'.6g'))==float(ds.attrib['timestep']);serr=abs(np.std(x[-1])-sig[j,2]);merr=abs(np.mean(x[-1])-sig[j,1])
  rec=dict(frame=fi,time_s=float(sig[j,0]),pvd_time_s=float(ds.attrib['timestep']),time_roundtrip=same_round,raw_path=str(f),raw_sha256=sha(b),bytes=len(b),rho_error=rhoerr,sw_error=swerr,sigma_error=serr,mean_error=merr,initial_field_error=float(np.max(abs(x-expected))) if fi==0 else None,x_max=float(x.max()),x_min=float(x.min()),top_sigma=float(np.std(x[-1])),top_mean=float(np.mean(x[-1])),max_layer_relative=float(np.max(sd/mean)),consistent=bool(rhoerr<1e-8 and swerr<1e-10 and serr<1e-10 and merr<1e-10 and same_round))
  records.append(rec);fields['x_'+str(j>0)]=x;fields['sw_'+str(j>0)]=sw;fields['mean_'+str(j>0)]=mean;fields['sd_'+str(j>0)]=sd;fields['xs']=xs;fields['z_m']=zs
 bio=io.BytesIO();np.savez_compressed(bio,**fields)
 # Frozen manifests contain neutral curves, not finite growth eigenvalues.
 manifest=json.loads((source/'case_manifest.json').read_text()) if (source/'case_manifest.json').exists() else {};lsa=manifest.get('lsa',{});runmeta=json.loads((source/'run_meta.json').read_text()) if (source/'run_meta.json').exists() else {}
 return dict(run_key=spec['run_key'],records=records,small_hashes=hashes,fields_base64=base64.b64encode(bio.getvalue()).decode(),lsa=lsa,run_meta=runmeta,log_tail=spec.get('log_tail'))
