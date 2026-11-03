"""Grid-consistent diagnostics from accepted VTU states; original ELF unchanged."""
import csv,hashlib,json,math,sys,time
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def savecsv(p,rows):
 if not rows:return
 keys=list(dict.fromkeys(k for r in rows for k in r))
 with p.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
def params(p):
 sec='';d={}
 for line in p.read_text().splitlines():
  line=line.split('#')[0].strip()
  if line.startswith('['):sec=line.strip('[]')
  elif '=' in line:k,v=line.split('=',1);d[sec+'.'+k.strip()]=v.strip()
 return d
def widths(a):return np.diff(np.r_[a[0],(a[1:]+a[:-1])/2,a[-1]])
def post(run,out):
 run=Path(run);out=Path(out);out.mkdir(exist_ok=True,parents=True);p=params(run/'params.input');scalar=lambda k:float(p['Problem.'+k]);profile=np.loadtxt(run/'eigenprofile.dat');sigma=np.loadtxt(next(run.glob('*sigmatop.csv')),ndmin=2)
 ds=list(ET.parse(next(run.glob('*.pvd'))).getroot().iter('DataSet'));rows=[];hashes=[];means=[];sds=[];times=[];representative={};grid=None;initial_error=None
 for fi,d in enumerate(ds):
  path=(run/d.attrib['file']).resolve();assert path.is_relative_to(run.resolve());raw=path.read_bytes();tree=ET.fromstring(raw);a={x.attrib['Name']:np.fromstring(x.text or '',sep=' ') for x in tree.iter('DataArray') if x.attrib.get('Name') in ['Coordinates','moleFraction','sw','rhoW','pw'] and x.attrib.get('format')=='ascii'};xyz=a['Coordinates'].reshape(-1,3)
  if grid is None:
   grid=xyz;xs=np.unique(xyz[:,0]);zs=np.unique(xyz[:,1]);nx=len(xs);nz=len(zs);W=xs[-1]-xs[0];H=zs[-1]-zs[0];wx=widths(xs);wz=widths(zs);order=np.argsort(np.searchsorted(zs,xyz[:,1])*nx+np.searchsorted(xs,xyz[:,0]));assert len(xyz)==nx*nz
   phase=2*np.pi*scalar('WaveNumber')*(xs-xs[0])/W;cos=np.cos(phase);sin=np.sin(phase);A=np.interp(zs,profile[:,0],profile[:,1]);basis=A[:,None]*cos[None,:];vol=wz[:,None]*wx[None,:]/W;denom=np.sum(basis*basis*vol)
   if nx<100:cm=np.cos(np.pi*np.arange(1,nx)[:,None]*(xs-xs[0])/W)
  else:assert np.array_equal(grid,xyz)
  x,sw,rho,pw=[a[k][order].reshape(nz,nx) for k in ['moleFraction','sw','rhoW','pw']];mean=x@wx/W;anom=x-mean[:,None];sd=np.sqrt(anom**2@wx/W);fourier=2*abs(np.sum(wx*anom[-1]*np.exp(-1j*phase)))/W;projection=np.sum(anom*basis*vol)/denom
  j=int(np.argmin(abs(sigma[:,0]-float(d.attrib['timestep']))));t=float(sigma[j,0]);assert float(format(t,'.6g'))==float(d.attrib['timestep']);assert abs(np.mean(x[-1])-sigma[j,1])<1e-10 and abs(np.std(x[-1])-sigma[j,2])<1e-10
  if fi==0:
   expected=scalar('SaltMoleFraction')+(scalar('Amplitude')*basis if p['Problem.PerturbationType']=='Eigenprofile' else 0);initial_error=float(np.max(abs(x-expected)))
  mass=x*scalar('MolarMassNaCl')/((1-x)*scalar('MolarMassH2O')+x*scalar('MolarMassNaCl'));rho_calc=scalar('PureWaterDensity')+1000*mass*(scalar('DensityLinearCoefficient')+scalar('DensityQuadraticCoefficient')*mass);se=(1+(scalar('PcAlpha')*np.maximum(1e5-pw,0))**scalar('PcN'))**(-(1-1/scalar('PcN')));swcalc=scalar('Swr')+(1-scalar('Swr'))*se
  odd=None
  if nx<100:
   coef=2*(anom*wx[None,:])@cm.T/W;energy=.5*coef**2;energy[:,-1]*=.5;total=np.sum(sd**2*wz);odd=float(np.sum(energy[:,np.arange(1,nx)%2==1]*wz[:,None])/total) if total>1e-26 else None
  row=dict(frame=fi,step=j,time_s=t,pvd_time_s=float(d.attrib['timestep']),weighted_top_mean=mean[-1],weighted_top_sigma=sd[-1],weighted_top_relative=sd[-1]/mean[-1],weighted_fourier=fourier,weighted_projection=projection,whole_salt_m=np.sum(scalar('Porosity')*sw*x*vol),whole_water_m=np.sum(scalar('Porosity')*sw*vol),top_max=x[-1].max(),domain_max=x.max(),max_mirror_error=np.max(abs(x-x[:,::-1])),max_endpoint_difference=np.max(abs(x[:,0]-x[:,-1])),odd_neumann_energy_fraction=odd,initial_field_error=initial_error,rho_error=np.max(abs(rho-rho_calc)),sw_error=np.max(abs(sw-swcalc)))
  rows.append({k:float(v) if isinstance(v,np.generic) else v for k,v in row.items()});hashes.append(dict(frame=fi,path=str(path),sha256=hashlib.sha256(raw).hexdigest(),bytes=len(raw)));times.append(t);means.append(mean);sds.append(sd)
  if fi in [0,len(ds)//2,len(ds)-1]:representative['x_'+str(fi)]=x;representative['sw_'+str(fi)]=sw
 savecsv(out/'weighted_diagnostics.csv',rows);np.savez_compressed(out/'profiles.npz',time_s=times,z_m=zs,x_mean=means,x_sigma=sds);np.savez_compressed(out/'representative_fields.npz',xs=xs,z_m=zs,**representative);(out/'field_hashes.json').write_text(json.dumps(hashes,indent=2))
 meta=dict(frames=len(ds),vertices=nx*nz,cells=[nx-1,nz-1],width_m=W,height_m=H,physical_wave_number=2*np.pi*scalar('WaveNumber')/W,initial_field_error=initial_error,all_rho_pass=all(r['rho_error']<1e-8 for r in rows),all_sw_pass=all(r['sw_error']<1e-10 for r in rows),profile_sha256=sha(run/'eigenprofile.dat'),input_sha256=sha(run/'params.input'),postprocessor_sha256=sha(Path(__file__)),npz_sha256={n:sha(out/n) for n in ['profiles.npz','representative_fields.npz']});(out/'field_audit.json').write_text(json.dumps(meta,indent=2));return meta
if __name__=='__main__':print(json.dumps(post(sys.argv[1],sys.argv[2])))
