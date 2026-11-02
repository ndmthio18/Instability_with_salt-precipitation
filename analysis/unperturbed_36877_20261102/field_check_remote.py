
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np,json,hashlib
r=Path('/temp/zhixin/paper2_transport_20261102/unperturbed_36877_20261102_r01')
finished=json.loads((r/'finished.json').read_text());records=[]
for ds in ET.parse(r/'unperturbed_36877.pvd').getroot().iter('DataSet'):
 p=r/ds.attrib['file'];assert hashlib.sha256(p.read_bytes()).hexdigest()==finished['files'][p.name]['sha256']
 arrays={a.attrib.get('Name'):np.fromstring(a.text or '',sep=' ') for a in ET.parse(p).getroot().iter('DataArray') if a.attrib.get('format')=='ascii' and a.attrib.get('Name') in ['Coordinates','moleFraction','sw','pw']}
 xyz=arrays['Coordinates'].reshape(-1,3);z,ix=np.unique(xyz[:,1],return_inverse=True);n=np.bincount(ix)
 row={'time_s':float(ds.attrib['timestep']),'file':p.name,'sha256':finished['files'][p.name]['sha256']}
 for name in ['moleFraction','sw','pw']:
  v=arrays[name];mean=np.bincount(ix,weights=v)/n
  lo=np.full(len(z),np.inf);hi=np.full(len(z),-np.inf);np.minimum.at(lo,ix,v);np.maximum.at(hi,ix,v)
  centered=v-lo[ix];mean_centered=np.bincount(ix,weights=centered)/n
  std=np.sqrt(np.bincount(ix,weights=(centered-mean_centered[ix])**2)/n)
  row[name+'_max_lateral_span']=float(np.max(hi-lo));row[name+'_max_lateral_std']=float(np.max(std))
 records.append(row)
print(json.dumps({'frames_checked':len(records),'records':records,'max_lateral_x_span':max(q['moleFraction_max_lateral_span'] for q in records),'max_lateral_x_std':max(q['moleFraction_max_lateral_std'] for q in records)},indent=2))
