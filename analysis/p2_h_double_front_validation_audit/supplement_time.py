from make_figures import *

def main():
 fig,ax=layout();handles=[];labels=[]
 pairs=[(f'fast_p0_ld_{lab}_r1e-09',f'fast_p0_ld_{lab}_r1e-10',lab+' (128)') for lab in ['u06','u09','x008']]
 pairs.append(('fast_grid_u06_n512','fast_production_u06_n512_tight','u = 0.6 (512)'))
 for i,(ca,cb,label) in enumerate(pairs):
  if not all((ROOT/'raw'/f'{c}_status.txt').exists() for c in [ca,cb]):continue
  a,b=map(read_raw,[ca,cb]);t=np.intersect1d(a.t,b.t);a=a.set_index('t').loc[t];b=b.set_index('t').loc[t]
  for axis,key in zip(ax.flat,['vt','vb','vnet','hc']):
   error=np.abs(a[key]-b[key])/(b.vref if key!='hc' else b.hc)*100;mask=error>0
   line,=axis.plot(t[mask]/3600,error[mask],color=colors[i]);axis.set_yscale('log');axis.set_xlabel('Time (h)')
  handles.append(line);labels.append(label)
 for a,label in zip(ax.flat,['Top time error / reference (%)','Bottom time error / reference (%)','Net time error / reference (%)','Thickness time error (%)']):a.set_ylabel(label)
 save(fig,'figure8b_temporal_convergence',handles,labels,2)
if __name__=='__main__':main()
