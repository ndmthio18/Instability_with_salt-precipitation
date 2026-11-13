from make_figures import *

def main():
 fig,ax=layout(3,2,175)
 for j,lab in enumerate(['u06','x008']):
  for k,n in enumerate([128,256,512]):
   h=read(f'setup_early_{lab}_n{n}_b1e-09');tight=read(f'setup_early_{lab}_n{n}_tight')
   h=h[h.t<=2];tight=tight[tight.t<=2]
   for a,key in zip(ax[:,j],['vt','vb','vnet']):
    a.plot(h.t,h[key]*1e9,color=colors[k]);a.plot(tight.t,tight[key]*1e9,color=colors[k],ls='--');a.set_xlim(0,2);a.set_xlabel('Time (s)')
  for a,label in zip(ax[:,j],['Top speed (nm s⁻¹)','Bottom speed (nm s⁻¹)','Net speed (nm s⁻¹)']):a.set_ylabel(label)
  ax[0,j].set_title('u = 0.6' if lab=='u06' else 'x008',fontsize=9)
 save(fig,'figure2b_first_seconds_velocities',[Line2D([],[],color=colors[k]) for k in range(3)]+[Line2D([],[],color='k',ls='--')],['128 cells','256 cells','512 cells','Refined time'],4)
 fig,ax=layout(2,2,130)
 for j,lab in enumerate(['u06','x008']):
  for k,gap in enumerate([1e-10,1e-9,1e-8]):
   h=read(f'setup_early_{lab}_n512_b{gap:g}');h=h[h.t>0]
   ax[0,j].plot(h.t,h.Fplus*1e6,color=colors[k]);ax[0,j].set_xscale('log')
  for k,suffix in enumerate(['b1e-09','consistent']):
   h=read(f'setup_early_{lab}_n512_{suffix}');h=h[h.t>0]
   ax[1,j].plot(h.t,h.Fplus*1e6,color=colors[1],ls='-' if k==0 else '--');ax[1,j].set_xscale('log')
  for a in ax[:,j]:a.set_xlabel('Time (s)');a.set_ylabel('Crust diffusion (mg m⁻² s⁻¹)')
  ax[0,j].set_title('u = 0.6' if lab=='u06' else 'x008',fontsize=9)
 save(fig,'figure2c_initial_condition_sensitivity',[Line2D([],[],color=colors[k]) for k in range(3)]+[Line2D([],[],color=colors[1],ls='--')],['0.1 nm (top)','1 nm / uniform','10 nm (top)','Consistent (bottom)'],2)
if __name__=='__main__':main()
