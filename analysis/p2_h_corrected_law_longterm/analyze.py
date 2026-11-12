from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.linalg import solve_banded
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results';OUT.mkdir(exist_ok=True)
FIG=ROOT/'figures';FIG.mkdir(exist_ok=True)
ceq=.264940869768;rho=1202.168577;cr=2165.;birth=.43;kr=.001;ds=1.3e-9;grain=25e-6

def law(h,j,B,eps,beta=2,avmult=1):
    md=eps**1.5*ds/(1-ceq);av=36*(1-eps)/grain*avmult
    ell=np.sqrt(md/(av*kr));X=beta*(1-birth)*kr*ell/md;H=h/ell
    sech=2*np.exp(-H)/(1+np.exp(-2*H));u=B/(ceq*j)
    v=ceq*j/(1-ceq)/(cr*(1-birth))*X*(1-u*sech)/(X+np.tanh(H))
    pe=j*ell/(rho*eps**1.5*ds)
    return v,H*H,X*X,pe,ell/h

def frozen_bvp(g,j,B,beta=2,avmult=1):
    """Conservative linear variable-coefficient RD BVP; no ALE top inputs."""
    eps=g.eps.to_numpy();h=g.h.iloc[0];w=g.width.to_numpy()*h
    z=g.xi.to_numpy()*h;faces=np.r_[0,np.cumsum(w)]
    a=rho*eps**1.5*ds/(1-ceq);r=rho*36*(1-eps)/grain*kr*avmult
    conduct=1/((faces[1:-1]-z[:-1])/a[:-1]+(z[1:]-faces[1:-1])/a[1:])
    kt=rho*beta*(1-birth)*kr;gt=a[-1]/(h-z[-1]);ge=kt*gt/(kt+gt)
    diagonal=r*w;diagonal[:-1]+=conduct;diagonal[1:]+=conduct;diagonal[-1]+=ge
    rhs=np.zeros(len(w));rhs[0]-=B/(1-ceq);rhs[-1]+=ceq*j/(1-ceq)*gt/(gt+kt)
    band=np.zeros((3,len(w)));band[1]=diagonal;band[0,1:]=-conduct;band[2,:-1]=-conduct
    c=solve_banded((1,1),band,rhs);ct=(ceq*j/(1-ceq)+gt*c[-1])/(gt+kt)
    return kt*ct/(cr*(1-birth))

def main():
    cases=pd.read_csv(ROOT/'cases.csv');hist=[];profiles=[];summary=[];frozen=[]
    if (ROOT/'cases_p1.csv').exists():cases=pd.concat([cases,pd.read_csv(ROOT/'cases_p1.csv')],ignore_index=True)
    if (ROOT/'cases_events.csv').exists():cases=pd.concat([cases,pd.read_csv(ROOT/'cases_events.csv')],ignore_index=True)
    if (ROOT/'cases_refine.csv').exists():cases=pd.concat([cases,pd.read_csv(ROOT/'cases_refine.csv')],ignore_index=True)
    pending=cases[~cases.id.map(lambda name:(ROOT/'raw'/f'{name}_status.txt').exists())]
    pending.to_csv(OUT/'pending_cases.csv',index=False)
    cases=cases[~cases.id.isin(pending.id)]
    for row in cases.itertuples(index=False):
        d=pd.read_csv(ROOT/'raw'/f'{row.id}_history.csv');d['id']=row.id;d['u_input']=row.u;d['n']=row.n
        d['beta']=row.beta;d['avmult']=row.av;d['jfactor']=row.jfactor
        for tag,eps in [('A',birth),('Bmean',d.epsmean)]:
            vals=law(d.h,d.j,d.B,eps,row.beta,row.av)
            for name,value in zip(['v','Da','Dastar','Pe','ell_over_h'],vals):d[f'{tag}_{name}']=value
        for tag in ['A','Bmean']:d[tag+'_error']=abs(d[tag+'_v']-d.v)/abs(d.v).clip(lower=1e-30)
        status=dict(line.split('=',1) for line in (ROOT/'raw'/f'{row.id}_status.txt').read_text().splitlines())
        noninitial=d[d.t>=1/row.jfactor]
        summary.append(dict(id=row.id,u=row.u,n=row.n,power=row.power,rtol=row.rtol,status=status['status'],
            previous_s=float(status['previous_s']),end_s=float(status['last_s']),beta=row.beta,avmult=row.av,jfactor=row.jfactor,final_h=d.h.iloc[-1],
            final_v=d.v.iloc[-1],final_A_error=d.A_error.iloc[-1],max_A_error=noninitial.A_error.max(),
            final_epsmin=d.epsmin.iloc[-1],final_epsmax=d.epsmax.iloc[-1],
            salt_res=d.salt_res.abs().max(),water_res=d.water_res.abs().max(),diss_res=d.diss_res.abs().max(),local_res=d.local_res.abs().max()))
        hist.append(d)
        if (row.n==64 and row.power==1.5 and row.rtol==1e-8) or (row.n>=128 and row.power==2):
            pp=pd.read_csv(ROOT/'raw'/f'{row.id}_profiles.csv');pp['id']=row.id
            for time,g in pp.groupby('t'):
                state=d.loc[np.isclose(d.t,time)].iloc[0];ell=law(state.h,state.j,state.B,birth,row.beta,row.av)[4]*state.h
                left=(g.xi-g.width/2)*state.h;right=(g.xi+g.width/2)*state.h
                weights=np.maximum(0,right-np.maximum(left,state.h-3*ell))
                epsupper=np.sum(g.eps*weights)/np.sum(weights)
                upper=law(state.h,state.j,state.B,epsupper,row.beta,row.av)[0]
                vb=frozen_bvp(g,state.j,state.B,row.beta,row.av)
                scale=state.j*ceq/(1-ceq)/ell
                terms={f'{term}_norm':float(np.sum(abs(g[term])*g.width)/scale) for term in ['storage','advection','diffusion','mesh','reaction']}
                # Pore-volume changes can dominate Q storage even with nearly
                # stationary concentration; reconstruct the same ALE E equation.
                faces=np.r_[0,np.cumsum(g.width.to_numpy())]
                epface=np.r_[g.eps.to_numpy()[0],g.eps.to_numpy()[1:],birth]
                H=-faces*state.v*epface
                de=(H[:-1]-H[1:])/(g.width.to_numpy()*state.h)-g.R.to_numpy()/cr
                concentration_storage=g.storage.to_numpy()-rho*g.C.to_numpy()*de
                terms['concentration_storage_norm']=float(np.sum(abs(concentration_storage)*g.width)/scale)
                terms['concentration_storage_to_reaction']=float(np.sum(abs(concentration_storage)*g.width)/max(np.sum(abs(g.R)*g.width),1e-30))
                frozen.append(dict(id=row.id,t=time,A=state.A_v,Bmean=state.Bmean_v,Bupper=upper,Cbvp=vb,ALE=state.v,epsupper=epsupper,**terms))
            profiles.append(pp)
    history=pd.concat(hist);summ=pd.DataFrame(summary);fr=pd.DataFrame(frozen)
    history[(history.t<=3600)|(history.t%600==0)].to_csv(OUT/'history.csv.gz',index=False,compression={'method':'gzip','mtime':0});summ.to_csv(OUT/'summary.csv',index=False);fr.to_csv(OUT/'frozen_bvp_terms.csv',index=False)
    pd.concat(profiles).to_csv(OUT/'profiles.csv.gz',index=False,compression={'method':'gzip','mtime':0})
    convergence=[]
    for (u,beta,av,jfactor),sub in cases.groupby(['u','beta','av','jfactor']):
        refid=sub.sort_values(['n','power']).id.iloc[-1];ref=history[history.id==refid]
        for row in sub.itertuples():
            d=history[(history.id==row.id)&(history.t>=1/jfactor)&(history.t<=ref.t.max())]
            if d.empty:continue
            rv=np.interp(d.t,ref.t,ref.v)
            extra={col+'_max_abs_diff':float(np.max(abs(d[col]-np.interp(d.t,ref.t,ref[col])))) for col in ['ct','cb','Rplus','Dminus','solid','dissolved']}
            convergence.append(dict(id=row.id,u=u,n=row.n,power=row.power,rtol=row.rtol,reference_id=refid,
                velocity_max_relative=float(np.max(abs(d.v-rv)/abs(rv))),common_end_s=d.t.max(),**extra))
    conv=pd.DataFrame(convergence);conv.to_csv(OUT/'convergence.csv',index=False)
    temporal=[]
    for u in cases.u.unique():
        pair=cases[(cases.u==u)&(cases.n==64)&(cases.power==1.5)&(~cases.id.str.startswith('p1'))]
        if len(pair)!=2:continue
        strict=pair.sort_values('rtol').iloc[0];normal=pair.sort_values('rtol').iloc[-1]
        a=history[history.id==normal.id];b=history[history.id==strict.id]
        a=a[(a.t>=1)&(a.t<=b.t.max())];bv=np.interp(a.t,b.t,b.v)
        temporal.append(dict(u=u,normal_id=normal.id,strict_id=strict.id,common_end_s=a.t.max(),
            max_velocity_relative=float(np.max(abs(a.v-bv)/abs(bv))),
            final_height_relative=float(abs(a.h.iloc[-1]/np.interp(a.t.iloc[-1],b.t,b.h)-1))))
    pd.DataFrame(temporal).to_csv(OUT/'temporal_convergence.csv',index=False)
    app=[]
    for keys,group in cases.groupby(['u','beta','av','jfactor']):
        unique=sorted(group.n.unique());finest=group[group.n==unique[-1]].sort_values(['power','id']).iloc[-1]
        lowergroup=group[(group.n==unique[-2])&(group.power==finest.power)&(group.rtol==finest.rtol)] if len(unique)>1 else group.iloc[0:0]
        lower=lowergroup.sort_values('id').iloc[-1] if len(lowergroup) else finest
        check=conv[conv.id==lower.id]
        error=float(check.velocity_max_relative.iloc[0]) if len(check) else float('nan')
        state=history[history.id==finest.id].iloc[-1]
        app.append(dict(id=finest.id,u=float(state.B/(ceq*state.j)),forcing_input=keys[0],beta=keys[1],avmult=keys[2],jfactor=keys[3],
          Da=state.A_Da,Dastar=state.A_Dastar,Pe=state.A_Pe,error=state.A_error,t=state.t,
          spatial_difference=error,spatial_supported=bool(error<.01),
          classification='within_1pct' if error<.01 and state.A_error<.01 else ('law_error_over_1pct' if error<.01 else 'unresolved_grid')))
    pd.DataFrame(app).to_csv(OUT/'applicability.csv',index=False)
    figures(history,summ,fr,pd.concat(profiles),conv,cases)
    print(summ.to_string(index=False));print(conv.to_string(index=False))

def figures(h,s,f,p,conv,cases):
    plt.rcParams.update({'font.family':'serif','font.serif':['Times New Roman'],'svg.fonttype':'none','pdf.fonttype':42,'font.size':9,
      'mathtext.fontset':'custom','mathtext.rm':'Times New Roman','mathtext.it':'Times New Roman:italic','mathtext.bf':'Times New Roman:bold'})
    base=cases[(cases.n==64)&(cases.power==1.5)&(cases.rtol==1e-8)].copy()
    refined=cases[cases.id=='u0_n512_p2_r1e-08_dt30']
    if len(refined):base=pd.concat([refined,base[base.u!=0]],ignore_index=True)
    def save(fig,name):
        for suffix in ['png','pdf','svg']:fig.savefig(FIG/f'{name}.{suffix}',dpi=250)
        svg=FIG/f'{name}.svg';svg.write_text('\n'.join(line.rstrip() for line in svg.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8')
        plt.close(fig)
    fig,axs=plt.subplots(2,3,figsize=(11,6),layout='constrained')
    for ax,row in zip(axs.flat,base.itertuples()):
        d=h[h.id==row.id];b=f[f.id==row.id];ax.plot(d.t/3600,d.v*1e9,label='ALE');ax.plot(d.t/3600,d.A_v*1e9,'--',label='A reference')
        ax.plot(b.t/3600,b.Bmean*1e9,'--',label='B mean');ax.plot(b.t/3600,b.Bupper*1e9,':',label='B upper');ax.plot(b.t/3600,b.Cbvp*1e9,'-.',label='C frozen BVP');ax.set(title=f'u={row.u:g}' if row.u>=0 else 'x008',xlabel='Time (h)',ylabel='Top speed (nm/s)');ax.legend(fontsize=7,ncol=2)
    for row in base.itertuples():
        d=h[(h.id==row.id)&(h.t>0)]
        axs.flat[-1].loglog(d.t/3600,d.h*1e6,label=f'u={row.u:g}' if row.u>=0 else 'x008')
    axs.flat[-1].set(xlabel='Time (h)',ylabel='Geometric height (um)',title='Actual ALE height');axs.flat[-1].legend(fontsize=7)
    save(fig,'figure1_velocity')
    fig,axs=plt.subplots(1,3,figsize=(11,3),layout='constrained')
    for row in base.itertuples():
        d=h[(h.id==row.id)&(h.t>=1)];label='x008' if row.u<0 else f'u={row.u:g}'
        for ax,col,ylabel in zip(axs,['A_error','Bmean_Da','Bmean_Dastar'],['Model A relative error','Da (mean porosity)','Da* (mean porosity)']):ax.semilogx(d.t/3600,d[col],label=label);ax.set(xlabel='Time (h)',ylabel=ylabel)
    axs[0].legend();save(fig,'figure2_error_groups')
    fig,axs=plt.subplots(2,2,figsize=(9,6),layout='constrained')
    for k,label in enumerate(['u0_n512_p2_r1e-08_dt30','x008_n64_p1.5_r1e-08_dt30']):
        q=p[p.id==label];tt=q.t.unique();targets=[0,3600,18000,72000,180000] if k==0 else [0,30,300,900]
        selected=sorted({tt[np.argmin(abs(tt-t))] for t in targets})
        for t in selected:
            g=q[q.t==t];axs[0,k].plot((g.C-ceq)*1e6,g.xi,label=f'{t/3600:.3g} h');axs[1,k].plot(g.eps,g.xi,label=f'{t/3600:.3g} h')
        axs[0,k].set(xlabel='Concentration excess (1e-6)',ylabel='z/h',title='u=0' if k==0 else 'x008');axs[1,k].set(xlabel='Porosity',ylabel='z/h');axs[0,k].legend()
    save(fig,'figure3_profiles')
    fig,axs=plt.subplots(2,3,figsize=(11,6),layout='constrained')
    for ax,row in zip(axs.flat,base.itertuples()):
        d=h[h.id==row.id]
        for col in ['G','Rplus','Dminus']:ax.plot(d.t/3600,d[col]*1e6,label=col)
        # Secants of independently reconstructed inventories, not same-RHS residual.
        ax.plot((d.t.to_numpy()[1:]+d.t.to_numpy()[:-1])/7200,np.diff(d.dissolved)/np.diff(d.t)*1e6,':',label='dissolved storage')
        ax.set(title=f'u={row.u:g}' if row.u>=0 else 'x008',xlabel='Time (h)',ylabel='Rate (1e-6 kg/m2/s)');ax.legend(fontsize=7)
    axs.flat[-1].axis('off');save(fig,'figure4_redistribution')
    fig,axs=plt.subplots(1,2,figsize=(9,3),layout='constrained')
    for ax,name in zip(axs,['u0_n512_p2_r1e-08_dt30','x008_n64_p1.5_r1e-08_dt30']):
        d=f[(f.id==name)&(f.t>0)]
        for col in ['storage','advection','diffusion','mesh','reaction','concentration_storage']:ax.loglog(d.t/3600,d[col+'_norm'].clip(lower=1e-15),label=col.replace('_',' '))
        ax.set(xlabel='Time (h)',ylabel='Volume-weighted |term| / (Sref/ell)',title=name.split('_')[0]);ax.legend(fontsize=7)
    save(fig,'figure5_terms')
    fig,axs=plt.subplots(1,3,figsize=(12,3.5),layout='constrained');app=pd.read_csv(OUT/'applicability.csv')
    for ax,xcol,ycol in zip(axs,['Da','Pe','u'],['Dastar','Dastar','Da']):
        for resolved,marker in [(True,'o'),(False,'x')]:
            d=app[app.spatial_supported==resolved]
            sc=ax.scatter(d[xcol],d[ycol],c=d.error*100,cmap='viridis',vmin=0,vmax=max(10,app.error.max()*100),marker=marker,s=45,label='grid checked' if resolved else 'grid unresolved')
        ax.set(xlabel='u at endpoint' if xcol=='u' else xcol,ylabel='Da*' if ycol=='Dastar' else ycol,yscale='log');
        if xcol!='u':ax.set_xscale('log')
    axs[0].legend(fontsize=7);bar=fig.colorbar(sc,ax=axs,label='Model A endpoint error (%)');bar.solids.set_rasterized(False);save(fig,'figure6_p0_applicability')
    fig,ax=plt.subplots(figsize=(7,4),layout='constrained');prior=pd.read_csv(ROOT.parent/'p2_h_corrected_vs_ale_x008_20261112/results/x008_forcing.csv')
    row=base[base.u<0].iloc[0];d=h[h.id==row.id];src=prior[prior.elapsed_time_s<=d.t.max()]
    ax.plot(src.elapsed_time_s/3600,src.soil_net_nacl_kg_m2_s*1e6,label='first internal soil SCVF')
    native=OUT/'soil_native_boundary_audit.csv'
    if native.exists():
        nd=pd.read_csv(native);nd=nd[nd.elapsed_time_s<=d.t.max()]
        ax.plot(nd.elapsed_time_s/3600,nd.native_soil_j_precip_kg_m2_s*1e6,':',label='native bare-soil surface')
    ax.plot(d.t/3600,(ceq*d.j-d.B)/(1-ceq)*1e6,'--',label='analytical reconstruction');ax.plot(d.t/3600,d.Jb*1e6,label='ALE bottom')
    ax.set(xlabel='Time after onset (h)',ylabel='Salt flux (1e-6 kg/m2/s)',title='Censored at framework threshold');ax.legend();save(fig,'figure7_forcing')
    fig,axs=plt.subplots(1,3,figsize=(12,3.5),layout='constrained')
    for u in cases.u.unique():
        d=conv[(conv.u==u)&(conv.power==(2 if u==0 else 1.5))&(conv.rtol==1e-8)&(~conv.id.str.startswith('p1'))].sort_values('n');axs[0].plot(d.n,d.velocity_max_relative*100,'o-',label=f'u={u:g}' if u>=0 else 'x008')
        d=s[(s.u==u)&(s.power==1.5)&(s.rtol==1e-8)];axs[1].plot(d.n,d.end_s/3600,'o-',label=f'u={u:g}')
    temporal=pd.read_csv(OUT/'temporal_convergence.csv')
    axs[2].bar(np.arange(len(temporal)),temporal.max_velocity_relative*100)
    axs[2].set(xticks=np.arange(len(temporal)),xticklabels=['x008' if x<0 else f'{x:g}' for x in temporal.u],xlabel='u (or forcing)',ylabel='Time-setting speed difference (%)',title='rtol 1e-8 / 1e-9; max step 30 / 15 s')
    axs[0].set(xlabel='ALE cells',ylabel='Max speed difference vs finest (%)');axs[1].set(xlabel='ALE cells',ylabel='Last time (h)',yscale='log');axs[0].legend(fontsize=7);save(fig,'figure8_convergence')

if __name__=='__main__':main()
