from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
RESULT=ROOT/'results';FIG=ROOT/'figures'
CEQ=.264940869768;RHO=1202.168577;CR=2165.;EPS=.43;DS=1.3e-9;KR=.001;BETA=2.;GRAIN=25e-6
MD=EPS**1.5*DS/(1-CEQ);AV=36*(1-EPS)/GRAIN;ELL=np.sqrt(MD/(AV*KR))
X=BETA*(1-EPS)*KR*ELL/MD;PHI=(1-EPS)*(1+X)/X

def top_law(h,j,B):
    H=np.asarray(h)/ELL
    sech=2*np.exp(-H)/(1+np.exp(-2*H))
    return CEQ*j/((1-CEQ)*CR*(1-EPS))*X*(1-B/(CEQ*j)*sech)/(X+np.tanh(H))

def cumulative(y,t):
    y=np.asarray(y);t=np.asarray(t)
    return np.r_[0,np.cumsum((y[1:]+y[:-1])*.5*np.diff(t))]

def main():
    RESULT.mkdir(exist_ok=True);FIG.mkdir(exist_ok=True)
    cases=pd.read_csv(ROOT/'cases.csv')
    if (ROOT/'cases_refine.csv').exists():cases=pd.concat([cases,pd.read_csv(ROOT/'cases_refine.csv')],ignore_index=True)
    histories=[];profiles=[];summaries=[];pending=[]
    for c in cases.itertuples():
        path=ROOT/'raw'/f'{c.id}_status.txt'
        if not path.exists():pending.append(c.id);continue
        status=dict(x.split('=',1) for x in path.read_text().splitlines())
        h=pd.read_csv(ROOT/'raw'/f'{c.id}_history.csv');h['id']=c.id
        h['u_input']=c.u;h['nc']=c.nc;h['nl']=c.nl;h['kappa']=c.kappa;h['gap']=c.gap;h['rtol']=c.rtol
        h['P_corr']=(CEQ*h.j-h.B)/(1-CEQ)
        h['vtop_corr']=top_law(h.hc,h.j,h.B)
        h['vtop_contact_diagnostic']=top_law(h.hc,h.j,h.Fplus)
        h['vnet_corr']=h.P_corr/(CR*PHI)
        h['h_corr']=h.hc.iloc[0]+cumulative(h.vnet_corr,h.t)
        h['h_inventory']=h.hc.iloc[0]+(h.solid-h.solid.iloc[0])/(CR*PHI)
        h['phi_mean']=h.solid/(CR*h.hc);h['dissolved']=h.crust_salt+h.liquid_salt
        h['solid_gain_rate']=h.G-h.Db+h.Rplus-h.Dminus
        h['dissolved_storage_rhs']=h.Jsoil-h.solid_gain_rate
        h['v_mass_equivalent']=h.solid_gain_rate/(CR*PHI)
        h['v_reference']=CEQ*h.j/((1-CEQ)*CR*PHI)
        h['top_error']=abs(h.vt-h.vtop_corr)/h.vt.abs().clip(lower=1e-30)
        h['net_error_scale']=abs(h.vnet-h.vnet_corr)/h.v_reference
        h['net_error_relative']=abs(h.vnet-h.vnet_corr)/h.vnet.abs().clip(lower=1e-30)
        # Exact rate decomposition, a diagnostic identity rather than a budget test.
        phi_dot=(h.solid_gain_rate/CR-h.phi_mean*h.vnet)/h.hc
        h['geometry_structure_term']=-h.hc*phi_dot/h.phi_mean
        h['mean_density_term']=h.solid_gain_rate/CR*(1/h.phi_mean-1/PHI)
        h['storage_term']=-h.dissolved_storage_rhs/(CR*PHI)
        h['supply_term']=(h.Jsoil-h.P_corr)/(CR*PHI)
        histories.append(h)
        p=pd.read_csv(ROOT/'raw'/f'{c.id}_profiles.csv');p['id']=c.id;profiles.append(p)
        last=h.iloc[-1];late=h[h.t>=.8*h.t.max()]
        row=dict(id=c.id,u=c.u,nc=c.nc,nl=c.nl,kappa=c.kappa,gap=c.gap,rtol=c.rtol,status=status['status'],
          end_s=float(status['last_s']),previous_s=float(status['previous_s']),transfers=int(status['transfers']))
        for key in ['zt','zb','hc','vt','vb','vnet','vtop_corr','vnet_corr','top_error','net_error_scale','net_error_relative','h_corr','h_inventory','solid','dissolved','phi_mean','core_phi','core_std','core_volume','phib','epsmin','epsmax']:
            row['final_'+key]=float(last[key])
        for key in ['salt_res','water_res','solid_res','liquid_salt_res','crust_salt_res','liquid_water_res','crust_water_res','contact_res','water_jump']:
            row['max_'+key]=float(h[key].abs().max())
        row['late_core_relative_std']=float((late.core_std/late.core_phi).max())
        row['late_core_relative_drift']=float((late.core_phi.max()-late.core_phi.min())/late.core_phi.mean())
        row['core_theory_error']=float(abs(last.core_phi/PHI-1))
        summaries.append(row)
    pd.DataFrame({'id':pending}).to_csv(RESULT/'pending.csv',index=False)
    if not histories:raise RuntimeError('no finished production histories')
    h=pd.concat(histories,ignore_index=True);p=pd.concat(profiles,ignore_index=True);s=pd.DataFrame(summaries)
    h.to_csv(RESULT/'history.csv.gz',index=False,compression={'method':'gzip','mtime':0})
    p.to_csv(RESULT/'profiles.csv.gz',index=False,compression={'method':'gzip','mtime':0})
    s.to_csv(RESULT/'summary.csv',index=False)
    conv=[];base=cases[(cases.kappa==1)&(cases.gap==1e-9)&(cases.rtol==1e-8)&(cases.nl==cases.nc//4)&cases.id.str.contains('_n')]
    for u,group in base.groupby('u'):
        group=group[group.id.isin(s.id)].sort_values('nc')
        for rank,case in enumerate(group.itertuples()):
            if rank+1==len(group):continue
            refid=group.iloc[rank+1].id;ref=h[h.id==refid]
            a=h[(h.id==case.id)&(h.t>=1)&(h.t<=ref.t.max())];row={'id':case.id,'reference_id':refid,'u':u,'nc':case.nc,'common_end_s':a.t.max()}
            scale=np.interp(a.t,ref.t,ref.v_reference)
            for key in ['vt','vb','vnet','zt','zb','hc','solid','dissolved','phib']:
                rv=np.interp(a.t,ref.t,ref[key]);delta=abs(a[key]-rv)
                row[key+'_max_abs']=float(delta.max())
                if key in ['vt','vb','vnet']:denom=scale
                elif key in ['zt','zb','hc']:denom=np.interp(a.t,ref.t,ref.hc)
                else:denom=np.maximum(abs(rv),1e-20)
                row[key+'_scaled_max']=float(np.max(delta/denom))
                if key=='vnet':
                    mask=abs(rv)>.05*scale
                    row['vnet_relative_max_away_zero']=float(np.max(np.asarray(delta)[mask]/abs(rv[mask]))) if np.any(mask) else np.nan
            for field in ['contact_res','water_jump','liquid_salt_res','crust_salt_res','liquid_water_res','crust_water_res']:
                row[field+'_max_abs']=float(a[field].abs().max())
                row['reference_'+field+'_max_abs']=float(ref[field].abs().max())
            conv.append(row)
    cv=pd.DataFrame(conv);cv.to_csv(RESULT/'convergence.csv',index=False)
    sensitivity=[]
    for case in cases[~cases.id.str.contains('_n')].itertuples():
        if case.id not in set(s.id):continue
        refid='x008_n128' if case.u<0 else f'u{case.u:g}_n128'
        a=h[(h.id==case.id)&(h.t>=1)];b=h[h.id==refid];a=a[a.t<=b.t.max()]
        row={'id':case.id,'reference_id':refid,'common_end_s':a.t.max()}
        for key in ['vt','vb','vnet','zt','zb','hc','solid','dissolved']:
            bv=np.interp(a.t,b.t,b[key]);den=np.interp(a.t,b.t,b.v_reference if key in ['vt','vb','vnet'] else b.hc) if key not in ['solid','dissolved'] else abs(bv)
            row[key+'_scaled_max']=float(np.max(abs(a[key]-bv)/den))
        sensitivity.append(row)
    pd.DataFrame(sensitivity).to_csv(RESULT/'sensitivity.csv',index=False)
    physical=s[(s.kappa==1)&(s.gap==1e-9)&(s.rtol==1e-8)&s.id.str.contains('_n')].sort_values('nc').groupby('u').tail(1)
    report={'accounted_cases':len(s),'completed_cases':int((s.status=='end').sum()),'pending':pending,'status_counts':s.status.value_counts().to_dict(),
      'theoretical_phi_core':PHI,'ell_m':ELL,'Dastar':X*X,
      'budget_max':float(s[[col for col in s if col.startswith('max_') and col not in ['max_contact_res','max_water_jump']]].max().max()),
      'finest_cases':physical.to_dict(orient='records')}
    (RESULT/'summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    plot(h,p,physical,cv)
    print(s[['id','status','end_s','final_hc','final_vt','final_vb','final_vnet','final_net_error_relative']].to_string(index=False))

def plot(h,p,base,conv):
    plt.rcParams.update({'font.family':'serif','font.serif':['Times New Roman'],'font.size':9,'svg.fonttype':'none','pdf.fonttype':42,
      'mathtext.fontset':'custom','mathtext.rm':'Times New Roman','mathtext.it':'Times New Roman:italic','mathtext.bf':'Times New Roman:bold'})
    def save(fig,name):
        for ext in ['png','pdf','svg']:fig.savefig(FIG/f'{name}.{ext}',dpi=250)
        path=FIG/f'{name}.svg';path.write_text('\n'.join(x.rstrip() for x in path.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8');plt.close(fig)
    def label(row):return 'x008' if row.u<0 else f'u={row.u:g}'
    base=base.sort_values('u',key=lambda x:np.where(x<0,2,x))
    for kind,name in [('geometry','figure1_interfaces'),('net','figure2_net_growth'),('top','figure3_top_growth')]:
        fig,axs=plt.subplots(2,3,figsize=(11,6),layout='constrained')
        for ax,row in zip(axs.flat,base.itertuples()):
            d=h[h.id==row.id]
            if kind=='geometry':
                for col,text,style in [('zt','Top','-'),('zb','Bottom','--'),('hc','Crust thickness',':')]:ax.plot(d.t/3600,d[col]*1e6,style,label=text)
                ylabel='Position / thickness (um)'
            elif kind=='net':
                ax.plot(d.t/3600,d.vnet*1e9,label='ALE: top minus bottom');ax.plot(d.t/3600,d.vnet_corr*1e9,'--',label='Corrected inventory law');ylabel='Net-thickening rate (nm/s)'
            else:
                ax.plot(d.t/3600,d.vt*1e9,label='ALE top');ax.plot(d.t/3600,d.vtop_corr*1e9,'--',label='Law with soil B');ax.plot(d.t/3600,d.vtop_contact_diagnostic*1e9,':',label='Law with contact flux');ylabel='Top speed (nm/s)'
            ax.set(xlabel='Time (h)',ylabel=ylabel,title=label(row));ax.legend(fontsize=7)
        axs.flat[-1].axis('off');axs.flat[-1].text(.05,.5,'Conditional moving-front model\nu=0.6: strict-tolerance runs stopped\nQuantitative validation remains unresolved',transform=axs.flat[-1].transAxes,fontsize=9,va='center');save(fig,name)
    row=base[base.u==.6].iloc[0];q=p[p.id==row.id];times=q.t.unique();selected=sorted({times[np.argmin(abs(times-t))] for t in [3600,18000,72000,180000]})
    fig,axs=plt.subplots(2,2,figsize=(9,6),layout='constrained')
    for ax,t in zip(axs.flat,selected):
        a=q[q.t==t].sort_values('zlo');state=h[(h.id==row.id)&np.isclose(h.t,t)].iloc[0]
        ax.stairs(a.phi.to_numpy(),np.r_[a.zlo.to_numpy(),a.zhi.iloc[-1]]*1e6,label='Computed solid fraction')
        ax.axvline(state.zb*1e6,color='C1',ls='--',label='Bottom');ax.axvline(state.zt*1e6,color='C2',ls=':',label='Top');ax.axhline(PHI,color='.5',ls='--',label='Theoretical core')
        ax.set(xlabel='Laboratory height (um)',ylabel='Solid fraction',title=f'u=0.6, {t/3600:g} h',ylim=(-.03,1));ax.legend(fontsize=7)
    fig.suptitle('u=0.6: time-tolerance verification unresolved',fontsize=10);save(fig,'figure4_solid_profiles')
    fig,axs=plt.subplots(2,3,figsize=(11,6),layout='constrained')
    for ax,row in zip(axs.flat,base.itertuples()):
        d=h[h.id==row.id]
        for col,text in [('G','Top precipitation'),('Db','Interface dissolution'),('Rplus','Internal precipitation'),('Dminus','Internal dissolution')]:ax.plot(d.t/3600,d[col]*1e6,label=text)
        ax.plot((d.t.to_numpy()[1:]+d.t.to_numpy()[:-1])/7200,np.diff(d.dissolved)/np.diff(d.t)*1e6,':',label='Dissolved storage secant')
        ax.set(xlabel='Time (h)',ylabel='Salt rate (1e-6 kg/m2/s)',title=label(row));ax.legend(fontsize=6.5)
    axs.flat[-1].axis('off');axs.flat[-1].text(.05,.5,'Conditional moving-front model\nu=0.6: strict-tolerance runs stopped\nQuantitative validation remains unresolved',transform=axs.flat[-1].transAxes,fontsize=9,va='center');save(fig,'figure5_redistribution')
    fig,axs=plt.subplots(2,3,figsize=(11,6),layout='constrained')
    for ax,row in zip(axs.flat,base.itertuples()):
        d=h[h.id==row.id];ax.plot(d.t/3600,d.hc*1e6,label='Geometric crust thickness');ax.plot(d.t/3600,d.h_inventory*1e6,'--',label='Inventory-equivalent thickness');ax.plot(d.t/3600,d.h_corr*1e6,':',label='Integrated corrected law')
        twin=ax.twinx();twin.plot(d.t/3600,d.solid,color='.6',ls='-.',label='Crystal mass');twin.set_ylabel('Crystal mass (kg/m2)',color='.4')
        ax.set(xlabel='Time (h)',ylabel='Thickness (um)',title=label(row));ax.legend(fontsize=6.5)
    axs.flat[-1].axis('off');axs.flat[-1].text(.05,.5,'Conditional moving-front model\nu=0.6: strict-tolerance runs stopped\nQuantitative validation remains unresolved',transform=axs.flat[-1].transAxes,fontsize=9,va='center');save(fig,'figure6_geometry_inventory')
    fig,axs=plt.subplots(2,3,figsize=(11,6),layout='constrained')
    for ax,col,text in zip(axs.flat,['vt','vb','vnet','zt','zb','hc'],['Top velocity','Bottom velocity','Net velocity','Top position','Bottom position','Crust thickness']):
        for u,d in conv.groupby('u'):
            d=d.sort_values('nc');ax.plot(d.nc,d[col+'_scaled_max']*100,'o-',label='x008' if u<0 else f'u={u:g}')
        ax.set(xlabel='Coarser grid cells (N to 2N)',ylabel='Max |delta v| / vref (%)' if col in ['vt','vb','vnet'] else 'Max |delta z| / h_fine (%)',title=text);ax.legend(fontsize=6.5)
    save(fig,'figure7_convergence')
    row=base[base.u<0].iloc[0];d=h[h.id==row.id];fig,axes=plt.subplots(2,2,figsize=(11,6.5),layout='constrained');axs=axes.flat
    axs[0].plot(d.t/3600,d.vnet*1e9,label='ALE net');axs[0].plot(d.t/3600,d.vnet_corr*1e9,'--',label='Corrected net');axs[0].plot(d.t/3600,d.vt*1e9,':',label='ALE top');axs[0].plot(d.t/3600,d.vtop_corr*1e9,'-.',label='Corrected top')
    axs[0].set(ylabel='Velocity (nm/s)',title='(a) Computed velocity histories')
    for col,text in [('hc','Geometric thickness'),('h_corr','Corrected thickness'),('h_inventory','Inventory equivalent')]:axs[1].plot(d.t/3600,d[col]*1e6,label=text)
    axs[1].set(ylabel='Thickness (um)',title='(b) Geometry and inventory')
    for col,text in [('B','Soil diffusion'),('Fplus','Crust-side diffusion'),('Jsoil','Actual salt supply'),('P_corr','Reconstructed supply')]:
        late=d[d.t>=3600];early=d[d.t<=3600]
        axs[2].plot(late.t/3600,late[col]*1e6,label=text)
        axs[3].plot(early.t,early[col]*1e6,label=text)
    axs[2].set(ylabel='Salt flux (1e-6 kg/m2/s)',title='(c) Fluxes: 1 h onward')
    axs[3].set(ylabel='Salt flux (1e-6 kg/m2/s)',title='(d) Startup fluxes: not converged',xlabel='Time since precipitation (s)')
    axs[3].set_xscale('symlog',linthresh=.01)
    for ax in list(axs)[:3]:ax.set_xlabel('Time since soil precipitation (h)')
    for ax in axes.flat:ax.legend(fontsize=7)
    save(fig,'figure8_x008')

if __name__=='__main__':main()
