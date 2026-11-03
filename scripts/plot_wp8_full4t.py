"""Command 15 figures. Reads only saved, hash-verified products; launches nothing."""
import json,re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import wp8_full4t_campaign as f
from analyze_wp8 import rows,pairs
from single_wavelength_post import params

ROOT=f.ROOT;OUT=f.OUT;CACHE=f.CACHE;STAGE=f.STAGE;T=f.T
HIST=ROOT/'analysis/wp8_strict_background_20261103'
HCACHE=ROOT/'tmp/wp8_strict_background_20261103'
NEW=OUT/STAGE
BLUE='#286f9e';RED='#b4442e';GREY='#999999';GREEN='#3f7d52';ORANGE='#d98c2b'


def iterseries(log):
    return np.array([int(x)+1 for x in re.findall(r'WP8_CONSTRAINT iterations=(\d+)',log)])


def figure_validation():
    a=dict(np.load(CACHE/STAGE/'post/profiles.npz'));h=dict(np.load(HCACHE/'4T/post/profiles.npz'))
    nev=rows(NEW/'surface_events.csv');hev=rows(HIST/'4T/surface_events.csv')
    nt=[float(r['time_s']) for r in nev];ht=[float(r['time_s']) for r in hev]
    pp=pairs(ht,nt);ti=np.array([nt[j] for i,j in pp])/T
    dx=np.array([np.abs(h['x_mean'][i]-a['x_mean'][j]).max() for i,j in pp])
    dp=np.array([np.abs(h['pw_mean'][i]-a['pw_mean'][j]).max() for i,j in pp])

    fig,axs=plt.subplots(2,2,figsize=(10.5,7),layout='constrained')

    ax=axs[0,0]
    ax.semilogy(ti,np.maximum(dx,1e-20),color=BLUE,lw=1,label=r'max$|\Delta x|$')
    ax.semilogy(ti,np.maximum(dp,1e-20),color=RED,lw=1,label=r'max$|\Delta p_w|$ (Pa)')
    ax.axvline(ht[-1]/T,color='black',ls=':',lw=.9)
    ax.text(ht[-1]/T,3e-10,' historical last\n accepted 3.6775T',fontsize=7,va='top')
    ax.set(xlabel='$t/T$',ylabel='Max absolute difference',
           title='A  Corrected vs historical strict background\n(exactly matched accepted times)')
    ax.legend(fontsize=8,loc='lower right')

    ax=axs[0,1]
    nit=iterseries((NEW/'solver.log').read_text());hit=iterseries((HIST/'4T/solver.log').read_text())
    ax.plot(np.arange(1,len(hit)+1)*113.68588647695695/T,hit,'.',ms=1.4,color=GREY,label=f'Historical (legacy), max {hit.max()}')
    ax.plot(np.arange(1,len(nit)+1)*113.68588647695695/T,nit,'.',ms=1.4,color=BLUE,label=f'Corrected (scaled), max {nit.max()}')
    ax.axvline(len(hit)*113.68588647695695/T,color='black',ls=':',lw=.9)
    ax.set(xlabel='$t/T$',ylabel='Newton iterations to acceptance',
           title='B  Accepted-step Newton cost',ylim=(0,19))
    ax.legend(fontsize=8,loc='upper left')

    ax=axs[1,0]
    fail=rows(HIST/'failed_step_iterations.csv')
    ax.semilogy([int(r['iteration']) for r in fail],[float(r['salt_residual']) for r in fail],'.-',color=RED,lw=.9,label='Rejected step 4414, salt residual')
    ax.axhline(2e-20,color='black',ls='--',lw=.9,label='Old fixed absolute limit $2\\times10^{-20}$')
    di=rows(NEW/'convergence_v2.csv');w=[r for r in di if int(r['step'])==4414]
    ax.semilogy([int(r['iteration']) for r in w],[float(r['salt_max']) for r in w],'o-',color=BLUE,ms=4,lw=.9,label='Corrected step 4414, accepted')
    ax.set(xlabel='Newton iteration',ylabel='Absolute salt residual (m)',
           title='C  The old failure window')
    ax.legend(fontsize=7.5,loc='upper right')

    ax=axs[1,1]
    bud=rows(NEW/'transport_budget.csv');reg=sorted({r['region'] for r in bud})
    for g,col in zip(['free_domain','upper_0.950000','upper_0.995000'],[BLUE,GREEN,ORANGE]):
        rr=[r for r in bud if r['region']==g]
        ax.semilogy([float(r['time_s'])/T for r in rr],[max(abs(float(r['salt_residual_m'])),1e-22) for r in rr],lw=.9,color=col,label=g)
    ax.axhline(1e-12,color='black',ls='--',lw=.9,label='Unchanged salt budget gate')
    ax.set(xlabel='$t/T$',ylabel='Regional salt closure residual (m)',
           title='D  Conservative closure, every accepted step',ylim=(1e-22,1e-11))
    ax.legend(fontsize=7.5,loc='lower right')

    for ax in axs.flat:ax.spines[['top','right']].set_visible(False);ax.grid(alpha=.15)
    fig.savefig(OUT/'Figure_A_background_validation.png',dpi=180)
    fig.savefig(OUT/'Figure_A_background_validation.pdf');plt.close(fig)


def figure_response():
    sci=json.loads((OUT/'scientific_analysis.json').read_text())
    resp=rows(OUT/'background_referenced_responses.csv')
    grid=rows(OUT/'grid_sensitivity_vs_background.csv')
    nev=rows(NEW/'surface_events.csv')
    runs={'R24_M0':'analysis/wp6_4T_exploratory_20261103/runs/WP6-M0-explore4T',
          'R24_M1':'analysis/wp6_4T_exploratory_20261103/runs/WP6-M1-explore4T',
          'R48_M0':'analysis/wp7_shortwave_validation_20261103/runs/WP7-R48-M0-4T',
          'R48_M1':'analysis/wp7_shortwave_validation_20261103/runs/WP7-R48-M1-4T'}

    fig,axs=plt.subplots(2,2,figsize=(10.5,7),layout='constrained')

    ax=axs[0,0]
    t=np.array([float(r['time_s']) for r in nev])/T
    ax.plot(t,[float(r['x_top_area_mean']) for r in nev],color='black',lw=1.4,label='Strict background (uniform)')
    for key,col in zip(['R24_M1','R48_M1'],[RED,BLUE]):
        ev=rows(ROOT/runs[key]/'surface_events.csv')
        ax.plot([float(r['time_s'])/T for r in ev],[float(r['x_top_max']) for r in ev],lw=1,color=col,label=key+' local max')
        ax.plot([float(r['time_s'])/T for r in ev],[float(r['x_top_area_mean']) for r in ev],lw=1,ls='--',color=col,alpha=.7,label=key+' surface mean')
    ax.axhline(.1,color=GREY,ls='--',lw=.9);ax.text(.1,.1015,'solubility threshold $x=0.1$',fontsize=7.5,color=GREY)
    ax.set(xlabel='$t/T$',ylabel='Top-row salt mole fraction',ylim=(0,.108),
           title='A  Surface mean versus local maximum')
    ax.legend(fontsize=7,loc='upper left')

    ax=axs[0,1]
    for name,col,ls in [('R24_B',RED,'--'),('R24_C',RED,'-'),('R48_B',BLUE,'--'),('R48_C',BLUE,'-')]:
        rr=[r for r in resp if r['comparison']==name]
        ax.semilogy([float(r['t_over_T']) for r in rr],[max(abs(float(r['mean_response_top'])),1e-20) for r in rr],
                    lw=1,color=col,ls=ls,label=name+(' (M0 - bg)' if name.endswith('B') else ' (M1 - bg)'))
    ax.axhline(9.02e-16,color='black',ls=':',lw=.9)
    ax.text(.1,1.3e-15,'solver-criterion floor',fontsize=7)
    ax.axvline(3.6775,color=GREY,ls=':',lw=.9)
    ax.set(xlabel='$t/T$',ylabel=r'$|\overline{x}_{\rm case}-\overline{x}_{\rm bg}|$ at the surface',
           ylim=(1e-18,1e-2),title='B  Background-referenced zero-mode response')
    ax.legend(fontsize=7.5,loc='upper left')

    ax=axs[1,0]
    for typ,col in [('B',RED),('C',BLUE)]:
        rr=[r for r in grid if r['comparison_type']==typ]
        y=[abs(float(r['R24_mean_response_top'])/float(r['R48_mean_response_top'])) if float(r['R48_mean_response_top']) else np.nan for r in rr]
        ax.semilogy([float(r['t_over_T']) for r in rr],y,lw=1,color=col,label='R24/R48, type '+typ)
    ax.axhline(1,color='black',ls='--',lw=.9,label='grid-insensitive')
    ax.axvline(3.,color=GREY,ls=':',lw=.9);ax.text(3.02,2e2,'3T',fontsize=7.5,color=GREY)
    ax.set(xlabel='$t/T$',ylabel='|R24 response| / |R48 response|',ylim=(1e-2,1e4),
           title='C  Grid sensitivity of the inferred response\n(separate from background validation)')
    ax.legend(fontsize=7.5,loc='upper left')

    ax=axs[1,1]
    red=sci['seven_plane_redistribution'];regs=['upper_0.500000','upper_0.750000','upper_0.900000','upper_0.950000','upper_0.975000','upper_0.990000','upper_0.995000']
    depth=[(1-float(g.split('_')[1]))*200 for g in regs]
    for name,col,m in [('R24_B',RED,'s'),('R24_C',RED,'o'),('R48_B',BLUE,'s'),('R48_C',BLUE,'o')]:
        y=[abs(red[name]['4T'][g]['delta_I_over_initial']) for g in regs]
        ax.semilogx(np.maximum(y,1e-16),depth,m+'-',ms=4,lw=.9,color=col,
                    markerfacecolor='none' if name.endswith('B') else col,label=name)
    ax.invert_yaxis()
    ax.set(xlabel=r'$|\Delta I_{\rm corrected}| / I_0$ at $4T$',ylabel='Plane depth below surface (mm)',
           xlim=(1e-15,1),title='D  Seven-plane initial-inventory-corrected\nsalt redistribution at $4T$')
    ax.legend(fontsize=7.5,loc='lower right')

    for ax in axs.flat:ax.spines[['top','right']].set_visible(False);ax.grid(alpha=.15)
    fig.savefig(OUT/'Figure_B_background_referenced_response.png',dpi=180)
    fig.savefig(OUT/'Figure_B_background_referenced_response.pdf');plt.close(fig)


if __name__=='__main__':
    figure_validation();figure_response();print('figures written to',OUT)
