"""WP9 figures. Reads only saved, hash-verified products; launches nothing."""
import collections,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import wp9_common as c

BLUE='#286f9e';RED='#b4442e';GREY='#999999';GREEN='#3f7d52';ORANGE='#d98c2b';PURPLE='#6c4f9e'
COL={'R24':RED,'R48':BLUE,'R96':GREEN}
CHECK=[('T',1200),('2T',2400),('3T',3600),('3.6775T',4413),('4T',4800)]


def save(fig,name):
    for ext in ('png','pdf'):fig.savefig(c.OUT/f'{name}.{ext}',dpi=150,bbox_inches='tight')
    plt.close(fig)


def series():
    by=collections.defaultdict(list)
    for r in c.rows(c.OUT/'resolution_series.csv'):by[r['leg']].append(r)
    return by


def arr(rows,key):
    return np.array([float(r[key]) if r[key] not in ('',None) else np.nan for r in rows])


def figure_a():
    lsa=json.loads((c.OUT/'lsa_audit.json').read_text(encoding='utf-8'))
    scan=c.rows(c.OUT/'lsa_wavenumber_scan.csv')
    vs=c.rows(c.OUT/'lsa_versus_dns.csv')
    fig,ax=plt.subplots(2,2,figsize=(12.6,8.8),constrained_layout=True)

    a=ax[0,0]
    for name,col in [('t_c',GREY),('2t_c','#c0c0c0'),('T',PURPLE),('2T',ORANGE),('3T',BLUE),('4T',RED)]:
        rows=[r for r in scan if r['checkpoint']==name]
        nx=arr(rows,'nx');s=arr(rows,'sigma_times_t_hat_c')
        a.plot(nx,s,color=col,lw=1.3,label=name)
    a.set_ylim(-2.0,4.0)
    for m,lab in [(2,'m=2 (A28)'),(6,'m=6'),(10,'m=10')]:
        a.axvline(14*m,color='k',ls=':',lw=.8)
        a.text(14*m,0.97,lab,rotation=90,va='top',ha='right',fontsize=7,
               transform=a.get_xaxis_transform())
    a.axhline(0,color='k',lw=.6);a.set_xscale('log');a.set_xlim(1,320)
    a.set_xlabel('LSA wavenumber index $n_x$');a.set_ylabel(r'$\sigma\,\hat t_c$')
    a.set_title('A  quasi-steady growth rate; the box admits only $n_x=14m$')
    a.legend(fontsize=7,ncol=2)

    b=ax[0,1]
    for name,col in [('t_c',GREY),('T',PURPLE),('2T',ORANGE),('3T',BLUE),('4T',RED)]:
        rows=[r for r in scan if r['checkpoint']==name]
        b.plot(arr(rows,'nx'),arr(rows,'critical_rayleigh_over_R_E'),color=col,lw=1.3,label=name)
    b.axhline(1,color='k',lw=.8,ls='--')
    b.set_xscale('log');b.set_yscale('log');b.set_xlim(1,320)
    b.set_xlabel('LSA wavenumber index $n_x$');b.set_ylabel(r'$R_c(n_x,t)/R_E$')
    b.set_title('B  neutral curve; unstable where the ratio is below one');b.legend(fontsize=7)

    d=ax[1,0]
    names=[k for k in ['t_c','2t_c','T','2T','3T','3.6775T','4T'] if k in lsa['per_checkpoint']]
    x=[lsa['per_checkpoint'][k]['t_over_tc'] for k in names]
    d.plot(x,[lsa['per_checkpoint'][k]['fastest_box_mode_m'] for k in names],'o-',color=BLUE,
           label='fastest box mode $m$')
    d.plot(x,[lsa['per_checkpoint'][k]['wide_scan_peak_nx']/14 for k in names],'s--',color=ORANGE,
           label='unconstrained peak, as $n_x/14$')
    for m,col in [(2,GREY),(6,'#c0c0c0'),(10,'#e0b0b0')]:d.axhline(m,color=col,lw=.8,ls=':')
    d.set_xlabel('$t/t_c$');d.set_ylabel('mode')
    d.set_title('C  the fastest mode drifts to short waves');d.legend(fontsize=7)

    e=ax[1,1]
    for leg in sorted({r['leg'] for r in vs}):
        rows=[r for r in vs if r['leg']==leg]
        k=leg.split('_')[0]
        adm=[r for r in rows if r['admitted_small_perturbation']=='True']
        out=[r for r in rows if r['admitted_small_perturbation']!='True']
        e.plot(arr(adm,'t_over_T'),arr(adm,'ratio'),'o-',color=COL[k],label=f'{k} admitted')
        if out:e.plot(arr(out,'t_over_T'),arr(out,'ratio'),'x',color=COL[k],alpha=.5,label=f'{k} outside window')
    e.axhline(1,color='k',lw=.8,ls='--');e.axhspan(.95,1.05,color=GREEN,alpha=.12)
    e.set_ylim(0.5,1.5);e.set_xlabel('$t/T$')
    off=[r for r in vs if r['admitted_small_perturbation']!='True' and abs(float(r['ratio']))>1.5]
    for r in off:
        e.annotate('%s %s: ratio %.1f, $m{=}2$ fraction %.3f'%(
                       r['leg'].split('_')[0],r['checkpoint'],float(r['ratio']),
                       float(r['canonical_energy_fraction'])),
                   (float(r['t_over_T']),0.56),fontsize=6,ha='right',
                   color=COL[r['leg'].split('_')[0]])
    e.set_ylabel(r'DuMuX $d\ln A_2/dt$  /  LSA $\sigma_{m=2}$')
    e.set_title('D  LSA against DuMuX, small-perturbation window only');e.legend(fontsize=7)
    fig.suptitle('Figure A  LSA wavenumber audit, case 36877',fontsize=11)
    save(fig,'Figure_A_lsa_wavenumber_audit')


def figure_b():
    by=series();spec=c.rows(c.OUT/'resolution_spectra.csv')
    comp=json.loads((c.OUT/'resolution_comparison.json').read_text(encoding='utf-8'))
    fig,ax=plt.subplots(2,2,figsize=(12.6,8.8),constrained_layout=True)

    a=ax[0,0]
    for leg in by:
        k=leg.split('_')[0];m1=leg.endswith('M1')
        t=arr(by[leg],'t_over_T')
        a.plot(t,arr(by[leg],'energy_m10'),color=COL[k],lw=1.4 if m1 else .9,
               ls='-' if m1 else ':',label=f'{leg}')
    a.set_yscale('log');a.set_ylim(1e-38,None)
    a.set_xlabel('$t/T$');a.set_ylabel('$m=10$ depth-integrated energy')
    a.set_title('A  short-wave energy; dotted lines are the unperturbed floors')
    a.legend(fontsize=7,ncol=2)

    b=ax[0,1]
    for leg in by:
        if not leg.endswith('M1'):continue
        k=leg.split('_')[0]
        rows=[r for r in spec if r['leg']==leg and r['checkpoint']=='4T']
        if not rows:continue
        b.plot(arr(rows,'k_rad_m'),arr(rows,'energy'),'o-',ms=3,color=COL[k],label=leg)
        res=[r for r in rows if r['resolvable']!='True']
        if res:b.plot(arr(res,'k_rad_m'),arr(res,'energy'),'x',color=COL[k],ms=5)
    for m,lab in [(2,'m=2'),(6,'m=6'),(10,'m=10')]:
        b.axvline(c.k_phys(m),color='k',ls=':',lw=.8)
    b.set_xscale('log');b.set_yscale('log')
    b.set_xlabel('physical wavenumber $k$ / rad m$^{-1}$');b.set_ylabel('energy at 4T')
    b.set_title('B  spectrum at 4T; crosses are under eight cells per wavelength');b.legend(fontsize=7)

    d=ax[1,0]
    for leg in by:
        if not leg.endswith('M1'):continue
        k=leg.split('_')[0];t=arr(by[leg],'t_over_T')
        d.plot(t,arr(by[leg],'peak_k_rad_m'),color=COL[k],lw=1.2,label=f'{leg} peak')
        d.plot(t,arr(by[leg],'energy_weighted_k_rad_m'),color=COL[k],lw=1.0,ls='--',
               label=f'{leg} energy weighted')
    for m in (2,6,10):d.axhline(c.k_phys(m),color=GREY,ls=':',lw=.8)
    d.set_yscale('log');d.set_xlabel('$t/T$');d.set_ylabel('$k$ / rad m$^{-1}$')
    d.set_title('C  peak and energy-weighted physical wavenumber');d.legend(fontsize=6,ncol=2)

    e=ax[1,1]
    disc=comp['discriminator']
    legs=[k for k in ['R24_M1','R48_M1','R96_M1'] if disc['m10_energy'].get(k) is not None]
    nx=[int(k[1:3].rstrip('_')) if False else int(k.split('_')[0][1:]) for k in legs]
    e.plot(nx,[disc['m10_energy'][k] for k in legs],'o-',color=RED,label='$m=10$ energy, seeded')
    floors=[(int(k.split('_')[0][1:]),v) for k,v in disc['unperturbed_total_energy'].items() if v]
    if floors:
        floors.sort();e.plot([f[0] for f in floors],[f[1] for f in floors],'s--',color=GREY,
                             label='unperturbed total energy')
    e.set_xscale('log',base=2);e.set_yscale('log');e.set_xticks(nx);e.set_xticklabels(nx)
    e.set_xlabel('horizontal cells');e.set_ylabel('energy at 4T')
    e.set_title('D  the registered discriminator');e.legend(fontsize=7)
    fig.suptitle('Figure B  R24 / R48 / R96 at shared accepted times',fontsize=11)
    save(fig,'Figure_B_resolution_comparison')


def figure_c():
    by=series();plane=c.rows(c.OUT/'resolution_seven_plane.csv')
    fig,ax=plt.subplots(2,2,figsize=(12.6,8.8),constrained_layout=True)

    a=ax[0,0]
    for leg in by:
        k=leg.split('_')[0];m1=leg.endswith('M1')
        a.plot(arr(by[leg],'t_over_T'),np.abs(arr(by[leg],'mean_response_top')),
               color=COL[k],lw=1.4 if m1 else .9,ls='-' if m1 else ':',label=leg)
    a.axhline(2.458970527197124e-15,color='k',ls='--',lw=.8)
    a.text(0.05,2.9e-15,'cross-executable resolvability floor',fontsize=6)
    a.set_yscale('log');a.set_xlabel('$t/T$')
    a.set_ylabel(r'$|\langle x\rangle_{top}-\langle x\rangle_{bg}|$')
    a.set_title('A  horizontal-mean response against the frozen strict background')
    a.legend(fontsize=7,ncol=2)

    b=ax[0,1]
    for leg in by:
        if not leg.endswith('M1'):continue
        k=leg.split('_')[0]
        b.plot(arr(by[leg],'t_over_T'),arr(by[leg],'surface_mean'),color=COL[k],lw=1.2,label=f'{leg} mean')
        b.plot(arr(by[leg],'t_over_T'),arr(by[leg],'domain_local_max'),color=COL[k],lw=1.0,ls='--',
               label=f'{leg} local max')
    b.axhline(c.SOLUBILITY,color='k',lw=1.0)
    b.text(0.05,c.SOLUBILITY*1.02,'solubility 0.1; right-censored at 4T',fontsize=7)
    b.set_xlabel('$t/T$');b.set_ylabel('salt mole fraction')
    b.set_title('B  surface mean and local maximum');b.legend(fontsize=6,ncol=2)

    order=['upper_0.995000','upper_0.990000','upper_0.975000','upper_0.950000',
           'upper_0.900000','upper_0.750000','upper_0.500000']
    d=ax[1,0]
    for leg in sorted({r['leg'] for r in plane if r['leg'].endswith('M1')}):
        k=leg.split('_')[0]
        rows=[r for r in plane if r['leg']==leg and r['checkpoint']=='4T']
        y=[abs(float(next(r for r in rows if r['region']==g)['relative_to_initial'])) for g in order if any(r['region']==g for r in rows)]
        d.plot(range(len(y)),y,'o-',color=COL[k],label=leg)
    d.set_yscale('log');d.set_xticks(range(len(order)))
    d.set_xticklabels(['1.0','1.9','5.0','9.9','20','50','100'])
    d.set_xlabel('plane depth / mm');d.set_ylabel(r'$|\Delta C_{corrected}|/I_0$ at 4T')
    d.set_title('C  seven-plane initial-corrected redistribution');d.legend(fontsize=7)

    e=ax[1,1]
    for leg in sorted({r['leg'] for r in plane if r['leg'].endswith('M1')}):
        k=leg.split('_')[0]
        rows=[r for r in plane if r['leg']==leg and r['region']=='upper_0.995000']
        e.plot([float(r['time_s'])/c.T_S for r in rows],
               [float(r['cumulative_cut_salt_m']) for r in rows],'o-',color=COL[k],label=f'{leg} total')
        e.plot([float(r['time_s'])/c.T_S for r in rows],
               [float(r['cumulative_cut_adv_m']) for r in rows],'^--',color=COL[k],ms=4,
               label=f'{leg} advective')
    e.axhline(0,color='k',lw=.6);e.set_xlabel('$t/T$')
    e.set_ylabel('signed cumulative downward salt / m')
    e.set_title('D  signed cumulative transport across the 1 mm plane');e.legend(fontsize=6)
    fig.suptitle('Figure C  background-referenced response and conservative transport',fontsize=11)
    save(fig,'Figure_C_background_and_transport')


def figure_d():
    corr=json.loads((c.OUT/'corrections_v2.json').read_text(encoding='utf-8'))
    import csv as _csv
    budget=collections.defaultdict(list)
    for r in c.rows(c.ROOT/'analysis/wp8_full4T_20261103/run/transport_budget.csv'):
        budget[r['region']].append(r)
    fd=budget['free_domain']
    t=np.array([float(r['time_s']) for r in fd])/c.T_S
    ws=np.array([float(r['whole_salt_inventory_m']) for r in fd])
    bs=np.array([float(r['bottom_salt_reaction_in_m_s']) for r in fd])
    dt=np.diff(np.array([float(r['time_s']) for r in fd]))
    cum=np.concatenate([[0.0],np.cumsum(bs[1:]*dt)])
    fig,ax=plt.subplots(2,2,figsize=(12.6,8.8),constrained_layout=True)

    a=ax[0,0]
    a.plot(t,ws-ws[0],color=BLUE,lw=1.6,label='whole-domain salt increase')
    a.plot(t,cum,color=RED,lw=1.0,ls='--',label='cumulative bottom-boundary input')
    a.set_xlabel('$t/T$');a.set_ylabel('m')
    a.set_title('A  the increase is bottom input, not evaporative concentration');a.legend(fontsize=7)

    b=ax[0,1]
    rows=corr['correction_1_inventory']['per_plane'][:-1]
    x=np.arange(len(rows));w=.38
    b.bar(x-w/2,[abs(r['advective_component_m']) for r in rows],w,color=BLUE,label='upward advective')
    b.bar(x+w/2,[abs(r['diffusive_component_m']) for r in rows],w,color=ORANGE,label='downward diffusive')
    b.plot(x,[r['delta_salt_m'] for r in rows],'ko-',ms=4,label='net stored')
    b.set_yscale('log');b.set_xticks(x)
    b.set_xticklabels([f"{r['depth_mm']:.1f}" for r in rows],fontsize=7)
    b.set_xlabel('plane depth / mm');b.set_ylabel('cumulative salt / m')
    b.set_title('B  near-surface enrichment is flux convergence');b.legend(fontsize=7)

    d=ax[1,0]
    rows=corr['correction_2_censoring']['trajectories']
    d.barh(range(len(rows)),[r['fraction_of_solubility'] for r in rows],color=GREY)
    d.axvline(1.0,color=RED,lw=1.2)
    d.set_yticks(range(len(rows)));d.set_yticklabels([r['run'] for r in rows],fontsize=7)
    d.set_xlim(0,1.1);d.set_xlabel('local maximum $x$ / solubility')
    d.set_title('C  every trajectory is right-censored at 4T')

    e=ax[1,1]
    gap=corr['elf_gap']
    keys=['max_abs_x_top_area_mean','max_abs_x_top_node_mean','max_abs_sigma_top_node',
          'max_abs_x_domain_max']
    e.barh(range(len(keys)),[gap[k] for k in keys],color=GREEN)
    e.set_xscale('log');e.set_yticks(range(len(keys)))
    e.set_yticklabels([k.replace('max_abs_','') for k in keys],fontsize=7)
    e.set_xlabel('maximum absolute difference over 0 to T')
    e.set_title('D  the two executables agree to round-off on the one-cell background')
    fig.suptitle('Figure D  v2 corrections to the command/15 deliverables',fontsize=11)
    save(fig,'Figure_D_corrections')


if __name__=='__main__':
    figure_a();figure_b();figure_c();figure_d()
    print(json.dumps({n:(c.OUT/f'{n}.png').stat().st_size for n in
        ['Figure_A_lsa_wavenumber_audit','Figure_B_resolution_comparison',
         'Figure_C_background_and_transport','Figure_D_corrections']},indent=2))
