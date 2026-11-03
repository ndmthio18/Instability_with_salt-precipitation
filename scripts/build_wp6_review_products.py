"""Local evidence summary and figures; no simulations, no frozen output writes."""
import csv,hashlib,json,re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from wp6_paired_stopping import validate_leg
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'analysis/wp6_review_20261103';BASE=ROOT/'analysis/wp6_execution_20261103';CACHE=ROOT/'tmp/wp6_review_20261103'
def load(p):return json.loads(p.read_text())
def rows(p):return list(csv.DictReader(p.open()))
def num(rr,k):return np.array([float(r[k]) for r in rr])
def save(name,v):(OUT/name).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def widths(x):return np.diff(np.r_[x[0],(x[1:]+x[:-1])/2,x[-1]])
def figsave(fig,name):
    fig.savefig(OUT/(name+'.png'),dpi=180);fig.savefig(OUT/(name+'.pdf'));plt.close(fig)


def main():
    summaries={};data={};series={};structure={}
    for name in ['baseline','half']:
        d=OUT/('background_'+name);summaries[name]=load(d/'noise_summary.json');data[name]=np.load(CACHE/('background_'+name)/'noise_fields.npz');series[name]=rows(d/'noise_timeseries.csv')
        a=data[name];rr=series[name];t=a['time_s'];late=t>=136423.0637723484
        sd=num(rr,'moleFraction_top_sigma');slope=np.polyfit(t[late],np.log(sd[late]),1)[0]
        wz=widths(a['zs']);energy=a['depth_sigma'][-1]**2*wz;depth=(a['zs'][-1]-a['zs'])[::-1]
        d90=depth[np.searchsorted(np.cumsum(energy[::-1])/sum(energy),.9)]
        modes=a['volume_mode_energy'];coeff=a['top_neumann_coeff'];first=int(np.argmax(modes[1])+1)
        wx=widths(a['xs']);end=a['top_anomaly'][-1]
        structure[name]=dict(late_log_noise_slope_per_s=float(slope),late_noise_gain=float(sd[-1]/sd[np.flatnonzero(late)[0]]),depth_containing_90percent_noise_energy_mm=float(d90*1000),first_step_dominant_mode=first,final_dominant_mode=int(np.argmax(modes[-1])+1),fraction_final_top_energy_in_two_endpoint_quadrature_weights=float(np.sum(end[[0,-1]]**2*wx[[0,-1]])/np.sum(end**2*wx)),final_top_neumann_m2=float(coeff[-1,1]),endpoint_first_exceeds_1e12_time_s=next((float(r['time_s']) for r in rr if float(r['moleFraction_max_endpoint'])>1e-12),None))
    m1=rows(BASE/'runs/WP6-M1-extend/post/weighted_diagnostics.csv');basew=rows(BASE/'runs/WP6-M0-extend/post/weighted_diagnostics.csv');halfw=rows(OUT/'runs/WP6-M0-dt_half/post/weighted_diagnostics.csv')
    base_events=rows(BASE/'runs/WP6-M0-extend/surface_events.csv');half_events=rows(OUT/'runs/WP6-M0-dt_half/surface_events.csv')
    control=load(OUT/'runs/WP6-M0-dt_half/finished.json');s=control['solver']
    leg=validate_leg(half_events,s['stop'],s['returncode'])
    checks=load(OUT/'runs/WP6-M0-dt_half/post/independent_checks.json')
    b=rows(OUT/'runs/WP6-M0-dt_half/transport_budget.csv')
    metrics=dict(structure=structure,
        endpoint_reduction_factor=summaries['baseline']['maxima']['moleFraction_max_endpoint']/summaries['half']['maxima']['moleFraction_max_endpoint'],
        sigma_reduction_factor=summaries['baseline']['max_top_noise_sigma']/summaries['half']['max_top_noise_sigma'],
        final_top_spatial_correlation=float(np.corrcoef(data['baseline']['top_anomaly'][-1],data['half']['top_anomaly'][-1])[0,1]),
        M0_sigma_over_M1_sigma_max=float(max(num(series['baseline'],'moleFraction_top_sigma')/num(m1,'weighted_top_sigma'))),
        M0_max_sigma_over_imposed_amplitude=summaries['baseline']['max_top_noise_sigma']/1e-6,
        M0_half_sigma_over_imposed_amplitude=summaries['half']['max_top_noise_sigma']/1e-6,
        control_validation=leg,independent_checks_passed=checks['passed'],
        control_max_salt_residual_m=float(max(abs(num(b,'salt_residual_m')))),control_max_water_residual_m=float(max(abs(num(b,'water_residual_m')))),
        original_background_screen='failed_retained',half_background_endpoint_screen='passed',symmetry_threshold_amended=False,
        final_mean_M0_base=float(basew[-1]['weighted_top_mean']),final_mean_M0_half=float(halfw[-1]['weighted_top_mean']),
        nominal_endpoint_mean_difference=float(halfw[-1]['weighted_top_mean'])-float(basew[-1]['weighted_top_mean']),
        actual_endpoint_time_difference_s=float(half_events[-1]['time_s'])-float(base_events[-1]['time_s']),
        endpoint_comparison_scope='nominal 2T endpoint sensitivity only, not a conservative M0/M1 pair',
        solver_max_relative_shift=1e-8,solver_min_newton_steps=2,solver_tolerance_varied=False,
        mechanism='timestep-dependent numerical background; solver, assembly and boundary contributions not uniquely separated',
        physical_numerical_instability_attribution='unresolved; growth of numerical noise observed, no clean eigenmode identified')
    for name,path in [('baseline',BASE/'runs/WP6-M0-extend/solver.log'),('half',OUT/'runs/WP6-M0-dt_half/solver.log')]:
        log=path.read_text();assert 'Newton.MaxRelativeShift = 1e-08' in log
        vals=re.findall(r'Newton iteration\s+\d+ done, maximum relative shift = ([0-9.e+-]+)',log)
        metrics[name+'_last_logged_newton_shift']=float(vals[-1])
    save('background_comparison.json',metrics)
    save('control_acceptance_and_cost.json',dict(accepted=bool(leg['valid'] and checks['passed']),original_gate_retained=True,
        steps=s['stop']['accepted_steps'],wall_s=s['wall_s'],cpu_s=s['cpu_s'],cpu_hours=s['cpu_s']/3600,peak_rss_kib=s['peak_rss_kib'],post_wall_s=control['post_wall_s'],raw_output_bytes=s['solver_output_bytes'],archive_bytes=control['archive']['bytes'],archive_sha256=control['archive']['sha256'],verified_vtu=2401,only_one_new_solver_run=True))
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(2,2,figsize=(11,7.5),layout='constrained')
    for name,label in [('baseline','M0 dt'),('half','M0 dt/2')]:
        rr=series[name];h=num(rr,'time_s')/3600
        ax[0,0].plot(h,num(rr,'moleFraction_max_endpoint'),label=label)
        ax[0,1].plot(h,num(rr,'moleFraction_top_sigma'),label=label)
        ax[1,0].plot(data[name]['xs']*1000,data[name]['top_anomaly'][-1]*1e12,'o-',ms=3,label=label)
        frac=np.array(summaries[name]['final_volume_mode_energy_fractions']);offset=-.18 if name=='baseline' else .18
        ax[1,1].bar(np.arange(1,25)+offset,frac,width=.36,label=label)
    ax[0,0].axhline(1e-12,color='k',ls='--',label='Frozen endpoint screen')
    for j in [0,1]:ax[0,j].set_xlabel('Accepted physical time (h)');ax[0,j].ticklabel_format(axis='y',style='sci',scilimits=(0,0));ax[0,j].legend(fontsize=8)
    ax[0,0].set_ylabel('Max endpoint difference (mol/mol)');ax[0,1].set_ylabel('Weighted top noise sigma (mol/mol)')
    ax[1,0].set(xlabel='Horizontal position (mm)',ylabel='Final top anomaly (10^-12 mol/mol)');ax[1,0].legend()
    ax[1,1].set(xlabel='Neumann cosine mode index m',ylabel='Final volume noise energy fraction');ax[1,1].legend()
    fig.suptitle('A  Background discrepancy: endpoint reduction is not full-field symmetry recovery')
    figsave(fig,'Figure_A_background')
    fig,ax=plt.subplots(1,3,figsize=(12,4),layout='constrained')
    peak=max(np.max(abs(data[k]['x_anomaly_'+str(len(data[k]['time_s'])-1)])) for k in data)*1e12
    for k,(name,label) in enumerate([('baseline','M0 dt'),('half','M0 dt/2')]):
        d=data[name];zdepth=(.2-d['zs'])*1000;field=d['x_anomaly_'+str(len(d['time_s'])-1)]*1e12
        im=ax[k].pcolormesh(d['xs']*1000,zdepth,field,cmap='RdBu_r',vmin=-peak,vmax=peak,shading='nearest');ax[k].set(xlabel='Horizontal position (mm)',ylabel='Depth below surface (mm)',ylim=(2,0),title=label)
    fig.colorbar(im,ax=list(ax[:2]),label='Concentration anomaly (10^-12 mol/mol)',shrink=.8)
    for name,label in [('baseline','M0 dt'),('half','M0 dt/2')]:
        d=data[name];ax[2].plot(d['depth_sigma'][-1]*1e12,(.2-d['zs'])*1000,label=label)
    ax[2].set(xlabel='Horizontal sigma (10^-12 mol/mol)',ylabel='Depth (mm)',ylim=(2,0));ax[2].legend()
    fig.suptitle('B  Final background noise: top 2 mm (full domain was analyzed)');figsave(fig,'Figure_B_spatial')
    pairs=rows(OUT/'paired_transport.csv');review=rows(OUT/'seven_plane_review.csv')
    fig,ax=plt.subplots(2,2,figsize=(11,7.5),layout='constrained')
    for reg,label,color in [('upper_0.995000','1.03 mm','C0'),('upper_0.990000','1.95 mm','C1')]:
        rr=[r for r in pairs if r['region']==reg];h=num(rr,'time_s')/3600
        ax[0,0].plot(h,num(rr,'delta_C_m'),label=label,color=color);ax[0,0].plot(h,num(rr,'residual_envelope_m'),ls=':',color=color)
        ax[0,1].plot(h,num(rr,'relative_delta_C_initial_inventory'),label=label,color=color)
    rr=[r for r in pairs if r['region']=='upper_0.995000'];h=num(rr,'time_s')/3600
    for key,label in [('delta_C_adv_m','Advective difference'),('delta_C_diff_m','Diffusive difference'),('delta_C_m','Total difference')]:ax[1,0].plot(h,num(rr,key),label=label)
    for key,label in [('M1_signed_cut_adv','M1 advective'),('M1_signed_cut_diff','M1 diffusive'),('M1_signed_cut_total','M1 total')]:ax[1,1].plot(h,num(rr,key),label=label)
    for a in ax.flat:a.axhline(0,color='.5',lw=.5);a.set_xlabel('Shared accepted physical time (h)');a.ticklabel_format(axis='y',style='sci',scilimits=(0,0));a.legend(fontsize=8)
    ax[0,0].set_ylabel('Cumulative difference / closure (m)');ax[0,0].set_title('Solid: delta C; dotted: residual envelope')
    ax[0,1].set_ylabel('Delta C / initial M0 inventory');ax[1,0].set_ylabel('1.03 mm cumulative difference (m)');ax[1,1].set_ylabel('1.03 mm signed downward flux (m/s)')
    fig.suptitle('C  Persistent shallow differences remain tiny; net flux is still upward');figsave(fig,'Figure_C_transport')
    print(json.dumps({k:v for k,v in metrics.items() if k not in ['control_validation']},indent=2))


if __name__=='__main__':main()
