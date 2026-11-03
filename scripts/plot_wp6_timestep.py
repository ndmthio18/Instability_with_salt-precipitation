"""Figures from accepted-state tables; line segments are display only."""
import csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'analysis/wp6_timestep_20261103';BASE=ROOT/'analysis/wp6_execution_20261103'
def rows(p):
    with p.open() as f:return list(csv.DictReader(f))
def nums(rr,k):return np.array([float(r[k]) for r in rr])
def save(fig,name):
    fig.savefig(OUT/(name+'.png'),dpi=180);fig.savefig(OUT/(name+'.pdf'));plt.close(fig)
def main():
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(2,2,figsize=(11,7.5),layout='constrained')
    for d,label,style in [(BASE/'runs/WP6-M1-extend','dt','-'),(OUT/'runs/WP6-M1-dt_half','dt/2','--')]:
        rr=rows(d/'post/weighted_diagnostics.csv');h=nums(rr,'time_s')/3600;s=nums(rr,'weighted_top_sigma')
        ax[0,0].plot(h,s/np.min(s),style,label=label);ax[0,1].plot(h,s,style,label=label)
        ax[1,0].plot(h,nums(rr,'weighted_top_mean'),style,label=label+' mean');ax[1,0].plot(h,nums(rr,'top_max'),style,lw=.7,label=label+' max')
    ax[1,0].axhline(.1,color='k',ls=':',label='Solubility applicability limit')
    profiles=rows(OUT/'vertical_profile_comparison.csv')
    for name,style in [('latest_exact',':'),('own_endpoints','-')]:
        rr=[r for r in profiles if r['report']==name]
        for k,c in [('coarse_sigma','C0'),('fine_sigma','C1')]:
            ax[1,1].plot(nums(rr,k)*1e6,(.2-nums(rr,'z_m'))*1000,style,color=c,label=('dt ' if k=='coarse_sigma' else 'dt/2 ')+('shared time' if name=='latest_exact' else 'own endpoint'))
    ax[0,0].set_ylabel('Weighted amplification from minimum');ax[0,1].set_ylabel('Weighted top sigma (mol/mol)')
    ax[1,0].set_ylabel('Surface mole fraction');ax[1,1].set(xlabel='Horizontal sigma (10^-6 mol/mol)',ylabel='Depth below surface (mm)',ylim=(5,0))
    for a in [ax[0,0],ax[0,1],ax[1,0]]:a.set_xlabel('Actual accepted physical time (h)')
    for a in ax.flat:a.legend(fontsize=8)
    fig.suptitle('A  Canonical M1 growth and profiles under timestep refinement')
    save(fig,'Figure_A_growth_profiles')
    fig,ax=plt.subplots(2,2,figsize=(11,7.5),layout='constrained')
    for reg,label,color in [('upper_0.995000','1.03 mm','C0'),('upper_0.990000','1.95 mm','C1')]:
        for d,style,protocol in [(ROOT/'analysis/wp6_review_20261103','-','dt'),(OUT/'refined_pair','--','dt/2')]:
            rr=[r for r in rows(d/'paired_transport.csv') if r['region']==reg]
            ax[0,0].plot(nums(rr,'time_s')/3600,nums(rr,'delta_C_m'),style,color=color,label=label+' '+protocol)
        rr=[r for r in rows(OUT/'paired_timestep_sensitivity.csv') if r['region']==reg]
        ax[0,1].plot(nums(rr,'common_report_time_for_bound_s')/3600,nums(rr,'U_dt_upper_m'),color=color,label=label+' U_dt upper')
        ax[0,1].plot(nums(rr,'common_report_time_for_bound_s')/3600,nums(rr,'timing_bound_m'),':',color=color,label=label+' timing bound')
    summary=rows(OUT/'seven_plane_timestep_summary.csv');summary.sort(key=lambda r:float(r['depth_mm']));x=np.arange(7)
    for off,key,label in [(-.25,'delta_C_fine_m','|refined delta C|'),(0,'U_dt_upper_m','U_dt upper'),(.25,'fine_residual_envelope_m','Refined closure envelope')]:
        ax[1,0].bar(x+off,abs(nums(summary,key)),width=.25,label=label)
    ax[1,0].set_yscale('log');ax[1,0].set_xticks(x,[f"{float(r['depth_mm']):.2f}" for r in summary],rotation=30);ax[1,0].set(xlabel='Depth (mm)',ylabel='Endpoint magnitude (m)')
    for d,label,style in [(ROOT/'analysis/wp6_review_20261103','dt','-'),(OUT/'refined_pair','dt/2','--')]:
        rr=[r for r in rows(d/'paired_transport.csv') if r['region']=='upper_0.995000'];h=nums(rr,'time_s')/3600
        for key,c in [('M1_signed_cut_adv','C0'),('M1_signed_cut_diff','C1'),('M1_signed_cut_total','C2')]:ax[1,1].plot(h,nums(rr,key),style,color=c,label=label+' '+key.split('_')[-1])
    ax[0,0].set(xlabel='Actual accepted physical time (h)',ylabel='Cumulative excess delta C (m)')
    ax[0,1].set(xlabel='Earlier accepted report time (h)',ylabel='Sensitivity / timing bound (m)',yscale='log')
    ax[1,1].set(xlabel='Actual accepted physical time (h)',ylabel='1.03 mm M1 downward flux (m/s)');ax[1,1].axhline(0,color='.5',lw=.5)
    for a in ax.flat:a.legend(fontsize=7)
    fig.suptitle('B  Paired transport sensitivity; spatial and domain errors remain unknown')
    save(fig,'Figure_B_transport_sensitivity')
if __name__=='__main__':main()
