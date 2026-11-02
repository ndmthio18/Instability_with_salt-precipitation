"""Summarize historical inventory changes and plot evidence, without solver calls."""
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'analysis/dense_history_20261102'
EVIDENCE=OUT/'evidence';FIG=OUT/'figures';FIG.mkdir(exist_ok=True)
PLAN=json.loads((OUT/'plan.json').read_text())
AUDITS={r['id']:r for r in json.loads((EVIDENCE/'summary.json').read_text())}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,
                     'savefig.dpi':180,'pdf.fonttype':42})


def readcsv(path):
    with path.open(newline='') as f:return list(csv.DictReader(f))


def writecsv(path,rows):
    if rows:
        fields=list(dict.fromkeys(k for row in rows for k in row))
        with path.open('w',encoding='utf-8',newline='') as f:
            w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)


def val(rows,key):return np.array([float(r[key]) for r in rows])


runs,inventory,events,allframes,time_review,initial_checks=[],[],[],[],[],[]
cache={}
for spec in PLAN['runs']:
    name=spec['id'];audit=AUDITS[name]
    if audit['status']=='analysis_failed':raise RuntimeError(name+' failed; review remote evidence')
    directory=EVIDENCE/'runs'/name
    f=readcsv(directory/'frame_metrics.csv');rr=readcsv(directory/'inventory_timeseries.csv')
    sig=np.loadtxt(next(directory.glob('*_sigmatop.csv')),comments='#',ndmin=2)
    original_times=val(f,'time_s');corrected=[];time_checks=[]
    for row,pvd_time in zip(f,original_times):
        k=int(np.argmin(abs(sig[:,0]-pvd_time)));exact=float(sig[k,0])
        rounding_match=float(format(exact,'.6g'))==pvd_time
        nearest_unique=np.sum(abs(sig[:,0]-pvd_time)<=abs(exact-pvd_time)+1e-9)==1
        field_match=float(row['sigmatop_mean_error'])<1e-10 and float(row['sigmatop_sigma_error'])<1e-10
        passed=bool(rounding_match and nearest_unique and field_match)
        time_checks.append(passed);corrected.append(exact if passed else pvd_time)
        time_review.append(dict(id=name,frame=row['frame'],pvd_time_s=pvd_time,matched_diagnostic_time_s=exact,
            absolute_difference_s=abs(exact-pvd_time),six_significant_digit_roundtrip=rounding_match,
            unique_nearest=bool(nearest_unique),top_statistics_match=field_match,accepted=passed))
    failed={k for k,v in audit['checks'].items() if not v}
    reconciled=all(time_checks)
    effective_status=audit['status']
    if failed=={'frame_time_matches_diagnostics'} and reconciled:
        effective_status='field_input_consistent_time_rounding_reconciled_binary_unverified'
    fields=np.load(directory/'representative_fields.npz')
    pr=spec['parameters'];tab=np.loadtxt(next(directory.glob('*_eigenprofile.dat')),comments='#')
    width=float(pr['Grid.Positions0'].split()[-1])-float(pr['Grid.Positions0'].split()[0])
    shape=np.interp(fields['z_m'],tab[:,0],tab[:,1])
    expected_x=float(pr['Problem.SaltMoleFraction'])+float(pr['Problem.Amplitude'])*shape[:,None]*np.cos(2*np.pi*float(pr['Problem.WaveNumber'])*fields['xs'][None,:]/width)
    initial_error=float(np.max(abs(fields['x_0']-expected_x)))
    initial_checks.append(dict(id=name,checked_points=int(expected_x.size),check_scope='all vertical rows, horizontal display sample',
        profile_mode_amplitude_max_x_error=initial_error,passed=initial_error<1e-12))
    if initial_error>=1e-12:effective_status='quarantine_initial_profile_mismatch'
    fields.close()
    t=np.array(corrected)
    for row,tnew in zip(f,t):row['time_s']=str(tnew)
    for row in rr:row['time_s']=str(t[int(row['frame'])])
    xmax=val(f,'domain_x_max');xtop=val(f,'top_x_max')
    first=np.flatnonzero(xmax>=.1);cut=int(first[0]) if len(first) else len(t)
    safe=np.arange(len(t))<cut
    if not np.any(safe):raise RuntimeError('No pre-observed-threshold frame: '+name)
    sf=sig[:,0] < (t[cut] if cut<len(t) else np.inf)
    threshold_time=t[cut] if cut<len(t) else None
    original=spec['frozen_audit']
    summary=dict(id=name,case_id=spec['case_id'],run_key=spec['run_key'],status=effective_status,
        raw_remote_status=audit['status'],pvd_time_rounding_reconciled=reconciled,
        initial_profile_sample_max_x_error=initial_error,
        mode_nx=int(original['seeded_nx']),cells_x=int(original['cells_x']),cells_z=int(original['cells_z']),
        canonical_seed010=original['realizes_seed_010'],onset_class=original['onset_class'],
        onset_s=original['onset_s'],onset_low_s=original['onset_low_s'],onset_high_s=original['onset_high_s'],
        t_seed_s=original['lsa_t_seed_010_s'],lsa_sat_reference_s=original['lsa_t_sat_reference_s'],
        final_time_day=t[-1]/86400,frames=len(f),last_pre_observed_threshold_day=t[safe][-1]/86400,
        threshold_low_s=t[cut-1] if 0<cut<len(t) else (0 if cut==0 else None),
        threshold_high_s=threshold_time,post_observed_threshold_frames=len(t)-cut,
        max_top_x_before_observed_threshold=float(max(xtop[safe])),max_top_x_all_frames=float(max(xtop)),
        top_x_max_last_pre_threshold=float(xtop[safe][-1]),
        top_sigma_last_pre_threshold=float(val(f,'top_sigma_node')[safe][-1]),
        top_sigma_over_mean_last_pre_threshold=float(val(f,'top_sigma_node')[safe][-1]/val(f,'top_x_mean_node')[safe][-1]),
        top_mean_last_pre_threshold=float(val(f,'top_x_mean_area')[safe][-1]),
        max_layer_sigma_pre_threshold=float(max(val(f,'max_layer_sigma')[safe])),
        max_layer_relative_sigma_pre_threshold=float(max(val(f,'max_layer_sigma_over_mean')[safe])),
        max_layer_relative_sigma_last_pre_threshold=float(val(f,'max_layer_sigma_over_mean')[safe][-1]),
        surface_mean_decreasing_intervals_pre_threshold=int(np.sum(np.diff(val(f,'top_x_mean_area')[safe]) < -1e-10)),
        whole_salt_relative_change_pre_threshold=float(val(f,'whole_salt_m')[safe][-1]/float(f[0]['whole_salt_m'])-1),
        sigma_amplification_full_record=float(sig[-1,2]/min(sig[1:,2])),
        sigma_amplification_pre_threshold=float(sig[sf][-1,2]/min(sig[sf][1:,2])) if sf.sum()>1 else None,
        density_input_max_error_kg_m3=audit['max_density_error_kg_m3'],sw_input_max_error=audit['max_sw_error'],
        historical_executable_verified=False,background_control_available=False,discrete_budget_verified=False)
    runs.append(summary)
    for kind,values in [('surface',xtop),('global',xmax)]:
        hit=np.flatnonzero(values>=.1);k=int(hit[0]) if len(hit) else None
        events.append(dict(id=name,event=kind+'_first_observed_ge_0p1',low_s=t[k-1] if k is not None and k>0 else (0 if k==0 else None),
            high_s=t[k] if k is not None else None,source='saved VTU only',hidden_earlier_crossing_excluded=False))
    mean_hit=np.flatnonzero(sig[:,1]>=.1)
    events.append(dict(id=name,event='node_mean_first_recorded_ge_0p1',low_s=None,
        high_s=float(sig[mean_hit[0],0]) if len(mean_hit) else None,source='every-step sigmatop; upper bound on local crossing only',hidden_earlier_crossing_excluded=False))
    for region in sorted(set(r['region'] for r in rr)):
        series=[r for r in rr if r['region']==region];values=val(series,'inventory_salt_m');initial=values[0]
        selected=values[safe];tt=t[safe]
        peak=np.maximum.accumulate(selected);drawdown=(peak-selected)/initial
        idraw=int(np.argmax(drawdown));ipeak=int(np.argmax(selected[:idraw+1]))
        differences=np.diff(selected)/initial
        inv=dict(id=name,case_id=spec['case_id'],region=region,actual_z_over_H=float(series[0]['actual_z_over_H']),
            initial_inventory_m=initial,last_pre_threshold_inventory_m=selected[-1],
            relative_inventory_change=selected[-1]/initial-1,inferred_net_downward_m=initial-selected[-1],
            inferred_net_downward_over_initial=1-selected[-1]/initial,
            min_interval_change_over_initial=float(min(differences)) if len(differences) else None,
            intervals_decreasing_over_1e_minus8_initial=int(np.sum(differences < -1e-8)),
            max_peak_drawdown_over_initial=float(drawdown[idraw]),drawdown_start_s=float(tt[ipeak]),drawdown_end_s=float(tt[idraw]),
            max_net_depletion_over_initial=float(max(0,1-min(selected)/initial)),
            whole_inventory_increase_not_removed=True,independent_residual_unknown=True)
        inventory.append(inv)
        for r,v in zip(series,values):
            r.update(inferred_net_downward_m=float(initial-v),relative_inventory_change=float(v/initial-1),
                before_first_observed_global_threshold=int(float(r['time_s']) < threshold_time) if threshold_time is not None else 1)
        allframes.extend(series)
    cache[name]=dict(spec=spec,summary=summary,frame=f,regions=rr,t=t,safe=safe,sig=sig)

# Preserve trajectory identity; repeated initial states are not duplicate whole trajectories.
trajectory_groups={};frame_hash_sets={}
for name in cache:
    hashes=json.loads((EVIDENCE/'runs'/name/'field_manifest.json').read_text())
    signature=tuple((r['time_s'],r['sha256']) for r in hashes)
    trajectory_groups.setdefault(signature,[]).append(name)
    frame_hash_sets[name]={r['sha256'] for r in hashes}
duplicates=[v for v in trajectory_groups.values() if len(v)>1]

# Compare only same-case / same-mode / same-grid / same-physics runs over shared times.
overlaps=[]
names=list(cache)
for ai,a in enumerate(names):
    for b in names[ai+1:]:
        ca,cb=cache[a],cache[b];pa,pb=ca['spec']['parameters'],cb['spec']['parameters']
        if ca['spec']['case_id']!=cb['spec']['case_id']:continue
        differing=[k for k in sorted(set(pa)|set(pb)) if pa.get(k)!=pb.get(k)]
        excluded={'Problem.Name','Problem.EigenprofilePath','TimeLoop.TEnd','TimeLoop.DtInitial','TimeLoop.MaxTimeStepSize','TimeLoop.OutputInterval'}
        if any(k not in excluded for k in differing):continue
        profa=next((EVIDENCE/'runs'/a).glob('*_eigenprofile.dat')).read_bytes()
        profb=next((EVIDENCE/'runs'/b).glob('*_eigenprofile.dat')).read_bytes()
        if profa!=profb:continue
        common=[]
        for i,t0 in enumerate(ca['t']):
            j=int(np.argmin(abs(cb['t']-t0)))
            if abs(cb['t'][j]-t0) <= max(1e-5,abs(t0)*1e-10) and ca['safe'][i] and cb['safe'][j]:common.append((i,j))
        if len(common)<2:continue
        for region in sorted(set(r['region'] for r in ca['regions'])):
            va=val([r for r in ca['regions'] if r['region']==region],'inventory_salt_m')
            vb=val([r for r in cb['regions'] if r['region']==region],'inventory_salt_m')
            diff=max(abs(va[i]-vb[j])/va[0] for i,j in common)
            overlaps.append(dict(a=a,b=b,case_id=ca['spec']['case_id'],region=region,common_saved_times=len(common),
                max_inventory_difference_over_initial=diff,last_shared_time_s=ca['t'][common[-1][0]],
                changed_parameters=';'.join(differing),interpretation='same physics/mode/grid/profile, differing time settings; limited shared-window sensitivity'))

writecsv(OUT/'run_assessment.csv',runs);writecsv(OUT/'inventory_summary.csv',inventory)
writecsv(OUT/'inventory_timeseries.csv',allframes);writecsv(OUT/'event_intervals.csv',events)
writecsv(OUT/'shared_window_comparison.csv',overlaps)
writecsv(OUT/'time_reconciliation.csv',time_review)
writecsv(OUT/'initial_protocol_check.csv',initial_checks)
aggregate=dict(runs=len(runs),physical_cases=len(set(r['case_id'] for r in runs)),frames=sum(r['frames'] for r in runs),
    field_input_consistent=sum(r['status'].startswith('field_input_consistent') for r in runs),
    quarantined=sum(not r['status'].startswith('field_input_consistent') for r in runs),
    pvd_time_rounding_reconciled_runs=sum(r['pvd_time_rounding_reconciled'] for r in runs),
    observed_threshold_runs=sum(r['threshold_high_s'] is not None for r in runs),
    full_trajectory_duplicate_groups=duplicates,
    max_relative_layer_sigma_pre_threshold=max(r['max_layer_relative_sigma_pre_threshold'] for r in runs),
    max_peak_drawdown_over_initial=max(r['max_peak_drawdown_over_initial'] for r in inventory),
    inventories_with_decrease_over_1e_minus8_initial=sum(r['intervals_decreasing_over_1e_minus8_initial']>0 for r in inventory),
    simulations_started=False,transport_effectiveness_established=False)
(OUT/'summary.json').write_text(json.dumps(aggregate,indent=2))

colors=['#0072B2','#D55E00','#009E73']
region_labels={'upper_0.9':'Top 10% (target)','upper_0.75':'Top 25% (target)','upper_0.5':'Top 50% (target)'}
with PdfPages(FIG/'dense_history_diagnostics.pdf') as pdf:
    fig,axes=plt.subplots(4,4,figsize=(15,11),layout='constrained')
    for ax,(name,c) in zip(axes.flat,cache.items()):
        for color,region in zip(colors,['upper_0.9','upper_0.75','upper_0.5']):
            q=[r for r in c['regions'] if r['region']==region]
            v=val(q,'inventory_salt_m');y=100*(v/v[0]-1)
            ax.plot(c['t']/86400,y,color=color,alpha=.25)
            ax.plot(c['t'][c['safe']]/86400,y[c['safe']],color=color,label=region_labels[region])
        if c['summary']['threshold_high_s'] is not None:ax.axvline(c['summary']['threshold_high_s']/86400,color='black',ls=':',lw=.8)
        ax.set_title(f"Case {c['summary']['case_id']}, mode {c['summary']['mode_nx']} ({name[:3]})",fontsize=9);ax.axhline(0,color='.6',lw=.5)
        ax.set_xlabel('Time (day)');ax.set_ylabel('Inventory change (%)')
    for ax in axes.flat[len(cache):]:ax.set_visible(False)
    axes.flat[0].legend(fontsize=7)
    fig.suptitle('Historical upper-layer salt inventory: positive means accumulation\nPale segments follow the first observed x >= 0.1; no unperturbed reference',fontsize=13)
    fig.savefig(FIG/'inventory_overview.png');pdf.savefig(fig);plt.close(fig)
    for name,c in cache.items():
        f=c['frame'];t=c['t']/86400;safe=c['safe'];s=c['summary'];sig=c['sig'];td=sig[:,0]/86400
        fig,axes=plt.subplots(2,2,figsize=(11,7),layout='constrained')
        for color,region in zip(colors,['upper_0.9','upper_0.75','upper_0.5']):
            q=[r for r in c['regions'] if r['region']==region];v=val(q,'inventory_salt_m');y=100*(v/v[0]-1)
            axes[0,0].plot(t,y,color=color,alpha=.25);axes[0,0].plot(t[safe],y[safe],color=color,label=region_labels[region])
        axes[0,0].set_ylabel('Upper inventory change (%)');axes[0,0].legend()
        axes[0,1].plot(t,val(f,'top_x_mean_area'),label='Surface area mean')
        axes[0,1].plot(t,val(f,'top_x_max'),ls='--',label='Surface max')
        axes[0,1].axhline(.1,color='black',ls=':',label='Threshold 0.1');axes[0,1].set_ylabel('Salt mole fraction');axes[0,1].legend()
        axes[1,0].plot(td,sig[:,2]/np.min(sig[1:,2]),label='Every-step surface sigma / recorded minimum')
        axes[1,0].set_ylabel('Amplitude relative to minimum');axes[1,0].legend(fontsize=8)
        axes[1,1].plot(t,100*val(f,'max_layer_sigma_over_mean'),label='Max over depth of lateral sigma / mean')
        axes[1,1].set_ylabel('Relative lateral variation (%)');axes[1,1].legend(fontsize=8)
        for ax in axes.flat:
            ax.set_xlabel('Time (day)');ax.grid(alpha=.15)
            if s['threshold_high_s'] is not None:ax.axvspan(s['threshold_high_s']/86400,t[-1],color='.8',alpha=.4)
            if s['onset_low_s'] and s['onset_high_s']:
                ax.axvspan(float(s['onset_low_s'])/86400,float(s['onset_high_s'])/86400,color='#CC79A7',alpha=.13)
        status_label='Input / field checks consistent' if s['status'].startswith('field_input_consistent') else 'Source review required'
        fig.suptitle(f"Case {s['case_id']} | mode {s['mode_nx']} | {status_label}\nPink: inherited onset interval; grey: after first observed threshold",fontsize=11)
        pdf.savefig(fig)
        if name=='r06_c36877_n60':fig.savefig(FIG/'case36877_n60_timeline.png')
        plt.close(fig)
    # Full-resolution horizontal means show the unscaled background concentration.
    name='r06_c36877_n60';p=np.load(EVIDENCE/'runs'/name/'profiles.npz')
    fig,ax=plt.subplots(figsize=(6,5),layout='constrained')
    for k in [0,25,50,75,100]:ax.plot(p['x_mean'][k],p['z_m'],label=f"{cache[name]['t'][k]/86400:.3f} d")
    ax.set(xlabel='Horizontally averaged salt mole fraction',ylabel='Height z (m)',title='Case 36877, mode 60: vertical concentration profiles')
    ax.legend();pdf.savefig(fig);fig.savefig(FIG/'case36877_n60_profiles.png');plt.close(fig)
    fields=np.load(EVIDENCE/'runs'/name/'representative_fields.npz')
    fig,axes=plt.subplots(2,2,figsize=(11,5.5),layout='constrained')
    vmax=max(float(fields['x_0'].max()),float(fields['x_100'].max()))
    amax=max(float(abs(fields['x_anomaly_0']).max()),float(abs(fields['x_anomaly_100']).max()))*1e6
    for col,k in enumerate([0,100]):
        im=axes[0,col].pcolormesh(fields['xs'],fields['z_m'],fields['x_'+str(k)],shading='nearest',vmin=.0035,vmax=vmax,cmap='viridis',rasterized=True)
        ia=axes[1,col].pcolormesh(fields['xs'],fields['z_m'],1e6*fields['x_anomaly_'+str(k)],shading='nearest',vmin=-amax,vmax=amax,cmap='RdBu_r',rasterized=True)
        axes[0,col].set_title(f"{cache[name]['t'][k]/86400:.3f} d")
        for ax in axes[:,col]:ax.set(xlabel='Horizontal coordinate (m)',ylabel='Height (m)');ax.set_aspect('equal')
    fig.colorbar(im,ax=axes[0,:],label='Salt mole fraction',shrink=.8)
    fig.colorbar(ia,ax=axes[1,:],label='Lateral anomaly (1e-6 mol/mol)',shrink=.8)
    fig.suptitle('Case 36877, mode 60: surface enrichment and small lateral perturbations\nCommon scales within rows; display subsampled, inventories use every node',fontsize=11)
    pdf.savefig(fig);fig.savefig(FIG/'case36877_n60_fields.png');plt.close(fig)

print(json.dumps(aggregate,indent=2))
print('Shared-window comparisons:',json.dumps(overlaps,indent=2))
