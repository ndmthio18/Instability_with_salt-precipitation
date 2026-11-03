"""Strict-reference contrasts at common accepted times, never interpolated."""
import csv,gzip,io,json,math
import numpy as np
from wp8_strict_campaign import ROOT,OUT,CACHE,T,sha,save
from analyze_wp8 import rows,pairs,conservative
from single_wavelength_post import widths,savecsv,params
from wp6_paired_stopping import same_time
from wp8_campaign import RUNS

def persistent_crossing(times,values,limit,n=5):
    mask=np.asarray(values)>limit
    for i in range(len(mask)-n+1):
        if mask[i:i+n].all():return dict(status='observed_diagnostic',low_s=float(times[max(0,i-1)]),high_s=float(times[i]),confirmed_s=float(times[i+n-1]),threshold=limit,consecutive_states=n)
    return dict(status='not_observed',lower_bound_s=float(times[-1]),threshold=limit,consecutive_states=n)
def main():
    gate=json.loads((OUT/'4T/validation.json').read_text());assert gate['accepted_prefix_passed'];b=dict(np.load(CACHE/'4T/post/profiles.npz'));be=rows(OUT/'4T/surface_events.csv');bt=[float(r['time_s']) for r in be];bb=rows(OUT/'4T/transport_budget.csv');phi=float(params(OUT/'4T/params.input')['Problem.Porosity']);w=widths(b['z_m']);response=[];profile=[];transport=[];summary={};depart=[];sources={};old=ROOT/'analysis/wp8_zero_mode_20261103'
    for key in ['R24_M0','R24_M1','R48_M0','R48_M1']:
        run=ROOT/RUNS[key];cache=ROOT/'tmp/wp8_zero_mode_20261103'/key/'profiles_spectra.npz';audit=json.loads((old/key/'audit.json').read_text());assert sha(cache)==audit['files']['profiles_spectra.npz'];d=dict(np.load(cache));assert np.array_equal(d['z_m'],b['z_m']);assert abs(float(d['phi'])-phi)<1e-16
        ev=rows(run/'surface_events.csv');budget=rows(run/'transport_budget.csv');pp=pairs(bt,d['time_s']);rr=[];checkpoint=0
        for i,j in pp:
            dm=d['x_mean'][j]-b['x_mean'][i];dq=d['salt_mean'][j]-phi*b['sx_mean'][i];delta_top_max=float(ev[j]['x_top_max'])-float(be[i]['x_top_max'])
            row=dict(run=key,time_s=bt[i],time_over_T=bt[i]/T,background_surface_mean=float(b['x_mean'][i,-1]),case_surface_mean=float(d['x_mean'][j,-1]),delta_mean_top=float(dm[-1]),delta_mean_Linf=float(abs(dm).max()),delta_mean_L2=float(np.sqrt(dm**2@w/w.sum())),delta_salt_mean_Linf=float(abs(dq).max()),background_surface_max=float(be[i]['x_top_max']),case_surface_max=float(ev[j]['x_top_max']),delta_surface_max=delta_top_max,background_domain_max=float(be[i]['x_domain_max']),case_domain_max=float(ev[j]['x_domain_max']))
            rr.append(row);response.append(row)
            if i%1200==0 or (i,j)==pp[-1]:
                for z,x,y,q in zip(b['z_m'],b['x_mean'][i],d['x_mean'][j],dq):profile.append(dict(run=key,time_s=bt[i],depth_mm=(.2-z)*1000,background_mean=x,case_mean=y,delta_mean=y-x,delta_salt_mean=q))
        for limit in [1e-12,1e-10,1e-8]:depart.append(dict(run=key,observable='max_depth_absolute_mean_difference',**persistent_crossing([r['time_s'] for r in rr],[r['delta_mean_Linf'] for r in rr],limit)))
        inventories={}
        for region in sorted({r['region'] for r in budget}-{'free_domain'}):
            br=[r for r in budget if r['region']==region];bp=pairs(bt,[float(r['time_s']) for r in br]);inv=[float(br[j]['inventory_salt_m']) for i,j in bp];inventories[region]=dict(initial_m=inv[0],final_m=inv[-1],last_time_s=bt[bp[-1][0]],minimum_accepted_increment_m=min(np.diff(inv)))
        summary[key]=dict(shared_states=len(pp),last=rr[-1],max_Linf=max(r['delta_mean_Linf'] for r in rr),own_inventory_history=inventories)
        for region in sorted({r['region'] for r in bb}-{'free_domain'}):
            tr=conservative([r for r in bb if r['region']==region],[r for r in budget if r['region']==region])
            for r in tr:r.update(run=key,region=region)
            transport+=tr
        for f in ['params.input','transport_budget.csv','surface_events.csv','post/field_hashes.json']:sources[str((run/f).relative_to(ROOT))]=sha(run/f)
    savecsv(OUT/'strict_responses.csv',response);savecsv(OUT/'strict_mean_profiles.csv',profile);savecsv(OUT/'departure_diagnostics.csv',depart)
    with gzip.GzipFile(filename=str(OUT/'strict_transport.csv.gz'),mode='wb',mtime=0) as g:
        with io.TextIOWrapper(g,newline='') as f:wri=csv.DictWriter(f,fieldnames=list(transport[0]));wri.writeheader();wri.writerows(transport)
    final=[r for r in transport if same_time(r['time_s'],summary[r['run']]['last']['time_s'])];savecsv(OUT/'strict_transport_final.csv',final)
    comparison=[];prior=rows(old/'grid_checkpoints.csv')
    for i,target in enumerate([T,2*T,3*T,bt[-1]]):
        if target>bt[-1]+1e-7:continue
        for role in ['M0','M1']:
            aa=min([r for r in response if r['run']=='R24_'+role],key=lambda r:abs(r['time_s']-target));cc=min([r for r in response if r['run']=='R48_'+role],key=lambda r:abs(r['time_s']-target));assert same_time(aa['time_s'],cc['time_s'])
            means={key:next(r['delta_mean_top'] for r in response if r['run']==key and same_time(r['time_s'],aa['time_s'])) for key in ['R24_M0','R24_M1','R48_M0','R48_M1']};agap=abs((means['R24_M1']-means['R24_M0'])-(means['R48_M1']-means['R48_M0']))
            comparison.append(dict(role=role,time_s=aa['time_s'],R24_delta_mean_top=aa['delta_mean_top'],R48_delta_mean_top=cc['delta_mean_top'],strict_grid_absolute_gap=abs(aa['delta_mean_top']-cc['delta_mean_top']),seeded_unseeded_grid_absolute_gap=agap,R24_surface_max=aa['case_surface_max'],R48_surface_max=cc['case_surface_max'],background_surface=aa['background_surface_mean']))
    savecsv(OUT/'grid_comparison.csv',comparison)
    save(OUT/'results.json',dict(background_4T_completed=gate['passed'],valid_background_end_s=bt[-1],summary=summary,departure_diagnostics=depart,max_budget_identity_error_m=max(abs(r['identity_error_m']) for r in transport),source_hashes=sources,t_eff=None,crystallization_delay=None,interpretation='accepted-prefix strict-reference contrasts only; full 4T validation failed, no background extrapolation'))
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
