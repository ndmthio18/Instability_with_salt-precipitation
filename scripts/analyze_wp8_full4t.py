"""Command 15 scientific analysis: R24/R48 M0/M1 against the corrected strict background.

Background validation (compare_wp8_full4t.py) and grid convergence are kept separate.
No t_eff, crystallization delay or surrogate label is assigned anywhere in this module.
"""
import csv,gzip,io,json,sys
from pathlib import Path
import numpy as np
import wp8_full4t_campaign as f
from analyze_wp8 import rows,pairs,conservative
from single_wavelength_post import widths,params,savecsv
from wp6_paired_stopping import same_time
from wp8_campaign import RUNS

ROOT=f.ROOT;OUT=f.OUT;CACHE=f.CACHE;STAGE=f.STAGE;sha=f.sha;save=f.save
ZERO=ROOT/'analysis/wp8_zero_mode_20261103'
ZCACHE=ROOT/'tmp/wp8_zero_mode_20261103'
T=f.T
TWO=['R24_M0','R24_M1','R48_M0','R48_M1']
CHECK=[('T',1200),('2T',2400),('3T',3600),('3.6775T',4413),('4T',4800)]


def background():
    d=OUT/STAGE
    v=json.loads((OUT/'validation.json').read_text())
    assert v['passed'],'the corrected background must pass its own gates first'
    a=dict(np.load(CACHE/STAGE/'post/profiles.npz'))
    ev=rows(d/'surface_events.csv')
    par=params(d/'params.input');phi=float(par['Problem.Porosity'])
    return dict(time_s=np.array([float(r['time_s']) for r in ev]),z_m=a['z_m'],x_mean=a['x_mean'],
                sw_mean=a['sw_mean'],salt_mean=phi*a['sx_mean'],covariance=np.zeros_like(a['x_mean']),
                phi=phi,ev=ev,bud=rows(d/'transport_budget.csv'),nx=1)


def main():
    bg=background()
    prod={k:dict(np.load(ZCACHE/k/'profiles_spectra.npz')) for k in TWO}
    bud={k:rows(ROOT/RUNS[k]/'transport_budget.csv') for k in TWO}
    ev2={k:rows(ROOT/RUNS[k]/'surface_events.csv') for k in TWO}
    wz=widths(bg['z_m'])

    response=[];profiles=[];transport=[];summary={};surface=[];events=[]

    for key in TWO:
        p=prod[key];assert np.array_equal(p['z_m'],bg['z_m']),'vertical nodes differ'
        pp=pairs(bg['time_s'],p['time_s'])
        name=key.split('_')[0]+('_B' if key.endswith('M0') else '_C')
        rr=[]
        for i,j in pp:
            delta=p['x_mean'][j]-bg['x_mean'][i]
            ds=p['salt_mean'][j]-bg['salt_mean'][i]
            r=dict(comparison=name,case=key,time_s=float(bg['time_s'][i]),t_over_T=float(bg['time_s'][i]/T),
                   mean_response_top=float(delta[-1]),mean_response_Linf=float(abs(delta).max()),
                   mean_response_L2=float(np.sqrt(delta**2@wz/wz.sum())),
                   salt_response_Linf=float(abs(ds).max()),
                   integrated_covariance_difference_m=float((p['phi']*p['covariance'][j]-bg['phi']*bg['covariance'][i])@wz))
            rr.append(r);response.append(r)
            if i%600==0 or (i,j)==pp[-1]:
                for z,xr,xc,s in zip(bg['z_m'],bg['x_mean'][i],p['x_mean'][j],ds):
                    profiles.append(dict(comparison=name,case=key,time_s=r['time_s'],z_m=float(z),
                        depth_mm=float((.2-z)*1000),background_x=float(xr),case_x=float(xc),
                        delta_mean=float(xc-xr),delta_salt_mean=float(s)))
        summary[name]=dict(case=key,shared_states=len(pp),end_s=rr[-1]['time_s'],end_over_T=rr[-1]['t_over_T'],
            final=rr[-1],max_mean_response_Linf=max(x['mean_response_Linf'] for x in rr),
            checkpoints={lab:rr[s] for lab,s in CHECK if s<len(rr)})
        for region in sorted({r['region'] for r in bg['bud']}-{'free_domain'}):
            tr=conservative([r for r in bg['bud'] if r['region']==region],[r for r in bud[key] if r['region']==region])
            for r in tr:r.update(comparison=name,case=key,region=region,t_over_T=r['time_s']/T)
            transport.extend(tr)

    # --- surface area mean versus local maximum ---------------------------
    for label,ev in [('strict_background',bg['ev'])]+[(k,ev2[k]) for k in TWO]:
        for lab,s in CHECK:
            if s>=len(ev):continue
            r=ev[s];m=float(r['x_top_area_mean']);mx=float(r['x_top_max'])
            surface.append(dict(run=label,checkpoint=lab,step=s,time_s=float(r['time_s']),
                t_over_T=float(r['time_s'])/T,x_top_area_mean=m,x_top_max=mx,
                local_excess=mx-m,local_excess_relative=(mx-m)/m,
                x_domain_max=float(r['x_domain_max']),
                fraction_of_solubility=mx/0.1))
        cross=[r for r in ev if int(r['surface_threshold_reached']) or int(r['global_threshold_reached'])]
        events.append(dict(run=label,accepted_steps=len(ev)-1,end_s=float(ev[-1]['time_s']),
            end_over_T=float(ev[-1]['time_s'])/T,max_x_top_max=max(float(r['x_top_max']) for r in ev),
            max_x_domain_max=max(float(r['x_domain_max']) for r in ev),
            threshold=0.1,threshold_observed=bool(cross),
            first_crossing_time_s=float(cross[0]['time_s']) if cross else None,
            status='no local solubility-threshold event within the covered horizon' if not cross else 'crossing observed'))

    # --- grid sensitivity of the inferred response, kept separate ---------
    grid=[]
    for typ in ['B','C']:
        a=[r for r in response if r['comparison']=='R24_'+typ];b=[r for r in response if r['comparison']=='R48_'+typ]
        for x,y in zip(a,b):
            assert same_time(x['time_s'],y['time_s'])
            grid.append(dict(comparison_type=typ,time_s=x['time_s'],t_over_T=x['t_over_T'],
                R24_mean_response_top=x['mean_response_top'],R48_mean_response_top=y['mean_response_top'],
                absolute_grid_difference=abs(x['mean_response_top']-y['mean_response_top']),
                grid_ratio=x['mean_response_top']/y['mean_response_top'] if y['mean_response_top'] else '',
                R24_Linf=x['mean_response_Linf'],R48_Linf=y['mean_response_Linf']))

    savecsv(OUT/'background_referenced_responses.csv',response)
    savecsv(OUT/'background_referenced_profiles.csv',profiles)
    savecsv(OUT/'surface_mean_versus_local_max.csv',surface)
    savecsv(OUT/'grid_sensitivity_vs_background.csv',grid)
    with gzip.GzipFile(filename=str(OUT/'seven_plane_redistribution.csv.gz'),mode='wb',mtime=0) as g:
        with io.TextIOWrapper(g,newline='') as fh:
            w=csv.DictWriter(fh,fieldnames=list(transport[0]));w.writeheader();w.writerows(transport)

    # seven-plane initial-inventory-corrected redistribution at the checkpoints
    red={}
    for name in ['R24_B','R24_C','R48_B','R48_C']:
        sub=[r for r in transport if r['comparison']==name]
        red[name]={}
        for lab,s in CHECK:
            t=bg['time_s'][s]
            red[name][lab]={r['region']:dict(time_s=r['time_s'],t_over_T=r['t_over_T'],
                delta_I_corrected_m=r['delta_I_corrected_m'],delta_C_m=r['delta_C_m'],
                delta_I_over_initial=r['delta_I_over_initial'],delta_C_over_initial=r['delta_C_over_initial'],
                identity_error_m=r['identity_error_m'])
                for r in sub if same_time(r['time_s'],t)}

    prev=json.loads((ZERO/'results.json').read_text())
    out=dict(command='command/15.md',
        background=dict(run=str((OUT/STAGE).relative_to(ROOT)),accepted_steps=len(bg['ev'])-1,
            end_over_T=float(bg['time_s'][-1]/T),criterion='scaled',nx=1,nz=len(bg['z_m'])-1,
            note='corrected strict horizontally uniform background, validated separately'),
        previous_state=dict(background_after_T=prev['background_after_T'],
            previously_shared_states_B_and_C=prev['comparisons']['R24_B']['shared_states']),
        comparisons=summary,
        seven_plane_redistribution=red,
        surface_mean_versus_local_max=surface,
        solubility_threshold_events=events,
        grid_sensitivity=dict(scope='difference between R24- and R48-inferred background-referenced response; NOT a convergence proof',
            final={typ:dict(time_s=grid[-1]['time_s'],
                R24=[r for r in grid if r['comparison_type']==typ][-1]['R24_mean_response_top'],
                R48=[r for r in grid if r['comparison_type']==typ][-1]['R48_mean_response_top'],
                absolute_difference=[r for r in grid if r['comparison_type']==typ][-1]['absolute_grid_difference'],
                ratio=[r for r in grid if r['comparison_type']==typ][-1]['grid_ratio']) for typ in ['B','C']},
            max_absolute_difference={typ:max(r['absolute_grid_difference'] for r in grid if r['comparison_type']==typ) for typ in ['B','C']}),
        max_identity_error_m=max(abs(r['identity_error_m']) for r in transport),
        solver_criterion_caveat=dict(
            background_elf=json.loads((OUT/STAGE/'plan.json').read_text())['elf_sha256'],
            two_d_runs_use_earlier_elf=True,
            bound_from_background_to_background_difference_max_abs_dx=json.loads((OUT/'historical_comparison.json').read_text())['references']['historical_4T_prefix']['max_abs_dx'],
            statement='the Newton criterion change alters accepted states at the 1e-15 level; responses below that level are not resolved'),
        t_eff=None,crystallization_delay=None,surrogate_label=None,
        interpretation_limits=['background validation is reported separately from grid convergence',
            'a difference exceeding round-off is not by itself validated physical convection',
            'no local solubility threshold was reached, so no precipitation-kinetics statement is made'])
    save(OUT/'scientific_analysis.json',out)
    print(json.dumps({k:v for k,v in out.items() if k not in ['surface_mean_versus_local_max','seven_plane_redistribution']},indent=2,default=str))
    return out


if __name__=='__main__':main()
