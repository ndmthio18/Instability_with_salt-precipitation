"""Command 15: corrected 4T background versus the frozen historical strict background.

Compares only at exactly matched accepted times (32 ULP), never by interpolation and
never by step index. Byte identity is not required after the Newton criterion changed;
the differences are measured and reported.
"""
import csv,json,re,sys
import numpy as np
import wp8_full4t_campaign as f
from analyze_wp8 import rows,pairs
from single_wavelength_post import widths,params
from wp6_paired_stopping import same_time

ROOT=f.ROOT;OUT=f.OUT;CACHE=f.CACHE;STAGE=f.STAGE;sha=f.sha;save=f.save
HIST=ROOT/'analysis/wp8_strict_background_20261103'
HCACHE=ROOT/'tmp/wp8_strict_background_20261103'
T=f.T
FLUX=['cut_salt_adv_down_m_s','cut_salt_diff_down_m_s','cut_salt_total_down_m_s','cut_water_down_m_s','other_salt_out_m_s','other_water_out_m_s']
INV=['inventory_salt_m','inventory_water_m']
CUM=['cumulative_cut_salt_m','cumulative_cut_adv_m','cumulative_cut_diff_m','cumulative_net_water_out_m']
CLOSE=['salt_residual_m','water_residual_m','whole_salt_residual_m','whole_water_residual_m']


def byregion(bud):
    out={}
    for r in bud:out.setdefault(r['region'],[]).append(r)
    return out


def iterations(log):
    return [int(x) for x in re.findall(r'WP8_CONSTRAINT iterations=(\d+)',log)]


def compare():
    new=OUT/STAGE
    plan=json.loads((new/'plan.json').read_text())
    nev=rows(new/'surface_events.csv');nbud=byregion(rows(new/'transport_budget.csv'))
    nprof=dict(np.load(CACHE/STAGE/'post/profiles.npz'))
    nff=json.loads((new/'post/field_manifest.json').read_text())
    nit=iterations((new/'solver.log').read_text())
    nt=[float(r['time_s']) for r in nev]

    refs={}
    for name,d,c in [('historical_4T_prefix',HIST/'4T',HCACHE/'4T'),('historical_T_run',HIST/'T',HCACHE/'T')]:
        ev=rows(d/'surface_events.csv')
        refs[name]=dict(dir=d,ev=ev,t=[float(r['time_s']) for r in ev],bud=byregion(rows(d/'transport_budget.csv')),
                        prof=dict(np.load(c/'post/profiles.npz')),ff=json.loads((d/'post/field_manifest.json').read_text()),
                        it=iterations((d/'solver.log').read_text()),
                        stop=json.loads((d/'stop_record.json').read_text()))

    regions=sorted(nbud)
    report=dict(command='command/15.md',new_run=dict(accepted_steps=len(nev)-1,last_accepted_time_s=nt[-1],
        last_accepted_over_T=nt[-1]/T,elf_sha256=plan['elf_sha256'],criterion='scaled'),references={},
        byte_identity_required=False,interpolation_used=False,matching='exactly matched accepted times, 32 ULP')

    rowsout=[]
    for name,R in refs.items():
        pp=pairs(R['t'],nt)
        assert pp and pp[0]==(0,0),'references must share the original initial state'
        assert np.array_equal(R['prof']['z_m'],nprof['z_m']),'vertical nodes differ'
        idx={i:j for i,j in pp}
        env=dict(shared_accepted_states=len(pp),first_shared_time_s=R['t'][pp[0][0]],last_shared_time_s=R['t'][pp[-1][0]],
                 last_shared_over_T=R['t'][pp[-1][0]]/T,
                 reference_accepted_steps=R['stop']['accepted_steps'],reference_stop_reason=R['stop']['reason'])
        env['max_abs_dx']=max(float(abs(R['prof']['x_mean'][i]-nprof['x_mean'][j]).max()) for i,j in pp)
        env['max_abs_dp_Pa']=max(float(abs(R['prof']['pw_mean'][i]-nprof['pw_mean'][j]).max()) for i,j in pp)
        env['max_abs_dsw']=max(float(abs(R['prof']['sw_mean'][i]-nprof['sw_mean'][j]).max()) for i,j in pp)
        env['max_rel_dx']=max(float((abs(R['prof']['x_mean'][i]-nprof['x_mean'][j])/abs(R['prof']['x_mean'][i])).max()) for i,j in pp)
        env['identical_field_sha256']=sum(R['ff'][i]['sha256']==nff[j]['sha256'] for i,j in pp)
        env['identical_field_fraction']=env['identical_field_sha256']/len(pp)
        for key in INV+FLUX+CUM+CLOSE:
            env['max_abs_d_'+key]=max(abs(float(R['bud'][g][i][key])-float(nbud[g][j][key])) for g in regions for i,j in pp)
        env['reference_iterations_max']=max(R['it'])+1;env['new_iterations_max']=max(nit[:len(pp)-1])+1
        env['reference_iterations_mean']=float(np.mean(R['it'])+1);env['new_iterations_mean']=float(np.mean(nit[:len(pp)-1])+1)
        report['references'][name]=env

        # --- requested checkpoints -------------------------------------
        cps=[('T',1200),('2T',2400),('3T',3600),('historical_last_accepted_3.6775T',R['stop']['accepted_steps']),('4T',4800)]
        for label,step in cps:
            if step>=len(nt) or step not in idx.values():
                hit=[(i,j) for i,j in pp if j==step]
            else:
                hit=[(i,j) for i,j in pp if j==step]
            if not hit:
                rowsout.append(dict(reference=name,checkpoint=label,new_step=step,
                    status='no matched historical accepted state',time_s=nt[step] if step<len(nt) else None,
                    region='',quantity='',reference_value='',new_value='',difference='',relative_difference=''))
                continue
            i,j=hit[0]
            rowsout.append(dict(reference=name,checkpoint=label,new_step=j,status='matched',time_s=nt[j],
                region='profile',quantity='max_abs_dx',reference_value='',new_value='',
                difference=float(abs(R['prof']['x_mean'][i]-nprof['x_mean'][j]).max()),relative_difference=''))
            rowsout.append(dict(reference=name,checkpoint=label,new_step=j,status='matched',time_s=nt[j],
                region='profile',quantity='max_abs_dp_Pa',reference_value='',new_value='',
                difference=float(abs(R['prof']['pw_mean'][i]-nprof['pw_mean'][j]).max()),relative_difference=''))
            rowsout.append(dict(reference=name,checkpoint=label,new_step=j,status='matched',time_s=nt[j],
                region='surface',quantity='x_top_area_mean',reference_value=float(R['ev'][i]['x_top_area_mean']),
                new_value=float(nev[j]['x_top_area_mean']),
                difference=float(nev[j]['x_top_area_mean'])-float(R['ev'][i]['x_top_area_mean']),
                relative_difference=(float(nev[j]['x_top_area_mean'])-float(R['ev'][i]['x_top_area_mean']))/float(R['ev'][i]['x_top_area_mean'])))
            for g in regions:
                for key in INV+FLUX+CUM+CLOSE:
                    a=float(R['bud'][g][i][key]);b=float(nbud[g][j][key])
                    rowsout.append(dict(reference=name,checkpoint=label,new_step=j,status='matched',time_s=nt[j],
                        region=g,quantity=key,reference_value=a,new_value=b,difference=b-a,
                        relative_difference=(b-a)/a if a!=0 else ''))

    # --- unmatched tail: the part no historical state covers -------------
    hlast=refs['historical_4T_prefix']['stop']['accepted_steps']
    report['beyond_historical'] = dict(
        first_step_without_historical_state=hlast+1,
        steps_beyond_history=len(nev)-1-hlast,
        time_span_s=nt[-1]-nt[hlast],time_span_over_T=(nt[-1]-nt[hlast])/T,
        new_surface_mean_at_history_end=float(nev[hlast]['x_top_area_mean']),
        new_surface_mean_at_4T=float(nev[-1]['x_top_area_mean']),
        new_domain_max_at_4T=float(nev[-1]['x_domain_max']))

    # --- old failure window ----------------------------------------------
    failed=rows(HIST/'failed_step_iterations.csv')
    ndi=rows(new/'convergence_v2.csv');nacc={int(r['step']):r for r in ndi if int(r['accepted'])}
    w=nacc[hlast+1]
    report['old_failure_window']=dict(
        historical_failed_attempt_step=hlast+1,historical_failed_iterations=len(failed),
        historical_failed_min_salt_residual=min(float(r['salt_residual']) for r in failed),
        historical_fixed_absolute_salt_threshold=2e-20,
        corrected_accepted=True,corrected_iterations=int(w['iteration'])+1,
        corrected_salt_max=float(w['salt_max']),corrected_water_max=float(w['water_max']),
        corrected_scaled_nonlinear=float(w['scaled_nonlinear']),
        corrected_relative_shift=float(w['relative_shift']),
        corrected_linear_eta_component=float(w['linear_eta_component']),
        corrected_regional_salt_error_m=float(w['regional_salt_error_m']),
        passed=True)

    # --- Newton iteration distribution, whole runs ------------------------
    hv,hc=np.unique(np.array(refs['historical_4T_prefix']['it'])+1,return_counts=True)
    nv,nc=np.unique(np.array(nit)+1,return_counts=True)
    report['newton_iteration_counts']=dict(
        historical_prefix={str(int(k)):int(v) for k,v in zip(hv,hc)},
        corrected_full_4T={str(int(k)):int(v) for k,v in zip(nv,nc)},
        historical_mean=float(np.mean(refs['historical_4T_prefix']['it'])+1),corrected_mean=float(np.mean(nit)+1),
        historical_max=int(max(refs['historical_4T_prefix']['it'])+1),corrected_max=int(max(nit)+1))

    with (OUT/'checkpoint_comparison.csv').open('w',newline='') as fh:
        wtr=csv.DictWriter(fh,fieldnames=list(rowsout[0]));wtr.writeheader();wtr.writerows(rowsout)
    report['checkpoint_comparison_rows']=len(rowsout)
    save(OUT/'historical_comparison.json',report)
    print(json.dumps(report,indent=2,default=str))
    return report


if __name__=='__main__':compare()
