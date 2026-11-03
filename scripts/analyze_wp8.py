"""WP8: zero-mode, conservative covariance, exact accepted-state comparisons."""
import csv,json,hashlib,gzip,io
from pathlib import Path
import numpy as np
from wp8_campaign import ROOT,OUT,CACHE,RUNS
from single_wavelength_post import widths,savecsv
from wp6_paired_stopping import same_time

def rows(p):
    opener=gzip.open if p.suffix=='.gz' else open
    with opener(p,'rt',newline='') as f:return list(csv.DictReader(f))
def pairs(a,b):
    i=j=0;out=[]
    while i<len(a) and j<len(b):
        if same_time(float(a[i]),float(b[j])):out.append((i,j));i+=1;j+=1
        elif a[i]<b[j]:i+=1
        else:j+=1
    return out
def conservative(ref,case):
    p=pairs([float(r['time_s']) for r in ref],[float(r['time_s']) for r in case]);out=[];initial=float(case[0]['inventory_salt_m'])-float(ref[0]['inventory_salt_m'])
    for i,j in p:
        a,b=ref[i],case[j];assert abs(float(a['actual_z_over_H'])-float(b['actual_z_over_H']))<1e-12
        di=float(b['inventory_salt_m'])-float(a['inventory_salt_m'])-initial;dc=float(b['cumulative_cut_salt_m'])-float(a['cumulative_cut_salt_m']);dr=float(b['salt_residual_m'])-float(a['salt_residual_m'])
        r=dict(time_s=float(a['time_s']),actual_z_over_H=float(a['actual_z_over_H']),initial_delta_I_m=initial,delta_I_corrected_m=di,delta_C_m=dc,delta_R_m=dr,identity_error_m=di+dc-dr,delta_I_over_initial=di/float(ref[0]['inventory_salt_m']),delta_I_over_current=di/float(a['inventory_salt_m']),delta_C_over_initial=dc/float(ref[0]['inventory_salt_m']),delta_C_over_current=dc/float(a['inventory_salt_m']))
        for name,key in [('I_m','inventory_salt_m'),('C_m','cumulative_cut_salt_m'),('adv_m_s','cut_salt_adv_down_m_s'),('diff_m_s','cut_salt_diff_down_m_s'),('total_m_s','cut_salt_total_down_m_s'),('R_m','salt_residual_m')]:r['ref_'+name]=float(a[key]);r['case_'+name]=float(b[key]);r['delta_'+name]=float(b[key])-float(a[key])
        assert abs(r['identity_error_m'])<1e-18
        out.append(r)
    return out
def main():
    products={k:dict(np.load(CACHE/k/'profiles_spectra.npz')) for k in RUNS};budget={k:rows(ROOT/v/'transport_budget.csv') for k,v in RUNS.items()};sources={};audit={};response=[];profiles=[];transport=[];summary={};gridrows=[]
    for key,rel in RUNS.items():
        run=ROOT/rel;d=products[key];a=json.loads((OUT/key/'audit.json').read_text());assert a['all_field_hashes_verified'];assert hashlib.sha256((CACHE/key/'profiles_spectra.npz').read_bytes()).hexdigest()==a['files']['profiles_spectra.npz']
        for f in ['params.input','finished.json','transport_budget.csv','transport_metadata.json','surface_events.csv','post/field_hashes.json']:sources[str(Path(rel)/f)] = hashlib.sha256((run/f).read_bytes()).hexdigest()
        audit[key]=dict(frames=len(d['time_s']),end_s=float(d['time_s'][-1]),nx=len(d['xs'])-1,nz=len(d['z_m'])-1)
    for grid in ['R24','R48']:
        for typ,ref,case in [('A',grid+'_M0',grid+'_M1'),('B','B',grid+'_M0'),('C','B',grid+'_M1')]:
            name=grid+'_'+typ;a,b=products[ref],products[case];pp=pairs(a['time_s'],b['time_s']);assert np.array_equal(a['z_m'],b['z_m']);wz=widths(a['z_m']);rr=[]
            for i,j in pp:
                delta=b['x_mean'][j]-a['x_mean'][i];ds=b['salt_mean'][j]-a['salt_mean'][i]
                r=dict(comparison=name,time_s=float(a['time_s'][i]),mean_response_top=float(delta[-1]),mean_response_Linf=float(abs(delta).max()),mean_response_L2=float(np.sqrt(delta**2@wz/wz.sum())),salt_response_Linf=float(abs(ds).max()),integrated_covariance_difference_m=float((b['phi']*b['covariance'][j]-a['phi']*a['covariance'][i])@wz))
                rr.append(r);response.append(r)
                if i%1200==0 or (i,j)==pp[-1]:
                    for z,x,y,s in zip(a['z_m'],a['x_mean'][i],b['x_mean'][j],ds):profiles.append(dict(comparison=name,time_s=r['time_s'],z_m=float(z),depth_mm=float((.2-z)*1000),ref_x=float(x),case_x=float(y),delta_mean=float(y-x),delta_salt_mean=float(s)))
            summary[name]=dict(shared_states=len(pp),end_s=rr[-1]['time_s'],final=rr[-1],max_mean_response_Linf=max(r['mean_response_Linf'] for r in rr))
            for region in sorted({r['region'] for r in budget[ref]}-{'free_domain'}):
                tr=conservative([r for r in budget[ref] if r['region']==region],[r for r in budget[case] if r['region']==region])
                for r in tr:r.update(comparison=name,region=region)
                transport.extend(tr)
    savecsv(OUT/'mean_responses.csv',response);savecsv(OUT/'selected_mean_profiles.csv',profiles)
    with gzip.GzipFile(filename=str(OUT/'conservative_comparisons.csv.gz'),mode='wb',mtime=0) as g:
        with io.TextIOWrapper(g,newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(transport[0]));w.writeheader();w.writerows(transport)
    for a,b in zip([r for r in response if r['comparison']=='R24_A'],[r for r in response if r['comparison']=='R48_A']):
        assert same_time(a['time_s'],b['time_s']);gridrows.append(dict(time_s=a['time_s'],R24_mean_response_top=a['mean_response_top'],R48_mean_response_top=b['mean_response_top'],absolute_grid_difference=abs(a['mean_response_top']-b['mean_response_top']),R24_Linf=a['mean_response_Linf'],R48_Linf=b['mean_response_Linf']))
    savecsv(OUT/'mean_grid_comparison.csv',gridrows)
    final=[r for r in transport if same_time(r['time_s'],summary[r['comparison']]['end_s'])];savecsv(OUT/'transport_final.csv',final)
    audit['background_gate']=json.loads((ROOT/'analysis/wp7_shortwave_validation_20261103/background_validation.json').read_text())
    result=dict(audit=audit,comparisons=summary,max_identity_error_m=max(abs(r['identity_error_m']) for r in transport),t_eff=None,background_after_T='unavailable',source_hashes=sources)
    (OUT/'results.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['source_hashes','audit']},indent=2))
if __name__=='__main__':main()
