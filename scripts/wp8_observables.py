"""Small reproducible tables from already hash-verified WP8 reductions."""
import json
import numpy as np
from analyze_wp8 import ROOT,OUT,CACHE,RUNS,rows
from single_wavelength_post import savecsv,widths

def main():
    modes=[];checkpoints=[];cov=[];check={}
    for key in RUNS:
        d=dict(np.load(CACHE/key/'profiles_spectra.npz'));ts=rows(OUT/key/'timeseries.csv');w=widths(d['z_m']);H=w.sum()
        for i,t in enumerate(d['time_s']):
            e=d['energy'][i];total=e.sum()
            for m in [2,4,6,10,14,20]:
                if m>len(e):continue
                modes.append(dict(run=key,time_s=t,mode=m,k_rad_m=d['k_rad_m'][m-1],depth_rms_mol_mol=np.sqrt(e[m-1]/H),signed_top_amplitude_mol_mol=d['top_coeff'][i,m-1],energy_mol_mol_squared_m=e[m-1]))
            if i%1200==0 or i==len(ts)-1:
                for m,ek in enumerate(e,1):checkpoints.append(dict(run=key,frame=i,time_s=t,mode=m,k_rad_m=d['k_rad_m'][m-1],energy=ek,fraction=ek/total if total>1e-26 else None))
                cov.append(dict(run=key,frame=i,time_s=t,max_covariance=float(abs(d['covariance'][i]).max()),integrated_phi_covariance_m=float(d['phi']*d['covariance'][i]@w),max_fraction_of_salt_mean=float(np.max(abs(d['phi']*d['covariance'][i]/d['salt_mean'][i])))))
        check[key]=[ts[i] for i in range(1200,len(ts),1200)]
    savecsv(OUT/'modal_amplitudes.csv',modes);savecsv(OUT/'modal_checkpoints.csv',checkpoints);savecsv(OUT/'covariance_checkpoints.csv',cov)
    profiles=rows(OUT/'selected_mean_profiles.csv');meanenergy=[]
    for label in sorted({r['comparison'] for r in profiles}):
        rr=[r for r in profiles if r['comparison']==label]
        for t in sorted({float(r['time_s']) for r in rr}):
            a=[r for r in rr if float(r['time_s'])==t];w=widths(np.array([float(r['z_m']) for r in a]));delta=np.array([float(r['delta_mean']) for r in a]);meanenergy.append(dict(comparison=label,time_s=t,mean_difference_energy_mol_mol_squared_m=float(delta**2@w)))
    savecsv(OUT/'mean_difference_energy.csv',meanenergy)
    comparison=[];r24=rows(OUT/'R24_M1/timeseries.csv');r48=rows(OUT/'R48_M1/timeseries.csv');means=rows(OUT/'mean_grid_comparison.csv')
    for i in [1200,2400,3600,4800]:
        r=dict(frame=i,time_s=float(r24[i]['time_s']))
        for key in ['canonical_rms','nonzero_energy','k_E_rad_m','k_peak_rad_m','mean_top','max_top']:
            a,b=float(r24[i][key]),float(r48[i][key]);r['R24_'+key]=a;r['R48_'+key]=b;r['abs_grid_gap_'+key]=abs(a-b)
        r.update({k:float(v) for k,v in means[i].items() if k!='time_s'});comparison.append(r)
    savecsv(OUT/'grid_checkpoints.csv',comparison)
    (OUT/'checkpoints.json').write_text(json.dumps(check,indent=2)+'\n')
    print(json.dumps(comparison,indent=2))

if __name__=='__main__':main()
