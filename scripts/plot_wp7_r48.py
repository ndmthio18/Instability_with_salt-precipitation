"""Publication-style WP7 R24/R48 evidence figures."""
import csv,json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'analysis/wp7_shortwave_validation_20261103';R24=ROOT/'analysis/wp6_4T_exploratory_20261103';T=136423.0637723484
def rows(p):
    with p.open() as f:return list(csv.DictReader(f))
def nums(r,k):return np.array([float(x[k]) for x in r])
def optional_nums(r,k):return np.array([float(x[k]) if x[k] else np.nan for x in r])
def save(fig,name):fig.tight_layout();fig.savefig(OUT/(name+'.png'),dpi=240);fig.savefig(OUT/(name+'.pdf'));plt.close(fig)
def main():
    colors={'R24':'#0072B2','R48':'#D55E00'}
    fig,ax=plt.subplots(1,3,figsize=(12,3.5))
    for grid,path in [('R24',R24/'runs/WP6-M1-explore4T/post/weighted_diagnostics.csv'),('R48',OUT/'runs/WP7-R48-M1-4T/post/weighted_diagnostics.csv')]:
        r=rows(path);t=nums(r,'time_s')/T
        for a,key,title in zip(ax,['weighted_top_sigma','weighted_fourier','weighted_projection'],['surface sigma','canonical Fourier','A28 projection']):a.semilogy(t,np.abs(nums(r,key)),color=colors[grid],label=grid);a.set(xlabel='t / T',title=title)
    ax[0].set_ylabel('amplitude (mol/mol)');ax[0].legend();fig.suptitle('A  R24-R48 growth at common physical protocol');save(fig,'Figure_A_R24_R48_growth')
    modes=rows(OUT/'final_mode_energy_R24_R48.csv');m=[r for r in modes if r['role']=='M1'];x=nums(m,'mode').astype(int);fig,ax=plt.subplots(figsize=(7,4));ax.bar(x-.18,optional_nums(m,'R24_energy_fraction'),.36,label='R24',color=colors['R24']);ax.bar(x+.18,optional_nums(m,'R48_energy_fraction'),.36,label='R48',color=colors['R48']);ax.axvline(12.5,color='.5',ls=':',label='R48 falls below 8 cells/wavelength');ax.set(xlabel='Neumann mode m (physical k=m pi/W)',ylabel='final volume anomaly energy fraction',yscale='log',ylim=(1e-8,1));ax.legend();fig.suptitle('B  Final physical-mode energy');save(fig,'Figure_B_R24_R48_modes')
    fig,ax=plt.subplots(1,2,figsize=(9,3.5))
    for grid,path in [('R24',R24/'runs/WP6-M0-explore4T/post/weighted_diagnostics.csv'),('R48',OUT/'runs/WP7-R48-M0-4T/post/weighted_diagnostics.csv')]:
        r=rows(path);t=nums(r,'time_s')/T;odd=optional_nums(r,'odd_neumann_energy_fraction');ax[0].semilogy(t,nums(r,'weighted_top_sigma'),color=colors[grid],label=grid);ax[1].semilogy(t,np.where(np.isnan(odd),np.nan,np.maximum(odd,1e-30)),color=colors[grid],label=grid)
    ax[0].set(xlabel='t / T',ylabel='M0 surface sigma',title='background growth');ax[1].set(xlabel='t / T',ylabel='odd-mode energy fraction',title='symmetry diagnostic');ax[0].legend();fig.suptitle('C  Two-dimensional M0 is not assumed one-dimensional');save(fig,'Figure_C_R24_R48_background')
    tr=rows(OUT/'transport_R24_R48_final.csv');depth=1-np.array([float(r['region'].split('_')[1]) for r in tr]);order=np.argsort(depth);depth=depth[order]*0.2*1000;fig,ax=plt.subplots(1,2,figsize=(9,3.5))
    for grid in ['R24','R48']:
        ax[0].plot(np.array([float(r[grid+'_delta_C_m']) for r in tr])[order],depth,'o-',label=grid,color=colors[grid]);ax[1].plot(np.array([float(r[grid+'_M1_total_down_m_s']) for r in tr])[order],depth,'o-',label=grid,color=colors[grid])
    ax[0].set(xlabel='M1-M0 cumulative extra salt, Delta C (m)',ylabel='depth below surface (mm)');ax[1].set(xlabel='M1 signed total salt flux (m/s)',ylabel='depth below surface (mm)');ax[0].invert_yaxis();ax[1].invert_yaxis();ax[0].legend();fig.suptitle('D  Seven-plane transport at the last common state');save(fig,'Figure_D_R24_R48_transport')
if __name__=='__main__':main()
