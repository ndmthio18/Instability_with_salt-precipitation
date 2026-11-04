"""WP9 Task D: LSA wavenumber audit for case 36877.

Uses the frozen Paper 1 solver `code/lsa_unsaturated.py` as a read-only reference; no
new dispersion relation is introduced. The quasi-steady growth rate is the one already
implemented and self-tested there (storage term restored in the chi-equation of paper
eq. 3.67, base state frozen at t). Its assumptions are restated in the output and its
stored anchors are reproduced before any new number is reported.

Run: python3 wp9_lsa_audit.py        (numpy only, seconds)
"""
import json,math,sys,time
from pathlib import Path
import numpy as np
import wp9_common as c

sys.path.insert(0,str(c.PAPER1))
import lsa_unsaturated as U                                    # noqa: E402

N_CHEB=96
CHI_MODE='paper_eq367'
NX_SET_PROTOTYPE=(1,2,3,5,8,12,20,30)       # the 2026-09-05 scan that produced the store
WIDE_NX=sorted({*range(1,31),*range(32,321,2)})
BOX_M=list(range(1,25))                      # box modes R24 can represent at all
CHECKPOINTS=[('t_c',1.0),('T',None),('2T',None),('3T',None),('3.6775T',None),('4T',None)]


def _state(t_hat,point):
    dt=t_hat/100.0
    z,Dz,pp=U._assemble(t_hat,point['n'],point['m'],point['L'],point['Swr'],
                        point['delta'],point['pbot'],point['beta'],N_CHEB,dt,0)
    S0=U.FSp(point['pbot']+z,point['delta'],point['n'],point['m'],point['Swr'])
    return Dz,pp,S0


def _sigma_and_rc(Dz,pp,S0,point,R_E,nx):
    """Quasi-steady sigma (1/t_hat) and neutral R_c at one quantized wavenumber."""
    B,Cm,D,E,F=U.build_matrices(Dz,pp,point['beta'],c.a_hat(nx),chi_coeff_mode=CHI_MODE)
    G=np.linalg.solve(B,Cm)
    M=np.diag(S0);M[0,:]=0.0;M[-1,:]=0.0
    with np.errstate(divide='ignore',invalid='ignore'):
        s=1.0/np.linalg.eigvals(np.linalg.solve(D-R_E*(E@G)-F@G,M))
        rc=1.0/np.linalg.eigvals(np.linalg.solve(D-F@G,E@G))
    ok=np.isfinite(s)&(np.abs(s.imag)<1e-8*np.abs(s.real)+1e-12)
    sigma=float(np.max(s.real[ok])) if ok.any() else float('nan')
    return sigma,float(U._select_R(rc,'real_positive'))


def scan(t_hat,point,R_E,nx_values):
    Dz,pp,S0=_state(t_hat,point)
    return {nx:_sigma_and_rc(Dz,pp,S0,point,R_E,nx) for nx in nx_values}


def main():
    start=time.time()
    point,R_E=c.lsa_point()
    case,case_sha=c.case_record()
    obs,obs_sha=c.observation_record()
    t_hat_c=float(obs['t_hat_onset'])
    t_c_s=float(obs['t_onset_s'])
    ref=float(case['reference_time_s'])
    assert abs(t_hat_c*ref-t_c_s)<1e-6*t_c_s,'stored dimensionless and dimensional onset disagree'

    # ---- reproduce the stored anchors before reporting anything new -----------------
    anchors={}
    got=scan(2.0*t_hat_c,point,R_E,sorted({*NX_SET_PROTOTYPE,int(float(obs['critical_nx']))}))
    sig_c=got[int(float(obs['critical_nx']))][0]*t_hat_c
    best_proto=max(NX_SET_PROTOTYPE,key=lambda k:got[k][0])
    anchors['sigma_c_tc_at2p0']=dict(recomputed=sig_c,stored=float(obs['sigma_c_tc_at2p0']),
                                     relative_difference=abs(sig_c/float(obs['sigma_c_tc_at2p0'])-1))
    anchors['sigma_max_tc_at2p0']=dict(recomputed=got[best_proto][0]*t_hat_c,
                                       stored=float(obs['sigma_max_tc_at2p0']),
                                       relative_difference=abs(got[best_proto][0]*t_hat_c/float(obs['sigma_max_tc_at2p0'])-1))
    wide2=scan(2.0*t_hat_c,point,R_E,WIDE_NX)
    wide2f={k:v[0] for k,v in wide2.items() if np.isfinite(v[0])}
    wide2_peak=max(wide2f,key=wide2f.get)
    anchors['nx_max_at2p0']=dict(recomputed=best_proto,stored=int(float(obs['nx_max_at2p0'])),
                                 scan_set=list(NX_SET_PROTOTYPE),
                                 wide_scan_peak_nx=wide2_peak,
                                 wide_scan_peak_sigma_times_t_hat_c=wide2f[wide2_peak]*t_hat_c,
                                 wide_scan_peak_is_interior=wide2_peak not in (min(WIDE_NX),max(WIDE_NX)),
                                 note=('the stored value is also the largest entry of the 2026-09-05 scan set '
                                       '(1,2,3,5,8,12,20,30); the wide scan confirms it is nevertheless the '
                                       'true interior maximum at 2 t_c, which is a coincidence that does not '
                                       'hold at later times'))
    rc_at_tc,nx_at_tc=U.critical_rayleigh(t_hat_c,point['n'],point['m'],point['L'],point['Swr'],
                                          point['delta'],point['pbot'],point['beta'],
                                          nx_max=30,ell_hat=c.ELL_HAT,N=N_CHEB,
                                          chi_coeff_mode=CHI_MODE)
    anchors['critical_rayleigh_at_tc']=dict(recomputed=rc_at_tc,target_R_E=R_E,
                                            relative_difference=abs(rc_at_tc/R_E-1),
                                            recomputed_critical_nx=nx_at_tc,
                                            stored_critical_nx=int(float(obs['critical_nx'])))
    anchors_ok=(anchors['sigma_c_tc_at2p0']['relative_difference']<1e-9
                and anchors['sigma_max_tc_at2p0']['relative_difference']<1e-9
                and anchors['critical_rayleigh_at_tc']['relative_difference']<2e-6
                and anchors['critical_rayleigh_at_tc']['recomputed_critical_nx']==anchors['critical_rayleigh_at_tc']['stored_critical_nx'])
    assert anchors_ok,'stored LSA anchors not reproduced; refusing to report new numbers'

    # ---- the registered times -------------------------------------------------------
    times=[('t_c',t_c_s),('2t_c',2.0*t_c_s)]
    for name,factor in [('T',1.0),('2T',2.0),('3T',3.0),('3.6775T',501695.8170227767/c.T_S),('4T',4.0)]:
        times.append((name,factor*c.T_S))
    t_amp_s=float(obs['t_hat_amp'])*ref
    times.append(('t_amp',t_amp_s))

    # ---- wide scan and box-representable scan ---------------------------------------
    wide=[];box=[]
    for name,t_s in times:
        t_hat=t_s/ref
        res=scan(t_hat,point,R_E,sorted({*WIDE_NX,*[int(c.lsa_nx(m)) for m in BOX_M]}))
        finite={k:v[0] for k,v in res.items() if np.isfinite(v[0])}
        peak=max(finite,key=finite.get) if finite else None
        for nx in sorted(res):
            sigma,rc=res[nx]
            m=nx/c.NX_PER_M
            wide.append(dict(checkpoint=name,time_s=t_s,t_hat=t_hat,t_over_tc=t_hat/t_hat_c,
                             nx=nx,a_hat=c.a_hat(nx),k_phys_rad_m=nx*2*math.pi/c.FULL_WIDTH_M,
                             wavelength_m=c.FULL_WIDTH_M/nx,
                             box_mode_m=m if abs(m-round(m))<1e-12 else None,
                             representable_in_box=abs(m-round(m))<1e-12 and round(m)>=1,
                             sigma_hat=sigma,sigma_times_t_hat_c=sigma*t_hat_c,
                             sigma_per_s=sigma/ref,critical_rayleigh=rc,
                             critical_rayleigh_over_R_E=rc/R_E,
                             is_wide_scan_peak=(nx==peak)))
        for m in BOX_M:
            nx=int(c.lsa_nx(m));sigma,rc=res[nx]
            box.append(dict(checkpoint=name,time_s=t_s,t_over_tc=t_hat/t_hat_c,box_mode_m=m,nx=nx,
                            wavelength_m=c.wavelength(m),k_phys_rad_m=c.k_phys(m),
                            cells_per_wavelength_R24=c.cells_per_wavelength(m,24),
                            cells_per_wavelength_R48=c.cells_per_wavelength(m,48),
                            cells_per_wavelength_R96=c.cells_per_wavelength(m,96),
                            sigma_hat=sigma,sigma_per_s=sigma/ref,sigma_times_t_hat_c=sigma*t_hat_c,
                            critical_rayleigh=rc,critical_rayleigh_over_R_E=rc/R_E))

    # fastest-growing box mode at each checkpoint, and the three modes the command names
    per_checkpoint={}
    for name,t_s in times:
        here=[r for r in box if r['checkpoint']==name]
        finite=[r for r in here if np.isfinite(r['sigma_hat'])]
        best=max(finite,key=lambda r:r['sigma_hat']) if finite else None
        wide_here=[r for r in wide if r['checkpoint']==name and r['is_wide_scan_peak']]
        per_checkpoint[name]=dict(time_s=t_s,t_over_tc=here[0]['t_over_tc'],
            fastest_box_mode_m=best['box_mode_m'] if best else None,
            fastest_box_sigma_per_s=best['sigma_per_s'] if best else None,
            wide_scan_peak_nx=wide_here[0]['nx'] if wide_here else None,
            wide_scan_peak_at_scan_edge=bool(wide_here and wide_here[0]['nx']>=max(WIDE_NX)),
            **{f'sigma_per_s_m{m}':next(r['sigma_per_s'] for r in here if r['box_mode_m']==m)
               for m in (2,6,10)},
            **{f'critical_rayleigh_over_R_E_m{m}':next(r['critical_rayleigh_over_R_E'] for r in here if r['box_mode_m']==m)
               for m in (2,6,10)})

    # ---- cumulative linear amplification of each named box mode over 0 -> 4T ---------
    grid=np.exp(np.linspace(math.log(0.2*t_hat_c),math.log(4.0*c.T_S/ref),60))
    amp_rows=[]
    sigmas={m:[] for m in (2,6,10)}
    for t_hat in grid:
        Dz,pp,S0=_state(t_hat,point)
        for m in (2,6,10):
            sigmas[m].append(_sigma_and_rc(Dz,pp,S0,point,R_E,int(c.lsa_nx(m)))[0])
    for m in (2,6,10):
        s=np.array(sigmas[m]);ln=np.concatenate([[0.0],np.cumsum(np.diff(grid)*(s[1:]+s[:-1])/2)])
        for t_hat,sv,lv in zip(grid,s,ln):
            amp_rows.append(dict(box_mode_m=m,t_hat=float(t_hat),time_s=float(t_hat*ref),
                                 t_over_T=float(t_hat*ref/c.T_S),sigma_hat=float(sv),
                                 ln_amplification_from_0p2tc=float(lv)))

    c.OUT.mkdir(parents=True,exist_ok=True)
    c.write_csv(c.OUT/'lsa_wavenumber_scan.csv',wide)
    c.write_csv(c.OUT/'lsa_box_modes.csv',box)
    c.write_csv(c.OUT/'lsa_mode_amplification.csv',amp_rows)
    report=dict(
        task='command/16 Task D: LSA wavenumber audit, case 36877',
        solver=dict(module=str(c.PAPER1/'lsa_unsaturated.py'),sha256=c.sha(c.PAPER1/'lsa_unsaturated.py'),
                    N_chebyshev=N_CHEB,chi_coeff_mode=CHI_MODE,ell_hat=c.ELL_HAT,
                    role='frozen Paper 1 read-only reference'),
        inputs=dict(cases_csv_sha256=case_sha,observations_csv_sha256=obs_sha,case_id=c.CASE),
        dimensionless_point=point,R_E=R_E,reference_time_s=ref,
        neutral_stability=dict(
            t_hat_onset=t_hat_c,t_onset_s=t_c_s,critical_nx=int(float(obs['critical_nx'])),
            critical_wavelength_m=c.FULL_WIDTH_M/float(obs['critical_nx']),
            crossing_class=obs['crossing_class'],stability_margin_log10=float(obs['stability_margin_log10']),
            bracket_sufficient=obs['bracket_sufficient']=='1',
            resolution_converged=obs['resolution_converged']=='1',
            resolution_relative_error=float(obs['resolution_relative_error']),
            seeded_nx28_t_hat=0.0017177719703013337,seeded_nx28_t_s=45474.35459078279,
            note='R_c(t) is the neutral sigma=0 curve; R_c at onset equals R_E by construction'),
        quasi_steady_assumptions=[
            'base state frozen at t; meaningful only where the base state changes slowly over 1/sigma',
            'single least-stable mode, no mode interaction and no nonlinear saturation',
            'laterally unbounded column with quantized a = nx*2*pi/ell_hat; the DuMuX box admits only nx = 14*m',
            'already implemented and self-tested in the frozen Paper 1 solver; nothing new is introduced here',
            'validated against DuMuX on cases 36545 and 36877 (2026-09-05) to within a factor of two near t_c'],
        reproduced_anchors=anchors,anchors_reproduced=bool(anchors_ok),
        box_mode_algebra=dict(width_m=c.WIDTH_M,full_width_m=c.FULL_WIDTH_M,nx_per_box_mode=c.NX_PER_M,
            canonical_box_mode=c.CANONICAL_M,canonical_nx=int(c.lsa_nx(2)),
            m6_nx=int(c.lsa_nx(6)),m10_nx=int(c.lsa_nx(10)),
            m2_wavelength_m=c.wavelength(2),m6_wavelength_m=c.wavelength(6),m10_wavelength_m=c.wavelength(10),
            m2_k_rad_m=c.k_phys(2),m6_k_rad_m=c.k_phys(6),m10_k_rad_m=c.k_phys(10),
            global_critical_nx_representable_in_box=False,
            note='the LSA global critical mode nx=5 has wavelength 0.12 m and does not fit the 0.0214 m box'),
        per_checkpoint=per_checkpoint,
        amplification_store=dict(ln_amplification_max=float(obs['ln_amplification_max']),
            nx_amp=int(float(obs['nx_amp'])),t_hat_amp=float(obs['t_hat_amp']),
            t_amp_s=t_amp_s,t_amp_over_T=t_amp_s/c.T_S,
            growth_stop_reason=obs['growth_stop_reason'],
            t_precipitation_s=float(obs['t_precipitation_s']),
            t_precipitation_over_T=float(obs['t_precipitation_s'])/c.T_S),
        wall_s=time.time()-start)
    c.save(c.OUT/'lsa_audit.json',report)
    print(json.dumps(dict(anchors_reproduced=report['anchors_reproduced'],
                          per_checkpoint={k:{kk:v[kk] for kk in ['t_over_tc','fastest_box_mode_m','wide_scan_peak_nx']}
                                          for k,v in per_checkpoint.items()},
                          wall_s=report['wall_s']),indent=2))


if __name__=='__main__':main()
