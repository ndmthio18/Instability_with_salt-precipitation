"""Shared WP9 constants, paths and mode algebra. No solver is invoked here."""
import csv,hashlib,json,math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'analysis/wp9_r96_convergence_20261104'
CACHE=ROOT/'tmp/wp9_r96_convergence_20261104'
PAPER1=ROOT.parents[1]/'code'
HANDOFF=ROOT/'data/external/paper1/158395f3d25a/handoff'
CASE='36877'

# Frozen physical grid of the WP6/WP7/WP8 campaign.
T_S=136423.0637723484                 # three seeded onset times; the campaign unit
HORIZON_S=545692.2550893936           # 4T, the registered requested horizon
DT_S=113.68588647695695               # R24/R48/background fixed accepted step
WIDTH_M=0.021428571428571432          # 0.6/28 m, one canonical A28 wavelength
HEIGHT_M=0.2
FULL_WIDTH_M=0.6                      # the LSA box width, ell_hat*H
ELL_HAT=3.0
REFERENCE_TIME_S=26472870.30932611    # phi*H/E
SOLUBILITY=0.1

# The 2D box carries no-flow lateral boundaries, so its anomaly field expands in the
# Neumann cosine basis cos(m*pi*x/W).  Mode m has wavelength 2W/m, which is the LSA
# wavelength FULL_WIDTH_M/nx at nx = m*FULL_WIDTH_M/(2*WIDTH_M) = 14*m.
# The box is exactly one wavelength of the 28-wavelength A28 mode of the 0.6 m LSA
# domain, so the ratio is the integer 14 by construction. It is written as an exact
# constant rather than derived from WIDTH_M, whose stored decimal (0.021428571428571432)
# differs from 0.6/28 in the last bit and would otherwise make the algebra inexact.
NX_PER_M=14.0
WIDTH_CONSISTENCY=abs(WIDTH_M-FULL_WIDTH_M/(2*NX_PER_M))/WIDTH_M   # 1.6e-16, one ulp
CANONICAL_M=2                         # the seeded A28 mode

def lsa_nx(m):
    """LSA full-box wavenumber index for Neumann box mode m (exact for integer m)."""
    return NX_PER_M*m

def a_hat(nx):
    """Dimensionless angular wavenumber used by lsa_unsaturated (paper eq. 3.54)."""
    return nx*2.0*math.pi/ELL_HAT

def k_phys(m):
    """Physical angular wavenumber of box mode m, rad/m."""
    return m*math.pi/WIDTH_M

def wavelength(m):
    return 2.0*WIDTH_M/m

def cells_per_wavelength(m,nx_cells):
    return 2.0*nx_cells/m

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def save(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(v,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')

def rows(p):
    with Path(p).open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f))

def write_csv(p,records,fields=None):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    fields=fields or list(records[0])
    with p.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(records)

def case_record():
    """The frozen Paper 1 row for the case, with its file hash."""
    path=HANDOFF/'cases.csv'
    match=[r for r in rows(path) if r['case_id']==CASE]
    assert len(match)==1,'case row must be unique'
    return match[0],sha(path)

def observation_record():
    path=HANDOFF/'observations.csv'
    match=[r for r in rows(path) if r['case_id']==CASE and r['fidelity']=='lsa']
    assert len(match)==1,'lsa observation must be unique'
    return match[0],sha(path)

def lsa_point():
    """(kwargs for lsa_unsaturated, R_E) in the convention of parameter_database.solvers."""
    rec,_=case_record()
    n=float(rec['vg_n'])
    return dict(n=n,m=1.0-1.0/n,L=float(rec['mualem_l']),Swr=float(rec['residual_saturation']),
                delta=1.0/float(rec['gamma']),pbot=float(rec['p_hat_b']),
                beta=float(rec['beta'])),float(rec['r_e'])
