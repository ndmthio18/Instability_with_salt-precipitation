"""Freeze WP7 R48, benchmark and one-cell background before any new result."""
import datetime,hashlib,json,shutil
from pathlib import Path
from single_wavelength_post import params

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'analysis/wp7_shortwave_validation_20261103'
BASE=ROOT/'analysis/wp6_4T_exploratory_20261103/runs'
T=136423.0637723484; DT=113.68588647695695
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,v): p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(v,indent=2)+'\n')
def replace(src,changes):
    lines=src.read_text().splitlines();seen=set();out=[]
    for line in lines:
        if '=' in line and not line.lstrip().startswith('#'):
            key=line.split('=',1)[0].strip()
            if key in changes: line=f'{key} = {changes[key]}';seen.add(key)
        out.append(line)
    assert seen==set(changes),(seen,changes)
    return '\n'.join(out)+'\n'
def main():
    assert not OUT.exists();OUT.mkdir(parents=True)
    command=ROOT/'command/11.md'; assert command.stat().st_size>1000
    source={r:BASE/f'WP6-{r}-explore4T' for r in ['M0','M1']}
    elf=json.loads((source['M1']/'plan.json').read_text())['elf_sha256']
    reg=dict(schema_version=1,id='WP7-E01',created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
      authority='command/11.md',base_commit='2443a3478c3146652980578e866e30382b73bf00',
      scope=['one R48 M1 ten-step cost benchmark','one R48 M0/M1 pair from t=0 to 4T if admitted','one one-cell extruded background from t=0 to T','isolated local-SCVF diagnostic build/bridge; no 8T'],
      fixed=dict(case_id=36877,width_m=.6/28,physical_mode='canonical nx28',wave_number=1,nz=200,dt_s=DT,T_s=T,R48_TEnd_s=4*T,amplitude_M1=1e-6,amplitude_M0=0,profile_sha256=sha(source['M1']/'eigenprofile.dat'),elf_sha256=elf,cuts=json.loads((source['M1']/'transport_metadata.json').read_text())['cuts'][:7]),
      cost=dict(existing_reduced_package_cap_core_h=28,R48_pair_registered_planning_cap_core_h=4.373333333333333,hard_wall_cap_s_per_leg=43200,new_tight_worker_limit_s=7200,benchmark_limit_s=600,background_limit_s=1200,parallel_physical_cores=2,minimum_free_disk_bytes=200*1024**3,minimum_available_memory_bytes=16*1024**3,
        admission='10-step M1 benchmark extrapolated solver + conservative post/archive estimate <=7200 s per leg and pair <=4.373333 core-h'),
      preregistered_metrics=dict(
        common_times='exact accepted physical times only; focus 3T-4T',
        spectral='Neumann physical k=m*pi/W; compare physical wavelength, not raw mode index',
        convergence=['dominant physical wavelength R48 within 10% of R24 and not shorter than 8 cells/wavelength','top three physical-mode energy ordering stable','weighted sigma/Fourier/projection within 20% over 3T-4T','M0/M1 background ratio and symmetry reported, never silently waived'],
        transport=['seven-plane signed advective/diffusive/total flux','initial-inventory-corrected Delta C','salt residual 1e-14 m; water residual 1e-11 m','candidate signal direction stable and change <=20% for quantitative GO'],
        background=['one-cell background lateral range <=1e-14','early-time mean profiles and seven inventories versus 2D M0 <=1e-6 relative or documented absolute near-zero criterion'],
        stopping=['x>=0.1 independent first crossing','no NaN/invalid state','no unchanged R24 8T','R96 only proposed after R48 cost and resolution evidence']),
      classification=dict(late_convection='grid-resolved only if wavelength and growth gates pass',local_downward_transport='pending unless conservative spatial face diagnostic is bridged and executed',t_eff=None,training_labels=False),
      retained_failures=['R24 M0 absolute symmetry failure','R24 M1 relative symmetry failure','R24 engineering background ratio failure'],
      forbidden=['unchanged R24 8T','R96 execution','multiwavelength execution','large campaign','restart','Paper 1 changes'])
    write(OUT/'preregistration.json',reg)
    specs=[
      ('WP7-R48-M0-4T','M0',48,4*T,1,7200,'r48pair'),('WP7-R48-M1-4T','M1',48,4*T,1,7200,'r48pair'),
      ('WP7-R48-M1-bench10','M1',48,10*DT,1,600,'preflight'),
      ('WP7-1D-background-T','M0',1,T,1,1200,'preflight')]
    checks={}
    for name,role,nx,tend,outint,wall,stage in specs:
        d=OUT/'runs'/name;d.mkdir(parents=True);src=source[role]
        text=replace(src/'params.input',{'Name':name,'Cells0':str(nx),'TEnd':repr(tend),'OutputInterval':str(outint)})
        (d/'params.input').write_text(text);shutil.copyfile(src/'eigenprofile.dat',d/'eigenprofile.dat')
        p=params(d/'params.input');assert abs(float(p['Grid.Positions0'].split()[-1])-.6/28)<1e-15
        assert int(p['Grid.Cells1'])==200 and float(p['TimeLoop.DtInitial'])==DT and float(p['TimeLoop.MaxTimeStepSize'])==DT
        plan=dict(id=name,role=role,stage=stage,case_id=36877,nx_cells=nx,nz_cells=200,T_s=tend,dt_s=DT,
          scalar_unknowns=2*(nx+1)*201,input_sha256=sha(d/'params.input'),profile_sha256=sha(d/'eigenprofile.dat'),elf_sha256=elf,
          preregistration_sha256=sha(OUT/'preregistration.json'),hard_wall_limit_s=wall,source_input_sha256=sha(src/'params.input'))
        write(d/'plan.json',plan)
        delta={k for k in params(src/'params.input') if params(src/'params.input')[k]!=p[k]}
        allowed={'Problem.Name','Grid.Cells0','TimeLoop.TEnd'}
        checks[name]=dict(allowed_input_diff=sorted(allowed),actual=sorted(delta),passed=delta<=allowed and {'Problem.Name','Grid.Cells0'}<=delta)
    assert all(v['passed'] for v in checks.values())
    write(OUT/'static_checks.json',dict(passed=True,legs=checks,command_sha256=sha(command),old_R24_manifest_unchanged=True))
    print('WP7-E01 frozen:',OUT)
if __name__=='__main__':main()
