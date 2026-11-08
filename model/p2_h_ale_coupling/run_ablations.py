"""Boundary-only test overlay, never linked into the DuMux production binary."""
import csv,json,os,pathlib,shutil,subprocess
assert os.environ.get('SLURM_JOB_ID')
root=pathlib.Path.cwd(); source=root/'repo/model/p2_h_ale_coupling'; build=root/'ablation_build'
build.mkdir(exist_ok=False); out=root/'ablations';out.mkdir(exist_ok=False)
s=(source/'ale.hh').read_text()
s=s.replace('    Parameters p;','    double forcedSalt=std::numeric_limits<double>::quiet_NaN();\n    Parameters p;').replace('#include <vector>','#include <vector>\n#include <limits>')
needle='        f.salt[n]=p.crystal*(1-p.birth)*f.speed;'
assert s.count(needle)==1
s=s.replace(needle,'        if(std::isfinite(forcedSalt)) f.salt[0]=forcedSalt;\n'+needle)
(build/'ale.hh').write_text(s)
for name in ('bdf.hh','ablation_driver.cc'): shutil.copyfile(source/name,build/name)
old=(root/'repo/model/p2_d_coupling/cpp/pratcrust.hh').read_text()
part=old.split('namespace prat {',1)[1].split('// ------------------------------------------------------ crust parameters',1)[0]
(build/'film.hh').write_text('#pragma once\n#include <cmath>\n#include <stdexcept>\n#include <string>\nnamespace p2gfilm {\n'+part+'\n}\n')
prefix=root/'dependencies/install'
cmd=['g++','-O2','-std=c++17',str(build/'ablation_driver.cc'),'-I'+str(prefix/'include'),'-o',str(build/'driver')]
cmd += [str(prefix/'lib64'/('libsundials_'+lib+'.a')) for lib in ['cvode','nvecserial','sunlinsoldense','sunmatrixdense']]
(build/'compile_command.json').write_text(json.dumps(cmd))
with (build/'compile.log').open('w') as f: subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
def read(p):
    with p.open() as f:return list(csv.DictReader(f))
for n in (160,320):
    path=root/f'ale_diagnostic/historical{n}_ale'
    rows=[r for r in read(path/'ale_diagnostics.csv') if r['status']=='accepted' and r['i']=='0']; start=float(rows[0]['time_s'])
    with (out/f'flux{n}.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['elapsed_s','bottom_c','adv','diff'])
        history=[r for r in read(path/'coupled_history.csv') if r['active']=='1']
        bytime={float(r['time_s']):r for r in rows}
        for before,after in zip(history,history[1:]):
            t=float(after['time_s']);dt=t-float(before['time_s']);r=bytime[t]
            # bottomC is constant inside each BDF/coupling trial. Therefore this
            # advective integral is exact from the accepted liquid ledger.
            adv=float(r['bottom_c'])*(float(after['interface_liquid_kg_m2'])-float(before['interface_liquid_kg_m2']))/dt
            total=(float(after['interface_salt_kg_m2'])-float(before['interface_salt_kg_m2']))/dt
            w.writerow([t-start,r['bottom_c'],adv,total-adv])
h=[r for r in read(root/'ale_diagnostic/historical320_ale/coupled_history.csv') if r['active']=='1'][0]
donor=[r for r in read(root/'ale_diagnostic/historical320_ale/donor_attempts.csv') if r['sufficient']=='1'][-1]
height=float(h['h_m']); salt=float(donor['salt_transfer_kg']);water=float(donor['water_transfer_kg'])
c=salt/(salt+water);rho=(salt+water)/(.001*.43*height)
results=[]
for eq in (0,1):
    for n,mode,end in [(320,'native',5),(320,'equilibrium_boundary',10),(320,'replay_total',5),(160,'replay_total',5),(320,'replay_advective',5),(320,'replay_diffusive',5)]:
        name=f'eq{eq}_n{n}_{mode}';d=out/name;d.mkdir()
        cmd=[str(build/'driver'),str(height),str(c),str(rho),str(eq),mode,str(out/f'flux{n}.csv'),str(end),'.025']
        (d/'command.json').write_text(json.dumps(cmd))
        with (d/'solver.out').open('w') as f:r=subprocess.run(cmd,cwd=d,stdout=f,stderr=subprocess.STDOUT,timeout=600)
        (d/'exit_status.txt').write_text(str(r.returncode));results.append(dict(case=name,exit_code=r.returncode))
(out/'outcomes.json').write_text(json.dumps(results,indent=2));print(json.dumps(results))
