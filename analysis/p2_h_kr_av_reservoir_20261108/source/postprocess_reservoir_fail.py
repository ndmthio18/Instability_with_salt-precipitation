"""Compact the failed 50 h pilot and produce controlled no-flux diagnostics."""
import csv, gzip, hashlib, json, pathlib, shutil
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

root=pathlib.Path.cwd(); case=root/'pilot/pilot_ref_ale64'; out=root/'analysis_fail'
res=out/'results'; figdir=out/'figures'; logs=out/'logs'; source=out/'source'
for p in (res,figdir,logs,source): p.mkdir(parents=True,exist_ok=True)

def compact(src,dst,interval=100.0):
    selected=[]; next_time=0.; first=last=None
    with src.open(newline='') as h:
        reader=csv.DictReader(h); names=reader.fieldnames
        for raw in reader:
            row={k:float(v) for k,v in raw.items()}; first=first or row; last=row
            if row['time_s']+1e-10>=next_time:
                selected.append(row)
                while next_time<=row['time_s']+1e-10: next_time+=interval
    if selected[-1]['time_s']!=last['time_s']: selected.append(last)
    with dst.open('w',newline='') as h:
        w=csv.DictWriter(h,fieldnames=names); w.writeheader(); w.writerows(selected)
    return first,last,selected

first,last,history=compact(case/'coupled_history.csv',res/'pilot_history.csv')
_,blast,boundary=compact(case/'open_boundary_budget.csv',res/'pilot_boundary.csv')
old_path=pathlib.Path('/WORK/hulm_work/zhixin/p2_h_kr_av_comparison_20261108/analysis/results/time_series/kr1_av1_ref.csv')
with old_path.open(newline='') as h: old=[{k:float(v) for k,v in r.items()} for r in csv.DictReader(h)]
profiles={}
with (case/'soil_fields.csv').open(newline='') as h:
    for raw in csv.DictReader(h):
        row={k:float(v) for k,v in raw.items()}; profiles.setdefault(row['time_s'],[]).append(row)
bottom_final=min(profiles[max(profiles)],key=lambda r:r['y_m'])
initial_salt=abs(first['soil_salt_kg']); initial_water=abs(first['soil_water_kg'])
summary={
  'status':'FAIL_NO_SURFACE_PRECIPITATION_NO_ALE_ACTIVATION','actual_end_s':last['time_s'],
  'surface_mass_fraction_initial':first['top_soil_mass_fraction'],'surface_mass_fraction_final':last['top_soil_mass_fraction'],
  'equilibrium_mass_fraction':0.264940870223385,'surface_sw_final':last['top_soil_sw'],'bottom_sw_final':bottom_final['saturation'],
  'final_evaporation_kg_m2_s':last['evaporation_kg_m2_s'],'cumulative_evaporation_kg_m2':last['vapour_kg']/0.001,
  'cumulative_bottom_water_kg_m2':blast['cumulative_bottom_water_kg']/0.001,
  'cumulative_bottom_salt_kg_m2':blast['cumulative_bottom_salt_kg']/0.001,
  'water_relative_residual':max(abs(r['water_residual_kg']) for r in history)/initial_water,
  'salt_relative_residual':max(abs(r['salt_residual_kg']) for r in history)/initial_salt,
  'surface_precipitated_solid_kg_m2':last['top_solid_kg_m2'],'ale_activated':int(last['active']),
  'old_no_flux_first_precipitation_s':108399.0,'old_no_flux_first_activation_attempt_s':108996.0,
  'old_no_flux_final_surface_sw':old[-1]['top_soil_sw'],'old_no_flux_final_surface_mass_fraction':old[-1]['top_soil_mass_fraction'],
}
(res/'pilot_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
with gzip.open(res/'soil_fields.csv.gz','wb',compresslevel=9) as z, (case/'soil_fields.csv').open('rb') as h: shutil.copyfileobj(h,z)

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8.5,'axes.labelsize':9,'legend.fontsize':8,'lines.linewidth':1.5})
def save(fig,name):
    for ext in ('png','pdf','svg'): fig.savefig(figdir/f'{name}.{ext}',dpi=300 if ext=='png' else None,bbox_inches='tight')
    plt.close(fig)

fig,ax=plt.subplots(2,1,figsize=(7.2,5.8),sharex=True,constrained_layout=True)
ax[0].plot([r['time_s']/3600 for r in history],[r['top_soil_mass_fraction'] for r in history],label='reservoir')
ax[0].plot([r['time_s']/3600 for r in old],[r['top_soil_mass_fraction'] for r in old],ls='--',label='old no_flux')
ax[0].axhline(summary['equilibrium_mass_fraction'],color='0.25',ls=':',label='equilibrium')
ax[1].plot([r['time_s']/3600 for r in history],[r['top_soil_sw'] for r in history])
ax[1].plot([r['time_s']/3600 for r in old],[r['top_soil_sw'] for r in old],ls='--')
ax[0].set_ylabel('Surface NaCl mass fraction'); ax[1].set_ylabel(r'Surface $S_w$'); ax[1].set_xlabel('Time (h)'); ax[0].legend(frameon=False)
save(fig,'figure1_reservoir_vs_no_flux')

fig,ax=plt.subplots(2,2,figsize=(8.5,6.2),sharex=True,constrained_layout=True)
t=[r['time_s']/3600 for r in boundary]; th=[r['time_s']/3600 for r in history]
ax[0,0].plot(t,[r['bottom_water_in_kg_m2_s'] for r in boundary]); ax[0,0].set_ylabel(r'Bottom water influx (kg m$^{-2}$ s$^{-1}$)')
ax[0,1].plot(th,[r['evaporation_kg_m2_s'] for r in history]); ax[0,1].set_ylabel(r'Top evaporation (kg m$^{-2}$ s$^{-1}$)')
ax[1,0].plot(t,[r['cumulative_bottom_water_kg']/0.001 for r in boundary]); ax[1,0].set_ylabel(r'Cumulative bottom water (kg m$^{-2}$)')
ax[1,1].plot(th,[r['top_soil_sw'] for r in history]); ax[1,1].set_ylabel(r'Surface $S_w$')
ax[1,1].ticklabel_format(axis='y',style='plain',useOffset=False)
for a in ax[1]: a.set_xlabel('Time (h)')
save(fig,'figure2_reservoir_supply')

times=(0.,36000.,90000.,180000.)
fig,ax=plt.subplots(1,2,figsize=(8.4,3.8),constrained_layout=True)
for time in times:
    p=min(profiles,key=lambda value:abs(value-time)); rows=sorted(profiles[p],key=lambda r:r['y_m'],reverse=True)
    depth=[(0.2-r['y_m'])*1000 for r in rows]
    ax[0].plot([r['mass_fraction'] for r in rows],depth,label=f'{p/3600:g} h')
    ax[1].plot([r['saturation'] for r in rows],depth,label=f'{p/3600:g} h')
for a in ax: a.set_ylim(200,0); a.set_ylabel('Depth below surface (mm)')
ax[0].set_xlabel('NaCl mass fraction'); ax[1].set_xlabel(r'$S_w$'); ax[0].legend(frameon=False)
ax[1].ticklabel_format(axis='x',style='plain',useOffset=False)
save(fig,'figure3_soil_profiles')

fig,ax=plt.subplots(1,2,figsize=(8.2,3.3),constrained_layout=True)
ax[0].plot(th,[r['water_residual_kg']/initial_water for r in history]); ax[1].plot(th,[r['salt_residual_kg']/initial_salt for r in history])
ax[0].set_ylabel('Relative water residual'); ax[1].set_ylabel('Relative salt residual')
for a in ax: a.set_xlabel('Time (h)')
save(fig,'figure4_conservation')

for name in ('solver.out','command.json','params.input','outcome.json','reservoir_config.csv','boundary_faces.csv','coupling_events.csv','donor_attempts.csv'):
    shutil.copy2(case/name,logs/name)
for name in ('verification_checks.json','pilot_gate.json'): shutil.copy2(root/name,logs/name)
for name in ('README.md','soil_methods.inc','RESERVOIR_FORMULATION.md','run_reservoir.py','verify_reservoir.py','postprocess_reservoir_fail.py'):
    shutil.copy2(root/'repo/model/p2_h_ale_coupling'/name,source/name)
checks=[]
from PIL import Image
for stem in ('figure1_reservoir_vs_no_flux','figure2_reservoir_supply','figure3_soil_profiles','figure4_conservation'):
    image=Image.open(figdir/f'{stem}.png').convert('RGB'); checks.append({'stem':stem,'pixels':image.size,'nonblank':any(a<b for a,b in image.getextrema()),'formats':all((figdir/f'{stem}.{e}').stat().st_size>1000 for e in ('png','pdf','svg'))})
(out/'figure_checks.json').write_text(json.dumps({'passed':all(c['nonblank'] and c['formats'] for c in checks),'checks':checks},indent=2)+'\n')
with (out/'CHECKSUMS.sha256').open('w') as manifest:
    for path in sorted(p for p in out.rglob('*') if p.is_file() and p.name!='CHECKSUMS.sha256'):
        manifest.write(hashlib.sha256(path.read_bytes()).hexdigest()+'  '+str(path.relative_to(out))+'\n')
print(json.dumps(summary,indent=2))
