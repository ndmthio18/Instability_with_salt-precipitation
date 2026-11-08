"""Read-only geometric attribution of the 320-cell preactivation water defect.

No model, boundary selection, ledger or acceptance criterion is corrected here.
"""
import csv,hashlib,json,pathlib,tarfile
root=pathlib.Path.cwd();out=root/'audit';out.mkdir(exist_ok=True)
source=(root/'repo/model/p1_r2_conservative_richards/source/appl/simpleRichardsnc/main.cc').read_text()
assert 'eps_ = 1.5e-7' in source and 'scvf.ipGlobal()' in source
assert 'return globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_;' in source
archive=pathlib.Path('/WORK/hulm_work/zhixin/Greal81_coupling_regime_20261106/staging/exact_stack_sources_r04_full.tar.gz')
upstream={}
with tarfile.open(archive) as tar:
    for name in ('subcontrolvolumeface.hh','boxgeometryhelper.hh'):
        path='dumux/dumux/discretization/box/'+name
        data=tar.extractfile(path).read();upstream[path]=hashlib.sha256(data).hexdigest()
        (out/('reference_'+name+'.txt')).write_bytes(data)
        if name=='subcontrolvolumeface.hh':
            assert 'return center_;' in data.decode() and 'center_(Dumux::center(corners))' in data.decode()
def read(p):
    with p.open() as f:return [{k:float(v) for k,v in r.items()} for r in csv.DictReader(f)]
geometry=[]
for n in (40,80,160,320):
    p=read(root/'controlled'/f'control_n{n}_h1_d5/donor_profiles.csv')
    depths=sorted(set(r['depth_m'] for r in p));bad=[]
    for upper,lower in zip(depths,depths[1:]):
        dz=lower-upper
        # Rectangular tensor BOX: each boundary line is split into two halves;
        # upstream ipGlobal is the half-face center (quarter/three-quarter).
        for q in (.25,.75):
            depth=upper+q*dz
            if depth<1.5e-7:
                for side in ('left','right'):bad.append(dict(side=side,ip_depth_m=depth,area_m2=dz/2))
    area=sum(f['area_m2'] for f in bad)
    geometry.append(dict(n=n,top_cell_height_m=depths[1]-depths[0],top_box_depth_m=max(r['box_deepest_m'] for r in p if r['depth_m']==0),
        minimum_spacing_m=min(b-a for a,b in zip(depths,depths[1:])),maximum_spacing_m=max(b-a for a,b in zip(depths,depths[1:])),
        misclassified_side_faces=bad,unintended_side_area_m2=area,intended_top_area_m2=.001,extra_area_fraction=area/.001))
ratio=geometry[-1]['extra_area_fraction'];comparisons=[]
for name in ('old_n320','pre_n320','new_n320_h1_d2'):
    rows=read(root/'matrix'/name/'coupled_history.csv');assert all(r['active']==0 for r in rows)
    initial=rows[0]['soil_water_kg']
    observed=max(abs(r['water_residual_kg']) for r in rows)/initial
    unexplained=max(abs(r['water_residual_kg']+ratio*r['vapour_kg']) for r in rows)/initial
    comparisons.append(dict(case=name,last_time_s=rows[-1]['time_s'],original_max_relative_water_defect=observed,
        max_relative_difference_from_predicted_side_loss=unexplained,
        last_observed_water_defect_kg=rows[-1]['water_residual_kg'],
        last_predicted_unrecorded_side_loss_kg=-ratio*rows[-1]['vapour_kg']))
result=dict(source_boundary_tolerance_m=1.5e-7,geometry=geometry,comparisons=comparisons,
    upstream_source_hashes=upstream,interpretation='Read-only attribution. Actual side faces satisfy the legacy height-only top predicate. Predicted extra vapor loss is compared with the observed ledger defect. No run is repaired or retroactively accepted; top-face selection requires a separate correction and rerun.')
(out/'boundary_geometry.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
