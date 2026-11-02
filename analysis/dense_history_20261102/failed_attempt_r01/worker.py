"""Read-only analysis of 14 frozen historical VTU sequences, executed on frodo.
Writes only into the separate analysis directory. Never invokes a solver.
"""
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import sys
import time
import traceback
import xml.etree.ElementTree as ET
import numpy as np

ROOT = Path(sys.argv[1]).resolve()
PLAN = json.loads((ROOT/'plan.json').read_text())


def sha(data):
    return hashlib.sha256(data).hexdigest()


def savecsv(path, data):
    if data:
        with path.open('w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(data[0])); w.writeheader(); w.writerows(data)


def parse_input(text):
    result, section = {}, None
    for line in text.splitlines():
        line = line.split('#', 1)[0].strip()
        if not line: continue
        if line.startswith('['): section = line.strip('[]')
        elif '=' in line:
            k, v = line.split('=', 1); result[section+'.'+k.strip()] = v.strip()
    return result


def analyze(spec):
    source = Path(spec['path']).resolve()
    target = ROOT/'runs'/spec['id']; target.mkdir(parents=True, exist_ok=False)
    checks, small_hashes = {}, {}
    for name, expected in spec['small_file_hashes'].items():
        data = (source/name).read_bytes(); small_hashes[name] = sha(data)
        checks['frozen_sha/'+name] = sha(data) == expected
        (target/name).write_bytes(data)
    inp = next(target.glob('*.input')); params = parse_input(inp.read_text())
    checks['input_matches_audit'] = params == spec['parameters']
    if not all(checks.values()): raise RuntimeError('Input/manifest/PVD/profile changed since frozen audit')
    sigma_file = next(source.glob('*_sigmatop.csv'))
    data = sigma_file.read_bytes(); checks['sigmatop_hash'] = sha(data) == spec['sigmatop_sha256']
    if not checks['sigmatop_hash']: raise RuntimeError('sigmatop changed since audit')
    (target/sigma_file.name).write_bytes(data)
    sig = np.loadtxt(sigma_file, comments='#', ndmin=2)
    datasets = list(ET.parse(next(target.glob('*.pvd'))).getroot().iter('DataSet'))
    times = [float(d.attrib['timestep']) for d in datasets]
    checks['101_monotonic_frames'] = len(times) == 101 and times[0] == 0 and np.all(np.diff(times) > 0)
    assert checks['101_monotonic_frames']
    phi = float(params['Problem.Porosity'])
    n, alpha, swr = [float(params['Problem.'+k]) for k in ['PcN', 'PcAlpha', 'Swr']]
    scalar = lambda k: float(params['Problem.'+k])
    mw, ms, r0, linear, quad = [scalar(k) for k in ['MolarMassH2O','MolarMassNaCl','PureWaterDensity','DensityLinearCoefficient','DensityQuadraticCoefficient']]
    frame_rows, region_rows, hashes, profile_mean, profile_sd, profile_sw, sample_fields = [], [], [], [], [], [], {}
    grid = None
    max_density_error, max_legacy_error, max_sw_error = 0., 0., 0.
    max_sigma_error, max_mean_error, max_time_match = 0., 0., 0.
    max_quadrature_error = 0.
    for fi, ds in enumerate(datasets):
        t = times[fi]; p = (source/ds.attrib['file']).resolve()
        assert p.is_relative_to(source)
        before = p.stat(); raw = p.read_bytes(); after = p.stat()
        assert (before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns)
        digest = sha(raw)
        tree = ET.fromstring(raw)
        piece = tree.find('.//Piece')
        arrays = {}
        wanted = {'Coordinates', 'pw', 'sw', 'moleFraction', 'rhoW'}
        for a in tree.iter('DataArray'):
            if a.attrib.get('Name') in wanted:
                assert a.attrib.get('format') == 'ascii'
                arrays[a.attrib['Name']] = np.fromstring(a.text or '', sep=' ')
        assert wanted == set(arrays)
        xyz = arrays['Coordinates'].reshape((-1,3))
        if grid is None:
            grid = xyz.copy(); xs = np.unique(xyz[:,0]); zs = np.unique(xyz[:,1])
            nx, nz = len(xs), len(zs)
            assert len(xyz) == nx*nz and np.unique(xyz[:,2]).size == 1
            assert int(piece.attrib['NumberOfCells']) == (nx-1)*(nz-1)
            ix, iz = np.searchsorted(xs, xyz[:,0]), np.searchsorted(zs, xyz[:,1])
            order = np.argsort(iz*nx+ix)
            assert np.array_equal((iz*nx+ix)[order], np.arange(nx*nz))
            width, height = xs[-1]-xs[0], zs[-1]-zs[0]
            wx = np.diff(np.r_[xs[0], .5*(xs[:-1]+xs[1:]), xs[-1]])
            wz = np.diff(np.r_[zs[0], .5*(zs[:-1]+zs[1:]), zs[-1]])
            weights = wz[:,None]*wx[None,:]/width
            checks['input_grid_shape'] = nx-1 == int(params['Grid.Cells0']) and nz-1 == int(params['Grid.Cells1'])
            checks['tensor_volume'] = abs(weights.sum()-height) < 1e-13
            mid = .5*(zs[:-1]+zs[1:]); cuts = []
            for q in [.9,.75,.5]:
                cut = mid[np.argmin(abs(mid-(zs[0]+q*height)))]
                cuts.append(dict(region='upper_'+str(q), requested_z_over_H=q, actual_z_over_H=(cut-zs[0])/height, mask=zs>cut))
            max_density_error = 0.
        else:
            assert np.array_equal(xyz, grid)
        x, sw, pw, rho = [arrays[k][order].reshape((nz,nx)) for k in ['moleFraction','sw','pw','rhoW']]
        assert all(np.all(np.isfinite(v)) for v in (x,sw,pw,rho))
        assert np.min(sw) >= 0 and np.max(sw) <= 1 and np.min(x) >= 0
        water, salt = phi*sw*weights, phi*sw*x*weights
        w = x*ms/((1-x)*mw+x*ms)
        err_density = float(np.max(abs(rho-r0-1000*w*(linear+quad*w))))
        wold = x*.0584428/((1-x)*.018+x*.0584428)
        err_legacy = float(np.max(abs(rho-1000-1000*wold*(.6708+.3739*wold))))
        pc = 1e5-pw
        sw_calc = swr+(1-swr)/(1+(alpha*np.maximum(pc,0))**n)**(1-1/n)
        err_sw = float(np.max(abs(sw-sw_calc)))
        max_density_error = max(max_density_error, err_density); max_legacy_error = max(max_legacy_error, err_legacy)
        max_sw_error = max(max_sw_error, err_sw)
        mean = (x*wx[None,:]).sum(axis=1)/width
        variance = ((x-mean[:,None])**2*wx[None,:]).sum(axis=1)/width
        sd = np.sqrt(variance)
        profile_mean.append(mean); profile_sd.append(sd); profile_sw.append((sw*wx[None,:]).sum(axis=1)/width)
        sig_index = int(np.argmin(abs(sig[:,0]-t)))
        time_error = float(abs(sig[sig_index,0]-t))
        max_time_match = max(max_time_match, time_error)
        sigma_node, mean_node = float(np.std(x[-1])), float(np.mean(x[-1]))
        sigma_error = abs(sigma_node-sig[sig_index,2]); mean_error = abs(mean_node-sig[sig_index,1])
        max_sigma_error = max(max_sigma_error, sigma_error); max_mean_error = max(max_mean_error, mean_error)
        fsalt = float(salt.sum()); fwater = float(water.sum())
        if fi in (0,len(times)-1):
            # Independent cell-corner quadrature, separate from nodal-width summation.
            area = np.diff(zs)[:,None]*np.diff(xs)[None,:]/width
            for v, total in [(phi*sw*x,fsalt),(phi*sw,fwater)]:
                independent = float(np.sum(area*(v[:-1,:-1]+v[:-1,1:]+v[1:,:-1]+v[1:,1:])/4))
                max_quadrature_error = max(max_quadrature_error, abs(independent-total)/abs(total))
        base = dict(id=spec['id'], case_id=spec['case_id'], frame=fi, time_s=t)
        frame_rows.append(dict(**base, whole_salt_m=fsalt, whole_water_m=fwater,
            salt_centroid_z_over_H=float(np.sum(salt*((zs-zs[0])/height)[:,None])/fsalt),
            top_x_mean_area=float(mean[-1]), top_x_mean_node=mean_node, top_x_max=float(x[-1].max()),
            domain_x_max=float(x.max()), domain_x_min=float(x.min()), top_sigma_node=sigma_node,
            max_layer_sigma=float(sd.max()), max_layer_sigma_over_mean=float(np.max(sd/mean)),
            x_lateral_rms=float(np.sqrt(np.sum(variance*wz)/height)),
            rho_input_max_error_kg_m3=err_density, rho_legacy_max_error_kg_m3=err_legacy,
            sw_input_max_error=err_sw, sigmatop_time_error_s=time_error,
            sigmatop_mean_error=mean_error, sigmatop_sigma_error=sigma_error))
        for cut in cuts:
            mask = cut['mask']
            region_rows.append(dict(**base, region=cut['region'], requested_z_over_H=cut['requested_z_over_H'],
                actual_z_over_H=cut['actual_z_over_H'], inventory_salt_m=float(salt[mask].sum()),
                inventory_water_m=float(water[mask].sum())))
        if fi in (0,25,50,75,100):
            stride = max(1,(nx-1)//300)
            sample_fields['x_'+str(fi)] = x[:,::stride].copy()
            sample_fields['sw_'+str(fi)] = sw[:,::stride].copy()
            sample_fields['x_anomaly_'+str(fi)] = (x-mean[:,None])[:,::stride].copy()
            sample_fields['xs'] = xs[::stride]
        hashes.append(dict(frame=fi, time_s=t, path=str(p), bytes=after.st_size, mtime_ns=after.st_mtime_ns, sha256=digest))
        if fi % 10 == 0:
            (ROOT/'progress.json').write_text(json.dumps(dict(id=spec['id'],frame=fi,total_frames=len(times),utc=datetime.now(timezone.utc).isoformat())))
        del tree, raw
    checks.update(density_matches_explicit_input=max_density_error<1e-8,
                  saturation_matches_input=max_sw_error<1e-10,
                  field_top_statistics_match_record=max_sigma_error<1e-10 and max_mean_error<1e-10,
                  frame_time_matches_diagnostics=max_time_match<max(1e-5, times[-1]*1e-10),
                  independent_quadrature=max_quadrature_error<1e-12)
    if any(not v for k,v in checks.items() if k not in ['density_matches_explicit_input','saturation_matches_input','field_top_statistics_match_record']):
        status = 'structural_review_required'
    elif not all(checks.values()): status = 'quarantine_source_or_constitutive_conflict'
    else: status = 'field_input_consistent_binary_unverified'
    savecsv(target/'frame_metrics.csv', frame_rows); savecsv(target/'inventory_timeseries.csv', region_rows)
    np.savez_compressed(target/'profiles.npz', time_s=times,z_m=zs,x_mean=profile_mean,x_sigma=profile_sd,sw_mean=profile_sw)
    np.savez_compressed(target/'representative_fields.npz', time_s=np.array(times)[[0,25,50,75,100]],z_m=zs,**sample_fields)
    (target/'field_manifest.json').write_text(json.dumps(hashes,indent=2))
    audit=dict(id=spec['id'], run_key=spec['run_key'], case_id=spec['case_id'], status=status, checks=checks,
        frames=len(times), points_per_frame=nx*nz, bytes_read=sum(h['bytes'] for h in hashes),
        actual_cuts=[{k:v for k,v in q.items() if k!='mask'} for q in cuts],
        max_density_error_kg_m3=max_density_error,max_legacy_density_error_kg_m3=max_legacy_error,
        max_sw_error=max_sw_error, max_sigmatop_mean_error=max_mean_error,max_sigmatop_sigma_error=max_sigma_error,
        max_quadrature_relative_error=max_quadrature_error,small_file_sha256=small_hashes,
        sigmatop_sha256=sha(data), historical_executable_identity='not established',
        physical_interpretation='inventory-derived net exchange; no independent face-flux closure, no unperturbed reference')
    (target/'audit.json').write_text(json.dumps(audit,indent=2))
    return audit


summary=[]
for spec in PLAN['runs']:
    t0=time.monotonic()
    try:
        result=analyze(spec); result['analysis_wall_s']=time.monotonic()-t0
    except Exception:
        result=dict(id=spec['id'],run_key=spec['run_key'],status='analysis_failed',error=traceback.format_exc())
        print(result['error'],flush=True)
    summary.append(result)
    (ROOT/'summary.json').write_text(json.dumps(summary,indent=2))
    print(spec['id'],result['status'],flush=True)
(ROOT/'finished.json').write_text(json.dumps(dict(completed=len(summary),planned=len(PLAN['runs']),
    analysis_failures=sum(r['status']=='analysis_failed' for r in summary),simulations_started=False)))
