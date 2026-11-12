"""Recheck delivered calculations; this does not establish physical validity."""
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd

R = Path(__file__).resolve().parent

def main():
    checks = []
    def require(condition, name):
        checks.append({'check': name, 'pass': bool(condition)})
        if not condition:
            raise AssertionError(name)
    cases = pd.read_csv(R/'cases.csv')
    if (R/'cases_refine.csv').exists():
        cases = pd.concat([cases, pd.read_csv(R/'cases_refine.csv')], ignore_index=True)
    s = pd.read_csv(R/'results/summary.csv')
    h = pd.read_csv(R/'results/history.csv.gz')
    require(set(s.id) == set(cases.id), 'all registered cases accounted for')
    failed=set(s.loc[s.status!='end','id'])
    require(failed == {'u0.6_strict','u0.6_strict_tighter'}, 'two registered temporal failures retained; all other endpoints reached')
    require(np.max(abs(h.hc-(h.zt-h.zb))) < 1e-15, 'geometric thickness from positions')
    require(np.max(abs(h.vnet-(h.vt-h.vb))) < 1e-18, 'net speed from independently evolved interfaces')
    require((h.solid >= 0).all(), 'nonnegative crystal inventory')
    require(((h.epsmin > .05) & (h.epsmax <= 1)).all(), 'accepted porosity within registered screens')
    require((h[['crust_salt','liquid_salt','water']] >= 0).all().all(), 'nonnegative liquid inventories')
    require(h.filter(regex='_res$').abs().max().max() < 2e-9, 'global and region inventory ledgers')
    require(h.water_jump.abs().max() < 1e-12, 'interface water jump')
    require(h.contact_res.abs().max() < 1e-12, 'interface concentration/flux matching')
    require(abs(h[h.u_input < 0].t.max()-160951.5) < 1e-8, 'no x008 forcing extrapolation')
    require((s[(s.u >= 0)&(s.status=='end')].end_s == 180000).all(), '50 h completed controlled forcing trajectories')
    require(np.max(abs(h.phi_mean*h.hc*2165-h.solid)) < 1e-12, 'computed mean solid fraction')
    decomposition = h.geometry_structure_term+h.mean_density_term+h.storage_term+h.supply_term
    require(np.max(abs((h.vnet-h.vnet_corr)-decomposition)) < 1e-18, 'exact net-growth discrepancy decomposition')
    # Raw files are kept locally/remotely; compressed published tables contain all histories.
    raw_manifests=[R/'provenance/raw_sha256.txt',R/'provenance/refine_raw_sha256.txt',R/'provenance/temporal_raw_sha256.txt']
    for line in '\n'.join(p.read_text() for p in raw_manifests if p.exists()).splitlines():
        if not line.strip(): continue
        digest, name = line.split(maxsplit=1)
        require(hashlib.sha256((R/name).read_bytes()).hexdigest() == digest, 'raw hash '+name)
    for line in (R/'provenance/final_source_sha256.txt').read_text().splitlines():
        if not line.strip(): continue
        digest, name = line.split(maxsplit=1)
        if name.startswith('source/'):
            path = R/'provenance'/name
            require(hashlib.sha256(path.read_bytes()).hexdigest() == digest, 'compiled source hash '+name)
            if path.name in ['two_front.hh','two_bdf.hh','verify.cc','run_case.cc']:
                require(path.read_bytes() == (R/path.name).read_bytes(), 'delivered source matches executed '+path.name)
    for k in range(1,9):
        paths = list((R/'figures').glob(f'figure{k}_*.svg'))
        require(len(paths) == 1, f'figure {k} SVG present')
        path = paths[0];root = ET.parse(path).getroot()
        require(not list(root.iter('{http://www.w3.org/2000/svg}image')), f'figure {k} remains vector')
        require('Times New Roman' in path.read_text(encoding='utf-8'), f'figure {k} Times New Roman')
        require(path.with_suffix('.png').exists() and path.with_suffix('.pdf').exists(), f'figure {k} exports')
    (R/'results/delivery_audit.json').write_text(json.dumps({'pass':True,'checks':checks},indent=2)+'\n')
    print(f'PASS {len(checks)} delivery checks; physical/model validity is evaluated separately.')

if __name__ == '__main__':
    main()
