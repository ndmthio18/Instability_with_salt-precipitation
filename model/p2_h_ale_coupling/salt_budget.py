"""Crust salt-budget decomposition for P2-H coupled runs (observer only).

Input is either the per-step ``salt_budget.csv`` written by the coupled solver
or a legacy ``coupled_history.csv``. A legacy history is differenced over
windows of at least ``--window`` seconds, so it can be subsampled.

The test is whether crust-base dissolution equals the gap between the
evaporation-limited crust precipitation ``C_eq j/(1-C_eq)`` and the soil's net
salt supply ``J_s = C_b J_L - D``. All rates are kg m^-2 s^-1 per unit soil area.

    python salt_budget.py RUN/salt_budget.csv [--out table.csv]
    python salt_budget.py RUN/coupled_history.csv --window 40 [--out table.csv]
"""
import argparse
import csv
import math
import sys

CEQ = 0.2649408697676898      # mass fraction at Transport.SolubilityMoleFraction
RHO_CRYSTAL = 2165.0
BIRTH = 0.43

def from_budget(rows):
    out = []
    for r in rows:
        f = lambda k: float(r[k])
        net, js = f('net_solid_kg_m2_s'), f('interface_salt_kg_m2_s')
        out.append({
            'time_s': f('time_s'), 'phase': r['phase'],
            'evaporation': f('evaporation_kg_m2_s'),
            'interface_liquid': f('interface_liquid_kg_m2_s'), 'interface_salt': js,
            'interface_minus_ceq': f('interface_minus_ceq'),
            'advective_salt': f('advective_salt_kg_m2_s'),
            'back_diffusion': f('implied_back_diffusion_kg_m2_s'),
            'evaporation_limited': f('evaporation_limited_precipitation_kg_m2_s'),
            'top_precipitation': f('top_precipitation_kg_m2_s'),
            'internal_positive': f('internal_positive_kg_m2_s'),
            'dissolution': f('dissolution_kg_m2_s'),
            'gross_precipitation': f('gross_precipitation_kg_m2_s'),
            'net_solid': net,
            'dissolution_implied': f('dissolution_implied_by_supply_kg_m2_s'),
            'predicted_interface_minus_ceq': f('predicted_interface_minus_ceq'),
            'salt_closure_residual': f('salt_closure_residual_kg_m2_s'),
            'water_closure_residual': f('water_closure_residual_kg_m2_s'),
        })
    return out


def from_history(rows, window, ceq):
    """Difference cumulative ledgers. Pore salt comes from J_s - net solid."""
    f = lambda r, k: float(r[k])
    out, a = [], rows[0]
    for b in rows[1:]:
        dt = f(b, 'time_s') - f(a, 'time_s')
        if dt < window:
            continue
        active_a, active_b = f(a, 'active') > 0, f(b, 'active') > 0
        if active_a != active_b:          # skip the window containing activation
            a = b
            continue
        j, cb = f(b, 'evaporation_kg_m2_s'), f(b, 'top_soil_mass_fraction')
        top = (f(b, 'top_solid_kg_m2') - f(a, 'top_solid_kg_m2')) / dt
        if active_b:
            jl = (f(b, 'interface_liquid_kg_m2') - f(a, 'interface_liquid_kg_m2')) / dt
            js = (f(b, 'interface_salt_kg_m2') - f(a, 'interface_salt_kg_m2')) / dt
            if 'internal_positive_kg_m2' in b:
                pos = (f(b, 'internal_positive_kg_m2') - f(a, 'internal_positive_kg_m2')) / dt
                diss = (f(b, 'internal_dissolution_kg_m2') - f(a, 'internal_dissolution_kg_m2')) / dt
            else:                         # older runs: only the net internal ledger exists
                pos = (f(b, 'internal_solid_kg_m2') - f(a, 'internal_solid_kg_m2')) / dt
                diss = math.nan
            net = (f(b, 'h_mass_m') - f(a, 'h_mass_m')) * RHO_CRYSTAL * (1 - BIRTH) / dt
        else:                             # history has no step fluxes before activation
            jl = js = math.nan
            pos = diss = 0.0
            net = top
        gross = top + pos
        out.append({
            'time_s': f(b, 'time_s'), 'phase': 'active' if active_b else 'precrust',
            'evaporation': j, 'interface_liquid': jl, 'interface_salt': js,
            'interface_minus_ceq': cb - ceq, 'advective_salt': cb * jl,
            'back_diffusion': cb * jl - js,
            'evaporation_limited': ceq * j / (1 - ceq),
            'top_precipitation': top, 'internal_positive': pos, 'dissolution': diss,
            'gross_precipitation': gross, 'net_solid': net,
            'dissolution_implied': gross - net,   # ledger vs inventory quadrature check
            'predicted_interface_minus_ceq': math.nan,
        })
        a = b
    return out


def summarize(table):
    pre = [r for r in table if r['phase'] == 'precrust']
    post = [r for r in table if r['phase'] == 'active']
    p_ref = pre[-1]['net_solid'] if pre else math.nan
    print(f"pre-crust precipitation P (last window)  {p_ref:.4e}")
    if not post:
        print('no active crust rows')
        return
    r = post[-1]
    g = lambda k: r[k]
    print(f"last active row t={g('time_s'):.1f} s")
    if math.isnan(g('dissolution')):
        print('  (no split reaction ledgers: "internal positive" is the net internal rate)')
    rows = [
        ('evaporation j', 'evaporation'),
        ('interface liquid J_L', 'interface_liquid'),
        ('advective salt C_b*J_L', 'advective_salt'),
        ('soil back-diffusion D', 'back_diffusion'),
        ('net soil salt supply J_s', 'interface_salt'),
        ('C_eq*j/(1-C_eq)', 'evaporation_limited'),
        ('gross crust precipitation', 'gross_precipitation'),
        ('  top', 'top_precipitation'),
        ('  internal positive', 'internal_positive'),
        ('dissolution (ledger)', 'dissolution'),
        ('dissolution implied', 'dissolution_implied'),
        ('net solid growth', 'net_solid'),
    ]
    for label, k in rows:
        print(f"  {label:28s} {g(k): .4e}")
    ratio = lambda x, y: x / y if y else math.nan
    print(f"  gross / evaporation limit    {ratio(g('gross_precipitation'), g('evaporation_limited')):.3f}")
    print(f"  dissolution / D              {ratio(g('dissolution'), g('back_diffusion')):.3f}")
    print(f"  net solid / pre-crust P      {ratio(g('net_solid'), p_ref):.3f}")
    print(f"  C_b - C_eq observed          {g('interface_minus_ceq'): .3e}")
    if not math.isnan(g('predicted_interface_minus_ceq')):
        print(f"  C_b - C_eq predicted         {g('predicted_interface_minus_ceq'): .3e}")
    if 'salt_closure_residual' in r:
        worst = lambda k: max(abs(x[k]) for x in post)
        print(f"  max |salt closure residual|  {worst('salt_closure_residual'):.3e}")
        print(f"  max |water closure residual| {worst('water_closure_residual'):.3e}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('path')
    ap.add_argument('--window', type=float, default=40.0, help='legacy differencing window [s]')
    ap.add_argument('--ceq', type=float, default=CEQ)
    ap.add_argument('--out')
    args = ap.parse_args()
    with open(args.path, newline='') as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        sys.exit('empty input')
    if 'phase' in rows[0]:
        table = from_budget(rows)
    elif 'interface_liquid_kg_m2' in rows[0]:
        table = from_history(rows, args.window, args.ceq)
    else:
        sys.exit('unrecognised input: expected salt_budget.csv or coupled_history.csv')
    summarize(table)
    if args.out:
        keys = list(table[0].keys())
        with open(args.out, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=keys)
            w.writeheader()
            w.writerows(table)


if __name__ == '__main__':
    main()
