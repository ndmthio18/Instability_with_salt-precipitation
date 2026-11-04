"""Evidence-bounded scalar/SCVF analysis for the single WP10 R24 trajectory."""
import csv, hashlib, json, sys
from pathlib import Path

T = 545692.2550893936 / 4.0

def read(path):
    with Path(path).open(newline="") as f: return list(csv.DictReader(f))

def prefix(old, new):
    a, b = Path(old).read_bytes().splitlines(), Path(new).read_bytes().splitlines()
    return {"old_lines": len(a), "new_lines": len(b), "bitwise_prefix": b[:len(a)] == a}

def nearest(rows, target):
    return min(rows, key=lambda r: abs(float(r["time_s"])-target))

def main(run, old, out):
    run, old, out = Path(run), Path(old), Path(out)
    out.mkdir(parents=True, exist_ok=True)
    stop = json.loads((run/"stop_record.json").read_text())
    surface, online, budget = [read(run/n) for n in
        ("surface_events.csv", "weighted_online.csv", "transport_budget.csv")]
    result = {
        "stop": dict(stop, time_over_T=stop["last_accepted_time_s"]/T),
        "through_4T": {
            n: prefix(old/n, run/n) for n in
            ("surface_events.csv", "weighted_online.csv", "transport_budget.csv")
        },
        "checkpoints": [], "planes": {},
        "sign_convention": "positive is downward; values are signed plane-integrated conservative SCVF fluxes",
        "local_split_available": False,
    }
    targets = [2,3,4,5] + list(range(10, int(stop["last_accepted_time_s"]//T)+1, 5))
    if not targets or abs(stop["last_accepted_time_s"]-targets[-1]*T)>1e-6:
        targets.append(stop["last_accepted_time_s"]/T)
    online_by = {int(r["step"]): r for r in online}
    for q in targets:
        r = nearest(surface, q*T); o = online_by[int(r["step"])]
        result["checkpoints"].append({
            "nominal_T": q, "step": int(r["step"]), "time_s": float(r["time_s"]),
            "actual_T": float(r["time_s"])/T,
            "surface_mean": float(r["x_top_area_mean"]),
            "surface_max": float(r["x_top_max"]), "domain_max": float(r["x_domain_max"]),
            "surface_sigma": float(r["sigma_top_node"]),
            "fundamental": float(o["weighted_fourier"]),
            "harmonic2": float(o["harmonic2_top"]), "harmonic3": float(o["harmonic3_top"]),
        })
    regions = sorted({r["region"] for r in budget})
    for region in regions:
        rr = [r for r in budget if r["region"] == region]
        adv = [float(r["cut_salt_adv_down_m_s"]) for r in rr]
        diff = [float(r["cut_salt_diff_down_m_s"]) for r in rr]
        total = [float(r["cut_salt_total_down_m_s"]) for r in rr]
        water = [float(r["cut_water_down_m_s"]) for r in rr]
        def first(vals, positive):
            for r,v in zip(rr, vals):
                if (v>0 if positive else v<0): return {"step":int(r["step"]),"time_s":float(r["time_s"]),"time_over_T":float(r["time_s"])/T,"value":v}
        def window(vals, threshold=1e-15):
            hit=[(r,v) for r,v in zip(rr,vals) if v>threshold]
            if not hit:return None
            a,b=hit[0],hit[-1]
            return {"threshold":threshold,"first_step":int(a[0]["step"]),"first_T":float(a[0]["time_s"])/T,
                    "last_step":int(b[0]["step"]),"last_T":float(b[0]["time_s"])/T,
                    "positive_step_fraction":len(hit)/len(vals),"maximum":max(vals)}
        result["planes"][region] = {
            "requested_z_over_H": float(rr[0]["requested_z_over_H"]),
            "actual_z_over_H": float(rr[0]["actual_z_over_H"]),
            "depth_below_surface_mm": (1-float(rr[0]["actual_z_over_H"]))*200,
            "adv_min": min(adv), "adv_max": max(adv), "diff_min": min(diff), "diff_max": max(diff),
            "total_min": min(total), "total_max": max(total),
            "water_min": min(water), "water_max": max(water),
            "adv_down_step_fraction": sum(v>0 for v in adv)/len(adv),
            "total_down_step_fraction": sum(v>0 for v in total)/len(total),
            "water_down_step_fraction": sum(v>0 for v in water)/len(water),
            "first_adv_down": first(adv, True), "first_adv_up": first(adv, False),
            "first_total_down": first(total, True), "first_total_up": first(total, False),
            "first_water_down": first(water, True),
            "adv_down_window_gt_1e-15":window(adv),
            "total_down_window_gt_1e-15":window(total),
            "final": {
                "adv": adv[-1], "diff": diff[-1], "total": total[-1],
                "inventory_salt_m": float(rr[-1]["inventory_salt_m"]),
                "inventory_change_m": float(rr[-1]["inventory_salt_m"])-float(rr[0]["inventory_salt_m"]),
                "cumulative_adv_m": float(rr[-1]["cumulative_cut_adv_m"]),
                "cumulative_diff_m": float(rr[-1]["cumulative_cut_diff_m"]),
            },
            "checkpoints": [{
                "nominal_T":q,
                "actual_T":float((r:=nearest(rr,q*T))["time_s"])/T,
                "adv":float(r["cut_salt_adv_down_m_s"]),
                "diff":float(r["cut_salt_diff_down_m_s"]),
                "total":float(r["cut_salt_total_down_m_s"]),
                "water":float(r["cut_water_down_m_s"]),
                "inventory_salt_m":float(r["inventory_salt_m"]),
            } for q in targets],
        }
    result["hashes"] = {n: hashlib.sha256((run/n).read_bytes()).hexdigest() for n in
        ("params.input","eigenprofile.dat","surface_events.csv","weighted_online.csv","transport_budget.csv","stop_record.json")}
    (out/"scalar_flux_metrics.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))

if __name__ == "__main__": main(*sys.argv[1:])
