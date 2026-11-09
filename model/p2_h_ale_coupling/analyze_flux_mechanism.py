"""Independent P2-H ALE face-flux reconstruction and mechanism audit.

The script deliberately re-implements the equations in ``ale.hh``.  It never
calls the production C++ flux function.  Production scalar diagnostics are
used only as boundary inputs (evaporation) and for independent comparisons.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import pathlib
from collections import defaultdict

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


CEQ = 0.2649408697676898
RHO_CRYSTAL = 2165.0
BIRTH = 0.43
DIFFUSION = 1.3e-9
REACTION = 1.0e-3
GRAIN = 25e-6
BETA = 2.0
TARGETS = [0, .5, 1, 2, 5, 10, 30, 100, 200, 300, 400, 500, 800, 1200]


def number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return value


def read_csv(path):
    with pathlib.Path(path).open(newline="") as handle:
        return [{key: number(value) for key, value in row.items()}
                for row in csv.DictReader(handle)]


def write_csv(path, rows, fields=None):
    rows = list(rows)
    if fields is None:
        fields = []
        for row in rows:
            for key in row:
                if key not in fields:
                    fields.append(key)
    with pathlib.Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path):
    digest = hashlib.sha256()
    with pathlib.Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def groups(path, statuses=("accepted",)):
    data = defaultdict(list)
    for row in read_csv(path):
        if row["status"] in statuses:
            data[(row["status"], row["time_s"])].append(row)
    return data


def density_from_group(rows):
    h = rows[0]["h_m"]
    denominator = sum(row["epsilon"] * row["c"] *
                      (row["z_top_m"] - row["z_bottom_m"])
                      for row in rows)
    return rows[0]["pore_salt_kg_m2"] / denominator


def reconstruct(rows, rho):
    """Reconstruct every face from stored primitive variables and geometry."""
    n = len(rows)
    h = rows[0]["h_m"]
    eps = np.array([r["epsilon"] for r in rows], dtype=float)
    concentration = np.array([r["c"] for r in rows], dtype=float)
    z0 = np.array([r["z_bottom_m"] for r in rows], dtype=float)
    z1 = np.array([r["z_top_m"] for r in rows], dtype=float)
    faces = np.concatenate((z0 / h, [z1[-1] / h]))
    centers = .5 * (faces[:-1] + faces[1:])
    widths = faces[1:] - faces[:-1]
    volume = h * widths
    deff = np.sqrt(eps) * DIFFUSION
    av = 36 * (1 - eps) / GRAIN
    reaction = rho * av * REACTION * (concentration - CEQ)
    evaporation = rows[0]["evap_kg_m2_s"]

    # Independent top quadratic, with evaporation supplied by the gas-film
    # observer. The resulting top concentration and speed are checked against
    # production diagnostics below.
    g = rho * eps[-1] * deff[-1] / (h * (1 - centers[-1]))
    a = BETA * rho * (1 - BIRTH) * REACTION
    b = evaporation - a * (1 - CEQ) - g
    d = evaporation * CEQ + g * (concentration[-1] - CEQ)
    discriminant = b * b - 4 * a * d
    if discriminant < 0:
        raise RuntimeError("independent top discriminant invalid")
    q = -.5 * (b + math.copysign(math.sqrt(discriminant), b))
    candidates = [d / q, q / a]
    valid = [x for x in candidates if math.isfinite(x) and 0 <= CEQ + x <= 1]
    if not valid:
        raise RuntimeError("independent top root absent")
    excess = min(valid, key=abs)
    top_c = CEQ + excess
    speed = BETA * rho / RHO_CRYSTAL * REACTION * excess

    liquid = np.zeros(n + 1)
    mesh = np.zeros(n + 1)
    advective = np.full(n + 1, np.nan)
    diffusive = np.full(n + 1, np.nan)
    salt = np.zeros(n + 1)
    relative = np.full(n + 1, np.nan)
    conductance = np.full(n, np.nan)
    delta_c = np.full(n, np.nan)
    liquid[n] = (rho * BIRTH + RHO_CRYSTAL * (1 - BIRTH)) * speed + evaporation
    for i in range(n - 1, -1, -1):
        liquid[i] = liquid[i + 1] + (1 - rho / RHO_CRYSTAL) * reaction[i] * volume[i]
    mesh[n] = -speed * BIRTH
    for i in range(1, n):
        eface = eps[i] if speed >= 0 else eps[i - 1]
        mesh[i] = -faces[i] * speed * eface
        resistance = h * ((faces[i] - centers[i - 1]) / (eps[i - 1] * deff[i - 1]) +
                          (centers[i] - faces[i]) / (eps[i] * deff[i]))
        relative[i] = liquid[i] + rho * mesh[i]
        advective[i] = relative[i] * (concentration[i - 1] if relative[i] >= 0
                                      else concentration[i])
        conductance[i] = rho / resistance
        delta_c[i] = concentration[i] - concentration[i - 1]
        diffusive[i] = -conductance[i] * delta_c[i]
        salt[i] = advective[i] + diffusive[i]
    bottom_c = rows[0]["bottom_c"]
    relative[0] = liquid[0]
    advective[0] = bottom_c * liquid[0]
    conductance[0] = rho * eps[0] * deff[0] / (h * centers[0])
    delta_c[0] = concentration[0] - bottom_c
    diffusive[0] = -conductance[0] * delta_c[0]
    salt[0] = advective[0] + diffusive[0]
    # The upper entry is a solid birth flux, not a liquid advective/diffusive face.
    salt[n] = RHO_CRYSTAL * (1 - BIRTH) * speed
    return {
        "n": n, "h": h, "eps": eps, "c": concentration, "faces": faces,
        "centers": centers, "widths": widths, "volume": volume, "deff": deff,
        "av": av, "reaction": reaction, "evap": evaporation, "top_c": top_c,
        "speed": speed, "liquid": liquid, "mesh": mesh, "relative": relative,
        "adv": advective, "diff": diffusive, "salt": salt,
        "conductance": conductance, "delta_c": delta_c, "bottom_c": bottom_c,
    }


def accepted_states(case_dir, case_name):
    raw = groups(pathlib.Path(case_dir) / "ale_diagnostics.csv", ("accepted", "bound_rejected"))
    accepted = []
    rejected = []
    for (status, time), values in sorted(raw.items(), key=lambda item: item[0][1]):
        n = int(max(row["i"] for row in values)) + 1
        if len(values) % n:
            raise RuntimeError("incomplete ALE diagnostic block at time %s" % time)
        for offset in range(0, len(values), n):
            rows = sorted(values[offset:offset + n], key=lambda row: int(row["i"]))
            rho = density_from_group(rows)
            state = reconstruct(rows, rho)
            state.update({"time": time, "status": status, "rows": rows,
                          "rho": rho, "case": case_name})
            (accepted if status == "accepted" else rejected).append(state)
    if not accepted:
        raise RuntimeError("no accepted ALE diagnostics in " + str(case_dir))
    activation = accepted[0]["time"]
    for state in accepted + rejected:
        state["elapsed"] = state["time"] - activation
    return accepted, rejected


def state_metrics(state):
    r = state
    reaction_mass = r["reaction"] * r["volume"]
    positive = np.maximum(reaction_mass, 0).sum()
    dissolution = np.maximum(-reaction_mass, 0).sum()
    net = reaction_mass.sum()
    ell = np.sqrt(r["eps"] * r["deff"] /
                  ((1 - CEQ) * r["av"] * REACTION))
    reverse = np.abs(r["diff"][:r["n"]])
    denominator = max(reverse[0], np.finfo(float).tiny)
    ratio = reverse / denominator
    k = int(np.nanargmin(ratio))
    negative_width = r["widths"][r["reaction"] < 0].sum() * r["h"]
    positive_width = r["widths"][r["reaction"] > 0].sum() * r["h"]
    crossings = []
    for i in range(r["n"] - 1):
        a, b = r["reaction"][i], r["reaction"][i + 1]
        if a == 0 or a * b < 0:
            fraction = abs(a) / (abs(a) + abs(b)) if a != b else 0
            crossings.append((r["centers"][i] + fraction *
                              (r["centers"][i + 1] - r["centers"][i])) * r["h"])
    sign_z = crossings[0] if crossings else math.nan
    pore = r["rho"] * np.sum(r["eps"] * r["c"] * r["volume"])
    return {
        "case": r["case"], "time_s": r["time"], "elapsed_s": r["elapsed"],
        "status": r["status"], "h_m": r["h"], "rho_kg_m3": r["rho"],
        "bottom_c": r["bottom_c"], "bottom_minus_ceq": r["bottom_c"] - CEQ,
        "top_c": r["top_c"], "speed_m_s": r["speed"],
        "evaporation_kg_m2_s": r["evap"], "g_top_kg_m2_s": r["salt"][-1],
        "bottom_liquid_kg_m2_s": r["liquid"][0],
        "bottom_salt_kg_m2_s": r["salt"][0],
        "bottom_advective_kg_m2_s": r["adv"][0],
        "bottom_diffusive_kg_m2_s": r["diff"][0],
        "r_positive_kg_m2_s": positive, "r_dissolution_kg_m2_s": dissolution,
        "r_net_kg_m2_s": net, "pore_salt_kg_m2": pore,
        "eps_min": float(r["eps"].min()), "eps_max": float(r["eps"].max()),
        "z_eps_min_m": float(r["centers"][np.argmin(r["eps"])] * r["h"]),
        "z_eps_max_m": float(r["centers"][np.argmax(r["eps"])] * r["h"]),
        "diff_ratio_min": float(ratio[k]), "z_diff_min_m": float(r["faces"][k] * r["h"]),
        "eps_at_diff_min": float(r["eps"][min(k, r["n"] - 1)]),
        "gradient_at_diff_min": float(r["delta_c"][k]),
        "conductance_at_diff_min": float(r["conductance"][k]),
        "gradient_ratio_at_diff_min": float(abs(r["delta_c"][k]) /
                                                  max(abs(r["delta_c"][0]), np.finfo(float).tiny)),
        "conductance_ratio_at_diff_min": float(r["conductance"][k] /
                                                     max(r["conductance"][0], np.finfo(float).tiny)),
        "epsilon_3over2_ratio_at_diff_min": float((r["eps"][min(k, r["n"] - 1)] /
                                                        r["eps"][0]) ** 1.5),
        "reaction_at_diff_min_kg_m3_s": float(r["reaction"][min(k, r["n"] - 1)]),
        "ell_min_m": float(ell.min()), "ell_median_m": float(np.median(ell)),
        "ell_bottom_m": float(ell[0]), "h_over_ell_bottom": float(r["h"] / ell[0]),
        "h_over_ell_median": float(r["h"] / np.median(ell)),
        "dissolution_zone_m": float(negative_width), "precipitation_zone_m": float(positive_width),
        "reaction_sign_change_m": sign_z,
        "g_top_over_saturated_evap_reference": r["salt"][-1] /
            (CEQ * r["evap"] / (1 - CEQ)),
    }


def validation(states):
    result = {"max_abs": defaultdict(float), "max_rel": defaultdict(float)}
    for state in states:
        rows, r = state["rows"], state
        observed = {
            "reaction": np.array([x["reaction_kg_m3_s"] for x in rows]),
            "deff": np.array([x["deff_m2_s"] for x in rows]),
            "av": np.array([x["av_m_inv"] for x in rows]),
            "speed": np.array([rows[0]["speed_m_s"]]),
            "bottom_liquid": np.array([rows[0]["bottom_liquid_kg_m2_s"]]),
            "bottom_salt": np.array([rows[0]["bottom_salt_kg_m2_s"]]),
            "bottom_adv": np.array([rows[0]["bottom_adv_kg_m2_s"]]),
            "bottom_diff": np.array([rows[0]["bottom_diff_kg_m2_s"]]),
        }
        rebuilt = {
            "reaction": r["reaction"], "deff": r["deff"], "av": r["av"],
            "speed": np.array([r["speed"]]), "bottom_liquid": np.array([r["liquid"][0]]),
            "bottom_salt": np.array([r["salt"][0]]), "bottom_adv": np.array([r["adv"][0]]),
            "bottom_diff": np.array([r["diff"][0]]),
        }
        for name in observed:
            delta = np.max(np.abs(observed[name] - rebuilt[name]))
            scale = max(np.max(np.abs(observed[name])), np.max(np.abs(rebuilt[name])), 1e-30)
            result["max_abs"][name] = max(result["max_abs"][name], float(delta))
            result["max_rel"][name] = max(result["max_rel"][name], float(delta / scale))
    return {key: dict(value) for key, value in result.items()}


def nearest(states, elapsed):
    return min(states, key=lambda state: abs(state["elapsed"] - elapsed))


def selected_states(states, rejected):
    selected = []
    seen = set()
    for target in TARGETS:
        if target <= states[-1]["elapsed"] + 1e-9:
            state = nearest(states, target)
            key = (state["status"], state["time"])
            if key not in seen:
                state["selection"] = "elapsed_%g_s" % target
                selected.append(state); seen.add(key)
    extras = [(states[-1], "last_accepted")]
    for threshold in (.42,):
        hit = next((s for s in states if s["eps"].min() < threshold), None)
        if hit: extras.append((hit, "first_eps_min_below_%g" % threshold))
    for threshold in (.60, .80, .90):
        hit = next((s for s in states if s["eps"][0] >= threshold), None)
        if hit: extras.append((hit, "first_bottom_eps_%g" % threshold))
    if rejected:
        extras.append((rejected[0], "first_rejected"))
    for state, label in extras:
        key = (state["status"], state["time"])
        if key not in seen:
            state["selection"] = label
            selected.append(state); seen.add(key)
    return sorted(selected, key=lambda state: state["time"])


def profile_rows(selected):
    face_rows, reaction_rows = [], []
    for s in selected:
        ell = np.sqrt(s["eps"] * s["deff"] /
                      ((1 - CEQ) * s["av"] * REACTION))
        for face in range(s["n"] + 1):
            face_rows.append({
                "case": s["case"], "selection": s["selection"], "status": s["status"],
                "time_s": s["time"], "elapsed_s": s["elapsed"], "face": face,
                "z_m": s["faces"][face] * s["h"], "z_over_h": s["faces"][face],
                "face_type": "bottom_liquid" if face == 0 else
                             ("top_solid_birth" if face == s["n"] else "interior_liquid"),
                "j_liquid_kg_m2_s": s["liquid"][face],
                "rho_j_mesh_kg_m2_s": s["rho"] * s["mesh"][face],
                "relative_liquid_kg_m2_s": s["relative"][face],
                "j_adv_kg_m2_s": s["adv"][face], "j_diff_kg_m2_s": s["diff"][face],
                "diff_ratio_to_bottom": abs(s["diff"][face]) /
                    max(abs(s["diff"][0]), np.finfo(float).tiny) if face < s["n"] else math.nan,
                "j_net_or_g_top_kg_m2_s": s["salt"][face],
                "g_top_kg_m2_s": s["salt"][face] if face == s["n"] else math.nan,
            })
        for i in range(s["n"]):
            reaction_rows.append({
                "case": s["case"], "selection": s["selection"], "status": s["status"],
                "time_s": s["time"], "elapsed_s": s["elapsed"], "cell": i,
                "z_bottom_m": s["faces"][i] * s["h"],
                "z_top_m": s["faces"][i + 1] * s["h"],
                "z_center_m": s["centers"][i] * s["h"], "z_over_h": s["centers"][i],
                "epsilon": s["eps"][i], "c": s["c"][i], "c_minus_ceq": s["c"][i] - CEQ,
                "reaction_kg_m3_s": s["reaction"][i], "reaction_kg_m2_s": s["reaction"][i] * s["volume"][i],
                "deff_m2_s": s["deff"][i], "av_m_inv": s["av"][i],
                "ell_m": ell[i], "h_over_ell": s["h"] / ell[i],
            })
    return face_rows, reaction_rows


def merge_step_budgets(metrics, case_dir):
    path = pathlib.Path(case_dir) / "salt_budget.csv"
    if not path.exists():
        return metrics
    budget = {round(row["time_s"], 9): row for row in read_csv(path) if row["phase"] == "active"}
    for row in metrics:
        b = budget.get(round(row["time_s"], 9))
        if b:
            for key in ("pore_salt_change_kg_m2_s", "pore_water_change_kg_m2_s",
                        "salt_closure_residual_kg_m2_s", "water_closure_residual_kg_m2_s",
                        "net_solid_kg_m2_s", "ledger_solid_residual_kg_m2_s"):
                row[key] = b[key]
    return metrics


def compact_local_balance(case_dir, selected_times):
    path = pathlib.Path(case_dir) / "ale_local_balance.csv"
    if not path.exists():
        return [], {"available": False}
    rows = read_csv(path)
    cells_by_time = defaultdict(list)
    for row in rows:
        if row["scope"] == "cell":
            cells_by_time[row["time_s"]].append(row)
    prefix_rows = []
    for time, cells in cells_by_time.items():
        cells.sort(key=lambda row: int(row["cell"]))
        q0 = q1 = reaction = 0.0
        for cell in cells:
            q0 += cell["q_before_kg_m2"]; q1 += cell["q_after_kg_m2"]
            reaction += cell["integrated_reaction_kg_m2"]
            incoming = cells[0]["integrated_salt_in_kg_m2"]
            outgoing = cell["integrated_salt_out_kg_m2"]
            residual = (q1 - q0) - (incoming - outgoing - reaction)
            scale = max(abs(q1 - q0), abs(incoming), abs(outgoing), abs(reaction), 1e-30)
            prefix_rows.append({
                "time_s": time, "dt_s": cell["dt_s"], "scope": "prefix",
                "cell": cell["cell"], "xi_bottom": 0.0, "xi_top": cell["xi_top"],
                "z_bottom_m": 0.0, "z_top_m": cell["z_top_m"],
                "q_before_kg_m2": q0, "q_after_kg_m2": q1,
                "integrated_salt_in_kg_m2": incoming,
                "integrated_salt_out_kg_m2": outgoing,
                "integrated_reaction_kg_m2": reaction,
                "balance_residual_kg_m2": residual, "relative_residual": residual / scale,
            })
    rows.extend(prefix_rows)
    selected_times = np.array(list(selected_times))
    keep, worst = [], {}
    for row in rows:
        time = row["time_s"]
        if row["scope"] == "whole":
            copy = dict(row); copy["selection"] = "all_whole_steps"; keep.append(copy)
        else:
            key = round(time, 9)
            if key not in worst or abs(row["relative_residual"]) > abs(worst[key]["relative_residual"]):
                worst[key] = row
            if selected_times.size and np.min(np.abs(selected_times - time)) < 1e-8:
                copy = dict(row); copy["selection"] = "selected_full_profile"; keep.append(copy)
    for row in worst.values():
        copy = dict(row); copy["selection"] = "worst_cell_each_step"; keep.append(copy)
    cell_rows = [r for r in rows if r["scope"] == "cell"]
    prefix_only = [r for r in rows if r["scope"] == "prefix"]
    whole_rows = [r for r in rows if r["scope"] == "whole"]
    summary = {
        "available": True,
        "max_abs_cell_residual_kg_m2": max(abs(r["balance_residual_kg_m2"]) for r in cell_rows),
        "max_abs_cell_relative": max(abs(r["relative_residual"]) for r in cell_rows),
        "max_abs_prefix_residual_kg_m2": max(abs(r["balance_residual_kg_m2"]) for r in prefix_only),
        "max_abs_prefix_relative": max(abs(r["relative_residual"]) for r in prefix_only),
        "max_abs_whole_residual_kg_m2": max(abs(r["balance_residual_kg_m2"]) for r in whole_rows),
        "max_abs_whole_relative": max(abs(r["relative_residual"]) for r in whole_rows),
    }
    return keep, summary


def history_summary(case_dir, states, salinity):
    history = read_csv(pathlib.Path(case_dir) / "coupled_history.csv")
    first_solid = next((r["time_s"] for r in history if r["top_solid_kg_m2"] > 0), math.nan)
    active = [r for r in history if r["active"] > .5]
    first, last = active[0], active[-1]
    metrics = [state_metrics(s) for s in states]
    return {
        "row_type": "case_summary", "case": states[0]["case"], "x_nacl": salinity,
        "first_surface_precipitation_s": first_solid, "activation_s": first["activation_s"],
        "activation_pore_c": states[0]["c"].mean(),
        "activation_bottom_minus_ceq": states[0]["bottom_c"] - CEQ,
        "last_accepted_s": last["time_s"], "post_activation_lifetime_s": last["time_s"] - first["activation_s"],
        "final_h_m": last["h_m"], "final_h_mass_m": last["h_mass_m"],
        "bottom_dissolution_kg_m2": last["internal_dissolution_kg_m2"] - first["internal_dissolution_kg_m2"],
        "positive_precipitation_kg_m2": last["internal_positive_kg_m2"] - first["internal_positive_kg_m2"],
        "top_growth_kg_m2": last["top_solid_kg_m2"] - first["top_solid_kg_m2"],
        "eps_min_final": last["epsilon_min"], "eps_max_final": max(states[-1]["eps"]),
        "bottom_minus_ceq_final": states[-1]["bottom_c"] - CEQ,
        "diff_ratio_min_final": metrics[-1]["diff_ratio_min"],
        "dissolution_zone_final_m": metrics[-1]["dissolution_zone_m"],
        "salt_residual_max_abs_kg": max(abs(r["salt_residual_kg"]) for r in history),
        "water_residual_max_abs_kg": max(abs(r["water_residual_kg"]) for r in history),
    }


def comparison_rows(a_states, b_states, a_dir, b_dir):
    out = [history_summary(a_dir, a_states, .08), history_summary(b_dir, b_states, .09)]
    common_elapsed = min(a_states[-1]["elapsed"], b_states[-1]["elapsed"])
    for elapsed in [10, 100, 300, common_elapsed]:
        if elapsed <= common_elapsed + 1e-9:
            for states in (a_states, b_states):
                s = nearest(states, elapsed); m = state_metrics(s)
                m.update({"row_type": "matched_elapsed", "alignment_target": elapsed})
                out.append(m)
    h_low = max(a_states[0]["h"], b_states[0]["h"])
    h_high = min(a_states[-1]["h"], b_states[-1]["h"])
    for height in np.linspace(h_low, h_high, 3):
        for states in (a_states, b_states):
            s = min(states, key=lambda state: abs(state["h"] - height)); m = state_metrics(s)
            m.update({"row_type": "matched_height", "alignment_target": height})
            out.append(m)
    return out


def heatmap(ax, states, key, title, cmap, centered=False):
    times = np.array([s["elapsed"] for s in states])
    zmax = max(s["h"] for s in states)
    z = np.linspace(0, zmax, 160)
    field = np.full((len(states), len(z)), np.nan)
    for j, s in enumerate(states):
        values = s[key]
        centers = s["centers"] * s["h"]
        inside = z <= s["h"]
        field[j, inside] = np.interp(z[inside], centers, values, left=values[0], right=values[-1])
    if centered:
        limit = np.nanpercentile(np.abs(field), 99)
        mesh = ax.pcolormesh(times, z * 1e6, field.T, shading="auto", cmap=cmap,
                             norm=mpl.colors.TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit))
    else:
        mesh = ax.pcolormesh(times, z * 1e6, field.T, shading="auto", cmap=cmap)
    ax.plot(times, [s["h"] * 1e6 for s in states], color="k", lw=.8)
    ax.set(xlabel="Time after activation (s)", ylabel="Physical height (µm)", title=title)
    plt.colorbar(mesh, ax=ax)


def figures(out, states, selected, metrics, local_rows, case_dir):
    figdir = out / "figures"; figdir.mkdir(exist_ok=True)
    mpl.rcParams.update({"font.size": 8, "axes.grid": True, "grid.alpha": .2,
                         "figure.dpi": 160, "savefig.dpi": 300})
    shown = [nearest(states, t) for t in (0, 10, 100, min(300, states[-1]["elapsed"]), states[-1]["elapsed"])]
    shown = list({s["time"]: s for s in shown}.values())
    colors = plt.cm.viridis(np.linspace(0, 1, len(shown)))
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.1))
    for s, color in zip(shown, colors):
        z = s["centers"] * s["h"] * 1e6; label = "%g s" % s["elapsed"]
        axes[0].plot((s["c"] - CEQ) * 1e3, z, color=color, label=label)
        axes[1].plot(s["eps"], z, color=color)
        axes[2].plot(s["reaction"], z, color=color)
    axes[0].set(xlabel=r"$(C-C_{eq})\times10^3$", ylabel="Physical height (µm)")
    axes[1].set(xlabel="Porosity")
    axes[2].set(xlabel=r"Reaction (kg m$^{-3}$ s$^{-1}$)")
    axes[0].legend(frameon=False, fontsize=7)
    fig.tight_layout(); fig.savefig(figdir / "figure_1_state_profiles.png"); plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(10, 3.1))
    for s, color in zip(shown, colors):
        z = s["faces"][:-1] * s["h"] * 1e6; label = "%g s" % s["elapsed"]
        axes[0].plot(s["adv"][:-1] * 1e6, z, color=color, label=label)
        axes[1].plot(s["diff"][:-1] * 1e6, z, color=color)
        axes[2].plot(s["salt"][:-1] * 1e6, z, color=color)
    for ax, title in zip(axes, ("Advective", "Diffusive", "Net liquid salt")):
        ax.set(xlabel=title + r" ($10^{-6}$ kg m$^{-2}$ s$^{-1}$)")
    axes[0].set_ylabel("Physical height (µm)"); axes[0].legend(frameon=False, fontsize=7)
    fig.tight_layout(); fig.savefig(figdir / "figure_2_face_flux_profiles.png"); plt.close(fig)

    for index, (key, title, cmap, centered) in enumerate((
            ("diff", r"Diffusive salt flux (kg m$^{-2}$ s$^{-1}$)", "coolwarm", True),
            ("reaction", r"Reaction (kg m$^{-3}$ s$^{-1}$)", "coolwarm", True),
            ("eps", "Porosity", "viridis", False)), start=3):
        # Diffusion is face based; omit the solid-birth top face and map it to cell centers.
        mapped = []
        for s in states:
            copy = dict(s)
            if key == "diff": copy[key] = s["diff"][:-1]
            mapped.append(copy)
        fig, ax = plt.subplots(figsize=(6.3, 3.8)); heatmap(ax, mapped, key, title, cmap, centered)
        fig.tight_layout(); fig.savefig(figdir / ("figure_%d_%s_heatmap.png" %
                                                  (index, "diffusion" if index == 3 else "reaction" if index == 4 else "porosity")))
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4))
    time = np.array([m["elapsed_s"] for m in metrics])
    for key, label in (("g_top_kg_m2_s", "Top solid birth"),
                       ("r_positive_kg_m2_s", "Internal precipitation"),
                       ("r_dissolution_kg_m2_s", "Internal dissolution magnitude"),
                       ("bottom_salt_kg_m2_s", "Soil net salt supply"),
                       ("pore_salt_change_kg_m2_s", "Pore-salt storage rate")):
        values = [m.get(key, math.nan) for m in metrics]
        ax.plot(time, np.array(values) * 1e6, label=label)
    ax.set(xlabel="Time after activation (s)", ylabel=r"Rate ($10^{-6}$ kg m$^{-2}$ s$^{-1}$)")
    ax.legend(frameon=False, ncol=2); fig.tight_layout()
    fig.savefig(figdir / "figure_6_salt_rate_decomposition.png"); plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9, 3.4))
    axes[0].plot(time, np.array([m["dissolution_zone_m"] for m in metrics]) * 1e6, label="Dissolution thickness")
    axes[0].plot(time, np.array([m["reaction_sign_change_m"] for m in metrics]) * 1e6, label="First sign change")
    axes[0].plot(time, np.array([m["z_diff_min_m"] for m in metrics]) * 1e6, label="Diffusion minimum")
    axes[0].set(xlabel="Time after activation (s)", ylabel="Physical height (µm)")
    axes[0].legend(frameon=False)
    axes[1].plot(time, [m["h_over_ell_bottom"] for m in metrics], label="bottom ell")
    axes[1].plot(time, [m["h_over_ell_median"] for m in metrics], label="median ell")
    axes[1].set(xlabel="Time after activation (s)", ylabel=r"$h/\ell$")
    axes[1].legend(frameon=False); fig.tight_layout()
    fig.savefig(figdir / "figure_7_zone_and_length_scales.png"); plt.close(fig)

    history = read_csv(pathlib.Path(case_dir) / "coupled_history.csv")
    active = [r for r in history if r["active"] > .5]
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.4))
    whole = [r for r in local_rows if r["scope"] == "whole"]
    worst = [r for r in local_rows if r.get("selection") == "worst_cell_each_step"]
    if whole: axes[0].semilogy([r["time_s"] - states[0]["time"] for r in whole],
                               np.maximum(np.abs([r["relative_residual"] for r in whole]), 1e-20), label="whole crust")
    if worst: axes[0].semilogy([r["time_s"] - states[0]["time"] for r in worst],
                               np.maximum(np.abs([r["relative_residual"] for r in worst]), 1e-20), label="worst cell")
    axes[0].set(xlabel="Time after activation (s)", ylabel="Accepted-step relative residual")
    axes[0].legend(frameon=False)
    axes[1].semilogy([r["time_s"] for r in active],
                     np.maximum(np.abs([r["salt_residual_kg"] for r in active]), 1e-20), label="salt")
    axes[1].semilogy([r["time_s"] for r in active],
                     np.maximum(np.abs([r["water_residual_kg"] for r in active]), 1e-20), label="water")
    axes[1].set(xlabel="Absolute time (s)", ylabel="Global residual (kg)")
    axes[1].legend(frameon=False); fig.tight_layout()
    fig.savefig(figdir / "figure_8_conservation_residuals.png"); plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--x008-progressive", required=True)
    parser.add_argument("--x008-adaptive", required=True)
    parser.add_argument("--x009", required=True)
    parser.add_argument("--x009-half", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    out = pathlib.Path(args.out); out.mkdir(parents=True, exist_ok=True)

    registered = {
        "x008_progressive": pathlib.Path(args.x008_progressive),
        "x008_adaptive_reference": pathlib.Path(args.x008_adaptive),
        "x009_progressive": pathlib.Path(args.x009),
        "x009_dt_half": pathlib.Path(args.x009_half),
    }
    states, rejected, verification = {}, {}, {"flux_reconstruction": {}, "local_balance": {}}
    for name, path in registered.items():
        states[name], rejected[name] = accepted_states(path, name)
        verification["flux_reconstruction"][name] = validation(states[name])
        check = verification["flux_reconstruction"][name]
        # Relative errors are ill-conditioned when the physical flux is near
        # zero at activation. Require either a strict absolute or relative gate
        # for each independently reconstructed quantity.
        check["pass_abs_1e-18_or_relative_1e-10"] = all(
            check["max_abs"][key] <= 1e-18 or check["max_rel"][key] <= 1e-10
            for key in check["max_abs"])

    selected = {name: selected_states(states[name], rejected[name]) for name in states}
    face_rows, reaction_rows = [], []
    for name in ("x008_adaptive_reference", "x008_progressive", "x009_progressive"):
        f, r = profile_rows(selected[name]); face_rows.extend(f); reaction_rows.extend(r)
    write_csv(out / "FACE_FLUX_PROFILES.csv", face_rows)
    write_csv(out / "REACTION_PROFILES.csv", reaction_rows)

    time_series = []
    for name in states:
        metrics = [state_metrics(state) for state in states[name]]
        merge_step_budgets(metrics, registered[name])
        time_series.extend(metrics)
    write_csv(out / "MECHANISM_TIME_SERIES.csv", time_series)

    compact, local_summary = compact_local_balance(
        registered["x009_progressive"], [s["time"] for s in selected["x009_progressive"]])
    verification["local_balance"]["x009_progressive"] = local_summary
    write_csv(out / "LOCAL_STORAGE_BALANCE.csv", compact)
    comparison = comparison_rows(states["x008_progressive"], states["x009_progressive"],
                                 registered["x008_progressive"], registered["x009_progressive"])
    write_csv(out / "X008_X009_COMPARISON.csv", comparison)

    main_metrics = [row for row in time_series if row["case"] == "x009_progressive"]
    figures(out, states["x009_progressive"], selected["x009_progressive"], main_metrics,
            compact, registered["x009_progressive"])

    # Timestep check at common elapsed times, without extrapolation.
    common = min(states["x009_progressive"][-1]["elapsed"], states["x009_dt_half"][-1]["elapsed"])
    checks = []
    for elapsed in (10, 100, 300, common):
        if elapsed <= common + 1e-9:
            a = state_metrics(nearest(states["x009_progressive"], elapsed))
            b = state_metrics(nearest(states["x009_dt_half"], elapsed))
            checks.append({"elapsed_s": elapsed,
                           "height_relative_difference": abs(a["h_m"] - b["h_m"]) / max(abs(b["h_m"]), 1e-30),
                           "g_top_relative_difference": abs(a["g_top_kg_m2_s"] - b["g_top_kg_m2_s"]) / max(abs(b["g_top_kg_m2_s"]), 1e-30),
                           "dissolution_relative_difference": abs(a["r_dissolution_kg_m2_s"] - b["r_dissolution_kg_m2_s"]) / max(abs(b["r_dissolution_kg_m2_s"]), 1e-30),
                           "eps_min_absolute_difference": abs(a["eps_min"] - b["eps_min"])})
    verification["timestep_check"] = checks
    case_spec = json.loads((registered["x009_progressive"] / "case_spec.json").read_text())
    reservoir = read_csv(registered["x009_progressive"] / "reservoir_config.csv")[0]
    verification["x009_input_checks"] = {
        "problem_salt_mole_fraction": case_spec["salt_mole_fraction"],
        "configured_reservoir_mole_fraction": reservoir["configured_salt_mole_fraction"],
        "reported_reservoir_mole_fraction": reservoir["reservoir_salt_mole_fraction"],
        "accepted_initial_min_mole_fraction": reservoir["initial_min_salt_mole_fraction"],
        "accepted_initial_max_mole_fraction": reservoir["initial_max_salt_mole_fraction"],
        "all_equal_0p09": all(abs(value - .09) <= 1e-14 for value in (
            case_spec["salt_mole_fraction"], reservoir["configured_salt_mole_fraction"],
            reservoir["reservoir_salt_mole_fraction"], reservoir["initial_min_salt_mole_fraction"],
            reservoir["initial_max_salt_mole_fraction"])),
    }
    verification["selected_states"] = {name: [{"selection": s["selection"], "status": s["status"],
                                                "time_s": s["time"], "elapsed_s": s["elapsed"]}
                                               for s in values] for name, values in selected.items()}
    verification["raw_files"] = {name: {filename: {"sha256": sha256(path / filename),
                                                     "bytes": (path / filename).stat().st_size}
                                                for filename in ("ale_diagnostics.csv", "coupled_history.csv",
                                                                 "ale_local_balance.csv", "salt_budget.csv",
                                                                 "outcome.json", "case_spec.json", "command.json",
                                                                 "params.input", "reservoir_config.csv")
                                                if (path / filename).exists()}
                                 for name, path in registered.items()}
    (out / "verification.json").write_text(json.dumps(verification, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"out": str(out), "cases": {name: len(value) for name, value in states.items()},
                      "local_balance": local_summary}, indent=2))


if __name__ == "__main__":
    main()
