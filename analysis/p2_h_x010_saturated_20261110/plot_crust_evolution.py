"""Plot accepted ALE fields from HPCC job 16481021."""

import hashlib
import json
from pathlib import Path

import matplotlib as mpl
import numpy as np

mpl.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent
RAW = ROOT / "raw"
CEQ = 0.26494086976769066


def verified_source(name):
    expected = {
        name: digest
        for digest, name in (
            line.split() for line in (RAW / "output_sha256.txt").read_text().splitlines()
        )
    }
    path = RAW / name
    with path.open("rb") as stream:
        digest = hashlib.file_digest(stream, "sha256").hexdigest()
    if digest != expected[name]:
        raise ValueError(f"SHA-256 mismatch for {path}")
    return path


def read_fields():
    spec = json.loads(verified_source("case_spec.json").read_text())
    outcome = json.loads(verified_source("outcome.json").read_text())
    data = np.loadtxt(verified_source("crust_fields.csv"), delimiter=",", skiprows=1)
    blocks = data.reshape(-1, spec["ale_cells"], 5)
    times = blocks[:, 0, 0]
    lower, upper, porosity, concentration = np.moveaxis(blocks[:, :, 1:], -1, 0)
    if not (
        np.all(np.isfinite(blocks))
        and np.all(blocks[:, :, 0] == times[:, None])
        and np.all(np.diff(times) > 0)
        and times[0] == outcome["activation_time_s"]
        and times[-1] == outcome["last_accepted_time_s"]
        and np.all(lower[:, 0] == 0)
        and np.array_equal(lower[:, 1:], upper[:, :-1])
        and np.all(upper > lower)
    ):
        raise ValueError("Incomplete or inconsistent ALE field sequence")

    history = np.genfromtxt(
        verified_source("coupled_history.csv"), delimiter=",", names=True
    )
    at_snapshot = history[np.searchsorted(history["time_s"], times)]
    np.testing.assert_array_equal(at_snapshot["time_s"], times)
    np.testing.assert_allclose(at_snapshot["h_m"], upper[:, -1], rtol=0, atol=1e-15)
    np.testing.assert_allclose(
        at_snapshot["epsilon_min"], porosity.min(axis=1), rtol=0, atol=1e-12
    )
    faces = np.column_stack([lower[:, 0], upper]) * 1e6
    return times, faces, porosity, (concentration - CEQ) * 1e6


def midpoint_edges(samples):
    """Build display edges from sample centers, keeping both endpoints exact."""
    return np.concatenate([samples[:1], (samples[:-1] + samples[1:]) / 2, samples[-1:]])


def main():
    time, faces, porosity, concentration = read_fields()
    time_edges = midpoint_edges(time)
    height_edges = midpoint_edges(faces)
    top_edges = height_edges[:, -1]
    depth_edges = top_edges[:, None] - height_edges
    mesh_time = np.broadcast_to(time_edges[:, None], height_edges.shape)
    blues = mpl.colors.ListedColormap(mpl.colormaps["Blues"](np.linspace(0.18, 0.94, 256)))

    style = {
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "figure.facecolor": "white",
    }
    with mpl.rc_context(style):
        fig = plt.figure(figsize=(11.5, 6.8))
        grid = fig.add_gridspec(
            2, 3, width_ratios=[1.2, 1, 0.04], hspace=0.40, wspace=0.28,
            left=0.09, right=0.85, bottom=0.16, top=0.84,
        )
        fig.text(0.09, 0.95, "Porosity and dissolved salt in the growing crust",
                 fontsize=17, weight="bold")
        fig.text(0.09, 0.902, "Saturated supply (x = 0.10)  |  80 soil layers  |  "
                 "ALE activation at 673 s", color="#454545", fontsize=10.5)
        settings = [
            (porosity, "viridis", 0.15, 0.43, [0.15, 0.22, 0.29, 0.36, 0.43],
             r"Porosity, $\epsilon$ (–)", "Porosity", "a", "b"),
            (concentration, blues, min(0, concentration.min()), 3.2, [0, 0.8, 1.6, 2.4, 3.2],
             r"$C-C_{\mathrm{eq}}$ (mg NaCl / kg solution)", "Dissolved NaCl", "c", "d"),
        ]
        for row, (values, cmap, vmin, vmax, ticks, label, title, left, right) in enumerate(settings):
            full = fig.add_subplot(grid[row, 0])
            zoom = fig.add_subplot(grid[row, 1], sharex=full)
            norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
            for ax, coordinates in ((full, height_edges), (zoom, depth_edges)):
                ax.pcolormesh(mesh_time, coordinates, values, cmap=cmap, norm=norm,
                              shading="flat", rasterized=True)
                ax.set_xlim(0, time[-1])
                ax.set_xticks([0, 2500, 5000, 7500, 10000, 12500])
                ax.axvline(time[0], color="#777777", linewidth=0.7, linestyle="--")
                if row == 1:
                    ax.set_xlabel("Time (s)")
                else:
                    ax.tick_params(labelbottom=False)
            full.plot(time_edges, top_edges, color="#242424", lw=0.85)
            full.set(ylim=(0, 48), ylabel=r"Height above soil ($\mu$m)")
            zoom.set(ylim=(2, 0), ylabel=r"Depth below moving top ($\mu$m)")
            full.set_title(f"{left}   {title} · full crust", loc="left", weight="bold")
            zoom.set_title(f"{right}   Top 2 µm", loc="left", weight="bold")
            fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
                         cax=fig.add_subplot(grid[row, 2]), ticks=ticks, label=label)

        fig.text(0.09, 0.071, "Each time band uses one saved ALE profile. White areas lie outside the active crust.",
                 fontsize=9, color="#454545")
        fig.text(0.09, 0.040, r"$C_{\mathrm{eq}}$ = 0.26494087 kg NaCl / kg solution. "
                 "The right panels follow the moving top.", fontsize=9, color="#454545")

        output = ROOT / "figures"
        output.mkdir(exist_ok=True)
        for extension in ("png", "pdf"):
            fig.savefig(output / f"crust_porosity_salinity_evolution.{extension}", dpi=300)
        plt.close(fig)
    print(f"{len(time)} profiles x {porosity.shape[1]} cells, {time[0]:g}–{time[-1]:g} s")
    print(f"Porosity {porosity.min():.6f}–{porosity.max():.6f}, "
          f"excess NaCl {concentration.min():.3g}–{concentration.max():.4f} mg/kg")


if __name__ == "__main__":
    main()
