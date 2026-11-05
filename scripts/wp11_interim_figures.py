"""Focused figures for the WP11 R96 interim three-grid assessment."""
import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from wp11_field_metrics import cell_grid, parse_vtu, point_grid

T_S = 545692.2550893936 / 4.0
COLORS = {"R24":"#0072B2", "R48":"#D55E00", "R96":"#009E73"}
STYLES = {"R24":"-", "R48":"--", "R96":"-."}


def rows(path):
    with Path(path).open(newline="") as handle: return list(csv.DictReader(handle))


def save(fig, out, stem):
    fig.savefig(out/f"{stem}.png", dpi=220, bbox_inches="tight")
    fig.savefig(out/f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def cumulative(records, field, positive=False):
    y=[0.0]
    for a,b in zip(records, records[1:]):
        value=float(b[field]); value=max(value,0.0) if positive else value
        y.append(y[-1]+value*(float(b["time_s"])-float(a["time_s"])))
    return np.asarray(y)


def load_panel(path):
    points, point_x, _, centers, _, velocity, _ = parse_vtu(path)
    xs,zs,x=point_grid(points,point_x)
    cxs,czs,down=cell_grid(centers,-velocity[:,1])
    return xs,zs,x,cxs,czs,down


def find_frame(root, frame):
    matches=list(Path(root).glob(f"*-{frame:05d}.vtu"))
    if len(matches)!=1: raise AssertionError(f"expected one frame {frame} under {root}, got {matches}")
    return matches[0]


def main(out, r24, r48, r96, raw24, raw48, raw96, end_T):
    out=Path(out); out.mkdir(parents=True,exist_ok=True)
    reduced={"R24":Path(r24),"R48":Path(r48),"R96":Path(r96)}
    surface={k:rows(v/"surface_series.csv") for k,v in reduced.items()}
    flux={k:rows(v/"flux_series.csv") for k,v in reduced.items()}
    modes={k:rows(v/"tracked_mode_series.csv") for k,v in reduced.items()}

    fig,axes=plt.subplots(2,1,figsize=(8.0,7.2),sharex=True)
    for grid in COLORS:
        q=[r for r in surface[grid] if float(r["time_over_T"])<=end_T+1e-9]
        t=[float(r["time_over_T"]) for r in q]
        axes[0].plot(t,[float(r["surface_mean"]) for r in q],color=COLORS[grid],ls=STYLES[grid],label=grid,lw=1.6)
        axes[1].plot(t,[float(r["surface_max"]) for r in q],color=COLORS[grid],ls=STYLES[grid],label=grid,lw=1.6)
    for ax in axes: ax.axvline(7.18,color="0.45",ls=":",lw=1,label="Paper 1 base-state context" if ax is axes[0] else None); ax.grid(alpha=.22)
    axes[1].axhline(.1,color="black",ls=":",lw=1,label="solubility threshold")
    axes[0].set_ylabel("surface mean mole fraction"); axes[1].set_ylabel("surface maximum mole fraction"); axes[1].set_xlabel("time / T")
    axes[0].legend(ncol=4,fontsize=8); axes[1].legend(fontsize=8); axes[0].set_title("A. Three-grid surface concentration (common interim window)")
    save(fig,out,"Figure_A_surface_concentration")

    depths=(1.,2.,5.,10.); fig,axes=plt.subplots(2,2,figsize=(10.2,7.0),sharex=True)
    for ax,depth in zip(axes.flat,depths):
        for grid in COLORS:
            q=[r for r in flux[grid] if float(r["depth_mm"])==depth and float(r["time_over_T"])<=end_T+1e-9]
            ax.plot([float(r["time_over_T"]) for r in q],[float(r["advective_down_m_s"]) for r in q],color=COLORS[grid],ls=STYLES[grid],lw=1.15,label=grid)
        ax.axhline(0,color="black",lw=.7); ax.axhline(1e-15,color="0.5",ls=":",lw=.7); ax.set_title(f"depth = {depth:g} mm"); ax.grid(alpha=.2)
    axes[0,0].legend(); axes[1,0].set_xlabel("time / T"); axes[1,1].set_xlabel("time / T"); axes[0,0].set_ylabel("advective flux (m s$^{-1}$)"); axes[1,0].set_ylabel("advective flux (m s$^{-1}$)")
    fig.suptitle("B. Signed conservative downward advective salt flux")
    save(fig,out,"Figure_B_advective_flux")

    fig,axes=plt.subplots(2,2,figsize=(10.2,7.0),sharex=True)
    for ax,depth in zip(axes.flat,depths):
        for grid in COLORS:
            q=[r for r in flux[grid] if float(r["depth_mm"])==depth and float(r["time_over_T"])<=end_T+1e-9]
            t=[float(r["time_over_T"]) for r in q]
            ax.plot(t,cumulative(q,"advective_down_m_s",True),color=COLORS[grid],ls=STYLES[grid],lw=1.4,label=f"{grid} adv+" )
            ax.plot(t,cumulative(q,"total_down_m_s",False),color=COLORS[grid],ls=":",lw=1.0,label=f"{grid} total signed")
        ax.axhline(0,color="black",lw=.7); ax.set_title(f"depth = {depth:g} mm"); ax.grid(alpha=.2)
    axes[0,0].legend(ncol=2,fontsize=7); axes[1,0].set_xlabel("time / T"); axes[1,1].set_xlabel("time / T"); axes[0,0].set_ylabel("cumulative transport (m)"); axes[1,0].set_ylabel("cumulative transport (m)")
    fig.suptitle("C. Positive-only advection (grid style) and signed total transport (dotted)")
    save(fig,out,"Figure_C_cumulative_transport")

    raw={"R24":Path(raw24),"R48":Path(raw48),"R96":Path(raw96)}; times=(5,6,8,20)
    panels={}
    for grid in COLORS:
        for time in times: panels[(grid,time)]=load_panel(find_frame(raw[grid],int(time*24)))
    cmin=min(p[2].min() for p in panels.values()); cmax=max(p[2].max() for p in panels.values())
    vmax=max(np.quantile(np.abs(p[5]),.995) for p in panels.values())
    fig,axes=plt.subplots(6,4,figsize=(12.5,13.5),sharex=False,sharey=True,constrained_layout=True)
    cim=vim=None
    for gi,grid in enumerate(COLORS):
        for j,time in enumerate(times):
            xs,zs,x,cxs,czs,down=panels[(grid,time)]
            a=axes[2*gi,j]; b=axes[2*gi+1,j]
            cim=a.pcolormesh(xs*1000,(zs.max()-zs)*1000,x,shading="auto",vmin=cmin,vmax=cmax,cmap="viridis")
            vim=b.pcolormesh(cxs*1000,(czs.max()-czs)*1000,down,shading="auto",vmin=-vmax,vmax=vmax,cmap="RdBu_r")
            # All resolved fingers remain within the upper 13 mm in this window;
            # a common 0--30 mm crop makes the three-grid comparison legible.
            a.set_ylim(30,0); b.set_ylim(30,0)
            if gi==0: a.set_title(f"{time}T")
            if j==0: a.set_ylabel(f"{grid}\nsalt depth (mm)"); b.set_ylabel(f"{grid}\nvelocity depth (mm)")
            b.set_xlabel("x (mm)")
    fig.colorbar(cim,ax=axes[0::2,:],shrink=.55,label="salt mole fraction"); fig.colorbar(vim,ax=axes[1::2,:],shrink=.55,label="downward cell-centre velocity (m/s; qualitative)")
    fig.suptitle("D. Matched-time fingering and qualitative reconstructed velocity",y=1.02)
    save(fig,out,"Figure_D_fingering_fields")

    fig,axes=plt.subplots(2,2,figsize=(10.0,7.0),sharex=True,sharey=True)
    for ax,mode in zip(axes.flat,(2,4,6,10)):
        for grid in COLORS:
            q=[r for r in modes[grid] if int(r["mode"])==mode and float(r["time_over_T"])<=end_T+1e-9]
            marker="o" if grid=="R96" else None
            ax.plot([float(r["time_over_T"]) for r in q],[float(r["energy_fraction"]) for r in q],color=COLORS[grid],ls=STYLES[grid],marker=marker,ms=3,lw=1.35,label=grid)
        ax.set_title(f"mode m={mode}"); ax.grid(alpha=.2)
    axes[0,0].legend(); axes[1,0].set_xlabel("time / T"); axes[1,1].set_xlabel("time / T"); axes[0,0].set_ylabel("energy fraction"); axes[1,0].set_ylabel("energy fraction")
    fig.suptitle("E. Modal development (R96 markers denote sparse saved fields)")
    save(fig,out,"Figure_E_modal_development")

    manifest={"common_end_T":end_T,"figures":["Figure_A_surface_concentration","Figure_B_advective_flux","Figure_C_cumulative_transport","Figure_D_fingering_fields","Figure_E_modal_development"],
        "styles":{"R24":"blue solid","R48":"vermillion dashed","R96":"green dash-dot"},
        "scope":"interim common-window figures; no smoothing; R96 spatial/modal series is sparse"}
    (out/"figure_manifest.json").write_text(json.dumps(manifest,indent=2))


if __name__=="__main__":
    p=argparse.ArgumentParser(); p.add_argument("out"); p.add_argument("r24"); p.add_argument("r48"); p.add_argument("r96"); p.add_argument("raw24"); p.add_argument("raw48"); p.add_argument("raw96"); p.add_argument("--end-T",type=float,required=True)
    a=p.parse_args(); main(a.out,a.r24,a.r48,a.r96,a.raw24,a.raw48,a.raw96,a.end_T)
