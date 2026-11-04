"""Sparse-field geometry and qualitative velocity diagnostics for WP10."""
import csv, hashlib, json, sys, xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
from single_wavelength_post import widths
from wp6_4t_metrics import energy_depth, neumann_spectrum

T = 545692.2550893936/4

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main(run, out):
    run, out = Path(run), Path(out); out.mkdir(parents=True, exist_ok=True)
    prof = np.load(run/"post/profiles.npz")
    times, zs, means, sigmas = [prof[k] for k in ("time_s","z_m","x_mean","x_sigma")]
    stop_T = times[-1]/T
    targets = [2,3,4,5] + list(range(10,int(stop_T)+1,5)) + [stop_T]
    profiles=[]
    for q in targets:
        i=int(np.argmin(abs(times-q*T)))
        profiles.append({"nominal_T":q,"frame":i,"actual_T":float(times[i]/T),
            "depth90_mm":energy_depth(zs,sigmas[i],.9),"depth99_mm":energy_depth(zs,sigmas[i],.99),
            "surface_mean":float(means[i,-1]),"bottom_mean":float(means[i,0]),
            "mean_at_1mm":float(means[i,np.argmin(abs((zs[-1]-zs)-.001))]),
            "mean_at_2mm":float(means[i,np.argmin(abs((zs[-1]-zs)-.002))]),
            "mean_at_5mm":float(means[i,np.argmin(abs((zs[-1]-zs)-.005))]),
            "mean_at_10mm":float(means[i,np.argmin(abs((zs[-1]-zs)-.010))]),
            "mean_at_20mm":float(means[i,np.argmin(abs((zs[-1]-zs)-.020))])})

    rep=np.load(run/"post/representative_fields.npz");xs=rep["xs"];wx=widths(xs);W=wx.sum()
    modes=[]
    for key in sorted(k for k in rep.files if k.startswith("x_")):
        i=int(key.split("_")[1]);x=rep[key];d=x-(x@wx/W)[:,None];coeff,en=neumann_spectrum(xs,zs,d)
        m=int(np.argmax(en)+1); total=float(en.sum())
        modes.append({"frame":i,"time_over_T":float(times[i]/T),"dominant_mode":m,
            "dominant_wavelength_m":float(2*W/m),"cells_per_wavelength":float(2*(len(xs)-1)/m),
            "dominant_energy_fraction":float(en[m-1]/total),
            "top_three_modes":[int(v+1) for v in np.argsort(en)[::-1][:3]]})

    result={"profile_checkpoints":profiles,"representative_modes":modes,
        "velocity":None,"velocity_scope":"pending field archive; cell-centre VTK reconstruction is qualitative only"}
    fields=sorted(run.glob("*.vtu"))
    if fields:
        diag=list(csv.DictReader((run/"post/weighted_diagnostics.csv").open()))
        manifest=json.loads((run/"post/field_hashes.json").read_text())
        assert len(fields)==len(diag)==len(manifest)
        mode_series=[]
        for i,path in enumerate(fields):
            assert sha(path)==manifest[i]["sha256"]
            root=ET.fromstring(path.read_bytes())
            arrays={a.attrib.get("Name"):np.fromstring(a.text or "",sep=" ") for a in root.iter("DataArray")}
            xyz=arrays["Coordinates"].reshape(-1,3);ux=np.unique(xyz[:,0]);uz=np.unique(xyz[:,1])
            order=np.argsort(np.searchsorted(uz,xyz[:,1])*len(ux)+np.searchsorted(ux,xyz[:,0]))
            x=arrays["moleFraction"][order].reshape(len(uz),len(ux));xw=widths(ux)
            _,en=neumann_spectrum(ux,uz,x-(x@xw/xw.sum())[:,None]);m=int(np.argmax(en)+1)
            mode_series.append({"frame":i,"step":int(diag[i]["step"]),"time_over_T":float(diag[i]["time_s"])/T,
                                "dominant_mode":m,"cells_per_wavelength":float(2*(len(ux)-1)/m)})
        windows=[]
        bounds=[(0,4),(4,5),(5,10),(10,15),(15,20),(20,stop_T+1e-9)]
        for a,b in bounds:
            mm=[r["dominant_mode"] for r in mode_series if a<=r["time_over_T"]<b]
            counts={str(m):mm.count(m) for m in sorted(set(mm))}
            windows.append({"start_T":a,"end_T":min(b,stop_T),"frames":len(mm),"mode_counts":counts})
        result["mode_timeseries"]=mode_series
        result["mode_windows"]=windows
        selected={int(np.argmin(abs(np.array([float(r['time_s']) for r in diag])-q*T))):q for q in targets}
        velocity=[]
        for i,q in selected.items():
            path=run/Path(manifest[i]["path"]).name
            root=ET.fromstring(path.read_bytes())
            arrays={a.attrib.get("Name"):(a,np.fromstring(a.text or "",sep=" ")) for a in root.iter("DataArray")}
            names=sorted(k for k in arrays if k)
            vname=next((k for k in names if "velocity" in k.lower()),None)
            if not vname: raise AssertionError(f"no velocity array in {names}")
            xyz=arrays["Coordinates"][1].reshape(-1,3);ux=np.unique(xyz[:,0]);uz=np.unique(xyz[:,1])
            order=np.argsort(np.searchsorted(uz,xyz[:,1])*len(ux)+np.searchsorted(ux,xyz[:,0]))
            x=arrays["moleFraction"][1][order].reshape(len(uz),len(ux))
            xw=widths(ux); xd=x-(x@xw/xw.sum())[:,None]
            _,mode_energy=neumann_spectrum(ux,uz,xd)
            dominant=int(np.argmax(mode_energy)+1)
            xc=.25*(x[:-1,:-1]+x[1:,:-1]+x[:-1,1:]+x[1:,1:])
            vel=arrays[vname][1].reshape(-1,int(arrays[vname][0].attrib.get("NumberOfComponents",3)))
            assert len(vel)==(len(ux)-1)*(len(uz)-1)
            down=(-vel[:,1]).reshape(len(uz)-1,len(ux)-1)
            high=xc>=np.quantile(xc,.9)
            xanom=xc-xc.mean(axis=1,keepdims=True)
            danom=down-down.mean(axis=1,keepdims=True)
            high_anom=xanom>=np.quantile(xanom,.9)
            positive_anom=xanom>0
            velocity.append({"nominal_T":q,"frame":i,"actual_T":float(diag[i]["time_s"])/T,
                "array_name":vname,"down_min_m_s":float(down.min()),"down_max_m_s":float(down.max()),
                "dominant_mode":dominant,"dominant_wavelength_m":float(2*xw.sum()/dominant),
                "dominant_cells_per_wavelength":float(2*(len(ux)-1)/dominant),
                "top_three_modes":[int(v+1) for v in np.argsort(mode_energy)[::-1][:3]],
                "down_cell_fraction":float(np.mean(down>0)),
                "high_salt_down_fraction":float(np.mean(down[high]>0)),
                "high_salt_mean_down_m_s":float(np.mean(down[high])),
                "other_mean_down_m_s":float(np.mean(down[~high])),
                "salt_down_correlation":float(np.corrcoef(xc.ravel(),down.ravel())[0,1]),
                "horizontal_anomaly_down_correlation":float(np.corrcoef(xanom.ravel(),danom.ravel())[0,1]),
                "top_anomaly_down_fraction":float(np.mean(down[high_anom]>0)),
                "positive_anomaly_down_fraction":float(np.mean(down[positive_anom]>0)),
                "source_sha256":sha(path)})
        result["velocity"]=velocity
        result["velocity_scope"]="qualitative cell-centre Darcy reconstruction only; not a conservative face flux"
    (out/"field_metrics.json").write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))

if __name__=="__main__":main(*sys.argv[1:])
