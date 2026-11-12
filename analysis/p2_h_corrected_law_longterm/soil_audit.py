"""Audit existing native boundary and internal-face exports without inventing hydraulics."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
R=Path(__file__).resolve().parent;P=R.parent/'p2_h_corrected_vs_ale_x008_20261112'
forcing=pd.read_csv(P/'results/x008_forcing.csv')
soil=pd.read_csv(P/'inputs/x008_timeseries.csv')
flux=pd.read_csv(P/'inputs/flux_decomposition.csv');flux=flux[np.isclose(flux.x0_mole_fraction,.08)]
t=forcing.elapsed_time_s.to_numpy();absolute=forcing.absolute_time_s.to_numpy()
out=forcing.copy()
for key in ['top_c_massfrac','top_sw','j_precip_kg_m2_s']:
    out['native_soil_'+key]=np.interp(absolute,soil.time_s,soil[key])
out['top_storage_rate']=np.interp(absolute,flux.time_s,flux.top_storage_rate_kg_m2_s)
out['internal_to_boundary_balance_residual']=out.soil_net_nacl_kg_m2_s-out.top_storage_rate-out.native_soil_j_precip_kg_m2_s
out['reconstructed_supply']=(.264940869768*out.evaporation_water_kg_m2_s-out.b_soil_kg_m2_s)/(1-.264940869768)
out['native_boundary_total_mass_loss']=out.evaporation_water_kg_m2_s+out.native_soil_j_precip_kg_m2_s
out.to_csv(R/'results/soil_native_boundary_audit.csv',index=False)
ale=pd.read_csv(R/'raw/x008_n64_p1.5_r1e-08_dt30_history.csv')
common=ale.copy()
for key in ['native_soil_top_c_massfrac','native_soil_top_sw','native_soil_j_precip_kg_m2_s','native_boundary_total_mass_loss','soil_net_nacl_kg_m2_s','soil_advective_nacl_kg_m2_s','b_soil_kg_m2_s','reconstructed_supply']:
    common[key]=np.interp(ale.t,t,out[key])
common['ALE_bottom_water']=common.Fb-common.Jb
common['native_soil_water_loss']=np.interp(ale.t,t,out.evaporation_water_kg_m2_s)
common.to_csv(R/'results/soil_ale_common_window.csv',index=False)
summary={key:float(np.trapezoid(out[col],t)) for key,col in {
 'internal_scvf_integral':'soil_net_nacl_kg_m2_s','native_surface_salt_integral':'native_soil_j_precip_kg_m2_s',
 'top_storage_integral':'top_storage_rate','analytical_reconstruction_integral':'reconstructed_supply'}.items()}
summary['pointwise_internal_boundary_balance_residual_max']=float(out.internal_to_boundary_balance_residual.abs().max())
summary['new_ALE_support_s']=float(ale.t.max())
summary['common_window_ALE_bottom_salt_accepted_integral']=float(ale.Qb.iloc[-1])
summary['common_window_ALE_water_demand_integral']=float(ale.Qliquid.iloc[-1]-ale.Qb.iloc[-1])
summary['common_window_evaporation_integral']=float(ale.Qevap.iloc[-1])
summary['interface_geometry']='export identifies top_internal_scvf; coordinates not provided'
summary['hydraulic_compatibility']='NOT_ESTABLISHED: native bare-soil boundary is not a coupled soil-crust interface'
(R/'results/soil_audit.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
print(json.dumps(summary,indent=2))
