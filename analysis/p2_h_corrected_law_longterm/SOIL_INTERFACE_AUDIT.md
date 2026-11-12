# P2: native boundary versus internal SCVF and ALE demand

This is an interface-consistency diagnostic using frozen source exports. It is
not a new DuMuX solution or a hydraulic coupling. The original CSV explicitly
labels its flux location `top_internal_scvf`. No physical face coordinates or
matching crust interface geometry are included in the frozen files. Therefore
that internal flux must not be relabeled as a soil/crust interface flux.

The independent native surface salt output is available as
`j_precip_kg_m2_s` in x008_timeseries.csv (also exported as
j_s_precip_boundary_kg_m2_s in flux_decomposition.csv). Its balance is

    J_internal - dM_top_control_volume/dt = J_native_surface_salt.

On the 44.70875 h source support, the independently tabulated terms give:

| Quantity | Integral (kg/m2) |
|---|---:|
| First internal SCVF net salt | 0.5949903634 |
| Native surface salt output | 0.5949832275 |
| Top control-volume dissolved storage | 0.0000071369 |
| Analytical reconstruction from j, B and Ceq | 0.6063993741 |

The pointwise internal/control-volume/native-boundary balance residual is at
most 2.56e-13 kg/m2/s. Thus the internal face is distinct from the boundary,
but the intervening storage accounts for only 7.14e-6 kg/m2. It does NOT explain
the approximately 0.011409 kg/m2 discrepancy between the analytical
reconstruction and the internal SCVF transport. Earlier wording that attributed
the entire 1.92% discrepancy to face location would be unsupported.

The original data provide internal advective/diffusive salt fluxes, native
surface water and salt outflows, and native top concentration/saturation.
They do not provide the actual advective donor concentration and liquid flux at
that internal face, nor a resolved common soil/crust face. Consequently the
remaining 1.92% mismatch cannot be uniquely assigned to concentration,
mass-average/Fick conventions or discretization from these exports alone.
The file soil_native_boundary_audit.csv retains the actual exported quantities.

For the new ALE Neumann experiment, the entire usable support is only 900 s,
because the bottom crystal fraction reaches the structural screen. Over that
common interval the accepted-step ALE salt influx is 0.00150662 kg/m2. ALE
bottom water demand is 0.0185280 kg/m2 while evaporation removes 0.0160342 kg/m2;
the additional 0.00249376 kg/m2 remains in the changing pore-water inventory.
The approximately 15.55% additional water demand cannot be ignored when passing
the bare-soil evaporation history to a saturated crust. It is not proof that
the soil cannot supply it; it shows that available water supply has not been
demonstrated by prescribing evaporation alone.

The prior 0.698607 kg/m2 ALE influx was integrated from the 44.70875 h
equilibrium-inflow, liquid-gap-permitting trajectory. It is not the cumulative
influx of the new structurally censored Neumann run. No extension of the latter
to 44.70875 h is claimed. Native top C and ALE reconstructed bottom-face C are
reported separately in soil_ale_common_window.csv; no shared interface
concentration has been solved.

P2 stops here as required by task 26: a quantitative hydraulic-compatibility
claim needs a common geometric interface, soil liquid flux and donor state,
crust pressure/permeability or a justified supply law, and a conservative joint
interface solve. A numerical correction factor would not supply those missing
conditions. No such factor or incompatible second water-flux condition was used.
