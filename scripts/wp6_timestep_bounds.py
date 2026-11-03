"""Bounds, not interpolated cumulative values, for asynchronous report times.

Applies to the backward-Euler accepted-step conservative quadrature already
independently reintegrated by wp6_post. It is not a bound on unknown continuous
physical flux fluctuations or on spatial discretization error.
"""
import math


def local_tail_bound(previous, current, report_time):
    t0=float(previous['time_s']);t1=float(current['time_s']);t=float(report_time)
    if not (t0<=t<=t1):raise ValueError('report time outside recorded accepted step; no extrapolation')
    values=[float(current[k]) for k in ['cut_salt_adv_down_m_s','cut_salt_diff_down_m_s','cut_salt_total_down_m_s','cumulative_cut_salt_m']]
    if not all(math.isfinite(v) for v in values):raise ValueError('nonfinite flux or integral')
    adv,diff,total,c=values;cp=float(previous['cumulative_cut_salt_m'])
    time_bound=(t1-t)*(abs(adv)+abs(diff))
    quadrature_defect=abs((c-cp)-total*(t1-t0))
    floating_bound=8*max(math.ulp(c),math.ulp(cp))
    return dict(bound_m=time_bound+quadrature_defect+floating_bound,
                time_flux_bound_m=time_bound,step_quadrature_defect_m=quadrature_defect,
                float_guard_m=floating_bound,tail_s=t1-t,actual_later_time_s=t1,report_time_s=t)


def pair_tail_bound(previous0,current0,previous1,current1,report_time):
    a=local_tail_bound(previous0,current0,report_time);b=local_tail_bound(previous1,current1,report_time)
    return dict(bound_m=a['bound_m']+b['bound_m'],M0=a,M1=b)


def sensitivity_interval(delta_coarse,delta_fine,timing_bound):
    if timing_bound<0 or not all(math.isfinite(v) for v in [delta_coarse,delta_fine,timing_bound]):raise ValueError('invalid sensitivity input')
    raw=abs(delta_fine-delta_coarse)
    return dict(U_dt_raw_m=raw,U_dt_lower_m=max(0.,raw-timing_bound),U_dt_upper_m=raw+timing_bound,timing_bound_m=timing_bound)
