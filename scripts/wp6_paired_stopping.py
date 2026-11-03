"""Accepted-state pairing with solubility censoring; no solver launches/interpolation.

Time equality allows only 32 binary64 ulps for serialization, not a time bin.
Observed brackets mean (low, high], except an initial crossing [0,0].
The global maximum includes the surface. An interior-first event is a global
crossing while the surface is still below threshold; it is not a surface event.
"""
import math


def same_time(a, b):
    return abs(a-b) <= 32*max(math.ulp(a), math.ulp(b))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def validate_leg(events, stop, returncode=0, threshold=0.1, horizon_tolerance_s=1e-7):
    require(returncode == 0, 'solver exit unsuccessful')
    require(bool(events), 'no accepted events')
    times = [float(r['time_s']) for r in events]
    require(all(math.isfinite(t) for t in times), 'nonfinite accepted time')
    require(times[0] == 0, 'missing initial state')
    require(all(b > a and not same_time(a,b) for a,b in zip(times,times[1:])), 'duplicate/nonmonotone time')
    crossings = {}
    for i, r in enumerate(events):
        require(int(r['step']) == i, 'missing accepted step')
        dt = float(r['dt_s'])
        require(math.isfinite(dt) and (dt == 0 if i == 0 else dt > 0), 'invalid dt')
        if i:
            require(same_time(times[i], times[i-1]+dt), 'accepted dt/time mismatch')
        top, domain = float(r['x_top_max']), float(r['x_domain_max'])
        require(math.isfinite(top) and math.isfinite(domain) and domain >= top, 'invalid maxima')
        for name, value, flag in [('surface',top,'surface_threshold_reached'),('global',domain,'global_threshold_reached')]:
            require(int(r[flag]) in (0,1) and bool(int(r[flag])) == (value >= threshold), 'event flag mismatch')
            if value >= threshold and name not in crossings:
                crossings[name] = dict(status='observed', low_s=times[max(0,i-1)], high_s=times[i], initial=(i==0))
    first = min((v['high_s'] for v in crossings.values()), default=None)
    require(first is None or same_time(first,times[-1]), 'accepted state after first threshold crossing')
    require(stop['initial_diagnostics_written'] and stop['accepted_steps'] == len(events)-1, 'stop step mismatch')
    require(same_time(float(stop['last_accepted_time_s']),times[-1]), 'stop time mismatch')
    require(times[-1] <= float(stop['requested_end_s'])+horizon_tolerance_s, 'past requested horizon')
    reason = ('surface_solubility_threshold' if 'surface' in crossings else
              'interior_solubility_threshold' if 'global' in crossings else 'requested_horizon')
    require(stop['reason'] == reason, 'stop reason mismatch')
    if not crossings:
        # DuMuX accumulated time may stop a few ulps short of requested end.
        # This terminal tolerance NEVER changes matching or event brackets.
        require(abs(times[-1]-float(stop['requested_end_s']))<=horizon_tolerance_s, 'premature uncensored stop')
    for name in ['surface','global']:
        for bound in ['low','high']:
            value = stop[name+'_threshold_'+bound+'_s']
            require((value is None if name not in crossings else value is not None and same_time(float(value),crossings[name][bound+'_s'])), 'stop event bracket mismatch')
    censored = dict(status='right_censored', lower_bound_s=times[-1], cause='applicability_stop' if crossings else 'requested_horizon')
    surface = crossings.get('surface',censored.copy())
    any_event = crossings.get('global',censored.copy())
    # Interior-specific maximum is not recorded separately. Do not invent its
    # crossing when a surface event occurred, even if global crossed too.
    interior = crossings['global'].copy() if reason == 'interior_solubility_threshold' else dict(status='not_separately_identifiable', cause=reason)
    return dict(times=times, stop_reason=reason, first_crossing_high_s=first,
                surface=surface, any_local=any_event, interior_first=interior,
                valid=True, threshold=threshold,horizon_tolerance_s=horizon_tolerance_s,
                accepted_minus_requested_end_s=times[-1]-float(stop['requested_end_s']))


def interval_difference(m0, m1):
    """Bounds for tau_M1 - tau_M0, never substitute censor time for event."""
    a,b = m0['status']=='observed', m1['status']=='observed'
    if a and b:
        return dict(status='interval_observed',lower_s=m1['low_s']-m0['high_s'],upper_s=m1['high_s']-m0['low_s'],bounds='conservative enclosure; crossing brackets are (low,high]')
    if a:
        return dict(status='lower_bound_only',lower_s=m1['lower_bound_s']-m0['high_s'],upper_s=None,bounds='strict lower bound')
    if b:
        return dict(status='upper_bound_only',lower_s=None,upper_s=m1['high_s']-m0['lower_bound_s'],bounds='strict upper bound')
    return dict(status='both_right_censored',lower_s=None,upper_s=None,bounds='difference unidentified')


def paired_times(m0, m1):
    a,b = m0['times'],m1['times']; i=j=0; shared=[]
    cutoff=min((v for v in [m0['first_crossing_high_s'],m1['first_crossing_high_s']] if v is not None),default=None)
    while i<len(a) and j<len(b):
        if same_time(a[i],b[j]):
            if cutoff is None or (max(a[i],b[j]) < cutoff and not same_time(max(a[i],b[j]),cutoff)):
                shared.append((i,j))
            i+=1; j+=1
        elif a[i]<b[j]: i+=1
        else: j+=1
    return dict(indices=shared,shared_times_s=[a[i] for i,j in shared],
                cutoff_crossing_high_s=cutoff, max_time_representation_error_s=max((abs(a[i]-b[j]) for i,j in shared),default=0),
                comparison_status='available' if len(shared)>1 else 'insufficient_common_prethreshold_states',
                surface_time_difference=interval_difference(m0['surface'],m1['surface']),
                any_local_time_difference=interval_difference(m0['any_local'],m1['any_local']),
                different_first_event_types=m0['stop_reason']!=m1['stop_reason'] and all(v['first_crossing_high_s'] is not None for v in [m0,m1]),
                interior_time_difference=dict(status='not_identifiable_without_separate_interior_maxima'))


def budget_rows_at_events(rows, leg, region):
    selected=[r for r in rows if r['region']==region]
    require(len(selected)==len(leg['times']), 'missing budget accepted state: '+region)
    for r,t in zip(selected,leg['times']):
        require(same_time(float(r['time_s']),t), 'budget/event time mismatch')
        require(all(math.isfinite(float(r[k])) for k in ['inventory_salt_m','inventory_initial_salt_m','cumulative_cut_salt_m','cumulative_cut_adv_m','cumulative_cut_diff_m','salt_residual_m']), 'nonfinite budget')
    return selected
