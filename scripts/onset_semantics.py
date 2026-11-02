"""Read-only, source-checked onset naming adapter. No solver or network imports.

Keep legacy fields unchanged. Use the explicit fields below for new analysis.
Missing / censored / unresolved values stay missing; no times are inferred.
"""
import csv
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data/external/paper1/158395f3d25a/handoff'
VERSION = 'paper2-onset-v1'


def rows(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


def number(value):
    if value is None or value == '':
        return None
    value = float(value)
    if not math.isfinite(value):
        raise ValueError('Non-finite frozen time or mode')
    return value


def integer(value):
    value = number(value)
    if value is not None and value != int(value):
        raise ValueError('Non-integer mode')
    return int(value) if value is not None else None


def normalize_case(observation, seed=None):
    seed = seed or {}
    accepted = seed.get('status') == 'ok'
    if seed and (number(seed['sharpness_min']) != .10 or number(seed['sharpness_floor']) != .03):
        raise ValueError('Expected theta=0.10, floor=0.03')
    if (accepted or seed.get('critical_nx') not in (None, '')) and integer(seed.get('critical_nx')) != integer(observation.get('critical_nx')):
        raise ValueError('Frozen global modes disagree')
    if accepted:
        selection = seed['selection']
        required = .10 if selection == 'preferred' else .03
        if selection not in ('preferred', 'fallback_best') or number(seed['sharpness']) < required:
            raise ValueError('Invalid frozen observable selection')
        if integer(seed['nx_seeded']) is None or number(seed['t_onset_seeded_s']) is None:
            raise ValueError('Accepted observable mode is missing its crossing')
    return dict(
        onset_semantics_version=VERSION,
        lsa_global_critical_mode_nx=integer(observation.get('critical_nx')),
        lsa_t_crit_s=number(observation.get('t_onset_s')),
        lsa_t_hat_crit=number(observation.get('t_hat_onset')),
        lsa_critical_status=observation.get('onset_status'),
        lsa_first_event_type=observation.get('first_event_type'),
        lsa_critical_censoring_type=observation.get('censoring_type'),
        lsa_observable_mode_nx_theta010=integer(seed.get('nx_seeded')) if accepted else None,
        lsa_t_seed_theta010_s=number(seed.get('t_onset_seeded_s')) if accepted else None,
        lsa_t_hat_seed_theta010=number(seed.get('t_hat_onset_seeded')) if accepted else None,
        lsa_observable_status_theta010=seed.get('status', 'not_available'),
        lsa_observable_selection_theta010=seed.get('selection') or None,
        lsa_observable_sharpness_theta010=number(seed.get('sharpness')),
        lsa_observable_theta=.10, lsa_observable_floor=.03,
        lsa_observable_fallback_enabled=True,
        lsa_observable_protocol_hash=seed.get('seeding_protocol_hash') or None,
        transport_t_eff_s=None, transport_t_eff_status='undefined_not_a_label')


class FrozenOnsets:
    def __init__(self):
        self.lock = json.loads((ROOT/'paper1_sources.lock.json').read_text())
        for name in ('observations.csv', 'seeded_modes.csv', 'dumux_onsets.csv', 'cases.csv'):
            expected = self.lock['release']['handoff_files'][name]
            if hashlib.sha256((DATA/name).read_bytes()).hexdigest() != expected:
                raise ValueError('Frozen input changed: '+name)
        self.observations = {r['case_id']: r for r in rows(DATA/'observations.csv')}
        seeds = [r for r in rows(DATA/'seeded_modes.csv') if number(r['sharpness_min']) == .10]
        self.seeds = {r['case_id']: r for r in seeds}
        if len(self.seeds) != len(seeds):
            raise ValueError('Duplicate theta=0.10 case')
        self.runs = {r['run_key']: r for r in rows(DATA/'dumux_onsets.csv')}
        self.cases = {r['case_id']: r for r in rows(DATA/'cases.csv')}

    def case(self, case_id):
        key = str(case_id)
        return normalize_case(self.observations[key], self.seeds.get(key))

    def run(self, raw):
        """Raw handoff or WP1 audit row; preserve every supplied legacy field."""
        result = dict(raw)
        standard = self.case(raw['case_id'])
        # A same-named existing field must agree, never silently overwrite it.
        for k, v in standard.items():
            if k in raw and raw[k] not in ('', None) and str(raw[k]) != str(v):
                if number(raw[k]) != number(v):
                    raise ValueError('Conflicting normalized field: '+k)
        result.update(standard)
        nx = integer(raw.get('seeded_nx'))
        t = number(raw.get('t_lsa_mode_s'))
        detected = number(raw.get('onset_s'))
        resolved = raw.get('onset_class') == 'resolved_onset'
        sample = number(raw.get('onset_sample_s'))
        if raw.get('local_manifest'):
            path = Path(raw['local_manifest'])
            if hashlib.sha256(path.read_bytes()).hexdigest() != raw['local_manifest_sha256']:
                raise ValueError('WP1 source manifest changed: '+str(path))
            manifest = json.loads(path.read_text(encoding='utf-8'))
            lsa, onset = manifest['lsa'], manifest['dumux']['onset']
            if (integer(lsa['requested_nx']) != nx or number(lsa['t_onset_lsa_s']) != t or
                    onset['onset_class'] != raw['onset_class']):
                raise ValueError('WP1 audit/manifest mapping mismatch')
            expected = onset.get('onset_refined_s')
            if expected is None:
                expected = onset.get('onset_s')
            if number(expected) != detected:
                raise ValueError('WP1 detection source mismatch')
            sample = number(onset.get('onset_s'))
        result.update(
            initialized_mode_nx=nx,
            lsa_t_initialized_mode_crossing_s=t,
            lsa_initialized_mode_source=raw.get('t_lsa_mode_source') or 'WP1 manifest requested_nx/t_onset_lsa_s',
            initializes_observable_theta010=(nx == standard['lsa_observable_mode_nx_theta010']) if nx is not None else None,
            dumux_t_detect_s=detected if resolved else None,
            dumux_t_detect_estimate_s=detected,
            dumux_t_detect_sample_s=sample,
            dumux_t_detect_low_s=number(raw.get('onset_low_s')),
            dumux_t_detect_high_s=number(raw.get('onset_high_s')),
            dumux_detection_class=raw.get('onset_class'),
            dumux_detection_point_valid=resolved and detected is not None,
            onset_case_source='frozen handoff observations.csv + seeded_modes.csv(theta=0.10)',
            onset_run_source=raw.get('local_manifest') or 'frozen handoff dumux_onsets.csv')
        return result

    def run_by_key(self, run_key):
        return self.run(self.runs[run_key])


def standard_fields(row):
    """Only new public fields, without legacy aliases or raw source columns."""
    return {k: v for k, v in row.items() if k.startswith((
        'lsa_global_', 'lsa_t_crit_', 'lsa_t_hat_', 'lsa_critical_', 'lsa_first_event_',
        'lsa_observable_', 'lsa_t_seed_theta', 'initialized_', 'initializes_',
        'lsa_t_initialized_', 'lsa_initialized_', 'dumux_', 'transport_t_eff', 'onset_'))
        and k not in ('onset_s', 'onset_class', 'onset_low_s', 'onset_high_s', 'onset_sample_s')}
