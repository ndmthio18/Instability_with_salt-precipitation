"""Regression checks for scientific meaning, source preservation and display labels."""
import ast
import csv
import hashlib
import io
import json
import subprocess
from pathlib import Path
from onset_semantics import ROOT, FrozenOnsets, normalize_case, rows
from audit_onset_semantics import OUT, save_csv, save_json

BASELINE = 'b3a960f8f230d6e6e96cc6a13f640ebce597e773'


def git(*args):
    return subprocess.check_output(['git','-C',str(ROOT),*args])


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    checks = {}; catalog = FrozenOnsets()
    provenance = json.loads((OUT/'source_provenance.json').read_text(encoding='utf-8'))
    checks['all_recorded_input_bytes_unchanged'] = all(sha(ROOT/p)==h for p,h in provenance['input_sha256'].items())
    c=catalog.case(36877)
    checks['36877_global_vs_observable'] = (c['lsa_global_critical_mode_nx']==5 and
        c['lsa_observable_mode_nx_theta010']==28 and c['lsa_t_crit_s']==31409.157792677746 and
        c['lsa_t_seed_theta010_s']==45474.35459078279)
    m=json.loads((OUT/'case36877_mapping.json').read_text())
    checks['manifest_retains_restricted_critical_nx'] = m['raw_manifest_lsa']['requested_nx']==m['raw_manifest_lsa']['critical_nx']==28
    checks['sample_minimum_not_refined_detection'] = m['normalized']['dumux_t_detect_sample_s'] != m['normalized']['dumux_t_detect_s']
    all_seeds=[catalog.case(k) for k in catalog.seeds]
    checks['fallback_accepted_without_requiring_theta'] = any(c['lsa_observable_selection_theta010']=='fallback_best' and
        .03<=c['lsa_observable_sharpness_theta010']<.1 and c['lsa_t_seed_theta010_s'] is not None for c in all_seeds)
    checks['no_observable_remains_missing'] = all(c['lsa_t_seed_theta010_s'] is None and c['lsa_observable_mode_nx_theta010'] is None
        for c in all_seeds if c['lsa_observable_status_theta010']=='no_sharp_mode')
    bad=dict(catalog.seeds['36877']);bad['critical_nx']='28'
    try: normalize_case(catalog.observations['36877'],bad)
    except ValueError: checks['inconsistent_global_source_rejected']=True
    else: checks['inconsistent_global_source_rejected']=False
    table_map=json.loads((OUT/'table_mapping.json').read_text())
    for source,spec in table_map.items():
        if not source.endswith('.csv'):continue
        old=rows(ROOT/source);new=rows(OUT/spec['view'])
        checks['raw_columns_preserved/'+source] = len(old)==len(new) and all(all(b[k]==v for k,v in a.items()) for a,b in zip(old,new))
    normalized=rows(OUT/'views/run_audit.csv')+rows(OUT/'views/wp1_supplement_audit.csv')
    checks['520_run_rows'] = len(normalized)==520
    checks['unresolved_not_promoted_to_points'] = all(not r['dumux_t_detect_s'] for r in normalized if r['onset_class']!='resolved_onset')
    checks['transport_time_never_inferred'] = all(not r['transport_t_eff_s'] and r['transport_t_eff_status']=='undefined_not_a_label' for r in normalized)
    assessment=rows(ROOT/'analysis/dense_history_20261102/run_assessment.csv')
    n60=next(r for r in assessment if r['id']=='r06_c36877_n60')
    checks['nx60_not_assigned_nx28_crossing'] = (n60['initialized_mode_nx']=='60' and n60['lsa_observable_mode_nx_theta010']=='28'
        and n60['lsa_t_initialized_mode_crossing_s']!=n60['lsa_t_seed_theta010_s'])
    old=list(csv.DictReader(io.StringIO(git('show',BASELINE+':analysis/dense_history_20261102/run_assessment.csv').decode())))
    checks['dense_original_columns_identical'] = len(old)==len(assessment) and all(all(b[k]==v for k,v in a.items()) for a,b in zip(old,assessment))
    candidates=json.loads((OUT/'views/pilot_candidates.json').read_text(encoding='utf-8'))
    old_candidates=json.loads((ROOT/'analysis/phase0_20261102/pilot_candidates.json').read_text(encoding='utf-8'))
    checks['candidate_raw_objects_identical'] = all({k:v for k,v in b.items() if k!='normalized_onset'}==a for a,b in zip(old_candidates,candidates))
    tracked=git('ls-tree','-r','--name-only',BASELINE).decode().splitlines()
    coverage=[]; protected=[]; identical_csv=[]
    for name in tracked:
        p=ROOT/name
        if name.startswith('analysis/'):
            # Every original numerical/evidence artifact stays byte-identical except the additive dense table and relabelled figures.
            if p.suffix in ('.json','.input','.dat','.csv','.cc','.hh') and name!='analysis/dense_history_20261102/run_assessment.csv':
                if name.endswith('verification.json'):continue
                equal=p.read_bytes()==git('show',BASELINE+':'+name)
                protected.append(equal)
                if p.suffix=='.csv':identical_csv.append(equal)
            if p.suffix in ('.csv','.json','.md','.png','.pdf'):
                if name in table_map: disposition='raw source retained; use '+table_map[name]['view']
                elif name=='analysis/dense_history_20261102/run_assessment.csv':disposition='original columns preserved; explicit standard fields appended'
                elif '/figures/' in name:disposition='time plots relabelled; inventory/field plots unchanged'
                elif p.suffix=='.md' and '/snapshots/' not in name:disposition='current report uses central terminology contract'
                else:disposition='frozen evidence or non-onset diagnostics; contextual mappings in CONTRACT_ZH.md; no event reassignment'
                coverage.append(dict(path=name,disposition=disposition))
    checks['original_evidence_and_registered_inputs_identical']=all(protected)
    checks['other_numeric_tables_byte_identical']=all(identical_csv)
    for p in (ROOT/'scripts').glob('*.py'):
        ast.parse(p.read_text(encoding='utf-8'),filename=str(p))
    checks['all_scripts_parse']=True
    import fitz
    pdf=fitz.open(ROOT/'analysis/dense_history_20261102/figures/dense_history_diagnostics.pdf')
    text='\n'.join(page.get_text() for page in pdf)
    checks['figure_labels_14_timelines']=len(pdf)==17 and text.count('Pink: DuMuX detection plateau')==14
    checks['figure_distinguishes_both_lsa_times']='LSA t_crit (global)' in text and 'LSA t_seed(0.10)' in text and 'LSA crossing of initialized mode' in text
    save_csv(OUT/'artifact_coverage.csv',coverage)
    result=dict(passed=all(checks.values()),checks=checks,baseline_commit=BASELINE,
        protected_artifacts=len(protected),artifact_coverage_rows=len(coverage),
        lsa_recomputed=False,simulations_started=False,
        scope='mapping and regression checks; no renewed physical or eigenvalue validation')
    save_json(OUT/'verification.json',result)
    outputs=[p for p in OUT.rglob('*') if p.is_file() and p.name!='output_manifest.json']
    outputs += [ROOT/'analysis/dense_history_20261102/run_assessment.csv',ROOT/'analysis/dense_history_20261102/REPORT_ZH.md']
    outputs += list((ROOT/'analysis/dense_history_20261102/figures').glob('*'))
    save_json(OUT/'output_manifest.json',dict(output_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in outputs},
        scripts_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in (ROOT/'scripts').glob('*.py')}))
    print(json.dumps(result,indent=2))
    raise SystemExit(0 if result['passed'] else 1)


if __name__=='__main__':main()
