"""Build additive onset views and provenance from frozen records; no solves."""
import argparse
import csv
import hashlib
import json
import subprocess
from collections import Counter
from pathlib import Path
from onset_semantics import ROOT, DATA, VERSION, FrozenOnsets, rows, number, standard_fields

OUT = ROOT/'analysis/onset_semantics_20261102'
PARENT = ROOT.parents[1]
PDF = Path(r'D:\3OutputFile\1Research\Article\22-1-surrogate_model\69e50893f5d9b55adfb56494\AWR_onstime\onstime.pdf')
PDF_SHA = 'e57fc093419d388ca9bdf005ac5a5b08bb73ae5a402ead7b08a33c1a39bfb498'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save_json(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')


def save_csv(path, records):
    fields = list(dict.fromkeys(k for r in records for k in r))
    with path.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader(); w.writerows(records)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pdf', type=Path, default=PDF)
    args = ap.parse_args()
    catalog = FrozenOnsets()
    OUT.mkdir(exist_ok=True)
    views = OUT/'views'; views.mkdir(exist_ok=True)
    inputs = {}
    for name in ('observations.csv', 'seeded_modes.csv', 'dumux_onsets.csv', 'cases.csv'):
        inputs[str((DATA/name).relative_to(ROOT))] = sha(DATA/name)
    mapping = {}
    tables = [
        ('analysis/phase0_20261102/run_audit.csv', 'run_audit.csv'),
        ('analysis/phase0_20261102/wp1_supplement_audit.csv', 'wp1_supplement_audit.csv'),
        ('analysis/phase0_20261102/pilot_screening.csv', 'pilot_screening.csv'),
        ('analysis/highk_unperturbed_42929_20261102/screening.csv', 'highk_screening.csv')]
    case_ids = set()
    for source, target in tables:
        records = rows(ROOT/source); normalized = []
        inputs[source] = sha(ROOT/source)
        for r in records:
            full = catalog.runs[r['run_key']] if 't_lsa_mode_s' not in r else r
            new = dict(r)
            new.update(standard_fields(catalog.run(full)))
            # Preserve original columns and exact CSV values, including legacy aliases.
            for k, v in r.items():
                if k in new and str(new[k]) != str(v):
                    assert number(new[k]) == number(v), (source, k)
                    new[k] = v
            normalized.append(new); case_ids.add(str(r['case_id']))
        save_csv(views/target, normalized)
        mapping[source] = dict(view='views/'+target, rows=len(records), raw_fields_preserved=True)

    candidates_path = ROOT/'analysis/phase0_20261102/pilot_candidates.json'
    candidates = json.loads(candidates_path.read_text())
    inputs[str(candidates_path.relative_to(ROOT))] = sha(candidates_path)
    for c in candidates:
        c['normalized_onset'] = catalog.case(c['case_id'])
        case_ids.add(str(c['case_id']))
    save_json(views/'pilot_candidates.json', candidates)
    mapping[str(candidates_path.relative_to(ROOT))] = dict(view='views/pilot_candidates.json', rows=2, raw_fields_preserved=True)
    save_csv(views/'case_onsets.csv', [dict(case_id=cid, **catalog.case(cid)) for cid in sorted(case_ids, key=int)])

    # All frozen theta=0.10 selections: audit acceptance and strict PDF time bound.
    accepted, equal_stop, beyond_stop = [], [], []
    for cid, seed in catalog.seeds.items():
        c = catalog.case(cid)
        if c['lsa_observable_status_theta010'] != 'ok':
            continue
        accepted.append(seed)
        obs = catalog.observations[cid]
        bounds = [number(catalog.cases[cid]['application_horizon_s']), number(obs['t_precipitation_s'])]
        stop = min(x for x in bounds if x is not None)
        t = c['lsa_t_seed_theta010_s']
        if t == stop: equal_stop.append(cid)
        if t > stop: beyond_stop.append(cid)
    assert not equal_stop and not beyond_stop, 'Frozen acceptance differs from strict PDF time bound'

    legacy_path = ROOT/'analysis/matched_pair_36877_20261102/legacy_nx28/case_manifest.json'
    legacy = json.loads(legacy_path.read_text())
    inputs[str(legacy_path.relative_to(ROOT))] = sha(legacy_path)
    lsa = legacy['lsa']; c = catalog.case(36877)
    assert lsa['requested_nx'] == lsa['critical_nx'] == c['lsa_observable_mode_nx_theta010'] == 28
    assert c['lsa_global_critical_mode_nx'] == 5
    assert lsa['t_onset_lsa_s'] == c['lsa_t_seed_theta010_s']
    native_key = 'frodo:c1b/runs_600/shard-001/0017--case-36877--native'
    native = catalog.run_by_key(native_key)
    detect_difference = native['dumux_t_detect_s'] - legacy['dumux']['onset']['onset_refined_s']
    assert abs(detect_difference) < 1e-6  # retain both source values, do not round / replace
    save_json(OUT/'case36877_mapping.json', dict(
        normalized=standard_fields(native),
        raw_manifest_lsa={k:lsa[k] for k in ('requested_nx', 'nx_max', 'critical_nx', 'least_stable_critical_nx',
            't_hat_onset_lsa', 't_onset_lsa_s', 'eigenprofile_t_hat', 'reference_time_s')},
        raw_manifest_detection=legacy['dumux']['onset'],
        handoff_minus_manifest_refined_s=detect_difference,
        source_manifest_sha256=sha(legacy_path),
        conclusion='critical_nx=28 is the result within requested_nx=28, not the global minimum over nx=1..30',
        limitation='Exact historical launcher argv and source revision are not recorded in this manifest; frozen code establishes the interface semantics, not the historical executable identity.'))

    # Git bytes, not potentially edited upstream working files; source snippets carry line numbers.
    commit = catalog.lock['paper1_repository']['commit']
    excerpts = {
        'code/lsa_unsaturated.py': [(571,605)],
        'code/parameter_database/solvers.py': [(286,335),(402,434),(730,823),(890,928),(950,975)],
        'code/parameter_database/pipeline.py': [(245,275),(280,345),(424,451)],
        'code/parameter_database/run_design_campaign_parallel.py': [(163,178)],
        'code/parameter_database/dominant_mode_cli.py': [(90,145)],
        'code/parameter_database/seeded_mode_cli.py': [(70,145)]}
    identities = []; text = []
    for path, ranges in excerpts.items():
        raw = subprocess.check_output(['git','-C',str(PARENT),'show',commit+':'+path])
        oid = subprocess.check_output(['git','-C',str(PARENT),'rev-parse',commit+':'+path],text=True).strip()
        identities.append(dict(path=path, commit=commit, git_blob_oid=oid, sha256=hashlib.sha256(raw).hexdigest(), lines=ranges))
        lines = raw.decode().splitlines()
        text += ['\n### '+path+' @ '+commit, '```text']
        for start,end in ranges:
            text += [f'{i}: {lines[i-1]}'.rstrip() for i in range(start,min(end,len(lines))+1)]
            text.append('...')
        text.append('```')
    (OUT/'SOURCE_EXCERPTS.md').write_text('\n'.join(text)+'\n', encoding='utf-8')
    if sha(args.pdf) != PDF_SHA:
        raise ValueError('User-selected PDF has changed; review before accepting new definitions')
    import fitz
    doc = fitz.open(args.pdf)
    extracted = '\n\n'.join(f'=== PDF page {i+1} (printed page {i}) ===\n'+doc[i].get_text() for i in range(7,14))
    # Escape PDF control glyphs so the extract remains reviewable text in Git.
    extracted = ''.join(f'\\x{ord(c):02x}' if ord(c)<32 and c not in '\n\t' else c for c in extracted)
    extracted = '\n'.join(line.rstrip() for line in extracted.splitlines())+'\n'
    (OUT/'paper1_sections_2_2_to_2_4.txt').write_text(extracted, encoding='utf-8')
    save_json(OUT/'source_provenance.json', dict(
        pdf=dict(path=str(args.pdf), sha256=PDF_SHA, bytes=args.pdf.stat().st_size, total_pages=len(doc),
                 relevant_pdf_pages=[8,9,10,11,12,13,14], selected_by_user=True),
        code=identities, input_sha256=inputs, paper1_read_only=True))
    save_json(OUT/'table_mapping.json', mapping)
    save_json(OUT/'selection_audit.json', dict(
        schema=VERSION, frozen_theta010_rows=len(catalog.seeds), accepted=len(accepted),
        accepted_selections=dict(Counter(r['selection'] for r in accepted)),
        statuses=dict(Counter(r['status'] for r in catalog.seeds.values())),
        equals_time_limit=equal_stop, beyond_time_limit=beyond_stop,
        boundary_note='PDF Eq.10 uses strict <; frozen selector excludes only >. No accepted frozen theta010 row is on or beyond the boundary.',
        verification_scope='Stored scalar acceptance and source trace only; candidate optimization and eigenvalues were not recomputed.',
        simulations_started=False, lsa_recomputed=False))
    print(json.dumps(dict(views=mapping, accepted=len(accepted), selections=dict(Counter(r['selection'] for r in accepted))), indent=2))


if __name__ == '__main__':
    main()
