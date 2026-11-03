"""Final input/artifact hashes and cross-table consistency; no solver invocation."""
import datetime,py_compile
import fitz
from analyze_post_onset import *
def main():
 checks={};rows=readcsv(OUT/'run_metrics.csv')
 compared=[r for r in rows if r['source_group']=='main' and r['dumux_detection_class']=='resolved_onset' and r['series_status']=='hash_verified']
 checks['234_recomputed_G_match_frozen']=len(compared)==234 and all(np.isclose(num(r['G_end']),num(r['frozen_G_end']),rtol=1e-12,atol=0) for r in compared)
 checks['404_series_sha_match']=all(hashlib.sha256((ROOT/r['local']).read_bytes()).hexdigest()==r['sha256'] for r in json.load(open(OUT/'series_source_manifest.json')))
 checks['t_eff_undefined']=all(r.get('transport_t_eff_s') in ('',None) for r in rows)
 checks['no_fake_resolved_detection']=all(r['dumux_t_detect_s']=='' for r in rows if r['dumux_detection_class']!='resolved_onset')
 checks['initial_profile_checks_pass']=all(r['initial_profile_passed']=='True' for r in readcsv(OUT/'sparse_field_metrics.csv'))
 checks['tests_passed']=json.load(open(OUT/'test_results.json'))['passed']
 pdf=[]
 for p in sorted(OUT.glob('figure*.pdf')):
  d=fitz.open(p);pdf.append(dict(path=p.name,pages=len(d),text_characters=sum(len(q.get_text()) for q in d),bytes=p.stat().st_size));d.close()
 checks['seven_single_page_pdfs']=len(pdf)==7 and all(p['pages']==1 and p['text_characters']>100 for p in pdf)
 checks['seven_pngs']=len(list(OUT.glob('figure*.png')))==7
 scripts=sorted((ROOT/'scripts').glob('*post_onset*.py'))
 for p in scripts:py_compile.compile(str(p),doraise=True)
 checks['scripts_compile']=True
 paths=[ROOT/'command/4.md',ROOT/'paper1_sources.lock.json',ROOT/'analysis/onset_semantics_20261102/CONTRACT_ZH.md',ROOT/'analysis/onset_semantics_20261102/views/run_audit.csv',ROOT/'analysis/onset_semantics_20261102/views/wp1_supplement_audit.csv',ROOT/'analysis/phase0_20261102/remote_frodo_inventory.json',ROOT/'analysis/dense_history_20261102/plan.json',ROOT/'analysis/dense_history_20261102/time_reconciliation.csv',ROOT/'analysis/dense_history_20261102/run_assessment.csv',ROOT/'analysis/dense_history_20261102/inventory_summary.csv',ROOT/'analysis/dense_history_20261102/evidence/product_manifest.json',ROOT/'analysis/pilot0_integrated_20261102/provenance/extract_onset_2025.py',*scripts]
 paths+=list((ROOT/'data/external/paper1/158395f3d25a/handoff').glob('*'))
 sources=[dict(path=str(p.relative_to(ROOT)).replace('\\','/'),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size) for p in paths if p.is_file()]
 savejson(OUT/'source_manifest.json',sources)
 savejson(OUT/'verification.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),checks=checks,pdf_checks=pdf,visual_QA='Figures 1-6 contact sheet and Figure7 individually inspected; explicit shallow crop and absolute shared anomaly scale; labels legible',simulations_started=0,lsa_solves=0,remote_original_writes=0,all_passed=all(checks.values())))
 manifest=[dict(path=p.name,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size) for p in sorted(OUT.iterdir()) if p.is_file() and p.name!='artifact_manifest.json'];savejson(OUT/'artifact_manifest.json',manifest)
 assert all(checks.values()),checks
 print(json.dumps(checks));print('artifacts',len(manifest),'total MB',sum(r['bytes'] for r in manifest)/1e6)
if __name__=='__main__':main()
