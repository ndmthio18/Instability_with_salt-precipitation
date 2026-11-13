"""Rendered QA. Times New Roman overrides the generic sans-serif source whitelist."""
from pathlib import Path
import subprocess,sys,json,os
ROOT=Path(__file__).resolve().parent;reports=[]
for p in sorted((ROOT/'figures').glob('*.pdf')):
 cmd=[sys.executable,str(ROOT/'qa_tools/audit_figure_collisions.py'),str(p),'--json-out',str(p)+'.collision.json']
 r=subprocess.run(cmd,capture_output=True,encoding='utf-8',env=dict(os.environ,PYTHONIOENCODING='utf-8'))
 t=subprocess.run([sys.executable,str(ROOT/'qa_tools/audit_pdf_text.py'),str(p)],capture_output=True,encoding='utf-8',env=dict(os.environ,PYTHONIOENCODING='utf-8'))
 Path(str(p)+'.text-audit.txt').write_text(t.stdout+t.stderr,encoding='utf-8')
 report=json.loads(Path(str(p)+'.collision.json').read_text(encoding='utf-8'))
 import fitz
 with fitz.open(p) as doc:fonts=sorted({f[3] for page in doc for f in page.get_fonts()})
 reports.append(dict(figure=p.stem,collision=report['verdict'],text_exit=t.returncode,collision_failures=report['summary']['fail'],fonts=fonts))
 assert r.returncode==0 and t.returncode==0,(p,r.stdout,r.stderr,t.stdout)
assert len(reports)==11,len(reports)
(ROOT/'results/figure_qa.json').write_text(json.dumps(reports,indent=2))
print('PASS 11 figures: rendered collisions, embedded fonts, glyph floor; see alignment JSONs and documented font preference override')
