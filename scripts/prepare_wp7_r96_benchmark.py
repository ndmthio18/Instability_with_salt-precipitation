"""Prepare the conditional R96 ten-step cost input without launching it."""
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'analysis/wp7_shortwave_validation_20261103';SOURCE=OUT/'runs/WP7-R48-M1-bench10';DEST=OUT/'runs/WP7-R96-M1-bench10';DT=113.68588647695695
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    assert not DEST.exists();DEST.mkdir()
    text=(SOURCE/'params.input').read_text().replace('Name = WP7-R48-M1-bench10','Name = WP7-R96-M1-bench10').replace('Cells0 = 48','Cells0 = 96')
    (DEST/'params.input').write_text(text);(DEST/'eigenprofile.dat').write_bytes((SOURCE/'eigenprofile.dat').read_bytes())
    src=json.loads((SOURCE/'plan.json').read_text());plan=dict(id='WP7-R96-M1-bench10',role='M1',stage='r96_cost_only',case_id=36877,nx_cells=96,nz_cells=200,T_s=10*DT,dt_s=DT,scalar_unknowns=2*97*201,input_sha256=sha(DEST/'params.input'),profile_sha256=sha(DEST/'eigenprofile.dat'),elf_sha256=src['elf_sha256'],preregistration_sha256=src['preregistration_sha256'],trigger_sha256=sha(OUT/'r96_cost_trigger_E06.json'),hard_wall_limit_s=600,minimum_free_disk_bytes=20*1024**3)
    (DEST/'plan.json').write_text(json.dumps(plan,indent=2)+'\n');print(json.dumps(plan,indent=2))
if __name__=='__main__':main()
