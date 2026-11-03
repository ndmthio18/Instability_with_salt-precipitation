"""Create isolated v2 source and build driver; preserve v1 and all old evidence."""
import json,hashlib,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'analysis/wp8_newton_correction_20261103';NEW=ROOT/'model/transport_wp8_background_v2/source'
def main(variant='candidate'):
    assert variant in ['candidate','probe','probe_fine']
    suffix='' if variant=='candidate' else '_'+variant
    new=ROOT/('model/transport_wp8_background_v2'+suffix+'/source')
    OUT.mkdir(exist_ok=True);cache=ROOT/'tmp/wp8_newton_correction_20261103';cache.mkdir(exist_ok=True)
    old=ROOT/'model/transport_wp8_background_v1/source';files={p.relative_to(old).as_posix():p.read_bytes() for p in old.rglob('*') if p.is_file()}
    s=files['appl/simpleRichardsnc/main.cc'].decode();s=s.replace('problem->applyInitialSolution(x);','problem->applyInitialSolution(x);\n    wp8DiagnosticReplay(x, *gridGeometry);').replace('paper2ConstrainedStep(*assembler, x, *gridGeometry);','paper2ConstrainedStep(*assembler, x, *gridGeometry, timeLoop->timeStepSize());');files['appl/simpleRichardsnc/main.cc']=s.encode();files['appl/simpleRichardsnc/wp8_equal_height.hh']=(ROOT/'scripts/wp8_convergence_v2.hh').read_bytes()
    key='appl/simpleRichardsnc/wp8_equal_height.hh'
    if variant!='candidate':files[key]=files[key].replace(b'&&iter==3)',b'&&iter==2)')
    if variant=='probe_fine':files[key]=files[key].replace(b'for(double rel:{1e-7,1e-6,1e-5})',b'for(double rel:{1e-10,1e-9,1e-8})')
    for n,b in files.items():
        p=new/n;p.parent.mkdir(parents=True,exist_ok=True)
        if p.exists():assert p.read_bytes()==b,n
        else:p.write_bytes(b)
    (new.parent/'source_manifest.json').write_text(json.dumps({n:dict(sha256=hashlib.sha256(b).hexdigest(),bytes=len(b)) for n,b in files.items()},indent=2)+'\n')
    s=(ROOT/'scripts/build_wp7_remote.py').read_text().replace('transport_wp7_localflux_v1','transport_wp8_background_v2').replace('paper2_wp7_localflux_build_20261103','paper2_wp8_background_v2_build_20261103').replace('analysis/wp7_shortwave_validation_20261103/build','analysis/wp8_newton_correction_20261103/build').replace("'--parallel','2'","'--parallel','1'").replace('os.getloadavg()[0]<8 and mem>8*1024**3','mem>4*1024**3').replace('ROOT=Path(__file__).resolve().parents[1]','ROOT=Path('+repr(str(ROOT))+')')
    if suffix:s=s.replace('transport_wp8_background_v2','transport_wp8_background_v2'+suffix).replace('paper2_wp8_background_v2_build','paper2_wp8_background_v2'+suffix+'_build').replace('analysis/wp8_newton_correction_20261103/build','analysis/wp8_newton_correction_20261103/'+variant+'_build')
    (cache/(('' if variant=='candidate' else variant+'_')+'build.py')).write_text(s);print('isolated v2 prepared; build script ready',variant)
if __name__=='__main__':main(sys.argv[1] if len(sys.argv)>1 else 'candidate')
