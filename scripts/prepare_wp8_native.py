"""Create an isolated WP8 source from frozen WP7; build staging is separate."""
import hashlib,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OLD=ROOT/'model/transport_wp7_localflux_v1/source';NEW=ROOT/'model/transport_wp8_background_v1/source';OUT=ROOT/'analysis/wp8_zero_mode_20261103'
def main():
    files={p.relative_to(OLD).as_posix():p.read_bytes() for p in OLD.rglob('*') if p.is_file()}
    mainpath='appl/simpleRichardsnc/main.cc';s=files[mainpath].decode()
    changes=[('#include "transportdiagnostics.hh"','#include "transportdiagnostics.hh"\n#include "wp8_equal_height.hh"'),('nonLinearSolver.solve(x, *timeLoop);','paper2ConstrainedStep(*assembler, x, *gridGeometry);'),('timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));','timeLoop->setTimeStepSize(timeLoop->timeStepSize()); // WP8 fixed dt; failure aborts')]
    for a,b in changes:assert s.count(a)==1;s=s.replace(a,b)
    files[mainpath]=s.encode();files['appl/simpleRichardsnc/wp8_equal_height.hh']=(ROOT/'scripts/wp8_equal_height.hh').read_bytes()
    for name,data in files.items():
        p=NEW/name;p.parent.mkdir(parents=True,exist_ok=True)
        if p.exists() and p.read_bytes()!=data:
            prior=json.loads((NEW.parent/'source_manifest.json').read_text());assert hashlib.sha256(p.read_bytes()).hexdigest()==prior[name]['sha256']
            history=OUT/'build/source_history'/prior[name]['sha256'];history.parent.mkdir(parents=True,exist_ok=True);history.write_bytes(p.read_bytes())
        p.write_bytes(data)
    manifest={n:dict(bytes=len(b),sha256=hashlib.sha256(b).hexdigest()) for n,b in files.items()}
    (NEW.parent/'source_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    original=(ROOT/'scripts/build_wp7_remote.py').read_text();generated=original.replace('transport_wp7_localflux_v1','transport_wp8_background_v1').replace('paper2_wp7_localflux_build_20261103','paper2_wp8_background_build_20261103').replace('analysis/wp7_shortwave_validation_20261103/build','analysis/wp8_zero_mode_20261103/build').replace("'--parallel','2'","'--parallel','1'").replace('os.getloadavg()[0]<8 and mem>8*1024**3','mem>4*1024**3')
    revision=sys.argv[1] if len(sys.argv)>1 else 'r03';assert revision in ['r02','r03']
    generated=generated.replace('paper2_wp8_background_build_20261103','paper2_wp8_background_build_20261103_'+revision).replace("out/(Path(j['attempt']).name+'_launch.json')","out/('"+revision+"_'+Path(j['attempt']).name+'_launch.json')")
    p=ROOT/'tmp/wp8_zero_mode_20261103/build_native.py';p.write_text(generated.replace("ROOT=Path(__file__).resolve().parents[1]",'ROOT=Path('+repr(str(ROOT))+')'))
    print('isolated source staged locally; generated build driver',p)
if __name__=='__main__':main()
