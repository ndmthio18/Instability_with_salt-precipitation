"""On compute node: archive exact selected rows; retain full remote outputs."""
import csv
import io
import tarfile
import sys
from pathlib import Path
R=Path.cwd()
prefix=sys.argv[1] if len(sys.argv)>1 else ''
stem=sys.argv[2] if len(sys.argv)>2 else 'compact'+('_'+prefix if prefix else '')
with tarfile.open(stem+'.tar.gz','w:gz') as tar:
    for path in sorted((R/'raw').glob(prefix+'*')):
        if not path.name.endswith('_profiles.csv'):
            tar.add(path,arcname='raw/'+path.name);continue
        with path.open() as f: rows=list(csv.reader(f))
        times=sorted({float(row[0]) for row in rows[1:]})
        selected=set(times[:12]+times[-1:])
        if times:
            for fraction in [.02,.05,.1,.2,.4,.6,.8,1]:
                selected.add(min(times,key=lambda t:abs(t-fraction*times[-1])))
        out=io.StringIO();writer=csv.writer(out,lineterminator='\n');writer.writerow(rows[0])
        writer.writerows(row for row in rows[1:] if float(row[0]) in selected)
        data=out.getvalue().encode();info=tarfile.TarInfo('raw/'+path.name);info.size=len(data);tar.addfile(info,io.BytesIO(data))
    for path in (R/'logs').glob('*'):tar.add(path,arcname='logs/'+path.name)
    for name in ['build_sha256.txt','build_p0_sha256.txt','build_p1_sha256.txt','accounting.txt']:
        if (R/name).exists():tar.add(R/name,arcname=name)
