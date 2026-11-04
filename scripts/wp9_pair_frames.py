"""Build a frame-aligned view of an OutputInterval=2 run, for the reviewed field passes.

`wp8_fields.py` and `wp7_spatial.py` both assume one archived frame per recorded accepted
step. The R96 leg advances with dt/2 and archives every second step, so the assumption
holds only on the even steps. Rather than edit either reviewed script, this builds a
sibling directory that presents exactly those even steps:

  * hard links to every .vtu and the .pvd, so no field byte is copied or rewritten;
    hard links rather than symlinks, because wp7_spatial.py resolves each frame path and
    requires it to stay inside the directory it was given;
  * params.input and eigenprofile.dat copied unchanged;
  * surface_events.csv restricted to the even steps, in order;
  * post/field_hashes.json with the same sha256 values and the paths retargeted.

The sha256 values are carried over untouched, so both passes still verify every frame
against the hash recorded when it was written.

Usage: python3 wp9_pair_frames.py <run_dir> <view_dir> <stride>
"""
import csv,hashlib,json,os,sys
from pathlib import Path


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main(run,view,stride='2'):
    run=Path(run).resolve();view=Path(view).resolve();stride=int(stride)
    assert run!=view and not view.exists(),'the view directory must be new'
    (view/'post').mkdir(parents=True)
    records=json.loads((run/'post/field_hashes.json').read_text())
    with (run/'surface_events.csv').open() as f:
        reader=csv.reader(f);header=next(reader);rows=list(reader)
    step=header.index('step')
    kept=[r for r in rows if int(r[step])%stride==0]
    assert len(kept)==len(records),f'{len(kept)} even steps but {len(records)} frames'
    assert [int(r[step]) for r in kept]==list(range(0,stride*len(records),stride))
    with (view/'surface_events.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(header);w.writerows(kept)
    for n in ['params.input','eigenprofile.dat']:
        (view/n).write_bytes((run/n).read_bytes())
    out=[]
    for rec in records:
        src=run/Path(rec['path']).name
        dst=view/src.name
        os.link(src,dst)
        assert dst.resolve()==dst,'the view must not resolve outside itself'
        out.append(dict(rec,path=str(dst)))
    for q in run.glob('*.pvd'):os.link(q,view/q.name)
    (view/'post/field_hashes.json').write_text(json.dumps(out,indent=2))
    receipt=dict(run=str(run),view=str(view),stride=stride,frames=len(out),
        recorded_steps=len(rows),kept_steps=len(kept),
        first_kept_step=int(kept[0][step]),last_kept_step=int(kept[-1][step]),
        hard_linked_fields=len(out),
        link_type='hard',copied_bytes=sum((view/n).stat().st_size for n in
            ['params.input','eigenprofile.dat','surface_events.csv']),
        sha256_values_unchanged=all(a['sha256']==b['sha256'] for a,b in zip(records,out)),
        source_manifest_sha256=sha(run/'post/field_hashes.json'),
        view_manifest_sha256=sha(view/'post/field_hashes.json'),
        script_sha256=sha(Path(__file__)))
    (view/'pairing_receipt.json').write_text(json.dumps(receipt,indent=2))
    print(json.dumps(receipt))


if __name__=='__main__':main(*sys.argv[1:])
