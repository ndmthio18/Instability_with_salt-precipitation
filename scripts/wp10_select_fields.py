"""Extract registered sparse-frame members from an immutable WP10 field archive."""
import json, sys, tarfile
from pathlib import Path

source, target = map(Path, sys.argv[1:3])
indices = [0,48,72,96,120,240,256,360,480,511]
with tarfile.open(source,"r:gz") as src:
    members=sorted((m for m in src.getmembers() if m.name.endswith(".vtu") and "/" not in m.name),key=lambda m:m.name)
    assert len(members)==512 and indices[-1]<len(members)
    chosen=[members[i] for i in indices]
    with tarfile.open(target,"w:gz",compresslevel=1) as dst:
        for m in chosen:
            f=src.extractfile(m); assert f is not None; dst.addfile(m,f)
manifest={"source":str(source),"source_members":len(members),"indices":indices,
          "members":[chosen[i].name for i in range(len(chosen))]}
target.with_suffix(target.suffix+".json").write_text(json.dumps(manifest,indent=2))
print(json.dumps(manifest))
