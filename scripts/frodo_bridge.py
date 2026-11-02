"""Run this inside WSL; isolate an existing SSH identity for one operation."""
import argparse, os, subprocess, tempfile
from pathlib import Path

ap=argparse.ArgumentParser()
ap.add_argument('operation',choices=['run','push','pull'])
ap.add_argument('items',nargs='+')
args=ap.parse_args()
with tempfile.TemporaryDirectory(prefix='lsa_audit_ssh_') as d:
    key=Path(d)/'key'
    fd=os.open(key,os.O_WRONLY|os.O_CREAT|os.O_EXCL,0o600)
    with os.fdopen(fd,'wb') as f:
        f.write(Path('/mnt/c/Users/12608/.ssh/zhixin_iws').read_bytes())
    options=['-o','HostKeyAlias=129.69.98.232',
        '-o','UserKnownHostsFile=/mnt/c/Users/12608/.ssh/known_hosts',
        '-o','StrictHostKeyChecking=yes','-o','BatchMode=yes',
        '-o','ConnectTimeout=8','-i',str(key)]
    host='zhixin@frodo.iws.uni-stuttgart.de'
    if args.operation=='run':
        assert len(args.items)==1
        cmd=['ssh',*options,host,args.items[0]]
    elif args.operation=='push':
        cmd=['scp',*options,*args.items[:-1],host+':'+args.items[-1]]
    else:
        assert len(args.items)==2
        cmd=['scp',*options,host+':'+args.items[0],args.items[1]]
    result=subprocess.run(cmd)
raise SystemExit(result.returncode)
