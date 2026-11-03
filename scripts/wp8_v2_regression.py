"""Command 14: bounded diagnostic runs only; no full-horizon launch interface."""
import json, sys
import numpy as np
import wp8_strict_campaign as c
c.OUT=c.ROOT/'analysis/wp8_newton_correction_20261103'
c.CACHE=c.ROOT/'tmp/wp8_newton_correction_20261103'
c.REMOTE='/temp/zhixin/paper2_wp8_newton_correction_20261103'
BUILD='/temp/zhixin/paper2_wp8_background_v2_build_20261103'
CASES={'late_legacy_probe':('4T',4413,1,'legacy',32,20,True),
       'late_legacy_window':('4T',4413,10,'legacy',32,20,False),
       'late_equation_probe':('4T',4413,1,'legacy',32,20,True),
       'late_equation_probe_fine':('4T',4413,1,'legacy',32,20,True),
       'early_scaled':(None,0,10,'scaled',32,20,False),
       'nearT_scaled':('T',1190,10,'scaled',32,20,False),
       'late_scaled':('4T',4413,10,'scaled',32,20,False),
       'late_scaled_16':('4T',4413,10,'scaled',16,20,False),
       'forced_failure':(None,0,1,'scaled',32,0,False)}

def prepare(name):
    old,index,steps,criterion,factor,maxiter,probe=CASES[name]
    d=c.OUT/name; assert not d.exists();d.mkdir()
    isprobe=name.startswith('late_equation_probe');fine=name.endswith('_fine')
    suffix=('_probe_fine' if fine else '_probe') if isprobe else ''
    model='transport_wp8_background_v2'+suffix
    b=json.loads((c.OUT/(suffix.lstrip('_')+'_build' if isprobe else 'build')/'attempt_01/result.json').read_text())
    m=json.loads((c.ROOT/'model'/model/'source_manifest.json').read_text())
    assert b['build_exit']==0 and b['source_sha256']=={k:v['sha256'] for k,v in m.items()}
    for k,v in m.items():assert c.sha(c.ROOT/'model'/model/'source'/k)==v['sha256']
    base=c.ROOT/'analysis/wp7_shortwave_validation_20261103/runs/WP7-1D-background-T/params.input'
    dt=113.68588647695695;offset=0.;replay=None
    text=base.read_text().replace('Name = WP7-1D-background-T','Name = WP8-v2-'+name).replace('TEnd = 136423.0637723484',f'TEnd = {steps*dt:.17g}')
    text+='\n[WP8]\nCriterion = '+criterion+f'\nRoundoffFactor = {factor}\nMaxIterations = {maxiter}\nProbe = '+str(probe).lower()+'\n'
    if old:
        from analyze_wp8 import rows
        cache=c.ROOT/f'tmp/wp8_strict_background_20261103/{old}/post/profiles.npz'
        a=dict(np.load(cache));ev=rows(c.ROOT/f'analysis/wp8_strict_background_20261103/{old}/surface_events.csv');offset=float(ev[index]['time_s'])
        np.savetxt(d/'replay.dat',np.c_[a['z_m'],a['pw_mean'][index],a['x_mean'][index]],fmt='%.17g')
        text+='AllowDiagnosticReplay = true\nReplayProfile = replay.dat\n'
        replay=dict(source_run=old,source_index=index,source_profiles_sha256=c.sha(cache),replay_sha256=c.sha(d/'replay.dat'),lossless_restart=False,physical_time_offset_s=offset)
    (d/'params.input').write_text(text)
    gates=json.loads((c.ROOT/'analysis/wp8_zero_mode_20261103/background_protocol.json').read_text())['gates']
    buildroot=BUILD.replace('v2_build','v2'+suffix+'_build')
    p=dict(stage=name,authorization_sha256=c.sha(c.ROOT/'command/14.md'),executable=b['executable']['path'],elf_sha256=b['executable']['sha256'],remote_source=buildroot+'/source',source_hashes=b['source_sha256'],input_sha256=c.sha(d/'params.input'),wall_cap_s=120,memory_gib=2,cpus=1,archive_cap_mib=10,horizon_s=steps*dt,steps_cap=steps,gates=gates,replay=replay,physical_time_offset_s=offset,criterion=criterion,roundoff_factor=factor,expected_failure=name=='forced_failure',full_4T_authorized=False)
    c.save(d/'plan.json',p)

def launch(name):
    # Reuse reviewed resource-aware worker/archiver, stage the diagnostic profile first
    # through a narrowly scoped wrapper around the original launch's remote payload.
    original=c.remote
    def send(code):
        p=c.OUT/name/'replay.dat'
        if p.exists():
            import base64
            marker="for n,b in FILES.items():"  # launch already substitutes FILES; inject before resource sampling
            code=code.replace('def snap():',"(r/'replay.dat').write_bytes(base64.b64decode("+repr(base64.b64encode(p.read_bytes()).decode())+"))\nassert hashlib.sha256((r/'replay.dat').read_bytes()).hexdigest()=="+repr(c.sha(p))+"\ndef snap():")
            code='import hashlib\n'+code
        return original(code)
    c.remote=send
    try:c.launch(name)
    finally:c.remote=original

if __name__=='__main__':
    action,name=sys.argv[1:];assert name in CASES
    if action=='prepare':prepare(name)
    elif action=='launch':launch(name)
    elif action=='collect':c.collect(name)
    elif action=='status':c.status(name)
    else:raise ValueError(action)
