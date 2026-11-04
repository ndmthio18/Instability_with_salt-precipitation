"""WP9 Task B: freeze the single R96-M1 validation protocol. Registration only.

Writes the pending input, the frozen protocol and the registered admission criteria.
It starts nothing: `wp9_r96_campaign.py` refuses to launch unless the authorization
record this script names exists and matches.
"""
import json,math,sys
from pathlib import Path
import wp9_common as c

R48=c.ROOT/'analysis/wp7_shortwave_validation_20261103/runs/WP7-R48-M1-4T'
PENDING=c.OUT/'pending_r96'
NAME='WP9-R96-M1-4T'
DT=c.DT_S/2.0
STEPS=9600
OUTPUT_INTERVAL=2
ELF_SHA='112076759fdff0942000307f729c2bd3c7e96963229b43c1b67c1f12bbf9df81'
PROFILE_SHA='d721e605d126e555f92b43200aa73b999ed4e2b61397edd1e98c0d9fddb36575'
REMOTE='/tmp/zhixin_paper2_wp9_r96_20261104'
WP7_REMOTE='/tmp/zhixin_paper2_wp7_shortwave_20261103_r01'
ARCHIVE_HOME='~/paper2_wp9_r96_20261104'

# changing only these keys relative to the frozen R48-M1 input is itself a gate
ALLOWED_CHANGES={'Grid.Cells0':('48','96'),
                 'TimeLoop.DtInitial':(repr(c.DT_S),repr(DT)),
                 'TimeLoop.MaxTimeStepSize':(repr(c.DT_S),repr(DT)),
                 'TimeLoop.OutputInterval':('1',str(OUTPUT_INTERVAL)),
                 'Problem.Name':('WP7-R48-M1-4T',NAME)}


def read_input(path):
    """Section-qualified key/value view of a DuMuX .input file, order preserved."""
    section='';out=[]
    for line in Path(path).read_text().splitlines():
        s=line.split('#')[0].strip()
        if not s:out.append(('raw',line));continue
        if s.startswith('[') and s.endswith(']'):section=s[1:-1];out.append(('raw',line));continue
        if '=' in s:
            k,v=s.split('=',1);out.append(('kv',(f'{section}.{k.strip()}' if section else k.strip()),v.strip(),line))
        else:out.append(('raw',line))
    return out


def flat(path):
    return {e[1]:e[2] for e in read_input(path) if e[0]=='kv'}


def build_input():
    text=[]
    for e in read_input(R48/'params.input'):
        if e[0]=='raw':text.append(e[3] if len(e)>3 else e[1]);continue
        key,value=e[1],e[2]
        if key=='Grid.Cells0':value='96'
        elif key in ('TimeLoop.DtInitial','TimeLoop.MaxTimeStepSize'):value=repr(DT)
        elif key=='TimeLoop.OutputInterval':value=str(OUTPUT_INTERVAL)
        elif key=='Problem.Name':value=NAME
        text.append(f'{key.split(".",1)[1]} = {value}')
    return '\n'.join(text)+'\n'


def time_axis_evidence():
    """How far the dt/2 accumulation can drift from the R24/R48 accumulation."""
    a=0.0;coarse=[]
    for _ in range(4801):coarse.append(a);a+=c.DT_S
    b=0.0;fine=[]
    for _ in range(STEPS+1):fine.append(b);b+=DT
    worst=max(abs(coarse[k]-fine[2*k]) for k in range(4801))
    return dict(exact_halving=2*DT==c.DT_S,
                coarse_endpoint_s=coarse[-1],fine_endpoint_s=fine[-1],
                max_paired_representation_error_s=worst,
                accumulation_budget_s=STEPS*math.ulp(c.HORIZON_S),
                pairing_rule='frame j of R96 pairs with frame j of R24/R48 (step 2j at dt/2 versus step j at dt)',
                tolerance_s=STEPS*math.ulp(c.HORIZON_S),
                note=('dt/2 is an exact binary halving, so the mathematical accepted times coincide '
                      'identically; only the double accumulation differs. The 32-ULP same_time rule '
                      'used inside one protocol is too tight across protocols and is replaced by the '
                      'accumulation budget plus the exact frame-index identity.'))


def main():
    PENDING.mkdir(parents=True,exist_ok=True)
    base=flat(R48/'params.input')
    (PENDING/'params.input').write_text(build_input(),encoding='utf-8')
    new=flat(PENDING/'params.input')
    changed={k:(base.get(k),new.get(k)) for k in set(base)|set(new) if base.get(k)!=new.get(k)}
    assert set(changed)==set(ALLOWED_CHANGES),f'unexpected input changes: {sorted(set(changed)^set(ALLOWED_CHANGES))}'
    assert abs(float(new['TimeLoop.DtInitial'])*STEPS-c.HORIZON_S)<1e-6
    import shutil;shutil.copyfile(R48/'eigenprofile.dat',PENDING/'eigenprofile.dat')
    assert c.sha(PENDING/'eigenprofile.dat')==PROFILE_SHA,'A28 eigenprofile must be the original file'

    protocol=dict(
        schema_version=1,id='WP9-R96-M1',created_date='2026-11-04',
        authority='command/16.md',authority_sha256=c.sha(c.ROOT/'command/16.md'),
        reviewed_baseline='fb4a91b2646d09394e428d408e2ffc005d4c7713',
        decision='FROZEN_AND_PENDING_EXPLICIT_USER_AUTHORIZATION',launched=False,
        scope=['exactly one R96 M1 leg from t=0 to 4T or the first local x>=0.1'],
        forbidden=['R96-M0','R192','multiwavelength','full-domain','8T','restart','retry',
                   'any change of the R24/R48 solver formulation','any Paper 1 change'],
        physics=dict(case_id=int(c.CASE),width_m=c.WIDTH_M,height_m=c.HEIGHT_M,
            canonical_mode='A28, Neumann m=2, one wavelength across the box',
            nx_cells=96,nz_cells=200,grading1=-0.98,
            perturbation=dict(type='Eigenprofile',amplitude=1e-06,wave_number=1,
                              profile_sha256=PROFILE_SHA,
                              note='the original A28 eigenprofile file used by R24 and R48, byte-identical'),
            constitutive='unchanged from the frozen R48-M1 input; KrwMualemExponent = true',
            stopping=dict(horizon_s=c.HORIZON_S,horizon_over_T=4.0,
                          solubility_mole_fraction=c.SOLUBILITY,stop_at_solubility=True,
                          first_local_crossing='whichever comes first')),
        numerics=dict(dt_s=DT,dt_relative_to_R24_R48=0.5,steps=STEPS,fixed_timestep=True,
            output_interval=OUTPUT_INTERVAL,expected_frames=STEPS//OUTPUT_INTERVAL+1,
            executable=dict(sha256=ELF_SHA,path=f'{WP7_REMOTE}/bundle/bin/paper2_transport',
                            provenance='the identical bundle that produced R24 and R48; no rebuild'),
            cuts='the seven original conservative planes, byte-identical CutCoordinates',
            time_axis=time_axis_evidence(),
            deviation=dict(key='TimeLoop.OutputInterval',from_='1',to_=str(OUTPUT_INTERVAL),
                reason=('the per-step conservative diagnostics and sigma_top are written every accepted '
                        'step irrespective of OutputInterval (main.cc: transport.write and writeSigmaTop '
                        'precede the interval test), so interval 2 loses no diagnostic and no comparison '
                        'time: the shared accepted times with R24/R48 are exactly the even steps'),
                effect='4801 frames instead of 9601; peak scratch 21.5 GiB instead of 43 GiB')),
        resources=dict(host='sal',partition='debug',ranks=1,cpus=1,
            hard_wall_s=43200,projected_solver_wall_s=24621.1,projected_post_archive_wall_s=1200.0,
            projected_worker_wall_s=25821.1,projection_method=(
                'the measured R96 ten-step benchmark 2.5646979548 s/step extrapolated to 9600 steps with '
                'no efficiency credit, plus twice the measured R48-M1 post/archive cost scaled by the '
                '2.02x per-frame size ratio; the same extrapolation overestimated the realised R48 leg '
                'by 32 percent, so the hard wall is 1.67x the projection'),
            core_h_cap=12.0,
            minimum_free_scratch_bytes=40*1024**3,projected_peak_scratch_bytes=23*1024**3,
            minimum_available_memory_bytes=16*1024**3,projected_peak_rss_bytes=256*1024**2,
            scratch=f'{REMOTE} on the sal tmpfs',
            archive_destination=ARCHIVE_HOME,
            archive_note=('sal /tmp is a RAM-backed tmpfs; the field archive is moved to the NFS home '
                          'directory so it survives a reboot. The workstation has 19 GiB free, so the '
                          'field archive is NOT pulled back; only hash-verified derived products are.'),
            minimum_local_free_bytes=4*1024**3),
        admission=['authorization record present and naming this protocol',
                   'no existing WP9 run directory, locally or remotely',
                   'remote executable hash equals the R24/R48 executable hash',
                   'pending input hash equals the registered hash',
                   'eigenprofile hash equals the original A28 profile hash',
                   'the pending input differs from the frozen R48-M1 input in exactly five registered keys',
                   'remote free scratch, available memory and local free disk above the registered minima',
                   'the R24/R48 reduced spectra and their recorded hashes are present'],
        physical_acceptance=[
            'solver exit 0 with stop reason requested_horizon or a recorded solubility crossing',
            'every accepted step uses the fixed dt; no Newton-driven reduction',
            'accepted steps equal 9600 unless the solubility stop fired first',
            'maximum regional salt residual <= 1e-14 m and water residual <= 1e-11 m (the WP7 thresholds)',
            'independent VTU inventory, flux reintegration and budget identity reproduce the online budget',
            'initial state matches the registered initial condition',
            'no NaN or non-finite state anywhere in the diagnostics',
            'frame count equals 4801 and every frame hash verifies at both ends',
            'the paired time axis satisfies the registered accumulation budget at every shared frame'],
        postprocessing=['wp8_fields.py reduction, identical to the R24/R48 products',
                        'wp7_spatial.py metrics, identical to the R48 products',
                        'both run against hash-verified frames and audited'],
        prediction=dict(
            registered_before_the_run=True,
            question='is the late short-wave energy seeded by discretization noise or generated nonlinearly',
            measured_R48_over_R24_m10_energy_ratio_at_4T=1.3121e-11/1.0085e-08,
            measured_R48_over_R24_unperturbed_total_energy_ratio_at_4T=7.2045e-14/2.2611e-11,
            noise_seeded_expectation=('m=10 energy at 4T falls from R48 to R96 by a comparable large '
                                      'factor (order 10^2 or more) and the m=2 fraction rises above the '
                                      'R48 value of 0.357'),
            nonlinear_cascade_expectation=('m=10 energy at 4T agrees between R48 and R96 within about one '
                                           'order of magnitude and its fraction settles near the R48 0.265'),
            discriminating_threshold='a factor of 10 in the m=10 energy ratio separates the two readings',
            caveat=('this is a two-point comparison with a simultaneous dt change; it discriminates the two '
                    'readings, it does not establish a convergence order')),
        known_confounds=[
            dict(name='simultaneous timestep change',
                 detail=('R96 halves dx and dt together, so an R48-R96 difference cannot be attributed to '
                         'horizontal resolution alone. The only existing bound is the WP6 R24 dt/2 pair to '
                         '2T: growth amplification -0.0999 percent, surface weighted sigma -0.1145 percent, '
                         'vertical sigma profile relative L2 0.0389 percent, detection time +0.0972 percent. '
                         'That bound does not extend to 3T-4T, to R96, or to the 1 mm transport difference, '
                         'where the same study found a 68 percent mid-window change under dt halving.'),
                 consequence='growth and modal conclusions carry a ~0.1 percent temporal prior; transport conclusions carry none'),
            dict(name='domain width',
                 detail=('the box is one A28 wavelength wide with no-flow sides, so only Neumann modes '
                         'm = 1,2,... exist and the LSA global critical mode nx=5 (wavelength 0.12 m) '
                         'cannot appear at all'),
                 consequence='no statement about the unbounded-domain wavelength selection follows'),
            dict(name='vertical resolution',
                 detail='nz=200 with Grading1=-0.98 is held fixed at every resolution and is not varied here',
                 consequence='vertical truncation is common-mode and untested'),
            dict(name='executable provenance',
                 detail='the strict one-cell background uses a different executable; only R24/R48/R96 share this one',
                 consequence='background-referenced differences retain the cross-executable gap recorded in command/15')],
        failure_policy=['preserve the last accepted state, the solver log and all provenance',
                        'no automatic retry and no restart',
                        'report the failure and stop'])
    c.save(c.OUT/'r96_protocol.json',protocol)

    admission=dict(schema_version=1,protocol_sha256=c.sha(c.OUT/'r96_protocol.json'),
        input_sha256=c.sha(PENDING/'params.input'),
        eigenprofile_sha256=c.sha(PENDING/'eigenprofile.dat'),
        reference_input=str(R48.relative_to(c.ROOT)/'params.input'),
        reference_input_sha256=c.sha(R48/'params.input'),
        input_key_changes={k:dict(frm=v[0],to=v[1]) for k,v in changed.items()},
        allowed_key_changes=sorted(ALLOWED_CHANGES),
        launch_requires='analysis/wp9_r96_convergence_20261104/user_authorization.json',
        status='PENDING')
    c.save(c.OUT/'r96_admission.json',admission)
    print(json.dumps(dict(input_sha256=admission['input_sha256'],
                          protocol_sha256=admission['protocol_sha256'],
                          changes=admission['input_key_changes'],
                          time_axis=protocol['numerics']['time_axis'],
                          projected_worker_wall_h=protocol['resources']['projected_worker_wall_s']/3600,
                          hard_wall_h=protocol['resources']['hard_wall_s']/3600),indent=2))


if __name__=='__main__':main()
