"""Prepare a future run's input and provenance. This script NEVER runs a solver.

Requires explicit paths; creates a fresh directory and refuses overwrites.
Use only after a separate decision to prepare numerical validation cases.
"""
import argparse,hashlib,json,math,re,shutil
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def parse(text):
    values={};section=''
    for line in text.splitlines():
        line=line.split('#',1)[0].strip()
        if not line:continue
        if line.startswith('[') and line.endswith(']'):section=line[1:-1]
        elif '=' in line:
            k,v=line.split('=',1);key=section+'.'+k.strip()
            if key in values:raise ValueError('Duplicate input key: '+key)
            values[key]=v.strip()
    return values
def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input',type=Path,required=True)
    ap.add_argument('--profile',type=Path)
    ap.add_argument('--executable',type=Path,required=True)
    ap.add_argument('--build-result',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--case-id',required=True)
    ap.add_argument('--role',choices=['uniform_2d','perturbed_2d'],required=True)
    args=ap.parse_args()
    if args.output.exists():raise ValueError('Run directory must not exist')
    text=args.input.read_text(encoding='utf-8');v=parse(text)
    build=json.loads(args.build_result.read_text())
    if build.get('build_exit')!=0 or sha(args.executable)!=build['executable']['sha256']:
        raise ValueError('Executable does not match successful build record')
    for name,h in build['source_sha256'].items():
        if sha(ROOT/'model/transport_v1/source'/name)!=h:raise ValueError('Source changed since build: '+name)
    def val(k):
        x=float(v['Problem.'+k])
        if not math.isfinite(x):raise ValueError(k)
        return x
    # Require explicit fluid constants: old defaults caused an audited conflict.
    rho0=val('PureWaterDensity');mw=val('MolarMassH2O');ms=val('MolarMassNaCl')
    x0=val('SaltMoleFraction');lin=val('DensityLinearCoefficient');quad=val('DensityQuadraticCoefficient')
    if v.get('Problem.KrwMualemExponent','').lower()!='true':raise ValueError('Exponent krw required')
    if args.role=='uniform_2d':
        if val('Amplitude')!=0 or v['Problem.PerturbationType'] not in ('Constant','Eigenprofile'):
            raise ValueError('Uniform reference requires zero amplitude and Constant or Eigenprofile initialization')
    elif v['Problem.PerturbationType']!='Eigenprofile':raise ValueError('Perturbed role requires Eigenprofile')
    if v['Problem.PerturbationType']=='Eigenprofile' and not args.profile:raise ValueError('Explicit --profile required')
    zs=[float(x) for x in v['Grid.Positions1'].split()];xs=[float(x) for x in v['Grid.Positions0'].split()]
    H=zs[-1]-zs[0];W=xs[-1]-xs[0]
    if len(zs)!=2 or len(xs)!=2 or zs[0]!=0 or not (H>0 and W>0):raise ValueError('Rectangular domain with bottom z=0 required')
    phi=val('Porosity');E=val('EvaporationRate');mu=val('ReferenceViscosity');g=val('GravitationalAcceleration')
    if not (0<phi<1 and E>0 and mu>0 and g>0 and 0<x0<.1):raise ValueError('Invalid reference physical parameters')
    mass_fraction=x0*ms/((1-x0)*mw+x0*ms);rho=rho0+1000*mass_fraction*(lin+quad*mass_fraction)
    groups={'R_E':val('Permeability')*rho*g/(mu*E),'beta':phi**1.5*val('SaltDiffusivity')/(E*H),
            'gamma':val('PcAlpha')*rho*g*H,'P_hat_B':val('WaterLevel')/H,'n':val('PcN'),
            'Swr':val('Swr'),'L':val('SweL'),'W_over_H':W/H,'reference_time_s':phi*H/E}
    if 'Transport.StopAtSolubility' in v and v['Transport.StopAtSolubility'].lower()!='true':
        raise ValueError('Preparation requires threshold stopping enabled')
    if 'Transport.SolubilityMoleFraction' in v and float(v['Transport.SolubilityMoleFraction'])!=.1:
        raise ValueError('Preparation protocol uses solubility threshold 0.1')
    # A fresh input contains the complete actual parameters, including defaults
    # chosen here, rather than relying on implicit values in a future executable.
    if args.profile:v['Problem.EigenprofilePath']='eigenprofile.dat'
    v['Transport.StopAtSolubility']='true';v['Transport.SolubilityMoleFraction']='0.1'
    sections={}
    for key,value in v.items():
        section,key=key.split('.',1);sections.setdefault(section,{})[key]=value
    normalized='\n\n'.join('['+s+']\n'+'\n'.join(k+' = '+v for k,v in entries.items()) for s,entries in sections.items())+'\n'
    args.output.mkdir(parents=True)
    (args.output/'params.input').write_text(normalized,encoding='utf-8')
    if args.profile:shutil.copyfile(args.profile,args.output/'eigenprofile.dat')
    manifest={'schema_version':1,'status':'prepared_not_executed','case_id':args.case_id,'role':args.role,
              'conserved_quantity':'phi*Sw*x normalized component inventory','dimensionless':groups,
              'reference_density_kg_m3':rho,'parameters':v,'initialization_time_s':0,
              'source_input_sha256':sha(args.input),'actual_input_sha256':sha(args.output/'params.input'),
              'eigenprofile_sha256':sha(args.profile) if args.profile else None,
              'executable_sha256':sha(args.executable),'executable_path_observed':str(args.executable.resolve()),
              'build_result_sha256':sha(args.build_result),'build_source_sha256':build['source_sha256'],
              'independent_1d_model_implemented':False,'solver_executed':False}
    (args.output/'run_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print('Prepared only; no solver executed:',args.output)
if __name__=='__main__':main()
