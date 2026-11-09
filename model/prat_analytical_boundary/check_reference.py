"""Independent scalar reference, no third-party packages; run on compute node."""
import csv,json,math,pathlib
r=pathlib.Path.cwd()
c=.1*.05844/(.9*.01801528+.1*.05844); e=.43; kr=.001; beta=2; ds=1.3e-9; dsub=25e-6
deff=math.sqrt(e)*ds; av=36*(1-e)/dsub
k=beta*(1-e)*(1-c)*kr/(e*deff)
a=c*1e-5/(2165*(1-e)*(1-c)); h0=3.5e-6
da=h0*h0*(1-c)*av*kr/(e*deff)
da_star=beta**2*(1-e)**2*(1-c)*kr/(e*deff*av)
da_top=k*h0
eff_m1=da_top/(1+da_top)
z=math.sqrt(da_star)*math.tanh(math.sqrt(da))
eff_m2=z/(1+z)
eff_m3=math.sqrt(da_star)/(1+math.sqrt(da_star))
lo=h0; hi=h0+a*1000
for _ in range(100):
    mid=(lo+hi)/2
    if mid-h0+math.log(mid/h0)/k>a*1000: hi=mid
    else: lo=mid
exact=(lo+hi)/2
rows=list(csv.DictReader((r/'unit/reference.csv').open()))
errors=[]
for row in rows:
    dt=float(row['dt']); h=h0
    for _ in range(round(1000/dt)): h+=dt*a*(k*h)/(1+k*h)
    cpp=float(row['h_final']); assert abs(cpp-h)<1e-17
    for key,value in [('da',da),('da_star',da_star),('da_top',da_top),('efficiency_m1',eff_m1),('efficiency_m2',eff_m2),('efficiency_m3',eff_m3)]:
        assert abs(float(row[key])-value)<=2e-15*max(1,abs(value)),(key,row[key],value)
    errors.append(abs(cpp-exact))
assert all(1.95<errors[i]/errors[i+1]<2.05 for i in range(3))
(r/'unit/reference.json').write_text(json.dumps({'status':'PASS','exact_h_m':exact,'errors_m':errors,'ratios':[errors[i]/errors[i+1] for i in range(3)],
    'activation_reference':{'h_m':h0,'da':da,'da_star':da_star,'da_top':da_top,'efficiency_m1':eff_m1,'efficiency_m2':eff_m2,'efficiency_m3':eff_m3,'m3_valid_da_ge_10':da>=10}},indent=2))
print('PASS independent M1/M2/M3 reference, first-order timestep convergence, C++/Python')
