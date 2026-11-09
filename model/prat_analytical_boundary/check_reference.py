"""Independent scalar reference, no third-party packages; run on compute node."""
import csv,json,math,pathlib
r=pathlib.Path.cwd()
c=.1*.05844/(.9*.01801528+.1*.05844); e=.43
k=2*(1-e)*(1-c)*.0023/(e*math.sqrt(e)*1.3e-9)
a=c*1e-5/(2165*(1-e)*(1-c)); h0=3.5e-6
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
    errors.append(abs(cpp-exact))
assert all(1.95<errors[i]/errors[i+1]<2.05 for i in range(3))
(r/'unit/reference.json').write_text(json.dumps({'status':'PASS','exact_h_m':exact,'errors_m':errors,'ratios':[errors[i]/errors[i+1] for i in range(3)]},indent=2))
print('PASS independent integral reference, first-order timestep convergence, C++/Python')
