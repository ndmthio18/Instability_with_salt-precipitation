"""Generate isolated DuMux application from pinned, immutable repository sources."""
import pathlib
import shutil
import sys

here=pathlib.Path(__file__).resolve().parent
repo=here.parents[1]
dest=pathlib.Path(sys.argv[1])
source=repo/'model/p1_r2_conservative_richards/source/appl'
for path in source.rglob('*'):
    target=dest/path.relative_to(source)
    if path.is_dir(): target.mkdir(parents=True,exist_ok=True)
    else:
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(str(path),str(target))
app=dest/'simpleRichardsnc'
old=(repo/'model/p2_h_ale_coupling/soil_methods.inc').read_text()
diag=(repo/'model/p2_h_ale_coupling/grid_diagnostics.inc').read_text()
def between(text,start,end):
    assert text.count(start)==1 and text.count(end)==1,(start,end)
    return text[text.index(start):text.index(end)]
helpers=between(old,'        bool topFace_','        struct SurfaceState')
helpers+=between(old,'        struct BottomFlux','        SurfaceState acceptedSurface_')
helpers+=between(old,'        double area_()','        double film_(')
helpers+=between(diag,'template<class Solution> std::pair<std::vector<double>','template<class Solution> void captureTransport_')
precip=between(diag,'template<class Solution,class State,class Solver,class TimeLoop>','std::ofstream gridBudget_')
precip=precip.replace('getParam<std::string>("P2H.DonorMode","top_only")=="progressive_top_boundary_filling"','true')
for a,b in [('ale_->p.rho','params_.rho()'),('ale_->p.birth','params_.porosity'),('ale_->p.crystal','params_.crystalDensity'),('ale_->p.ceq','params_.equilibriumMassFraction'),('P2H.MinimumSupplySaturation','PratBoundary.MinimumSupplySaturation')]:
    precip=precip.replace(a,b)
methods=(here/'methods.inc').read_text().replace('// GENERATED_HELPERS',helpers+precip)
src=(app/'main.cc').read_text()
def replace(a,b):
    global src
    assert src.count(a)==1,a
    src=src.replace(a,b)
replace('#include "model.hh"','#include "model.hh"\n#include "boundary.hh"\n#include "film.hh"')
replace('if (onUpperBoundary_(globalPos))','if (topFace_(scvf))')
replace('else if (mode == "no_flux")','else if (mode == "prat_analytical") { values[conti0EqIdx]=heldLiquid_; values[conti1EqIdx]=heldSalt_; }\n                else if (mode == "no_flux")')
replace('if (mode != "prescribed_components")','if (mode != "prescribed_components" && mode != "prat_analytical")')
replace('        Scalar nonWettingReferencePressure() const',methods+'\n        Scalar nonWettingReferencePressure() const')
replace('    auto xOld = x;','    problem->initializeSurface(x);\n    auto xOld = x;')
replace('        nonLinearSolver.solve(x, *timeLoop);','        problem->advanceCoupled(x,nonLinearSolver,*timeLoop,*gridVariables);')
(app/'main.cc').write_text(src)
for origin,name in [(here/'boundary.hh','boundary.hh'),(repo/'model/p2_f_analytical_growth/cpp/analyticalgrowth.hh','analyticalgrowth.hh'),(repo/'model/p2_h_ale_coupling/progressive.hh','progressive.hh')]: shutil.copy2(str(origin),str(app/name))
film=(repo/'model/p2_d_coupling/cpp/pratcrust.hh').read_text().split('namespace prat {',1)[1].split('// ------------------------------------------------------ crust parameters',1)[0]
(app/'film.hh').write_text('#pragma once\n#include <cmath>\n#include <stdexcept>\n#include <string>\nnamespace p2gfilm {\n'+film+'\n}\n')
assert 'bdf.hh' not in src and 'ale_->' not in src
print('Generated analytical application without crust grid or SUNDIALS')
