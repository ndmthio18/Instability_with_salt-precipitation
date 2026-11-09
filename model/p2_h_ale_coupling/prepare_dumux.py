"""Create a separate P1-R2 application; frozen P2-G and upstream stay read-only."""
import pathlib
import shutil
import sys

here=pathlib.Path(__file__).resolve().parent
root=here.parents[1]
target=pathlib.Path(sys.argv[1])
prefix=pathlib.Path(sys.argv[2])
source=root/'model/p1_r2_conservative_richards/source/appl'
# The HPCC production image provides Python 3.6, before copytree gained
# dirs_exist_ok. Merge the frozen application tree explicitly and preserve
# metadata instead of depending on the client-side Python version.
for path in source.rglob('*'):
    destination=target/path.relative_to(source)
    if path.is_dir():
        destination.mkdir(parents=True,exist_ok=True)
    else:
        destination.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(str(path),str(destination))
app=target/'simpleRichardsnc'
src=(app/'main.cc').read_text()
def replace(old,new):
    global src
    assert src.count(old)==1,(old,src.count(old))
    src=src.replace(old,new)
replace('#include "model.hh"','#include "model.hh"\n#include "bdf.hh"\n#include "surface.hh"\n#include "film.hh"\n#include "donor.hh"\n#include "progressive.hh"')
replace('if (onUpperBoundary_(globalPos))','if (topFace_(scvf))')
replace('else if (mode == "no_flux")','else if (mode == "p2h_ale")\n                { values[conti0EqIdx]=heldLiquid_; values[conti1EqIdx]=heldSalt_; }\n                else if (mode == "no_flux")')
replace('if (mode != "prescribed_components")','if (mode != "prescribed_components" && mode != "p2h_ale")')
replace('        Scalar nonWettingReferencePressure() const',(here/'soil_methods.inc').read_text()+'\n        Scalar nonWettingReferencePressure() const')
replace('    auto xOld = x;','    problem->initializeSurface(x);\n    if (problem->controlledDonorComplete()) return 0;\n    auto xOld = x;')
replace('        nonLinearSolver.solve(x, *timeLoop);','        problem->advanceCoupled(x,nonLinearSolver,*timeLoop,*gridVariables);')
replace('            || timeLoop->finished() || transport.shouldStop())',
        '            || timeLoop->finished() || transport.shouldStop() || problem->activationOnlyComplete())')
replace('    } while (!timeLoop->finished() && !transport.shouldStop());',
        '    } while (!timeLoop->finished() && !transport.shouldStop() && !problem->activationOnlyComplete());')
replace('transport.finish(transport.shouldStop() ? transport.thresholdReason() : "requested_horizon");',
        'transport.finish(problem->activationOnlyComplete() ? "conservative_activation_complete" : (transport.shouldStop() ? transport.thresholdReason() : "requested_horizon"));')
(app/'main.cc').write_text(src)
for name in ('ale.hh','bdf.hh','donor.hh','progressive.hh','grid_diagnostics.inc'):
    shutil.copyfile(here/name,app/name)
shutil.copyfile(root/'model/p2_g_surface_precipitation/surface.hh',app/'surface.hh')
old=(root/'model/p2_d_coupling/cpp/pratcrust.hh').read_text()
part=old.split('namespace prat {',1)[1].split('// ------------------------------------------------------ crust parameters',1)[0]
(app/'film.hh').write_text('#pragma once\n#include <cmath>\n#include <stdexcept>\n#include <string>\nnamespace p2gfilm {\n'+part+'\n}\n')
cmake=app/'CMakeLists.txt'
libs=['cvode','nvecserial','sunlinsoldense','sunmatrixdense']
with cmake.open('a') as f:
    f.write('\ntarget_include_directories(paper2_transport PRIVATE "'+str(prefix/'include')+'")\n')
    f.write('target_link_libraries(paper2_transport PRIVATE '+
            ' '.join('"'+str(prefix/'lib64'/('libsundials_'+lib+'.a'))+'"' for lib in libs)+')\n')
