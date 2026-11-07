"""Build an isolated application by deterministic edits of frozen P1-R2.
Run on the HPCC compute node; only the supplied destination is written.
"""
import pathlib,sys,shutil
root=pathlib.Path(__file__).resolve().parents[2]
here=pathlib.Path(__file__).resolve().parent
target=pathlib.Path(sys.argv[1])
shutil.copytree(root/'model/p1_r2_conservative_richards/source/appl',target,dirs_exist_ok=True)
app=target/'simpleRichardsnc'
src=(app/'main.cc').read_text()
def replace(old,new):
    global src
    assert src.count(old)==1,(old,src.count(old))
    src=src.replace(old,new)
replace('#include "model.hh"','#include "model.hh"\n#include "surface.hh"\n#include "film.hh"')
replace('else if (mode == "no_flux")','else if (mode == "surface_film")\n                { values[conti0EqIdx]=filmFlux_; values[conti1EqIdx]=0; }\n                else if (mode == "no_flux")')
replace('        Scalar nonWettingReferencePressure() const',(here/'methods.inc').read_text()+'\n        Scalar nonWettingReferencePressure() const')
replace('    auto xOld = x;','    problem->prepareSurfaceFlux(x);\n    problem->surfaceStep(x,0.0,0.0);\n    auto xOld = x;')
replace('        gridVariables->update(x);','        problem->surfaceStep(x,timeBefore+acceptedDt,acceptedDt);\n        gridVariables->update(x);')
replace('        // set new dt as suggested by the newton solver','        problem->prepareSurfaceFlux(x);\n        // set new dt as suggested by the newton solver')
(app/'main.cc').write_text(src)
shutil.copy(here/'surface.hh',app/'surface.hh')
# Extract only thermodynamics/gas-film code; exclude all reaction/profile code.
old=(root/'model/p2_d_coupling/cpp/pratcrust.hh').read_text()
part=old.split('namespace prat {',1)[1].split('// ------------------------------------------------------ crust parameters',1)[0]
(app/'film.hh').write_text('#pragma once\n#include <cmath>\n#include <stdexcept>\n#include <string>\nnamespace p2gfilm {\n'+part+'\n}\n')
