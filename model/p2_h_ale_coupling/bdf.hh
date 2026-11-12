#pragma once
#include "ale.hh"
#include <cvode/cvode.h>
#include <cvode/cvode_ls.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sundials/sundials_context.h>

namespace p2h {
// CVODE BDF/Newton normally advances E=eps*volume and Q=rho*C*E directly.
// In the explicit complete-dissolution diagnostic it instead advances the
// non-negative solid volume S=volume-E, preventing epsilon from crossing 1.
// A diagonal
// error-weight transform uses M3's tolerance coefficients; its norm is not
// identical to SciPy's primitive-variable norm. Optional diagnostic states
// integrate every salt-face flux and every cell reaction with the same BDF
// history as E and Q, so accepted-step local balances need no endpoint-flux
// approximation. Every trial owns its solver.
class BDF {
    SUNContext context=nullptr;
    N_Vector state=nullptr, tolerance=nullptr, constraints=nullptr;
    SUNMatrix matrix=nullptr;
    SUNLinearSolver linear=nullptr;
    void* memory=nullptr;
    ALE model;
    double relativeTolerance;
    bool integrateFluxes=false;
    bool solidMode=false;
    int stateSize=0;
    Vec integratedSaltFaces_, integratedReaction_;
    std::string error;
    static int evaluate(realtype t, N_Vector y, N_Vector dy, void* pointer) {
        auto& self=*static_cast<BDF*>(pointer);
        try {
            Vec input(N_VGetArrayPointer(y),N_VGetArrayPointer(y)+self.model.size());
            auto primitive=self.primitive(input);
            auto f=self.model.rhs(primitive,t);
            const auto flux=self.model.flux(primitive,t);
            for(int i=0;i<self.model.n;++i) {
                const double dE=flux.mesh[i]-flux.mesh[i+1]-flux.reaction[i]*flux.volume[i]/self.model.p.crystal;
                f[i]=self.solidMode ? self.model.widths[i]*flux.speed-dE : dE;
                f[self.model.n+i]=flux.salt[i]-flux.salt[i+1]-flux.reaction[i]*flux.volume[i];
            }
            std::copy(f.begin(),f.end(),N_VGetArrayPointer(dy));
            if(self.integrateFluxes) {
                const int offset=self.model.size();
                for(int face=0;face<=self.model.n;++face)
                    NV_Ith_S(dy,offset+face)=flux.salt[face];
                for(int cell=0;cell<self.model.n;++cell)
                    NV_Ith_S(dy,offset+self.model.n+1+cell)=flux.reaction[cell]*flux.volume[cell];
            }
            return 0;
        } catch(const std::exception& e) { self.error=e.what(); return 1; }
    }
    Vec primitive(const Vec& conservative) const {
        auto y=conservative;
        for(int i=0;i<model.n;++i) {
            const double volume=conservative[2*model.n]*model.widths[i];
            // CVODE can evaluate an unaccepted Newton iterate infinitesimally
            // outside a one-sided constraint.  Project only that evaluation
            // back into the primitive-variable domain; accepted states remain
            // governed by CVodeSetConstraints and are never modified here.
            const double solid=solidMode
                ? std::max(0.0,std::min(conservative[i],volume*(1-1e-14)))
                : volume-conservative[i];
            const double pore=solidMode ? volume-solid : conservative[i];
            y[i]=pore/volume;
            y[model.n+i]=conservative[model.n+i]/(model.p.rho*pore);
        }
        return y;
    }
    static int weights(N_Vector y,N_Vector w,void* pointer) {
        auto& self=*static_cast<BDF*>(pointer);
        const int n=self.model.n; const double h=NV_Ith_S(y,2*n);
        if(!(h>0)) return -1;
        for(int i=0;i<n;++i) {
            const double primary=NV_Ith_S(y,i),q=NV_Ith_S(y,n+i),vol=h*self.model.widths[i];
            if(!(vol>0) || !std::isfinite(primary) || !std::isfinite(q)) return -1;
            const double boundedPrimary=self.solidMode
                ? std::max(0.0,std::min(primary,vol*(1-1e-14))) : primary;
            const double pore=self.solidMode ? vol-boundedPrimary : boundedPrimary;
            if(!(pore>0)) return -1;
            NV_Ith_S(w,i)=1/(1e-12*vol+self.relativeTolerance*std::abs(boundedPrimary));
            NV_Ith_S(w,n+i)=1/(1e-12*self.model.p.rho*pore+self.relativeTolerance*std::abs(q));
        }
        for(int i=2*n;i<self.model.size();++i)
            NV_Ith_S(w,i)=1/((i==2*n?1e-14:1e-12)+self.relativeTolerance*std::abs(NV_Ith_S(y,i)));
        for(int i=self.model.size();i<self.stateSize;++i)
            NV_Ith_S(w,i)=1/(1e-16+self.relativeTolerance*std::abs(NV_Ith_S(y,i)));
        return 0;
    }
    static void checked(int flag,const char* name) {
        if(flag<0) throw std::runtime_error(std::string(name)+" flag="+std::to_string(flag));
    }
public:
    std::function<void(const Vec&,double)> rejectedState;
    BDF(const ALE& source,const Vec& initial,double maxStep=60,double relative=1e-8,
        bool integrateAcceptedFluxes=false)
        : model(source), relativeTolerance(relative), integrateFluxes(integrateAcceptedFluxes),
          solidMode(source.p.allowCompleteDissolution),
          stateSize(model.size()+(integrateFluxes ? 2*model.n+1 : 0)) {
        checked(SUNContext_Create(nullptr,&context),"context");
        state=N_VNew_Serial(stateSize,context); tolerance=N_VClone(state); constraints=N_VClone(state);
        N_VConst(0.0,state);
        N_VConst(0.0,constraints);
        std::copy(initial.begin(),initial.end(),N_VGetArrayPointer(state));
        for(int i=0;i<model.n;++i) {
            const double volume=initial[2*model.n]*model.widths[i];
            const double pore=initial[i]*volume;
            NV_Ith_S(state,i)=solidMode ? volume-pore : pore;
            NV_Ith_S(state,model.n+i)=model.p.rho*initial[model.n+i]*pore;
            if(solidMode) NV_Ith_S(constraints,i)=1.0;
        }
        for(int i=0;i<stateSize;++i)
            NV_Ith_S(tolerance,i)=i==2*model.n?1e-14:(i<model.size()?1e-12:1e-16);
        memory=CVodeCreate(CV_BDF,context);
        checked(CVodeInit(memory,evaluate,0,state),"init");
        checked(CVodeSetUserData(memory,this),"user data");
        checked(CVodeWFtolerances(memory,weights),"conservative error weights");
        if(solidMode) checked(CVodeSetConstraints(memory,constraints),"solid non-negativity");
        checked(CVodeSetMaxNumSteps(memory,200000),"step limit");
        checked(CVodeSetMaxStep(memory,maxStep),"maximum step");
        matrix=SUNDenseMatrix(stateSize,stateSize,context);
        linear=SUNLinSol_Dense(state,matrix,context);
        checked(CVodeSetLinearSolver(memory,linear,matrix),"linear solver");
    }
    BDF(const BDF&)=delete;
    ~BDF() {
        if(memory) CVodeFree(&memory);
        if(linear) SUNLinSolFree(linear);
        if(matrix) SUNMatDestroy(matrix);
        if(state) N_VDestroy(state);
        if(tolerance) N_VDestroy(tolerance);
        if(constraints) N_VDestroy(constraints);
        if(context) SUNContext_Free(&context);
    }
    Vec advance(double time) {
        double reached=0;
        const int flag=CVode(memory,time,state,&reached,CV_NORMAL);
        if(flag<0) throw std::runtime_error("BDF failed "+std::to_string(flag)+": "+error);
        Vec y(N_VGetArrayPointer(state),N_VGetArrayPointer(state)+model.size());
        y=primitive(y);
        try { model.validate(y); } catch(...) { if(rejectedState) rejectedState(y,reached); throw; }
        if(integrateFluxes) {
            const auto* q=N_VGetArrayPointer(state)+model.size();
            integratedSaltFaces_.assign(q,q+model.n+1);
            integratedReaction_.assign(q+model.n+1,q+2*model.n+1);
        }
        return y;
    }
    const Vec& integratedSaltFaces() const { return integratedSaltFaces_; }
    const Vec& integratedReaction() const { return integratedReaction_; }
    long steps() const { long v; CVodeGetNumSteps(memory,&v); return v; }
};
}
