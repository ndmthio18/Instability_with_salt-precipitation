#pragma once
#include "two_front.hh"
#include <cvode/cvode.h>
#include <cvode/cvode_ls.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sundials/sundials_context.h>

namespace twofront {
class Solver {
    SUNContext ctx=nullptr;N_Vector state=nullptr,constraints=nullptr;
    SUNMatrix matrix=nullptr;SUNLinearSolver linear=nullptr;void* mem=nullptr;
    const Model& model;double rtol;std::string error;
    static void check(int value,const char* name){if(value<0)throw std::runtime_error(std::string(name)+":"+std::to_string(value));}
    static int eval(realtype t,N_Vector a,N_Vector b,void* ptr){
        auto& s=*static_cast<Solver*>(ptr);
        try {Vec y(N_VGetArrayPointer(a),N_VGetArrayPointer(a)+s.model.size());auto dy=s.model.rhs(y,t);
             std::copy(dy.begin(),dy.end(),N_VGetArrayPointer(b));return 0;
        }catch(const std::exception& e){s.error=e.what();return 1;}
    }
    static int weights(N_Vector a,N_Vector w,void* ptr){
        auto& s=*static_cast<Solver*>(ptr);const auto& m=s.model;
        double b=NV_Ith_S(a,m.bottomIndex()),h=NV_Ith_S(a,m.topIndex())-b;
        if(!(h>0) || (m.l && !(b>0)))return -1;
        for(int i=0;i<m.n;++i){double v=h*m.widths[i],solid=std::clamp(NV_Ith_S(a,i),0.,v*(1-1e-12)),pore=v-solid;
            NV_Ith_S(w,i)=1/(1e-12*v+s.rtol*solid);
            NV_Ith_S(w,m.n+i)=1/(1e-12*m.p.rho*pore+s.rtol*std::abs(NV_Ith_S(a,m.n+i)));
        }
        for(int i=0;i<m.l;++i)NV_Ith_S(w,2*m.n+i)=1/(1e-12*m.p.rho*b*m.lw[i]+s.rtol*std::abs(NV_Ith_S(a,2*m.n+i)));
        for(int i=m.topIndex();i<m.size();++i)NV_Ith_S(w,i)=1/((i<m.ledger()?1e-14:1e-14)+s.rtol*std::abs(NV_Ith_S(a,i)));
        return 0;
    }
public:
    Solver(const Model& m,const Vec& initial,double tolerance=1e-8,double maxstep=30,double start=0):model(m),rtol(tolerance){
        check(SUNContext_Create(nullptr,&ctx),"context");state=N_VNew_Serial(m.size(),ctx);constraints=N_VClone(state);
        std::copy(initial.begin(),initial.end(),N_VGetArrayPointer(state));N_VConst(0.,constraints);
        for(int i=0;i<m.n;++i)NV_Ith_S(constraints,i)=1;
        NV_Ith_S(constraints,m.topIndex())=1;NV_Ith_S(constraints,m.bottomIndex())=1;
        mem=CVodeCreate(CV_BDF,ctx);check(CVodeInit(mem,eval,start,state),"init");check(CVodeSetUserData(mem,this),"data");
        check(CVodeWFtolerances(mem,weights),"weights");check(CVodeSetConstraints(mem,constraints),"constraints");
        check(CVodeSetMaxNumSteps(mem,300000),"steps");check(CVodeSetMaxStep(mem,maxstep),"maxstep");
        matrix=SUNDenseMatrix(m.size(),m.size(),ctx);linear=SUNLinSol_Dense(state,matrix,ctx);
        check(CVodeSetLinearSolver(mem,linear,matrix),"linear");
    }
    ~Solver(){if(mem)CVodeFree(&mem);if(linear)SUNLinSolFree(linear);if(matrix)SUNMatDestroy(matrix);if(state)N_VDestroy(state);if(constraints)N_VDestroy(constraints);if(ctx)SUNContext_Free(&ctx);}
    Vec advance(double t){double reached=0;int flag=CVode(mem,t,state,&reached,CV_NORMAL);
        if(flag<0)throw std::runtime_error("CVODE "+std::to_string(flag)+":"+error);
        Vec y(N_VGetArrayPointer(state),N_VGetArrayPointer(state)+model.size());model.validate(y,t);return y;}
    long steps()const{long n;CVodeGetNumSteps(mem,&n);return n;}
};
}
