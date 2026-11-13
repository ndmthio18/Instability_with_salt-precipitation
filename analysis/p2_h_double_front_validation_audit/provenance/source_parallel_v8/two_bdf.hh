#pragma once
#include "two_front.hh"
#include <cvode/cvode.h>
#include <cvode/cvode_ls.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sundials/sundials_context.h>
#include <limits>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <atomic>

namespace audit28 {
class Solver {
    SUNContext ctx=nullptr;N_Vector state=nullptr,constraints=nullptr;
    SUNMatrix matrix=nullptr;SUNLinearSolver linear=nullptr;void* mem=nullptr;
    const Model& model;double rtol;std::string error;std::ofstream diagnostics;Vec scale;bool physicalWeights=false;bool fdJac=false;
    static void check(int value,const char* name){if(value<0)throw std::runtime_error(std::string(name)+":"+std::to_string(value));}
    static int eval(realtype t,N_Vector a,N_Vector b,void* ptr){
        auto& s=*static_cast<Solver*>(ptr);
        try {std::vector<long double> y(s.model.size());
             for(int i=0;i<s.model.size();++i)y[i]=static_cast<long double>(NV_Ith_S(a,i))*s.scale[i];
             auto dy=s.model.rhs(y,t);
             for(int i=0;i<s.model.size();++i)NV_Ith_S(b,i)=static_cast<double>(dy[i]/s.scale[i]);
             return 0;
        }catch(const std::exception& e){s.error=e.what();return 1;}
    }
    static int weights(N_Vector a,N_Vector w,void* ptr){
        auto& s=*static_cast<Solver*>(ptr);const auto& m=s.model;
        auto physical=[&](int i){return NV_Ith_S(a,i)*s.scale[i];};
        double b=physical(m.bottomIndex()),h=physical(m.topIndex())-b;
        if(!(h>0) || (m.l && !(b>0)))return -1;
        for(int i=0;i<m.n;++i){double v=h*m.widths[i],solid=std::clamp(physical(i),0.,v*(1-1e-12)),pore=v-solid;
            NV_Ith_S(w,i)=s.scale[i]/(1e-12*v+s.rtol*solid);
            double T=physical(m.n+i),Q=T+m.p.rho*m.referenceC*pore;
            double den=s.physicalWeights?1e-12*m.p.rho*pore+s.rtol*std::abs(Q):s.rtol*(m.p.rho*pore*1e-5+std::abs(T));
            NV_Ith_S(w,m.n+i)=s.scale[m.n+i]/den;
        }
        for(int i=0;i<m.l;++i){int k=2*m.n+i;double pore=b*m.lw[i],T=physical(k),Q=T+m.p.rho*m.referenceC*pore;
            double den=s.physicalWeights?1e-12*m.p.rho*pore+s.rtol*std::abs(Q):s.rtol*(m.p.rho*pore*1e-5+std::abs(T));
            NV_Ith_S(w,k)=s.scale[k]/den;}
        for(int i=m.topIndex();i<m.size();++i){double den=s.physicalWeights?1e-14+s.rtol*std::abs(physical(i)):s.rtol*(std::abs(physical(i))+(i<m.ledger()?1e-6:m.p.crystal*m.p.initialHeight*1e-4));
            NV_Ith_S(w,i)=s.scale[i]/den;}
        return 0;
    }
    static int jacobian(realtype t,N_Vector state,N_Vector fy,SUNMatrix J,void* ptr,N_Vector,N_Vector,N_Vector){
        auto& s=*static_cast<Solver*>(ptr);const auto& m=s.model;
        try{
            std::vector<Dual> dual(m.size());Vec y(m.size());
            for(int i=0;i<m.size();++i){y[i]=NV_Ith_S(state,i)*s.scale[i];dual[i]=Dual(static_cast<long double>(NV_Ith_S(state,i))*s.scale[i]);}
            std::atomic<int> failures{0};
            #pragma omp parallel for schedule(static)
            for(int j=0;j<m.size();++j){
              try {
                if(j>=m.ledger()){for(int i=0;i<m.size();++i)SM_ELEMENT_D(J,i,j)=0;continue;}
                if(!s.fdJac){auto column=dual;column[j].d=s.scale[j];auto f=m.rhs(column,t);
                    for(int i=0;i<m.size();++i)SM_ELEMENT_D(J,i,j)=f[i].d/s.scale[i];
                }else{
                    double dx=1e-7*std::max(1.,std::abs(NV_Ith_S(state,j))),d=dx*s.scale[j];
                    if(j==m.bottomIndex() && m.l){d=std::min(d,.01*y[j]);dx=d/s.scale[j];}
                    Vec plus=y;plus[j]+=d;auto fp=m.rhs(plus,t);
                    bool central=!(j<m.n && y[j]<=d) && !(j==m.bottomIndex() && y[j]<=d);
                    if(central){Vec minus=y;minus[j]-=d;auto fm=m.rhs(minus,t);
                        for(int i=0;i<m.size();++i)SM_ELEMENT_D(J,i,j)=(fp[i]-fm[i])/(2*dx*s.scale[i]);
                    }else for(int i=0;i<m.size();++i)SM_ELEMENT_D(J,i,j)=(fp[i]/s.scale[i]-NV_Ith_S(fy,i))/dx;
                }
              } catch(const std::exception& e){
                failures.fetch_add(1);
                #pragma omp critical
                {s.error=e.what();}
              }
            }
            if(failures.load())return 1;

            return 0;
        }catch(const std::exception& e){s.error=e.what();return 1;}
    }
    void diagnose(double requested,int flag){
        double now=0,hs=0,hl=0;long ni=0,nf=0,ne=0,nj=0,ns=0;
        CVodeGetCurrentTime(mem,&now);CVodeGetCurrentStep(mem,&hs);CVodeGetLastStep(mem,&hl);
        CVodeGetNonlinSolvStats(mem,&ni,&nf);CVodeGetNumErrTestFails(mem,&ne);CVodeGetNumJacEvals(mem,&nj);CVodeGetNumSteps(mem,&ns);
        N_Vector internal=nullptr;CVodeGetCurrentState(mem,&internal);
        auto ew=N_VClone(state),er=N_VClone(state);CVodeGetErrWeights(mem,ew);CVodeGetEstLocalErrors(mem,er);
        Vec y(N_VGetArrayPointer(internal),N_VGetArrayPointer(internal)+model.size());for(int i=0;i<model.size();++i)y[i]*=scale[i];auto f=model.flux(y,now);
        std::vector<int> rank(model.size());std::iota(rank.begin(),rank.end(),0);
        auto score=[&](int i){return std::abs(NV_Ith_S(ew,i)*NV_Ith_S(er,i));};
        std::sort(rank.begin(),rank.end(),[&](int a,int b){return score(a)>score(b);});
        double cmin=*std::min_element(f.cc.begin(),f.cc.end()),cmax=*std::max_element(f.cc.begin(),f.cc.end());
        double pore=1,solid=1,lv=1;for(int i=0;i<model.n;++i){pore=std::min(pore,f.vol[i]-y[i]);solid=std::min(solid,y[i]);}
        for(int i=0;i<model.l;++i){cmin=std::min(cmin,f.cl[i]);cmax=std::max(cmax,f.cl[i]);lv=std::min(lv,f.vl[i]);}
        diagnostics<<std::setprecision(17)<<requested<<','<<now<<','<<flag<<','<<hs<<','<<hl<<','<<ns<<','<<ni<<','<<nf<<','<<ne<<','<<nj<<','<<rank[0]<<','<<score(rank[0])<<','<<rank[1]<<','<<score(rank[1])<<','<<cmin<<','<<cmax<<','<<pore<<','<<solid<<','<<lv<<','<<f.vt<<','<<f.vb<<'\n';diagnostics.flush();
        if(flag<0){std::ofstream snap(std::string("raw/")+(std::getenv("DIAG_ID")?std::getenv("DIAG_ID"):"audit")+"_failure_state.csv");snap<<std::setprecision(17)<<"index,state,weight,estimated_error,weighted_error\n";for(int i:rank)snap<<i<<','<<y[i]<<','<<NV_Ith_S(ew,i)/scale[i]<<','<<NV_Ith_S(er,i)*scale[i]<<','<<score(i)<<'\n';}
        N_VDestroy(ew);N_VDestroy(er);
    }

public:
    Solver(const Model& m,const Vec& initial,double tolerance=1e-8,double maxstep=30,double start=0):model(m),rtol(tolerance){
        diagnostics.open(std::string("raw/")+(std::getenv("DIAG_ID")?std::getenv("DIAG_ID"):"audit")+"_diagnostics.csv");diagnostics<<"requested,internal_t,flag,current_step,last_step,steps,nonlinear_iters,nonlinear_failures,error_failures,jacobians,worst_index,worst_weighted_error,second_index,second_weighted_error,Cmin,Cmax,pore_min,solid_min,liquid_cell_min,vt,vb\n";
        check(SUNContext_Create(nullptr,&ctx),"context");state=N_VNew_Serial(m.size(),ctx);constraints=N_VClone(state);
        physicalWeights=std::getenv("WEIGHT_MODE") && std::string(std::getenv("WEIGHT_MODE"))=="physical";
        fdJac=std::getenv("JAC_MODE") && std::string(std::getenv("JAC_MODE"))=="fd";
        scale.assign(m.size(),m.p.crystal*m.p.initialHeight);
        for(int i=0;i<m.n;++i){scale[i]=m.p.initialHeight*m.widths[i];scale[m.n+i]=m.p.rho*scale[i]*1e-5;}
        for(int i=0;i<m.l;++i)scale[2*m.n+i]=m.p.rho*m.gap0*m.lw[i]*1e-5;
        scale[m.topIndex()]=scale[m.bottomIndex()]=m.p.initialHeight;
        if(std::getenv("SCALE_NONE"))std::fill(scale.begin(),scale.end(),1.);
        for(int i=0;i<m.size();++i){NV_Ith_S(state,i)=initial[i]/scale[i];}
        N_VConst(0.,constraints);
        for(int i=0;i<m.n;++i)NV_Ith_S(constraints,i)=1;
        NV_Ith_S(constraints,m.topIndex())=NV_Ith_S(constraints,m.bottomIndex())=1;
        mem=CVodeCreate(CV_BDF,ctx);check(CVodeInit(mem,eval,start,state),"init");check(CVodeSetUserData(mem,this),"data");
        check(CVodeWFtolerances(mem,weights),"weights");check(CVodeSetConstraints(mem,constraints),"constraints");
        check(CVodeSetMaxNumSteps(mem,300000),"steps");check(CVodeSetMaxStep(mem,maxstep),"maxstep");
        matrix=SUNDenseMatrix(m.size(),m.size(),ctx);linear=SUNLinSol_Dense(state,matrix,ctx);
        check(CVodeSetLinearSolver(mem,linear,matrix),"linear");
        check(CVodeSetJacFn(mem,jacobian),"physical-scale Jacobian");
        if(std::getenv("FRESH_JAC")){
          if(std::string(std::getenv("FRESH_JAC"))!="nonlinear"){
          if(std::string(std::getenv("FRESH_JAC"))!="setup")check(CVodeSetLSetupFrequency(mem,1),"fresh linear setup");
          check(CVodeSetJacEvalFrequency(mem,1),"fresh Jacobian");
          }
          check(CVodeSetNonlinConvCoef(mem,.01),"tight nonlinear residual");
          check(CVodeSetMaxNonlinIters(mem,8),"nonlinear iterations");
        }
    }
    ~Solver(){if(mem)CVodeFree(&mem);if(linear)SUNLinSolFree(linear);if(matrix)SUNMatDestroy(matrix);if(state)N_VDestroy(state);if(constraints)N_VDestroy(constraints);if(ctx)SUNContext_Free(&ctx);}
    Vec advance(double t){double reached=0;int flag=CVode(mem,t,state,&reached,CV_NORMAL);
        diagnose(t,flag);
        if(flag<0)throw std::runtime_error("CVODE "+std::to_string(flag)+":"+error);
        Vec y(N_VGetArrayPointer(state),N_VGetArrayPointer(state)+model.size());for(int i=0;i<model.size();++i)y[i]*=scale[i];model.validate(y,t);return y;}
    void maxStep(double h){check(CVodeSetMaxStep(mem,h),"set maxstep");}
    long steps()const{long n;CVodeGetNumSteps(mem,&n);return n;}
};
}
