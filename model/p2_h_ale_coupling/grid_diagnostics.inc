// Included inside the isolated Problem class. Observers never modify x.
std::string precipitationMode_() const {
    const auto mode=getParam<std::string>("P2H.PrecipitationMode","projection");
    if(mode!="projection" && mode!="implicit_surface") throw std::runtime_error("unknown precipitation mode");
    return mode;
}
template<class Solution,class State,class Solver,class TimeLoop>
void solvePrecipitation_(Solution& x,const Solution& base,const State& original,Solver& solver,TimeLoop& loop,GridVariables& variables,double& dt,double evaporation) {
    double flux=original.precipFlux,previousFlux=0,previousResidual=0;
    for(int iteration=0;iteration<24;++iteration) {
        heldSalt_=flux; heldLiquid_=evaporation+flux;
        x=base; variables.update(x); loop.setTimeStepSize(dt); solver.solve(x,loop);
        if(loop.timeStepSize()<dt*(1-1e-12)) {
            dt=loop.timeStepSize(); throw std::runtime_error("soil reduced precipitation trial dt");
        }
        const double residual=topC_(x)-ale_->p.ceq;
        if((flux==0 && residual<=0) || (flux>0 && std::abs(residual)<5e-13)) return;
        double next=flux+std::max(1e-8,evaporation*ale_->p.ceq/(1-ale_->p.ceq));
        if(iteration>0 && flux!=previousFlux) {
            const double slope=(residual-previousResidual)/(flux-previousFlux);
            if(!std::isfinite(slope) || slope>=0) throw std::runtime_error("invalid surface precipitation root slope");
            next=std::max(0.,flux-residual/slope);
        }
        previousFlux=flux; previousResidual=residual; flux=next;
    }
    throw std::runtime_error("surface precipitation complementarity not converged");
}
std::ofstream gridBudget_,gridProfile_,aleDiagnostic_;
std::vector<double> cutDepths_{.0001,.0005,.001,.002,.005,.01};
std::vector<std::vector<double>> cutWeights_;
std::vector<double> transportSalt_,transportWater_,cutAdv_,cutDiff_,cutLiquid_;
std::vector<double> cutBoundarySalt_;
double diagnosticBottomC_=0;
bool gridDiagnostics_() const { return getParam<bool>("P2H.GridDiagnostics",false); }
template<class Solution> void initializeGridDiagnostics_(const Solution& x) {
    if(!gridDiagnostics_()) return;
    const auto vol=volumes_(x); const double top=this->gridGeometry().bBoxMax()[1];
    cutWeights_.assign(cutDepths_.size(),std::vector<double>(x.size(),0));
    for(const auto& e:elements(this->gridGeometry().gridView())) {
        auto g=localView(this->gridGeometry()); g.bind(e);
        for(const auto& scv:scvs(g)) {
            const double a=scv.dofPosition()[1],b=2*scv.center()[1]-a;
            const double lo=std::min(a,b),hi=std::max(a,b);
            for(std::size_t k=0;k<cutDepths_.size();++k)
                cutWeights_[k][scv.dofIndex()]+=scv.volume()/vol[scv.dofIndex()]*
                    std::clamp((hi-std::max(lo,top-cutDepths_[k]))/(hi-lo),0.,1.);
        }
    }
    gridBudget_.open("grid_budgets.csv"); gridProfile_.open("grid_profiles.csv");
    aleDiagnostic_.open("ale_diagnostics.csv");
    gridBudget_<<std::setprecision(17)<<"time_s,dt_s,active,depth_m,salt_before_kg,salt_transport_kg,salt_after_kg,water_before_kg,water_transport_kg,water_after_kg,adv_in_kg_s,diff_in_kg_s,liquid_in_kg_s,projection_salt_kg,surface_solid_kg,soil_salt_kg,soil_water_kg,evap_kg_m2_s,supersaturation_excess_kg,boundary_reconstruction_salt_kg_s,boundary_precip_kg_m2_s\n";
    gridProfile_<<std::setprecision(17)<<"time_s,dof,depth_m,box_volume_m3,pressure_pa,sw,c,rho\n";
    aleDiagnostic_<<std::setprecision(17)<<"time_s,status,i,z_bottom_m,z_top_m,epsilon,c,reaction_kg_m3_s,deff_m2_s,av_m_inv,h_m,speed_m_s,bottom_c,bottom_liquid_kg_m2_s,bottom_salt_kg_m2_s,bottom_adv_kg_m2_s,bottom_diff_kg_m2_s,evap_kg_m2_s,pore_salt_kg_m2,pore_water_kg_m2,solid_kg_m2\n";
    writeGridProfile_(x,0);
}
template<class Solution> std::pair<std::vector<double>,std::vector<double>> dofMasses_(const Solution& x) const {
    const auto v=volumes_(x); std::vector<double> s(x.size()),w(x.size());
    for(std::size_t i=0;i<x.size();++i) {
        const double c=massFraction(x[i][moleFractionIdx]);
        const double m=porosity_*sw(x[i][pressureIdx])*rhoW(x[i][moleFractionIdx])*v[i];
        s[i]=m*c; w[i]=m*(1-c);
    }
    return {s,w};
}
template<class Solution> void captureTransport_(const Solution& x,GridVariables& gv) {
    if(!gridDiagnostics_()) return;
    const auto masses=dofMasses_(x); transportSalt_=masses.first; transportWater_=masses.second;
    cutAdv_.assign(cutDepths_.size(),0); cutDiff_=cutAdv_; cutLiquid_=cutAdv_; cutBoundarySalt_=cutAdv_;
    using Residual=RichardsSaltModelLocalResidual<TypeTag>;
    // bind from x explicitly: never use cached variables from a previous trial.
    gv.update(x);
    for(const auto& e:elements(this->gridGeometry().gridView())) {
        auto g=localView(this->gridGeometry()); g.bind(e);
        auto vv=localView(gv.curGridVolVars()); vv.bind(e,g,x);
        auto cache=localView(gv.gridFluxVarsCache()); cache.bind(e,g,vv);
        for(const auto& face:scvfs(g)) {
            if(face.boundary()) {
                if(topFace_(face)) {
                    const auto i=g.scv(face.insideScvIdx()).dofIndex();
                    for(std::size_t k=0;k<cutDepths_.size();++k) {
                        cutLiquid_[k]+=(1-cutWeights_[k][i])*heldLiquid_*face.area();
                        cutBoundarySalt_[k]+=(1-cutWeights_[k][i])*heldSalt_*face.area();
                    }
                }
                continue;
            }
            const auto i=g.scv(face.insideScvIdx()).dofIndex(),j=g.scv(face.outsideScvIdx()).dofIndex();
            const auto f=Residual::transportFlux(*this,e,g,vv,face,cache);
            for(std::size_t k=0;k<cutDepths_.size();++k) {
                const double weight=cutWeights_[k][j]-cutWeights_[k][i];
                cutAdv_[k]+=weight*f.saltAdvective; cutDiff_[k]+=weight*f.saltDiffusive;
                cutLiquid_[k]+=weight*f.totalMass;
            }
        }
    }
}
template<class Solution> void writeGridProfile_(const Solution& x,double time) {
    const auto v=volumes_(x); const double top=this->gridGeometry().bBoxMax()[1];
    for(const auto& vertex:vertices(this->gridGeometry().gridView())) {
        const auto i=this->gridGeometry().vertexMapper().index(vertex);
        gridProfile_<<time<<','<<i<<','<<top-vertex.geometry().center()[1]<<','<<v[i]<<','
            <<x[i][pressureIdx]<<','<<sw(x[i][pressureIdx])<<','<<massFraction(x[i][moleFractionIdx])<<','<<rhoW(x[i][moleFractionIdx])<<'\n';
    }
}
void writeAleDiagnostic_(const p2h::ALE& model,const p2h::Vec& y,double time,const char* status,double bottomC) {
    if(!gridDiagnostics_()) return;
    const auto f=model.flux(y); const auto m=model.inventory(y); const double h=y[2*model.n];
    for(int i=0;i<model.n;++i)
        aleDiagnostic_<<time<<','<<status<<','<<i<<','<<model.faces[i]*h<<','<<model.faces[i+1]*h<<','<<y[i]<<','<<y[model.n+i]<<','<<f.reaction[i]<<','<<std::sqrt(y[i])*model.p.diffusion<<','<<36*(1-y[i])/model.p.grain<<','<<h<<','<<f.speed<<','<<bottomC<<','<<f.liquid[0]<<','<<f.salt[0]<<','<<bottomC*f.liquid[0]<<','<<f.salt[0]-bottomC*f.liquid[0]<<','<<f.evaporation<<','<<m.salt<<','<<m.water<<','<<m.solid<<'\n';
    aleDiagnostic_.flush();
}
template<class Solution,class State> void writeGridDiagnostics_(const Solution& base,const Solution& x,const State& old,const State& now,double time,double dt,double j) {
    if(!gridDiagnostics_()) return;
    const auto before=dofMasses_(base),after=dofMasses_(x); const auto all=inventory_(x);
    if(!old.active) for(std::size_t k=0;k<cutDepths_.size();++k) {
        double bs=0,bw=0,ts=0,tw=0,as=0,aw=0,excess=0;
        for(std::size_t i=0;i<x.size();++i) {
            const double w=cutWeights_[k][i];
            bs+=w*before.first[i]; bw+=w*before.second[i]; ts+=w*transportSalt_[i]; tw+=w*transportWater_[i];
            as+=w*after.first[i]; aw+=w*after.second[i];
            excess+=w*std::max(0.,after.first[i]-ale_->p.ceq/(1-ale_->p.ceq)*after.second[i]);
        }
        gridBudget_<<time<<','<<dt<<','<<now.active<<','<<cutDepths_[k]<<','<<bs<<','<<ts<<','<<as<<','<<bw<<','<<tw<<','<<aw<<','<<cutAdv_[k]<<','<<cutDiff_[k]<<','<<cutLiquid_[k]<<','<<ts-as<<','<<(now.active?now.initialSolid:now.loose)*area_()<<','<<all.first<<','<<all.second<<','<<j<<','<<excess<<','<<cutBoundarySalt_[k]<<','<<heldSalt_<<'\n';
    }
    // Full near-activation fields on every accepted step; elsewhere every 100 s.
    if(time<=10 || (now.active && time-now.activation<=10) || std::fmod(time,100)<1e-8)
        writeGridProfile_(x,time);
    if(now.active) {
        auto model=*ale_; model.p.bottomC=old.active?diagnosticBottomC_:topC_(x);
        writeAleDiagnostic_(model,now.crust,time,"accepted",model.p.bottomC);
    }
    gridBudget_.flush(); gridProfile_.flush();
}
