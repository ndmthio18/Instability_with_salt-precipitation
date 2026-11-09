        struct SurfaceState {
            pratboundary::State growth;
            p2h::PrecursorInventory precursor;
            double loose=0,precipFlux=0,vapour=0,activation=-1,bottomSalt=0,bottomLiquid=0;
        };
        pratboundary::Parameters params_;
        SurfaceState acceptedSurface_;
        double heldLiquid_=0,heldSalt_=0,initialSalt_=0,initialWater_=0;
        double evaporationFlux_=0,precipitationFlux_=0,fillLiquidFlux_=0,fillSaltFlux_=0,fillDensity_=0;
        double nextOutput_=0;
        pratboundary::GrowthLaw growthLaw_=pratboundary::GrowthLaw::M1;
        double substrateGrainSize_=25e-6;
        bool injected_=false;
        std::ofstream history_,events_,profiles_,rejected_;
        // GENERATED_HELPERS
        double film_(double c) const {
            if(getParam<bool>("PratBoundary.ZeroEvaporation",false)) return 0;
            p2gfilm::GasFilm f;
            f.thickness=getParam<double>("PratBoundary.FilmThickness",.0083);
            f.relativeHumidity=getParam<double>("PratBoundary.RelativeHumidity",.4);
            return f.fluxFromMassFraction(c);
        }
        template<class Solution> void initializeSurface(const Solution& x) {
            if(getParam<int>("Grid.Cells0")!=1) throw std::runtime_error("analytical column requires one horizontal cell");
            params_.equilibriumMassFraction=massFraction(getParam<double>("Transport.SolubilityMoleFraction"));
            params_.solutionDensity=rhoW(getParam<double>("Transport.SolubilityMoleFraction"));
            params_.reactionCoefficient=getParam<double>("PratBoundary.ReactionCoefficient",.0023);
            params_.minimumInterfaceSaturation=getParam<double>("PratBoundary.MinimumSupplySaturation",.99);
            growthLaw_=pratboundary::parseGrowthLaw(getParam<std::string>("PratBoundary.GrowthLaw","M1"));
            substrateGrainSize_=getParam<double>("PratBoundary.SubstrateGrainSize",25e-6);
            if(!(substrateGrainSize_>0 && std::isfinite(substrateGrainSize_))) throw std::runtime_error("invalid substrate grain size");
            params_.check(); auditBoundaryFaces_();
            const auto m=inventory_(x); initialSalt_=m.first; initialWater_=m.second;
            history_.open("history.csv"); events_.open("events.csv"); profiles_.open("soil_profiles.csv"); rejected_.open("rejected_trials.csv");
            history_<<std::setprecision(17)<<"time_s,dt_s,active,activation_s,growth_law,h_m,solid_kg_m2,pore_salt_kg_m2,pore_water_kg_m2,j_kg_m2_s,v_m_s,da,da_star,da_top,specific_area_m_inv,efficiency_m1,efficiency_m2,efficiency_m3,selected_efficiency,m3_relative_error,m3_valid,pe_u,pe_move,c_top,c_soil,c_soil_relative_to_ceq,sw_soil,total_out_kg_m2_s,salt_out_kg_m2_s,solid_rate_kg_m2_s,pore_salt_rate_kg_m2_s,pore_water_rate_kg_m2_s,soil_salt_kg,soil_water_kg,vapour_kg,bottom_salt_kg,bottom_water_kg,salt_residual_rel,water_residual_rel,precip_flux,fill_liquid_flux\n";
            events_<<std::setprecision(17)<<"time_s,event,value1,value2\n";
            rejected_<<std::setprecision(17)<<"time_s,dt_s,retry,growth_law,reason,c_soil,sw_soil,total_out_kg_m2_s,salt_out_kg_m2_s,da,da_star,da_top,selected_efficiency\n";
            profiles_<<std::setprecision(17)<<"time_s,z_m,c,sw,pressure_pa\n";
            write_(x,0,0,0,0,0,0,0,{});
        }
        template<class Solution> void budget_(const Solution& x,const SurfaceState& s) const {
            for(const auto& value:x) {
                const double c=massFraction(value[moleFractionIdx]),sat=sw(value[pressureIdx]);
                if(!(c>=0 && c<1 && sat>=0 && sat<=1 && rhoW(value[moleFractionIdx])>0))
                    throw std::runtime_error("invalid local soil inventory");
            }
            const auto m=inventory_(x);
            const double solid=s.growth.active?s.growth.solidSalt:s.loose;
            const double salt=s.growth.active?s.growth.dissolvedSalt:s.precursor.poreSalt;
            const double water=s.growth.active?s.growth.water:s.precursor.poreWater;
            if(!(solid>=0 && salt>=0 && water>=0)) throw std::runtime_error("negative surface inventory");
            if(std::abs(m.first+(solid+salt)*area_()-s.bottomSalt-initialSalt_)>1e-8*initialSalt_ ||
               std::abs(m.second+water*area_()+s.vapour-(s.bottomLiquid-s.bottomSalt)-initialWater_)>1e-8*initialWater_)
                throw std::runtime_error("combined budget failure");
        }
        template<class Solution,class Solver,class TimeLoop,class Variables>
        void advanceCoupled(Solution& x,Solver& solver,TimeLoop& loop,Variables& variables) {
            const auto base=x; const auto original=acceptedSurface_;
            const double savedLiquid=heldLiquid_,savedSalt=heldSalt_;
            double dt=loop.timeStepSize();
            for(int retry=0;retry<9;++retry) {
                auto next=original; x=base; variables.update(x); loop.setTimeStepSize(dt);
                double j=0,v=0,pu=0,pm=0,ct=0;
                pratboundary::GrowthDiagnostics growth;
                try {
                    precipitationFlux_=fillLiquidFlux_=fillSaltFlux_=0;
                    if(!original.growth.active) {
                        j=film_(topC_(base));
                        solvePrecipitation_(x,base,original,solver,loop,variables,dt,j);
                        auto inc=p2h::progressiveIncrement(precipitationFlux_,dt,params_.porosity,
                            params_.crystalDensity,params_.equilibriumMassFraction,fillDensity_);
                        inc.poreSalt=fillSaltFlux_*dt;
                        inc.poreWater=(fillLiquidFlux_-fillSaltFlux_)*dt;
                        // Store exactly the actual accepted interface transfer.
                        next.loose+=inc.solid; next.precipFlux=precipitationFlux_;
                        p2h::acceptProgressiveIncrement(next.precursor,inc);
                        const double h=next.loose/(params_.crystalDensity*(1-params_.porosity));
                        if(h>=getParam<double>("PratBoundary.ActivationHeight",3.5e-6) && !getParam<bool>("PratBoundary.PrecrustOnly",false)) {
                            pratboundary::validateSupply(params_,topC_(x),topSw_(x));
                            next.growth=pratboundary::activate(params_,next.loose,next.precursor.poreSalt,next.precursor.poreWater);
                            next.growth.time=loop.time()+dt; next.activation=loop.time()+dt;
                            const double rs=(next.growth.totalSalt()-next.loose-next.precursor.poreSalt)*area_()/initialSalt_;
                            const double rw=(next.growth.water-next.precursor.poreWater)*area_()/initialWater_;
                            if(std::abs(rs)>1e-12 || std::abs(rw)>1e-12) throw std::runtime_error("activation budget failed");
                            events_<<next.activation<<",activation_trial,"<<rs<<','<<rw<<'\n';
                        }
                    } else {
                        auto trial=pratboundary::evaluate(original.growth,params_,dt,[this](double c){return film_(c);},growthLaw_,substrateGrainSize_);
                        growth=trial.growth;
                        heldLiquid_=trial.totalFlux; heldSalt_=trial.saltFlux;
                        solver.solve(x,loop);
                        if(loop.timeStepSize()<dt*(1-1e-12)) { dt=loop.timeStepSize(); throw std::runtime_error("soil reduced active trial dt"); }
                        if(topSw_(x)<params_.minimumInterfaceSaturation ||
                           std::abs(topC_(x)/params_.equilibriumMassFraction-1)>params_.interfaceRelativeTolerance)
                            events_<<loop.time()+dt<<",invalid_trial_interface,"<<topC_(x)<<','<<topSw_(x)<<'\n';
                        pratboundary::validateSupply(params_,topC_(x),topSw_(x));
                        p2f::commit(next.growth,trial);
                        j=trial.evaporation; v=trial.velocity; pu=trial.solutionPe; pm=trial.movingPe; ct=trial.topMassFraction;
                    }
                    next.vapour+=j*dt*area_();
                    const auto b=reconstructBottomFlux_(base,x,dt,variables);
                    next.bottomSalt+=b.salt*dt; next.bottomLiquid+=b.liquid*dt;
                    budget_(x,next);
                    if(next.growth.active && !injected_ && getParam<bool>("PratBoundary.InjectRejectedTrial",false)) {
                        injected_=true; throw std::runtime_error("injected rollback");
                    }
                    acceptedSurface_=next;
                    if(original.loose==0 && next.loose>0) events_<<loop.time()+dt<<",accepted_first_precipitation,"<<next.loose<<",0\n";
                    if(!original.growth.active && next.growth.active) events_<<next.activation<<",accepted_activation,"<<next.growth.height<<','<<topC_(x)<<'\n';
                    if(next.growth.active && growth.damkohler==0)
                        growth=pratboundary::growthDiagnostics(params_,next.growth.height,substrateGrainSize_,growthLaw_);
                    write_(x,loop.time()+dt,dt,j,v,pu,pm,ct,growth); return;
                } catch(const std::exception& e) {
                    double rejectedC=std::numeric_limits<double>::quiet_NaN(),rejectedSw=rejectedC;
                    try { rejectedC=topC_(x); rejectedSw=topSw_(x); } catch(...) {}
                    rejected_<<loop.time()+dt<<','<<dt<<','<<retry<<','<<pratboundary::growthLawName(growthLaw_)<<','<<e.what()<<','
                        <<rejectedC<<','<<rejectedSw<<','<<heldLiquid_<<','<<heldSalt_<<','<<growth.damkohler<<','
                        <<growth.damkohlerStar<<','<<growth.damkohlerTop<<','<<growth.selectedEfficiency<<'\n'; rejected_.flush();
                    events_<<loop.time()<<",rejected_"<<e.what()<<','<<dt<<','<<retry<<'\n'; events_.flush();
                    x=base; variables.update(x); heldLiquid_=savedLiquid; heldSalt_=savedSalt;
                    precipitationFlux_=fillLiquidFlux_=fillSaltFlux_=0;
                    dt*=.5;
                    if(retry==8 || dt<1e-6) {
                        std::ofstream last("last_accepted_soil.csv");
                        last<<std::setprecision(17)<<"time_s,z_m,c,sw,pressure_pa\n";
                        for(const auto& vertex:vertices(this->gridGeometry().gridView())) {
                            const auto i=this->gridGeometry().vertexMapper().index(vertex);
                            last<<loop.time()<<','<<vertex.geometry().center()[1]<<','<<massFraction(x[i][moleFractionIdx])<<','<<sw(x[i][pressureIdx])<<','<<x[i][pressureIdx]<<'\n';
                        }
                        throw;
                    }
                }
            }
        }
        template<class Solution> void write_(const Solution& x,double t,double dt,double j,double v,double pu,double pm,double ct,
                                             const pratboundary::GrowthDiagnostics& growth) {
            const auto& s=acceptedSurface_; const auto m=inventory_(x); const auto a=area_();
            const bool active=s.growth.active;
            const double solid=active?s.growth.solidSalt:s.loose;
            const double salt=active?s.growth.dissolvedSalt:s.precursor.poreSalt;
            const double water=active?s.growth.water:s.precursor.poreWater;
            const double h=active?s.growth.height:solid/(params_.crystalDensity*(1-params_.porosity));
            const double solidRate=params_.crystalDensity*(1-params_.porosity)*v;
            const double poreSaltRate=params_.rho()*params_.porosity*params_.equilibriumMassFraction*v;
            const double poreWaterRate=params_.rho()*params_.porosity*(1-params_.equilibriumMassFraction)*v;
            history_<<t<<','<<dt<<','<<active<<','<<s.activation<<','<<pratboundary::growthLawName(growthLaw_)<<','<<h<<','<<solid<<','<<salt<<','<<water<<','<<j<<','<<v<<','
                <<growth.damkohler<<','<<growth.damkohlerStar<<','<<growth.damkohlerTop<<','<<growth.specificArea<<','
                <<growth.efficiencyM1<<','<<growth.efficiencyM2<<','<<growth.efficiencyM3<<','<<growth.selectedEfficiency<<','
                <<growth.m3RelativeError<<','<<growth.m3Valid<<','<<pu<<','<<pm<<','<<ct<<','<<topC_(x)<<','
                <<topC_(x)/params_.equilibriumMassFraction-1<<','<<topSw_(x)<<','<<heldLiquid_<<','<<heldSalt_<<','
                <<solidRate<<','<<poreSaltRate<<','<<poreWaterRate<<','<<m.first<<','<<m.second<<','<<s.vapour<<','<<s.bottomSalt<<','<<s.bottomLiquid-s.bottomSalt<<','
                <<(m.first+(solid+salt)*a-s.bottomSalt-initialSalt_)/initialSalt_<<','
                <<(m.second+water*a+s.vapour-(s.bottomLiquid-s.bottomSalt)-initialWater_)/initialWater_<<','<<precipitationFlux_<<','<<fillLiquidFlux_<<'\n';
            if(t==0 || t>=nextOutput_ || t==s.activation) {
                for(const auto& vertex:vertices(this->gridGeometry().gridView())) {
                    const auto i=this->gridGeometry().vertexMapper().index(vertex);
                    profiles_<<t<<','<<vertex.geometry().center()[1]<<','<<massFraction(x[i][moleFractionIdx])<<','<<sw(x[i][pressureIdx])<<','<<x[i][pressureIdx]<<'\n';
                }
                nextOutput_=t+getParam<double>("PratBoundary.FieldOutputInterval",100.);
            }
            history_.flush(); events_.flush(); profiles_.flush();
        }
