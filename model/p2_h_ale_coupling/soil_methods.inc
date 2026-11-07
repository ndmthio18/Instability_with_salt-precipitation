        struct SurfaceState {
            bool active=false;
            p2h::Vec crust;
            double loose=0, initialSolid=0, vapour=0, activation=-1, pratHeight=0;
            double poreSaltTransfer=0, poreWaterTransfer=0;
        };
        SurfaceState acceptedSurface_;
        std::unique_ptr<p2h::ALE> ale_;
        double heldLiquid_=0,heldSalt_=0,initialSalt_=0,initialWater_=0;
        unsigned long coupledTrials_=0,rejectedTrials_=0;
        bool injectedFailure_=false;
        std::ofstream surfaceLog_,eventLog_,crustFields_,soilFields_;

        double area_() const { return this->gridGeometry().bBoxMax()[0]-this->gridGeometry().bBoxMin()[0]; }
        template<class Solution> std::vector<double> volumes_(const Solution& x) const {
            std::vector<double> v(x.size(),0);
            for(const auto& e:elements(this->gridGeometry().gridView())) {
                auto g=localView(this->gridGeometry()); g.bind(e);
                for(const auto& s:scvs(g)) v[s.dofIndex()]+=s.volume();
            }
            return v;
        }
        std::vector<std::size_t> topDofs_() const {
            std::vector<std::size_t> indices;
            for(const auto& v:vertices(this->gridGeometry().gridView()))
                if(std::abs(v.geometry().center()[1]-this->gridGeometry().bBoxMax()[1])<1e-12)
                    indices.push_back(this->gridGeometry().vertexMapper().index(v));
            return indices;
        }
        template<class Solution> std::pair<double,double> inventory_(const Solution& x) const {
            auto v=volumes_(x); double s=0,w=0;
            for(std::size_t i=0;i<x.size();++i) {
                const double c=massFraction(x[i][moleFractionIdx]);
                const double m=porosity_*sw(x[i][pressureIdx])*rhoW(x[i][moleFractionIdx])*v[i];
                s+=m*c; w+=m*(1-c);
            }
            return {s,w};
        }
        template<class Solution> double topC_(const Solution& x) const {
            double c=0; auto ids=topDofs_();
            for(auto i:ids) c+=massFraction(x[i][moleFractionIdx]);
            return c/ids.size();
        }
        template<class Solution> double topSw_(const Solution& x) const {
            double s=1;
            for(auto i:topDofs_()) s=std::min(s,double(sw(x[i][pressureIdx])));
            return s;
        }
        template<class Solution> void assignMass_(Solution& x,std::size_t i,double salt,double water,double vol) const {
            if(!(salt>=0 && water>0)) throw std::runtime_error("negative donor inventory");
            const double w=salt/(salt+water);
            const double mole=w*molarMassH2O()/(molarMassNaCl()-w*(molarMassNaCl()-molarMassH2O()));
            const double sat=(salt+water)/(porosity_*rhoW(mole)*vol);
            const double swr=getParam<double>("Problem.Swr"),n=getParam<double>("Problem.PcN"),m=1-1/n;
            if(!(sat>swr && sat<1)) throw std::runtime_error("donor saturation outside inversion range");
            const double se=(sat-swr)/(1-swr);
            x[i][pressureIdx]=nonWettingReferencePressure()-std::pow(std::pow(se,-1/m)-1,1/n)/getParam<double>("Problem.PcAlpha");
            x[i][moleFractionIdx]=mole;
        }
        double film_(double c) const {
            p2gfilm::GasFilm f;
            f.thickness=getParam<double>("P2G.FilmThickness",.0083);
            f.relativeHumidity=getParam<double>("P2G.RelativeHumidity",.4);
            return f.fluxFromMassFraction(c);
        }
        template<class Solution> void initializeSurface(const Solution& x) {
            if(getParam<int>("Grid.Cells0")!=1) throw std::runtime_error("P2-H supports one horizontal cell only");
            p2h::Parameters p;
            p.ceq=massFraction(getParam<double>("Transport.SolubilityMoleFraction"));
            p.bottomC=p.ceq; p.crystal=2165;
            p.rho=rhoW(getParam<double>("Transport.SolubilityMoleFraction"));
            ale_=std::make_unique<p2h::ALE>(getParam<int>("P2H.CrustCells",64),2,p);
            ale_->gasFilm=[this](double c){return film_(c);};
            const auto m=inventory_(x); initialSalt_=m.first; initialWater_=m.second;
            heldLiquid_=film_(topC_(x)); heldSalt_=0;
            surfaceLog_.open("coupled_history.csv"); eventLog_.open("coupling_events.csv"); crustFields_.open("crust_fields.csv");
            soilFields_.open("soil_fields.csv");
            surfaceLog_<<std::setprecision(17)<<"time_s,dt_s,active,activation_s,h_m,h_mass_m,h_prat_m,top_solid_kg_m2,internal_solid_kg_m2,epsilon_min,top_soil_mass_fraction,top_soil_sw,evaporation_kg_m2_s,soil_salt_kg,crust_salt_kg,soil_water_kg,crust_water_kg,vapour_kg,salt_residual_kg,water_residual_kg,solid_ledger_residual_kg,pore_salt_transfer_kg,pore_water_transfer_kg,coupled_trials,rejected_trials,interface_liquid_kg_m2,interface_salt_kg_m2\n";
            eventLog_<<std::setprecision(17)<<"time_s,event,value1,value2\n";
            crustFields_<<std::setprecision(17)<<"time_s,z_bottom_m,z_top_m,epsilon,concentration\n";
            soilFields_<<std::setprecision(17)<<"time_s,y_m,mass_fraction,saturation,density_kg_m3\n";
            writeSurface_(x,0,0,0);
        }
        template<class Solution> void precrust_(Solution& x,SurfaceState& next,double time,double dt) {
            const double ceq=ale_->p.ceq,area=area_(); auto vol=volumes_(x);
            for(auto i:topDofs_()) {
                const double c=massFraction(x[i][moleFractionIdx]);
                const double total=porosity_*sw(x[i][pressureIdx])*rhoW(x[i][moleFractionIdx])*vol[i];
                const double s=total*c,w=total-s,p=std::max(0.,s-ceq/(1-ceq)*w);
                if(p>0) { assignMass_(x,i,s-p,w,vol[i]); next.loose+=p/area; }
            }
            next.vapour+=heldLiquid_*dt*area;
            const double h=next.loose/(ale_->p.crystal*(1-ale_->p.birth));
            if(h<getParam<double>("P2H.ActivationHeight",1e-6)) return;
            double available=0;
            for(auto i:topDofs_()) available+=porosity_*sw(x[i][pressureIdx])*rhoW(x[i][moleFractionIdx])*vol[i];
            const double required=ale_->p.rho*ale_->p.birth*h*area;
            if(!(available>required)) { eventLog_<<time<<",activation_blocked_inventory,"<<available<<','<<required<<'\n'; return; }
            auto donor=x; double saltTransfer=0,waterTransfer=0;
            const auto before=inventory_(x);
            for(auto i:topDofs_()) {
                const double c=massFraction(x[i][moleFractionIdx]);
                const double total=porosity_*sw(x[i][pressureIdx])*rhoW(x[i][moleFractionIdx])*vol[i];
                const double fraction=required/available;
                saltTransfer+=total*c*fraction; waterTransfer+=total*(1-c)*fraction;
                assignMass_(donor,i,total*c*(1-fraction),total*(1-c)*(1-fraction),vol[i]);
            }
            const double minSw=getParam<double>("P2H.MinimumSupplySaturation",.99);
            if(topSw_(donor)<minSw) { eventLog_<<time<<",activation_blocked_saturation,"<<topSw_(donor)<<','<<minSw<<'\n'; return; }
            auto crust=ale_->initial(h,saltTransfer/(saltTransfer+waterTransfer));
            auto cm=ale_->inventory(crust); const auto after=inventory_(donor);
            if(std::abs(after.first+cm.salt*area-before.first)>1e-12*initialSalt_ ||
               std::abs(after.second+cm.water*area-before.second)>1e-12*initialWater_ ||
               std::abs(cm.solid-next.loose)>1e-12*initialSalt_/area)
                throw std::runtime_error("activation transfer budget failed");
            next.active=true; next.crust=std::move(crust); next.initialSolid=next.loose; next.loose=0;
            next.poreSaltTransfer=saltTransfer; next.poreWaterTransfer=waterTransfer;
            next.activation=time; next.pratHeight=h; x=donor;
            eventLog_<<time<<",activation_trial,"<<h<<','<<topSw_(x)<<'\n';
        }
        template<class Solution> void checkBudgets_(const Solution& x,const SurfaceState& state) const {
            const auto sm=inventory_(x); double salt=state.loose*area_(),water=0;
            if(state.active) { auto cm=ale_->inventory(state.crust); salt=(cm.salt+cm.solid)*area_(); water=cm.water*area_(); }
            if(std::abs(sm.first+salt-initialSalt_)>1e-8*initialSalt_ ||
               std::abs(sm.second+water+state.vapour-initialWater_)>1e-8*initialWater_)
                throw std::runtime_error("global coupled budget failed");
        }
        template<class Solution,class Solver,class TimeLoop,class GridVariables>
        void advanceCoupled(Solution& x,Solver& solver,TimeLoop& loop,GridVariables& variables) {
            const Solution base=x; const SurfaceState original=acceptedSurface_;
            double dt=loop.timeStepSize();
            for(int retry=0;retry<9;++retry) {
                SurfaceState next=original; x=base; variables.update(x); loop.setTimeStepSize(dt);
                double j=heldLiquid_;
                try {
                    if(!original.active) {
                        heldLiquid_=film_(topC_(base)); heldSalt_=0; j=heldLiquid_;
                        solver.solve(x,loop); dt=loop.timeStepSize();
                        precrust_(x,next,loop.time()+dt,dt);
                    } else {
                        double c=topC_(base),previousC=0,previousR=0; bool solved=false;
                        for(int iteration=0;iteration<16;++iteration) {
                            ++coupledTrials_;
                            auto model=*ale_; model.p.bottomC=c;
                            p2h::BDF integrator(model,original.crust,std::min(dt,60.));
                            auto trial=integrator.advance(dt); const int k=model.inv();
                            heldLiquid_=(trial[k+3]-original.crust[k+3])/dt;
                            heldSalt_=(trial[k+2]-original.crust[k+2])/dt;
                            if(heldLiquid_<0) throw std::runtime_error("reverse interface supply outside demonstration");
                            x=base; variables.update(x); loop.setTimeStepSize(dt); solver.solve(x,loop);
                            if(loop.timeStepSize()<dt*(1-1e-12)) { dt=loop.timeStepSize(); throw std::runtime_error("soil reduced trial dt"); }
                            const double residual=topC_(x)-c;
                            if(std::abs(residual)<1e-10) {
                                if(topSw_(x)<getParam<double>("P2H.MinimumSupplySaturation",.99))
                                    throw std::runtime_error("saturated crust supply screen failed");
                                j=(trial[k+4]-original.crust[k+4])/dt;
                                next.crust=std::move(trial); next.vapour+=j*dt*area_();
                                next.pratHeight+=p2g::velocity(next.pratHeight,j,ale_->p.ceq)*dt;
                                solved=true; break;
                            }
                            double newC=c+1e-7;
                            if(iteration>0) {
                                const double slope=(residual-previousR)/(c-previousC);
                                if(!std::isfinite(slope) || std::abs(slope)<1e-8) throw std::runtime_error("interface root slope invalid");
                                newC=c-std::clamp(residual/slope,-1e-3,1e-3);
                            }
                            if(!(newC>0 && newC<.95)) throw std::runtime_error("interface root outside composition range");
                            previousC=c; previousR=residual; c=newC;
                        }
                        if(!solved) throw std::runtime_error("coupled interface root not converged");
                    }
                    checkBudgets_(x,next);
                    if(next.active && (!getParam<bool>("P2H.RejectAfterActive",false) || original.active) &&
                       !injectedFailure_ && getParam<bool>("P2H.InjectRejectedTrial",false)) {
                        injectedFailure_=true; throw std::runtime_error("injected rollback test");
                    }
                    acceptedSurface_=std::move(next);
                    writeSurface_(x,loop.time()+dt,dt,j);
                    return;
                } catch(const std::exception& e) {
                    ++rejectedTrials_;
                    std::string reason=e.what();
                    std::replace(reason.begin(),reason.end(),',',';');
                    std::replace(reason.begin(),reason.end(),'\n',' ');
                    eventLog_<<loop.time()<<",rejected_"<<reason<<','<<dt<<','<<retry<<'\n'; eventLog_.flush();
                    // Both states remain at the old accepted time. No compensation.
                    x=base; variables.update(x); dt*=.5;
                    if(retry==8 || dt<1e-6) throw;
                }
            }
        }
        bool activationOnlyComplete() const { return acceptedSurface_.active && getParam<bool>("P2H.ActivationOnly",false); }
        template<class Solution> void writeSurface_(const Solution& x,double time,double dt,double j) {
            const auto& st=acceptedSurface_; const auto sm=inventory_(x); const double a=area_();
            double cs=0,cw=0,top=st.loose,inside=0,eps=ale_->p.birth,h=0,solid=st.loose,defect=0;
            if(st.active) {
                auto m=ale_->inventory(st.crust); const int k=ale_->inv();
                cs=m.salt*a; cw=m.water*a; solid=m.solid; top=st.initialSolid+st.crust[k]; inside=st.crust[k+1];
                h=st.crust[2*ale_->n]; eps=*std::min_element(st.crust.begin(),st.crust.begin()+ale_->n);
                defect=(solid-top-inside)*a;
                if(time==st.activation || std::fmod(time,100)<1e-8 || activationOnlyComplete())
                    for(int i=0;i<ale_->n;++i) crustFields_<<time<<','<<ale_->faces[i]*h<<','<<ale_->faces[i+1]*h<<','<<st.crust[i]<<','<<st.crust[ale_->n+i]<<'\n';
            }
            surfaceLog_<<time<<','<<dt<<','<<st.active<<','<<st.activation<<','<<h<<','<<solid/(ale_->p.crystal*(1-ale_->p.birth))<<','<<st.pratHeight<<','
                <<top<<','<<inside<<','<<eps<<','<<topC_(x)<<','<<topSw_(x)<<','<<j<<','<<sm.first<<','<<cs<<','<<sm.second<<','<<cw<<','<<st.vapour<<','
                <<sm.first+cs+solid*a-initialSalt_<<','<<sm.second+cw+st.vapour-initialWater_<<','<<defect<<','<<st.poreSaltTransfer<<','<<st.poreWaterTransfer<<','<<coupledTrials_<<','<<rejectedTrials_<<','
                <<(st.active?st.crust[ale_->inv()+3]:0)<<','<<(st.active?st.crust[ale_->inv()+2]:0)<<'\n';
            if(time==0 || time==st.activation || std::fmod(time,100)<1e-8 || activationOnlyComplete())
                for(const auto& vertex:vertices(this->gridGeometry().gridView())) {
                    const auto i=this->gridGeometry().vertexMapper().index(vertex);
                    soilFields_<<time<<','<<vertex.geometry().center()[1]<<','<<massFraction(x[i][moleFractionIdx])<<','<<sw(x[i][pressureIdx])<<','<<rhoW(x[i][moleFractionIdx])<<'\n';
                }
            surfaceLog_.flush(); eventLog_.flush(); crustFields_.flush(); soilFields_.flush();
        }
