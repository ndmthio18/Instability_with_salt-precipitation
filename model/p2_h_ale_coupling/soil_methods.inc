        #include "grid_diagnostics.inc"
        bool topFace_(const SubControlVolumeFace& face) const {
            // This one-column rectangular mesh has an upward top and lateral
            // faces which can lie inside the legacy coordinate tolerance.
            return onUpperBoundary_(face.ipGlobal()) && face.unitOuterNormal()[1] > 0.5;
        }
        void auditBoundaryFaces_() const {
            std::ofstream log("boundary_faces.csv");
            log<<std::setprecision(17)<<"id,x_m,y_m,nx,ny,area_m2,legacy_top,selected_top\n";
            double selectedArea=0,sideArea=0; unsigned id=0;
            for(const auto& e:elements(this->gridGeometry().gridView())) {
                auto g=localView(this->gridGeometry()); g.bind(e);
                for(const auto& face:scvfs(g)) if(face.boundary()) {
                    const auto x=face.ipGlobal(); const auto n=face.unitOuterNormal();
                    const bool selected=topFace_(face);
                    log<<id++<<','<<x[0]<<','<<x[1]<<','<<n[0]<<','<<n[1]<<','<<face.area()<<','
                       <<onUpperBoundary_(x)<<','<<selected<<'\n';
                    if(selected) { selectedArea+=face.area(); if(n[1]<.5) sideArea+=face.area(); }
                }
            }
            if(std::abs(selectedArea-area_())>1e-12*area_() || sideArea!=0)
                throw std::runtime_error("top boundary area audit failed");
        }
        struct SurfaceState {
            bool active=false;
            p2h::Vec crust;
            double loose=0, initialSolid=0, vapour=0, activation=-1, pratHeight=0;
            double poreSaltTransfer=0, poreWaterTransfer=0, precipFlux=0;
            p2h::PrecursorInventory precursor;
            double internalPositive=0, internalDissolution=0;
            double bottomLiquidFlux=0, bottomSaltFlux=0;
            double cumulativeBottomLiquid=0, cumulativeBottomSalt=0;
        };
        struct BottomFlux {
            double liquid=0, salt=0;
            double bottomLiquidBefore=0, bottomLiquidAfter=0;
            double bottomSaltBefore=0, bottomSaltAfter=0;
            double cutLiquidDown=0, cutSaltDown=0;
        };
        SurfaceState acceptedSurface_;
        std::unique_ptr<p2h::ALE> ale_;
        double heldLiquid_=0,heldSalt_=0,initialSalt_=0,initialWater_=0;
        double evaporationFlux_=0,precipitationFlux_=0,fillLiquidFlux_=0,fillSaltFlux_=0,fillDensity_=0;
        double nextFieldOutputTime_=0;
        unsigned long coupledTrials_=0,rejectedTrials_=0;
        bool injectedFailure_=false;
        std::ofstream surfaceLog_,eventLog_,crustFields_,soilFields_,donorLog_,allocationLog_,boundaryLog_;
        unsigned long allocationAttempt_=0, candidateAttempt_=0;
        bool controlledDonorComplete() const { return getParam<bool>("P2H.ControlledDonorTest",false); }

        template<class Solution> std::vector<p2h::Donor> donors_(const Solution& x) const {
            const auto vol=volumes_(x); std::vector<p2h::Donor> donors(x.size());
            const double top=this->gridGeometry().bBoxMax()[1];
            for(const auto& vertex:vertices(this->gridGeometry().gridView())) {
                const auto i=this->gridGeometry().vertexMapper().index(vertex);
                donors[i]={i,top-vertex.geometry().center()[1],0,vol[i],porosity_,
                    sw(x[i][pressureIdx]),rhoW(x[i][moleFractionIdx]),massFraction(x[i][moleFractionIdx])};
            }
            // Sum each element SCV once into its global BOX DOF; rectangular
            // tensor meshes may be graded. The opposite corner follows from
            // SCV center and DOF position, not an assumed uniform spacing.
            for(const auto& e:elements(this->gridGeometry().gridView())) {
                auto g=localView(this->gridGeometry()); g.bind(e);
                for(const auto& s:scvs(g)) {
                    const auto i=s.dofIndex();
                    const double low=std::min(s.dofPosition()[1],2*s.center()[1]-s.dofPosition()[1]);
                    donors[i].deepest=std::max(donors[i].deepest,top-low);
                }
            }
            return donors;
        }
        template<class Solution> bool adaptiveDonors_(Solution& x,SurfaceState& next,double time,double h) {
            const auto before=inventory_(x); const auto original=x;
            const double area=area_(),required=ale_->p.rho*ale_->p.birth*h*area;
            const auto plan=p2h::allocate(donors_(x),required,getParam<double>("P2H.MaxDonorDepth",.005),
                                        getParam<double>("P2H.MinimumSupplySaturation",.99));
            const auto attempt=++allocationAttempt_; auto candidate=x;
            double minSw=1,saltResidual=0,waterResidual=0;
            unsigned count=0;
            if(plan.sufficient) {
                for(const auto& d:plan.donors) if(d.extracted>0) {
                    ++count;
                    const double mass=d.porosity*d.saturation*d.density*d.volume;
                    assignMass_(candidate,d.id,(mass-d.extracted)*d.concentration,
                                (mass-d.extracted)*(1-d.concentration),d.volume);
                    minSw=std::min(minSw,double(sw(candidate[d.id][pressureIdx])));
                }
                if(minSw<getParam<double>("P2H.MinimumSupplySaturation",.99))
                    throw std::runtime_error("reconstructed donor saturation below retained minimum");
                const double mix=plan.salt/(plan.salt+plan.water);
                auto crust=ale_->initial(h,mix); const auto cm=ale_->inventory(crust);
                const auto after=inventory_(candidate);
                saltResidual=after.first+cm.salt*area-before.first;
                waterResidual=after.second+cm.water*area-before.second;
                if(std::abs(saltResidual)>1e-12*initialSalt_ || std::abs(waterResidual)>1e-12*initialWater_ ||
                   std::abs(cm.solid-next.loose)>1e-12*initialSalt_/area)
                    throw std::runtime_error("adaptive activation transfer budget failed");
                next.active=true; next.crust=std::move(crust); next.initialSolid=next.loose; next.loose=0;
                next.poreSaltTransfer=plan.salt; next.poreWaterTransfer=plan.water;
                next.activation=time; next.pratHeight=h; x=candidate; candidateAttempt_=attempt;
            }
            allocationLog_<<attempt<<','<<time<<','<<h<<','<<required<<','<<plan.available<<','<<plan.sufficient<<','
                <<count<<','<<plan.usedDepth<<','<<(plan.sufficient?plan.salt:0)<<','<<(plan.sufficient?plan.water:0)<<','
                <<saltResidual<<','<<waterResidual<<','<<initialSalt_<<','<<initialWater_<<'\n';
            // First refusal and every successful trial keep complete profiles;
            // accepted_attempt events distinguish rejected trials from acceptance.
            if(attempt==1 || plan.sufficient || controlledDonorComplete()) {
                const auto all=donors_(original);
                for(const auto& d:all) {
                    double extracted=0,available=0;
                    for(const auto& p:plan.donors) if(p.id==d.id) {
                        available=p.available; if(plan.sufficient) extracted=p.extracted;
                    }
                    donorLog_<<attempt<<','<<time<<','<<d.id<<','<<d.depth<<','<<d.deepest<<','<<d.volume<<','
                        <<available<<','<<extracted<<','<<d.saturation<<','<<sw(candidate[d.id][pressureIdx])<<','
                        <<d.concentration<<','<<massFraction(candidate[d.id][moleFractionIdx])<<','
                        <<original[d.id][pressureIdx]<<','<<candidate[d.id][pressureIdx]<<','<<d.density<<','
                        <<rhoW(candidate[d.id][moleFractionIdx])<<'\n';
                }
            }
            eventLog_<<time<<(plan.sufficient?",activation_trial,":",activation_blocked_donor_depth,")
                <<attempt<<','<<(plan.sufficient?plan.usedDepth:required-plan.available)<<'\n';
            allocationLog_.flush(); donorLog_.flush(); eventLog_.flush();
            return plan.sufficient;
        }

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
        template<class Solution,class GridVariables>
        BottomFlux reconstructBottomFlux_(const Solution& before,const Solution& after,double dt,GridVariables& gv) const {
            BottomFlux result;
            if(!bottomIsDirichlet()) return result;
            if(!(dt>0)) throw std::runtime_error("invalid dt for bottom-flux reconstruction");
            const double bottom=this->gridGeometry().bBoxMin()[1];
            const double tol=1e-12*std::max(1.,this->gridGeometry().bBoxMax()[1]-bottom);
            std::vector<double> y(after.size()),levels;
            for(const auto& vertex:vertices(this->gridGeometry().gridView())) {
                const auto i=this->gridGeometry().vertexMapper().index(vertex);
                y[i]=vertex.geometry().center()[1]; levels.push_back(y[i]);
            }
            std::sort(levels.begin(),levels.end());
            levels.erase(std::unique(levels.begin(),levels.end(),[&](double a,double b){return std::abs(a-b)<tol;}),levels.end());
            if(levels.size()<2) throw std::runtime_error("bottom flux needs two vertical DOF levels");
            const double cut=.5*(levels[0]+levels[1]);
            const auto mb=dofMasses_(before),ma=dofMasses_(after);
            for(std::size_t i=0;i<after.size();++i) if(std::abs(y[i]-bottom)<tol) {
                result.bottomSaltBefore+=mb.first[i]; result.bottomSaltAfter+=ma.first[i];
                result.bottomLiquidBefore+=mb.first[i]+mb.second[i];
                result.bottomLiquidAfter+=ma.first[i]+ma.second[i];
            }
            using Residual=RichardsSaltModelLocalResidual<TypeTag>;
            gv.update(after);
            for(const auto& e:elements(this->gridGeometry().gridView())) {
                auto g=localView(this->gridGeometry()); g.bind(e);
                auto vv=localView(gv.curGridVolVars()); vv.bind(e,g,after);
                auto cache=localView(gv.gridFluxVarsCache()); cache.bind(e,g,vv);
                for(const auto& face:scvfs(g)) if(!face.boundary()) {
                    const auto zi=y[g.scv(face.insideScvIdx()).dofIndex()];
                    const auto zo=y[g.scv(face.outsideScvIdx()).dofIndex()];
                    const int sign=int(zi>cut)-int(zo>cut);
                    if(!sign) continue;
                    const auto f=Residual::transportFlux(*this,e,g,vv,face,cache);
                    result.cutLiquidDown+=sign*f.totalMass;
                    result.cutSaltDown+=sign*f.saltTotal;
                }
            }
            // Backward-Euler balance over the bottom BOX control volumes:
            // dM_bottom/dt = B_in + F_cut,down.
            result.liquid=(result.bottomLiquidAfter-result.bottomLiquidBefore)/dt-result.cutLiquidDown;
            result.salt=(result.bottomSaltAfter-result.bottomSaltBefore)/dt-result.cutSaltDown;
            return result;
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
        template<class Solution> double topRho_(const Solution& x) const {
            double rho=0; const auto ids=topDofs_();
            for(auto i:ids) rho+=rhoW(x[i][moleFractionIdx]);
            return rho/ids.size();
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
            if(getParam<bool>("P2H.VerificationZeroEvaporation",false)) return 0;
            p2gfilm::GasFilm f;
            f.thickness=getParam<double>("P2G.FilmThickness",.0083);
            f.relativeHumidity=getParam<double>("P2G.RelativeHumidity",.4);
            return f.fluxFromMassFraction(c);
        }
        template<class Solution> void initializeSurface(const Solution& x) {
            if(getParam<int>("Grid.Cells0")!=1) throw std::runtime_error("P2-H supports one horizontal cell only");
            auditBoundaryFaces_();
            precipitationMode_();
            std::ofstream config("activation_config.csv");
            config<<std::setprecision(17)<<"threshold_m,birth_porosity,crust_cells,reaction_coefficient_m_s,av_multiplier,initial_height_m,initial_solid_kg_m2\n"
                  <<getParam<double>("P2H.ActivationHeight",3.5e-6)<<",0.43,"
                  <<getParam<int>("P2H.CrustCells",64)<<','
                  <<getParam<double>("P2H.ReactionCoefficient",2.3e-3)<<','
                  <<getParam<double>("P2H.AvMultiplier",1.0)<<",0,0\n";
            p2h::Parameters p;
            p.ceq=massFraction(getParam<double>("Transport.SolubilityMoleFraction"));
            p.bottomC=p.ceq; p.crystal=2165;
            p.rho=rhoW(getParam<double>("Transport.SolubilityMoleFraction"));
            p.reaction=getParam<double>("P2H.ReactionCoefficient",p.reaction);
            p.avMultiplier=getParam<double>("P2H.AvMultiplier",1.0);
            p.splitReactionLedgers=true;
            if(!(std::isfinite(p.reaction) && p.reaction>0 &&
                 std::isfinite(p.avMultiplier) && p.avMultiplier>0))
                throw std::runtime_error("invalid P2-H reaction parameter");
            ale_=std::make_unique<p2h::ALE>(getParam<int>("P2H.CrustCells",64),2,p);
            ale_->gasFilm=[this](double c){return film_(c);};
            const auto m=inventory_(x); initialSalt_=m.first; initialWater_=m.second;
            heldLiquid_=film_(topC_(x)); heldSalt_=0;
            evaporationFlux_=heldLiquid_; precipitationFlux_=fillLiquidFlux_=fillSaltFlux_=0;
            surfaceLog_.open("coupled_history.csv"); eventLog_.open("coupling_events.csv"); crustFields_.open("crust_fields.csv");
            soilFields_.open("soil_fields.csv");
            donorLog_.open("donor_profiles.csv"); allocationLog_.open("donor_attempts.csv");
            boundaryLog_.open("open_boundary_budget.csv");
            nextFieldOutputTime_=getParam<double>("P2H.FieldOutputInterval",100.0);
            if(!(std::isfinite(nextFieldOutputTime_) && nextFieldOutputTime_>0))
                throw std::runtime_error("invalid P2-H field output interval");
            donorLog_<<std::setprecision(17)<<"attempt,time_s,dof,depth_m,box_deepest_m,volume_m3,available_kg,extracted_kg,sw_before,sw_after,c_before,c_after,p_before_pa,p_after_pa,rho_before,rho_after\n";
            allocationLog_<<std::setprecision(17)<<"attempt,time_s,h_m,required_kg,available_kg,sufficient,donor_count,used_depth_m,salt_transfer_kg,water_transfer_kg,salt_transfer_residual_kg,water_transfer_residual_kg,initial_salt_kg,initial_water_kg\n";
            boundaryLog_<<std::setprecision(17)<<"time_s,dt_s,bottom_liquid_in_kg_s,bottom_salt_in_kg_s,bottom_water_in_kg_s,bottom_liquid_in_kg_m2_s,bottom_salt_in_kg_m2_s,bottom_water_in_kg_m2_s,cumulative_bottom_liquid_kg,cumulative_bottom_salt_kg,cumulative_bottom_water_kg,bottom_liquid_before_kg,bottom_liquid_after_kg,bottom_salt_before_kg,bottom_salt_after_kg,cut_liquid_down_kg_s,cut_salt_down_kg_s,soil_water_kg,crust_water_kg,vapour_kg,soil_salt_kg,crust_dissolved_salt_kg,solid_salt_kg,water_residual_kg,salt_residual_kg\n";
            const double x0=getParam<double>("Problem.SaltMoleFraction");
            double xmin=std::numeric_limits<double>::infinity(),xmax=-xmin;
            for(const auto& value:x) { xmin=std::min(xmin,double(value[moleFractionIdx])); xmax=std::max(xmax,double(value[moleFractionIdx])); }
            const double yBottom=this->gridGeometry().bBoxMin()[1];
            const double reservoirPressure=nonWettingReferencePressure()-(getParam<double>("Problem.WaterLevel")+yBottom)*rhoW(x0)*gravitationalAcceleration();
            std::ofstream reservoir("reservoir_config.csv");
            reservoir<<std::setprecision(17)<<"bottom_mode,domain_bottom_y_m,domain_top_y_m,water_level_parameter_m,water_table_y_m,reservoir_pressure_pa,configured_salt_mole_fraction,reservoir_salt_mole_fraction,initial_min_salt_mole_fraction,initial_max_salt_mole_fraction\n"
                     <<(bottomIsDirichlet()?"reservoir":"no_flux")<<','<<yBottom<<','<<this->gridGeometry().bBoxMax()[1]<<','
                     <<getParam<double>("Problem.WaterLevel")<<','<<-getParam<double>("Problem.WaterLevel")<<','<<reservoirPressure<<','
                     <<x0<<','<<x0<<','<<xmin<<','<<xmax<<'\n';
            if(std::max(std::abs(xmin-x0),std::abs(xmax-x0))>1e-14)
                throw std::runtime_error("accepted initial salt field differs from Problem.SaltMoleFraction");
            const auto mode=getParam<std::string>("P2H.DonorMode","top_only");
            if(mode!="top_only" && mode!="adaptive_multilayer" && mode!="progressive_top_boundary_filling")
                throw std::runtime_error("unknown donor mode");
            if(mode=="progressive_top_boundary_filling" && precipitationMode_()!="implicit_surface")
                throw std::runtime_error("progressive filling requires implicit surface precipitation");
            surfaceLog_<<std::setprecision(17)<<"time_s,dt_s,active,activation_s,h_m,h_mass_m,h_prat_m,top_solid_kg_m2,internal_solid_kg_m2,internal_positive_kg_m2,internal_dissolution_kg_m2,internal_quadrature_residual_kg_m2,top_new_solid_kg_m2,internal_fraction_net,internal_fraction_positive,positive_reaction_kg_m2_s,dissolution_kg_m2_s,top_growth_velocity_m_s,crust_top_mass_fraction,epsilon_min,top_soil_mass_fraction,top_soil_sw,evaporation_kg_m2_s,soil_salt_kg,crust_salt_kg,soil_water_kg,crust_water_kg,vapour_kg,salt_residual_kg,water_residual_kg,solid_ledger_residual_kg,pore_salt_transfer_kg,pore_water_transfer_kg,coupled_trials,rejected_trials,interface_liquid_kg_m2,interface_salt_kg_m2,reaction_coefficient_m_s,av_multiplier,precipitation_flux_kg_m2_s,fill_liquid_flux_kg_m2_s,fill_salt_flux_kg_m2_s,precursor_pore_salt_kg_m2,precursor_pore_water_kg_m2,precursor_required_volume_m3_m2,precursor_filled_volume_m3_m2,precursor_unfilled_volume_m3_m2,precursor_concentration\n";
            eventLog_<<std::setprecision(17)<<"time_s,event,value1,value2\n";
            crustFields_<<std::setprecision(17)<<"time_s,z_bottom_m,z_top_m,epsilon,concentration\n";
            soilFields_<<std::setprecision(17)<<"time_s,y_m,mass_fraction,saturation,density_kg_m3\n";
            initializeGridDiagnostics_(x);
            writeSurface_(x,0,0,0);
            if(controlledDonorComplete()) {
                // Separate prescribed-inventory allocation harness only. Its
                // synthetic solid inventory never enters a coupled simulation.
                auto candidate=x; SurfaceState test;
                const double h=getParam<double>("P2H.ActivationHeight",3.5e-6);
                test.loose=h*p.crystal*(1-p.birth);
                const bool ok=adaptiveDonors_(candidate,test,0,h);
                eventLog_<<0<<(ok?",controlled_allocation_pass,":",controlled_allocation_refused,")<<h<<",0\n";
            }
        }
        template<class Solution> void precrust_(Solution& x,SurfaceState& next,double time,double dt) {
            const double ceq=ale_->p.ceq,area=area_(); auto vol=volumes_(x);
            const auto donorMode=getParam<std::string>("P2H.DonorMode","top_only");
            if(precipitationMode_()=="projection") for(auto i:topDofs_()) {
                const double c=massFraction(x[i][moleFractionIdx]);
                const double total=porosity_*sw(x[i][pressureIdx])*rhoW(x[i][moleFractionIdx])*vol[i];
                const double s=total*c,w=total-s,p=std::max(0.,s-ceq/(1-ceq)*w);
                if(p>0) { assignMass_(x,i,s-p,w,vol[i]); next.loose+=p/area; }
            }
            if(precipitationMode_()=="implicit_surface") {
                next.loose+=precipitationFlux_*dt; next.precipFlux=precipitationFlux_;
                if(donorMode=="progressive_top_boundary_filling") {
                    const auto increment=p2h::progressiveIncrement(precipitationFlux_,dt,ale_->p.birth,
                        ale_->p.crystal,topC_(x),topRho_(x));
                    const double fluxTolerance=1e-12*std::max(1e-12,fillLiquidFlux_);
                    if(std::abs((increment.poreSalt+increment.poreWater)/dt-fillLiquidFlux_)>fluxTolerance ||
                       std::abs(increment.poreSalt/dt-fillSaltFlux_)>fluxTolerance)
                        throw std::runtime_error("progressive accepted interface flux mismatch");
                    p2h::acceptProgressiveIncrement(next.precursor,increment);
                }
            }
            next.vapour+=evaporationFlux_*dt*area;
            const double h=next.loose/(ale_->p.crystal*(1-ale_->p.birth));
            if(h<getParam<double>("P2H.ActivationHeight",3.5e-6)) return;
            if(getParam<bool>("P2H.PrecrustOnly",false)) return;
            if(donorMode=="progressive_top_boundary_filling") {
                const double volumeTolerance=1e-12*std::max(h,1e-12);
                const double volumeResidual=p2h::precursorVolumeResidual(next.precursor);
                if(std::abs(volumeResidual)>volumeTolerance) {
                    eventLog_<<time<<",activation_blocked_unfilled_precursor,"<<volumeResidual<<','<<next.precursor.requiredPoreVolume<<'\n';
                    return;
                }
                if(topSw_(x)<getParam<double>("P2H.MinimumSupplySaturation",.99)) {
                    eventLog_<<time<<",activation_blocked_top_state,"<<topSw_(x)<<','<<getParam<double>("P2H.MinimumSupplySaturation",.99)<<'\n';
                    return;
                }
                const double c=p2h::precursorConcentration(next.precursor);
                auto crust=ale_->initial(h,c); const auto cm=ale_->inventory(crust);
                const double saltResidual=cm.salt-next.precursor.poreSalt;
                const double waterResidual=cm.water-next.precursor.poreWater;
                const double solidResidual=cm.solid-next.loose;
                const double scale=std::max(initialSalt_/area,1e-30);
                if(std::abs(saltResidual)>1e-12*scale || std::abs(waterResidual)>1e-12*scale ||
                   std::abs(solidResidual)>1e-12*scale || std::abs(cm.pore-next.precursor.filledPoreVolume)>volumeTolerance)
                    throw std::runtime_error("progressive activation transfer budget failed");
                next.active=true; next.crust=std::move(crust); next.initialSolid=next.loose; next.loose=0;
                next.poreSaltTransfer=next.precursor.poreSalt*area;
                next.poreWaterTransfer=next.precursor.poreWater*area;
                next.activation=time; next.pratHeight=h; candidateAttempt_=++allocationAttempt_;
                eventLog_<<time<<",activation_trial_progressive,"<<h<<','<<c<<'\n';
                return;
            }
            if(donorMode=="adaptive_multilayer") {
                adaptiveDonors_(x,next,time,h); return;
            }
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
            const auto sm=inventory_(x); double salt=(state.loose+state.precursor.poreSalt)*area_();
            double water=state.precursor.poreWater*area_();
            if(state.active) { auto cm=ale_->inventory(state.crust); salt=(cm.salt+cm.solid)*area_(); water=cm.water*area_(); }
            const double bottomWater=state.cumulativeBottomLiquid-state.cumulativeBottomSalt;
            if(std::abs(sm.first+salt-state.cumulativeBottomSalt-initialSalt_)>1e-8*initialSalt_ ||
               std::abs(sm.second+water+state.vapour-bottomWater-initialWater_)>1e-8*initialWater_)
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
                        evaporationFlux_=film_(topC_(base)); heldLiquid_=evaporationFlux_; heldSalt_=0;
                        precipitationFlux_=fillLiquidFlux_=fillSaltFlux_=0; j=evaporationFlux_;
                        if(precipitationMode_()=="implicit_surface")
                            solvePrecipitation_(x,base,original,solver,loop,variables,dt,j);
                        else { solver.solve(x,loop); dt=loop.timeStepSize(); }
                        captureTransport_(x,variables);
                        precrust_(x,next,loop.time()+dt,dt);
                    } else {
                        precipitationFlux_=fillLiquidFlux_=fillSaltFlux_=0;
                        double c=topC_(base),previousC=0,previousR=0; bool solved=false;
                        for(int iteration=0;iteration<16;++iteration) {
                            ++coupledTrials_;
                            auto model=*ale_; model.p.bottomC=c;
                            p2h::BDF integrator(model,original.crust,std::min(dt,60.));
                            integrator.rejectedState=[&](const p2h::Vec& y,double t){ writeAleDiagnostic_(model,y,loop.time()+t,"bound_rejected",c); };
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
                                next.internalPositive=trial[model.positiveInv()];
                                next.internalDissolution=trial[model.dissolutionInv()];
                                diagnosticBottomC_=c;
                                next.crust=std::move(trial); next.vapour+=j*dt*area_();
                                evaporationFlux_=j;
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
                    const auto bottom=reconstructBottomFlux_(base,x,dt,variables);
                    next.bottomLiquidFlux=bottom.liquid; next.bottomSaltFlux=bottom.salt;
                    next.cumulativeBottomLiquid+=bottom.liquid*dt;
                    next.cumulativeBottomSalt+=bottom.salt*dt;
                    checkBudgets_(x,next);
                    if(next.active && (!getParam<bool>("P2H.RejectAfterActive",false) || original.active) &&
                       !injectedFailure_ && getParam<bool>("P2H.InjectRejectedTrial",false)) {
                        injectedFailure_=true; throw std::runtime_error("injected rollback test");
                    }
                    acceptedSurface_=std::move(next);
                    if(!original.active && acceptedSurface_.active)
                        eventLog_<<loop.time()+dt<<",accepted_activation_attempt,"<<candidateAttempt_<<','<<acceptedSurface_.activation<<'\n';
                    const double topEvaporation=original.active?j:evaporationFlux_;
                    writeGridDiagnostics_(base,x,original,acceptedSurface_,loop.time()+dt,dt,topEvaporation);
                    writeSurface_(x,loop.time()+dt,dt,topEvaporation,&bottom);
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
        template<class Solution> void writeSurface_(const Solution& x,double time,double dt,double j,const BottomFlux* bottom=nullptr) {
            const auto& st=acceptedSurface_; const auto sm=inventory_(x); const double a=area_();
            const double outputInterval=getParam<double>("P2H.FieldOutputInterval",100.0);
            const double outputTolerance=1e-10*std::max(1.0,time);
            const bool periodicOutput=time==0 || time+outputTolerance>=nextFieldOutputTime_;
            const bool fieldOutput=periodicOutput || time==st.activation || activationOnlyComplete();
            if(periodicOutput && time>0)
                while(nextFieldOutputTime_<=time+outputTolerance) nextFieldOutputTime_+=outputInterval;
            double cs=st.precursor.poreSalt*a,cw=st.precursor.poreWater*a,top=st.loose,inside=0,eps=ale_->p.birth,h=0,solid=st.loose,defect=0;
            double topNew=0,netFraction=0,positiveFraction=0,positiveRate=0,dissolutionRate=0,growth=0,crustTopC=0;
            if(st.active) {
                auto m=ale_->inventory(st.crust); const int k=ale_->inv();
                cs=m.salt*a; cw=m.water*a; solid=m.solid; top=st.initialSolid+st.crust[k]; inside=st.crust[k+1];
                h=st.crust[2*ale_->n]; eps=*std::min_element(st.crust.begin(),st.crust.begin()+ale_->n);
                topNew=st.crust[k];
                const auto flux=ale_->flux(st.crust);
                growth=flux.speed; crustTopC=flux.topC;
                for(int i=0;i<ale_->n;++i) {
                    positiveRate+=std::max(0.,flux.reaction[i])*flux.volume[i];
                    dissolutionRate+=std::max(0.,-flux.reaction[i])*flux.volume[i];
                }
                if(inside+topNew>0) netFraction=inside/(inside+topNew);
                if(st.internalPositive+topNew>0) positiveFraction=st.internalPositive/(st.internalPositive+topNew);
                defect=(solid-top-inside)*a;
                if(fieldOutput)
                    for(int i=0;i<ale_->n;++i) crustFields_<<time<<','<<ale_->faces[i]*h<<','<<ale_->faces[i+1]*h<<','<<st.crust[i]<<','<<st.crust[ale_->n+i]<<'\n';
            }
            const double bottomWater=st.cumulativeBottomLiquid-st.cumulativeBottomSalt;
            const double saltResidual=sm.first+cs+solid*a-st.cumulativeBottomSalt-initialSalt_;
            const double waterResidual=sm.second+cw+st.vapour-bottomWater-initialWater_;
            surfaceLog_<<time<<','<<dt<<','<<st.active<<','<<st.activation<<','<<h<<','<<solid/(ale_->p.crystal*(1-ale_->p.birth))<<','<<st.pratHeight<<','
                <<top<<','<<inside<<','<<st.internalPositive<<','<<st.internalDissolution<<','
                <<st.internalPositive-st.internalDissolution-inside<<','<<topNew<<','<<netFraction<<','<<positiveFraction<<','
                <<positiveRate<<','<<dissolutionRate<<','<<growth<<','<<crustTopC<<','<<eps<<','<<topC_(x)<<','<<topSw_(x)<<','<<j<<','<<sm.first<<','<<cs<<','<<sm.second<<','<<cw<<','<<st.vapour<<','
                <<saltResidual<<','<<waterResidual<<','<<defect<<','<<st.poreSaltTransfer<<','<<st.poreWaterTransfer<<','<<coupledTrials_<<','<<rejectedTrials_<<','
                <<(st.active?st.crust[ale_->inv()+3]:0)<<','<<(st.active?st.crust[ale_->inv()+2]:0)<<','
                <<ale_->p.reaction<<','<<ale_->p.avMultiplier<<','<<precipitationFlux_<<','<<fillLiquidFlux_<<','<<fillSaltFlux_<<','
                <<st.precursor.poreSalt<<','<<st.precursor.poreWater<<','<<st.precursor.requiredPoreVolume<<','<<st.precursor.filledPoreVolume<<','
                <<p2h::precursorVolumeResidual(st.precursor)<<','<<((st.precursor.poreSalt+st.precursor.poreWater)>0?p2h::precursorConcentration(st.precursor):0)<<'\n';
            const BottomFlux zero; const auto& bf=bottom?*bottom:zero;
            boundaryLog_<<time<<','<<dt<<','<<st.bottomLiquidFlux<<','<<st.bottomSaltFlux<<','<<st.bottomLiquidFlux-st.bottomSaltFlux<<','
                <<st.bottomLiquidFlux/a<<','<<st.bottomSaltFlux/a<<','<<(st.bottomLiquidFlux-st.bottomSaltFlux)/a<<','
                <<st.cumulativeBottomLiquid<<','<<st.cumulativeBottomSalt<<','<<bottomWater<<','
                <<bf.bottomLiquidBefore<<','<<bf.bottomLiquidAfter<<','<<bf.bottomSaltBefore<<','<<bf.bottomSaltAfter<<','
                <<bf.cutLiquidDown<<','<<bf.cutSaltDown<<','<<sm.second<<','<<cw<<','<<st.vapour<<','<<sm.first<<','<<cs<<','<<solid*a<<','
                <<waterResidual<<','<<saltResidual<<'\n';
            if(fieldOutput)
                for(const auto& vertex:vertices(this->gridGeometry().gridView())) {
                    const auto i=this->gridGeometry().vertexMapper().index(vertex);
                    soilFields_<<time<<','<<vertex.geometry().center()[1]<<','<<massFraction(x[i][moleFractionIdx])<<','<<sw(x[i][pressureIdx])<<','<<rhoW(x[i][moleFractionIdx])<<'\n';
                }
            surfaceLog_.flush(); eventLog_.flush(); crustFields_.flush(); soilFields_.flush(); boundaryLog_.flush();
        }
