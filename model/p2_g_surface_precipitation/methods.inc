        // P2-G: called once after a converged transport step, never by Newton.
        template<class Solution> void surfaceStep(Solution& x, Scalar time, Scalar dt)
        {
            const auto& gg=this->gridGeometry();
            const Scalar area=gg.bBoxMax()[0]-gg.bBoxMin()[0];
            const Scalar top=gg.bBoxMax()[1];
            const Scalar ceq=massFraction(getParam<Scalar>("Transport.SolubilityMoleFraction"));
            std::vector<Scalar> vol(x.size(),0), yy(x.size(),0);
            for(const auto& e:elements(gg.gridView())) {
                auto g=localView(gg); g.bind(e);
                for(const auto& s:scvs(g)) vol[s.dofIndex()]+=s.volume();
            }
            for(const auto& v:vertices(gg.gridView())) yy[gg.vertexMapper().index(v)]=v.geometry().center()[1];
            Scalar beforeS=0,beforeW=0,afterS=0,afterW=0,removed=0,topW=0,minSw=1;
            int count=0;
            auto trial=x;
            for(std::size_t i=0;i<x.size();++i) {
                const Scalar mole=x[i][moleFractionIdx], w=massFraction(mole);
                const Scalar total=porosity_*sw(x[i][pressureIdx])*rhoW(mole)*vol[i];
                const Scalar salt=total*w, water=total-salt;
                beforeS+=salt; beforeW+=water;
                Scalar p=0;
                if(std::abs(yy[i]-top)<1e-12 && dt>0 && getParam<bool>("P2G.EnablePrecipitation",true))
                    p=p2g::precipitate(salt,water,ceq);
                if(p>0) {
                    const Scalar ws=(salt-p)/(salt-p+water);
                    const Scalar xs=ws*molarMassH2O()/(molarMassNaCl()-ws*(molarMassNaCl()-molarMassH2O()));
                    const Scalar sat=(salt-p+water)/(porosity_*rhoW(xs)*vol[i]);
                    const Scalar swr=getParam<Scalar>("Problem.Swr"), n=getParam<Scalar>("Problem.PcN"), m=1-1/n;
                    if(!(sat>swr && sat<1)) throw std::runtime_error("surface projection saturation outside invertible range");
                    const Scalar se=(sat-swr)/(1-swr);
                    trial[i][pressureIdx]=nonWettingReferencePressure()-std::pow(std::pow(se,-1/m)-1,1/n)/getParam<Scalar>("Problem.PcAlpha");
                    trial[i][moleFractionIdx]=xs;
                }
                const Scalar newTotal=porosity_*sw(trial[i][pressureIdx])*rhoW(trial[i][moleFractionIdx])*vol[i];
                afterS+=newTotal*massFraction(trial[i][moleFractionIdx]);
                afterW+=newTotal*(1-massFraction(trial[i][moleFractionIdx]));
                removed+=p;
                if(std::abs(yy[i]-top)<1e-12) {
                    topW+=massFraction(trial[i][moleFractionIdx]); ++count;
                    minSw=std::min(minSw,sw(trial[i][pressureIdx]));
                }
            }
            if(dt==0) {
                initialSalt_=beforeS; initialWater_=beforeW;
                surfaceLog_.open("surface_precipitation.csv");
                surfaceLog_<<std::setprecision(17)<<"time_s,dt_s,evaporation_kg_m2_s,precipitated_step_kg_m2,solid_kg_m2,h_mass_m,h_prat_m,prat_demand_kg_m2,prat_zero_reaction_demand_kg_m2,top_mass_fraction,min_surface_sw,soil_salt_kg,soil_water_kg,evaporated_water_kg,salt_residual_kg,water_residual_kg,projection_salt_residual_kg,projection_water_residual_kg,prat_active\n";
            }
            Scalar demand=0,zeroDemand=0;
            // Prediction is independent of available salt, and never debits soil.
            if(pratActive_) {
                const Scalar dh=p2g::velocity(predictedHeight_,filmFlux_,ceq)*dt;
                zeroDemand=2165*.57*p2g::velocity(predictedHeight_,filmFlux_,ceq,false)*dt;
                demand=2165*.57*dh; predictedHeight_+=dh;
            }
            solid_+=removed/area;
            const Scalar hmass=solid_/(2165*.57);
            if(!pratActive_ && hmass>=getParam<Scalar>("P2G.ContinuumHeight",1e-6)) {
                pratActive_=true; predictedHeight_=hmass;
            }
            evaporated_+=filmFlux_*area*dt;
            surfaceLog_<<time<<','<<dt<<','<<filmFlux_<<','<<removed/area<<','<<solid_<<','<<hmass<<','<<predictedHeight_<<','<<demand<<','<<zeroDemand<<','<<topW/count<<','<<minSw<<','<<afterS<<','<<afterW<<','<<evaporated_<<','<<afterS+solid_*area-initialSalt_<<','<<afterW+evaporated_-initialWater_<<','<<afterS+removed-beforeS<<','<<afterW-beforeW<<','<<pratActive_<<'\n';
            surfaceLog_.flush();
            x=trial;
        }
        template<class Solution> void prepareSurfaceFlux(const Solution& x) {
            Scalar w=0; int count=0;
            for(const auto& v:vertices(this->gridGeometry().gridView())) {
                if(std::abs(v.geometry().center()[1]-this->gridGeometry().bBoxMax()[1])<1e-12) {
                    w+=massFraction(x[this->gridGeometry().vertexMapper().index(v)][moleFractionIdx]); ++count;
                }
            }
            p2gfilm::GasFilm film;
            film.thickness=getParam<Scalar>("P2G.FilmThickness",.0083);
            film.relativeHumidity=getParam<Scalar>("P2G.RelativeHumidity",.4);
            filmFlux_=film.fluxFromMassFraction(w/count);
            if(filmFlux_<0) throw std::runtime_error("condensation outside P2-G scope");
        }
        Scalar filmFlux_=0,solid_=0,predictedHeight_=0,evaporated_=0,initialSalt_=0,initialWater_=0;
        bool pratActive_=false;
        std::ofstream surfaceLog_;
