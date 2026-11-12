#include "bdf.hh"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct Arguments {
    std::string id, suite, mode, output;
    double da=0, u=0, power=0, rtol=0, maxStepDiv=0, endTau=0, outputsPerTau=0;
    int cells=0;
};

Arguments arguments(int argc,char** argv) {
    if(argc!=13) throw std::runtime_error(
        "usage: driver id suite Da u cells power rtol max_step_div end_tau "
        "outputs_per_tau mode output_dir");
    Arguments a;
    a.id=argv[1]; a.suite=argv[2]; a.da=std::stod(argv[3]); a.u=std::stod(argv[4]);
    a.cells=std::stoi(argv[5]); a.power=std::stod(argv[6]); a.rtol=std::stod(argv[7]);
    a.maxStepDiv=std::stod(argv[8]); a.endTau=std::stod(argv[9]);
    a.outputsPerTau=std::stod(argv[10]); a.mode=argv[11]; a.output=argv[12];
    if(!(a.da>0 && a.u>=0 && a.u<=1 && a.cells>=4 && a.power>=1 &&
         a.rtol>0 && a.maxStepDiv>0 && a.endTau>0 && a.outputsPerTau>0))
        throw std::runtime_error("invalid driver arguments");
    if(a.mode!="neumann" && a.mode!="equilibrium_inflow" && a.mode!="dirichlet")
        throw std::runtime_error("unknown bottom boundary mode");
    return a;
}
bool snapshot(double theta) {
    const double requested[]={0,1,2,5,10,15,20};
    for(double value:requested)
        if(std::abs(theta-value)<1e-8) return true;
    return false;
}
}

int main(int argc,char** argv) {
    try {
        const Arguments a=arguments(argc,argv);
        std::filesystem::create_directories(a.output);
        p2h::Parameters p;
        p.ceq=0.264940869768; p.rho=1202.168577; p.crystal=2165.0; p.birth=0.43;
        p.diffusion=1.3e-9; p.reaction=1.0e-3; p.avMultiplier=1.0; p.beta=2.0;
        p.grain=25.0e-6; p.bottomC=p.ceq; p.splitReactionLedgers=true;
        p.allowCompleteDissolution=true;
        p.prescribedDiffusionUsesEquilibriumInflow=(a.mode=="equilibrium_inflow");
        const double deff0=std::sqrt(p.birth)*p.diffusion;
        const double av0=p.avMultiplier*36*(1-p.birth)/p.grain;
        const double mathcalD0=p.birth*deff0/(1-p.ceq);
        const double ell0=std::sqrt(mathcalD0/(av0*p.reaction));
        p.initialHeight=ell0*std::sqrt(a.da);
        const double diffusionTime=p.initialHeight*p.initialHeight/deff0;
        const double endTime=a.endTau*diffusionTime;
        const double outputDt=diffusionTime/a.outputsPerTau;
        const double maxStep=diffusionTime/a.maxStepDiv;
        const double back=a.u*p.ceq*1.78e-5;

        p2h::ALE model(a.cells,a.power,p);
        model.evaporationHistory=[](double){return 1.78e-5;};
        if(a.mode!="dirichlet")
            model.bottomBackDiffusion=[=](double){return back;};
        auto initial=model.initial(p.initialHeight,p.ceq), y=initial;
        const auto start=model.inventory(initial);
        p2h::BDF solver(model,initial,maxStep,a.rtol,true);

        const std::string prefix=a.output+"/"+a.id;
        std::ofstream history(prefix+"_history.csv"), profiles(prefix+"_profiles.csv");
        history<<std::setprecision(17)
          <<"case_id,suite,boundary_mode,time_s,time_over_tau,height_m,"
          <<"target_da,da_birth_properties,da_volume_average,da_star_volume_average,"
          <<"velocity_m_s,top_concentration,bottom_cell_concentration,bottom_face_concentration,"
          <<"evaporation_kg_m2_s,b_soil_kg_m2_s,bottom_liquid_flux_kg_m2_s,"
          <<"bottom_advective_salt_kg_m2_s,bottom_diffusive_salt_kg_m2_s,"
          <<"bottom_total_salt_kg_m2_s,top_solid_flux_kg_m2_s,"
          <<"internal_signed_kg_m2_s,internal_positive_kg_m2_s,"
          <<"internal_dissolution_kg_m2_s,upper_positive_2ell_kg_m2_s,"
          <<"upper_positive_3ell_kg_m2_s,upper_positive_4ell_kg_m2_s,"
          <<"dissolved_storage_rate_kg_m2_s,epsilon_min,epsilon_max,epsilon_mean,"
          <<"solid_kg_m2,dissolved_salt_kg_m2,pore_water_kg_m2,"
          <<"cumulative_bottom_salt_kg_m2,cumulative_top_solid_kg_m2,"
          <<"cumulative_internal_signed_kg_m2,cumulative_positive_kg_m2,"
          <<"cumulative_dissolution_kg_m2,salt_residual_kg_m2,water_residual_kg_m2\n";
        profiles<<std::setprecision(17)
          <<"case_id,suite,boundary_mode,time_s,time_over_tau,cell,xi,z_m,"
          <<"epsilon,concentration,reaction_kg_m3_s,cell_volume_m\n";

        auto write=[&](double time) {
            const auto inventory=model.inventory(y);
            const auto f=model.flux(y,time);
            const int k=model.inv();
            double signedRate=0,positive=0,dissolution=0,meanEps=0;
            double upper2=0,upper3=0,upper4=0;
            for(int i=0;i<a.cells;++i) {
                const double rate=f.reaction[i]*f.volume[i];
                signedRate+=rate; positive+=std::max(0.0,rate);
                dissolution+=std::max(0.0,-rate); meanEps+=y[i]*model.widths[i];
            }
            const double deff=std::sqrt(meanEps)*p.diffusion;
            const double av=p.avMultiplier*36*(1-meanEps)/p.grain;
            const double md=meanEps*deff/(1-p.ceq);
            const double ell=std::sqrt(md/(av*p.reaction));
            const double daAverage=y[2*a.cells]*y[2*a.cells]/(ell*ell);
            const double ks=p.beta*(1-p.birth)*p.reaction;
            const double daStarAverage=ks*ks/(md*av*p.reaction);
            for(int i=0;i<a.cells;++i) {
                const double z=y[2*a.cells]*model.centers[i];
                const double positiveRate=std::max(0.0,f.reaction[i]*f.volume[i]);
                if(z>=y[2*a.cells]-2*ell) upper2+=positiveRate;
                if(z>=y[2*a.cells]-3*ell) upper3+=positiveRate;
                if(z>=y[2*a.cells]-4*ell) upper4+=positiveRate;
            }
            const double dissolvedRate=f.salt[0]-f.salt[a.cells]-signedRate;
            const double saltResidual=inventory.salt+inventory.solid
                -start.salt-start.solid-y[k+2];
            const double waterResidual=inventory.water-start.water
                -(y[k+3]-y[k+2])+y[k+4];
            const double theta=time/diffusionTime;
            history<<a.id<<','<<a.suite<<','<<a.mode<<','<<time<<','<<theta<<','
              <<y[2*a.cells]<<','<<a.da<<','
              <<std::pow(y[2*a.cells]/ell0,2)<<','<<daAverage<<','<<daStarAverage<<','
              <<f.speed<<','<<f.topC<<','<<y[a.cells]<<','<<f.bottomC<<','
              <<f.evaporation<<','<<(a.mode=="dirichlet"?-f.bottomDiffusiveSalt:back)<<','
              <<f.liquid[0]<<','<<f.bottomAdvectiveSalt<<','<<f.bottomDiffusiveSalt<<','
              <<f.salt[0]<<','<<f.salt[a.cells]<<','<<signedRate<<','<<positive<<','
              <<dissolution<<','<<upper2<<','<<upper3<<','<<upper4<<','
              <<dissolvedRate<<','<<*std::min_element(y.begin(),y.begin()+a.cells)<<','
              <<*std::max_element(y.begin(),y.begin()+a.cells)<<','<<meanEps<<','
              <<inventory.solid<<','<<inventory.salt<<','<<inventory.water<<','
              <<y[k+2]<<','<<y[k]<<','<<y[k+1]<<','<<y[model.positiveInv()]<<','
              <<y[model.dissolutionInv()]<<','<<saltResidual<<','<<waterResidual<<'\n';
            if(snapshot(theta)) {
                for(int i=0;i<a.cells;++i)
                    profiles<<a.id<<','<<a.suite<<','<<a.mode<<','<<time<<','<<theta<<','
                      <<i<<','<<model.centers[i]<<','<<y[2*a.cells]*model.centers[i]<<','
                      <<y[i]<<','<<y[a.cells+i]<<','<<f.reaction[i]<<','<<f.volume[i]<<'\n';
            }
            if(std::abs(saltResidual)>2e-9 || std::abs(waterResidual)>2e-9)
                throw std::runtime_error("mass budget failure");
        };

        write(0.0);
        double last=0; std::string stop="forcing_end";
        const int outputs=int(std::llround(endTime/outputDt));
        for(int step=1;step<=outputs;++step) {
            const double target=(step==outputs)?endTime:step*outputDt;
            try { y=solver.advance(target); write(target); last=target; }
            catch(const std::exception& e) { stop=std::string("ALE_stop: ")+e.what(); break; }
        }
        std::ofstream status(prefix+"_status.txt");
        status<<"case_id="<<a.id<<"\nstatus="<<stop<<"\nlast_valid_time_s="
          <<std::setprecision(17)<<last<<"\nend_time_s="<<endTime<<"\n"
          <<"diffusion_time_s="<<diffusionTime<<"\ncells="<<a.cells<<"\npower="<<a.power
          <<"\nrtol="<<a.rtol<<"\nmax_step_s="<<maxStep<<"\nbdf_steps="<<solver.steps()<<"\n";
        std::cout<<a.id<<" status="<<stop<<" last="<<last<<" steps="<<solver.steps()<<'\n';
        return stop=="forcing_end"?0:2;
    } catch(const std::exception& e) {
        std::cerr<<e.what()<<'\n'; return 1;
    }
}
