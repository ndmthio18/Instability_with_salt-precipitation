#include "bdf.hh"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

int main(int argc,char** argv) {
    if(argc!=4) {
        std::cerr<<"usage: av_mechanism AV_MULTIPLIER CELLS END_S\n";
        return 2;
    }
    double lastAccepted=0;
    try {
        p2h::Parameters p;
        p.ceq=.264940870;
        p.reaction=1e-3;
        p.avMultiplier=std::stod(argv[1]);
        p.splitReactionLedgers=true;
        p.bottomC=p.ceq+2e-4;
        const int cells=std::stoi(argv[2]);
        const double end=std::stod(argv[3]);
        p2h::ALE model(cells,2,p);
        auto y=model.initial(3.5e-6,p.ceq+1e-4);
        p2h::BDF solver(model,y,.05);
        std::ofstream history("mechanism_history.csv"),profiles("mechanism_profiles.csv");
        solver.rejectedState=[&](const p2h::Vec& rejected,double time) {
            std::ofstream output("mechanism_first_rejected.csv");
            output<<std::setprecision(17)<<"time_s,i,z_bottom_m,z_top_m,epsilon,c,av_m_inv,reaction_kg_m3_s,deff_m2_s\n";
            const auto flux=model.flux(rejected);
            for(int i=0;i<cells;++i) output<<time<<','<<i<<','<<model.faces[i]*rejected[2*cells]<<','
                <<model.faces[i+1]*rejected[2*cells]<<','<<rejected[i]<<','<<rejected[cells+i]<<','
                <<model.specificArea(rejected[i])<<','<<flux.reaction[i]<<','<<std::sqrt(rejected[i])*model.p.diffusion<<'\n';
        };
        history<<std::setprecision(17)<<"time_s,av_multiplier,height_m,top_velocity_m_s,top_c,positive_rate_kg_m2_s,dissolution_rate_kg_m2_s,cumulative_positive_kg_m2,cumulative_dissolution_kg_m2,net_internal_solid_kg_m2,epsilon_min,salt_residual_kg_m2,water_residual_kg_m2\n";
        profiles<<std::setprecision(17)<<"time_s,i,z_bottom_m,z_top_m,epsilon,c,av_m_inv,reaction_kg_m3_s,deff_m2_s\n";
        const auto initial=model.inventory(y);
        auto write=[&](double time,bool profile) {
            const auto f=model.flux(y); const auto inventory=model.inventory(y); const int k=model.inv();
            double positive=0,dissolution=0;
            for(int i=0;i<cells;++i) {
                positive+=std::max(0.,f.reaction[i])*f.volume[i];
                dissolution+=std::max(0.,-f.reaction[i])*f.volume[i];
                if(profile) profiles<<time<<','<<i<<','<<model.faces[i]*y[2*cells]<<','<<model.faces[i+1]*y[2*cells]<<','
                    <<y[i]<<','<<y[cells+i]<<','<<model.specificArea(y[i])<<','<<f.reaction[i]<<','<<std::sqrt(y[i])*model.p.diffusion<<'\n';
            }
            const double cumulativePositive=y[model.positiveInv()];
            const double cumulativeDissolution=y[model.dissolutionInv()];
            const double saltResidual=inventory.salt+inventory.solid-initial.salt-initial.solid-y[k+2];
            const double waterResidual=inventory.water-initial.water-(y[k+3]-y[k+2])+y[k+4];
            history<<time<<','<<p.avMultiplier<<','<<y[2*cells]<<','<<f.speed<<','<<f.topC<<','<<positive<<','<<dissolution<<','
                <<cumulativePositive<<','<<cumulativeDissolution<<','<<y[k+1]<<','
                <<*std::min_element(y.begin(),y.begin()+cells)<<','<<saltResidual<<','<<waterResidual<<'\n';
        };
        write(0,true);
        for(double time=.1;time<=end+1e-12;time+=.1) {
            const double target=std::min(time,end);
            y=solver.advance(target);
            lastAccepted=target;
            write(target,std::abs(std::round(target)-target)<1e-10 || target==end);
            if(target==end) break;
        }
        std::ofstream("mechanism_status.txt")<<"COMPLETED "<<end<<"\n";
        std::cout<<"PASS mechanism av="<<p.avMultiplier<<" cells="<<cells<<" end="<<end<<'\n';
    } catch(const std::exception& e) {
        std::ofstream("mechanism_status.txt")<<"FAILED last_accepted_s="<<lastAccepted<<' '<<e.what()<<"\n";
        std::cerr<<e.what()<<'\n';
        return 1;
    }
}
