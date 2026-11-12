#include "bdf.hh"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

struct Forcing {
    std::vector<double> t,j,b;
    static std::vector<std::string> split(const std::string& line) {
        std::vector<std::string> out; std::stringstream stream(line); std::string item;
        while(std::getline(stream,item,',')) out.push_back(item);
        return out;
    }
    explicit Forcing(const std::string& path) {
        std::ifstream input(path); if(!input) throw std::runtime_error("cannot open forcing CSV");
        std::string line; std::getline(input,line); const auto header=split(line);
        auto index=[&](const std::string& name){
            auto it=std::find(header.begin(),header.end(),name);
            if(it==header.end()) throw std::runtime_error("forcing column missing: "+name);
            return int(it-header.begin());
        };
        const int itime=index("elapsed_time_s"), ij=index("evaporation_water_kg_m2_s"), ib=index("b_soil_kg_m2_s");
        while(std::getline(input,line)) {
            if(line.empty()) continue;
            const auto row=split(line);
            t.push_back(std::stod(row[itime])); j.push_back(std::stod(row[ij])); b.push_back(std::stod(row[ib]));
        }
        if(t.size()<2 || t.front()!=0.0) throw std::runtime_error("invalid forcing support");
        for(std::size_t i=1;i<t.size();++i) if(!(t[i]>t[i-1])) throw std::runtime_error("forcing times not increasing");
    }
    double value(const std::vector<double>& y,double time) const {
        if(time<=t.front()) return y.front();
        if(time>=t.back()) return y.back();
        const auto upper=std::upper_bound(t.begin(),t.end(),time); const std::size_t i=std::size_t(upper-t.begin()-1);
        const double w=(time-t[i])/(t[i+1]-t[i]); return y[i]+w*(y[i+1]-y[i]);
    }
    double evaporation(double time) const { return value(j,time); }
    double backDiffusion(double time) const { return value(b,time); }
};

int main(int argc,char** argv) {
    try {
        if(argc<3) throw std::runtime_error("usage: x008_prescribed_flux_ale forcing.csv history.csv [cells] [output_dt_s]");
        const Forcing forcing(argv[1]); const int n=argc>3?std::stoi(argv[3]):64;
        const double outputDt=argc>4?std::stod(argv[4]):10.0;
        p2h::Parameters p;
        p.ceq=0.264940869768; p.rho=1202.168577; p.crystal=2165.0; p.birth=0.43;
        p.diffusion=1.3e-9; p.reaction=1.0e-3; p.avMultiplier=1.0; p.beta=2.0;
        p.grain=25.0e-6; p.initialHeight=3.5e-6; p.bottomC=p.ceq; p.splitReactionLedgers=true;
        p.allowCompleteDissolution=true;
        p2h::ALE model(n,3,p);
        model.evaporationHistory=[&](double time){return forcing.evaporation(time);};
        model.bottomBackDiffusion=[&](double time){return forcing.backDiffusion(time);};
        auto initial=model.initial(p.initialHeight,p.ceq), y=initial; const auto start=model.inventory(initial);
        p2h::BDF solver(model,initial,2.0,1e-8,true);
        std::ofstream log(argv[2]);
        log<<std::setprecision(17)<<"time_s,height_m,epsilon_min,epsilon_bottom,epsilon_top,concentration_bottom,concentration_top_cell,top_boundary_concentration,growth_velocity_m_s,evaporation_kg_m2_s,b_soil_kg_m2_s,bottom_total_salt_flux_kg_m2_s,top_solid_rate_kg_m2_s,internal_positive_rate_kg_m2_s,internal_dissolution_rate_kg_m2_s,solid_kg_m2,dissolved_salt_kg_m2,pore_water_kg_m2,cumulative_bottom_salt_kg_m2,cumulative_evaporation_kg_m2,salt_residual_kg_m2,water_residual_kg_m2\n";
        auto write=[&](double time) {
            const auto m=model.inventory(y);
            const auto f=model.flux(y,time);
            const int k=model.inv();
            double positive=0.0,dissolution=0.0;
            for(int i=0;i<n;++i) {
                const double rate=f.reaction[i]*f.volume[i]; positive+=std::max(0.0,rate); dissolution+=std::max(0.0,-rate);
            }
            const double saltResidual=m.salt+m.solid-start.salt-start.solid-y[k+2];
            const double waterResidual=m.water-start.water-(y[k+3]-y[k+2])+y[k+4];
            log<<time<<','<<y[2*n]<<','<<*std::min_element(y.begin(),y.begin()+n)<<','<<y[0]<<','<<y[n-1]<<','
               <<y[n]<<','<<y[2*n-1]<<','<<f.topC<<','<<f.speed<<','<<f.evaporation<<','<<forcing.backDiffusion(time)<<','
               <<f.salt[0]<<','<<f.salt[n]<<','<<positive<<','<<dissolution<<','<<m.solid<<','<<m.salt<<','<<m.water<<','
               <<y[k+2]<<','<<y[k+4]<<','<<saltResidual<<','<<waterResidual<<'\n'; log.flush();
            if(std::abs(saltResidual)>2e-8 || std::abs(waterResidual)>2e-8) throw std::runtime_error("mass budget failure");
        };
        write(0.0); double last=0.0; std::string stop="forcing_end";
        for(double target=std::min(outputDt,forcing.t.back()); target<=forcing.t.back()+1e-12; target=std::min(target+outputDt,forcing.t.back())) {
            if(target<=last) break;
            try { y=solver.advance(target); write(target); last=target; }
            catch(const std::exception& e) { stop=std::string("ALE_stop: ")+e.what(); break; }
            if(target==forcing.t.back()) break;
        }
        std::ofstream status("ale_status.txt");
        status<<"status="<<stop<<"\nlast_valid_time_s="<<std::setprecision(17)<<last<<"\nforcing_end_s="<<forcing.t.back()<<"\ncells="<<n<<"\noutput_dt_s="<<outputDt<<"\nbdf_steps="<<solver.steps()<<"\n";
        std::cout<<"ALE_DONE status="<<stop<<" last_valid_time_s="<<last<<" steps="<<solver.steps()<<'\n';
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
