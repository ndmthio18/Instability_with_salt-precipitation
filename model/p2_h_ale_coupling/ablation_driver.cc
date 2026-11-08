// Built ONLY against an explicitly generated test boundary overlay.
// Production ALE equations and production initialization are not altered.
#include "bdf.hh"
#include "film.hh"
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>
struct Record { double t,c,adv,diff; };
int main(int argc,char** argv) {
    if(argc!=9) return 2;
    const double h=std::stod(argv[1]),c=std::stod(argv[2]),rho=std::stod(argv[3]);
    const bool equilibrium=std::stoi(argv[4]); const std::string mode=argv[5];
    std::ifstream input(argv[6]); std::string line; std::vector<Record> records;
    std::getline(input,line);
    while(std::getline(input,line)) { std::replace(line.begin(),line.end(),',',' '); std::istringstream in(line); Record r; in>>r.t>>r.c>>r.adv>>r.diff; records.push_back(r); }
    const double end=std::stod(argv[7]),step=std::stod(argv[8]);
    p2h::Parameters p; p.ceq=.2649408697676898; p.rho=rho; p.crystal=2165;
    p2h::ALE model(64,2,p); p2gfilm::GasFilm film;
    film.thickness=.0083; film.relativeHumidity=.4;
    model.gasFilm=[&](double x){return film.fluxFromMassFraction(x);};
    auto y=model.initial(h,c); const auto original=model.inventory(y);
    if(equilibrium) {
        // Conserve BOTH components including skeleton. No external reservoir.
        const double pore=original.water/(p.rho*(1-p.ceq));
        const double solid=original.solid+original.salt-p.ceq/(1-p.ceq)*original.water;
        const double newHeight=pore+solid/p.crystal;
        y=model.initial(newHeight,p.ceq);
        std::fill(y.begin(),y.begin()+model.n,pore/newHeight);
    }
    const auto initial=model.inventory(y);
    if(std::abs(initial.water-original.water)>1e-12*original.water ||
       std::abs(initial.salt+initial.solid-original.salt-original.solid)>1e-12*(original.salt+original.solid)) return 3;
    std::ofstream out("history.csv"),fields("fields.csv"),meta("initial.csv");
    out<<std::setprecision(17)<<"time_s,status,h_m,epsilon_min,top_added_kg_m2,internal_added_kg_m2,salt_input_kg_m2,liquid_input_kg_m2,water_residual_kg_m2,salt_residual_kg_m2\n";
    fields<<std::setprecision(17)<<"time_s,status,i,z_m,epsilon,c,reaction\n";
    meta<<std::setprecision(17)<<"h_original,h_initial,epsilon_initial,c_initial,initial_internal_conversion,salt_residual,water_residual\n"
        <<h<<','<<y[2*model.n]<<','<<y[0]<<','<<y[model.n]<<','<<initial.solid-original.solid<<','
        <<initial.salt+initial.solid-original.salt-original.solid<<','<<initial.water-original.water<<'\n';
    auto write=[&](double t,const p2h::Vec& state,const char* status){
        const auto m=model.inventory(state); const int k=model.inv();
        out<<t<<','<<status<<','<<state[2*model.n]<<','<<*std::min_element(state.begin(),state.begin()+model.n)<<','
            <<state[k]<<','<<state[k+1]<<','<<state[k+2]<<','<<state[k+3]<<','
            <<m.water-initial.water-state[k+3]+state[k+2]+state[k+4]<<','
            <<m.salt+m.solid-initial.salt-initial.solid-state[k+2]<<'\n';
        const auto f=model.flux(state);
        for(int i=0;i<model.n;++i) fields<<t<<','<<status<<','<<i<<','<<model.centers[i]*state[2*model.n]<<','<<state[i]<<','<<state[model.n+i]<<','<<f.reaction[i]<<'\n';
    };
    double time=0; write(0,y,"accepted");
    while(time<end-1e-12) {
        double dt=std::min(step,end-time); Record record=records.front();
        // Exact step-integrated input divided by dt. Split at each recorded
        // endpoint, so replay never smears a flux across a temporal boundary.
        for(const auto& r:records) { record=r; if(r.t>time+1e-12) break; }
        if(mode!="equilibrium_boundary") dt=std::min(dt,record.t-time);
        if(!(dt>0)) throw std::runtime_error("recorded input horizon exhausted");
        model.p.bottomC=mode=="equilibrium_boundary"?p.ceq:record.c;
        if(mode=="replay_total") model.forcedSalt=record.adv+record.diff;
        if(mode=="replay_advective") model.forcedSalt=record.adv;
        if(mode=="replay_diffusive") model.forcedSalt=record.diff;
        try {
            p2h::BDF solver(model,y,dt);
            solver.rejectedState=[&](const p2h::Vec& state,double t){write(time+t,state,"rejected");};
            auto next=solver.advance(dt); time+=dt; y=std::move(next); write(time,y,"accepted");
        } catch(const std::exception& e) { std::cerr<<"t="<<time<<" "<<e.what()<<'\n'; return 1; }
    }
    return 0;
}
