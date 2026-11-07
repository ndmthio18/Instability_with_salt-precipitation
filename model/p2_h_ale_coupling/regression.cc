#include "bdf.hh"
#include <fstream>
#include <iomanip>
#include <iostream>

int main(int argc,char** argv) {
    try {
        int n=argc>1?std::stoi(argv[1]):640;
        double end=argc>2?std::stod(argv[2]):21600;
        double maxStep=argc>3?std::stod(argv[3]):60;
        p2h::ALE model(n,3);
        auto initial=model.initial(), y=initial;
        auto start=model.inventory(initial);
        p2h::BDF solver(model,initial,maxStep);
        std::ofstream log("history.csv"), fields("fields.csv");
        log<<std::setprecision(17)<<"time_s,height_m,epsilon_min,velocity_m_s,internal_fraction,top_solid_kg_m2,internal_solid_kg_m2,salt_residual_kg_m2,water_residual_kg_m2,evaporation_kg_m2_s,h_mass_m,dissolved_salt_kg_m2,pore_water_kg_m2\n";
        fields<<std::setprecision(17)<<"time_s,z_bottom_m,z_top_m,epsilon,concentration\n";
        auto write=[&](double t) {
            const auto m=model.inventory(y); const auto f=model.flux(y); const int k=model.inv();
            double sr=m.salt+m.solid-start.salt-start.solid-y[k+2];
            double wr=m.water-start.water-(y[k+3]-y[k+2])+y[k+4];
            log<<t<<','<<y[2*n]<<','<<*std::min_element(y.begin(),y.begin()+n)<<','<<f.speed<<','
               <<(y[k]+y[k+1]!=0?y[k+1]/(y[k]+y[k+1]):0)<<','<<y[k]<<','<<y[k+1]<<','<<sr<<','<<wr<<','
               <<f.evaporation<<','<<m.solid/(model.p.crystal*(1-model.p.birth))<<','<<m.salt<<','<<m.water<<'\n';
            log.flush();
            if(t==0 || t==300 || t==3600 || t==10800 || t==21600 || t==end)
                for(int i=0;i<n;++i) fields<<t<<','<<model.faces[i]*y[2*n]<<','<<model.faces[i+1]*y[2*n]<<','<<y[i]<<','<<y[n+i]<<'\n';
            if(std::abs(sr)>1e-7 || std::abs(wr)>1e-7) throw std::runtime_error("standalone mass budget failure");
        };
        write(0);
        for(double t=std::min(300.,end);;t=std::min(t+300,end)) {
            y=solver.advance(t); write(t); if(t==end) break;
        }
        std::cout<<"PASS integration n="<<n<<" t="<<end<<" BDF_steps="<<solver.steps()<<'\n';
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
