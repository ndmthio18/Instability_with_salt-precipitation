#include "ale.hh"
#include "bdf.hh"
#include "pratcrust.hh"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace fs = std::filesystem;

int main(int argc, char** argv)
{
    if (argc != 8)
        throw std::runtime_error("usage: counterfactual OUT CELLS C_INIT C_BOTTOM HORIZON DT MODE");
    const fs::path out = argv[1];
    const int cells = std::stoi(argv[2]);
    const double cInitial = std::stod(argv[3]);
    const double cBottom = std::stod(argv[4]);
    const double horizon = std::stod(argv[5]);
    const double outputDt = std::stod(argv[6]);
    const std::string mode = argv[7];
    if (mode != "actual_mix" && mode != "equilibrium_reconstruction")
        throw std::runtime_error("unknown reconstruction mode");
    fs::create_directories(out);

    p2h::Parameters p;
    p.ceq = 0.2649408697676898;
    p.rho = 1202.1685774416314;
    p.crystal = 2165.0;
    p.birth = 0.43;
    p.diffusion = 1.3e-9;
    p.reaction = 1e-3;
    p.avMultiplier = 1.0;
    p.beta = 2.0;
    p.grain = 25e-6;
    p.bottomC = cBottom;
    p.splitReactionLedgers = true;
    const double height = 3.5001270242627885e-6;
    p2h::ALE model(cells, 2.0, p);
    prat::GasFilm film;
    film.thickness = 0.0083;
    film.relativeHumidity = 0.4;
    model.gasFilm = [&film](double c) { return film.fluxFromMassFraction(c); };
    auto state = model.initial(height, cInitial);
    const auto initialInventory = model.inventory(state);

    // Diagnostic conservative reconstruction: at fixed pore-liquid mass, the
    // liquid salt increase is supplied by a declared interfacial solid ledger
    // and the same mass of pore water is returned to the soil ledger. Geometry
    // and porosity are intentionally held fixed to isolate concentration.
    const double originalMix = 0.2640671415653032;
    double solidAdjustment = 0.0;
    double soilWaterReturn = 0.0;
    if (mode == "equilibrium_reconstruction") {
        const double poreLiquid = p.rho * p.birth * height;
        const double delta = poreLiquid * (p.ceq - originalMix);
        solidAdjustment = -delta;
        soilWaterReturn = delta;
    }

    std::ofstream metadata(out / "metadata.csv");
    metadata << std::setprecision(17)
             << "cells,c_initial,c_bottom,ceq,height_m,mode,initial_pore_salt_kg_m2,"
                "initial_pore_water_kg_m2,geometric_solid_kg_m2,solid_adjustment_kg_m2,"
                "soil_water_return_kg_m2\n"
             << cells << ',' << cInitial << ',' << cBottom << ',' << p.ceq << ',' << height << ','
             << mode << ',' << initialInventory.salt << ',' << initialInventory.water << ','
             << initialInventory.solid << ',' << solidAdjustment << ',' << soilWaterReturn << '\n';

    std::ofstream history(out / "history.csv");
    history << std::setprecision(17)
            << "time_s,status,h_m,epsilon_min,epsilon_max,c_bottom_cell,c_top_cell,"
               "reaction_positive_kg_m2_s,reaction_dissolution_kg_m2_s,bottom_liquid_kg_m2_s,"
               "bottom_salt_kg_m2_s,bottom_advective_salt_kg_m2_s,bottom_diffusive_salt_kg_m2_s,"
               "top_c,top_speed_m_s,evaporation_kg_m2_s,pore_salt_kg_m2,pore_water_kg_m2,"
               "geometric_solid_kg_m2,adjusted_system_salt_kg_m2,adjusted_system_water_kg_m2,"
               "cumulative_bottom_salt_kg_m2,cumulative_bottom_liquid_kg_m2,"
               "cumulative_internal_net_kg_m2,cumulative_internal_positive_kg_m2,"
               "cumulative_internal_dissolution_kg_m2\n";

    auto writeState = [&](double time, const char* status, const p2h::Vec& y) {
        const auto flux = model.flux(y);
        const auto inv = model.inventory(y);
        double positive = 0.0, dissolution = 0.0;
        for (int i = 0; i < model.n; ++i) {
            positive += std::max(0.0, flux.reaction[i]) * flux.volume[i];
            dissolution += std::max(0.0, -flux.reaction[i]) * flux.volume[i];
        }
        const int k = model.inv();
        history << time << ',' << status << ',' << y[2*model.n] << ','
                << *std::min_element(y.begin(), y.begin()+model.n) << ','
                << *std::max_element(y.begin(), y.begin()+model.n) << ','
                << y[model.n] << ',' << y[2*model.n-1] << ',' << positive << ',' << dissolution << ','
                << flux.liquid[0] << ',' << flux.salt[0] << ',' << p.bottomC*flux.liquid[0] << ','
                << flux.salt[0]-p.bottomC*flux.liquid[0] << ',' << flux.topC << ',' << flux.speed << ','
                << flux.evaporation << ',' << inv.salt << ',' << inv.water << ',' << inv.solid << ','
                << inv.salt+inv.solid+solidAdjustment << ',' << inv.water+soilWaterReturn << ','
                << y[k+2] << ',' << y[k+3] << ',' << y[k+1] << ',' << y[model.positiveInv()] << ','
                << y[model.dissolutionInv()] << '\n';
    };
    writeState(0.0, "initial", state);

    p2h::BDF integrator(model, state, 0.5);
    std::ofstream rejected(out / "first_rejected.csv");
    rejected << std::setprecision(17)
             << "time_s,cell,z_bottom_m,z_top_m,epsilon,c,reaction_kg_m3_s\n";
    bool wroteRejected = false;
    integrator.rejectedState = [&](const p2h::Vec& y, double time) {
        if (wroteRejected) return;
        const auto flux = model.flux(y);
        const double h = y[2*model.n];
        for (int i = 0; i < model.n; ++i)
            rejected << time << ',' << i << ',' << model.faces[i]*h << ',' << model.faces[i+1]*h << ','
                     << y[i] << ',' << y[model.n+i] << ',' << flux.reaction[i] << '\n';
        wroteRejected = true;
    };

    std::string finalStatus = "COMPLETED";
    double acceptedTime = 0.0;
    try {
        for (double target = outputDt; target <= horizon + 1e-12; target += outputDt) {
            state = integrator.advance(std::min(target, horizon));
            acceptedTime = std::min(target, horizon);
            writeState(acceptedTime, "accepted", state);
        }
    } catch (const std::exception& e) {
        finalStatus = std::string("STOPPED: ") + e.what();
    }
    std::ofstream status(out / "status.txt");
    status << finalStatus << '\n' << std::setprecision(17) << acceptedTime << '\n';
    std::cout << out.string() << ',' << finalStatus << ',' << acceptedTime << '\n';
    return 0;
}
