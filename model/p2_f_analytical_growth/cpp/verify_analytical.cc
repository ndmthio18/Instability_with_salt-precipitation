#include "analyticalgrowth.hh"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool value, const char* message)
{ if (!value) throw std::runtime_error(message); }

double exactHeight(const p2f::Parameters& p, double h0, double j, double time)
{
    const double g = p2f::growthPrefactor(p, j);
    const double k = p.heightCoefficient();
    auto residual = [&](double h) {
        return h-h0 + std::log(h/h0)/k - g*time;
    };
    double low = h0, high = h0 + g*time + 1.0/k;
    for (int i = 0; i < 160; ++i)
    {
        const double mid = 0.5*(low+high);
        if (residual(mid) > 0.0) high = mid; else low = mid;
    }
    return 0.5*(low+high);
}
}

int main()
{
    try
    {
        p2f::Parameters p;
        p.check();
        require(std::abs(p.equilibriumMassFraction-0.2649408697676898) < 1e-15,
                "mole-to-mass conversion drift");
        require(p2f::growthVelocity(p, 1e-5, 0.0) == 0.0,
                "zero evaporation grew crust");
        const double j = 1.8e-5;
        const double h0 = 1e-5;
        const double exact = exactHeight(p, h0, j, 2000.0);
        double previous = 0.0;
        for (double dt : {20.0, 10.0, 5.0})
        {
            auto state = p2f::activate(p, h0);
            const int steps = static_cast<int>(std::llround(2000.0/dt));
            for (int i = 0; i < steps; ++i)
            {
                auto trial = p2f::evaluate(state, p, j, dt);
                p2f::commit(state, trial);
            }
            const double error = std::abs(state.height-exact);
            if (previous > 0.0) require(previous/error > 1.9, "Euler did not converge");
            previous = error;
        }

        auto state = p2f::activate(p, h0);
        auto a = p2f::evaluate(state, p, j, 1.0);
        auto b = p2f::evaluate(state, p, j, 1.0);
        require(a.proposed.height == b.proposed.height, "trial not repeatable");
        p2f::rollback(state, a);
        bool refused = false;
        try { p2f::commit(state, a); } catch (const p2f::DomainExit&) { refused = true; }
        require(refused, "discarded trial committed");
        p2f::commit(state, b);
        refused = false;
        try { p2f::commit(state, b); } catch (const p2f::DomainExit&) { refused = true; }
        require(refused, "trial committed twice");

        const auto c = p2f::evaluate(state, p, j, 1.0);
        const auto good = p2f::validity(p, c, p.equilibriumMassFraction, 1.0);
        require(good.valid, "declared benchmark outside validity envelope");
        const auto bad = p2f::validity(p, c, 0.95*p.equilibriumMassFraction, 0.5);
        require(!bad.valid, "insufficient supply was not flagged");

        std::cout << std::setprecision(17)
                  << "PASS analytical growth, inventories, convergence, lifecycle, validity\n"
                  << "VECTOR equilibrium_mass_fraction=" << p.equilibriumMassFraction
                  << " solution_density=" << p.rho()
                  << " effective_diffusivity=" << p.effectiveDiffusivity()
                  << " velocity=" << p2f::growthVelocity(p, h0, j)
                  << " reference_height=" << exact << '\n';
    }
    catch (const std::exception& error)
    {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
