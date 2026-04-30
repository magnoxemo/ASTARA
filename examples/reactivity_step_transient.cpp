/**
 * @file   reactivity_step_transient.cpp
 * @brief  Step-reactivity transient demonstration.
 *
 * Initialises a Westinghouse-like 4-loop PWR at full power and applies a
 * step in external reactivity, then integrates for 60 s.  Writes the time
 * trace of normalised power, average fuel/moderator temperatures, and
 * total reactivity (external + feedback) to `reactivity_step_transient.csv`.
 *
 * This reproduces qualitatively the kind of transient analysed in
 * Naghedolfeizi (1990) Section 5.1, Fig. 3.3 (positive reactivity insertion
 * with negative temperature feedback).
 *
 * Usage:
 *     ./reactivity_step_transient                # default +50 pcm step
 *     ./reactivity_step_transient -- -100        # -100 pcm step
 */

#include "astara/reactor/Reactor.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>

using astara::reactor::DelayedGroupConstants;
using astara::reactor::Reactor;
using astara::reactor::ReactivityModel;
using astara::reactor::ReactorThermalParameters;

namespace {

Reactor makeWestinghouseLikeReactor() {
    ReactorThermalParameters tp;
    tp.num_fuel_nodes        = 3;
    tp.num_moderator_nodes   = 6;
    tp.fuel_mass_total_kg    = 1.0e5;
    tp.fuel_cp_J_per_kgK     = 300.0;
    tp.fission_power_in_fuel = 0.974;
    tp.moderator_mass_total_kg = 1.5e4;
    tp.moderator_cp_J_per_kgK  = 5400.0;
    tp.mass_flow_rate_kg_s     = 17600.0;
    tp.overall_h_W_per_m2K     = 28000.0;
    tp.heat_transfer_area_m2   = 5400.0;
    tp.lower_plenum_mass_kg    = 5000.0;
    tp.upper_plenum_mass_kg    = 5000.0;
    tp.hot_leg_mass_kg         = 2500.0;
    tp.cold_leg_mass_kg        = 2500.0;

    auto groups = DelayedGroupConstants::u235SixGroup();
    constexpr double P_rated = 3.4e9;       // 3400 MWth

    ReactivityModel rm;
    rm.alpha_fuel_per_K      = -2.5e-5;     // Doppler
    rm.alpha_moderator_per_K = -3.0e-4;     // moderator (negative)
    return Reactor(groups, tp, P_rated, rm);
}

}  // namespace

int main(int argc, char** argv) {
    // Reactivity step in pcm (1 pcm = 1e-5).  Default = +50 pcm.
    double rho_pcm = 50.0;
    if (argc >= 2) {
        try {
            rho_pcm = std::stod(argv[1]);
        } catch (const std::exception& e) {
            std::cerr << "Failed to parse reactivity from '" << argv[1]
                      << "': " << e.what() << "\n";
            return 2;
        }
    }
    const double rho_step = rho_pcm * 1.0e-5;

    auto reactor = makeWestinghouseLikeReactor();
    constexpr double T_inlet = 559.0;       // K, ~ 547 F
    reactor.initialiseSteadyState(/*n0=*/1.0, T_inlet);

    // Apply the step.
    reactor.reactivity().rho_external = rho_step;

    constexpr double dt        = 1.0e-3;     // 1 ms
    constexpr double t_end_s   = 60.0;
    constexpr int    out_every = 50;         // log every 50 ms
    const int        N_steps   = static_cast<int>(t_end_s / dt + 0.5);

    const std::string filename = "reactivity_step_transient.csv";
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Cannot open " << filename << " for writing\n";
        return 1;
    }
    out << "# t [s], n/n0, T_fuel_avg [K], T_mod_avg [K], rho_total [pcm], "
           "T_hot [K], T_cold [K]\n";
    out << std::fixed << std::setprecision(6);

    auto write_row = [&](const Reactor& r) {
        const auto& s = r.state();
        const double Tf = s.thermal.averageFuelTemperature();
        const double Tm = s.thermal.averageModeratorTemperature();
        const double rho_total = r.reactivity().evaluate(Tf, Tm);
        out << s.t_s << ", "
            << s.kinetics.power() << ", "
            << Tf << ", " << Tm << ", "
            << rho_total / 1e-5 << ", "
            << s.thermal.T_hot_leg << ", "
            << s.thermal.T_cold_leg << "\n";
    };
    write_row(reactor);

    for (int i = 1; i <= N_steps; ++i) {
        reactor.timeStep(dt);
        if (i % out_every == 0) write_row(reactor);
    }

    std::cout << "Step reactivity = " << rho_pcm << " pcm; integrated "
              << t_end_s << " s in " << N_steps << " steps.\n";
    std::cout << "Final n/n0 = " << reactor.state().kinetics.power() << "\n";
    std::cout << "Final T_hot_leg = " << reactor.state().thermal.T_hot_leg << " K\n";
    std::cout << "CSV trace written to " << filename << "\n";
    return 0;
}
