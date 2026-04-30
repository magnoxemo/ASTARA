/**
 * @file   pressurizer_heatup.cpp
 * @brief  Pressurizer heater-on / spray-on demonstration.
 *
 * Initialises the pressurizer at 15.5 MPa, 8 m water level.  At t = 30 s
 * turns on the heater (1.4 MW); at t = 90 s turns on the spray (5 kg/s of
 * cold cold-leg water).  Writes pressure and water-level trace to
 * `pressurizer_heatup.csv`.
 */

#include "astara/pressurizer/Pressurizer.hpp"
#include "astara/props/IF97Water.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>

using astara::pressurizer::Pressurizer;
using astara::pressurizer::PressurizerParameters;
using astara::props::IF97Water;

int main() {
    IF97Water water;
    PressurizerParameters geom;
    geom.cross_section_area_m2 = 4.0;
    geom.total_height_m        = 13.0;

    Pressurizer pz(geom, &water);
    pz.initialiseSteadyState(15.5e6, 8.0);

    constexpr double dt        = 0.05;
    constexpr double t_end_s   = 150.0;
    constexpr int    out_every = 20;          // every 1 s
    const int        N_steps   = static_cast<int>(t_end_s / dt + 0.5);

    std::ofstream out("pressurizer_heatup.csv");
    out << "# t [s], P [Pa], Lw [m], heater [W], spray [kg/s]\n";
    out << std::fixed << std::setprecision(6);
    auto write = [&]() {
        const auto& s = pz.state();
        out << s.t_s << ", " << s.pressure_Pa << ", " << s.water_level_m << ", "
            << pz.inputs().heater_power_W << ", "
            << pz.inputs().spray_mass_flow_kg_s << "\n";
    };
    write();

    for (int i = 1; i <= N_steps; ++i) {
        const double t = pz.state().t_s;
        if (t >= 30.0 && pz.inputs().heater_power_W == 0.0) {
            pz.inputs().heater_power_W = 1.4e6;
        }
        if (t >= 90.0 && pz.inputs().spray_mass_flow_kg_s == 0.0) {
            pz.inputs().spray_mass_flow_kg_s = 5.0;
            pz.inputs().spray_enthalpy_J_kg  = 1.3e6;
        }
        pz.timeStep(dt);
        if (i % out_every == 0) write();
    }
    std::cout << "Final P  = " << pz.state().pressure_Pa / 1e6 << " MPa\n";
    std::cout << "Final Lw = " << pz.state().water_level_m << " m\n";
    std::cout << "CSV trace written to pressurizer_heatup.csv\n";
    return 0;
}
