/**
 * @file   pump_trip.cpp
 * @brief  Reactor coolant pump (RCP) trip transient.
 *
 * Initialises the pump at rated speed with rated input power, then sets
 * P_d = 0 ("trip") and integrates for 60 s.  Writes the speed and
 * volumetric-flow trace to `pump_trip.csv`.  Reproduces the qualitative
 * shape of thesis Fig. 3.23 (homologous-pump coastdown).
 */

#include "astara/pump/HomologousPump.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>

using astara::pump::HomologousPump;
using astara::pump::HomologousPumpParameters;
using astara::pump::PumpCurveCoefficients;

int main() {
    HomologousPumpParameters p;
    p.curve.A0 = 90.0;
    p.curve.A1 = 1.0;
    p.curve.A2 = -1.0;
    p.rated_speed_rev_s          = 20.0;
    p.rated_volumetric_flow_m3_s = 6.0;
    p.loop_resistance_K_s2_m5    = 60.0 / 36.0;
    p.effective_flow_area_m2     = 0.4;
    p.loop_length_m              = 70.0;
    p.fluid_density_kg_m3        = 720.0;
    p.moment_of_inertia_kg_m2    = 1500.0;
    p.rated_input_power_W        = 720.0 * 9.80665 * 6.0 * 60.0;

    HomologousPump pump(p);
    pump.initialiseAtRated();
    pump.setInputPowerW(0.0);                 // pump trip

    constexpr double dt        = 0.01;
    constexpr double t_end_s   = 60.0;
    constexpr int    out_every = 10;          // every 100 ms
    const int        N_steps   = static_cast<int>(t_end_s / dt + 0.5);

    std::ofstream out("pump_trip.csv");
    out << "# t [s], N [rev/s], Q [m^3/s], H [m]\n";
    out << std::fixed << std::setprecision(6);
    auto write = [&]() {
        const auto& s = pump.state();
        out << s.t_s << ", " << s.speed_rev_s << ", " << s.volumetric_flow << ", "
            << pump.developedHead_m() << "\n";
    };
    write();

    for (int i = 1; i <= N_steps; ++i) {
        pump.timeStep(dt);
        if (i % out_every == 0) write();
    }
    std::cout << "Final speed = " << pump.state().speed_rev_s << " rev/s\n";
    std::cout << "Final flow  = " << pump.state().volumetric_flow << " m^3/s\n";
    std::cout << "CSV trace written to pump_trip.csv\n";
    return 0;
}
