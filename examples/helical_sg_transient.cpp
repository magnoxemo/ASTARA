/**
 * @file   helical_sg_transient.cpp
 * @brief  Helical-coil SG transient demo, reproducing Arda 2015 Fig. 11/14
 *         (feedwater-temperature step and feedwater-flow step).
 *
 * Initialises one NuScale-like helical-coil SG to its 100% steady state,
 * settles for 30 s, then applies the same 5 % perturbations described in
 * Section 7.2 of Arda (2015):
 *
 *     - feedwater inlet enthalpy bumped by +30 kJ/kg (~ +7 K)
 *     - feedwater mass flow bumped by +5 %
 *
 * Two 400-second runs are written as CSV files in the working directory
 * so the user can compare the time series visually with the published
 * figures.
 */

#include "astara/sg/HelicalCoilSteamGenerator.hpp"
#include "astara/props/IF97Water.hpp"

#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>

using astara::sg::HelicalCoilSteamGenerator;
using astara::sg::HelicalCoilSteamGeneratorParameters;
using astara::props::IF97Water;

namespace {

void runTransient(const std::string& name,
                  std::function<void(HelicalCoilSteamGenerator&)> perturb)
{
    IF97Water w;
    auto params = HelicalCoilSteamGeneratorParameters::nuscaleSMRTwoSG();
    HelicalCoilSteamGenerator sg(params, &w);
    sg.initialiseSteadyState(/*T_pi*/ 291.0 + 273.15,
                             /*W_p */ 354.0,                 // 708/2 kg/s per SG
                             /*W_fw*/ 35.6,                   // 71.25/2 kg/s per SG
                             /*T_fw*/ 148.5 + 273.15);

    std::ofstream out(name);
    out << "# t_s, L1_m, L2_m, L3_m, p_S_MPa, T_steam_C, T_P1_C\n";
    out << std::fixed << std::setprecision(4);

    constexpr double dt   = 0.05;
    constexpr double t_pre  = 200.0;    // settle (Arda's region-length time-constant ~ 100 s)
    constexpr double t_total = 600.0;   // 400 s post-perturbation
    bool perturbed = false;

    for (int i = 0; i * dt <= t_total; ++i) {
        const double t = i * dt;
        if (!perturbed && t >= t_pre) {
            perturb(sg);
            perturbed = true;
        }
        if (i % 20 == 0) {              // log every 1.0 s
            const auto& s = sg.state();
            out << t << ", "
                << s.L_1 << ", "
                << s.L_2 << ", "
                << sg.superheatedRegionLength_m() << ", "
                << s.p_S / 1.0e6 << ", "
                << sg.steamOutletTemperatureK() - 273.15 << ", "
                << s.T_P1 - 273.15 << "\n";
        }
        sg.timeStep(dt);
    }
    std::cout << "Wrote " << name << "\n";
}

}  // namespace

int main() {
    std::cout << "ASTARA helical-coil SG transient demo (Arda 2015 Fig. 11/14)\n";

    runTransient("helical_sg_feedwater_T_step.csv",
                 [](HelicalCoilSteamGenerator& sg) {
                     sg.inputs().feedwater_enthalpy_J_kg += 30.0e3;     // +7 K
                 });

    runTransient("helical_sg_feedwater_W_step.csv",
                 [](HelicalCoilSteamGenerator& sg) {
                     sg.inputs().feedwater_mass_flow_kg_s *= 1.05;
                     sg.inputs().steam_outlet_mass_flow_kg_s *= 1.05;
                 });

    return 0;
}
