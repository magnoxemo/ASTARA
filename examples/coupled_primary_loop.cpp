/**
 * @file   primary_loop_load_change.cpp
 * @brief  Coupled-loop demonstration: a load-following transient.
 *
 * Sets up a complete PWR primary loop (reactor + Ali Model D U-tube SG +
 * homologous pump + pressurizer) with all three controllers attached
 * (reactor turbine-following, pressurizer heater+spray, SG three-element
 * feedwater).  Initialises everything at full-power steady state and then
 * drops the turbine load fraction from 100% to 80% at t = 30 s.
 *
 * trace:
 *
 *   - Reactor power: ramps from 1.0 down toward ~ 0.80 over ~ 100 s as the
 *     reactor controller withdraws positive reactivity (rate-limited at
 *     ~ 1.4 pcm/s -> 100 pcm in ~ 70 s).
 *   - Primary T_avg: tracks the load-dependent setpoint (562 K at zero
 *     load, 583 K at full load -> 578.8 K at 80% load).
 *   - SG pressure: rises briefly as the steam demand falls faster than
 *     the reactor can ramp down, then settles at a new equilibrium.
 *   - Pressurizer pressure: held within ~ 0.05 MPa of the 15.5 MPa
 *     setpoint by heater/spray action.
 *   - Drum level: rides on the three-element controller; level error
 *     should stay well under 0.1 m throughout.
 */

#include "astara/primary/PrimaryLoop.hpp"
#include "astara/props/IF97Water.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>

using astara::primary::PrimaryLoop;
using astara::reactor::Reactor;
using astara::reactor::ReactorThermalParameters;
using astara::reactor::DelayedGroupConstants;
using astara::reactor::ReactivityModel;
using astara::sg::AliSteamGenerator;
using astara::sg::AliSteamGeneratorParameters;
using astara::pump::HomologousPump;
using astara::pump::HomologousPumpParameters;
using astara::pressurizer::Pressurizer;
using astara::pressurizer::PressurizerParameters;
using astara::props::IF97Water;
using astara::control::ReactorController;
using astara::control::PressurizerController;
using astara::control::ThreeElementController;
using astara::control::PIDConfig;

namespace {

    IF97Water& sharedWater() {
        static IF97Water w;
        return w;
    }


    std::unique_ptr<Reactor> makeReactor() {
        ReactorThermalParameters tp;
        tp.num_fuel_nodes          = 3;
        tp.num_moderator_nodes     = 6;
        tp.fuel_mass_total_kg      = 1.0e5;
        tp.fuel_cp_J_per_kgK       = 300.0;
        tp.fission_power_in_fuel   = 0.974;
        tp.moderator_mass_total_kg = 1.5e4;
        tp.moderator_cp_J_per_kgK  = 5400.0;
        tp.mass_flow_rate_kg_s     = 17600.0;
        tp.overall_h_W_per_m2K     = 28000.0;
        tp.heat_transfer_area_m2   = 5400.0;
        tp.lower_plenum_mass_kg    = 5000.0;
        tp.upper_plenum_mass_kg    = 5000.0;
        tp.hot_leg_mass_kg         = 2500.0;
        tp.cold_leg_mass_kg        = 2500.0;

        ReactivityModel rm;
        rm.alpha_fuel_per_K      = -2.5e-5;     // Doppler
        rm.alpha_moderator_per_K = -3.0e-4;     // moderator-temperature feedback

        auto r = std::make_unique<Reactor>(DelayedGroupConstants::u235SixGroup(), tp,
                /*P_rated_W=*/3.4e9, rm);
        r->initialiseSteadyState(/*n0=*/1.0, /*T_inlet=*/559.0);
        return r;
    }

    std::unique_ptr<AliSteamGenerator> makeSG() {
        auto p  = AliSteamGeneratorParameters::westinghouseModelD5();
        auto sg = std::make_unique<AliSteamGenerator>(p, &sharedWater());
        sg->initialiseSteadyState(/*T_pi_in*/ 597.0,        // hot-leg out
                /*W_p   */ 4400.0,        // 17600 / 4 loops
                /*P_sg  */ 6.9e6,
                /*L_dw  */ 3.0);
        return sg;
    }

    std::unique_ptr<HomologousPump> makePump() {
        HomologousPumpParameters p;
        p.curve.A0                   = 90.0;
        p.curve.A1                   = 1.0;
        p.curve.A2                   = -1.0;
        p.rated_speed_rev_s          = 20.0;
        p.rated_volumetric_flow_m3_s = 6.0;
        p.loop_resistance_K_s2_m5    = 60.0 / 36.0;
        p.effective_flow_area_m2     = 0.4;
        p.loop_length_m              = 70.0;
        p.fluid_density_kg_m3        = 720.0;
        p.moment_of_inertia_kg_m2    = 1500.0;
        p.rated_input_power_W        = 720.0 * 9.80665 * 6.0 * 60.0;
        auto pump = std::make_unique<HomologousPump>(p);
        pump->initialiseAtRated();
        return pump;
    }

    std::unique_ptr<Pressurizer> makePressurizer() {
        PressurizerParameters geom;
        geom.cross_section_area_m2 = 4.0;
        geom.total_height_m        = 13.0;
        auto pz = std::make_unique<Pressurizer>(geom, &sharedWater());
        pz->initialiseSteadyState(/*P*/ 15.5e6, /*Lw*/ 8.0);
        return pz;
    }

// ---------------------------------------------------------------------------
// Controllers.
// ---------------------------------------------------------------------------

    std::unique_ptr<ReactorController> makeReactorController() {
        ReactorController::Config cfg;
        // Modest gains: the reactor controller is slow on purpose -- it should
        // not chase fast disturbances, temperature feedbacks will just take care of that.
        cfg.pid.Kp    = 1.0e-6;        // (1/s) of demanded reactivity per K of T_avg error
        cfg.pid.Ki    = 1.0e-7;
        cfg.pid.u_min = -1.0;
        cfg.pid.u_max =  1.0;
        cfg.T_avg_setpoint_no_load_K   = 562.0;
        cfg.T_avg_setpoint_full_load_K = 583.0;
        cfg.max_reactivity_rate        = 1.4e-5;     // ~ 1.4 pcm/s rod-bank speed
        cfg.rho_min                    = -5.0e-3;
        cfg.rho_max                    =  5.0e-3;
        auto c = std::make_unique<ReactorController>(cfg);
        c->setTurbineLoadFraction(1.0);              // start at 100% load
        return c;
    }

    std::unique_ptr<PressurizerController> makePressurizerController() {
        PressurizerController::Config cfg;
        cfg.pressure_setpoint_Pa     = 15.5e6;
        cfg.heater_max_W             = 1.4e6;
        cfg.heater_steady_state_W    = 0.2e6;
        cfg.heater_P_gain_W_per_Pa   = 0.5;
        cfg.spray_threshold_Pa       = 1.0e5;
        cfg.spray_max_kg_s           = 5.0;
        cfg.spray_gain_kg_per_Pa     = 5.0e-5;
        cfg.spray_dead_band_Pa       = 1.0e4;
        return std::make_unique<PressurizerController>(cfg);
    }

    std::unique_ptr<ThreeElementController> makeFeedwaterController() {
        PIDConfig pi;
        pi.Kp    = 50.0;     // kg/s of feedwater trim per metre of level error
        pi.Ki    = 5.0;
        pi.u_min = -200.0;
        pi.u_max =  200.0;
        return std::make_unique<ThreeElementController>(/*L_setpoint=*/3.0, pi);
    }

}  // namespace


int main() {
    std::cout << "ASTARA primary-loop demo: 100% -> 80% load reduction\n";

    PrimaryLoop loop(makeReactor(), makeSG(), makePump(), makePressurizer());
    loop.setReactorController(makeReactorController());
    loop.setPressurizerController(makePressurizerController());
    loop.setFeedwaterController(makeFeedwaterController());

    // We need a handle on the reactor controller to drop the load
    // fraction at the right moment.  The PrimaryLoop owns it but exposes
    // accessors via a small helper -- in this simple demo we'll re-issue
    // the entire controller.  (A future API revision could add a setter.)
    // For now the simplest path is to grab it before transferring
    // ownership, but the demo's structure is clear enough as written:
    // we'll re-attach a fresh controller at t = 30 s with the new load.
    auto rebuildReactorControllerAtLoad = [](double load_fraction) {
        auto c = makeReactorController();
        c->setTurbineLoadFraction(load_fraction);
        // Preload the integrator so we don't reset the rod position when
        // we swap the controller.  In a real plant this is handled by
        // bumpless transfer; here we simply build a fresh one.
        return c;
    };

    std::ofstream out("primary_loop_load_change.csv");
    out << "# t_s, P_reactor_MW, T_avg_K, T_hot_K, T_cold_K, "
        << "rho_external_pcm, mdot_kg_s, "
        << "P_pzr_MPa, L_pzr_m, Q_heater_kW, W_spray_kg_s, "
        << "P_sg_MPa, L_drum_m, W_steam_kg_s, W_fw_kg_s\n";
    out << std::fixed << std::setprecision(4);

    constexpr double dt        = 1.0e-3;       // 1 ms; making it bigger will cause instability
    constexpr double t_total   = 300.0;
    constexpr double t_step    = 30.0;
    constexpr double log_every = 1.0;          // log once per simulated second

    bool stepped_down  = false;
    double next_log_t  = 0.0;

    // For a quick summary at the end.
    double P0 = loop.reactor().state().kinetics.power();
    double Psg0 = loop.steamGenerator().state().P;
    double Ppzr0 = loop.pressurizer().state().pressure_Pa;

    const long n_steps = static_cast<long>(t_total / dt);
    for (long i = 0; i <= n_steps; ++i) {
        const double t = i * dt;

        if (!stepped_down && t >= t_step) {
            loop.setReactorController(rebuildReactorControllerAtLoad(0.80));
            stepped_down = true;
            std::cout << "  t = " << t << " s: turbine load -> 80%\n";
        }

        if (t + 0.5 * dt >= next_log_t) {
            const auto& rs  = loop.reactor().state();
            const auto& sg  = loop.steamGenerator();
            const auto& pz  = loop.pressurizer();

            const double T_avg = 0.5 * (rs.thermal.T_hot_leg + rs.thermal.T_cold_leg);
            const double P_MW  = rs.kinetics.power()
                                 * loop.reactor().ratedPowerW() * 1.0e-6;

            out << t << ", "
                << P_MW << ", "
                << T_avg << ", "
                << rs.thermal.T_hot_leg << ", "
                << rs.thermal.T_cold_leg << ", "
                << loop.reactor().reactivity().rho_external * 1.0e5 << ", "
                << loop.pump().massFlow_kg_s() << ", "
                << pz.state().pressure_Pa * 1.0e-6 << ", "
                << pz.state().water_level_m << ", "
                << pz.inputs().heater_power_W * 1.0e-3 << ", "
                << pz.inputs().spray_mass_flow_kg_s << ", "
                << sg.state().P * 1.0e-6 << ", "
                << sg.state().L_dw << ", "
                << sg.steamMassFlow_kg_s() << ", "
                << sg.inputs().feedwater_mass_flow_kg_s << "\n";

            next_log_t += log_every;
        }

        loop.timeStep(dt);
    }

    const auto& rs  = loop.reactor().state();
    const auto& sg  = loop.steamGenerator();
    const auto& pz  = loop.pressurizer();
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "\nFinal state at t = " << t_total << " s:\n";
    std::cout << "  Reactor power       : "
              << rs.kinetics.power() * 100.0 << " % rated  "
              << "(was " << P0 * 100.0 << " %)\n";
    std::cout << "  Primary T_avg       : "
              << 0.5 * (rs.thermal.T_hot_leg + rs.thermal.T_cold_leg) << " K  "
              << "(setpoint at 80% = "
              << 562.0 + 0.80 * (583.0 - 562.0) << " K)\n";
    std::cout << "  External reactivity : "
              << loop.reactor().reactivity().rho_external * 1.0e5 << " pcm\n";
    std::cout << "  SG pressure         : "
              << sg.state().P * 1.0e-6 << " MPa  "
              << "(was " << Psg0 * 1.0e-6 << " MPa)\n";
    std::cout << "  SG drum level       : "
              << sg.state().L_dw << " m  (setpoint 3.00 m)\n";
    std::cout << "  Pressurizer P       : "
              << pz.state().pressure_Pa * 1.0e-6 << " MPa  "
              << "(was " << Ppzr0 * 1.0e-6 << " MPa, setpoint 15.50 MPa)\n";
    std::cout << "\nTrace written to primary_loop_load_change.csv\n";
    return 0;
}