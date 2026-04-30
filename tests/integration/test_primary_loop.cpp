/**
 * @file   test_primary_loop.cpp
 * @brief  Integration test for the full primary loop (reactor + SG + pump
 *         + pressurizer, optionally with controllers).
 */

#include "astara/primary/PrimaryLoop.hpp"
#include "astara/props/IF97Water.hpp"

#include <gtest/gtest.h>
#include <cmath>
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

namespace {

// IF97Water needs static lifetime; tests that build a primary loop must
// hand it a reference that outlives the loop.  We provide a function-local
// static for that.
IF97Water& sharedWater() {
    static IF97Water w;
    return w;
}

std::unique_ptr<Reactor> makeReactor() {
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
    ReactivityModel rm;
    rm.alpha_fuel_per_K      = -2.5e-5;
    rm.alpha_moderator_per_K = -3.0e-4;
    auto r = std::make_unique<Reactor>(groups, tp, /*P_rated*/ 3.4e9, rm);
    r->initialiseSteadyState(/*n0=*/1.0, /*T_inlet=*/559.0);
    return r;
}

std::unique_ptr<AliSteamGenerator> makeSG() {
    auto p = AliSteamGeneratorParameters::westinghouseModelD5();
    auto sg = std::make_unique<AliSteamGenerator>(p, &sharedWater());
    sg->initialiseSteadyState(/*T_pi_in*/ 597.0,
                              /*W_p   */ 4400.0,
                              /*P_sg  */ 6.9e6,
                              /*L_dw  */ 3.0);
    return sg;
}

std::unique_ptr<HomologousPump> makePump() {
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

PrimaryLoop makeFullLoop() {
    return PrimaryLoop(makeReactor(), makeSG(), makePump(), makePressurizer());
}

}  // namespace

// -----------------------------------------------------------------------------
// 1. Constructor accepts well-initialised components.
// -----------------------------------------------------------------------------
TEST(PrimaryLoop, ConstructsFromValidComponents) {
    auto loop = makeFullLoop();
    EXPECT_GT(loop.reactor().state().kinetics.power(), 0.0);
    EXPECT_GT(loop.steamGenerator().state().P, 0.0);
    EXPECT_GT(loop.pump().state().speed_rev_s, 0.0);
    EXPECT_GT(loop.pressurizer().state().pressure_Pa, 0.0);
}

// -----------------------------------------------------------------------------
// 2. NULL component is rejected.
// -----------------------------------------------------------------------------
TEST(PrimaryLoop, NullComponentThrows) {
    EXPECT_THROW({
        PrimaryLoop(nullptr, makeSG(), makePump(), makePressurizer());
    }, std::invalid_argument);
}

// -----------------------------------------------------------------------------
// 3. Time-stepping the loop is stable for at least 30 s.
// -----------------------------------------------------------------------------
TEST(PrimaryLoop, IntegratesStablyForThirtySeconds) {
    auto loop = makeFullLoop();
    constexpr double dt = 1.0e-3;
    for (int i = 0; i < 30000; ++i) {
        loop.timeStep(dt);
    }
    EXPECT_TRUE(std::isfinite(loop.reactor().state().kinetics.power()));
    EXPECT_TRUE(std::isfinite(loop.steamGenerator().state().P));
    EXPECT_TRUE(std::isfinite(loop.pump().state().speed_rev_s));
    EXPECT_TRUE(std::isfinite(loop.pressurizer().state().pressure_Pa));
    // Power should still be in a sensible range
    EXPECT_GT(loop.reactor().state().kinetics.power(), 0.5);
    EXPECT_LT(loop.reactor().state().kinetics.power(), 2.0);
}

// -----------------------------------------------------------------------------
// 4. Attaching all controllers and stepping yields a stable loop.
// -----------------------------------------------------------------------------
TEST(PrimaryLoop, ControlledLoopIsStable) {
    auto loop = makeFullLoop();

    ReactorController::Config rcfg;
    rcfg.pid.Kp = 1.0e-6;
    rcfg.pid.Ki = 1.0e-7;
    rcfg.pid.u_min = -1.0;
    rcfg.pid.u_max =  1.0;
    rcfg.T_avg_setpoint_no_load_K   = 562.0;
    rcfg.T_avg_setpoint_full_load_K = 583.0;
    rcfg.max_reactivity_rate = 1.4e-5;
    auto rctl = std::make_unique<ReactorController>(rcfg);
    rctl->setTurbineLoadFraction(1.0);
    loop.setReactorController(std::move(rctl));

    PressurizerController::Config pcfg;
    auto pctl = std::make_unique<PressurizerController>(pcfg);
    loop.setPressurizerController(std::move(pctl));

    EXPECT_TRUE(loop.hasReactorController());
    EXPECT_TRUE(loop.hasPressurizerController());

    constexpr double dt = 1.0e-3;
    for (int i = 0; i < 5000; ++i) loop.timeStep(dt);
    EXPECT_TRUE(std::isfinite(loop.reactor().state().kinetics.power()));
}
