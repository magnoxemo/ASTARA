/**
 * @file   test_pressurizer_analytical.cpp
 * @brief  Analytical tests for the pressurizer model.
 *
 * Five tests:
 *
 *   1. Initialisation at saturation pressure with zero inputs -> dP/dt and
 *      dLw/dt are both zero (steady state).
 *
 *   2. Heater on, no surge, no spray -> pressure rises monotonically.
 *
 *   3. Spray on, heater off, no surge -> pressure falls monotonically.
 *
 *   4. Surge in (water entering) at saturation enthalpy -> water level
 *      rises monotonically; pressure may rise slightly due to compression.
 *
 *   5. Surge out (water leaving) -> water level falls monotonically.
 */

#include "astara/pressurizer/Pressurizer.hpp"
#include "astara/props/IF97Water.hpp"
#include <gtest/gtest.h>

using astara::pressurizer::Pressurizer;
using astara::pressurizer::PressurizerInputs;
using astara::pressurizer::PressurizerParameters;
using astara::props::IF97Water;

namespace {

PressurizerParameters makePressurizerParams() {
    PressurizerParameters p;
    p.cross_section_area_m2 = 4.0;       // ~3 m diameter cylinder
    p.total_height_m        = 13.0;      // typical PWR pressurizer
    return p;
}

}  // namespace

// -----------------------------------------------------------------------------
// 1. Steady state holds.
// -----------------------------------------------------------------------------
TEST(Pressurizer, SteadyStateHoldsWithZeroInputs) {
    IF97Water w;
    Pressurizer pz(makePressurizerParams(), &w);
    pz.initialiseSteadyState(15.5e6, /*Lw=*/8.0);
    pz.inputs() = {};   // everything zero

    auto dy = pz.evaluateDerivative(pz.state());
    // Pressure: derivative might be slightly non-zero because the surge_enthalpy
    // and spray_enthalpy default to h_f, and W_sg = W_sp = 0 means those terms
    // vanish.  So dP/dt should be exactly zero (Q=0, NUM=0).
    EXPECT_NEAR(dy.pressure_Pa,   0.0, 1.0e-3);
    EXPECT_NEAR(dy.water_level_m, 0.0, 1.0e-9);
}

// -----------------------------------------------------------------------------
// 2. Heater on, isolated (no surge, no spray) -> pressure rises.
// -----------------------------------------------------------------------------
TEST(Pressurizer, HeaterRaisesPressure) {
    IF97Water w;
    Pressurizer pz(makePressurizerParams(), &w);
    pz.initialiseSteadyState(15.5e6, 8.0);
    pz.inputs().heater_power_W = 1.4e6;   // 1.4 MW (typical PWR proportional heater bank)

    const double P0 = pz.state().pressure_Pa;
    constexpr double dt = 0.1;
    for (int i = 0; i < 200; ++i) pz.timeStep(dt);   // 20 s
    EXPECT_GT(pz.state().pressure_Pa, P0);
    // After 20 s the rise should be measurable but not crazy (< 5 bar).
    EXPECT_LT(pz.state().pressure_Pa - P0, 5.0e5);
}

// -----------------------------------------------------------------------------
// 3. Spray on, isolated -> pressure falls.
// -----------------------------------------------------------------------------
TEST(Pressurizer, SprayLowersPressure) {
    IF97Water w;
    Pressurizer pz(makePressurizerParams(), &w);
    pz.initialiseSteadyState(15.5e6, 8.0);
    pz.inputs().spray_mass_flow_kg_s = 5.0;
    // Cold spray water from cold leg, h ~ 1300 kJ/kg at 290 C.
    pz.inputs().spray_enthalpy_J_kg  = 1.3e6;

    const double P0 = pz.state().pressure_Pa;
    constexpr double dt = 0.1;
    for (int i = 0; i < 200; ++i) pz.timeStep(dt);   // 20 s
    EXPECT_LT(pz.state().pressure_Pa, P0);
}

// -----------------------------------------------------------------------------
// 4. In-surge raises water level.
// -----------------------------------------------------------------------------
TEST(Pressurizer, InSurgeRaisesLevel) {
    IF97Water w;
    Pressurizer pz(makePressurizerParams(), &w);
    pz.initialiseSteadyState(15.5e6, 8.0);
    pz.inputs().surge_mass_flow_kg_s = 20.0;          // 20 kg/s into pressurizer
    pz.inputs().surge_enthalpy_J_kg  = w.satLiquidEnthalpy_P(15.5e6);  // hot leg ~ saturated

    const double Lw0 = pz.state().water_level_m;
    constexpr double dt = 0.1;
    for (int i = 0; i < 100; ++i) pz.timeStep(dt);   // 10 s
    EXPECT_GT(pz.state().water_level_m, Lw0);
}

// -----------------------------------------------------------------------------
// 5. Out-surge lowers water level.
// -----------------------------------------------------------------------------
TEST(Pressurizer, OutSurgeLowersLevel) {
    IF97Water w;
    Pressurizer pz(makePressurizerParams(), &w);
    pz.initialiseSteadyState(15.5e6, 8.0);
    pz.inputs().surge_mass_flow_kg_s = -20.0;         // 20 kg/s out
    pz.inputs().surge_enthalpy_J_kg  = 0.0;            // unused on out-surge

    const double Lw0 = pz.state().water_level_m;
    constexpr double dt = 0.1;
    for (int i = 0; i < 100; ++i) pz.timeStep(dt);   // 10 s
    EXPECT_LT(pz.state().water_level_m, Lw0);
}
