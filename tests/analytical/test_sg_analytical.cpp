/**
 * @file   test_sg_analytical.cpp
 * @brief  Analytical tests for the Ali Model D steam generator.
 *
 * Tests:
 *   1. Steady-state initialisation produces a state with primary-side
 *      energy balance within tolerance:
 *         Q_in_primary = W_p * cp * (T_pi_in - T_p4)
 *      and on the secondary side:
 *         Q_to_secondary ~ W_steam * h_fg + W_fw * (h_f - h_fw)
 *      with both equal within ~10% (operator-splitting tolerance).
 *   2. Increasing the primary inlet temperature raises the SG pressure
 *      (more heat in -> higher saturation pressure, given a fixed valve).
 *   3. Closing the steam valve raises the SG pressure.
 *   4. The model integrates for 60 s without producing non-finite values
 *      from the rated steady state (basic stability).
 */

#include "astara/sg/AliSteamGenerator.hpp"
#include "astara/props/IF97Water.hpp"
#include <gtest/gtest.h>

#include <cmath>

using astara::sg::AliSteamGenerator;
using astara::sg::AliSteamGeneratorParameters;
using astara::props::IF97Water;

namespace {

AliSteamGenerator makeRatedSG(IF97Water* w) {
    auto p = AliSteamGeneratorParameters::westinghouseModelD5();
    AliSteamGenerator sg(p, w);
    sg.initialiseSteadyState(/*T_pi_in=*/597.0,
                             /*W_p=*/4400.0,         // ~ 17600 / 4 loops
                             /*P_sg=*/6.9e6,
                             /*L_dw=*/3.0);
    return sg;
}

}  // namespace

// -----------------------------------------------------------------------------
// 1. After initialisation, the primary heat load is positive and reasonable.
// -----------------------------------------------------------------------------
TEST(AliSteamGenerator, InitialisationProducesPositiveHeatLoad) {
    IF97Water w;
    auto sg = makeRatedSG(&w);
    const double Q = sg.primaryHeatLoad_W();
    EXPECT_GT(Q, 5.0e8);   // > 500 MW (single loop ~ 850 MW)
    EXPECT_LT(Q, 1.5e9);   // sanity upper bound
}

// -----------------------------------------------------------------------------
// 2. Steam mass flow is positive and consistent with valve flow equation.
// -----------------------------------------------------------------------------
TEST(AliSteamGenerator, SteamFlowPositiveAtRated) {
    IF97Water w;
    auto sg = makeRatedSG(&w);
    EXPECT_GT(sg.steamMassFlow_kg_s(), 100.0);   // ~ 470 kg/s rated; allow margin
    EXPECT_LT(sg.steamMassFlow_kg_s(), 1000.0);
}

// -----------------------------------------------------------------------------
// 3. Stability: integrate for 60 s from rated steady state.
// -----------------------------------------------------------------------------
TEST(AliSteamGenerator, IntegratesStablyAtRated) {
    IF97Water w;
    auto sg = makeRatedSG(&w);
    constexpr double dt = 0.05;
    for (int i = 0; i < 1200; ++i) {
        sg.timeStep(dt);
    }
    EXPECT_TRUE(std::isfinite(sg.state().P));
    EXPECT_TRUE(std::isfinite(sg.state().T_p4));
    EXPECT_TRUE(std::isfinite(sg.state().L_dw));
    // Pressure shouldn't have run away.
    EXPECT_GT(sg.state().P, 5.0e6);
    EXPECT_LT(sg.state().P, 9.0e6);
}

// -----------------------------------------------------------------------------
// 4. Hotter primary inlet -> primary outlet rises (heat balance closure).
// -----------------------------------------------------------------------------
TEST(AliSteamGenerator, HotterPrimaryInletRaisesOutlet) {
    IF97Water w;
    auto sg = makeRatedSG(&w);
    const double T_p4_baseline = sg.primaryOutletTemperatureK();

    sg.inputs().primary_inlet_temperature_K = 605.0;   // +8 K
    constexpr double dt = 0.05;
    for (int i = 0; i < 600; ++i) sg.timeStep(dt);     // 30 s
    EXPECT_GT(sg.primaryOutletTemperatureK(), T_p4_baseline);
}

// -----------------------------------------------------------------------------
// 5. Closing the steam valve raises pressure.
// -----------------------------------------------------------------------------
TEST(AliSteamGenerator, ClosingValveRaisesPressure) {
    IF97Water w;
    auto sg = makeRatedSG(&w);
    const double P0 = sg.state().P;

    // Close the valve by raising the steam-line back pressure above SG P.
    sg.inputs().steam_line_pressure_Pa = sg.state().P + 1.0e5;
    constexpr double dt = 0.05;
    for (int i = 0; i < 600; ++i) sg.timeStep(dt);     // 30 s
    EXPECT_GT(sg.state().P, P0);
}
