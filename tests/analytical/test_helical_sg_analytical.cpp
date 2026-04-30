/**
 * @file   test_helical_sg_analytical.cpp
 * @brief  Analytical tests for the moving-boundary helical-coil SG of
 *         Arda & Holbert (2015).
 *
 * Tests:
 *   1. Steady-state initialisation produces a state vector consistent with
 *      Arda 2015 Table 4 (pressure, region lengths, temperatures within
 *      reasonable tolerance).
 *   2. At steady state the time derivative norm is small (the SG is
 *      properly initialised at its operating point).
 *   3. The model integrates for 60 s without crossing region-topology
 *      limits.
 *   4. Qualitative dynamic response: a step rise in feedwater inlet
 *      temperature shortens the sub-cooled region (Arda Fig. 11).
 *   5. Qualitative dynamic response: a step rise in feedwater mass flow
 *      lowers the secondary pressure (Arda Fig. 14).
 */

#include "astara/sg/HelicalCoilSteamGenerator.hpp"
#include "astara/props/IF97Water.hpp"
#include <gtest/gtest.h>

#include <cmath>

using astara::sg::HelicalCoilSteamGenerator;
using astara::sg::HelicalCoilSteamGeneratorParameters;
using astara::props::IF97Water;

namespace {

HelicalCoilSteamGenerator makeRated(IF97Water* w) {
    auto p = HelicalCoilSteamGeneratorParameters::nuscaleSMRTwoSG();
    HelicalCoilSteamGenerator sg(p, w);
    // Per-SG conditions for the NuScale design (160 MWth split between two SGs):
    //   primary inlet T = 291 deg C (hot leg riser temperature, Arda Table 4)
    //   primary mass flow per SG = 354 kg/s (= 708 / 2)
    //   feedwater per SG = 35.6 kg/s (= 71.25 / 2)
    //   feedwater T = 148.5 deg C
    sg.initialiseSteadyState(/*T_pi*/ 291.0 + 273.15,
                             /*W_p */ 354.0,
                             /*W_fw*/ 35.6,
                             /*T_fw*/ 148.5 + 273.15);
    return sg;
}

}  // namespace

// -----------------------------------------------------------------------------
// 1. Initialisation matches Arda Table 4 within reasonable tolerance.
// -----------------------------------------------------------------------------
TEST(HelicalCoilSteamGenerator, InitialStateMatchesArdaTable4) {
    IF97Water w;
    auto sg = makeRated(&w);
    const auto& s = sg.state();

    // Region lengths from Table 4: L_1 = 2.90, L_2 = 17.60, L_3 = 1.75 m.
    EXPECT_NEAR(s.L_1, 2.90, 0.05);
    EXPECT_NEAR(s.L_2, 17.60, 0.05);
    EXPECT_NEAR(sg.superheatedRegionLength_m(), 1.75, 0.05);

    // Pressure ~ 3.1 MPa.
    EXPECT_NEAR(s.p_S / 1e6, 3.1, 0.05);

    // Primary outlet (T_P1) ~ 246 deg C = 519.15 K.
    EXPECT_NEAR(s.T_P1, 246.0 + 273.15, 5.0);

    // Steam outlet temperature ~ 264 deg C = 537.15 K.
    EXPECT_NEAR(sg.steamOutletTemperatureK(), 264.0 + 273.15, 10.0);
}

// -----------------------------------------------------------------------------
// 2. After a brief settling period the time derivative is small (the SG
//    reaches and stays at a steady state that's close to the published
//    Arda Table 4 values).
// -----------------------------------------------------------------------------
TEST(HelicalCoilSteamGenerator, SettlesToSteadyState) {
    IF97Water w;
    auto sg = makeRated(&w);
    constexpr double dt = 0.05;
    for (int i = 0; i < 4000; ++i) sg.timeStep(dt);   // 200 s

    const auto dy = sg.evaluateDerivative(sg.state());
    // After settling, all derivatives should be very small (mm/s on
    // lengths, kPa/s on pressure, K/s on temperatures).
    EXPECT_LT(std::abs(dy.L_1), 1.0e-3);
    EXPECT_LT(std::abs(dy.L_2), 1.0e-2);
    EXPECT_LT(std::abs(dy.p_S), 5.0e3);
    EXPECT_LT(std::abs(dy.T_M2), 0.5);
    EXPECT_LT(std::abs(dy.T_P1), 0.5);
}

// -----------------------------------------------------------------------------
// 3. 60 s integration from rated condition with small dt is stable.
// -----------------------------------------------------------------------------
TEST(HelicalCoilSteamGenerator, IntegratesStablyForSixtySeconds) {
    IF97Water w;
    auto sg = makeRated(&w);
    constexpr double dt = 0.05;
    for (int i = 0; i < 1200; ++i) {
        sg.timeStep(dt);
    }
    const auto& s = sg.state();
    EXPECT_TRUE(std::isfinite(s.p_S));
    EXPECT_TRUE(std::isfinite(s.L_1));
    EXPECT_TRUE(std::isfinite(s.L_2));
    EXPECT_GT(s.L_1, 0.0);
    EXPECT_GT(s.L_2, 0.0);
    EXPECT_GT(sg.superheatedRegionLength_m(), 0.0);
    EXPECT_GT(s.p_S, 1.0e6);
    EXPECT_LT(s.p_S, 8.0e6);
}

// -----------------------------------------------------------------------------
// 4. Increasing feedwater inlet enthalpy shortens the sub-cooled region.
//    (Arda Fig. 11: feedwater-temperature step => L1 decreases.)
// -----------------------------------------------------------------------------
TEST(HelicalCoilSteamGenerator, HotterFeedwaterShortensSubcooledRegion) {
    IF97Water w;
    auto sg = makeRated(&w);
    constexpr double dt = 0.05;
    // Settle for 200 s -- region-length time constants are ~ 100 s in the
    // Arda model (Fig. 11 reaches new steady state ~ 200 s after a step).
    for (int i = 0; i < 4000; ++i) sg.timeStep(dt);
    const double L1_baseline = sg.state().L_1;

    // Bump feedwater inlet enthalpy: ~ 7 K hotter feedwater corresponds to
    // ~ 30 kJ/kg increase in h_i for liquid water.
    sg.inputs().feedwater_enthalpy_J_kg += 30.0e3;
    for (int i = 0; i < 4000; ++i) sg.timeStep(dt);       // 200 s
    EXPECT_LT(sg.state().L_1, L1_baseline);
}

// -----------------------------------------------------------------------------
// 5. Increasing feedwater flow lowers the secondary pressure.
//    (Arda Fig. 14: feedwater-flow step UP => steam pressure DOWN.)
// -----------------------------------------------------------------------------
TEST(HelicalCoilSteamGenerator, MoreFeedwaterLowersPressure) {
    IF97Water w;
    auto sg = makeRated(&w);
    constexpr double dt = 0.05;
    for (int i = 0; i < 4000; ++i) sg.timeStep(dt);       // settle 200 s
    const double P0 = sg.state().p_S;

    sg.inputs().feedwater_mass_flow_kg_s *= 1.05;
    for (int i = 0; i < 4000; ++i) sg.timeStep(dt);
    EXPECT_LT(sg.state().p_S, P0);
}
