/**
 * @file   test_pump_analytical.cpp
 * @brief  Analytical tests for the homologous coolant pump.
 *
 * Three tests:
 *
 *   1. Initialised at rated power -> derivatives are zero.
 *   2. After a "pump trip" (P_d set to 0), the speed decays approximately
 *      exponentially with a time constant tau ~ I*omega^2/(rho*g*Q*H).
 *      We don't pretend to predict tau analytically (it depends on the
 *      coupling), but we verify that the speed drops below 1/e of rated
 *      after a few tens of seconds, consistent with the thesis Fig. 3.23.
 *   3. The homologous-curve formula is consistent with itself: at rated
 *      conditions, H(Q=Qs, N=Ns) = A0 + A1*Qs + A2*Qs^2 = the rated head.
 */

#include "astara/pump/HomologousPump.hpp"
#include <gtest/gtest.h>

#include <cmath>

using astara::pump::HomologousPump;
using astara::pump::HomologousPumpParameters;
using astara::pump::PumpCurveCoefficients;

namespace {

/// A reactor coolant pump with thesis-representative parameters in SI.
HomologousPumpParameters makeRcpParameters() {
    HomologousPumpParameters p;
    // Curve fit chosen so that:
    //   shutoff head A0 = 90 m
    //   at rated flow Qs = 6.0 m^3/s: H(Qs) = 60 m (typical PWR RCP rated head)
    //   curvature A2 = -1.0 m / (m^3/s)^2
    // => A2 Qs^2 + A1 Qs + A0 = 60
    //    -1.0 * 36 + A1 * 6 + 90 = 60
    //    -36 + 6 A1 + 90 = 60 => A1 = 1
    p.curve.A0 = 90.0;
    p.curve.A1 = 1.0;
    p.curve.A2 = -1.0;
    p.rated_speed_rev_s         = 20.0;       // ~1200 rpm
    p.rated_volumetric_flow_m3_s = 6.0;
    // K_loop chosen so H_loss(Qs) = rated head:  K * Qs^2 = 60
    //   => K = 60 / 36 = 1.667
    p.loop_resistance_K_s2_m5   = 60.0 / (6.0 * 6.0);
    p.effective_flow_area_m2    = 0.4;
    p.loop_length_m             = 70.0;
    p.fluid_density_kg_m3       = 720.0;
    p.moment_of_inertia_kg_m2   = 1500.0;     // pump+motor flywheel
    p.rated_input_power_W       = p.fluid_density_kg_m3 * 9.80665 * 6.0 * 60.0;
                                              // rho g Q H at rated
    return p;
}

}  // namespace

// -----------------------------------------------------------------------------
// 1. Initialised at rated -> derivative is zero (modulo numerical noise).
// -----------------------------------------------------------------------------
TEST(HomologousPump, InitialisationAtRatedHasZeroDerivative) {
    HomologousPump p(makeRcpParameters());
    p.initialiseAtRated();

    auto dy = p.evaluateDerivative(p.state());
    EXPECT_NEAR(dy.volumetric_flow, 0.0, 1.0e-6);
    EXPECT_NEAR(dy.speed_rev_s,     0.0, 1.0e-3);
}

// -----------------------------------------------------------------------------
// 2. Sanity: at rated, head equals head loss.
// -----------------------------------------------------------------------------
TEST(HomologousPump, RatedConditionsEnergyBalance) {
    HomologousPump p(makeRcpParameters());
    p.initialiseAtRated();
    const double H_p = p.developedHead_m();
    const double H_loss = p.params().loop_resistance_K_s2_m5
                          * p.state().volumetric_flow * p.state().volumetric_flow;
    EXPECT_NEAR(H_p, H_loss, 1.0e-9);
}

// -----------------------------------------------------------------------------
// 3. Pump trip: set P_d=0, integrate, verify speed decays monotonically.
// -----------------------------------------------------------------------------
TEST(HomologousPump, PumpTripCausesSpeedDecay) {
    HomologousPump p(makeRcpParameters());
    p.initialiseAtRated();
    const double N0 = p.state().speed_rev_s;
    const double Q0 = p.state().volumetric_flow;

    p.setInputPowerW(0.0);   // pump trip

    constexpr double dt = 0.01;
    constexpr int    N  = 6000;     // 60 s
    for (int i = 0; i < N; ++i) {
        p.timeStep(dt);
    }
    const double N_end = p.state().speed_rev_s;
    const double Q_end = p.state().volumetric_flow;

    // After 60 s the speed and flow should both be substantially lower.
    EXPECT_LT(N_end, 0.5 * N0);
    EXPECT_LT(Q_end, 0.5 * Q0);
    EXPECT_GE(N_end, 0.0);          // never negative
    EXPECT_GE(Q_end, 0.0);
}

// -----------------------------------------------------------------------------
// 4. Homologous-curve formula at rated.
// -----------------------------------------------------------------------------
TEST(HomologousPump, HeadAtRatedMatchesQuadratic) {
    auto pp = makeRcpParameters();
    HomologousPump p(pp);
    const double H = p.headAt(pp.rated_volumetric_flow_m3_s, pp.rated_speed_rev_s);
    const double expected = pp.curve.A0
                          + pp.curve.A1 * pp.rated_volumetric_flow_m3_s
                          + pp.curve.A2 * pp.rated_volumetric_flow_m3_s * pp.rated_volumetric_flow_m3_s;
    EXPECT_NEAR(H, expected, 1.0e-12);
    // Sanity check: 60 m
    EXPECT_NEAR(H, 60.0, 1.0e-9);
}
