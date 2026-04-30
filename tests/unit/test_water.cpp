/**
 * @file   test_water.cpp
 * @brief  Validate water-property implementations against published values.
 *
 * Reference values are taken from the NIST/ASME steam tables (1997
 * formulation) at PWR-relevant operating points.  Tolerances are tight for
 * IF97 (it *is* the reference) and looser for LinearizedWater (the linear
 * fit only matches around its centre point).
 */

#include "astara/props/IF97Water.hpp"
#include "astara/props/LinearizedWater.hpp"
#include <gtest/gtest.h>

#include <cmath>

using astara::props::IF97Water;
using astara::props::LinearizedWater;
using astara::props::PropertyError;

// -----------------------------------------------------------------------------
// IF97 -- spot checks at PWR-relevant operating points.
// -----------------------------------------------------------------------------

TEST(IF97Water, SaturationAt6MPa) {
    IF97Water w;
    // ASME steam tables (rounded): Tsat(6 MPa) = 275.59 C = 548.74 K
    EXPECT_NEAR(w.saturationTemperature(6.0e6), 548.74, 0.05);
    // Saturation densities at 6 MPa: rho_l ~ 758.0, rho_g ~ 30.82 kg/m^3
    EXPECT_NEAR(w.satLiquidDensity_P(6.0e6), 758.0, 0.5);
    EXPECT_NEAR(w.satVapourDensity_P(6.0e6),  30.82, 0.05);
    // h_fg(6 MPa) ~ 1571 kJ/kg
    EXPECT_NEAR(w.latentHeat_P(6.0e6), 1571.0e3, 1.0e3);
}

TEST(IF97Water, SubcooledLiquidAtPWRPrimaryConditions) {
    IF97Water w;
    // 575 K at 15.5 MPa: subcooled liquid (well below Tsat~617.9 K).
    // IF97 reference value: rho ~ 722.7 kg/m^3.  Bracket loosely so the test
    // still has signal if IF97 is replaced or rebuilt.
    const double rho = w.density_TP(575.0, 15.5e6);
    EXPECT_GT(rho, 700.0);
    EXPECT_LT(rho, 750.0);
    const double cp = w.cp_TP(575.0, 15.5e6);
    EXPECT_GT(cp, 5000.0);
    EXPECT_LT(cp, 6500.0);
}

TEST(IF97Water, SaturationCurveRoundTrip) {
    IF97Water w;
    for (double P : {1.0e6, 5.0e6, 7.0e6, 10.0e6, 15.0e6}) {
        const double T = w.saturationTemperature(P);
        const double P_back = w.saturationPressure(T);
        EXPECT_NEAR(P_back, P, P * 1.0e-6);
    }
}

TEST(IF97Water, OutOfRangeThrowsInStrictMode) {
    IF97Water w(/*strict=*/true);
    EXPECT_THROW({ w.saturationTemperature(0.5); }, PropertyError);     // 0.5 Pa
    EXPECT_THROW({ w.saturationPressure(900.0); }, PropertyError);      // > Tcrit
    EXPECT_THROW({ w.density_TP(200.0, 1.0e6); }, PropertyError);       // T < Ttrip
}

// -----------------------------------------------------------------------------
// LinearizedWater -- coefficients fit around 7 MPa should agree with IF97 at
// 7 MPa to high precision and degrade with distance.
// -----------------------------------------------------------------------------

TEST(LinearizedWater, MatchesIF97AtCentrePoint) {
    LinearizedWater lin;
    IF97Water if97;
    constexpr double P0 = 7.0e6;
    // Central-difference fits over dP = 0.5 MPa pick up curvature, so the
    // linear value at the centre is offset by ~ (1/6) * d^2 phi/dP^2 * dP^2.
    // For Tsat that's ~0.15 K; for the densities and enthalpies it's much
    // smaller (rho_l, rho_g, h_f, h_g are nearly linear in P at PWR
    // pressures).  Loosen the temperature tolerance accordingly.
    EXPECT_NEAR(lin.saturationTemperature(P0), if97.saturationTemperature(P0), 0.2);
    EXPECT_NEAR(lin.satLiquidDensity_P(P0),    if97.satLiquidDensity_P(P0),    1.0);
    EXPECT_NEAR(lin.satVapourDensity_P(P0),    if97.satVapourDensity_P(P0),    0.05);
    EXPECT_NEAR(lin.satLiquidEnthalpy_P(P0),   if97.satLiquidEnthalpy_P(P0),   1.0e3);
    EXPECT_NEAR(lin.satVapourEnthalpy_P(P0),   if97.satVapourEnthalpy_P(P0),   1.0e3);
}

TEST(LinearizedWater, CustomFitAroundOperatingPoint) {
    IF97Water if97;
    constexpr double P0 = 5.0e6;
    auto lin = LinearizedWater::fitAround(if97, P0, 2.0e5);
    EXPECT_NEAR(lin.saturationTemperature(P0), if97.saturationTemperature(P0), 0.1);
    EXPECT_NEAR(lin.satLiquidDensity_P(P0),    if97.satLiquidDensity_P(P0),    1.0);
}

TEST(LinearizedWater, SatPressureRoundTrip) {
    LinearizedWater lin;
    for (double P : {6.5e6, 7.0e6, 7.5e6}) {
        const double T = lin.saturationTemperature(P);
        const double P_back = lin.saturationPressure(T);
        EXPECT_NEAR(P_back, P, P * 1.0e-12);  // exact within roundoff (linear)
    }
}
