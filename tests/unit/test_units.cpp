/**
 * @file   test_units.cpp
 * @brief  Tests for unit-conversion helpers.
 *
 * Reference values are taken from NIST SP-811 and verified by hand to at
 * least 6 significant figures.
 */

#include "astara/core/Units.hpp"
#include <gtest/gtest.h>

namespace u = astara::units;

TEST(Units, FahrenheitConversionsAreInverses) {
    for (double T_F : {-40.0, 32.0, 70.0, 212.0, 547.0, 1000.0}) {
        const double round_trip = u::kelvinToFahrenheit(u::fahrenheitToKelvin(T_F));
        EXPECT_NEAR(round_trip, T_F, 1.0e-12);
    }
}

TEST(Units, FahrenheitToKelvinKeyValues) {
    // -40 F = 233.15 K (and -40 C, classic reference)
    EXPECT_NEAR(u::fahrenheitToKelvin(-40.0), 233.15, 1.0e-12);
    // Water freezing: 32 F = 273.15 K
    EXPECT_NEAR(u::fahrenheitToKelvin(32.0), 273.15, 1.0e-12);
    // Water boiling at 1 atm: 212 F = 373.15 K
    EXPECT_NEAR(u::fahrenheitToKelvin(212.0), 373.15, 1.0e-12);
}

TEST(Units, PsiToPaIsExact) {
    // 1 psi = 6894.757293168361 Pa (international pound, 1959).
    EXPECT_NEAR(1.0 * u::kPsiToPa, 6894.757293168361, 1.0e-9);
    // PWR primary: 2250 psia = 15.51 MPa
    EXPECT_NEAR(2250.0 * u::kPsiToPa, 1.5513e7, 1.0e3);
}

TEST(Units, BtuPerLbmFConversion) {
    // 1 Btu/(lbm*F) = 4186.8 J/(kg*K) (within 0.01%).
    EXPECT_NEAR(1.0 * u::kBtuPerLbmFToJperKgK, 4186.8, 1.0);
}

TEST(Units, LbmPerHrToKgS) {
    // 100 lbm/hr = 0.012600 kg/s
    EXPECT_NEAR(100.0 * u::kLbmPerHrToKgPerS, 0.0126, 1.0e-6);
    // PWR primary mass flow ~ 1.4e8 lbm/hr ~ 17630 kg/s
    EXPECT_NEAR(1.4e8 * u::kLbmPerHrToKgPerS, 17640.0, 10.0);
}
