/**
 * @file   test_integrator.cpp
 * @brief  Verify RK4 and Euler against closed-form ODE solutions.
 *
 * Two test problems with known exact solutions:
 *
 *   1. Linear decay:  dy/dt = -k y,  y(0) = y0.   Exact: y(t) = y0 e^{-kt}.
 *   2. Harmonic osc.: dx/dt = v, dv/dt = -w^2 x.  Exact: x(t) = x0 cos(w t).
 *
 * The harmonic oscillator is a stress test for the integrator's handling of
 * vector-valued states; the linear decay tests scalar accuracy.
 */

#include "astara/core/Integrator.hpp"
#include "astara/core/State.hpp"
#include <gtest/gtest.h>

#include <cmath>

using astara::core::FixedState;
using astara::core::eulerStep;
using astara::core::rk4Step;

TEST(RK4, LinearDecayHitsExactToHighPrecision) {
    // y' = -k y;  y(0) = 1
    constexpr double k    = 0.5;
    constexpr double t_end = 4.0;
    constexpr double dt    = 1.0e-2;
    constexpr int    N     = static_cast<int>(t_end / dt + 0.5);

    auto f = [k](double, double y) { return -k * y; };

    double y = 1.0;
    double t = 0.0;
    for (int i = 0; i < N; ++i) {
        y = rk4Step(t, y, dt, f);
        t += dt;
    }
    const double exact = std::exp(-k * t_end);
    EXPECT_NEAR(y, exact, 1.0e-9);
}

TEST(Euler, LinearDecayConvergesAtFirstOrder) {
    constexpr double k     = 0.5;
    constexpr double t_end = 4.0;
    auto f = [k](double, double y) { return -k * y; };

    auto err = [&](double dt) {
        const int N = static_cast<int>(t_end / dt + 0.5);
        double y = 1.0;
        double t = 0.0;
        for (int i = 0; i < N; ++i) {
            y = eulerStep(t, y, dt, f);
            t += dt;
        }
        return std::abs(y - std::exp(-k * t_end));
    };
    const double e1 = err(1.0e-2);
    const double e2 = err(5.0e-3);
    // Halving dt should roughly halve the error for a first-order scheme.
    EXPECT_LT(e2, 0.55 * e1);
    EXPECT_GT(e2, 0.45 * e1);
}

TEST(RK4, HarmonicOscillatorPreservesAmplitude) {
    // x' = v;  v' = -w^2 x.  Exact period T = 2*pi/w.
    constexpr double w     = 2.0;
    constexpr double dt    = 1.0e-3;

    using S = FixedState<2>;
    auto f = [w](double, const S& y) {
        return S{{ y[1], -w * w * y[0] }};
    };

    constexpr int    N  = static_cast<int>(2.0 * M_PI / w / dt + 0.5);
    constexpr double dt_exact = (2.0 * M_PI / w) / static_cast<double>(N);
    S y{{1.0, 0.0}};
    double t = 0.0;
    for (int i = 0; i < N; ++i) {
        y = rk4Step(t, y, dt_exact, f);
        t += dt_exact;
    }
    // After exactly one period, position should return to 1.0 and velocity to 0.
    EXPECT_NEAR(y[0], 1.0, 1.0e-6);
    EXPECT_NEAR(y[1], 0.0, 1.0e-6);
}
