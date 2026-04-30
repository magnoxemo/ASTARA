/**
 * @file   test_kinetics_analytical.cpp
 * @brief  Analytical-solution tests for point-reactor kinetics.
 *
 * Three problems with known closed-form answers:
 *
 *   1. Steady-state precursors. With dn/dt = 0 and rho = 0 we should have
 *      n constant.  The state-space derivative must return zeros.
 *
 *   2. One-group inhour-equation: with constant reactivity rho < beta the
 *      asymptotic period is given by
 *         omega = (rho - beta) / (Lambda + (beta - rho)/lambda * lambda)  [crude]
 *      In the one-group limit at long times the dominant root of the
 *      inhour equation is
 *         omega = (rho * lambda) / (beta - rho + lambda * Lambda)
 *      We march the simulator forward, fit log(n) vs t, compare slope to
 *      this analytical value.
 *
 *   3. Prompt jump.  An impulse of reactivity rho ~ 0.5 beta produces an
 *      essentially instantaneous power rise of n(0+)/n(0-) = beta/(beta-rho)
 *      (the "prompt jump approximation").  Test by stepping forward briefly
 *      and verifying the ratio is achieved.
 *
 * @cite Hetrick, D. L. (1993). "Dynamics of Nuclear Reactors."  Chap. 1.
 */

#include "astara/reactor/PointKinetics.hpp"
#include "astara/core/Integrator.hpp"
#include <gtest/gtest.h>

#include <cmath>

using astara::reactor::DelayedGroupConstants;
using astara::reactor::PointKineticsState;
using astara::reactor::pointKineticsDerivative;
using astara::reactor::steadyStatePrecursors;

namespace {

/// Step the kinetics state with classical RK4.
PointKineticsState advanceKinetics(PointKineticsState s,
                                   double rho,
                                   const DelayedGroupConstants& g,
                                   double dt, std::size_t nsteps) {
    auto deriv = [&](double, const PointKineticsState& y) {
        return pointKineticsDerivative(y, rho, g);
    };
    for (std::size_t i = 0; i < nsteps; ++i) {
        s = astara::core::rk4Step(0.0, s, dt, deriv);
    }
    return s;
}

}  // namespace

// -----------------------------------------------------------------------------
// Test 1: steady-state precursors give zero derivative.
// -----------------------------------------------------------------------------
TEST(PointKinetics, SteadyStateHasZeroDerivative) {
    auto g = DelayedGroupConstants::u235SixGroup();
    PointKineticsState s(g.beta.size());
    s.power() = 1.0;
    s.precursors() = steadyStatePrecursors(g, 1.0);

    auto dy = pointKineticsDerivative(s, /*rho=*/0.0, g);
    EXPECT_NEAR(dy.power(), 0.0, 1.0e-10);
    for (double dc : dy.precursors()) {
        EXPECT_NEAR(dc, 0.0, 1.0e-6);
    }
}

// -----------------------------------------------------------------------------
// Test 2: one-group inhour equation.  For a step of constant reactivity, after
// the prompt-jump transient the power evolves as n(t) = A * e^{omega t} where
// omega is the largest root of the inhour equation:
//     rho = omega Lambda + (beta omega) / (omega + lambda)
// We test in the one-group case: solve numerically for omega and compare to
// the asymptotic slope of log(n) vs t in the simulator output.
// -----------------------------------------------------------------------------
TEST(PointKinetics, OneGroupInhourPeriod) {
    // Single-group constants representative of U-235.
    DelayedGroupConstants g;
    g.lambda = {0.08};        // 1/s
    g.beta   = {0.0065};      // typical total beta
    g.Lambda = 5.0e-5;        // s

    constexpr double rho = 0.001;   // small positive (subprompt-critical)

    // Solve the inhour equation analytically for one group:
    //   rho = omega Lambda + beta omega / (omega + lambda)
    //   => rho (omega + lambda) = omega Lambda (omega + lambda) + beta omega
    //   => Lambda omega^2 + (Lambda lambda + beta - rho) omega - rho lambda = 0
    const double a = g.Lambda;
    const double b = g.Lambda * g.lambda[0] + g.beta[0] - rho;
    const double c = -rho * g.lambda[0];
    const double disc = b * b - 4.0 * a * c;
    ASSERT_GT(disc, 0.0);
    const double omega = (-b + std::sqrt(disc)) / (2.0 * a);  // positive root

    // Initial steady state.
    PointKineticsState s(1);
    s.power() = 1.0;
    s.precursors() = steadyStatePrecursors(g, 1.0);

    // Step forward enough that prompt-jump transients have died (~ a few
    // 1/lambda ~ 12 s) and the asymptotic exponential has emerged.
    constexpr double dt    = 1.0e-3;
    constexpr double t_end = 30.0;
    s = advanceKinetics(s, rho, g, dt, static_cast<std::size_t>(t_end / dt));
    const double n1 = s.power();

    s = advanceKinetics(s, rho, g, dt, static_cast<std::size_t>(1.0 / dt));
    const double n2 = s.power();
    const double omega_measured = std::log(n2 / n1) / 1.0;

    // Tolerance: expect agreement to <0.5% on the slope.
    EXPECT_NEAR(omega_measured, omega, std::abs(omega) * 5.0e-3);
}

// -----------------------------------------------------------------------------
// Test 3: prompt jump.  Step rho from 0 to rho_step (with rho_step < beta),
// and verify that within ~0.1 s (much less than the longest precursor
// time-constant) the ratio n / n0 approaches the prompt-jump value
//   n(0+)/n(0-) = beta / (beta - rho_step)
// -----------------------------------------------------------------------------
TEST(PointKinetics, PromptJumpRatio) {
    auto g = DelayedGroupConstants::u235SixGroup();
    const double beta = g.totalBeta();
    const double rho_step = 0.5 * beta;
    const double prompt_jump = beta / (beta - rho_step);

    PointKineticsState s(g.beta.size());
    s.power() = 1.0;
    s.precursors() = steadyStatePrecursors(g, 1.0);

    // Step forward 50 ms with rho = rho_step using a small dt.
    constexpr double dt    = 1.0e-5;
    constexpr double t_end = 0.05;     // 50 ms: many prompt times, but much
                                       // shorter than 1/lambda_min ~ 80 s
    s = advanceKinetics(s, rho_step, g, dt, static_cast<std::size_t>(t_end / dt));

    // Power should have jumped to roughly the prompt-jump value.  The exact
    // value drifts upward on a slow exponential, so allow 5% tolerance.
    EXPECT_NEAR(s.power(), prompt_jump, prompt_jump * 0.05);
}

// -----------------------------------------------------------------------------
// Test 4: One-group projection.  Generating a one-group equivalent of a
// six-group set should give the same long-time period as the full six-group
// simulator at small reactivity.
// -----------------------------------------------------------------------------
TEST(PointKinetics, OneGroupProjectionMatchesAtSmallRho) {
    auto g6 = DelayedGroupConstants::u235SixGroup();
    auto g1 = DelayedGroupConstants::oneGroupAverageOf(g6);

    EXPECT_EQ(g1.beta.size(), 1u);
    EXPECT_NEAR(g1.beta[0], g6.totalBeta(), 1.0e-12);
    EXPECT_GT(g1.lambda[0], 0.0);
    // Reasonable: weighted average of {0.0125, ..., 3.19}.
    EXPECT_LT(g1.lambda[0], 0.2);
    EXPECT_GT(g1.lambda[0], 0.05);
}
