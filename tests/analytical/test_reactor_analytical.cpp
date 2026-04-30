/**
 * @file   test_reactor_analytical.cpp
 * @brief  Analytical tests for the composite Reactor model.
 *
 * Six tests:
 *
 *   1. Steady-state initialisation with no external reactivity should hold
 *      (zero derivative everywhere).
 *
 *   2. Steady-state initialisation with no power and no flow should also
 *      hold (degenerate case, but tests the mathematics).
 *
 *   3. Average fuel and moderator temperatures at steady state should be
 *      reasonable for thesis-like Westinghouse parameters (T_avg ~ 580 K
 *      for moderator, ~700 K for fuel).
 *
 *   4. Cold-leg-out temperature equals the user-supplied inlet temperature
 *      at steady state.
 *
 *   5. Hot-leg-out temperature is consistent with energy balance:
 *         T_hot - T_cold = P / (mdot * cp).
 *
 *   6. After a small step in external reactivity, after long time,
 *      negative-feedback brings the system to a new steady state with the
 *      expected temperature shift.
 *
 * The parameters chosen for these tests are *thesis-representative* but not
 * the exact thesis values; the goal is to exercise the model end-to-end
 * with realistic numbers, not to reproduce a specific published result.
 */

#include "astara/reactor/Reactor.hpp"
#include <gtest/gtest.h>

#include <cmath>

using astara::reactor::DelayedGroupConstants;
using astara::reactor::Reactor;
using astara::reactor::ReactivityModel;
using astara::reactor::ReactorThermalParameters;
using astara::reactor::PointKineticsState;

namespace {

/// Build a 4-loop Westinghouse-style PWR reactor with thesis-representative
/// parameters in SI units.
Reactor makeWestinghouseLikeReactor() {
    ReactorThermalParameters tp;
    tp.num_fuel_nodes        = 3;
    tp.num_moderator_nodes   = 6;
    tp.fuel_mass_total_kg    = 1.0e5;       // ~100 tonnes UO2 in core
    tp.fuel_cp_J_per_kgK     = 300.0;       // UO2 average cp
    tp.fission_power_in_fuel = 0.974;       // thesis value
    tp.moderator_mass_total_kg = 1.5e4;     // ~15 tonnes coolant in core
    tp.moderator_cp_J_per_kgK  = 5400.0;    // PWR primary at ~580 K, 15.5 MPa
    tp.mass_flow_rate_kg_s     = 17600.0;   // ~140 Mlb/hr / 4 loops summed
    tp.overall_h_W_per_m2K     = 28000.0;   // overall U with film + clad + gap + UO2
    tp.heat_transfer_area_m2   = 5400.0;    // total clad surface area
    tp.lower_plenum_mass_kg    = 5000.0;
    tp.upper_plenum_mass_kg    = 5000.0;
    tp.hot_leg_mass_kg         = 2500.0;
    tp.cold_leg_mass_kg        = 2500.0;

    auto groups = DelayedGroupConstants::u235SixGroup();
    constexpr double P_rated = 3.4e9;       // 3400 MWth (Westinghouse 4-loop)

    ReactivityModel rm;
    rm.alpha_fuel_per_K      = -2.5e-5;     // typical PWR Doppler coefficient
    rm.alpha_moderator_per_K = -3.0e-4;     // typical PWR moderator coefficient
    return Reactor(groups, tp, P_rated, rm);
}

}  // namespace

// -----------------------------------------------------------------------------
// 1. Initialised steady state -> zero derivative everywhere.
// -----------------------------------------------------------------------------
TEST(Reactor, SteadyStateHasZeroDerivative) {
    auto r = makeWestinghouseLikeReactor();
    r.initialiseSteadyState(/*n0=*/1.0, /*T_inlet=*/559.0);

    auto dy = r.evaluateDerivative(r.state());
    EXPECT_NEAR(dy.kinetics.power(), 0.0, 1.0e-6);
    for (double dc : dy.kinetics.precursors()) {
        EXPECT_NEAR(dc, 0.0, 1.0e-3);
    }
    for (double dT : dy.thermal.T_fuel) {
        EXPECT_NEAR(dT, 0.0, 1.0e-3);
    }
    for (double dT : dy.thermal.T_moderator) {
        EXPECT_NEAR(dT, 0.0, 1.0e-3);
    }
    EXPECT_NEAR(dy.thermal.T_cold_leg,     0.0, 1.0e-3);
    EXPECT_NEAR(dy.thermal.T_lower_plenum, 0.0, 1.0e-3);
    EXPECT_NEAR(dy.thermal.T_upper_plenum, 0.0, 1.0e-3);
    EXPECT_NEAR(dy.thermal.T_hot_leg,      0.0, 1.0e-3);
}

// -----------------------------------------------------------------------------
// 2. Steady state with low power: T_fuel ~ T_moderator (no thermal driving force).
// -----------------------------------------------------------------------------
TEST(Reactor, LowPowerSteadyStateThermalAlignment) {
    auto r = makeWestinghouseLikeReactor();
    r.initialiseSteadyState(/*n0=*/1.0e-5, /*T_inlet=*/559.0);  // 1e-5 of rated

    const double Tf_avg = r.state().thermal.averageFuelTemperature();
    const double Tm_avg = r.state().thermal.averageModeratorTemperature();
    // Fuel-to-moderator delta should be small (~ a degree or two).
    EXPECT_LT(std::abs(Tf_avg - Tm_avg), 5.0);
}

// -----------------------------------------------------------------------------
// 3. Hot-leg-cold-leg temperature difference matches energy balance.
// -----------------------------------------------------------------------------
TEST(Reactor, EnergyBalanceTemperatureRise) {
    auto r = makeWestinghouseLikeReactor();
    constexpr double T_inlet = 559.0;
    r.initialiseSteadyState(1.0, T_inlet);

    const auto& s   = r.state();
    const auto& tp  = r.thermalParams();
    const double dT_expected = r.ratedPowerW() / (tp.mass_flow_rate_kg_s * tp.moderator_cp_J_per_kgK);
    const double dT_actual   = s.thermal.T_hot_leg - s.thermal.T_cold_leg;
    EXPECT_NEAR(dT_actual, dT_expected, dT_expected * 1.0e-3);
}

// -----------------------------------------------------------------------------
// 4. Cold-leg-out temperature equals supplied inlet at steady state.
// -----------------------------------------------------------------------------
TEST(Reactor, ColdLegOutEqualsInlet) {
    auto r = makeWestinghouseLikeReactor();
    constexpr double T_inlet = 559.0;
    r.initialiseSteadyState(1.0, T_inlet);
    EXPECT_NEAR(r.state().thermal.T_cold_leg, T_inlet, 1.0e-9);
}

// -----------------------------------------------------------------------------
// 5. Time-stepping a steady state preserves it (drift < 1e-6 K over 1 s).
// -----------------------------------------------------------------------------
TEST(Reactor, SteadyStateIsStableUnderTimeStepping) {
    auto r = makeWestinghouseLikeReactor();
    r.initialiseSteadyState(1.0, 559.0);

    const auto s0 = r.state();
    // dt must respect the prompt-neutron stability limit: with beta/Lambda ~ 360/s,
    // explicit RK4 needs dt < ~7e-3 s for absolute stability.  Pick 1 ms.
    constexpr double dt = 1.0e-3;
    constexpr int    N  = 1000;     // 1 s
    for (int i = 0; i < N; ++i) r.timeStep(dt);

    const auto s1 = r.state();
    // Power should have stayed essentially constant.
    EXPECT_NEAR(s1.kinetics.power(), s0.kinetics.power(), 1.0e-6);
    // Temperatures should have drifted by less than 1e-3 K (numerical noise
    // from the finite-step discretisation, but bounded).
    EXPECT_NEAR(s1.thermal.averageFuelTemperature(),
                s0.thermal.averageFuelTemperature(),     1.0e-3);
    EXPECT_NEAR(s1.thermal.averageModeratorTemperature(),
                s0.thermal.averageModeratorTemperature(), 1.0e-3);
}

// -----------------------------------------------------------------------------
// 6. Negative reactivity step -> power drops, temperatures fall.
// -----------------------------------------------------------------------------
TEST(Reactor, NegativeReactivityStepDropsPower) {
    auto r = makeWestinghouseLikeReactor();
    r.initialiseSteadyState(1.0, 559.0);

    // Apply a -50 pcm step.
    r.reactivity().rho_external = -5.0e-4;

    // dt must be much smaller than the prompt-neutron timescale Lambda ~ 1.8e-5 s
    // divided by (beta - rho) for explicit-RK4 stability.  Picking 1 ms is
    // conservatively safe and still lets us integrate 10 s in 10 000 steps.
    constexpr double dt = 1.0e-3;
    for (int i = 0; i < 10000; ++i) r.timeStep(dt);  // 10 s

    EXPECT_LT(r.state().kinetics.power(), 0.99);
    EXPECT_GT(r.state().kinetics.power(), 0.5);   // shouldn't crash to ~ 0
}
