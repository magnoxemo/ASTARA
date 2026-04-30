/**
 * @file   test_controllers.cpp
 * @brief  Unit tests for PID and the higher-level controller classes.
 */

#include "astara/control/PIDController.hpp"
#include "astara/control/ThreeElementController.hpp"
#include "astara/control/ReactorController.hpp"
#include "astara/control/PressurizerController.hpp"

#include <gtest/gtest.h>

using astara::control::PIDController;
using astara::control::PIDConfig;
using astara::control::ThreeElementController;
using astara::control::ReactorController;
using astara::control::PressurizerController;

// -----------------------------------------------------------------------------
// PID
// -----------------------------------------------------------------------------

TEST(PID, PureProportionalProducesProportionalOutput) {
    PIDController pid({2.0, 0.0, 0.0, -10.0, 10.0});
    pid.setSetpoint(5.0);
    EXPECT_DOUBLE_EQ(pid.update(/*meas=*/3.0, /*dt=*/0.1), 4.0);  // 2*(5-3)=4
}

TEST(PID, IntegratorAccumulates) {
    PIDController pid({0.0, 1.0, 0.0, -10.0, 10.0});
    pid.setSetpoint(1.0);
    // With a constant 1.0 error, integrator grows by error*dt each step.
    pid.update(0.0, 0.1);
    pid.update(0.0, 0.1);
    EXPECT_NEAR(pid.integrator(), 0.2, 1.0e-12);
}

TEST(PID, OutputIsClampedToLimits) {
    PIDController pid({100.0, 0.0, 0.0, -1.0, 1.0});
    pid.setSetpoint(10.0);
    EXPECT_DOUBLE_EQ(pid.update(0.0, 0.1), 1.0);   // saturates at u_max
}

TEST(PID, AntiWindupHoldsIntegratorWhenSaturated) {
    PIDController pid({100.0, 1.0, 0.0, -1.0, 1.0});
    pid.setSetpoint(10.0);
    pid.update(0.0, 1.0);   // big error, output saturated
    pid.update(0.0, 1.0);
    pid.update(0.0, 1.0);
    // Integrator should not have grown -- anti-windup held it.
    EXPECT_NEAR(pid.integrator(), 0.0, 1.0e-12);
}

TEST(PID, ResetClearsState) {
    PIDController pid({1.0, 1.0, 0.0, -10.0, 10.0});
    pid.setSetpoint(5.0);
    pid.update(0.0, 0.1);
    pid.update(0.0, 0.1);
    pid.reset();
    EXPECT_DOUBLE_EQ(pid.integrator(), 0.0);
    EXPECT_DOUBLE_EQ(pid.lastError(),  0.0);
}

// -----------------------------------------------------------------------------
// ThreeElementController
// -----------------------------------------------------------------------------

TEST(ThreeElementController, FeedforwardEqualsSteamFlowAtZeroError) {
    ThreeElementController c(/*Lsp=*/3.0, {5.0, 0.5, 0.0, -50, 50});
    const double demand = c.update(/*L=*/3.0, /*Wsteam=*/100.0, /*dt=*/0.1);
    EXPECT_NEAR(demand, 100.0, 1.0e-12);  // pure feedforward when level error = 0
}

TEST(ThreeElementController, LowLevelDemandsMoreFeedwater) {
    ThreeElementController c(/*Lsp=*/3.0, {5.0, 0.0, 0.0, -50, 50});
    const double demand = c.update(/*L=*/2.5, /*Wsteam=*/100.0, /*dt=*/0.1);
    EXPECT_GT(demand, 100.0);
}

// -----------------------------------------------------------------------------
// ReactorController
// -----------------------------------------------------------------------------

TEST(ReactorController, BelowSetpointInsertsPositiveReactivity) {
    ReactorController::Config cfg;
    cfg.pid = {1.0e-5, 0.0, 0.0, -1, 1};
    cfg.T_avg_setpoint_full_load_K = 583.0;
    cfg.T_avg_setpoint_no_load_K   = 562.0;
    ReactorController c(cfg);
    c.setTurbineLoadFraction(1.0);
    // Measured T_avg below 583 K -> demand positive reactivity rate.
    c.update(575.0, 1.0);
    EXPECT_GT(c.externalReactivity(), 0.0);
}

TEST(ReactorController, ReactivityIsRateLimited) {
    ReactorController::Config cfg;
    cfg.pid = {1e3, 0.0, 0.0, -1.0, 1.0};   // huge gain
    cfg.max_reactivity_rate = 1.0e-5;                                 // 1 pcm/s
    ReactorController c(cfg);
    c.setTurbineLoadFraction(1.0);
    c.update(/*T_avg=*/500.0, /*dt=*/1.0);  // huge negative error
    // After 1 s rate-limited at 1 pcm/s, reactivity should be exactly +1 pcm.
    EXPECT_NEAR(c.externalReactivity(), 1.0e-5, 1.0e-12);
}

// -----------------------------------------------------------------------------
// PressurizerController
// -----------------------------------------------------------------------------

TEST(PressurizerController, HeaterAtSetpointEqualsBias) {
    PressurizerController::Config cfg;
    PressurizerController c(cfg);
    EXPECT_DOUBLE_EQ(c.heaterDemand(cfg.pressure_setpoint_Pa),
                     cfg.heater_steady_state_W);
}

TEST(PressurizerController, LowPressureRaisesHeater) {
    PressurizerController::Config cfg;
    PressurizerController c(cfg);
    EXPECT_GT(c.heaterDemand(cfg.pressure_setpoint_Pa - 1.0e5),
              cfg.heater_steady_state_W);
}

TEST(PressurizerController, SprayOpensAboveThreshold) {
    PressurizerController::Config cfg;
    PressurizerController c(cfg);
    EXPECT_DOUBLE_EQ(c.sprayDemand(cfg.pressure_setpoint_Pa), 0.0);
    const double P_high = cfg.pressure_setpoint_Pa + cfg.spray_threshold_Pa + 1e4;
    EXPECT_GT(c.sprayDemand(P_high), 0.0);
    EXPECT_TRUE(c.spraysOpen());
}

TEST(PressurizerController, SprayDeadBandKeepsItOpen) {
    PressurizerController::Config cfg;
    PressurizerController c(cfg);
    // Open it.
    c.sprayDemand(cfg.pressure_setpoint_Pa + cfg.spray_threshold_Pa + 5e4);
    EXPECT_TRUE(c.spraysOpen());
    // Drop pressure inside the dead-band: should stay open.
    const double P_in_band = cfg.pressure_setpoint_Pa + cfg.spray_threshold_Pa
                           - 0.5 * cfg.spray_dead_band_Pa;
    c.sprayDemand(P_in_band);
    EXPECT_TRUE(c.spraysOpen());
    // Drop fully below threshold - dead-band: should close.
    c.sprayDemand(cfg.pressure_setpoint_Pa);
    EXPECT_FALSE(c.spraysOpen());
}
