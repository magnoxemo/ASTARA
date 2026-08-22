/**
 * @file   test_control_rod_system.cpp
 * @brief  Unit tests for the manual control-rod reactivity-rate source.
 */

#include "astara/control/ControlRodSystem.hpp"

#include <gtest/gtest.h>

using astara::control::ControlRodSystem;
using astara::control::ControlRodSystemConfig;
using astara::control::RodMode;

namespace {
ControlRodSystemConfig makeConfig() {
    ControlRodSystemConfig cfg;
    cfg.insertion_speed  = 0.5;
    cfg.withdrawal_speed = 0.4;
    cfg.step_reactivity  = 0.002;
    return cfg;
}
}  // namespace

TEST(ControlRodSystem, InsertProducesNegativeReactivityRate) {
    ControlRodSystem rods(makeConfig());
    EXPECT_NEAR(rods.reactivityRate(RodMode::Insert), -0.5 * 0.002, 1.0e-12);
}

TEST(ControlRodSystem, WithdrawProducesPositiveReactivityRate) {
    ControlRodSystem rods(makeConfig());
    EXPECT_NEAR(rods.reactivityRate(RodMode::Withdraw), 0.4 * 0.002, 1.0e-12);
}

TEST(ControlRodSystem, HoldProducesZeroRate) {
    ControlRodSystem rods(makeConfig());
    EXPECT_DOUBLE_EQ(rods.reactivityRate(RodMode::Hold), 0.0);
}

TEST(ControlRodSystem, ConfigIsMutable) {
    ControlRodSystem rods(makeConfig());
    rods.config().insertion_speed = 1.0;
    EXPECT_NEAR(rods.reactivityRate(RodMode::Insert), -1.0 * 0.002, 1.0e-12);
}
