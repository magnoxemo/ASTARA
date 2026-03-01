#include "Reactor.h"
#include "Function.h"
#include <gtest/gtest.h>
#include <memory>

using namespace astara;

// ===== Test Fixture for Common Setup =====
class ReactorTest : public ::testing::Test {
protected:
    std::vector<double> decay_constants = {0.0127, 0.0317, 0.115,
                                           0.311,  1.40,   3.87};
    std::vector<double> beta_values = {0.000215, 0.001424, 0.001274,
                                       0.002568, 0.000748, 0.000273};
    double Lambda = 5.0e-5;

    std::unique_ptr<Reactor> reactor;

    void SetUp() override {
        reactor = std::make_unique<Reactor>(6, decay_constants, beta_values, Lambda);
    }

    void TearDown() override {
        reactor.reset();
    }
};

// ===== Tests =====

TEST_F(ReactorTest, AllInitialConditionsIndependent) {
// Set each initial condition independently
reactor->setInitialPower(100.0);
reactor->setInitialReactivity(0.0);
reactor->setInitialFuelTemperature(565.0);
reactor->setInitialModeratorTemperature(300.0);

// Verify they're all set correctly
EXPECT_DOUBLE_EQ(reactor->getPower(), 100.0);
EXPECT_DOUBLE_EQ(reactor->getReactivity(), 0.0);
EXPECT_DOUBLE_EQ(reactor->getFuelTemperature(), 565.0);
EXPECT_DOUBLE_EQ(reactor->getModeratorTemperature(), 300.0);
}

TEST_F(ReactorTest, SetHeatTransferCoefficientWithConstantFunction) {
// Use helper to create constant function
Function* h_const = Reactor::createConstantFunction(5000.0);
reactor->setHeatTransferCoefficient(h_const);

reactor->setInitialPower(100.0);
reactor->setInitialFuelTemperature(565.0);
reactor->setInitialModeratorTemperature(300.0);

// Should not segfault - execute one time step
EXPECT_NO_THROW(reactor->timeStep(0.01));

// Verify reactor still works
EXPECT_GT(reactor->getPower(), 0.0);
}

TEST_F(ReactorTest, SetHeatTransferCoefficientWithExpressionFunction) {
// Create temperature-dependent heat transfer coefficient
Function* h_func = new Function("5000 + 10 * (x - 300)", {"x"});
reactor->setHeatTransferCoefficient(h_func);

reactor->setInitialPower(100.0);
reactor->setInitialFuelTemperature(565.0);
reactor->setInitialModeratorTemperature(300.0);

// Should handle temperature-dependent function
EXPECT_NO_THROW(reactor->timeStep(0.01));
EXPECT_GT(reactor->getPower(), 0.0);
}

TEST_F(ReactorTest, NoFunctionsDefaultBehavior) {
// Don't set any functions - should use defaults
reactor->setInitialPower(100.0);
reactor->setInitialReactivity(0.0);
reactor->setInitialFuelTemperature(565.0);
reactor->setInitialModeratorTemperature(300.0);

// Should work without any functions
EXPECT_NO_THROW(reactor->timeStep(0.01));
EXPECT_GT(reactor->getPower(), 0.0);
}

TEST_F(ReactorTest, ConfigurablePhysicalParameters) {
// Set all physical parameters
reactor->setFuelMass(2000.0);
reactor->setModeratorMass(50000.0);
reactor->setFuelThermalCapacity(500.0);
reactor->setModeratorThermalCapacity(4186.0);
reactor->setHeatTransferArea(500.0);
reactor->setFissionEnergy(3.2e-11);
reactor->setControlRodEffectiveness(-0.0001);
reactor->setBoronEffectiveness(-1.0e-5);

// Verify getters work
EXPECT_DOUBLE_EQ(reactor->getFuelMass(), 2000.0);
EXPECT_DOUBLE_EQ(reactor->getHeatTransferArea(), 500.0);

reactor->setInitialPower(100.0);
EXPECT_NO_THROW(reactor->timeStep(0.01));
EXPECT_GT(reactor->getPower(), 0.0);
}

TEST_F(ReactorTest, ControlRodInsertion) {
reactor->setControlRodEffectiveness(-0.0001);

reactor->setInitialPower(100.0);
reactor->setInitialReactivity(0.0);
reactor->setInitialFuelTemperature(565.0);

double rho_before = reactor->getReactivity();
reactor->insertControlRod(50.0); // Insert 50 cm
double rho_after = reactor->getReactivity();

// Reactivity should decrease
EXPECT_LT(rho_after, rho_before);
EXPECT_DOUBLE_EQ(rho_after, rho_before + (-0.0001 * 50.0));
}

TEST_F(ReactorTest, BoronInjection) {
reactor->setBoronEffectiveness(-1.0e-5);

reactor->setInitialReactivity(0.0);

double rho_before = reactor->getReactivity();
reactor->injectBoron(1000.0); // Inject 1000 ppm
double rho_after = reactor->getReactivity();

// Reactivity should decrease
EXPECT_LT(rho_after, rho_before);
EXPECT_DOUBLE_EQ(rho_after, rho_before + (-1.0e-5 * 1000.0));
}

TEST_F(ReactorTest, MultipleTimeSteps) {
reactor->setInitialPower(100.0);
reactor->setInitialFuelTemperature(565.0);
reactor->setInitialModeratorTemperature(300.0);

double dt = 0.01;
for (int i = 0; i < 100; ++i) {
EXPECT_NO_THROW(reactor->timeStep(dt));
}

// Should complete 1 second of simulation without issues
EXPECT_NEAR(reactor->getState().time, 1.0, 1e-10);
EXPECT_GT(reactor->getPower(), 0.0);
}

// ===== Parametrized Test Example =====
class ReactorParameterizedTest : public ReactorTest,
                                 public ::testing::WithParamInterface<double> {
};

TEST_P(ReactorParameterizedTest, VariousInitialPowers) {
double initial_power = GetParam();

reactor->setInitialPower(initial_power);
reactor->setInitialFuelTemperature(565.0);
reactor->setInitialModeratorTemperature(300.0);

EXPECT_NO_THROW(reactor->timeStep(0.01));

// Power should remain positive
EXPECT_GE(reactor->getPower(), 0.0);
}

// Test with various power levels: 1 MW, 10 MW, 100 MW, 1000 MW
INSTANTIATE_TEST_SUITE_P(DifferentPowerLevels, ReactorParameterizedTest,
        ::testing::Values(1.0, 10.0, 100.0, 1000.0));

// ===== Edge Case Tests =====
class ReactorEdgeCaseTest : public ReactorTest {};

TEST_F(ReactorEdgeCaseTest, VeryLowPower) {
reactor->setInitialPower(1e-6); // Extremely low power
reactor->setInitialFuelTemperature(565.0);
reactor->setInitialModeratorTemperature(300.0);

// Should handle without numerical issues
EXPECT_NO_THROW(reactor->timeStep(0.01));
EXPECT_GE(reactor->getPower(), 0.0);
}

TEST_F(ReactorEdgeCaseTest, HighInitialReactivity) {
reactor->setInitialPower(100.0);
reactor->setInitialReactivity(0.01); // 1000 pcm - very high
reactor->setInitialFuelTemperature(565.0);
reactor->setInitialModeratorTemperature(300.0);

// Should handle without crashing
EXPECT_NO_THROW(reactor->timeStep(0.01));
}

TEST_F(ReactorEdgeCaseTest, LargeTemperatureGradient) {
reactor->setInitialPower(100.0);
reactor->setInitialFuelTemperature(800.0);      // Very hot fuel
reactor->setInitialModeratorTemperature(250.0); // Cool moderator

// Should handle large temperature differences
EXPECT_NO_THROW(reactor->timeStep(0.01));
}

TEST_F(ReactorEdgeCaseTest, InvalidControlRodLength) {
reactor->setInitialReactivity(0.0);

// Negative length should throw
EXPECT_THROW(reactor->insertControlRod(-10.0), std::runtime_error);
}

TEST_F(ReactorEdgeCaseTest, InvalidBoronConcentration) {
reactor->setInitialReactivity(0.0);

// Negative concentration should throw
EXPECT_THROW(reactor->injectBoron(-100.0), std::runtime_error);
}

// ===== Main function for running tests =====
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}