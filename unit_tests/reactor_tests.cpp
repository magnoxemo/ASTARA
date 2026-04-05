// test_reactor.cpp
#include "Function.h"
#include "Reactor.h"
#include <gtest/gtest.h>
#include <memory>

using namespace astara;

class ReactorTest : public ::testing::Test {
protected:
  // From thesis Table A.1
  std::vector<double> decay_constants = {0.0125, 0.0308, 0.114,
                                         0.307,  1.19,   3.19};
  std::vector<double> beta_values = {0.000209, 0.001414, 0.001309,
                                     0.002727, 0.000925, 0.000314};
  double Lambda = 1.79e-5; // From thesis

  std::unique_ptr<Reactor> reactor;

  void SetUp() override {
    reactor =
        std::make_unique<Reactor>(6, decay_constants, beta_values, Lambda);

    // Set basic parameters
    reactor->setFuelMass(222739.0 * 0.453592); // Convert lbm to kg
    reactor->setModeratorMass(50000.0);
    reactor->setFuelThermalCapacity(0.059 *
                                    4186.8); // Convert BTU/lbm-°F to J/kg-K
    reactor->setModeratorThermalCapacity(1.39 * 4186.8);
    reactor->setHeatTransferArea(59900.0 * 0.092903); // Convert ft² to m²
    reactor->setRatedPower(3436.0);                   // MWt
    reactor->setCoolantFlowRate(1.5e8 * 0.000125998); // Convert lbm/hr to kg/s
    reactor->setCoolantInletTemperature(296.96);      // °C

    // Set plenum and leg masses
    reactor->setUpperPlenumMass(38.96 * 700.0); // kg
    reactor->setLowerPlenumMass(50.72 * 700.0); // kg
    reactor->setHotLegMass(28.32 * 700.0);      // kg
    reactor->setColdLegMass(56.63 * 700.0);     // kg

    reactor->setFractionPowerInFuel(0.974); // From thesis

    // Initialize precursor concentrations for steady state
    reactor->setInitialPower(1.0); // Fraction of full power
    for (unsigned int i = 0; i < 6; ++i) {
      double C_i = (beta_values[i] * reactor->getPower()) /
                   (decay_constants[i] * Lambda);
      reactor->setInitialPrecursorConcentration(i, C_i);
    }
  }

  void TearDown() override { reactor.reset(); }
};

TEST_F(ReactorTest, AllInitialConditionsIndependent) {
  reactor->setInitialPower(1.0);
  reactor->setInitialReactivity(0.0);
  reactor->setInitialFuelTemperature(565.0);
  reactor->setInitialModeratorTemperature(300.0);
  reactor->setInitialUpperPlenumTemperature(325.0);
  reactor->setInitialLowerPlenumTemperature(300.0);
  reactor->setInitialHotLegTemperature(325.0);
  reactor->setInitialColdLegTemperature(300.0);

  EXPECT_DOUBLE_EQ(reactor->getPower(), 1.0);
  EXPECT_DOUBLE_EQ(reactor->getReactivity(), 0.0);
  EXPECT_DOUBLE_EQ(reactor->getFuelTemperature(), 565.0);
  EXPECT_DOUBLE_EQ(reactor->getModeratorTemperature(), 300.0);
  EXPECT_DOUBLE_EQ(reactor->getUpperPlenumTemperature(), 325.0);
  EXPECT_DOUBLE_EQ(reactor->getLowerPlenumTemperature(), 300.0);
  EXPECT_DOUBLE_EQ(reactor->getHotLegTemperature(), 325.0);
  EXPECT_DOUBLE_EQ(reactor->getColdLegTemperature(), 300.0);
}

TEST_F(ReactorTest, NodeTemperaturesInitialized) {
  reactor->setInitialFuelTemperature(565.0);
  reactor->setInitialModeratorTemperature(300.0);

  const auto &fuel_temps = reactor->getFuelTemperatures();
  const auto &coolant_temps = reactor->getCoolantTemperatures();

  EXPECT_EQ(fuel_temps.size(), 3);
  EXPECT_EQ(coolant_temps.size(), 6);

  for (double t : fuel_temps) {
    EXPECT_DOUBLE_EQ(t, 565.0);
  }

  for (double t : coolant_temps) {
    EXPECT_DOUBLE_EQ(t, 300.0);
  }
}

TEST_F(ReactorTest, SetHeatTransferCoefficientWithConstantFunction) {
  Function *h_const = Reactor::createConstantFunction(5000.0);
  reactor->setHeatTransferCoefficient(h_const);

  reactor->setInitialPower(1.0);
  reactor->setInitialFuelTemperature(565.0);
  reactor->setInitialModeratorTemperature(300.0);

  EXPECT_NO_THROW(reactor->timeStep(0.01));
  EXPECT_GT(reactor->getPower(), 0.0);
}

TEST_F(ReactorTest, SetHeatTransferCoefficientWithExpressionFunction) {
  Function *h_func = new Function("5000 + 10 * (x - 300)", {"x"});
  reactor->setHeatTransferCoefficient(h_func);

  reactor->setInitialPower(1.0);
  reactor->setInitialFuelTemperature(565.0);
  reactor->setInitialModeratorTemperature(300.0);

  EXPECT_NO_THROW(reactor->timeStep(0.01));
  EXPECT_GT(reactor->getPower(), 0.0);
}

TEST_F(ReactorTest, NoFunctionsDefaultBehavior) {
  reactor->setInitialPower(1.0);
  reactor->setInitialReactivity(0.0);
  reactor->setInitialFuelTemperature(565.0);
  reactor->setInitialModeratorTemperature(300.0);

  EXPECT_NO_THROW(reactor->timeStep(0.01));
  EXPECT_GT(reactor->getPower(), 0.0);
}

TEST_F(ReactorTest, ConfigurablePhysicalParameters) {
  reactor->setFuelMass(2000.0);
  reactor->setModeratorMass(50000.0);
  reactor->setFuelThermalCapacity(500.0);
  reactor->setModeratorThermalCapacity(4186.0);
  reactor->setHeatTransferArea(500.0);
  reactor->setFissionEnergy(3.2e-11);
  reactor->setControlRodEffectiveness(-0.0001);
  reactor->setBoronEffectiveness(-1.0e-5);
  reactor->setRatedPower(3400.0);
  reactor->setCoolantFlowRate(17700.0);
  reactor->setCoolantInletTemperature(296.96);
  reactor->setUpperPlenumMass(30000.0);
  reactor->setLowerPlenumMass(35000.0);
  reactor->setHotLegMass(20000.0);
  reactor->setColdLegMass(40000.0);
  reactor->setFractionPowerInFuel(0.97);

  // Since these getters don't exist, we'll test through behavior
  reactor->setInitialPower(1.0);
  EXPECT_NO_THROW(reactor->timeStep(0.01));
  EXPECT_GT(reactor->getPower(), 0.0);
}

TEST_F(ReactorTest, ControlRodInsertion) {
  reactor->setControlRodEffectiveness(-0.0001);

  reactor->setInitialPower(1.0);
  reactor->setInitialReactivity(0.0);
  reactor->setInitialFuelTemperature(565.0);

  double rho_before = reactor->getReactivity();
  reactor->insertControlRod(0.5); // 50% insertion
  double rho_after = reactor->getReactivity();

  EXPECT_LT(rho_after, rho_before);
  EXPECT_DOUBLE_EQ(rho_after, rho_before + (-0.0001 * 0.5));
}

TEST_F(ReactorTest, BoronInjection) {
  reactor->setBoronEffectiveness(-1.0e-5);

  reactor->setInitialReactivity(0.0);

  double rho_before = reactor->getReactivity();
  reactor->injectBoron(1000.0);
  double rho_after = reactor->getReactivity();

  EXPECT_LT(rho_after, rho_before);
  EXPECT_DOUBLE_EQ(rho_after, rho_before + (-1.0e-5 * 1000.0));
}

TEST_F(ReactorTest, MultipleTimeSteps) {
  reactor->setInitialPower(1.0);
  reactor->setInitialFuelTemperature(565.0);
  reactor->setInitialModeratorTemperature(300.0);
  reactor->setInitialUpperPlenumTemperature(325.0);
  reactor->setInitialLowerPlenumTemperature(300.0);
  reactor->setInitialHotLegTemperature(325.0);
  reactor->setInitialColdLegTemperature(300.0);

  double dt = 0.01;
  for (int i = 0; i < 100; ++i) {
    EXPECT_NO_THROW(reactor->timeStep(dt));
  }

  EXPECT_NEAR(reactor->getState().time, 1.0, 1e-10);
  EXPECT_GT(reactor->getPower(), 0.0);

  // Verify temperature gradients developed
  const auto &fuel_temps = reactor->getFuelTemperatures();
  const auto &coolant_temps = reactor->getCoolantTemperatures();

  // Fuel nodes should have temperature gradient (higher at top)
  EXPECT_GT(fuel_temps[2], fuel_temps[0]);

  // Coolant nodes should show heating along flow path
  EXPECT_GT(coolant_temps[5], coolant_temps[0]);

  // Plenum temperatures should be reasonable
  EXPECT_GT(reactor->getUpperPlenumTemperature(),
            reactor->getLowerPlenumTemperature());
  EXPECT_GT(reactor->getHotLegTemperature(), reactor->getColdLegTemperature());
}

TEST_F(ReactorTest, TemperatureFeedback) {
  // Set reactivity coefficients
  Function *alpha_f_func = new Function("-1.98e-5 * (x - 565)", {"x"});
  Function *alpha_m_func = new Function("-3.6e-4 * (x - 300)", {"x"});
  reactor->setFuelTemperatureCoEfficientFunction(alpha_f_func);
  reactor->setModeratorTemperatureCoEfficientFunction(alpha_m_func);

  reactor->setInitialPower(1.0);
  reactor->setInitialReactivity(0.0);
  reactor->setInitialFuelTemperature(565.0);
  reactor->setInitialModeratorTemperature(300.0);

  // Run a few steps
  for (int i = 0; i < 10; ++i) {
    reactor->timeStep(0.01);
  }

  // Reactivity should be negative due to temperature feedback
  EXPECT_LT(reactor->getReactivity(), 0.0);
}

class ReactorParameterizedTest : public ReactorTest,
                                 public ::testing::WithParamInterface<double> {
};

TEST_P(ReactorParameterizedTest, VariousInitialPowers) {
  double initial_power = GetParam();

  reactor->setInitialPower(initial_power);
  reactor->setInitialFuelTemperature(565.0);
  reactor->setInitialModeratorTemperature(300.0);

  // Update precursor concentrations for new power level
  for (unsigned int i = 0; i < 6; ++i) {
    double C_i =
        (beta_values[i] * reactor->getPower()) / (decay_constants[i] * Lambda);
    reactor->setInitialPrecursorConcentration(i, C_i);
  }

  EXPECT_NO_THROW(reactor->timeStep(0.01));
  EXPECT_GE(reactor->getPower(), 0.0);
}

INSTANTIATE_TEST_SUITE_P(DifferentPowerLevels, ReactorParameterizedTest,
                         ::testing::Values(0.1, 0.5, 1.0, 1.5));

class ReactorEdgeCaseTest : public ReactorTest {};

TEST_F(ReactorEdgeCaseTest, VeryLowPower) {
  reactor->setInitialPower(1e-6);
  reactor->setInitialFuelTemperature(565.0);
  reactor->setInitialModeratorTemperature(300.0);

  // Update precursor concentrations
  for (unsigned int i = 0; i < 6; ++i) {
    double C_i =
        (beta_values[i] * reactor->getPower()) / (decay_constants[i] * Lambda);
    reactor->setInitialPrecursorConcentration(i, C_i);
  }

  EXPECT_NO_THROW(reactor->timeStep(0.01));
  EXPECT_GE(reactor->getPower(), 0.0);
}

TEST_F(ReactorEdgeCaseTest, HighInitialReactivity) {
  reactor->setInitialPower(1.0);
  reactor->setInitialReactivity(0.01); // 1000 pcm
  reactor->setInitialFuelTemperature(565.0);
  reactor->setInitialModeratorTemperature(300.0);

  EXPECT_NO_THROW(reactor->timeStep(0.01));
  EXPECT_GT(reactor->getPower(), 1.0); // Power should increase
}

TEST_F(ReactorEdgeCaseTest, LargeTemperatureGradient) {
  reactor->setInitialPower(1.0);
  reactor->setInitialFuelTemperature(800.0);      // Very hot fuel
  reactor->setInitialModeratorTemperature(250.0); // Cool moderator

  EXPECT_NO_THROW(reactor->timeStep(0.01));
}

TEST_F(ReactorEdgeCaseTest, InvalidControlRodLength) {
  reactor->setInitialReactivity(0.0);
  EXPECT_THROW(reactor->insertControlRod(-10.0), std::runtime_error);
}

TEST_F(ReactorEdgeCaseTest, InvalidBoronConcentration) {
  reactor->setInitialReactivity(0.0);
  EXPECT_THROW(reactor->injectBoron(-100.0), std::runtime_error);
}

TEST_F(ReactorEdgeCaseTest, NegativeTimeStep) {
  EXPECT_THROW(reactor->timeStep(-0.01), std::runtime_error);
}

TEST_F(ReactorEdgeCaseTest, ZeroTimeStep) {
  EXPECT_THROW(reactor->timeStep(0.0), std::runtime_error);
}

TEST_F(ReactorEdgeCaseTest, VeryLargeTimeStep) {
  reactor->setInitialPower(1.0);
  reactor->setInitialFuelTemperature(565.0);
  reactor->setInitialModeratorTemperature(300.0);

  // Should still work but might be inaccurate
  EXPECT_NO_THROW(reactor->timeStep(1.0));
}

TEST_F(ReactorEdgeCaseTest, TemperatureClamping) {
  reactor->setInitialPower(1.0);
  reactor->setInitialFuelTemperature(-300.0);      // Below absolute zero
  reactor->setInitialModeratorTemperature(2000.0); // Extremely hot

  EXPECT_NO_THROW(reactor->timeStep(0.01));
  EXPECT_GE(reactor->getFuelTemperature(), -273.15);
  EXPECT_LE(reactor->getModeratorTemperature(), 2000.0);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}