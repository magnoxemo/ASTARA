// test_reactor.cpp
#include "Function.h"
#include "Reactor.h"
#include <gtest/gtest.h>
#include <memory>

using namespace astara;

class ReactorTest : public ::testing::Test {
protected:
  std::vector<double> decay_constants = {0.0127, 0.0317, 0.115,
                                         0.311,  1.40,   3.87};
  std::vector<double> beta_values = {0.000215, 0.001424, 0.001274,
                                     0.002568, 0.000748, 0.000273};
  double Lambda = 5.0e-5;

  std::unique_ptr<Reactor> reactor;

  void SetUp() override {
    reactor =
        std::make_unique<Reactor>(6, decay_constants, beta_values, Lambda);
    reactor->setFuelMass(2000.0);
    reactor->setModeratorMass(50000.0);
    reactor->setFuelThermalCapacity(500.0);
    reactor->setModeratorThermalCapacity(4186.0);
    reactor->setHeatTransferArea(500.0);
    reactor->setRatedPower(100.0);
    reactor->setCoolantFlowRate(17700.0);
    reactor->setCoolantInletTemperature(296.96);
  }

  void TearDown() override { reactor.reset(); }
};

TEST_F(ReactorTest, AllInitialConditionsIndependent) {
  reactor->setInitialPower(100.0);
  reactor->setInitialReactivity(0.0);
  reactor->setInitialFuelTemperature(565.0);
  reactor->setInitialModeratorTemperature(300.0);

  EXPECT_DOUBLE_EQ(reactor->getPower(), 100.0);
  EXPECT_DOUBLE_EQ(reactor->getReactivity(), 0.0);
  EXPECT_DOUBLE_EQ(reactor->getFuelTemperature(), 565.0);
  EXPECT_DOUBLE_EQ(reactor->getModeratorTemperature(), 300.0);
}

TEST_F(ReactorTest, SetHeatTransferCoefficientWithConstantFunction) {
  Function *h_const = Reactor::createConstantFunction(5000.0);
  reactor->setHeatTransferCoefficient(h_const);

  reactor->setInitialPower(100.0);
  reactor->setInitialFuelTemperature(565.0);
  reactor->setInitialModeratorTemperature(300.0);

  EXPECT_NO_THROW(reactor->timeStep(0.01));
  EXPECT_GT(reactor->getPower(), 0.0);
}

TEST_F(ReactorTest, SetHeatTransferCoefficientWithExpressionFunction) {
  Function *h_func = new Function("5000 + 10 * (x - 300)", {"x"});
  reactor->setHeatTransferCoefficient(h_func);

  reactor->setInitialPower(100.0);
  reactor->setInitialFuelTemperature(565.0);
  reactor->setInitialModeratorTemperature(300.0);

  EXPECT_NO_THROW(reactor->timeStep(0.01));
  EXPECT_GT(reactor->getPower(), 0.0);
}

TEST_F(ReactorTest, NoFunctionsDefaultBehavior) {
  reactor->setInitialPower(100.0);
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

  EXPECT_DOUBLE_EQ(reactor->getFuelMass(), 2000.0);
  EXPECT_DOUBLE_EQ(reactor->getHeatTransferArea(), 500.0);
  EXPECT_DOUBLE_EQ(reactor->getRatedPower(), 3400.0);
  EXPECT_DOUBLE_EQ(reactor->getCoolantFlowRate(), 17700.0);

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
  reactor->insertControlRod(50.0);
  double rho_after = reactor->getReactivity();

  EXPECT_LT(rho_after, rho_before);
  EXPECT_DOUBLE_EQ(rho_after, rho_before + (-0.0001 * 50.0));
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
  reactor->setInitialPower(100.0);
  reactor->setInitialFuelTemperature(565.0);
  reactor->setInitialModeratorTemperature(300.0);

  double dt = 0.01;
  for (int i = 0; i < 100; ++i) {
    EXPECT_NO_THROW(reactor->timeStep(dt));
  }

  EXPECT_NEAR(reactor->getState().time, 1.0, 1e-10);
  EXPECT_GT(reactor->getPower(), 0.0);
}

class ReactorParameterizedTest : public ReactorTest,
                                 public ::testing::WithParamInterface<double> {
};

TEST_P(ReactorParameterizedTest, VariousInitialPowers) {
  double initial_power = GetParam();

  reactor->setInitialPower(initial_power);
  reactor->setInitialFuelTemperature(565.0);
  reactor->setInitialModeratorTemperature(300.0);

  EXPECT_NO_THROW(reactor->timeStep(0.01));
  EXPECT_GE(reactor->getPower(), 0.0);
}

INSTANTIATE_TEST_SUITE_P(DifferentPowerLevels, ReactorParameterizedTest,
                         ::testing::Values(1.0, 10.0, 100.0, 1000.0));

class ReactorEdgeCaseTest : public ReactorTest {};

TEST_F(ReactorEdgeCaseTest, VeryLowPower) {
  reactor->setInitialPower(1e-6);
  reactor->setInitialFuelTemperature(565.0);
  reactor->setInitialModeratorTemperature(300.0);

  EXPECT_NO_THROW(reactor->timeStep(0.01));
  EXPECT_GE(reactor->getPower(), 0.0);
}

TEST_F(ReactorEdgeCaseTest, HighInitialReactivity) {
  reactor->setInitialPower(100.0);
  reactor->setInitialReactivity(0.01);
  reactor->setInitialFuelTemperature(565.0);
  reactor->setInitialModeratorTemperature(300.0);

  EXPECT_NO_THROW(reactor->timeStep(0.01));
}

TEST_F(ReactorEdgeCaseTest, LargeTemperatureGradient) {
  reactor->setInitialPower(100.0);
  reactor->setInitialFuelTemperature(800.0);
  reactor->setInitialModeratorTemperature(250.0);

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

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}