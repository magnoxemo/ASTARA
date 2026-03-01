#include <gtest/gtest.h>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <memory>
#include "Function.h"
#include "Reactor.h"


/*most of the tests are written by Claude*/

using namespace astara;

class ReactorTest : public ::testing::Test {
protected:
  void SetUp() override {
    // Standard 6-group delayed neutron parameters for PWR
    decay_constants = {0.0127, 0.0317, 0.115, 0.311, 1.40, 3.87};
    beta_values = {0.000215, 0.001424, 0.001274, 0.002568, 0.000748, 0.000273};
    Lambda = 5.0e-5;
    n_groups = 6;
  }

  std::vector<double> decay_constants;
  std::vector<double> beta_values;
  double Lambda;
  unsigned int n_groups;

  // Helper to calculate total beta
  double calculateTotalBeta(const std::vector<double>& betas) {
    double total = 0.0;
    for (double b : betas) {
      total += b;
    }
    return total;
  }
};

// ===== Constructor Tests =====

TEST_F(ReactorTest, ConstructorThrowsWhenGroupCountDoesNotMatchDecayConstants) {
  EXPECT_THROW(Reactor reactor(3, {0.1, 0.2}, beta_values, Lambda),
               std::runtime_error);
}

TEST_F(ReactorTest, ConstructorThrowsWhenGroupCountDoesNotMatchBetaValues) {
  EXPECT_THROW(Reactor reactor(3, decay_constants, {0.1, 0.2}, Lambda),
               std::runtime_error);
}

TEST_F(ReactorTest, ConstructorThrowsWhenNeutronGenerationTimeIsZero) {
  EXPECT_THROW(Reactor reactor(n_groups, decay_constants, beta_values, 0.0),
               std::runtime_error);
}

TEST_F(ReactorTest, ConstructorThrowsWhenNeutronGenerationTimeIsNegative) {
  EXPECT_THROW(Reactor reactor(n_groups, decay_constants, beta_values, -1e-5),
               std::runtime_error);
}

TEST_F(ReactorTest, ConstructorThrowsWhenTotalBetaIsNegative) {
  std::vector<double> negative_beta = {-0.1, 0.01};
  EXPECT_THROW(Reactor reactor(2, {0.1, 0.2}, negative_beta, 1e-5),
               std::runtime_error);
}

TEST_F(ReactorTest, ConstructorThrowsWhenTotalBetaEqualsOrExceedsOne) {
  std::vector<double> invalid_beta = {0.5, 0.5};  // Total = 1.0
  EXPECT_THROW(Reactor reactor(2, {0.1, 0.2}, invalid_beta, 1e-5),
               std::runtime_error);
}

TEST_F(ReactorTest, ConstructorSucceedsWithValidParameters) {
  EXPECT_NO_THROW(Reactor reactor(n_groups, decay_constants, beta_values, Lambda));
}

TEST_F(ReactorTest, ConstructorInitializesPrecursorConcentrations) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  
  for (unsigned int i = 0; i < n_groups; ++i) {
    EXPECT_EQ(reactor.getPrecursorConcentration(i), 0.0);
  }
}

TEST_F(ReactorTest, ConstructorThrowsOnInvalidPrecursorAccess) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  
  EXPECT_THROW(reactor.getPrecursorConcentration(n_groups), std::runtime_error);
  EXPECT_THROW(reactor.getPrecursorConcentration(n_groups + 5), std::runtime_error);
}

// ===== Initial Condition Setters =====

TEST_F(ReactorTest, SetInitialFuelTemperature) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  double T_fuel = 565.0;
  
  reactor.setInitialFuelTemperature(T_fuel);
  EXPECT_DOUBLE_EQ(reactor.getFuelTemperature(), T_fuel);
}

TEST_F(ReactorTest, SetInitialModeratorTemperature) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  double T_mod = 300.0;
  
  reactor.setInitialModeratorTemperature(T_mod);
  EXPECT_DOUBLE_EQ(reactor.getModeratorTemperature(), T_mod);
}

TEST_F(ReactorTest, SetInitialReactivity) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  double rho = 0.001;
  
  reactor.setInitialReactivity(rho);
  EXPECT_DOUBLE_EQ(reactor.getReactivity(), rho);
}

TEST_F(ReactorTest, SetInitialPower) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  double power = 100.0;
  
  reactor.setInitialPower(power);
  EXPECT_DOUBLE_EQ(reactor.getPower(), power);
}

TEST_F(ReactorTest, SetInitialPrecursorConcentration) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  double C_0 = 1.5e-5;
  
  reactor.setInitialPrecursorConcentration(0, C_0);
  EXPECT_DOUBLE_EQ(reactor.getPrecursorConcentration(0), C_0);
}

TEST_F(ReactorTest, SetInitialPrecursorConcentrationThrowsOnInvalidIndex) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  
  EXPECT_THROW(reactor.setInitialPrecursorConcentration(n_groups, 1.0),
               std::runtime_error);
}

TEST_F(ReactorTest, AllInitialConditionsIndependent) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  
  reactor.setInitialFuelTemperature(565.0);
  reactor.setInitialModeratorTemperature(300.0);
  reactor.setInitialReactivity(0.001);
  reactor.setInitialPower(100.0);
  
  EXPECT_DOUBLE_EQ(reactor.getFuelTemperature(), 565.0);
  EXPECT_DOUBLE_EQ(reactor.getModeratorTemperature(), 300.0);
  EXPECT_DOUBLE_EQ(reactor.getReactivity(), 0.001);
  EXPECT_DOUBLE_EQ(reactor.getPower(), 100.0);
}

// ===== Function Setter Tests =====


TEST_F(ReactorTest, SetHeatTransferCoefficientWithLinearFunction) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  
  std::unique_ptr<Function> h_func(
      new Function("5000 + 10 * (x - 300)", {"x"}));
  EXPECT_NO_THROW(reactor.setHeatTransferCoefficient(h_func.get()));
}

TEST_F(ReactorTest, SetFuelSpecificHeatFunction) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  
  std::unique_ptr<Function> c_func(new Function("500 + 0.2 * (x - 300)", {"x"}));
  EXPECT_NO_THROW(reactor.setFuelSpecificHeatFunction(c_func.get()));
}

TEST_F(ReactorTest, SetFuelTemperatureCoEfficientFunction) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  
  std::unique_ptr<Function> alpha_f(
      new Function("-0.00001 * (x - 565)", {"x"}));
  EXPECT_NO_THROW(reactor.setFuelTemperatureCoEfficientFunction(alpha_f.get()));
}

TEST_F(ReactorTest, SetModeratorTemperatureCoEfficientFunction) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  
  std::unique_ptr<Function> alpha_m(
      new Function("-0.00005 * (x - 300)", {"x"}));
  EXPECT_NO_THROW(
      reactor.setModeratorTemperatureCoEfficientFunction(alpha_m.get()));
}

TEST_F(ReactorTest, SetAllFunctionsDoesNotThrow) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  
  std::unique_ptr<Function> h_func(new Function("5000.0", {"x"}));
  std::unique_ptr<Function> c_func(new Function("500.0", {"x"}));
  std::unique_ptr<Function> alpha_f(new Function("-0.00001", {"x"}));
  std::unique_ptr<Function> alpha_m(new Function("-0.00005", {"x"}));
  std::unique_ptr<Function> alpha_b(new Function("0.0", {"x"}));
  
  EXPECT_NO_THROW({
    reactor.setHeatTransferCoefficient(h_func.get());
    reactor.setFuelSpecificHeatFunction(c_func.get());
    reactor.setFuelTemperatureCoEfficientFunction(alpha_f.get());
    reactor.setModeratorTemperatureCoEfficientFunction(alpha_m.get());
    reactor.setBoronTemperatureCoEfficientFunction(alpha_b.get());
  });
}

// ===== Physical Parameter Configuration Tests =====

TEST_F(ReactorTest, SetAndGetFuelMass) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  double fuel_mass = 2000.0;
  
  reactor.setFuelMass(fuel_mass);
  EXPECT_DOUBLE_EQ(reactor.getFuelMass(), fuel_mass);
}

TEST_F(ReactorTest, SetAndGetModeratorMass) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  double moderator_mass = 50000.0;
  
  reactor.setModeratorMass(moderator_mass);
  EXPECT_DOUBLE_EQ(reactor.getModeratorMass(), moderator_mass);
}

TEST_F(ReactorTest, SetAndGetFuelThermalCapacity) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  double capacity = 600.0;
  
  reactor.setFuelThermalCapacity(capacity);
  EXPECT_DOUBLE_EQ(reactor.getFuelThermalCapacity(), capacity);
}

TEST_F(ReactorTest, SetAndGetModeratorThermalCapacity) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  double capacity = 4200.0;
  
  reactor.setModeratorThermalCapacity(capacity);
  EXPECT_DOUBLE_EQ(reactor.getModeratorThermalCapacity(), capacity);
}

TEST_F(ReactorTest, SetAndGetHeatTransferArea) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  double area = 750.0;
  
  reactor.setHeatTransferArea(area);
  EXPECT_DOUBLE_EQ(reactor.getHeatTransferArea(), area);
}

TEST_F(ReactorTest, SetAndGetFissionEnergy) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  double energy = 3.2e-11;
  
  reactor.setFissionEnergy(energy);
  EXPECT_DOUBLE_EQ(reactor.getFissionEnergy(), energy);
}

TEST_F(ReactorTest, DefaultPhysicalParametersAreReasonable) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  
  EXPECT_GT(reactor.getFuelMass(), 0.0);
  EXPECT_GT(reactor.getModeratorMass(), 0.0);
  EXPECT_GT(reactor.getFuelThermalCapacity(), 0.0);
  EXPECT_GT(reactor.getModeratorThermalCapacity(), 0.0);
  EXPECT_GT(reactor.getHeatTransferArea(), 0.0);
  EXPECT_GT(reactor.getFissionEnergy(), 0.0);
}

// ===== Time Stepping Tests =====

TEST_F(ReactorTest, TimeStepThrowsOnZeroTimeStep) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  reactor.setInitialPower(100.0);
  reactor.setInitialReactivity(0.0);
  reactor.setInitialFuelTemperature(565.0);
  reactor.setInitialModeratorTemperature(300.0);
  
  EXPECT_THROW(reactor.timeStep(0.0), std::runtime_error);
}

TEST_F(ReactorTest, TimeStepThrowsOnNegativeTimeStep) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  reactor.setInitialPower(100.0);
  reactor.setInitialReactivity(0.0);
  reactor.setInitialFuelTemperature(565.0);
  reactor.setInitialModeratorTemperature(300.0);
  
  EXPECT_THROW(reactor.timeStep(-0.01), std::runtime_error);
}

TEST_F(ReactorTest, TimeStepAdvancesSimulationTime) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  reactor.setInitialPower(100.0);
  reactor.setInitialReactivity(0.0);
  reactor.setInitialFuelTemperature(565.0);
  reactor.setInitialModeratorTemperature(300.0);
  
  double initial_time = reactor.getState().time;
  double dt = 0.01;
  
  reactor.timeStep(dt);
  
  EXPECT_NEAR(reactor.getState().time, initial_time + dt, 1e-10);
}

TEST_F(ReactorTest, TimeStepMultipleTimes) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  reactor.setInitialPower(100.0);
  reactor.setInitialReactivity(0.0);
  reactor.setInitialFuelTemperature(565.0);
  reactor.setInitialModeratorTemperature(300.0);
  
  double dt = 0.01;
  int steps = 100;
  
  for (int i = 0; i < steps; ++i) {
    EXPECT_NO_THROW(reactor.timeStep(dt));
  }
  
  EXPECT_NEAR(reactor.getState().time, steps * dt, 1e-8);
}

TEST_F(ReactorTest, CriticalReactorHasSteadyPower) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  
  reactor.setInitialPower(100.0);
  reactor.setInitialReactivity(0.0);  // Critical
  reactor.setInitialFuelTemperature(565.0);
  reactor.setInitialModeratorTemperature(300.0);
  
  // Initialize precursors at steady state
  for (unsigned int i = 0; i < n_groups; ++i) {
    double C_ss = (beta_values[i] / decay_constants[i]) * 
                  (reactor.getPower() / Lambda);
    reactor.setInitialPrecursorConcentration(i, C_ss);
  }
  
  // Disable feedback effects
  std::unique_ptr<Function> zero_feedback(new Function("0.0", {"x"}));
  reactor.setFuelTemperatureCoEfficientFunction(zero_feedback.get());
  reactor.setModeratorTemperatureCoEfficientFunction(zero_feedback.get());
  
  double initial_power = reactor.getPower();
  double dt = 0.01;
  
  // Run 10 steps
  for (int i = 0; i < 10; ++i) {
    reactor.timeStep(dt);
  }
  
  // Power should remain approximately constant
  EXPECT_NEAR(reactor.getPower(), initial_power, initial_power * 0.01);
}

// ===== Control Rod Tests =====

TEST_F(ReactorTest, InsertControlRodReducesReactivity) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  reactor.setInitialReactivity(0.0);
  reactor.setControlRodEffectiveness(-0.0001);
  
  double initial_rho = reactor.getReactivity();
  reactor.insertControlRod(50.0);  // Insert 50 cm
  
  EXPECT_LT(reactor.getReactivity(), initial_rho);
}

TEST_F(ReactorTest, InsertControlRodChangeIsProportional) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  reactor.setInitialReactivity(0.0);
  reactor.setControlRodEffectiveness(-0.0001);
  
  reactor.insertControlRod(50.0);
  double rho_after_50 = reactor.getReactivity();
  
  // Reset and try 100 cm
  reactor.setInitialReactivity(0.0);
  reactor.insertControlRod(100.0);
  double rho_after_100 = reactor.getReactivity();
  
  // Change should be proportional
  EXPECT_NEAR(rho_after_100, 2.0 * rho_after_50, 1e-10);
}

TEST_F(ReactorTest, InsertControlRodThrowsOnNegativeLength) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  
  EXPECT_THROW(reactor.insertControlRod(-10.0), std::runtime_error);
}

TEST_F(ReactorTest, InsertControlRodAcceptsZeroLength) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  reactor.setInitialReactivity(0.001);
  
  double initial_rho = reactor.getReactivity();
  EXPECT_NO_THROW(reactor.insertControlRod(0.0));
  
  EXPECT_DOUBLE_EQ(reactor.getReactivity(), initial_rho);
}

// ===== Boron Injection Tests =====

TEST_F(ReactorTest, InjectBoronReducesReactivity) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  reactor.setInitialReactivity(0.0);
  reactor.setBoronEffectiveness(-1.0e-5);
  
  double initial_rho = reactor.getReactivity();
  reactor.injectBoron(1000.0);  // 1000 ppm
  
  EXPECT_LT(reactor.getReactivity(), initial_rho);
}

TEST_F(ReactorTest, InjectBoronChangeIsProportional) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  reactor.setInitialReactivity(0.0);
  reactor.setBoronEffectiveness(-1.0e-5);
  
  reactor.injectBoron(500.0);
  double rho_after_500 = reactor.getReactivity();
  
  // Reset and try 1000 ppm
  reactor.setInitialReactivity(0.0);
  reactor.injectBoron(1000.0);
  double rho_after_1000 = reactor.getReactivity();
  
  // Change should be proportional
  EXPECT_NEAR(rho_after_1000, 2.0 * rho_after_500, 1e-10);
}

TEST_F(ReactorTest, InjectBoronThrowsOnNegativeConcentration) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  
  EXPECT_THROW(reactor.injectBoron(-100.0), std::runtime_error);
}

TEST_F(ReactorTest, InjectBoronAcceptsZeroConcentration) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  reactor.setInitialReactivity(0.001);
  
  double initial_rho = reactor.getReactivity();
  EXPECT_NO_THROW(reactor.injectBoron(0.0));
  
  EXPECT_DOUBLE_EQ(reactor.getReactivity(), initial_rho);
}

// ===== State Access Tests =====

TEST_F(ReactorTest, GetStateReturnsConstReference) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  reactor.setInitialPower(100.0);
  reactor.setInitialFuelTemperature(565.0);
  
  const ReactorState& state = reactor.getState();
  
  EXPECT_DOUBLE_EQ(state.power, 100.0);
  EXPECT_DOUBLE_EQ(state.fuel_temperature, 565.0);
}

TEST_F(ReactorTest, GetMutableStateAllowsModification) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  
  reactor.getMutableState().power = 150.0;
  
  EXPECT_DOUBLE_EQ(reactor.getPower(), 150.0);
}

TEST_F(ReactorTest, StateConsistencyAcrossAccess) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  
  reactor.setInitialPower(125.5);
  reactor.setInitialReactivity(0.00123);
  reactor.setInitialFuelTemperature(575.3);
  reactor.setInitialModeratorTemperature(310.8);
  
  EXPECT_DOUBLE_EQ(reactor.getPower(), reactor.getState().power);
  EXPECT_DOUBLE_EQ(reactor.getReactivity(), reactor.getState().reactivity);
  EXPECT_DOUBLE_EQ(reactor.getFuelTemperature(), 
                   reactor.getState().fuel_temperature);
  EXPECT_DOUBLE_EQ(reactor.getModeratorTemperature(), 
                   reactor.getState().moderator_temperature);
}

// ===== Transient Analysis Tests =====

TEST_F(ReactorTest, PositiveReactivityIncreasePower) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  
  reactor.setInitialPower(100.0);
  reactor.setInitialReactivity(0.001);  // Positive reactivity
  reactor.setInitialFuelTemperature(565.0);
  reactor.setInitialModeratorTemperature(300.0);
  
  // Initialize precursors
  for (unsigned int i = 0; i < n_groups; ++i) {
    double C_ss = (beta_values[i] / decay_constants[i]) * 
                  (reactor.getPower() / Lambda);
    reactor.setInitialPrecursorConcentration(i, C_ss);
  }
  
  // Disable feedback to isolate power increase
  std::unique_ptr<Function> zero_feedback(new Function("0.0", {"x"}));
  reactor.setFuelTemperatureCoEfficientFunction(zero_feedback.get());
  reactor.setModeratorTemperatureCoEfficientFunction(zero_feedback.get());
  
  double initial_power = reactor.getPower();
  
  reactor.timeStep(0.01);
  
  EXPECT_GT(reactor.getPower(), initial_power);
}

TEST_F(ReactorTest, NegativeReactivityDecreasesPower) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  
  reactor.setInitialPower(100.0);
  reactor.setInitialReactivity(-0.001);  // Negative reactivity
  reactor.setInitialFuelTemperature(565.0);
  reactor.setInitialModeratorTemperature(300.0);
  
  // Initialize precursors
  for (unsigned int i = 0; i < n_groups; ++i) {
    double C_ss = (beta_values[i] / decay_constants[i]) * 
                  (reactor.getPower() / Lambda);
    reactor.setInitialPrecursorConcentration(i, C_ss);
  }
  
  // Disable feedback
  std::unique_ptr<Function> zero_feedback(new Function("0.0", {"x"}));
  reactor.setFuelTemperatureCoEfficientFunction(zero_feedback.get());
  reactor.setModeratorTemperatureCoEfficientFunction(zero_feedback.get());
  
  double initial_power = reactor.getPower();
  
  reactor.timeStep(0.01);
  
  EXPECT_LT(reactor.getPower(), initial_power);
}

TEST_F(ReactorTest, PowerEnforceNonNegativity) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  
  reactor.setInitialPower(0.001);  // Very low power
  reactor.setInitialReactivity(-0.1);  // Strong negative reactivity
  reactor.setInitialFuelTemperature(565.0);
  reactor.setInitialModeratorTemperature(300.0);
  
  // Initialize precursors
  for (unsigned int i = 0; i < n_groups; ++i) {
    reactor.setInitialPrecursorConcentration(i, 1e-10);
  }
  
  // Disable feedback
  std::unique_ptr<Function> zero_feedback(new Function("0.0", {"x"}));
  reactor.setFuelTemperatureCoEfficientFunction(zero_feedback.get());
  reactor.setModeratorTemperatureCoEfficientFunction(zero_feedback.get());
  
  for (int i = 0; i < 100; ++i) {
    reactor.timeStep(0.01);
    EXPECT_GE(reactor.getPower(), 0.0);
  }
}

TEST_F(ReactorTest, PrecursorConcentrationsNonNegative) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  
  reactor.setInitialPower(100.0);
  reactor.setInitialReactivity(-0.01);
  reactor.setInitialFuelTemperature(565.0);
  reactor.setInitialModeratorTemperature(300.0);
  
  // Initialize with very small precursors
  for (unsigned int i = 0; i < n_groups; ++i) {
    reactor.setInitialPrecursorConcentration(i, 1e-12);
  }
  
  for (int i = 0; i < 50; ++i) {
    reactor.timeStep(0.01);
    for (unsigned int j = 0; j < n_groups; ++j) {
      EXPECT_GE(reactor.getPrecursorConcentration(j), 0.0);
    }
  }
}

// ===== Temperature Feedback Tests =====

TEST_F(ReactorTest, FuelTemperatureIncreasesWithPower) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  
  reactor.setInitialPower(100.0);
  reactor.setInitialReactivity(0.0);
  reactor.setInitialFuelTemperature(565.0);
  reactor.setInitialModeratorTemperature(300.0);
  
  // Configure heat transfer
  std::unique_ptr<Function> h_func(new Function("5000.0", {"x"}));
  reactor.setHeatTransferCoefficient(h_func.get());
  
  std::unique_ptr<Function> c_func(new Function("500.0", {"x"}));
  reactor.setFuelSpecificHeatFunction(c_func.get());
  
  // Disable feedback
  std::unique_ptr<Function> zero_feedback(new Function("0.0", {"x"}));
  reactor.setFuelTemperatureCoEfficientFunction(zero_feedback.get());
  reactor.setModeratorTemperatureCoEfficientFunction(zero_feedback.get());
  
  double initial_T_fuel = reactor.getFuelTemperature();
  
  reactor.timeStep(0.1);
  
  EXPECT_GT(reactor.getFuelTemperature(), initial_T_fuel);
}

TEST_F(ReactorTest, NegativeFuelFeedbackStabilizesPower) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  
  reactor.setInitialPower(100.0);
  reactor.setInitialReactivity(0.002);  // Positive step
  reactor.setInitialFuelTemperature(565.0);
  reactor.setInitialModeratorTemperature(300.0);
  
  // Configure functions
  std::unique_ptr<Function> h_func(new Function("5000.0", {"x"}));
  reactor.setHeatTransferCoefficient(h_func.get());
  
  std::unique_ptr<Function> c_func(new Function("500.0", {"x"}));
  reactor.setFuelSpecificHeatFunction(c_func.get());
  
  // Negative fuel feedback
  std::unique_ptr<Function> alpha_f(
      new Function("-0.00001 * (x - 565)", {"x"}));
  reactor.setFuelTemperatureCoEfficientFunction(alpha_f.get());
  
  std::unique_ptr<Function> zero_feedback(new Function("0.0", {"x"}));
  reactor.setModeratorTemperatureCoEfficientFunction(zero_feedback.get());
  
  // Initialize precursors
  for (unsigned int i = 0; i < n_groups; ++i) {
    double C_ss = (beta_values[i] / decay_constants[i]) * 
                  (reactor.getPower() / Lambda);
    reactor.setInitialPrecursorConcentration(i, C_ss);
  }
  
  double power_0 = reactor.getPower();
  reactor.timeStep(0.01);
  double power_1 = reactor.getPower();
  reactor.timeStep(0.01);
  double power_2 = reactor.getPower();
  
  // With negative feedback, power should not increase indefinitely
  EXPECT_LT(power_2 - power_1, power_1 - power_0);
}

// ===== Integration Tests =====

TEST_F(ReactorTest, FullTransientSimulation) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  
  // Set up realistic reactor state
  reactor.setInitialPower(100.0);
  reactor.setInitialReactivity(0.0);
  reactor.setInitialFuelTemperature(565.0);
  reactor.setInitialModeratorTemperature(300.0);
  
  // Configure all functions
  std::unique_ptr<Function> h_func(
      new Function("5000 + 10 * (x - 300)", {"x"}));
  reactor.setHeatTransferCoefficient(h_func.get());
  
  std::unique_ptr<Function> c_func(
      new Function("500 + 0.2 * (x - 300)", {"x"}));
  reactor.setFuelSpecificHeatFunction(c_func.get());
  
  std::unique_ptr<Function> alpha_f(
      new Function("-0.00001 * (x - 565)", {"x"}));
  reactor.setFuelTemperatureCoEfficientFunction(alpha_f.get());
  
  std::unique_ptr<Function> alpha_m(
      new Function("-0.00005 * (x - 300)", {"x"}));
  reactor.setModeratorTemperatureCoEfficientFunction(alpha_m.get());
  
  // Set physical parameters
  reactor.setFuelMass(2000.0);
  reactor.setModeratorMass(50000.0);
  reactor.setHeatTransferArea(500.0);
  
  // Initialize precursors
  for (unsigned int i = 0; i < n_groups; ++i) {
    double C_ss = (beta_values[i] / decay_constants[i]) * 
                  (reactor.getPower() / Lambda);
    reactor.setInitialPrecursorConcentration(i, C_ss);
  }
  
  // Run 100 steps
  double dt = 0.01;
  for (int i = 0; i < 100; ++i) {
    EXPECT_NO_THROW(reactor.timeStep(dt));
  }
  
  // Verify state is physically reasonable
  EXPECT_GT(reactor.getPower(), 0.0);
  EXPECT_GT(reactor.getFuelTemperature(), 300.0);
  EXPECT_GT(reactor.getModeratorTemperature(), 200.0);
  EXPECT_EQ(reactor.getState().time, 1.0);  // 100 * 0.01 = 1.0
  
  // All precursor concentrations should be non-negative
  for (unsigned int i = 0; i < n_groups; ++i) {
    EXPECT_GE(reactor.getPrecursorConcentration(i), 0.0);
  }
}

// ===== Edge Cases =====

TEST_F(ReactorTest, VerySmallTimeStep) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  reactor.setInitialPower(100.0);
  reactor.setInitialReactivity(0.001);
  reactor.setInitialFuelTemperature(565.0);
  reactor.setInitialModeratorTemperature(300.0);
  
  EXPECT_NO_THROW(reactor.timeStep(1e-6));
}

TEST_F(ReactorTest, VeryLargeTimeStep) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  reactor.setInitialPower(100.0);
  reactor.setInitialReactivity(0.0);
  reactor.setInitialFuelTemperature(565.0);
  reactor.setInitialModeratorTemperature(300.0);
  
  // Should complete without error (though accuracy may suffer)
  EXPECT_NO_THROW(reactor.timeStep(10.0));
}

TEST_F(ReactorTest, ZeroInitialPower) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  reactor.setInitialPower(0.0);
  reactor.setInitialReactivity(0.5);
  reactor.setInitialFuelTemperature(565.0);
  reactor.setInitialModeratorTemperature(300.0);
  
  EXPECT_NO_THROW(reactor.timeStep(0.01));
  EXPECT_GE(reactor.getPower(), 0.0);
}

TEST_F(ReactorTest, HighTemperatureFeedback) {
  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);
  reactor.setInitialPower(1000.0);
  reactor.setInitialReactivity(0.005);
  reactor.setInitialFuelTemperature(800.0);
  reactor.setInitialModeratorTemperature(350.0);
  
  // Strong negative feedback
  std::unique_ptr<Function> alpha_f(
      new Function("-0.0001 * (x - 565)", {"x"}));
  reactor.setFuelTemperatureCoEfficientFunction(alpha_f.get());
  
  EXPECT_NO_THROW(reactor.timeStep(0.01));
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
