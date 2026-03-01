#ifndef ASTARA_REACTOR_H
#define ASTARA_REACTOR_H

#include "Function.h"
#include <memory>
#include <vector>
#include <cmath>

namespace astara {

/**
 * Represents the current state of the reactor core.
 * All state variables are stored here for access by other components.
 * All the documenation strings are geneated using Claude
 */
struct ReactorState {
  double reactivity = 0.0;
  double power = 0.0;
  double fuel_temperature = 0.0;
  double moderator_temperature = 0.0;
  std::vector<double> precursor_concentrations;
  double time = 0.0;
};

/**
 * Nuclear reactor simulator with thermal feedback.
 * Models neutron kinetics with delayed neutrons, thermal transport,
 * and reactivity feedback from temperature changes.
 */
class Reactor {
public:
  /**
   * Constructor initializes reactor parameters.
   * @param n_groups Number of delayed neutron groups
   * @param neutron_group_const Decay constants (1/s) for each precursor group
   * @param delayed_neutron_constants Delayed neutron fractions (beta) for each group
   * @param neutron_generation_time Prompt neutron generation time (s)
   */
  Reactor(unsigned int n_groups,
          const std::vector<double>& neutron_group_const,
          const std::vector<double>& delayed_neutron_constants,
          double neutron_generation_time);

  ~Reactor() = default;

  /* ===== Initial Condition Setters ===== */
  
  void setInitialFuelTemperature(double fuel_temperature) {
    _state.fuel_temperature = fuel_temperature;
  }

  void setInitialModeratorTemperature(double moderator_temperature) {
    _state.moderator_temperature = moderator_temperature;
  }

  void setInitialReactivity(double reactivity) {
    _state.reactivity = reactivity;
  }

  void setInitialPower(double power) {
    _state.power = power;
  }

  void setInitialPrecursorConcentration(unsigned int group_index, 
                                       double concentration) {
    if (group_index >= _state.precursor_concentrations.size()) {
      throw std::runtime_error("Precursor group index out of range");
    }
    _state.precursor_concentrations[group_index] = concentration;
  }

  /* ===== Function Setters for Temperature/Time Dependencies ===== */
  
  void setHeatTransferCoefficient(Function* heat_transfer_coeff_function) {
    _convective_heat_transfer_coefficient.reset(heat_transfer_coeff_function);
  }

  void setFuelSpecificHeatFunction(Function* fuel_specific_function) {
    _fuel_specific_heat.reset(fuel_specific_function);
  }

  void setFuelTemperatureCoEfficientFunction(Function* fuel_temp_feedback_func) {
    _fuel_temperature_coefficient.reset(fuel_temp_feedback_func);
  }

  void setModeratorTemperatureCoEfficientFunction(Function* moderator_temp_feedback_func) {
    _moderator_temperature_coefficient.reset(moderator_temp_feedback_func);
  }

  void setBoronTemperatureCoEfficientFunction(Function* boron_temp_feedback_func) {
    _boron_temperature_coefficient.reset(boron_temp_feedback_func);
  }

  /* ===== Simulation Control ===== */
  
  /**
   * Advance the reactor simulation by one time step using RK4 integration.
   * @param dt Time step size (seconds)
   */
  void timeStep(double dt);

  /**
   * Get the current reactor state
   */
  const ReactorState& getState() const {
    return _state;
  }

  /**
   * Get mutable reference to reactor state
   */
  ReactorState& getMutableState() {
    return _state;
  }

  /* ===== Reactivity Control ===== */
  
  /**
   * Insert control rod into reactor core.
   * @param length Length of inserted control rod (cm)
   */
  void insertControlRod(double length);

  /**
   * Inject boron into cooling system for absorption.
   * @param boron_concentration Boron concentration (ppm)
   */
  void injectBoron(double boron_concentration);

  /* ===== State Access ===== */
  
  double getReactivity() const { return _state.reactivity; }
  double getPower() const { return _state.power; }
  double getFuelTemperature() const { return _state.fuel_temperature; }
  double getModeratorTemperature() const { return _state.moderator_temperature; }
  double getPrecursorConcentration(unsigned int group) const {
    if (group >= _state.precursor_concentrations.size()) {
      throw std::runtime_error("Precursor group index out of range");
    }
    return _state.precursor_concentrations[group];
  }

private:
  /* ===== Differential Equations (ODEs) ===== */
  
  /**
   * Calculate drho/dt - reactivity change due to feedback and control
   */
  double calculateDRhoDt() const;

  /**
   * Calculate dP/dt - power change from neutron kinetics
   */
  double calculateDPowerDt() const;

  /**
   * Calculate dT_fuel/dt - fuel temperature change from heat generation and transfer
   */
  double calculateDFuelTempDt() const;

  /**
   * Calculate dC_i/dt - precursor concentration for group i
   * @param group_index Index of the precursor group
   */
  double calculateDCDt(unsigned int group_index) const;

  /**
   * Apply Runge-Kutta 4th order integration step
   */
  void integrateRK4(double dt);

  /**
   * Create a copy of current state for RK4 calculations
   */
  ReactorState stateSnapshot() const;

  /**
   * Restore state from a snapshot
   */
  void restoreState(const ReactorState& snapshot);

  /* ===== Parameter Storage ===== */
  
  // Physical constants
  const unsigned int _number_of_neutron_groups;
  const std::vector<double> _neutron_group_decay_constants;  // lambda_i (1/s)
  const std::vector<double> _delayed_neutron_fractions;      // beta_i
  const double _neutron_generation_time;                      // Lambda (s)
  
  // Derived constants
  const double _total_delayed_neutron_fraction;  // beta_total = sum(beta_i)
  const double _prompt_neutron_fraction;         // (1 - beta)
  
  // Temperature feedback functions
  std::unique_ptr<Function> _convective_heat_transfer_coefficient;
  std::unique_ptr<Function> _fuel_specific_heat;
  std::unique_ptr<Function> _fuel_temperature_coefficient;      // alpha_f
  std::unique_ptr<Function> _moderator_temperature_coefficient; // alpha_m
  std::unique_ptr<Function> _boron_temperature_coefficient;     // alpha_b
  
  // Reactor core state
  ReactorState _state;
  
  // Physical parameters (configurable, no magic numbers)
  double _fuel_mass = 1.0;                          // kg
  double _moderator_mass = 100.0;                   // kg
  double _thermal_capacity_fuel = 500.0;            // J/(kg·K)
  double _thermal_capacity_moderator = 4186.0;      // J/(kg·K)
  double _heat_transfer_area = 10.0;                // m^2
  double _fission_energy = 3.2e-11;                 // J/fission
  
  // Control parameters
  double _control_rod_effectiveness = -0.001;       // (1/cm) - reactivity change per cm
  double _boron_effectiveness = -1.0e-4;            // (1/ppm) - reactivity change per ppm

public:
  /* ===== Configuration Methods ===== */
  
  void setFuelMass(double mass) { _fuel_mass = mass; }
  void setModeratorMass(double mass) { _moderator_mass = mass; }
  void setFuelThermalCapacity(double capacity) { _thermal_capacity_fuel = capacity; }
  void setModeratorThermalCapacity(double capacity) { _thermal_capacity_moderator = capacity; }
  void setHeatTransferArea(double area) { _heat_transfer_area = area; }
  void setFissionEnergy(double energy) { _fission_energy = energy; }
  void setControlRodEffectiveness(double effectiveness) { _control_rod_effectiveness = effectiveness; }
  void setBoronEffectiveness(double effectiveness) { _boron_effectiveness = effectiveness; }

  /* Getters */
  
  double getFuelMass() const { return _fuel_mass; }
  double getModeratorMass() const { return _moderator_mass; }
  double getFuelThermalCapacity() const { return _thermal_capacity_fuel; }
  double getModeratorThermalCapacity() const { return _thermal_capacity_moderator; }
  double getHeatTransferArea() const { return _heat_transfer_area; }
  double getFissionEnergy() const { return _fission_energy; }
};

} // namespace astara

#endif // ASTARA_REACTOR_H
