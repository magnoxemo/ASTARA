#ifndef REACTOR_H
#define REACTOR_H

#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

namespace astara {

// Forward declaration
class Function;

// Structure to hold reactor state variables
struct ReactorState {
  double time = 0.0;
  double power = 0.0;                           // Thermal power [MW]
  double reactivity = 0.0;                      // Reactivity [delta-k/k]
  double fuel_temperature = 0.0;                // Fuel temperature [C]
  double moderator_temperature = 0.0;           // Moderator temperature [C]
  std::vector<double> precursor_concentrations; // Delayed neutron precursors
};

class Reactor {
public:
  // Constructor
  Reactor(unsigned int n_groups, const std::vector<double> &neutron_group_const,
          const std::vector<double> &delayed_neutron_constants,
          double neutron_generation_time);

  // Destructor
  ~Reactor() = default;

  // Disable copy constructor and assignment
  Reactor(const Reactor &) = delete;
  Reactor &operator=(const Reactor &) = delete;

  // Allow move operations
  Reactor(Reactor &&) = default;
  Reactor &operator=(Reactor &&) = default;

  // Time integration
  void timeStep(double dt);

  // State getters
  ReactorState stateSnapshot() const;
  void restoreState(const ReactorState &snapshot);

  // Individual state getters
  double getPower() const { return _state.power; }
  double getReactivity() const { return _state.reactivity; }
  double getFuelTemperature() const { return _state.fuel_temperature; }
  double getModeratorTemperature() const {
    return _state.moderator_temperature;
  }
  const ReactorState &getState() const { return _state; }
  ReactorState &getMutableState() { return _state; }

  // Physical parameter getters
  double getFuelMass() const { return _fuel_mass; }
  double getModeratorMass() const { return _moderator_mass; }
  double getFuelThermalCapacity() const { return _thermal_capacity_fuel; }
  double getModeratorThermalCapacity() const {
    return _thermal_capacity_moderator;
  }
  double getHeatTransferArea() const { return _heat_transfer_area; }
  double getFissionEnergy() const { return _fission_energy; }
  double getControlRodEffectiveness() const {
    return _control_rod_effectiveness;
  }
  double getBoronEffectiveness() const { return _boron_effectiveness; }

  // Initial condition setters
  void setInitialPower(double power) { _state.power = power; }
  void setInitialReactivity(double reactivity) {
    _state.reactivity = reactivity;
  }
  void setInitialFuelTemperature(double temp) {
    _state.fuel_temperature = temp;
  }
  void setInitialModeratorTemperature(double temp) {
    _state.moderator_temperature = temp;
  }
  void setInitialPrecursorConcentration(unsigned int group,
                                        double concentration) {
    if (group < _number_of_neutron_groups) {
      _state.precursor_concentrations[group] = concentration;
    }
  }

  // Physical parameter setters
  void setFuelMass(double mass) { _fuel_mass = mass; }
  void setModeratorMass(double mass) { _moderator_mass = mass; }
  void setFuelThermalCapacity(double capacity) {
    _thermal_capacity_fuel = capacity;
  }
  void setModeratorThermalCapacity(double capacity) {
    _thermal_capacity_moderator = capacity;
  }
  void setHeatTransferArea(double area) { _heat_transfer_area = area; }
  void setFissionEnergy(double energy) { _fission_energy = energy; }

  // Control effectiveness setters
  void setControlRodEffectiveness(double effectiveness) {
    _control_rod_effectiveness = effectiveness;
  }
  void setBoronEffectiveness(double effectiveness) {
    _boron_effectiveness = effectiveness;
  }

  // Function pointer setters - taking raw pointers (for test compatibility)
  // The Reactor will take ownership and wrap in unique_ptr
  void setHeatTransferCoefficient(Function *func) {
    _convective_heat_transfer_coefficient.reset(func);
  }
  void setFuelSpecificHeatFunction(Function *func) {
    _fuel_specific_heat.reset(func);
  }
  void setFuelTemperatureCoEfficientFunction(Function *func) {
    _fuel_temperature_coefficient.reset(func);
  }
  void setModeratorTemperatureCoEfficientFunction(Function *func) {
    _moderator_temperature_coefficient.reset(func);
  }
  void setBoronTemperatureCoEfficientFunction(Function *func) {
    _boron_temperature_coefficient.reset(func);
  }

  // unique_ptr versions (for modern code)
  void setHeatTransferCoefficient(std::unique_ptr<Function> func) {
    _convective_heat_transfer_coefficient = std::move(func);
  }
  void setFuelSpecificHeatFunction(std::unique_ptr<Function> func) {
    _fuel_specific_heat = std::move(func);
  }
  void setFuelTemperatureCoEfficientFunction(std::unique_ptr<Function> func) {
    _fuel_temperature_coefficient = std::move(func);
  }
  void
  setModeratorTemperatureCoEfficientFunction(std::unique_ptr<Function> func) {
    _moderator_temperature_coefficient = std::move(func);
  }
  void setBoronTemperatureCoEfficientFunction(std::unique_ptr<Function> func) {
    _boron_temperature_coefficient = std::move(func);
  }

  // Static helper for creating constant functions
  static Function *createConstantFunction(double value);

  // Control actions
  void insertControlRod(double length);
  void injectBoron(double boron_concentration);

  void print_states() {
    std::cout << std::setw(10) << this->getState().time << std::setw(14)
              << this->getPower() << std::setw(14)
              << this->getReactivity() * 1e5 // Convert to pcm
              << std::setw(14) << this->getFuelTemperature() << std::setw(14)
              << this->getModeratorTemperature() << std::endl;
  };

private:
  // Core reactor parameters
  unsigned int _number_of_neutron_groups;
  std::vector<double> _neutron_group_decay_constants;
  std::vector<double> _delayed_neutron_fractions;
  double _neutron_generation_time;
  double _total_delayed_neutron_fraction;
  double _prompt_neutron_fraction;

  // Physical parameters
  double _fuel_mass = 0.0;                  // [kg]
  double _moderator_mass = 0.0;             // [kg]
  double _thermal_capacity_fuel = 0.0;      // [J/(kg·K)]
  double _thermal_capacity_moderator = 0.0; // [J/(kg·K)]
  double _heat_transfer_area = 0.0;         // [m^2]
  double _fission_energy = 3.2e-11;         // [J/fission]

  // Control parameters
  double _control_rod_effectiveness = 0.0; // [delta-k/k per meter]
  double _boron_effectiveness = 0.0;       // [delta-k/k per ppm]

  // Function pointers for temperature-dependent properties
  std::unique_ptr<Function> _convective_heat_transfer_coefficient;
  std::unique_ptr<Function> _fuel_specific_heat;
  std::unique_ptr<Function> _fuel_temperature_coefficient;
  std::unique_ptr<Function> _moderator_temperature_coefficient;
  std::unique_ptr<Function> _boron_temperature_coefficient;

  // Current state
  ReactorState _state;

  // Derivative calculations
  double calculateDRhoDt() const;
  double calculateDPowerDt() const;
  double calculateDFuelTempDt() const;
  double calculateDModeratorTempDt() const;
  double calculateDCDt(unsigned int group_index) const;

  ReactorState evaluateDerivatives(const ReactorState &s) {
    ReactorState saved = _state;
    _state = s;

    ReactorState k;
    k.reactivity = calculateDRhoDt();
    k.power = calculateDPowerDt();
    k.fuel_temperature = calculateDFuelTempDt();
    k.moderator_temperature = calculateDModeratorTempDt();

    k.precursor_concentrations.resize(_number_of_neutron_groups);
    for (unsigned int i = 0; i < _number_of_neutron_groups; ++i) {
      k.precursor_concentrations[i] = calculateDCDt(i);
    }

    _state = saved;
    return k;
  }

  void clampState() {
    // Power cannot be negative
    if (_state.power < 0.0) {
      _state.power = 0.0;
    }

    // Temperatures cannot be below absolute zero (Celsius scale)
    const double absolute_zero_c = -273.15;

    if (_state.fuel_temperature < absolute_zero_c) {
      _state.fuel_temperature = absolute_zero_c;
    }

    if (_state.moderator_temperature < absolute_zero_c) {
      _state.moderator_temperature = absolute_zero_c;
    }

    // Precursor concentrations must remain non-negative
    for (unsigned int i = 0; i < _number_of_neutron_groups; ++i) {
      if (_state.precursor_concentrations[i] < 0.0) {
        _state.precursor_concentrations[i] = 0.0;
      }
    }
  }

  // Integration methods
  void integrateRK4(double dt);
};

    ReactorState operator+(const ReactorState &a, const ReactorState &b) {
        ReactorState r = a;

        r.reactivity += b.reactivity;
        r.power += b.power;
        r.fuel_temperature += b.fuel_temperature;
        r.moderator_temperature += b.moderator_temperature;

        for (size_t i = 0; i < r.precursor_concentrations.size(); ++i) {
            r.precursor_concentrations[i] += b.precursor_concentrations[i];
        }
        return r;
    }

    ReactorState operator*(const ReactorState &a, double scalar) {
        ReactorState r;

        r.reactivity = a.reactivity * scalar;
        r.power = a.power * scalar;
        r.fuel_temperature = a.fuel_temperature * scalar;
        r.moderator_temperature = a.moderator_temperature * scalar;

        r.precursor_concentrations.resize(a.precursor_concentrations.size());
        for (size_t i = 0; i < a.precursor_concentrations.size(); ++i) {
            r.precursor_concentrations[i] = a.precursor_concentrations[i] * scalar;
        }
        return r;
    }
} // namespace astara

#endif // REACTOR_H