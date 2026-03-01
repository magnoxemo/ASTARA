#include "Reactor.h"
#include <stdexcept>
#include <numeric>
#include <algorithm>

namespace astara {

Reactor::Reactor(unsigned int n_groups,
                 const std::vector<double>& neutron_group_const,
                 const std::vector<double>& delayed_neutron_constants,
                 double neutron_generation_time)
    : _number_of_neutron_groups(n_groups),
      _neutron_group_decay_constants(neutron_group_const),
      _delayed_neutron_fractions(delayed_neutron_constants),
      _neutron_generation_time(neutron_generation_time),
      _total_delayed_neutron_fraction(
          std::accumulate(delayed_neutron_constants.begin(),
                         delayed_neutron_constants.end(), 0.0)),
      _prompt_neutron_fraction(1.0 - 
          std::accumulate(delayed_neutron_constants.begin(),
                         delayed_neutron_constants.end(), 0.0)) {
  
  // Validate input consistency
  if (neutron_group_const.size() != n_groups) {
    throw std::runtime_error(
        "Number of decay constants does not match number of groups");
  }
  
  if (delayed_neutron_constants.size() != n_groups) {
    throw std::runtime_error(
        "Number of delayed neutron fractions does not match number of groups");
  }
  
  if (neutron_generation_time <= 0.0) {
    throw std::runtime_error("Neutron generation time must be positive");
  }
  
  if (_total_delayed_neutron_fraction < 0.0 || _total_delayed_neutron_fraction >= 1.0) {
    throw std::runtime_error("Total delayed neutron fraction must be in [0, 1)");
  }
  
  // Initialize precursor concentrations
  _state.precursor_concentrations.resize(n_groups, 0.0);
  
  // Create default (constant) functions - users can override
  // Default: heat transfer = 5000 W/(m²·K), fuel specific heat = 500 J/(kg·K)
  _convective_heat_transfer_coefficient = std::make_unique<Function>("5000.0", std::vector<std::string>{"x"});
  _fuel_specific_heat = std::make_unique<Function>("500.0", std::vector<std::string>{"x"});
  _fuel_temperature_coefficient = std::make_unique<Function>("0.0", std::vector<std::string>{"x"});
  _moderator_temperature_coefficient = std::make_unique<Function>("0.0", std::vector<std::string>{"x"});
  _boron_temperature_coefficient = std::make_unique<Function>("0.0", std::vector<std::string>{"x"});
}

void Reactor::timeStep(double dt) {
  if (dt <= 0.0) {
    throw std::runtime_error("Time step must be positive");
  }
  integrateRK4(dt);
  _state.time += dt;
}

double Reactor::calculateDRhoDt() const {
  // Reactivity feedback from fuel temperature
  double fuel_feedback = 0.0;
  if (_fuel_temperature_coefficient) {
    fuel_feedback = (*_fuel_temperature_coefficient)(_state.fuel_temperature);
  }
  
  // Reactivity feedback from moderator temperature
  double moderator_feedback = 0.0;
  if (_moderator_temperature_coefficient) {
    moderator_feedback = (*_moderator_temperature_coefficient)(_state.moderator_temperature);
  }
  
  // Reactivity feedback from boron absorption
  double boron_feedback = 0.0;
  if (_boron_temperature_coefficient) {
    boron_feedback = (*_boron_temperature_coefficient)(_state.moderator_temperature);
  }
  
  // Total feedback reactivity rate (1/s)
  double drho_dt = fuel_feedback + moderator_feedback + boron_feedback;
  
  return drho_dt;
}

double Reactor::calculateDPowerDt() const {
  // Point kinetics equation for power
  // dP/dt = [(rho - beta) / Lambda] * P + sum(lambda_i * C_i) * P
  
  if (_state.power < 1e-12) {
    // Avoid numerical issues at very low power
    return 0.0;
  }
  
  const double rho_minus_beta = _state.reactivity - _total_delayed_neutron_fraction;
  const double inv_generation_time = 1.0 / _neutron_generation_time;
  
  // First term: prompt neutron contribution
  double dp_dt = (rho_minus_beta * inv_generation_time) * _state.power;
  
  // Second term: delayed neutron contribution from precursors
  for (unsigned int i = 0; i < _number_of_neutron_groups; ++i) {
    double lambda_i = _neutron_group_decay_constants[i];
    double c_i = _state.precursor_concentrations[i];
    dp_dt += lambda_i * c_i * _state.power;
  }
  
  return dp_dt;
}

double Reactor::calculateDFuelTempDt() const {
  // Energy balance for fuel:
  // dT_f/dt = (Power - Heat_transfer) / (m_f * c_f)
  
  if (_state.power < 1e-12 || _fuel_mass <= 0.0 || _thermal_capacity_fuel <= 0.0) {
    return 0.0;
  }
  
  // Fission heat generation (converted from power in MW to watts)
  const double fission_power_watts = _state.power * 1e6; // Assuming power is in MW
  const double fission_heat = fission_power_watts / _fission_energy; // Fissions per second
  const double heat_generated = fission_power_watts; // (J/s)
  
  // Heat transfer coefficient (W/(m²·K))
  double h = 1.0;
  if (_convective_heat_transfer_coefficient) {
    h = (*_convective_heat_transfer_coefficient)(_state.fuel_temperature);
  }
  
  // Heat transfer to moderator
  const double temperature_difference = _state.fuel_temperature - _state.moderator_temperature;
  const double heat_transferred = h * _heat_transfer_area * temperature_difference;
  
  // Fuel specific heat capacity (J/(kg·K))
  double c_f = _thermal_capacity_fuel;
  if (_fuel_specific_heat) {
    c_f = (*_fuel_specific_heat)(_state.fuel_temperature);
  }
  
  // Temperature rate of change
  double dt_f_dt = (heat_generated - heat_transferred) / (_fuel_mass * c_f);
  
  return dt_f_dt;
}

double Reactor::calculateDCDt(unsigned int group_index) const {
  // Precursor concentration equation:
  // dC_i/dt = (beta_i / Lambda) * P - lambda_i * C_i
  
  if (group_index >= _number_of_neutron_groups) {
    throw std::runtime_error("Precursor group index out of range");
  }
  
  const double beta_i = _delayed_neutron_fractions[group_index];
  const double lambda_i = _neutron_group_decay_constants[group_index];
  const double inv_generation_time = 1.0 / _neutron_generation_time;
  
  // Production term from fission
  double dc_dt = beta_i * inv_generation_time * _state.power;
  
  // Decay term
  dc_dt -= lambda_i * _state.precursor_concentrations[group_index];
  
  return dc_dt;
}

ReactorState Reactor::stateSnapshot() const {
  return _state;
}

void Reactor::restoreState(const ReactorState& snapshot) {
  _state = snapshot;
}

void Reactor::integrateRK4(double dt) {
  // Stage 1: Evaluate at current point
  ReactorState k1_state;
  
  k1_state.reactivity = calculateDRhoDt();
  k1_state.power = calculateDPowerDt();
  k1_state.fuel_temperature = calculateDFuelTempDt();
  k1_state.precursor_concentrations.resize(_number_of_neutron_groups);
  for (unsigned int i = 0; i < _number_of_neutron_groups; ++i) {
    k1_state.precursor_concentrations[i] = calculateDCDt(i);
  }
  
  // Stage 2: Evaluate at t + dt/2 using k1
  ReactorState saved_state = stateSnapshot();
  _state.reactivity += k1_state.reactivity * dt / 2.0;
  _state.power += k1_state.power * dt / 2.0;
  _state.fuel_temperature += k1_state.fuel_temperature * dt / 2.0;
  _state.moderator_temperature += 0.0 * dt / 2.0; // Moderator doesn't change directly (could be added)
  for (unsigned int i = 0; i < _number_of_neutron_groups; ++i) {
    _state.precursor_concentrations[i] += k1_state.precursor_concentrations[i] * dt / 2.0;
  }
  
  ReactorState k2_state;
  k2_state.reactivity = calculateDRhoDt();
  k2_state.power = calculateDPowerDt();
  k2_state.fuel_temperature = calculateDFuelTempDt();
  k2_state.precursor_concentrations.resize(_number_of_neutron_groups);
  for (unsigned int i = 0; i < _number_of_neutron_groups; ++i) {
    k2_state.precursor_concentrations[i] = calculateDCDt(i);
  }
  
  // Stage 3: Evaluate at t + dt/2 using k2
  restoreState(saved_state);
  _state.reactivity += k2_state.reactivity * dt / 2.0;
  _state.power += k2_state.power * dt / 2.0;
  _state.fuel_temperature += k2_state.fuel_temperature * dt / 2.0;
  for (unsigned int i = 0; i < _number_of_neutron_groups; ++i) {
    _state.precursor_concentrations[i] += k2_state.precursor_concentrations[i] * dt / 2.0;
  }
  
  ReactorState k3_state;
  k3_state.reactivity = calculateDRhoDt();
  k3_state.power = calculateDPowerDt();
  k3_state.fuel_temperature = calculateDFuelTempDt();
  k3_state.precursor_concentrations.resize(_number_of_neutron_groups);
  for (unsigned int i = 0; i < _number_of_neutron_groups; ++i) {
    k3_state.precursor_concentrations[i] = calculateDCDt(i);
  }
  
  // Stage 4: Evaluate at t + dt using k3
  restoreState(saved_state);
  _state.reactivity += k3_state.reactivity * dt;
  _state.power += k3_state.power * dt;
  _state.fuel_temperature += k3_state.fuel_temperature * dt;
  for (unsigned int i = 0; i < _number_of_neutron_groups; ++i) {
    _state.precursor_concentrations[i] += k3_state.precursor_concentrations[i] * dt;
  }
  
  ReactorState k4_state;
  k4_state.reactivity = calculateDRhoDt();
  k4_state.power = calculateDPowerDt();
  k4_state.fuel_temperature = calculateDFuelTempDt();
  k4_state.precursor_concentrations.resize(_number_of_neutron_groups);
  for (unsigned int i = 0; i < _number_of_neutron_groups; ++i) {
    k4_state.precursor_concentrations[i] = calculateDCDt(i);
  }
  
  // Combine stages with RK4 weights
  restoreState(saved_state);
  
  _state.reactivity += (dt / 6.0) * (k1_state.reactivity + 2.0 * k2_state.reactivity + 
                                     2.0 * k3_state.reactivity + k4_state.reactivity);
  _state.power += (dt / 6.0) * (k1_state.power + 2.0 * k2_state.power + 
                               2.0 * k3_state.power + k4_state.power);
  _state.fuel_temperature += (dt / 6.0) * (k1_state.fuel_temperature + 2.0 * k2_state.fuel_temperature + 
                                           2.0 * k3_state.fuel_temperature + k4_state.fuel_temperature);
  
  for (unsigned int i = 0; i < _number_of_neutron_groups; ++i) {
    _state.precursor_concentrations[i] += 
        (dt / 6.0) * (k1_state.precursor_concentrations[i] + 
                     2.0 * k2_state.precursor_concentrations[i] + 
                     2.0 * k3_state.precursor_concentrations[i] + 
                     k4_state.precursor_concentrations[i]);
  }
  
  // Ensure physical constraints
  if (_state.power < 0.0) {
    _state.power = 0.0;
  }
  
  for (unsigned int i = 0; i < _number_of_neutron_groups; ++i) {
    if (_state.precursor_concentrations[i] < 0.0) {
      _state.precursor_concentrations[i] = 0.0;
    }
  }
}

void Reactor::insertControlRod(double length) {
  if (length < 0.0) {
    throw std::runtime_error("Control rod length cannot be negative");
  }
  
  // Control rod changes reactivity proportionally to insertion length
  double reactivity_change = _control_rod_effectiveness * length;
  _state.reactivity += reactivity_change;
}

void Reactor::injectBoron(double boron_concentration) {
  if (boron_concentration < 0.0) {
    throw std::runtime_error("Boron concentration cannot be negative");
  }
  
  // Boron absorption reduces reactivity proportionally to concentration
  double reactivity_change = _boron_effectiveness * boron_concentration;
  _state.reactivity += reactivity_change;
}

} // namespace astara
