#include "Reactor.h"
#include "Function.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <string>

namespace astara {

Reactor::Reactor(unsigned int n_groups,
                 const std::vector<double> &neutron_group_const,
                 const std::vector<double> &delayed_neutron_constants,
                 double neutron_generation_time)
    : _number_of_neutron_groups(n_groups),
      _neutron_group_decay_constants(neutron_group_const),
      _delayed_neutron_fractions(delayed_neutron_constants),
      _neutron_generation_time(neutron_generation_time),
      _total_delayed_neutron_fraction(
          std::accumulate(delayed_neutron_constants.begin(),
                          delayed_neutron_constants.end(), 0.0)),
      _prompt_neutron_fraction(1.0 - _total_delayed_neutron_fraction) {

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

  if (_total_delayed_neutron_fraction < 0.0 ||
      _total_delayed_neutron_fraction >= 1.0) {
    throw std::runtime_error(
        "Total delayed neutron fraction must be in [0, 1)");
  }

  // Initialize precursor concentrations
  _state.precursor_concentrations.resize(n_groups, 0.0);
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
    try {
      fuel_feedback = (*_fuel_temperature_coefficient)(_state.fuel_temperature);
      if (std::isnan(fuel_feedback) || std::isinf(fuel_feedback)) {
        fuel_feedback = 0.0;
      }
    } catch (const std::exception &e) {
      fuel_feedback = 0.0;
    }
  }

  // Reactivity feedback from moderator temperature
  double moderator_feedback = 0.0;
  if (_moderator_temperature_coefficient) {
    try {
      moderator_feedback =
          (*_moderator_temperature_coefficient)(_state.moderator_temperature);
      if (std::isnan(moderator_feedback) || std::isinf(moderator_feedback)) {
        moderator_feedback = 0.0;
      }
    } catch (const std::exception &e) {
      moderator_feedback = 0.0;
    }
  }

  // Return rate of change of reactivity (dρ/dt)
  return fuel_feedback + moderator_feedback;
}

double Reactor::calculateDPowerDt() const {
  // Point kinetics equation for power
  // dP/dt = [(rho - beta) / Lambda] * P + sum(lambda_i * C_i)

  if (_state.power < 1e-12) {
    return 0.0;
  }

  const double rho_minus_beta =
      _state.reactivity - _total_delayed_neutron_fraction;
  const double inv_generation_time = 1.0 / _neutron_generation_time;

  // Prompt neutron contribution
  double dp_dt = (rho_minus_beta * inv_generation_time) * _state.power;

  // Delayed neutron contribution from precursors
  for (unsigned int i = 0; i < _number_of_neutron_groups; ++i) {
    dp_dt +=
        _neutron_group_decay_constants[i] * _state.precursor_concentrations[i];
  }

  return dp_dt;
}

double Reactor::calculateDFuelTempDt() const {
  // Energy balance for fuel: dT_f/dt = (Power - Heat_transfer) / (m_f * c_f)

  if (_state.power < 1e-12 || _fuel_mass <= 0.0 ||
      _thermal_capacity_fuel <= 0.0) {
    return 0.0;
  }

  // Fission heat generation (MW to W)
  const double fission_power_watts = _state.power * 1e6;

  // Heat transfer coefficient
  double h = 5000.0; // Default value
  if (_convective_heat_transfer_coefficient) {
    try {
      h = (*_convective_heat_transfer_coefficient)(_state.fuel_temperature);
      if (std::isnan(h) || std::isinf(h) || h < 0.0) {
        h = 5000.0;
      }
    } catch (const std::exception &e) {
      h = 5000.0;
    }
  }

  // Heat transfer to moderator
  const double temperature_difference =
      _state.fuel_temperature - _state.moderator_temperature;
  const double heat_transferred =
      h * _heat_transfer_area * temperature_difference;

  // Fuel specific heat capacity
  double c_f = _thermal_capacity_fuel;
  if (_fuel_specific_heat) {
    try {
      c_f = (*_fuel_specific_heat)(_state.fuel_temperature);
      if (std::isnan(c_f) || std::isinf(c_f) || c_f <= 0.0) {
        c_f = _thermal_capacity_fuel;
      }
    } catch (const std::exception &e) {
      c_f = _thermal_capacity_fuel;
    }
  }

  // Temperature rate of change
  double dt_f_dt =
      (fission_power_watts - heat_transferred) / (_fuel_mass * c_f);

  return dt_f_dt;
}

double Reactor::calculateDModeratorTempDt() const {
  // Energy balance for moderator: dT_m/dt = Heat_transfer / (m_m * c_m)

  if (_moderator_mass <= 0.0 || _thermal_capacity_moderator <= 0.0) {
    return 0.0;
  }

  // Heat transfer coefficient
  double h = 5000.0;
  if (_convective_heat_transfer_coefficient) {
    try {
      h = (*_convective_heat_transfer_coefficient)(_state.fuel_temperature);
      if (std::isnan(h) || std::isinf(h) || h < 0.0) {
        h = 5000.0;
      }
    } catch (const std::exception &e) {
      h = 5000.0;
    }
  }

  // Heat received from fuel
  const double temperature_difference =
      _state.fuel_temperature - _state.moderator_temperature;
  const double heat_transferred =
      h * _heat_transfer_area * temperature_difference;

  // Moderator temperature rate of change
  double dt_m_dt =
      heat_transferred / (_moderator_mass * _thermal_capacity_moderator);

  return dt_m_dt;
}

double Reactor::calculateDCDt(unsigned int group_index) const {
  // Precursor concentration equation: dC_i/dt = (beta_i / Lambda) * P -
  // lambda_i * C_i

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

ReactorState Reactor::stateSnapshot() const { return _state; }

void Reactor::restoreState(const ReactorState &snapshot) { _state = snapshot; }

void Reactor::integrateRK4(double dt) {
  ReactorState y0 = _state;

  ReactorState k1 = evaluateDerivatives(y0);
  ReactorState k2 = evaluateDerivatives(y0 + k1 * (dt * 0.5));
  ReactorState k3 = evaluateDerivatives(y0 + k2 * (dt * 0.5));
  ReactorState k4 = evaluateDerivatives(y0 + k3 * dt);

  _state = y0 + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (dt / 6.0);

  clampState();
}

void Reactor::insertControlRod(double length) {
  if (length < 0.0) {
    throw std::runtime_error("Control rod length cannot be negative");
  }

  double reactivity_change = _control_rod_effectiveness * length;
  _state.reactivity += reactivity_change;
}

void Reactor::injectBoron(double boron_concentration) {
  if (boron_concentration < 0.0) {
    throw std::runtime_error("Boron concentration cannot be negative");
  }

  double reactivity_change = _boron_effectiveness * boron_concentration;
  _state.reactivity += reactivity_change;
}

// Static helper for creating constant functions
Function *Reactor::createConstantFunction(double value) {
  // Create a function that returns a constant value
  // This assumes Function can be constructed with a string expression
  return new Function(std::to_string(value), std::vector<std::string>{});
}

} // namespace astara