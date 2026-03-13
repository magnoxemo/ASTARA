// Reactor.cpp
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
              _total_delayed_neutron_fraction(std::accumulate(delayed_neutron_constants.begin(), delayed_neutron_constants.end(), 0.0)),
              _prompt_neutron_fraction(1.0 - _total_delayed_neutron_fraction) {

        if (neutron_group_const.size() != n_groups)
            throw std::runtime_error("Number of decay constants does not match number of groups");

        if (delayed_neutron_constants.size() != n_groups)
            throw std::runtime_error("Number of delayed neutron fractions does not match number of groups");

        if (neutron_generation_time <= 0.0)
            throw std::runtime_error("Neutron generation time must be positive");

        if (_total_delayed_neutron_fraction < 0.0 || _total_delayed_neutron_fraction >= 1.0)
            throw std::runtime_error("Total delayed neutron fraction must be in [0, 1)");

        _state.precursor_concentrations.resize(n_groups, 0.0);
    }

    void Reactor::timeStep(double dt) {
        if (dt <= 0.0)
            throw std::runtime_error("Time step must be positive");

        integrateRK4(dt);
        _state.time += dt;
    }

    double Reactor::calculateDRhoDt() const {
        return 0.0;
    }

    double Reactor::calculateTotalReactivity() const {
        double rho_total = _state.reactivity;

        if (_fuel_temperature_coefficient) {
            try {
                double alpha_f = (*_fuel_temperature_coefficient)(_state.fuel_temperature);
                if (!std::isnan(alpha_f) && !std::isinf(alpha_f)) {
                    rho_total += alpha_f;
                }
            } catch (const std::exception &e) {
            }
        }

        if (_moderator_temperature_coefficient) {
            try {
                double alpha_m = (*_moderator_temperature_coefficient)(_state.moderator_temperature);
                if (!std::isnan(alpha_m) && !std::isinf(alpha_m)) {
                    rho_total += alpha_m;
                }
            } catch (const std::exception &e) {
            }
        }

        return rho_total;
    }

    double Reactor::calculateDPowerDt() const {
        if (_state.power < 1e-12) {
            return 0.0;
        }

        double total_reactivity = calculateTotalReactivity();

        const double rho_minus_beta = total_reactivity - _total_delayed_neutron_fraction;
        const double inv_generation_time = 1.0 / _neutron_generation_time;

        double dp_dt = (rho_minus_beta * inv_generation_time) * _state.power;

        for (unsigned int i = 0; i < _number_of_neutron_groups; ++i) {
            dp_dt += _neutron_group_decay_constants[i] * _state.precursor_concentrations[i];
        }

        return dp_dt;
    }

    double Reactor::calculateDFuelTempDt() const {
        if (_state.power < 1e-12 || _fuel_mass <= 0.0 || _thermal_capacity_fuel <= 0.0) {
            return 0.0;
        }

        const double thermal_power_watts = _state.power * _rated_thermal_power_MW * 1e6;

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

        const double temperature_difference = _state.fuel_temperature - _state.moderator_temperature;
        const double heat_transferred = h * _heat_transfer_area * temperature_difference;

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

        double dt_f_dt = (thermal_power_watts - heat_transferred) / (_fuel_mass * c_f);

        return dt_f_dt;
    }

    double Reactor::calculateDModeratorTempDt() const {
        if (_moderator_mass <= 0.0 || _thermal_capacity_moderator <= 0.0) {
            return 0.0;
        }

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

        const double temperature_difference = _state.fuel_temperature - _state.moderator_temperature;
        const double heat_from_fuel = h * _heat_transfer_area * temperature_difference;

        const double cp_c = _thermal_capacity_moderator;
        const double heat_from_flow = _coolant_mass_flow_rate * cp_c *
                                      (_coolant_inlet_temperature - _state.moderator_temperature);

        double dt_m_dt = (heat_from_fuel + heat_from_flow) / (_moderator_mass * cp_c);

        return dt_m_dt;
    }

    double Reactor::calculateDCDt(unsigned int group_index) const {
        if (group_index >= _number_of_neutron_groups) {
            throw std::runtime_error("Precursor group index out of range");
        }

        const double beta_i = _delayed_neutron_fractions[group_index];
        const double lambda_i = _neutron_group_decay_constants[group_index];
        const double inv_generation_time = 1.0 / _neutron_generation_time;

        double dc_dt = beta_i * inv_generation_time * _state.power;
        dc_dt -= lambda_i * _state.precursor_concentrations[group_index];

        return dc_dt;
    }

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

        double reactivity_change = _control_rod_worth * length;
        _state.reactivity += reactivity_change;
    }

    void Reactor::injectBoron(double boron_concentration) {
        if (boron_concentration < 0.0) {
            throw std::runtime_error("Boron concentration cannot be negative");
        }

        double reactivity_change = _boron_worth * boron_concentration;
        _state.reactivity += reactivity_change;
    }

    Function *Reactor::createConstantFunction(double value) {
        return new Function(std::to_string(value), std::vector<std::string>{});
    }

    ReactorState Reactor::evaluateDerivatives(const ReactorState &s) {
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

    void Reactor::clampState() {
        if (_state.power < 0.0) {
            _state.power = 0.0;
        }

        constexpr double absolute_zero_c = -273.15;

        if (_state.fuel_temperature < absolute_zero_c) {
            _state.fuel_temperature = absolute_zero_c;
        }

        if (_state.moderator_temperature < absolute_zero_c) {
            _state.moderator_temperature = absolute_zero_c;
        }

        for (unsigned int i = 0; i < _number_of_neutron_groups; ++i) {
            if (_state.precursor_concentrations[i] < 0.0) {
                _state.precursor_concentrations[i] = 0.0;
            }
        }
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

} // namespace astara