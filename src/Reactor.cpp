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

        // Initialize additional state variables
        _state.upper_plenum_temperature = 300.0;
        _state.lower_plenum_temperature = 300.0;
        _state.hot_leg_temperature = 300.0;
        _state.cold_leg_temperature = 300.0;

        // Initialize fuel and coolant node temperatures
        _state.fuel_temperatures.resize(3, 565.0);
        _state.coolant_temperatures.resize(6, 300.0);
    }

    void Reactor::timeStep(double dt) {
        if (dt <= 0.0)
            throw std::runtime_error("Time step must be positive");

        integrateRK4(dt);
        _state.time += dt;

        // Update average fuel temperature for compatibility with existing getters
        if (!_state.fuel_temperatures.empty()) {
            double sum = 0.0;
            for (double t : _state.fuel_temperatures) {
                sum += t;
            }
            _state.fuel_temperature = sum / _state.fuel_temperatures.size();
        }

        // Update average moderator temperature
        if (!_state.coolant_temperatures.empty()) {
            double sum = 0.0;
            for (double t : _state.coolant_temperatures) {
                sum += t;
            }
            _state.moderator_temperature = sum / _state.coolant_temperatures.size();
        }
    }

    double Reactor::calculateDRhoDt() const {
        return 0.0;
    }

    double Reactor::calculateTotalReactivity() const {
        double rho_total = _state.reactivity;

        // Calculate average fuel temperature from nodes
        double avg_fuel_temp = 0.0;
        for (double t : _state.fuel_temperatures) {
            avg_fuel_temp += t;
        }
        avg_fuel_temp /= _state.fuel_temperatures.size();

        // Calculate average coolant temperature from nodes
        double avg_coolant_temp = 0.0;
        for (double t : _state.coolant_temperatures) {
            avg_coolant_temp += t;
        }
        avg_coolant_temp /= _state.coolant_temperatures.size();

        if (_fuel_temperature_coefficient) {
            try {
                double alpha_f = (*_fuel_temperature_coefficient)(avg_fuel_temp);
                if (!std::isnan(alpha_f) && !std::isinf(alpha_f)) {
                    rho_total += alpha_f;
                }
            } catch (const std::exception &e) {
            }
        }

        if (_moderator_temperature_coefficient) {
            try {
                double alpha_m = (*_moderator_temperature_coefficient)(avg_coolant_temp);
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

    double Reactor::calculateFuelNodeDt(int node) const {
        if (node < 0 || node >= 3 || _fuel_mass_per_node <= 0.0 || _thermal_capacity_fuel <= 0.0) {
            return 0.0;
        }

        const double thermal_power_watts = _state.power * _rated_thermal_power_MW * 1e6;
        const double power_fraction = _fraction_power_in_fuel / 3.0; // Distribute equally among 3 nodes

        double h = 5000.0;
        if (_convective_heat_transfer_coefficient) {
            try {
                h = (*_convective_heat_transfer_coefficient)(_state.fuel_temperatures[node]);
                if (std::isnan(h) || std::isinf(h) || h < 0.0) {
                    h = 5000.0;
                }
            } catch (const std::exception &e) {
                h = 5000.0;
            }
        }

        // Each fuel node transfers heat to two adjacent coolant nodes
        // For node 0: transfers to coolant nodes 0 and 1
        // For node 1: transfers to coolant nodes 2 and 3
        // For node 2: transfers to coolant nodes 4 and 5
        int coolant_node1 = node * 2;
        int coolant_node2 = node * 2 + 1;

        double avg_coolant_temp = (_state.coolant_temperatures[coolant_node1] +
                                   _state.coolant_temperatures[coolant_node2]) / 2.0;

        const double temperature_difference = _state.fuel_temperatures[node] - avg_coolant_temp;
        const double heat_transferred = h * _heat_transfer_area_per_node * temperature_difference;

        double c_f = _thermal_capacity_fuel;
        if (_fuel_specific_heat) {
            try {
                c_f = (*_fuel_specific_heat)(_state.fuel_temperatures[node]);
                if (std::isnan(c_f) || std::isinf(c_f) || c_f <= 0.0) {
                    c_f = _thermal_capacity_fuel;
                }
            } catch (const std::exception &e) {
                c_f = _thermal_capacity_fuel;
            }
        }

        double dt_f_dt = (power_fraction * thermal_power_watts - heat_transferred) / (_fuel_mass_per_node * c_f);
        return dt_f_dt;
    }

    double Reactor::calculateCoolantNodeDt(int node) const {
        if (node < 0 || node >= 6 || _coolant_mass_per_node <= 0.0 || _thermal_capacity_moderator <= 0.0) {
            return 0.0;
        }

        double h = 5000.0;
        if (_convective_heat_transfer_coefficient) {
            try {
                h = (*_convective_heat_transfer_coefficient)(_state.fuel_temperatures[node/2]);
                if (std::isnan(h) || std::isinf(h) || h < 0.0) {
                    h = 5000.0;
                }
            } catch (const std::exception &e) {
                h = 5000.0;
            }
        }

        // Determine which fuel node heats this coolant node
        int fuel_node = node / 2;
        const double temperature_difference = _state.fuel_temperatures[fuel_node] - _state.coolant_temperatures[node];
        const double heat_from_fuel = h * _heat_transfer_area_per_node * temperature_difference;

        const double cp_c = _thermal_capacity_moderator;

        // Flow between coolant nodes
        double inflow_temp = 0.0;
        if (node == 0) {
            // First coolant node receives from lower plenum
            inflow_temp = _state.lower_plenum_temperature;
        } else {
            // Other nodes receive from previous node
            inflow_temp = _state.coolant_temperatures[node - 1];
        }

        const double heat_from_flow = _coolant_mass_flow_rate_per_node * cp_c * (inflow_temp - _state.coolant_temperatures[node]);

        double dt_m_dt = (heat_from_fuel + heat_from_flow) / (_coolant_mass_per_node * cp_c);
        return dt_m_dt;
    }

    double Reactor::calculateDUpperPlenumDt() const {
        if (_upper_plenum_mass <= 0.0 || _thermal_capacity_moderator <= 0.0) {
            return 0.0;
        }

        // Upper plenum receives flow from last coolant node (node 5)
        const double cp_c = _thermal_capacity_moderator;
        const double heat_from_flow = _coolant_mass_flow_rate_per_node * cp_c *
                                      (_state.coolant_temperatures[5] - _state.upper_plenum_temperature);

        return heat_from_flow / (_upper_plenum_mass * cp_c);
    }

    double Reactor::calculateDLowerPlenumDt() const {
        if (_lower_plenum_mass <= 0.0 || _thermal_capacity_moderator <= 0.0) {
            return 0.0;
        }

        // Lower plenum receives from cold leg
        const double cp_c = _thermal_capacity_moderator;
        const double heat_from_flow = _coolant_mass_flow_rate_per_node * cp_c *
                                      (_state.cold_leg_temperature - _state.lower_plenum_temperature);

        return heat_from_flow / (_lower_plenum_mass * cp_c);
    }

    double Reactor::calculateDHotLegDt() const {
        if (_hot_leg_mass <= 0.0 || _thermal_capacity_moderator <= 0.0) {
            return 0.0;
        }

        // Hot leg receives from upper plenum
        const double cp_c = _thermal_capacity_moderator;
        const double heat_from_flow = _coolant_mass_flow_rate_per_node * cp_c *
                                      (_state.upper_plenum_temperature - _state.hot_leg_temperature);

        return heat_from_flow / (_hot_leg_mass * cp_c);
    }

    double Reactor::calculateDColdLegDt() const {
        if (_cold_leg_mass <= 0.0 || _thermal_capacity_moderator <= 0.0) {
            return 0.0;
        }

        // Cold leg receives from steam generator outlet
        const double cp_c = _thermal_capacity_moderator;
        const double heat_from_flow = _coolant_mass_flow_rate_per_node * cp_c *
                                      (_coolant_inlet_temperature - _state.cold_leg_temperature);

        return heat_from_flow / (_cold_leg_mass * cp_c);
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

    void Reactor::setFuelMass(double mass) {
        _fuel_mass = mass;
        _fuel_mass_per_node = mass / 3.0;
    }

    void Reactor::setModeratorMass(double mass) {
        _moderator_mass = mass;
        _coolant_mass_per_node = mass / 6.0;
    }

    void Reactor::setCoolantFlowRate(double flow_rate) {
        _coolant_mass_flow_rate = flow_rate;
        _coolant_mass_flow_rate_per_node = flow_rate / 6.0;
    }

    void Reactor::setHeatTransferArea(double area) {
        _heat_transfer_area = area;
        _heat_transfer_area_per_node = area / 3.0;
    }

    void Reactor::setUpperPlenumMass(double mass) {
        _upper_plenum_mass = mass;
    }

    void Reactor::setLowerPlenumMass(double mass) {
        _lower_plenum_mass = mass;
    }

    void Reactor::setHotLegMass(double mass) {
        _hot_leg_mass = mass;
    }

    void Reactor::setColdLegMass(double mass) {
        _cold_leg_mass = mass;
    }

    void Reactor::setFractionPowerInFuel(double fraction) {
        _fraction_power_in_fuel = fraction;
    }

    ReactorState Reactor::evaluateDerivatives(const ReactorState &s) {
        ReactorState saved = _state;
        _state = s;

        ReactorState k;
        k.reactivity = calculateDRhoDt();
        k.power = calculateDPowerDt();

        // Calculate derivatives for each fuel node
        k.fuel_temperatures.resize(3);
        for (int i = 0; i < 3; ++i) {
            k.fuel_temperatures[i] = calculateFuelNodeDt(i);
        }

        // Calculate derivatives for each coolant node
        k.coolant_temperatures.resize(6);
        for (int i = 0; i < 6; ++i) {
            k.coolant_temperatures[i] = calculateCoolantNodeDt(i);
        }

        // Calculate plenum and leg derivatives
        k.upper_plenum_temperature = calculateDUpperPlenumDt();
        k.lower_plenum_temperature = calculateDLowerPlenumDt();
        k.hot_leg_temperature = calculateDHotLegDt();
        k.cold_leg_temperature = calculateDColdLegDt();

        // For compatibility with existing code
        if (!k.fuel_temperatures.empty()) {
            double sum = 0.0;
            for (double t : k.fuel_temperatures) {
                sum += t;
            }
            k.fuel_temperature = sum / k.fuel_temperatures.size();
        }

        if (!k.coolant_temperatures.empty()) {
            double sum = 0.0;
            for (double t : k.coolant_temperatures) {
                sum += t;
            }
            k.moderator_temperature = sum / k.coolant_temperatures.size();
        }

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
        constexpr double max_temp = 2000.0;

        // Clamp fuel node temperatures
        for (double& temp : _state.fuel_temperatures) {
            if (temp < absolute_zero_c) temp = absolute_zero_c;
            if (temp > max_temp) temp = max_temp;
        }

        // Clamp coolant node temperatures
        for (double& temp : _state.coolant_temperatures) {
            if (temp < absolute_zero_c) temp = absolute_zero_c;
            if (temp > max_temp) temp = max_temp;
        }

        // Clamp plenum and leg temperatures
        if (_state.upper_plenum_temperature < absolute_zero_c) _state.upper_plenum_temperature = absolute_zero_c;
        if (_state.upper_plenum_temperature > max_temp) _state.upper_plenum_temperature = max_temp;

        if (_state.lower_plenum_temperature < absolute_zero_c) _state.lower_plenum_temperature = absolute_zero_c;
        if (_state.lower_plenum_temperature > max_temp) _state.lower_plenum_temperature = max_temp;

        if (_state.hot_leg_temperature < absolute_zero_c) _state.hot_leg_temperature = absolute_zero_c;
        if (_state.hot_leg_temperature > max_temp) _state.hot_leg_temperature = max_temp;

        if (_state.cold_leg_temperature < absolute_zero_c) _state.cold_leg_temperature = absolute_zero_c;
        if (_state.cold_leg_temperature > max_temp) _state.cold_leg_temperature = max_temp;


        if (!_state.fuel_temperatures.empty()) {
            double sum = 0.0;
            for (double t : _state.fuel_temperatures) {
                sum += t;
            }
            _state.fuel_temperature = sum / _state.fuel_temperatures.size();
        }

        if (!_state.coolant_temperatures.empty()) {
            double sum = 0.0;
            for (double t : _state.coolant_temperatures) {
                sum += t;
            }
            _state.moderator_temperature = sum / _state.coolant_temperatures.size();
        }

        for (unsigned int i = 0; i < _number_of_neutron_groups; ++i) {
            if (_state.precursor_concentrations[i] < 0.0) {
                _state.precursor_concentrations[i] = 0.0;
            }
        }
    }




    void Reactor::recordState(std::ostream& csv_file) {
        csv_file << std::fixed << std::setprecision(6)
                 << _state.time << ","
                 << _state.power * _rated_thermal_power_MW << ","
                 << _state.reactivity * 1e5 << ","  // Convert to pcm
                 << _state.fuel_temperature << ","
                 << _state.moderator_temperature << ","
                 << _state.upper_plenum_temperature << ","
                 << _state.lower_plenum_temperature << ","
                 << _state.hot_leg_temperature << ","
                 << _state.cold_leg_temperature;

        // Write individual fuel node temperatures
        for (double temp : _state.fuel_temperatures) {
            csv_file << "," << temp;
        }

        // Write individual coolant node temperatures
        for (double temp : _state.coolant_temperatures) {
            csv_file << "," << temp;
        }

        // Write precursor concentrations
        for (double conc : _state.precursor_concentrations) {
            csv_file << "," << conc;
        }

        csv_file << "\n";
    }

    void Reactor::writeCSVHeader(std::ostream& csv_file) {
        csv_file << "Time(s),Power(MW),Rho(pcm),T_fuel(C),T_mod(C),"
                 << "T_up(C),T_lp(C),T_hot(C),T_cold(C)";

        // Add fuel node headers
        for (size_t i = 0; i < _state.fuel_temperatures.size(); ++i) {
            csv_file << ",Fuel_Node" << i << "(C)";
        }

        // Add coolant node headers
        for (size_t i = 0; i < _state.coolant_temperatures.size(); ++i) {
            csv_file << ",Coolant_Node" << i << "(C)";
        }

        // Add precursor headers
        for (size_t i = 0; i < _state.precursor_concentrations.size(); ++i) {
            csv_file << ",Precursor_" << i;
        }

        csv_file << "\n";
    }

    ReactorState operator*(const ReactorState &a, double scalar) {
        ReactorState r;

        r.reactivity = a.reactivity * scalar;
        r.power = a.power * scalar;
        r.fuel_temperature = a.fuel_temperature * scalar;
        r.moderator_temperature = a.moderator_temperature * scalar;

        r.upper_plenum_temperature = a.upper_plenum_temperature * scalar;
        r.lower_plenum_temperature = a.lower_plenum_temperature * scalar;
        r.hot_leg_temperature = a.hot_leg_temperature * scalar;
        r.cold_leg_temperature = a.cold_leg_temperature * scalar;

        r.fuel_temperatures.resize(a.fuel_temperatures.size());
        for (size_t i = 0; i < a.fuel_temperatures.size(); ++i) {
            r.fuel_temperatures[i] = a.fuel_temperatures[i] * scalar;
        }

        r.coolant_temperatures.resize(a.coolant_temperatures.size());
        for (size_t i = 0; i < a.coolant_temperatures.size(); ++i) {
            r.coolant_temperatures[i] = a.coolant_temperatures[i] * scalar;
        }

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

        r.upper_plenum_temperature += b.upper_plenum_temperature;
        r.lower_plenum_temperature += b.lower_plenum_temperature;
        r.hot_leg_temperature += b.hot_leg_temperature;
        r.cold_leg_temperature += b.cold_leg_temperature;

        for (size_t i = 0; i < r.fuel_temperatures.size() && i < b.fuel_temperatures.size(); ++i) {
            r.fuel_temperatures[i] += b.fuel_temperatures[i];
        }

        for (size_t i = 0; i < r.coolant_temperatures.size() && i < b.coolant_temperatures.size(); ++i) {
            r.coolant_temperatures[i] += b.coolant_temperatures[i];
        }

        for (size_t i = 0; i < r.precursor_concentrations.size(); ++i) {
            r.precursor_concentrations[i] += b.precursor_concentrations[i];
        }
        return r;
    }

}