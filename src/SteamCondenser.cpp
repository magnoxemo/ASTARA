// Condenser.cpp
#include "SteamCondenser.h"
#include "CoolProp.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace astara {

    Condenser::Condenser() {
        _state.steam_pressure = 10.0e4;
        _state.air_pressure = 9.0e3;
        _state.steam_temperature = getSaturationTemperature(_state.steam_pressure);
        _state.steam_mass = 100.0;
        _state.air_mass = 10.0;
        _state.hotwell_water_mass = 5000.0;
        _state.hotwell_temperature = 273.15 + 80.0;
        updateDerivedQuantities();
    }

    void Condenser::timeStep(double dt) {
        if (dt <= 0.0)
            throw std::runtime_error("Time step must be positive");

        integrateRK4(dt);
        _state.time += dt;
        updateDerivedQuantities();
    }

    double Condenser::calculateDSteamMassDt() const {
        const double R_mixture = (_state.air_pressure * _steam_gas_constant) /
                                 (_state.air_pressure * _air_gas_constant +
                                  _state.steam_pressure * _steam_gas_constant);

        const double steam_in_air_outlet = _steam_air_outlet_flow * (1.0 - R_mixture);

        const double condensation_rate = getCondensationRate();

        const double dm_steam_dt = _turbine_inlet_flow +
                                   _other_steam_inlet_flow -
                                   condensation_rate -
                                   steam_in_air_outlet;

        return dm_steam_dt;
    }

    double Condenser::calculateDAirMassDt() const {
        const double dm_air_dt = _vacuum_breaker_flow +
                                 _drain_condenser_flow +
                                 _steam_gas_outlet_flow -
                                 _steam_air_outlet_flow;

        return dm_air_dt;
    }

    double Condenser::calculateDSteamPressureDt() const {
        const double dm_steam_dt = calculateDSteamMassDt();

        const double dp_steam_dt = (_steam_gas_constant * _state.steam_temperature * dm_steam_dt) /
                                   _volume;

        return dp_steam_dt;
    }

    double Condenser::calculateDAirPressureDt() const {
        const double dm_air_dt = calculateDAirMassDt();

        const double dp_air_dt = (_air_gas_constant * _state.steam_temperature * dm_air_dt) /
                                 _volume;

        return dp_air_dt;
    }

    double Condenser::calculateDHotwellMassDt() const {
        const double condensation_rate = getCondensationRate();

        const double dm_hotwell_dt = condensation_rate +
                                     _oxygen_bubbling_flow -
                                     _hotwell_water_outlet_flow;

        return dm_hotwell_dt;
    }

    double Condenser::calculateDHotwellTemperatureDt() const {
        const double h_steam = getSteamEnthalpy(_turbine_inlet_temperature,
                                                _state.steam_pressure);
        const double h_cold_water = getWaterEnthalpy(_cold_water_temperature, 101325.0);
        const double h_hot_water = getWaterEnthalpy(_state.hotwell_temperature,
                                                    _state.steam_pressure);

        const double T_steam_sat = getSaturationTemperature(_state.steam_pressure);

        const double delta_T_lmtd = (_hot_water_outlet_temperature - _cold_water_temperature) /
                                    std::log((T_steam_sat - _cold_water_temperature) /
                                             (T_steam_sat - _hot_water_outlet_temperature));

        const double Q_heat_transfer = _UA * delta_T_lmtd;

        const double Q_from_condensate = _state.condensate_flow_rate * h_steam;
        const double Q_from_oxygen = _oxygen_bubbling_flow * h_cold_water;
        const double Q_to_outlet = _hotwell_water_outlet_flow * h_hot_water;
        const double Q_to_cold_water = _cold_water_flow * _water_specific_heat *
                                       (_hot_water_outlet_temperature - _cold_water_temperature);

        const double dT_dt = (Q_from_condensate + Q_from_oxygen - Q_to_outlet - Q_to_cold_water) /
                             (_state.hotwell_water_mass * _water_specific_heat);

        return dT_dt;
    }

    double Condenser::getCondensationRate() const {
        const double h_steam = getSteamEnthalpy(_turbine_inlet_temperature,
                                                _state.steam_pressure);
        const double h_cold_water = getWaterEnthalpy(_cold_water_temperature, 101325.0);

        const double T_steam_sat = getSaturationTemperature(_state.steam_pressure);

        const double delta_T_lmtd = (_hot_water_outlet_temperature - _cold_water_temperature) /
                                    std::log((T_steam_sat - _cold_water_temperature) /
                                             (T_steam_sat - _hot_water_outlet_temperature));

        const double condensation_rate = (_UA * delta_T_lmtd) / (h_steam - h_cold_water);

        return std::max(0.0, condensation_rate);
    }

    double Condenser::getSteamEnthalpy(double temperature, double pressure) const {
        try {
            return CoolProp::PropsSI("H", "T", temperature, "P", pressure, "Water");
        } catch (...) {
            return 2.7e6;
        }
    }

    double Condenser::getWaterEnthalpy(double temperature, double pressure) const {
        try {
            return CoolProp::PropsSI("H", "T", temperature, "P", pressure, "Water");
        } catch (...) {
            return 4.2e3 * (temperature - 273.15);
        }
    }

    double Condenser::getSaturationTemperature(double pressure) const {
        try {
            return CoolProp::PropsSI("T", "P", pressure, "Q", 0.0, "Water");
        } catch (...) {
            return 273.15 + 100.0;
        }
    }

    double Condenser::getWaterDensity(double temperature, double pressure) const {
        try {
            return CoolProp::PropsSI("D", "T", temperature, "P", pressure, "Water");
        } catch (...) {
            return 1000.0;
        }
    }

    void Condenser::integrateRK4(double dt) {
        CondenserState y0 = _state;

        CondenserState k1 = evaluateDerivatives(y0);
        CondenserState k2 = evaluateDerivatives(y0 + k1 * (dt * 0.5));
        CondenserState k3 = evaluateDerivatives(y0 + k2 * (dt * 0.5));
        CondenserState k4 = evaluateDerivatives(y0 + k3 * dt);

        _state = y0 + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (dt / 6.0);

        clampState();
    }

    CondenserState Condenser::evaluateDerivatives(const CondenserState &s) {
        CondenserState saved = _state;
        _state = s;
        _state.steam_temperature = getSaturationTemperature(_state.steam_pressure);
        _state.condensate_flow_rate = getCondensationRate();

        CondenserState k;
        k.steam_mass = calculateDSteamMassDt();
        k.air_mass = calculateDAirMassDt();
        k.steam_pressure = calculateDSteamPressureDt();
        k.air_pressure = calculateDAirPressureDt();
        k.hotwell_water_mass = calculateDHotwellMassDt();
        k.hotwell_temperature = calculateDHotwellTemperatureDt();

        _state = saved;
        return k;
    }

    void Condenser::clampState() {
        if (_state.steam_pressure < 1000.0) {
            _state.steam_pressure = 1000.0;
        }
        if (_state.steam_pressure > 500.0e3) {
            _state.steam_pressure = 500.0e3;
        }

        if (_state.air_pressure < 100.0) {
            _state.air_pressure = 100.0;
        }
        if (_state.air_pressure > 100.0e3) {
            _state.air_pressure = 100.0e3;
        }

        if (_state.steam_mass < 0.1) {
            _state.steam_mass = 0.1;
        }

        if (_state.air_mass < 0.01) {
            _state.air_mass = 0.01;
        }

        if (_state.hotwell_water_mass < 100.0) {
            _state.hotwell_water_mass = 100.0;
        }

        if (_state.hotwell_temperature < 273.15) {
            _state.hotwell_temperature = 273.15;
        }
        if (_state.hotwell_temperature > 473.15) {
            _state.hotwell_temperature = 473.15;
        }
    }

    void Condenser::updateDerivedQuantities() {
        _state.steam_temperature = getSaturationTemperature(_state.steam_pressure);
        _state.total_pressure = _state.steam_pressure + _state.air_pressure;
        _state.condensate_flow_rate = getCondensationRate();

        const double water_density = getWaterDensity(_state.hotwell_temperature,
                                                     _state.steam_pressure);
        _state.water_level = _state.hotwell_water_mass / (_hotwell_area * water_density);
    }

    CondenserState operator*(const CondenserState &a, double scalar) {
        CondenserState r;

        r.steam_pressure = a.steam_pressure * scalar;
        r.air_pressure = a.air_pressure * scalar;
        r.steam_temperature = a.steam_temperature * scalar;
        r.steam_mass = a.steam_mass * scalar;
        r.air_mass = a.air_mass * scalar;
        r.hotwell_water_mass = a.hotwell_water_mass * scalar;
        r.hotwell_temperature = a.hotwell_temperature * scalar;
        r.condensate_flow_rate = a.condensate_flow_rate * scalar;

        return r;
    }

    CondenserState operator+(const CondenserState &a, const CondenserState &b) {
        CondenserState r = a;

        r.steam_pressure += b.steam_pressure;
        r.air_pressure += b.air_pressure;
        r.steam_temperature += b.steam_temperature;
        r.steam_mass += b.steam_mass;
        r.air_mass += b.air_mass;
        r.hotwell_water_mass += b.hotwell_water_mass;
        r.hotwell_temperature += b.hotwell_temperature;
        r.condensate_flow_rate += b.condensate_flow_rate;

        return r;
    }

} // namespace astara