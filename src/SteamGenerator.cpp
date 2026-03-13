#include "SteamGenerator.h"
#include "CoolProp.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace astara {

    SteamGenerator::SteamGenerator() {
        _state.primary_node1_temperature = 320.0;
        _state.primary_node2_temperature = 310.0;
        _state.metal_node1_temperature = 300.0;
        _state.metal_node2_temperature = 295.0;
        _state.secondary_pressure = 7.0e6;
        _state.secondary_temperature = getSaturationTemperature(_state.secondary_pressure);
        _state.secondary_water_mass = 50000.0;
        _state.secondary_steam_mass = 2000.0;
    }

    void SteamGenerator::timeStep(double dt) {
        if (dt <= 0.0)
            throw std::runtime_error("Time step must be positive");

        integrateRK4(dt);
        _state.time += dt;
        updateDerivedQuantities();
    }

    double SteamGenerator::calculateDTp1Dt() const {
        const double mdot = _primary_mass_flow_rate;
        const double cp = _primary_specific_heat;
        const double m1 = _primary_node1_mass;

        const double T_in = _state.primary_inlet_temperature;
        const double T_p1 = _state.primary_node1_temperature;
        const double T_m1 = _state.metal_node1_temperature;

        const double h_pm = _heat_transfer_coeff_pm;
        const double A_pm1 = _heat_transfer_area_pm1;

        const double flow_term = (mdot * cp / m1) * (T_in - T_p1);
        const double heat_term = (h_pm * A_pm1 / (m1 * cp)) * (T_p1 - T_m1);

        return flow_term - heat_term;
    }

    double SteamGenerator::calculateDTp2Dt() const {
        const double mdot = _primary_mass_flow_rate;
        const double cp = _primary_specific_heat;
        const double m2 = _primary_node2_mass;

        const double T_p1 = _state.primary_node1_temperature;
        const double T_p2 = _state.primary_node2_temperature;
        const double T_m2 = _state.metal_node2_temperature;

        const double h_pm = _heat_transfer_coeff_pm;
        const double A_pm2 = _heat_transfer_area_pm2;

        const double flow_term = (mdot * cp / m2) * (T_p1 - T_p2);
        const double heat_term = (h_pm * A_pm2 / (m2 * cp)) * (T_p2 - T_m2);

        return flow_term - heat_term;
    }

    double SteamGenerator::calculateDTm1Dt() const {
        const double m_m1 = _metal_node1_mass;
        const double cp_m = _metal_specific_heat;

        const double T_p1 = _state.primary_node1_temperature;
        const double T_m1 = _state.metal_node1_temperature;
        const double T_s = _state.secondary_temperature;

        const double h_pm = _heat_transfer_coeff_pm;
        const double h_ms = _heat_transfer_coeff_ms;
        const double A_pm1 = _heat_transfer_area_pm1;
        const double A_ms1 = _heat_transfer_area_ms1;

        const double heat_from_primary = (h_pm * A_pm1 / (m_m1 * cp_m)) * (T_p1 - T_m1);
        const double heat_to_secondary = (h_ms * A_ms1 / (m_m1 * cp_m)) * (T_m1 - T_s);

        return heat_from_primary - heat_to_secondary;
    }

    double SteamGenerator::calculateDTm2Dt() const {
        const double m_m2 = _metal_node2_mass;
        const double cp_m = _metal_specific_heat;

        const double T_p2 = _state.primary_node2_temperature;
        const double T_m2 = _state.metal_node2_temperature;
        const double T_s = _state.secondary_temperature;

        const double h_pm = _heat_transfer_coeff_pm;
        const double h_ms = _heat_transfer_coeff_ms;
        const double A_pm2 = _heat_transfer_area_pm2;
        const double A_ms2 = _heat_transfer_area_ms2;

        const double heat_from_primary = (h_pm * A_pm2 / (m_m2 * cp_m)) * (T_p2 - T_m2);
        const double heat_to_secondary = (h_ms * A_ms2 / (m_m2 * cp_m)) * (T_m2 - T_s);

        return heat_from_primary - heat_to_secondary;
    }

    double SteamGenerator::calculateDPsDt() const {
        const double p_s = _state.secondary_pressure;
        const double T_s = _state.secondary_temperature;
        const double m_ws = _state.secondary_water_mass;
        const double m_ss = _state.secondary_steam_mass;

        const double T_m1 = _state.metal_node1_temperature;
        const double T_m2 = _state.metal_node2_temperature;

        const double h_ms = _heat_transfer_coeff_ms;
        const double A_ms1 = _heat_transfer_area_ms1;
        const double A_ms2 = _heat_transfer_area_ms2;

        const double Q1 = h_ms * A_ms1 * (T_m1 - T_s);
        const double Q2 = h_ms * A_ms2 * (T_m2 - T_s);
        const double Q_total = Q1 + Q2;

        const double mdot_fw = _feedwater_flow_rate;
        const double T_fw = _feedwater_temperature;
        const double mdot_steam = _state.steam_outlet_flow_rate;

        const double h_fw = getWaterEnthalpy(T_fw, p_s);
        const double h_s = getSteamEnthalpy(p_s);
        const double h_w = getWaterEnthalpy(T_s, p_s);

        const double dh_w_dp = 1.0 / getWaterDensity(T_s, p_s);
        const double dh_s_dp = 1.0 / getSteamDensity(p_s);

        const double nu_w = getWaterSpecificVolume(T_s, p_s);
        const double nu_s = getSteamSpecificVolume(p_s);
        const double dnu_s_dp = -nu_s * nu_s / (p_s * 1e-6);

        const double K_s = m_ws * dh_w_dp + m_ss * dh_s_dp -
                           m_ws * ((h_w - h_s) / (nu_w - nu_s)) * dnu_s_dp;

        if (std::abs(K_s) < 1e-10) {
            return 0.0;
        }

        const double numerator = Q_total - mdot_steam * (h_s - h_fw);

        return numerator / K_s;
    }

    double SteamGenerator::calculateDMwsDt() const {
        const double mdot_fw = _feedwater_flow_rate;
        const double mdot_steam = _state.steam_outlet_flow_rate;

        const double T_s = _state.secondary_temperature;
        const double p_s = _state.secondary_pressure;

        const double T_m1 = _state.metal_node1_temperature;
        const double T_m2 = _state.metal_node2_temperature;

        const double h_ms = _heat_transfer_coeff_ms;
        const double A_ms1 = _heat_transfer_area_ms1;
        const double A_ms2 = _heat_transfer_area_ms2;

        const double Q1 = h_ms * A_ms1 * (T_m1 - T_s);
        const double Q2 = h_ms * A_ms2 * (T_m2 - T_s);

        const double h_fg = getSteamEnthalpy(p_s) - getWaterEnthalpy(T_s, p_s);

        const double evaporation_rate = (Q1 + Q2) / h_fg;

        return mdot_fw - evaporation_rate;
    }

    double SteamGenerator::calculateDMssDt() const {
        const double mdot_steam = _state.steam_outlet_flow_rate;

        const double T_s = _state.secondary_temperature;
        const double p_s = _state.secondary_pressure;

        const double T_m1 = _state.metal_node1_temperature;
        const double T_m2 = _state.metal_node2_temperature;

        const double h_ms = _heat_transfer_coeff_ms;
        const double A_ms1 = _heat_transfer_area_ms1;
        const double A_ms2 = _heat_transfer_area_ms2;

        const double Q1 = h_ms * A_ms1 * (T_m1 - T_s);
        const double Q2 = h_ms * A_ms2 * (T_m2 - T_s);

        const double h_fg = getSteamEnthalpy(p_s) - getWaterEnthalpy(T_s, p_s);

        const double evaporation_rate = (Q1 + Q2) / h_fg;

        return evaporation_rate - mdot_steam;
    }

    double SteamGenerator::getSaturationTemperature(double pressure) const {
        try {
            double T_sat = CoolProp::PropsSI("T", "P", pressure, "Q", 0.0, "Water");
            return T_sat - 273.15;
        } catch (...) {
            return 285.0;
        }
    }

    double SteamGenerator::getWaterEnthalpy(double temperature, double pressure) const {
        try {
            double h = CoolProp::PropsSI("H", "T", temperature + 273.15, "P", pressure, "Water");
            return h;
        } catch (...) {
            return 1.0e6;
        }
    }

    double SteamGenerator::getSteamEnthalpy(double pressure) const {
        try {
            double h = CoolProp::PropsSI("H", "P", pressure, "Q", 1.0, "Water");
            return h;
        } catch (...) {
            return 2.7e6;
        }
    }

    double SteamGenerator::getWaterSpecificVolume(double temperature, double pressure) const {
        try {
            double rho = CoolProp::PropsSI("D", "T", temperature + 273.15, "P", pressure, "Water");
            return 1.0 / rho;
        } catch (...) {
            return 1.0e-3;
        }
    }

    double SteamGenerator::getSteamSpecificVolume(double pressure) const {
        try {
            double rho = CoolProp::PropsSI("D", "P", pressure, "Q", 1.0, "Water");
            return 1.0 / rho;
        } catch (...) {
            return 0.01;
        }
    }

    double SteamGenerator::getWaterDensity(double temperature, double pressure) const {
        try {
            return CoolProp::PropsSI("D", "T", temperature + 273.15, "P", pressure, "Water");
        } catch (...) {
            return 700.0;
        }
    }

    double SteamGenerator::getSteamDensity(double pressure) const {
        try {
            return CoolProp::PropsSI("D", "P", pressure, "Q", 1.0, "Water");
        } catch (...) {
            return 40.0;
        }
    }

    void SteamGenerator::integrateRK4(double dt) {
        SteamGeneratorState y0 = _state;

        SteamGeneratorState k1 = evaluateDerivatives(y0);
        SteamGeneratorState k2 = evaluateDerivatives(y0 + k1 * (dt * 0.5));
        SteamGeneratorState k3 = evaluateDerivatives(y0 + k2 * (dt * 0.5));
        SteamGeneratorState k4 = evaluateDerivatives(y0 + k3 * dt);

        _state = y0 + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (dt / 6.0);

        clampState();
    }

    SteamGeneratorState SteamGenerator::evaluateDerivatives(const SteamGeneratorState &s) {
        SteamGeneratorState saved = _state;
        _state = s;
        _state.secondary_temperature = getSaturationTemperature(_state.secondary_pressure);
        _state.steam_outlet_flow_rate = _steam_valve_coefficient * _state.secondary_pressure / 1e6;

        SteamGeneratorState k;
        k.primary_node1_temperature = calculateDTp1Dt();
        k.primary_node2_temperature = calculateDTp2Dt();
        k.metal_node1_temperature = calculateDTm1Dt();
        k.metal_node2_temperature = calculateDTm2Dt();
        k.secondary_pressure = calculateDPsDt();
        k.secondary_water_mass = calculateDMwsDt();
        k.secondary_steam_mass = calculateDMssDt();

        _state = saved;
        return k;
    }

    void SteamGenerator::clampState() {
        if (_state.secondary_pressure < 1.0e5) {
            _state.secondary_pressure = 1.0e5;
        }
        if (_state.secondary_pressure > 20.0e6) {
            _state.secondary_pressure = 20.0e6;
        }

        if (_state.secondary_water_mass < 1000.0) {
            _state.secondary_water_mass = 1000.0;
        }

        if (_state.secondary_steam_mass < 100.0) {
            _state.secondary_steam_mass = 100.0;
        }

        if (_state.primary_node1_temperature < -273.15) {
            _state.primary_node1_temperature = -273.15;
        }
        if (_state.primary_node2_temperature < -273.15) {
            _state.primary_node2_temperature = -273.15;
        }
        if (_state.metal_node1_temperature < -273.15) {
            _state.metal_node1_temperature = -273.15;
        }
        if (_state.metal_node2_temperature < -273.15) {
            _state.metal_node2_temperature = -273.15;
        }
    }

    void SteamGenerator::updateDerivedQuantities() {
        _state.secondary_temperature = getSaturationTemperature(_state.secondary_pressure);
        _state.primary_outlet_temperature = _state.primary_node2_temperature;
        _state.steam_outlet_flow_rate = _steam_valve_coefficient * _state.secondary_pressure / 1e6;
    }

    SteamGeneratorState operator*(const SteamGeneratorState &a, double scalar) {
        SteamGeneratorState r;

        r.primary_node1_temperature = a.primary_node1_temperature * scalar;
        r.primary_node2_temperature = a.primary_node2_temperature * scalar;
        r.metal_node1_temperature = a.metal_node1_temperature * scalar;
        r.metal_node2_temperature = a.metal_node2_temperature * scalar;
        r.secondary_temperature = a.secondary_temperature * scalar;
        r.secondary_pressure = a.secondary_pressure * scalar;
        r.secondary_water_mass = a.secondary_water_mass * scalar;
        r.secondary_steam_mass = a.secondary_steam_mass * scalar;
        r.steam_outlet_flow_rate = a.steam_outlet_flow_rate * scalar;

        return r;
    }

    SteamGeneratorState operator+(const SteamGeneratorState &a, const SteamGeneratorState &b) {
        SteamGeneratorState r = a;

        r.primary_node1_temperature += b.primary_node1_temperature;
        r.primary_node2_temperature += b.primary_node2_temperature;
        r.metal_node1_temperature += b.metal_node1_temperature;
        r.metal_node2_temperature += b.metal_node2_temperature;
        r.secondary_temperature += b.secondary_temperature;
        r.secondary_pressure += b.secondary_pressure;
        r.secondary_water_mass += b.secondary_water_mass;
        r.secondary_steam_mass += b.secondary_steam_mass;
        r.steam_outlet_flow_rate += b.steam_outlet_flow_rate;

        return r;
    }

} // namespace astara