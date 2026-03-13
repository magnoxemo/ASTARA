#ifndef ASTARA_STEAM_GENERATOR_H
#define ASTARA_STEAM_GENERATOR_H


#include <vector>

namespace astara {

    struct SteamGeneratorState {
        double time = 0.0;
        double primary_inlet_temperature = 0.0;
        double primary_outlet_temperature = 0.0;
        double primary_node1_temperature = 0.0;
        double primary_node2_temperature = 0.0;
        double metal_node1_temperature = 0.0;
        double metal_node2_temperature = 0.0;
        double secondary_temperature = 0.0;
        double secondary_pressure = 0.0;
        double secondary_water_mass = 0.0;
        double secondary_steam_mass = 0.0;
        double steam_outlet_flow_rate = 0.0;
    };


    SteamGeneratorState operator+(const SteamGeneratorState &a, const SteamGeneratorState &b);
    SteamGeneratorState operator*(const SteamGeneratorState &a, double scalar);

    class SteamGenerator {

    public:
        SteamGenerator();
        ~SteamGenerator() = default;

        SteamGenerator(const SteamGenerator &) = delete;
        SteamGenerator &operator=(const SteamGenerator &) = delete;

        SteamGenerator(SteamGenerator &&) = default;
        SteamGenerator &operator=(SteamGenerator &&) = default;

        void timeStep(double dt);

        const SteamGeneratorState &getState() const { return _state; }
        SteamGeneratorState &getMutableState() { return _state; }

        void setPrimaryInletTemperature(double temp) {
            _state.primary_inlet_temperature = temp;
        }
        void setPrimaryMassFlowRate(double flow) {
            _primary_mass_flow_rate = flow;
        }
        void setFeedwaterTemperature(double temp) {
            _feedwater_temperature = temp;
        }
        void setFeedwaterFlowRate(double flow) {
            _feedwater_flow_rate = flow;
        }
        void setSteamValveCoefficient(double coeff) {
            _steam_valve_coefficient = coeff;
        }

        double getPrimaryOutletTemperature() const {
            return _state.primary_outlet_temperature;
        }
        double getSecondaryPressure() const {
            return _state.secondary_pressure;
        }
        double getSecondaryTemperature() const {
            return _state.secondary_temperature;
        }
        double getSteamOutletFlowRate() const {
            return _state.steam_outlet_flow_rate;
        }

        void setInitialPrimaryNode1Temperature(double temp) {
            _state.primary_node1_temperature = temp;
        }
        void setInitialPrimaryNode2Temperature(double temp) {
            _state.primary_node2_temperature = temp;
        }
        void setInitialMetalNode1Temperature(double temp) {
            _state.metal_node1_temperature = temp;
        }
        void setInitialMetalNode2Temperature(double temp) {
            _state.metal_node2_temperature = temp;
        }
        void setInitialSecondaryPressure(double pressure) {
            _state.secondary_pressure = pressure;
        }
        void setInitialSecondaryWaterMass(double mass) {
            _state.secondary_water_mass = mass;
        }
        void setInitialSecondarySteamMass(double mass) {
            _state.secondary_steam_mass = mass;
        }

        void setPrimaryNode1Mass(double mass) { _primary_node1_mass = mass; }
        void setPrimaryNode2Mass(double mass) { _primary_node2_mass = mass; }
        void setPrimarySpecificHeat(double cp) { _primary_specific_heat = cp; }

        void setMetalNode1Mass(double mass) { _metal_node1_mass = mass; }
        void setMetalNode2Mass(double mass) { _metal_node2_mass = mass; }
        void setMetalSpecificHeat(double cp) { _metal_specific_heat = cp; }

        void setHeatTransferAreaPM1(double area) { _heat_transfer_area_pm1 = area; }
        void setHeatTransferAreaPM2(double area) { _heat_transfer_area_pm2 = area; }
        void setHeatTransferAreaMS1(double area) { _heat_transfer_area_ms1 = area; }
        void setHeatTransferAreaMS2(double area) { _heat_transfer_area_ms2 = area; }

        void setHeatTransferCoefficientPM(double h) { _heat_transfer_coeff_pm = h; }
        void setHeatTransferCoefficientMS(double h) { _heat_transfer_coeff_ms = h; }

        void setSecondaryVolume(double volume) { _secondary_volume = volume; }

    private:
        SteamGeneratorState _state;


        // for now I am going with magic numbers
        // I need to come back and implement methods to calcualte these
        // to make it more generalized
        double _primary_mass_flow_rate = 17700.0;
        double _feedwater_flow_rate = 2164.0;
        double _feedwater_temperature = 232.2;
        double _steam_valve_coefficient = 2.0;

        double _primary_node1_mass = 3000.0;
        double _primary_node2_mass = 3000.0;
        double _primary_specific_heat = 5400.0;

        double _metal_node1_mass = 50000.0;
        double _metal_node2_mass = 50000.0;
        double _metal_specific_heat = 500.0;

        double _heat_transfer_area_pm1 = 5000.0;
        double _heat_transfer_area_pm2 = 5000.0;
        double _heat_transfer_area_ms1 = 5000.0;
        double _heat_transfer_area_ms2 = 5000.0;

        double _heat_transfer_coeff_pm = 10000.0;
        double _heat_transfer_coeff_ms = 5000.0;

        double _secondary_volume = 100.0;

        double calculateDTp1Dt() const;
        double calculateDTp2Dt() const;
        double calculateDTm1Dt() const;
        double calculateDTm2Dt() const;
        double calculateDPsDt() const;
        double calculateDMwsDt() const;
        double calculateDMssDt() const;

        double getSaturationTemperature(double pressure) const;
        double getWaterEnthalpy(double temperature, double pressure) const;
        double getSteamEnthalpy(double pressure) const;
        double getWaterSpecificVolume(double temperature, double pressure) const;
        double getSteamSpecificVolume(double pressure) const;
        double getWaterDensity(double temperature, double pressure) const;
        double getSteamDensity(double pressure) const;

        SteamGeneratorState evaluateDerivatives(const SteamGeneratorState &s);
        void integrateRK4(double dt);
        void clampState();
        void updateDerivedQuantities();
    };

} // namespace astara


#endif //ASTARA_STEAM_GENERATOR_H
