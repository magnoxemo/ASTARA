#ifndef ASTARA_REACTOR_H
#define ASTARA_REACTOR_H

#include <vector>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace astara {

    class Function;

    struct ReactorState {
        double time = 0.0;
        double reactivity = 0.0;
        double power = 0.0;
        double fuel_temperature = 0.0;
        double moderator_temperature = 0.0;

        double upper_plenum_temperature = 0.0;
        double lower_plenum_temperature = 0.0;
        double hot_leg_temperature = 0.0;
        double cold_leg_temperature = 0.0;

        std::vector<double> fuel_temperatures;      // 3 nodes
        std::vector<double> coolant_temperatures;   // 6 nodes
        std::vector<double> precursor_concentrations;
    };

    ReactorState operator+(const ReactorState &a, const ReactorState &b);
    ReactorState operator*(const ReactorState &a, double scalar);

    class Reactor {
    public:
        Reactor(unsigned int n_groups,
                const std::vector<double> &neutron_group_const,
                const std::vector<double> &delayed_neutron_constants,
                double neutron_generation_time);
        ~Reactor() = default;

        Reactor(const Reactor &) = delete;
        Reactor &operator=(const Reactor &) = delete;

        Reactor(Reactor &&) = default;
        Reactor &operator=(Reactor &&) = default;

        void timeStep(double dt);

        const ReactorState &getState() const { return _state; }
        ReactorState &getMutableState() { return _state; }

        double getPower() const { return _state.power; }
        double getReactivity() const { return _state.reactivity; }
        double getFuelTemperature() const { return _state.fuel_temperature; }
        double getModeratorTemperature() const { return _state.moderator_temperature; }

        // New getters
        double getRatedPower()const {return _state.power;}
        double getUpperPlenumTemperature() const { return _state.upper_plenum_temperature; }
        double getLowerPlenumTemperature() const { return _state.lower_plenum_temperature; }
        double getHotLegTemperature() const { return _state.hot_leg_temperature; }
        double getColdLegTemperature() const { return _state.cold_leg_temperature; }
        const std::vector<double>& getFuelTemperatures() const { return _state.fuel_temperatures; }
        const std::vector<double>& getCoolantTemperatures() const { return _state.coolant_temperatures; }

        void setInitialPower(double power) { _state.power = power; }
        void setInitialReactivity(double reactivity) { _state.reactivity = reactivity; }
        void setInitialFuelTemperature(double temp) {
            _state.fuel_temperature = temp;
            // Initialize all fuel nodes to same temperature
            for (double& t : _state.fuel_temperatures) {
                t = temp;
            }
        }
        void setInitialModeratorTemperature(double temp) {
            _state.moderator_temperature = temp;
            // Initialize all coolant nodes to same temperature
            for (double& t : _state.coolant_temperatures) {
                t = temp;
            }
            _state.upper_plenum_temperature = temp;
            _state.lower_plenum_temperature = temp;
            _state.hot_leg_temperature = temp;
            _state.cold_leg_temperature = temp;
        }

        void setInitialUpperPlenumTemperature(double temp) { _state.upper_plenum_temperature = temp; }
        void setInitialLowerPlenumTemperature(double temp) { _state.lower_plenum_temperature = temp; }
        void setInitialHotLegTemperature(double temp) { _state.hot_leg_temperature = temp; }
        void setInitialColdLegTemperature(double temp) { _state.cold_leg_temperature = temp; }

        void setInitialPrecursorConcentration(unsigned int group, double value) {
            if (group < _state.precursor_concentrations.size()) {
                _state.precursor_concentrations[group] = value;
            }
        }

        void setFuelMass(double mass);
        void setModeratorMass(double mass);
        void setFuelThermalCapacity(double cp) { _thermal_capacity_fuel = cp; }
        void setModeratorThermalCapacity(double cp) { _thermal_capacity_moderator = cp; }
        void setHeatTransferArea(double area);
        void setFissionEnergy(double energy) { _fission_energy = energy; }
        void setRatedPower(double power) { _rated_thermal_power_MW = power; }
        void setCoolantFlowRate(double flow_rate);
        void setCoolantInletTemperature(double temp) { _coolant_inlet_temperature = temp; }

        // New setters for plenum and leg masses
        void setUpperPlenumMass(double mass);
        void setLowerPlenumMass(double mass);
        void setHotLegMass(double mass);
        void setColdLegMass(double mass);
        void setFractionPowerInFuel(double fraction);

        void setControlRodEffectiveness(double worth) { _control_rod_worth = worth; }
        void setBoronEffectiveness(double worth) { _boron_worth = worth; }

        void setHeatTransferCoefficient(Function *h) { _convective_heat_transfer_coefficient = h; }
        void setFuelSpecificHeatFunction(Function *c_f) { _fuel_specific_heat = c_f; }
        void setFuelTemperatureCoEfficientFunction(Function *alpha_f) { _fuel_temperature_coefficient = alpha_f; }
        void setModeratorTemperatureCoEfficientFunction(Function *alpha_m) { _moderator_temperature_coefficient = alpha_m; }
        void setBoronTemperatureCoEfficientFunction(Function *alpha_b) { _boron_temperature_coefficient = alpha_b; }

        void insertControlRod(double length);
        void injectBoron(double boron_concentration);

        static Function *createConstantFunction(double value);

        void printStates() const;

        void recordState(std::ostream& csv_file);
        void writeCSVHeader(std::ostream& csv_file) ;

    private:
        unsigned int _number_of_neutron_groups;
        std::vector<double> _neutron_group_decay_constants;
        std::vector<double> _delayed_neutron_fractions;
        double _neutron_generation_time;
        double _total_delayed_neutron_fraction;
        double _prompt_neutron_fraction;

        ReactorState _state;

        // Reactor thermal parameters
        double _fuel_mass = 0.0;
        double _fuel_mass_per_node = 0.0;
        double _moderator_mass = 0.0;
        double _coolant_mass_per_node = 0.0;
        double _thermal_capacity_fuel = 0.0;
        double _thermal_capacity_moderator = 0.0;
        double _heat_transfer_area = 0.0;
        double _heat_transfer_area_per_node = 0.0;
        double _fission_energy = 0.0;
        double _rated_thermal_power_MW = 0.0;
        double _coolant_mass_flow_rate = 0.0;
        double _coolant_mass_flow_rate_per_node = 0.0;
        double _coolant_inlet_temperature = 0.0;

        // New parameters for plenum and legs
        double _upper_plenum_mass = 0.0;
        double _lower_plenum_mass = 0.0;
        double _hot_leg_mass = 0.0;
        double _cold_leg_mass = 0.0;
        double _fraction_power_in_fuel = 0.97;  // Default from thesis

        // Reactivity coefficients
        double _control_rod_worth = 0.0;
        double _boron_worth = 0.0;

        // Temperature-dependent functions
        Function *_convective_heat_transfer_coefficient = nullptr;
        Function *_fuel_specific_heat = nullptr;
        Function *_fuel_temperature_coefficient = nullptr;
        Function *_moderator_temperature_coefficient = nullptr;
        Function *_boron_temperature_coefficient = nullptr;

        double calculateDRhoDt() const;
        double calculateTotalReactivity() const;
        double calculateDPowerDt() const;
        double calculateFuelNodeDt(int node) const;
        double calculateCoolantNodeDt(int node) const;
        double calculateDUpperPlenumDt() const;
        double calculateDLowerPlenumDt() const;
        double calculateDHotLegDt() const;
        double calculateDColdLegDt() const;
        double calculateDFuelTempDt() const;
        double calculateDModeratorTempDt() const;
        double calculateDCDt(unsigned int group_index) const;

        ReactorState evaluateDerivatives(const ReactorState &s);
        void integrateRK4(double dt);
        void clampState();
    };

}

#endif // ASTARA_REACTOR_H