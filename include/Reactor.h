// Reactor.h
#ifndef REACTOR_H
#define REACTOR_H

#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

namespace astara {

    class Function;

    struct ReactorState {
        double time = 0.0;
        double power = 0.0;
        double reactivity = 0.0;
        double fuel_temperature = 0.0;
        double moderator_temperature = 0.0;
        std::vector<double> precursor_concentrations;
    };

    ReactorState operator+(const ReactorState &a, const ReactorState &b);
    ReactorState operator*(const ReactorState &a, double scalar);

    class Reactor {
    public:
        Reactor(unsigned int n_groups, const std::vector<double> &neutron_group_const,
                const std::vector<double> &delayed_neutron_constants,
                double neutron_generation_time);

        ~Reactor() = default;

        Reactor(const Reactor &) = delete;
        Reactor &operator=(const Reactor &) = delete;

        Reactor(Reactor &&) = default;
        Reactor &operator=(Reactor &&) = default;

        void timeStep(double dt);

        ReactorState stateSnapshot() const { return _state; }
        void restoreState(const ReactorState &snapshot) { _state = snapshot; }

        double getPower() const { return _state.power; }
        double getReactivity() const { return _state.reactivity; }
        double getFuelTemperature() const { return _state.fuel_temperature; }
        double getModeratorTemperature() const {
            return _state.moderator_temperature;
        }
        const ReactorState &getState() const { return _state; }
        ReactorState &getMutableState() { return _state; }

        double getFuelMass() const { return _fuel_mass; }
        double getModeratorMass() const { return _moderator_mass; }
        double getFuelThermalCapacity() const { return _thermal_capacity_fuel; }
        double getModeratorThermalCapacity() const {
            return _thermal_capacity_moderator;
        }
        double getHeatTransferArea() const { return _heat_transfer_area; }
        double getFissionEnergy() const { return _fission_energy; }
        double getControlRodEffectiveness() const { return _control_rod_worth; }
        double getBoronEffectiveness() const { return _boron_worth; }
        double getRatedPower() const { return _rated_thermal_power_MW; }
        double getCoolantFlowRate() const { return _coolant_mass_flow_rate; }
        double getCoolantInletTemperature() const { return _coolant_inlet_temperature; }

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
        void setRatedPower(double power_MW) { _rated_thermal_power_MW = power_MW; }
        void setCoolantFlowRate(double flow_rate_kg_s) {
            _coolant_mass_flow_rate = flow_rate_kg_s;
        }
        void setCoolantInletTemperature(double temp_C) {
            _coolant_inlet_temperature = temp_C;
        }

        void setControlRodEffectiveness(double effectiveness) {
            _control_rod_worth = effectiveness;
        }
        void setBoronEffectiveness(double effectiveness) {
            _boron_worth = effectiveness;
        }

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

        static Function *createConstantFunction(double value);

        void insertControlRod(double length);
        void injectBoron(double boron_concentration);

        void print_states() {
            std::cout << std::setw(25) << this->getState().time << std::setw(14)
                      << this->getPower() << std::setw(14)
                      << this->getReactivity() * 1e5 << std::setw(14)
                      << this->getFuelTemperature() << std::setw(14)
                      << this->getModeratorTemperature() << std::endl;
        };

    private:
        unsigned int _number_of_neutron_groups;
        std::vector<double> _neutron_group_decay_constants;
        std::vector<double> _delayed_neutron_fractions;
        double _neutron_generation_time;
        double _total_delayed_neutron_fraction;
        double _prompt_neutron_fraction;

        double _fuel_mass = 0.0;
        double _moderator_mass = 0.0;
        double _thermal_capacity_fuel = 0.0;
        double _thermal_capacity_moderator = 0.0;
        double _heat_transfer_area = 0.0;
        double _fission_energy = 3.2e-11;

        double _control_rod_worth = 0.0;
        double _boron_worth = 0.0;

        double _rated_thermal_power_MW = 3400.0;
        double _coolant_mass_flow_rate = 17700.0;
        double _coolant_inlet_temperature = 296.96;

        std::unique_ptr<Function> _convective_heat_transfer_coefficient;
        std::unique_ptr<Function> _fuel_specific_heat;
        std::unique_ptr<Function> _fuel_temperature_coefficient;
        std::unique_ptr<Function> _moderator_temperature_coefficient;
        std::unique_ptr<Function> _boron_temperature_coefficient;

        ReactorState _state;

        double calculateDRhoDt() const;
        double calculateDPowerDt() const;
        double calculateDFuelTempDt() const;
        double calculateDModeratorTempDt() const;
        double calculateDCDt(unsigned int group_index) const;
        double calculateTotalReactivity() const;

        ReactorState evaluateDerivatives(const ReactorState &s);

        void clampState();

        void integrateRK4(double dt);
    };

} // namespace astara

#endif // REACTOR_H