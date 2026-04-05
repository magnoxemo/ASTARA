#ifndef ASTARA_STEAM_CONDENSER_H
#define ASTARA_STEAM_CONDENSER_H


namespace astara {

    struct CondenserState {
        double time = 0.0;
        double steam_pressure = 0.0;
        double air_pressure = 0.0;
        double steam_temperature = 0.0;
        double steam_mass = 0.0;
        double air_mass = 0.0;
        double hotwell_water_mass = 0.0;
        double hotwell_temperature = 0.0;
        double condensate_flow_rate = 0.0;
        double total_pressure = 0.0;
        double water_level = 0.0;
    };

    CondenserState operator+(const CondenserState &a, const CondenserState &b);
    CondenserState operator*(const CondenserState &a, double scalar);

    class Condenser {
    public:
        Condenser();
        Condenser(CondenserState& condenser_state);
        ~Condenser() = default;


        Condenser(const Condenser &) = delete;
        Condenser &operator=(const Condenser &) = delete;

        Condenser(Condenser &&) = default;
        Condenser &operator=(Condenser &&) = default;

        void timeStep(double dt);

        const CondenserState &getState() const { return _state; }
        CondenserState &getMutableState() { return _state; }

        void setTurbineInletFlow(double flow) { _turbine_inlet_flow = flow; }
        void setTurbineInletTemperature(double temp) { _turbine_inlet_temperature = temp; }
        void setColdWaterFlow(double flow) { _cold_water_flow = flow; }
        void setColdWaterTemperature(double temp) { _cold_water_temperature = temp; }
        void setVacuumBreakerFlow(double flow) { _vacuum_breaker_flow = flow; }
        void setOxygenBubblingFlow(double flow) { _oxygen_bubbling_flow = flow; }

        double getHotwellWaterOutletFlow() const { return _hotwell_water_outlet_flow; }
        double getHotwellTemperature() const { return _state.hotwell_temperature; }
        double getSteamPressure() const { return _state.steam_pressure; }
        double getTotalPressure() const { return _state.total_pressure; }
        double getWaterLevel() const { return _state.water_level; }
        double getCondensateFlow() const { return _state.condensate_flow_rate; }

        void setInitialSteamPressure(double pressure) {
            _state.steam_pressure = pressure;
        }
        void setInitialAirPressure(double pressure) {
            _state.air_pressure = pressure;
        }
        void setInitialSteamMass(double mass) {
            _state.steam_mass = mass;
        }
        void setInitialAirMass(double mass) {
            _state.air_mass = mass;
        }
        void setInitialHotwellWaterMass(double mass) {
            _state.hotwell_water_mass = mass;
        }
        void setInitialHotwellTemperature(double temp) {
            _state.hotwell_temperature = temp;
        }

        void setVolume(double volume) { _volume = volume; }
        void setHeatTransferUA(double ua) { _UA = ua; }
        void setHotwellArea(double area) { _hotwell_area = area; }
        void setWaterSpecificHeat(double cp) { _water_specific_heat = cp; }

    private:

        CondenserState _state;

        double _volume = 3.0;
        double _UA = 356.972e3;
        double _water_specific_heat = 4.2e3;
        double _steam_gas_constant = 0.4615e3;
        double _air_gas_constant = 0.287e3;

        double _turbine_inlet_flow = 4.0;
        double _turbine_inlet_temperature = 600.0 + 273.15;
        double _other_steam_inlet_flow = 0.0;
        double _other_steam_temperature = 400.0 + 273.15;
        double _steam_air_outlet_flow = 0.0;

        double _vacuum_breaker_flow = 0.0;
        double _drain_condenser_flow = 0.0;
        double _steam_gas_outlet_flow = 0.0;

        double _hotwell_area = 0.2;
        double _hotwell_water_outlet_flow = 110.0;
        double _oxygen_bubbling_flow = 10.0;

        double _cold_water_flow = 107.881;
        double _cold_water_temperature = 273.15 + 60.0;
        double _hot_water_outlet_temperature = 273.15 + 80.0;

        double calculateDSteamMassDt() const;
        double calculateDAirMassDt() const;
        double calculateDSteamPressureDt() const;
        double calculateDAirPressureDt() const;
        double calculateDHotwellMassDt() const;
        double calculateDHotwellTemperatureDt() const;

        double getSteamEnthalpy(double temperature, double pressure) const;
        double getWaterEnthalpy(double temperature, double pressure) const;
        double getSaturationTemperature(double pressure) const;
        double getWaterDensity(double temperature, double pressure) const;
        double getCondensationRate() const;

        CondenserState evaluateDerivatives(const CondenserState &s);
        void integrateRK4(double dt);
        void clampState();
        void updateDerivedQuantities();
    };

} // namespace astara

#endif //ASTARA_STEAM_CONDENSER_H
