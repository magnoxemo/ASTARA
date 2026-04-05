#include "Function.h"
#include "Reactor.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>

using namespace astara;

int main() {
    unsigned int n_groups = 6;

    std::vector<double> decay_constants = {0.0125, 0.0308, 0.114,
                                           0.307,  1.19,   3.19};
    std::vector<double> beta_values = {0.000209, 0.001414, 0.001309,
                                       0.002727, 0.000925, 0.000314};

    double total_beta = 0.0;
    for (double beta : beta_values) {
        total_beta += beta;
    }

    std::cout << "Total delayed neutron fraction: " << total_beta << std::endl;

    double Lambda = 1.79e-5;

    Reactor reactor(n_groups, decay_constants, beta_values, Lambda);

    reactor.setInitialPower(1.0);
    reactor.setInitialReactivity(0.0);
    reactor.setInitialFuelTemperature(565.0);
    reactor.setInitialModeratorTemperature(300.0);
    reactor.setInitialUpperPlenumTemperature(325.0);
    reactor.setInitialLowerPlenumTemperature(300.0);
    reactor.setInitialHotLegTemperature(325.0);
    reactor.setInitialColdLegTemperature(300.0);
    reactor.setFuelMass(222739.0 * 0.453592);
    reactor.setModeratorMass(50000.0);
    reactor.setFuelThermalCapacity(0.059 * 4186.8);
    reactor.setModeratorThermalCapacity(1.39 * 4186.8);
    reactor.setHeatTransferArea(59900.0 * 0.092903);
    reactor.setFissionEnergy(3.2e-11);
    reactor.setRatedPower(3436.0);
    reactor.setCoolantFlowRate(1.5e8 * 0.000125998);
    reactor.setCoolantInletTemperature(296.96);
    reactor.setUpperPlenumMass(38.96 * 700.0);
    reactor.setLowerPlenumMass(50.72 * 700.0);
    reactor.setHotLegMass(28.32 * 700.0);
    reactor.setColdLegMass(56.63 * 700.0);
    reactor.setFractionPowerInFuel(0.974);
    reactor.setControlRodEffectiveness(-0.001);
    reactor.setBoronEffectiveness(-1.0e-5);


    // Set temperature-dependent functions
    Function *h_func = new Function("5000 + 10 * (x - 300)", {"x"});
    Function *alpha_m_func = new Function("-3.6e-4 * (x - 300)", {"x"});
    Function *alpha_f_func = new Function("-1.98e-5 * (x - 600)", {"x"});
    Function *c_f_func = new Function("5000 + 0.2 * (x - 300)", {"x"});

    reactor.setHeatTransferCoefficient(h_func);
    reactor.setFuelSpecificHeatFunction(c_f_func);
    reactor.setFuelTemperatureCoEfficientFunction(alpha_f_func);
    reactor.setModeratorTemperatureCoEfficientFunction(alpha_m_func);

    for (unsigned int i = 0; i < n_groups; ++i) {
        double C_i = (beta_values[i] * reactor.getPower()) / (decay_constants[i] * Lambda);
        reactor.setInitialPrecursorConcentration(i, C_i);
    }


    // Open CSV file for writing
    std::ofstream csv_file("reactor_transient.csv");
    if (!csv_file.is_open()) {
        std::cerr << "Error: Could not open CSV file for writing" << std::endl;
        return 1;
    }

    reactor.writeCSVHeader(csv_file);

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "\nReactor Transient Analysis\n";
    std::cout << std::string(120, '=') << std::endl;
    std::cout << std::setw(12) << "Time (s)"
              << std::setw(14) << "Power (MW)"
              << std::setw(14) << "Rho (pcm)"
              << std::setw(14) << "T_fuel (C)"
              << std::setw(14) << "T_mod (C)"
              << std::setw(14) << "T_up (C)"
              << std::setw(14) << "T_lp (C)"
              << std::setw(14) << "T_hot (C)"
              << std::setw(14) << "T_cold (C)" << std::endl;
    std::cout << std::string(120, '=') << std::endl;


    reactor.recordState(csv_file);

    auto print_state=[&reactor] (double time){
            std::cout << std::setw(12) << std::fixed << std::setprecision(2) << time
          << std::setw(14) << std::setprecision(2)
          << reactor.getPower() * reactor.getRatedPower()
          << std::setw(14) << std::setprecision(2) << reactor.getReactivity() * 1e5
          << std::setw(14) << std::setprecision(2) << reactor.getFuelTemperature()
          << std::setw(14) << std::setprecision(2) << reactor.getModeratorTemperature()
          << std::setw(14) << std::setprecision(2) << reactor.getUpperPlenumTemperature()
          << std::setw(14) << std::setprecision(2) << reactor.getLowerPlenumTemperature()
          << std::setw(14) << std::setprecision(2) << reactor.getHotLegTemperature()
          << std::setw(14) << std::setprecision(2) << reactor.getColdLegTemperature()
          << std::endl;
    };

    double dt = 0.01;
    double t_end = 1000.0;
    int steps = static_cast<int>(t_end / dt);
    int output_freq = 100;
    int csv_freq = 10;

    double insert_time = 200.0;
    double withdraw_time = 400.0;

    bool rod_inserted = false;
    bool rod_withdrawn = false;

    print_state(0);

    for (int i = 1; i <= steps; ++i) {
        double current_time = i * dt;

        if (!rod_inserted && current_time >= insert_time) {
            reactor.setCoolantInletTemperature(340);
            reactor.setHeatTransferArea(59900);
            rod_inserted = true;
            std::cout << "\n--- Rod inserted at t = " << current_time << " s ---\n" << std::endl;
        }

        if (!rod_withdrawn && current_time >= withdraw_time) {
            reactor.setCoolantInletTemperature(3000);
            reactor.setHeatTransferArea(200);
            reactor.setInitialReactivity(0.00001);
            rod_withdrawn = false;
        }

        reactor.timeStep(dt);

        if (i % csv_freq == 0) {
            reactor.recordState(csv_file);
        }

        if (i % output_freq == 0)
            print_state(current_time);
    }

    csv_file.close();

    std::cout << "\n" << std::string(120, '=') << std::endl;
    std::cout << "Simulation Complete\n";
    std::cout << "Data saved to reactor_transient.csv\n";
    std::cout << std::string(120, '=') << std::endl;

    return 0;
}