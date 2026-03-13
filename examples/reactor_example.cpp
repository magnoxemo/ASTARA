// main.cpp
#include "Function.h"
#include "Reactor.h"
#include <cmath>
#include <iomanip>
#include <iostream>

using namespace astara;

int main() {

  unsigned int n_groups = 6;

  std::vector<double> decay_constants = {0.0127, 0.0317, 0.115,
                                         0.311,  1.40,   3.87};
  std::vector<double> beta_values = {0.000215, 0.001424, 0.001274,
                                     0.002568, 0.000748, 0.000273};

  double total_beta = 0.0;
  for (double beta : beta_values) {
    total_beta += beta;
  }

  std::cout << "Total delayed neutron fraction: " << total_beta << std::endl;

  double Lambda = 5e-5;

  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);

  reactor.setInitialPower(1.0);
  reactor.setInitialReactivity(0.0);
  reactor.setInitialFuelTemperature(565.0);
  reactor.setInitialModeratorTemperature(300.0);

  for (unsigned int i = 0; i < n_groups; ++i) {
    double C_i =
            (beta_values[i] * reactor.getPower()) / (decay_constants[i] * Lambda);
    reactor.setInitialPrecursorConcentration(i, C_i);
    std::cout << "Precursor " << i << " initial concentration: " << C_i
              << std::endl;
  }

  reactor.setFuelMass(2000.0);
  reactor.setModeratorMass(50000.0);
  reactor.setFuelThermalCapacity(500.0);
  reactor.setModeratorThermalCapacity(4186.0);
  reactor.setHeatTransferArea(500.0);
  reactor.setFissionEnergy(3.2e-11);
  reactor.setRatedPower(100.0);
  reactor.setCoolantFlowRate(17700.0);
  reactor.setCoolantInletTemperature(296.96);

  reactor.setControlRodEffectiveness(-0.0001);
  reactor.setBoronEffectiveness(-1.0e-5);

  Function *h_func = new Function("5000 + 10 * (x - 300)", {"x"});
  reactor.setHeatTransferCoefficient(h_func);

  Function *c_f_func = new Function("500 + 0.2 * (x - 300)", {"x"});
  reactor.setFuelSpecificHeatFunction(c_f_func);

  Function *alpha_f_func = new Function("-0.00001 * (x - 565)", {"x"});
  reactor.setFuelTemperatureCoEfficientFunction(alpha_f_func);

  Function *alpha_m_func = new Function("-0.00005 * (x - 300)", {"x"});
  reactor.setModeratorTemperatureCoEfficientFunction(alpha_m_func);

  Function *alpha_boron_func = new Function("-0.005 * (x - 300)", {"x"});
  reactor.setBoronTemperatureCoEfficientFunction(alpha_boron_func);

  std::cout << std::fixed << std::setprecision(6);
  std::cout << "\nReactor Transient Analysis\n";
  std::cout << std::setw(25) << "Time (s)" << std::setw(14) << "Power (MW)"
            << std::setw(14) << "Rho (pcm)" << std::setw(14) << "T_fuel (C)"
            << std::setw(14) << "T_mod (C)" << std::endl;

  reactor.print_states();

    double dt = 0.01;
    double t_end = 10000.0;
    int steps = static_cast<int>(t_end / dt);

    double insert_time = 5000.0;
    double withdraw_time = 7000.0;  // time when rod is removed
    double boron_inject_time = 8000.0;

    bool rod_inserted = false;
    bool rod_withdrawn = false;
    bool boron_injected = false;

    for (int i = 1; i <= steps; ++i) {
    double current_time = i * dt;

    if (!rod_inserted && current_time >= insert_time) {
        reactor.insertControlRod(1.0);
        rod_inserted = true;
        // std::cout << "--- Rod inserted at t = " << current_time << " s ---" << std::endl;
    }

    if (!rod_withdrawn && current_time >= withdraw_time) {
        reactor.setInitialReactivity(0);
        rod_withdrawn = true;
        // std::cout << "--- Rod withdrawn at t = " << current_time << " s ---" << std::endl;
    }

    if (!boron_injected and current_time>=boron_inject_time) {
      reactor.injectBoron(100);
      boron_injected = true;
    }

    reactor.timeStep(dt);

    if (i % 10000 == 0) {
        reactor.print_states();
    }
    }

  std::cout << "\n\nSimulation Complete\n";
  std::cout << "Final State:\n";
  std::cout << "  Power: " << reactor.getPower() << " MW\n";
  std::cout << "  Fuel Temperature: " << reactor.getFuelTemperature() << " C\n";
  std::cout << "  Moderator Temperature: " << reactor.getModeratorTemperature() << " C\n";
  std::cout << "  Reactivity: " << reactor.getReactivity() * 1e5 << " pcm\n";

  return 0;
}