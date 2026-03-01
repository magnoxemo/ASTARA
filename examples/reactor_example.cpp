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

  // Total delayed neutron fraction (should sum to ~0.0065)
  double total_beta = 0.0;
  for (double beta : beta_values) {
    total_beta += beta;
  }
  std::cout << "Total delayed neutron fraction: " << total_beta << std::endl;

  // Prompt neutron generation time (s)
  double Lambda = 5.0e-5; // 50 microseconds

  Reactor reactor(n_groups, decay_constants, beta_values, Lambda);

  // ===== Set Initial Conditions =====
  reactor.setInitialPower(100.0);                // 100 MW
  reactor.setInitialReactivity(0.0);             // Critical
  reactor.setInitialFuelTemperature(565.0);      // Celsius
  reactor.setInitialModeratorTemperature(300.0); // Celsius

  for (unsigned int i = 0; i < n_groups; ++i) {
    double C_i =
        (beta_values[i] * reactor.getPower()) / (decay_constants[i] * Lambda);
    reactor.setInitialPrecursorConcentration(i, C_i);
    std::cout << "Precursor " << i << " initial concentration: " << C_i
              << std::endl;
  }
  reactor.setFuelMass(2000.0);                 // kg
  reactor.setModeratorMass(50000.0);           // kg
  reactor.setFuelThermalCapacity(500.0);       // J/(kg·K)
  reactor.setModeratorThermalCapacity(4186.0); // J/(kg·K)
  reactor.setHeatTransferArea(500.0);          // m^2
  reactor.setFissionEnergy(3.2e-11);           // J/fission

  reactor.setControlRodEffectiveness(-0.0001); // -10 pcm/cm
  reactor.setBoronEffectiveness(-1.0e-5);      // -1 pcm/ppm

  Function *h_func = new Function("5000 + 10 * (x - 300)", {"x"}
                                  // x represents fuel temperature
  );
  reactor.setHeatTransferCoefficient(h_func);

  // Fuel specific heat: c_f = 500 + 0.2*(T - 300)
  Function *c_f_func = new Function("500 + 0.2 * (x - 300)", {"x"});
  reactor.setFuelSpecificHeatFunction(c_f_func);

  Function *alpha_f_func = new Function("-0.01 * (x - 56)", {"x"});
  reactor.setFuelTemperatureCoEfficientFunction(alpha_f_func);
  Function *alpha_m_func = new Function("-0.000005 * (x - 300)", {"x"});
  reactor.setModeratorTemperatureCoEfficientFunction(alpha_m_func);

  std::cout << std::fixed << std::setprecision(6);
  std::cout << "\nReactor Transient Analysis\n";
  std::cout << "===========================\n\n";
  std::cout << std::setw(10) << "Time (s)" << std::setw(14) << "Power (MW)"
            << std::setw(14) << "Rho (pcm)" << std::setw(14) << "T_fuel (C)"
            << std::setw(14) << "T_mod (C)" << std::endl;
  std::cout << std::string(70, '-') << std::endl;

  reactor.print_states();

  double dt = 0.0001;
  double t_end = 100.0;
  int steps = static_cast<int>(t_end / dt);
  double step_time = 1.0;
  double step_magnitude = 1;
  bool step_applied = false;

  for (int i = 1; i <= steps; ++i) {
    double current_time = i * dt;
    reactor.timeStep(dt);

    if (i % 100 == 0) {
      reactor.print_states();
    }
  }

  std::cout << "\n\nSimulation Complete\n";
  std::cout << "Final State:\n";
  std::cout << "  Power: " << reactor.getPower() << " MW\n";
  std::cout << "  Fuel Temperature: " << reactor.getFuelTemperature() << " C\n";
  std::cout << "  Moderator Temperature: " << reactor.getModeratorTemperature()
            << " C\n";
  std::cout << "  Reactivity: " << reactor.getReactivity() * 1e5 << " pcm\n";

  return 0;
}