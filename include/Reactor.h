#ifndef ASTARA_REACTOR_H
#define ASTARA_REACTOR_H

#include <memory>
#include <vector>

namespace astara {

class Function;

class Reactor {

public:
  Reactor(unsigned int n_groups,
          std::vector<double> neutron_group_const,
          std::vector<double> delayed_neutron_constants,
          double neutron_generation_time);
  ~Reactor();

  /* initial condition setters*/

  void setInitialFuelTemperature(double fuel_temperature) {
    _fuel_temperature = fuel_temperature;
  };
  void setInitialModeratorTemperature(double moderator_temperature) {
    _moderator_temperature = moderator_temperature;
  };

  /* Methods for setting functional expression */
  void setHeatTransferCoefficient(Function *heat_transfer_co_eff_function);
  void setFuelSpecificHeatFunction(Function *fuel_specific_function);
  void
  setFuelTemperatureCoEfficientFunction(Function *fuel_temp_feed_back_func);
  void setModeratorTemperatureCoEfficientFunction(
      Function *moderator_temp_feed_back_func);

protected:
  // transients
  /** differential eq for transient reactivity calculation */
  void dRhoDt();

  /*differential eq for transient power calculation */
  void dPowerDt();

  /*differential eq for transient fuel temperature  calculation*/
  void dFuelTempDt();

  /*differential eq for transient precursor calculation */
  void dCDt();

  // reactivity control via external medium
  void insertControlRod(double length);
  void injectBoron(double boron_concentration);

  /* calculate the over all states */
  void updateCurrentState();

  /* Method for broadcasting current state and data to other components*/
  void broadCastState();

private:
  /**
   * Most of the time heat transfer coefficient and fuel specific heat can be
   * time and temperature depended. So making it a function with arbitrary
   * number of inputs is a better idea
   */
  std::unique_ptr<Function> _convective_heat_transfer_co_efficient;
  std::unique_ptr<Function> _fuel_specific_heat;
  std::unique_ptr<Function> _fuel_temperature_co_efficient;
  std::unique_ptr<Function> _moderator_temperature_co_efficient;
  std::unique_ptr<Function> _boron_temperature_co_efficient;

  //neutron generation time
  const double _neutron_generation_time;

  //variables
  /* There is an initial reactivity which can be calculated from the startup
   * condition. Like how much boron in it, fuel moderator temperature, initial xenon presence
   * fuel inventory presence etc */
  double _reactivity = 0;

  // temperatures
  double _fuel_temperature;
  double _moderator_temperature;

  // delayed neutron fractions
  const unsigned int _number_of_neutron_groups;
  const std::vector<double> _neutron_group_const;
  const std::vector<double> _decay_constants;
  const double _total_decay_constant;
  const double _total_group_const;

};
} // namespace astara

#endif // ASTARA_REACTOR_H
