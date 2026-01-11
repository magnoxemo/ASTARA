#include <iostream>
#include <numeric>

#include "Parser.h"
#include "Reactor.h"

astara::Reactor::Reactor(unsigned int n_groups,
                         std::vector<double> delayed_neutron_constants)
    : _number_of_neutron_groups(n_groups),
      _delayed_neutron_constants(delayed_neutron_constants),
      _sum_delayed_neutron_constants(
          std::accumulate(_delayed_neutron_constants.begin(),
                          _delayed_neutron_constants.end(), 0.0)) {

  if (_number_of_neutron_groups != _delayed_neutron_constants.size())
    throw std::runtime_error("Number of neutron group must be equal to the "
                             "number of delayed consts");
}

void astara::Reactor::setHeatTransferCoefficient(
    Function *heat_transfer_co_eff_function) {
  if (heat_transfer_co_eff_function != nullptr)
    _convective_heat_transfer_co_efficient.reset(heat_transfer_co_eff_function);
  else
    throw std::runtime_error("heat_transfer_co_eff_function is null_ptr");
}

void astara::Reactor::setFuelSpecificHeatFunction(
    Function *fuel_specific_function) {
  if (fuel_specific_function != nullptr)
    _fuel_specific_heat.reset(fuel_specific_function);
  else
    throw std::runtime_error("fuel_specific_function is null_ptr");
}

void astara::Reactor::setFuelTemperatureCoEfficientFunction(
    Function *fuel_temp_feed_back_func) {
  if (fuel_temp_feed_back_func != nullptr)
    _fuel_temperature_co_efficient.reset(fuel_temp_feed_back_func);
  else
    throw std::runtime_error("fuel_temp_feed_back_func is null_ptr");
}

void astara::Reactor::setModeratorTemperatureCoEfficientFunction(
    Function *moderator_temp_feed_back_func) {
  if (moderator_temp_feed_back_func != nullptr)
    _moderator_temperature_co_efficient.reset(moderator_temp_feed_back_func);
  else
    throw std::runtime_error("moderator_temp_feed_back_func is null_ptr");
}

astara::Reactor::Reactor::~Reactor() = default;