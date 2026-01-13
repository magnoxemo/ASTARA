#include <iostream>
#include <numeric>

#include "Function.h"
#include "Reactor.h"

astara::Reactor::Reactor(unsigned int n_groups,
                         std::vector<double> neutron_group_const,
                         std::vector<double> delayed_neutron_constants,
                         double neutron_generation_time)
    : _number_of_neutron_groups(n_groups),
      _neutron_group_const(neutron_group_const),
      _decay_constants(delayed_neutron_constants),
      _total_decay_constant(std::accumulate(_decay_constants.begin(),
                                            _decay_constants.end(), 0.0)),
      _total_group_const(std::accumulate(_neutron_group_const.begin(),
                                         _neutron_group_const.end(), 0)),
      _neutron_generation_time(neutron_generation_time) {

  if (_number_of_neutron_groups != _decay_constants.size() and
      _number_of_neutron_groups != _neutron_group_const.size())
    throw std::runtime_error("Number of neutron group must be equal to the "
                             "number of delayed and group consts");
}

astara::Reactor::Reactor::~Reactor() = default;

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

void astara::Reactor::dRhoDt() {}