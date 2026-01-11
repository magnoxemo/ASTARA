#include "iostream"
#include "Reactor.h"
#include "Parser.h"

astara::Reactor::Reactor(unsigned int n_groups, std::vector<double> delayed_neutron_constants)
               : _number_of_neutron_groups(n_groups),_delayed_neutron_constants(delayed_neutron_constants){

    if (_number_of_neutron_groups != _delayed_neutron_constants.size())
        std::runtime_error("Number of neutron group must be equal to the number of delayed consts");
}
