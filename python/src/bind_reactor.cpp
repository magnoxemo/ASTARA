#include "astara/reactor/PointKinetics.hpp"
#include "astara/reactor/ReactorThermal.hpp"
#include "astara/reactor/Reactor.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>

namespace nb = nanobind;
using namespace astara::reactor;

void bind_reactor(nb::module_& m) {
    nb::class_<DelayedGroupConstants>(m, "DelayedGroupConstants")
        .def(nb::init<>())
        .def_rw("decay_constants", &DelayedGroupConstants::lambda)
        .def_rw("beta",            &DelayedGroupConstants::beta)
        .def_rw("Lambda",          &DelayedGroupConstants::Lambda)
        .def("total_beta",         &DelayedGroupConstants::totalBeta)
        .def("validate",           &DelayedGroupConstants::validate)
        .def_static("u235_six_group",       &DelayedGroupConstants::u235SixGroup)
        .def_static("one_group_average_of", &DelayedGroupConstants::oneGroupAverageOf);

    nb::class_<PointKineticsState>(m, "PointKineticsState")
        .def(nb::init<>())
        .def(nb::init<std::size_t>(), nb::arg("num_groups"))
        .def_prop_rw("power",
            [](const PointKineticsState& s) { return s.power(); },
            [](PointKineticsState& s, double v) { s.power() = v; })
        .def_prop_rw("precursors",
            [](const PointKineticsState& s) { return s.precursors(); },
            [](PointKineticsState& s, std::vector<double> v) { s.precursors() = std::move(v); })
        .def("num_groups", &PointKineticsState::numGroups);

    m.def("steady_state_precursors", &steadyStatePrecursors, nb::arg("g"), nb::arg("n0"));
    m.def("point_kinetics_derivative", &pointKineticsDerivative,
          nb::arg("state"), nb::arg("rho"), nb::arg("g"));

    nb::class_<ReactivityModel>(m, "ReactivityModel")
        .def(nb::init<>())
        .def_rw("alpha_fuel_per_K",         &ReactivityModel::alpha_fuel_per_K)
        .def_rw("alpha_moderator_per_K",    &ReactivityModel::alpha_moderator_per_K)
        .def_rw("T_fuel_reference_K",       &ReactivityModel::T_fuel_reference_K)
        .def_rw("T_moderator_reference_K",  &ReactivityModel::T_moderator_reference_K)
        .def_rw("rho_external",             &ReactivityModel::rho_external)
        .def("evaluate", &ReactivityModel::evaluate);

    nb::class_<ReactorThermalParameters>(m, "ReactorThermalParameters")
        .def(nb::init<>())
        .def_rw("num_fuel_nodes",       &ReactorThermalParameters::num_fuel_nodes)
        .def_rw("num_moderator_nodes",  &ReactorThermalParameters::num_moderator_nodes)
        .def_rw("fuel_mass_total_kg",   &ReactorThermalParameters::fuel_mass_total_kg)
        .def_rw("fuel_cp_J_per_kgK",    &ReactorThermalParameters::fuel_cp_J_per_kgK)
        .def_rw("fission_power_in_fuel", &ReactorThermalParameters::fission_power_in_fuel)
        .def_rw("moderator_mass_total_kg", &ReactorThermalParameters::moderator_mass_total_kg)
        .def_rw("moderator_cp_J_per_kgK",  &ReactorThermalParameters::moderator_cp_J_per_kgK)
        .def_rw("mass_flow_rate_kg_s",  &ReactorThermalParameters::mass_flow_rate_kg_s)
        .def_rw("overall_h_W_per_m2K",  &ReactorThermalParameters::overall_h_W_per_m2K)
        .def_rw("heat_transfer_area_m2", &ReactorThermalParameters::heat_transfer_area_m2)
        .def_rw("lower_plenum_mass_kg", &ReactorThermalParameters::lower_plenum_mass_kg)
        .def_rw("upper_plenum_mass_kg", &ReactorThermalParameters::upper_plenum_mass_kg)
        .def_rw("hot_leg_mass_kg",      &ReactorThermalParameters::hot_leg_mass_kg)
        .def_rw("cold_leg_mass_kg",     &ReactorThermalParameters::cold_leg_mass_kg)
        .def("validate", &ReactorThermalParameters::validate);

    nb::class_<ReactorThermalState>(m, "ReactorThermalState")
        .def(nb::init<>())
        .def_rw("T_fuel",          &ReactorThermalState::T_fuel)
        .def_rw("T_moderator",     &ReactorThermalState::T_moderator)
        .def_rw("T_cold_leg",      &ReactorThermalState::T_cold_leg)
        .def_rw("T_lower_plenum",  &ReactorThermalState::T_lower_plenum)
        .def_rw("T_upper_plenum",  &ReactorThermalState::T_upper_plenum)
        .def_rw("T_hot_leg",       &ReactorThermalState::T_hot_leg)
        .def("size",  &ReactorThermalState::size)
        .def("average_fuel_temperature",      &ReactorThermalState::averageFuelTemperature)
        .def("average_moderator_temperature", &ReactorThermalState::averageModeratorTemperature);

    m.def("steady_state_thermal", &steadyStateThermal,
          nb::arg("p"), nb::arg("P_thermal"), nb::arg("T_inlet_K"));
    m.def("reactor_thermal_derivative", &reactorThermalDerivative,
          nb::arg("state"), nb::arg("p"), nb::arg("P_thermal_W"), nb::arg("T_inlet_K"));

    nb::class_<ReactorState>(m, "ReactorState")
        .def(nb::init<>())
        .def_rw("t_s",      &ReactorState::t_s)
        .def_rw("kinetics", &ReactorState::kinetics)
        .def_rw("thermal",  &ReactorState::thermal);

    nb::class_<Reactor>(m, "Reactor")
        .def(nb::init<DelayedGroupConstants, ReactorThermalParameters, double, ReactivityModel>(),
             nb::arg("groups"), nb::arg("thermal_parameters"),
             nb::arg("rated_thermal_power_W"), nb::arg("reactivity") = ReactivityModel{})
        .def("initialise_steady_state", &Reactor::initialiseSteadyState,
             nb::arg("n0"), nb::arg("T_inlet_K"))
        .def("set_cold_leg_inlet_temperature", &Reactor::setColdLegInletTemperature)
        .def("cold_leg_inlet_temperature",     &Reactor::coldLegInletTemperature)
        .def("reactivity",
             static_cast<ReactivityModel& (Reactor::*)()>(&Reactor::reactivity),
             nb::rv_policy::reference_internal)
        .def("time_step", &Reactor::timeStep)
        .def("state",         &Reactor::state)
        .def("groups",        &Reactor::groups)
        .def("thermal_params", &Reactor::thermalParams)
        .def("rated_power_W", &Reactor::ratedPowerW)
        .def("thermal_power_W", &Reactor::thermalPowerW)
        .def("hot_leg_outlet_temperature_K", &Reactor::hotLegOutletTemperatureK)
        .def("evaluate_derivative", &Reactor::evaluateDerivative);
}
