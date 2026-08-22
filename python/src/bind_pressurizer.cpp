#include "astara/pressurizer/Pressurizer.hpp"
#include "astara/props/WaterProperties.hpp"

#include <nanobind/nanobind.h>

namespace nb = nanobind;
using namespace astara::pressurizer;
using astara::props::WaterProperties;

void bind_pressurizer(nb::module_& m) {
    nb::class_<PressurizerParameters>(m, "PressurizerParameters")
        .def(nb::init<>())
        .def_rw("cross_section_area_m2", &PressurizerParameters::cross_section_area_m2)
        .def_rw("total_height_m",        &PressurizerParameters::total_height_m)
        .def("validate", &PressurizerParameters::validate);

    nb::class_<PressurizerState>(m, "PressurizerState")
        .def(nb::init<>())
        .def_rw("t_s",           &PressurizerState::t_s)
        .def_rw("pressure_Pa",   &PressurizerState::pressure_Pa)
        .def_rw("water_level_m", &PressurizerState::water_level_m);

    nb::class_<PressurizerInputs>(m, "PressurizerInputs")
        .def(nb::init<>())
        .def_rw("surge_mass_flow_kg_s", &PressurizerInputs::surge_mass_flow_kg_s)
        .def_rw("spray_mass_flow_kg_s", &PressurizerInputs::spray_mass_flow_kg_s)
        .def_rw("heater_power_W",       &PressurizerInputs::heater_power_W)
        .def_rw("spray_enthalpy_J_kg",  &PressurizerInputs::spray_enthalpy_J_kg)
        .def_rw("surge_enthalpy_J_kg",  &PressurizerInputs::surge_enthalpy_J_kg);

    nb::class_<Pressurizer>(m, "Pressurizer")
        .def(nb::init<PressurizerParameters, const WaterProperties*>(),
             nb::arg("params"), nb::arg("props"), nb::keep_alive<1, 3>())
        .def("initialise_steady_state", &Pressurizer::initialiseSteadyState,
             nb::arg("P_Pa"), nb::arg("L_w_m"))
        .def("inputs",
             static_cast<PressurizerInputs& (Pressurizer::*)()>(&Pressurizer::inputs),
             nb::rv_policy::reference_internal)
        .def("time_step", &Pressurizer::timeStep)
        .def("state",  &Pressurizer::state)
        .def("params", &Pressurizer::params)
        .def("water_mass_kg", &Pressurizer::waterMass_kg)
        .def("steam_mass_kg", &Pressurizer::steamMass_kg)
        .def("evaluate_derivative", &Pressurizer::evaluateDerivative);
}
