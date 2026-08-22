#include "astara/core/Port.hpp"
#include "astara/core/Integrator.hpp"

#include <nanobind/nanobind.h>

namespace nb = nanobind;
using namespace astara::core;

void bind_core(nb::module_& m) {
    nb::class_<FlowPort>(m, "FlowPort")
        .def(nb::init<>())
        .def_rw("mass_flow_kg_s", &FlowPort::mass_flow_kg_s)
        .def_rw("enthalpy_J_kg",  &FlowPort::enthalpy_J_kg)
        .def_rw("pressure_Pa",    &FlowPort::pressure_Pa)
        .def_rw("temperature_K",  &FlowPort::temperature_K);

    nb::enum_<IntegratorKind>(m, "IntegratorKind")
        .value("Euler", IntegratorKind::Euler)
        .value("RK4",   IntegratorKind::RK4);
}
