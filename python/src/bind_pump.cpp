#include "astara/pump/HomologousPump.hpp"

#include <nanobind/nanobind.h>

namespace nb = nanobind;
using namespace astara::pump;

void bind_pump(nb::module_& m) {
    nb::class_<PumpCurveCoefficients>(m, "PumpCurveCoefficients")
        .def(nb::init<>())
        .def_rw("A0", &PumpCurveCoefficients::A0)
        .def_rw("A1", &PumpCurveCoefficients::A1)
        .def_rw("A2", &PumpCurveCoefficients::A2);

    nb::class_<HomologousPumpParameters>(m, "HomologousPumpParameters")
        .def(nb::init<>())
        .def_rw("curve",                        &HomologousPumpParameters::curve)
        .def_rw("rated_speed_rev_s",             &HomologousPumpParameters::rated_speed_rev_s)
        .def_rw("rated_volumetric_flow_m3_s",    &HomologousPumpParameters::rated_volumetric_flow_m3_s)
        .def_rw("loop_resistance_K_s2_m5",       &HomologousPumpParameters::loop_resistance_K_s2_m5)
        .def_rw("effective_flow_area_m2",        &HomologousPumpParameters::effective_flow_area_m2)
        .def_rw("loop_length_m",                 &HomologousPumpParameters::loop_length_m)
        .def_rw("fluid_density_kg_m3",           &HomologousPumpParameters::fluid_density_kg_m3)
        .def_rw("moment_of_inertia_kg_m2",       &HomologousPumpParameters::moment_of_inertia_kg_m2)
        .def_rw("rated_input_power_W",           &HomologousPumpParameters::rated_input_power_W)
        .def_rw("gravity_m_s2",                  &HomologousPumpParameters::gravity_m_s2)
        .def("validate", &HomologousPumpParameters::validate);

    nb::class_<HomologousPumpState>(m, "HomologousPumpState")
        .def(nb::init<>())
        .def_rw("t_s",             &HomologousPumpState::t_s)
        .def_rw("volumetric_flow", &HomologousPumpState::volumetric_flow)
        .def_rw("speed_rev_s",     &HomologousPumpState::speed_rev_s);

    nb::class_<HomologousPump>(m, "HomologousPump")
        .def(nb::init<HomologousPumpParameters>())
        .def("initialise_at_rated", &HomologousPump::initialiseAtRated)
        .def("set_input_power_W",   &HomologousPump::setInputPowerW)
        .def("time_step",           &HomologousPump::timeStep)
        .def("state",               &HomologousPump::state)
        .def("params",              &HomologousPump::params)
        .def("developed_head_m",    &HomologousPump::developedHead_m)
        .def("mass_flow_kg_s",      &HomologousPump::massFlow_kg_s)
        .def("evaluate_derivative", &HomologousPump::evaluateDerivative)
        .def("head_at",             &HomologousPump::headAt);
}
