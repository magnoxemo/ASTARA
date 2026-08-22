#include "astara/props/WaterProperties.hpp"
#include "astara/props/IF97Water.hpp"
#include "astara/props/LinearizedWater.hpp"

#include <nanobind/nanobind.h>

namespace nb = nanobind;
using namespace astara::props;

void bind_props(nb::module_& m) {
    nb::exception<PropertyError>(m, "PropertyError", PyExc_RuntimeError);

    nb::class_<WaterProperties>(m, "WaterProperties")
        .def("density_TP",              &WaterProperties::density_TP)
        .def("enthalpy_TP",             &WaterProperties::enthalpy_TP)
        .def("cp_TP",                   &WaterProperties::cp_TP)
        .def("saturation_temperature",  &WaterProperties::saturationTemperature)
        .def("saturation_pressure",     &WaterProperties::saturationPressure)
        .def("sat_liquid_density_P",    &WaterProperties::satLiquidDensity_P)
        .def("sat_vapour_density_P",    &WaterProperties::satVapourDensity_P)
        .def("sat_liquid_enthalpy_P",   &WaterProperties::satLiquidEnthalpy_P)
        .def("sat_vapour_enthalpy_P",   &WaterProperties::satVapourEnthalpy_P)
        .def("latent_heat_P",           &WaterProperties::latentHeat_P)
        .def("sat_liquid_specific_volume_P", &WaterProperties::satLiquidSpecificVolume_P)
        .def("sat_vapour_specific_volume_P", &WaterProperties::satVapourSpecificVolume_P)
        .def("two_phase_enthalpy_PQ",   &WaterProperties::twoPhaseEnthalpy_PQ);

    nb::class_<IF97Water, WaterProperties>(m, "IF97Water")
        .def(nb::init<bool>(), nb::arg("strict") = true);

    auto lw = nb::class_<LinearizedWater, WaterProperties>(m, "LinearizedWater");

    nb::class_<LinearizedWater::Coefficients>(lw, "Coefficients")
        .def(nb::init<>())
        .def_rw("X_Tsat", &LinearizedWater::Coefficients::X_Tsat)
        .def_rw("K_Tsat", &LinearizedWater::Coefficients::K_Tsat)
        .def_rw("X_hf",   &LinearizedWater::Coefficients::X_hf)
        .def_rw("K_hf",   &LinearizedWater::Coefficients::K_hf)
        .def_rw("X_hg",   &LinearizedWater::Coefficients::X_hg)
        .def_rw("K_hg",   &LinearizedWater::Coefficients::K_hg)
        .def_rw("X_vf",   &LinearizedWater::Coefficients::X_vf)
        .def_rw("K_vf",   &LinearizedWater::Coefficients::K_vf)
        .def_rw("X_vg",   &LinearizedWater::Coefficients::X_vg)
        .def_rw("K_vg",   &LinearizedWater::Coefficients::K_vg)
        .def_rw("X_rhof", &LinearizedWater::Coefficients::X_rhof)
        .def_rw("K_rhof", &LinearizedWater::Coefficients::K_rhof)
        .def_rw("X_rhog", &LinearizedWater::Coefficients::X_rhog)
        .def_rw("K_rhog", &LinearizedWater::Coefficients::K_rhog);

    lw.def(nb::init<>())
      .def(nb::init<const LinearizedWater::Coefficients&>())
      .def_static("fit_around", &LinearizedWater::fitAround,
                  nb::arg("ref"), nb::arg("P0_Pa"), nb::arg("dP_Pa") = 5.0e5)
      .def("coefficients", &LinearizedWater::coefficients);
}
