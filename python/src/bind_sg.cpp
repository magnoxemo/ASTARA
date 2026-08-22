#include "astara/sg/AliSteamGenerator.hpp"
#include "astara/sg/HelicalCoilSteamGenerator.hpp"
#include "astara/props/WaterProperties.hpp"

#include <nanobind/nanobind.h>

namespace nb = nanobind;
using namespace astara::sg;
using astara::props::WaterProperties;

void bind_sg(nb::module_& m) {
    // ---------- AliSteamGenerator ----------
    nb::class_<AliSteamGeneratorParameters>(m, "AliSteamGeneratorParameters")
        .def(nb::init<>())
        .def_rw("primary_mass_per_node_kg",     &AliSteamGeneratorParameters::primary_mass_per_node_kg)
        .def_rw("primary_mass_inlet_plenum_kg", &AliSteamGeneratorParameters::primary_mass_inlet_plenum_kg)
        .def_rw("primary_cp_J_per_kgK",         &AliSteamGeneratorParameters::primary_cp_J_per_kgK)
        .def_rw("metal_mass_per_node_kg",       &AliSteamGeneratorParameters::metal_mass_per_node_kg)
        .def_rw("metal_cp_J_per_kgK",           &AliSteamGeneratorParameters::metal_cp_J_per_kgK)
        .def_rw("area_pm_per_node_m2",          &AliSteamGeneratorParameters::area_pm_per_node_m2)
        .def_rw("overall_h_pm_W_per_m2K",       &AliSteamGeneratorParameters::overall_h_pm_W_per_m2K)
        .def_rw("area_ms_per_node_m2",          &AliSteamGeneratorParameters::area_ms_per_node_m2)
        .def_rw("overall_h_ms_W_per_m2K",       &AliSteamGeneratorParameters::overall_h_ms_W_per_m2K)
        .def_rw("sec_flow_area_m2",             &AliSteamGeneratorParameters::sec_flow_area_m2)
        .def_rw("drum_water_area_m2",           &AliSteamGeneratorParameters::drum_water_area_m2)
        .def_rw("tube_bundle_height_m",         &AliSteamGeneratorParameters::tube_bundle_height_m)
        .def_rw("downcomer_length_m",           &AliSteamGeneratorParameters::downcomer_length_m)
        .def_rw("downcomer_volume_m3",          &AliSteamGeneratorParameters::downcomer_volume_m3)
        .def_rw("recirc_pressure_drop_coeff",   &AliSteamGeneratorParameters::recirc_pressure_drop_coeff)
        .def_rw("steam_valve_coefficient",      &AliSteamGeneratorParameters::steam_valve_coefficient)
        .def("validate", &AliSteamGeneratorParameters::validate)
        .def_static("westinghouse_model_D5", &AliSteamGeneratorParameters::westinghouseModelD5);

    nb::class_<AliSteamGeneratorState>(m, "AliSteamGeneratorState")
        .def(nb::init<>())
        .def_rw("t_s",   &AliSteamGeneratorState::t_s)
        .def_rw("T_pi",  &AliSteamGeneratorState::T_pi)
        .def_rw("T_p1",  &AliSteamGeneratorState::T_p1)
        .def_rw("T_p2",  &AliSteamGeneratorState::T_p2)
        .def_rw("T_p3",  &AliSteamGeneratorState::T_p3)
        .def_rw("T_p4",  &AliSteamGeneratorState::T_p4)
        .def_rw("T_m1",  &AliSteamGeneratorState::T_m1)
        .def_rw("T_m2",  &AliSteamGeneratorState::T_m2)
        .def_rw("T_m3",  &AliSteamGeneratorState::T_m3)
        .def_rw("T_m4",  &AliSteamGeneratorState::T_m4)
        .def_rw("T_sub", &AliSteamGeneratorState::T_sub)
        .def_rw("L_s1",  &AliSteamGeneratorState::L_s1)
        .def_rw("h_b",   &AliSteamGeneratorState::h_b)
        .def_rw("x_e",   &AliSteamGeneratorState::x_e)
        .def_rw("P",     &AliSteamGeneratorState::P)
        .def_rw("L_dw",  &AliSteamGeneratorState::L_dw)
        .def_rw("T_dw",  &AliSteamGeneratorState::T_dw)
        .def_rw("rho_r", &AliSteamGeneratorState::rho_r)
        .def_rw("T_dc",  &AliSteamGeneratorState::T_dc);

    nb::class_<AliSteamGeneratorInputs>(m, "AliSteamGeneratorInputs")
        .def(nb::init<>())
        .def_rw("primary_inlet_temperature_K", &AliSteamGeneratorInputs::primary_inlet_temperature_K)
        .def_rw("primary_mass_flow_kg_s",      &AliSteamGeneratorInputs::primary_mass_flow_kg_s)
        .def_rw("feedwater_mass_flow_kg_s",    &AliSteamGeneratorInputs::feedwater_mass_flow_kg_s)
        .def_rw("feedwater_enthalpy_J_kg",     &AliSteamGeneratorInputs::feedwater_enthalpy_J_kg)
        .def_rw("steam_line_pressure_Pa",      &AliSteamGeneratorInputs::steam_line_pressure_Pa);

    nb::class_<AliSteamGenerator>(m, "AliSteamGenerator")
        .def(nb::init<AliSteamGeneratorParameters, const WaterProperties*>(),
             nb::arg("params"), nb::arg("props"), nb::keep_alive<1, 3>())
        .def("initialise_steady_state", &AliSteamGenerator::initialiseSteadyState,
             nb::arg("T_pi_in_K"), nb::arg("primary_mass_flow_kg_s"),
             nb::arg("sg_pressure_Pa"), nb::arg("drum_level_m"))
        .def("inputs",
             static_cast<AliSteamGeneratorInputs& (AliSteamGenerator::*)()>(&AliSteamGenerator::inputs),
             nb::rv_policy::reference_internal)
        .def("time_step", &AliSteamGenerator::timeStep)
        .def("state",  &AliSteamGenerator::state)
        .def("params", &AliSteamGenerator::params)
        .def("primary_outlet_temperature_K", &AliSteamGenerator::primaryOutletTemperatureK)
        .def("steam_mass_flow_kg_s",         &AliSteamGenerator::steamMassFlow_kg_s)
        .def("primary_heat_load_W",          &AliSteamGenerator::primaryHeatLoad_W)
        .def("evaluate_derivative", &AliSteamGenerator::evaluateDerivative);

    // ---------- HelicalCoilSteamGenerator ----------
    nb::class_<HelicalCoilSteamGeneratorParameters>(m, "HelicalCoilSteamGeneratorParameters")
        .def(nb::init<>())
        .def_rw("tube_total_length_m",         &HelicalCoilSteamGeneratorParameters::tube_total_length_m)
        .def_rw("tube_inner_diameter_m",       &HelicalCoilSteamGeneratorParameters::tube_inner_diameter_m)
        .def_rw("tube_outer_diameter_m",       &HelicalCoilSteamGeneratorParameters::tube_outer_diameter_m)
        .def_rw("tube_count",                  &HelicalCoilSteamGeneratorParameters::tube_count)
        .def_rw("secondary_flow_area_m2",      &HelicalCoilSteamGeneratorParameters::secondary_flow_area_m2)
        .def_rw("primary_flow_area_m2",        &HelicalCoilSteamGeneratorParameters::primary_flow_area_m2)
        .def_rw("tube_metal_area_m2",          &HelicalCoilSteamGeneratorParameters::tube_metal_area_m2)
        .def_rw("tube_metal_density_kg_m3",    &HelicalCoilSteamGeneratorParameters::tube_metal_density_kg_m3)
        .def_rw("tube_metal_cp_J_per_kgK",     &HelicalCoilSteamGeneratorParameters::tube_metal_cp_J_per_kgK)
        .def_rw("primary_density_kg_m3",       &HelicalCoilSteamGeneratorParameters::primary_density_kg_m3)
        .def_rw("primary_cp_J_per_kgK",        &HelicalCoilSteamGeneratorParameters::primary_cp_J_per_kgK)
        .def_rw("alpha_o_subcooled_W_per_m2K", &HelicalCoilSteamGeneratorParameters::alpha_o_subcooled_W_per_m2K)
        .def_rw("alpha_o_twophase_W_per_m2K",  &HelicalCoilSteamGeneratorParameters::alpha_o_twophase_W_per_m2K)
        .def_rw("alpha_o_superheat_W_per_m2K", &HelicalCoilSteamGeneratorParameters::alpha_o_superheat_W_per_m2K)
        .def_rw("alpha_i_subcooled_W_per_m2K", &HelicalCoilSteamGeneratorParameters::alpha_i_subcooled_W_per_m2K)
        .def_rw("alpha_i_twophase_W_per_m2K",  &HelicalCoilSteamGeneratorParameters::alpha_i_twophase_W_per_m2K)
        .def_rw("alpha_i_superheat_W_per_m2K", &HelicalCoilSteamGeneratorParameters::alpha_i_superheat_W_per_m2K)
        .def("validate", &HelicalCoilSteamGeneratorParameters::validate)
        .def_static("nuscale_smr_two_sg", &HelicalCoilSteamGeneratorParameters::nuscaleSMRTwoSG);

    nb::class_<HelicalCoilSteamGeneratorState>(m, "HelicalCoilSteamGeneratorState")
        .def(nb::init<>())
        .def_rw("t_s",  &HelicalCoilSteamGeneratorState::t_s)
        .def_rw("L_1",  &HelicalCoilSteamGeneratorState::L_1)
        .def_rw("L_2",  &HelicalCoilSteamGeneratorState::L_2)
        .def_rw("p_S",  &HelicalCoilSteamGeneratorState::p_S)
        .def_rw("h_o",  &HelicalCoilSteamGeneratorState::h_o)
        .def_rw("T_M1", &HelicalCoilSteamGeneratorState::T_M1)
        .def_rw("T_M2", &HelicalCoilSteamGeneratorState::T_M2)
        .def_rw("T_M3", &HelicalCoilSteamGeneratorState::T_M3)
        .def_rw("T_P1", &HelicalCoilSteamGeneratorState::T_P1)
        .def_rw("T_P2", &HelicalCoilSteamGeneratorState::T_P2)
        .def_rw("T_P3", &HelicalCoilSteamGeneratorState::T_P3);

    nb::class_<HelicalCoilSteamGeneratorInputs>(m, "HelicalCoilSteamGeneratorInputs")
        .def(nb::init<>())
        .def_rw("primary_inlet_temperature_K", &HelicalCoilSteamGeneratorInputs::primary_inlet_temperature_K)
        .def_rw("primary_mass_flow_kg_s",      &HelicalCoilSteamGeneratorInputs::primary_mass_flow_kg_s)
        .def_rw("feedwater_mass_flow_kg_s",    &HelicalCoilSteamGeneratorInputs::feedwater_mass_flow_kg_s)
        .def_rw("feedwater_enthalpy_J_kg",     &HelicalCoilSteamGeneratorInputs::feedwater_enthalpy_J_kg)
        .def_rw("steam_outlet_mass_flow_kg_s", &HelicalCoilSteamGeneratorInputs::steam_outlet_mass_flow_kg_s);

    nb::class_<HelicalCoilSteamGenerator>(m, "HelicalCoilSteamGenerator")
        .def(nb::init<HelicalCoilSteamGeneratorParameters, const WaterProperties*>(),
             nb::arg("params"), nb::arg("props"), nb::keep_alive<1, 3>())
        .def("initialise_steady_state", &HelicalCoilSteamGenerator::initialiseSteadyState,
             nb::arg("T_pi_in_K"), nb::arg("primary_mass_flow_kg_s"),
             nb::arg("feedwater_flow_kg_s"), nb::arg("feedwater_temperature_K"))
        .def("inputs",
             static_cast<HelicalCoilSteamGeneratorInputs& (HelicalCoilSteamGenerator::*)()>(
                 &HelicalCoilSteamGenerator::inputs),
             nb::rv_policy::reference_internal)
        .def("time_step", &HelicalCoilSteamGenerator::timeStep)
        .def("state",  &HelicalCoilSteamGenerator::state)
        .def("params", &HelicalCoilSteamGenerator::params)
        .def("primary_outlet_temperature_K", &HelicalCoilSteamGenerator::primaryOutletTemperatureK)
        .def("steam_outlet_temperature_K",   &HelicalCoilSteamGenerator::steamOutletTemperatureK)
        .def("primary_heat_load_W",          &HelicalCoilSteamGenerator::primaryHeatLoad_W)
        .def("superheated_region_length_m",  &HelicalCoilSteamGenerator::superheatedRegionLength_m)
        .def("evaluate_derivative", &HelicalCoilSteamGenerator::evaluateDerivative);
}
