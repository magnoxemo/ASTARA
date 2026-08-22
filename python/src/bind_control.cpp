#include "astara/control/PIDController.hpp"
#include "astara/control/ReactorController.hpp"
#include "astara/control/PressurizerController.hpp"
#include "astara/control/ThreeElementController.hpp"
#include "astara/control/ControlRodSystem.hpp"

#include <nanobind/nanobind.h>

namespace nb = nanobind;
using namespace astara::control;

void bind_control(nb::module_& m) {
    nb::class_<PIDConfig>(m, "PIDConfig")
        .def(nb::init<>())
        .def_rw("Kp",    &PIDConfig::Kp)
        .def_rw("Ki",    &PIDConfig::Ki)
        .def_rw("Kd",    &PIDConfig::Kd)
        .def_rw("u_min", &PIDConfig::u_min)
        .def_rw("u_max", &PIDConfig::u_max);

    nb::class_<PIDController>(m, "PIDController")
        .def(nb::init<PIDConfig>())
        .def("config", static_cast<PIDConfig& (PIDController::*)()>(&PIDController::config),
             nb::rv_policy::reference_internal)
        .def("set_setpoint", &PIDController::setSetpoint)
        .def("setpoint",     &PIDController::setpoint)
        .def("reset",        &PIDController::reset)
        .def("update",       &PIDController::update)
        .def("integrator",   &PIDController::integrator)
        .def("last_error",   &PIDController::lastError);

    auto rc = nb::class_<ReactorController>(m, "ReactorController");
    nb::class_<ReactorController::Config>(rc, "Config")
        .def(nb::init<>())
        .def_rw("pid",                          &ReactorController::Config::pid)
        .def_rw("T_avg_setpoint_no_load_K",     &ReactorController::Config::T_avg_setpoint_no_load_K)
        .def_rw("T_avg_setpoint_full_load_K",   &ReactorController::Config::T_avg_setpoint_full_load_K)
        .def_rw("max_reactivity_rate",          &ReactorController::Config::max_reactivity_rate)
        .def_rw("rho_min",                      &ReactorController::Config::rho_min)
        .def_rw("rho_max",                      &ReactorController::Config::rho_max);
    rc.def(nb::init<ReactorController::Config>())
      .def("set_turbine_load_fraction", &ReactorController::setTurbineLoadFraction)
      .def("turbine_load_fraction",     &ReactorController::turbineLoadFraction)
      .def("reset",  &ReactorController::reset, nb::arg("rho0") = 0.0)
      .def("update", &ReactorController::update)
      .def("external_reactivity", &ReactorController::externalReactivity)
      .def("config",
           static_cast<ReactorController::Config& (ReactorController::*)()>(&ReactorController::config),
           nb::rv_policy::reference_internal);

    auto pc = nb::class_<PressurizerController>(m, "PressurizerController");
    nb::class_<PressurizerController::Config>(pc, "Config")
        .def(nb::init<>())
        .def_rw("pressure_setpoint_Pa",   &PressurizerController::Config::pressure_setpoint_Pa)
        .def_rw("heater_P_gain_W_per_Pa", &PressurizerController::Config::heater_P_gain_W_per_Pa)
        .def_rw("heater_min_W",           &PressurizerController::Config::heater_min_W)
        .def_rw("heater_max_W",           &PressurizerController::Config::heater_max_W)
        .def_rw("heater_steady_state_W",  &PressurizerController::Config::heater_steady_state_W)
        .def_rw("spray_threshold_Pa",     &PressurizerController::Config::spray_threshold_Pa)
        .def_rw("spray_max_kg_s",         &PressurizerController::Config::spray_max_kg_s)
        .def_rw("spray_gain_kg_per_Pa",   &PressurizerController::Config::spray_gain_kg_per_Pa)
        .def_rw("spray_dead_band_Pa",     &PressurizerController::Config::spray_dead_band_Pa);
    pc.def(nb::init<PressurizerController::Config>())
      .def("heater_demand", &PressurizerController::heaterDemand)
      .def("spray_demand",  &PressurizerController::sprayDemand)
      .def("reset",         &PressurizerController::reset)
      .def("config",
           static_cast<PressurizerController::Config& (PressurizerController::*)()>(
               &PressurizerController::config),
           nb::rv_policy::reference_internal)
      .def("sprays_open", &PressurizerController::spraysOpen);

    nb::class_<ThreeElementController>(m, "ThreeElementController")
        .def(nb::init<double, PIDConfig>(), nb::arg("level_setpoint_m"), nb::arg("level_pi_cfg") = PIDConfig{})
        .def("set_level_setpoint", &ThreeElementController::setLevelSetpoint)
        .def("level_setpoint",     &ThreeElementController::levelSetpoint)
        .def("reset",  &ThreeElementController::reset)
        .def("update", &ThreeElementController::update)
        .def("config",
             static_cast<PIDConfig& (ThreeElementController::*)()>(&ThreeElementController::config),
             nb::rv_policy::reference_internal)
        .def("level_integrator", &ThreeElementController::levelIntegrator);

    nb::enum_<RodMode>(m, "RodMode")
        .value("Insert",   RodMode::Insert)
        .value("Withdraw", RodMode::Withdraw)
        .value("Hold",     RodMode::Hold);

    nb::class_<ControlRodSystemConfig>(m, "ControlRodSystemConfig")
        .def(nb::init<>())
        .def_rw("insertion_speed",  &ControlRodSystemConfig::insertion_speed)
        .def_rw("withdrawal_speed", &ControlRodSystemConfig::withdrawal_speed)
        .def_rw("step_reactivity",  &ControlRodSystemConfig::step_reactivity);

    nb::class_<ControlRodSystem>(m, "ControlRodSystem")
        .def(nb::init<ControlRodSystemConfig>())
        .def("config",
             static_cast<ControlRodSystemConfig& (ControlRodSystem::*)()>(&ControlRodSystem::config),
             nb::rv_policy::reference_internal)
        .def("reactivity_rate", &ControlRodSystem::reactivityRate);
}
