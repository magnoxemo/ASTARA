/**
 * @file   PrimaryLoop.cpp
 * @brief  Wiring + sequential time-stepping of the integrated primary loop.
 */

#include "astara/primary/PrimaryLoop.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <cmath>
#include <stdexcept>
#include <utility>

namespace astara::primary {

PrimaryLoop::PrimaryLoop(std::unique_ptr<reactor::Reactor>       reactor,
                         std::unique_ptr<sg::AliSteamGenerator>  sg,
                         std::unique_ptr<pump::HomologousPump>   pump,
                         std::unique_ptr<pressurizer::Pressurizer> pressurizer)
        : reactor_(std::move(reactor)),
          sg_(std::move(sg)),
          pump_(std::move(pump)),
          pressurizer_(std::move(pressurizer)) {
    if (!reactor_)     throw std::invalid_argument("PrimaryLoop: reactor null");
    if (!sg_)          throw std::invalid_argument("PrimaryLoop: SG null");
    if (!pump_)        throw std::invalid_argument("PrimaryLoop: pump null");
    if (!pressurizer_) throw std::invalid_argument("PrimaryLoop: pressurizer null");
}

bool PrimaryLoop::isConsistent(double tol_K) const noexcept {
    // SG primary-in should track reactor hot-leg-out; SG primary-out should
    // be in the right ballpark for reactor cold-leg-in.
    const double T_hot      = reactor_->state().thermal.T_hot_leg;
    const double T_sg_in    = sg_->inputs().primary_inlet_temperature_K;
    const double T_sg_out   = sg_->primaryOutletTemperatureK();
    const double T_cold_set = reactor_->coldLegInletTemperature();
    return std::abs(T_hot - T_sg_in) < tol_K
        && std::abs(T_sg_out - T_cold_set) < tol_K;
}

void PrimaryLoop::timeStep(double dt) {
    if (!(dt > 0.0)) {
        throw std::invalid_argument("PrimaryLoop::timeStep: dt must be positive");
    }

    // -------------------------------------------------------------------------
    // 1. Update boundary conditions before the step (operator-splitting).
    // -------------------------------------------------------------------------

    // Pump output -> reactor + SG primary mass flow.  We only read the
    // current pump state; any update from a controller would have been
    // applied via setInputPowerW() prior to this call.
    const double mdot = pump_->massFlow_kg_s();

    // Reactor sees SG primary outlet as cold-leg-in (after the pump, but
    // we treat the pump as instantaneous on the thermal timescale).
    reactor_->setColdLegInletTemperature(sg_->primaryOutletTemperatureK());

    // SG sees reactor hot-leg-out as primary inlet.
    sg_->inputs().primary_inlet_temperature_K = reactor_->state().thermal.T_hot_leg;
    sg_->inputs().primary_mass_flow_kg_s      = mdot;

    // Pressurizer surge: a simple proxy is to drive surge mass flow
    // proportional to dP/dt of cold-leg pressure, which we don't model
    // here.  In the absence of explicit surge dynamics, set surge = 0
    // (the pressurizer remains an isolated heater + spray volume).
    // The user can override this between calls.
    pressurizer_->inputs().surge_enthalpy_J_kg = 0.0;  // h_f at current P (default)

    // -------------------------------------------------------------------------
    // 2. Run controllers if attached.
    // -------------------------------------------------------------------------
    if (reactor_controller_) {
        const double T_avg = 0.5 * (reactor_->state().thermal.T_hot_leg
                                  + reactor_->state().thermal.T_cold_leg);
        const double rho   = reactor_controller_->update(T_avg, dt);
        reactor_->reactivity().rho_external = rho;
    }
    if (pressurizer_controller_) {
        pressurizer_->inputs().heater_power_W
                = pressurizer_controller_->heaterDemand(pressurizer_->state().pressure_Pa);
        pressurizer_->inputs().spray_mass_flow_kg_s
                = pressurizer_controller_->sprayDemand(pressurizer_->state().pressure_Pa);
    }
    if (feedwater_controller_) {
        const double L_drum = sg_->state().L_dw;
        const double W_st   = sg_->steamMassFlow_kg_s();
        sg_->inputs().feedwater_mass_flow_kg_s
                = feedwater_controller_->update(L_drum, W_st, dt);
    }

    // -------------------------------------------------------------------------
    // 3. Step each component (in parallel where it's safe and helpful).
    // -------------------------------------------------------------------------

#ifdef _OPENMP
    #pragma omp parallel sections default(none) shared(dt)
    {
        #pragma omp section
        { reactor_->timeStep(dt); }
        #pragma omp section
        { sg_->timeStep(dt); }
        #pragma omp section
        { pump_->timeStep(dt); }
        #pragma omp section
        { pressurizer_->timeStep(dt); }
    }
#else
    reactor_->timeStep(dt);
    sg_->timeStep(dt);
    pump_->timeStep(dt);
    pressurizer_->timeStep(dt);
#endif
}

}  // namespace astara::primary
