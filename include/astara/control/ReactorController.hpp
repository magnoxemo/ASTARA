#ifndef ASTARA_CONTROL_REACTOR_CONTROLLER_HPP
#define ASTARA_CONTROL_REACTOR_CONTROLLER_HPP

/**
 * @file   ReactorController.hpp
 * @brief  Reactor "turbine-following" controller.
 *
 * In a Westinghouse PWR the reactor is normally operated in turbine-following
 * mode: the turbine demand is an upstream variable, and the reactor adjusts
 * its power (via control-rod motion -> external reactivity) to match.  The
 * average primary temperature `T_avg = (T_hot + T_cold)/2` is the controlled
 * variable, with a setpoint that varies linearly with turbine load.
 *
 * This class implements that loop:
 *
 *   - input:  measured T_avg (K), turbine load fraction `0..1`
 *   - output: demanded external reactivity (absolute, not pcm) bounded by
 *             user-specified rod-bank speed limits
 *
 * The reactivity rate is rate-limited to mimic finite control-rod speed
 * (typical 72 steps/min ~ 0.6 mm/s ~ 1.4 pcm/s of insertion at full power
 * worth).  A PID controller computes the desired reactivity *rate*, which
 * is integrated to produce the absolute external reactivity actually
 * applied to the reactor.
 *
 * @cite Naghedolfeizi (1990), Chapter 4, Fig. 4.1 (reactor control system).
 */

#include "astara/control/PIDController.hpp"

#include <algorithm>

namespace astara::control {

/**
 * @brief Turbine-following reactor controller (T_avg control).
 */
class ReactorController {
public:
    struct Config {
        /// PID gains acting on (T_avg_setpoint - T_avg_measured); output in
        /// units of d(rho)/dt (1/s).
        PIDConfig pid;
        /// T_avg setpoint at zero load, K.
        double T_avg_setpoint_no_load_K = 562.0;
        /// T_avg setpoint at full load, K.
        double T_avg_setpoint_full_load_K = 583.0;
        /// Maximum rate of change of external reactivity (1/s).
        double max_reactivity_rate = 1.4e-5;   // ~ 1.4 pcm/s
        /// Hard limits on integrated external reactivity.
        double rho_min = -5.0e-3;              // -500 pcm
        double rho_max =  5.0e-3;              // +500 pcm
    };

    explicit ReactorController(Config cfg)
            : cfg_(cfg), pid_(cfg.pid) {}

    /// Update the turbine load fraction (0 = no load, 1 = full load).
    void setTurbineLoadFraction(double frac) noexcept {
        load_fraction_ = std::clamp(frac, 0.0, 1.0);
    }
    double turbineLoadFraction() const noexcept { return load_fraction_; }

    /// Reset accumulated reactivity and integrator.
    void reset(double rho0 = 0.0) noexcept {
        rho_external_ = rho0;
        pid_.reset();
    }

    /**
     * @brief Compute the external reactivity to apply this step.
     * @param T_avg_measured_K  primary average temperature (T_hot + T_cold)/2
     * @param dt_s              time since last update, s
     * @return                  external reactivity, dimensionless
     */
    double update(double T_avg_measured_K, double dt_s) noexcept {
        const double T_sp = cfg_.T_avg_setpoint_no_load_K
                          + load_fraction_ * (cfg_.T_avg_setpoint_full_load_K
                                              - cfg_.T_avg_setpoint_no_load_K);
        pid_.setSetpoint(T_sp);
        // PID output is the demanded reactivity *rate* (positive when T_avg
        // is below setpoint, requiring more power -> positive rho rate).
        const double rate_demanded = pid_.update(T_avg_measured_K, dt_s);
        const double rate = std::clamp(rate_demanded,
                                       -cfg_.max_reactivity_rate,
                                        cfg_.max_reactivity_rate);
        rho_external_ = std::clamp(rho_external_ + rate * dt_s,
                                   cfg_.rho_min, cfg_.rho_max);
        return rho_external_;
    }

    double externalReactivity() const noexcept { return rho_external_; }
    Config&       config()       noexcept { return cfg_; }
    const Config& config() const noexcept { return cfg_; }

private:
    Config         cfg_;
    PIDController  pid_;
    double         load_fraction_ = 1.0;
    double         rho_external_  = 0.0;
};

}  // namespace astara::control

#endif  // ASTARA_CONTROL_REACTOR_CONTROLLER_HPP
