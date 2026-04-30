#ifndef ASTARA_CONTROL_PID_CONTROLLER_HPP
#define ASTARA_CONTROL_PID_CONTROLLER_HPP

/**
 * @file   PIDController.hpp
 * @brief  Generic PID controller with output limits and clamping anti-windup.
 *
 * Implements the classical "ideal" PID form
 *
 * @f[
 *    u(t) \;=\; K_p\, e(t)
 *           \;+\; K_i \int_0^t e(\tau)\, d\tau
 *           \;+\; K_d\, \frac{de(t)}{dt}
 * @f]
 *
 * with discrete-time backward-difference derivative (avoids derivative kick
 * on a setpoint change) and clamping anti-windup on the integrator (the
 * integrator is frozen whenever the unsaturated output saturates and the
 * sign of the error matches the saturation direction).
 *
 * Output is clamped to `[u_min, u_max]`.  `setpoint`, `Kp`, `Ki`, `Kd`,
 * `u_min`, `u_max` may all be changed at run time.
 *
 * Threading: this class is **not** thread-safe.  Each PID instance is
 * intended to be owned by a single controller object.
 *
 * Usage:
 * @code
 *   PIDController pid({.Kp = 0.5, .Ki = 0.05, .Kd = 0.0,
 *                      .u_min = 0.0, .u_max = 1.0});
 *   pid.setSetpoint(target);
 *   const double u = pid.update(measurement, dt_s);
 * @endcode
 */

namespace astara::control {

/**
 * @brief Configuration for a `PIDController`.
 */
struct PIDConfig {
    double Kp = 0.0;     ///< proportional gain
    double Ki = 0.0;     ///< integral gain     (units: 1/s)
    double Kd = 0.0;     ///< derivative gain   (units: s)
    double u_min = 0.0;  ///< lower output limit
    double u_max = 0.0;  ///< upper output limit
};

/**
 * @brief Discrete-time PID controller with clamping anti-windup.
 */
class PIDController {
public:
    explicit PIDController(PIDConfig cfg) : cfg_(cfg) {}

    /// Read/write configuration access.
    PIDConfig&       config()       noexcept { return cfg_; }
    const PIDConfig& config() const noexcept { return cfg_; }

    void setSetpoint(double sp) noexcept { setpoint_ = sp; }
    double setpoint() const noexcept    { return setpoint_; }

    /// Reset integrator and last-error memory to zero.
    void reset() noexcept {
        integrator_   = 0.0;
        last_error_   = 0.0;
        has_last_     = false;
    }

    /**
     * @brief Compute one control step.
     * @param measurement  current process value
     * @param dt_s         time elapsed since previous call, s
     * @return             clamped control output in [u_min, u_max]
     */
    double update(double measurement, double dt_s) noexcept {
        const double error = setpoint_ - measurement;

        // Tentative integrator and derivative.
        const double integrator_unsaturated = integrator_ + error * dt_s;
        const double derivative = (has_last_ && dt_s > 0.0)
                ? (error - last_error_) / dt_s
                : 0.0;

        const double u_unsat = cfg_.Kp * error
                             + cfg_.Ki * integrator_unsaturated
                             + cfg_.Kd * derivative;

        // Clamp.
        double u = u_unsat;
        if (u > cfg_.u_max)      u = cfg_.u_max;
        else if (u < cfg_.u_min) u = cfg_.u_min;

        // Anti-windup: only commit the integrator update if we are not
        // saturating in the same direction as the error would push us.
        const bool saturated_high = (u_unsat > cfg_.u_max);
        const bool saturated_low  = (u_unsat < cfg_.u_min);
        if ((saturated_high && error > 0.0) || (saturated_low && error < 0.0)) {
            // Hold the integrator (do not commit `integrator_unsaturated`).
        } else {
            integrator_ = integrator_unsaturated;
        }
        last_error_ = error;
        has_last_   = true;

        return u;
    }

    /// Inspect internal state (mostly for tests).
    double integrator() const noexcept { return integrator_; }
    double lastError()  const noexcept { return last_error_; }

private:
    PIDConfig cfg_;
    double    setpoint_   = 0.0;
    double    integrator_ = 0.0;
    double    last_error_ = 0.0;
    bool      has_last_   = false;
};

}  // namespace astara::control

#endif  // ASTARA_CONTROL_PID_CONTROLLER_HPP
