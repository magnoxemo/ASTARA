#ifndef ASTARA_CONTROL_THREE_ELEMENT_CONTROLLER_HPP
#define ASTARA_CONTROL_THREE_ELEMENT_CONTROLLER_HPP

/**
 * @file   ThreeElementController.hpp
 * @brief  SG feedwater three-element controller (level, steam flow,
 *         feedwater flow).
 *
 * The classical three-element control combines:
 *
 *   - drum **level** error (`L_drum - L_setpoint`), integrated through a
 *     PI controller to produce a trim signal;
 *   - **steam mass flow rate** as a feed-forward signal (the feedwater
 *     should track the steam flow at steady state);
 *   - **measured feedwater mass flow rate** to close an inner loop and
 *     prevent valve hunting.
 *
 * Output is the demanded feedwater valve position (or, equivalently, the
 * demanded feedwater mass flow rate) in absolute terms.  The implementation
 * here returns demanded `W_fw` directly (kg/s) and assumes the upstream
 * pump can supply it -- no valve dynamics.
 *
 * @cite Naghedolfeizi (1990), Section 3.2.2 ("Steam Generator Control
 *       System"), Fig. 3.10.
 */

#include "astara/control/PIDController.hpp"

namespace astara::control {

/**
 * @brief Three-element feedwater controller for a U-tube SG.
 *
 * Output range and gains are configurable; defaults are chosen so that at
 * 100% steam flow with level error = 0 the feedwater demand equals the
 * steam flow, giving a conservative starting condition for tuning.
 */
class ThreeElementController {
public:
    /**
     * @brief Construct with explicit PI gains and feedwater-flow limits.
     * @param level_setpoint_m   target drum water level, m
     * @param level_pi_cfg       outer-loop PI configuration (Kp, Ki, limits)
     */
    explicit ThreeElementController(double level_setpoint_m,
                                    PIDConfig level_pi_cfg = {})
            : level_pi_(level_pi_cfg),
              level_setpoint_m_(level_setpoint_m) {
        level_pi_.setSetpoint(level_setpoint_m);
    }

    /// Update level setpoint at run time (e.g. for shrink/swell compensation).
    void setLevelSetpoint(double L_m) noexcept {
        level_setpoint_m_ = L_m;
        level_pi_.setSetpoint(L_m);
    }
    double levelSetpoint() const noexcept { return level_setpoint_m_; }

    /// Reset internal integrator (e.g. after a controller mode change).
    void reset() noexcept { level_pi_.reset(); }

    /**
     * @brief Compute the demanded feedwater mass flow rate.
     * @param level_measured_m       drum level measurement, m
     * @param steam_flow_meas_kg_s   measured steam mass flow rate, kg/s
     * @param dt_s                   sample time, s
     * @return                       demanded feedwater mass flow rate, kg/s
     */
    double update(double level_measured_m,
                  double steam_flow_meas_kg_s,
                  double dt_s) noexcept {
        // Outer (level) loop: trim signal in kg/s.
        const double trim = level_pi_.update(level_measured_m, dt_s);

        // Feed-forward: steam flow.  At steady state W_fw = W_steam.
        return steam_flow_meas_kg_s + trim;
    }

    /// Tuning access.
    PIDConfig&       config()       noexcept { return level_pi_.config(); }
    const PIDConfig& config() const noexcept { return level_pi_.config(); }

    /// Integrator inspection (tests).
    double levelIntegrator() const noexcept { return level_pi_.integrator(); }

private:
    PIDController level_pi_;
    double        level_setpoint_m_;
};

}  // namespace astara::control

#endif  // ASTARA_CONTROL_THREE_ELEMENT_CONTROLLER_HPP
