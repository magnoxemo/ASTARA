#ifndef ASTARA_CONTROL_CONTROL_ROD_SYSTEM_HPP
#define ASTARA_CONTROL_CONTROL_ROD_SYSTEM_HPP

/**
 * @file   ControlRodSystem.hpp
 * @brief  Open-loop manual control-rod reactivity-rate source.
 *
 * Models a bank of control rods driven at a fixed insertion/withdrawal speed,
 * producing a constant external-reactivity rate while in that mode:
 *
 * @f[
 *    \frac{d\rho_{ex}}{dt} =
 *      \begin{cases}
 *        -v_{ins}\, \Delta\rho_{step} & \text{Insert} \\
 *        +v_{wd}\,  \Delta\rho_{step} & \text{Withdraw} \\
 *        0                            & \text{Hold}
 *      \end{cases}
 * @f]
 *
 * Sign convention: control rods are neutron absorbers, so driving them
 * *into* the core removes reactivity (negative rate) and withdrawing them
 * adds reactivity (positive rate). The legacy Python prototype this class
 * replaces returned a positive rate for both modes -- a sign bug, fixed
 * here.
 *
 * This is a simple bang-bang rod-speed model (no rod-worth curve, no axial
 * position tracking) -- the caller is responsible for integrating the
 * returned rate into the reactor's total external reactivity, the same
 * contract as the legacy prototype this replaces. It complements the
 * closed-loop `ReactorController` (T_avg-tracking automatic control): use
 * this class instead when a transient calls for scripted/manual rod motion
 * (e.g. a manual rod withdrawal test, or an ATWS scenario with rods held).
 *
 * @cite Naghedolfeizi (1990), control-rod reactivity model (VVER-1200
 *       default step-reactivity/speed values).
 */

namespace astara::control {

/// Commanded control-rod motion.
enum class RodMode { Insert, Withdraw, Hold };

/**
 * @brief Configuration for a `ControlRodSystem`.
 */
struct ControlRodSystemConfig {
    double insertion_speed  = 0.5;    ///< normalised insertion rate, 1/s
    double withdrawal_speed = 0.4;    ///< normalised withdrawal rate, 1/s
    double step_reactivity  = 0.002;  ///< reactivity worth per unit rate, dimensionless
};

/**
 * @brief Manual control-rod reactivity-rate source.
 */
class ControlRodSystem {
public:
    explicit ControlRodSystem(ControlRodSystemConfig cfg) : cfg_(cfg) {}

    /// Read/write configuration access.
    ControlRodSystemConfig&       config()       noexcept { return cfg_; }
    const ControlRodSystemConfig& config() const noexcept { return cfg_; }

    /**
     * @brief External-reactivity rate demanded by the given rod motion.
     * @param mode  Insert, Withdraw, or Hold
     * @return      d(rho_ex)/dt, dimensionless / s
     */
    double reactivityRate(RodMode mode) const noexcept {
        switch (mode) {
            case RodMode::Insert:   return -cfg_.insertion_speed  * cfg_.step_reactivity;
            case RodMode::Withdraw: return  cfg_.withdrawal_speed * cfg_.step_reactivity;
            case RodMode::Hold:     return 0.0;
        }
        return 0.0;
    }

private:
    ControlRodSystemConfig cfg_;
};

}  // namespace astara::control

#endif  // ASTARA_CONTROL_CONTROL_ROD_SYSTEM_HPP
