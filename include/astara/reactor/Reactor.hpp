#ifndef ASTARA_REACTOR_REACTOR_HPP
#define ASTARA_REACTOR_REACTOR_HPP

/**
 * @file   Reactor.hpp
 * @brief  Composite reactor: point kinetics + Mann-style thermal hydraulics
 *         + temperature-feedback reactivity coupling.
 *
 * Composes `PointKinetics` and `ReactorThermal` and supplies the temperature
 * feedback that closes the loop:
 *
 * @f[
 *   \rho(t) \;=\; \rho_{ext}(t)
 *               \;+\; \alpha_F \,\bigl(\bar T_F(t) - \bar T_{F,0}\bigr)
 *               \;+\; \alpha_M \,\bigl(\bar T_M(t) - \bar T_{M,0}\bigr)
 * @f]
 *
 * where `bar T_F` and `bar T_M` are the spatial averages of the fuel and
 * moderator temperatures, and `alpha_F`, `alpha_M` are the fuel and moderator
 * temperature coefficients of reactivity (typically negative for a PWR).
 *
 * The reactor exposes `timeStep(dt)` for advancing one classical RK4 step.
 * Within the step, point kinetics and thermal hydraulics are advanced
 * **together** using a single shared RK4 -- this is essential for stability
 * during fast transients where the fuel-temperature feedback couples back
 * into the prompt jump.
 *
 * @cite Naghedolfeizi (1990), Section 3.1, eqs. (3.1)-(3.16).
 */

#include "astara/reactor/PointKinetics.hpp"
#include "astara/reactor/ReactorThermal.hpp"

namespace astara::reactor {

/**
 * @brief Temperature-feedback (and external) reactivity model.
 */
struct ReactivityModel {
    double alpha_fuel_per_K      = 0.0;   ///< fuel temp coefficient,      1/K (typically negative)
    double alpha_moderator_per_K = 0.0;   ///< moderator temp coefficient, 1/K (typically negative)
    double T_fuel_reference_K      = 0.0; ///< bar T_F at which feedback = 0
    double T_moderator_reference_K = 0.0; ///< bar T_M at which feedback = 0
    double rho_external          = 0.0;   ///< externally-imposed reactivity (e.g. control rods)

    /// Total feedback reactivity at given averaged fuel and moderator T.
    double evaluate(double T_fuel_avg_K, double T_mod_avg_K) const noexcept {
        return rho_external
             + alpha_fuel_per_K      * (T_fuel_avg_K - T_fuel_reference_K)
             + alpha_moderator_per_K * (T_mod_avg_K  - T_moderator_reference_K);
    }
};

/**
 * @brief Combined reactor state: kinetics + thermal.
 */
struct ReactorState {
    double t_s = 0.0;                  ///< simulation time, s
    PointKineticsState   kinetics;     ///< (n, c_1, ..., c_G)
    ReactorThermalState  thermal;      ///< T_f, T_m, T_cl, T_lp, T_up, T_hl

    /// Vector-space arithmetic used by RK4 (excludes time, which advances
    /// monotonically and is owned by the integrator).
    ReactorState& operator+=(const ReactorState& rhs);
    ReactorState& operator-=(const ReactorState& rhs);
    ReactorState& operator*=(double s) noexcept;
};

inline ReactorState operator+(ReactorState a, const ReactorState& b) { a += b; return a; }
inline ReactorState operator-(ReactorState a, const ReactorState& b) { a -= b; return a; }
inline ReactorState operator*(ReactorState a, double s) noexcept     { a *= s; return a; }
inline ReactorState operator*(double s, ReactorState a) noexcept     { a *= s; return a; }

/**
 * @brief Top-level Reactor model.
 *
 * Holds parameters and current state.  The state is advanced in place by
 * `timeStep(dt)`.  External reactivity is set on the `reactivity()` model
 * (e.g. by a controller) between time steps.
 *
 * Power scaling: the point-kinetics state `n` is **normalised** (n = 1 at
 * full rated power).  The instantaneous thermal power deposited in the core
 * is `n * rated_thermal_power_W` and is what `ReactorThermal` consumes.
 */
class Reactor {
public:
    /**
     * @brief Construct from kinetics constants and thermal parameters.
     * @param  groups                six-group constants (or one-group average)
     * @param  thermal_parameters    Mann nodalisation + masses + heat transfer
     * @param  rated_thermal_power_W full-power thermal output, W
     * @param  reactivity            initial feedback model (often left default)
     */
    Reactor(DelayedGroupConstants    groups,
            ReactorThermalParameters thermal_parameters,
            double                   rated_thermal_power_W,
            ReactivityModel          reactivity = {});

    // ---------- Initialisation ----------

    /**
     * @brief Initialise the reactor at steady state at given normalised power
     *        and inlet temperature.
     *
     * - precursors set to their steady-state values for the given power
     * - thermal field set by `steadyStateThermal`
     * - reactivity reference temperatures set to the resulting averages so
     *   that initial feedback reactivity = rho_external (= 0 by default)
     *
     * @param n0           normalised power (1.0 = full power)
     * @param T_inlet_K    cold-leg-in temperature, K
     */
    void initialiseSteadyState(double n0, double T_inlet_K);

    // ---------- Inputs ----------

    /// Set the temperature of the fluid returning from the SG into the cold leg.
    void setColdLegInletTemperature(double T_K) noexcept { T_inlet_K_ = T_K; }
    /// Read the cold-leg-in temperature currently being applied as boundary.
    double coldLegInletTemperature() const noexcept { return T_inlet_K_; }

    /// External reactivity (in absolute units, not pcm).  Apply +1 cent of
    /// reactivity by passing 0.01 * beta_total.
    ReactivityModel&       reactivity()       noexcept { return reactivity_; }
    const ReactivityModel& reactivity() const noexcept { return reactivity_; }

    // ---------- Time stepping ----------

    /**
     * @brief Advance the reactor state one RK4 time step of size `dt`.
     * @throws std::runtime_error if `dt <= 0` or the resulting state is non-finite.
     */
    void timeStep(double dt);

    // ---------- Outputs ----------

    const ReactorState&             state()      const noexcept { return state_; }
    const DelayedGroupConstants&    groups()     const noexcept { return groups_; }
    const ReactorThermalParameters& thermalParams() const noexcept { return tp_; }
    double                          ratedPowerW() const noexcept { return P_rated_W_; }

    /// Instantaneous thermal power, W.
    double thermalPowerW() const noexcept { return state_.kinetics.power() * P_rated_W_; }
    /// Hot-leg-out temperature (the port that feeds the SG primary inlet), K.
    double hotLegOutletTemperatureK() const noexcept { return state_.thermal.T_hot_leg; }

    /// Compute the right-hand side of the full reactor ODE system at the given
    /// state.  Public so test code can inspect derivatives without stepping.
    ReactorState evaluateDerivative(const ReactorState& s) const;

private:
    DelayedGroupConstants     groups_;
    ReactorThermalParameters  tp_;
    double                    P_rated_W_;
    ReactivityModel           reactivity_;
    ReactorState              state_;
    double                    T_inlet_K_ = 0.0;
};

}  // namespace astara::reactor

#endif  // ASTARA_REACTOR_REACTOR_HPP
