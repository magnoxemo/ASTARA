#ifndef ASTARA_CORE_INTEGRATOR_HPP
#define ASTARA_CORE_INTEGRATOR_HPP

/**
 * @file   Integrator.hpp
 * @brief  Explicit time integrators (forward Euler, classical RK4).
 *
 * The integrators are *function templates*, not classes -- they take a state
 * `y0`, a step size `dt`, and a derivative-evaluating callable `f` such that
 * `f(t, y) -> y_dot`, and return the advanced state.  This decouples the
 * integrator from any particular component class and makes them trivially
 * unit-testable against analytical solutions:
 *
 * @code
 * // dy/dt = -k*y, exact: y(t) = y0*exp(-k*t)
 * auto f = [k](double, double y){ return -k * y; };
 * double y = rk4Step(0.0, 1.0, dt, f);
 * @endcode
 *
 * The state type `Y` must support the vector-space operations
 * `Y + Y`, `Y - Y`, and `Y * double`.  Both `FixedState<N>` and `DynamicState`
 * (see State.hpp) satisfy these.  Plain `double` works too.
 */

#include <type_traits>
#include <utility>

namespace astara::core {

/**
 * @brief Forward Euler step.  First-order accurate; useful only for testing.
 *
 * @tparam Y      state type (must support `Y + Y` and `Y * double`)
 * @tparam F      callable with signature `Y(double t, const Y& y)`
 * @param  t     current time
 * @param  y     current state
 * @param  dt    step size
 * @param  f     derivative function (returns dy/dt at given t,y)
 * @return       state at time `t + dt`
 */
template <class Y, class F>
inline Y eulerStep(double t, const Y& y, double dt, F&& f) {
    const Y k1 = f(t, y);
    return y + k1 * dt;
}

/**
 * @brief Classical fourth-order Runge-Kutta step.
 *
 * Standard "RK4" with weights (1, 2, 2, 1)/6.  Fourth-order accurate, but not
 * adaptive; the user is responsible for picking a `dt` small enough for the
 * fastest dynamics in the system.
 *
 * For the PWR primary loop, the prompt-neutron timescale (Λ ~ 1e-5 s) sets
 * the upper bound on `dt` during fast transients; the moderator thermal-lag
 * (τ ~ 1 s) sets it during slow transients.  A typical safe choice for full
 * primary-loop transients is `dt = 0.01 s`.
 *
 * @tparam Y     state type (must support `+`, `*double`)
 * @tparam F     callable with signature `Y(double t, const Y& y)`
 * @param  t    current time
 * @param  y    current state
 * @param  dt   step size
 * @param  f    derivative function
 * @return      state at time `t + dt`
 */
template <class Y, class F>
inline Y rk4Step(double t, const Y& y, double dt, F&& f) {
    const double half_dt = 0.5 * dt;

    const Y k1 = f(t,            y);
    const Y k2 = f(t + half_dt,  y + k1 * half_dt);
    const Y k3 = f(t + half_dt,  y + k2 * half_dt);
    const Y k4 = f(t + dt,       y + k3 * dt);

    // y_{n+1} = y_n + dt/6 * (k1 + 2 k2 + 2 k3 + k4)
    return y + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (dt / 6.0);
}

/**
 * @brief Enumeration of supported integrators (used for runtime selection).
 */
enum class IntegratorKind {
    Euler,        ///< First-order forward Euler -- testing only
    RK4,          ///< Classical fourth-order Runge-Kutta -- production default
};

/**
 * @brief Dispatch helper: take a single step using the chosen scheme.
 *
 * Equivalent to calling `eulerStep` or `rk4Step` directly, but lets test code
 * sweep schemes by passing an enum.
 */
template <class Y, class F>
inline Y integrateStep(IntegratorKind kind, double t, const Y& y, double dt, F&& f) {
    switch (kind) {
        case IntegratorKind::Euler: return eulerStep(t, y, dt, std::forward<F>(f));
        case IntegratorKind::RK4:   return rk4Step  (t, y, dt, std::forward<F>(f));
    }
    // Unreachable; placate compilers that warn on missing return.
    return y;
}

}  // namespace astara::core

#endif  // ASTARA_CORE_INTEGRATOR_HPP
