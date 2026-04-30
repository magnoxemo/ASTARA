#ifndef ASTARA_PUMP_HOMOLOGOUS_PUMP_HPP
#define ASTARA_PUMP_HOMOLOGOUS_PUMP_HPP

/**
 * @file   HomologousPump.hpp
 * @brief  Centrifugal coolant-pump model using homologous theory.
 *
 * Reproduces the model in Naghedolfeizi (1990), Section 3.4: a single
 * centrifugal pump with a quadratic head-flow characteristic at rated speed,
 * scaled to other speeds via the pump similarity (homologous) laws:
 *
 * @f[
 *   \frac{H_p}{H_s} = \left(\frac{N_p}{N_s}\right)^2,
 *   \qquad
 *   \frac{Q_p}{Q_s} = \frac{N_p}{N_s}.
 * @f]
 *
 * The pump characteristic is fitted as
 * @f[
 *   H_s(Q_s) \;=\; A_2 Q_s^2 + A_1 Q_s + A_0
 * @f]
 * over the operating range, with the coefficients (A_0, A_1, A_2) supplied
 * at construction.  Combined with the homologous laws this gives
 * @f[
 *   H_p(Q_p, N_p) \;=\; A_2 Q_p^2
 *                       + A_1 Q_p \,(N_p / N_s)
 *                       + A_0    \,(N_p / N_s)^2.
 * @f]
 *
 * The momentum equation around the primary loop, assuming a single effective
 * loop-pressure-drop coefficient `K_loop` such that `H_loss = K_loop * Q^2`:
 *
 * @f[
 *    \rho g L A_{ef}\, \frac{d}{dt}\!\left(\frac{Q}{A_{ef}}\right)
 *      \;=\; \rho g A_{ef}\,(H_p(Q, N_p) - K_{loop} Q^2),
 * @f]
 *
 * which simplifies to
 *
 * @f[
 *   \frac{dQ}{dt} \;=\; \frac{g A_{ef}}{L}\,(H_p - K_{loop} Q^2).
 * @f]
 *
 * The torque equation:
 *
 * @f[
 *   2 \pi I \, \frac{dN_p}{dt}
 *      \;=\; T_d - T_h,
 *      \quad\text{with}\quad
 *   T_d = \frac{P_d}{2 \pi N_p},\quad
 *   T_h = \frac{\rho g Q_p H_p}{2 \pi N_p}.
 * @f]
 *
 * Substituting:
 *
 * @f[
 *   \frac{dN_p}{dt}
 *      \;=\; \frac{1}{(2\pi)^2 I N_p}\,
 *            \bigl(P_d - \rho g Q_p H_p\bigr).
 * @f]
 *
 * State: `(Q, N)` with Q in m^3/s (volumetric flow) and N in rev/s
 * (rotational frequency, **not** rad/s).  Mass flow rate is recovered as
 * `mdot = rho * Q`.
 *
 * @cite Naghedolfeizi (1990), Section 3.4, eqs. (3.50)-(3.62), Fig. 3.22-3.23.
 * @cite Streeter, V. L., and Wylie, E. B. (1967). "Hydraulic Transients."
 *       McGraw-Hill.  (Source of the homologous pump theory.)
 */

#include <stdexcept>

namespace astara::pump {

/**
 * @brief Quadratic fit to the pump characteristic at rated speed.
 *
 * H = A2 * Q^2 + A1 * Q + A0   (Q in m^3/s, H in m).  For a centrifugal
 * pump A2 < 0 (head decreases as Q^2), A1 has either sign, and A0 > 0
 * is the shutoff head.
 */
struct PumpCurveCoefficients {
    double A0 = 0.0;   ///< shutoff (zero-flow) head, m
    double A1 = 0.0;   ///< linear coefficient, m / (m^3/s)
    double A2 = 0.0;   ///< quadratic coefficient, m / (m^3/s)^2
};

/**
 * @brief Geometric, inertial, and rated parameters of the loop + pump.
 */
struct HomologousPumpParameters {
    PumpCurveCoefficients curve;          ///< rated-speed head-flow fit

    double rated_speed_rev_s    = 0.0;    ///< N_s, rated rotational frequency (rev/s)
    double rated_volumetric_flow_m3_s = 0.0;  ///< Q_s, used for sanity checks

    double loop_resistance_K_s2_m5 = 0.0; ///< K_loop, head-loss coefficient
                                          ///< such that H_loss = K * Q^2

    double effective_flow_area_m2  = 0.0; ///< A_ef, loop cross-section for inertia
    double loop_length_m           = 0.0; ///< L, hydraulic length of the loop
    double fluid_density_kg_m3     = 0.0; ///< rho, fluid density used in head/torque

    double moment_of_inertia_kg_m2 = 0.0; ///< I, pump+motor combined inertia

    double rated_input_power_W     = 0.0; ///< P_d at rated conditions
                                          ///< (driven externally between time steps)

    double gravity_m_s2 = 9.80665;        ///< g (default = standard gravity)

    /// Throw on any non-positive required parameter.
    void validate() const;
};

/**
 * @brief State vector of the pump + loop hydraulic dynamics.
 */
struct HomologousPumpState {
    double t_s             = 0.0;  ///< time, s
    double volumetric_flow = 0.0;  ///< Q, m^3/s
    double speed_rev_s     = 0.0;  ///< N_p, rev/s

    /// Vector-space arithmetic used by RK4 (excludes time).
    HomologousPumpState& operator+=(const HomologousPumpState& rhs) noexcept {
        volumetric_flow += rhs.volumetric_flow;
        speed_rev_s     += rhs.speed_rev_s;
        return *this;
    }
    HomologousPumpState& operator-=(const HomologousPumpState& rhs) noexcept {
        volumetric_flow -= rhs.volumetric_flow;
        speed_rev_s     -= rhs.speed_rev_s;
        return *this;
    }
    HomologousPumpState& operator*=(double s) noexcept {
        volumetric_flow *= s;
        speed_rev_s     *= s;
        return *this;
    }
};

inline HomologousPumpState operator+(HomologousPumpState a, const HomologousPumpState& b) noexcept { a += b; return a; }
inline HomologousPumpState operator-(HomologousPumpState a, const HomologousPumpState& b) noexcept { a -= b; return a; }
inline HomologousPumpState operator*(HomologousPumpState a, double s) noexcept                      { a *= s; return a; }
inline HomologousPumpState operator*(double s, HomologousPumpState a) noexcept                      { a *= s; return a; }

/**
 * @brief Centrifugal coolant pump with quadratic curve and homologous scaling.
 */
class HomologousPump {
public:
    explicit HomologousPump(HomologousPumpParameters params);

    /**
     * @brief Initialise to steady state for the rated input power.
     *
     * Solves H_p(Q, N_s) = K_loop * Q^2 (=> the operating point on the rated
     * curve) and sets `state.speed = rated_speed`.  Throws if no positive root
     * exists (configuration error).
     */
    void initialiseAtRated();

    /**
     * @brief Set the input shaft power supplied to the pump (W).
     *
     * The user (or the controller) drives this between time steps.  Setting
     * P_d = 0 represents a pump trip.
     */
    void setInputPowerW(double Pd) noexcept { Pd_W_ = Pd; }

    /// Single RK4 step.
    void timeStep(double dt);

    // ---------- Outputs ----------
    const HomologousPumpState&    state()  const noexcept { return state_; }
    const HomologousPumpParameters& params() const noexcept { return p_; }

    /// Pump head at the current operating point (m).
    double developedHead_m() const noexcept;

    /// Mass flow rate through the loop (kg/s).
    double massFlow_kg_s() const noexcept { return p_.fluid_density_kg_m3 * state_.volumetric_flow; }

    /// Compute the right-hand side at given state (used by tests + integrator).
    HomologousPumpState evaluateDerivative(const HomologousPumpState& s) const;

    /// Compute the developed head at arbitrary (Q, N) using the homologous laws.
    double headAt(double Q_m3_s, double N_rev_s) const noexcept;

private:
    HomologousPumpParameters p_;
    HomologousPumpState      state_;
    double                   Pd_W_ = 0.0;  ///< current input power
};

}  // namespace astara::pump

#endif  // ASTARA_PUMP_HOMOLOGOUS_PUMP_HPP
