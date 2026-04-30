#ifndef ASTARA_REACTOR_POINT_KINETICS_HPP
#define ASTARA_REACTOR_POINT_KINETICS_HPP

/**
 * @file   PointKinetics.hpp
 * @brief  Point-reactor kinetics with N delayed-neutron precursor groups.
 *
 * Solves the standard point-kinetics equations
 *
 * @f[
 *   \frac{d n}{d t}
 *     = \frac{\rho - \beta}{\Lambda}\, n
 *       + \sum_{i=1}^{G} \lambda_i\, c_i
 * @f]
 * @f[
 *   \frac{d c_i}{d t}
 *     = \frac{\beta_i}{\Lambda}\, n
 *       - \lambda_i\, c_i
 *       \qquad i = 1, \dots, G
 * @f]
 *
 * with `n` the (normalised) neutron population (or thermal power), `c_i` the
 * i-th precursor concentration in equivalent power units, `beta_i` the i-th
 * delayed-neutron fraction (`beta = sum beta_i`), `lambda_i` the i-th decay
 * constant, `Lambda` the prompt-neutron generation time, and `rho` the total
 * (external + feedback) reactivity.
 *
 * The reactor temperature feedback (`alpha_f * dT_f`, `alpha_m * dT_m`) is
 * **not** evaluated here -- this class is purely the neutronics.  The
 * `Reactor` class composes `PointKinetics` with `ReactorThermal` and supplies
 * the total reactivity to `derivative()`.
 *
 * @cite Lamarsh, J. R. (1966). "Introduction to Nuclear Reactor Theory."
 *       Addison-Wesley.  Sect. 7-1.
 * @cite Naghedolfeizi (1990), eq. (3.1)-(3.3).
 */

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace astara::reactor {

/**
 * @brief A pre-computed set of group constants (lambda_i, beta_i, Lambda).
 *
 * Provided as a value type rather than constructor arguments so the same
 * constants can be shared between an instance and its tests.
 */
struct DelayedGroupConstants {
    std::vector<double> lambda;  ///< decay constants of each group, 1/s
    std::vector<double> beta;    ///< delayed-neutron fractions of each group, dimensionless
    double Lambda = 0.0;         ///< prompt-neutron generation time, s

    /// Total delayed-neutron fraction beta = sum(beta_i).
    double totalBeta() const noexcept;

    /// Throw `std::invalid_argument` if the constants are inconsistent
    /// (size mismatch, non-positive lambda, beta out of [0,1), Lambda <= 0).
    void validate() const;

    /// Six-group constants for U-235 thermal fission, ENDF/B-VI averages.
    /// These match the values used in the thesis Appendix A (Table A.1) and
    /// in the existing example program of the legacy code.
    static DelayedGroupConstants u235SixGroup() noexcept;

    /// Single-group equivalent (the thesis Section 3.1.1 uses this form).
    /// `lambda_eff = beta / sum_i (beta_i / lambda_i)`.
    static DelayedGroupConstants oneGroupAverageOf(const DelayedGroupConstants& g);
};

/**
 * @brief Point-kinetics state: power and one precursor per group.
 *
 * The state vector is laid out as `(n, c_1, c_2, ..., c_G)`.
 */
class PointKineticsState {
public:
    PointKineticsState() = default;
    explicit PointKineticsState(std::size_t G) : c_(G, 0.0) {}

    double                power() const noexcept { return n_; }
    double&               power()       noexcept { return n_; }

    /// Full vector of precursor concentrations.
    const std::vector<double>& precursors() const noexcept { return c_; }
    std::vector<double>&       precursors()       noexcept { return c_; }

    /// Number of delayed-neutron groups.
    std::size_t numGroups() const noexcept { return c_.size(); }

    /// Vector-space arithmetic required by `rk4Step`.
    PointKineticsState& operator+=(const PointKineticsState& rhs) {
        if (c_.size() != rhs.c_.size()) {
            throw std::invalid_argument("PointKineticsState size mismatch in +=");
        }
        n_ += rhs.n_;
        for (std::size_t i = 0; i < c_.size(); ++i) c_[i] += rhs.c_[i];
        return *this;
    }
    PointKineticsState& operator-=(const PointKineticsState& rhs) {
        if (c_.size() != rhs.c_.size()) {
            throw std::invalid_argument("PointKineticsState size mismatch in -=");
        }
        n_ -= rhs.n_;
        for (std::size_t i = 0; i < c_.size(); ++i) c_[i] -= rhs.c_[i];
        return *this;
    }
    PointKineticsState& operator*=(double s) noexcept {
        n_ *= s;
        for (auto& v : c_) v *= s;
        return *this;
    }

private:
    double n_ = 0.0;
    std::vector<double> c_;
};

inline PointKineticsState operator+(PointKineticsState a, const PointKineticsState& b) { a += b; return a; }
inline PointKineticsState operator-(PointKineticsState a, const PointKineticsState& b) { a -= b; return a; }
inline PointKineticsState operator*(PointKineticsState a, double s) noexcept           { a *= s; return a; }
inline PointKineticsState operator*(double s, PointKineticsState a) noexcept           { a *= s; return a; }

/**
 * @brief Compute steady-state precursor concentrations for a given power.
 *
 * At steady state, dc_i/dt = 0 implies c_i = (beta_i / (lambda_i * Lambda)) * n.
 * This is the standard initial condition used at the start of a transient
 * calculation -- otherwise the simulation will exhibit a non-physical
 * "kick" as the precursors equilibrate.
 *
 * @param  g   group constants
 * @param  n0  initial power (any unit; the precursors come out in the same)
 * @return     precursor concentrations
 */
std::vector<double> steadyStatePrecursors(const DelayedGroupConstants& g,
                                          double n0);

/**
 * @brief Pure derivative function for the point-kinetics equations.
 *
 * Stateless: takes the current state, the total reactivity rho, and the
 * group constants; returns dy/dt.  This is the function the time integrator
 * calls.
 *
 * @param  state   current (n, c_1, ..., c_G)
 * @param  rho     total reactivity (external + feedback) at this evaluation
 * @param  g       group constants
 * @return         dy/dt
 */
PointKineticsState pointKineticsDerivative(const PointKineticsState& state,
                                           double rho,
                                           const DelayedGroupConstants& g);

}  // namespace astara::reactor

#endif  // ASTARA_REACTOR_POINT_KINETICS_HPP
