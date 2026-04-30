/**
 * @file   HomologousPump.cpp
 * @brief  Implementation of the homologous coolant-pump model.
 */

#include "astara/pump/HomologousPump.hpp"
#include "astara/core/Integrator.hpp"

#include <cmath>
#include <sstream>

namespace astara::pump {

void HomologousPumpParameters::validate() const {
    auto check_pos = [](double v, const char* name) {
        if (!(v > 0.0)) {
            std::ostringstream os;
            os << "HomologousPumpParameters: " << name
               << " must be positive, got " << v;
            throw std::invalid_argument(os.str());
        }
    };
    check_pos(rated_speed_rev_s,            "rated_speed_rev_s");
    check_pos(rated_volumetric_flow_m3_s,   "rated_volumetric_flow_m3_s");
    check_pos(loop_resistance_K_s2_m5,      "loop_resistance_K_s2_m5");
    check_pos(effective_flow_area_m2,       "effective_flow_area_m2");
    check_pos(loop_length_m,                "loop_length_m");
    check_pos(fluid_density_kg_m3,          "fluid_density_kg_m3");
    check_pos(moment_of_inertia_kg_m2,      "moment_of_inertia_kg_m2");
    check_pos(rated_input_power_W,          "rated_input_power_W");
    if (!(curve.A0 > 0.0)) {
        std::ostringstream os;
        os << "HomologousPumpParameters: shutoff head curve.A0 must be positive, got "
           << curve.A0;
        throw std::invalid_argument(os.str());
    }
    if (gravity_m_s2 <= 0.0) {
        throw std::invalid_argument("HomologousPumpParameters: gravity_m_s2 must be positive");
    }
}

HomologousPump::HomologousPump(HomologousPumpParameters params) : p_(params) {
    p_.validate();
    Pd_W_ = p_.rated_input_power_W;
}

double HomologousPump::headAt(double Q, double N) const noexcept {
    // Homologous scaling of the rated curve to operating speed:
    //   H(Q, N) = A2 Q^2 + A1 Q (N/Ns) + A0 (N/Ns)^2
    const double r = N / p_.rated_speed_rev_s;
    return p_.curve.A2 * Q * Q
         + p_.curve.A1 * Q * r
         + p_.curve.A0 * r * r;
}

double HomologousPump::developedHead_m() const noexcept {
    return headAt(state_.volumetric_flow, state_.speed_rev_s);
}

HomologousPumpState HomologousPump::evaluateDerivative(const HomologousPumpState& s) const {
    HomologousPumpState dy;

    const double H_p = headAt(s.volumetric_flow, s.speed_rev_s);

    // dQ/dt = (g A_ef / L) * (H_p - K_loop * Q^2)
    const double H_loss = p_.loop_resistance_K_s2_m5 * s.volumetric_flow * s.volumetric_flow;
    dy.volumetric_flow = (p_.gravity_m_s2 * p_.effective_flow_area_m2 / p_.loop_length_m)
                       * (H_p - H_loss);

    // dN/dt = (1 / ((2 pi)^2 I N)) * (P_d - rho g Q H_p).
    // Guard against N == 0 (divide-by-zero).  At zero speed the equation
    // becomes singular -- in practice we never run at exactly zero, but a
    // pump trip transient drives N very low, so we floor at a tiny positive
    // value (1e-3 rev/s = 0.06 RPM) below which the simulation is meaningless
    // anyway.
    constexpr double kMinN = 1.0e-3;
    const double N_safe = (s.speed_rev_s > kMinN) ? s.speed_rev_s : kMinN;
    const double inv_denom = 1.0 / ((2.0 * M_PI) * (2.0 * M_PI) * p_.moment_of_inertia_kg_m2 * N_safe);
    dy.speed_rev_s = inv_denom * (Pd_W_
                                  - p_.fluid_density_kg_m3 * p_.gravity_m_s2
                                    * s.volumetric_flow * H_p);
    return dy;
}

void HomologousPump::initialiseAtRated() {
    // At rated conditions the pump head equals the loop friction head:
    //   A2 Q^2 + A1 Q + A0 = K_loop Q^2
    // i.e. (A2 - K_loop) Q^2 + A1 Q + A0 = 0
    const double a = p_.curve.A2 - p_.loop_resistance_K_s2_m5;
    const double b = p_.curve.A1;
    const double c = p_.curve.A0;
    const double disc = b * b - 4.0 * a * c;
    if (disc < 0.0) {
        throw std::runtime_error("HomologousPump::initialiseAtRated: pump curve "
                                 "and loop resistance have no real intersection");
    }
    // Two roots; pick the positive one closest to the rated flow.
    const double sqrtD = std::sqrt(disc);
    const double q1 = (-b + sqrtD) / (2.0 * a);
    const double q2 = (-b - sqrtD) / (2.0 * a);
    double Q = -1.0;
    if (q1 > 0.0 && q2 > 0.0) {
        Q = (std::abs(q1 - p_.rated_volumetric_flow_m3_s)
           < std::abs(q2 - p_.rated_volumetric_flow_m3_s)) ? q1 : q2;
    } else if (q1 > 0.0) {
        Q = q1;
    } else if (q2 > 0.0) {
        Q = q2;
    } else {
        throw std::runtime_error("HomologousPump::initialiseAtRated: no positive "
                                 "root of pump-curve / loop-resistance intersection");
    }

    state_.t_s             = 0.0;
    state_.volumetric_flow = Q;
    state_.speed_rev_s     = p_.rated_speed_rev_s;
    Pd_W_                  = p_.rated_input_power_W;
}

void HomologousPump::timeStep(double dt) {
    if (!(dt > 0.0)) {
        throw std::invalid_argument("HomologousPump::timeStep: dt must be positive");
    }
    auto deriv = [this](double /*t*/, const HomologousPumpState& s) {
        return this->evaluateDerivative(s);
    };
    state_ = astara::core::rk4Step(state_.t_s, state_, dt, deriv);
    state_.t_s += dt;

    if (!std::isfinite(state_.volumetric_flow) || !std::isfinite(state_.speed_rev_s)) {
        throw std::runtime_error("HomologousPump::timeStep produced non-finite state");
    }
    // Floor the speed at the same minimum used in the derivative; otherwise
    // very small negative undershoots can occur with RK4 on stiff spindown.
    if (state_.speed_rev_s < 0.0) state_.speed_rev_s = 0.0;
    if (state_.volumetric_flow < 0.0) state_.volumetric_flow = 0.0;
}

}  // namespace astara::pump
