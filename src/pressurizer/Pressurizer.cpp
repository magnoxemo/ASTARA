/**
 * @file   Pressurizer.cpp
 * @brief  Implementation of the two-region pressurizer model.
 */

#include "astara/pressurizer/Pressurizer.hpp"
#include "astara/props/WaterProperties.hpp"
#include "astara/core/Integrator.hpp"

#include <cmath>
#include <sstream>
#include <stdexcept>

namespace astara::pressurizer {

Pressurizer::Pressurizer(PressurizerParameters params,
                         const props::WaterProperties* p)
        : p_(params), props_(p) {
    p_.validate();
    if (props_ == nullptr) {
        throw std::invalid_argument("Pressurizer: water properties service may not be null");
    }
}

void Pressurizer::initialiseSteadyState(double P_Pa, double L_w_m) {
    if (!(P_Pa > 0.0))           throw std::invalid_argument("Pressurizer::initialiseSteadyState: P > 0 required");
    if (!(L_w_m > 0.0))          throw std::invalid_argument("Pressurizer::initialiseSteadyState: Lw > 0 required");
    if (!(L_w_m < p_.total_height_m)) {
        throw std::invalid_argument("Pressurizer::initialiseSteadyState: Lw must be less than total height");
    }
    state_.t_s            = 0.0;
    state_.pressure_Pa    = P_Pa;
    state_.water_level_m  = L_w_m;
    inputs_ = {};                       // no surge, no spray, heater off
    inputs_.spray_enthalpy_J_kg = props_->satLiquidEnthalpy_P(P_Pa);  // sane default
    inputs_.surge_enthalpy_J_kg = props_->satLiquidEnthalpy_P(P_Pa);  // sane default
}

double Pressurizer::waterMass_kg() const noexcept {
    const double rho_w = props_->satLiquidDensity_P(state_.pressure_Pa);
    return rho_w * p_.cross_section_area_m2 * state_.water_level_m;
}

double Pressurizer::steamMass_kg() const noexcept {
    const double rho_s = props_->satVapourDensity_P(state_.pressure_Pa);
    return rho_s * p_.cross_section_area_m2 * (p_.total_height_m - state_.water_level_m);
}

PressurizerState Pressurizer::evaluateDerivative(const PressurizerState& s) const {
    const double P    = s.pressure_Pa;
    const double Lw   = s.water_level_m;
    const double Apr  = p_.cross_section_area_m2;
    const double L    = p_.total_height_m;
    const double Vw   = Apr * Lw;
    const double Vs   = Apr * (L - Lw);

    if (Vs <= 0.0) {
        throw std::runtime_error(
                "Pressurizer: water level has reached the top of the vessel "
                "(Vs <= 0); model is no longer valid -- this means the entire "
                "vessel filled with water.  Either model relief valve action or "
                "stop the simulation here.");
    }
    if (Vw <= 0.0) {
        throw std::runtime_error(
                "Pressurizer: water level has reached zero (vessel emptied); "
                "model is no longer valid.");
    }

    // Saturation properties at current pressure.  We use central differences
    // to approximate the dphi/dP partial derivatives that close the
    // compressibility equation (eqs. C.13-C.18 of the thesis use these as
    // K_p* coefficients of the linearised property model; with IF97 we
    // compute them numerically).
    constexpr double dP = 5.0e3;     // 5 kPa: small enough that the slope is
                                     // accurate, large enough that IF97
                                     // numerical noise is averaged out.
    const double rho_w = props_->satLiquidDensity_P(P);
    const double rho_s = props_->satVapourDensity_P(P);
    const double h_f   = props_->satLiquidEnthalpy_P(P);
    const double h_g   = props_->satVapourEnthalpy_P(P);
    const double h_fg  = h_g - h_f;
    const double v_f   = 1.0 / rho_w;
    const double v_g   = 1.0 / rho_s;

    auto slope = [&](auto fn) {
        const double v_lo = fn(P - dP);
        const double v_hi = fn(P + dP);
        return (v_hi - v_lo) / (2.0 * dP);
    };
    const double dvg_dP   = slope([&](double pp){ return 1.0 / props_->satVapourDensity_P(pp); });
    const double dvf_dP   = slope([&](double pp){ return 1.0 / props_->satLiquidDensity_P(pp); });
    const double dhg_dP   = slope([&](double pp){ return props_->satVapourEnthalpy_P(pp); });
    const double dhf_dP   = slope([&](double pp){ return props_->satLiquidEnthalpy_P(pp); });

    // Inputs.
    const double Wsg = inputs_.surge_mass_flow_kg_s;
    const double Wsp = inputs_.spray_mass_flow_kg_s;
    const double Q   = inputs_.heater_power_W;
    const double hsp = inputs_.spray_enthalpy_J_kg;
    // Surge enthalpy: when surge is *into* the pressurizer the water carries
    // the hot-leg enthalpy supplied by the user; when it surges *out* the
    // mass leaving the pressurizer carries the saturated-liquid enthalpy.
    const double h_surge_in = (Wsg >= 0.0) ? inputs_.surge_enthalpy_J_kg : h_f;

    // -----------------------------------------------------------------------
    // Pressure derivative.
    //
    // The model treats both phases as saturated.  The simplification we use
    // (consistent with the thesis once the linearised properties are replaced
    // by IF97 finite differences) is the energy balance over the steam
    // volume:
    //
    //   M_s * dh_g/dt + M_w * dh_f/dt - V_total * dP/dt
    //         = Q + W_sp * (h_sp - h_f) + W_sg * (h_surge - h_f) + W_co * h_fg
    //
    // Combined with the mass-balance evaporation rate W_co
    // (eq. C.8 of the thesis), and with d(rho_x)/dt = (drho_x/dP) dP/dt for
    // saturated phases, this collapses to:
    //
    //   dP/dt = NUM / DEN
    //
    // where NUM contains the source terms (Q, surge, spray) and DEN contains
    // the compressibility "stiffness" of the two-phase system.  We assemble
    // it term-by-term below and keep them in named locals for readability
    // and for the analytical-test verification (the steady-state pressure
    // controller test compares dP/dt computed here to a finite difference
    // of the simulator).
    // -----------------------------------------------------------------------

    // Mass balance:  W_co (condensation/evaporation in the steam dome).
    // From eqs. (C.7)-(C.8), and using the notation C1 = rho_w/rho_s - 1:
    const double C1 = (rho_w / rho_s) - 1.0;
    if (!(std::abs(C1) > 1e-9)) {
        throw std::runtime_error("Pressurizer: rho_w == rho_s (critical point); "
                                 "model invalid");
    }
    // Compressibility coefficient (eq. C.11 with linearised properties; here
    // computed from finite differences):
    //    C2 = A_pr * (L - L_w) * (rho_w / rho_s) * d(rho_s)/dP
    //       + A_pr *  L_w     *                       d(rho_w)/dP
    // where d(rho_s)/dP = -(1/v_g^2) * dv_g/dP, and similarly for rho_w.
    const double drhog_dP = -(1.0 / (v_g * v_g)) * dvg_dP;
    const double drhof_dP = -(1.0 / (v_f * v_f)) * dvf_dP;
    const double C2 = Apr * (L - Lw) * (rho_w / rho_s) * drhog_dP
                    + Apr *      Lw  *                   drhof_dP;

    // The pressure ODE in compact form, using a numerical compressibility
    // coefficient kappa = M_s * dv_g/dP + M_w * dv_f/dP that captures the
    // total volumetric "softness" of the contents:
    const double M_w = rho_w * Vw;
    const double M_s = rho_s * Vs;
    const double kappa = M_s * dvg_dP + M_w * dvf_dP;
    if (!(std::abs(kappa) > 1e-30)) {
        throw std::runtime_error("Pressurizer: degenerate compressibility (kappa==0)");
    }

    // Volumetric balance:  V_w + V_s = V_total = const
    //   d(V_w)/dt + d(V_s)/dt = 0
    //   => A_pr * dLw/dt + d(M_s v_g)/dt = 0   (since V_s = M_s v_g)
    //   => A_pr * dLw/dt + (W_co) v_g + M_s dv_g/dP * dP/dt = 0
    //
    // Independently from mass balance:
    //   dM_w/dt = W_sg + W_sp - W_co = d(rho_w V_w)/dt
    //   = drho_w/dP * V_w * dP/dt + rho_w * A_pr * dLw/dt
    // Solving these together for dP/dt gives:
    //
    //   dP/dt = (W_sg + W_sp - rho_w (Wsg+Wsp-W_co)/rho_w * 0 ... )  -- expansion below
    //
    // For clarity we use the energy balance directly.  Energy balance on the
    // entire pressurizer (water + steam):
    //
    //   d/dt (M_w h_f + M_s h_g - P V_total) = Q + W_sg h_surge + W_sp h_sp
    //                                          - W_lo h_f
    //
    // where W_lo is mass leaving the surge line in out-surge (which we
    // already accounted for by using h_surge that switches at Wsg sign).
    //
    // Expanding the LHS (using h_f, h_g, M_w, M_s as functions of P only,
    // since the phases are saturated, plus dependence on Lw via M_w, M_s):
    //
    //   = (M_w dh_f/dP + M_s dh_g/dP) dP/dt
    //     + h_f dM_w/dt + h_g dM_s/dt - V_total dP/dt
    //
    //   = [M_w dh_f/dP + M_s dh_g/dP - V_total] dP/dt
    //     + h_f (W_sg + W_sp - W_co) + h_g W_co
    //
    //   = [M_w dh_f/dP + M_s dh_g/dP - V_total] dP/dt
    //     + h_f (W_sg + W_sp) + (h_g - h_f) W_co
    //
    // Source terms = Q + W_sg h_surge + W_sp h_sp - W_sg h_f (when out-surge)
    // (covered by h_surge_in switching).  We rearrange:
    //
    //   [M_w dh_f/dP + M_s dh_g/dP - V_total] dP/dt
    //     = Q + W_sg (h_surge_in - h_f) + W_sp (h_sp - h_f) - h_fg * W_co
    //
    // We still need W_co.  Substituting from mass balance into the volumetric
    // closure yields, after some algebra:
    //
    //   W_co = ((1/C1) * dP/dt * C2_vol - W_sg - W_sp) / 1     (linearised)
    //
    // To keep the implementation closed-form *with* IF97 properties we use
    // the simpler approximation that on the timescale of the pressure
    // dynamics in a PWR (~minutes), W_co is small except during fast
    // transients, and we can approximate W_co by the volumetric-closure
    // expression evaluated at the current dP/dt = 0 estimate; then we
    // iterate once.  In practice for the two-state pressurizer, the
    // following non-iterative form is used by Naghedolfeizi (eqs. C.19,
    // C.20) and is what we implement:
    //
    //    DEN = M_w dh_f/dP + M_s dh_g/dP - V_total
    //          + (h_fg / C1) * (- C2)        (sign per thesis eq. C.20)
    //    NUM = Q + W_sp (h_sp - h_f) + W_sg (h_surge - h_f)
    //          + (h_fg / C1) * (W_sg + W_sp)
    //
    // Cross-checked against the dimensions: every term is W/(Pa) or W,
    // giving dP/dt in Pa/s.  Good.
    //
    // (See thesis eqs. C.19-C.20 for the original derivation.)
    const double DEN = M_w * dhf_dP + M_s * dhg_dP - (Vw + Vs)
                       - (h_fg / C1) * C2;
    const double NUM = Q
                       + Wsp * (hsp        - h_f)
                       + Wsg * (h_surge_in - h_f)
                       + (h_fg / C1) * (Wsg + Wsp);

    if (!(std::abs(DEN) > 1e-30)) {
        throw std::runtime_error("Pressurizer: degenerate energy-balance denominator");
    }
    const double dP_dt = NUM / DEN;

    // Water level.  From the volumetric closure:
    //   A_pr * dLw/dt = -V_s * (1/v_g) * dv_g/dP * dP/dt - W_co * v_g
    // and W_co = (1/C1) * (C2 * dP/dt - W_sg - W_sp)   (eq. C.8)
    const double W_co = (1.0 / C1) * (C2 * dP_dt - Wsg - Wsp);
    // dLw/dt:
    //   d/dt (rho_w V_w) = W_sg + W_sp - W_co
    //   drho_w/dP * V_w * dP/dt + rho_w * A_pr * dLw/dt = W_sg + W_sp - W_co
    const double dLw_dt = ((Wsg + Wsp - W_co) - drhof_dP * Vw * dP_dt)
                          / (rho_w * Apr);

    PressurizerState dy;
    dy.t_s           = 0.0;
    dy.pressure_Pa   = dP_dt;
    dy.water_level_m = dLw_dt;
    return dy;
}

void Pressurizer::timeStep(double dt) {
    if (!(dt > 0.0)) {
        throw std::invalid_argument("Pressurizer::timeStep: dt must be positive");
    }
    auto deriv = [this](double /*t*/, const PressurizerState& s) {
        return this->evaluateDerivative(s);
    };
    state_ = astara::core::rk4Step(state_.t_s, state_, dt, deriv);
    state_.t_s += dt;

    if (!std::isfinite(state_.pressure_Pa) || !std::isfinite(state_.water_level_m)) {
        throw std::runtime_error("Pressurizer::timeStep produced non-finite state");
    }
    if (state_.water_level_m <= 0.0 || state_.water_level_m >= p_.total_height_m) {
        std::ostringstream os;
        os << "Pressurizer: water level out of valid range (0, " << p_.total_height_m
           << ") m, got " << state_.water_level_m
           << " -- physical limit reached, model no longer valid.";
        throw std::runtime_error(os.str());
    }
}

}  // namespace astara::pressurizer
