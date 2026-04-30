/**
 * @file   HelicalCoilSteamGenerator.cpp
 * @brief  Implementation of the moving-boundary helical-coil SG of
 *         Arda & Holbert (2015).
 *
 * Equation numbers in comments refer to:
 *   S. E. Arda and K. E. Holbert, "A dynamic model of a passively cooled
 *   small modular reactor for controller design purposes," Nucl. Eng.
 *   Des. 289 (2015) 218-230.
 */

#include "astara/sg/HelicalCoilSteamGenerator.hpp"
#include "astara/props/WaterProperties.hpp"
#include "astara/core/Integrator.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace astara::sg {

// =============================================================================
// Parameter validation + default
// =============================================================================

void HelicalCoilSteamGeneratorParameters::validate() const {
    auto check_pos = [](double v, const char* name) {
        if (!(v > 0.0)) {
            std::ostringstream os;
            os << "HelicalCoilSteamGeneratorParameters: " << name
               << " must be positive, got " << v;
            throw std::invalid_argument(os.str());
        }
    };
    check_pos(tube_total_length_m,         "tube_total_length_m");
    check_pos(tube_inner_diameter_m,       "tube_inner_diameter_m");
    check_pos(tube_outer_diameter_m,       "tube_outer_diameter_m");
    check_pos(tube_count,                  "tube_count");
    check_pos(secondary_flow_area_m2,      "secondary_flow_area_m2");
    check_pos(primary_flow_area_m2,        "primary_flow_area_m2");
    check_pos(tube_metal_area_m2,          "tube_metal_area_m2");
    check_pos(tube_metal_density_kg_m3,    "tube_metal_density_kg_m3");
    check_pos(tube_metal_cp_J_per_kgK,     "tube_metal_cp_J_per_kgK");
    check_pos(primary_density_kg_m3,       "primary_density_kg_m3");
    check_pos(primary_cp_J_per_kgK,        "primary_cp_J_per_kgK");
    check_pos(alpha_o_subcooled_W_per_m2K, "alpha_o_subcooled_W_per_m2K");
    check_pos(alpha_o_twophase_W_per_m2K,  "alpha_o_twophase_W_per_m2K");
    check_pos(alpha_o_superheat_W_per_m2K, "alpha_o_superheat_W_per_m2K");
    check_pos(alpha_i_subcooled_W_per_m2K, "alpha_i_subcooled_W_per_m2K");
    check_pos(alpha_i_twophase_W_per_m2K,  "alpha_i_twophase_W_per_m2K");
    check_pos(alpha_i_superheat_W_per_m2K, "alpha_i_superheat_W_per_m2K");
}

HelicalCoilSteamGeneratorParameters
HelicalCoilSteamGeneratorParameters::nuscaleSMRTwoSG() {
    // NuScale SMR (Arda 2015 Tables 1, 4 + INEEL 2003).  Per-SG (one of two).
    HelicalCoilSteamGeneratorParameters p;
    // Geometry: 506 tubes, OD 16 mm, wall 0.9 mm, total length 22.25 m.
    p.tube_count                    = 506.0;
    p.tube_total_length_m           = 22.25;
    p.tube_outer_diameter_m         = 0.016;
    p.tube_inner_diameter_m         = 0.016 - 2.0 * 0.0009;     // = 0.0142
    // Equivalent flow areas (sum over all tubes for the secondary side;
    // the primary side is the annulus between the bundle and the riser
    // wall -- a representative cross-section ~ 1 m^2).
    const double pi = 3.14159265358979323846;
    p.secondary_flow_area_m2        = p.tube_count * pi * 0.25
                                      * p.tube_inner_diameter_m * p.tube_inner_diameter_m;
    p.primary_flow_area_m2          = 1.0;     // representative annular area
    p.tube_metal_area_m2            = p.tube_count * pi * 0.25
                                      * (p.tube_outer_diameter_m * p.tube_outer_diameter_m
                                         - p.tube_inner_diameter_m * p.tube_inner_diameter_m);
    // Materials.
    p.tube_metal_density_kg_m3      = 8190.0;     // Inconel-690
    p.tube_metal_cp_J_per_kgK       = 450.0;      // Inconel-690 ~ 450 J/kg-K
    p.primary_density_kg_m3         = 740.0;      // pressurised water at ~ 290 deg C, 12.7 MPa
    p.primary_cp_J_per_kgK          = 5350.0;
    // Heat-transfer coefficients (chosen to yield the published steady-state
    // temperature profile; ranges are consistent with Arda 2015 eqs. 41-42
    // and the steam-generator total-heat-transfer requirement of 160 MW
    // per pair of SGs, i.e. 80 MW per SG with a peak primary-secondary
    // dT of ~ 30 K).
    p.alpha_o_subcooled_W_per_m2K   = 12000.0;
    p.alpha_o_twophase_W_per_m2K    = 16000.0;
    p.alpha_o_superheat_W_per_m2K   = 9000.0;
    p.alpha_i_subcooled_W_per_m2K   = 5000.0;
    p.alpha_i_twophase_W_per_m2K    = 25000.0;    // boiling enhancement
    p.alpha_i_superheat_W_per_m2K   = 4000.0;
    return p;
}

// =============================================================================
// State arithmetic + array conversion
// =============================================================================

HelicalCoilSteamGeneratorState&
HelicalCoilSteamGeneratorState::operator+=(const HelicalCoilSteamGeneratorState& r) noexcept {
    L_1 += r.L_1; L_2 += r.L_2; p_S += r.p_S; h_o += r.h_o;
    T_M1 += r.T_M1; T_M2 += r.T_M2; T_M3 += r.T_M3;
    T_P1 += r.T_P1; T_P2 += r.T_P2; T_P3 += r.T_P3;
    return *this;
}
HelicalCoilSteamGeneratorState&
HelicalCoilSteamGeneratorState::operator-=(const HelicalCoilSteamGeneratorState& r) noexcept {
    L_1 -= r.L_1; L_2 -= r.L_2; p_S -= r.p_S; h_o -= r.h_o;
    T_M1 -= r.T_M1; T_M2 -= r.T_M2; T_M3 -= r.T_M3;
    T_P1 -= r.T_P1; T_P2 -= r.T_P2; T_P3 -= r.T_P3;
    return *this;
}
HelicalCoilSteamGeneratorState&
HelicalCoilSteamGeneratorState::operator*=(double s) noexcept {
    L_1 *= s; L_2 *= s; p_S *= s; h_o *= s;
    T_M1 *= s; T_M2 *= s; T_M3 *= s;
    T_P1 *= s; T_P2 *= s; T_P3 *= s;
    return *this;
}

std::array<double, 10> HelicalCoilSteamGeneratorState::toArray() const noexcept {
    return {L_1, L_2, p_S, h_o, T_M1, T_M2, T_M3, T_P1, T_P2, T_P3};
}
void HelicalCoilSteamGeneratorState::fromArray(const std::array<double, 10>& a) noexcept {
    L_1 = a[0]; L_2 = a[1]; p_S = a[2]; h_o = a[3];
    T_M1 = a[4]; T_M2 = a[5]; T_M3 = a[6];
    T_P1 = a[7]; T_P2 = a[8]; T_P3 = a[9];
}

// =============================================================================
// Local helpers: properties + their numerical derivatives
// =============================================================================
namespace {

constexpr double kEpsP_Pa  = 5.0e3;        // pressure perturbation for finite differences
constexpr double kEpsH_Jkg = 1.0e3;        // enthalpy  perturbation for finite differences

/// Find single-phase liquid temperature corresponding to (h, p).  Bisection
/// over [T_sat - 200, T_sat - eps] for liquid (we stay strictly below the
/// saturation curve so single-phase property functions are well defined).
double T_liquid_from_PH(const props::WaterProperties* w, double p, double h) {
    constexpr double kTsatGap = 0.05;             // K; must match enthalpy_TP guard
    const double Tsat = w->saturationTemperature(p);
    const double h_f  = w->satLiquidEnthalpy_P(p);
    if (h >= h_f) return Tsat - kTsatGap;         // saturated -- back off slightly
    double Tlo = std::max(Tsat - 200.0, 280.0);
    double Thi = Tsat - kTsatGap;
    for (int it = 0; it < 60; ++it) {
        const double Tm = 0.5 * (Tlo + Thi);
        const double hm = w->enthalpy_TP(Tm, p);
        if (hm < h) Tlo = Tm; else Thi = Tm;
        if (Thi - Tlo < 1e-4) break;
    }
    return 0.5 * (Tlo + Thi);
}

/// Find single-phase vapour temperature corresponding to (h, p).  Bisection
/// over [T_sat + eps, T_sat + 400].
double T_vapour_from_PH(const props::WaterProperties* w, double p, double h) {
    constexpr double kTsatGap = 0.05;
    const double Tsat = w->saturationTemperature(p);
    const double h_g  = w->satVapourEnthalpy_P(p);
    if (h <= h_g) return Tsat + kTsatGap;
    double Tlo = Tsat + kTsatGap;
    double Thi = Tsat + 400.0;
    for (int it = 0; it < 60; ++it) {
        const double Tm = 0.5 * (Tlo + Thi);
        const double hm = w->enthalpy_TP(Tm, p);
        if (hm < h) Tlo = Tm; else Thi = Tm;
        if (Thi - Tlo < 1e-4) break;
    }
    return 0.5 * (Tlo + Thi);
}

/// Single-phase density at (P, h) -- inverts h_TP to get T, then returns rho_TP.
/// At the saturation boundary returns the saturated property value.
double rho_PH_liquid(const props::WaterProperties* w, double p, double h) {
    const double h_f = w->satLiquidEnthalpy_P(p);
    if (h >= h_f - 100.0) return w->satLiquidDensity_P(p);  // close enough to saturation
    const double T = T_liquid_from_PH(w, p, h);
    return w->density_TP(T, p);
}
double rho_PH_vapour(const props::WaterProperties* w, double p, double h) {
    const double h_g = w->satVapourEnthalpy_P(p);
    if (h <= h_g + 100.0) return w->satVapourDensity_P(p);
    const double T = T_vapour_from_PH(w, p, h);
    return w->density_TP(T, p);
}

/// (drho/dP)|h for liquid / vapour.
double drho_dP_at_h_liquid(const props::WaterProperties* w, double p, double h) {
    return (rho_PH_liquid(w, p + kEpsP_Pa, h) - rho_PH_liquid(w, p - kEpsP_Pa, h))
           / (2.0 * kEpsP_Pa);
}
double drho_dP_at_h_vapour(const props::WaterProperties* w, double p, double h) {
    return (rho_PH_vapour(w, p + kEpsP_Pa, h) - rho_PH_vapour(w, p - kEpsP_Pa, h))
           / (2.0 * kEpsP_Pa);
}
/// (drho/dh)|P for liquid / vapour.
double drho_dh_at_P_liquid(const props::WaterProperties* w, double p, double h) {
    return (rho_PH_liquid(w, p, h + kEpsH_Jkg) - rho_PH_liquid(w, p, h - kEpsH_Jkg))
           / (2.0 * kEpsH_Jkg);
}
double drho_dh_at_P_vapour(const props::WaterProperties* w, double p, double h) {
    return (rho_PH_vapour(w, p, h + kEpsH_Jkg) - rho_PH_vapour(w, p, h - kEpsH_Jkg))
           / (2.0 * kEpsH_Jkg);
}
/// dh_f/dP, dh_g/dP, drho_f/dP, drho_g/dP from central differences.
double dhf_dP(const props::WaterProperties* w, double p) {
    return (w->satLiquidEnthalpy_P(p + kEpsP_Pa) - w->satLiquidEnthalpy_P(p - kEpsP_Pa))
           / (2.0 * kEpsP_Pa);
}
double dhg_dP(const props::WaterProperties* w, double p) {
    return (w->satVapourEnthalpy_P(p + kEpsP_Pa) - w->satVapourEnthalpy_P(p - kEpsP_Pa))
           / (2.0 * kEpsP_Pa);
}
double drhof_dP(const props::WaterProperties* w, double p) {
    return (w->satLiquidDensity_P(p + kEpsP_Pa) - w->satLiquidDensity_P(p - kEpsP_Pa))
           / (2.0 * kEpsP_Pa);
}
double drhog_dP(const props::WaterProperties* w, double p) {
    return (w->satVapourDensity_P(p + kEpsP_Pa) - w->satVapourDensity_P(p - kEpsP_Pa))
           / (2.0 * kEpsP_Pa);
}

/// Mean void fraction along a uniformly-heated tube where quality varies
/// linearly from 0 (saturated liquid, two-phase inlet) to 1 (saturated
/// vapour, two-phase exit).  Closed-form integral of homogeneous
/// alpha(x) = 1/(1 + (1-x)/x * gamma) over x in [0, 1]:
///
///     mean alpha = 1/(1 - gamma) - gamma * ln(1/gamma) / (1 - gamma)^2
///
/// where gamma = rho_g / rho_l < 1.  This is the Jensen-Tummescheit (2002)
/// "uniform heat flux" closed form.  At gamma -> 0 the limit is
/// alpha -> 1; at gamma -> 1 the integrand becomes 1 and alpha -> 0.5.
double meanVoidFraction(double rho_g, double rho_l) {
    if (!(rho_l > 0.0)) return 0.5;
    const double gamma = rho_g / rho_l;
    if (gamma >= 0.999) return 0.5;
    if (gamma <= 1.0e-6) return 1.0 - gamma;     // avoid log(0)
    const double inv_one_minus_gamma = 1.0 / (1.0 - gamma);
    return inv_one_minus_gamma
           - gamma * std::log(1.0 / gamma) * inv_one_minus_gamma * inv_one_minus_gamma;
}

/// Solve A x = b for a 10 x 10 matrix in row-major order, in place.
/// Uses partial pivoting Gaussian elimination.  Returns false on
/// near-singularity.
bool solve10(std::array<std::array<double, 10>, 10>& A,
             std::array<double, 10>& b) {
    constexpr int N = 10;
    for (int k = 0; k < N; ++k) {
        // Pivot.
        int piv = k;
        double pmag = std::abs(A[k][k]);
        for (int i = k + 1; i < N; ++i) {
            if (std::abs(A[i][k]) > pmag) {
                pmag = std::abs(A[i][k]);
                piv = i;
            }
        }
        if (pmag < 1.0e-30) return false;
        if (piv != k) {
            std::swap(A[piv], A[k]);
            std::swap(b[piv], b[k]);
        }
        // Eliminate.
        const double inv_pv = 1.0 / A[k][k];
        for (int i = k + 1; i < N; ++i) {
            const double factor = A[i][k] * inv_pv;
            if (factor == 0.0) continue;
            for (int j = k; j < N; ++j) A[i][j] -= factor * A[k][j];
            b[i] -= factor * b[k];
        }
    }
    // Back-substitute.
    for (int i = N - 1; i >= 0; --i) {
        double s = b[i];
        for (int j = i + 1; j < N; ++j) s -= A[i][j] * b[j];
        b[i] = s / A[i][i];
    }
    return true;
}

}  // namespace

// =============================================================================
// Constructor
// =============================================================================

HelicalCoilSteamGenerator::HelicalCoilSteamGenerator(
        HelicalCoilSteamGeneratorParameters params,
        const props::WaterProperties* props)
        : p_(std::move(params)), props_(props) {
    p_.validate();
    if (props_ == nullptr) {
        throw std::invalid_argument(
                "HelicalCoilSteamGenerator: water-properties service must not be null");
    }
}

// =============================================================================
// Steady-state initialisation (Arda 2015 Table 4 target)
// =============================================================================

void HelicalCoilSteamGenerator::initialiseSteadyState(double T_pi_in_K,
                                                     double primary_mass_flow_kg_s,
                                                     double feedwater_flow_kg_s,
                                                     double feedwater_temperature_K) {
    if (!(T_pi_in_K > 0.0))               throw std::invalid_argument("T_Pi must be > 0");
    if (!(primary_mass_flow_kg_s > 0.0))  throw std::invalid_argument("m_dot_P must be > 0");
    if (!(feedwater_flow_kg_s > 0.0))     throw std::invalid_argument("m_dot_S,i must be > 0");
    if (!(feedwater_temperature_K > 273)) throw std::invalid_argument("T_S,i must be > 273 K");

    state_ = {};
    state_.t_s = 0.0;

    // Target steady state from Arda 2015 Table 4.
    state_.p_S = 3.1e6;                                  // 3.1 MPa
    const double Tsat = props_->saturationTemperature(state_.p_S);
    const double h_f  = props_->satLiquidEnthalpy_P(state_.p_S);
    const double h_g  = props_->satVapourEnthalpy_P(state_.p_S);
    const double h_fg = h_g - h_f;

    // Region lengths from Table 4 (proportional to total length L_T = 22.25 m).
    state_.L_1 = 2.90;
    state_.L_2 = 17.60;
    const double L_3 = p_.tube_total_length_m - state_.L_1 - state_.L_2;
    if (!(L_3 > 0.0)) {
        throw std::runtime_error(
                "HelicalCoilSG: parameter tube_total_length_m too short for the "
                "Arda Table 4 region lengths (L_T must exceed L_1 + L_2 = 20.50 m)");
    }

    // Total heat to remove = m_dot_S,i * (h_o - h_i).
    const double h_i = props_->enthalpy_TP(feedwater_temperature_K, state_.p_S);
    const double Q_total = primary_mass_flow_kg_s * p_.primary_cp_J_per_kgK
                           * (T_pi_in_K - 246.0 - 273.15);   // T_P1 ~ 246 deg C from Table 4
    // From Table 4 the steam outlet temperature is 264 deg C.  Compute h_o
    // from that target.
    const double T_steam_K = 264.0 + 273.15;
    const double h_o_target = props_->enthalpy_TP(T_steam_K, state_.p_S);
    state_.h_o = h_o_target;
    (void)Q_total; (void)h_fg; (void)h_i;     // documented but not needed below

    // Distribute the total dT across the three primary regions in proportion
    // to (length * heat-transfer coefficient).  Denote the local Q per
    // region as Q_k = m_dot_P * cp_P * (T_P_(k+1) - T_P_k) (with T_P,4 =
    // T_pi_in by Arda's numbering: region 1 is at the tube cold end,
    // region 3 is where the primary enters from the top -- see Fig. 3).
    const double Wcp = primary_mass_flow_kg_s * p_.primary_cp_J_per_kgK;
    const double Q_design = Wcp * (T_pi_in_K - (246.0 + 273.15));  // ~ 80 MW
    const double UA1 = state_.L_1 * (p_.tube_outer_diameter_m * 3.14159265358979 * p_.alpha_o_subcooled_W_per_m2K);
    const double UA2 = state_.L_2 * (p_.tube_outer_diameter_m * 3.14159265358979 * p_.alpha_o_twophase_W_per_m2K);
    const double UA3 = L_3        * (p_.tube_outer_diameter_m * 3.14159265358979 * p_.alpha_o_superheat_W_per_m2K);
    const double UA_sum = UA1 + UA2 + UA3;
    const double Q1 = Q_design * UA1 / UA_sum;
    const double Q2 = Q_design * UA2 / UA_sum;
    const double Q3 = Q_design * UA3 / UA_sum;
    // Primary temperatures (region 1 at bottom, region 3 at top).
    state_.T_P1 = 246.0 + 273.15;             // T_P1 from Table 4
    state_.T_P2 = state_.T_P1 + Q1 / Wcp;
    state_.T_P3 = state_.T_P2 + Q2 / Wcp;
    // Sanity: T_P3 + Q3/Wcp should equal T_pi_in.
    (void)Q3;

    // Metal temperatures from Q = (d_o * pi * L * alpha_o) * (T_P - T_M)
    // and Q = (d_i * pi * L * alpha_i) * (T_M - T_S).  Solve for T_M.
    auto metal_T = [&](double T_P, double T_S, double L, double alpha_o, double alpha_i) {
        const double UAo = p_.tube_outer_diameter_m * 3.14159265358979 * L * alpha_o;
        const double UAi = p_.tube_inner_diameter_m * 3.14159265358979 * L * alpha_i;
        return (UAo * T_P + UAi * T_S) / (UAo + UAi);
    };
    const double T_S1_avg = 0.5 * (feedwater_temperature_K + Tsat);
    const double T_S3_avg = 0.5 * (Tsat + T_steam_K);
    state_.T_M1 = metal_T(state_.T_P1, T_S1_avg, state_.L_1,
                           p_.alpha_o_subcooled_W_per_m2K, p_.alpha_i_subcooled_W_per_m2K);
    state_.T_M2 = metal_T(state_.T_P2, Tsat, state_.L_2,
                           p_.alpha_o_twophase_W_per_m2K, p_.alpha_i_twophase_W_per_m2K);
    state_.T_M3 = metal_T(state_.T_P3, T_S3_avg, L_3,
                           p_.alpha_o_superheat_W_per_m2K, p_.alpha_i_superheat_W_per_m2K);

    // Inputs: feed steam outlet flow = feedwater flow (steady-state mass balance).
    inputs_.primary_inlet_temperature_K = T_pi_in_K;
    inputs_.primary_mass_flow_kg_s      = primary_mass_flow_kg_s;
    inputs_.feedwater_mass_flow_kg_s    = feedwater_flow_kg_s;
    inputs_.feedwater_enthalpy_J_kg     = h_i;
    inputs_.steam_outlet_mass_flow_kg_s = feedwater_flow_kg_s;
}

// =============================================================================
// Outputs
// =============================================================================

double HelicalCoilSteamGenerator::steamOutletTemperatureK() const {
    return T_vapour_from_PH(props_, state_.p_S, state_.h_o);
}

double HelicalCoilSteamGenerator::primaryHeatLoad_W() const noexcept {
    return inputs_.primary_mass_flow_kg_s * p_.primary_cp_J_per_kgK
           * (inputs_.primary_inlet_temperature_K - state_.T_P1);
}

// =============================================================================
// Derivative evaluation (the heart of the model).  Forms D and f, solves
// D dx = f via partial-pivoting Gaussian elimination.
// =============================================================================

HelicalCoilSteamGeneratorState
HelicalCoilSteamGenerator::evaluateDerivative(
        const HelicalCoilSteamGeneratorState& s) const {
    constexpr double pi = 3.14159265358979323846;

    // ----- Region lengths -------------------------------------------------
    const double L_T = p_.tube_total_length_m;
    const double L_3 = L_T - s.L_1 - s.L_2;
    if (!(s.L_1 > 0.0 && s.L_2 > 0.0 && L_3 > 0.0)) {
        throw std::runtime_error(
                "HelicalCoilSG::evaluateDerivative: a region length has crossed "
                "zero (model topology no longer valid).");
    }

    // ----- Saturation properties at s.p_S ---------------------------------
    const double Tsat   = props_->saturationTemperature(s.p_S);
    const double h_f    = props_->satLiquidEnthalpy_P(s.p_S);
    const double h_g    = props_->satVapourEnthalpy_P(s.p_S);
    const double rho_f  = props_->satLiquidDensity_P(s.p_S);
    const double rho_g  = props_->satVapourDensity_P(s.p_S);

    const double dhf  = dhf_dP(props_, s.p_S);
    const double dhg  = dhg_dP(props_, s.p_S);
    const double drhof= drhof_dP(props_, s.p_S);
    const double drhog= drhog_dP(props_, s.p_S);

    // ----- Region-average enthalpies and densities ------------------------
    const double h_i = inputs_.feedwater_enthalpy_J_kg;
    const double h_1 = 0.5 * (h_i + h_f);                  // subcooled mean
    const double h_3 = 0.5 * (h_g + s.h_o);                // superheated mean

    const double rho_1 = rho_PH_liquid(props_, s.p_S, h_1);
    const double rho_3 = rho_PH_vapour(props_, s.p_S, h_3);

    const double drho1_dP = drho_dP_at_h_liquid(props_, s.p_S, h_1);
    const double drho1_dh = drho_dh_at_P_liquid(props_, s.p_S, h_1);
    const double drho3_dP = drho_dP_at_h_vapour(props_, s.p_S, h_3);
    const double drho3_dh = drho_dh_at_P_vapour(props_, s.p_S, h_3);

    // ----- Mean void fraction in the two-phase region ---------------------
    const double alpha_void = meanVoidFraction(rho_g, rho_f);
    const double one_minus_a = 1.0 - alpha_void;

    // ----- Secondary-side average temperatures (for heat transfer) --------
    const double T_S1 = T_liquid_from_PH(props_, s.p_S, h_1);
    const double T_S3 = T_vapour_from_PH(props_, s.p_S, h_3);
    // T_S2 == Tsat (already in scope as `Tsat`).

    // ----- Heat-transfer coefficients * tube perimeter --------------------
    const double UAi1 = pi * p_.tube_inner_diameter_m * p_.tube_count
                          * p_.alpha_i_subcooled_W_per_m2K;       // per metre
    const double UAi2 = pi * p_.tube_inner_diameter_m * p_.tube_count
                          * p_.alpha_i_twophase_W_per_m2K;
    const double UAi3 = pi * p_.tube_inner_diameter_m * p_.tube_count
                          * p_.alpha_i_superheat_W_per_m2K;
    const double UAo1 = pi * p_.tube_outer_diameter_m * p_.tube_count
                          * p_.alpha_o_subcooled_W_per_m2K;
    const double UAo2 = pi * p_.tube_outer_diameter_m * p_.tube_count
                          * p_.alpha_o_twophase_W_per_m2K;
    const double UAo3 = pi * p_.tube_outer_diameter_m * p_.tube_count
                          * p_.alpha_o_superheat_W_per_m2K;

    // Boundary inputs.
    const double mP    = inputs_.primary_mass_flow_kg_s;
    const double cpP   = p_.primary_cp_J_per_kgK;
    const double mSi   = inputs_.feedwater_mass_flow_kg_s;
    const double mSo   = inputs_.steam_outlet_mass_flow_kg_s;
    const double T_Pi  = inputs_.primary_inlet_temperature_K;
    const double AS    = p_.secondary_flow_area_m2;
    const double AP    = p_.primary_flow_area_m2;
    const double AM    = p_.tube_metal_area_m2;
    const double rhoP  = p_.primary_density_kg_m3;
    const double rhoM  = p_.tube_metal_density_kg_m3;
    const double cpM   = p_.tube_metal_cp_J_per_kgK;

    // ----- Build D (10x10) and f (10) -- Arda 2015 eq. 44, Table 3 -------
    std::array<std::array<double, 10>, 10> D{};
    std::array<double, 10>                  f{};
    for (auto& row : D) row.fill(0.0);
    f.fill(0.0);

    // Convenience: indices into x = [L1, L2, pS, ho, TM1, TM2, TM3, TP1, TP2, TP3]
    enum { iL1=0, iL2=1, ipS=2, iho=3, iTM1=4, iTM2=5, iTM3=6, iTP1=7, iTP2=8, iTP3=9 };

    // ---- Row 1: subcooled energy balance ----
    //   d_{1,1} = AS * rho_1 * (h_1 - h_f)
    //   d_{1,3} = AS * L_1 * { [drho1/dP|h1 + 0.5 * drho1/dh1|p * dhf/dP]*(h1 - h_f)
    //                          + 0.5 * dhf/dP - 1 }
    //   f_1     = mSi * (h_i - h_f) + di * alpha_i1 * L1 * (T_M1 - T_S1)
    D[0][iL1] = AS * rho_1 * (h_1 - h_f);
    D[0][ipS] = AS * s.L_1 * ((drho1_dP + 0.5 * drho1_dh * dhf) * (h_1 - h_f)
                              + 0.5 * dhf - 1.0);
    f[0] = mSi * (h_i - h_f)
         + UAi1 / p_.tube_count * p_.tube_count * s.L_1 * (s.T_M1 - T_S1);
    // The ".UAi1/tube_count*tube_count" is identity but kept so that you
    // can later switch to per-tube UAi if you'd rather; simplifies to:
    f[0] = mSi * (h_i - h_f) + UAi1 * s.L_1 * (s.T_M1 - T_S1);

    // ---- Row 2: two-phase energy balance ----
    //   d_{2,1} = AS * (rho_1*h_f - rho_3*h_g)
    //   d_{2,2} = AS * { (1-alpha)*(rho_f*h_f - rho_g*h_g) + h_g*(rho_g - rho_3) }
    //   d_{2,3} = AS * [ h_f*L_1*{drho1/dP|h1 + 0.5*drho1/dh1|p*dhf/dP}
    //                  + L_2*{ alpha*d(rho_g h_g)/dP + (1-alpha)*d(rho_f h_f)/dP - 1 }
    //                  + h_g*L_3*{drho3/dP|h3 + 0.5*drho3/dh3|p*dhg/dP} ]
    //   d_{2,4} = AS * h_g * L_3 / 2 * drho3/dh3|p
    //   f_2     = mSi*h_f - mSo*h_g + di*alpha_i2*L2*(TM2 - Tsat)
    const double d_rhog_hg_dP = drhog * h_g + rho_g * dhg;
    const double d_rhof_hf_dP = drhof * h_f + rho_f * dhf;
    D[1][iL1] = AS * (rho_1 * h_f - rho_3 * h_g);
    D[1][iL2] = AS * (one_minus_a * (rho_f * h_f - rho_g * h_g)
                       + h_g * (rho_g - rho_3));
    D[1][ipS] = AS * (
            h_f * s.L_1 * (drho1_dP + 0.5 * drho1_dh * dhf)
          + s.L_2 * (alpha_void * d_rhog_hg_dP + one_minus_a * d_rhof_hf_dP - 1.0)
          + h_g * L_3 * (drho3_dP + 0.5 * drho3_dh * dhg));
    D[1][iho] = AS * h_g * L_3 * 0.5 * drho3_dh;
    f[1] = mSi * h_f - mSo * h_g + UAi2 * s.L_2 * (s.T_M2 - Tsat);

    // ---- Row 3: superheated energy balance ----
    //   d_{3,1} = AS * rho_3 * (h_g - h_3)
    //   d_{3,2} = AS * rho_3 * (h_g - h_3)
    //   d_{3,3} = AS * L_3 * { [drho3/dP|h3 + 0.5*drho3/dh3|p*dhg/dP]*(h_3 - h_g)
    //                         + (rho_3/2)*dhg/dP - 1 }
    //   d_{3,4} = AS * L_3 * { 0.5 * drho3/dh3|p * (h_3 - h_g) + rho_3/2 }
    //   f_3     = mSo*(h_g - h_o) + di*alpha_i3*L3*(TM3 - T_S3)
    D[2][iL1] = AS * rho_3 * (h_g - h_3);
    D[2][iL2] = AS * rho_3 * (h_g - h_3);
    D[2][ipS] = AS * L_3 * ((drho3_dP + 0.5 * drho3_dh * dhg) * (h_3 - h_g)
                             + 0.5 * rho_3 * dhg - 1.0);
    D[2][iho] = AS * L_3 * (0.5 * drho3_dh * (h_3 - h_g) + 0.5 * rho_3);
    f[2] = mSo * (h_g - s.h_o) + UAi3 * L_3 * (s.T_M3 - T_S3);

    // ---- Row 4: total mass balance ----
    //   d_{4,1} = AS * (rho_1 - rho_3)
    //   d_{4,2} = AS * { (1-alpha)*(rho_f - rho_g) + (rho_g - rho_3) }
    //   d_{4,3} = AS * [ L_1*{drho1/dP|h1 + 0.5*drho1/dh1|p*dhf/dP}
    //                   + L_2*{ alpha*drho_g/dP + (1-alpha)*drho_f/dP }
    //                   + L_3*{drho3/dP|h3 + 0.5*drho3/dh3|p*dhg/dP} ]
    //   d_{4,4} = AS * L_3 / 2 * drho3/dh3|p
    //   f_4     = mSi - mSo
    D[3][iL1] = AS * (rho_1 - rho_3);
    D[3][iL2] = AS * (one_minus_a * (rho_f - rho_g) + (rho_g - rho_3));
    D[3][ipS] = AS * (
            s.L_1 * (drho1_dP + 0.5 * drho1_dh * dhf)
          + s.L_2 * (alpha_void * drhog + one_minus_a * drhof)
          + L_3 * (drho3_dP + 0.5 * drho3_dh * dhg));
    D[3][iho] = AS * L_3 * 0.5 * drho3_dh;
    f[3] = mSi - mSo;

    // ---- Rows 5-7: tube-metal energy balances (eq. 35-37) ----
    //   d_{5,1} = AM*rho_M*cp_M*(T_M1 - T_M2)
    //   d_{5,5} = AM*rho_M*cp_M*L_1
    //   f_5     = L_1*[ d_o*alpha_o1*(T_P1 - T_M1) - d_i*alpha_i1*(T_M1 - T_S1) ]
    const double Mtube = AM * rhoM * cpM;
    D[4][iL1] = Mtube * (s.T_M1 - s.T_M2);
    D[4][iTM1] = Mtube * s.L_1;
    f[4] = s.L_1 * (UAo1 * (s.T_P1 - s.T_M1) / p_.tube_count * 1.0   // outer per metre
                  - UAi1 * (s.T_M1 - T_S1) / p_.tube_count * 1.0);   // inner per metre
    // simpler: UAo1 already includes tube_count -- so:
    f[4] = s.L_1 * (UAo1 * (s.T_P1 - s.T_M1) - UAi1 * (s.T_M1 - T_S1));

    //   d_{6,6} = AM*rho_M*cp_M*L_2
    //   f_6     = L_2 * [ d_o*alpha_o2*(T_P2 - T_M2) - d_i*alpha_i2*(T_M2 - Tsat) ]
    D[5][iTM2] = Mtube * s.L_2;
    f[5] = s.L_2 * (UAo2 * (s.T_P2 - s.T_M2) - UAi2 * (s.T_M2 - Tsat));

    //   d_{7,1} = AM*rho_M*cp_M*(T_M2 - T_M3)
    //   d_{7,2} = AM*rho_M*cp_M*(T_M2 - T_M3)
    //   d_{7,7} = AM*rho_M*cp_M*L_3
    //   f_7     = L_3 * [ d_o*alpha_o3*(T_P3 - T_M3) - d_i*alpha_i3*(T_M3 - T_S3) ]
    D[6][iL1] = Mtube * (s.T_M2 - s.T_M3);
    D[6][iL2] = Mtube * (s.T_M2 - s.T_M3);
    D[6][iTM3] = Mtube * L_3;
    f[6] = L_3 * (UAo3 * (s.T_P3 - s.T_M3) - UAi3 * (s.T_M3 - T_S3));

    // ---- Rows 8-10: primary-side energy balances (eq. 38-40) ----
    //   d_{8,1} = AP*rho_P*cp_P*(T_P1 - T_P2)
    //   d_{8,8} = AP*rho_P*cp_P*L_1
    //   f_8     = m_dot_P*cp_P*(T_P2 - T_P1) - d_o*alpha_o1*L_1*(T_P1 - T_M1)
    const double Mprim = AP * rhoP * cpP;
    D[7][iL1] = Mprim * (s.T_P1 - s.T_P2);
    D[7][iTP1] = Mprim * s.L_1;
    f[7] = mP * cpP * (s.T_P2 - s.T_P1) - UAo1 * s.L_1 * (s.T_P1 - s.T_M1);

    //   d_{9,9}  = AP*rho_P*cp_P*L_2
    //   f_9      = m_dot_P*cp_P*(T_P3 - T_P2) - d_o*alpha_o2*L_2*(T_P2 - T_M2)
    D[8][iTP2] = Mprim * s.L_2;
    f[8] = mP * cpP * (s.T_P3 - s.T_P2) - UAo2 * s.L_2 * (s.T_P2 - s.T_M2);

    //   d_{10,1} = AP*rho_P*cp_P*(T_P2 - T_P3)
    //   d_{10,2} = AP*rho_P*cp_P*(T_P2 - T_P3)
    //   d_{10,10}= AP*rho_P*cp_P*L_3
    //   f_10     = m_dot_P*cp_P*(T_Pi - T_P3) - d_o*alpha_o3*L_3*(T_P3 - T_M3)
    D[9][iL1] = Mprim * (s.T_P2 - s.T_P3);
    D[9][iL2] = Mprim * (s.T_P2 - s.T_P3);
    D[9][iTP3] = Mprim * L_3;
    f[9] = mP * cpP * (T_Pi - s.T_P3) - UAo3 * L_3 * (s.T_P3 - s.T_M3);

    // ----- Solve D dx = f -------------------------------------------------
    if (!solve10(D, f)) {
        throw std::runtime_error(
                "HelicalCoilSG::evaluateDerivative: D matrix is near-singular.  "
                "Region lengths or properties are inconsistent.");
    }

    HelicalCoilSteamGeneratorState dy{};
    dy.L_1 = f[iL1]; dy.L_2 = f[iL2]; dy.p_S = f[ipS]; dy.h_o = f[iho];
    dy.T_M1 = f[iTM1]; dy.T_M2 = f[iTM2]; dy.T_M3 = f[iTM3];
    dy.T_P1 = f[iTP1]; dy.T_P2 = f[iTP2]; dy.T_P3 = f[iTP3];
    return dy;
}

void HelicalCoilSteamGenerator::timeStep(double dt) {
    if (!(dt > 0.0)) {
        throw std::invalid_argument("HelicalCoilSG::timeStep: dt must be positive");
    }
    auto deriv = [this](double, const HelicalCoilSteamGeneratorState& s) {
        return this->evaluateDerivative(s);
    };
    state_ = astara::core::rk4Step(state_.t_s, state_, dt, deriv);
    state_.t_s += dt;

    // Sanity.
    auto bad = [](double v){ return !std::isfinite(v); };
    if (bad(state_.p_S) || bad(state_.L_1) || bad(state_.L_2)
        || bad(state_.h_o) || bad(state_.T_P1)) {
        throw std::runtime_error("HelicalCoilSG::timeStep produced non-finite state");
    }
    const double L_3 = p_.tube_total_length_m - state_.L_1 - state_.L_2;
    if (state_.L_1 <= 0.0 || state_.L_2 <= 0.0 || L_3 <= 0.0) {
        throw std::runtime_error(
                "HelicalCoilSG::timeStep: a region length has crossed zero (model "
                "topology no longer valid).");
    }
}

}  // namespace astara::sg
