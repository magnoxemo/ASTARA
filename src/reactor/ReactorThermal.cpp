/**
 * @file   ReactorThermal.cpp
 * @brief  Implementation of the lumped-parameter core thermal model.
 */

#include "astara/reactor/ReactorThermal.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>

namespace astara::reactor {

void ReactorThermalParameters::validate() const {
    auto check_pos = [](double v, const char* name) {
        if (!(v > 0.0)) {
            std::ostringstream os;
            os << "ReactorThermalParameters: " << name << " must be positive, got " << v;
            throw std::invalid_argument(os.str());
        }
    };
    if (num_fuel_nodes < 1) {
        throw std::invalid_argument("ReactorThermalParameters: num_fuel_nodes >= 1 required");
    }
    if (num_moderator_nodes != 2 * num_fuel_nodes) {
        std::ostringstream os;
        os << "ReactorThermalParameters: num_moderator_nodes (" << num_moderator_nodes
           << ") must be exactly 2 * num_fuel_nodes (" << 2*num_fuel_nodes << ") "
           << "for Mann nodalisation.";
        throw std::invalid_argument(os.str());
    }
    check_pos(fuel_mass_total_kg,      "fuel_mass_total_kg");
    check_pos(fuel_cp_J_per_kgK,       "fuel_cp_J_per_kgK");
    check_pos(moderator_mass_total_kg, "moderator_mass_total_kg");
    check_pos(moderator_cp_J_per_kgK,  "moderator_cp_J_per_kgK");
    check_pos(mass_flow_rate_kg_s,     "mass_flow_rate_kg_s");
    check_pos(overall_h_W_per_m2K,     "overall_h_W_per_m2K");
    check_pos(heat_transfer_area_m2,   "heat_transfer_area_m2");
    check_pos(lower_plenum_mass_kg,    "lower_plenum_mass_kg");
    check_pos(upper_plenum_mass_kg,    "upper_plenum_mass_kg");
    check_pos(hot_leg_mass_kg,         "hot_leg_mass_kg");
    check_pos(cold_leg_mass_kg,        "cold_leg_mass_kg");
    if (!(fission_power_in_fuel > 0.0 && fission_power_in_fuel <= 1.0)) {
        throw std::invalid_argument(
                "ReactorThermalParameters: fission_power_in_fuel must lie in (0, 1].");
    }
}

// =============================================================================
// ReactorThermalState arithmetic
// =============================================================================

std::size_t ReactorThermalState::size() const noexcept {
    return T_fuel.size() + T_moderator.size() + 4;
}

double ReactorThermalState::averageFuelTemperature() const noexcept {
    if (T_fuel.empty()) return 0.0;
    double s = 0.0;
    for (double v : T_fuel) s += v;
    return s / static_cast<double>(T_fuel.size());
}

double ReactorThermalState::averageModeratorTemperature() const noexcept {
    if (T_moderator.empty()) return 0.0;
    double s = 0.0;
    for (double v : T_moderator) s += v;
    return s / static_cast<double>(T_moderator.size());
}

ReactorThermalState& ReactorThermalState::operator+=(const ReactorThermalState& rhs) {
    if (T_fuel.size() != rhs.T_fuel.size() || T_moderator.size() != rhs.T_moderator.size()) {
        throw std::invalid_argument("ReactorThermalState size mismatch in +=");
    }
    for (std::size_t i = 0; i < T_fuel.size(); ++i)      T_fuel[i]      += rhs.T_fuel[i];
    for (std::size_t i = 0; i < T_moderator.size(); ++i) T_moderator[i] += rhs.T_moderator[i];
    T_cold_leg     += rhs.T_cold_leg;
    T_lower_plenum += rhs.T_lower_plenum;
    T_upper_plenum += rhs.T_upper_plenum;
    T_hot_leg      += rhs.T_hot_leg;
    return *this;
}

ReactorThermalState& ReactorThermalState::operator-=(const ReactorThermalState& rhs) {
    if (T_fuel.size() != rhs.T_fuel.size() || T_moderator.size() != rhs.T_moderator.size()) {
        throw std::invalid_argument("ReactorThermalState size mismatch in -=");
    }
    for (std::size_t i = 0; i < T_fuel.size(); ++i)      T_fuel[i]      -= rhs.T_fuel[i];
    for (std::size_t i = 0; i < T_moderator.size(); ++i) T_moderator[i] -= rhs.T_moderator[i];
    T_cold_leg     -= rhs.T_cold_leg;
    T_lower_plenum -= rhs.T_lower_plenum;
    T_upper_plenum -= rhs.T_upper_plenum;
    T_hot_leg      -= rhs.T_hot_leg;
    return *this;
}

ReactorThermalState& ReactorThermalState::operator*=(double s) noexcept {
    for (double& v : T_fuel)      v *= s;
    for (double& v : T_moderator) v *= s;
    T_cold_leg     *= s;
    T_lower_plenum *= s;
    T_upper_plenum *= s;
    T_hot_leg      *= s;
    return *this;
}

// =============================================================================
// Steady state
// =============================================================================

ReactorThermalState steadyStateThermal(const ReactorThermalParameters& p,
                                       double P_thermal,
                                       double T_inlet_K) {
    p.validate();
    if (!(P_thermal >= 0.0)) {
        throw std::invalid_argument("steadyStateThermal: P_thermal must be non-negative");
    }
    if (!(T_inlet_K > 0.0)) {
        throw std::invalid_argument("steadyStateThermal: T_inlet_K must be positive (K)");
    }

    const std::size_t Nf = p.num_fuel_nodes;
    const std::size_t Nm = p.num_moderator_nodes;
    const double      mdot   = p.mass_flow_rate_kg_s;
    const double      cp_c   = p.moderator_cp_J_per_kgK;
    const double      h      = p.overall_h_W_per_m2K;
    const double      A_per  = p.heat_transfer_area_m2 / static_cast<double>(Nf);
    const double      P_per_fuel = p.fission_power_in_fuel * P_thermal / static_cast<double>(Nf);

    ReactorThermalState s;
    s.T_fuel.assign(Nf, 0.0);
    s.T_moderator.assign(Nm, 0.0);

    // Lower plenum sees the cold-leg-in coolant.
    s.T_lower_plenum = T_inlet_K;

    // Closed-form steady-state distribution within each fuel node.
    //
    // Derivation: setting dT_m,2k/dt = 0 and dT_m,2k+1/dt = 0 (using the
    // derivative function below) and subtracting yields
    //
    //     T_m,2k+1 - T_m,2k  =  S / (X + Y)
    //
    // where
    //     S = P / (2 N_f)            (total energy added across fuel node k)
    //     X = h * A_per / 4          (per-half fuel-side heat-transfer conductance)
    //     Y = mdot * cp_c            (advective conductance per moderator node;
    //                                 the full primary mass flow passes through
    //                                 every node in series)
    //
    // and substituting into either node-equation gives
    //
    //     T_m,2k - T_up_k    =  S * (2X + Y) / (Y * (X + Y))
    //
    // The total rise across fuel node k is then 2 S / Y = (P/N_f) / (mdot*cp_c),
    // and the total rise from inlet to outlet is P/(mdot*cp_c) as required by
    // global energy balance.  See Naghedolfeizi (1990) Section 3.1.1; the same
    // closed form is implicit in the thesis but never written out explicitly.
    const double X = 0.25 * h * A_per;
    const double Y = mdot * cp_c;
    const double S = P_thermal / (2.0 * static_cast<double>(Nf));
    const double dT_within = S / (X + Y);
    const double dT_first  = S * (2.0 * X + Y) / (Y * (X + Y));

    double T_up = T_inlet_K;
    for (std::size_t k = 0; k < Nf; ++k) {
        const double Tm_lo = T_up + dT_first;
        const double Tm_hi = Tm_lo + dT_within;
        s.T_moderator[2 * k]     = Tm_lo;
        s.T_moderator[2 * k + 1] = Tm_hi;
        const double Tm_avg = 0.5 * (Tm_lo + Tm_hi);
        s.T_fuel[k] = Tm_avg + P_per_fuel / (h * A_per);
        T_up = Tm_hi;
    }

    // Upper plenum sees the last moderator node; legs follow.  At steady
    // state these reduce to identity transformations because they are pure
    // first-order lags driven by their inlets.
    s.T_upper_plenum = s.T_moderator.back();
    s.T_hot_leg      = s.T_upper_plenum;
    s.T_cold_leg     = T_inlet_K;  // the user-specified inlet -- this is what
                                    // the upstream SG/pump produces in steady state

    return s;
}

// =============================================================================
// Derivatives
// =============================================================================

ReactorThermalState reactorThermalDerivative(const ReactorThermalState& s,
                                             const ReactorThermalParameters& p,
                                             double P_thermal_W,
                                             double T_inlet_K) {
    if (s.T_fuel.size() != p.num_fuel_nodes ||
        s.T_moderator.size() != p.num_moderator_nodes) {
        throw std::invalid_argument(
                "reactorThermalDerivative: state node counts disagree with parameters");
    }

    const std::size_t Nf = p.num_fuel_nodes;
    const std::size_t Nm = p.num_moderator_nodes;
    const double      mdot     = p.mass_flow_rate_kg_s;
    const double      cp_c     = p.moderator_cp_J_per_kgK;
    const double      h        = p.overall_h_W_per_m2K;
    const double      A_per    = p.heat_transfer_area_m2 / static_cast<double>(Nf);
    const double      MF_per_node = p.fuel_mass_total_kg      / static_cast<double>(Nf);
    const double      MC_per_node = p.moderator_mass_total_kg / static_cast<double>(Nm);
    const double      cp_F     = p.fuel_cp_J_per_kgK;
    const double      P_per_fuel = p.fission_power_in_fuel * P_thermal_W / static_cast<double>(Nf);
    const double      P_per_mod  = (1.0 - p.fission_power_in_fuel) * P_thermal_W / static_cast<double>(Nm);

    ReactorThermalState dy;
    dy.T_fuel.assign(Nf, 0.0);
    dy.T_moderator.assign(Nm, 0.0);

    // Fuel nodes (parallel-safe; each only depends on its own and its two
    // moderator neighbours).  OpenMP earns its keep here only when Nf is
    // large; for the thesis Nf=3 we leave it serial.
    for (std::size_t k = 0; k < Nf; ++k) {
        const double Tm_avg = 0.5 * (s.T_moderator[2*k] + s.T_moderator[2*k + 1]);
        const double q_to_mod = h * A_per * (s.T_fuel[k] - Tm_avg);
        dy.T_fuel[k] = (P_per_fuel - q_to_mod) / (MF_per_node * cp_F);
    }

    // Moderator nodes.  Node j is bridged by fuel node k = j/2.
    // Upstream advection comes from node j-1 (or T_lower_plenum for j=0).
    // The full primary mass flow `mdot` passes through every moderator node
    // in series, so the advection rate is `mdot * cp * (T_up - T_node)`
    // *per node*, not divided by Nm.
    for (std::size_t j = 0; j < Nm; ++j) {
        const std::size_t k    = j / 2;
        const double T_up      = (j == 0) ? s.T_lower_plenum : s.T_moderator[j - 1];
        const double q_from_f  = 0.5 * h * A_per * (s.T_fuel[k] - s.T_moderator[j]);
        const double q_advect  = mdot * cp_c * (T_up - s.T_moderator[j]);
        dy.T_moderator[j] = (P_per_mod + q_from_f + q_advect)
                            / (MC_per_node * cp_c);
    }

    // Plena and legs are pure first-order lags driven by their inlets.
    // Time constants: tau = M / mdot.  The thesis writes them as M*cp/(mdot*cp)
    // = M/mdot with cp cancelling.
    auto firstOrder = [&](double T_in, double T_self, double M) {
        const double tau = M / mdot;
        return (T_in - T_self) / tau;
    };
    dy.T_lower_plenum = firstOrder(s.T_cold_leg,         s.T_lower_plenum, p.lower_plenum_mass_kg);
    dy.T_upper_plenum = firstOrder(s.T_moderator.back(), s.T_upper_plenum, p.upper_plenum_mass_kg);
    dy.T_hot_leg      = firstOrder(s.T_upper_plenum,     s.T_hot_leg,      p.hot_leg_mass_kg);
    dy.T_cold_leg     = firstOrder(T_inlet_K,            s.T_cold_leg,     p.cold_leg_mass_kg);

    return dy;
}

}  // namespace astara::reactor
