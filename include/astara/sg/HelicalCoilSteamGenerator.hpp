#ifndef ASTARA_SG_HELICAL_COIL_STEAM_GENERATOR_HPP
#define ASTARA_SG_HELICAL_COIL_STEAM_GENERATOR_HPP

/**
 * @file   HelicalCoilSteamGenerator.hpp
 * @brief  Once-through helical-coil steam generator with a moving-boundary
 *         secondary-side model (subcooled / two-phase / superheated).
 *
 * Implements the model published in
 *
 *   S. E. Arda and K. E. Holbert,
 *   "A dynamic model of a passively cooled small modular reactor for
 *   controller design purposes,"
 *   Nuclear Engineering and Design 289 (2015) 218-230.
 *   https://doi.org/10.1016/j.nucengdes.2015.04.026
 *
 * Hereafter cited as Arda (2015).  The default parameter set reproduces the
 * NuScale-like 160 MWth helical-coil SG modelled in that paper.
 *
 * # Model summary (Arda 2015 sections 5.1 - 5.6)
 *
 * The secondary side of the SG is divided into three regions that move
 * along the tube axis as operating conditions change:
 *
 *      L_1  : sub-cooled liquid (feedwater entry)
 *      L_2  : two-phase saturated mixture
 *      L_3  : superheated steam (= L_T - L_1 - L_2; algebraic)
 *
 * Each region exchanges heat with a corresponding metal-tube node and a
 * primary-coolant node.  The 10 dynamic states are
 *
 *      x = [ L_1, L_2, p_S, h_o,
 *            T_M1, T_M2, T_M3,
 *            T_P1, T_P2, T_P3 ]^T
 *
 * with `p_S` the secondary pressure, `h_o` the steam-outlet enthalpy,
 * `T_Mk` and `T_Pk` the metal and primary average temperatures of the
 * three regions (Arda Fig. 3).
 *
 * The governing equations have the form
 *
 *      D(x, u) * dx/dt = f(x, u)
 *
 * (Arda eq. 44) where D is a 10x10 sparse matrix whose entries are listed
 * in Arda Table 3.  At each time step we form D and f, solve the linear
 * system D dx = f for the time derivative, and advance with classical
 * RK4.  The matrix is well-conditioned at typical SMR operating points
 * but does become near-singular if a region length collapses to zero;
 * the integrator throws in that case rather than producing garbage.
 *
 * # Inputs
 *
 *      - primary inlet temperature        T_Pi       [K]
 *      - primary mass flow rate           m_dot_P    [kg/s]
 *      - secondary inlet (feedwater) flow m_dot_Si   [kg/s]
 *      - feedwater enthalpy               h_i        [J/kg]
 *      - steam outlet flow                m_dot_So   [kg/s]
 *
 * `m_dot_So` is treated as a *boundary condition* (set by an external
 * valve or load demand), as in the paper.  In a coupled simulation it
 * may be obtained from a turbine model; the helper `setOutletFromValve()`
 * provides a simple pressure-driven valve closure for stand-alone use.
 *
 * # Outputs
 *
 *      - primary outlet temperature       primaryOutletTemperatureK()
 *      - steam pressure                   state().p_S
 *      - steam outlet enthalpy            state().h_o
 *      - region lengths                   state().L_1, L_2, L_T - L_1 - L_2
 *
 * # Validity
 *
 * The Arda model assumes:
 *
 *      - All three regions exist (L_1, L_2, L_3 all > 0).  If any region
 *        length crosses zero the moving-boundary topology has changed and
 *        the model is no longer valid (would need re-derivation with
 *        merged regions).  The integrator throws on this condition.
 *      - Primary and secondary pressures are uniform across each region.
 *      - The two-phase region is at thermal equilibrium.
 *      - Heat-transfer coefficients are evaluated at steady-state
 *        Reynolds numbers for the typical operating point and held
 *        constant during transients (a common simplification for
 *        controller-design models -- see Arda section 5.5).
 *
 * For pressures and temperatures that violate any of the above (e.g. a
 * scram that completely empties the superheated region), use a different
 * topology -- e.g. the Ali Model D recirculation-loop SG.
 */

#include <array>
#include <stdexcept>

namespace astara::props { class WaterProperties; }

namespace astara::sg {

/**
 * @brief Geometry, mass, and heat-transfer parameters for the helical-coil SG.
 *
 * Defaults reproduce a NuScale-like 160 MWth, two-SG-per-vessel design
 * (Arda 2015 Table 1).  Each SG models the aggregate thermal behaviour of
 * 506 tubes in parallel as a single equivalent tube.
 */
struct HelicalCoilSteamGeneratorParameters {
    // Tube-bundle geometry (single equivalent tube)
    double tube_total_length_m         = 0.0;   ///< L_T  total tube length
    double tube_inner_diameter_m       = 0.0;   ///< d_i
    double tube_outer_diameter_m       = 0.0;   ///< d_o
    double tube_count                  = 0.0;   ///< number of tubes (for scaling)

    // Cross-sectional flow areas, summed over all tubes
    double secondary_flow_area_m2      = 0.0;   ///< A_S  (sum over tubes)
    double primary_flow_area_m2        = 0.0;   ///< A_P  outside the bundle
    double tube_metal_area_m2          = 0.0;   ///< A_M  metal cross-section

    // Material properties
    double tube_metal_density_kg_m3    = 0.0;   ///< rho_M
    double tube_metal_cp_J_per_kgK     = 0.0;   ///< c_pM (Inconel-690)
    double primary_density_kg_m3       = 0.0;   ///< rho_P
    double primary_cp_J_per_kgK        = 0.0;   ///< c_pP

    // Heat-transfer coefficients (constant during transients)
    double alpha_o_subcooled_W_per_m2K = 0.0;   ///< alpha_o,1
    double alpha_o_twophase_W_per_m2K  = 0.0;   ///< alpha_o,2
    double alpha_o_superheat_W_per_m2K = 0.0;   ///< alpha_o,3
    double alpha_i_subcooled_W_per_m2K = 0.0;   ///< alpha_i,1
    double alpha_i_twophase_W_per_m2K  = 0.0;   ///< alpha_i,2
    double alpha_i_superheat_W_per_m2K = 0.0;   ///< alpha_i,3

    /// Throw on any non-positive required parameter.
    void validate() const;

    /// NuScale-like default parameters (Arda 2015 Table 1 + Table 4).
    static HelicalCoilSteamGeneratorParameters nuscaleSMRTwoSG();
};

/**
 * @brief Dynamic state of the helical-coil SG (10 components).
 *
 * Order matches Arda (2015) eq. 44:
 *      x = [L_1, L_2, p_S, h_o, T_M1, T_M2, T_M3, T_P1, T_P2, T_P3]
 */
struct HelicalCoilSteamGeneratorState {
    double t_s = 0.0;               ///< simulation time, s
    double L_1  = 0.0;              ///< sub-cooled region length, m
    double L_2  = 0.0;              ///< two-phase   region length, m
    double p_S  = 0.0;              ///< secondary pressure, Pa
    double h_o  = 0.0;              ///< steam outlet enthalpy, J/kg
    double T_M1 = 0.0;              ///< metal node 1 average temperature, K
    double T_M2 = 0.0;              ///< metal node 2 average temperature, K
    double T_M3 = 0.0;              ///< metal node 3 average temperature, K
    double T_P1 = 0.0;              ///< primary node 1 average temperature, K
    double T_P2 = 0.0;              ///< primary node 2 average temperature, K
    double T_P3 = 0.0;              ///< primary node 3 average temperature, K

    static constexpr std::size_t kSize = 10;

    /// Vector-space arithmetic for RK4 (excludes time).
    HelicalCoilSteamGeneratorState& operator+=(const HelicalCoilSteamGeneratorState&) noexcept;
    HelicalCoilSteamGeneratorState& operator-=(const HelicalCoilSteamGeneratorState&) noexcept;
    HelicalCoilSteamGeneratorState& operator*=(double) noexcept;

    /// Pack/unpack to a length-10 array (used to form D dx = f).
    std::array<double, 10> toArray() const noexcept;
    void                    fromArray(const std::array<double, 10>& a) noexcept;
};

inline HelicalCoilSteamGeneratorState
operator+(HelicalCoilSteamGeneratorState a, const HelicalCoilSteamGeneratorState& b) noexcept { a += b; return a; }
inline HelicalCoilSteamGeneratorState
operator-(HelicalCoilSteamGeneratorState a, const HelicalCoilSteamGeneratorState& b) noexcept { a -= b; return a; }
inline HelicalCoilSteamGeneratorState
operator*(HelicalCoilSteamGeneratorState a, double s) noexcept { a *= s; return a; }
inline HelicalCoilSteamGeneratorState
operator*(double s, HelicalCoilSteamGeneratorState a) noexcept { a *= s; return a; }

/**
 * @brief External boundary conditions (Arda 2015 eq. 44 input vector u).
 */
struct HelicalCoilSteamGeneratorInputs {
    double primary_inlet_temperature_K = 0.0;   ///< T_Pi
    double primary_mass_flow_kg_s      = 0.0;   ///< m_dot_P
    double feedwater_mass_flow_kg_s    = 0.0;   ///< m_dot_S,i
    double feedwater_enthalpy_J_kg     = 0.0;   ///< h_i  (= h(T_S,i, p_S))
    double steam_outlet_mass_flow_kg_s = 0.0;   ///< m_dot_S,o (boundary)
};

/**
 * @brief Helical-coil moving-boundary SG (Arda 2015).
 */
class HelicalCoilSteamGenerator {
public:
    HelicalCoilSteamGenerator(HelicalCoilSteamGeneratorParameters params,
                              const props::WaterProperties* props);

    /**
     * @brief Initialise to the published 100% steady-state values from
     *        Arda (2015) Table 4.
     *
     * Sets the 10 dynamic states + the inputs to a self-consistent
     * 100%-power operating point: L_1 = 2.90 m, L_2 = 17.60 m,
     * p_S = 3.1 MPa, T_steam = 264 deg C, h_o evaluated from those.
     * Adjusts heat-transfer coefficients and tube-bundle parameters in
     * the order (1) compute steady-state primary-side temperature
     * profile, (2) compute metal-tube profile so that primary-to-metal
     * and metal-to-secondary heat fluxes balance, (3) set h_o from the
     * superheated-region energy balance.
     */
    void initialiseSteadyState(double T_pi_in_K,
                               double primary_mass_flow_kg_s,
                               double feedwater_flow_kg_s,
                               double feedwater_temperature_K);

    HelicalCoilSteamGeneratorInputs&       inputs()       noexcept { return inputs_; }
    const HelicalCoilSteamGeneratorInputs& inputs() const noexcept { return inputs_; }

    /// Single RK4 step.  Each stage forms D, f and solves D dx = f.
    void timeStep(double dt);

    const HelicalCoilSteamGeneratorState&    state()  const noexcept { return state_; }
    const HelicalCoilSteamGeneratorParameters& params() const noexcept { return p_; }

    /// Primary outlet temperature, K (= T_P1, the region nearest the
    /// feedwater inlet at the bottom of the tubes).
    double primaryOutletTemperatureK() const noexcept { return state_.T_P1; }

    /// Steam outlet temperature, K (computed from h_o at p_S).
    double steamOutletTemperatureK() const;

    /// Total heat removed from the primary side, W.
    double primaryHeatLoad_W() const noexcept;

    /// Length of the superheated region (algebraic = L_T - L_1 - L_2).
    double superheatedRegionLength_m() const noexcept {
        return std::max(p_.tube_total_length_m - state_.L_1 - state_.L_2, 0.0);
    }

    /// Public derivative evaluation (forms and solves D dx = f).
    HelicalCoilSteamGeneratorState
    evaluateDerivative(const HelicalCoilSteamGeneratorState& s) const;

private:
    HelicalCoilSteamGeneratorParameters p_;
    const props::WaterProperties*       props_;
    HelicalCoilSteamGeneratorState      state_;
    HelicalCoilSteamGeneratorInputs     inputs_;
};

}  // namespace astara::sg

#endif  // ASTARA_SG_HELICAL_COIL_STEAM_GENERATOR_HPP
