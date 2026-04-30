#ifndef ASTARA_SG_ALI_STEAM_GENERATOR_HPP
#define ASTARA_SG_ALI_STEAM_GENERATOR_HPP

/**
 * @file   AliSteamGenerator.hpp
 * @brief  U-tube recirculating steam generator after Ali (1985), as adapted
 *         in Naghedolfeizi (1990) Appendix B and Section 3.2.
 *
 * The model lumps the steam generator into four primary nodes (water inside
 * the U-tubes), four metal-tube nodes, and four secondary lumps (subcooled,
 * boiling, riser/separator, drum/downcomer water + steam dome).  The state
 * vector has 18 components:
 *
 *   primary side (5) :  T_pi (inlet plenum), T_p1..T_p4
 *   metal     (4) :     T_m1, T_m2, T_m3, T_m4
 *   secondary  (9):
 *      T_sub   sub-cooled water temperature
 *      L_s1    sub-cooled region length
 *      h_b     boiling-region average enthalpy
 *      x_e     exit quality of the boiling region
 *      P       SG secondary pressure
 *      L_dw    drum water level
 *      T_dw    drum water temperature
 *      rho_r   riser/separator average density
 *      T_dc    downcomer temperature (recirculation loop)
 *
 * # Inputs (set externally, e.g. by the reactor + controller)
 *   - primary inlet temperature  T_pi_in    [K]   (= reactor hot-leg-out)
 *   - primary mass flow rate     W_p        [kg/s]
 *   - feedwater mass flow rate   W_fw       [kg/s]
 *   - feedwater enthalpy         h_fw       [J/kg]
 *   - steam-line pressure (back) P_steam_bc [Pa]  (sets the steam valve)
 *
 * # Outputs the loop integrator consumes
 *   - primary outlet temperature `primaryOutletTemperatureK()` (= cold-leg-in)
 *   - steam mass flow rate
 *   - drum level (for control)
 *   - SG pressure
 *
 * # Modelling notes
 *
 * The thesis equations are written in Ali's notation with linearised water
 * properties (eq. B.22-B.34).  The implementation here uses the same
 * structural form but evaluates the constitutive properties through the
 * `WaterProperties` abstraction, defaulting to `IF97Water`.  This deviates
 * from the thesis in the *values* of the constitutive coefficients but not
 * in the *structure* of the dynamics, and is consistent with the spirit of
 * Section 3.2.1 (where the author replaces the linearised properties with
 * a more accurate steam-table at the cost of comparability with Ali).
 *
 * For users who need to reproduce Ali's published frequency response
 * exactly, construct the SG with a `LinearizedWater` instance fit around
 * the operating pressure (typically 6.9 MPa = 1000 psia).
 *
 * # Limitations
 *
 *   - The model is valid for pressures 5-8 MPa and feedwater temperatures
 *     above 470 K (typical PWR secondary side).
 *   - The recirculation ratio is treated as a slowly-varying algebraic
 *     variable (eq. B.34); the dynamics of the riser are quasi-static.
 *   - Single tube-bundle representation: no separate hot-leg / cold-leg
 *     U-tube columns.  Adequate for full-power / low-power transients but
 *     not for asymmetric perturbations.
 *
 * @cite Ali, M. R. A. (1985).  "Lumped Parameter, State Variable Dynamic
 *       Models for U-Tube Recirculation Type Nuclear Steam Generators."
 *       PhD diss., Univ. of Tennessee, Knoxville.
 * @cite Naghedolfeizi (1990), Section 3.2 + Appendix B.
 */

#include <stdexcept>
#include <vector>

namespace astara::props { class WaterProperties; }

namespace astara::sg {

/**
 * @brief Geometry, mass, and heat-transfer parameters of an Ali Model D SG.
 *
 * Defaults reflect a Westinghouse Model D5 U-tube SG (Naghedolfeizi
 * Table A.4): four-loop plant, 3411 MW thermal, ~ 100 m^2 metal area per
 * tube bundle, ~ 7 MPa secondary, ~ 1700 t total water inventory.
 */
struct AliSteamGeneratorParameters {
    // Primary side
    double primary_mass_per_node_kg = 0.0;     ///< M_pi  primary inlet plenum mass (proxy)
    double primary_mass_inlet_plenum_kg = 0.0; ///< M_pi  inlet-plenum water mass
    double primary_cp_J_per_kgK     = 0.0;     ///< C_p1, primary fluid cp (constant)

    // Metal tubes
    double metal_mass_per_node_kg   = 0.0;     ///< M_m,  per metal node
    double metal_cp_J_per_kgK       = 0.0;     ///< C_m,  metal cp
    double area_pm_per_node_m2      = 0.0;     ///< S_pm, prim->metal area per node
    double overall_h_pm_W_per_m2K   = 0.0;     ///< U_pm, prim->metal coefficient
    double area_ms_per_node_m2      = 0.0;     ///< S_ms, metal->secondary area per node
    double overall_h_ms_W_per_m2K   = 0.0;     ///< U_ms, metal->secondary coefficient

    // Secondary side
    double sec_flow_area_m2         = 0.0;     ///< A_fs, secondary flow area in U-tube region
    double drum_water_area_m2       = 0.0;     ///< A_dw, drum water effective area
    double tube_bundle_height_m     = 0.0;     ///< L,    effective height of U-tubes
    double downcomer_length_m       = 0.0;     ///< L_d,  downcomer length
    double downcomer_volume_m3      = 0.0;     ///< V_dw, downcomer volume

    // Recirculation
    double recirc_pressure_drop_coeff = 0.0;   ///< C_l,  recirculation loop dP coefficient

    // Steam valve
    double steam_valve_coefficient    = 0.0;   ///< C_v,  steam valve admittance

    /// Throw on any non-positive required parameter.
    void validate() const;

    /// Westinghouse Model D5 four-loop default parameters in SI.
    static AliSteamGeneratorParameters westinghouseModelD5();
};

/**
 * @brief State of the Ali Model D SG (18 components).
 */
struct AliSteamGeneratorState {
    double t_s = 0.0;                          ///< simulation time, s

    // Primary side (5 components)
    double T_pi = 0.0;                         ///< inlet plenum temperature, K
    double T_p1 = 0.0;                         ///< primary node 1 (hot leg side, top), K
    double T_p2 = 0.0;                         ///< primary node 2, K
    double T_p3 = 0.0;                         ///< primary node 3, K
    double T_p4 = 0.0;                         ///< primary node 4 (cold leg side, top), K

    // Metal (4 components)
    double T_m1 = 0.0;                         ///< metal node 1 (sub-cooled boundary), K
    double T_m2 = 0.0;                         ///< metal node 2 (sub-cooled boundary), K
    double T_m3 = 0.0;                         ///< metal node 3 (boiling region),   K
    double T_m4 = 0.0;                         ///< metal node 4 (boiling region),   K

    // Secondary (9 components)
    double T_sub  = 0.0;                       ///< sub-cooled liquid temperature, K
    double L_s1   = 0.0;                       ///< sub-cooled region length, m
    double h_b    = 0.0;                       ///< boiling-region average enthalpy, J/kg
    double x_e    = 0.0;                       ///< exit quality of the boiling region, -
    double P      = 0.0;                       ///< SG secondary pressure, Pa
    double L_dw   = 0.0;                       ///< drum water level, m
    double T_dw   = 0.0;                       ///< drum water temperature, K
    double rho_r  = 0.0;                       ///< riser / separator density, kg/m^3
    double T_dc   = 0.0;                       ///< downcomer temperature, K

    /// Number of state variables (18).
    static constexpr std::size_t kSize = 18;

    /// Vector-space arithmetic for RK4 (excludes time).
    AliSteamGeneratorState& operator+=(const AliSteamGeneratorState& rhs) noexcept;
    AliSteamGeneratorState& operator-=(const AliSteamGeneratorState& rhs) noexcept;
    AliSteamGeneratorState& operator*=(double s) noexcept;
};

inline AliSteamGeneratorState operator+(AliSteamGeneratorState a, const AliSteamGeneratorState& b) noexcept { a += b; return a; }
inline AliSteamGeneratorState operator-(AliSteamGeneratorState a, const AliSteamGeneratorState& b) noexcept { a -= b; return a; }
inline AliSteamGeneratorState operator*(AliSteamGeneratorState a, double s) noexcept                       { a *= s; return a; }
inline AliSteamGeneratorState operator*(double s, AliSteamGeneratorState a) noexcept                       { a *= s; return a; }

/**
 * @brief External boundary inputs to the steam generator.
 */
struct AliSteamGeneratorInputs {
    double primary_inlet_temperature_K = 0.0;  ///< T_pi_in (= reactor hot-leg-out)
    double primary_mass_flow_kg_s      = 0.0;  ///< W_p
    double feedwater_mass_flow_kg_s    = 0.0;  ///< W_fw
    double feedwater_enthalpy_J_kg     = 0.0;  ///< h_fw
    double steam_line_pressure_Pa      = 0.0;  ///< downstream pressure for valve flow
};

/**
 * @brief U-tube recirculating steam generator (Ali Model D).
 */
class AliSteamGenerator {
public:
    /**
     * @brief Construct with given parameters and a water-properties service.
     * @param params  geometry, mass, and heat-transfer parameters
     * @param props   water-properties service (must outlive the SG object)
     */
    AliSteamGenerator(AliSteamGeneratorParameters params,
                      const props::WaterProperties* props);

    /**
     * @brief Initialise to a representative full-power steady state.
     *
     * Sets all 18 states to values consistent with rated full-power
     * operation (PWR ~ 7 MPa secondary, T_pi_in ~ 597 K).  Useful as a
     * starting condition for transient calculations.  The exact values
     * solve a simplified energy-balance system; transients should converge
     * to a slightly different operating point if the parameters differ
     * from the nominal ones.
     */
    void initialiseSteadyState(double T_pi_in_K,
                               double primary_mass_flow_kg_s,
                               double sg_pressure_Pa,
                               double drum_level_m);

    /// Read/write external inputs.
    AliSteamGeneratorInputs&       inputs()       noexcept { return inputs_; }
    const AliSteamGeneratorInputs& inputs() const noexcept { return inputs_; }

    /// Single RK4 step.
    void timeStep(double dt);

    // ---------- Outputs ----------
    const AliSteamGeneratorState&    state()  const noexcept { return state_; }
    const AliSteamGeneratorParameters& params() const noexcept { return p_; }

    /// Primary outlet temperature (= cold-leg-in for the reactor).
    double primaryOutletTemperatureK() const noexcept { return state_.T_p4; }

    /// Steam mass flow rate through the steam valve, kg/s.
    double steamMassFlow_kg_s() const noexcept;

    /// Total heat removed from the primary side, W.
    double primaryHeatLoad_W() const noexcept;

    /// Right-hand side at given state.  Public so tests can inspect the
    /// derivative without stepping.
    AliSteamGeneratorState evaluateDerivative(const AliSteamGeneratorState& s) const;

private:
    AliSteamGeneratorParameters    p_;
    const props::WaterProperties*  props_;
    AliSteamGeneratorState         state_;
    AliSteamGeneratorInputs        inputs_;
};

}  // namespace astara::sg

#endif  // ASTARA_SG_ALI_STEAM_GENERATOR_HPP
