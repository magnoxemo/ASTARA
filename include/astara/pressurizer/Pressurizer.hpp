#ifndef ASTARA_PRESSURIZER_PRESSURIZER_HPP
#define ASTARA_PRESSURIZER_PRESSURIZER_HPP

/**
 * @file   Pressurizer.hpp
 * @brief  Two-region (water + steam) pressurizer with heater and spray.
 *
 * Implements the pressurizer model in Naghedolfeizi (1990), Appendix C
 * (which itself adapts an earlier model with explicit steam-compressibility
 * coupling).  The pressurizer is treated as a vertical cylindrical vessel
 * with a saturated water region at the bottom and a saturated steam dome
 * at the top.  Heat is added by the electric heater, removed by the
 * pressurizer spray (cold water sprayed into the steam space, condensing
 * it), and mass is exchanged with the primary loop via the surge line.
 *
 * State (2 ODEs): pressurizer pressure `P` and water level `L_w`.
 *
 * # Mass balances
 * Water:    dM_w/dt = W_sg + W_sp - W_co              (surge in + spray - condensation)
 * Steam:    dM_s/dt =                W_co              (only condensation/evaporation)
 *
 * with `M_w = rho_w * A_pr * L_w` and `M_s = rho_s * A_pr * (L - L_w)`.
 *
 * # Pressure (energy balance on the steam volume + compressibility)
 * Eq. (C.18)-(C.20) of the thesis give a closed-form expression for
 * `dP/dt = NUM/DEN`, with NUM containing heater power, spray-cooling, and
 * surge enthalpy; DEN containing the steam-compressibility (M_s K_p1 +
 * M_w K_p2) terms.
 *
 * The implementation below uses the simpler but equivalent reduced form
 * obtained when the linearised property model is replaced by IF97 lookups:
 * (`v_g`, `h_f`, `h_g`, `rho_w`, `rho_s` come from `WaterProperties`).
 *
 * # Inputs (set externally between time steps)
 *   - `surge_mass_flow`  W_sg   (positive into pressurizer = water from hot leg)
 *   - `spray_mass_flow`  W_sp   (positive = cold water injected into steam)
 *   - `heater_power_W`   Q      (positive = heat added by heaters)
 *   - `spray_enthalpy`   h_sp   (enthalpy of the spray water, J/kg)
 *
 * @cite Naghedolfeizi (1990), Section 3.3, Appendix C.
 */

#include <memory>
#include <stdexcept>

namespace astara::props { class WaterProperties; }

namespace astara::pressurizer {

/**
 * @brief Geometric parameters of the pressurizer.
 */
struct PressurizerParameters {
    double cross_section_area_m2 = 0.0;   ///< A_pr, effective horizontal area
    double total_height_m        = 0.0;   ///< L, total internal height of vessel

    void validate() const {
        if (!(cross_section_area_m2 > 0.0)) {
            throw std::invalid_argument(
                    "PressurizerParameters: cross_section_area_m2 must be positive");
        }
        if (!(total_height_m > 0.0)) {
            throw std::invalid_argument(
                    "PressurizerParameters: total_height_m must be positive");
        }
    }
};

/**
 * @brief State of the pressurizer.
 */
struct PressurizerState {
    double t_s            = 0.0;   ///< time, s
    double pressure_Pa    = 0.0;   ///< P, Pa
    double water_level_m  = 0.0;   ///< L_w, m

    /// Vector-space arithmetic used by RK4 (excludes time).
    PressurizerState& operator+=(const PressurizerState& rhs) noexcept {
        pressure_Pa   += rhs.pressure_Pa;
        water_level_m += rhs.water_level_m;
        return *this;
    }
    PressurizerState& operator-=(const PressurizerState& rhs) noexcept {
        pressure_Pa   -= rhs.pressure_Pa;
        water_level_m -= rhs.water_level_m;
        return *this;
    }
    PressurizerState& operator*=(double s) noexcept {
        pressure_Pa   *= s;
        water_level_m *= s;
        return *this;
    }
};

inline PressurizerState operator+(PressurizerState a, const PressurizerState& b) noexcept { a += b; return a; }
inline PressurizerState operator-(PressurizerState a, const PressurizerState& b) noexcept { a -= b; return a; }
inline PressurizerState operator*(PressurizerState a, double s) noexcept                  { a *= s; return a; }
inline PressurizerState operator*(double s, PressurizerState a) noexcept                  { a *= s; return a; }

/**
 * @brief External boundary inputs.  Set between time steps by the caller
 *        (or by the controller).
 */
struct PressurizerInputs {
    double surge_mass_flow_kg_s = 0.0;     ///< W_sg, positive = water in (out-surge < 0)
    double spray_mass_flow_kg_s = 0.0;     ///< W_sp, positive = cold spray flow on
    double heater_power_W       = 0.0;     ///< Q,    positive = heat added
    double spray_enthalpy_J_kg  = 0.0;     ///< h_sp, enthalpy of the spray water
    double surge_enthalpy_J_kg  = 0.0;     ///< h_sg, enthalpy of the surge stream
                                           ///<       (used only when W_sg > 0; on out-surge
                                           ///<        the water leaving carries the saturated
                                           ///<        liquid enthalpy at current P)
};

/**
 * @brief Two-state lumped-parameter pressurizer.
 *
 * Properties are obtained from a `WaterProperties` service, which the
 * caller supplies at construction.  Use `IF97Water` for production runs
 * (default) or `LinearizedWater` to reproduce the Ali-style derivation.
 */
class Pressurizer {
public:
    /**
     * @brief Construct with given geometry and a property service.
     * @param params   pressurizer geometry
     * @param props    water properties (must outlive the Pressurizer object)
     */
    Pressurizer(PressurizerParameters params,
                const props::WaterProperties* props);

    /**
     * @brief Initialise to steady state at given pressure and water level.
     */
    void initialiseSteadyState(double P_Pa, double L_w_m);

    PressurizerInputs&       inputs()       noexcept { return inputs_; }
    const PressurizerInputs& inputs() const noexcept { return inputs_; }

    /// Single RK4 step.
    void timeStep(double dt);

    const PressurizerState&    state()  const noexcept { return state_; }
    const PressurizerParameters& params() const noexcept { return p_; }

    /// Current mass of water in the pressurizer, kg.
    double waterMass_kg() const noexcept;
    /// Current mass of steam in the pressurizer, kg.
    double steamMass_kg() const noexcept;

    /// Right-hand side at given state (used by tests + integrator).
    /// Returns dP/dt in Pa/s and dLw/dt in m/s as the two-component state.
    PressurizerState evaluateDerivative(const PressurizerState& s) const;

private:
    PressurizerParameters         p_;
    const props::WaterProperties* props_;
    PressurizerState              state_;
    PressurizerInputs             inputs_;
};

}  // namespace astara::pressurizer

#endif  // ASTARA_PRESSURIZER_PRESSURIZER_HPP
