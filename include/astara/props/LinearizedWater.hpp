#ifndef ASTARA_PROPS_LINEARIZED_WATER_HPP
#define ASTARA_PROPS_LINEARIZED_WATER_HPP

/**
 * @file   LinearizedWater.hpp
 * @brief  Pressure-linearised water/steam properties (thesis-style).
 *
 * Implements eq. 3.32 of Naghedolfeizi (1990) and the associated
 * Appendix-B constitutive equations:
 *
 * @f[
 *      \phi(P) \;=\; X_\phi \;+\; K_\phi \, P
 * @f]
 *
 * for each saturation property `phi`.  This is a deliberate, faithful
 * implementation of the Ali Model D constitutive layer used in the thesis
 * Steam Generator: it is **not** physically accurate on its own, but it is
 * the property model the Ali model was derived from.  Replacing it with IF97
 * inside the Ali SG would change the dynamics in ways that are not in the
 * literature.
 *
 * The slope/intercept coefficients are fit to IAPWS data over an operating
 * window (default: 5 to 8 MPa around 7 MPa, the typical Westinghouse SG
 * secondary pressure).  An alternative window may be specified at
 * construction.  Single-phase properties are returned at saturation values
 * for backward compatibility -- the Ali model never queries them at non-
 * saturation states.
 *
 * @cite Naghedolfeizi, M. (1990). "Dynamic Modeling of a Pressurized Water
 *       Reactor Plant for Diagnostics and Control." MSc thesis, Univ. of
 *       Tennessee, Knoxville.  Eq. (3.32) and Appendix B.
 * @cite Ali, M. R. A. (1985). "Lumped Parameter, State Variable Dynamic
 *       Models for U-Tube Recirculation Type Nuclear Steam Generators." PhD
 *       dissertation, Univ. of Tennessee, Knoxville.
 */

#include "astara/props/WaterProperties.hpp"

namespace astara::props {

/**
 * @brief Linear-in-pressure saturation properties around an operating point.
 *
 * Constants are stored as `(X, K)` pairs so that property = X + K*P, with
 * P in pascals and the property in SI units described in `WaterProperties`.
 *
 * The default constructor produces coefficients fit to IF97 around 7 MPa.
 * A user can also pass a `WaterProperties` reference to compute the
 * coefficients automatically by central-difference around a chosen P0.
 */
class LinearizedWater final : public WaterProperties {
public:
    /**
     * @brief Coefficients of the linearised property model.  Each property
     *        is computed as `X + K * P` (P in Pa).
     */
    struct Coefficients {
        // Saturation curve
        double X_Tsat = 0.0;  double K_Tsat = 0.0;   ///< Tsat(P)  [K, K/Pa]
        double X_hf   = 0.0;  double K_hf   = 0.0;   ///< h_f(P)   [J/kg, J/(kg*Pa)]
        double X_hg   = 0.0;  double K_hg   = 0.0;   ///< h_g(P)   [J/kg, J/(kg*Pa)]
        double X_vf   = 0.0;  double K_vf   = 0.0;   ///< v_f(P)   [m^3/kg, m^3/(kg*Pa)]
        double X_vg   = 0.0;  double K_vg   = 0.0;   ///< v_g(P)   [m^3/kg, m^3/(kg*Pa)]
        double X_rhof = 0.0;  double K_rhof = 0.0;   ///< rho_f(P) [kg/m^3, kg/(m^3*Pa)]
        double X_rhog = 0.0;  double K_rhog = 0.0;   ///< rho_g(P) [kg/m^3, kg/(m^3*Pa)]
    };

    /**
     * @brief Construct with default coefficients fit around 7 MPa
     *        (typical Westinghouse SG secondary pressure).
     */
    LinearizedWater();

    /**
     * @brief Construct with user-supplied coefficients.
     */
    explicit LinearizedWater(const Coefficients& c);

    /**
     * @brief Build a `LinearizedWater` by central-differencing another
     *        property service around the given operating pressure.
     * @param ref     reference property service (e.g. an IF97Water)
     * @param P0_Pa   centre of the linearisation window
     * @param dP_Pa   half-width used for the central difference
     */
    static LinearizedWater fitAround(const WaterProperties& ref,
                                     double P0_Pa,
                                     double dP_Pa = 5.0e5);

    const Coefficients& coefficients() const noexcept { return c_; }

    // Single-phase properties: in this thesis-style model we return saturation
    // values regardless of T (the Ali model never evaluates them off the
    // saturation curve).  Tests should use `IF97Water` if they need true
    // single-phase properties.
    double density_TP   (double T_K, double P_Pa) const override;
    double enthalpy_TP  (double T_K, double P_Pa) const override;
    double cp_TP        (double T_K, double P_Pa) const override;

    // Saturation
    double saturationTemperature(double P_Pa) const override;
    double saturationPressure   (double T_K)  const override;

    double satLiquidDensity_P (double P_Pa) const override;
    double satVapourDensity_P (double P_Pa) const override;
    double satLiquidEnthalpy_P(double P_Pa) const override;
    double satVapourEnthalpy_P(double P_Pa) const override;

private:
    Coefficients c_;
};

}  // namespace astara::props

#endif  // ASTARA_PROPS_LINEARIZED_WATER_HPP
