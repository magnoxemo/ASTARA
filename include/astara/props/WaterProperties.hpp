#ifndef ASTARA_PROPS_WATER_PROPERTIES_HPP
#define ASTARA_PROPS_WATER_PROPERTIES_HPP

/**
 * @file   WaterProperties.hpp
 * @brief  Abstract interface to water/steam thermodynamic properties.
 *
 * Two implementations are provided:
 *
 *   - `IF97Water`        -- IAPWS Industrial Formulation 1997 (default).  This
 *                          is the same engine CoolProp uses for water and is
 *                          appropriate for production runs.  See IF97Water.hpp.
 *   - `LinearizedWater`  -- pressure-linearised properties as used in the
 *                          thesis (Naghedolfeizi 1990, eq. 3.32 and Ali 1985).
 *                          Used by `AliSteamGenerator` to faithfully reproduce
 *                          the thesis Model D, and as a unit-test reference.
 *
 * Why an interface?  Two reasons:
 *
 *   1. Tests can substitute a deterministic `LinearizedWater` for a
 *      thread-safe `IF97Water` to make analytical comparisons exact.
 *   2. The Ali SG model deliberately uses linearised properties (its
 *      equations are derived from them).  Using IF97 inside that model would
 *      change its dynamics in ways the thesis did not analyse.
 *
 * All quantities are in SI: pressures in Pa, temperatures in K, enthalpies in
 * J/kg, densities in kg/m^3, specific heats in J/(kg*K), specific volumes in
 * m^3/kg.
 */

#include <stdexcept>
#include <string>

namespace astara::props {

/**
 * @brief Thrown when a property query is outside the implementation's valid
 *        range.
 */
class PropertyError : public std::runtime_error {
public:
    explicit PropertyError(const std::string& msg) : std::runtime_error(msg) {}
};

/**
 * @brief Abstract interface to water/steam properties.
 *
 * Subclasses must be **thread-safe for read-only queries** (the integrator
 * may call them from multiple OpenMP threads simultaneously when evaluating
 * derivatives across components).  Both supplied implementations are stateless
 * and therefore trivially thread-safe.
 */
class WaterProperties {
public:
    virtual ~WaterProperties() = default;

    // ---------- Single-phase, by (T, P) ----------

    /// Density at given temperature and pressure, kg/m^3.
    virtual double density_TP    (double T_K, double P_Pa) const = 0;
    /// Specific enthalpy at given temperature and pressure, J/kg.
    virtual double enthalpy_TP   (double T_K, double P_Pa) const = 0;
    /// Isobaric specific heat capacity at given temperature and pressure, J/(kg*K).
    virtual double cp_TP         (double T_K, double P_Pa) const = 0;

    // ---------- Saturation curve ----------

    /// Saturation temperature at given pressure, K.
    virtual double saturationTemperature(double P_Pa) const = 0;
    /// Saturation pressure at given temperature, Pa.
    virtual double saturationPressure   (double T_K)  const = 0;

    /// Saturated-liquid density at given pressure, kg/m^3.
    virtual double satLiquidDensity_P (double P_Pa) const = 0;
    /// Saturated-vapour  density at given pressure, kg/m^3.
    virtual double satVapourDensity_P (double P_Pa) const = 0;
    /// Saturated-liquid enthalpy at given pressure, J/kg.
    virtual double satLiquidEnthalpy_P(double P_Pa) const = 0;
    /// Saturated-vapour  enthalpy at given pressure, J/kg.
    virtual double satVapourEnthalpy_P(double P_Pa) const = 0;

    /// Latent heat of vaporisation at given pressure, J/kg.
    /// (Default implementation: hvap_p - hliq_p.)
    virtual double latentHeat_P(double P_Pa) const {
        return satVapourEnthalpy_P(P_Pa) - satLiquidEnthalpy_P(P_Pa);
    }

    /// Saturated-liquid specific volume at given pressure, m^3/kg.
    virtual double satLiquidSpecificVolume_P(double P_Pa) const {
        return 1.0 / satLiquidDensity_P(P_Pa);
    }
    /// Saturated-vapour  specific volume at given pressure, m^3/kg.
    virtual double satVapourSpecificVolume_P(double P_Pa) const {
        return 1.0 / satVapourDensity_P(P_Pa);
    }

    /// Two-phase enthalpy at given (P, quality x), J/kg.
    virtual double twoPhaseEnthalpy_PQ(double P_Pa, double x) const {
        return satLiquidEnthalpy_P(P_Pa) + x * latentHeat_P(P_Pa);
    }
};

}  // namespace astara::props

#endif  // ASTARA_PROPS_WATER_PROPERTIES_HPP
