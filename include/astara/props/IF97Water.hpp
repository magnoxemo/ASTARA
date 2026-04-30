#ifndef ASTARA_PROPS_IF97_WATER_HPP
#define ASTARA_PROPS_IF97_WATER_HPP

/**
 * @file   IF97Water.hpp
 * @brief  IAPWS Industrial Formulation 1997 water/steam properties.
 *
 * Wraps the IF97 header-only library (the same engine used by CoolProp for
 * water) behind the `WaterProperties` interface.  Stateless and trivially
 * thread-safe.
 *
 * @cite IAPWS R7-97(2012). Revised Release on the IAPWS Industrial Formulation
 *       1997 for the Thermodynamic Properties of Water and Steam.
 */

#include "astara/props/WaterProperties.hpp"

namespace astara::props {

/**
 * @brief Production water/steam properties via IAPWS-IF97.
 *
 * Valid range: 273.15 K <= T <= 1073.15 K (Region 1, 2, 4); pressure up to
 * 100 MPa.  PWR primary (~155 bar, ~575 K) and secondary (~70 bar, ~565 K)
 * conditions are well inside Region 1 and the saturation curve respectively.
 *
 * Throws `PropertyError` on out-of-range queries (typed error, not silent
 * NaN -- this is what kept the user's original Python prototype "untrusted",
 * per the user's own description of the legacy code).
 */
class IF97Water final : public WaterProperties {
public:
    /**
     * @brief Construct an IF97 property service.
     * @param strict  if true, raise PropertyError on out-of-range queries;
     *                if false, allow IF97 to extrapolate (not recommended).
     */
    explicit IF97Water(bool strict = true) noexcept;

    double density_TP    (double T_K, double P_Pa) const override;
    double enthalpy_TP   (double T_K, double P_Pa) const override;
    double cp_TP         (double T_K, double P_Pa) const override;

    double saturationTemperature(double P_Pa) const override;
    double saturationPressure   (double T_K)  const override;

    double satLiquidDensity_P (double P_Pa) const override;
    double satVapourDensity_P (double P_Pa) const override;
    double satLiquidEnthalpy_P(double P_Pa) const override;
    double satVapourEnthalpy_P(double P_Pa) const override;

private:
    bool strict_;

    void checkSinglePhase (double T_K, double P_Pa, const char* fn) const;
    void checkPressure    (double P_Pa, const char* fn) const;
    void checkTemperature (double T_K,  const char* fn) const;
};

}  // namespace astara::props

#endif  // ASTARA_PROPS_IF97_WATER_HPP
