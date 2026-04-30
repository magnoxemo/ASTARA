#ifndef ASTARA_CORE_UNITS_HPP
#define ASTARA_CORE_UNITS_HPP

/**
 * @file   Units.hpp
 * @brief  SI unit conventions and imperial-to-SI conversions.
 *
 * ASTARA uses SI internally **without exception**:
 *
 * | Quantity            | Unit              | Symbol      |
 * |---------------------|-------------------|-------------|
 * | length              | metre             | m           |
 * | mass                | kilogram          | kg          |
 * | time                | second            | s           |
 * | temperature         | kelvin            | K           |
 * | pressure            | pascal            | Pa          |
 * | mass flow rate      | kg/s              | kg/s        |
 * | enthalpy (specific) | J/kg              | J/kg        |
 * | heat capacity       | J/(kg*K)          | J/(kg*K)    |
 * | density             | kg/m^3            | kg/m^3      |
 * | thermal power       | watt              | W           |
 * | reactivity          | dimensionless     | 1           |
 *
 * The thesis (Naghedolfeizi 1990) and many other PWR references quote design
 * values in mixed imperial-engineering units (psia, lbm, ft, btu, °F).  All
 * such values must be converted at the **boundary** of the simulator (e.g. in
 * the `WestinghousePWR.hpp` parameter table).  Helpers below do the
 * conversions in a way that is `constexpr` and self-documenting.
 */

namespace astara::units {

// ---------- Length ----------
inline constexpr double kFootToMetre   = 0.3048;
inline constexpr double kInchToMetre   = 0.0254;

// ---------- Mass ----------
inline constexpr double kPoundToKg     = 0.45359237;

// ---------- Pressure ----------
/// 1 psi = 6894.757... Pa (exact by definition since 1959 international pound).
inline constexpr double kPsiToPa       = 6894.757293168361;
/// 1 bar = 1e5 Pa.
inline constexpr double kBarToPa       = 1.0e5;

// ---------- Energy / heat ----------
/// IT (International Table) British thermal unit: 1 Btu_IT = 1055.05585262 J.
inline constexpr double kBtuToJoule    = 1055.05585262;

// ---------- Temperature ----------
/**
 * @brief Convert Fahrenheit to Kelvin.
 * @param T_F  temperature in degrees Fahrenheit
 * @return     temperature in kelvin
 */
inline constexpr double fahrenheitToKelvin(double T_F) noexcept {
    return (T_F - 32.0) * 5.0 / 9.0 + 273.15;
}

/**
 * @brief Convert Kelvin to Fahrenheit.
 */
inline constexpr double kelvinToFahrenheit(double T_K) noexcept {
    return (T_K - 273.15) * 9.0 / 5.0 + 32.0;
}

/**
 * @brief Convert Celsius to Kelvin.
 */
inline constexpr double celsiusToKelvin(double T_C) noexcept {
    return T_C + 273.15;
}

/**
 * @brief Convert Kelvin to Celsius.
 */
inline constexpr double kelvinToCelsius(double T_K) noexcept {
    return T_K - 273.15;
}

// ---------- Composite conversions used a lot in the thesis ----------

/// Btu/(lbm * °F) -> J/(kg * K).  Note: a temperature *interval* of 1 °F = 5/9 K,
/// and 1 lbm = 0.45359237 kg.
inline constexpr double kBtuPerLbmFToJperKgK =
        kBtuToJoule / (kPoundToKg * (5.0 / 9.0));

/// Btu/(hr * ft^2 * °F) -> W/(m^2 * K).  1 hr = 3600 s, ft^2 = (0.3048)^2 m^2.
inline constexpr double kBtuPerHrFt2FToWperM2K =
        kBtuToJoule / (3600.0 * (kFootToMetre * kFootToMetre) * (5.0 / 9.0));

/// lbm/hr -> kg/s.
inline constexpr double kLbmPerHrToKgPerS = kPoundToKg / 3600.0;

/// lbm/sec -> kg/s.
inline constexpr double kLbmPerSecToKgPerS = kPoundToKg;

/// ft^3 -> m^3.
inline constexpr double kFt3ToM3 = kFootToMetre * kFootToMetre * kFootToMetre;

/// lbm/ft^3 -> kg/m^3.
inline constexpr double kLbmPerFt3ToKgPerM3 = kPoundToKg / kFt3ToM3;

/// gpm (US gallon per minute) -> m^3/s.  1 US gal = 0.003785411784 m^3.
inline constexpr double kGpmToM3PerS = 0.003785411784 / 60.0;

}  // namespace astara::units

#endif  // ASTARA_CORE_UNITS_HPP
