/**
 * @file   LinearizedWater.cpp
 * @brief  Implementation of pressure-linearised water/steam properties.
 */

#include "astara/props/LinearizedWater.hpp"
#include "astara/props/IF97Water.hpp"

#include <stdexcept>

namespace astara::props {

namespace {
/**
 * @brief Default coefficients fit by central difference of IF97 around
 *        P0 = 7 MPa, dP = 0.5 MPa.  Recomputed at static-init time so they
 *        track any changes in the IF97 backend.
 */
LinearizedWater::Coefficients defaultCoefficients() {
    IF97Water if97;
    constexpr double P0 = 7.0e6;
    constexpr double dP = 5.0e5;

    auto fit = [&](auto fn) {
        const double v_lo = fn(P0 - dP);
        const double v_hi = fn(P0 + dP);
        const double K = (v_hi - v_lo) / (2.0 * dP);
        const double v0 = 0.5 * (v_lo + v_hi);
        const double X = v0 - K * P0;
        return std::pair<double, double>{X, K};
    };

    LinearizedWater::Coefficients c;
    {
        auto [X, K] = fit([&](double P){ return if97.saturationTemperature(P); });
        c.X_Tsat = X; c.K_Tsat = K;
    }
    {
        auto [X, K] = fit([&](double P){ return if97.satLiquidEnthalpy_P(P); });
        c.X_hf = X; c.K_hf = K;
    }
    {
        auto [X, K] = fit([&](double P){ return if97.satVapourEnthalpy_P(P); });
        c.X_hg = X; c.K_hg = K;
    }
    {
        auto [X, K] = fit([&](double P){ return 1.0 / if97.satLiquidDensity_P(P); });
        c.X_vf = X; c.K_vf = K;
    }
    {
        auto [X, K] = fit([&](double P){ return 1.0 / if97.satVapourDensity_P(P); });
        c.X_vg = X; c.K_vg = K;
    }
    {
        auto [X, K] = fit([&](double P){ return if97.satLiquidDensity_P(P); });
        c.X_rhof = X; c.K_rhof = K;
    }
    {
        auto [X, K] = fit([&](double P){ return if97.satVapourDensity_P(P); });
        c.X_rhog = X; c.K_rhog = K;
    }
    return c;
}
}  // namespace

LinearizedWater::LinearizedWater() : c_(defaultCoefficients()) {}

LinearizedWater::LinearizedWater(const Coefficients& c) : c_(c) {}

LinearizedWater LinearizedWater::fitAround(const WaterProperties& ref,
                                           double P0_Pa, double dP_Pa) {
    if (dP_Pa <= 0.0) {
        throw std::invalid_argument("LinearizedWater::fitAround: dP must be positive");
    }
    auto fit = [&](auto fn) {
        const double v_lo = fn(P0_Pa - dP_Pa);
        const double v_hi = fn(P0_Pa + dP_Pa);
        const double K = (v_hi - v_lo) / (2.0 * dP_Pa);
        const double v0 = 0.5 * (v_lo + v_hi);
        const double X = v0 - K * P0_Pa;
        return std::pair<double, double>{X, K};
    };
    Coefficients c;
    { auto [X,K] = fit([&](double P){ return ref.saturationTemperature(P); });   c.X_Tsat=X; c.K_Tsat=K; }
    { auto [X,K] = fit([&](double P){ return ref.satLiquidEnthalpy_P(P);    });   c.X_hf  =X; c.K_hf  =K; }
    { auto [X,K] = fit([&](double P){ return ref.satVapourEnthalpy_P(P);    });   c.X_hg  =X; c.K_hg  =K; }
    { auto [X,K] = fit([&](double P){ return 1.0/ref.satLiquidDensity_P(P); });   c.X_vf  =X; c.K_vf  =K; }
    { auto [X,K] = fit([&](double P){ return 1.0/ref.satVapourDensity_P(P); });   c.X_vg  =X; c.K_vg  =K; }
    { auto [X,K] = fit([&](double P){ return ref.satLiquidDensity_P(P);     });   c.X_rhof=X; c.K_rhof=K; }
    { auto [X,K] = fit([&](double P){ return ref.satVapourDensity_P(P);     });   c.X_rhog=X; c.K_rhog=K; }
    return LinearizedWater(c);
}

// ---------- single-phase: return saturation values (Ali model assumption) ----

double LinearizedWater::density_TP(double /*T_K*/, double P_Pa) const {
    return c_.X_rhof + c_.K_rhof * P_Pa;
}
double LinearizedWater::enthalpy_TP(double /*T_K*/, double P_Pa) const {
    return c_.X_hf + c_.K_hf * P_Pa;
}
double LinearizedWater::cp_TP(double /*T_K*/, double /*P_Pa*/) const {
    // Constant cp for liquid water at PWR-secondary conditions.  The Ali
    // model assumes constant cp for the primary fluid and the secondary
    // sub-cooled region; this matches the thesis assumption.
    return 5400.0;  // J/(kg*K), thesis value: 1.39 Btu/(lbm*F) ~= 5821 J/(kg*K).
                    // We use 5400 to match the IAPWS value at the operating point.
}

// ---------- saturation curve ----------

double LinearizedWater::saturationTemperature(double P_Pa) const {
    return c_.X_Tsat + c_.K_Tsat * P_Pa;
}
double LinearizedWater::saturationPressure(double T_K) const {
    if (c_.K_Tsat == 0.0) {
        throw PropertyError("LinearizedWater::saturationPressure: Tsat slope is zero");
    }
    return (T_K - c_.X_Tsat) / c_.K_Tsat;
}
double LinearizedWater::satLiquidDensity_P (double P_Pa) const { return c_.X_rhof + c_.K_rhof * P_Pa; }
double LinearizedWater::satVapourDensity_P (double P_Pa) const { return c_.X_rhog + c_.K_rhog * P_Pa; }
double LinearizedWater::satLiquidEnthalpy_P(double P_Pa) const { return c_.X_hf   + c_.K_hf   * P_Pa; }
double LinearizedWater::satVapourEnthalpy_P(double P_Pa) const { return c_.X_hg   + c_.K_hg   * P_Pa; }

}  // namespace astara::props
