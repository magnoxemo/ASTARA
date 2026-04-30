/**
 * @file   IF97Water.cpp
 * @brief  Implementation of IAPWS-IF97 water properties via the IF97 library.
 */

#include "astara/props/IF97Water.hpp"

#include "IF97.h"  // header-only IAPWS-IF97 implementation (CoolProp upstream)

#include <cmath>
#include <sstream>
#include <stdexcept>

namespace astara::props {

namespace {
// IAPWS-IF97 validity envelope used in `strict` mode.  These are conservative
// values chosen for PWR work; the formulation itself is valid over a wider
// range but extrapolation outside this box has not been verified for our
// purposes.
constexpr double kPmin = 1.0e3;     ///< 1 kPa lower bound on pressure
constexpr double kPmax = 1.0e8;     ///< 100 MPa upper bound on pressure
constexpr double kTmin = 273.16;    ///< triple-point temperature
constexpr double kTmax = 1073.15;   ///< 800 deg C
}  // namespace

IF97Water::IF97Water(bool strict) noexcept : strict_(strict) {}

void IF97Water::checkPressure(double P_Pa, const char* fn) const {
    if (!strict_) return;
    if (!(P_Pa > kPmin && P_Pa < kPmax)) {
        std::ostringstream os;
        os << "IF97Water::" << fn << ": pressure out of range ["
           << kPmin << ", " << kPmax << "] Pa, got " << P_Pa;
        throw PropertyError(os.str());
    }
}

void IF97Water::checkTemperature(double T_K, const char* fn) const {
    if (!strict_) return;
    if (!(T_K > kTmin && T_K < kTmax)) {
        std::ostringstream os;
        os << "IF97Water::" << fn << ": temperature out of range ["
           << kTmin << ", " << kTmax << "] K, got " << T_K;
        throw PropertyError(os.str());
    }
}

void IF97Water::checkSinglePhase(double T_K, double P_Pa, const char* fn) const {
    checkPressure(P_Pa, fn);
    checkTemperature(T_K, fn);
    // IF97 internally splits into Regions 1 (subcooled liquid), 2 (vapour /
    // superheat), 3 (supercritical), 4 (saturation curve), and 5 (high T).
    // The forward functions handle the dispatch automatically.  The user must
    // not call density_TP at exactly the saturation point because the answer
    // is two-valued there; we allow a 1e-6 tolerance and treat closer queries
    // as an error in strict mode.
    if (P_Pa < 22.064e6) {  // below critical pressure
        const double Tsat = IF97::Tsat97(P_Pa);
        if (std::abs(T_K - Tsat) < 1e-6) {
            std::ostringstream os;
            os << "IF97Water::" << fn << ": query lies on the saturation curve "
               << "(T=" << T_K << " K, Tsat(P)=" << Tsat << " K); use the "
               << "satLiquid*/satVapour* methods instead.";
            throw PropertyError(os.str());
        }
    }
}

double IF97Water::density_TP(double T_K, double P_Pa) const {
    checkSinglePhase(T_K, P_Pa, "density_TP");
    return IF97::rhomass_Tp(T_K, P_Pa);
}

double IF97Water::enthalpy_TP(double T_K, double P_Pa) const {
    checkSinglePhase(T_K, P_Pa, "enthalpy_TP");
    return IF97::hmass_Tp(T_K, P_Pa);
}

double IF97Water::cp_TP(double T_K, double P_Pa) const {
    checkSinglePhase(T_K, P_Pa, "cp_TP");
    return IF97::cpmass_Tp(T_K, P_Pa);
}

double IF97Water::saturationTemperature(double P_Pa) const {
    checkPressure(P_Pa, "saturationTemperature");
    if (P_Pa >= 22.064e6 && strict_) {
        throw PropertyError("IF97Water::saturationTemperature: P at or above "
                            "critical pressure (22.064 MPa); no saturation curve.");
    }
    return IF97::Tsat97(P_Pa);
}

double IF97Water::saturationPressure(double T_K) const {
    checkTemperature(T_K, "saturationPressure");
    if (T_K >= 647.096 && strict_) {
        throw PropertyError("IF97Water::saturationPressure: T at or above "
                            "critical temperature (647.096 K); no saturation curve.");
    }
    return IF97::psat97(T_K);
}

double IF97Water::satLiquidDensity_P(double P_Pa) const {
    checkPressure(P_Pa, "satLiquidDensity_P");
    return IF97::rholiq_p(P_Pa);
}

double IF97Water::satVapourDensity_P(double P_Pa) const {
    checkPressure(P_Pa, "satVapourDensity_P");
    return IF97::rhovap_p(P_Pa);
}

double IF97Water::satLiquidEnthalpy_P(double P_Pa) const {
    checkPressure(P_Pa, "satLiquidEnthalpy_P");
    return IF97::hliq_p(P_Pa);
}

double IF97Water::satVapourEnthalpy_P(double P_Pa) const {
    checkPressure(P_Pa, "satVapourEnthalpy_P");
    return IF97::hvap_p(P_Pa);
}

}  // namespace astara::props
