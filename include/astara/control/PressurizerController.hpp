#ifndef ASTARA_CONTROL_PRESSURIZER_CONTROLLER_HPP
#define ASTARA_CONTROL_PRESSURIZER_CONTROLLER_HPP

/**
 * @file   PressurizerController.hpp
 * @brief  Pressurizer pressure controller: proportional heaters + spray.
 *
 * Industrial PWR pressurizer control uses two distinct actuators:
 *
 *   - **Proportional heaters** (~ 1.4 MW): always somewhat on at full load
 *     to compensate for ambient losses; ramped up when pressure falls.
 *     Modelled here as proportional to the negative pressure deviation,
 *     bounded by `[Q_min, Q_max]`.
 *
 *   - **Spray valve** (~ 25 kg/s rated): opens only when pressure exceeds
 *     a threshold by more than a small dead-band; provides large cooling
 *     authority via condensation in the steam dome.
 *
 * @cite Naghedolfeizi (1990), Section 3.3.2, Fig. 3.16.
 */

#include <algorithm>

namespace astara::control {

/**
 * @brief Pressurizer pressure controller (proportional heaters + spray).
 */
class PressurizerController {
public:
    struct Config {
        double pressure_setpoint_Pa = 15.5e6;   ///< nominal PWR primary pressure
        // Proportional heaters
        double heater_P_gain_W_per_Pa = 0.5;     ///< W of heater per Pa of *low* deviation
        double heater_min_W           = 0.0;
        double heater_max_W           = 1.4e6;
        double heater_steady_state_W  = 0.2e6;   ///< small bias to balance ambient losses
        // Spray
        double spray_threshold_Pa     = 1.0e5;   ///< P > setpoint + threshold to open spray
        double spray_max_kg_s         = 5.0;     ///< rated spray mass flow
        double spray_gain_kg_per_Pa   = 5.0e-5;  ///< proportional spray opening
        double spray_dead_band_Pa     = 1.0e4;   ///< hysteresis dead-band
    };

    explicit PressurizerController(Config cfg) : cfg_(cfg) {}

    /// Demanded heater power (W) for given measured pressure.
    double heaterDemand(double P_measured_Pa) const noexcept {
        const double err = cfg_.pressure_setpoint_Pa - P_measured_Pa;
        const double Q_p = cfg_.heater_steady_state_W
                         + cfg_.heater_P_gain_W_per_Pa * std::max(err, 0.0);
        return std::clamp(Q_p, cfg_.heater_min_W, cfg_.heater_max_W);
    }

    /// Demanded spray mass flow (kg/s) for given measured pressure.
    /// Latching dead-band: spray turns on at threshold, stays on until
    /// pressure has dropped by `spray_dead_band_Pa` below the threshold.
    double sprayDemand(double P_measured_Pa) noexcept {
        const double dev = P_measured_Pa - cfg_.pressure_setpoint_Pa;
        if (!spray_open_ && dev > cfg_.spray_threshold_Pa) {
            spray_open_ = true;
        } else if (spray_open_
                   && dev < cfg_.spray_threshold_Pa - cfg_.spray_dead_band_Pa) {
            spray_open_ = false;
        }
        if (!spray_open_) return 0.0;
        const double overshoot = std::max(dev - cfg_.spray_threshold_Pa, 0.0);
        return std::min(cfg_.spray_gain_kg_per_Pa * overshoot, cfg_.spray_max_kg_s);
    }

    void reset() noexcept { spray_open_ = false; }

    Config&       config()       noexcept { return cfg_; }
    const Config& config() const noexcept { return cfg_; }
    bool spraysOpen() const noexcept { return spray_open_; }

private:
    Config cfg_;
    bool   spray_open_ = false;
};

}  // namespace astara::control

#endif  // ASTARA_CONTROL_PRESSURIZER_CONTROLLER_HPP
