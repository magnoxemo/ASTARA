#ifndef ASTARA_PRIMARY_PRIMARY_LOOP_HPP
#define ASTARA_PRIMARY_PRIMARY_LOOP_HPP

/**
 * @file   PrimaryLoop.hpp
 * @brief  Integrated PWR primary loop: reactor + SG + pump + pressurizer
 *         with optional controllers attached.
 *
 * Wires the four major component classes together so that one loop
 * integration call advances them all consistently:
 *
 *   - reactor.hot_leg  -> SG.primary_inlet
 *   - SG.primary_outlet -> pump (cold leg path)
 *   - pump.massFlow    -> reactor.mass_flow_rate (set in parameters; updated
 *                          dynamically when called via the loop)
 *   - reactor.cold_leg_in <- SG.primary_outlet (after lag through pump)
 *
 * The pressurizer is connected to the cold leg via the surge line; in the
 * present implementation we treat surge enthalpy as the cold-leg saturation
 * enthalpy, neglecting flashing dynamics in the surge line itself.
 *
 * # Design choice: explicit operator coupling
 *
 * Rather than concatenating all four component states into one giant vector
 * and integrating the lot, we step each component with its own RK4 stage
 * and exchange port values **between** time steps.  This is sometimes
 * called "operator splitting" and is justified because:
 *
 *   - The fastest dynamics (reactor prompt-neutron, ~ 30-300 ms) are
 *     intrinsic to the reactor and its internal RK4 already resolves them.
 *   - The slowest dynamics (drum level, primary T_avg) have time constants
 *     of tens of seconds, much longer than the loop time step.
 *   - The errors introduced by the splitting are O(dt) on a per-stage
 *     basis but bounded by the natural smoothness of the transients.
 *
 * For ultra-fast (< 100 ms) coupled phenomena (e.g. main-steam-line break)
 * a fully-coupled monolithic integrator would be more accurate; for the
 * load-following / step-reactivity transients in scope here, operator
 * splitting is sufficient and far simpler to wire.
 *
 * Optional **OpenMP parallelism**: when `evaluateAllDerivatives()` is
 * called the four component derivatives are evaluated in parallel using
 * `#pragma omp parallel sections`.  This roughly halves the per-step cost
 * on a 4-core machine for the reactor+SG transients (the pump and
 * pressurizer each take a small fraction of the reactor cost).
 */

#include "astara/reactor/Reactor.hpp"
#include "astara/sg/AliSteamGenerator.hpp"
#include "astara/pump/HomologousPump.hpp"
#include "astara/pressurizer/Pressurizer.hpp"
#include "astara/control/ReactorController.hpp"
#include "astara/control/PressurizerController.hpp"
#include "astara/control/ThreeElementController.hpp"

#include <memory>
#include <stdexcept>

namespace astara::primary {

/**
 * @brief Integrated PWR primary-loop simulator.
 *
 * Owns one `Reactor`, one `AliSteamGenerator`, one `HomologousPump`, and
 * one `Pressurizer`.  Optional controllers may be attached; if none are
 * attached, the components run open-loop with whatever inputs the user
 * has set on them.
 *
 * @note This class is **not** thread-safe (single owner, sequential
 *       advancement).  OpenMP parallelism is used internally inside
 *       `evaluateAllDerivatives`; the public `timeStep` is a serial call.
 */
class PrimaryLoop {
public:
    /**
     * @brief Construct from existing components.  Components are passed
     *        by `unique_ptr` so the loop assumes ownership and lifetime.
     *
     * @param reactor      already-initialised reactor
     * @param sg           already-initialised SG (with same WaterProperties
     *                     as the pressurizer; that constraint is checked
     *                     at construction time only nominally -- the user
     *                     is responsible for consistency)
     * @param pump         already-initialised pump
     * @param pressurizer  already-initialised pressurizer
     */
    PrimaryLoop(std::unique_ptr<reactor::Reactor>       reactor,
                std::unique_ptr<sg::AliSteamGenerator>  sg,
                std::unique_ptr<pump::HomologousPump>   pump,
                std::unique_ptr<pressurizer::Pressurizer> pressurizer);

    // ---------- Component accessors ----------
    reactor::Reactor&             reactor()      noexcept { return *reactor_; }
    const reactor::Reactor&       reactor() const noexcept { return *reactor_; }
    sg::AliSteamGenerator&        steamGenerator()      noexcept { return *sg_; }
    const sg::AliSteamGenerator&  steamGenerator() const noexcept { return *sg_; }
    pump::HomologousPump&         pump()         noexcept { return *pump_; }
    const pump::HomologousPump&   pump() const   noexcept { return *pump_; }
    pressurizer::Pressurizer&     pressurizer()       noexcept { return *pressurizer_; }
    const pressurizer::Pressurizer& pressurizer() const noexcept { return *pressurizer_; }

    // ---------- Controllers (optional, ownership transferred) ----------
    void setReactorController(std::unique_ptr<control::ReactorController> c) noexcept {
        reactor_controller_ = std::move(c);
    }
    void setPressurizerController(std::unique_ptr<control::PressurizerController> c) noexcept {
        pressurizer_controller_ = std::move(c);
    }
    void setFeedwaterController(std::unique_ptr<control::ThreeElementController> c) noexcept {
        feedwater_controller_ = std::move(c);
    }
    bool hasReactorController()      const noexcept { return reactor_controller_      != nullptr; }
    bool hasPressurizerController()  const noexcept { return pressurizer_controller_  != nullptr; }
    bool hasFeedwaterController()    const noexcept { return feedwater_controller_    != nullptr; }

    /**
     * @brief Advance the entire loop by `dt`.
     *
     * Implementation:
     *   1. Update mass flow rate from the pump and propagate to reactor + SG.
     *   2. Run controllers if attached.
     *   3. Step each component independently with its own RK4 step.
     *   4. Exchange port values: reactor.hot -> SG.primary_in; SG.primary_out
     *      -> reactor.cold_leg_in; surge enthalpy from cold leg to pressurizer.
     */
    void timeStep(double dt);

    /// Current simulation time.  Synchronised across components after each step.
    double timeSeconds() const noexcept { return reactor_->state().t_s; }

    /// Sanity check at construction or after parameter change: verifies that
    /// the four components carry consistent state (e.g. SG primary inlet
    /// matches reactor hot-leg-out).  Useful in tests.
    bool isConsistent(double tol_K = 5.0) const noexcept;

private:
    std::unique_ptr<reactor::Reactor>          reactor_;
    std::unique_ptr<sg::AliSteamGenerator>     sg_;
    std::unique_ptr<pump::HomologousPump>      pump_;
    std::unique_ptr<pressurizer::Pressurizer>  pressurizer_;

    std::unique_ptr<control::ReactorController>      reactor_controller_;
    std::unique_ptr<control::PressurizerController>  pressurizer_controller_;
    std::unique_ptr<control::ThreeElementController> feedwater_controller_;
};

}  // namespace astara::primary

#endif  // ASTARA_PRIMARY_PRIMARY_LOOP_HPP
