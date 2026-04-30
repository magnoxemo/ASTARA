#ifndef ASTARA_REACTOR_THERMAL_HPP
#define ASTARA_REACTOR_THERMAL_HPP

/**
 * @file   ReactorThermal.hpp
 * @brief  Lumped-parameter core thermal-hydraulics with Mann nodalisation.
 *
 * This implements the thermal model of the Naghedolfeizi (1990) thesis,
 * Section 3.1 and Figure 3.2: the core is divided into `Nf` fuel nodes and
 * `Nm = 2*Nf` moderator nodes (the moderator runs at half the spatial node
 * size of the fuel, classic Mann nodalisation).  In addition to the in-core
 * nodes, the model carries lower plenum, upper plenum, hot leg, and cold leg
 * lumps, each as a single first-order lag with a fluid residence-time
 * constant.
 *
 * For the thesis defaults (`Nf = 3`, `Nm = 6`) the state vector contains:
 *
 *   - 3 fuel temperatures    T_f1, T_f2, T_f3        [K]
 *   - 6 moderator temperatures T_m1, ..., T_m6       [K]
 *   - 1 cold leg temperature  T_cl                   [K]
 *   - 1 lower plenum          T_lp                   [K]
 *   - 1 upper plenum          T_up                   [K]
 *   - 1 hot leg               T_hl                   [K]
 *
 * Total: **13 thermal states**.  Combined with `PointKinetics` (1 + G
 * states), the full Reactor carries 14..20 states for typical G = 1..6.
 *
 * # Energy balances
 *
 * Fuel node i (fraction F_R of the core power deposited in fuel):
 * @f[
 *   (M_F c_{p,F})\, \frac{dT_{f,i}}{dt}
 *      = F_R\, \frac{P}{N_f} \;-\; h A_i \,(T_{f,i} - \bar T_{m,i})
 * @f]
 * where `P` is the total reactor thermal power, `N_f` is the number of fuel
 * nodes, and `bar T_{m,i}` is the average moderator temperature surrounding
 * fuel node `i`.
 *
 * Moderator node j (Mann's model: heat from fuel + advection from upstream):
 * @f[
 *   (M_C c_{p,C})\, \frac{dT_{m,j}}{dt}
 *      = (1 - F_R)\, \frac{P}{N_m}
 *      \;+\; \frac{h A_j}{2}\,(T_{f,k} - T_{m,j})
 *      \;+\; \frac{2 \dot m c_{p,C}}{N_m}\,(T_{m,j-1} - T_{m,j})
 * @f]
 * (with `T_{m,0} = T_{lp}`, the lower-plenum temperature feeding the bottom
 * of the core).
 *
 * Plenum / leg j (single first-order lag):
 * @f[
 *   \frac{dT_j}{dt} = \frac{T_{j,in} - T_j}{\tau_j}, \qquad
 *   \tau_j = \frac{M_j}{\dot m}.
 * @f]
 *
 * # Outputs
 *
 *   - **Cold-leg-out temperature** (= `T_cl`) is what the SG sees.
 *   - **Average fuel temperature** and **average moderator temperature**
 *     drive the temperature feedback in `Reactor`.
 *
 * @cite Naghedolfeizi (1990), Section 3.1.1, Figure 3.2, eqs. (3.4)-(3.15).
 * @cite Mann, W. R. (1949). "Transient temperature distribution in a
 *       counter-flow heat exchanger." Trans. ASME 71, 663.
 */

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace astara::reactor {

/**
 * @brief Geometry, mass, and heat-transfer parameters of the core thermal model.
 */
struct ReactorThermalParameters {
    // Spatial nodalisation -----------------------------------------------------
    std::size_t num_fuel_nodes      = 3;    ///< N_f, default thesis value
    std::size_t num_moderator_nodes = 6;    ///< N_m = 2*N_f for Mann nodalisation

    // Fuel ---------------------------------------------------------------------
    double fuel_mass_total_kg       = 0.0;  ///< total mass of fuel in the core
    double fuel_cp_J_per_kgK        = 0.0;  ///< specific heat capacity of fuel
    double fission_power_in_fuel    = 0.974;///< F_R, fraction of P deposited in fuel
                                            ///<      (rest goes directly to coolant via gamma)

    // Moderator (coolant in core) ---------------------------------------------
    double moderator_mass_total_kg  = 0.0;  ///< total coolant mass in core (all M nodes)
    double moderator_cp_J_per_kgK   = 0.0;  ///< specific heat capacity of moderator
    double mass_flow_rate_kg_s      = 0.0;  ///< primary coolant mass flow rate

    // Heat transfer -----------------------------------------------------------
    double overall_h_W_per_m2K      = 0.0;  ///< average overall heat-transfer coefficient
    double heat_transfer_area_m2    = 0.0;  ///< total fuel-to-coolant surface area

    // Plenum / piping lumps ---------------------------------------------------
    double lower_plenum_mass_kg     = 0.0;
    double upper_plenum_mass_kg     = 0.0;
    double hot_leg_mass_kg          = 0.0;
    double cold_leg_mass_kg         = 0.0;

    /// Throw `std::invalid_argument` if any required parameter is non-positive
    /// or `num_moderator_nodes != 2 * num_fuel_nodes`.
    void validate() const;
};

/**
 * @brief Snapshot of the thermal-hydraulic state.
 *
 * Sized to the parameters supplied at construction; do not resize at runtime.
 */
struct ReactorThermalState {
    std::vector<double> T_fuel;           ///< fuel-node temperatures, K (size N_f)
    std::vector<double> T_moderator;      ///< moderator-node temperatures, K (size N_m)
    double T_cold_leg     = 0.0;          ///< cold-leg outlet temperature, K
    double T_lower_plenum = 0.0;          ///< lower-plenum temperature,  K
    double T_upper_plenum = 0.0;          ///< upper-plenum temperature,  K
    double T_hot_leg      = 0.0;          ///< hot-leg temperature,       K

    /// Total number of state variables.
    std::size_t size() const noexcept;

    /// Average fuel temperature: (1/N_f) sum_i T_{f,i}.
    double averageFuelTemperature() const noexcept;
    /// Average moderator temperature: (1/N_m) sum_j T_{m,j}.
    double averageModeratorTemperature() const noexcept;

    /// Element-wise vector-space arithmetic, required by the integrator.
    ReactorThermalState& operator+=(const ReactorThermalState& rhs);
    ReactorThermalState& operator-=(const ReactorThermalState& rhs);
    ReactorThermalState& operator*=(double s) noexcept;
};

inline ReactorThermalState operator+(ReactorThermalState a, const ReactorThermalState& b) { a += b; return a; }
inline ReactorThermalState operator-(ReactorThermalState a, const ReactorThermalState& b) { a -= b; return a; }
inline ReactorThermalState operator*(ReactorThermalState a, double s) noexcept           { a *= s; return a; }
inline ReactorThermalState operator*(double s, ReactorThermalState a) noexcept           { a *= s; return a; }

/**
 * @brief Construct an initial `ReactorThermalState` representing steady state
 *        for given inlet temperature and total power.
 *
 * Solves the algebraic system that results from setting all derivatives to
 * zero with the given boundary conditions.  Used to start a transient at a
 * meaningful steady state rather than at zero everywhere.
 *
 * @param  p          thermal parameters
 * @param  P_thermal  total reactor thermal power, W
 * @param  T_inlet_K  cold-leg-in (= primary-side SG outlet) temperature, K
 * @return            steady-state thermal field
 */
ReactorThermalState steadyStateThermal(const ReactorThermalParameters& p,
                                       double P_thermal,
                                       double T_inlet_K);

/**
 * @brief Compute dT/dt for every thermal state variable.
 *
 * The reactor cold-leg-in temperature (the temperature of the fluid returning
 * from the SG / pump) enters as `T_inlet_K`.  The thermal power deposited in
 * the core enters as `P_thermal_W` and is computed by the caller from the
 * point-kinetics state and the rated power.
 *
 * @param  state         current thermal state
 * @param  p             thermal parameters
 * @param  P_thermal_W   instantaneous reactor thermal power, W
 * @param  T_inlet_K     cold-leg-in (lower-plenum-feeding) coolant temperature, K
 * @return               time derivatives, same layout as `state`
 */
ReactorThermalState reactorThermalDerivative(const ReactorThermalState& state,
                                             const ReactorThermalParameters& p,
                                             double P_thermal_W,
                                             double T_inlet_K);

}  // namespace astara::reactor

#endif  // ASTARA_REACTOR_THERMAL_HPP
