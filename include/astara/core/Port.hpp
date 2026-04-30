#ifndef ASTARA_CORE_PORT_HPP
#define ASTARA_CORE_PORT_HPP

/**
 * @file   Port.hpp
 * @brief  Boundary-condition struct for inter-component coupling.
 *
 * In a primary-loop transient, components pass each other a small bundle of
 * thermohydraulic data at their interfaces:
 *
 *   - mass flow rate (kg/s)
 *   - specific enthalpy (J/kg)         OR temperature, depending on context
 *   - pressure (Pa)
 *
 * `FlowPort` is the canonical struct.  When wiring components in the primary
 * loop, the integrator copies the *output* port of one component into the
 * *input* port of the next, after each RK4 stage.  This keeps the coupling
 * explicit and inspectable in tests.
 *
 * Direction convention: a `FlowPort` always describes flow *out of* a
 * component.  If `mass_flow > 0`, fluid moves out of the component; this is
 * the normal forward direction in steady operation.  Reverse flow uses
 * `mass_flow < 0` and is supported but rare.
 */

namespace astara::core {

/**
 * @brief Thermohydraulic boundary condition at a component port.
 *
 * Used to wire reactor-core outlet -> hot leg -> SG primary inlet, etc.
 */
struct FlowPort {
    double mass_flow_kg_s = 0.0;   ///< mass flow rate, kg/s, positive in forward flow
    double enthalpy_J_kg  = 0.0;   ///< specific enthalpy of the stream, J/kg
    double pressure_Pa    = 0.0;   ///< stream pressure, Pa
    double temperature_K  = 0.0;   ///< stream temperature, K (informational; redundant with (P,h))
};

}  // namespace astara::core

#endif  // ASTARA_CORE_PORT_HPP
