#!/usr/bin/env python3
"""Coupled primary-loop demonstration: reactor coolant pump trip (Python API).

Builds two identical PWR primary loops (reactor + Ali Model D U-tube SG +
homologous pump + pressurizer), started from the same initial conditions.
At t = 20 s, one of them ("tripped") has its reactor coolant pump's input
shaft power set to zero; the other ("baseline") keeps running undisturbed.

Why compare two runs instead of one: these components are each initialised
at their own individually-consistent steady state, but the four together are
only *approximately* a joint equilibrium -- `PrimaryLoop::timeStep` only
propagates the pump's mass flow into the SG's primary-side input (see
`include/astara/primary/PrimaryLoop.hpp`), not into the reactor's own core
thermal-hydraulics parameters, so even an untouched loop drifts slowly for
tens of seconds while the four components settle against each other. A
single tripped run can't distinguish "caused by the trip" from "the model
settling on its own" -- the baseline run subtracts that out, via
`astara.postprocess.compare()`.

Uses `astara.Simulation` for the run loop/scheduling/recording and
`astara.postprocess` for the isolated-trip-effect comparison and CSV export,
rather than the hand-rolled stepping loop this demo used to write.

Writes both traces to `pump_trip_primary_loop.csv` and prints the trip's
isolated effect (tripped - baseline) at the end of the run.

Usage:
    PYTHONPATH=<build>/python python3 demo/pump_trip_primary_loop.py
Requires `astara[analysis]` (pandas) for the DataFrame recording/CSV export.
"""

import astara


def make_reactor():
    tp = astara.ReactorThermalParameters()
    tp.num_fuel_nodes = 3
    tp.num_moderator_nodes = 6
    tp.fuel_mass_total_kg = 1.0e5
    tp.fuel_cp_J_per_kgK = 300.0
    tp.fission_power_in_fuel = 0.974
    tp.moderator_mass_total_kg = 1.5e4
    tp.moderator_cp_J_per_kgK = 5400.0
    tp.mass_flow_rate_kg_s = 17600.0
    tp.overall_h_W_per_m2K = 28000.0
    tp.heat_transfer_area_m2 = 5400.0
    tp.lower_plenum_mass_kg = 5000.0
    tp.upper_plenum_mass_kg = 5000.0
    tp.hot_leg_mass_kg = 2500.0
    tp.cold_leg_mass_kg = 2500.0

    rm = astara.ReactivityModel()
    rm.alpha_fuel_per_K = -2.5e-5       # Doppler
    rm.alpha_moderator_per_K = -3.0e-4  # moderator (negative)

    reactor = astara.Reactor(astara.DelayedGroupConstants.u235_six_group(), tp, 3.4e9, rm)
    reactor.initialise_steady_state(n0=1.0, T_inlet_K=559.0)
    return reactor


def make_sg(props):
    p = astara.AliSteamGeneratorParameters.westinghouse_model_D5()
    sg = astara.AliSteamGenerator(p, props)
    sg.initialise_steady_state(597.0, 4400.0, 6.9e6, 3.0)
    return sg


def make_pump():
    p = astara.HomologousPumpParameters()
    p.curve.A0, p.curve.A1, p.curve.A2 = 90.0, 1.0, -1.0
    p.rated_speed_rev_s = 20.0
    p.rated_volumetric_flow_m3_s = 6.0
    p.loop_resistance_K_s2_m5 = 60.0 / 36.0
    p.effective_flow_area_m2 = 0.4
    p.loop_length_m = 70.0
    p.fluid_density_kg_m3 = 720.0
    p.moment_of_inertia_kg_m2 = 1500.0
    p.rated_input_power_W = 720.0 * 9.80665 * 6.0 * 60.0
    pump = astara.HomologousPump(p)
    pump.initialise_at_rated()
    return pump


def make_pressurizer(props):
    geom = astara.PressurizerParameters()
    geom.cross_section_area_m2 = 4.0
    geom.total_height_m = 13.0
    pz = astara.Pressurizer(geom, props)
    pz.initialise_steady_state(15.5e6, 8.0)
    return pz


def make_loop(props):
    return astara.PrimaryLoop(make_reactor(), make_sg(props), make_pump(), make_pressurizer(props))


def main():
    print("ASTARA primary-loop demo: reactor coolant pump trip")
    print("  (no controllers, no scram; compared against an untripped baseline)\n")

    # One shared, stateless IF97Water service for every component in both loops.
    props = astara.IF97Water()
    baseline_loop = make_loop(props)
    tripped_loop = make_loop(props)

    def make_sim(loop, trip_at_s):
        sim = astara.Simulation(loop, dt=1.0e-3)
        if trip_at_s is not None:
            sim.schedule(t=trip_at_s, action=lambda: loop.pump().set_input_power_W(0.0))
        sim.record(
            lambda: {
                "power_pct": loop.reactor().state().kinetics.power * 100.0,
                "pump_speed_rev_s": loop.pump().state().speed_rev_s,
                "T_hot_K": loop.reactor().state().thermal.T_hot_leg,
            },
            every=0.1,
        )
        return sim

    print("  t = 20.0 s: PUMP TRIP (tripped loop only)")
    baseline_df = make_sim(baseline_loop, trip_at_s=None).run(t_end=60.0)
    tripped_df = make_sim(tripped_loop, trip_at_s=20.0).run(t_end=60.0)

    filename = "pump_trip_primary_loop.csv"
    astara.postprocess.to_csv(tripped_df.merge(baseline_df, on="t_s", suffixes=("_tripped", "_baseline")), filename)

    effect = astara.postprocess.compare(
        baseline_df, tripped_df, columns=["power_pct", "pump_speed_rev_s", "T_hot_K"]
    )

    pb = baseline_df["power_pct"].iloc[-1]
    pt = tripped_df["power_pct"].iloc[-1]
    Tb = baseline_df["T_hot_K"].iloc[-1]
    Tt = tripped_df["T_hot_K"].iloc[-1]

    print("\nFinal state at t = 60.0 s:")
    print(f"  Reactor power   : baseline {pb:.2f} %   tripped {pt:.2f} %   "
          f"(trip effect: {effect['power_pct'].iloc[-1]:+.2f} pp)")
    print(f"  Pump speed      : baseline {baseline_df['pump_speed_rev_s'].iloc[-1]:.2f} rev/s   "
          f"tripped {tripped_df['pump_speed_rev_s'].iloc[-1]:.2f} rev/s")
    print(f"  Primary T_hot   : baseline {Tb:.2f} K   tripped {Tt:.2f} K   "
          f"(trip effect: {effect['T_hot_K'].iloc[-1]:+.2f} K)")
    print(f"\nTrace: {filename}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
