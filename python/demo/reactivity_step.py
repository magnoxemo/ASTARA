#!/usr/bin/env python3
"""Step-reactivity transient demonstration (Python API).

Initialises a Westinghouse-like 4-loop PWR at full power and applies a step
in external reactivity, then integrates for 600 s. Writes the time trace of
normalised power, average fuel/moderator temperatures, and total reactivity
(external + feedback) to `reactivity_step_transient.csv`.

Python re-implementation of ``examples/reactivity_step_transient.cpp`` --
same model, same defaults, same output columns -- to show the C++ and Python
APIs are interchangeable for this kind of transient study. Uses
`astara.Simulation` for the run loop/recording and `astara.postprocess` for
the CSV write, rather than a hand-rolled logging loop.

Usage:
    python3 reactivity_step.py            # default +50 pcm step
    python3 reactivity_step.py -100       # -100 pcm step

Requires the `astara` extension to be importable, e.g.:
    PYTHONPATH=<build>/python python3 demo/reactivity_step.py
Requires `astara[analysis]` (pandas) for the CSV export.
"""

import sys

import astara


def make_westinghouse_like_reactor():
    tp = astara.ReactorThermalParameters()
    tp.num_fuel_nodes = 10
    tp.num_moderator_nodes = 2 * tp.num_fuel_nodes
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

    groups = astara.DelayedGroupConstants.u235_six_group()
    P_rated = 3.4e9  # 3400 MWth

    rm = astara.ReactivityModel()
    rm.alpha_fuel_per_K = -2.5e-5       # Doppler
    rm.alpha_moderator_per_K = -3.0e-4  # moderator (negative)
    return astara.Reactor(groups, tp, P_rated, rm)


def main():
    rho_pcm = 50.0
    if len(sys.argv) >= 2:
        try:
            rho_pcm = float(sys.argv[1])
        except ValueError:
            print(f"Failed to parse reactivity from '{sys.argv[1]}'", file=sys.stderr)
            return 2
    rho_step = rho_pcm * 1.0e-5

    reactor = make_westinghouse_like_reactor()
    reactor.initialise_steady_state(n0=1.0, T_inlet_K=559.0)
    reactor.reactivity().rho_external = rho_step

    def sample():
        s = reactor.state()
        Tf = s.thermal.average_fuel_temperature()
        Tm = s.thermal.average_moderator_temperature()
        return {
            "n_over_n0": s.kinetics.power,
            "T_fuel_avg_K": Tf,
            "T_mod_avg_K": Tm,
            "rho_total_pcm": reactor.reactivity().evaluate(Tf, Tm) / 1e-5,
            "T_hot_K": s.thermal.T_hot_leg,
            "T_cold_K": s.thermal.T_cold_leg,
        }

    sim = astara.Simulation(reactor, dt=1.0e-3)
    sim.record(sample, every=0.05)  # every 50 ms, matching the C++ example
    df = sim.run(t_end=600.0)

    filename = "reactivity_step_transient.csv"
    astara.postprocess.to_csv(df, filename)

    print(f"Step reactivity = {rho_pcm} pcm; integrated 600 s.")
    print(f"Final n/n0 = {reactor.state().kinetics.power}")
    print(f"Final T_hot_leg = {reactor.state().thermal.T_hot_leg} K")
    print(f"CSV trace written to {filename}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
