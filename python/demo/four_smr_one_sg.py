#!/usr/bin/env python3
"""Four SMR reactor+pump loops feeding one shared steam generator.

The motivating example for `astara.Plant`/`astara.Loop`: the compiled
`astara.PrimaryLoop` couples exactly one reactor to one SG to one pump to
one pressurizer, so it can't represent a small-modular-reactor plant where
several reactor+pump loops share a single steam generator. `Plant` builds
that topology instead, out of the same bound `Reactor`/`HomologousPump`/
`AliSteamGenerator` classes `PrimaryLoop` itself uses -- see
`python/astara/network.py` for the flow-weighted mixing rule at the shared
SG's inlet (fan-in) and how its single outlet temperature is broadcast back
to every reactor's cold leg (fan-out).

At t = 10 s, loop 0's control rods are driven in (`astara.Rod`) for 5 s and
then held, dropping its power; at t = 20 s loop 1's rods do the same.
`astara.Simulation` runs and records the transient; `astara.postprocess`
plots it -- this is the "import astara in a notebook, run, post-process"
workflow end to end.

Meant to be run as a plain script (`python3 demo/four_smr_one_sg.py`, writes
`four_smr_one_sg.png`) or pasted cell-by-cell into a Jupyter notebook, where
`sim.run(...)` and `postprocess.plot_traces(...)` display inline instead.

Requires the `astara` extension to be importable and `astara[analysis]`
(pandas + matplotlib):
    PYTHONPATH=<build>/python python3 demo/four_smr_one_sg.py
"""

import astara

N_LOOPS = 4


def make_reactor():

    tp = astara.ReactorThermalParameters()
    tp.num_fuel_nodes = 20
    tp.num_moderator_nodes = 40  # fixed 2x coupling: 2 moderator nodes bridge each fuel node
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
    rm.alpha_fuel_per_K = -2.5e-5
    rm.alpha_moderator_per_K = -3.0e-4

    reactor = astara.Reactor(astara.DelayedGroupConstants.u235_six_group(), tp, 3.4e9 / N_LOOPS, rm)
    reactor.initialise_steady_state(n0=1.0, T_inlet_K=559.0)
    return reactor


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


def make_shared_sg(props):
    p = astara.AliSteamGeneratorParameters.westinghouse_model_D5()
    sg = astara.AliSteamGenerator(p, props)
    sg.initialise_steady_state(597.0, 4400.0, 6.9e6, 3.0)
    return sg


def animate_axial_temperatures(df, n_loops, out_path, target_frames=150, fps=20):
    """GIF of every loop's fuel + moderator axial temperature profile over time.

    `df` must carry `fuel_profile{i}`/`moderator_profile{i}` columns (each
    cell a list of per-node temperatures, one column per loop -- see
    `loop_fuel_profile`/`loop_moderator_profile` in `main()`) plus `t_s`.
    Subsamples to ~`target_frames` frames so the file stays a reasonable
    size regardless of how finely the simulation was recorded.
    """
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation

    stride = max(1, len(df) // target_frames)
    sub = df.iloc[::stride].reset_index(drop=True)

    fuel_cols = [f"fuel_profile{i}" for i in range(n_loops)]
    mod_cols = [f"moderator_profile{i}" for i in range(n_loops)]
    x_fuel = list(range(1, len(sub[fuel_cols[0]].iloc[0]) + 1))
    x_mod = list(range(1, len(sub[mod_cols[0]].iloc[0]) + 1))

    def bounds(cols):
        lo = min(min(row) for col in cols for row in sub[col])
        hi = max(max(row) for col in cols for row in sub[col])
        pad = 0.05 * (hi - lo) if hi > lo else 1.0
        return lo - pad, hi + pad

    fuel_lo, fuel_hi = bounds(fuel_cols)
    mod_lo, mod_hi = bounds(mod_cols)

    fig, (ax_fuel, ax_mod) = plt.subplots(1, 2, figsize=(11, 4.5))
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    fuel_lines = [ax_fuel.plot([], [], color=colors[i % len(colors)], label=f"loop{i}")[0] for i in range(n_loops)]
    mod_lines = [ax_mod.plot([], [], color=colors[i % len(colors)], label=f"loop{i}")[0] for i in range(n_loops)]

    for ax, x, ylo, yhi, ylabel, title in (
        (ax_fuel, x_fuel, fuel_lo, fuel_hi, "fuel temperature [K]", "Fuel axial profile"),
        (ax_mod, x_mod, mod_lo, mod_hi, "coolant/moderator temperature [K]", "Coolant axial profile"),
    ):
        ax.set_xlim(x[0], x[-1])
        ax.set_ylim(ylo, yhi)
        ax.set_xlabel("axial node (inlet -> outlet)")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend(loc="upper left", fontsize="small")

    time_title = fig.suptitle("")

    def update(frame_idx):
        row = sub.iloc[frame_idx]
        for i in range(n_loops):
            fuel_lines[i].set_data(x_fuel, row[f"fuel_profile{i}"])
            mod_lines[i].set_data(x_mod, row[f"moderator_profile{i}"])
        time_title.set_text(f"t = {row['t_s']:.1f} s")
        return [*fuel_lines, *mod_lines, time_title]

    anim = animation.FuncAnimation(fig, update, frames=len(sub), interval=1000.0 / fps, blit=False)
    anim.save(out_path, writer=animation.PillowWriter(fps=fps))
    plt.close(fig)


def main():
    props = astara.IF97Water()

    plant = astara.Plant()
    loops = [plant.add(f"loop{i}", astara.Loop(make_reactor(), make_pump())) for i in range(N_LOOPS)]
    sg = plant.add("sg", make_shared_sg(props))

    plant.connect(loops, sg)  # fan-in:  4 loops -> shared SG inlet
    plant.connect(sg, loops)  # fan-out: shared SG outlet -> every loop's cold leg

    rod1 = astara.Rod(loops[0].reactor)
    rod2 = astara.Rod(loops[1].reactor)

    plant.add("rod0", rod1)
    plant.add("rod1", rod2)

    sim = astara.Simulation(plant, dt=1.0e-3)

    # Initial power increase
    sim.schedule(t=5.0,   action=rod1.insert)
    sim.schedule(t=8.0,   action=rod1.hold)

    sim.schedule(t=15.0,  action=rod2.insert)
    sim.schedule(t=18.0,  action=rod2.hold)

    # Partial withdrawal / power increase
    sim.schedule(t=30.0,  action=rod1.withdraw)
    sim.schedule(t=35.0,  action=rod1.hold)

    # Second rod movement
    sim.schedule(t=45.0,  action=rod2.withdraw)
    sim.schedule(t=50.0,  action=rod2.hold)

    # Larger transient
    sim.schedule(t=65.0,  action=rod1.insert)
    sim.schedule(t=70.0,  action=rod2.insert)

    # Hold both rods
    sim.schedule(t=80.0,  action=rod1.hold)
    sim.schedule(t=82.0,  action=rod2.hold)

    # Oscillatory control sequence
    sim.schedule(t=95.0,  action=rod1.withdraw)
    sim.schedule(t=100.0, action=rod1.insert)

    sim.schedule(t=105.0, action=rod2.withdraw)
    sim.schedule(t=110.0, action=rod2.insert)

    # Long-term withdrawal
    sim.schedule(t=120.0, action=rod1.withdraw)
    sim.schedule(t=125.0, action=rod1.hold)

    sim.schedule(t=135.0, action=rod2.withdraw)
    sim.schedule(t=140.0, action=rod2.hold)

    # Final reactivity excursion
    sim.schedule(t=155.0, action=rod1.insert)
    sim.schedule(t=160.0, action=rod2.insert)

    sim.schedule(t=170.0, action=rod1.hold)
    sim.schedule(t=172.0, action=rod2.hold)

    def loop_power_pct(loop):
        return lambda: loop.reactor.state().kinetics.power * 100.0

    def loop_temp_K(loop):
        return lambda: loop.get_outlet()["temperature_K"]

    def loop_fuel_temp_K(loop):
        return lambda: loop.reactor.state().thermal.average_fuel_temperature()

    def loop_reactivity_pcm(loop):
        def get():
            state = loop.reactor.state()
            Tf = state.thermal.average_fuel_temperature()
            Tm = state.thermal.average_moderator_temperature()
            return loop.reactor.reactivity().evaluate(Tf, Tm) / 1.0e-5
        return get

    def loop_fuel_profile(loop):
        return lambda: list(loop.reactor.state().thermal.T_fuel)

    def loop_moderator_profile(loop):
        return lambda: list(loop.reactor.state().thermal.T_moderator)

    sim.record(
        {
            "sg_primary_outlet_T_K": sg.primary_outlet_temperature_K,
            **{f"power{i}_pct": loop_power_pct(loop) for i, loop in enumerate(loops)},
            **{f"temp{i}_K": loop_temp_K(loop) for i, loop in enumerate(loops)},
            **{f"fuel_temp{i}_K": loop_fuel_temp_K(loop) for i, loop in enumerate(loops)},
            **{f"reactivity{i}_pcm": loop_reactivity_pcm(loop) for i, loop in enumerate(loops)},
            **{f"fuel_profile{i}": loop_fuel_profile(loop) for i, loop in enumerate(loops)},
            **{f"moderator_profile{i}": loop_moderator_profile(loop) for i, loop in enumerate(loops)},
            "rod0_position": lambda: rod1.position,
            "rod1_position": lambda: rod2.position,
        },
        every=0.1,
    )
    df = sim.run(t_end=300.0)

    scalar_cols = [c for c in df.columns if not c.startswith(("fuel_profile", "moderator_profile"))]
    print("Four-SMR-into-one-SG demo: rods driven into loop 0 at t=10-15 s, loop 1 at t=20-25 s\n")
    print(df[scalar_cols].iloc[[0, -1]].to_string(index=False))

    ax = astara.postprocess.plot_traces(
        df, [f"power{i}_pct" for i in range(N_LOOPS)],
        title="Reactor power per loop (rods inserted into loop 0 at t=10-15 s, loop 1 at t=20-25 s)",
    )
    ax.set_ylabel("power [% rated]")
    ax.figure.savefig("four_smr_one_sg.png", dpi=100)
    print("\nPlot written to four_smr_one_sg.png")

    ax1 = astara.postprocess.plot_traces(df, [f"temp{i}_K" for i in range(N_LOOPS)])

    ax1.set_ylabel("outlet temperature [K]")
    ax1.figure.savefig("four_smr_one_sg_temp.png", dpi=100)
    print("\nPlot written to four_smr_one_sg_temp.png")

    ax2 = astara.postprocess.plot_traces(df, [f"fuel_temp{i}_K" for i in range(N_LOOPS)])
    ax2.set_ylabel("average fuel temperature [K]")
    ax2.figure.savefig("four_smr_one_sg_fuel_temp.png", dpi=100)
    print("\nPlot written to four_smr_one_sg_fuel_temp.png")

    ax3 = astara.postprocess.plot_traces(df, [f"reactivity{i}_pcm" for i in range(N_LOOPS)])
    ax3.set_ylabel("total reactivity [pcm]")
    ax3.figure.savefig("four_smr_one_sg_reactivity.png", dpi=100)
    print("\nPlot written to four_smr_one_sg_reactivity.png")

    animate_axial_temperatures(df, N_LOOPS, "four_smr_one_sg_axial_temp.gif")
    print("\nAnimation written to four_smr_one_sg_axial_temp.gif")
    #
    #
    # ax4 = astara.postprocess.plot_traces(df, [f"mass_flow_kg{i}" for i in range(N_LOOPS)])
    #
    # ax4.set_ylabel("outlet temperature [K]")
    # ax4.figure.savefig("four_smr_ms_f.png", dpi=100)
    # print("\nPlot written to four_smr_one_sg_temp.png")
    # return 0


if __name__ == "__main__":
    raise SystemExit(main())
