"""Smoke tests for the astara Python bindings.

One test per bound module: construct at (or near) a physically sensible
operating point, run at least one time_step, and assert the outputs are
finite and in a sane range. These are not physics-validation tests (the C++
analytical/integration test suites already cover that) -- they exist to
catch binding regressions: wrong signatures, broken ownership/keep_alive
wiring, wrong default arguments.
"""

import math
import astara


def _finite(*values):
    return all(math.isfinite(v) for v in values)


def test_core_flow_port():
    port = astara.FlowPort()
    port.mass_flow_kg_s = 4.0

    port.pressure_Pa = 7.0e6
    assert port.mass_flow_kg_s == 4.0
    assert port.pressure_Pa == 7.0e6


def test_props_if97_water():
    props = astara.IF97Water()
    T_sat = props.saturation_temperature(7.0e6)
    assert _finite(T_sat)
    assert 500.0 < T_sat < 600.0
    h = props.sat_vapour_enthalpy_P(7.0e6)
    assert _finite(h) and h > 0.0


def test_props_linearized_water():
    lw = astara.LinearizedWater()
    T_sat = lw.saturation_temperature(7.0e6)
    assert _finite(T_sat) and T_sat > 0.0


def test_reactor():
    groups = astara.DelayedGroupConstants.u235_six_group()
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
    rm.alpha_fuel_per_K = -2.5e-5
    rm.alpha_moderator_per_K = -3.0e-4

    reactor = astara.Reactor(groups, tp, 3.4e9, rm)
    reactor.initialise_steady_state(1.0, 559.0)
    reactor.time_step(1.0e-3)
    p = reactor.state().kinetics.power
    assert _finite(p) and 0.9 < p < 1.1


def test_pump():
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
    pump.time_step(0.01)
    s = pump.state()
    assert _finite(s.volumetric_flow, s.speed_rev_s)
    assert s.speed_rev_s > 0.0


def test_pressurizer():
    props = astara.IF97Water()
    p = astara.PressurizerParameters()
    p.cross_section_area_m2 = 4.0
    p.total_height_m = 13.0
    pz = astara.Pressurizer(p, props)
    pz.initialise_steady_state(15.5e6, 8.0)
    pz.time_step(0.01)
    s = pz.state()
    assert _finite(s.pressure_Pa, s.water_level_m)
    assert s.pressure_Pa > 0.0


def test_ali_steam_generator():
    props = astara.IF97Water()
    p = astara.AliSteamGeneratorParameters.westinghouse_model_D5()
    sg = astara.AliSteamGenerator(p, props)
    sg.initialise_steady_state(597.0, 4400.0, 6.9e6, 3.0)
    sg.time_step(0.01)
    assert _finite(sg.primary_outlet_temperature_K(), sg.state().P)
    assert sg.state().P > 0.0


def test_helical_coil_steam_generator():
    props = astara.IF97Water()
    p = astara.HelicalCoilSteamGeneratorParameters.nuscale_smr_two_sg()
    sg = astara.HelicalCoilSteamGenerator(p, props)
    sg.initialise_steady_state(600.0, 550.0, 67.0, 460.0)
    sg.time_step(0.01)
    assert _finite(sg.state().p_S, sg.state().L_1, sg.state().L_2)


def test_pid_controller():
    cfg = astara.PIDConfig()
    cfg.Kp, cfg.Ki, cfg.Kd, cfg.u_min, cfg.u_max = 2.0, 0.0, 0.0, -10.0, 10.0
    pid = astara.PIDController(cfg)
    pid.set_setpoint(5.0)
    u = pid.update(3.0, 0.1)
    assert abs(u - 4.0) < 1.0e-9


def test_control_rod_system():
    rods = astara.ControlRodSystem(astara.ControlRodSystemConfig())
    assert rods.reactivity_rate(astara.RodMode.Insert) < 0.0
    assert rods.reactivity_rate(astara.RodMode.Withdraw) > 0.0
    assert rods.reactivity_rate(astara.RodMode.Hold) == 0.0


def test_primary_loop_integrates_stably():
    groups = astara.DelayedGroupConstants.u235_six_group()
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
    rm.alpha_fuel_per_K = -2.5e-5
    rm.alpha_moderator_per_K = -3.0e-4
    reactor = astara.Reactor(groups, tp, 3.4e9, rm)
    reactor.initialise_steady_state(1.0, 559.0)

    props = astara.IF97Water()
    sgp = astara.AliSteamGeneratorParameters.westinghouse_model_D5()
    sg = astara.AliSteamGenerator(sgp, props)
    sg.initialise_steady_state(597.0, 4400.0, 6.9e6, 3.0)

    pumpp = astara.HomologousPumpParameters()
    pumpp.curve.A0, pumpp.curve.A1, pumpp.curve.A2 = 90.0, 1.0, -1.0
    pumpp.rated_speed_rev_s = 20.0
    pumpp.rated_volumetric_flow_m3_s = 6.0
    pumpp.loop_resistance_K_s2_m5 = 60.0 / 36.0
    pumpp.effective_flow_area_m2 = 0.4
    pumpp.loop_length_m = 70.0
    pumpp.fluid_density_kg_m3 = 720.0
    pumpp.moment_of_inertia_kg_m2 = 1500.0
    pumpp.rated_input_power_W = 720.0 * 9.80665 * 6.0 * 60.0
    pump = astara.HomologousPump(pumpp)
    pump.initialise_at_rated()

    ppz = astara.PressurizerParameters()
    ppz.cross_section_area_m2 = 4.0
    ppz.total_height_m = 13.0
    pz = astara.Pressurizer(ppz, props)
    pz.initialise_steady_state(15.5e6, 8.0)

    loop = astara.PrimaryLoop(reactor, sg, pump, pz)
    dt = 1.0e-3
    for _ in range(2000):  # 2 s, enough to catch a broken binding without a slow test
        loop.time_step(dt)

    p = loop.reactor().state().kinetics.power
    assert _finite(p) and 0.5 < p < 2.0
