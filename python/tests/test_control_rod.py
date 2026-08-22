"""Tests for the `astara.Rod` control-rod wrapper.

`ControlRodSystem.reactivity_rate()` itself is exercised by
`tests/unit/test_control_rod_system.cpp` (sign convention) and
`test_smoke.py::test_control_rod_system`. This checks the layer `Rod` adds:
automatic integration of that rate into the bound reactor's
`reactivity().rho_external`, gated by insert/withdraw/hold, and clamped to
the rod's `position` travel range.
"""

import warnings

import astara


def _make_reactor():
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
    reactor = astara.Reactor(astara.DelayedGroupConstants.u235_six_group(), tp, 3.4e9, rm)
    reactor.initialise_steady_state(n0=1.0, T_inlet_K=559.0)
    return reactor


def test_insert_drives_reactivity_negative():
    reactor = _make_reactor()
    rod = astara.Rod(reactor)
    rod.insert()
    rho0 = reactor.reactivity().rho_external
    for _ in range(100):
        rod.time_step(1.0e-3)
    assert reactor.reactivity().rho_external < rho0


def test_withdraw_drives_reactivity_positive():
    reactor = _make_reactor()
    rod = astara.Rod(reactor)
    rod.withdraw()
    rho0 = reactor.reactivity().rho_external
    for _ in range(100):
        rod.time_step(1.0e-3)
    assert reactor.reactivity().rho_external > rho0


def test_hold_leaves_reactivity_unchanged():
    reactor = _make_reactor()
    rod = astara.Rod(reactor)
    rod.hold()
    rho0 = reactor.reactivity().rho_external
    for _ in range(100):
        rod.time_step(1.0e-3)
    assert reactor.reactivity().rho_external == rho0


def test_mode_switch_stops_the_drift():
    reactor = _make_reactor()
    rod = astara.Rod(reactor)
    rod.insert()
    for _ in range(50):
        rod.time_step(1.0e-3)
    rod.hold()
    rho_after_insert = reactor.reactivity().rho_external
    for _ in range(50):
        rod.time_step(1.0e-3)
    assert reactor.reactivity().rho_external == rho_after_insert


def test_position_starts_at_midpoint():
    rod = astara.Rod(_make_reactor())
    assert rod.position == 0.5


def test_insert_clamps_position_at_zero_and_warns_once():
    reactor = _make_reactor()
    rod = astara.Rod(reactor)
    rod.insert()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        for _ in range(2000):  # far more than enough to bottom out from 0.5
            rod.time_step(1.0e-3)
    assert rod.position == 0.0
    limit_warnings = [w for w in caught if issubclass(w.category, RuntimeWarning)]
    assert len(limit_warnings) == 1
    assert "fully inserted" in str(limit_warnings[0].message)


def test_reactivity_stops_changing_once_fully_inserted():
    reactor = _make_reactor()
    rod = astara.Rod(reactor)
    rod.insert()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for _ in range(1000):
            rod.time_step(1.0e-3)
        rho_at_limit = reactor.reactivity().rho_external
        for _ in range(1000):
            rod.time_step(1.0e-3)
    assert reactor.reactivity().rho_external == rho_at_limit


def test_withdraw_clamps_position_at_one_and_warns_once():
    reactor = _make_reactor()
    rod = astara.Rod(reactor)
    rod.withdraw()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        for _ in range(2000):
            rod.time_step(1.0e-3)
    assert rod.position == 1.0
    limit_warnings = [w for w in caught if issubclass(w.category, RuntimeWarning)]
    assert len(limit_warnings) == 1
    assert "fully withdrawn" in str(limit_warnings[0].message)
