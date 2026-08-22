"""Tests for the pure-Python `astara.network` plant/topology layer.

These don't re-validate reactor/SG physics (the C++ analytical test suites
already do that) -- they check that `Plant.connect()`/`step()` correctly
implements the documented flow-weighted mixing rule and moves values
between the right accessors.
"""

import astara


class _FakeSource:
    """Minimal stand-in with a fixed outlet, registered via register_adapter."""

    def __init__(self, temperature_K, mass_flow_kg_s):
        self.temperature_K = temperature_K
        self.mass_flow_kg_s = mass_flow_kg_s

    def time_step(self, dt):
        pass


class _FakeSink:
    def __init__(self):
        self.last_inlet = None

    def time_step(self, dt):
        pass


def _register_fakes():
    astara.network.register_adapter(
        "_FakeSource",
        lambda src: {"temperature_K": src.temperature_K, "mass_flow_kg_s": src.mass_flow_kg_s},
        lambda src, **fields: None,
    )
    astara.network.register_adapter(
        "_FakeSink",
        lambda sink: {},
        lambda sink, **fields: setattr(sink, "last_inlet", fields),
    )


def test_fan_in_mixing_is_flow_weighted():
    _register_fakes()
    a = _FakeSource(temperature_K=600.0, mass_flow_kg_s=10.0)
    b = _FakeSource(temperature_K=500.0, mass_flow_kg_s=30.0)
    sink = _FakeSink()

    plant = astara.Plant()
    plant.add("a", a)
    plant.add("b", b)
    plant.add("sink", sink)
    plant.connect([a, b], sink)

    plant.step(dt=1.0e-3)

    assert sink.last_inlet["mass_flow_kg_s"] == 40.0
    expected_T = (600.0 * 10.0 + 500.0 * 30.0) / 40.0
    assert abs(sink.last_inlet["temperature_K"] - expected_T) < 1.0e-9


def test_fan_out_broadcasts_same_value():
    _register_fakes()
    astara.network.register_adapter(
        "_FakeSource2",
        lambda src: {"temperature_K": src.temperature_K},
        lambda src, **fields: None,
    )

    class _FakeSource2(_FakeSource):
        pass

    source = _FakeSource(temperature_K=590.0, mass_flow_kg_s=100.0)
    sink_a, sink_b = _FakeSink(), _FakeSink()

    plant = astara.Plant()
    plant.add("source", source)
    plant.add("sink_a", sink_a)
    plant.add("sink_b", sink_b)
    plant.connect(source, [sink_a, sink_b])

    plant.step(dt=1.0e-3)

    assert sink_a.last_inlet["temperature_K"] == 590.0
    assert sink_b.last_inlet["temperature_K"] == 590.0


def test_unknown_type_raises_on_connect():
    plant = astara.Plant()
    try:
        plant.connect(object(), object())
    except TypeError:
        pass
    else:
        raise AssertionError("expected TypeError for an unadapted type")
