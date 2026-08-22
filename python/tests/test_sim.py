"""Tests for `astara.Simulation`'s scheduling/recording, decoupled from any
particular plant/component -- a tiny fake plant is enough to exercise the
run-loop logic itself."""

import astara


class _CountingPlant:
    def __init__(self):
        self.steps = 0

    def step(self, dt):
        self.steps += 1


def test_schedule_fires_exactly_once_at_the_right_time():
    plant = _CountingPlant()
    sim = astara.Simulation(plant, dt=0.1)
    calls = []
    sim.schedule(t=0.5, action=lambda: calls.append(plant.steps))
    sim.run(t_end=1.0)
    assert len(calls) == 1
    # dt=0.1 -> 10 steps to reach t=1.0; t=0.5 is step index 5.
    assert calls[0] == 5


def test_record_produces_expected_row_count():
    plant = _CountingPlant()
    sim = astara.Simulation(plant, dt=0.1)
    sim.record(lambda: {"steps": plant.steps}, every=0.2)
    df = sim.run(t_end=1.0)
    assert len(df) == 6  # t = 0.0, 0.2, ..., 1.0 inclusive


def test_run_returns_none_with_no_recorders():
    plant = _CountingPlant()
    sim = astara.Simulation(plant, dt=0.1)
    result = sim.run(t_end=0.5)
    assert result is None
    assert plant.steps == 5


def test_multiple_recorders_return_a_list_in_order():
    plant = _CountingPlant()
    sim = astara.Simulation(plant, dt=0.1)
    sim.record(lambda: {"a": 1}, every=0.5)
    sim.record(lambda: {"b": 2}, every=0.25)
    dfs = sim.run(t_end=1.0)
    assert isinstance(dfs, list) and len(dfs) == 2
    assert list(dfs[0]["a"]) == [1, 1, 1]
    assert list(dfs[1]["b"]) == [2, 2, 2, 2, 2]
