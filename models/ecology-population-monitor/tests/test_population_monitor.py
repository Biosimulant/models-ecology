from __future__ import annotations


def test_emits_population_summary(biosim):
    from biosim.signals import BioSignal, SignalMetadata
    from src.population_monitor import PopulationMonitor

    mon = PopulationMonitor(max_points=10, min_dt=1.0)
    mon.set_inputs(
        {
            "population_state": BioSignal(
                source="rabbits",
                name="population_state",
                value={"species": "Rabbits", "count": 50, "t": 0.0},
                time=0.0,
                metadata=SignalMetadata(description="test", kind="state"),
            )
        }
    )
    mon.advance_to(1.0)
    out = mon.get_outputs()
    assert "population_summary" in out
    payload = out["population_summary"].value
    assert payload["latest_counts"]["Rabbits"] == 50

