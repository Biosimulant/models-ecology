from __future__ import annotations


def test_emits_metrics_signal(biosim):
    from biosim.signals import BioSignal, SignalMetadata
    from src.ecology_metrics import EcologyMetrics

    mod = EcologyMetrics(min_dt=1.0)
    mod.set_inputs(
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
    mod.set_inputs(
        {
            "population_state": BioSignal(
                source="foxes",
                name="population_state",
                value={"species": "Foxes", "count": 5, "t": 0.0},
                time=0.0,
                metadata=SignalMetadata(description="test", kind="state"),
            )
        }
    )
    mod.advance_to(1.0)
    out = mod.get_outputs()
    assert "metrics" in out
    payload = out["metrics"].value
    assert payload["n_species"] >= 2

