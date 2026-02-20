from __future__ import annotations


def test_predation_emits_events(biosim):
    from biosim.signals import BioSignal, SignalMetadata
    from src.predator_prey import PredatorPreyInteraction

    mod = PredatorPreyInteraction(predation_rate=1.0, conversion_efficiency=0.5, seed=1, min_dt=1.0)
    mod.set_inputs(
        {
            "prey_state": BioSignal(
                source="prey",
                name="prey_state",
                value={"species": "Rabbits", "count": 200, "t": 0.0},
                time=0.0,
                metadata=SignalMetadata(description="test", kind="state"),
            ),
            "predator_state": BioSignal(
                source="pred",
                name="predator_state",
                value={"species": "Foxes", "count": 10, "t": 0.0},
                time=0.0,
                metadata=SignalMetadata(description="test", kind="state"),
            ),
        }
    )
    mod.advance_to(1.0)
    out = mod.get_outputs()
    # With a large expected kill count, we should emit at least predation.
    assert "predation" in out
    assert out["predation"].value["kills"] > 0

