from __future__ import annotations


def test_emits_phase_point(bsim):
    from bsim.signals import BioSignal, SignalMetadata
    from src.phase_space import PhaseSpaceMonitor

    mon = PhaseSpaceMonitor(x_species="Rabbits", y_species="Foxes", max_points=10, min_dt=1.0)
    mon.set_inputs(
        {
            "population_state": BioSignal(
                source="rabbits",
                name="population_state",
                value={"species": "Rabbits", "count": 80, "t": 0.0},
                time=0.0,
                metadata=SignalMetadata(description="test", kind="state"),
            )
        }
    )
    mon.set_inputs(
        {
            "population_state": BioSignal(
                source="foxes",
                name="population_state",
                value={"species": "Foxes", "count": 12, "t": 0.0},
                time=0.0,
                metadata=SignalMetadata(description="test", kind="state"),
            )
        }
    )
    mon.advance_to(1.0)
    out = mon.get_outputs()
    assert "phase_point" in out
    payload = out["phase_point"].value
    assert payload["x"] == 80
    assert payload["y"] == 12

