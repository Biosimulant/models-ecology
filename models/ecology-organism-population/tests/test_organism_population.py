from __future__ import annotations


def test_emits_population_state(bsim):
    from bsim.signals import BioSignal, SignalMetadata
    from src.organism_population import OrganismPopulation

    pop = OrganismPopulation(name="Rabbits", initial_count=100, seed=1, min_dt=1.0)
    pop.set_inputs(
        {
            "conditions": BioSignal(
                source="env",
                name="conditions",
                value={"temperature": 20.0, "water": 100.0, "food": 1.0, "sunlight": 1.0, "t": 0.0},
                time=0.0,
                metadata=SignalMetadata(description="test", kind="state"),
            )
        }
    )
    pop.advance_to(1.0)
    out = pop.get_outputs()
    assert "population_state" in out
    payload = out["population_state"].value
    assert payload["species"] == "Rabbits"
    assert isinstance(payload["count"], int)

