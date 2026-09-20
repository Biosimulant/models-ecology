from biosim import BioWorld
from src.lotka_volterra import LotkaVolterraSystem

def test_world_emits_typed_terminal_signals():
    world=BioWorld(communication_step=.01)
    m=LotkaVolterraSystem()
    world.add_biomodule("core",m)
    world.setup()
    world.run(.03)
    out=world.get_outputs("core")
    assert set(out)==set(m.outputs())
    assert out["prey_population_state"].spec.emitted_unit=="count"
    assert out["prey_population_state"].value["t"]==.03
    assert len(out["visualisation_payload"].value["payload"]["history"])==4
    assert m.visualize() is None
