import importlib.util
from pathlib import Path
from types import SimpleNamespace
import pytest
from biosim import BioWorld, ExecutionPolicy
from src.lotka_volterra import LotkaVolterraSystem
p=Path(__file__).resolve().parents[1]/"src/ecology_visualisation.py"
spec=importlib.util.spec_from_file_location("lotka_presenter",p)
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)

@pytest.mark.parametrize("prey,predator",[(10,5),(0,5),(10,0)])
def test_full_graph_includes_first_and_last_points(prey,predator):
    world=BioWorld(communication_step=.01)
    world.add_biomodule("core",LotkaVolterraSystem(prey_initial=prey,predator_initial=predator))
    v=mod.EcologyVisualisationModel(source_alias="core")
    world.add_biomodule("visualisation",v)
    world.connect("core.visualisation_payload","visualisation.core_visualisation_payload")
    world.setup();world.run(2)
    assert v.execution_policy is ExecutionPolicy.ONCE_AFTER_RUN
    history=world.get_outputs("core")["visualisation_payload"].value["payload"]["history"]
    visuals=v.visualize()
    population=next(x for x in visuals if x["data"]["title"]=="Population Trajectories")
    assert population["data"]["series"][0]["points"]==[[p["t"],p["prey"]] for p in history]
    assert history[0]["t"]==0 and history[-1]["t"]==2
    audit=next(x for x in visuals if "Invariant" in x["data"]["title"])
    assert audit["render"] == ("timeseries" if prey and predator else "table")

@pytest.mark.parametrize("times", [[0,1],[0,float("nan"),2],[0,2,1,2],[.1,2]])
def test_presenter_rejects_incomplete_or_invalid_time_axis(times):
    v=mod.EcologyVisualisationModel(source_alias="core")
    payload={"status":"completed","history":[{"t":t,"prey":10,"predator":5} for t in times]}
    with pytest.raises(ValueError):v.execute({"core_visualisation_payload":{"payload":payload}},context=SimpleNamespace(run_start=0,run_end=2))
