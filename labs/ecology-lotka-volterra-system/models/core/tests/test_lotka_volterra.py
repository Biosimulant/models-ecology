from types import SimpleNamespace
from pathlib import Path
import math
import json
import numpy as np
import pytest
import yaml
from scipy.integrate import solve_ivp
from biosim import BioWorld, ExecutionPolicy
from src.lotka_volterra import LotkaVolterraSystem


def run(module, end=20.0, step=0.01):
    world = BioWorld(communication_step=step)
    world.add_biomodule("core", module)
    world.setup()
    world.run(end)
    return world.get_outputs("core")["visualisation_payload"].value["payload"]


def test_default_full_trajectory_against_independent_dop853():
    payload = run(LotkaVolterraSystem(), end=250)
    history = payload["history"]
    t = np.array([p["t"] for p in history])
    actual = np.array([[p["prey"],p["predator"]] for p in history])
    ref = solve_ivp(lambda t,y: [1.1*y[0]-.4*y[0]*y[1],.1*y[0]*y[1]-.4*y[1]], (0,250), [10,5], method="DOP853", rtol=1e-12, atol=1e-13, dense_output=True)
    assert ref.success
    expected = ref.sol(t).T
    np.testing.assert_allclose(actual,expected,rtol=1e-5,atol=1e-6)
    assert t[0] == 0 and t[-1] == pytest.approx(250)
    assert len(t) >= 25001
    assert np.all(np.diff(t) > 0)
    drift = max(abs(p["drift"]) for p in history)
    assert drift < 1e-6
    assert payload["invariant_applicable"] is True
    assert payload["status"] == "completed"
    print(json.dumps({"protocol":"250-day default","points":len(t),"max_absolute_error":float(np.max(np.abs(actual-expected))),"max_invariant_drift":drift}))


def test_equilibrium_is_stationary():
    h = run(LotkaVolterraSystem(prey_initial=4,predator_initial=2.75))["history"]
    assert all(p["prey"] == 4 and p["predator"] == 2.75 for p in h)


@pytest.mark.parametrize("prey,predator", [(10,0),(0,5),(0,0)])
def test_boundary_cases_against_analytic_solution(prey,predator):
    payload = run(LotkaVolterraSystem(prey_initial=prey,predator_initial=predator),end=2)
    for p in payload["history"]:
        assert p["prey"] == pytest.approx(prey*math.exp(1.1*p["t"]), rel=1e-8)
        assert p["predator"] == pytest.approx(predator*math.exp(-.4*p["t"]), rel=1e-8)
        assert p["invariant"] is None and p["drift"] is None
    assert not payload["invariant_applicable"]


def test_step_halving_converges_against_dop853():
    ref=solve_ivp(lambda t,y:[1.1*y[0]-.4*y[0]*y[1],.1*y[0]*y[1]-.4*y[1]],(0,20),[10,5],method="DOP853",rtol=1e-12,atol=1e-13,dense_output=True)
    errors=[]
    for step in [.1,.05,.025]:
        history=run(LotkaVolterraSystem(integration_step=step),step=step)["history"]
        expected=ref.sol([p["t"] for p in history]).T
        actual=np.array([[p["prey"],p["predator"]] for p in history])
        errors.append(float(np.max(np.abs(actual-expected))))
    assert errors[0]/errors[1] > 12
    assert errors[1]/errors[2] > 12
    print(json.dumps({"protocol":"RK4 step halving","steps":[.1,.05,.025],"errors":errors}))


@pytest.mark.parametrize("bad", [-1,float("nan"),float("inf"),True,None,"invalid"])
def test_invalid_port_update_is_atomic(bad):
    m=LotkaVolterraSystem()
    with pytest.raises(ValueError):
        m.set_inputs({"prey_growth_rate":2,"predation_rate":bad})
    assert m.alpha == 1.1 and m.beta == .4


@pytest.mark.parametrize("bad", [-1,float("nan"),float("inf"),True])
@pytest.mark.parametrize("arg", ["alpha","beta","gamma","delta","prey_initial","predator_initial","integration_step"])
def test_constructor_rejects_invalid_values(arg,bad):
    with pytest.raises(ValueError): LotkaVolterraSystem(**{arg:bad})


def test_initial_overrides_reset_invariant_baseline():
    m=LotkaVolterraSystem()
    m.set_inputs({"prey_initial_population":12,"prey_growth_rate":.8})
    p=run(m,end=2)
    assert p["history"][0]["prey"] == 12
    assert p["history"][0]["drift"] == 0
    assert max(abs(x["drift"]) for x in p["history"]) < 1e-6


def test_changing_initial_population_after_start_is_rejected():
    m=LotkaVolterraSystem()
    m.execute({},context=SimpleNamespace(window_end=1,run_end=2))
    with pytest.raises(ValueError):m.set_inputs({"prey_initial_population":11})
    m.set_inputs({"prey_growth_rate":.8})
    p=m.execute({},context=SimpleNamespace(window_end=2,run_end=2))["visualisation_payload"]["payload"]
    assert not p["invariant_applicable"]


def test_unstable_step_raises_instead_of_clipping():
    m=LotkaVolterraSystem(integration_step=10)
    with pytest.raises(FloatingPointError):m.execute({},context=SimpleNamespace(window_end=10,run_end=10))


def test_reset_repeats_exact_results():
    m=LotkaVolterraSystem()
    a=run(m,end=2)
    b=run(m,end=2)
    assert a==b


def test_manifest_contract():
    root=Path(__file__).resolve().parents[3]
    lab=yaml.safe_load((root/"lab.yaml").read_text())
    m=LotkaVolterraSystem()
    assert m.execution_policy is ExecutionPolicy.EACH_WINDOW
    assert {p["maps_to"].split(".")[1] for p in lab["io"]["inputs"]} == set(m.inputs())
    assert {p["name"] for p in lab["io"]["outputs"]} == {"prey_population_state","predator_population_state","trajectory"}
