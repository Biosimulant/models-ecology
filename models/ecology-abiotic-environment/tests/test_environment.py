from __future__ import annotations


def test_emits_conditions(bsim):
    from src.environment import Environment

    env = Environment(temperature=22.0, water=80.0, food_availability=1.2, min_dt=1.0)
    env.advance_to(1.0)
    out = env.get_outputs()
    assert "conditions" in out
    sig = out["conditions"]
    assert isinstance(sig.value, dict)
    assert sig.value["temperature"] == 22.0
    assert sig.value["water"] == 80.0
    assert sig.value["food"] == 1.2

