from __future__ import annotations

import sys
from pathlib import Path


MODEL_ROOT = Path(__file__).resolve().parents[1]
MONOREPO_ROOT = MODEL_ROOT.parents[3]
BSIM_SRC = MONOREPO_ROOT / "bsim-active" / "biosim" / "src"

for path in (str(MODEL_ROOT), str(BSIM_SRC)):
    if path not in sys.path:
        sys.path.insert(0, path)


from src.leibovich2022_multispecies_eco_competition_descr_model2212080001_model import (  # noqa: E402
    Leibovich2022MultispeciesEcoCompetitionDescrModel2212080001Model,
)


def test_history_and_outputs_accumulate() -> None:
    module = Leibovich2022MultispeciesEcoCompetitionDescrModel2212080001Model(min_dt=0.1, rng_seed=3)
    module.advance_to(2.0)

    outputs = module.get_outputs()
    assert set(outputs) == {"community_state", "diversity_metrics"}
    assert len(module._history) == 20
    assert outputs["community_state"].value["total_abundance"] >= 0


def test_diversity_metrics_are_well_formed() -> None:
    module = Leibovich2022MultispeciesEcoCompetitionDescrModel2212080001Model(min_dt=0.1, rng_seed=4)
    module.advance_to(3.0)

    metrics = module.get_outputs()["diversity_metrics"].value
    assert 0 <= metrics["richness"] <= module.species_count
    assert metrics["shannon_diversity"] >= 0.0
    assert 0.0 <= metrics["evenness"] <= 1.0
    assert metrics["dominant_species"].startswith("species_")


def test_visuals_contain_multi_point_trajectories() -> None:
    module = Leibovich2022MultispeciesEcoCompetitionDescrModel2212080001Model(min_dt=0.2, rng_seed=1)
    module.advance_to(4.0)
    visuals = module.visualize()

    assert isinstance(visuals, list)
    assert [visual["render"] for visual in visuals] == ["timeseries", "timeseries", "table"]
    for visual in visuals[:2]:
        for series in visual["data"]["series"]:
            assert len(series["points"]) > 1


def test_higher_immigration_supports_richness() -> None:
    low = Leibovich2022MultispeciesEcoCompetitionDescrModel2212080001Model(immigration_rate=0.05, rng_seed=2, min_dt=0.1)
    high = Leibovich2022MultispeciesEcoCompetitionDescrModel2212080001Model(immigration_rate=1.2, rng_seed=2, min_dt=0.1)

    low.advance_to(5.0)
    high.advance_to(5.0)

    low_richness = low.get_outputs()["diversity_metrics"].value["richness"]
    high_richness = high.get_outputs()["diversity_metrics"].value["richness"]
    assert high_richness >= low_richness
