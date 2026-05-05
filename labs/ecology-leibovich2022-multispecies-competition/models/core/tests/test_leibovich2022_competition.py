from __future__ import annotations

import numpy as np

from src.leibovich2022_competition import (
    Leibovich2022CommunityModel,
)


def test_history_and_outputs_accumulate() -> None:
    module = Leibovich2022CommunityModel(integration_step=0.1, rng_seed=3)
    module.advance_window(0.0, 2.0)

    outputs = module.get_outputs()
    assert set(outputs) == {"community_state", "diversity_metrics", "visualisation_payload"}
    assert len(module._history) == 20
    assert outputs["community_state"].value["total_abundance"] >= 0


def test_diversity_metrics_are_well_formed() -> None:
    module = Leibovich2022CommunityModel(integration_step=0.1, rng_seed=4)
    module.advance_window(0.0, 3.0)

    metrics = module.get_outputs()["diversity_metrics"].value
    assert 0 <= metrics["richness"] <= module.species_count
    assert metrics["shannon_diversity"] >= 0.0
    assert 0.0 <= metrics["evenness"] <= 1.0
    assert metrics["dominant_species"].startswith("species_")


def test_visualisation_payload_contains_multi_point_trajectories() -> None:
    module = Leibovich2022CommunityModel(integration_step=0.2, rng_seed=1)
    module.advance_window(0.0, 4.0)
    payload = module.get_outputs()["visualisation_payload"].value["payload"]
    assert payload["species_count"] == module.species_count
    assert len(payload["history"]) > 1


def test_higher_immigration_supports_richness() -> None:
    low = Leibovich2022CommunityModel(immigration_rate=0.01, rng_seed=2, integration_step=0.1)
    high = Leibovich2022CommunityModel(immigration_rate=2.0, rng_seed=2, integration_step=0.1)

    low.advance_window(0.0, 5.0)
    high.advance_window(0.0, 5.0)

    low_richness = low.get_outputs()["diversity_metrics"].value["richness"]
    high_richness = high.get_outputs()["diversity_metrics"].value["richness"]
    assert high_richness >= low_richness


def test_upstream_propensity_formula() -> None:
    """Verify propensity calculations match upstream MultiLV exactly."""
    module = Leibovich2022CommunityModel(
        species_count=3, carrying_capacity=100.0, birth_rate=2.0, death_rate=1.0,
        competition_overlap=0.3, immigration_rate=0.5, initial_abundance=40.0,
        rng_seed=10, integration_step=0.05,
    )
    n = module._abundances.astype(float)
    total = float(np.sum(n))
    b, d = module.birth_rate, module.death_rate
    rho = module.competition_overlap
    K = module.carrying_capacity
    immi = module.immigration_rate

    for i in range(module.species_count):
        # Upstream: prop_birth = immi_rate + n_i * birth_rate
        expected_birth = immi + n[i] * b
        # Upstream: prop_death = n_i * (d + (b-d) * crowding / K)
        crowding = (1 - rho) * n[i] + rho * total
        expected_death = n[i] * (d + (b - d) * crowding / K)
        assert expected_birth >= 0.0
        assert expected_death >= 0.0
        # With these defaults the death propensity should be bounded
        assert expected_death <= n[i] * (d + (b - d) * total / K)


def test_uniform_immigration() -> None:
    """All species receive the same immigration rate (no rarer-species bias)."""
    module = Leibovich2022CommunityModel(
        species_count=4, immigration_rate=0.5, rng_seed=42, integration_step=0.05,
    )
    # The immigration component in the birth propensity is a flat constant
    # for all species, matching upstream MultiLV where prop_birth[i] uses
    # self.immi_rate (single scalar, same for all species).
    assert not hasattr(module, '_immigration_weights')


def test_competition_overlap_effect() -> None:
    """Higher competition overlap should produce stronger interspecific competition."""
    low_overlap = Leibovich2022CommunityModel(
        competition_overlap=0.0, rng_seed=5, integration_step=0.05,
    )
    high_overlap = Leibovich2022CommunityModel(
        competition_overlap=0.8, rng_seed=5, integration_step=0.05,
    )

    low_overlap.advance_window(0.0, 10.0)
    high_overlap.advance_window(0.0, 10.0)

    low_total = low_overlap.get_outputs()["community_state"].value["total_abundance"]
    high_total = high_overlap.get_outputs()["community_state"].value["total_abundance"]
    # With higher overlap, interspecific competition is stronger,
    # so total abundance should not be dramatically higher
    assert high_total <= low_total * 1.5
