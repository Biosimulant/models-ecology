from __future__ import annotations

import numpy as np
import pytest

from src.geci2022_gene_drive import (
    Geci2022GeneDriveModel,
    GENOTYPES, GAMETES, N_GENOTYPES, N_GAMETES,
    _build_selection_vector, _build_homing_matrix, _build_editing_matrix,
    _build_recombination_matrix, _build_mutation_matrix, _build_gamete_matrices,
)


# ---------------------------------------------------------------------------
# Genotype enumeration tests
# ---------------------------------------------------------------------------

def test_genotype_count() -> None:
    """1071 genotypes = 51 sex genotypes x 21 autosome pairs."""
    assert N_GENOTYPES == 1071


def test_gamete_count() -> None:
    """66 gametes = 11 sex chromosomes x 6 autosomes."""
    assert N_GAMETES == 66


# ---------------------------------------------------------------------------
# Matrix conservation tests
# ---------------------------------------------------------------------------

_DEFAULT_PARAMS = {
    "homing_efficiency": 0.95, "editing_efficiency": 0.95, "shredding_efficiency": 0.95,
    "copy_mutation_rate": 0.01, "background_mutation_rate": 0.01,
    "editing_resistance_rate": 0.05, "shredding_resistance_rate": 0.05, "homing_resistance_rate": 0.05,
    "fitness_cost_cas9": 0.0, "fitness_cost_grna": 0.0, "fitness_cost_shredder": 0.0,
    "fitness_cost_nuclease": 0.0, "fitness_cost_shredder_activity": 0.0,
    "fitness_cost_edited_female": 0.0, "fitness_cost_edited_male": 0.0,
    "dominance_editing": 0.5, "dominance_shredder": 0.5, "dominance_shredder_gamete": 0.5,
    "cas9_cofactor": 1.0, "recombination_rate": 0.05,
}


def test_homing_matrix_row_sums() -> None:
    """Homing matrix rows must sum to 1 (conservative)."""
    mat = _build_homing_matrix(_DEFAULT_PARAMS)
    row_sums = mat.sum(axis=1)
    np.testing.assert_allclose(row_sums, 1.0, atol=1e-12)


def test_editing_matrix_row_sums() -> None:
    """Editing matrix rows must sum to 1 (conservative)."""
    mat = _build_editing_matrix(_DEFAULT_PARAMS)
    row_sums = mat.sum(axis=1)
    np.testing.assert_allclose(row_sums, 1.0, atol=1e-12)


def test_recombination_matrix_row_sums() -> None:
    """Recombination matrix rows must sum to 1 (conservative)."""
    mat = _build_recombination_matrix(_DEFAULT_PARAMS)
    row_sums = mat.sum(axis=1)
    np.testing.assert_allclose(row_sums, 1.0, atol=1e-12)


def test_gamete_matrix_row_sums() -> None:
    """Gamete matrix rows sum to 1 for genotypes that produce gametes, 0 otherwise."""
    sperm, egg = _build_gamete_matrices(_DEFAULT_PARAMS)
    for mat, label in [(sperm, "sperm"), (egg, "egg")]:
        for i in range(N_GENOTYPES):
            rs = mat[i, :].sum()
            assert rs == pytest.approx(0.0, abs=1e-12) or rs == pytest.approx(1.0, abs=1e-12), (
                f"{label} matrix row {i} sums to {rs}")


# ---------------------------------------------------------------------------
# Ecology behavior tests
# ---------------------------------------------------------------------------

def test_wildtype_equilibrium() -> None:
    """Without transgenic release, population should remain stable."""
    module = Geci2022GeneDriveModel(release_size=0.0)
    module.advance_window(0.0, 10.0)

    eco = module._compute_ecology()
    # All wild-type females, no males (no Y chromosomes without release)
    assert eco["total_adults"] > 0
    # Drive frequency should be 0
    assert eco["drive_frequency"] == pytest.approx(0.0, abs=1e-12)


def test_gene_drive_suppresses_population() -> None:
    """With transgenic release, population should decline over generations."""
    module = Geci2022GeneDriveModel(release_size=0.1)
    initial_total = module._initial_total
    module.advance_window(0.0, 30.0)

    final_total = module._compute_ecology()["total_adults"]
    assert final_total < initial_total


def test_sex_ratio_distortion() -> None:
    """X-shredding should bias sex ratio toward males."""
    module = Geci2022GeneDriveModel(release_size=0.1, shredding_efficiency=0.95)
    module.advance_window(0.0, 10.0)

    eco = module._compute_ecology()
    if eco["total_adults"] > 0:
        assert eco["male_fraction"] > 0.5 or eco["drive_frequency"] < 0.01


def test_outputs_and_history_accumulate() -> None:
    module = Geci2022GeneDriveModel()
    module.advance_window(0.0, 5.0)

    outputs = module.get_outputs()
    assert set(outputs) == {"population_state", "gene_drive_metrics", "visualisation_payload"}
    assert len(module._history) == 5


def test_drive_metrics_stay_in_physical_bounds() -> None:
    module = Geci2022GeneDriveModel()
    module.advance_window(0.0, 20.0)

    for record in module._history:
        assert 0.0 <= record["drive_frequency"] <= 1.0
        assert 0.0 <= record["resistance_frequency"] <= 1.0
        assert record["total_adults"] >= 0.0


def test_visualisation_payload_uses_multi_point_trajectories() -> None:
    module = Geci2022GeneDriveModel()
    module.advance_window(0.0, 5.0)
    payload = module.get_outputs()["visualisation_payload"].value["payload"]
    assert payload["initial_total"] > 0.0
    assert len(payload["history"]) > 1
