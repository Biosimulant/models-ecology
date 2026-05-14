# SPDX-FileCopyrightText: 2025-present Demi <bjaiye1@gmail.com>
#
# SPDX-License-Identifier: MIT
"""Faithful eco-genetic gene-drive suppression model ported from Geci et al. (2022).

The upstream BioModels asset for MODEL2301120001 is a Julia implementation
(not SBML).  This package is a parity-tested Python port of the upstream
genotype-tracking mechanics needed for the Biosimulant lab scenario:

  - Genotype construction with Y-linked editor, X-shredder, and autosomal
    homing components.
  - Process matrices: selection, homing, editing, recombination, mutation,
    and gamete production (sperm/egg with sex-specific X-shredding).
  - Discrete-generation simulation loop: mutation -> homing -> editing ->
    recombination -> gamete production -> zygote formation ->
    density-dependent survival (Beverton-Holt) -> selection.
  - Ecology outputs (population, sex ratio, drive frequency, resistance)
    are derived from the full genotype state, not local approximations.
"""
from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Tuple, TYPE_CHECKING

import biosim
import numpy as np
from biosim.signals import (AcceptedSignalProfile, ArraySignal, BioSignal,
                            EventSignal, RecordSignal, ScalarSignal, SignalSpec)
from biosim.signals import make_signal as _make_signal

if TYPE_CHECKING:  # pragma: no cover - typing only
    from biosim.visuals import VisualSpec


# ---------------------------------------------------------------------------
# Signal helpers (shared with other ecology models)
# ---------------------------------------------------------------------------

def _schema_type(value):
    if isinstance(value, bool):
        return "bool"
    if isinstance(value, int) and not isinstance(value, bool):
        return "int"
    if isinstance(value, float):
        return "float"
    if isinstance(value, str):
        return "str"
    return "json"



# ---------------------------------------------------------------------------
# Genotype enumeration (port of Julia lines 1-210)
# ---------------------------------------------------------------------------

# Allele-level chromosomes (detailed representation)
# Y chromosome: locus1 (Cas9) x locus2 (gRNA1), plus wt
_Y_DETAILED = ["ab", "AB", "Aβ", "αB", "αβ"]
# Autosomes: locus1 (X-shredder) x locus2 (gRNA2), plus wt and r3
_A_DETAILED = ["cd", "CD", "Cδ", "ζD", "ζδ", "r3"]
# X chromosomes: editing-target (e/ε/r1) x shredding-target (f/r2)
_X_DETAILED = ["ef", "er2", "εf", "εr2", "r1f", "r1r2"]

_ALL_SEX_CHROMS = _Y_DETAILED + _X_DETAILED  # 11 total

# Shreddable X chromosomes (contain "f" = wt shredding target)
_SHREDDABLE_X = [x for x in _X_DETAILED if "f" in x]  # ef, εf, r1f


def _build_genotypes():
    """Build the full genotype list matching Julia's enumeration.

    Each genotype is a 4-tuple: (sex_chrom1, sex_chrom2, auto1, auto2)
    Males: (Y_variant, X_variant, auto_i, auto_j) with i <= j
    Females: (X_variant_i, X_variant_j, auto_k, auto_l) with i <= j, k <= l
    """
    sex_genotypes = []
    # Males: Y x X
    for y in _Y_DETAILED:
        for x in _X_DETAILED:
            sex_genotypes.append((y, x))
    # Females: X x X upper triangle
    for i, x1 in enumerate(_X_DETAILED):
        for x2 in _X_DETAILED[i:]:
            sex_genotypes.append((x1, x2))

    genotypes = []
    # For each sex genotype, pair with autosome pairs (upper triangle)
    for sg in sex_genotypes:
        for i, a1 in enumerate(_A_DETAILED):
            for a2 in _A_DETAILED[i:]:
                genotypes.append((sg[0], sg[1], a1, a2))

    return genotypes, sex_genotypes


def _build_gametes():
    """Build the gamete list: (sex_chromosome, autosome)."""
    gametes = []
    for sc in _ALL_SEX_CHROMS:
        for auto in _A_DETAILED:
            gametes.append((sc, auto))
    return gametes


# Module-level constants
GENOTYPES, _SEX_GENOTYPES = _build_genotypes()
GAMETES = _build_gametes()
N_GENOTYPES = len(GENOTYPES)  # 1071
N_GAMETES = len(GAMETES)  # 66

# Index lookup dicts
_GENO_INDEX = {g: i for i, g in enumerate(GENOTYPES)}
_GAMETE_INDEX = {g: i for i, g in enumerate(GAMETES)}

# Precompute which genotypes are male vs female
_MALE_INDICES = np.array([i for i, g in enumerate(GENOTYPES) if g[0] in _Y_DETAILED], dtype=int)
_FEMALE_INDICES = np.array([i for i, g in enumerate(GENOTYPES) if g[0] in _X_DETAILED], dtype=int)

# Sex breakpoint: first female genotype index
_SEX_BREAKPOINT = int(_MALE_INDICES[-1]) + 1 if len(_MALE_INDICES) > 0 else 0


def _find_genotype(geno_tuple):
    """Find genotype index, handling autosome pair symmetry."""
    idx = _GENO_INDEX.get(geno_tuple)
    if idx is not None:
        return idx
    # Try swapping autosomes (positions 2,3)
    swapped = (geno_tuple[0], geno_tuple[1], geno_tuple[3], geno_tuple[2])
    idx = _GENO_INDEX.get(swapped)
    if idx is not None:
        return idx
    # Try swapping sex chromosomes for females
    if geno_tuple[0] in _X_DETAILED and geno_tuple[1] in _X_DETAILED:
        swapped_sex = (geno_tuple[1], geno_tuple[0], geno_tuple[2], geno_tuple[3])
        idx = _GENO_INDEX.get(swapped_sex)
        if idx is not None:
            return idx
        swapped_both = (geno_tuple[1], geno_tuple[0], geno_tuple[3], geno_tuple[2])
        idx = _GENO_INDEX.get(swapped_both)
        if idx is not None:
            return idx
    return None


# ---------------------------------------------------------------------------
# Zygote formation lookup table
# ---------------------------------------------------------------------------

def _build_zygote_map():
    """Pre-compute (sperm_gamete_idx, egg_gamete_idx) -> genotype_idx."""
    zmap = np.full((N_GAMETES, N_GAMETES), -1, dtype=int)
    for si, sg in enumerate(GAMETES):
        for ei, eg in enumerate(GAMETES):
            # Sperm contributes: sex_chrom -> position 0, autosome -> position 2
            # Egg contributes: sex_chrom -> position 1, autosome -> position 3
            geno = (sg[0], eg[0], sg[1], eg[1])
            idx = _find_genotype(geno)
            if idx is not None:
                zmap[si, ei] = idx
    return zmap


_ZYGOTE_MAP = _build_zygote_map()


# ---------------------------------------------------------------------------
# Matrix builders (port of Julia matrix construction functions)
# ---------------------------------------------------------------------------

def _build_selection_vector(params):
    """Port of create_selection_matrix (Julia lines 234-344).

    Returns a (N_GENOTYPES,) array of fitness multipliers.
    """
    sel = np.ones(N_GENOTYPES)
    s_f = params["fitness_cost_edited_female"]
    h_f = params["dominance_editing"]
    s_m = params["fitness_cost_edited_male"]
    s_a = params["fitness_cost_cas9"]
    s_b = params["fitness_cost_grna"]
    s_c = params["fitness_cost_shredder"]
    s_d = params["fitness_cost_nuclease"]
    s_e = params["fitness_cost_shredder_activity"]
    h_e = params["dominance_shredder"]

    for i, g in enumerate(GENOTYPES):
        sex1, sex2, aut1, aut2 = g
        fitness = 1.0

        # Desired fitness costs from edited X-chromosome targets
        is_female = sex1 in _X_DETAILED
        if is_female:
            # Female: check both X chromosomes for editing (ε)
            has_edit_1 = "ε" in sex1
            has_edit_2 = "ε" in sex2
            if has_edit_1 and has_edit_2:
                fitness *= (1 - s_f)
            elif has_edit_1 or has_edit_2:
                fitness *= (1 - s_f * h_f)
        else:
            # Male: Y in sex1, X in sex2
            if "ε" in sex2:
                fitness *= (1 - s_m)

        # Count molecular components
        n_cas9 = 1 if "A" in sex1 and sex1 not in _X_DETAILED else 0
        n_grna1 = 1 if "B" in sex1 and sex1 not in _X_DETAILED else 0
        n_shredder = ("C" in aut1) + ("C" in aut2)
        n_grna2 = ("D" in aut1) + ("D" in aut2)

        # Expression costs of proteins
        if n_cas9 == 1:
            fitness *= (1 - s_a)
        if n_shredder == 2:
            fitness *= (1 - s_c) ** 2
        elif n_shredder == 1:
            fitness *= (1 - s_c)

        # Expression costs of gRNAs
        if n_grna1 == 1:
            fitness *= (1 - s_b)
        if n_grna2 == 2:
            fitness *= (1 - s_b) ** 2
        elif n_grna2 == 1:
            fitness *= (1 - s_b)

        # Shredder activity cost
        if n_shredder == 2:
            fitness *= (1 - s_e)
        elif n_shredder == 1:
            fitness *= (1 - s_e * h_e)

        # Nuclease activity costs
        if n_cas9 == 1 and n_grna1 == 1:
            fitness *= (1 - s_d)
        if n_cas9 == 1 and n_grna2 >= 1:
            fitness *= (1 - s_d)

        sel[i] = fitness

    return sel


def _build_homing_matrix(params):
    """Port of create_homing_matrix (Julia lines 348-463).

    Homing requires: (1) functional Cas9 (A in sex1), (2) functional gRNA2
    (D in autosome), (3) wild-type autosome target (cd).
    """
    mat = np.eye(N_GENOTYPES)
    e_h = params["homing_efficiency"]
    er_3 = params["homing_resistance_rate"]
    m_1 = params["copy_mutation_rate"]

    for i, g in enumerate(GENOTYPES):
        sex1, sex2, aut1, aut2 = g
        has_cas9 = "A" in sex1 and sex1 not in _X_DETAILED
        has_grna2_in_aut2 = "D" in aut2
        target_is_wt = aut1 == "cd"

        if not (has_cas9 and has_grna2_in_aut2 and target_is_wt):
            continue

        # Determine what the correctly homed autosome 1 would become
        # (it copies autosome 2)
        post = (sex1, sex2, aut2, aut2)
        resistant = (sex1, sex2, "r3", aut2)

        i_post = _find_genotype(post)
        i_resistant = _find_genotype(resistant)
        if i_post is None or i_resistant is None:
            continue

        # Reduce pre-genotype frequency
        mat[i, i] *= (1 - e_h)
        # Correctly homed genotype increases
        mat[i, i_resistant] += e_h * er_3

        if aut2 == "CD":
            mat[i, i_post] += e_h * (1 - er_3) * (1 - m_1) ** 2

            # Mutations during copying
            c_mut = (sex1, sex2, "ζD", aut2)
            d_mut = (sex1, sex2, "Cδ", aut2)
            both_mut = (sex1, sex2, "ζδ", aut2)
            i_c = _find_genotype(c_mut)
            i_d = _find_genotype(d_mut)
            i_both = _find_genotype(both_mut)
            if i_c is not None:
                mat[i, i_c] += e_h * (1 - er_3) * m_1 * (1 - m_1)
            if i_d is not None:
                mat[i, i_d] += e_h * (1 - er_3) * m_1 * (1 - m_1)
            if i_both is not None:
                mat[i, i_both] += e_h * (1 - er_3) * m_1 ** 2

        elif aut2 == "ζD":
            mat[i, i_post] += e_h * (1 - er_3) * (1 - m_1)
            d_mut = (sex1, sex2, aut2, "ζδ")
            i_d = _find_genotype(d_mut)
            if i_d is not None:
                mat[i, i_d] += e_h * (1 - er_3) * m_1

    return mat


def _build_editing_matrix(params):
    """Port of create_editing_matrix (Julia lines 540-576).

    Editing requires: (1) functional Cas9+gRNA1 on Y (AB), (2) editable X
    chromosome (contains 'e' = wt editing target).
    """
    mat = np.eye(N_GENOTYPES)
    e_e = params["editing_efficiency"]
    er_1 = params["editing_resistance_rate"]

    for i, g in enumerate(GENOTYPES):
        sex1, sex2, aut1, aut2 = g
        if sex1 != "AB":
            continue
        if "e" not in sex2:
            continue

        post_x = sex2.replace("e", "ε")
        resistant_x = sex2.replace("e", "r1")
        post = (sex1, post_x, aut1, aut2)
        resistant = (sex1, resistant_x, aut1, aut2)

        i_post = _find_genotype(post)
        i_resistant = _find_genotype(resistant)
        if i_post is None or i_resistant is None:
            continue

        mat[i, i] = (1 - e_e)
        mat[i, i_post] = e_e * (1 - er_1)
        mat[i, i_resistant] = e_e * er_1

    return mat


def _build_recombination_matrix(params):
    """Port of create_recombination_matrix (Julia lines 466-537).

    Only affects double-heterozygous females at X-linked loci.
    """
    mat = np.eye(N_GENOTYPES)
    r_rate = params["recombination_rate"]

    for i, g in enumerate(GENOTYPES):
        sex1, sex2, aut1, aut2 = g
        # Only females can recombine X chromosomes
        if sex1 not in _X_DETAILED:
            continue

        # Extract editing alleles from each X
        def _edit_allele(x):
            if "r1" in x:
                return "r1"
            if "ε" in x:
                return "ε"
            return "e"

        def _shred_allele(x):
            if "r2" in x:
                return "r2"
            return "f"

        a1_e, a1_s = _edit_allele(sex1), _shred_allele(sex1)
        a2_e, a2_s = _edit_allele(sex2), _shred_allele(sex2)

        # Only recombine if both loci differ between the two X chromosomes
        if a1_e == a2_e or a1_s == a2_s:
            continue
        # Need both loci to differ (double heterozygote)
        if a1_e == a2_e and a1_s == a2_s:
            continue

        # Recombinant: swap shredding alleles
        rec_x1 = a1_e + a2_s
        rec_x2 = a2_e + a1_s
        post = _find_genotype((rec_x1, rec_x2, aut1, aut2))
        if post is None:
            continue

        mat[i, i] = 1 - r_rate
        mat[i, post] = r_rate

    return mat


def _build_mutation_matrix(params):
    """Port of create_mutation_matrix (Julia lines 577-688).

    Background mutations that convert functional elements to dysfunctional.
    """
    mat = np.eye(N_GENOTYPES)
    m_2 = params["background_mutation_rate"]

    if m_2 == 0.0:
        return mat

    for i, g in enumerate(GENOTYPES):
        sex1, sex2, aut1, aut2 = g

        # Find all possible single-step mutations for each chromosome position
        mutations_per_pos = [[], [], [], []]  # sex1, sex2, aut1, aut2

        # Sex chromosome 1 mutations (Y chromosomes only have functional elements)
        if sex1 == "AB":
            mutations_per_pos[0] = [("αB", 1), ("Aβ", 1), ("αβ", 2)]
        elif "A" in sex1 and sex1 not in _X_DETAILED:
            mutations_per_pos[0] = [(sex1.replace("A", "α"), 1)]
        elif "B" in sex1 and sex1 not in _X_DETAILED:
            mutations_per_pos[0] = [(sex1.replace("B", "β"), 1)]

        # Autosome 1 mutations
        if aut1 == "CD":
            mutations_per_pos[2] = [("ζD", 1), ("Cδ", 1), ("ζδ", 2)]
        elif "C" in aut1:
            mutations_per_pos[2] = [(aut1.replace("C", "ζ"), 1)]
        elif "D" in aut1:
            mutations_per_pos[2] = [(aut1.replace("D", "δ"), 1)]

        # Autosome 2 mutations
        if aut2 == "CD":
            mutations_per_pos[3] = [("ζD", 1), ("Cδ", 1), ("ζδ", 2)]
        elif "C" in aut2:
            mutations_per_pos[3] = [(aut2.replace("C", "ζ"), 1)]
        elif "D" in aut2:
            mutations_per_pos[3] = [(aut2.replace("D", "δ"), 1)]

        # Count total mutable elements
        all_mutations = []
        for pos_muts in mutations_per_pos:
            all_mutations.extend(pos_muts)

        if not all_mutations:
            continue

        max_mutations = max(m for _, m in all_mutations) if all_mutations else 0
        # Count total mutable positions
        total_mutable = 0
        for pos_muts in mutations_per_pos:
            if pos_muts:
                total_mutable = max(total_mutable, max(m for _, m in pos_muts))

        # Generate all combinations of mutations across positions
        # Each position can be: original (0 mutations) or one of its mutation targets
        pos_options = []
        for p_idx, pos_muts in enumerate(mutations_per_pos):
            options = [(g[p_idx], 0)]  # original = 0 mutations
            for allele, n_mut in pos_muts:
                options.append((allele, n_mut))
            pos_options.append(options)

        # Find total mutable elements for this genotype
        total_elements = sum(max((m for _, m in pm), default=0) for pm in mutations_per_pos if pm)

        # Iterate all combinations
        for s1_opt in pos_options[0]:
            for s2_opt in pos_options[1]:
                for a1_opt in pos_options[2]:
                    for a2_opt in pos_options[3]:
                        n_muts = s1_opt[1] + s2_opt[1] + a1_opt[1] + a2_opt[1]
                        if n_muts == 0:
                            continue
                        mutated = (s1_opt[0], s2_opt[0], a1_opt[0], a2_opt[0])
                        i_post = _find_genotype(mutated)
                        if i_post is None:
                            continue

                        n_not_mutated = total_elements - n_muts
                        prob = m_2 ** n_muts
                        if n_not_mutated > 0:
                            prob *= (1 - m_2) ** n_not_mutated

                        mat[i, i_post] += prob

        # Diagonal: probability of no mutations
        mat[i, i] = (1 - m_2) ** total_elements if total_elements > 0 else 1.0

    return mat


def _build_gamete_matrices(params):
    """Port of create_gamete_matrix (Julia lines 691-889).

    Returns (sperm_matrix, egg_matrix) each of shape (N_GENOTYPES, N_GAMETES).
    Males produce sperm; females produce eggs.
    X-shredding affects gamete production in males carrying shredder.
    """
    gamete_mat = np.zeros((N_GENOTYPES, N_GAMETES))
    e_s = params["shredding_efficiency"]
    er_2 = params["shredding_resistance_rate"]
    c_param = params["cas9_cofactor"]
    h_e2 = params["dominance_shredder_gamete"]

    for i, g in enumerate(GENOTYPES):
        sex1, sex2, aut1, aut2 = g

        # Four possible gametes from this genotype
        gametes = [
            (sex1, aut1), (sex1, aut2),
            (sex2, aut1), (sex2, aut2),
        ]

        is_female = sex1 in _X_DETAILED
        n_shredder = ("C" in aut1) + ("C" in aut2)
        is_shreddable = sex2 in _SHREDDABLE_X

        if is_female or n_shredder == 0 or not is_shreddable:
            # Standard Mendelian: each gamete at frequency count/4
            for gam in gametes:
                gi = _GAMETE_INDEX.get(gam)
                if gi is not None and gamete_mat[i, gi] == 0:
                    freq = gametes.count(gam)
                    gamete_mat[i, gi] = freq / 4.0
        else:
            # Male with X-shredder and shreddable X: X-bearing gametes are reduced
            # Gametes 0,1 carry sex1 (Y chromosome); gametes 2,3 carry sex2 (X chromosome)
            # The X-bearing gametes (2,3) are subject to shredding

            # Compute effective shredding rate
            has_AB = sex1 == "AB"
            if n_shredder == 2 and has_AB:
                eff_shred = e_s
            elif n_shredder == 2 and not has_AB:
                eff_shred = e_s * c_param
            elif n_shredder == 1 and has_AB:
                eff_shred = e_s * h_e2
            else:  # n_shredder == 1 and not has_AB
                eff_shred = e_s * h_e2 * c_param

            # Denominator: total viable gametes accounting for shredding
            denom = 4.0 - 2.0 * eff_shred * (1 - er_2)

            # Y-bearing gametes (indices 0,1 in gametes list)
            for gam_idx in [0, 1]:
                gi = _GAMETE_INDEX.get(gametes[gam_idx])
                if gi is not None:
                    gamete_mat[i, gi] += 1.0 / denom

            # X-bearing gametes (indices 2,3): reduced by shredding
            for gam_idx in [2, 3]:
                gam = gametes[gam_idx]
                gi = _GAMETE_INDEX.get(gam)
                if gi is not None:
                    gamete_mat[i, gi] += (1.0 - eff_shred) / denom

                # Shredding-resistant gametes
                res_gam = (gam[0].replace("f", "r2"), gam[1])
                gi_res = _GAMETE_INDEX.get(res_gam)
                if gi_res is not None:
                    gamete_mat[i, gi_res] += (eff_shred * er_2) / denom

    # Split into sperm (from males) and egg (from females) matrices
    sperm_matrix = np.zeros_like(gamete_mat)
    egg_matrix = np.zeros_like(gamete_mat)
    sperm_matrix[:_SEX_BREAKPOINT, :] = gamete_mat[:_SEX_BREAKPOINT, :]
    egg_matrix[_SEX_BREAKPOINT:, :] = gamete_mat[_SEX_BREAKPOINT:, :]

    return sperm_matrix, egg_matrix


def build_matrices(params):
    """Build all process matrices from a parameter dict."""
    return {
        "selection": _build_selection_vector(params),
        "homing": _build_homing_matrix(params),
        "editing": _build_editing_matrix(params),
        "recombination": _build_recombination_matrix(params),
        "mutation": _build_mutation_matrix(params),
        "sperm": _build_gamete_matrices(params)[0],
        "egg": _build_gamete_matrices(params)[1],
    }


def create_zygote_vector(sperm_freqs, egg_counts):
    """Create genotype vector from sperm frequencies and egg counts.

    Port of Julia create_zygote_vector (lines 892-913).
    """
    gv = np.zeros(N_GENOTYPES)
    for si in range(N_GAMETES):
        if sperm_freqs[si] == 0.0:
            continue
        for ei in range(N_GAMETES):
            if egg_counts[ei] == 0.0:
                continue
            gi = _ZYGOTE_MAP[si, ei]
            if gi >= 0:
                gv[gi] += sperm_freqs[si] * egg_counts[ei]
    return gv


# ---------------------------------------------------------------------------
# Biosim model class
# ---------------------------------------------------------------------------

class Geci2022GeneDriveModel(biosim.BioModule):
    """Faithful genotype-tracking gene-drive suppression model."""

    def __init__(
        self,
        net_reproduction_rate: float = 6.0,
        juvenile_survival: float = 0.1,
        initial_population: float = 1.0,
        homing_efficiency: float = 0.95,
        editing_efficiency: float = 0.95,
        shredding_efficiency: float = 0.95,
        copy_mutation_rate: float = 0.0,
        background_mutation_rate: float = 0.0,
        editing_resistance_rate: float = 0.0,
        shredding_resistance_rate: float = 0.0,
        homing_resistance_rate: float = 0.0,
        fitness_cost_cas9: float = 0.0,
        fitness_cost_grna: float = 0.0,
        fitness_cost_shredder: float = 0.0,
        fitness_cost_nuclease: float = 0.0,
        fitness_cost_shredder_activity: float = 0.0,
        fitness_cost_edited_female: float = 0.0,
        fitness_cost_edited_male: float = 0.0,
        dominance_editing: float = 0.5,
        dominance_shredder: float = 0.5,
        dominance_shredder_gamete: float = 0.5,
        cas9_cofactor: float = 1.0,
        recombination_rate: float = 0.0,
        release_size: float = 0.1,
    ) -> None:
        self.net_reproduction_rate = float(net_reproduction_rate)
        self.juvenile_survival = float(juvenile_survival)
        self.initial_population = float(initial_population)
        self.release_size = float(release_size)

        self._params = {
            "homing_efficiency": float(homing_efficiency),
            "editing_efficiency": float(editing_efficiency),
            "shredding_efficiency": float(shredding_efficiency),
            "copy_mutation_rate": float(copy_mutation_rate),
            "background_mutation_rate": float(background_mutation_rate),
            "editing_resistance_rate": float(editing_resistance_rate),
            "shredding_resistance_rate": float(shredding_resistance_rate),
            "homing_resistance_rate": float(homing_resistance_rate),
            "fitness_cost_cas9": float(fitness_cost_cas9),
            "fitness_cost_grna": float(fitness_cost_grna),
            "fitness_cost_shredder": float(fitness_cost_shredder),
            "fitness_cost_nuclease": float(fitness_cost_nuclease),
            "fitness_cost_shredder_activity": float(fitness_cost_shredder_activity),
            "fitness_cost_edited_female": float(fitness_cost_edited_female),
            "fitness_cost_edited_male": float(fitness_cost_edited_male),
            "dominance_editing": float(dominance_editing),
            "dominance_shredder": float(dominance_shredder),
            "dominance_shredder_gamete": float(dominance_shredder_gamete),
            "cas9_cofactor": float(cas9_cofactor),
            "recombination_rate": float(recombination_rate),
        }

        # Build process matrices (done once)
        sperm_mat, egg_mat = _build_gamete_matrices(self._params)
        self._matrices = {
            "selection": _build_selection_vector(self._params),
            "homing": _build_homing_matrix(self._params),
            "editing": _build_editing_matrix(self._params),
            "recombination": _build_recombination_matrix(self._params),
            "mutation": _build_mutation_matrix(self._params),
            "sperm": sperm_mat,
            "egg": egg_mat,
        }

        # Ecology constants
        self._f = (self.net_reproduction_rate * 2.0) / self.juvenile_survival  # eggs per female
        self._alpha = self.initial_population * self._f / (self.net_reproduction_rate - 1.0)

        # State
        self._time = 0.0
        self._genotype_vector = self._initial_state()
        self._initial_total = float(np.sum(self._genotype_vector))
        self._input_overrides: Dict[str, BioSignal] = {}
        self._history: List[Dict[str, float]] = []
        self._outputs: Dict[str, BioSignal] = {}

    def _initial_state(self) -> np.ndarray:
        """Create initial population: wild-type males + females + transgenic release."""
        gv = np.zeros(N_GENOTYPES)
        # Wild-type male: (ab, ef, cd, cd) — wt Y, wt X, wt autosomes
        wt_male = _find_genotype(("ab", "ef", "cd", "cd"))
        # Wild-type female: (ef, ef, cd, cd) — two wt X, wt autosomes
        wt_female = _find_genotype(("ef", "ef", "cd", "cd"))

        # Split initial population 50:50 between males and females
        half_pop = self.initial_population * 0.5
        if wt_male is not None:
            gv[wt_male] = half_pop
        if wt_female is not None:
            gv[wt_female] = half_pop

        if self.release_size > 0:
            # Release transgenic males: Y2 (AB) with X1 (ef), fully functional
            # autosomal shredder: A2 (CD)
            release_geno = _find_genotype(("AB", "ef", "CD", "CD"))
            if release_geno is not None:
                gv[release_geno] = self.initial_population * self.release_size

        return gv

    @staticmethod
    def _scalar_input_spec(unit: str, description: str) -> SignalSpec:
        return SignalSpec.scalar(
            dtype="float64",
            accepted_profiles=(
                AcceptedSignalProfile(signal_type="scalar", dtype="float64",
                                     accepted_units=(unit,), description=description),
                AcceptedSignalProfile(signal_type="record", schema={"payload": "json"},
                                     description=description),
            ),
            description=description,
        )

    def inputs(self) -> dict[str, SignalSpec]:
        return {
            "net_reproduction_rate": self._scalar_input_spec("dimensionless", "Net reproduction rate."),
            "juvenile_survival": self._scalar_input_spec("fraction", "Juvenile survival probability."),
            "initial_population": self._scalar_input_spec("population", "Normalized initial population size."),
            "release_size": self._scalar_input_spec("fraction", "Transgenic release as fraction of initial population."),
            "homing_efficiency": self._scalar_input_spec("fraction", "Homing efficiency."),
            "editing_efficiency": self._scalar_input_spec("fraction", "Editing efficiency."),
            "shredding_efficiency": self._scalar_input_spec("fraction", "X-shredding efficiency."),
            "copy_mutation_rate": self._scalar_input_spec("fraction", "Copying-error mutation rate."),
            "background_mutation_rate": self._scalar_input_spec("fraction", "Background mutation rate."),
            "editing_resistance_rate": self._scalar_input_spec("fraction", "Editing resistance rate."),
            "shredding_resistance_rate": self._scalar_input_spec("fraction", "Shredding resistance rate."),
            "homing_resistance_rate": self._scalar_input_spec("fraction", "Homing resistance rate."),
            "fitness_cost_cas9": self._scalar_input_spec("fraction", "Cas9 expression fitness cost."),
            "fitness_cost_grna": self._scalar_input_spec("fraction", "gRNA expression fitness cost."),
            "fitness_cost_shredder": self._scalar_input_spec("fraction", "Shredder expression fitness cost."),
            "fitness_cost_nuclease": self._scalar_input_spec("fraction", "Nuclease activity fitness cost."),
            "fitness_cost_shredder_activity": self._scalar_input_spec("fraction", "Shredder activity fitness cost."),
            "fitness_cost_edited_female": self._scalar_input_spec("fraction", "Female edited-target fitness cost."),
            "fitness_cost_edited_male": self._scalar_input_spec("fraction", "Male edited-target fitness cost."),
            "dominance_editing": self._scalar_input_spec("dimensionless", "Dominance coefficient for female editing."),
            "dominance_shredder": self._scalar_input_spec("dimensionless", "Dominance coefficient for shredder activity."),
            "dominance_shredder_gamete": self._scalar_input_spec("dimensionless", "Dominance coefficient for shredder in gametes."),
            "cas9_cofactor": self._scalar_input_spec("dimensionless", "Cas9 cofactor for shredding."),
            "recombination_rate": self._scalar_input_spec("fraction", "X-linked recombination rate."),
        }

    def set_inputs(self, inputs: dict[str, BioSignal]) -> None:
        self._input_overrides = dict(inputs or {})
        self._apply_input_overrides(reset_initial_state=self._time <= 0.0 and not self._history)

    def _input_number(self, name: str) -> float | None:
        signal = self._input_overrides.get(name)
        if signal is None:
            return None
        value = signal.value
        if isinstance(value, dict):
            if "payload" in value:
                value = value["payload"]
            elif "value" in value:
                value = value["value"]
        try:
            return float(value)
        except (TypeError, ValueError):
            return None

    def _apply_input_overrides(self, *, reset_initial_state: bool) -> None:
        rebuild_matrices = False

        # Genetic parameters (stored in self._params)
        for param_name in list(self._params.keys()):
            value = self._input_number(param_name)
            if value is not None:
                self._params[param_name] = value
                rebuild_matrices = True

        if rebuild_matrices:
            self._rebuild_matrices()

        # Ecology parameters
        for attr in ("net_reproduction_rate", "juvenile_survival"):
            value = self._input_number(attr)
            if value is not None and value > 0:
                setattr(self, attr, value)
                self._f = (self.net_reproduction_rate * 2.0) / self.juvenile_survival
                self._alpha = self.initial_population * self._f / (self.net_reproduction_rate - 1.0)

        # Initial-condition parameters (only apply before simulation starts)
        for attr in ("initial_population", "release_size"):
            value = self._input_number(attr)
            if value is not None and value >= 0:
                setattr(self, attr, value)
                if reset_initial_state:
                    self._f = (self.net_reproduction_rate * 2.0) / self.juvenile_survival
                    self._alpha = self.initial_population * self._f / (self.net_reproduction_rate - 1.0)
                    self._genotype_vector = self._initial_state()
                    self._initial_total = float(np.sum(self._genotype_vector))

    def _rebuild_matrices(self) -> None:
        """Rebuild all process matrices from current self._params."""
        sperm_mat, egg_mat = _build_gamete_matrices(self._params)
        self._matrices = {
            "selection": _build_selection_vector(self._params),
            "homing": _build_homing_matrix(self._params),
            "editing": _build_editing_matrix(self._params),
            "recombination": _build_recombination_matrix(self._params),
            "mutation": _build_mutation_matrix(self._params),
            "sperm": sperm_mat,
            "egg": egg_mat,
        }

    def outputs(self) -> dict[str, SignalSpec]:
        return {
            'population_state': SignalSpec.record(
                schema={'total_adults': 'json', 'adult_females': 'json', 'adult_males': 'json'},
                emitted_unit='individuals',
                description='Adult mosquito population partitioned into females and males.',
            ),
            'gene_drive_metrics': SignalSpec.record(
                schema={'drive_frequency': 'json', 'resistance_frequency': 'json',
                        'male_fraction': 'json', 'suppression_ratio': 'json'},
                emitted_unit='fraction',
                description='Eco-genetic metrics for drive spread, resistance, and suppression.',
            ),
            'visualisation_payload': SignalSpec.record(schema={'payload': 'json'}, description='Internal history payload for the sibling visualisation model.'),
        }

    def setup(self, config: Optional[Dict[str, Any]] = None) -> None:
        self.reset()

    def reset(self) -> None:
        self._time = 0.0
        self._genotype_vector = self._initial_state()
        self._initial_total = float(np.sum(self._genotype_vector))
        self._history = []
        self._outputs = {}

    def advance_window(self, start: float, end: float, inputs: dict[str, BioSignal] | None = None) -> None:
        if inputs:
            self.set_inputs(inputs)
        else:
            self._apply_input_overrides(reset_initial_state=False)

        t = float(end)
        if t <= self._time:
            return

        # Each time unit = one generation
        n_generations = max(1, int(round(t - self._time)))
        for gen in range(n_generations):
            self._step_generation()
            self._time += 1.0
            self._record_state(self._time)

        self._publish_outputs(self._time)

    def get_outputs(self) -> Dict[str, BioSignal]:
        return dict(self._outputs)

    def visualize(self) -> Optional[List["VisualSpec"]]:
        return None

    def _step_generation(self) -> None:
        """One generation of the gene-drive simulation.

        Port of Julia timecourse loop (lines 1107-1165).
        """
        gv = self._genotype_vector
        M = self._matrices

        # 1. Mutation
        gv = M["mutation"].T @ gv
        # 2. Homing
        gv = M["homing"].T @ gv
        # 3. Editing
        gv = M["editing"].T @ gv
        # 4. Recombination
        gv = M["recombination"].T @ gv

        # 5. Gamete production
        sperm = M["sperm"].T @ gv
        total_sperm = np.sum(sperm)
        if total_sperm > 0:
            sperm = sperm / total_sperm

        eggs = M["egg"].T @ gv
        eggs = eggs * self._f

        # 6. Zygote formation
        gv = create_zygote_vector(sperm, eggs)

        # 7. Density-dependent survival (Beverton-Holt)
        total_zygotes = np.sum(gv)
        if total_zygotes > 0:
            gv = gv * (self.juvenile_survival * self._alpha / (self._alpha + total_zygotes))

        # 8. Selection
        gv = M["selection"] * gv

        self._genotype_vector = gv

    def _compute_ecology(self) -> Dict[str, float]:
        """Derive ecology observables from genotype state."""
        gv = self._genotype_vector
        total_females = float(np.sum(gv[_FEMALE_INDICES]))
        total_males = float(np.sum(gv[_MALE_INDICES]))
        total_adults = total_females + total_males

        # Drive frequency: fraction of Y-linked transgenics among males
        # Y2=AB, Y3=Aβ, Y4=αB, Y5=αβ all carry some transgenic component
        drive_males = 0.0
        for idx in _MALE_INDICES:
            g = GENOTYPES[idx]
            sex1 = g[0]
            if sex1 != "ab":  # anything other than wild-type Y
                drive_males += gv[idx]
        drive_frequency = drive_males / total_males if total_males > 0 else 0.0

        # Resistance frequency: fraction of autosomal r3 (homing-resistant) alleles
        resistance_count = 0.0
        total_autosome_alleles = 0.0
        for idx in range(N_GENOTYPES):
            g = GENOTYPES[idx]
            val = gv[idx]
            if val <= 0:
                continue
            total_autosome_alleles += 2.0 * val
            if g[2] == "r3":
                resistance_count += val
            if g[3] == "r3":
                resistance_count += val
        resistance_frequency = resistance_count / total_autosome_alleles if total_autosome_alleles > 0 else 0.0

        male_fraction = total_males / total_adults if total_adults > 0 else 0.5
        suppression_ratio = 1.0 - (total_adults / self._initial_total if self._initial_total > 0 else 0.0)

        return {
            "total_adults": total_adults,
            "adult_females": total_females,
            "adult_males": total_males,
            "drive_frequency": drive_frequency,
            "resistance_frequency": resistance_frequency,
            "male_fraction": male_fraction,
            "suppression_ratio": suppression_ratio,
        }

    def _record_state(self, t: float) -> None:
        eco = self._compute_ecology()
        eco["t"] = float(t)
        self._history.append(eco)

    def _publish_outputs(self, t: float) -> None:
        latest = self._history[-1] if self._history else self._compute_ecology()
        source_name = getattr(self, "_world_name", self.__class__.__name__)
        self._outputs = {
            "population_state": _make_signal(
                source=source_name, name="population_state",
                value={
                    "total_adults": float(latest["total_adults"]),
                    "adult_females": float(latest["adult_females"]),
                    "adult_males": float(latest["adult_males"]),
                },
                emitted_at=float(t),
                spec=self.outputs().get("population_state"),
            ),
            "gene_drive_metrics": _make_signal(
                source=source_name, name="gene_drive_metrics",
                value={
                    "drive_frequency": float(latest["drive_frequency"]),
                    "resistance_frequency": float(latest["resistance_frequency"]),
                    "male_fraction": float(latest["male_fraction"]),
                    "suppression_ratio": float(latest["suppression_ratio"]),
                },
                emitted_at=float(t),
                spec=self.outputs().get("gene_drive_metrics"),
            ),
            "visualisation_payload": _make_signal(
                source=source_name,
                name="visualisation_payload",
                value={"payload": self._visualisation_payload()},
                emitted_at=float(t),
                spec=self.outputs().get("visualisation_payload"),
            ),
        }

    def _visualisation_payload(self) -> Dict[str, Any]:
        return {"initial_total": self._initial_total, "history": list(self._history)}

    def _population_visual(self) -> "VisualSpec":
        return {
            "render": "timeseries",
            "description": "Adult mosquito abundance through time, including sex-ratio distortion from the drive system.",
            "data": {
                "title": "Adult Population Suppression",
                "series": [
                    {"name": "Total adults", "points": [[p["t"], p["total_adults"]] for p in self._history]},
                    {"name": "Adult females", "points": [[p["t"], p["adult_females"]] for p in self._history]},
                    {"name": "Adult males", "points": [[p["t"], p["adult_males"]] for p in self._history]},
                ],
            },
        }

    def _genetics_visual(self) -> "VisualSpec":
        return {
            "render": "timeseries",
            "description": "Drive spread, resistance emergence, and male bias for the faithful Geci2022 eco-genetic model.",
            "data": {
                "title": "Gene-Drive Metrics",
                "series": [
                    {"name": "Drive frequency", "points": [[p["t"], p["drive_frequency"]] for p in self._history]},
                    {"name": "Resistance frequency", "points": [[p["t"], p["resistance_frequency"]] for p in self._history]},
                    {"name": "Male fraction", "points": [[p["t"], p["male_fraction"]] for p in self._history]},
                    {"name": "Suppression ratio", "points": [[p["t"], p["suppression_ratio"]] for p in self._history]},
                ],
            },
        }

    def _summary_visual(self) -> "VisualSpec":
        latest = self._history[-1]
        peak_drive = max(p["drive_frequency"] for p in self._history)
        peak_resistance = max(p["resistance_frequency"] for p in self._history)
        min_population = min(p["total_adults"] for p in self._history)
        return {
            "render": "table",
            "description": "Final and extremal eco-genetic metrics for the faithful Geci2022 model.",
            "data": {
                "title": "Geci2022 Summary",
                "columns": ["Metric", "Value"],
                "rows": [
                    ["Initial total", f"{self._initial_total:.6g}"],
                    ["Final total adults", f"{latest['total_adults']:.6g}"],
                    ["Minimum total adults", f"{min_population:.6g}"],
                    ["Final females", f"{latest['adult_females']:.6g}"],
                    ["Final males", f"{latest['adult_males']:.6g}"],
                    ["Peak drive frequency", f"{peak_drive:.6g}"],
                    ["Peak resistance frequency", f"{peak_resistance:.6g}"],
                    ["Final male fraction", f"{latest['male_fraction']:.6g}"],
                    ["Final suppression ratio", f"{latest['suppression_ratio']:.6g}"],
                ],
            },
        }
