# models-ecology

> Storage-only repo: each former root model now lives in `labs/<slug>/models/core/` and is wrapped by
> `labs/<slug>/lab.yaml`. This repo has no repo-level import catalog and no composed labs at the root.

Focused research-grade ecology and population-dynamics model pack for the `biosim` platform.

## What's Inside

### Models (5 packages)

#### Native Ecology Model
- `ecology-lotka-volterra-system`
  Canonical deterministic predator-prey dynamics with RK4 integration, a phase portrait, summary diagnostics, and invariant-drift auditing.

#### Publication-Grounded Ecology Models
- `ecology-sbml-geci2022-model2301120001-model`
  Curated deterministic gene-drive population suppression model inspired by Geci et al. (2022).
- `ecology-sbml-leibovich2022-multispecies-eco-competition-descr-model2212080001-model`
  Curated stochastic multispecies competition model with immigration and demographic noise inspired by Leibovich et al. (2022).
- `ecology-sbml-pfeiffer2001-atp-producingpathways-cooperationco-biomd0000000337-model`
  Curated pathway-competition model tracking resource use and pathway dominance.
- `ecology-sbml-turner2015-human-mosquito-elp-model-biomd0000000922-model`
  Curated mosquito early-life-stage population dynamics model tracking eggs, larvae, and pupae.

### Spaces (1 package)

- `ecology-predator-prey`
  Canonical Lotka-Volterra rabbit/fox dynamics using the native ecology model.

## Principles

- Keep only ecology and population-dynamics models.
- Use visuals that match the modeled observables.
- Prefer explicit assumptions over hidden heuristics.
- Treat imported source material as publication-grounded only when the packaging is scientifically coherent and visually interpretable.

## Layout

```text
models-ecology/
├── models/<model-slug>/     # One model package per folder, each with model.yaml
├── labs/<space-slug>/     # Composed simulation spaces with lab.yaml
├── scripts/                 # Validation scripts
├── templates/model-pack/    # Starter template for new model packs
├── docs/                    # Governance documentation
└── .github/workflows/       # CI/CD pipeline
```

## Running Validation

```bash
python scripts/validate_manifests.py
python scripts/check_entrypoints.py
```

## Notes

- This repository is intentionally small and thematic.
- Non-ecology imports should live in domain-appropriate repositories instead of being mixed into this catalog.
