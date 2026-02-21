# models-ecology

Curated collection of **ecology** and **population dynamics** simulation models for the **biosim** platform. This repository contains computational models of ecosystem interactions, predator-prey dynamics, population regulation, and environmental coupling, plus related biochemical and cellular models.

## What's Inside

### Models (30 packages)

Each model is a self-contained simulation component with a `model.yaml` manifest.

**Ecology & Population Dynamics** — ecosystem modeling, species interactions, and environmental regulation:

#### Core Ecology Models (Custom-Built)
- `ecology-abiotic-environment` — Broadcasts environmental conditions (temperature, water, food, sunlight)
- `ecology-organism-population` — Population dynamics with birth, death, and predation
- `ecology-predator-prey-interaction` — Predation rates and functional response
- `ecology-population-monitor` — Population size tracking over time
- `ecology-phase-space-monitor` — Predator vs prey phase-space visualization
- `ecology-population-metrics` — Ecosystem summary statistics

#### Ecological & Biological Systems Models (SBML)
- `ecology-sbml-leibovich2022-multispecies-eco-competition-descr` — Multi-species ecological competition
- `ecology-sbml-edelstein1996-epsp-ach-species` — Acetylcholine species dynamics

#### Blood Coagulation & Hemodynamics
- `ecology-sbml-hansen2019-nine-species-reduced-model-of-blood-c` — Nine-species blood coagulation model
- `ecology-sbml-hansen2019-seven-species-reduced-model-of-blood` — Seven-species reduced coagulation model
- `ecology-sbml-makin2013-blood-coagulation-cascade-model` — Blood coagulation cascade
- `ecology-sbml-mitrophanov2014-extended-hockin-blood-coagulatio` — Extended Hockin coagulation model

#### Cellular & Molecular Systems
- `ecology-sbml-novak2001-fissionyeast-cellcycle` — Fission yeast cell cycle
- `ecology-sbml-obeyesekere1999-cellcycle` — Cell cycle dynamics
- `ecology-sbml-coggins2014-cxcl12-dependent-recruitment-of-beta` — CXCL12-dependent cell recruitment

#### Pharmacokinetics & Distribution
- `ecology-sbml-e2-pbpk` — PBPK model for E2 distribution
- `ecology-sbml-mouse-iron-distribution-adequate-iron-diet-trace` — Mouse iron distribution modeling

#### Protein Aggregation & Signaling
- `ecology-sbml-masel2000-drugs-to-stop-prion-aggregates-and-oth` — Prion aggregation and drug intervention
- `ecology-sbml-nik-dependent-p100-processing-into-p52-with-relb` — NIK-dependent NF-κB processing
- `ecology-sbml-geci2022` — Genetically encoded calcium indicators

**Note:** This repository contains 30 models total, including 6 custom-built ecology models and 24 SBML models from various biological domains. For a complete list, see the `models/` directory.

## Layout

```
models-ecology/
├── models/<model-slug>/     # One model package per folder, each with model.yaml
├── libs/                    # Shared helper code for curated models
├── templates/model-pack/    # Starter template for new model packs
├── scripts/                 # Manifest and entrypoint validation scripts
├── docs/                    # Governance documentation
└── .github/workflows/       # CI/CD pipeline
```

## How It Works

### Model Interface

Every model implements the `biosim.BioModule` interface:

- **`inputs()`** — declares named input signals the module consumes
- **`outputs()`** — declares named output signals the module produces
- **`advance_to(t)`** — advances the model's internal state to time `t`

Custom ecology models include Python source under `src/` and can be wired together via `space.yaml` without additional code.

### Model Standards

Models in this repository include:
- **Custom ecology models**: Native Python implementations with composable interfaces
- **SBML models**: Use SBML format with tellurium runtime for execution
- All provide `state` or domain-specific outputs for monitoring
- Support configurable timesteps via `min_dt` parameter

### Running Models

Models are loaded and executed by the `biosim-platform`. The platform reads `model.yaml`, instantiates the model from its entrypoint, and runs the simulation loop at the configured timestep for the specified duration.

Models can be wired together into composed simulations (spaces) for multi-scale ecological or biological modeling.

## Getting Started

### Prerequisites

- Python 3.11+
- `biosim` framework

### Install biosim

```bash
pip install "biosim @ git+https://github.com/BioSimulant/biosim.git@main"
```

### Create a New Model

1. Copy `templates/model-pack/` to `models/<your-model-slug>/`
2. Edit `model.yaml` with metadata, entrypoint, and pinned dependencies
3. Implement your module (subclass `biosim.BioModule`)
4. Add ecology-specific tags and categorization
5. Validate: `python scripts/validate_manifests.py && python scripts/check_entrypoints.py`

### Creating Ecology Spaces

Example predator-prey space configuration:

```yaml
wiring:
  - from: environment.temperature
    to: [prey_pop.temperature, predator_pop.temperature]
  - from: prey_pop.population
    to: [predator_prey_interaction.prey_count]
  - from: predator_pop.population
    to: [predator_prey_interaction.predator_count]
  - from: predator_prey_interaction.predation_rate
    to: [prey_pop.mortality, predator_pop.food_intake]
```

## Linking in biosim-platform

- Models can be linked with explicit paths:
  - `models/ecology-predator-prey-interaction/model.yaml`
- Ecology models can be composed with other domain models in multi-scale simulations

## External Repos

External authors can keep models in independent repositories and link them directly in `biosim-platform`. This repository is curated, not exclusive.

## Validation & CI

Three scripts enforce repository integrity on every push:

| Script | Purpose |
|--------|---------|
| `scripts/validate_manifests.py` | Schema validation for all model.yaml files |
| `scripts/check_entrypoints.py` | Verifies Python entrypoints are importable and callable |
| `scripts/check_public_boundary.sh` | Prevents business-sensitive content in this public repo |

The CI pipeline (`.github/workflows/ci.yml`) runs: **secret scan** → **manifest validation** → **smoke sandbox** (Docker).

## Contributing

- All dependencies must use exact version pinning (`==`)
- Model slugs use kebab-case with domain prefix (`ecology-` or `ecology-sbml-`)
- Custom modules must follow the `biosim.BioModule` interface
- SBML models use tellurium runtime for execution
- Pre-commit hooks enforce trailing whitespace, EOF newlines, YAML syntax, and secret detection
- See [docs/PUBLIC_INTERNAL_BOUNDARY.md](docs/PUBLIC_INTERNAL_BOUNDARY.md) for content policy

## Domain-Specific Notes

**Ecology & Population Dynamics Focus Areas:**
- **Predator-Prey Interactions**: Classic Lotka-Volterra dynamics, functional responses
- **Environmental Coupling**: Abiotic factors affecting population dynamics
- **Population Monitoring**: Time series, phase space, and summary metrics
- **Multi-Species Competition**: Species coexistence and competitive exclusion

**Other Biological Systems** (included in this repository):
- Blood coagulation cascades
- Cell cycle regulation
- Pharmacokinetic distribution
- Protein aggregation dynamics

## License

This repository is dual-licensed:

- **Code** (scripts, templates, Python modules): Apache-2.0 (`LICENSE-CODE.txt`)
- **Model/content** (manifests, docs, wiring/config): CC BY 4.0 (`LICENSE-CONTENT.txt`)

Attribution guidance: `ATTRIBUTION.md`
