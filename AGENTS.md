# AGENTS.md

Instructions for AI agents working in `models-ecology`.

## Repository Purpose

This repository is a public, tightly scoped ecology and population-dynamics pack
for Biosimulant/BioSim. Keep it focused on runnable ecology labs with clear
scientific assumptions, explicit units, and visual outputs that match the
modeled quantities.

The current layout is lab-first: every runnable item lives under `labs/<slug>/`
and contains its own `lab.yaml`, `wiring-layout.json`, embedded model folder,
tests, and lab-level README when appropriate. Do not recreate the old root
`models/` or `spaces/` layout.

## Repository Structure

```text
models-ecology/
├── labs/
│   ├── ecology-lotka-volterra-system/
│   │   ├── README.md
│   │   ├── lab.yaml
│   │   ├── wiring-layout.json
│   │   └── model/
│   │       ├── model.yaml
│   │       ├── src/
│   │       └── tests/
│   ├── ecology-rosenzweig-macarthur-predator-prey-system/
│   └── ecology-sbml-.../
├── docs/
├── .github/
├── README.md
├── LICENSE
└── ATTRIBUTION.md
```

There may be legacy references to `templates/`, `models/`, or `spaces/` in old
docs or scripts. Treat those as stale unless the path exists and is actively
used by the current lab workflow.

## Current Native Labs

- `ecology-lotka-volterra-system`
  - Baseline deterministic Lotka-Volterra predator-prey system.
  - Keep the equations unchanged:
    - `dN/dt = alpha*N - beta*N*P`
    - `dP/dt = delta*N*P - gamma*P`
  - Keep this lab simple: no carrying capacity, seasonality, disease,
    migration, spatial patches, or age structure.

- `ecology-rosenzweig-macarthur-predator-prey-system`
  - Applied predator-prey lab using logistic prey growth and Holling type II
    predation.
  - This is the place for carrying capacity, food limits, seasonality, noise,
    disease, migration, patch structure, stage bookkeeping, risk summaries, and
    richer ecology diagnostics.

## Working Rules

- Keep the repository ecology-only.
- Prefer generic scientific names such as `prey` and `predator` in ports and
  manifests. Species names can be display labels, not wire-facing contracts.
- Use explicit units for inputs and outputs.
- Keep lab-level `io` mappings in sync with model ports so labs can be used as
  sublabs.
- Use visuals that explain the modeled system: trajectories, phase portraits,
  functional responses, resource/patch tables, risk summaries, and diagnostics
  that are actually tied to the equations.
- Do not silently change the scientific identity of a lab. If a model needs new
  assumptions, create or use the appropriate lab rather than stretching a
  baseline model beyond its scope.
- Add or update tests when changing model behavior, ports, manifests, or
  visualization contracts.

## Validation

When the validation scripts are present, run:

```bash
python scripts/validate_manifests.py
python scripts/check_entrypoints.py
```

For focused model tests, run the relevant lab test folders, for example:

```bash
python -m pytest labs/ecology-lotka-volterra-system/model/tests
python -m pytest labs/ecology-rosenzweig-macarthur-predator-prey-system/model/tests
```

In the local monorepo, the working Python environment is usually:

```bash
/Volumes/dem-ssd/imp/projects/Nitoons/Biosimulant/bsim-active/biosim/.venv/bin/python
```

Use that interpreter if the system `python3` lacks `pytest`, `yaml`, or
`biosim`.

## Public Boundary

This repository is public. Do not add private operational details, credentials,
customer data, production URLs, or business-only notes. Keep examples and docs
portable.
