# AGENTS.md

Instructions for AI agents working in `models-ecology`.

## Repository Purpose

This repository is a public, tightly scoped ecology and population-dynamics pack for `biosim`.
Only research-grade ecology models and ecology-focused composed spaces should remain here.

## Repository Structure

```text
models-ecology/
├── models/          # 5 model packages
├── spaces/          # 1 composed ecology space
├── scripts/         # Validation scripts
├── templates/       # Starter template
├── docs/            # Governance docs
└── .github/         # CI workflows
```

## Working Rules

- Keep the repository ecology-only.
- Remove stale catalog references when models are deleted.
- Prefer explicit scientific assumptions in manifests and READMEs.
- Visual outputs must reflect the real scientific quantities exposed by the model.

## Validation

Before merging:

```bash
python scripts/validate_manifests.py
python scripts/check_entrypoints.py
```

## Public Boundary

This repository is public. Do not add private operational or business material.
