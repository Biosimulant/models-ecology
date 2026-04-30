# Leibovich2022 Multispecies Competition Lab

This lab runs a stochastic multispecies community model. It asks: when several species compete for the same resources and new individuals arrive by immigration, what determines how many species coexist and how abundant each one is?

The model uses tau-leaping to approximate Gillespie dynamics. Each species has a birth propensity driven by its current abundance and a flat immigration rate, and a death propensity that increases with crowding. The crowding term mixes intraspecific and interspecific competition through a single overlap parameter. Higher overlap means species interact more strongly, which tends to reduce total abundance and can lead to competitive exclusion.

This is a faithful Python port of the upstream MultiLV Gillespie model from Leibovich et al. (2022), published as [BioModels MODEL2212080001](https://www.ebi.ac.uk/biomodels/MODEL2212080001). The upstream asset is Python Gillespie code, not SBML.

## What You'll See

The lab opens as a canvas with one community model node and a run-results panel. After running, you will see three visualizations: species-resolved abundance trajectories, community-level metrics over time, and a summary table.

## How to Read the Visualizations

The abundance trajectory plot shows each species as a separate line. The x-axis is time and the y-axis is individual count. In a typical run with six species, you will see stochastic fluctuations around a shared carrying capacity. Some species may go extinct if competition is strong and immigration is low.

The community metrics plot tracks three quantities:

- **Total abundance**: sum of all species. This should settle near the carrying capacity.
- **Species richness**: how many species have at least one individual. Immigration prevents extinction, so higher immigration rates support higher richness.
- **Shannon diversity**: an information-theoretic measure of evenness. Higher values mean more equal abundances across species.

The summary table gives the parameter settings, final and extremal abundances, peak richness, and the dominant species.

## What This Lab Contains

- `lab.yaml` describes the lab and exposes its outputs.
- `wiring-layout.json` places the model on the canvas.
- `model/model.yaml` describes the model package, parameters, and ports.
- `model/src/leibovich2022_competition.py` contains the stochastic simulation and propensity formulas.
- `model/tests/` checks output accumulation, diversity metrics, propensity formulas, immigration uniformity, and competition effects.
- `model/upstream/` contains the original Python Gillespie code from BioModels.

## Inputs

All parameters are set through `model/model.yaml` init_kwargs. The six runtime-tunable parameters below are also exposed as input ports for workflow wiring. `species_count` and `rng_seed` are structural and not available as input ports.

- `species_count` (`count`): number of competing species (default 6).
- `carrying_capacity` (`individuals`): habitat carrying capacity (default 100).
- `birth_rate` (`1/time`): per-capita birth rate (default 2.0).
- `death_rate` (`1/time`): per-capita baseline death rate (default 1.0).
- `competition_overlap` (`dimensionless`): fraction of competition that is interspecific (default 0.2). At 0, species only compete with themselves. At 1, all species compete equally.
- `immigration_rate` (`individuals/time`): constant immigration rate per species (default 0.1).
- `initial_abundance` (`individuals`): starting abundance per species (default 50).
- `rng_seed` (`integer`): random number generator seed for reproducibility.

## Outputs

- `community_state` (`individuals`): species-resolved abundances and total abundance.
- `diversity_metrics`: species richness, Shannon diversity, evenness, and dominant species identity.

## Recreate and Run with the Biosim CLI

From this lab folder:

```bash
cd /path/to/models-ecology/labs/ecology-leibovich2022-multispecies-competition
mkdir -p dist
python -m biosim pack build . --out dist/leibovich2022-competition.bsilab
python -m biosim pack run dist/leibovich2022-competition.bsilab
```

If you are working from this monorepo without installing `biosim`, use the local package environment instead:

```bash
mkdir -p dist
/path/to/bsim-active/biosim/.venv/bin/python -m biosim pack build . --out dist/leibovich2022-competition.bsilab
/path/to/bsim-active/biosim/.venv/bin/python -m biosim pack run dist/leibovich2022-competition.bsilab
```

You can also validate the package before running:

```bash
python -m biosim pack validate dist/leibovich2022-competition.bsilab
```

## Run in the Desktop App

1. Open Biosimulant Desktop.
2. Go to Projects or Labs.
3. Choose the option to open or import an existing lab.
4. Select this folder's `lab.yaml`.
5. Open the lab and press Run.

The right side of the app should show the run result and the visualizations.

## How to Edit It

For scenario changes, start with `model/model.yaml` and `lab.yaml`.

- Change `runtime.duration` in `lab.yaml` for a longer or shorter run.
- Change `runtime.communication_step` if you want more or fewer reported points.
- Change defaults in `model/model.yaml` for species count, carrying capacity, birth/death rates, competition overlap, or immigration rate.
- Try `competition_overlap: 0.0` for pure intraspecific competition, or `competition_overlap: 1.0` for full interspecific competition.
- Try higher `immigration_rate` values to see how immigration rescues species from extinction.

Edit `model/src/leibovich2022_competition.py` only if you are changing the model mechanics. Good code-level edits include adding a new diversity metric, changing the density-dependence formula, or adding a new visualization. If you only want a different scenario, prefer changing parameters rather than changing the Python code.
