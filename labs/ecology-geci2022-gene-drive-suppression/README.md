# Geci2022 Gene-Drive Suppression Lab 1.1.0

Explore deterministic genotype and population dynamics with 24 controls for reproduction, survival, initial population, release ratio, molecular efficiencies, fitness costs and resistance. The model tracks 1,071 genotypes and 66 gametes. Outputs describe the specified mathematical scenario; they are not a field prediction or deployment recommendation.

The source is the Julia implementation distributed as [BioModels MODEL2301120001](https://www.ebi.ac.uk/biomodels/MODEL2301120001), associated with Geci, Willis and Burt, [PLOS Genetics (2022)](https://doi.org/10.1371/journal.pgen.1010550). BioModels identifies the record as non-curated. The exact bundled source and native verification artifacts are pinned by checksums in `verification/`.

## Read the outputs correctly

Population values are normalized continuous quantities, not counts of individual animals. The existing wrapper initializes a wild-type total I, split equally between sexes, and adds a transgenic male release of I × release_size. Its retained density-dependence scale implies a wild-type equilibrium total of 2I. Thus the no-release example starts below equilibrium and grows toward 2I. This initialization differs from the upstream convenience function, which starts with I of each sex. Both the equations and the wrapper's existing initialization are preserved in this repair.

`drive_frequency` is the fraction of males carrying a non-wild-type Y chromosome, including dysfunctional transgenic variants. `resistance_frequency` counts autosomal r3 homing-resistance alleles; it does not summarize all resistance mechanisms. `suppression_ratio` retains the legacy formula 1 − N(t)/N(0), with the initial total including release. It can be negative. It is not an effect relative to a matched no-release control. The plots and summary state these definitions.

The example runs for 100 whole generations in 10-generation communication windows. It includes generation zero and every generation through the final state. The final presenter consumes that same full history. Rate and initial-condition controls are fixed after the run starts; use separate runs for separate scenarios.

## Verification and limitations

44 tests pass with biosimulant==0.0.34 and numpy==1.26.4. Four 100-generation scenarios match native Julia across all genotype abundances and seven process matrices. The largest absolute trajectory difference was 1.23e-15. Native scenarios use exactly matched parameters and initial vectors; one adopts upstream baseline parameters, but this is not a claim to reproduce a published figure. See MTS.md and verification/comparison.json.

Invalid probabilities and fractional-generation windows fail explicitly. A successful calculation does not establish empirical ecological accuracy, stochastic extinction, spatial spread or release efficacy. The Lab implements a single deterministic population, not every spatial or multiple-release workflow in the Julia file. Handwritten scientific logic requires manual review before promotion. The repaired version has not been publicly released or assigned a new managed Experiment Passport.

Run locally: `biosimulant labs run . --results-file results.json --report-file report.html`.
Tests: install requirements-test.txt and run `OPENBLAS_NUM_THREADS=1 python -m pytest . -q`.
