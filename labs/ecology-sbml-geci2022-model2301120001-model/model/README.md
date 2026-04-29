# Geci2022 Eco-Genetic Gene Drive Suppression Model

**Upstream reference**: [BioModels MODEL2301120001](https://www.ebi.ac.uk/biomodels/MODEL2301120001)  
**Original asset type**: Julia source (`new model v9.jl`), not SBML  
**Package standard**: `other`

## Scientific Scope

This package keeps the ecology-facing question from Geci et al. (2022): how a
gene-drive construct changes sex ratio, spreads through a mosquito population,
and suppresses adult abundance through time.

Because the BioModels record is a Julia implementation rather than an SBML
model, Biosimulant now exposes an explicit reduced deterministic
eco-genetic reimplementation instead of a misleading SBML wrapper.

## Observables

- Adult population abundance, split into females and males
- Drive allele frequency
- Resistance frequency
- Population suppression ratio relative to the starting population

## Visuals

- Adult population time series
- Gene-drive metric time series
- Summary table with peak drive, peak resistance, and minimum abundance

## Limits

- This is a reduced ecology package, not a line-by-line execution of the
  upstream Julia code.
- Climate forcing, stage structure, and genotype-resolved bookkeeping are
  intentionally compressed into aggregate eco-genetic observables.
