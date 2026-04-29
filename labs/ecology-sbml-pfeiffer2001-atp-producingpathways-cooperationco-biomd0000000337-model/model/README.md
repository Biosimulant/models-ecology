# Pfeiffer2001 ATP-Pathway Cooperation-Competition Model

**Upstream reference**: [BioModels BIOMD0000000337](https://www.ebi.ac.uk/biomodels/BIOMD0000000337)  
**Original asset type**: SBML  
**Package standard**: `other`

## Scientific Scope

This package keeps the published SBML dynamics and surfaces the ecological
observables that matter for interpretation:

- Shared substrate resource (`S`)
- Population using strategy `N1`
- Population using strategy `N2`
- Relative dominance of the two ATP-producing strategies

## Visuals

- Resource and biomass time series
- Strategy-fraction time series
- Summary table for depletion and competitive outcome

## Packaging Notes

- Units are propagated from SBML only when the source model declares them.
- Time-series visuals are generated from the accumulated simulation history, not
  a one-point snapshot.
