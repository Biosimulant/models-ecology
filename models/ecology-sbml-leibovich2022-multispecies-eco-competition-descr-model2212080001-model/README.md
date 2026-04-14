# Leibovich2022 Multispecies Competition-Immigration Community Model

**Upstream reference**: [BioModels MODEL2212080001](https://www.ebi.ac.uk/biomodels/MODEL2212080001)  
**Original asset type**: Python Gillespie implementation  
**Package standard**: `other`

## Scientific Scope

This package represents the ecological core of Leibovich et al. (2022): a
multispecies community shaped by immigration, competition overlap, and
demographic noise.

The upstream BioModels record is Python simulation code rather than SBML, so
this package now exposes an explicit stochastic ecology model instead of a
misleading SBML wrapper.

## Observables

- Species-resolved abundance trajectories
- Total community abundance
- Species richness
- Shannon diversity and evenness
- Dominant species identity

## Visuals

- Multi-species abundance time series
- Community metric time series
- Summary table with diversity and dominance diagnostics

## Limits

- This implementation is a curated stochastic competition-immigration model,
  not a verbatim execution of the original research code.
- Environmental forcing and parameter sweeps are not included in the package
  outputs.
