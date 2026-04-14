# Turner2015 Mosquito Life-Stage Ecology Model

**Upstream reference**: [BioModels BIOMD0000000922](https://www.ebi.ac.uk/biomodels/BIOMD0000000922)  
**Original asset type**: SBML  
**Package standard**: `other`

## Scientific Scope

This package keeps the original SBML dynamics but interprets them through the
ecology observables that matter for this model:

- Eggs
- Larvae
- Pupae
- Total immature mosquito abundance
- Stage composition through time

## Visuals

- Life-stage abundance time series
- Life-stage fraction time series
- Summary table for immature population size

## Packaging Notes

- Units are emitted only when declared by the SBML source.
- Visual payloads are built from accumulated trajectories rather than a single
  endpoint sample.
