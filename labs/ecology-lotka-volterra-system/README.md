# Lotka–Volterra System Lab 1.1.0

Explore how prey growth, predation, predator mortality and reproduction shape idealized predator–prey cycles. This is an interpretable ecological baseline, not a species-calibrated population forecast.

The six public controls set initial prey and predator populations and the four rate coefficients. Defaults are prey 10, predator 5, alpha 1.1/day, beta 0.4/(count day), gamma 0.4/day and delta 0.1/(count day). Time is in days; population counts are continuous quantities, not rounded individual animals.

The default run spans 250 days at 0.01-day communication intervals. Outputs include both final populations, the complete trajectory from time zero through the final state, population and phase plots, a summary and a conserved-quantity diagnostic. A near-zero numerical threshold does not establish biological extinction. At zero populations, or if rates change during a run, the invariant audit is marked not applicable.

Assumptions: closed, homogeneous populations; unlimited prey resources; mass-action encounters; constant default rates; no carrying capacity, age structure, seasons, migration, demographic noise or fitted species observations. Conclusions apply to these equations and inputs. Handwritten scientific logic remains subject to manual scientific review.

Local verification against biosimulant 0.0.34: 54 checks pass, including an independent SciPy DOP853 comparison over the complete default trajectory, fourth-order convergence, exact equilibrium and boundary solutions, invalid-input rejection, reproducibility and final visual coverage. See MTS.md and verification.json for numerical tolerances and measured evidence. This repair has not been published or given a new managed Experiment Passport.

Run locally with `biosimulant labs run . --results-file results.json --report-file report.html`. The model runtime is exactly pinned in each embedded model manifest. For verification, install requirements-test.txt and run `python -m pytest . -q`.
