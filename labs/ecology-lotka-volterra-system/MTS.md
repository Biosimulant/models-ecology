# Modeling technical specification

## Scope and source

A transparent deterministic ecological baseline for exploring the canonical Lotka–Volterra equations. The original Hub package is pinned in source-provenance.json. This revision retains its equations and example parameters; it changes validation, execution timing, trajectory retention and presentation. No source dataset or empirical validation is claimed.

## Equations and units

For prey N and predator P: dN/dt = alpha N − beta N P; dP/dt = delta N P − gamma P. N and P are continuous counts, t is days, alpha and gamma have units day⁻¹, beta and delta have units (count day)⁻¹. Positive-population conserved quantity H = delta N − gamma ln N + beta P − alpha ln P uses population numbers expressed in the declared count unit. H is a numerical diagnostic of the fixed-rate equations, not an ecological measurement.

Defaults: N(0)=10, P(0)=5, alpha=1.1, beta=0.4, gamma=0.4, delta=0.1. No stochastic seed is needed. Initial populations can be overridden before integration and cannot be changed during a run. Rate inputs can change at communication boundaries; doing so disables the conserved-quantity audit. The public example uses constant rates.

## Integration and runtime

Python 3.12 and biosimulant==0.0.34. Classical explicit RK4 uses h=min(integration_step, remaining communication window). Constructor integration_step defaults to 0.1 day; the default Lab communication interval is 0.01 day, so the effective default step is 0.01 day. Duration is 250 days. Explicit RK4 is not positivity-preserving for arbitrary rates/steps: nonfinite or negative candidate populations raise an error with advice to reduce the step; they are not clipped to zero.

The core uses EACH_WINDOW; the presenter uses ONCE_AFTER_RUN. The core records time zero and each internal step, sending the full trajectory on the final boundary. Earlier payloads are marked pending. Runtime floating-point window accumulation may add a very short final interval; those computed states are retained. More than 100000 stored points causes an explicit error rather than silent history truncation. Both model manifests pin the same runtime version.

## Contracts and failures

Six scalar-or-record inputs use the exact names and units in model.yaml and lab.yaml. Nonnegative finite numerical values are required; booleans, unknown ports and malformed values fail. All values in a multi-port update are validated before mutation. Population outputs are typed records containing role, label, count and t. The public trajectory output maps to a mixed terminal record containing status, time_unit, parameters, history, threshold times and invariant applicability. History rows contain t, prey, predator, invariant and drift. Invariant/drift are null where not applicable.

The presenter requires a complete finite, strictly increasing time axis covering the requested run and nonnegative finite populations. It reports the complete population and phase trajectories, final populations, scope and the invariant audit. The ≤1e-9 count threshold is a numerical diagnostic, not an extinction forecast. Visuals and terminal result records use the same final history.

## Acceptance and evidence

The full 250-day default is compared with independently implemented equations integrated by SciPy DOP853 at rtol=1e-12 and atol=1e-13. Every population sample must satisfy rtol=1e-5 and atol=1e-6; maximum absolute invariant drift must be below 1e-6. Convergence compares the maximum error over the whole 20-day trajectory at independently calculated RK4 steps 0.1, 0.05 and 0.025 day; each error reduction must exceed 12. A whole-trajectory norm avoids interpreting accidental endpoint cancellation as convergence order.

Additional checks cover the stationary equilibrium N=gamma/delta, P=alpha/beta, analytic exponential boundary solutions, input rejection, initial-condition overrides, reset reproducibility, manifest contracts, full graph timing and malformed presenter input. Execution results and computed errors are in verification.json; tests are executable and reference solver values are not hardcoded.

## Scientific limits

No empirical calibration, field-data comparison, uncertainty estimation, intervention optimization or conservation-management validation is supplied. Rates and initial conditions outside the example remain exploratory and may require smaller steps. A successful local run is not a managed Passport or authorization to publish. Source-faithful custom equations and changed presentation need manual scientific review before promotion.
