# Ecology: Predator-Prey

## Scientific Question
How do prey and predator populations co-evolve in the canonical Lotka-Volterra system?

## Biological Context
This space runs a single deterministic Lotka-Volterra model with built-in population trajectories, phase portrait, summary metrics, and invariant-drift diagnostics.

## Mechanistic Assumptions
- The dynamics are exactly:
  - `dX/dt = alpha * X - beta * X * Y`
  - `dY/dt = delta * X * Y - gamma * Y`
- State variables are continuous real-valued populations.
- RK4 is used with a fixed simulation step.

## Model Rationale
- The implementation is intentionally strict and omits environmental forcing, carrying capacity, satiation, and stochasticity.
- Visualizations are owned directly by the Lotka-Volterra model instead of separate observer modules.

## Expected Behaviors
- Oscillatory prey and predator trajectories.
- Closed phase-space orbits around the non-trivial equilibrium.
- Small conserved-quantity drift under RK4 integration.

## Known Limitations
- Two-species canonical system only.
- No environmental forcing or density dependence.
- No parameter sweep or batch-analysis visuals are included in this space.

## How to Run
```bash
python spaces/ecology-predator-prey/run_local.py --duration auto --tick-dt auto
python spaces/ecology-predator-prey/simui_local.py --port 8765
```

## How to Interpret Outputs
- Use the population trajectories panel to compare prey and predator timing.
- Use the phase portrait to inspect orbit geometry, nullclines, and the equilibrium point.
- Use the summary table and invariant audit to verify that the run behaves like a faithful numerical Lotka-Volterra simulation.
