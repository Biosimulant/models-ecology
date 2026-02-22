# Ecology: Predator-Prey

## Scientific Question
How do rabbit and fox populations co-evolve under fixed abiotic conditions?

## Biological Context
This space composes a classic two-species predator-prey system with a shared environment module, explicit predation coupling, and monitor modules for timeseries, phase-space, and summary metrics.

## Mechanistic Assumptions
- Environment conditions are broadcast to both populations at each simulation tick.
- Predation is captured by a Lotka-Volterra-style interaction module.
- Predation removes prey and returns food gain to predators.
- Monitoring modules are passive observers and do not affect dynamics.

## Wiring Rationale
- `environment.conditions` drives both `rabbits` and `foxes`.
- Population states feed into `predation` for causal interaction.
- Both population streams are fanned out to monitor/metrics modules.
- Predation outputs are routed back to species-specific input ports.

## Expected Behaviors
- Oscillation or damped cycles in prey/predator counts.
- Structured trajectories in rabbit-vs-fox phase-space.
- Stability/extinction metrics available in the metrics table.

## Known Limitations
- Two-species system only; no additional trophic levels.
- No stochastic external forcing besides module-internal randomness.
- Fixed baseline environment in default configuration.

## How to Run
```bash
python spaces/ecology-predator-prey/run_local.py --duration auto --tick-dt auto
python spaces/ecology-predator-prey/simui_local.py --port 8765
```

## How to Interpret Outputs
- Use `PopulationMonitor` for direct trajectory comparison.
- Use `PhaseSpaceMonitor` to inspect cycle geometry and attractor behavior.
- Use `EcologyMetrics` for aggregate indicators (extinctions, diversity, stability).
