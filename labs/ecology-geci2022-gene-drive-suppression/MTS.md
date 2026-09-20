# Modeling technical specification

## Scientific scope and source

Single-population deterministic genotype tracking from the exact bundled BioModels MODEL2301120001 Julia source. Source SHA-256 and executed source ranges are recorded in verification/native-provenance.toml. The Python process order is mutation → homing → editing → recombination → sperm/egg production → zygote formation → density-dependent survival → selection. The 1,071 genotype and 66 gamete enumerations are compared elementwise with native Julia. This repair does not change any matrix builder, initial vector or generation equation.

## State, parameters and units

Genotype abundances are normalized continuous population quantities (dimensionless signal unit 1), not individual-animal counts. One model time unit is one discrete generation. All 24 controls and their preserved defaults are exposed in lab.yaml and the core model.yaml.

For initial wild-type total I, initialize I/2 wild-type males and I/2 wild-type females, plus I × release_size males of genotype (AB, ef, CD, CD). Eggs per female f=2Rm/theta; density coefficient alpha=I f/(Rm−1). If Z is total zygote abundance, density survival is theta alpha/(alpha+Z), followed by genotype-specific selection. This preserved wrapper normalization means a no-release wild-type equilibrium of 2I, not I. The source convenience initializer instead assigns I to each sex. Native comparisons pass the exact wrapper initial vector to the original timecourse and retain initial_size=I. Parameters named upstream_parameters match the native baseline parameter dictionary, while the explicit initial vector still matches the wrapper.

Domains: Rm>1; theta in (0,1]; I>0; release_size≥0 (a ratio, allowing more released individuals than the initial population); recombination rate in [0,0.5]; other efficiencies, costs, resistance rates, dominance factors and cofactor in [0,1]. All must be finite; booleans are rejected. Whole input updates are validated before mutation. Scenario parameters cannot change after integration starts. Identical repeated inputs do not rebuild matrices.

## Runtime and outputs

Python 3.12; core numpy==1.26.4 and biosimulant==0.0.34, presenter biosimulant==0.0.34. The core uses the supported temporal advance_window hook with EACH_WINDOW. Windows must be contiguous, ordered and bounded by nonnegative whole generations. Fractional windows fail instead of rounding or advancing beyond the requested time. Every generation is calculated, independent of the communication-window size. The presenter uses execute with ONCE_AFTER_RUN, checks full generation coverage and uses the actual final state.

The core emits population_state, gene_drive_metrics and visualisation_payload; the Lab exposes the payload as trajectory. Population records contain total_adults, adult_females and adult_males. Metric definitions: drive_frequency is all non-wild-type Y chromosomes among males; resistance_frequency is autosomal r3 alleles divided by all autosomal alleles; male_fraction is males/total (legacy fallback 0.5 for an empty population); suppression_ratio is 1−N(t)/N(0), with N(0) including the release. The latter is a signed change from initial abundance, not a causal treatment effect versus a no-release counterfactual. Zero-denominator drive/resistance frequencies retain the source wrapper's zero fallback. All definitions and normalization appear in the terminal trajectory and summary.

Runtime rejects nonfinite or negative genotype abundances. No state is clipped or rounded to an extinction threshold. Extremely small positive deterministic populations remain mathematical abundances, not evidence that individual organisms persist. Matrix and genotype arithmetic remains unchanged from the verified wrapper.

## Native reference and acceptance

The headless Julia harness includes the original scientific source ranges unchanged, omitting plotting imports and unrelated plotting/spatial workflows. It runs Julia 1.10.12, SymPy.jl from its pinned Project/Manifest and Python SymPy 1.14.0, with the source's 120-bit BigFloat setting. Native floating-point conversions follow the original source; high-precision reference is not claimed independently of those conversions.

Four scenarios: preserved published defaults; upstream standard parameters with matched wrapper initial state; nonzero costs, mutation and resistance; no release. Compare all 1,071 genotype values at generations 0..100 with rtol=1e-9, atol=1e-11; compare all seven process matrices with rtol=1e-10, atol=1e-12. Native sparse matrix exports retain original row/column order; the selection matrix is an N×1 column and is compared with the Python vector. Exact genotype/gamete ordering and compressed fixture checksums are checked separately. Maximum observed trajectory error is about 1.23e-15. The nonzero-cost scenario prevents default zero costs from masking selection defects.

Native fixtures were produced by verification/generate_reference.jl and stored as checksummed compressed TSV, not fitted to the Python outputs. Re-running that harness writes geci-native-reference beside it; tests consume the pinned gzip fixtures. Four initial vectors and both native-parameter and Python-keyword mappings are included. The test suite also verifies probability domains, invalid-window atomicity, parameter overrides, communication-window independence, reset reproducibility and final visual coverage.

## Limits and review

Native code parity establishes implementation consistency for the tested scenarios, not empirical validation. This Lab is not calibrated to an observed species population and does not implement the Julia file's spatial, migration or multiple-release studies. The default wrapper parameters differ from the original standard parameter set and are labeled as such. No stochastic extinction, field efficacy or biosafety assurance is implied. Manual scientific review and managed validation remain necessary before public promotion.
