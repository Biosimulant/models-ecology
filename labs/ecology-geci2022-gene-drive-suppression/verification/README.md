# Native reference fixtures

The gzip TSV files are exact compressed exports from the original Julia scientific functions. `fixture-sha256.json` pins their bytes; `native-provenance.toml` pins the source, source ranges and runtime. The four scenarios have explicit matched initial vectors and parameters in the scenario files. Tests compare every genotype at every recorded generation and every process-matrix element, rather than relying on summary population plots.

To regenerate, use Julia 1.10.12 and a Python 3.12 environment with `sympy==1.14.0`. Point PyCall at that Python interpreter via its `PYTHON` environment variable, instantiate the pinned Julia environment and build PyCall:

```
julia --project=verification/julia -e 'using Pkg; Pkg.instantiate(); Pkg.build("PyCall")'
julia --project=verification/julia verification/generate_reference.jl
```

PyCall requires a discoverable shared Python library. Run from the Lab directory. The harness writes raw TSV files under `verification/geci-native-reference`; it does not replace the checked-in fixtures automatically. Compare regenerated files and source checksums before updating the fixtures and their hashes. The core Python tests use compressed native fixtures without requiring Julia.

`upstream_parameters` uses the original standard parameter values but retains the explicitly matched wrapper initial vector. It is not an asserted published-figure reproduction. The no-release scenario captures the existing wrapper's growth from initial total I toward equilibrium total 2I; that normalization is preserved and documented, not silently changed.
# Portable reference retrieval

Source checkouts include the compressed fixtures. Environments that transfer
only text files may omit those binaries: the native-reference tests then fetch
each missing file from exact repository commit
`ad82fc61a6a67843ebca8911a4bf6ffb8627e407`, bounded to 2 MiB per file, and require
the original SHA-256 in `fixture-sha256.json` before loading it. Network access
is needed only for missing fixtures. Existing corrupted files fail the check;
they are never silently replaced. Every genotype and matrix comparison and its
tolerance are unchanged. The tests never regenerate reference expectations from
the implementation under test.
