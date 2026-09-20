using SymPy
using LinearAlgebra
using DelimitedFiles
using SHA
using TOML
const ROOT = @__DIR__
const SOURCE = joinpath(ROOT, "../models/core/upstream/MODEL2301120001.jl")
const SOURCE_BYTES = read(SOURCE)
const SOURCE_LINES = split(String(copy(SOURCE_BYTES)), '\n')
const OUT = joinpath(ROOT, "geci-native-reference")
mkpath(OUT)
given_precision = 120
setprecision(given_precision)
function load_exact(first, last)
    println("Loading original source lines ", first, ":", last); flush(stdout)
    include_string(Main, join(SOURCE_LINES[first:last], "\n"), SOURCE * ":" * string(first))
end
load_exact(45, 1176)
load_exact(2549, 2688)
load_exact(3679, 3730)
load_exact(4580, 4581)
for line in 4583:4588
    load_exact(line, line)
end
for line in 4591:4596
    load_exact(line, line)
end
for line in 4598:4603
    load_exact(line, line)
end
load_exact(4607,4664)
load_exact(5274,5335)
writedlm(joinpath(OUT,"genotypes.tsv"),reduce(vcat,genotypes_detailed),'\t')
writedlm(joinpath(OUT,"gametes.tsv"),reduce(vcat,gametes_detailed),'\t')
println("Native scientific functions and symbolic matrices ready");flush(stdout)
# The default Python scenario is supplied as an explicit matched-parameter fixture.
fixtures=TOML.parsefile(joinpath(ROOT,"geci-reference-scenarios.toml"))
for case in fixtures["cases"]
    name=case["name"]
    println("Evaluating ",name);flush(stdout)
    parameters=Dict(k=>BigFloat(v) for (k,v) in case["parameters"])
    apply_parameters_set(parameters)
    # Explicit initial state, including its normalization, is supplied to the
    # original native timecourse. Do not silently change the published wrapper.
    start=BigFloat.(readdlm(joinpath(ROOT,case["initial_vector"]),Float64)[:,1])
    native=timecourse(case["generations"],start;initial_size=case["initial_size"])
    writedlm(joinpath(OUT,name*"-genotypes.tsv"),Float64.(native["genotypes"]),'\t')
    for (key,matrix) in current_matrices
        open(joinpath(OUT,name*"-"*key*".tsv"),"w") do io
            for index in CartesianIndices(matrix)
                value=Float64(matrix[index])
                if value != 0
                    indices=Tuple(index)
                    println(io,join(indices,'\t'),'\t',repr(value))
                end
            end
        end
    end
    println("Completed ",name);flush(stdout)
end
open(joinpath(OUT,"provenance.toml"),"w") do io
    TOML.print(io,Dict("source_sha256"=>bytes2hex(sha256(SOURCE_BYTES)),"julia_version"=>string(VERSION),"sympy_version"=>string(SymPy.sympy.__version__),"precision_bits"=>given_precision,"source_ranges"=>["45:1176","2549:2688","3679:3730","4580:4603","4607:4664","5274:5335"],"science_source_modified"=>false))
end
println("All reference exports complete")
