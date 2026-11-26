# Simulation scenarios by name, with JSON parameters.
#
# `rssn_sim_run` takes a scenario name and a JSON parameter object; `{}` selects
# the defaults of the scenarios that have them. `scenarios.jl` lists every name
# and the parameters its request type reads.

"""
    run_scenario(name, parameters = Dict{String, Any}()) -> Dict

Run the simulation scenario `name` and return its JSON result as a `Dict`.

```julia
run_scenario("ising", Dict("width" => 8, "height" => 8, "temperature" => 2.0,
                           "mc_steps" => 10, "seed" => 3))
```

`{}` selects the defaults where a scenario has them. See `SCENARIO_NAMES` and
`SCENARIO_PARAMS` for the catalog.
"""
function run_scenario(name::AbstractString, parameters = Dict{String, Any}())
    payload = json_encode(parameters)
    out = Ref{Ptr{UInt8}}(C_NULL)
    _check(sim_run(name, payload, out))
    return json_decode(take_string(out[]))
end
