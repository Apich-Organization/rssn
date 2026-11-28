"""
    RSSN

Julia bindings for [rssn](https://github.com/Apich-Organization/rssn), a
scientific computing library: a computer algebra system, numerical kernels and
physics simulations behind one small API.

The package binds the C ABI in `rssn.h`, so it needs a shared library: set
`RSSN_LIB` to `rssn.dll` / `librssn.so` / `librssn.dylib`, or build one in a
source checkout with `bash scripts/build_ffi.sh`.

```julia
using RSSN

session = Session()
x = sym(session, "x")

answer = compute(diff(sin(x) * cos(x), x))
answer.text                       # "cos(2*x)"

cfg = numeric!(Config(), 1e-12)
bind!(cfg, "x", 0.0)
compute(diff(sin(x) * cos(x), x), cfg).value    # 1.0
```

Operators are requests, not results: `diff(f, x)` returns an unevaluated term
and `compute` reduces it. Every registered operator has a function here, named
as the manifest names it; `OPS`, `NAMES` and `BY_RULE_SET` in `ops.jl` hold the
table, and `apply(session, "name", args...)` reaches one by name at runtime.

Names that `Base` already exports (`sin`, `cos`, `diff`, `sum`, ...) are not
re-exported — that would make `using RSSN` ambiguous — but they are defined on
`Term`, so `sin(x)` and `diff(f, x)` work as written. The remaining 5300-odd
operators are exported by `exports.jl`.
"""
module RSSN

include("ffi.jl")
include("json.jl")
include("config.jl")
include("session.jl")
include("value.jl")
include("callback.jl")
include("backend.jl")
include("npy.jl")
include("sim.jl")
include("ops.jl")
include("base_ext.jl")
include("exports.jl")
include("units.jl")
include("constants.jl")
include("scenarios.jl")
include("docs.jl")

export
    # lifetime and configuration
    Session, Config,
    symbolic!, numeric!, bind!, assume!, budget!, dimension_mode!,
    bypass_dimensions!, dimension!, define_base_dimension!,
    Assumption, DimensionMode,
    # terms
    Term, Coercible, NAMES, sym, integer, number, rational, parse_term, term, apply, text, latex,
    typst, pretty, value, tensor, to_float, handle,
    # computation
    Answer, compute, compute_all, simplify, compute_json, compute_bincode,
    compute_all_bincode, eval_term, dimension, last_error,
    # dynamic values (operators with native-object arguments)
    Value, call_value, call_key, kind, json, value_count, VALUE_KINDS,
    with_coherent_dimension, tensor_cells,
    available_backends, current_backend, open_backend, release_backend!, set_default_backend!,
    encode_npy, decode_npy,
    # host callbacks (operators that take a function)
    Callback, callback_signatures, parse_signature,
    # sessions and metadata
    default_session, with_session, libpath, set_library!, version,
    manifest_json_string, manifest_hash,
    # simulations
    run_scenario, SCENARIO_NAMES, SCENARIO_PARAMS,
    # the operator reference
    op_info, op_help, op_search, op_names, op_categories, op_rule_sets, BY_CATEGORY,
    OPERATOR_DOCS,
    # units and constants
    UnitDef, UNITS, UNIT_BY_ALIAS, unit, unit_names,
    PhysicalConstant, CONSTANTS, CONSTANT_BY_SYMBOL, constant,
    # errors
    RssnError

end # module RSSN
