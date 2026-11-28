# rssn for Julia

`ccall` bindings for the rssn C ABI, plus one function per registered operator.
No package dependencies: the JSON codec the simulation entry point needs is
included.

```julia
using RSSN

session = Session()
x = sym(session, "x")

compute(diff(sin(x) * cos(x), x)).text                 # "cos(2*x)"
cfg = bind!(numeric!(Config(), 1e-12), "x", 0.0)
compute(diff(sin(x) * cos(x), x), cfg).value           # 1.0
```

## Install

```sh
julia --project=bindings/julia -e 'using Pkg; Pkg.instantiate()'
export RSSN_LIB=/path/to/library-directory  # a full shared-library path also works
```

Or develop it into an environment:

```julia
# not run: the installation step
using Pkg; Pkg.develop(path = "bindings/julia"); using RSSN
```

Without `RSSN_LIB` the package looks in `lib/` and at the extracted
distribution root when running from its `bindings/julia` project. It also
searches a source checkout's Bazel FFI output directories, including
`bazel-bin/ffi_distribution/ffi`, `bazel-bin/ffi_distribution_portable/ffi` and
`bazel-bin/ffi_dev/ffi`, then its Cargo output directories and the platform
loader path. `RSSN.libpath()` reports what it bound
and `RSSN.set_library!(path)` accepts a file or directory before you create
sessions.

## Shape of the API

| Object | Role |
|---|---|
| `Session()` | owns every term; single-threaded; `close(session)` or a finalizer |
| `Term` | a handle; `text(t)`, `latex(t)`, `typst(t)`, `pretty(t)`, `value(t)`, `tensor(t)` |
| `Config()` | target, bindings, assumptions, budget, dimension mode |
| `Answer` | `term`, `text`, `value`, `error`, `reduced`, `iterations`, `dimension`, `diagnostic_codes`, plus `.latex` / `.typst` / `.pretty` / `.tensor` |
| `diff(f, x)` | a generated operator function; returns an **unevaluated** term |
| `OPS` | `name => (min_arity, max_arity, ffi_reachable)` |
| `NAMES` | `name => binding name` |
| `BY_RULE_SET` | discovery |
| `UNITS`, `CONSTANTS` | the manifest's tables, exact rationals as strings |
| `Value` | a dynamic value: `kind`, `value`, `text`, `json`, `term`, `tensor`, `length`, `[i]` |
| `call_value(session, op, args...)` | call an operator whose arguments are native objects |
| `run_scenario(name, params)` | a simulation, JSON in and a `Dict` out |
| `manifest_json_string()`, `manifest_hash()` | the kernel manifest |

Operators are *requests*: `diff(f, x)` builds the derivative, `compute(...)`
reduces it. `apply(session, "name", args...)` reaches an operator by name.

## Workflow: symbolic calculus

```julia
x = sym(session, "x")

t = diff(sin(x) * cos(x), x)      # a request, not a result
text(t)                           # "diff(cos(x)*sin(x), x)"
compute(t).text                   # "cos(2*x)"

compute(integral(x^2, x)).text                       # "1/3*x^3"
compute(parse_term(session, "defint(x^2, x, 0, 1)")).text
compute(parse_term(session, "taylor(sin(x), x, 0, 5)")).text
compute(parse_term(session, "limit(sin(x)/x, x, 0)")).text
simplify(session, "(x^2 - 1)/(x - 1)").text          # "x + 1"
```

A `Term` overloads the arithmetic operators and the `Base` functions, so a
request reads like the mathematics:

```julia
x^2 + 2*x + 1        # a request; compute() gives "x^2 + 2*x + 1"
sin(x); sqrt(x); exp(x); log(x); log(2, x); abs(x)
```

## Workflow: a numeric target

```julia
cfg = bind!(numeric!(Config(), 1e-12), "x", 0.5)
answer = compute(diff(sin(x), x), cfg)
answer.value          # 0.877582561890373
answer.error          # the estimated absolute error

eval_term(parse_term(session, "sin(x) + 1"), Dict("x" => 0.0))   # 1.0
```

`compute_all` runs several requests in parallel inside the library:

```julia
answers = compute_all(session, ["diff(x^3, x)", "integral(x^2, x)", "gcd(4, 6)"])
compute_all(session, ["1 + 1", "diff(x^2, x)"], config = cfg, errors = :collect)
```

## Workflow: linear algebra and tensors

```julia
m = parse_term(session, "matmul(list(list(1,2),list(3,4)), list(list(0,1),list(1,0)))")
tensor(compute(m).term)          # a 2x2 Matrix, row-major reading order

compute(parse_term(session, "det(list(list(2,0),list(0,3)))")).text   # "6"
compute(parse_term(session, "inverse(list(list(1,2),list(3,4)))")).text
```

## Workflow: NumPy arrays

The NumPy codec accepts Julia arrays and returns standard NumPy bytes. The
decoded tensor preserves the original shape.

```julia
if get(RSSN.json_decode(manifest_json_string())["features"], "npy", false)
    payload = encode_npy(session, [1.0 2.0; 3.0 4.0])
    restored = decode_npy(session, payload)
    println(tensor(restored))
    close(restored)
end
```

## Workflow: a simulation

```julia
SCENARIO_NAMES                          # the 29 scenarios
run_scenario("poisson_multigrid")       # the defaults
result = run_scenario("ising", Dict("width" => 8, "height" => 8, "temperature" => 2.0,
                                    "mc_steps" => 10, "seed" => 3))
result["magnetization"]
```

## Workflow: dynamic values

```julia
tape = call_value(session, "num_ad_tape_new", 16)
v = call_value(session, "num_ad_tape_var", tape, 1.5)
z = call_value(session, "num_ad_var_mul", call_value(session, "num_ad_var_sin", v), v)
value(call_value(session, "num_ad_var_value", z))      # sin(1.5) * 1.5

record = value(session, (estimate = 2.5, method = "fit"))
value(record[1])                                       # 2.5
text(record["method"])                                 # "fit"

moments = call_value(session, "num_timeseries_rolling_moments", [1.0, 2.0, 3.0], 2, 1)
tensor(moments["variance"])                             # [0.5, 0.5]
```

A record accepts exact field names and one-based positions. An unknown field
name raises `KeyError`.

`value(session, x)` accepts `BigInt`. Read an integer cell as decimal digits
with `text`, then parse those digits as a Julia `BigInt`.

```julia
exact = value(session, big(2)^100 + 7)
parse(BigInt, text(exact))
close(exact)
```

## Workflow: units, constants and dimensions

```julia
unit("ft").dimension        # "Length"
unit("metre").name          # "m"
value(constant("c"))        # 2.99792458e8

cfg = dimension!(dimension!(Config(), "x", "m"), "t", "s")
dimension(diff(x^2, x), cfg)     # "Length"

strict = dimension_mode!(dimension!(dimension!(Config(), "x", "m"), "y", "s"), "strict")
try
    compute(parse_term(session, "x + y"), strict)
catch error
    error isa RssnError || rethrow()
    error.status                 # 12, the RssnStatus dimension code
    sprint(showerror, error)     # starts with "RssnError(dimension):"
end
```

## Workflow: choosing a backend

The loaded library uses automatic selection until the process default is set explicitly.
Query the host's usable backends and change the process default when needed.

```julia
println(available_backends(session))
println(current_backend(session))
try
    @assert set_default_backend!(session, "interpreter") == "interpreter"
    @assert current_backend(session) == "interpreter"
finally
    set_default_backend!(session, "auto")
end
```

For explicit native compilation, pass the handle from
`open_backend(session, name)` to `call_value(session, "num_backend_compile", ...)`.
Release the handle with `release_backend!(handle)` and close its value when
finished.

## Workflow: passing measured floats

The binary64 values below are already in coherent units. The dimension
object follows the kernel manifest's `DimExpr` JSON format. Attaching a
dimension preserves tensor shape and does not apply a unit scale or offset.

```julia
distance = Dict("base" => "Length")
differences = with_coherent_dimension(RSSN.value(session, [2.0, 3.0]), distance)
initial = with_coherent_dimension(RSSN.value(session, [0.0]), distance)
positions = call_value(session, "num_timeseries_integrate", differences, initial)
@assert tensor(positions) == [0.0, 2.0, 5.0]
```

## Workflow: a checked dynamic operator

Pass `checked=true` on a dynamic operator to send native inputs directly to its
declared signature. The result is a session-owned `Value`.

```julia
lerch_value = num_special_ext_lerch(-1234.0, 2.0, 1.0; session, checked=true)
println(value(lerch_value))  # 0.0218615893044645
close(lerch_value)
```

## Workflow: finding an operator

```julia
print(op_help("diff"))            # the signature, arity, category and the engine's docs
op_search("brent")                # [(name, first line of the doc), ...]
op_search(category = "calculus", limit = 5)
length(op_names())                # count for this binding release
```

The documentation of every operator is generated from the kernel manifest into
`OPERATOR_DOCS`, which is also what each function's docstring interpolates, so
`?diff` and `op_help("diff")` read the same text.

## Names

`using RSSN` brings in every operator whose name `Base` does not already export.
The 46 that are mathematical `Base` functions (`sin`, `cos`, `diff`, `sum`,
`sqrt`, ...) are added as `Base` methods on `Term`, so they are called
unqualified and nothing is ambiguous:

```julia
sin(x)            # Base.sin(::Term)
diff(x^2, x)      # Base.diff(::Term, ::Any)
sum(x, 1)         # Base.sum(::Term, ...)
```

The rest of the `Base` collisions keep their own meaning, so the operator is
mangled with a trailing underscore: `RSSN.collect_`, `RSSN.convert_`,
`RSSN.pi_`, `RSSN.zeros_`. `NAMES["collect"] == "collect_"`. The same suffix
applies to names the hand-written API uses (`apply_`, `dimension_`, `sym_`) and
to Julia keywords.

Arithmetic on `Term` goes through the engine's own `sub`, `neg` and `divide`
operators (the quotient is `divide` because `div` is the divergence).

## Notes

- Operators whose arguments are native objects are called through the value ABI:

      current_backend(session)                                # effective backend

  `OPS[name][3]` is `true` for every operator: the value ABI carries native
  objects and host callbacks alike.
- A higher-order operator takes a function. Any Julia function works, and the
  `do` form puts it last:

      call_value(session, "num_brent", 1.0, 2.0, 1e-12, 64) do x
          x * x - 2
      end

  `Callback` registers one explicitly; `callback_signatures` says what an
  operator accepts at a position.
- `tensor(term)` returns an `Array` in row-major reading order (the C side is
  row-major; the reshape permutes the axes back).
- `default_session()` is per task, and `with_session(session) do ... end` scopes
  it, which is what the operator functions use when no argument carries a
  session.

## Tests

```sh
export RSSN_LIB=/path/to/rssn.dll
julia --project=bindings/julia bindings/julia/test/runtests.jl
python bindings/tests/run_readme_examples.py julia   # every example in this file
```

`test/runtests.jl` runs the cross-language table in `bindings/tests/cases.tsv`
plus Julia-specific tests (the JSON codec, `Base` extensions, the row-major
tensor, unit and constant tables).
