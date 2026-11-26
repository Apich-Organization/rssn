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
export RSSN_LIB=/path/to/rssn.dll       # or librssn.so / librssn.dylib
```

Or develop it into an environment:

```julia
using Pkg; Pkg.develop(path = "bindings/julia"); using RSSN
```

Without `RSSN_LIB` the package looks in `src/../lib/`, then in a source
checkout's `target/debug`, `target/ffi` and `target/release`, then on the
platform loader path. `RSSN.libpath()` reports what it bound and
`RSSN.set_library!(path)` rebinds.

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

## Names

`using RSSN` brings in every operator whose name `Base` does not already export.
The 46 that are mathematical `Base` functions (`sin`, `cos`, `diff`, `sum`,
`sqrt`, ...) are added as `Base` methods on `Term`, so they are called
unqualified and nothing is ambiguous:

```julia
sin(x)            # Base.sin(::Term)
diff(f, x)        # Base.diff(::Term, ::Any)
sum(...)          # Base.sum(::Term, ...)
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

      backend = call_value(session, "num_backend_current")
      text(call_value(session, "num_backend_name", backend))   # "interpreter"

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
```

`test/runtests.jl` runs the cross-language table in `bindings/tests/cases.tsv`
plus Julia-specific tests (the JSON codec, `Base` extensions, the row-major
tensor, unit and constant tables).
