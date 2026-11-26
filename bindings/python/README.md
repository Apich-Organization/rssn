# rssn for Python

`ctypes` bindings for the rssn C ABI, plus one function per registered
operator. No dependencies; numpy is used for tensors when it is installed.

```python
import rssn

s = rssn.Session()
x = s.sym("x")

s.compute(rssn.diff(rssn.sin(x) * rssn.cos(x), x)).text      # 'cos(2*x)'
s.compute(rssn.diff(rssn.sin(x) * rssn.cos(x), x),
          rssn.Config().numeric(1e-12).bind("x", 0.0)).value  # 1.0
```

## Install

```sh
pip install -e bindings/python          # or: pip install bindings/python
export RSSN_LIB=/path/to/rssn.dll       # or librssn.so / librssn.dylib
```

Without `RSSN_LIB` the package looks in `rssn/_lib/`, then in a source
checkout's `target/debug`, `target/ffi` and `target/release`, then on the
platform loader path. `rssn.library_path()` reports what it bound.

## Shape of the API

| Object | Role |
|---|---|
| `rssn.Session()` | owns every term; single-threaded; `with` block or `.close()` |
| `rssn.Term` | a handle; `str()` is infix text, `.latex` / `.typst` / `.pretty`, `.value`, `.tensor()` |
| `rssn.Config()` | target, bindings, assumptions, budget, dimension mode |
| `rssn.Answer` | `.term`, `.text`, `.value`, `.error`, `.reduced`, `.iterations`, `.dimension`, `.diagnostic_codes` |
| `rssn.diff(f, x)` | a generated operator function; returns an **unevaluated** term |
| `rssn.OPS` | `name -> (min_arity, max_arity, ffi_reachable)` |
| `rssn.NAMES` | `name -> Python identifier` (differs only for keywords) |
| `rssn.BY_RULE_SET`, `rssn.BY_CATEGORY` | discovery |
| `rssn.units`, `rssn.constants` | the manifest's unit and constant tables, exact rationals |
| `rssn.Value` | a dynamic value: `kind()`, `float()`, `text()`, `json()`, `to_term()`, `tensor()`, `len()`, `[i]` |
| `rssn.sim` | `run(name, params)`, `names()`, `params(name)`, `defaults(name)` |
| `rssn.manifest()`, `rssn.manifest_hash()` | the kernel manifest |

Operators are *requests*: `rssn.diff(f, x)` builds the derivative,
`session.compute(...)` reduces it. `Session.apply("name", *args)` reaches an
operator by name; `Session.call` is the same function, matching the kernel ABI's
vocabulary.

## Notes

- `import rssn` is the recommended form. `from rssn import *` also brings in the
  operators named after builtins (`abs`, `sum`, `pow`, `round`, `list`) and
  would shadow them.
- Operator names that are Python keywords get a trailing underscore: `and_`,
  `not_`, `or_`. The original spelling stays in `OPS` and works with
  `Session.apply`.
- Arithmetic on `Term` goes through the engine's own `sub`, `neg` and `divide`
  operators (the quotient is `divide` because `div` is the divergence). The
  engine rewrites them to `add(a, mul(-1, b))` and `mul(a, pow(b, -1))` and
  prints both spellings the same way.
- Operators whose arguments are native objects — AD tapes, compiled backends,
  solver handles — are called with `Session.call_value`, which uses the value
  ABI instead of the term ABI:

      backend = s.call_value("num_backend_current")
      s.call_value("num_backend_name", backend).text()      # 'interpreter'

  `OPS[name][2]` is `True` for every operator reachable through one of the two
  ABIs — which is all 5411 of them.
- Higher-order operators (`num_bfgs`, `num_brent`, the adaptive quadrature
  family, ...) take a *function*; pass any Python callable and it becomes a host
  callback:

      s.call_value("num_brent", 1.0, 2.0, 1e-12, 64, lambda x: x * x - 2)

  See `rssn.Callback` for a registration that outlives one call, and
  `rssn.callback_signatures` for what an operator accepts at a position.
- `Answer.diagnostic_codes` holds the stable codes of the advisory dimension
  findings; the prose of each finding is in `Session.compute_json(...)`.

## Tests

```sh
export RSSN_LIB=/path/to/rssn.dll
python -m unittest discover -s bindings/python/tests -v
python bindings/tests/run_conformance.py     # the cross-language table
```
