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

## Workflow: symbolic calculus

```python
x = s.sym("x")

t = rssn.diff(rssn.sin(x) * rssn.cos(x), x)   # a request, not a result
str(t)                                        # 'diff(cos(x)*sin(x), x)'
s.compute(t).text                             # 'cos(2*x)'

s.compute(rssn.integral(x**2, x)).text                        # '1/3*x^3'
s.compute("defint(x^2, x, 0, 1)").text                        # '1/3'
s.compute("taylor(sin(x), x, 0, 5)").text                     # 'x - 1/6*x^3 + 1/120*x^5'
s.compute("limit(sin(x)/x, x, 0)").text                       # '1'
s.compute("sum(1/(k*(k+1)), k, 1, n)").text                   # '1 - 1/(n + 1)'
s.simplify("(x^2 - 1)/(x - 1)").text                          # 'x + 1'
```

A `Term` overloads the arithmetic operators, so a request reads like the
mathematics:

```python
x**2 + 2*x + 1        # a request; compute() gives 'x^2 + 2*x + 1'
abs(x); rssn.sin(x); rssn.sqrt(x); rssn.exp(x); rssn.ln(x); rssn.log(2, x)
```

## Workflow: a numeric target

```python
cfg = rssn.Config().numeric(1e-12).bind("x", 0.5)
answer = s.compute(rssn.diff(rssn.sin(x), x), cfg)
answer.value          # 0.877582561890373
answer.error          # the estimated absolute error

s.eval("sin(x) + 1", {"x": 0.0})     # evaluation without running any rule
```

`compute_all` runs several requests in parallel inside the library:

```python
answers = s.compute_all(["diff(x^3, x)", "integral(x^2, x)", "gcd(4, 6)"])
s.compute_all(["1 + 1", "diff(x^2, x)"], cfg, errors="collect")   # per-request failures
```

## Workflow: linear algebra and tensors

```python
m = s.parse("matmul(list(list(1,2),list(3,4)), list(list(0,1),list(1,0)))")
s.compute(m).tensor()                  # a numpy array, row-major reading order

s.compute("det(list(list(2,0),list(0,3)))").text        # '6'
s.compute("inverse(list(list(1,2),list(3,4)))").text
s.compute("eigenvals(list(list(2,0),list(0,3)))").text
```

## Workflow: a simulation

```python
rssn.sim.names()                                   # the 29 scenarios
rssn.sim.defaults("poisson_multigrid")             # a starting point for the parameters

result = rssn.sim.run("ising", {"width": 8, "height": 8, "temperature": 2.0,
                                "mc_steps": 10, "seed": 3})
result["magnetization"]
```

## Workflow: dynamic values

Operators whose arguments or results are native objects go through the value
ABI:

```python
backend = s.call_value("num_backend_current")
backend.kind()                                     # 'opaque'
s.call_value("num_backend_name", backend).text()   # 'interpreter'

tape = s.call_value("num_ad_tape_new", 16)
v = s.call_value("num_ad_tape_var", tape, 1.5)
z = s.call_value("num_ad_var_mul", s.call_value("num_ad_var_sin", v), v)
s.call_value("num_ad_var_value", z).float()        # sin(1.5) * 1.5

grad = s.call_value("num_ad_var_backward", z)
s.call_value("num_ad_var_grad_in", v, grad).float()   # cos(1.5)*1.5 + sin(1.5)
```

## Workflow: units, constants and dimensions

```python
rssn.units.get("ft").scale             # Fraction(381, 1250) -- exact
rssn.units.get("m").dimension          # 'Length'
rssn.constants.get("c").value()        # 299792458.0

cfg = rssn.Config().dimension("x", "m").dimension("t", "s")
s.dimension(rssn.diff(x**2, x), cfg)   # 'Length'

strict = (rssn.Config().dimension("x", "m").dimension("y", "s")
          .dimension_mode("strict"))
try:
    s.compute("x + y", strict)
except rssn.RssnError as error:
    error.status.name                  # 'dimension'
    error                              # dimension mismatch in `x + y`: Length vs Time
```

## Workflow: finding an operator

```python
print(rssn.op_help("diff"))            # the signature, arity, category and docs
rssn.op_search("brent")                # [(name, first line of the doc), ...]
rssn.op_search(category="calculus", limit=5)
len(rssn.op_names())                   # 5411
```

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
python bindings/tests/run_readme_examples.py python   # every example in this file
```
