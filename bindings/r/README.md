# rssn for R

An R package over the rssn C ABI: a computer algebra system, numerical kernels
and physics simulations behind one small API.  No package dependencies -- the
JSON codec the value and simulation entry points need is included, and the
shared library is loaded at run time.

```r
library(rssn)

s <- rssn_session()
x <- sym("x", s)

compute(rssn_diff(sin(x) * cos(x), x))$text                    # "cos(2*x)"
compute(rssn_diff(sin(x) * cos(x), x),
        rssn_config() |> cfg_numeric(1e-12) |> cfg_bind("x", 0))$value   # 1
```

## Install

```sh
R CMD INSTALL bindings/r
export RSSN_LIB=/path/to/extracted-rssn  # a full shared-library path also works
```

The search order starts with `options(rssn.library =)` and `RSSN_LIB`, then
checks the installed package and the root of an extracted SDK. From a source
checkout the search prefers the Bazel FFI output directories, including
`bazel-bin/ffi_distribution/ffi`, `bazel-bin/ffi_distribution_portable/ffi` and
`bazel-bin/ffi_dev/ffi`. It then checks Cargo output directories and the
platform loader path.
`rssn_library()` reports the bound library. You can also call
`rssn_load(path)` with a file or directory before creating a session.

The install takes a minute or two: the package carries a generated function for each operator
and the operator reference, which is the price of not spelling operator names
as strings.

## Shape of the API

| Object | Role |
|---|---|
| `rssn_session()` | owns every term; single-threaded; `close()` or the garbage collector |
| `sym()`, `parse_term()`, `term()`, `rational()` | build a term: a symbol, infix text, a coercion, an exact fraction |
| `rssn_term` | a handle; `term_text()`, `term_latex()`, `term_typst()`, `term_pretty()`, `term_value()`, `to_float()` |
| `rssn_config()` + `cfg_*()` | target, bindings, assumptions, budget, dimension mode |
| `rssn_answer` | `$term`, `$text`, `$value`, `$error`, `$reduced`, `$iterations`, `$dimension`, `$diagnostic_codes`; `tensor(answer)` |
| `rssn_diff(f, x)` | a generated operator function; returns an **unevaluated** request |
| `OPS`, `OPERATORS` | the arity table and the searchable operator reference |
| `op_help()`, `op_table()`, `op_names()`, `op_apply()` | documentation and lookup |
| `call_value()`, `rssn_value` | the value ABI: operators whose arguments are native objects |
| `rssn_callback()` | a host callback the engine calls |
| `sim_run()`, `SCENARIOS` | the simulation catalog |
| `UNITS`, `CONSTANTS` | the manifest's tables, exact rationals as text |
| `manifest()`, `manifest_json()`, `manifest_hash()` | the kernel manifest |

Operators are *requests*: `rssn_diff(f, x)` builds the derivative,
`compute()` reduces it.  `op_apply("name", ...)` reaches an operator by name at
run time.

## Workflow: symbolic calculus

```r
x <- sym("x")

t <- rssn_diff(sin(x) * cos(x), x)   # a request
term_text(t)                         # "diff(cos(x)*sin(x), x)"
compute(t)$text                      # "cos(2*x)"

compute(rssn_integral(x^2, x))$text                 # "1/3*x^3"
compute(parse_term("defint(x^2, x, 0, 1)"))$text    # "1/3"
compute(parse_term("taylor(sin(x), x, 0, 5)"))$text # "x - 1/6*x^3 + 1/120*x^5"
compute(parse_term("limit(sin(x)/x, x, 0)"))$text   # "1"
compute(parse_term("sum(1/(k*(k+1)), k, 1, n)"))$text
simplify("(x^2 - 1)/(x - 1)")$text                  # "x + 1"
```

Arithmetic and the `Math` group work directly on a term, through the engine's
own operators:

```r
x^2 + 2*x + 1        # a request; compute() gives "x^2 + 2*x + 1"
sin(x); sqrt(x); exp(x); log(x); log(x, 2); abs(x); floor(x); ceiling(x); gamma(x)
x == 0               # the engine's eq, as a term
```

## Workflow: a numeric target

```r
cfg <- rssn_config() |> cfg_numeric(1e-12) |> cfg_bind("x", 0.5)
compute(rssn_diff(sin(x), x), cfg)$value      # 0.877582561890373
compute(rssn_diff(sin(x), x), cfg)$error      # the estimated absolute error

# Without running any rule at all:
eval_term(parse_term("sin(x) + 1"), x = 0)    # 1
```

## Workflow: linear algebra and tensors

```r
m <- parse_term("matmul(list(list(1,2),list(3,4)), list(list(0,1),list(1,0)))")
tensor(compute(m))                 # a 2x2 matrix in reading order
value_tensor(value_new(compute(m)$term))   # the same, through the value ABI

compute(parse_term("det(list(list(2,0),list(0,3)))"))$text   # "6"
compute(parse_term("inverse(list(list(1,2),list(3,4)))"))$text
compute(parse_term("eigenvals(list(list(2,0),list(0,3)))"))$text
```

The engine stores tensors row-major and R stores arrays column-major, so the
extents are permuted on the way in and out: what you read is the matrix the
expression describes.

## Workflow: NumPy arrays

Use `rssn_features()` to query the loaded library's capabilities. The NumPy
codec takes a numeric array or a flat row-major vector and returns raw bytes.

```r
if (isTRUE(rssn_features()$npy)) {
    samples <- matrix(c(1.0, 2.0, 3.0, 4.0), nrow = 2L, byrow = TRUE)
    payload <- rssn_encode_npy(samples, session = s)
    restored <- rssn_decode_npy(payload, session = s)
    print(value_tensor(restored))
    close(restored)
}
```

## Workflow: a simulation

```r
sim_names()                                  # the 29 scenarios
sim_defaults("poisson_multigrid")            # a starting point for the parameters

result <- sim_run("ising", width = 8, height = 8, temperature = 2.0,
                  mc_steps = 10, seed = 3)
result$magnetization
image(matrix(unlist(result$spins), nrow = 8))
```

## Workflow: an R function as a host callback

189 operators take a function and call it once per quadrature node,
optimisation step or root-finder probe.  Any R function works:

```r
# the root of x^2 - 2
answer <- call_value("num_brent", 1, 2, 1e-12, 64, function(t) t * t - 2)
value_float(answer[[1]])                      # 1.4142135623730951

# an objective and its gradient, at argument positions 5 and 6
point <- call_value("num_bfgs", c(3, -2), 1e-10, 100, 1000,
                    function(p) sum((p - 1)^2),
                    function(p) 2 * (p - 1))
value_float(point[[1]][[1]])                  # 1

# an explicit registration, for a function that goes to several calls
f <- rssn_callback(function(t) t * t - 2, signature = "(float) -> float")
call_value("num_brent", 1, 2, 1e-12, 64, f)
close(f)
```

The function receives native R objects -- numbers, integers, logicals, complex
numbers, strings, arrays for tensors and named lists for records -- and its
return value is converted back.  An error inside it aborts the operator and its
message becomes the `rssn_error` message.

## Workflow: dynamic values

Operators whose arguments or results are native objects -- an AD tape, a
compiled backend, a solver handle -- go through the value ABI:

```r
tape <- call_value("num_ad_tape_new", 16)
v <- call_value("num_ad_tape_var", tape, 1.5)
z <- call_value("num_ad_var_mul", call_value("num_ad_var_sin", v), v)
value_float(call_value("num_ad_var_value", z))         # sin(1.5) * 1.5

grad <- call_value("num_ad_var_backward", z)
value_float(call_value("num_ad_var_grad_in", v, grad)) # cos(1.5)*1.5 + sin(1.5)

moments <- call_value("num_timeseries_rolling_moments", c(1, 2, 3), 2L, 1L)
value_tensor(moments[["variance"]])                    # 0.5, 0.5
```

`[[` reads a record field by its exact name or one-based position. An unknown
field name raises an R error.

Use `array(x, dim = 1L)` when an operator needs a one-element vector. R treats
`x` and `c(x)` as the same scalar when `x` is a single number. The dimension
makes `value_new()` and `call_value()` preserve the vector's shape.

```r
single <- array(3, dim = 1L)
value_float(call_value("num_dot", single, single))   # 9
```

Pass signed decimal text to `value_integer_decimal()` when you need an exact
integer beyond R's numeric precision. Use `value_text()` to read it back.

```r
exact <- value_integer_decimal("1267650600228229401496703205383", session = s)
value_text(exact)
close(exact)
```

## Workflow: units, constants and dimensions

```r
unit_get("ft")$dimension          # "Length"
unit_scale(unit_get("ft"))        # 0.3048
constant_value(constant_get("c")) # 299792458

cfg <- rssn_config() |> cfg_dimension("x", "m") |> cfg_dimension("t", "s")
dimension(rssn_diff(x^2, x), cfg) # "Length"

strict <- rssn_config() |>
    cfg_dimension("x", "m") |>
    cfg_dimension("y", "s") |>
    cfg_dimension_mode("strict")
try(compute(parse_term("x + y"), strict))   # an rssn_error with status "dimension"
```

## Workflow: choosing a backend

The loaded library uses automatic selection until the process default is set explicitly.
Query the host's usable backends and change the process default when needed.

```r
rssn_available_backends(s)
rssn_current_backend(s)
tryCatch({
    stopifnot(rssn_set_default_backend("interpreter", s) == "interpreter")
    stopifnot(rssn_current_backend(s) == "interpreter")
}, finally = rssn_set_default_backend("auto", s))
```

For explicit native compilation, pass a handle from
`rssn_open_backend(name, s)` to `call_value("num_backend_compile", ...)`.
Release the handle with `rssn_release_backend(handle)` before closing its value.

## Workflow: passing measured floats

The binary64 values below are already in coherent units. The dimension
object follows the kernel manifest's `DimExpr` JSON format. Attaching a
dimension preserves tensor shape and does not apply a unit scale or offset.

```r
distance <- list(base = "Length")
differences <- value_with_coherent_dimension(value_new(array(c(2.0, 3.0), 2L), s), distance)
initial <- value_with_coherent_dimension(value_new(array(0.0, 1L), s), distance)
positions <- call_value("num_timeseries_integrate", differences, initial, session = s)
stopifnot(identical(as.numeric(value_tensor(positions)), c(0.0, 2.0, 5.0)))
```

## Workflow: a checked dynamic operator

Pass `checked = TRUE` on a dynamic operator to send native inputs directly to
its declared signature. The result is a session-owned `rssn_value`.

```r
lerch_value <- rssn_num_special_ext_lerch(-1234.0, 2.0, 1.0, session = s, checked = TRUE)
value_float(lerch_value)  # 0.0218615893044645
close(lerch_value)
```

## Workflow: finding an operator

```r
op_help("diff")                  # the signature, arity, category and documentation
op_table(pattern = "brent")      # a search over names and documentation
op_table(category = "calculus")
length(op_names())               # count for this binding release
```

## Notes

- **Names.** The generated facade is prefixed: `rssn_diff()` rather than `diff`,
  so `library(rssn)` masks nothing in base R.  Arithmetic and the `Math` group
  are the exception, by design: `x^2`, `sin(x)` and `sqrt(x)` read naturally
  because they are methods on `rssn_term`, not new functions.
- **Sessions.** A session is single-threaded (its store sits behind a `RefCell`
  in Rust).  Use one per thread, or `compute_all()` for parallel requests inside
  one session.  The operator functions take the session from their first term
  argument and fall back to `default_session()`; `with_session(s, ...)` scopes
  that fallback.
- **Engine errors.** A failed engine call raises an `rssn_error` with
  `$status` (the `RssnStatus` code) and `$status_name`.
  `tryCatch(parse_term("sin("), rssn_error = function(e) e$status_name)` reads
  the status name.
- **Integers.** R has no integer literal: `1` is a double and `1L` an integer,
  and the value ABI is strict about the difference.  `call_value()` flips an
  argument's numeric kind when the engine says it wanted the other one, so `64`
  reaches an `Int` parameter without `64L`.
- **Partial matching.** `sim_run()` names its first argument `scenario`, not
  `name`, so that `sim_run("ising", n = 8)` cannot be read as the scenario name
  by R's partial matching.
- **No `text()`/`pretty()`.** They would mask `graphics::text` and
  `base::pretty`; the renderings are `term_text()`, `term_latex()`,
  `term_typst()` and `term_pretty()`, and `as.character()` is the text.
- **The JSON codec** is `json_decode()` / `json_encode()`: an object becomes a
  named list, an array an unnamed list, and `null` `NULL`.

## Tests

```sh
export RSSN_LIB=/path/to/rssn.dll
Rscript bindings/tests/run_conformance.R       # the cross-language table
python bindings/tests/run_readme_examples.py r   # every example in this file
Rscript bindings/r/tests/test-rssn.R           # the R-specific half
```

`run_conformance.R` checks the same `bindings/tests/cases.tsv` as the Python,
Julia and C++ runners.  `test-rssn.R` covers what the table cannot see: the R
coercions, the `Ops`/`Math` groups, the JSON codec, the callback trampoline and
the value objects.
