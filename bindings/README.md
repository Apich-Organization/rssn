# rssn language bindings: Python, Julia, C++

Complete wrappers of rssn's embedding surface for three host languages, all
built on one ABI and held to one conformance table.

| Directory | Package | Entry point |
|---|---|---|
| [`python/`](python/) | `rssn` (ctypes) | `import rssn` |
| [`julia/`](julia/) | `RSSN` (ccall) | `using RSSN` |
| [`cpp/`](cpp/) | `rssn::client` (header-only) | `#include "rssn/client.hpp"` |
| [`api/`](api/) | generated API description | `rssn_api.json` |
| [`tools/`](tools/) | the two generators | `extract_api.py`, `generate.py` |
| [`tests/`](tests/) | shared conformance suite | `cases.tsv` |
| [`examples/`](examples/) | host callbacks, one program per language | `python/callbacks.py`, `julia/callbacks.jl`, `cpp/callbacks.cpp` |

## What "complete" covers

rssn exposes one small C ABI (`rssn.h`, ~50 entry points) and one large
operator table (5411 registered operators). A wrapper that only forwarded the
C functions would force every caller to spell operator names as strings, so each
binding covers both layers:

**The C ABI, one to one.**

- sessions and configuration: creation, release, symbolic/numeric target,
  symbol bindings, assumptions, search budget, dimension mode, dimension
  declarations, user base dimensions;
- terms: parse, symbol, integer, float, rational, apply-by-name, numeric
  evaluation, literal value, dense tensor extraction, and the four renderings
  (infix text, LaTeX, Typst, Unicode pretty);
- computation: single, batch, JSON record, `bincode-next` record, batch record;
- the value ABI: session-owned cells for dynamic values, `call_value` for the
  operators whose arguments have no term form, and readers for every value kind;
- the callback ABI: a host function registered as a `Value::Function` and
  passed to a higher-order operator, with the signature the operator declares
  or one the host spells out;
- the dimension report: inferred dimension text and the stable codes of the
  advisory findings;
- the embedding kernel: manifest JSON, manifest hash, version;
- simulations: scenario by name with JSON parameters;
- ownership: every string, byte buffer and tensor rssn hands out is released
  through the matching `rssn_*_free`, and every failure becomes an exception
  carrying the `RssnStatus` code.

**The operator table, one function per operator.** `python/rssn/_ops.py`,
`julia/src/ops.jl` and `cpp/include/rssn/ops.hpp` are generated from the kernel
manifest: 5411 functions, each with the operator's own documentation and arity.
Aliases are generated as forwarders. A name that is a keyword in the target
language is mangled (`and_` in Python, `and_op` in C++); the original spelling
stays reachable through `apply("and", ...)` and the `OPS` table.

`sub`, `neg` and `divide` are engine operators (the arithmetic quotient is
`divide` because `div` is the divergence), so `x - 1`, `-x` and `x/2` are real
requests rather than sugar the wrappers assemble.

**The data that travels with the manifest.** All 82 units with their exact
rational scale (`numer/denom * pi^pi`), affine offset, aliases and category; all
30 physical constants with exact value, dimension, unit and uncertainty; and the
29 simulation scenarios with the parameters their request types read.

## Quick start

Build the shared library once (this is the slow step; see
[the build note](#build-notes)):

```sh
bash scripts/build_ffi.sh                       # release profile, ~1h
RSSN_FFI_PROFILE=dev bash scripts/build_ffi.sh  # development profile, faster
```

Point a wrapper at the result with `RSSN_LIB` (or let it find a checkout's
`target/debug`, `target/ffi` or `target/release` on its own).

```python
import rssn

s = rssn.Session()
x = s.sym("x")

print(s.compute(rssn.diff(rssn.sin(x) * rssn.cos(x), x)).text)          # cos(2*x)
print(s.compute(rssn.diff(rssn.sin(x) * rssn.cos(x), x),
                rssn.Config().numeric(1e-12).bind("x", 0.0)).value)     # 1.0
```

```julia
using RSSN

session = Session()
x = sym(session, "x")

compute(diff(sin(x) * cos(x), x)).text                                  # "cos(2*x)"
cfg = bind!(numeric!(Config(), 1e-12), "x", 0.0)
compute(diff(sin(x) * cos(x), x), cfg).value                            # 1.0
```

```cpp
#include "rssn/client.hpp"
#include "rssn/ops.hpp"

rssn::client::Session session;
const auto x = session.sym("x");
const auto answer = session.compute(diff(sin(x) * cos(x), x));
std::cout << answer.text << "\n";                                       // cos(2*x)
```

## Host callbacks

189 operators are *higher-order*: they take a function and call it once per
quadrature node, optimisation step or root-finder probe (`num_bfgs`,
`num_brent`, `num_analysis_integrate_adaptive_scaled`, `num_bootstrap`, ...).
Each binding turns a function in its own language into the C function pointer
the engine calls, so these are ordinary calls:

```python
s.call_value("num_brent", 1.0, 2.0, 1e-12, 64, lambda x: x * x - 2.0)
rssn.num_brent(1.0, 2.0, 1e-12, 64, lambda x: x * x - 2.0)     # generated facade
```

```julia
call_value(session, "num_brent", 1.0, 2.0, 1e-12, 64, x -> x * x - 2)
call_value(session, "num_brent", 1.0, 2.0, 1e-12, 64) do x     # do-block sugar
    x * x - 2
end
```

```cpp
auto f = rssn::client::Callback::for_argument(session.raw(), "num_brent", 4,
    [](rssn::client::SessionHandle s, const std::vector<rssn::client::Value>& args) {
        const double x = args[0].as_f64();
        return rssn::client::value_number(s, x * x - 2.0);
    });
auto answer = num_brent(1.0, 2.0, 1e-12, 64, f);
```

Three pieces make it work:

- `rssn_callback_new_for(session, op, index, fn, user, out)` copies the
  signature the operator declares for that argument position, so a binding
  never guesses the host's types;
- `rssn_callback_new(session, signature_json, fn, user, out)` takes the
  signature explicitly (`{"params": ["real"], "ret": "real"}`);
- `rssn_callback_signatures(session, op, index, out_json)` lists what an
  operator accepts there, for the few that declare more than one function type
  at the same position (an objective and a gradient, say).

A callback receives its arguments as cells and reads them with the ordinary
`rssn_value_*` readers; it answers with a cell it built the same way. The
wrappers do that conversion for you, so a Python callback gets `float` /
`numpy` arrays / `dict`s and a Julia one gets `Float64` / arrays / `Dict`s.

## Tests

One table, three runners. `tests/cases.tsv` holds 27 cases — calculus, algebra,
number theory, linear algebra, tensors, numeric binding, the three failure
statuses, two simulations, the manifest hash and the unit count. The Python
runner records the expectations; the Julia and C++ runners verify the same rows,
so a wrapper that disagrees with the engine *or with another wrapper* fails.

```sh
export RSSN_LIB=$PWD/target-bindings/debug/rssn.dll     # or librssn.so

python bindings/tests/run_conformance.py                # verify
python bindings/tests/run_conformance.py --record       # rewrite the expectations
python -m unittest discover -s bindings/python/tests    # Python unit tests
julia --project=bindings/julia bindings/julia/test/runtests.jl
cmake -S bindings/cpp -B build/cpp -DRSSN_ROOT=$PWD       -DRSSN_LIBRARY=$PWD/target-bindings/debug/rssn.dll.lib
cmake --build build/cpp && ctest --test-dir build/cpp --output-on-failure
```

`bindings/tests/run_all.sh` runs all four in order.

Each language also has its own unit tests for the parts the table cannot see:
`ctypes`/`ccall` marshalling, the exception mapping, the Julia JSON codec,
Python operator overloading, C++ RAII and move semantics.

## Regenerating the generated files

```
target/manifest.json  ──extract_api.py──▶  api/rssn_api.json  ──generate.py──▶  10 files
```

```sh
python bindings/tools/extract_api.py    # refresh api/rssn_api.json from the manifest
python bindings/tools/generate.py       # rewrite the operator facade, units, constants
python bindings/tools/extract_api.py --check   # exit 1 when stale
python bindings/tools/generate.py --check      # exit 1 when stale
```

`bindings/api/julia_base_names.txt` and `julia_base_functions.txt` are one-off
scans of `names(Base)` (1.13) that the generator uses to decide which operator
names Julia can export as-is, which it extends `Base` with, and which it must
mangle. Regenerate them with

```sh
julia -e 'print(join(sort(string.(names(Base))), "
"))' > bindings/api/julia_base_names.txt
julia -e 'for n in names(Base); v = try getproperty(Base, n) catch; continue end; v isa Function && println(n); end' > bindings/api/julia_base_functions.txt
```

`extract_api.py` reads `bindings/api/manifest.json` when it exists, else
`target/manifest.json`. The manifest is the engine's own description of itself,
so the wrappers cannot drift from the operator table without the check failing.
The scenario catalog is the one thing the manifest does not carry: it is scanned
from `src/sim/scenario.rs` and `src/sim/scenario_engines.rs`, including the
parameter names and the literal defaults their `Default` implementations write.

## Design notes

**Operators are requests.** `diff(f, x)` returns an unevaluated term; `compute`
reduces it. This is the engine's own model (`call_value` vs `compute`), and it
is why the wrappers never "evaluate" behind your back: `str(term)` shows the
request, `session.compute(term).text` shows the answer.

**Handles, not objects.** A term is a `uint32_t` in one session. The wrappers
keep the session alive for as long as a term can reach it (Python keeps a
reference, Julia a field, C++ a raw pointer plus documented lifetime) and turn
the C ABI's `RSSN_NO_TERM` into an exception.

**Sessions are single-threaded.** The Rust store sits behind a `RefCell`, so two
threads in one session is a borrow conflict, which aborts in a release build.
Python and Julia give each thread/task its own default session; C++ has a
`thread_local` one. `compute_all` is the supported way to parallelise inside a
session.

**Failures are codes.** Every wrapper raises an exception carrying the
`RssnStatus` code, and never parses the message to decide what happened — the
engine does not promise wording.

## Known limits

These are the C ABI's limits, documented in `src/ffi/mod.rs`, not the wrappers':

1. **Every one of the 5411 operators is callable from C**, and
   `backend_support.ffi` is `true` for all of them. An operator is reached by
   the ABI its arguments need: `apply` when they are terms, `call_value` when
   they are native objects, records or other values with no term form, and a
   host callback when they are functions. The 189 higher-order operators are
   the last group.
2. **A value with no term form cannot become a term.** `to_term` refuses a
   native handle, a record and a host callback; `value_kinds` in the manifest
   says which kinds do have a graph form. Sparse tensors and complex entries are
   likewise read through the value readers, not through `rssn_term_tensor_data`.
3. **`bincode-next` buffers are opaque.** The wrappers return the bytes and
   document the record layout; they do not decode them. The JSON record carries
   the same fields and is the convenient path.
4. **The manifest is the standard snapshot.** It describes the process-wide rule
   sets, not a session's plugins (a C host cannot install plugins anyway).
5. **Panic containment is a build property.** The `ffi` profile inherits the
   release profile, whose `panic = "abort"` means a Rust panic aborts the host
   process. A `dev` build unwinds and reports `Panic` instead.

## Build notes

`cargo build` produces an `rlib` only; the C libraries come from
`scripts/build_ffi.sh` (`cargo rustc --lib --crate-type staticlib,cdylib`).
The default `ffi` profile inherits release optimizations with LTO and takes
about an hour on a development machine; `RSSN_FFI_PROFILE=dev` is much faster
and is what the tests here use. Set `CARGO_TARGET_DIR` to keep that build out of
a `target/` another process is using.

The wrappers never call Cargo themselves. They load `rssn.dll` / `librssn.so` /
`librssn.dylib` from `RSSN_LIB`, from `rssn/_lib/`, from a checkout's build
directories, or from the platform loader path.
