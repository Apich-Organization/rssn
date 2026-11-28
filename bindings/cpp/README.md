# rssn for C++

[English](README.md) | [简体中文](README.zh-CN.md)

A header-only RAII layer over the rssn C ABI, plus one function per registered
operator. C++17, no dependencies beyond the standard library and `rssn.h`.

```cpp
#include <iostream>
#include "rssn/client.hpp"
#include "rssn/ops.hpp"
#include "rssn/value.hpp"    // the value ABI

int main() {
    rssn::client::Session session;
    const auto x = session.sym("x");

    const auto answer = session.compute(diff(sin(x) * cos(x), x));
    std::cout << answer.text << "\n";                     // cos(2*x)

    rssn::client::Config numeric;
    numeric.numeric(1e-12).bind("x", 0.0);
    std::cout << *session.compute(diff(sin(x) * cos(x), x), &numeric).value << "\n";
}
```

## Build

The client needs `rssn.h` on the include path (the repository root) and the C
library:

```sh
bazel build //:ffi_dev
bazel test //:ffi_native_clients_test
```

Or by hand:

```sh
g++ -std=c++17 -I bindings/cpp/include -I . \
    bindings/cpp/tests/smoke.cpp -L<dir> -lrssn -o smoke
```

CMake accepts a full shared-library path or a directory via `RSSN_LIB`.
Pass `-DRSSN_LIBRARY=/path/to/link/library` to choose a link library explicitly.

| Stock archive | Shared linkage | Static archive |
| --- | --- | --- |
| Linux x86-64 and ARM64 | `librssn.so` | `librssn.a` |
| macOS x86-64 and ARM64 | `librssn.dylib` | `librssn.a` |
| Windows x86-64 and ARM64 (MSVC) | `rssn.dll.lib` (import library) and `rssn.dll` (runtime) | `rssn.lib` |

On Linux, select the static archive with
`cmake -S bindings/cpp -B build-cpp-static "-DRSSN_LIBRARY=$PWD/librssn.a"`
from the extracted root. The `rssn_client` CMake target supplies the
platform math, thread and dynamic-loader libraries. The package tests exercise
static linking on Linux and shared linking on all three operating systems.
For a macOS or Windows static link, supply the target's native system libraries
from the RSSN Rust build's `--print native-static-libs` output.

For MinGW, generate an import library from `rssn.dll`.

```sh
gendef rssn.dll && dlltool -d rssn.def -l librssn.dll.a
```

## Shape of the API

| Type | Role |
|---|---|
| `rssn::client::Session` | owns every term; move-only; single-threaded |
| `rssn::client::Term` | a handle; `.text()`, `.latex()`, `.typst()`, `.pretty()`, `.value()`, `.tensor()` |
| `rssn::client::Config` | target, bindings, assumptions, budget, dimension mode (fluent) |
| `rssn::client::Answer` | `.term`, `.text`, `.value`, `.error`, `.reduced`, `.iterations`, `.dimension`, `.diagnostic_codes`, `.latex()`, `.tensor()` |
| `rssn::client::diff(f, x)` | a generated operator function; returns an **unevaluated** term |
| `rssn::client::kOps` | the operator table: name, arity, `ffi` reachability |
| `rssn::client::op_spec(name)` | one lookup into that table |
| `rssn::client::Value` | a dynamic value (`rssn/value.hpp`): `kind()`, `as_f64()`, `text()`, `json()`, `to_term()`, `tensor()`, `size()`, `at()` |
| `rssn::client::kUnits`, `lookup_unit(name)` | the unit table |
| `rssn::client::kConstants`, `lookup_constant(name)` | the constant table |
| `rssn::client::kScenarios`, `lookup_scenario(name)` | the simulation catalog |
| `rssn::client::Session::run_scenario(name, params)` | run a simulation |
| `rssn::client::Error` | an exception with `.status()` |

Operators are *requests*: `diff(f, x)` builds the derivative,
`session.compute(...)` reduces it. `session.apply("name", args, count)` reaches
an operator by name at run time.

## Workflow: symbolic calculus

```cpp
rssn::client::Session session;
const auto x = session.sym("x");

const auto request = diff(sin(x) * cos(x), x);   // a request, not a result
std::cout << request.text() << "\n";            // diff(cos(x)*sin(x), x)
std::cout << session.compute(request).text << "\n";   // cos(2*x)

std::cout << session.compute(integral(x * x, x)).text << "\n";          // 1/3*x^3
std::cout << session.compute("defint(x^2, x, 0, 1)").text << "\n";      // 1/3
std::cout << session.compute("taylor(sin(x), x, 0, 5)").text << "\n";
std::cout << session.compute("limit(sin(x)/x, x, 0)").text << "\n";     // 1
```

A `Term` overloads the arithmetic operators, and a number or a string converts
into one implicitly, so a request reads like the mathematics:

```cpp
x * x + 2 * x + 1                  // a request; compute() gives "x^2 + 2*x + 1"
session.compute("diff(x^2, x)")    // text is parsed in the session that computes it
```

## Workflow: a numeric target

```cpp
rssn::client::Config cfg;
cfg.numeric(1e-12).bind("x", 0.5);
const auto answer = session.compute(diff(sin(x), x), &cfg);
answer.value      // std::optional<double>: 0.877582561890373
answer.error      // the estimated absolute error

session.compute_all({diff(x * x * x, x), integral(x * x, x), session.parse("gcd(4, 6)")});
```

## Workflow: linear algebra and tensors

```cpp
const auto matrix = session.parse("matmul(list(list(1,2),list(3,4)), list(list(0,1),list(1,0)))");
std::vector<std::size_t> shape;
const auto entries = session.compute(matrix).tensor(&shape);   // row-major, 4 doubles
```

## Workflow: NumPy arrays

With NumPy support in the native build, pass the extents and flat row-major
samples to the codec. The decoded value is a native tensor.

```cpp
const auto payload = session.encode_npy({2, 2}, {1.0, 2.0, 3.0, 4.0});
std::vector<std::size_t> extents;
const auto restored = session.decode_npy(payload).tensor(&extents);
```

## Workflow: a simulation

```cpp
#include "rssn/client.hpp"
#include "rssn/scenarios.hpp"

const auto result = rssn::client::Session::run_scenario("ising",
    R"({"width":8,"height":8,"temperature":2.0,"mc_steps":10,"seed":3})");
// `result` is the scenario's JSON text.

const rssn::client::ScenarioDef* scenario = rssn::client::lookup_scenario("poisson_multigrid");
scenario->params    // "dim:int=2; n:int=33; sigma:float=0.0; cycle:str=v; tol:float=1e-10"
```

## Workflow: units and constants

```cpp
#include "rssn/units.hpp"
#include "rssn/constants.hpp"

const auto* foot = rssn::client::lookup_unit("ft");
foot->dimension                     // "Length"
rssn::client::constant_value(*rssn::client::lookup_constant("c"))   // 299792458
```

For an integer beyond the signed 64-bit range, construct a value cell from its
decimal spelling and read it with `Value::text()`.

```cpp
const auto exact = rssn::client::value_integer_decimal(
    session.raw(), "1267650600228229401496703205383");
exact.text();
```

## Workflow: choosing a backend

The loaded library uses automatic selection until the process default is set explicitly.
Query the host's usable backends before changing the process default.

```cpp
const auto available = session.available_backends();
const auto active = session.current_backend();
session.set_default_backend("interpreter");
session.set_default_backend("auto");   // restore the engine's automatic policy
```

For explicit native compilation, pass the handle from `session.open_backend(name)`
to `session.call_value("num_backend_compile", ...)`. Use
`session.release_backend(handle)` to release the native handle when finished.

## Workflow: passing measured floats

The binary64 values below are already in coherent units. The dimension
object follows the kernel manifest's `DimExpr` JSON format. Attaching a
dimension preserves tensor shape and does not apply a unit scale or offset.

```cpp
const auto dimension = R"({"base":"Length"})";
const auto differences = session.value(std::vector<double>{2.0, 3.0}).with_coherent_dimension(dimension);
const auto initial = session.value(std::vector<double>{0.0}).with_coherent_dimension(dimension);
const auto positions = session.call_value("num_timeseries_integrate", differences, initial);
// positions.tensor() contains 0.0, 2.0 and 5.0 in metres.
```

## Workflow: a checked dynamic operator

Use the operator in `rssn::client::checked` with the session as its first
argument. It accepts the native inputs of that signature and returns a
session-owned `Value`.

```cpp
const auto lerch_value = checked::num_special_ext_lerch(session, -1234.0, 2.0, 1.0);
lerch_value.as_f64();  // 0.0218615893044645
```

## Workflow: finding an operator

The documentation of all registered operators is generated from the kernel manifest
into `rssn/ops.hpp`, as a comment above each function, together with the
`kOps` table (`name`, arity, reachability) that `op_spec()` searches:

```cpp
const rssn::client::OpSpec* spec = rssn::client::op_spec("diff");
spec->min_arity;    // 2
spec->max_arity;    // 2
spec->ffi;          // true
```

## Notes

- `Term` is a `(RssnSession*, uint32_t)` pair, so it is copyable and cheap, and
  it stays valid as long as the session lives — exactly the C contract. Moving a
  `Session` does not invalidate terms.
- A number or a string converts implicitly into a `Term` in the thread's default
  session, which is what makes `x * 2.0` and `diff("x^2", "x")` read naturally.
  A string argument to `Session::compute` is parsed in *that* session instead,
  so `session.compute("x^2")` never mixes two sessions.
  When an argument comes from a *different* session, `apply` copies it in
  through its text form (exact, because the engine's text round-trips).
- C++ keywords are mangled with an `_op` suffix: `and_op`, `not_op`, `or_op`,
  `xor_op`, `true_op`, `false_op`. `apply` takes the original spelling, and
  `op_spec("and")` finds it.
- Arithmetic on `Term` goes through the engine's own `sub`, `neg` and `divide`
  operators (the quotient is `divide` because `div` is the divergence).
- Operators whose arguments are native objects are called through the value ABI,
  in `rssn/value.hpp`:

      std::cout << session.current_backend() << '\n';

  `op_spec(name)->ffi` is `true` for every operator: the value ABI carries
  native objects and host callbacks alike. For record outputs, `Value::field(name)`
  retrieves a single field and `Value::field_name(index)` reads its label.

      const auto moments = session.call_value(
          "num_timeseries_rolling_moments", std::vector<double>{1.0, 3.0, 5.0}, 2, 0);
      const auto variance = moments.field("variance").tensor();  // {1.0, 1.0}

- A higher-order operator takes a function; `rssn/callback.hpp` wraps a
  `std::function` in the C function pointer the engine calls:

      auto f = rssn::client::Callback::for_argument(session.raw(), "num_brent", 4,
          [](SessionHandle s, const std::vector<Value>& args) {
              const double x = args[0].as_f64();
              return rssn::client::value_number(s, x * x - 2.0);
          });
      // `f` must outlive the call: destroying a Callback retires the
      // registration, and the cell it handed over stops being live.
- `#include "rssn.h"` puts the C declarations in `namespace rssn`; the client
  lives in `rssn::client`. `rssn.hpp` and `rssn.h` share an include guard, so
  including both is harmless but redundant.

## Tests

```sh
ctest --test-dir build/cpp --output-on-failure
```

`tests/smoke.cpp` is the tour (and a self-check); `tests/workflow.cpp` is the
workflows this README shows, compiled and run; `../tests/conformance.cpp` is the
cross-language table runner.
