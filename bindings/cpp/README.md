# rssn for C++

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
bash scripts/build_ffi.sh                       # produces rssn.dll + rssn.dll.lib

cmake -S bindings/cpp -B build/cpp \
      -DRSSN_ROOT=$PWD \
      -DRSSN_LIBRARY=$PWD/target/ffi/rssn.dll.lib
cmake --build build/cpp
ctest --test-dir build/cpp --output-on-failure
```

Or by hand:

```sh
g++ -std=c++17 -I bindings/cpp/include -I . \
    bindings/cpp/tests/smoke.cpp -L<dir> -lrssn -o smoke
```

`RSSN_LIBRARY` is whichever library matches the toolchain: the MSVC import
library `rssn.dll.lib`, the static archive `rssn.lib`, or a MinGW import library
built from the DLL:

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

session.compute_all({diff(x * x * x, x), integral(x * x, x), gcd(4, 6)});
```

## Workflow: linear algebra and tensors

```cpp
const auto matrix = session.parse("matmul(list(list(1,2),list(3,4)), list(list(0,1),list(1,0)))");
std::vector<std::size_t> shape;
const auto entries = session.compute(matrix).tensor(&shape);   // row-major, 4 doubles
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

## Workflow: finding an operator

The documentation of all 5411 operators is generated from the kernel manifest
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

      const auto backend = session.call_value("num_backend_current");
      std::cout << session.call_value("num_backend_name", backend).text() << "
";

  `op_spec(name)->ffi` is `true` for every operator: the value ABI carries
  native objects and host callbacks alike.
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
