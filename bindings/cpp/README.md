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
| `rssn::client::Error` | an exception with `.status()` |

Operators are *requests*: `diff(f, x)` builds the derivative,
`session.compute(...)` reduces it. `session.apply("name", args, count)` reaches
an operator by name at run time.

## Notes

- `Term` is a `(RssnSession*, uint32_t)` pair, so it is copyable and cheap, and
  it stays valid as long as the session lives — exactly the C contract. Moving a
  `Session` does not invalidate terms.
- A number or a string converts implicitly into a `Term` in the thread's default
  session, which is what makes `x * 2.0` and `diff("x^2", "x")` read naturally.
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
- `#include "rssn.h"` puts the C declarations in `namespace rssn`; the client
  lives in `rssn::client`. `rssn.hpp` and `rssn.h` share an include guard, so
  including both is harmless but redundant.

## Tests

```sh
ctest --test-dir build/cpp --output-on-failure
```

`tests/smoke.cpp` is the tour (and a self-check); `../tests/conformance.cpp` is
the cross-language table runner.
