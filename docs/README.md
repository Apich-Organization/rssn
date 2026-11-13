# Documentation index

Start with [../ARCHITECTURE.md](../ARCHITECTURE.md) for the overall design and
[../README.md](../README.md) for features and how to run the tests.

| Document | Subject |
|---|---|
| [BACKENDS.md](BACKENDS.md) | Execution backends and devices: tape, JIT, CPU and wgpu devices, source emitters, double-single `f64` on the GPU, launch fusion, how to add a backend, known gaps |
| [BUILDING.md](BUILDING.md) | Building the Rust library and the C artifacts, generated headers, required source files (`signatures.overlay.json`, `signatures.witnessed.json`), build metadata |
| [CALLBACKS.md](CALLBACKS.md) | Checked host callbacks through the dynamic ABI |
| [DIMENSIONS.md](DIMENSIONS.md) | The dimension and unit layer: design, algebra, API, coverage protocol, known limits |
| [KERNEL.md](KERNEL.md) | The embedding kernel: wire format, manifest, typed and witnessed signatures, dynamic calls, plugins, units, structured results |
| [MATH_COVERAGE.md](MATH_COVERAGE.md) | Edge-case conventions of the numerical kernels, where each area is tested, known gaps, space-dimension naming and coverage, kernels exposed as operators |
| [NUMERICAL_CALLBACKS.md](NUMERICAL_CALLBACKS.md) | The `numerical_callback` operators (quadrature, Brent, least squares, ODE integrators) |

Outside this directory: [../CONTRIBUTING.md](../CONTRIBUTING.md) (workflow and
checks), [../CHANGELOG.md](../CHANGELOG.md), [../SECURITY.md](../SECURITY.md),
and the generator documentation in
[../tools/rssn-api-gen/README.md](../tools/rssn-api-gen/README.md).
