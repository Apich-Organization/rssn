# Documentation index

[中文](README.zh.md)

Start with [../ARCHITECTURE.md](../ARCHITECTURE.md) for the overall design and
[../README.md](../README.md) for features and how to run the tests.

| Document | Subject |
|---|---|
| [BACKENDS.md](BACKENDS.md) · [中文](BACKENDS.zh.md) | Execution backends and devices: tape, JIT, CPU and wgpu devices, source emitters, double-single `f64` on the GPU, launch fusion, how to add a backend, known gaps |
| [BUILDING.md](BUILDING.md) · [中文](BUILDING.zh.md) | Building the Rust library and the C artifacts, generated headers, required source files (`signatures.overlay.json`, `signatures.witnessed.json`), build metadata |
| [CALLBACKS.md](CALLBACKS.md) · [中文](CALLBACKS.zh.md) | Checked host callbacks through the dynamic ABI |
| [Dimensions and units](DIMENSIONS.md) · [中文](DIMENSIONS.zh.md) | Dimension algebra, conversion, inference and embedding contracts |
| [Embedding kernel](KERNEL.md) · [中文](KERNEL.zh.md) | The embedding kernel: wire format, manifest, typed and witnessed signatures, dynamic calls, plugins, units, structured results |
| [MATH_COVERAGE.md](MATH_COVERAGE.md) · [中文](MATH_COVERAGE.zh.md) | Edge-case conventions of the numerical kernels, where each area is tested, known gaps, space-dimension naming and coverage, kernels exposed as operators |
| [NUMERICAL_CALLBACKS.md](NUMERICAL_CALLBACKS.md) · [中文](NUMERICAL_CALLBACKS.zh.md) | The `numerical_callback` operators (quadrature, Brent, least squares, ODE integrators) |

Outside this directory: [../CONTRIBUTING.md](../CONTRIBUTING.md) (workflow and
checks), [../CHANGELOG.md](../CHANGELOG.md), [../SECURITY.md](../SECURITY.md),
and the generator documentation in
[../tools/rssn-api-gen/README.md](../tools/rssn-api-gen/README.md).
