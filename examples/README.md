# RSSN Examples

[English](README.md) | [简体中文](README.zh-CN.md)

Run a Rust example from the repository root with `bazel run //:example -- NAME`.

| Example | What it exercises |
|---|---|
| [quickstart](quickstart.rs) | Sessions, symbolic and numeric targets |
| [symbolic_calculus](symbolic_calculus.rs), [solve_equations](solve_equations.rs) | Rewrites, calculus, and equations |
| [units_and_dimensions](units_and_dimensions.rs) | Physical units and dimension checks |
| [linear_algebra](linear_algebra.rs), [numeric_kernels](numeric_kernels.rs) | Matrix factorizations, solvers, and real FFT reconstruction |
| [ode_pde](ode_pde.rs), [simulation_field](simulation_field.rs) | Evolution models and fields |
| [simulation](simulation.rs), [scenario_catalog](scenario_catalog.rs) | Named simulation dispatch and typed scenario metadata |
| [compile_expression](compile_expression.rs) | Reusable compiled evaluation |
| [weighted_calibration](weighted_calibration.rs) | Weighted quantiles and threshold confidence |
| [engine_timing](engine_timing.rs), [backend_timing](backend_timing.rs) | Symbolic workload and CPU backend timing |
| [io_formats](io_formats.rs) | HDF5 (chunked/filters/vlen/soft-links), Parquet (Brotli/Bloom/PageIndex), VTK, WAV, PNG, GIF & CSV |
| [sim_postprocess_vis](sim_postprocess_vis.rs) | 2D simulation post-processing to HDF5, VTK, PNG, SVG, GIF & WAV |
| [C smoke example](c/smoke.c) | Session handles and ownership across the C ABI |

Use the `tools` profile for timing; backend comparison also requires `jit`:

```sh
bazel run //:example -- --profile tools engine_timing
bazel run //:example -- --profile tools backend_timing
```

See the [build guide](../docs/BUILDING.md) for features and the C toolchain.
