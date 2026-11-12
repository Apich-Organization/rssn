# Execution backends and devices

rssn runs compiled numeric code along two independent axes. Keep them apart
when you add something.

| Axis | Trait | Question it answers | Implementations |
|---|---|---|---|
| **How is one term turned into a function?** | `backend::Backend` (`Compiled`, `CompiledMulti`) | "evaluate this closed form at a point" | `Interpreter`, the tape `TapeFunction`, Cranelift JIT (`jit`), `TieredBackend` |
| **Where do many rows run?** | `backend::device::Device` | "evaluate this program on a million rows / integrate ten thousand initial conditions" | `cpu`, `cpu-parallel`, `gpu:wgpu:N` (feature `gpu-wgpu`) |

`device::DeviceBackend` joins them: a `Backend` whose compiled functions send
batches to a `Device` and fall back to the host tape for programs a device cannot
run.

## Audit: what existed, what was missing

Honest state before the device layer:

* One scalar abstraction (`Backend::compile -> Compiled: f64^n -> f64`). Its
  batch methods (`call_batch`, `call_batch_par`) were a loop over rows, or rayon
  chunks over that loop. There was no notion of memory that lives somewhere else,
  of a launch, or of a reduction.
* The tape IR was the right seed for a kernel IR (SSA, hashed, constants folded)
  but it also carries instructions that only make sense on the host: `Call` (an
  operator's `EvalFn` pointer) and `Extern1`/`Extern2` (C function pointers).
  Nothing marked which tapes were portable.
* The JIT (Cranelift) is a scalar CPU code generator. It is not a template for
  GPUs: its `unsafe` is about calling machine code, not about buffers.
* "Add a backend" in `ARCHITECTURE.md` meant "implement `Backend`", which for a
  GPU means nothing: a GPU is not a better way to evaluate one point.

## Architecture

```
Term --lower--> Tape --KernelIr::from_tape--> KernelIr --Device::compile--> CompiledKernel
                  |                              |                              |
           (host: Call/Extern,           (portable: arithmetic,        Device::launch(rows,
            autodiff, intervals)          intrinsics, select, gather)   input, params, output)
                                                 |
                                                 +--> wgsl::emit_wgsl  (source emitters)
```

Module layout (`src/backend/device/`):

| File | Contents |
|---|---|
| `mod.rs` | `Device` (`launch`, `launch_cols`, `copy`, ...), `DeviceInfo`/`DeviceKind`/`DeviceOps`, `Buffer`, `DeviceRegistry`, `registry()`, `build_devices()`, `map_columns`, `ensemble_rk4`, `DeviceBackend` |
| `kernel.rs` | `KernelIr` (elementwise core, index values, loads, one fold), `KernelOp`, `Fold`, `Analysis`, `KernelBuilder`, `rk4_step`, the reference evaluators `eval_row` / `eval_at` |
| `array.rs` | array programs as kernel builders and their drivers: reductions, fused map-reduce, matmul, CSR product, stencils, all-pairs, `ensemble_dopri5` |
| `cpu.rs` | `CpuDevice`: lane-block interpreter for elementwise kernels, the reference evaluator for structured ones, rayon splitting |
| `emit.rs` | one emitter, five targets: WGSL, GLSL, OpenCL C, CUDA C++, Metal |
| `wgsl.rs` | `KernelEmitter` trait, `emit_wgsl`, `SourceEmitter` |
| `gpu.rs` | `WgpuDevice` (feature `gpu-wgpu` only) |
| `golden/*` | expected emitter output per target |
| `tests.rs` | all tests above |

### Capability-described devices

`DeviceInfo` is data, so a caller can choose without trying:

* `kind`: `Cpu`, `CpuParallel`, `Gpu(GpuApi)`, `Custom(name)`;
* `dtypes`: `F64` and/or `F32`. **The kernel IR has `f64` semantics;** a device
  that only has `f32` says so and its results are `f32`-rounded;
* `ops`: `elementwise`, `reduce`, `matmul`, `structured` flags; an optional
  method a device does not implement returns `DeviceError::Unsupported`
  (nothing is emulated silently). `structured` is the one that matters: a
  device with it runs every program in `array`;
* `parallelism`, `max_buffer_elems`, `unified_memory`, `deterministic`.

`registry()` lists them (`list()`), resolves by id (`get("cpu-parallel")`), and
picks a default (`default_device()`: first GPU, else `cpu-parallel`). A build's
*possible* kinds are `build_devices()`, which is what the kernel manifest
publishes in `features.devices`; whether a GPU is *present* is a run-time
question for the registry.

### Memory model

* `Buffer` is a typed array owned by one device. Its contents are visible to
  the host only through `Device::download`; `upload` copies in. Nothing
  converts implicitly, and a buffer handed to a device that does not own it is a
  `DeviceError::Mismatch`.
* Layout is column-major (structure of arrays): `input[i * rows + row]`,
  `output[k * rows + row]`. One flat buffer carries any number of columns, which
  is also one storage binding on a GPU.
* `params` are a host slice delivered with each launch (a uniform).

### Launch interface

```rust
let dev = backend::device::registry().get("cpu-parallel").unwrap();
let ir  = KernelIr::from_tape(&tape)?;            // portable?  else Unsupported
let k   = dev.compile(&ir)?;                       // cached by structural hash
let inp = dev.upload(&flat_columns, Dtype::F64)?;
let mut out = dev.alloc(ir.outputs.len() * rows, Dtype::F64)?;
dev.launch(&*k, rows, &inp, &params, &mut out)?;
dev.download(&out, &mut host)?;
```

Convenience layers that use only the trait:

* `map_columns(device, ir, columns, params)`: upload, launch, download.
* `dev.reduce(ReduceOp::Sum | Prod | Min | Max, &buf)`: fixed 4096-element
  chunks combined in order, so the result does not depend on thread count.
* `dev.matmul(a, b, c, m, k, n)`: row-major.
* `ensemble_rk4(device, rhs, params, y0, t0, t1, steps)`: the "same program,
  many initial conditions" pattern. `rk4_step` *inlines the right-hand side four
  times into one kernel*, so each step is one launch on any device; trajectories
  are independent, so a row's answer does not depend on the batch.

### Structured kernels: one primitive instead of a list of operations

The first device layer could run an elementwise formula per row and nothing
else; `reduce` and `matmul` were separate trait methods that only the CPU
implemented, and every further pattern (a stencil, a sparse product) would
have meant another method and another hand-written shader per device. That
does not scale to "broader GPU support", so the IR was extended instead of
the trait.

`KernelIr` keeps the tape as its elementwise core and adds three things
(`device/kernel.rs`):

* **index values**: `Row`, `Rows`, `Len(col)`, integer constants and integer
  arithmetic. Registers are typed (`Ty::Num` / `Ty::Idx`) and never mixed
  without `ToNum` / `ToIdx`;
* **addressed reads**: `Load(col, idx)` reads any element of an input column.
  Columns have independent lengths, so one launch can read a grid, a table and
  a CSR structure;
* **one fold**: `Fold { count, accs }` runs a body `count` times with the
  inner index `Inner` and combines values into accumulators (`Sum`, `Prod`,
  `Min`, `Max`); `Acc(k)` reads one afterwards.

Phases are *derived*, not written: `KernelIr::analyse` assigns every register
to `Pre` (once per row), `Body` (once per inner index) or `Post` (after the
loop) from its data flow, and rejects a register that needs both the inner
index and an accumulator. Loop-invariant code is therefore hoisted for free.

A device implements **one** entry point for all of it:

```rust
fn launch_cols(&self, kernel, rows, input: &Buffer, lens: &[usize], params, output: &mut Buffer)
```

`launch` is the special case "every column has one element per row". The
default `launch_cols` forwards exactly that case and returns `Unsupported`
otherwise, so a device written against the old trait keeps compiling and says
what it cannot do (`DeviceOps::structured`).

Everything else is a **library of kernel builders** (`device/array.rs`) that
uses only the trait:

| pattern | builder | rows | the fold |
|---|---|---|---|
| reduction | `reduce_pass` | `ceil(n / 256)` | one block of the column |
| map, then reduce, fused | `map_reduce_pass` | `ceil(n / 256)` | evaluates the map inside the block: the mapped column is never stored |
| matrix product | `matmul_kernel` | `m * n` | a row times a column |
| CSR matrix times vector | `spmv_csr_kernel` | matrix rows | the stored entries of a row |
| stencils (zero / clamp / periodic) | `stencil1d_kernel`, `stencil2d_kernel` | grid points | none: a fixed set of loads |
| all-pairs interaction | `pairwise_kernel` | bodies | every partner |
| adaptive ODE ensemble | `dopri5_step` | trajectories | none: one embedded step attempt with a per-row step size |

with drivers `array::reduce`, `reduce_column` (one column of a wider buffer,
in place), `map_reduce`, `matmul`, `run_cols` and `ensemble_dopri5`. The
`gpu-wgpu` device's `reduce` and `matmul` *are* these drivers: no shader was
written for either. `array::reduce` combines in a fixed tree (blocks of 256,
level by level), so two devices with the same float type agree bit for bit.

`ensemble_dopri5` integrates many initial conditions with adaptive
Dormand-Prince 5(4) steps. Time and step size are columns, so a stiff row
takes small steps while the others take large ones; the host asks only for
the minimum of the time column every few launches (one number per check).

`Device::copy` moves a range between two buffers of one device (wgpu: a
buffer-to-buffer copy; the default goes through the host and is still
correct).

### What is portable, and what is not

`KernelIr::from_tape` rejects `Call`, `Extern1` and `Extern2` with
`DeviceError::Unsupported`. Special functions that rssn attaches as `Lowering`
externs (gamma, erf, Bessel, ...) therefore stay on the host; `DeviceBackend`
falls back to the tape interpreter for such terms and still compiles. Moving them
onto a device means writing them in kernel instructions or in the target
language (see "Adding intrinsics" below).

### The CPU reference devices

`cpu` and `cpu-parallel` are bit-identical to each other and to
`KernelIr::eval_row` (which equals `Tape::eval`), including NaN, `select`
and out-of-range `gather`. Registers are 128-row blocks; each instruction is one
loop over a block, which the compiler vectorises for arithmetic. Transcendentals
call `libm` per lane (no SIMD math library is assumed). `cpu-parallel` halves the
rows recursively with `rayon::join` over disjoint `&mut` slices down to 2048
rows: no `unsafe`, no atomics.

### Source emitters and the `gpu-wgpu` feature

`emit::emit(&KernelIr, &EmitOptions) -> String` walks a program once and asks
a `Target` how to spell each construct; the hard part (phases, the fold,
bounds-checked loads, NaN conventions) exists once.

| target | language | precisions | status |
|---|---|---|---|
| `Target::Wgsl` | WGSL | `f32` | executed by `gpu-wgpu`; validated by naga in the tests |
| `Target::Glsl` | GLSL 4.50 compute | `f32`, `f64` | emitter, golden-tested |
| `Target::OpenCl` | OpenCL C | `f32`, `f64` | emitter, golden-tested |
| `Target::Cuda` | CUDA C++ | `f32`, `f64` | emitter, golden-tested |
| `Target::Metal` | Metal Shading Language | `f32` | emitter, golden-tested |

"Emitter" means the text is generated and pinned by golden files
(`golden/*.{glsl,cl,cu,metal}`); this crate has no runtime for those APIs and
their output has **not** been compiled by the vendors' compilers here. A device
for one of them supplies the launch half and reuses the emitter.

All targets share one interface: `rssn_in` (every column, concatenated),
`rssn_params`, `rssn_out` (column-major), `rssn_cols` (offset and length of
each input column) and the row and column counts; bindings 0 to 4 for WGSL,
GLSL and Metal, kernel arguments in that order for OpenCL and CUDA.
`wgsl::emit_wgsl` remains as the WGSL-only spelling.

`backend::device::gpu::WgpuDevice` (feature `gpu-wgpu`, off by default) compiles
the WGSL with wgpu, one cached pipeline per kernel, `f32` buffers, a
two-dimensional dispatch for more than `65535 * 64` rows, and runs structured
kernels (`ops.structured`, `ops.reduce`, `ops.matmul`). `gpu::probe()` returns
one device per adapter; none is not an error.

## Adding a backend, step by step

The same recipe covers CUDA, Metal, Vulkan (raw) and OpenCL. Everything that is
API-independent already exists; you add an *emitter* and a *device*.

1. **Pick the feature name** and add it to `Cargo.toml` as an optional
   dependency plus feature, additive only and off by default, e.g.
   `gpu-cuda = ["dep:cudarc"]`. Add the kind to `GpuApi` if it is not there
   (`Cuda`, `Metal`, `Vulkan`, `OpenCl` already are) and to `build_devices()` under
   `cfg!(feature = "...")` so the manifest reports it.
2. **Write the emitter** `src/backend/device/<target>.rs` implementing
   `wgsl::KernelEmitter` (PTX text, SPIR-V assembly, MSL, OpenCL C). Walk
   `KernelIr::ops` in order; op `d` defines value `v{d}`. Reuse the WGSL file as
   the template: one thread per row, column-major `in[i*rows+row]`, parameters
   from a constant buffer, helpers only when used. Decide `f32` or `f64`
   deliberately and write the differences from the tape in the module docs, as
   `wgsl.rs` does.
3. **Golden-test the emitter** before any hardware is involved: build small
   `Tape`s by hand, emit, compare with `golden/*.<ext>` (copy the
   `golden()` helper; the `RSSN_BLESS_GOLDEN` switch regenerates). If the target
   has an offline validator (naga for WGSL, `ptxas`/`spirv-val` where installed),
   run it in a test that skips when the tool is missing.
4. **Implement `Device`** in `src/backend/device/<api>.rs` behind the feature:
   * `info()`: fill `DeviceInfo` truthfully (`dtypes`, `ops`, memory limits,
     `unified_memory`). Claim only what is implemented.
   * `uid()`: `next_uid()` once per instance; every `Buffer` you create carries
     it.
   * `alloc`/`upload`/`download`: wrap the API buffer in a type implementing
     `Storage` and recover it with `as_any().downcast_ref`. Convert `f64` to your
     dtype on upload, back on download. Check `buffer.owner() == self.uid()`.
   * `compile`: `ir.validate()?`, emit source, build the pipeline/module, cache
     by `ir.structural_hash()` and confirm `cached.ir == *ir` (hash collisions).
     Pipeline-creation errors become `DeviceError::Backend`, never a panic.
   * `launch`: call `check_launch(self.uid(), ir, rows, input, output)?` first,
     then bind `input`, `params`, `output` and a dims record, dispatch
     `ceil(rows / group)` groups (spread over two dimensions past the API's
     per-dimension limit, as `gpu.rs` does), and submit. Launch may be
     asynchronous.
   * `download` and `synchronize` are the sync points; surface errors the API
     reported late there.
   * `reduce`, `matmul`: leave the defaults (`Unsupported`) until you implement
     them and set the `ops` flags to match.
5. **Discovery**: write `fn probe() -> Vec<Arc<dyn Device>>` that never fails (no
   adapter means an empty list) and call it from `DeviceRegistry::new()` under
   the feature, as `gpu::probe()` is. Ids are `gpu:<api>:<ordinal>`.
6. **Test with the shared contract**, not with new expectations:
   * `ensemble_rk4` of a decay equation against the closed form, tolerance
     chosen for your dtype;
   * `map_columns` of the kernel in `tests.rs::everything()` against
     `KernelIr::eval_row` within dtype tolerance (skip NaN rows on devices that
     do not guarantee them);
   * foreign-buffer and wrong-size launches return `Mismatch`;
   * gate any test that needs hardware on `probe()` returning a device so CI
     machines without one pass.
7. **Document** the precision, the unsupported ops and the feature in this file's
   table, and keep `unsafe` out of the device (if the API forces it, confine it
   to one small module with a safety comment per block, as the JIT does in its
   `abi` module).

### Adding intrinsics

A new elementary function enters at `backend::Intrinsic` (tape and CPU), then
each emitter's name table (`intrinsic_name` in `wgsl.rs`). Functions that the
target lacks (WGSL has no `cbrt`) get a helper emitted only when used. A special
function with no closed instruction sequence stays a host `Call`/`Extern` and
therefore host-only until someone writes its kernel.

## Known gaps

* GPU arithmetic is `f32`. There is no `f64` emulation (double-single) and
  WGSL has no `f64`; the GLSL, OpenCL and CUDA emitters produce `f64` source,
  but no device in this crate runs it.
* The GLSL, OpenCL, CUDA and Metal emitters are golden-tested text. Nothing
  here compiles them with the vendor toolchains; expect to fix details when
  the first device for one of them is written.
* Index data that a kernel reads from a column (CSR row pointers and column
  indices) is stored as floats: exact below `2^24` on an `f32` device. Typed
  integer columns are the fix and are not implemented.
* A kernel has one fold, not nested loops, and reads one input buffer. A
  product of two buffers is joined with `Device::copy` first.
* `array::reduce` is a tree of launches with a download of one number at the
  end; there are no on-device scans, sorts or compactions.
* Special functions that are host externs do not run on devices.
* Single-row calls through `DeviceBackend` pay a full launch; use batches.
* `ensemble_dopri5` has no dense output and no event detection, and the CPU
  path for structured kernels is the reference evaluator (correct, parallel,
  not vectorised).
* The `gpu-wgpu` path is validated with naga; execution is tested only where
  an adapter is present (see the test's skip message).
* Kernel launches are fused only where a builder does it (`map_reduce_pass`);
  there is no general planner that fuses a chain of launches.
