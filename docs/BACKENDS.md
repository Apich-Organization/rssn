# RSSN execution backends and devices

[简体中文](BACKENDS.zh.md)

RSSN separates compilation of a numeric expression from execution of a kernel
over a batch. Choose a `Backend` to turn a term into a callable function. Choose
a `Device` to place data, select arithmetic precision and launch a program over
rows. A row can be a data point, a grid cell or the initial state of an
independent trajectory.

| Interface | Implementations and uses |
|---|---|
| `backend::Backend` | `Interpreter` is the default scalar backend. The shared tape implements `Compiled` and `CompiledMulti` through `TapeFunction`. With `jit`, `CraneliftBackend` compiles to CPU machine code and `TieredBackend` compiles scalar functions in the background after a configurable call or row threshold. |
| `backend::device::Device` | `cpu` and `cpu-parallel` run on the host. With `gpu-wgpu`, each usable wgpu adapter can register a `gpu:wgpu:N` device. |
| `backend::device::DeviceBackend` | Connects the interfaces. A batch uploads its inputs, launches the compiled device program and downloads its outputs. A single-point call uses a one-row launch. |

`Backend::compile` produces a scalar `Compiled` function of its declared
inputs. `compile_multi` produces a `CompiledMulti` program with shared
subexpressions, parameters constant across a batch, and multiple outputs.
`backend::current()` uses the innermost thread-local `with_backend` scope, then
the process default set by `set_default`, then the interpreter.
[Backend API](../src/backend/mod.rs) · [Cranelift and tiering](../src/backend/jit.rs)

A term lowers to a straight-line tape with shared subexpressions and
constant folding. The tape supports host operator calls, differentiation
and interval evaluation. `KernelIr::from_tape` translates portable arithmetic,
intrinsics, comparisons, selection and indexed gathering. It returns
`DeviceError::Unsupported` for host function pointers in `Call`,
`Extern1` and `Extern2`. A term that lowers successfully to a tape can still
run through `DeviceBackend` when kernel translation or device compilation fails.
In that case the backend uses `TapeFunction`. During device execution, the
infallible `Compiled` interfaces represent device errors as NaNs. Use the
explicit `Device` API for structured error handling.
[Tape](../src/backend/tape.rs) · [Kernel translation](../src/backend/device/kernel.rs) · [Device adapter](../src/backend/device/mod.rs)

## Device selection and buffers

`registry().list()` describes devices present at run time. Use `get(id)` for an
explicit choice or `default_device()` for the first registered GPU, then
`cpu-parallel` when the active Rayon pool has multiple workers, then `cpu` or
the first available device. wgpu may expose software adapters
and marks them in their names. An absent adapter leaves both CPU devices
available. The kernel manifest's `features.devices` uses `build_devices()` to
describe capabilities compiled into a build, independently of adapter discovery.
[Registry](../src/backend/device/mod.rs) · [Adapter discovery](../src/backend/device/gpu.rs) · [Kernel manifest](../src/kernel/manifest.rs)

The global Rayon pool detects available CPUs. `RAYON_NUM_THREADS` sets its
size, and a host can install a dedicated Rayon pool around its kernel work.
Parallel helpers use the active pool and skip dispatch below their work limits
or with one worker. Embedders with their own pools can query
`execution::default_cpu_threads()` to apply the host default without starting
the global pool. [CPU execution policy](../src/execution.rs)

Compiled CPU batches lease exclusive temporary `f64` buffers. Numerical
kernels can use the same reuse module for typed scratch and shared immutable
plans. Sized scratch requests can reuse an idle allocation that already fits.
Concurrent requests for the same cached plan share its construction. Returned
results and active scratch leases remain owned by their callers.
[Scratch and plan APIs](../src/kernels/reuse.rs)

`DeviceInfo` describes device kind, supported `Dtype` values and operations,
parallelism, buffer capacity, unified memory and determinism. Use `ops.structured`
to select programs with independent column lengths, addressed reads and a
fold. Optional operations report `DeviceError::Unsupported` when the device
offers no implementation. Use `emulated_f64` to distinguish the implementation
of the advertised `F64`, and `preferred_dtype` to choose `F64` when available or
`F32` otherwise.

A `Buffer` belongs to one device instance. `upload` converts a host `f64` slice
to the requested device type, `download` converts it back, and `copy` moves a
range between buffers of the same device and type. The default copy stages
through the host. The CPU and wgpu devices use their own memory copies.
Foreign buffers, incompatible types for a copy and invalid shapes report
`DeviceError::Mismatch`.

Input columns are concatenated in one buffer. An equal-length input column
`i` starts at `i * rows`. Output column `k` starts at `k * rows`. In
`launch_cols`, `lens` supplies the length of each input column and its prefix
sum gives the offset. All outputs still have `rows` elements. `launch` is the
equal-length form. Device implementations can support that form without
implementing structured launches. The common launch checks limit `rows` and
the total input size to `2^31 - 1` elements. Adapter limits can impose smaller
buffers. Parameters are supplied as a separate host slice for each launch.
[Device and buffer contracts](../src/backend/device/mod.rs)

For a portable tape and equal-length host columns, the complete transfer and
launch can use `map_columns`.

```rust
use rssn::backend::{device, tape::Tape};
use rssn::backend::device::{DeviceError, KernelIr};

fn evaluate(
    tape: &Tape,
    columns: &[&[f64]],
    params: &[f64],
) -> Result<Vec<Vec<f64>>, DeviceError> {
    let dev = device::registry()
        .get("cpu-parallel")
        .ok_or_else(|| DeviceError::Unsupported("cpu-parallel is unavailable".into()))?;
    let ir = KernelIr::from_tape(tape)?;
    device::map_columns(dev.as_ref(), &ir, columns, params)
}
```

For repeated launches, retain the compiled kernel and buffers. Allocate the
output with `ir.outputs.len() * rows` elements. Call `compile`, `upload`, `alloc`
and `launch` or `launch_cols`, then `download` when host code needs the result.
On wgpu, submission can be asynchronous. `download` and `synchronize` wait
for work and can report deferred errors. [wgpu device](../src/backend/device/gpu.rs)

## Kernel IR and array programs

`KernelIr` extends the elementwise tape with integer index registers, addressed
column reads and one optional fold. `Row`, `Rows`, `Len` and integer arithmetic
operate on `Ty::Idx`. Numeric values have type `Ty::Num`. `ToIdx` and `ToNum`
convert between the two. `Load(col, idx)` reads a cell of one column and
`Gather` selects an input column for the current row. A `Fold` repeats its
body using `Inner` and combines values into `Sum`, `Prod`, `Min` or `Max`
accumulators, which are subsequently read with `Acc`.

`KernelIr::analyse` validates operand types and derives the three execution
phases from data dependencies. `Pre` runs once per row and includes
loop-invariant work and the fold count. `Body` runs per inner index. `Post`
reads accumulated values. A register that depends on both the inner index and
an accumulator is invalid. `eval_row` and `eval_at` provide the host
reference semantics. [Kernel IR and validation](../src/backend/device/kernel.rs)

The array library expresses these algorithms using the same launch interface.

| Computation | Builder | Rows and work per row |
|---|---|---|
| Reduction | `reduce_pass` | One row per block of up to 256 elements, with a fold over the block. |
| Fused map and reduction | `map_reduce_pass` | The fold evaluates the map for each element and stores block partials. |
| Row-major matrix product | `matmul_kernel` | One row per result element, with a fold over its dot product. |
| CSR matrix-vector product | `spmv_csr_kernel` | One row per matrix row, with a fold over its stored entries. |
| One- and two-dimensional stencil | `stencil1d_kernel`, `stencil2d_kernel` | One row per grid point, with fixed loads and zero, clamp or periodic boundaries. |
| All-pairs interaction | `pairwise_kernel` | One row per body, with a fold over matching partners. |
| Adaptive ODE ensemble | `dopri5_step` | One row per trajectory, with one embedded step attempt using that row's time and step size. |

`array::run_cols` runs a structured kernel on host columns. `reduce` and
`reduce_column` combine 256-way blocks level by level and download a single
value. `map_reduce` combines map and first reduction pass in one launch.
`array::matmul` packs the two row-major operands into a single input buffer
using `Device::copy`. wgpu implements its reduction and matrix-product
methods through these drivers. The CPU device's dedicated `Device::reduce`
uses fixed 4,096-element chunks and combines them in order, independent of
thread count. Reduction order and floating-point behavior can differ across
devices and entry points. [Array kernels and drivers](../src/backend/device/array.rs) · [CPU reductions](../src/backend/device/cpu.rs)

`ensemble_rk4` runs a fixed number of classical fourth-order Runge–Kutta
steps. `rk4_step` inlines four right-hand-side evaluations into one kernel
launch per step. `ensemble_dopri5` integrates forward with adaptive
Dormand–Prince 5(4), keeping time and step size on the device for each
trajectory. The host reduces the time column to its minimum every four
step launches and downloads the final state columns. The return value
includes the number of step launches. An ensemble that fails to finish
within its launch budget reports `DeviceError::Backend`.
[Fixed-step driver](../src/backend/device/mod.rs) · [Adaptive driver](../src/backend/device/array.rs)

## CPU execution and LIR

The sequential and parallel CPU devices use native `f64`. For elementwise
kernels, the interpreter evaluates each instruction over blocks of 128
rows, exposing contiguous arithmetic loops to compiler vectorization.
Transcendental functions use scalar math operations. Structured kernels use
the reference evaluator. The parallel device splits disjoint output ranges
with Rayon. The grain is 2,048 rows for the elementwise interpreter and
native path and 128 rows for the structured interpreter.

With `jit`, CPU device launches from 256 rows try to compile a row loop on
first use and cache the result. `Gather` uses the interpreter. A native
bail, such as index overflow, reruns the launch with the interpreter.
The native and reference paths follow the same index and numeric conventions
in ordinary cases. An equal `Min` or `Max` fold over signed zeros can choose a
different zero sign. [CPU execution](../src/backend/device/cpu.rs) · [Native lowering](../src/backend/device/native.rs)

LIR is a CPU intermediate representation shared by tape batching, CPU kernels
and language front ends. It has basic blocks, jumps, branches and direct
module-local calls. Registers of type `F64`, `I64`, `Bool` and `Slice` are
mutable and zero-initialized. A jump assigns the target block parameters in
parallel. Host slices have bounded loads and stores. An out-of-range load
produces NaN and an out-of-range store has no effect. Validation traces
transitive writes so execution can reject read-only slice arguments before
entering a function that writes to them.

Checked integer arithmetic and float-to-integer conversion can return
`Outcome::Bailed`, as can a `Term::Bail` or a module call exceeding
`MAX_CALL_DEPTH` (512). Hosts can rerun a bailed computation when its reads
remain valid. `Module::interpret` provides a reference implementation and
`lir_jit` compiles modules with Cranelift when `jit` is enabled.
[LIR](../src/backend/lir.rs) · [LIR code generation](../src/backend/lir_jit.rs)

| Front end | LIR use |
|---|---|
| Tape | `lir::batch_from_tape` places eligible instructions in a loop over column slices. The Cranelift batch path attempts this from 64 rows. |
| `KernelIr` | `device::native::lower` lowers the row loop, fold and kernel index semantics to LIR for the CPU devices. |
| Aevia Core IR | Eligible pure scalar `@jit` functions use LIR for branches, checked `Int` arithmetic, recursion and self tail calls with fuel accounting. |

[Backend timing example](../examples/backend_timing.rs) compares reference,
tape row and column execution, and both Cranelift modes over 100,000 rows.
Run it with `cargo run --profile tools --features jit --example backend_timing`
to measure the executors on your machine.

## GPU source targets and precision

`emit::emit(&KernelIr, &EmitOptions)` generates source for the requested
target. The emitter shares validation, phase assignment, fold structure and
bounds-checking conventions across targets. `wgsl::emit_wgsl` is the
single-precision WGSL convenience entry point.

| Target | Arithmetic | Integration |
|---|---|---|
| WGSL | `f32` or double-single | wgpu execution with `gpu-wgpu`, naga validation tests and adapter-dependent execution tests. |
| GLSL 4.50 compute | `f32` or native `f64` source | Source export and golden tests. |
| OpenCL C 1.2 | `f32` or native `f64` source | Source export and golden tests. |
| CUDA C++ | `f32` or native `f64` source | Source export and golden tests. |
| Metal Shading Language | `f32` | Source export and golden tests. |

The latter four integrations supply source for an external runtime and
compiler. The golden tests compare generated text. WGSL validation tests
parse and validate emitted modules with naga without requiring an adapter.
[Emitter](../src/backend/device/emit.rs) · [Emitter and device tests](../src/backend/device/tests.rs)

The generator assigns one invocation to each row. WGSL and GLSL use five
bindings in group or set 0 for input, parameters, output, dimensions
(`rows, inputs, outputs, pad`) and column metadata (`offset, length` pairs),
in that order.
Metal uses the corresponding buffer indices. OpenCL and CUDA take input,
parameters, output, column metadata, row count and input column count as
kernel arguments. wgpu uses storage buffers for parameters and a uniform for
dimensions. It dispatches groups of 64 rows, splitting into two dimensions
beyond 65,535 groups on the first axis.

The wgpu device caches compiled kernels by structural hash and checks the
original IR on a cache hit. A compiled kernel starts with a single-precision
pipeline and builds its double-single pipeline on the first `F64` launch.
On wgpu, use one element type for the input and output of a nonempty launch.
The reported `max_buffer_elems` is conservatively counted in eight-byte
elements. [wgpu device](../src/backend/device/gpu.rs)

### Double-single arithmetic

`Dtype::F64` on wgpu uses two `f32` components with value `hi + lo`. Upload
splits host values and download joins them. Shader inputs, parameters and
outputs are arrays of eight-byte pairs. The format carries about 48 significant
bits and retains the exponent range of `f32`. `DeviceInfo::emulated_f64` marks
this representation. `dtype_epsilon` supplies `2^-23` for `F32`, `2^-52`
for native `F64` and `2^-44` for the pair implementation. The adaptive
ensemble uses these values to set its endpoint slack.

The pair arithmetic covers basic arithmetic, comparisons, square root,
powers, rounding and elementary functions. Range reduction and series or
Newton corrections evaluate transcendental functions in pair precision.
The adapter-backed tests compare representative elementary results with the
CPU at absolute tolerance `5e-12 * max(abs(reference), 1e-3)`. Treat the
precision as workload- and driver-dependent. For double-single data, keep
nonzero magnitudes and intermediate values within roughly `1e-30` to
`1e30`. Trigonometric range reduction loses absolute accuracy as the argument
grows and returns NaN for `sin`, `cos` and `tan` beyond `|x| = 1e6`.
GPU NaN and infinity behavior is implementation-dependent.

The pair library hides intermediate rounding points with `ds_hide`, a
bitwise operation whose zero comes from the dimensions uniform. This
preserves the rounding required by error-free transformations when shader
optimizers would otherwise reassociate or contract expressions. Preserve
that barrier when changing the arithmetic.

Kernel index registers on shaders use signed 32-bit integers. Index data
carried in numeric columns, such as CSR row pointers and column indices,
also has to fit the buffer's representation. `F32` loses unit resolution
above `2^24`. The double-single representation allows a wider exact
integer range, subject to the shader's index conversions and device bounds.
[Pair arithmetic](../src/backend/device/ds.rs) · [Precision tests](../src/backend/device/tests.rs)

## Launch fusion

`device::plan::Plan` constructs a graph over equal-length external columns and
outputs of earlier kernels. Mark the desired results with `result`, call
`fuse` to combine eligible launches, and call `run` on the chosen device. In
a chain such as `y = f(x); z = g(y, w)`, fusion can replace the intermediate
column `y` with a register. An output needed by another consumer or requested
as a final result can remain an output of the merged kernel.

Fusion substitutes same-row `Input` reads, renumbers parameters and registers,
and shares repeated input sources. A producer and consumer with two folds,
a consumer that uses `Gather`, or an addressed `Load` or `Len` on the
producer's column keeps its launch boundary. The merged kernel must validate
and preserve execution order. An earlier reader keeps the producer in its
original slot. For equal-length workloads the plan handles device copies
between launches. Use `launch_cols` for programs with different input
column lengths. [Planner and fusion tests](../src/backend/device/plan.rs)

## Implementing a device

1. Use a `DeviceRegistry` entry with a unique ID and report implemented
   precision, operations and limits through `DeviceInfo`. The built-in
   feature-gated devices use `gpu:<api>:<ordinal>` and publish build capabilities
   through `build_devices()`. Discovery returns an empty list when no device
   can be initialized.
2. Wrap storage in `Storage`. Use `next_uid()` to identify a device instance.
   Check owner, shape and type at the transfer and launch boundaries. Override
   `copy` when the device can move data without a host transfer.
3. Validate and compile `KernelIr`, reusing a source emitter where the target
   matches. Cache by structural hash and check IR equality before reuse.
   Surface validation and backend failures through `DeviceError`. For
   asynchronous runtimes, collect deferred failures in `download` and
   `synchronize`.
4. Implement and advertise the supported launch and reduction capabilities.
   Use `check_launch` or `check_launch_cols` to enforce the shared buffer
   contract. A structured device can implement reductions and matrix products
   through the array drivers. Compare numerical results with the reference at
   the declared precision and test foreign buffers and wrong sizes.
5. Keep source-emitter tests independent of hardware. Golden outputs live in
   the device directory and can be regenerated with `RSSN_BLESS_GOLDEN=1`.
   Validate generated code with the target compiler where available and
   run execution tests when discovery finds a suitable adapter. A new
   elementary `Intrinsic` needs tape and CPU semantics plus a mapping or
   helper per source target.
