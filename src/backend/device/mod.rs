//! # Devices: data-parallel execution
//!
//! [`Backend`](super::Backend) compiles one closed-form term to a *scalar*
//! function. A [`Device`] is the other axis: it runs one program over *many
//! rows* (the rows of a data set, the cells of a grid, the initial conditions
//! of an ensemble), on whatever hardware it represents. The two compose:
//! [`DeviceBackend`] is a `Backend` whose compiled functions evaluate batches
//! on a device.
//!
//! The layers, from the user down:
//!
//! 1. [`DeviceRegistry`] / [`registry`]: discovery. Lists [`DeviceInfo`]
//!    (kind, capabilities, precision) without touching device memory.
//! 2. [`Device`]: allocate, [`upload`](Device::upload) /
//!    [`download`](Device::download) ([`Buffer`] transfers are explicit and
//!    never implicit), [`compile`](Device::compile) a [`KernelIr`],
//!    [`launch`](Device::launch) it over `rows`, [`reduce`](Device::reduce),
//!    [`matmul`](Device::matmul).
//! 3. [`KernelIr`]: the backend-neutral program, lowered from a
//!    [`Tape`](super::tape::Tape). Everything a device compiles, whether Rust
//!    loops, WGSL ([`wgsl`]) or a future PTX/SPIR-V/MSL emitter, starts here.
//! 4. [`ensemble_rk4`]: "the same program over many initial conditions",
//!    expressed *only* through `launch`, so any device gets it for free.
//!
//! ## Data layout
//!
//! Kernel inputs and outputs are **column-major (structure of arrays)** in one
//! flat buffer: `input[i * rows + row]`, `output[k * rows + row]`. This is
//! what a GPU wants (coalesced loads), what a CPU auto-vectoriser wants, and
//! what lets one storage binding carry any number of columns. Parameters
//! (`Param`) are a small host slice delivered with every launch (a uniform).
//!
//! ## Built-in devices
//!
//! | id | kind | notes |
//! |---|---|---|
//! | `cpu` | [`DeviceKind::Cpu`] | single thread, lane blocks the compiler vectorises; always available; the reference |
//! | `cpu-parallel` | [`DeviceKind::CpuParallel`] | same loops, rows split over the rayon pool; bit-identical to `cpu` |
//! | `gpu:wgpu:N` | [`DeviceKind::Gpu`] | feature `gpu-wgpu`; WGSL, `f32` |
//!
//! New devices: [`DeviceRegistry::register`]; see `docs/BACKENDS.md`.
//!
//! ## Unsafe code
//!
//! None in this module. The only `unsafe` a GPU device would need is buffer
//! byte-casting, which the `gpu-wgpu` implementation delegates to `bytemuck`.

#![allow(
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::missing_docs_in_private_items,
    clippy::module_name_repetitions,
    clippy::single_call_fn,
    clippy::too_many_arguments,
    clippy::needless_range_loop,
    reason = "numerical kernels index column-major buffers whose extents are \
              validated once per launch (see `check_launch`)"
)]

use std::any::Any;
use std::fmt;
use std::sync::Arc;
use std::sync::OnceLock;
use std::sync::RwLock;
use std::sync::atomic::AtomicU64;
use std::sync::atomic::Ordering;

use serde::Deserialize;
use serde::Serialize;

pub mod cpu;
#[cfg(feature = "gpu-wgpu")]
pub mod gpu;
pub mod kernel;
pub mod wgsl;

#[cfg(test)]
mod tests;

pub use cpu::CpuDevice;
pub use kernel::KernelBuilder;
pub use kernel::KernelIr;
pub use kernel::KernelOp;
pub use kernel::rk4_step;

use super::Backend;
use super::BackendError;
use super::Compiled;
use super::CompiledMulti;
use super::TapeFunction;
use super::tape;
use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::SymbolId;

/// Why a device operation failed.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum DeviceError {
    /// The device cannot do this (an operation, a data type, a host call).
    Unsupported(String),
    /// The arguments do not fit: wrong sizes, a buffer of another device or
    /// another data type.
    Mismatch(String),
    /// The kernel is malformed.
    InvalidKernel(String),
    /// A buffer could not be allocated.
    OutOfMemory(String),
    /// The underlying API reported an error.
    Backend(String),
}

impl fmt::Display for DeviceError {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        let (what, why) = match self {
            | Self::Unsupported(w) => ("unsupported on this device", w),
            | Self::Mismatch(w) => ("argument mismatch", w),
            | Self::InvalidKernel(w) => ("invalid kernel", w),
            | Self::OutOfMemory(w) => ("out of device memory", w),
            | Self::Backend(w) => ("device backend error", w),
        };
        write!(f, "{what}: {why}")
    }
}

impl std::error::Error for DeviceError {}

impl From<DeviceError> for BackendError {
    fn from(e: DeviceError) -> Self {
        Self::Codegen(e.to_string())
    }
}

/// Element type of a [`Buffer`] on the device. Host transfers always use
/// `f64`; the device converts.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Dtype {
    /// 32-bit float (what WGSL and most GPUs compute in).
    F32,
    /// 64-bit float (IEEE double; the semantics of the tape).
    F64,
}

/// Graphics/compute API of a GPU device.
#[derive(Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum GpuApi {
    /// Through the `wgpu` crate (Vulkan, Metal, DX12 or GL underneath).
    Wgpu,
    /// NVIDIA CUDA.
    Cuda,
    /// Apple Metal.
    Metal,
    /// Vulkan directly.
    Vulkan,
    /// OpenCL.
    OpenCl,
}

/// What a device is.
#[derive(Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum DeviceKind {
    /// One host thread.
    Cpu,
    /// The host's thread pool (rayon) over auto-vectorised loops.
    CpuParallel,
    /// A GPU through the given API.
    Gpu(GpuApi),
    /// Anything else (FPGA, remote worker, test double); the string names it.
    Custom(String),
}

impl DeviceKind {
    /// A short stable tag: `cpu`, `cpu-parallel`, `gpu-wgpu`, `gpu-cuda`...
    #[must_use]
    pub fn tag(&self) -> String {
        match self {
            | Self::Cpu => "cpu".into(),
            | Self::CpuParallel => "cpu-parallel".into(),
            | Self::Gpu(api) => {
                format!("gpu-{}", format!("{api:?}").to_ascii_lowercase())
            },
            | Self::Custom(name) => format!("custom-{name}"),
        }
    }
}

/// The operations a device implements. A launch of an unsupported one
/// returns [`DeviceError::Unsupported`]; nothing is emulated silently.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct DeviceOps {
    /// [`Device::launch`] of elementwise kernels.
    pub elementwise: bool,
    /// [`Device::reduce`].
    pub reduce: bool,
    /// [`Device::matmul`].
    pub matmul: bool,
}

/// Description of a device, fixed for its lifetime.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct DeviceInfo {
    /// Registry id (`cpu`, `cpu-parallel`, `gpu:wgpu:0`); unique.
    pub id: String,
    /// Human name (adapter name for a GPU).
    pub name: String,
    /// What it is.
    pub kind: DeviceKind,
    /// Data types its buffers and kernels support.
    pub dtypes: Vec<Dtype>,
    /// Operations it implements.
    pub ops: DeviceOps,
    /// Independent execution units (threads, or a GPU's reported maximum
    /// workgroup invocations; informational).
    pub parallelism: usize,
    /// Largest buffer, in elements of the narrowest supported type.
    pub max_buffer_elems: usize,
    /// Device memory is host memory (transfers are copies within RAM).
    pub unified_memory: bool,
    /// Results do not depend on thread count or scheduling.
    pub deterministic: bool,
}

/// Device-private buffer contents.
pub trait Storage: Any + Send + Sync {
    /// Upcast for the owning device to recover its concrete type.
    fn as_any(&self) -> &dyn Any;
    /// Mutable upcast.
    fn as_any_mut(&mut self) -> &mut dyn Any;
}

/// A typed array that lives on a device. Created only by that device.
/// Contents cross to the host only through [`Device::download`].
pub struct Buffer {
    owner: u64,
    dtype: Dtype,
    len: usize,
    storage: Box<dyn Storage>,
}

impl Buffer {
    /// Wraps device storage. For [`Device`] implementors.
    #[must_use]
    pub fn new(
        owner: u64,
        dtype: Dtype,
        len: usize,
        storage: Box<dyn Storage>,
    ) -> Self {
        Self {
            owner,
            dtype,
            len,
            storage,
        }
    }

    /// Number of elements.
    #[must_use]
    pub const fn len(&self) -> usize {
        self.len
    }

    /// Whether the buffer has no elements.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Element type.
    #[must_use]
    pub const fn dtype(&self) -> Dtype {
        self.dtype
    }

    /// The [`DeviceInfo`]-independent unique number of the owning device.
    #[must_use]
    pub const fn owner(&self) -> u64 {
        self.owner
    }

    /// The device-private storage.
    #[must_use]
    pub fn storage(&self) -> &dyn Storage {
        &*self.storage
    }

    /// The device-private storage, mutably.
    pub fn storage_mut(&mut self) -> &mut dyn Storage {
        &mut *self.storage
    }
}

impl fmt::Debug for Buffer {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        f.debug_struct("Buffer")
            .field("owner", &self.owner)
            .field("dtype", &self.dtype)
            .field("len", &self.len)
            .finish_non_exhaustive()
    }
}

/// A kernel compiled by a device.
pub trait CompiledKernel: Send + Sync {
    /// The program it was compiled from.
    fn ir(&self) -> &KernelIr;
    /// Upcast for the owning device.
    fn as_any(&self) -> &dyn Any;
}

/// A reduction.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ReduceOp {
    /// Sum.
    Sum,
    /// Product.
    Prod,
    /// Minimum (NaN if any element is NaN).
    Min,
    /// Maximum (NaN if any element is NaN).
    Max,
}

impl ReduceOp {
    /// The neutral element.
    #[must_use]
    pub const fn identity(self) -> f64 {
        match self {
            | Self::Sum => 0.0,
            | Self::Prod => 1.0,
            | Self::Min => f64::INFINITY,
            | Self::Max => f64::NEG_INFINITY,
        }
    }

    /// `acc ∘ x`.
    #[inline]
    #[must_use]
    pub fn combine(
        self,
        acc: f64,
        x: f64,
    ) -> f64 {
        match self {
            | Self::Sum => acc + x,
            | Self::Prod => acc * x,
            | Self::Min => {
                if acc.is_nan() || x.is_nan() {
                    f64::NAN
                } else {
                    acc.min(x)
                }
            },
            | Self::Max => {
                if acc.is_nan() || x.is_nan() {
                    f64::NAN
                } else {
                    acc.max(x)
                }
            },
        }
    }
}

static NEXT_UID: AtomicU64 = AtomicU64::new(1);

/// A fresh owner number for a [`Buffer`]'s device. Each device instance takes
/// one at construction.
#[must_use]
pub fn next_uid() -> u64 {
    NEXT_UID.fetch_add(1, Ordering::Relaxed)
}

/// An execution device.
///
/// Implementations must be usable from several threads. Operations may be
/// asynchronous; [`Device::download`] and [`Device::synchronize`] are the
/// points at which their effects are visible.
pub trait Device: Send + Sync {
    /// Static description.
    fn info(&self) -> &DeviceInfo;

    /// The number identifying this instance as the owner of its buffers.
    fn uid(&self) -> u64;

    /// A zero-filled buffer of `len` elements.
    ///
    /// # Errors
    /// [`DeviceError::Unsupported`] for a data type the device lacks,
    /// [`DeviceError::OutOfMemory`].
    fn alloc(
        &self,
        len: usize,
        dtype: Dtype,
    ) -> Result<Buffer, DeviceError>;

    /// Host to device transfer (converting to `dtype`).
    ///
    /// # Errors
    /// As [`Device::alloc`].
    fn upload(
        &self,
        host: &[f64],
        dtype: Dtype,
    ) -> Result<Buffer, DeviceError>;

    /// Device to host transfer; `out.len()` must equal `buffer.len()`.
    ///
    /// # Errors
    /// [`DeviceError::Mismatch`] on a foreign buffer or wrong length.
    fn download(
        &self,
        buffer: &Buffer,
        out: &mut [f64],
    ) -> Result<(), DeviceError>;

    /// Compiles a kernel. Implementations cache by
    /// [`KernelIr::structural_hash`].
    ///
    /// # Errors
    /// [`DeviceError::InvalidKernel`], [`DeviceError::Unsupported`].
    fn compile(
        &self,
        ir: &KernelIr,
    ) -> Result<Arc<dyn CompiledKernel>, DeviceError>;

    /// Runs `kernel` for `rows` rows. `input` holds `ir.inputs` columns of
    /// `rows` elements (column-major), `output` `ir.outputs.len()` columns.
    ///
    /// # Errors
    /// [`DeviceError::Mismatch`] on foreign buffers or wrong sizes.
    fn launch(
        &self,
        kernel: &dyn CompiledKernel,
        rows: usize,
        input: &Buffer,
        params: &[f64],
        output: &mut Buffer,
    ) -> Result<(), DeviceError>;

    /// Reduces a whole buffer to one number.
    ///
    /// # Errors
    /// [`DeviceError::Unsupported`] unless `info().ops.reduce`.
    fn reduce(
        &self,
        _op: ReduceOp,
        _input: &Buffer,
    ) -> Result<f64, DeviceError> {
        Err(DeviceError::Unsupported("reduce".into()))
    }

    /// `c = a * b` for row-major `a` (`m x k`), `b` (`k x n`), `c` (`m x n`).
    ///
    /// # Errors
    /// [`DeviceError::Unsupported`] unless `info().ops.matmul`.
    fn matmul(
        &self,
        _a: &Buffer,
        _b: &Buffer,
        _c: &mut Buffer,
        _m: usize,
        _k: usize,
        _n: usize,
    ) -> Result<(), DeviceError> {
        Err(DeviceError::Unsupported("matmul".into()))
    }

    /// Waits for submitted work.
    ///
    /// # Errors
    /// [`DeviceError::Backend`].
    fn synchronize(&self) -> Result<(), DeviceError> {
        Ok(())
    }
}

/// Checks the buffers and sizes of a launch. For implementors.
///
/// # Errors
/// [`DeviceError::Mismatch`].
pub fn check_launch(
    uid: u64,
    ir: &KernelIr,
    rows: usize,
    input: &Buffer,
    output: &Buffer,
) -> Result<(), DeviceError> {
    let bad = |m: String| Err(DeviceError::Mismatch(m));
    if input.owner != uid || output.owner != uid {
        return bad("a buffer belongs to another device; transfer it explicitly".into());
    }
    let want_in = ir.inputs.checked_mul(rows);
    let want_out = ir.outputs.len().checked_mul(rows);
    if want_in != Some(input.len) {
        return bad(format!(
            "input has {} elements, {} columns of {rows} rows need {:?}",
            input.len, ir.inputs, want_in
        ));
    }
    if want_out != Some(output.len) {
        return bad(format!(
            "output has {} elements, {} columns of {rows} rows need {:?}",
            output.len,
            ir.outputs.len(),
            want_out
        ));
    }
    Ok(())
}

/// The registered devices.
pub struct DeviceRegistry {
    devices: RwLock<Vec<Arc<dyn Device>>>,
}

impl Default for DeviceRegistry {
    fn default() -> Self {
        Self::new()
    }
}

impl DeviceRegistry {
    /// A registry holding the built-in CPU devices (and, with feature
    /// `gpu-wgpu`, every adapter wgpu finds; none is not an error).
    #[must_use]
    pub fn new() -> Self {
        let r = Self::empty();
        let _ = r.register(Arc::new(CpuDevice::sequential()));
        let _ = r.register(Arc::new(CpuDevice::parallel()));
        #[cfg(feature = "gpu-wgpu")]
        for device in gpu::probe() {
            let _ = r.register(device);
        }
        r
    }

    /// A registry with no devices.
    #[must_use]
    pub const fn empty() -> Self {
        Self {
            devices: RwLock::new(Vec::new()),
        }
    }

    /// Adds a device.
    ///
    /// # Errors
    /// [`DeviceError::Mismatch`] when the id is taken.
    pub fn register(
        &self,
        device: Arc<dyn Device>,
    ) -> Result<(), DeviceError> {
        let mut guard = self
            .devices
            .write()
            .map_err(|_| DeviceError::Backend("registry poisoned".into()))?;
        if guard.iter().any(|d| d.info().id == device.info().id) {
            return Err(DeviceError::Mismatch(format!(
                "device `{}` is already registered",
                device.info().id
            )));
        }
        guard.push(device);
        Ok(())
    }

    /// Descriptions of every device, in registration order.
    #[must_use]
    pub fn list(&self) -> Vec<DeviceInfo> {
        self.devices
            .read()
            .map(|g| g.iter().map(|d| d.info().clone()).collect())
            .unwrap_or_default()
    }

    /// The device with this id.
    #[must_use]
    pub fn get(
        &self,
        id: &str,
    ) -> Option<Arc<dyn Device>> {
        self.devices
            .read()
            .ok()?
            .iter()
            .find(|d| d.info().id == id)
            .cloned()
    }

    /// The first device of a kind.
    #[must_use]
    pub fn first_of(
        &self,
        kind: &DeviceKind,
    ) -> Option<Arc<dyn Device>> {
        self.devices
            .read()
            .ok()?
            .iter()
            .find(|d| &d.info().kind == kind)
            .cloned()
    }

    /// The device to use when the caller names none: the first GPU, else the
    /// parallel CPU, else any.
    #[must_use]
    pub fn default_device(&self) -> Option<Arc<dyn Device>> {
        let guard = self.devices.read().ok()?;
        guard
            .iter()
            .find(|d| matches!(d.info().kind, DeviceKind::Gpu(_)))
            .or_else(|| {
                guard
                    .iter()
                    .find(|d| d.info().kind == DeviceKind::CpuParallel)
            })
            .or_else(|| guard.first())
            .cloned()
    }
}

/// The process-wide registry (created on first use).
#[must_use]
pub fn registry() -> &'static DeviceRegistry {
    static REGISTRY: OnceLock<DeviceRegistry> = OnceLock::new();
    REGISTRY.get_or_init(DeviceRegistry::new)
}

/// The device kinds this *build* can offer, as [`DeviceKind::tag`]s: a
/// compile-time fact (it does not probe hardware), so it can go in the kernel
/// manifest. `gpu-wgpu` here means the build can drive GPUs; whether one is
/// present is answered by [`registry`]`().list()`.
#[must_use]
pub fn build_devices() -> Vec<String> {
    let mut v = vec![DeviceKind::Cpu.tag(), DeviceKind::CpuParallel.tag()];
    if cfg!(feature = "gpu-wgpu") {
        v.push(DeviceKind::Gpu(GpuApi::Wgpu).tag());
    }
    v
}

// --- host-side conveniences built on the trait --------------------------------

/// Packs columns into one column-major host vector.
///
/// # Errors
/// [`DeviceError::Mismatch`] when the columns differ in length.
pub fn pack_columns(columns: &[&[f64]]) -> Result<(Vec<f64>, usize), DeviceError> {
    let rows = columns.first().map_or(0, |c| c.len());
    if columns.iter().any(|c| c.len() != rows) {
        return Err(DeviceError::Mismatch("columns differ in length".into()));
    }
    let mut flat = Vec::with_capacity(rows * columns.len());
    for c in columns {
        flat.extend_from_slice(c);
    }
    Ok((flat, rows))
}

/// Upload, launch, download: `columns[i][row]` in, one vector per output out.
/// The simple path; hold [`Buffer`]s yourself to keep data on the device
/// between launches.
///
/// # Errors
/// Any [`DeviceError`].
pub fn map_columns(
    device: &dyn Device,
    ir: &KernelIr,
    columns: &[&[f64]],
    params: &[f64],
) -> Result<Vec<Vec<f64>>, DeviceError> {
    if columns.len() != ir.inputs {
        return Err(DeviceError::Mismatch(format!(
            "kernel has {} inputs, {} columns given",
            ir.inputs,
            columns.len()
        )));
    }
    let (flat, rows) = if columns.is_empty() {
        (Vec::new(), 0)
    } else {
        pack_columns(columns)?
    };
    let dtype = preferred_dtype(device.info());
    let kernel = device.compile(ir)?;
    let input = device.upload(&flat, dtype)?;
    let mut output = device.alloc(ir.outputs.len() * rows, dtype)?;
    device.launch(&*kernel, rows, &input, params, &mut output)?;
    let mut host = vec![0.0; output.len()];
    device.download(&output, &mut host)?;
    Ok(if rows == 0 {
        vec![Vec::new(); ir.outputs.len()]
    } else {
        host.chunks(rows).map(<[f64]>::to_vec).collect()
    })
}

/// `F64` when the device has it, else `F32`.
#[must_use]
pub fn preferred_dtype(info: &DeviceInfo) -> Dtype {
    if info.dtypes.contains(&Dtype::F64) {
        Dtype::F64
    } else {
        Dtype::F32
    }
}

/// Integrates `y' = f([t,] y; params)` for many initial conditions at once
/// with fixed-step RK4, entirely through [`Device::launch`].
///
/// `y0[i][row]` is state variable `i` of trajectory `row`. Returns the states
/// at `t1` in the same layout. The step kernel ([`rk4_step`]) is compiled once;
/// each of the `steps` launches advances every trajectory by
/// `(t1 - t0) / steps`. Trajectories are independent, so the result for a row
/// does not depend on the rest of the batch.
///
/// # Errors
/// Any [`DeviceError`]; [`DeviceError::Mismatch`] for `steps == 0`.
pub fn ensemble_rk4(
    device: &dyn Device,
    rhs: &KernelIr,
    params: &[f64],
    y0: &[&[f64]],
    t0: f64,
    t1: f64,
    steps: usize,
) -> Result<Vec<Vec<f64>>, DeviceError> {
    if steps == 0 {
        return Err(DeviceError::Mismatch("steps must be positive".into()));
    }
    if params.len() != rhs.params {
        return Err(DeviceError::Mismatch(format!(
            "rhs has {} parameters, {} given",
            rhs.params,
            params.len()
        )));
    }
    let n = y0.len();
    let step = rk4_step(rhs, n)?;
    let (flat, rows) = pack_columns(y0)?;
    let dtype = preferred_dtype(device.info());
    let kernel = device.compile(&step)?;
    let mut a = device.upload(&flat, dtype)?;
    let mut b = device.alloc(n * rows, dtype)?;
    let h = (t1 - t0) / steps as f64;
    let mut p = params.to_vec();
    p.extend([t0, h]);
    let (t_at, h_at) = (params.len(), params.len() + 1);
    for s in 0..steps {
        // t is recomputed from the step number, not accumulated, so the
        // schedule does not drift with the number of steps.
        p[t_at] = t0 + h * s as f64;
        p[h_at] = h;
        device.launch(&*kernel, rows, &a, &p, &mut b)?;
        std::mem::swap(&mut a, &mut b);
    }
    let mut host = vec![0.0; n * rows];
    device.download(&a, &mut host)?;
    Ok(if rows == 0 {
        vec![Vec::new(); n]
    } else {
        host.chunks(rows).map(<[f64]>::to_vec).collect()
    })
}

// --- Backend integration --------------------------------------------------------

/// A [`Backend`] whose compiled functions run batches on a [`Device`].
///
/// Single-point calls still work (a one-row launch) but are slow on a GPU;
/// the benefit is in `call_batch` / `eval_batch`. Programs the device cannot
/// run (operator calls, extern functions: [`KernelIr::from_tape`]) fall back
/// to the host tape interpreter, so compilation never fails because of the
/// device.
#[derive(Clone)]
pub struct DeviceBackend {
    device: Arc<dyn Device>,
}

impl DeviceBackend {
    /// A backend executing on `device`.
    #[must_use]
    pub fn new(device: Arc<dyn Device>) -> Self {
        Self { device }
    }
}

impl fmt::Debug for DeviceBackend {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        write!(f, "DeviceBackend({})", self.device.info().id)
    }
}

struct DeviceProgram {
    device: Arc<dyn Device>,
    kernel: Arc<dyn CompiledKernel>,
}

impl DeviceProgram {
    fn run(
        &self,
        columns: &[&[f64]],
        params: &[f64],
        rows: usize,
    ) -> Result<Vec<Vec<f64>>, DeviceError> {
        let ir = self.kernel.ir();
        let dtype = preferred_dtype(self.device.info());
        let mut flat = Vec::with_capacity(rows * ir.inputs);
        for i in 0..ir.inputs {
            let column = columns.get(i).copied().unwrap_or(&[]);
            for row in 0..rows {
                flat.push(column.get(row).copied().unwrap_or(f64::NAN));
            }
        }
        let input = self.device.upload(&flat, dtype)?;
        let mut output = self.device.alloc(ir.outputs.len() * rows, dtype)?;
        self.device
            .launch(&*self.kernel, rows, &input, params, &mut output)?;
        let mut host = vec![0.0; output.len()];
        self.device.download(&output, &mut host)?;
        Ok(if rows == 0 {
            vec![Vec::new(); ir.outputs.len()]
        } else {
            host.chunks(rows).map(<[f64]>::to_vec).collect()
        })
    }
}

impl CompiledMulti for DeviceProgram {
    fn inputs(&self) -> usize {
        self.kernel.ir().inputs
    }

    fn params(&self) -> usize {
        self.kernel.ir().params
    }

    fn outputs(&self) -> usize {
        self.kernel.ir().outputs.len()
    }

    fn eval(
        &self,
        inputs: &[f64],
        params: &[f64],
        out: &mut [f64],
    ) {
        let cols: Vec<&[f64]> = inputs.iter().map(std::slice::from_ref).collect();
        let result = self.run(&cols, params, 1);
        for (k, slot) in out.iter_mut().enumerate() {
            *slot = result
                .as_ref()
                .ok()
                .and_then(|r| r.get(k))
                .and_then(|c| c.first())
                .copied()
                .unwrap_or(f64::NAN);
        }
    }

    fn eval_batch(
        &self,
        columns: &[&[f64]],
        params: &[f64],
        out: &mut [&mut [f64]],
    ) {
        let rows = out.iter().map(|o| o.len()).min().unwrap_or(0);
        match self.run(columns, params, rows) {
            | Ok(result) => {
                for (dst, src) in out.iter_mut().zip(&result) {
                    for (d, s) in dst.iter_mut().zip(src) {
                        *d = *s;
                    }
                }
            },
            | Err(_) => {
                for dst in out.iter_mut() {
                    dst.fill(f64::NAN);
                }
            },
        }
    }
}

struct SingleOutput(Box<dyn CompiledMulti>);

impl Compiled for SingleOutput {
    fn arity(&self) -> usize {
        self.0.inputs()
    }

    fn call(
        &self,
        args: &[f64],
    ) -> f64 {
        let mut out = [f64::NAN];
        self.0.eval(args, &[], &mut out);
        out[0]
    }

    fn call_batch(
        &self,
        columns: &[&[f64]],
        out: &mut [f64],
    ) {
        self.0.eval_batch(columns, &[], &mut [out]);
    }
}

impl Backend for DeviceBackend {
    fn name(&self) -> &'static str {
        "device"
    }

    fn compile(
        &self,
        graph: &Graph,
        root: NodeId,
        inputs: &[SymbolId],
    ) -> Result<Box<dyn Compiled>, BackendError> {
        let multi = self.compile_multi(graph, &[root], inputs, &[])?;
        Ok(Box::new(SingleOutput(multi)))
    }

    fn compile_multi(
        &self,
        graph: &Graph,
        roots: &[NodeId],
        inputs: &[SymbolId],
        params: &[SymbolId],
    ) -> Result<Box<dyn CompiledMulti>, BackendError> {
        let t = tape::lower(graph, roots, inputs, params)?.optimise();
        let program = KernelIr::from_tape(&t)
            .and_then(|ir| self.device.compile(&ir))
            .map(|kernel| {
                DeviceProgram {
                    device: self.device.clone(),
                    kernel,
                }
            });
        let boxed: Box<dyn CompiledMulti> = match program {
            | Ok(p) => Box::new(p),
            | Err(_) => Box::new(TapeFunction::from_tape(t)),
        };
        Ok(boxed)
    }
}
