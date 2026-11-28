//! # Execution backends
//!
//! Once a request has been reduced to a closed form, evaluating it many
//! times — over the rows of a data set, inside an ODE stepper, at every
//! cell of a simulation grid — is a different problem from finding it. A
//! [`Backend`] turns a closed-form term into a [`Compiled`] function of its
//! free symbols.
//!
//! The contract is deliberately small so that very different executors fit
//! behind it:
//!
//! * [`Interpreter`] (this module) flattens the term's DAG into a tape and
//!   walks it; always available, no code generation;
//! * a JIT backend lowers the same DAG to machine code. Each operator's
//!   [`EvalFn`] is a plain `fn` pointer it can
//!   emit a direct call to, and richer lowerings can be attached to
//!   operators as typed attributes
//!   ([`OpTable::set_attr`](crate::graph::OpTable::set_attr)) without the
//!   kernel knowing about them.
//!
//! * [`device`] is the data-parallel axis: it runs one program over many rows
//!   on a CPU thread pool or a GPU, behind a capability-described `Device`
//!   trait (see `docs/BACKENDS.md`).
//!
//! Numeric reduction kernels and the simulation steppers only see
//! `dyn Compiled`, so switching executors never touches them.

use std::cell::RefCell;
use std::fmt;
use std::sync::Arc;
use std::sync::OnceLock;
use std::sync::RwLock;

pub(crate) mod access;
pub use access::SourceAccess;
pub use access::gpu_allowed;
pub use access::restrict_to_source;
/// Checked forward/reverse differentiation of elementary expressions.
pub mod autodiff;
pub mod device;
#[cfg(feature = "jit")]
pub mod jit;
pub mod lir;
#[cfg(feature = "jit")]
pub mod lir_jit;
pub mod tape;

use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::SymbolId;
use crate::graph::op::EvalFn;

/// Why a term could not be compiled.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum BackendError {
    /// The term depends on a symbol that is not among the inputs.
    UnboundSymbol(String),
    /// The term uses an operator without scalar semantics (for example an
    /// unevaluated heavy operator).
    NoSemantics(String),
    /// The term contains a non-numeric literal.
    NotNumeric,
    /// The code generator failed.
    Codegen(String),
}

impl fmt::Display for BackendError {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        match self {
            | Self::UnboundSymbol(name) => {
                write!(
                    f,
                    "symbol `{name}` is not an input of the compiled function"
                )
            },
            | Self::NoSemantics(op) => write!(f, "operator `{op}` has no numeric semantics"),
            | Self::NotNumeric => write!(f, "the term contains a non-numeric literal"),
            | Self::Codegen(reason) => write!(f, "code generation failed: {reason}"),
        }
    }
}

impl std::error::Error for BackendError {}

/// A term compiled to a function `f64^n -> f64`.
pub trait Compiled: Send + Sync {
    /// The executable's retained allocation, including storage it owns.
    ///
    /// Backends without a bound return `None`. A detached program needs a
    /// bound before another native object can reserve storage for it.
    fn retained_bytes(&self) -> Option<usize> {
        None
    }

    /// Number of inputs.
    fn arity(&self) -> usize;

    /// Evaluates at one point. `args` must hold [`Compiled::arity`] values;
    /// missing ones are read as NaN.
    fn call(
        &self,
        args: &[f64],
    ) -> f64;

    /// Evaluates at many points: `columns[i][row]` is input `i` of `row`.
    /// Writes one result per element of `out`.
    fn call_batch(
        &self,
        columns: &[&[f64]],
        out: &mut [f64],
    ) {
        let arity = self.arity();
        let mut args = crate::kernels::reuse::f64_scratch(arity);
        args.resize(arity, f64::NAN);
        for (row, slot) in out.iter_mut().enumerate() {
            for (arg, column) in args.iter_mut().zip(columns) {
                *arg = column.get(row).copied().unwrap_or(f64::NAN);
            }
            *slot = self.call(&args);
        }
    }
}

/// Rows from which [`call_batch_par`] splits a batch over the thread pool.
const BATCH_PAR_ROWS: usize = 4096;

/// [`Compiled::call_batch`] split into fixed blocks of `BATCH_PAR_ROWS` rows.
///
/// The blocks are
/// evaluated concurrently. Every row is an independent evaluation, so the output is identical to
/// `call_batch`; batches of at most two blocks take the sequential call.
pub fn call_batch_par(
    f: &dyn Compiled,
    columns: &[&[f64]],
    out: &mut [f64],
) {
    if !crate::execution::should_split(out.len().div_ceil(BATCH_PAR_ROWS), 3) {
        f.call_batch(columns, out);
        return;
    }
    crate::kernels::par::for_chunks_mut(out, BATCH_PAR_ROWS, 3, |k, chunk| {
        call_batch_chunk(f, columns, k * BATCH_PAR_ROWS, chunk);
    });
}

fn call_batch_chunk(
    f: &dyn Compiled,
    columns: &[&[f64]],
    start: usize,
    out: &mut [f64],
) {
    let cols = batch_columns(columns, start, out.len());
    f.call_batch(&cols, out);
}

fn batch_columns<'a>(
    columns: &[&'a [f64]],
    start: usize,
    rows: usize,
) -> Vec<&'a [f64]> {
    let end = start + rows;
    columns
        .iter()
        .map(|c| &c[start.min(c.len())..end.min(c.len())])
        .collect()
}

/// Run a checked batch in fixed blocks sized by estimated executable work.
pub(crate) fn call_batch_checked(
    f: &dyn Compiled,
    columns: &[&[f64]],
    out: &mut [f64],
    work_per_row: usize,
) -> Result<(), crate::cancellation::Cancelled> {
    crate::cancellation::checkpoint()?;
    let rows_per_block = (65_536 / work_per_row.max(1)).clamp(256, BATCH_PAR_ROWS);
    if out.len() <= rows_per_block {
        f.call_batch(columns, out);
    } else {
        crate::kernels::par::for_chunks_mut(out, rows_per_block, 2, |k, chunk| {
            if !crate::cancellation::is_cancelled() {
                call_batch_chunk(f, columns, k * rows_per_block, chunk);
            }
        });
    }
    crate::cancellation::checkpoint()
}

/// A compiled program with several inputs, parameters and outputs (an ODE
/// right-hand side, a gradient, a set of expressions sharing
/// subexpressions).
pub trait CompiledMulti: Send + Sync {
    /// Number of inputs.
    fn inputs(&self) -> usize;
    /// Number of parameters (constant during a batch).
    fn params(&self) -> usize;
    /// Number of outputs.
    fn outputs(&self) -> usize;
    /// Evaluates every output at one point.
    fn eval(
        &self,
        inputs: &[f64],
        params: &[f64],
        out: &mut [f64],
    );
    /// Evaluates at many points: `columns[i][row]` is input `i` of `row`,
    /// `out[k][row]` output `k`.
    fn eval_batch(
        &self,
        columns: &[&[f64]],
        params: &[f64],
        out: &mut [&mut [f64]],
    ) {
        let rows = out.iter().map(|o| o.len()).min().unwrap_or(0);
        let (inputs, outputs) = (self.inputs(), self.outputs());
        let mut point = crate::kernels::reuse::f64_scratch(inputs);
        point.resize(inputs, f64::NAN);
        let mut values = crate::kernels::reuse::f64_scratch(outputs);
        values.resize(outputs, f64::NAN);
        for row in 0..rows {
            for (arg, column) in point.iter_mut().zip(columns) {
                *arg = column.get(row).copied().unwrap_or(f64::NAN);
            }
            self.eval(&point, params, &mut values);
            for (column, &v) in out.iter_mut().zip(values.iter()) {
                if let Some(slot) = column.get_mut(row) {
                    *slot = v;
                }
            }
        }
    }
}

/// Elementary functions a code generator implements natively or through a
/// direct call to the platform's math library.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Intrinsic {
    /// `exp`.
    Exp,
    /// `ln`.
    Ln,
    /// `sin`.
    Sin,
    /// `cos`.
    Cos,
    /// `tan`.
    Tan,
    /// `asin`.
    Asin,
    /// `acos`.
    Acos,
    /// `atan`.
    Atan,
    /// `sinh`.
    Sinh,
    /// `cosh`.
    Cosh,
    /// `tanh`.
    Tanh,
    /// `asinh`.
    Asinh,
    /// `acosh`.
    Acosh,
    /// `atanh`.
    Atanh,
    /// `sqrt`.
    Sqrt,
    /// `abs`.
    Abs,
    /// `floor`.
    Floor,
    /// `ceil`.
    Ceil,
    /// `cbrt`.
    Cbrt,
}

/// How a compiling backend should lower an operator.
///
/// It is attached with
/// [`OpTable::set_attr`](crate::graph::OpTable::set_attr). Operators
/// without it are matched by name and otherwise called through their
/// [`EvalFn`].
#[derive(Clone, Copy, Debug)]
pub enum Lowering {
    /// A known elementary function.
    Intrinsic(Intrinsic),
    /// A scalar C function of one argument.
    Extern1(extern "C" fn(f64) -> f64),
    /// A scalar C function of two arguments.
    Extern2(extern "C" fn(f64, f64) -> f64),
    /// `gather(tensor, offset)`: read `inputs[base + offset]`, where `base`
    /// is the position of the tensor's first element in the caller's
    /// flattened input buffer.
    ///
    /// A tensor has no tape representation of its own — the tape is scalar —
    /// so indexing is an offset into the input buffer. An index outside the
    /// tensor reads NaN, never memory outside the buffer.
    Gather,
}

impl PartialEq for Lowering {
    /// Compares *kinds*, not function addresses.
    ///
    /// Two function pointers are not reliably comparable — the compiler may
    /// merge or duplicate a function — so comparing them would give an answer
    /// that depends on the build. Callers here ask "is this the same kind of
    /// lowering", which the discriminant answers.
    fn eq(
        &self,
        other: &Self,
    ) -> bool {
        match (self, other) {
            | (Self::Intrinsic(a), Self::Intrinsic(b)) => a == b,
            | (Self::Extern1(_), Self::Extern1(_))
            | (Self::Extern2(_), Self::Extern2(_))
            | (Self::Gather, Self::Gather) => true,
            | _ => false,
        }
    }
}

impl Eq for Lowering {}

/// A [`tape::Tape`] evaluated by the reference tape interpreter.
#[derive(Clone, Debug)]
pub struct TapeFunction {
    tape: tape::Tape,
}

impl TapeFunction {
    /// Wraps a lowered program (the host fallback of [`device::DeviceBackend`]).
    #[must_use]
    pub const fn from_tape(tape: tape::Tape) -> Self {
        Self { tape }
    }

    /// The underlying program.
    #[must_use]
    pub const fn tape(&self) -> &tape::Tape {
        &self.tape
    }
}

impl CompiledMulti for TapeFunction {
    fn inputs(&self) -> usize {
        self.tape.inputs
    }

    fn params(&self) -> usize {
        self.tape.params
    }

    fn outputs(&self) -> usize {
        self.tape.outputs.len()
    }

    fn eval(
        &self,
        inputs: &[f64],
        params: &[f64],
        out: &mut [f64],
    ) {
        self.tape.eval(inputs, params, out);
    }

    fn eval_batch(
        &self,
        columns: &[&[f64]],
        params: &[f64],
        out: &mut [&mut [f64]],
    ) {
        self.tape.eval_columns(columns, params, out);
    }
}

impl Compiled for TapeFunction {
    fn retained_bytes(&self) -> Option<usize> {
        let instructions = self
            .tape
            .insts
            .capacity()
            .checked_mul(size_of::<tape::Inst>())?;
        let arguments = self
            .tape
            .args
            .capacity()
            .checked_mul(size_of::<tape::Reg>())?;
        let outputs = self
            .tape
            .outputs
            .capacity()
            .checked_mul(size_of::<tape::Reg>())?;
        size_of::<Self>()
            .checked_add(instructions)?
            .checked_add(arguments)?
            .checked_add(outputs)
    }

    fn arity(&self) -> usize {
        self.tape.inputs
    }

    fn call(
        &self,
        args: &[f64],
    ) -> f64 {
        let mut out = [f64::NAN];
        self.tape.eval(args, &[], &mut out);
        out[0]
    }

    fn call_batch(
        &self,
        columns: &[&[f64]],
        out: &mut [f64],
    ) {
        self.tape.eval_columns(columns, &[], &mut [out]);
    }
}

/// Turns closed-form terms into executable functions.
pub trait Backend: Send + Sync {
    /// The backend instance's retained allocation, when bounded.
    ///
    /// This excludes allocations retained individually by compiled functions.
    fn retained_bytes(&self) -> Option<usize> {
        None
    }

    /// Name for diagnostics.
    fn name(&self) -> &'static str;

    /// Canonical selector for this backend, including a concrete device ID.
    fn selector(&self) -> &str {
        self.name()
    }

    /// Whether this backend may execute on a GPU.
    ///
    /// Unclassified host backends are treated conservatively. Implementations
    /// that execute solely on CPUs can return `false`.
    fn uses_gpu(&self) -> bool {
        true
    }

    /// Whether session-owned compilation may reserve memory for this executor.
    /// Unbounded custom backends are kept out of the portable value ABI.
    #[doc(hidden)]
    fn native_program_budget(&self) -> Option<NativeProgramBudget> {
        None
    }

    /// Compiles a scalar subject to the limits of [`Self::native_program_budget`].
    #[doc(hidden)]
    fn compile_native(
        &self,
        graph: &Graph,
        root: NodeId,
        inputs: &[SymbolId],
    ) -> Result<Box<dyn Compiled>, BackendError> {
        self.compile(graph, root, inputs)
    }

    /// Compiles the concrete term `root` as a function of `inputs`, in that
    /// order.
    ///
    /// # Errors
    /// See [`BackendError`].
    fn compile(
        &self,
        graph: &Graph,
        root: NodeId,
        inputs: &[SymbolId],
    ) -> Result<Box<dyn Compiled>, BackendError>;

    /// Compiles several terms into one program with shared
    /// subexpressions, as functions of `inputs` and `params`.
    ///
    /// The default lowers to the [Tape IR](tape) and interprets it.
    ///
    /// # Errors
    /// See [`BackendError`].
    fn compile_multi(
        &self,
        graph: &Graph,
        roots: &[NodeId],
        inputs: &[SymbolId],
        params: &[SymbolId],
    ) -> Result<Box<dyn CompiledMulti>, BackendError> {
        let tape = tape::lower(graph, roots, inputs, params)?.optimise();
        Ok(Box::new(TapeFunction { tape }))
    }
}

/// Limits for compiled expressions retained by the portable value ABI.
#[doc(hidden)]
#[derive(Clone, Copy)]
pub struct NativeProgramBudget {
    /// Maximum number of lowered instructions.
    pub instructions: usize,
    /// Reservation for the executor beyond the retained graph and tape.
    pub extra_bytes: usize,
}

pub(crate) const NATIVE_TAPE_BUDGET: NativeProgramBudget = NativeProgramBudget {
    instructions: 4096,
    extra_bytes: 0,
};
pub(crate) const NATIVE_CODE_BUDGET: NativeProgramBudget = NativeProgramBudget {
    instructions: 512,
    extra_bytes: 4_194_304,
};

type Shared = Arc<dyn Backend>;

static DEFAULT_BACKEND: RwLock<Option<Shared>> = RwLock::new(None);

thread_local! {
    static SCOPED_BACKEND: RefCell<Option<Shared>> = const { RefCell::new(None) };
}

fn interpreter() -> Shared {
    static BACKEND: OnceLock<Shared> = OnceLock::new();
    BACKEND.get_or_init(|| Arc::new(Interpreter)).clone()
}

#[cfg(feature = "jit")]
fn jit_available() -> bool {
    static AVAILABLE: OnceLock<bool> = OnceLock::new();
    *AVAILABLE.get_or_init(|| cranelift_native::builder().is_ok())
}

fn automatic() -> Shared {
    static AUTOMATIC: OnceLock<Shared> = OnceLock::new();
    AUTOMATIC
        .get_or_init(|| {
            #[cfg(feature = "jit")]
            if jit_available() {
                return Arc::new(jit::TieredBackend::default());
            }
            interpreter()
        })
        .clone()
}

/// Returns the backend numeric kernels compile with.
///
/// A thread scope takes priority over a process override, then the automatic
/// choice. On supported JIT hosts this starts with the CPU tape and compiles
/// hot programs in the background.
#[must_use]
pub fn current() -> Shared {
    let allow_gpu = access::gpu_allowed();
    let allowed = |backend: &&Shared| allow_gpu || !backend.uses_gpu();
    SCOPED_BACKEND
        .with(|scope| scope.borrow().as_ref().filter(allowed).cloned())
        .or_else(|| {
            DEFAULT_BACKEND
                .read()
                .ok()
                .and_then(|slot| slot.as_ref().filter(allowed).cloned())
        })
        .unwrap_or_else(automatic)
}

fn cpu_backend(parallel: bool) -> Shared {
    static SEQUENTIAL: OnceLock<Shared> = OnceLock::new();
    static PARALLEL: OnceLock<Shared> = OnceLock::new();
    let slot = if parallel {
        &PARALLEL
    } else {
        &SEQUENTIAL
    };
    slot.get_or_init(|| {
        let device = if parallel {
            device::CpuDevice::parallel()
        } else {
            device::CpuDevice::sequential()
        };
        Arc::new(device::DeviceBackend::requiring_device(Arc::new(device)))
    })
    .clone()
}

/// Lists selectors supported by this build and its discovered devices.
/// The first call may initialize the GPU device registry.
#[must_use]
pub fn available() -> Vec<String> {
    let mut names = vec!["auto".into(), "interpreter".into()];
    #[cfg(feature = "jit")]
    if jit_available() {
        names.extend(["tiered".into(), "cranelift".into()]);
    }
    names.extend(
        device::registry_for_source()
            .list()
            .into_iter()
            .map(|info| info.id)
            .filter(|id| {
                !matches!(
                    id.as_str(),
                    "auto" | "interpreter" | "tiered" | "cranelift" | "jit" | "gpu"
                )
            }),
    );
    names
}

/// Resolves `auto`, `interpreter`, `tiered`, `cranelift`, a CPU selector,
/// or a device ID from [`available`]. `jit` is an alias for `cranelift`.
/// `gpu` picks the first discovered `wgpu` adapter.
///
/// # Errors
/// Returns [`BackendError::Codegen`] when the requested backend is unavailable.
pub fn named(selector: &str) -> Result<Shared, BackendError> {
    match selector {
        | "auto" => return Ok(automatic()),
        | "interpreter" => return Ok(interpreter()),
        | "cpu" => return Ok(cpu_backend(false)),
        | "cpu-parallel" => return Ok(cpu_backend(true)),
        #[cfg(feature = "jit")]
        | "tiered" if jit_available() => return Ok(Arc::new(jit::TieredBackend::default())),
        #[cfg(feature = "jit")]
        | "cranelift" | "jit" if jit_available() => {
            return Ok(Arc::new(jit::CraneliftBackend::uncached()));
        },
        | _ => {},
    }
    if !access::gpu_allowed() && !access::is_host_selector(selector) {
        return Err(BackendError::Codegen(
            "permission denied: the `gpu` capability is not granted here".into(),
        ));
    }
    let selected: Option<Arc<dyn device::Device>> = match selector {
        | "gpu" => device::registry().first_of(&device::DeviceKind::Gpu(device::GpuApi::Wgpu)),
        | "tiered" | "cranelift" | "jit" => None,
        | id => device::registry().get(id),
    };
    selected
        .map(|d| Arc::new(device::DeviceBackend::requiring_device(d)) as Shared)
        .ok_or_else(|| BackendError::Codegen(format!("backend `{selector}` is unavailable")))
}

/// Sets the process-wide default backend. `None` restores automatic selection.
pub fn set_default(backend: Option<Shared>) {
    if let Ok(mut slot) = DEFAULT_BACKEND.write() {
        *slot = backend;
    }
}

/// Sets the process default to a named backend and returns its canonical selector.
/// `auto` clears the process override. Thread scopes retain priority.
///
/// # Errors
/// Returns [`BackendError::Codegen`] when the requested backend is unavailable.
pub fn select_default(selector: &str) -> Result<String, BackendError> {
    let selected = named(selector)?;
    let name = selected.selector().to_owned();
    set_default((selector != "auto").then_some(selected));
    Ok(name)
}

struct ScopedBackend(Option<Shared>);

impl ScopedBackend {
    fn install(backend: Option<Shared>) -> Self {
        Self(SCOPED_BACKEND.with(|scope| scope.replace(backend)))
    }
}

impl Drop for ScopedBackend {
    fn drop(&mut self) {
        SCOPED_BACKEND.with(|scope| *scope.borrow_mut() = self.0.take());
    }
}

/// Runs `f` with `backend` as the current backend on this thread.
pub fn with_backend<R>(
    backend: Shared,
    f: impl FnOnce() -> R,
) -> R {
    let _scope = ScopedBackend::install(Some(backend));
    f()
}

/// The caller's explicit override and source grant, carried into worker jobs.
pub(crate) struct WorkerContext {
    backend: Option<Shared>,
    source_gpu: Option<bool>,
}

/// Restores the worker's previous backend and source scopes after a job.
pub(crate) struct WorkerScope {
    _backend: ScopedBackend,
    _source: Option<SourceAccess>,
}

impl WorkerContext {
    pub(crate) fn enter(&self) -> WorkerScope {
        WorkerScope {
            _backend: ScopedBackend::install(self.backend.clone()),
            _source: self.source_gpu.map(restrict_to_source),
        }
    }
}

pub(crate) fn worker_context() -> WorkerContext {
    WorkerContext {
        backend: SCOPED_BACKEND.with(|scope| scope.borrow().clone()),
        source_gpu: access::source_gpu_grant(),
    }
}

/// Compiles `root` with the [current] backend.
///
/// # Errors
/// See [`BackendError`].
pub fn compile(
    graph: &Graph,
    root: NodeId,
    inputs: &[SymbolId],
) -> Result<Box<dyn Compiled>, BackendError> {
    current().compile(graph, root, inputs)
}

#[derive(Clone, Copy)]
enum Instr {
    Const(f64),
    Input(u32),
    Call { f: EvalFn, start: u32, len: u32 },
}

/// The reference backend: a register tape in topological order.
///
/// Every distinct subterm becomes one instruction, so shared subterms are
/// computed once per call, exactly as a JIT would keep them in a register.
#[derive(Copy, Clone, Debug, Default)]
pub struct Interpreter;

struct Tape {
    arity: usize,
    instrs: Vec<Instr>,
    args: Vec<u32>,
}

impl Backend for Interpreter {
    fn uses_gpu(&self) -> bool {
        false
    }

    fn native_program_budget(&self) -> Option<NativeProgramBudget> {
        Some(NATIVE_TAPE_BUDGET)
    }

    fn retained_bytes(&self) -> Option<usize> {
        Some(size_of::<Self>())
    }

    fn name(&self) -> &'static str {
        "interpreter"
    }

    fn compile(
        &self,
        graph: &Graph,
        root: NodeId,
        inputs: &[SymbolId],
    ) -> Result<Box<dyn Compiled>, BackendError> {
        let mut tape = Tape {
            arity: inputs.len(),
            instrs: Vec::new(),
            args: Vec::new(),
        };
        let mut slot: crate::graph::id::FastMap<NodeId, u32> = crate::graph::id::FastMap::default();
        let mut stack = vec![(root, false)];
        while let Some((node, expanded)) = stack.pop() {
            if slot.contains_key(&node) {
                continue;
            }
            let children = graph.children(node);
            if !expanded && !children.is_empty() {
                stack.push((node, true));
                stack.extend(
                    children
                        .iter()
                        .filter(|c| !slot.contains_key(c))
                        .map(|&c| (c, false)),
                );
                continue;
            }
            let instr = if let Some(n) = graph.as_number(node) {
                Instr::Const(n.to_f64())
            } else if let Some(symbol) = graph.as_symbol(node) {
                let index = inputs.iter().position(|&s| s == symbol).ok_or_else(|| {
                    BackendError::UnboundSymbol(graph.interner().symbol_name(symbol).to_owned())
                })?;
                Instr::Input(u32::try_from(index).unwrap_or(u32::MAX))
            } else if graph.payload(node).is_some() {
                return Err(BackendError::NotNumeric);
            } else {
                let desc = graph.ops().get(graph.op(node));
                let f = desc
                    .eval
                    .ok_or_else(|| BackendError::NoSemantics(desc.name.to_string()))?;
                let start = u32::try_from(tape.args.len()).unwrap_or(u32::MAX);
                for child in children {
                    tape.args.push(slot.get(child).copied().unwrap_or(u32::MAX));
                }
                Instr::Call {
                    f,
                    start,
                    len: u32::try_from(children.len()).unwrap_or(u32::MAX),
                }
            };
            slot.insert(node, u32::try_from(tape.instrs.len()).unwrap_or(u32::MAX));
            tape.instrs.push(instr);
        }
        Ok(Box::new(tape))
    }
}

impl Tape {
    fn run(
        &self,
        args: &[f64],
        values: &mut Vec<f64>,
        scratch: &mut Vec<f64>,
    ) -> f64 {
        values.clear();
        for instr in &self.instrs {
            let value = match *instr {
                | Instr::Const(c) => c,
                | Instr::Input(i) => args.get(i as usize).copied().unwrap_or(f64::NAN),
                | Instr::Call { f, start, len } => {
                    scratch.clear();
                    let range = start as usize..(start as usize).saturating_add(len as usize);
                    for &arg in self.args.get(range).unwrap_or(&[]) {
                        scratch.push(values.get(arg as usize).copied().unwrap_or(f64::NAN));
                    }
                    f(scratch)
                },
            };
            values.push(value);
        }
        values.last().copied().unwrap_or(f64::NAN)
    }
}

impl Compiled for Tape {
    fn retained_bytes(&self) -> Option<usize> {
        let instructions = self.instrs.capacity().checked_mul(size_of::<Instr>())?;
        let arguments = self.args.capacity().checked_mul(size_of::<u32>())?;
        size_of::<Self>()
            .checked_add(instructions)?
            .checked_add(arguments)
    }

    fn arity(&self) -> usize {
        self.arity
    }

    fn call(
        &self,
        args: &[f64],
    ) -> f64 {
        thread_local! {
            static INTERP_SCRATCH: std::cell::RefCell<(Vec<f64>, Vec<f64>)> = const {
                std::cell::RefCell::new((Vec::new(), Vec::new()))
            };
        }
        INTERP_SCRATCH.with(|buffers| {
            match buffers.try_borrow_mut() {
                | Ok(mut pair) => {
                    let (values, scratch) = &mut *pair;
                    self.run(args, values, scratch)
                },
                | Err(_) => {
                    self.run(
                        args,
                        &mut Vec::with_capacity(self.instrs.len()),
                        &mut Vec::new(),
                    )
                },
            }
        })
    }

    fn call_batch(
        &self,
        columns: &[&[f64]],
        out: &mut [f64],
    ) {
        let mut values = crate::kernels::reuse::f64_scratch(self.instrs.len());
        let mut scratch = crate::kernels::reuse::f64_scratch(8);
        let mut args = crate::kernels::reuse::f64_scratch(self.arity);
        args.resize(self.arity, f64::NAN);
        for (row, slot) in out.iter_mut().enumerate() {
            for (arg, column) in args.iter_mut().zip(columns) {
                *arg = column.get(row).copied().unwrap_or(f64::NAN);
            }
            *slot = self.run(&args, &mut values, &mut scratch);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Engine;
    use crate::graph::Env;
    use crate::rules;

    fn setup(src: &str) -> (Graph, NodeId, Vec<SymbolId>) {
        let mut g = Graph::new();
        assert!(Engine::install(&mut g, &rules::standard()).is_ok());
        let root = g.parse(src).unwrap_or_else(|e| panic!("{e}"));
        let inputs = vec![g.interner_mut().symbol("x"), g.interner_mut().symbol("y")];
        (g, root, inputs)
    }

    #[test]
    fn restricted_source_skips_a_preselected_gpu_without_probing() {
        struct TestGpu;

        impl Backend for TestGpu {
            fn name(&self) -> &'static str {
                "accelerator-with-an-arbitrary-id"
            }

            fn uses_gpu(&self) -> bool {
                true
            }

            fn compile(
                &self,
                _: &Graph,
                _: NodeId,
                _: &[SymbolId],
            ) -> Result<Box<dyn Compiled>, BackendError> {
                Err(BackendError::Codegen(
                    "test device is discovery-only".into(),
                ))
            }
        }

        with_backend(Arc::new(TestGpu), || {
            assert!(current().uses_gpu());
            {
                let _source = restrict_to_source(false);
                assert!(!current().uses_gpu());
                assert!(named("accelerator-with-an-arbitrary-id").is_err());
                assert!(named("cpu").is_ok());
            }
            assert!(current().uses_gpu());
        });
    }

    #[test]
    fn worker_context_carries_override_and_keeps_the_narrower_grant() {
        let selected: Shared = Arc::new(Interpreter);
        let context = with_backend(Arc::clone(&selected), || {
            let _source = restrict_to_source(true);
            worker_context()
        });
        let worker = std::thread::spawn(move || {
            let _restricted = restrict_to_source(false);
            let _worker = context.enter();
            assert!(!gpu_allowed());
            assert!(Arc::ptr_eq(&current(), &selected));
        });
        assert!(
            worker.join().is_ok(),
            "worker should keep its narrowed scope"
        );
    }

    #[test]
    fn agrees_with_the_reference_evaluator() {
        let (g, root, inputs) = setup("sin(x)^2 * exp(-y) + (x - y)/(1 + x^2) + pi");
        let f = Interpreter
            .compile(&g, root, &inputs)
            .unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(f.arity(), 2);
        for (x, y) in [(0.0, 0.0), (1.5, -2.0), (-0.3, 4.0)] {
            let mut env = Env::numeric(0.0);
            env.bind(inputs[0], x);
            env.bind(inputs[1], y);
            let want = g.eval(root, &env).unwrap_or(f64::NAN);
            assert_eq!(f.call(&[x, y]).to_bits(), want.to_bits());
        }
    }

    #[test]
    fn batch_evaluation() {
        let (g, root, inputs) = setup("x^2 + y");
        let f = Interpreter
            .compile(&g, root, &inputs)
            .unwrap_or_else(|e| panic!("{e}"));
        let xs = [1.0, 2.0, 3.0];
        let ys = [10.0, 20.0, 30.0];
        let mut out = [0.0; 3];
        f.call_batch(&[&xs, &ys], &mut out);
        for (got, want) in out.iter().zip([11.0, 24.0, 39.0]) {
            assert!((got - want).abs() < 1e-9);
        }
    }

    #[test]
    fn checked_batch_sizes_blocks_by_work_without_reordering_rows() {
        use std::sync::atomic::AtomicUsize;
        use std::sync::atomic::Ordering;

        #[derive(Default)]
        struct CountingBatch(AtomicUsize);

        impl Compiled for CountingBatch {
            fn arity(&self) -> usize {
                1
            }

            fn call(
                &self,
                args: &[f64],
            ) -> f64 {
                args[0] * 2.0
            }

            fn call_batch(
                &self,
                columns: &[&[f64]],
                out: &mut [f64],
            ) {
                self.0.fetch_add(1, Ordering::Relaxed);
                for (slot, &value) in out.iter_mut().zip(columns[0]) {
                    *slot = self.call(&[value]);
                }
            }
        }

        let executor = CountingBatch::default();
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .unwrap_or_else(|error| panic!("{error}"));
        let input = (0..4096).map(f64::from).collect::<Vec<_>>();
        let mut out = vec![0.0; input.len()];
        pool.install(|| {
            for (work_per_row, expected_batches) in [(1, 1), (64, 4)] {
                out.fill(f64::NAN);
                executor.0.store(0, Ordering::Relaxed);
                call_batch_checked(&executor, &[&input], &mut out, work_per_row).unwrap();
                assert_eq!(executor.0.load(Ordering::Relaxed), expected_batches);
                for (&value, &actual) in input.iter().zip(&out) {
                    assert_eq!(actual.to_bits(), executor.call(&[value]).to_bits());
                }
            }
        });
    }

    #[test]
    fn checked_batch_stops_between_blocks_in_a_single_worker_pool() {
        use std::sync::atomic::AtomicBool;
        use std::sync::atomic::Ordering;

        struct CancelAfterBlock(Arc<AtomicBool>);

        impl Compiled for CancelAfterBlock {
            fn arity(&self) -> usize {
                0
            }

            fn call(
                &self,
                _: &[f64],
            ) -> f64 {
                1.0
            }

            fn call_batch(
                &self,
                _: &[&[f64]],
                out: &mut [f64],
            ) {
                out.fill(1.0);
                self.0.store(true, Ordering::Relaxed);
            }
        }

        let token = Arc::new(AtomicBool::new(false));
        let executor = CancelAfterBlock(Arc::clone(&token));
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .unwrap_or_else(|error| panic!("{error}"));
        let mut out = [0.0; 1024];
        let result = pool.install(|| {
            let _scope = crate::cancellation::install(Arc::clone(&token));
            call_batch_checked(&executor, &[], &mut out, 512)
        });
        assert_eq!(result, Err(crate::cancellation::Cancelled));
        let (first, rest) = out.split_at(256);
        assert!(first.iter().all(|v| f64::to_bits(*v) == 1.0f64.to_bits()));
        assert!(rest.iter().all(|v| f64::to_bits(*v) == 0.0f64.to_bits()));
    }

    #[test]
    fn shared_subterms_become_one_instruction() {
        let (mut g, _, inputs) = setup("x");
        let mut node = g.sym("x");
        for _ in 0..30 {
            node = g.node(crate::graph::op::core::POW, &[node, node]);
        }
        // 2^30 paths, 31 distinct subterms: must compile and run instantly.
        let f = Interpreter
            .compile(&g, node, &inputs)
            .unwrap_or_else(|e| panic!("{e}"));
        assert!((f.call(&[1.0, 0.0]) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn errors() {
        let (g, root, inputs) = setup("x + z");
        assert_eq!(
            Interpreter.compile(&g, root, &inputs).err(),
            Some(BackendError::UnboundSymbol("z".to_owned()))
        );
        let (g, root, inputs) = setup("diff(x^2, x)");
        assert_eq!(
            Interpreter.compile(&g, root, &inputs).err(),
            Some(BackendError::NoSemantics("diff".to_owned()))
        );
        assert_eq!(Interpreter.name(), "interpreter");
    }

    #[test]
    fn tape_function_batches_preserve_bits_across_worker_counts() {
        let (g, root, inputs) = setup("sin(x)^2 * exp(-y) + (x - y)/(1 + x^2)");
        let compiled = Interpreter
            .compile_multi(&g, &[root], &inputs, &[])
            .unwrap_or_else(|e| panic!("{e}"));
        let rows = 9_001_usize;
        let xs: Vec<f64> = (0..rows)
            .map(|i| f64::from(u32::try_from(i).unwrap_or(0)) * 0.013 - 4.0)
            .collect();
        let ys: Vec<f64> = (0..rows)
            .map(|i| f64::from(u32::try_from(i % 17).unwrap_or(0)) * 0.3 - 2.0)
            .collect();
        let mut reference = Vec::with_capacity(rows);
        let mut pt = [0.0];
        for (&x, &y) in xs.iter().zip(&ys) {
            compiled.eval(&[x, y], &[], &mut pt);
            reference.push(pt[0]);
        }
        for workers in [1_usize, 4] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(workers)
                .build()
                .unwrap_or_else(|e| panic!("{e}"));
            pool.install(|| {
                let mut out = vec![0.0; rows];
                compiled.eval_batch(&[&xs, &ys], &[], &mut [&mut out]);
                for (row, (&got, &want)) in out.iter().zip(&reference).enumerate() {
                    assert_eq!(got.to_bits(), want.to_bits(), "workers={workers} row={row}");
                }
            });
        }
    }
}
