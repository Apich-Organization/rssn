//! Session-owned compiled backends with explicit physical normalization.
//! Values and reference scales use scalar handles so heterogeneous dimensions
//! never masquerade as a homogeneous physical numeric tensor.
use super::numeric;
use super::output_scalar;
use super::runtime_dimensions;
use crate::api::Session;
use crate::backend::Backend;
use crate::backend::Compiled;
use crate::backend::Interpreter;
use crate::backend::{
    self,
};
use crate::graph::Arity;
use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::OpDescriptor;
use crate::graph::RuleError;
use crate::graph::SymbolId;
use crate::graph::op::EvalDyn;
use crate::graph::rule::Installer;
use crate::kernel::call::CallCtx;
use crate::kernel::call::OpError;
use crate::kernel::native_objects::NativeObjects;
use crate::kernel::signature::ArityInfo;
use crate::kernel::signature::BackendSupport;
use crate::kernel::signature::Determinism;
use crate::kernel::signature::DimExpr;
use crate::kernel::signature::Effect;
use crate::kernel::signature::EffectRow;
use crate::kernel::signature::NatExpr;
use crate::kernel::signature::Purity;
use crate::kernel::signature::RepExpr;
use crate::kernel::signature::Signature;
use crate::kernel::signature::Stability;
use crate::kernel::signature::ThreadSafety;
use crate::kernel::signature::TypeExpr;
use crate::kernel::signature::TypeParam;
use crate::kernel::signature::TypeParamKind;
use crate::kernel::value::DynValue;
use crate::kernel::value::TensorValue;
use crate::kernel::value::Value;
use num_traits::ToPrimitive;
use std::collections::BTreeSet;
use std::sync::Arc;
const BACKEND: &str = "rssn/native/backend/v1";
const SCALAR: &str = "rssn/native/backend-scalar/v1";
const PROGRAM: &str = "rssn/native/compiled/v1";
const MAX_INPUTS: usize = 63;
const MAX_ROWS: usize = 16_384;
const MAX_CELLS: usize = 262_144;
const MAX_NODES: usize = 4096;
const MAX_WORK: usize = 100_000_000;
const PROGRAM_SESSION_RESERVATION: usize = 4_194_304;
#[derive(Clone)]
struct Scalar {
    value: f64,
    // Stored in canonical base order when the scalar handle is created.
    dimension: DimExpr,
}
struct Program {
    payload_bytes: usize,
    retained_bytes: usize,
    compiled: Arc<dyn Compiled>,
    bounded_point_tape: Option<Arc<backend::TapeFunction>>,
    // Retain every object owning graph/plugin function pointers for at least
    // as long as the executable, including a snapshot of its actual backend.
    _session: Session,
    backend: Arc<dyn Backend>,
    references: Vec<Scalar>,
    output: Scalar,
    work_per_row: usize,
    work_budget: usize,
}

impl Program {
    fn check_batch(
        &self,
        columns: usize,
        rows: usize,
    ) -> Result<(), OpError> {
        if columns != self.compiled.arity()
            || rows > MAX_ROWS
            || self
                .work_per_row
                .checked_mul(rows)
                .is_none_or(|n| n > self.work_budget)
        {
            Err(error("compiled column shape, row or work budget exceeded"))
        } else {
            Ok(())
        }
    }

    fn batch_result(
        &self,
        out: Vec<f64>,
    ) -> Result<Value, OpError> {
        let rows = out.len();
        let data = out
            .into_iter()
            .map(|value| output_scalar(value * self.output.value))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Value::Tensor(TensorValue {
            shape: vec![rows],
            elem: TypeExpr::Float.with_dim(self.output.dimension.clone()),
            data,
        }))
    }
}

struct SelectedCompile(Arc<dyn Backend>, bool);
impl Backend for SelectedCompile {
    fn uses_gpu(&self) -> bool {
        self.0.uses_gpu()
    }

    fn name(&self) -> &'static str {
        self.0.name()
    }

    fn compile(
        &self,
        graph: &Graph,
        root: NodeId,
        inputs: &[SymbolId],
    ) -> Result<Box<dyn Compiled>, backend::BackendError> {
        let compile = || self.0.compile_native(graph, root, inputs);
        if self.1 {
            backend::with_backend(self.0.clone(), compile)
        } else {
            compile()
        }
    }
}
fn error(message: impl std::fmt::Display) -> OpError {
    OpError::Domain(message.to_string())
}
fn unavailable(
    selector: &str,
    error: &backend::BackendError,
) -> OpError {
    OpError::Unsupported {
        backend: selector.into(),
        feature: error.to_string(),
    }
}
fn objects<'a>(ctx: &'a CallCtx<'_>) -> Result<&'a NativeObjects, OpError> {
    ctx.native_objects
        .as_ref()
        .ok_or_else(|| error("backend objects require a session call context"))
}
fn canonical(dim: &DimExpr) -> Result<DimExpr, OpError> {
    runtime_dimensions::expression(&runtime_dimensions::canonical(dim)?)
}
fn dimension(value: &Value) -> Result<DimExpr, OpError> {
    canonical(match value {
        | Value::Dyn(d) => &d.dim,
        | Value::Quantity(q) => &q.unit.dimension,
        | _ => &DimExpr::Dimless,
    })
}
fn dyn_type() -> TypeExpr {
    TypeExpr::Dyn(RepExpr {
        repr: "term".into(),
        dim: None,
    })
}
fn tensor_type(
    rank: u64,
    elem: TypeExpr,
) -> TypeExpr {
    TypeExpr::Tensor {
        rank: NatExpr::Lit(rank),
        shape: None,
        elem: Box::new(elem),
    }
}
fn opaque(kind: &str) -> TypeExpr {
    TypeExpr::Opaque(kind.into())
}
fn integer(value: &Value) -> Result<usize, OpError> {
    let Value::Int(value) = value else {
        return Err(error("expected a nonnegative integer"));
    };
    value
        .to_usize()
        .ok_or_else(|| error("integer does not fit platform size"))
}
fn text(value: &Value) -> Result<&str, OpError> {
    let Value::String(value) = value else {
        return Err(error("expected text"));
    };
    Ok(value)
}
fn handles(
    value: &Value,
    rank: usize,
) -> Result<&TensorValue, OpError> {
    let Value::Tensor(t) = value else {
        return Err(error("expected a tensor of native scalar handles"));
    };
    if t.shape.len() != rank
        || t.shape.iter().try_fold(1usize, |a, b| a.checked_mul(*b)) != Some(t.data.len())
        || t.data.len() > MAX_CELLS
    {
        return Err(error(
            "malformed scalar-handle tensor or allocation budget exceeded",
        ));
    }
    Ok(t)
}
fn scalar(
    bank: &NativeObjects,
    value: &Value,
) -> Result<Scalar, OpError> {
    bank.with::<Scalar, _>(value, SCALAR, |v| Ok(v.clone()))
}
fn normalized_scalar(
    bank: &NativeObjects,
    handle: &Value,
    reference: &Scalar,
    nonfinite: &'static str,
) -> Result<f64, OpError> {
    bank.with::<Scalar, _>(handle, SCALAR, |value| {
        if value.dimension != reference.dimension {
            return Err(error(
                "compiled input dimension does not match its declared reference",
            ));
        }
        let normalized = value.value / reference.value;
        if normalized.is_finite() {
            Ok(normalized)
        } else {
            Err(error(nonfinite))
        }
    })
}
fn measured(
    value: f64,
    dimension: DimExpr,
) -> Result<Value, OpError> {
    let value = output_scalar(value)?;
    if dimension == DimExpr::Dimless {
        Ok(value)
    } else {
        Ok(Value::Dyn(DynValue {
            value: Box::new(value),
            dim: dimension,
        }))
    }
}
fn normalized_point(
    program: &Program,
    bank: &NativeObjects,
    values: &[Value],
) -> Result<Vec<f64>, OpError> {
    if values.len() != program.references.len() {
        return Err(error("compiled input count differs from declared arity"));
    }
    values
        .iter()
        .zip(&program.references)
        .map(|(handle, reference)| {
            normalized_scalar(
                bank,
                handle,
                reference,
                "normalized compiled input is not finite",
            )
        })
        .collect()
}
fn strings(value: &Value) -> Result<Vec<String>, OpError> {
    let t = handles(value, 1)?;
    if t.data.len() > MAX_INPUTS {
        return Err(error("too many compiled input symbols"));
    }
    let mut names = BTreeSet::new();
    let mut result = Vec::with_capacity(t.data.len());
    for value in &t.data {
        let name = text(value)?;
        let mut characters = name.chars();
        if name.len() > 128
            || !characters
                .next()
                .is_some_and(|c| c == '_' || c.is_ascii_alphabetic())
            || !characters.all(|c| c == '_' || c.is_ascii_alphanumeric())
            || !names.insert(name)
        {
            return Err(error(
                "compiled input names must be unique ASCII identifiers of at most 128 bytes",
            ));
        }
        result.push(name.into());
    }
    Ok(result)
}
fn source(text: &str) -> Result<(), OpError> {
    if text.is_empty() || text.len() > 16_384 || !text.is_ascii() {
        return Err(error("compiled expressions require 1..=16384 ASCII bytes"));
    }
    let mut depth = 0usize;
    let mut operations = 0usize;
    for c in text.chars() {
        match c {
            | '(' | '[' => {
                depth += 1;
                if depth > 128 {
                    return Err(error("compiled expression nesting exceeds 128"));
                }
            },
            | ')' | ']' => depth = depth.saturating_sub(1),
            | '+' | '-' | '*' | '/' | '^' => operations += 1,
            | _ => {},
        }
    }
    if operations > 1024 {
        return Err(error("compiled expression operation budget exceeded"));
    }
    Ok(())
}
#[derive(Clone, Copy)]
enum Method {
    Current,
    CurrentName,
    Available,
    Named,
    SelectDefault,
    Interpreter,
    Name,
    Reference,
    Value,
    Compile,
    CompileCurrent,
    CompileScoped,
    Arity,
    Call,
    Batch,
    ParallelBatch,
    PhysicalBatch,
    SetDefault,
    Release,
}
struct Adapter(Method);
impl EvalDyn for Adapter {
    fn signature(&self) -> Signature {
        use Method as M;
        let scalar = || opaque(SCALAR);
        let program = || opaque(PROGRAM);
        let selected = || opaque(BACKEND);
        let (params, ret) = match self.0 {
            | M::Current | M::Interpreter => (vec![], selected()),
            | M::CurrentName => (vec![], TypeExpr::String),
            | M::Available => (vec![], tensor_type(1, TypeExpr::String)),
            | M::Named => (vec![TypeExpr::String], selected()),
            | M::SelectDefault => (vec![TypeExpr::String], TypeExpr::String),
            | M::Name => (vec![selected()], TypeExpr::String),
            | M::Reference | M::Value => {
                (
                    vec![TypeExpr::Float.with_dim(DimExpr::Var("D".into()))],
                    scalar(),
                )
            },
            | M::Compile | M::CompileScoped => {
                (
                    vec![
                        selected(),
                        TypeExpr::String,
                        tensor_type(1, TypeExpr::String),
                        tensor_type(1, scalar()),
                        scalar(),
                        TypeExpr::Int,
                    ],
                    program(),
                )
            },
            | M::CompileCurrent => {
                (
                    vec![
                        TypeExpr::String,
                        tensor_type(1, TypeExpr::String),
                        tensor_type(1, scalar()),
                        scalar(),
                        TypeExpr::Int,
                    ],
                    program(),
                )
            },
            | M::Arity => (vec![program()], TypeExpr::Int),
            | M::Call => (vec![program(), tensor_type(1, scalar())], dyn_type()),
            | M::Batch | M::ParallelBatch => {
                (vec![program(), tensor_type(2, scalar())], dyn_type())
            },
            | M::PhysicalBatch => {
                (
                    vec![program(), TypeExpr::Int, tensor_type(1, dyn_type())],
                    dyn_type(),
                )
            },
            | M::SetDefault => (vec![selected(), TypeExpr::String], TypeExpr::String),
            | M::Release => (vec![dyn_type()], TypeExpr::Bool),
        };
        let mut signature = Signature::new(params, ret);
        signature.effects = EffectRow {
            effects: vec![Effect::Stateful, Effect::Diverges],
        };
        if matches!(self.0, M::Reference | M::Value) {
            signature.type_params.push(TypeParam {
                name: "D".into(),
                kind: TypeParamKind::Dim,
                bounds: vec![],
            });
        }
        signature
    }

    #[allow(clippy::too_many_lines)] // Each native family keeps its preflight checks beside dispatch.
    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        use Method as M;
        let count = self.signature().params.len();
        if args.len() != count {
            return Err(OpError::Arity {
                expected: ArityInfo::fixed(count as u8),
                got: args.len(),
            });
        }
        let bank = objects(ctx)?;
        match self.0 {
            | M::Current | M::Interpreter | M::Named => {
                let selected: Arc<dyn Backend> = match self.0 {
                    | M::Named => {
                        let selector = text(&args[0])?;
                        backend::access::require_gpu_if(!backend::access::is_host_selector(
                            selector,
                        ))?;
                        backend::named(selector).map_err(|err| unavailable(selector, &err))?
                    },
                    | M::Current => backend::current(),
                    | _ => Arc::new(Interpreter),
                };
                bank.insert(BACKEND, vec![], selected, 1024)
            },
            | M::CurrentName => Ok(Value::String(backend::current().selector().into())),
            | M::Available => {
                let names = backend::available();
                Ok(Value::Tensor(TensorValue {
                    shape: vec![names.len()],
                    elem: TypeExpr::String,
                    data: names.into_iter().map(Value::String).collect(),
                }))
            },
            | M::SelectDefault => {
                backend::access::require_process_default()?;
                let selector = text(&args[0])?;
                let active =
                    backend::select_default(selector).map_err(|err| unavailable(selector, &err))?;
                Ok(Value::String(active))
            },
            | M::Name => {
                bank.with::<Arc<dyn Backend>, _>(&args[0], BACKEND, |backend| {
                    backend::access::require_gpu_if(backend.uses_gpu())?;
                    Ok(Value::String(backend.selector().into()))
                })
            },
            | M::Reference | M::Value => {
                let value = numeric(&args[0], 0)?;
                if matches!(self.0, M::Reference) && value <= 0.0 {
                    return Err(error("physical reference values must be positive"));
                }
                let dimension = dimension(&args[0])?;
                bank.insert(
                    SCALAR,
                    vec![dimension.clone()],
                    Scalar { value, dimension },
                    256,
                )
            },
            | M::Compile | M::CompileCurrent | M::CompileScoped => {
                let offset = usize::from(!matches!(self.0, M::CompileCurrent));
                let chosen = if offset == 0 {
                    backend::current()
                } else {
                    bank.with::<Arc<dyn Backend>, _>(&args[0], BACKEND, |value| Ok(value.clone()))?
                };
                backend::access::require_gpu_if(chosen.uses_gpu())?;
                let budget = chosen.native_program_budget().ok_or_else(|| {
                    OpError::Unsupported {
                        backend: chosen.selector().into(),
                        feature: "bounded native expression compilation".into(),
                    }
                })?;
                let expression = text(&args[offset])?;
                source(expression)?;
                let names = strings(&args[offset + 1])?;
                let refs = handles(&args[offset + 2], 1)?;
                if refs.data.len() != names.len() {
                    return Err(error(
                        "one physical reference is required per declared input",
                    ));
                }
                let references = refs
                    .data
                    .iter()
                    .map(|value| scalar(bank, value))
                    .collect::<Result<Vec<_>, _>>()?;
                if references.iter().any(|r| r.value <= 0.0) {
                    return Err(error("compiled coordinate references must be positive"));
                }
                let output = scalar(bank, &args[offset + 3])?;
                if output.value <= 0.0 {
                    return Err(error("compiled output reference must be positive"));
                }
                let work_budget = integer(&args[offset + 4])?;
                if work_budget == 0 || work_budget > MAX_WORK {
                    return Err(error("compiled work budget must be 1..=100000000"));
                }
                let session = Session::with_rules(&[crate::rules::elementary()]).map_err(error)?;
                let term = session.parse(expression).map_err(error)?;
                let borrowed = names.iter().map(String::as_str).collect::<Vec<_>>();
                let tape = term.tape(&borrowed).map_err(error)?;
                let nodes = tape.insts.len();
                if nodes > MAX_NODES.min(budget.instructions) {
                    return Err(error(format!(
                        "compiled tape exceeds {} nodes for backend `{}`",
                        MAX_NODES.min(budget.instructions),
                        chosen.selector()
                    )));
                }
                let work_per_row = nodes
                    .checked_add(tape.args.len())
                    .ok_or_else(|| error("compiled work count overflows"))?
                    .max(1);
                if work_per_row > work_budget {
                    return Err(error("work budget cannot evaluate even one compiled row"));
                }
                let payload_bytes = (expression.len() * 64)
                    .checked_add(work_per_row * 512)
                    .ok_or_else(|| error("compiled payload-memory charge overflows"))?;
                let retained_bytes = PROGRAM_SESSION_RESERVATION
                    .checked_add(payload_bytes)
                    .and_then(|n| n.checked_add(budget.extra_bytes))
                    .ok_or_else(|| error("compiled retained-memory charge overflows"))?;
                let compiler = SelectedCompile(chosen.clone(), matches!(self.0, M::CompileScoped));
                let compiled = term.compile_with(&compiler, &borrowed).map_err(error)?;
                if compiled.arity() != references.len() {
                    return Err(error("backend returned an inconsistent compiled arity"));
                }
                // A solver that detaches a point callback needs an executor with
                // a retained-memory bound. Reuse the preflight tape when the
                // selected executor has no such bound.
                let bounded_point_tape =
                    if compiled.retained_bytes().is_none() || chosen.retained_bytes().is_none() {
                        Some(Arc::new(backend::TapeFunction::from_tape(tape)))
                    } else {
                        None
                    };
                let reserved_bytes = if let Some(tape) = &bounded_point_tape {
                    tape.retained_bytes()
                        .and_then(|bytes| retained_bytes.checked_add(bytes))
                        .ok_or_else(|| error("compiled point tape memory charge overflows"))?
                } else {
                    retained_bytes
                };
                let mut dimensions = references
                    .iter()
                    .map(|r| r.dimension.clone())
                    .collect::<Vec<_>>();
                dimensions.push(output.dimension.clone());
                bank.insert(
                    PROGRAM,
                    dimensions,
                    Program {
                        payload_bytes,
                        retained_bytes: reserved_bytes,
                        compiled: Arc::from(compiled),
                        bounded_point_tape,
                        _session: session,
                        backend: chosen,
                        references,
                        output,
                        work_per_row,
                        work_budget,
                    },
                    reserved_bytes,
                )
            },
            | M::Arity => {
                bank.with::<Program, _>(&args[0], PROGRAM, |p| {
                    backend::access::require_gpu_if(p.backend.uses_gpu())?;
                    Ok(Value::int(p.compiled.arity()))
                })
            },
            | M::Call => {
                let values = handles(&args[1], 1)?;
                bank.with::<Program, _>(&args[0], PROGRAM, |p| {
                    backend::access::require_gpu_if(p.backend.uses_gpu())?;
                    let point = normalized_point(p, bank, &values.data)?;
                    measured(
                        p.compiled.call(&point) * p.output.value,
                        p.output.dimension.clone(),
                    )
                })
            },
            | M::Batch | M::ParallelBatch => {
                let values = handles(&args[1], 2)?;
                bank.with::<Program, _>(&args[0], PROGRAM, |p| {
                    backend::access::require_gpu_if(p.backend.uses_gpu())?;
                    let rows = values.shape[1];
                    p.check_batch(values.shape[0], rows)?;
                    let parallel = matches!(self.0, M::ParallelBatch);
                    if parallel {
                        crate::cancellation::checkpoint()?;
                    }
                    let mut columns = Vec::with_capacity(values.shape[0]);
                    for (i, reference) in p.references.iter().enumerate() {
                        let column = values.data[i * rows..(i + 1) * rows]
                            .iter()
                            .map(|handle| {
                                normalized_scalar(
                                    bank,
                                    handle,
                                    reference,
                                    "normalized batch input is not finite",
                                )
                            })
                            .collect::<Result<Vec<_>, _>>()?;
                        columns.push(column);
                    }
                    let columns = columns.iter().map(Vec::as_slice).collect::<Vec<_>>();
                    let mut out = vec![0.0; rows];
                    if parallel {
                        backend::call_batch_checked(
                            &*p.compiled,
                            &columns,
                            &mut out,
                            p.work_per_row,
                        )?;
                    } else {
                        p.compiled.call_batch(&columns, &mut out);
                    }
                    p.batch_result(out)
                })
            },
            | M::PhysicalBatch => {
                let rows = integer(&args[1])?;
                bank.with::<Program, _>(&args[0], PROGRAM, |p| {
                    backend::access::require_gpu_if(p.backend.uses_gpu())?;
                    p.check_batch(p.references.len(), rows)?;
                    crate::cancellation::checkpoint()?;
                    let scales = p.references.iter().map(|s| (s.value, &s.dimension));
                    let columns = super::native_backend_extended::physical_columns(
                        ctx, &args[2], rows, scales,
                    )?;
                    let slices = columns.iter().map(Vec::as_slice).collect::<Vec<_>>();
                    let mut out = vec![0.0; rows];
                    backend::call_batch_checked(&*p.compiled, &slices, &mut out, p.work_per_row)?;
                    p.batch_result(out)
                })
            },
            | M::SetDefault => {
                backend::access::require_process_default()?;
                let chosen =
                    bank.with::<Arc<dyn Backend>, _>(&args[0], BACKEND, |value| Ok(value.clone()))?;
                match text(&args[1])? {
                    | "set_process_default" => {
                        let name = chosen.selector().to_owned();
                        backend::set_default(Some(chosen));
                        Ok(Value::String(name))
                    },
                    | "restore_interpreter" => {
                        backend::select_default("interpreter")
                            .map(Value::String)
                            .map_err(error)
                    },
                    | _ => {
                        Err(error(
                            "explicit action must be set_process_default or restore_interpreter",
                        ))
                    },
                }
            },
            | M::Release => {
                bank.release(&args[0])?;
                Ok(Value::Bool(true))
            },
        }
    }
}
pub(super) fn install(installer: &mut Installer<'_>) -> Result<(), RuleError> {
    for (name, path, method) in METHODS {
        let adapter = Adapter(*method);
        let signature = adapter.signature();
        let documentation = if matches!(method, Method::PhysicalBatch) {
            "Evaluate a retained scalar program over measured rank-one columns in declared symbol order. Each column has its own physical unit and matches the explicit row count. Returns a single measured vector in input order. The active CPU pool splits sufficiently large work into fixed independent blocks. Cancellation and the program's row and work limits apply."
        } else if matches!(method, Method::CurrentName) {
            "Return the canonical active backend selector without allocating a handle."
        } else if matches!(method, Method::Available) {
            "List valid backend selectors discovered on this host. The first device query can initialize GPU adapters."
        } else if matches!(method, Method::Named) {
            "Resolve a named backend into a session-owned handle for measured compilation. An unavailable selector returns Unsupported."
        } else if matches!(method, Method::SelectDefault) {
            "Set a named process default and return its resolved selector. Use auto to clear the override. Thread-local overrides keep priority."
        } else if matches!(method, Method::ParallelBatch) {
            "Evaluate a retained scalar program over a tensor of session-owned input handles with shape `[arity, rows]`. Columns follow declared symbol order and retain their individual physical dimensions.\n\nThe executor estimates the work per row from the compiled program and uses the active CPU pool for independent blocks when there is enough work.\n\nThe measured result preserves row order. Cancellation is checked between blocks, and shape, units, and the program’s row and work limits apply."
        } else {
            "Session-owned native backend/compiled/scalar handles. Compile an elementary ASCII expression in declared symbol order, normalized by positive physical coordinate references and a positive output reference. Scalar-value handles preserve mixed dimensions individually; call columns are [arity,rows] handles, not mixed-unit numeric tensors. The compiled object retains its private elementary RSSN Session and actual backend. Compilation enforces source, allocation and per-call work bounds for the selected backend. Current snapshots environment; scoped compile selects its backend explicitly. Process-default changes require an explicit action string. Kind/session/lifetime and units are checked; release invalidates every handle copy."
        };
        let mut descriptor = OpDescriptor::new(name, Arity::Fixed(signature.params.len() as u8))
            .category("numerical_matrix")
            .doc(documentation)
            .stability(Stability::Experimental)
            .purity(Purity::Effectful)
            .determinism(Determinism::EnvDependent)
            .thread_safety(ThreadSafety::SessionBound)
            // Reachable from C: the session-owned handles these methods take
            // and return cross the value ABI as opaque cells
            // (`rssn_value_call`), which is exactly what they are.
            .backend(BackendSupport {
                interpreter: true,
                jit: false,
                ffi: true,
                gpu: false,
            })
            .signature(signature)
            .eval_dyn(Arc::new(adapter));
        if !path.is_empty() {
            descriptor = descriptor.rust_path(path);
        }
        let op = installer.op(descriptor)?;
        installer.dim(op, crate::graph::dim::DimRule::Container);
    }
    Ok(())
}
const METHODS: &[(&str, &str, Method)] = &[
    (
        "num_backend_current",
        "crate::backend::current",
        Method::Current,
    ),
    ("num_backend_current_name", "", Method::CurrentName),
    (
        "num_backend_available",
        "crate::backend::available",
        Method::Available,
    ),
    ("num_backend_named", "crate::backend::named", Method::Named),
    (
        "num_backend_select_default",
        "crate::backend::select_default",
        Method::SelectDefault,
    ),
    ("num_backend_interpreter", "", Method::Interpreter),
    (
        "num_backend_name",
        "crate::backend::Backend::name",
        Method::Name,
    ),
    ("num_backend_reference", "", Method::Reference),
    ("num_backend_value", "", Method::Value),
    (
        "num_backend_compile",
        "crate::backend::Backend::compile",
        Method::Compile,
    ),
    (
        "num_backend_compile_current",
        "crate::backend::compile",
        Method::CompileCurrent,
    ),
    (
        "num_backend_with_backend_compile",
        "crate::backend::with_backend",
        Method::CompileScoped,
    ),
    (
        "num_backend_arity",
        "crate::backend::Compiled::arity",
        Method::Arity,
    ),
    (
        "num_backend_call",
        "crate::backend::Compiled::call",
        Method::Call,
    ),
    (
        "num_backend_call_batch",
        "crate::backend::Compiled::call_batch",
        Method::Batch,
    ),
    (
        "num_backend_call_batch_par",
        "crate::backend::call_batch_par",
        Method::ParallelBatch,
    ),
    ("num_backend_call_batch_columns", "", Method::PhysicalBatch),
    (
        "num_backend_set_default",
        "crate::backend::set_default",
        Method::SetDefault,
    ),
    ("num_backend_release", "", Method::Release),
];

/// Immutable snapshot of an elementary-only compiled function. The executable
/// and backend own their code; the constructor session has only static kernels.
#[derive(Clone)]
pub(in crate::rules) struct CompiledSnapshot {
    payload_bytes: usize,
    pub retained_bytes: usize,
    compiled: Arc<dyn Compiled>,
    bounded_point_tape: Option<Arc<backend::TapeFunction>>,
    backend: Arc<dyn Backend>,
    pub input_dimensions: Vec<DimExpr>,
    pub output_dimension: DimExpr,
    pub input_scales: Vec<f64>,
    pub output_scale: f64,
    remaining: Arc<std::sync::atomic::AtomicUsize>,
    failure: Arc<std::sync::Mutex<Option<OpError>>>,
}
impl CompiledSnapshot {
    /// The source and instruction reservation for a program kept alive by a
    /// native parent handle. The parent separately reserves the session and
    /// backend storage for the entire lifetime of its dependents.
    #[must_use]
    pub const fn program_payload_bytes(&self) -> usize {
        self.payload_bytes
    }

    fn charge(&self) -> Result<(), OpError> {
        use std::sync::atomic::Ordering;
        let mut n = self.remaining.load(Ordering::Relaxed);
        loop {
            if n == 0 {
                return Err(error("compiled function evaluation budget exhausted"));
            }
            match self.remaining.compare_exchange_weak(
                n,
                n - 1,
                Ordering::Relaxed,
                Ordering::Relaxed,
            ) {
                | Ok(_) => return Ok(()),
                | Err(next) => n = next,
            }
        }
    }

    pub fn evaluate_normalized(
        &self,
        x: &[f64],
    ) -> Result<f64, OpError> {
        self.charge()?;
        if x.len() != self.input_scales.len() || x.iter().any(|v| !v.is_finite()) {
            return Err(error(
                "compiled function inputs have wrong arity or nonfinite values",
            ));
        }
        let y = self.compiled.call(x);
        if !y.is_finite() {
            return Err(error("compiled function result is not finite"));
        }
        Ok(y)
    }

    pub fn evaluate_coherent(
        &self,
        x: &[f64],
    ) -> Result<f64, OpError> {
        if x.len() != self.input_scales.len() {
            return Err(error("compiled function arity mismatch"));
        }
        let evaluate = |args: &mut [f64]| {
            for ((arg, &value), &scale) in args.iter_mut().zip(x).zip(&self.input_scales) {
                *arg = value / scale;
            }
            self.evaluate_normalized(args)
        };
        // Avoid global scratch-pool traffic for short, independently invoked callbacks.
        let y = match x.len() {
            | 0 => self.evaluate_normalized(x),
            | 1..=8 => evaluate(&mut [0.0; 8][..x.len()]),
            | 9..=MAX_INPUTS => evaluate(&mut [0.0; MAX_INPUTS][..x.len()]),
            | _ => return Err(error("compiled function arity exceeds the snapshot bound")),
        }? * self.output_scale;
        if !y.is_finite() {
            return Err(error("scaled compiled function result is not finite"));
        }
        Ok(y)
    }

    pub fn sample(
        &self,
        x: &[f64],
    ) -> f64 {
        match self.evaluate_coherent(x) {
            | Ok(y) => y,
            | Err(e) => {
                if let Ok(mut failure) = self.failure.lock()
                    && failure.is_none()
                {
                    *failure = Some(e);
                }
                f64::NAN
            },
        }
    }

    /// Switch an unbounded executor to its shared CPU point tape.
    ///
    /// The physical scales and the evaluation ledger are kept intact.
    pub fn bound_detached_point_storage(&mut self) -> bool {
        let Some(tape) = self.bounded_point_tape.take() else {
            return false;
        };
        self.compiled = tape;
        self.backend = Arc::new(Interpreter);
        true
    }

    /// Return the identity and charge of the shared executable plus the cost
    /// of each retained snapshot. The program's constructor Session is absent.
    pub fn detached_storage(&self) -> Option<(usize, usize, usize)> {
        fn dimension_bytes(dimension: &DimExpr) -> usize {
            match dimension {
                | DimExpr::Base(name) | DimExpr::Var(name) => name.capacity(),
                | DimExpr::Mul(factors) => {
                    factors
                        .capacity()
                        .saturating_mul(size_of::<DimExpr>())
                        .saturating_add(
                            factors
                                .iter()
                                .map(dimension_bytes)
                                .fold(0, usize::saturating_add),
                        )
                },
                | DimExpr::Inverse(inner) => {
                    size_of::<DimExpr>().saturating_add(dimension_bytes(inner))
                },
                | DimExpr::Pow(inner, exponent) => {
                    let symbol = match exponent {
                        | crate::kernel::signature::Exp::Symbol(symbol) => symbol.capacity(),
                        | crate::kernel::signature::Exp::Rat(_) => 0,
                    };
                    size_of::<DimExpr>()
                        .saturating_add(dimension_bytes(inner))
                        .saturating_add(symbol)
                },
                | _ => 0,
            }
        }

        let key = Arc::as_ptr(&self.compiled).cast::<()>() as usize;
        let executable = self
            .compiled
            .retained_bytes()?
            .checked_add(self.backend.retained_bytes()?)?;
        let dimensions = self
            .input_dimensions
            .iter()
            .chain(std::iter::once(&self.output_dimension))
            .map(dimension_bytes)
            .fold(0, usize::saturating_add);
        // Covers Arc headers, the per-use failure message, bookkeeping,
        // and allocation slack for a solver's Arc and callback entries.
        let snapshot = size_of::<Self>()
            .saturating_add(size_of::<std::sync::Mutex<Option<OpError>>>())
            .saturating_add(512)
            .saturating_add(
                self.input_scales
                    .capacity()
                    .saturating_mul(size_of::<f64>()),
            )
            .saturating_add(
                self.input_dimensions
                    .capacity()
                    .saturating_mul(size_of::<DimExpr>()),
            )
            .saturating_add(dimensions);
        Some((key, executable, snapshot))
    }

    pub fn finish(&self) -> Result<(), OpError> {
        let failure = self
            .failure
            .lock()
            .map_err(|_| error("compiled function failure state poisoned"))?;
        failure.as_ref().map_or(Ok(()), |e| Err(e.clone()))
    }
}
pub(in crate::rules) fn compiled_snapshot(
    ctx: &CallCtx<'_>,
    handle: &Value,
    max_evaluations: usize,
) -> Result<CompiledSnapshot, OpError> {
    if max_evaluations == 0 || max_evaluations > 10_000_000 {
        return Err(error(
            "compiled function evaluation limit must be 1..=10000000",
        ));
    }
    objects(ctx)?.with::<Program, _>(handle, PROGRAM, |p| {
        Ok(CompiledSnapshot {
            payload_bytes: p.payload_bytes,
            retained_bytes: p.retained_bytes,
            compiled: p.compiled.clone(),
            bounded_point_tape: p.bounded_point_tape.clone(),
            backend: p.backend.clone(),
            input_dimensions: p.references.iter().map(|r| r.dimension.clone()).collect(),
            output_dimension: p.output.dimension.clone(),
            input_scales: p.references.iter().map(|r| r.value).collect(),
            output_scale: p.output.value,
            remaining: Arc::new(std::sync::atomic::AtomicUsize::new(
                max_evaluations.min(p.work_budget / p.work_per_row),
            )),
            failure: Arc::new(std::sync::Mutex::new(None)),
        })
    })
}

pub(in crate::rules) fn scalar_value(
    ctx: &CallCtx<'_>,
    handle: &Value,
) -> Result<(f64, DimExpr), OpError> {
    objects(ctx)?.with::<Scalar, _>(handle, SCALAR, |s| Ok((s.value, s.dimension.clone())))
}
pub(in crate::rules) fn scalar_reference(
    ctx: &CallCtx<'_>,
    handle: &Value,
) -> Result<(f64, DimExpr), OpError> {
    let (value, dim) = scalar_value(ctx, handle)?;
    if value <= 0.0 {
        return Err(error("reference scalar must be positive"));
    }
    Ok((value, dim))
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Sum(usize);

    impl Compiled for Sum {
        fn arity(&self) -> usize {
            self.0
        }

        fn call(
            &self,
            args: &[f64],
        ) -> f64 {
            assert_eq!(self.0, args.len());
            args.iter().sum()
        }
    }

    #[test]
    fn coherent_callbacks_share_normalization_and_budget_across_arities() {
        for arity in [0, 2, 8, 9, MAX_INPUTS] {
            let scales: Vec<_> = (0..arity).map(|i| (i + 1) as f64).collect();
            let values: Vec<_> = scales.iter().map(|scale| 0.5 * scale).collect();
            let snapshot = CompiledSnapshot {
                payload_bytes: 0,
                retained_bytes: 0,
                compiled: Arc::new(Sum(arity)),
                bounded_point_tape: None,
                backend: Arc::new(Interpreter),
                input_dimensions: vec![DimExpr::Dimless; arity],
                output_dimension: DimExpr::Dimless,
                input_scales: scales,
                output_scale: 4.0,
                remaining: Arc::new(std::sync::atomic::AtomicUsize::new(2)),
                failure: Arc::new(std::sync::Mutex::new(None)),
            };
            assert!(snapshot.evaluate_coherent(&vec![0.0; arity + 1]).is_err());
            let expected = 4.0 * std::iter::repeat_n(0.5, arity).sum::<f64>();
            for _ in 0..2 {
                assert_eq!(
                    snapshot.evaluate_coherent(&values).unwrap().to_bits(),
                    expected.to_bits()
                );
            }
            assert!(snapshot.evaluate_coherent(&values).is_err());
        }
    }
}
