//! Checked session objects for reverse AD. Handles never serialize native state.
use super::numeric;
use super::output_scalar;
use super::runtime_dimensions;
use crate::graph::Arity;
use crate::graph::OpDescriptor;
use crate::graph::RuleError;
use crate::graph::op::EvalDyn;
use crate::graph::rule::Installer;
use crate::kernel::call::CallCtx;
use crate::kernel::call::OpError;
use crate::kernel::native_objects::NativeObjects;
use crate::kernel::signature::ArityInfo;
use crate::kernel::signature::DimExpr;
use crate::kernel::signature::Signature;
use crate::kernel::signature::Stability;
use crate::kernel::signature::TypeExpr;
use crate::kernel::signature::TypeParam;
use crate::kernel::signature::TypeParamKind;
use crate::kernel::value::DynValue;
use crate::kernel::value::Value;
use crate::kernels::autodiff::Tape;
use crate::kernels::autodiff::Var;
use num_traits::ToPrimitive;
use std::sync::Arc;

const TAPE: &str = "rssn/native/ad/tape/v1";
const VAR: &str = "rssn/native/ad/var/v1";
const GRAD: &str = "rssn/native/ad/gradient/v1";
const MAX_NODES: usize = 65_536;
struct TapeObject {
    tape: Tape,
    limit: usize,
}
#[derive(Clone)]
struct AdVariable {
    value: Var,
    tape: Option<Value>,
    limit: usize,
}
struct Gradient {
    values: Vec<f64>,
    origin: Var,
}
fn error(message: &str) -> OpError {
    OpError::Domain(message.into())
}
fn objects<'a>(ctx: &'a CallCtx<'_>) -> Result<&'a NativeObjects, OpError> {
    ctx.native_objects
        .as_ref()
        .ok_or_else(|| OpError::Host("native objects require a session call context".into()))
}
fn canonical(dim: &DimExpr) -> Result<DimExpr, OpError> {
    runtime_dimensions::expression(&runtime_dimensions::canonical(dim)?)
}
fn dimension(value: &Value) -> Result<DimExpr, OpError> {
    canonical(&match value {
        | Value::Dyn(d) => d.dim.clone(),
        | Value::Quantity(q) => q.unit.dimension.clone(),
        | _ => DimExpr::Dimless,
    })
}
fn result(
    value: f64,
    dim: DimExpr,
) -> Result<Value, OpError> {
    let value = output_scalar(value)?;
    if dim == DimExpr::Dimless {
        Ok(value)
    } else {
        Ok(Value::Dyn(DynValue {
            dim,
            value: Box::new(value),
        }))
    }
}
fn read(
    ctx: &CallCtx<'_>,
    handle: &Value,
) -> Result<(AdVariable, DimExpr), OpError> {
    let bank = objects(ctx)?;
    let value = bank.with::<AdVariable, _>(handle, VAR, |v| Ok(v.clone()))?;
    let dim = bank
        .dimensions(handle)?
        .into_iter()
        .next()
        .ok_or_else(|| error("AD variable has no dimension"))?;
    Ok((value, dim))
}
fn store(
    ctx: &CallCtx<'_>,
    value: AdVariable,
    dim: DimExpr,
) -> Result<Value, OpError> {
    if !value.value.val.is_finite() || value.value.node_count() > value.limit {
        return Err(error("AD domain or tape budget exceeded"));
    }
    let parents = value.tape.clone().into_iter().collect::<Vec<_>>();
    let mut domains = vec![dim];
    domains[0] = canonical(&domains[0])?;
    objects(ctx)?.insert_dependent(VAR, domains, value, 256, &parents)
}
// A native elementary operation appends exactly one node. Roll it back if
// validation or handle retention fails, so an error does not consume tape space.
fn store_new(
    ctx: &CallCtx<'_>,
    value: AdVariable,
    dim: DimExpr,
) -> Result<Value, OpError> {
    let rollback = value.value.clone();
    let result = store(ctx, value, dim);
    if result.is_err() {
        rollback.discard_newest();
    }
    result
}
/// Type of a checked native reverse AD variable; its physical domain is stored with the handle.
pub(in crate::rules) fn ad_var_type(_: DimExpr) -> TypeExpr {
    TypeExpr::Opaque(VAR.into())
}
/// Retain a callback-local AD variable. Callers hold a scoped object guard and cap tape nodes.
pub(in crate::rules) fn encode_ad_var(
    ctx: &CallCtx<'_>,
    value: Var,
    dim: DimExpr,
) -> Result<Value, OpError> {
    store(
        ctx,
        AdVariable {
            value,
            tape: None,
            limit: MAX_NODES,
        },
        dim,
    )
}
/// Recover an AD variable only from the current session and exact native kind.
pub(in crate::rules) fn decode_ad_var(
    ctx: &CallCtx<'_>,
    handle: &Value,
) -> Result<(Var, DimExpr), OpError> {
    let (v, d) = read(ctx, handle)?;
    Ok((v.value, d))
}
#[derive(Clone, Copy)]
enum Method {
    TapeNew,
    TapeLen,
    TapeEmpty,
    Variable,
    Primal,
    Nodes,
    Same,
    Backward,
    Grad,
    Sin,
    Cos,
    Exp,
    Ln,
    Sqrt,
    Pow,
    Tanh,
    Add,
    Sub,
    Mul,
    Div,
    Release,
    Usage,
}
struct Adapter(Method);
impl EvalDyn for Adapter {
    fn signature(&self) -> Signature {
        use Method::Add;
        use Method::Backward;
        use Method::Cos;
        use Method::Div;
        use Method::Exp;
        use Method::Grad;
        use Method::Ln;
        use Method::Mul;
        use Method::Nodes;
        use Method::Pow;
        use Method::Primal;
        use Method::Release;
        use Method::Same;
        use Method::Sin;
        use Method::Sqrt;
        use Method::Sub;
        use Method::Tanh;
        use Method::TapeEmpty;
        use Method::TapeLen;
        use Method::TapeNew;
        use Method::Usage;
        use Method::Variable;
        let tape = || TypeExpr::Opaque(TAPE.into());
        let var = || ad_var_type(DimExpr::Dimless);
        let grad = || TypeExpr::Opaque(GRAD.into());
        let (params, ret) = match self.0 {
            | TapeNew => (vec![TypeExpr::Int], tape()),
            | TapeLen => (vec![tape()], TypeExpr::Int),
            | TapeEmpty => (vec![tape()], TypeExpr::Bool),
            | Variable => {
                (
                    vec![tape(), TypeExpr::Float.with_dim(DimExpr::Var("D".into()))],
                    var(),
                )
            },
            | Primal => {
                (
                    vec![var()],
                    TypeExpr::Dyn(crate::kernel::signature::RepExpr {
                        repr: "term".into(),
                        dim: None,
                    }),
                )
            },
            | Nodes => (vec![var()], TypeExpr::Int),
            | Same => (vec![var(), var()], TypeExpr::Bool),
            | Backward => (vec![var()], grad()),
            | Grad => {
                (
                    vec![var(), grad()],
                    TypeExpr::Dyn(crate::kernel::signature::RepExpr {
                        repr: "term".into(),
                        dim: None,
                    }),
                )
            },
            | Sin | Cos | Exp | Ln | Sqrt | Tanh => (vec![var()], var()),
            | Pow => (vec![var(), TypeExpr::Int], var()),
            | Add | Sub | Mul | Div => (vec![var(), var()], var()),
            | Release => {
                (
                    vec![TypeExpr::Dyn(crate::kernel::signature::RepExpr {
                        repr: "term".into(),
                        dim: None,
                    })],
                    TypeExpr::Bool,
                )
            },
            | Usage => {
                (
                    vec![],
                    TypeExpr::Record(vec![
                        ("objects".into(), TypeExpr::Int),
                        ("bytes".into(), TypeExpr::Int),
                    ]),
                )
            },
        };
        let mut signature = Signature::new(params, ret);
        signature.effects = crate::kernel::signature::EffectRow {
            effects: vec![crate::kernel::signature::Effect::Stateful],
        };
        if matches!(self.0, Variable) {
            signature.type_params.push(TypeParam {
                name: "D".into(),
                kind: TypeParamKind::Dim,
                bounds: vec![],
            });
        }
        signature
    }

    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        use Method::Add;
        use Method::Backward;
        use Method::Cos;
        use Method::Div;
        use Method::Exp;
        use Method::Grad;
        use Method::Ln;
        use Method::Mul;
        use Method::Nodes;
        use Method::Pow;
        use Method::Primal;
        use Method::Release;
        use Method::Same;
        use Method::Sin;
        use Method::Sqrt;
        use Method::Sub;
        use Method::Tanh;
        use Method::TapeEmpty;
        use Method::TapeLen;
        use Method::TapeNew;
        use Method::Usage;
        use Method::Variable;
        let count = self.signature().params.len();
        if args.len() != count {
            return Err(OpError::Arity {
                expected: ArityInfo::fixed(count as u8),
                got: args.len(),
            });
        }
        let bank = objects(ctx)?;
        match self.0 {
            | TapeNew => {
                let limit = match &args[0] {
                    | Value::Int(n) => n.to_usize().filter(|n| *n > 0 && *n <= MAX_NODES),
                    | _ => None,
                }
                .ok_or_else(|| error("AD tape node limit must be 1..=65536"))?;
                bank.insert(
                    TAPE,
                    vec![],
                    TapeObject {
                        tape: Tape::new(),
                        limit,
                    },
                    limit * 128,
                )
            },
            | TapeLen | TapeEmpty => {
                bank.with::<TapeObject, _>(&args[0], TAPE, |t| {
                    Ok(if matches!(self.0, TapeLen) {
                        Value::int(t.tape.len())
                    } else {
                        Value::Bool(t.tape.is_empty())
                    })
                })
            },
            | Variable => {
                let val = numeric(&args[1], 1)?;
                let dim = dimension(&args[1])?;
                let variable = bank.with::<TapeObject, _>(&args[0], TAPE, |t| {
                    if t.tape.len() >= t.limit {
                        return Err(error("AD tape node limit exceeded"));
                    }
                    Ok(AdVariable {
                        value: t.tape.var(val),
                        tape: Some(args[0].clone()),
                        limit: t.limit,
                    })
                })?;
                store_new(ctx, variable, dim)
            },
            | Release => {
                let handle = if let Value::Dyn(d) = &args[0] {
                    &d.value
                } else {
                    &args[0]
                };
                bank.release(handle)?;
                Ok(Value::Bool(true))
            },
            | Usage => {
                let (n, b) = bank.usage()?;
                Ok(Value::Record(crate::kernel::value::RecordValue {
                    fields: vec![
                        ("objects".into(), Value::int(n)),
                        ("bytes".into(), Value::int(b)),
                    ],
                }))
            },
            | _ => {
                let (v, d) = read(ctx, &args[0])?;
                match self.0 {
                    | Primal => result(v.value.val, d),
                    | Nodes => Ok(Value::int(v.value.node_count())),
                    | Same => {
                        let (w, _) = read(ctx, &args[1])?;
                        Ok(Value::Bool(v.value.shares_tape(&w.value)))
                    },
                    | Backward => {
                        let values = v.value.backward();
                        if values.iter().any(|n| !n.is_finite()) {
                            return Err(error("non-finite reverse AD gradient"));
                        }
                        let bytes = values.len() * 8 + 128;
                        let parents = v.tape.into_iter().collect::<Vec<_>>();
                        bank.insert_dependent(
                            GRAD,
                            vec![d],
                            Gradient {
                                values,
                                origin: v.value,
                            },
                            bytes,
                            &parents,
                        )
                    },
                    | Grad => {
                        bank.with::<Gradient, _>(&args[1], GRAD, |g| {
                            if !v.value.shares_tape(&g.origin) {
                                return Err(error("gradient belongs to another tape"));
                            }
                            // Later tape nodes cannot affect the earlier output. Native
                            // grad_in returns zero for those later independent variables.
                            let output = bank
                                .dimensions(&args[1])?
                                .into_iter()
                                .next()
                                .ok_or_else(|| error("gradient dimension missing"))?;
                            result(
                                v.value.grad_in(&g.values),
                                canonical(&DimExpr::Mul(vec![
                                    output,
                                    DimExpr::Inverse(Box::new(d)),
                                ]))?,
                            )
                        })
                    },
                    | Add | Sub | Mul | Div => {
                        let (w, e) = read(ctx, &args[1])?;
                        if !v.value.shares_tape(&w.value) {
                            return Err(error("AD operands belong to different tapes"));
                        }
                        if v.value.node_count() >= v.limit {
                            return Err(error("AD tape node limit exceeded"));
                        }
                        let dim = match self.0 {
                            | Add | Sub => {
                                if runtime_dimensions::canonical(&d)?
                                    != runtime_dimensions::canonical(&e)?
                                {
                                    return Err(error("AD sum requires equal dimensions"));
                                }
                                d
                            },
                            | Mul => canonical(&DimExpr::Mul(vec![d, e]))?,
                            | _ => {
                                canonical(&DimExpr::Mul(vec![d, DimExpr::Inverse(Box::new(e))]))?
                            },
                        };
                        let value = match self.0 {
                            | Add => &v.value + &w.value,
                            | Sub => &v.value - &w.value,
                            | Mul => &v.value * &w.value,
                            | _ => &v.value / &w.value,
                        };
                        store_new(ctx, AdVariable { value, ..v }, dim)
                    },
                    | Sin | Cos | Exp | Ln | Sqrt | Pow | Tanh => {
                        if v.value.node_count() >= v.limit {
                            return Err(error("AD tape node limit exceeded"));
                        }
                        let power = if matches!(self.0, Pow) {
                            match &args[1] {
                                | Value::Int(n) => n.to_i32().filter(|n| n.unsigned_abs() <= 1024),
                                | _ => None,
                            }
                            .ok_or_else(|| error("AD integer exponent must be -1024..=1024"))?
                        } else {
                            1
                        };
                        let dim = match self.0 {
                            | Sqrt => {
                                canonical(&d.pow(crate::kernel::signature::Exp::Rat(
                                    crate::kernel::signature::Rational::new(1, 2),
                                )))?
                            },
                            | Pow => {
                                canonical(&d.pow(crate::kernel::signature::Exp::Rat(
                                    crate::kernel::signature::Rational::int(i64::from(power)),
                                )))?
                            },
                            | _ => {
                                if d != DimExpr::Dimless {
                                    return Err(error(
                                        "AD transcendental requires a dimensionless argument",
                                    ));
                                }
                                d
                            },
                        };
                        let value = match self.0 {
                            | Sin => v.value.sin(),
                            | Cos => v.value.cos(),
                            | Exp => v.value.exp(),
                            | Ln => v.value.ln(),
                            | Sqrt => v.value.sqrt(),
                            | Pow => v.value.powf(f64::from(power)),
                            | _ => v.value.tanh(),
                        };
                        store_new(ctx, AdVariable { value, ..v }, dim)
                    },
                    | _ => unreachable!(),
                }
            },
        }
    }
}
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    use Method::Add;
    use Method::Backward;
    use Method::Cos;
    use Method::Div;
    use Method::Exp;
    use Method::Grad;
    use Method::Ln;
    use Method::Mul;
    use Method::Nodes;
    use Method::Pow;
    use Method::Primal;
    use Method::Release;
    use Method::Same;
    use Method::Sin;
    use Method::Sqrt;
    use Method::Sub;
    use Method::Tanh;
    use Method::TapeEmpty;
    use Method::TapeLen;
    use Method::TapeNew;
    use Method::Usage;
    use Method::Variable;
    for (name, path, method) in [
        (
            "num_ad_tape_new",
            "crate::kernels::autodiff::Tape::new",
            TapeNew,
        ),
        (
            "num_ad_tape_len",
            "crate::kernels::autodiff::Tape::len",
            TapeLen,
        ),
        (
            "num_ad_tape_is_empty",
            "crate::kernels::autodiff::Tape::is_empty",
            TapeEmpty,
        ),
        (
            "num_ad_tape_var",
            "crate::kernels::autodiff::Tape::var",
            Variable,
        ),
        ("num_ad_var_value", "", Primal),
        (
            "num_ad_var_node_count",
            "crate::kernels::autodiff::Var::node_count",
            Nodes,
        ),
        (
            "num_ad_var_shares_tape",
            "crate::kernels::autodiff::Var::shares_tape",
            Same,
        ),
        (
            "num_ad_var_backward",
            "crate::kernels::autodiff::Var::backward",
            Backward,
        ),
        (
            "num_ad_var_grad_in",
            "crate::kernels::autodiff::Var::grad_in",
            Grad,
        ),
        ("num_ad_var_sin", "crate::kernels::autodiff::Var::sin", Sin),
        ("num_ad_var_cos", "crate::kernels::autodiff::Var::cos", Cos),
        ("num_ad_var_exp", "crate::kernels::autodiff::Var::exp", Exp),
        ("num_ad_var_ln", "crate::kernels::autodiff::Var::ln", Ln),
        (
            "num_ad_var_sqrt",
            "crate::kernels::autodiff::Var::sqrt",
            Sqrt,
        ),
        ("num_ad_var_pow", "crate::kernels::autodiff::Var::powf", Pow),
        (
            "num_ad_var_tanh",
            "crate::kernels::autodiff::Var::tanh",
            Tanh,
        ),
        ("num_ad_var_add", "", Add),
        ("num_ad_var_sub", "", Sub),
        ("num_ad_var_mul", "", Mul),
        ("num_ad_var_div", "", Div),
        ("num_native_release", "", Release),
        ("num_native_usage", "", Usage),
    ] {
        let adapter = Adapter(method);
        let signature = adapter.signature();
        let mut descriptor=OpDescriptor::new(name,Arity::Fixed(signature.params.len()as u8)).category("numerical_vector").purity(crate::kernel::signature::Purity::Effectful).determinism(crate::kernel::signature::Determinism::EnvDependent).thread_safety(crate::kernel::signature::ThreadSafety::SessionBound).backend(crate::kernel::signature::BackendSupport{interpreter:true,jit:false,ffi:false,gpu:false}).doc("Session-local reverse AD object. Tape constructor reserves 1..65536 nodes. Linear quantity values use coherent units; addition requires matching dimensions, products/quotients and gradients propagate physical powers; transcendental arguments must be dimensionless. Opaque handles are kind checked, invalid across sessions and after release. Release dependent variables/gradients before their tape. Callbacks scope temporary handles. No native pointer/state is serialized.").stability(Stability::Experimental).signature(signature).eval_dyn(Arc::new(adapter));
        if !path.is_empty() {
            descriptor = descriptor.rust_path(path);
        }
        let op = i.op(descriptor)?;
        i.dim(op, crate::graph::dim::DimRule::Container);
    }
    Ok(())
}
