//! Normalized sparse SQP with checked callback packets and a shared call budget.

use super::CallCtx;
use super::DimExpr;
use super::EvalDyn;
use super::Installer;
use super::OpError;
use super::RecordValue;
use super::RuleError;
use super::Signature;
use super::TypeExpr;
use super::Value;
use super::common::CallbackRun;
use super::common::effectful;
use super::common::positive_budget;
use super::common::register;
use super::evolution::output_tensor;
use super::evolution::tensor_type;
use super::evolution::vector;
use super::numeric;
use crate::kernel::signature::ArityInfo;
use crate::kernel::signature::NatExpr;
use crate::kernel::value::TensorValue;
use crate::kernels::milp::sqp::NlpOracle;
use crate::kernels::milp::sqp::SparseMat;
use crate::kernels::milp::sqp::SparseSym;
use crate::kernels::milp::sqp::SqpOptions;
use crate::kernels::milp::sqp::solve_sqp;
use num_traits::ToPrimitive;
use std::cell::Cell;
use std::cell::RefCell;

const MAX_SIDE: usize = 128;
const MAX_TRIPLETS: usize = 16_384;

fn floats() -> TypeExpr {
    tensor_type(1, DimExpr::Dimless)
}

fn integers() -> TypeExpr {
    TypeExpr::Tensor {
        rank: NatExpr::Lit(1),
        shape: None,
        elem: Box::new(TypeExpr::Int),
    }
}

fn record_type(fields: &[(&str, TypeExpr)]) -> TypeExpr {
    TypeExpr::Record(
        fields
            .iter()
            .map(|(name, ty)| ((*name).into(), ty.clone()))
            .collect(),
    )
}

fn typed_record<const N: usize>(fields: [(&str, Value); N]) -> Value {
    Value::Record(RecordValue {
        fields: fields
            .into_iter()
            .map(|(name, value)| (name.into(), value))
            .collect(),
    })
}

fn sparse_type() -> TypeExpr {
    record_type(&[
        ("rows", integers()),
        ("columns", integers()),
        ("values", floats()),
    ])
}

fn shape_type() -> TypeExpr {
    record_type(&[
        ("equalities", TypeExpr::Int),
        ("inequalities", TypeExpr::Int),
    ])
}

fn evaluation_type() -> TypeExpr {
    TypeExpr::Function(
        vec![floats()],
        Box::new(record_type(&[
            ("objective", TypeExpr::Float),
            ("gradient", floats()),
            ("constraints", floats()),
            ("jacobian", sparse_type()),
        ])),
    )
}

fn hessian_type() -> TypeExpr {
    TypeExpr::Function(vec![floats(), floats(), floats()], Box::new(sparse_type()))
}

fn checked_record<'a>(
    value: &'a Value,
    names: &[&str],
) -> Result<&'a RecordValue, OpError> {
    crate::rules::numerical_abi::native_object_api::bounded(value)?;
    let Value::Record(record) = value else {
        return Err(OpError::Domain(
            "SQP callback must return the declared record".into(),
        ));
    };
    if record.fields.len() != names.len() || names.iter().any(|name| record.get(name).is_none()) {
        return Err(OpError::Domain(
            "SQP record fields do not match the declared fields".into(),
        ));
    }
    Ok(record)
}

fn field<'a>(
    record: &'a RecordValue,
    name: &str,
) -> Result<&'a Value, OpError> {
    record
        .get(name)
        .ok_or_else(|| OpError::Domain(format!("missing SQP field {name}")))
}

fn integer(value: &Value) -> Result<usize, OpError> {
    let Value::Int(value) = value else {
        return Err(OpError::Domain(
            "SQP index or count must be an integer".into(),
        ));
    };
    value
        .to_usize()
        .ok_or_else(|| OpError::Domain("SQP index or count is out of range".into()))
}

fn index_vector(value: &Value) -> Result<Vec<usize>, OpError> {
    let Value::Tensor(vector) = value else {
        return Err(OpError::Domain(
            "SQP sparse indices must be a vector".into(),
        ));
    };
    if vector.shape != [vector.data.len()] || vector.data.len() > MAX_TRIPLETS {
        return Err(OpError::Domain(
            "SQP sparse indices exceed their shape or size budget".into(),
        ));
    }
    vector.data.iter().map(integer).collect()
}

fn sparse(
    value: &Value,
    rows_n: usize,
    cols_n: usize,
    symmetric: bool,
    ctx: &CallCtx<'_>,
) -> Result<SparseMat, OpError> {
    let record = checked_record(value, &["rows", "columns", "values"])?;
    let rows = index_vector(field(record, "rows")?)?;
    let cols = index_vector(field(record, "columns")?)?;
    let vals = vector(
        field(record, "values")?,
        &DimExpr::Dimless,
        Some(rows.len()),
        0,
        ctx,
    )?;
    if cols.len() != rows.len()
        || rows
            .iter()
            .zip(&cols)
            .any(|(&r, &c)| r >= rows_n || c >= cols_n || (symmetric && r > c))
    {
        return Err(OpError::Domain(
            "SQP sparse triplet lengths or bounds do not match".into(),
        ));
    }
    Ok(SparseMat {
        rows_n,
        cols_n,
        rows,
        cols,
        vals,
    })
}

fn wire_sparse(
    rows: &[usize],
    cols: &[usize],
    vals: Vec<f64>,
) -> Result<Value, OpError> {
    let indices = |values: &[usize]| {
        Value::Tensor(TensorValue {
            shape: vec![values.len()],
            data: values.iter().copied().map(Value::int).collect(),
            elem: TypeExpr::Int,
        })
    };
    Ok(typed_record([
        ("rows", indices(rows)),
        ("columns", indices(cols)),
        (
            "values",
            output_tensor(vec![vals.len()], vals, DimExpr::Dimless)?,
        ),
    ]))
}

#[derive(Clone, Copy)]
struct Shape {
    eq: usize,
    ineq: usize,
}

impl Shape {
    fn decode(value: &Value) -> Result<Self, OpError> {
        let value = checked_record(value, &["equalities", "inequalities"])?;
        let eq = integer(field(value, "equalities")?)?;
        let ineq = integer(field(value, "inequalities")?)?;
        if eq.saturating_add(ineq) > MAX_SIDE {
            return Err(OpError::Domain(
                "SQP constraint count exceeds its budget".into(),
            ));
        }
        Ok(Self { eq, ineq })
    }

    const fn total(self) -> usize {
        self.eq + self.ineq
    }
}

struct Evaluation {
    objective: f64,
    gradient: Vec<f64>,
    constraints: Vec<f64>,
    jacobian: SparseMat,
}

impl Evaluation {
    fn decode(
        value: &Value,
        shape: Shape,
        n: usize,
        ctx: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        let record = checked_record(value, &["objective", "gradient", "constraints", "jacobian"])?;
        let objective = field(record, "objective")?;
        super::check_dimension(ctx, objective, &DimExpr::Dimless, 0)?;
        Ok(Self {
            objective: numeric(objective, 0)?,
            gradient: vector(
                field(record, "gradient")?,
                &DimExpr::Dimless,
                Some(n),
                0,
                ctx,
            )?,
            constraints: vector(
                field(record, "constraints")?,
                &DimExpr::Dimless,
                Some(shape.total()),
                0,
                ctx,
            )?,
            jacobian: sparse(field(record, "jacobian")?, shape.total(), n, false, ctx)?,
        })
    }
}

struct Oracle<'a, 'ctx> {
    context: &'a CallCtx<'ctx>,
    shape: Shape,
    evaluator: CallbackRun<'a, 'ctx>,
    hessian: CallbackRun<'a, 'ctx>,
    remaining: Cell<usize>,
    failure: RefCell<Option<OpError>>,
    last: RefCell<Option<(Vec<f64>, Evaluation)>>,
}

impl Oracle<'_, '_> {
    fn permit(&self) -> bool {
        if self.failure.borrow().is_some() || self.evaluator.failed() || self.hessian.failed() {
            return false;
        }
        let remaining = self.remaining.get();
        if remaining == 0 {
            *self.failure.borrow_mut() =
                Some(OpError::Domain("SQP callback budget exhausted".into()));
            return false;
        }
        if let Err(error) = crate::cancellation::checkpoint() {
            *self.failure.borrow_mut() = Some(error.into());
            return false;
        }
        self.remaining.set(remaining - 1);
        true
    }

    fn with_point<T>(
        &self,
        x: &[f64],
        get: impl FnOnce(&Evaluation) -> T,
    ) -> Option<T> {
        {
            let cached = self.last.borrow();
            if let Some((point, value)) = cached.as_ref()
                && point.len() == x.len()
                && point.iter().zip(x).all(|(a, b)| a.to_bits() == b.to_bits())
            {
                return Some(get(value));
            }
        }
        if !self.permit() {
            return None;
        }
        let arg = output_tensor(vec![x.len()], x.to_vec(), DimExpr::Dimless).ok()?;
        let sample = self.evaluator.evaluate(&[arg], |value| {
            Evaluation::decode(&value, self.shape, x.len(), self.context)
        })?;
        let output = get(&sample);
        *self.last.borrow_mut() = Some((x.to_vec(), sample));
        Some(output)
    }
}

impl NlpOracle for Oracle<'_, '_> {
    fn objective(
        &self,
        x: &[f64],
    ) -> f64 {
        self.with_point(x, |value| value.objective)
            .unwrap_or(f64::NAN)
    }

    fn gradient(
        &self,
        x: &[f64],
    ) -> Vec<f64> {
        self.with_point(x, |value| value.gradient.clone())
            .unwrap_or_default()
    }

    fn hessian(
        &self,
        x: &[f64],
        eq: &[f64],
        ineq: &[f64],
    ) -> SparseSym {
        let output = || -> Option<SparseSym> {
            if !self.permit() {
                return None;
            }
            let a = [x, eq, ineq]
                .map(|vector| output_tensor(vec![vector.len()], vector.to_vec(), DimExpr::Dimless));
            let [x_value, eq_value, ineq_value] = a;
            self.hessian.evaluate(
                &[x_value.ok()?, eq_value.ok()?, ineq_value.ok()?],
                |value| {
                    let m = sparse(&value, x.len(), x.len(), true, self.context)?;
                    Ok(SparseSym {
                        n: x.len(),
                        rows: m.rows,
                        cols: m.cols,
                        vals: m.vals,
                    })
                },
            )
        };
        output().unwrap_or_else(|| {
            SparseSym {
                n: 0,
                rows: vec![],
                cols: vec![],
                vals: vec![],
            }
        })
    }

    fn constraints(
        &self,
        x: &[f64],
    ) -> Vec<f64> {
        self.with_point(x, |value| value.constraints.clone())
            .unwrap_or_default()
    }

    fn jacobian(
        &self,
        x: &[f64],
    ) -> SparseMat {
        self.with_point(x, |value| value.jacobian.clone())
            .unwrap_or_else(|| {
                SparseMat {
                    rows_n: 0,
                    cols_n: 0,
                    rows: vec![],
                    cols: vec![],
                    vals: vec![],
                }
            })
    }

    fn n_eq(&self) -> usize {
        self.shape.eq
    }

    fn n_ineq(&self) -> usize {
        self.shape.ineq
    }
}

fn point(
    value: &Value,
    ctx: &CallCtx<'_>,
) -> Result<Vec<f64>, OpError> {
    crate::rules::numerical_abi::native_object_api::bounded(value)?;
    let x = vector(value, &DimExpr::Dimless, None, 0, ctx)?;
    if x.is_empty() || x.len() > MAX_SIDE {
        return Err(OpError::Domain("SQP point exceeds its side budget".into()));
    }
    Ok(x)
}

struct SolveSqp;
impl EvalDyn for SolveSqp {
    fn signature(&self) -> Signature {
        effectful(Signature::new(
            vec![
                floats(),
                shape_type(),
                TypeExpr::Float,
                TypeExpr::Float,
                TypeExpr::Int,
                TypeExpr::Int,
                evaluation_type(),
                hessian_type(),
            ],
            record_type(&[
                ("x", floats()),
                ("objective", TypeExpr::Float),
                ("equality_multipliers", floats()),
                ("inequality_multipliers", floats()),
                ("iterations", TypeExpr::Int),
                ("evaluations", TypeExpr::Int),
                ("kkt_residual", TypeExpr::Float),
                ("converged", TypeExpr::Bool),
            ]),
        ))
    }

    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        if args.len() != 8 {
            return Err(OpError::Arity {
                expected: ArityInfo::fixed(8),
                got: args.len(),
            });
        }
        let x = point(&args[0], ctx)?;
        let shape = Shape::decode(&args[1])?;
        super::check_dimension(ctx, &args[2], &DimExpr::Dimless, 2)?;
        super::check_dimension(ctx, &args[3], &DimExpr::Dimless, 3)?;
        let tol = numeric(&args[2], 2)?;
        let trust_radius = numeric(&args[3], 3)?;
        let max_iter = positive_budget(&args[4], 4, "SQP iteration limit")?;
        let budget = positive_budget(&args[5], 5, "SQP callback budget")?;
        if tol <= 0.0 || trust_radius <= 0.0 || max_iter > 256 || budget > 10_000 {
            return Err(OpError::Domain(
                "SQP needs positive tolerances within its iteration and callback budgets".into(),
            ));
        }
        let (evaluator, _) = CallbackRun::new(ctx, &args[6], 6, 1, evaluation_type())?;
        let (hessian, _) = CallbackRun::new(ctx, &args[7], 7, 3, hessian_type())?;
        let oracle = Oracle {
            context: ctx,
            shape,
            evaluator,
            hessian,
            remaining: Cell::new(budget),
            failure: RefCell::new(None),
            last: RefCell::new(None),
        };
        let result = solve_sqp(
            &oracle,
            &x,
            &SqpOptions {
                max_iter,
                tol,
                trust_radius,
            },
        );
        let evaluations = budget - oracle.remaining.get();
        if let Some(error) = oracle.failure.into_inner() {
            return Err(error);
        }
        oracle.evaluator.finish()?;
        oracle.hessian.finish()?;
        let result = result.map_err(OpError::Domain)?;
        let dense = |values: Vec<f64>| output_tensor(vec![values.len()], values, DimExpr::Dimless);
        Ok(typed_record([
            ("x", dense(result.x)?),
            ("objective", Value::Float(result.objective)),
            ("equality_multipliers", dense(result.mult_eq)?),
            ("inequality_multipliers", dense(result.mult_ineq)?),
            ("iterations", Value::int(result.iterations)),
            ("evaluations", Value::int(evaluations)),
            ("kkt_residual", Value::Float(result.kkt_residual)),
            ("converged", Value::Bool(result.kkt_residual <= tol)),
        ]))
    }
}

#[derive(Clone, Copy)]
enum Projected {
    Objective,
    Gradient,
    Constraints,
    Jacobian,
}

impl EvalDyn for Projected {
    fn signature(&self) -> Signature {
        let result = match self {
            | Self::Objective => TypeExpr::Float,
            | Self::Gradient | Self::Constraints => floats(),
            | Self::Jacobian => sparse_type(),
        };
        effectful(Signature::new(
            vec![floats(), shape_type(), evaluation_type()],
            result,
        ))
    }

    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        if args.len() != 3 {
            return Err(OpError::Arity {
                expected: ArityInfo::fixed(3),
                got: args.len(),
            });
        }
        let x = point(&args[0], ctx)?;
        let shape = Shape::decode(&args[1])?;
        let (callback, _) = CallbackRun::new(ctx, &args[2], 2, 1, evaluation_type())?;
        crate::cancellation::checkpoint()?;
        let side = x.len();
        let argument = output_tensor(vec![side], x, DimExpr::Dimless)?;
        let sample = callback.evaluate(&[argument], |value| {
            Evaluation::decode(&value, shape, side, ctx)
        });
        callback.finish()?;
        let sample =
            sample.ok_or_else(|| OpError::Domain("SQP callback produced no evaluation".into()))?;
        match self {
            | Self::Objective => Ok(Value::Float(sample.objective)),
            | Self::Gradient => {
                output_tensor(
                    vec![sample.gradient.len()],
                    sample.gradient,
                    DimExpr::Dimless,
                )
            },
            | Self::Constraints => {
                output_tensor(
                    vec![sample.constraints.len()],
                    sample.constraints,
                    DimExpr::Dimless,
                )
            },
            | Self::Jacobian => {
                wire_sparse(
                    &sample.jacobian.rows,
                    &sample.jacobian.cols,
                    sample.jacobian.vals,
                )
            },
        }
    }
}

struct HessianProjection;
impl EvalDyn for HessianProjection {
    fn signature(&self) -> Signature {
        effectful(Signature::new(
            vec![floats(), floats(), floats(), hessian_type()],
            sparse_type(),
        ))
    }

    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        if args.len() != 4 {
            return Err(OpError::Arity {
                expected: ArityInfo::fixed(4),
                got: args.len(),
            });
        }
        let x = point(&args[0], ctx)?;
        let eq = vector(&args[1], &DimExpr::Dimless, None, 1, ctx)?;
        let ineq = vector(&args[2], &DimExpr::Dimless, None, 2, ctx)?;
        if eq.len().saturating_add(ineq.len()) > MAX_SIDE {
            return Err(OpError::Domain(
                "SQP multiplier count exceeds its budget".into(),
            ));
        }
        let (callback, _) = CallbackRun::new(ctx, &args[3], 3, 3, hessian_type())?;
        crate::cancellation::checkpoint()?;
        let size = x.len();
        let dense = |values: Vec<f64>| output_tensor(vec![values.len()], values, DimExpr::Dimless);
        let evaluated = callback.evaluate(&[dense(x)?, dense(eq)?, dense(ineq)?], |value| {
            sparse(&value, size, size, true, ctx)
        });
        callback.finish()?;
        let matrix =
            evaluated.ok_or_else(|| OpError::Domain("SQP callback produced no Hessian".into()))?;
        wire_sparse(&matrix.rows, &matrix.cols, matrix.vals)
    }
}

struct CountProjection {
    equality: bool,
}
impl EvalDyn for CountProjection {
    fn signature(&self) -> Signature {
        Signature::new(vec![shape_type()], TypeExpr::Int)
    }

    fn eval(
        &self,
        args: &[Value],
        _: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        if args.len() != 1 {
            return Err(OpError::Arity {
                expected: ArityInfo::fixed(1),
                got: args.len(),
            });
        }
        let shape = Shape::decode(&args[0])?;
        Ok(Value::int(if self.equality {
            shape.eq
        } else {
            shape.ineq
        }))
    }
}

fn install_count(
    installer: &mut Installer<'_>,
    name: &'static str,
    path: &'static str,
    doc: &'static str,
    equality: bool,
) -> Result<(), RuleError> {
    let evaluator = CountProjection { equality };
    let op = installer.op(super::OpDescriptor::new(name, super::Arity::Fixed(1))
        .rust_path(path)
        .doc(doc)
        .category("numerical_callback")
        .stability(super::Stability::Experimental)
        .thread_safety(super::ThreadSafety::Send)
        .signature(evaluator.signature())
        .eval_dyn(super::Arc::new(evaluator)))?;
    installer.dim(op, crate::graph::dim::DimRule::Dimensionless);
    Ok(())
}

pub(super) fn install(installer: &mut Installer<'_>) -> Result<(), RuleError> {
    register(
        installer,
        "num_sqp_solve",
        "crate::kernels::milp::sqp::solve_sqp",
        "Solves a normalized nonlinear program with bounded typed evaluation and Lagrangian-Hessian callbacks. The evaluation callback receives each consecutive point once. Equality rows precede inequalities g(x) <= 0. Returns KKT residual, convergence, multipliers and callback counts.",
        8,
        SolveSqp,
    )?;
    for (name, path, docs, selection) in [
        (
            "num_sqp_oracle_objective",
            "crate::kernels::milp::sqp::NlpOracle::objective",
            "Reads the objective from the checked SQP evaluation callback.",
            Projected::Objective,
        ),
        (
            "num_sqp_oracle_gradient",
            "crate::kernels::milp::sqp::NlpOracle::gradient",
            "Reads the objective gradient from the checked SQP evaluation callback.",
            Projected::Gradient,
        ),
        (
            "num_sqp_oracle_constraints",
            "crate::kernels::milp::sqp::NlpOracle::constraints",
            "Reads equality-first values from the checked SQP evaluation callback.",
            Projected::Constraints,
        ),
        (
            "num_sqp_oracle_jacobian",
            "crate::kernels::milp::sqp::NlpOracle::jacobian",
            "Reads the zero-based sparse constraint Jacobian from the checked SQP evaluation callback.",
            Projected::Jacobian,
        ),
    ] {
        register(installer, name, path, docs, 3, selection)?;
    }
    register(
        installer,
        "num_sqp_oracle_hessian",
        "crate::kernels::milp::sqp::NlpOracle::hessian",
        "Evaluates the callback for the sparse upper-triangular Lagrangian Hessian.",
        4,
        HessianProjection,
    )?;
    install_count(
        installer,
        "num_sqp_oracle_equality_count",
        "crate::kernels::milp::sqp::NlpOracle::n_eq",
        "Returns the declared equality row count of the normalized SQP oracle.",
        true,
    )?;
    install_count(
        installer,
        "num_sqp_oracle_inequality_count",
        "crate::kernels::milp::sqp::NlpOracle::n_ineq",
        "Returns the declared inequality row count of the normalized SQP oracle.",
        false,
    )
}
