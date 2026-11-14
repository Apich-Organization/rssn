//! Serial, dimensionless adapters for advanced numeric callback algorithms.
//!
//! Every call has an explicit host-evaluation budget shared by all callbacks.
//! Numerical work stays in `kernels`; this module only checks and marshals ABI values.
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
use super::check_dimension;
use super::common::CallbackRun;
use super::common::effectful;
use super::common::positive_budget;
use super::common::register;
use super::evolution::output_tensor;
use super::evolution::tensor_type;
use super::evolution::vector;
use super::evolution::vector_dimension;
use super::numeric;
use super::scalar_dimension;
use crate::kernel::signature::ArityInfo;
use crate::kernel::value::TensorValue;
use crate::kernels::optim;
use crate::kernels::quadrature;
use crate::kernels::rootfind;
use std::cell::Cell;

struct Budget {
    maximum: usize,
    used: Cell<usize>,
    failed: Cell<bool>,
    exhausted: Cell<bool>,
}
impl Budget {
    const fn new(maximum: usize) -> Self {
        Self {
            maximum,
            used: Cell::new(0),
            failed: Cell::new(false),
            exhausted: Cell::new(false),
        }
    }

    fn evaluate<T>(
        &self,
        callback: &CallbackRun<'_, '_>,
        args: &[Value],
        decode: impl FnOnce(Value) -> Result<T, OpError>,
    ) -> Option<T> {
        if self.failed.get() {
            return None;
        }
        if self.used.get() >= self.maximum {
            self.failed.set(true);
            self.exhausted.set(true);
            return None;
        }
        self.used.set(self.used.get() + 1);
        let result = callback.evaluate(args, decode);
        if result.is_none() {
            self.failed.set(true);
        }
        result
    }

    fn finish(&self) -> Result<(), OpError> {
        if self.exhausted.get() {
            Err(OpError::NotReduced(
                "numerical callback evaluation budget exhausted".into(),
            ))
        } else {
            Ok(())
        }
    }
}
fn arity(
    args: &[Value],
    count: u8,
) -> Result<(), OpError> {
    if args.len() == usize::from(count) {
        Ok(())
    } else {
        Err(OpError::Arity {
            expected: ArityInfo::fixed(count),
            got: args.len(),
        })
    }
}
fn scalar(
    args: &[Value],
    index: usize,
    ctx: &CallCtx<'_>,
) -> Result<f64, OpError> {
    check_dimension(ctx, &args[index], &DimExpr::Dimless, index)?;
    numeric(&args[index], index)
}
fn tolerance(
    args: &[Value],
    index: usize,
    ctx: &CallCtx<'_>,
) -> Result<f64, OpError> {
    let tol = scalar(args, index, ctx)?;
    if tol > 0.0 {
        Ok(tol)
    } else {
        Err(OpError::Domain(
            "tolerance must be finite and positive".into(),
        ))
    }
}
fn dimensionless(d: &DimExpr) -> Result<(), OpError> {
    if super::canonical(d)? == super::canonical(&DimExpr::Dimless)? {
        Ok(())
    } else {
        Err(OpError::Domain(
            "this raw numerical callback interface requires dimensionless values".into(),
        ))
    }
}
fn scalar_callback(signature: &Signature) -> Result<(), OpError> {
    dimensionless(&scalar_dimension(&signature.params[0])?)?;
    dimensionless(&scalar_dimension(&signature.ret)?)
}
fn vector_callback(
    signature: &Signature,
    length: usize,
    vector_result: bool,
) -> Result<(), OpError> {
    dimensionless(&vector_dimension(&signature.params[0], Some(length))?)?;
    if vector_result {
        dimensionless(&vector_dimension(&signature.ret, Some(length))?)
    } else {
        dimensionless(&scalar_dimension(&signature.ret)?)
    }
}
fn vec_value(x: &[f64]) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![x.len()],
        data: x.iter().copied().map(Value::Float).collect(),
        elem: TypeExpr::Float,
    })
}
fn record(fields: Vec<(&str, Value)>) -> Value {
    Value::Record(RecordValue {
        fields: fields
            .into_iter()
            .map(|(name, value)| (name.into(), value))
            .collect(),
    })
}
fn function(
    input: TypeExpr,
    output: TypeExpr,
) -> TypeExpr {
    TypeExpr::Function(vec![input], Box::new(output))
}
fn vector_type() -> TypeExpr {
    tensor_type(1, DimExpr::Dimless)
}
fn vector_result() -> TypeExpr {
    TypeExpr::Record(vec![
        ("x".into(), vector_type()),
        ("evaluations".into(), TypeExpr::Int),
    ])
}
fn root_error(error: rootfind::RootError) -> OpError {
    match error {
        | rootfind::RootError::NoConvergence => {
            OpError::NotReduced(
                "nonlinear solver exhausted its iteration budget or stagnated".into(),
            )
        },
        | other => OpError::Domain(format!("nonlinear solver failed: {other:?}")),
    }
}

#[derive(Clone, Copy)]
enum Quad {
    GaussLegendre,
    ClenshawCurtis,
    Romberg,
    TanhSinh,
    RightInfinite,
    LeftInfinite,
    Infinite,
}
impl Quad {
    const fn fixed(self) -> bool {
        matches!(self, Self::GaussLegendre | Self::ClenshawCurtis)
    }

    const fn finite(self) -> bool {
        matches!(
            self,
            Self::GaussLegendre | Self::ClenshawCurtis | Self::Romberg | Self::TanhSinh
        )
    }

    const fn count(self) -> u8 {
        if self.fixed() {
            5
        } else if self.finite() {
            6
        } else if matches!(self, Self::Infinite) {
            3
        } else {
            4
        }
    }
}
impl EvalDyn for Quad {
    fn signature(&self) -> Signature {
        let mut params = if self.fixed() {
            vec![
                TypeExpr::Float,
                TypeExpr::Float,
                TypeExpr::Int,
                TypeExpr::Int,
            ]
        } else if self.finite() {
            vec![
                TypeExpr::Float,
                TypeExpr::Float,
                TypeExpr::Float,
                TypeExpr::Int,
                TypeExpr::Int,
            ]
        } else if matches!(self, Self::Infinite) {
            vec![TypeExpr::Float, TypeExpr::Int]
        } else {
            vec![TypeExpr::Float, TypeExpr::Float, TypeExpr::Int]
        };
        params.push(function(TypeExpr::Float, TypeExpr::Float));
        let mut fields = vec![
            ("value".into(), TypeExpr::Float),
            ("evaluations".into(), TypeExpr::Int),
        ];
        if !self.fixed() {
            fields.extend([
                ("error".into(), TypeExpr::Float),
                ("converged".into(), TypeExpr::Bool),
            ]);
        }
        effectful(Signature::new(params, TypeExpr::Record(fields)))
    }

    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        arity(args, self.count())?;
        let last = args.len() - 1;
        let (callback, signature) = CallbackRun::new(
            ctx,
            &args[last],
            last,
            1,
            self.signature().params[last].clone(),
        )?;
        scalar_callback(&signature)?;
        let budget = Budget::new(positive_budget(
            &args[last - 1],
            last - 1,
            "quadrature evaluation budget",
        )?);
        let evaluate = |x: f64| {
            budget
                .evaluate(&callback, &[Value::Float(x)], |value| numeric(&value, last))
                .unwrap_or(f64::NAN)
        };
        let result = if self.fixed() {
            let a = scalar(args, 0, ctx)?;
            let b = scalar(args, 1, ctx)?;
            if b < a {
                return Err(OpError::Domain("quadrature requires a <= b".into()));
            }
            let order = positive_budget(&args[2], 2, "quadrature order")?;
            if order > 4096 {
                return Err(OpError::Domain(
                    "fixed quadrature order must be at most 4096".into(),
                ));
            }
            let value = if crate::kernels::exact_eq(a, b) {
                0.0
            } else if matches!(self, Self::GaussLegendre) {
                quadrature::gauss_legendre(evaluate, a, b, order)
            } else {
                quadrature::clenshaw_curtis(evaluate, a, b, order)
            };
            budget.finish()?;
            callback.finish()?;
            if !value.is_finite() {
                return Err(OpError::Domain("quadrature value is not finite".into()));
            }
            return Ok(record(vec![
                ("value", Value::Float(value)),
                ("evaluations", Value::int(budget.used.get())),
            ]));
        } else if self.finite() {
            let a = scalar(args, 0, ctx)?;
            let b = scalar(args, 1, ctx)?;
            if b < a {
                return Err(OpError::Domain("quadrature requires a <= b".into()));
            }
            let tol = tolerance(args, 2, ctx)?;
            let levels = positive_budget(&args[3], 3, "quadrature refinement budget")?;
            if matches!(self, Self::Romberg) {
                if !(2..=24).contains(&levels) {
                    return Err(OpError::Domain(
                        "Romberg requires 2 to 24 refinement levels".into(),
                    ));
                }
                quadrature::romberg(evaluate, a, b, tol, levels)
            } else {
                if levels > 20 {
                    return Err(OpError::Domain(
                        "tanh-sinh supports at most 20 refinement levels".into(),
                    ));
                }
                quadrature::tanh_sinh_with_levels(evaluate, a, b, tol, levels)
            }
        } else if matches!(self, Self::Infinite) {
            quadrature::integrate_infinite(evaluate, tolerance(args, 0, ctx)?)
        } else {
            let bound = scalar(args, 0, ctx)?;
            let tol = tolerance(args, 1, ctx)?;
            if matches!(self, Self::RightInfinite) {
                quadrature::integrate_semi_infinite(evaluate, bound, tol)
            } else {
                quadrature::integrate_semi_infinite(|x| evaluate(-x), -bound, tol)
            }
        };
        budget.finish()?;
        callback.finish()?;
        if !result.value.is_finite() || !result.error.is_finite() {
            return Err(OpError::NotReduced(
                "quadrature did not produce a finite value and error estimate".into(),
            ));
        }
        Ok(record(vec![
            ("value", Value::Float(result.value)),
            ("evaluations", Value::int(budget.used.get())),
            ("error", Value::Float(result.error)),
            ("converged", Value::Bool(result.converged)),
        ]))
    }
}

#[derive(Clone, Copy)]
enum System {
    Broyden,
    Dogleg,
}
impl EvalDyn for System {
    fn signature(&self) -> Signature {
        effectful(Signature::new(
            vec![
                vector_type(),
                TypeExpr::Float,
                TypeExpr::Int,
                TypeExpr::Int,
                function(vector_type(), vector_type()),
            ],
            vector_result(),
        ))
    }

    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        arity(args, 5)?;
        let initial = vector(&args[0], &DimExpr::Dimless, None, 0, ctx)?;
        if initial.is_empty() {
            return Err(OpError::Domain("initial vector must be nonempty".into()));
        }
        let tol = tolerance(args, 1, ctx)?;
        let iterations = positive_budget(&args[2], 2, "system iteration budget")?;
        let budget = Budget::new(positive_budget(&args[3], 3, "system evaluation budget")?);
        let (callback, signature) =
            CallbackRun::new(ctx, &args[4], 4, 1, self.signature().params[4].clone())?;
        vector_callback(&signature, initial.len(), true)?;
        let evaluate = |x: &[f64]| {
            budget
                .evaluate(&callback, &[vec_value(x)], |value| {
                    vector(&value, &DimExpr::Dimless, Some(initial.len()), 4, ctx)
                })
                .unwrap_or_default()
        };
        let result = match self {
            | Self::Broyden => rootfind::broyden(evaluate, &initial, tol, iterations),
            | Self::Dogleg => rootfind::dogleg(evaluate, &initial, tol, iterations),
        };
        budget.finish()?;
        callback.finish()?;
        let result = result.map_err(root_error)?;
        Ok(record(vec![
            (
                "x",
                output_tensor(vec![result.x.len()], result.x, DimExpr::Dimless)?,
            ),
            ("evaluations", Value::int(budget.used.get())),
        ]))
    }
}

struct Halley;
impl EvalDyn for Halley {
    fn signature(&self) -> Signature {
        let f = function(TypeExpr::Float, TypeExpr::Float);
        effectful(Signature::new(
            vec![
                TypeExpr::Float,
                TypeExpr::Float,
                TypeExpr::Int,
                TypeExpr::Int,
                f.clone(),
                f.clone(),
                f,
            ],
            TypeExpr::Record(vec![
                ("root".into(), TypeExpr::Float),
                ("evaluations".into(), TypeExpr::Int),
            ]),
        ))
    }

    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        arity(args, 7)?;
        let initial = scalar(args, 0, ctx)?;
        let tol = tolerance(args, 1, ctx)?;
        let iterations = positive_budget(&args[2], 2, "Halley iteration budget")?;
        let budget = Budget::new(positive_budget(&args[3], 3, "Halley evaluation budget")?);
        let make = |index| -> Result<_, OpError> {
            let (callback, signature) = CallbackRun::new(
                ctx,
                &args[index],
                index,
                1,
                self.signature().params[index].clone(),
            )?;
            scalar_callback(&signature)?;
            Ok(callback)
        };
        let f = make(4)?;
        let df = make(5)?;
        let d2f = make(6)?;
        let evaluate = |callback: &CallbackRun<'_, '_>, x| {
            budget
                .evaluate(callback, &[Value::Float(x)], |value| numeric(&value, 4))
                .unwrap_or(f64::NAN)
        };
        let result = rootfind::halley(
            |x| evaluate(&f, x),
            |x| evaluate(&df, x),
            |x| evaluate(&d2f, x),
            initial,
            tol,
            iterations,
        );
        budget.finish()?;
        f.finish()?;
        df.finish()?;
        d2f.finish()?;
        let result = result.map_err(root_error)?;
        Ok(record(vec![
            ("root", Value::Float(result)),
            ("evaluations", Value::int(budget.used.get())),
        ]))
    }
}

#[derive(Clone, Copy)]
enum Optimizer {
    Bfgs,
    Lbfgs,
    ConjugateGradient,
    Box,
}
impl Optimizer {
    const fn bounded(self) -> bool {
        matches!(self, Self::Box)
    }
}
impl EvalDyn for Optimizer {
    fn signature(&self) -> Signature {
        let mut params = vec![vector_type()];
        if self.bounded() {
            params.extend([vector_type(), vector_type()]);
        }
        params.extend([TypeExpr::Float, TypeExpr::Int, TypeExpr::Int]);
        params.extend([
            function(vector_type(), TypeExpr::Float),
            function(vector_type(), vector_type()),
        ]);
        effectful(Signature::new(
            params,
            TypeExpr::Record(vec![
                ("x".into(), vector_type()),
                ("value".into(), TypeExpr::Float),
                ("iterations".into(), TypeExpr::Int),
                ("evaluations".into(), TypeExpr::Int),
                ("converged".into(), TypeExpr::Bool),
            ]),
        ))
    }

    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        arity(args, if self.bounded() { 8 } else { 6 })?;
        let initial = vector(&args[0], &DimExpr::Dimless, None, 0, ctx)?;
        let n = initial.len();
        if n == 0 {
            return Err(OpError::Domain("initial vector must be nonempty".into()));
        }
        let shift = if self.bounded() { 2 } else { 0 };
        let (lower, upper) = if self.bounded() {
            let lo = vector(&args[1], &DimExpr::Dimless, Some(n), 1, ctx)?;
            let hi = vector(&args[2], &DimExpr::Dimless, Some(n), 2, ctx)?;
            if lo.iter().zip(&hi).any(|(a, b)| a > b) {
                return Err(OpError::Domain(
                    "lower bounds must not exceed upper bounds".into(),
                ));
            }
            (lo, hi)
        } else {
            (Vec::new(), Vec::new())
        };
        let tol = tolerance(args, 1 + shift, ctx)?;
        let iterations =
            positive_budget(&args[2 + shift], 2 + shift, "optimizer iteration budget")?;
        let budget = Budget::new(positive_budget(
            &args[3 + shift],
            3 + shift,
            "optimizer evaluation budget",
        )?);
        let signature = self.signature();
        let (objective, fsig) = CallbackRun::new(
            ctx,
            &args[4 + shift],
            4 + shift,
            1,
            signature.params[4 + shift].clone(),
        )?;
        let (gradient, gsig) = CallbackRun::new(
            ctx,
            &args[5 + shift],
            5 + shift,
            1,
            signature.params[5 + shift].clone(),
        )?;
        vector_callback(&fsig, n, false)?;
        vector_callback(&gsig, n, true)?;
        let f = |x: &[f64]| {
            budget
                .evaluate(&objective, &[vec_value(x)], |value| {
                    numeric(&value, 4 + shift)
                })
                .unwrap_or(f64::NAN)
        };
        let g = |x: &[f64]| {
            budget
                .evaluate(&gradient, &[vec_value(x)], |value| {
                    vector(&value, &DimExpr::Dimless, Some(n), 5 + shift, ctx)
                })
                .unwrap_or_else(|| vec![f64::NAN; n])
        };
        let result = match self {
            | Self::Bfgs => optim::bfgs(f, Some(&g), &initial, tol, iterations),
            | Self::Lbfgs => optim::lbfgs(f, Some(&g), &initial, 10, tol, iterations),
            | Self::ConjugateGradient => {
                optim::nonlinear_cg(f, Some(&g), &initial, tol, iterations)
            },
            | Self::Box => {
                optim::lbfgs_b(f, Some(&g), &initial, &lower, &upper, 10, tol, iterations)
            },
        };
        budget.finish()?;
        objective.finish()?;
        gradient.finish()?;
        if !result.fx.is_finite() {
            return Err(OpError::Domain(
                "optimizer returned a nonfinite objective".into(),
            ));
        }
        Ok(record(vec![
            (
                "x",
                output_tensor(vec![result.x.len()], result.x, DimExpr::Dimless)?,
            ),
            ("value", Value::Float(result.fx)),
            ("iterations", Value::int(result.iterations)),
            ("evaluations", Value::int(budget.used.get())),
            ("converged", Value::Bool(result.converged)),
        ]))
    }
}

pub(super) fn install(installer: &mut Installer<'_>) -> Result<(), RuleError> {
    for (name, path, method) in [
        (
            "num_gauss_legendre",
            "crate::kernels::quadrature::gauss_legendre",
            Quad::GaussLegendre,
        ),
        (
            "num_clenshaw_curtis",
            "crate::kernels::quadrature::clenshaw_curtis",
            Quad::ClenshawCurtis,
        ),
        (
            "num_romberg",
            "crate::kernels::quadrature::romberg",
            Quad::Romberg,
        ),
        (
            "num_tanh_sinh",
            "crate::kernels::quadrature::tanh_sinh_with_levels",
            Quad::TanhSinh,
        ),
        (
            "num_integrate_semi_infinite",
            "crate::kernels::quadrature::integrate_semi_infinite",
            Quad::RightInfinite,
        ),
        (
            "num_integrate_left_infinite",
            "crate::kernels::quadrature::integrate_semi_infinite",
            Quad::LeftInfinite,
        ),
        (
            "num_integrate_infinite",
            "crate::kernels::quadrature::integrate_infinite",
            Quad::Infinite,
        ),
    ] {
        register(
            installer,
            name,
            path,
            "Dimensionless serial quadrature with a strict callback evaluation budget. Fixed rules: a,b,order,max_evaluations,callback. Romberg/tanh-sinh: a,b,tolerance,max_levels,max_evaluations,callback. Semi-infinite rules: finite_bound,tolerance,max_evaluations,callback; left-infinite integrates (-infinity,bound], right-infinite [bound,infinity). Full infinite rule: tolerance,max_evaluations,callback. Adaptive results include value,error,evaluations,converged; fixed rules report value,evaluations. Callback failures propagate unchanged.",
            method.count(),
            method,
        )?;
    }
    for (name, path, method) in [
        (
            "num_broyden",
            "crate::kernels::rootfind::broyden",
            System::Broyden,
        ),
        (
            "num_dogleg",
            "crate::kernels::rootfind::dogleg",
            System::Dogleg,
        ),
    ] {
        register(
            installer,
            name,
            path,
            "Dimensionless nonlinear system solver: initial vector,positive tolerance,max_iterations,max_evaluations,residual callback. Callback input and output retain the initial vector length. Returns x,evaluations; callback failures and iteration exhaustion are explicit errors.",
            5,
            method,
        )?;
    }
    register(
        installer,
        "num_halley",
        "crate::kernels::rootfind::halley",
        "Dimensionless Halley root solver: initial guess,positive tolerance,max_iterations,max_evaluations,f,first derivative,second derivative. The evaluation budget is shared across all three callbacks. Returns root,evaluations; callback failures propagate unchanged.",
        7,
        Halley,
    )?;
    for (name, path, method) in [
        ("num_bfgs", "crate::kernels::optim::bfgs", Optimizer::Bfgs),
        (
            "num_lbfgs",
            "crate::kernels::optim::lbfgs",
            Optimizer::Lbfgs,
        ),
        (
            "num_nonlinear_cg",
            "crate::kernels::optim::nonlinear_cg",
            Optimizer::ConjugateGradient,
        ),
        (
            "num_lbfgs_b",
            "crate::kernels::optim::lbfgs_b",
            Optimizer::Box,
        ),
    ] {
        register(
            installer,
            name,
            path,
            "Dimensionless serial optimization: initial vector,[lower bounds,upper bounds for lbfgs_b],positive tolerance,max_iterations,max_evaluations,objective callback,gradient callback. Both callbacks share a strict evaluation budget. Gradient length matches initial vector; L-BFGS methods retain ten correction pairs. Returns x,value,iterations,evaluations,converged; callback failures propagate unchanged.",
            if method.bounded() { 8 } else { 6 },
            method,
        )?;
    }
    Ok(())
}
