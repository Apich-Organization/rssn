//! Dimension-aware adaptive geometric ODE, residual DAE, and retarded DDE callbacks.
use super::CallCtx;
use super::DimExpr;
use super::EvalDyn;
use super::Installer;
use super::OpError;
use super::RecordValue;
use super::RuleError;
use super::Signature;
use super::TypeExpr;
use super::TypeParam;
use super::TypeParamKind;
use super::Value;
use super::canonical;
use super::check_dimension;
use super::common::CallbackRun;
use super::common::effectful;
use super::common::positive_budget;
use super::common::register;
use super::dimensional_float;
use super::dimensional_type;
use super::evolution::output_tensor;
use super::evolution::state_type;
use super::evolution::tensor_type;
use super::evolution::vector;
use super::evolution::vector_dimension;
use super::expression;
use super::numeric;
use super::scalar_dimension;
use crate::kernel::signature::ArityInfo;
use crate::kernel::signature::NatExpr;
use crate::kernel::value::TensorValue;
use crate::kernels::ode_geom::AdaptiveGeomOptions;
use crate::kernels::ode_geom::DaeOptions;
use crate::kernels::ode_geom::DdeOptions;
use crate::kernels::ode_geom::GeometricMethod;
use crate::kernels::ode_geom::{
    self,
};
use std::cell::Cell;

fn rate_dim(
    state: &DimExpr,
    time: &DimExpr,
) -> Result<DimExpr, OpError> {
    expression(&canonical(&DimExpr::Mul(vec![
        state.clone(),
        DimExpr::Inverse(Box::new(time.clone())),
    ]))?)
}
fn same(
    a: &DimExpr,
    b: &DimExpr,
) -> Result<(), OpError> {
    if canonical(a)? == canonical(b)? {
        Ok(())
    } else {
        Err(OpError::Domain(
            "callback dimensions do not match the declared differential equation".into(),
        ))
    }
}
fn vector_value(
    x: &[f64],
    dim: &DimExpr,
) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![x.len()],
        data: x.iter().copied().map(Value::Float).collect(),
        elem: dimensional_type(dim.clone()),
    })
}
fn solution_type(
    time: DimExpr,
    state: DimExpr,
    rate: DimExpr,
) -> TypeExpr {
    TypeExpr::Record(vec![
        ("times".into(), tensor_type(1, time)),
        ("states".into(), tensor_type(2, state)),
        ("rates".into(), tensor_type(2, rate)),
        ("evaluations".into(), TypeExpr::Int),
        ("rejected".into(), TypeExpr::Int),
    ])
}
fn solution_value(
    times: Vec<f64>,
    states: Vec<Vec<f64>>,
    rates: Vec<Vec<f64>>,
    evaluations: usize,
    rejected: usize,
    dimensions: (DimExpr, DimExpr, DimExpr),
) -> Result<Value, OpError> {
    let (time, state, rate) = dimensions;
    let rows = times.len();
    let cols = states.first().map_or(0, Vec::len);
    if states.len() != rows
        || rates.len() != rows
        || states.iter().chain(&rates).any(|v| v.len() != cols)
    {
        return Err(OpError::Domain(
            "malformed adaptive solution storage".into(),
        ));
    }
    Ok(Value::Record(RecordValue {
        fields: vec![
            ("times".into(), output_tensor(vec![rows], times, time)?),
            (
                "states".into(),
                output_tensor(
                    vec![rows, cols],
                    states.into_iter().flatten().collect(),
                    state,
                )?,
            ),
            (
                "rates".into(),
                output_tensor(
                    vec![rows, cols],
                    rates.into_iter().flatten().collect(),
                    rate,
                )?,
            ),
            ("evaluations".into(), Value::int(evaluations)),
            ("rejected".into(), Value::int(rejected)),
        ],
    }))
}
fn generic(
    mut signature: Signature,
    dimensions: &[&str],
    delays: bool,
) -> Signature {
    signature.type_params = dimensions
        .iter()
        .map(|name| {
            TypeParam {
                name: (*name).into(),
                kind: TypeParamKind::Dim,
                bounds: Vec::new(),
            }
        })
        .collect();
    signature.type_params.push(TypeParam {
        name: "N".into(),
        kind: TypeParamKind::Nat,
        bounds: Vec::new(),
    });
    if delays {
        signature.type_params.push(TypeParam {
            name: "K".into(),
            kind: TypeParamKind::Nat,
            bounds: Vec::new(),
        });
    }
    effectful(signature)
}
fn arity(
    args: &[Value],
    n: u8,
) -> Result<(), OpError> {
    if args.len() == usize::from(n) {
        Ok(())
    } else {
        Err(OpError::Arity {
            expected: ArityInfo::fixed(n),
            got: args.len(),
        })
    }
}
fn error(e: ode_geom::OdeGeomError) -> OpError {
    match e {
        | ode_geom::OdeGeomError::NotConverged => {
            OpError::NotReduced(
                "ODE/DAE/DDE tolerance, implicit iteration, or work limit could not be satisfied"
                    .into(),
            )
        },
        | other => OpError::Domain(other.to_string()),
    }
}
fn controls(
    args: &[Value],
    time: &DimExpr,
    state: &DimExpr,
    ctx: &CallCtx<'_>,
) -> Result<AdaptiveGeomOptions, OpError> {
    check_dimension(ctx, &args[4], &DimExpr::Dimless, 4)?;
    check_dimension(ctx, &args[3], state, 3)?;
    for i in [0, 1, 5, 6, 7] {
        check_dimension(ctx, &args[i], time, i)?;
    }
    Ok(AdaptiveGeomOptions {
        atol: numeric(&args[3], 3)?,
        rtol: numeric(&args[4], 4)?,
        h0: numeric(&args[5], 5)?,
        h_min: numeric(&args[6], 6)?,
        h_max: numeric(&args[7], 7)?,
        max_steps: positive_budget(&args[8], 8, "ODE attempted-step limit")?,
        max_evaluations: positive_budget(&args[9], 9, "ODE callback limit")?,
    })
}
struct Geometric(GeometricMethod);
impl EvalDyn for Geometric {
    fn signature(&self) -> Signature {
        let d = DimExpr::Var("D".into());
        let e = DimExpr::Var("E".into());
        let r = DimExpr::Mul(vec![e.clone(), DimExpr::Inverse(Box::new(d.clone()))]);
        generic(
            Signature::new(
                vec![
                    dimensional_type(d.clone()),
                    dimensional_type(d.clone()),
                    state_type(e.clone()),
                    dimensional_type(e.clone()),
                    TypeExpr::Float,
                    dimensional_type(d.clone()),
                    dimensional_type(d.clone()),
                    dimensional_type(d.clone()),
                    TypeExpr::Int,
                    TypeExpr::Int,
                    TypeExpr::Function(
                        vec![dimensional_type(d.clone()), state_type(e.clone())],
                        Box::new(state_type(r.clone())),
                    ),
                ],
                solution_type(d, e, r),
            ),
            &["D", "E"],
            false,
        )
    }

    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        arity(args, 11)?;
        let (callback, signature) =
            CallbackRun::new(ctx, &args[10], 10, 2, self.signature().params[10].clone())?;
        let d = scalar_dimension(&signature.params[0])?;
        let e = vector_dimension(&signature.params[1], None)?;
        let rate = rate_dim(&e, &d)?;
        let initial = vector(&args[2], &e, None, 2, ctx)?;
        vector_dimension(&signature.params[1], Some(initial.len()))?;
        same(
            &vector_dimension(&signature.ret, Some(initial.len()))?,
            &rate,
        )?;
        let options = controls(args, &d, &e, ctx)?;
        let mut rhs = |t: f64, y: &[f64]| {
            callback
                .evaluate(
                    &[dimensional_float(t, d.clone()), vector_value(y, &e)],
                    |v| vector(&v, &rate, Some(initial.len()), 10, ctx),
                )
                .unwrap_or_else(|| vec![f64::NAN; initial.len()])
        };
        let result = ode_geom::solve_adaptive_geometric(
            &mut rhs,
            &initial,
            (numeric(&args[0], 0)?, numeric(&args[1], 1)?),
            self.0,
            options,
            &[],
            &mut |_, _, _| 1.0,
            &Cell::new(false),
        );
        callback.finish()?;
        let result = result.map_err(error)?;
        solution_value(
            result.t,
            result.y,
            result.rates,
            result.nfev,
            result.rejected,
            (d, e, rate),
        )
    }
}

struct ImplicitDae;
impl EvalDyn for ImplicitDae {
    fn signature(&self) -> Signature {
        let d = DimExpr::Var("D".into());
        let e = DimExpr::Var("E".into());
        let residual = DimExpr::Var("R".into());
        let rate = DimExpr::Mul(vec![e.clone(), DimExpr::Inverse(Box::new(d.clone()))]);
        generic(
            Signature::new(
                vec![
                    dimensional_type(d.clone()),
                    dimensional_type(d.clone()),
                    state_type(e.clone()),
                    state_type(rate.clone()),
                    dimensional_type(e.clone()),
                    TypeExpr::Float,
                    dimensional_type(d.clone()),
                    dimensional_type(d.clone()),
                    dimensional_type(d.clone()),
                    TypeExpr::Int,
                    TypeExpr::Int,
                    dimensional_type(residual.clone()),
                    TypeExpr::Int,
                    TypeExpr::Function(
                        vec![
                            dimensional_type(d.clone()),
                            state_type(e.clone()),
                            state_type(rate.clone()),
                        ],
                        Box::new(state_type(residual)),
                    ),
                ],
                solution_type(d, e, rate),
            ),
            &["D", "E", "R"],
            false,
        )
    }

    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        arity(args, 14)?;
        let (callback, signature) =
            CallbackRun::new(ctx, &args[13], 13, 3, self.signature().params[13].clone())?;
        let d = scalar_dimension(&signature.params[0])?;
        let e = vector_dimension(&signature.params[1], None)?;
        let rate = rate_dim(&e, &d)?;
        let residual = vector_dimension(&signature.ret, None)?;
        let initial = vector(&args[2], &e, None, 2, ctx)?;
        let n = initial.len();
        vector_dimension(&signature.params[1], Some(n))?;
        same(&vector_dimension(&signature.params[2], Some(n))?, &rate)?;
        vector_dimension(&signature.ret, Some(n))?;
        let yp = vector(&args[3], &rate, Some(n), 3, ctx)?;
        for i in [0, 1, 6, 7, 8] {
            check_dimension(ctx, &args[i], &d, i)?;
        }
        check_dimension(ctx, &args[4], &e, 4)?;
        check_dimension(ctx, &args[5], &DimExpr::Dimless, 5)?;
        check_dimension(ctx, &args[11], &residual, 11)?;
        let options = DaeOptions {
            rtol: numeric(&args[5], 5)?,
            atol: numeric(&args[4], 4)?,
            h0: numeric(&args[6], 6)?,
            h_min: numeric(&args[7], 7)?,
            h_max: numeric(&args[8], 8)?,
            max_steps: positive_budget(&args[9], 9, "DAE attempted-step limit")?,
            max_evaluations: positive_budget(&args[10], 10, "DAE callback limit")?,
            newton_tol: numeric(&args[11], 11)?,
            max_newton_iter: positive_budget(&args[12], 12, "DAE Newton iteration limit")?,
        };
        let mut f = |t: f64, y: &[f64], yp: &[f64]| {
            callback
                .evaluate(
                    &[
                        dimensional_float(t, d.clone()),
                        vector_value(y, &e),
                        vector_value(yp, &rate),
                    ],
                    |v| vector(&v, &residual, Some(n), 13, ctx),
                )
                .unwrap_or_else(|| vec![f64::NAN; n])
        };
        let result = ode_geom::solve_dae_implicit(
            &mut f,
            &initial,
            &yp,
            (numeric(&args[0], 0)?, numeric(&args[1], 1)?),
            options,
        );
        callback.finish()?;
        let result = result.map_err(error)?;
        solution_value(
            result.t,
            result.y,
            result.rates,
            result.nfev,
            result.rejected,
            (d, e, rate),
        )
    }
}

struct ConstantDde;
impl EvalDyn for ConstantDde {
    fn signature(&self) -> Signature {
        let d = DimExpr::Var("D".into());
        let e = DimExpr::Var("E".into());
        let rate = DimExpr::Mul(vec![e.clone(), DimExpr::Inverse(Box::new(d.clone()))]);
        let delays = TypeExpr::Tensor {
            rank: NatExpr::Lit(1),
            shape: Some(vec![NatExpr::Var("K".into())]),
            elem: Box::new(dimensional_type(d.clone())),
        };
        let history = TypeExpr::Tensor {
            rank: NatExpr::Lit(2),
            shape: Some(vec![NatExpr::Var("K".into()), NatExpr::Var("N".into())]),
            elem: Box::new(dimensional_type(e.clone())),
        };
        generic(
            Signature::new(
                vec![
                    dimensional_type(d.clone()),
                    dimensional_type(d.clone()),
                    state_type(e.clone()),
                    delays,
                    dimensional_type(e.clone()),
                    TypeExpr::Float,
                    dimensional_type(d.clone()),
                    dimensional_type(d.clone()),
                    dimensional_type(d.clone()),
                    TypeExpr::Int,
                    TypeExpr::Int,
                    TypeExpr::Function(
                        vec![dimensional_type(d.clone()), state_type(e.clone()), history],
                        Box::new(state_type(rate.clone())),
                    ),
                    TypeExpr::Function(
                        vec![dimensional_type(d.clone())],
                        Box::new(state_type(e.clone())),
                    ),
                ],
                solution_type(d, e, rate),
            ),
            &["D", "E"],
            true,
        )
    }

    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        arity(args, 13)?;
        let schema = self.signature();
        let (rhs, signature) = CallbackRun::new(ctx, &args[11], 11, 3, schema.params[11].clone())?;
        let (history, hsig) = CallbackRun::new(ctx, &args[12], 12, 1, schema.params[12].clone())?;
        let d = scalar_dimension(&signature.params[0])?;
        let e = vector_dimension(&signature.params[1], None)?;
        let rate = rate_dim(&e, &d)?;
        let initial = vector(&args[2], &e, None, 2, ctx)?;
        let n = initial.len();
        vector_dimension(&signature.params[1], Some(n))?;
        same(&vector_dimension(&signature.ret, Some(n))?, &rate)?;
        same(&scalar_dimension(&hsig.params[0])?, &d)?;
        same(&vector_dimension(&hsig.ret, Some(n))?, &e)?;
        let delays = vector(&args[3], &d, None, 3, ctx)?;
        match &signature.params[2] {
            | TypeExpr::Tensor {
                rank: NatExpr::Lit(2),
                shape,
                elem,
            } => {
                same(&scalar_dimension(elem)?, &e)?;
                if let Some(shape) = shape
                    && shape.as_slice()
                        != [NatExpr::Lit(delays.len() as u64), NatExpr::Lit(n as u64)]
                {
                    return Err(OpError::Domain(
                        "DDE callback history matrix must have delay-count by state-count shape"
                            .into(),
                    ));
                }
            },
            | _ => {
                return Err(OpError::Domain(
                    "DDE callback history argument must be a rank-two state matrix".into(),
                ));
            },
        }

        for i in [0, 1, 6, 7, 8] {
            check_dimension(ctx, &args[i], &d, i)?;
        }
        check_dimension(ctx, &args[4], &e, 4)?;
        check_dimension(ctx, &args[5], &DimExpr::Dimless, 5)?;
        let options = DdeOptions {
            rtol: numeric(&args[5], 5)?,
            atol: numeric(&args[4], 4)?,
            h0: numeric(&args[6], 6)?,
            h_min: numeric(&args[7], 7)?,
            h_max: numeric(&args[8], 8)?,
            max_steps: positive_budget(&args[9], 9, "DDE attempted-step limit")?,
            max_evaluations: positive_budget(&args[10], 10, "DDE callback limit")?,
        };
        let mut f = |t: f64, y: &[f64], past: &[Vec<f64>]| {
            if rhs.failed() || history.failed() {
                return vec![f64::NAN; n];
            }
            let matrix = Value::Tensor(TensorValue {
                shape: vec![past.len(), n],
                data: past.iter().flatten().copied().map(Value::Float).collect(),
                elem: dimensional_type(e.clone()),
            });
            rhs.evaluate(
                &[dimensional_float(t, d.clone()), vector_value(y, &e), matrix],
                |v| vector(&v, &rate, Some(n), 11, ctx),
            )
            .unwrap_or_else(|| vec![f64::NAN; n])
        };
        let mut hist = |t: f64| {
            if rhs.failed() || history.failed() {
                return vec![f64::NAN; n];
            }
            history
                .evaluate(&[dimensional_float(t, d.clone())], |v| {
                    vector(&v, &e, Some(n), 12, ctx)
                })
                .unwrap_or_else(|| vec![f64::NAN; n])
        };
        let result = ode_geom::solve_dde(
            &mut f,
            &mut hist,
            &initial,
            &delays,
            (numeric(&args[0], 0)?, numeric(&args[1], 1)?),
            options,
        );
        rhs.finish()?;
        history.finish()?;
        let result = result.map_err(error)?;
        solution_value(
            result.t,
            result.y,
            result.rates,
            result.nfev,
            result.rejected,
            (d, e, rate),
        )
    }
}

pub(super) fn install(installer: &mut Installer<'_>) -> Result<(), RuleError> {
    for (name, method) in [
        ("num_ode_rk4_adaptive", GeometricMethod::Rk4),
        ("num_ode_verlet_adaptive", GeometricMethod::StormerVerlet),
        ("num_ode_ruth4_adaptive", GeometricMethod::Ruth4),
        (
            "num_ode_midpoint_adaptive",
            GeometricMethod::ImplicitMidpoint,
        ),
        ("num_ode_gauss4_adaptive", GeometricMethod::GaussLegendre4),
        ("num_ode_gauss6_adaptive", GeometricMethod::GaussLegendre6),
    ] {
        register(
            installer,
            name,
            "crate::kernels::ode_geom::solve_adaptive_geometric",
            "Dimension-aware adaptive requested-method integration: t0,t1,y0,atol,rtol,h0,h_min,h_max,max_steps,max_evaluations,rhs. Callback derivative has state/time dimensions. Returns accepted times/states/Hermite rates and actual work counters; no method substitution. Adaptive state-dependent stepping does not preserve fixed-step symplecticity.",
            11,
            Geometric(method),
        )?;
    }
    for index2 in [false, true] {
        register(
            installer,
            if index2 {
                "num_dae_index2_adaptive"
            } else {
                "num_dae_index1_adaptive"
            },
            if index2 {
                "crate::kernels::ode_geom::solve_dae_index2"
            } else {
                "crate::kernels::ode_geom::solve_dae_index1"
            },
            "Dimension-aware semi-explicit DAE with explicit nondimensionalization: t0,t1,y0,z0,time_scale,state_scale,algebraic_scale,constraint_scale,normalized_atol,rtol,normalized_h0,normalized_h_min,normalized_h_max,max_steps,max_evaluations,normalized_newton_tolerance,max_newton_iterations,rhs,constraint. Four positive physical reference scales make mixed state/constraint Newton blocks dimensionless. Index1 constraint receives t,y,z; index2 receives t,y. Returns physical times,states,algebraic,rates and actual work counts.",
            19,
            SemiDae { index2 },
        )?;
    }
    register(
        installer,
        "num_dae_implicit_adaptive",
        "crate::kernels::ode_geom::solve_dae_implicit",
        "Dimension-aware residual DAE BDF1/BDF2 with genuine step-doubling control: t0,t1,y0,yp0,atol,rtol,h0,h_min,h_max,max_steps,max_evaluations,newton_tolerance,max_newton_iterations,residual. State has arbitrary E dimensions, derivatives E/D, residuals independent R; Newton tolerance carries R. Algebraic initial values must be consistent. Returns times/states/rates/evaluations/rejected.",
        14,
        ImplicitDae,
    )?;
    register(
        installer,
        "num_dde_constant_adaptive",
        "crate::kernels::ode_geom::solve_dde",
        "Dimension-aware constant-delay method-of-steps DP5(4): t0,t1,y0,positive_delays,atol,rtol,h0,h_min,h_max,max_steps,max_evaluations,rhs,history. RHS receives t,y,delay-by-state matrix and returns state/time. History receives time and returns state. Strict total RHS/prehistory callback budget; no fabricated intra-step history. Returns times/states/Hermite rates/evaluations/rejected.",
        13,
        ConstantDde,
    )
}

struct SemiDae {
    index2: bool,
}
fn algebraic_type(dim: DimExpr) -> TypeExpr {
    TypeExpr::Tensor {
        rank: NatExpr::Lit(1),
        shape: Some(vec![NatExpr::Var("M".into())]),
        elem: Box::new(dimensional_type(dim)),
    }
}
impl EvalDyn for SemiDae {
    fn signature(&self) -> Signature {
        let d = DimExpr::Var("D".into());
        let e = DimExpr::Var("E".into());
        let z = DimExpr::Var("Z".into());
        let r = DimExpr::Var("R".into());
        let rate = DimExpr::Mul(vec![e.clone(), DimExpr::Inverse(Box::new(d.clone()))]);
        let mut constraint_args = vec![dimensional_type(d.clone()), state_type(e.clone())];
        if !self.index2 {
            constraint_args.push(algebraic_type(z.clone()));
        }
        let mut params = vec![
            dimensional_type(d.clone()),
            dimensional_type(d.clone()),
            state_type(e.clone()),
            algebraic_type(z.clone()),
            dimensional_type(d.clone()),
            dimensional_type(e.clone()),
            dimensional_type(z.clone()),
            dimensional_type(r.clone()),
        ];
        params.extend([
            TypeExpr::Float,
            TypeExpr::Float,
            TypeExpr::Float,
            TypeExpr::Float,
            TypeExpr::Float,
            TypeExpr::Int,
            TypeExpr::Int,
            TypeExpr::Float,
            TypeExpr::Int,
        ]);
        params.extend([
            TypeExpr::Function(
                vec![
                    dimensional_type(d.clone()),
                    state_type(e.clone()),
                    algebraic_type(z.clone()),
                ],
                Box::new(state_type(rate.clone())),
            ),
            TypeExpr::Function(constraint_args, Box::new(algebraic_type(r))),
        ]);
        let result = TypeExpr::Record(vec![
            ("times".into(), tensor_type(1, d)),
            ("states".into(), tensor_type(2, e)),
            ("algebraic".into(), tensor_type(2, z)),
            ("rates".into(), tensor_type(2, rate)),
            ("evaluations".into(), TypeExpr::Int),
            ("rejected".into(), TypeExpr::Int),
        ]);
        let mut sig = generic(Signature::new(params, result), &["D", "E", "Z", "R"], false);
        sig.type_params.push(TypeParam {
            name: "M".into(),
            kind: TypeParamKind::Nat,
            bounds: Vec::new(),
        });
        sig
    }

    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        arity(args, 19)?;
        let schema = self.signature();
        let (rhs, fsig) = CallbackRun::new(ctx, &args[17], 17, 3, schema.params[17].clone())?;
        let (constraint, gsig) = CallbackRun::new(
            ctx,
            &args[18],
            18,
            if self.index2 { 2 } else { 3 },
            schema.params[18].clone(),
        )?;
        let d = scalar_dimension(&fsig.params[0])?;
        let e = vector_dimension(&fsig.params[1], None)?;
        let z = vector_dimension(&fsig.params[2], None)?;
        let r = vector_dimension(&gsig.ret, None)?;
        let rate = rate_dim(&e, &d)?;
        let initial = vector(&args[2], &e, None, 2, ctx)?;
        let algebraic = vector(&args[3], &z, None, 3, ctx)?;
        let n = initial.len();
        let m = algebraic.len();
        vector_dimension(&fsig.params[1], Some(n))?;
        vector_dimension(&fsig.params[2], Some(m))?;
        same(&vector_dimension(&fsig.ret, Some(n))?, &rate)?;
        same(&scalar_dimension(&gsig.params[0])?, &d)?;
        same(&vector_dimension(&gsig.params[1], Some(n))?, &e)?;
        vector_dimension(&gsig.ret, Some(m))?;
        if !self.index2 {
            same(&vector_dimension(&gsig.params[2], Some(m))?, &z)?;
        }
        for (i, dim) in [(0, &d), (1, &d), (4, &d), (5, &e), (6, &z), (7, &r)] {
            check_dimension(ctx, &args[i], dim, i)?;
        }
        for i in [8, 9, 10, 11, 12, 15] {
            check_dimension(ctx, &args[i], &DimExpr::Dimless, i)?;
        }
        let t0 = numeric(&args[0], 0)?;
        let t1 = numeric(&args[1], 1)?;
        let scales = [
            numeric(&args[4], 4)?,
            numeric(&args[5], 5)?,
            numeric(&args[6], 6)?,
            numeric(&args[7], 7)?,
        ];
        if scales.iter().any(|v| *v <= 0.0) {
            return Err(OpError::Domain(
                "DAE reference scales must be finite and positive".into(),
            ));
        }
        let [td, yd, zd, gd] = scales;
        let options = DaeOptions {
            atol: numeric(&args[8], 8)?,
            rtol: numeric(&args[9], 9)?,
            h0: numeric(&args[10], 10)?,
            h_min: numeric(&args[11], 11)?,
            h_max: numeric(&args[12], 12)?,
            max_steps: positive_budget(&args[13], 13, "DAE attempted-step limit")?,
            max_evaluations: positive_budget(&args[14], 14, "DAE callback limit")?,
            newton_tol: numeric(&args[15], 15)?,
            max_newton_iter: positive_budget(&args[16], 16, "DAE Newton iteration limit")?,
        };
        let y0: Vec<f64> = initial.iter().map(|x| x / yd).collect();
        let z0: Vec<f64> = algebraic.iter().map(|x| x / zd).collect();
        let mut f = |t: f64, y: &[f64], zv: &[f64]| {
            if rhs.failed() || constraint.failed() {
                return vec![f64::NAN; n];
            }
            let y: Vec<f64> = y.iter().map(|x| x * yd).collect();
            let zv: Vec<f64> = zv.iter().map(|x| x * zd).collect();
            rhs.evaluate(
                &[
                    dimensional_float(t0 + t * td, d.clone()),
                    vector_value(&y, &e),
                    vector_value(&zv, &z),
                ],
                |v| vector(&v, &rate, Some(n), 17, ctx),
            )
            .map_or_else(
                || vec![f64::NAN; n],
                |out| out.into_iter().map(|x| x * td / yd).collect(),
            )
        };
        let eval_g = |t: f64, y: &[f64], zv: Option<&[f64]>| {
            if rhs.failed() || constraint.failed() {
                return vec![f64::NAN; m];
            }
            let y: Vec<f64> = y.iter().map(|x| x * yd).collect();
            let mut values = vec![
                dimensional_float(t0 + t * td, d.clone()),
                vector_value(&y, &e),
            ];
            if let Some(zv) = zv {
                values.push(vector_value(
                    &zv.iter().map(|x| x * zd).collect::<Vec<_>>(),
                    &z,
                ));
            }
            constraint
                .evaluate(&values, |v| vector(&v, &r, Some(m), 18, ctx))
                .map_or_else(
                    || vec![f64::NAN; m],
                    |out| out.into_iter().map(|x| x / gd).collect(),
                )
        };
        let result = if self.index2 {
            ode_geom::solve_dae_index2(
                &mut f,
                &mut |t, y| eval_g(t, y, None),
                &y0,
                &z0,
                (0.0, (t1 - t0) / td),
                options,
            )
        } else {
            ode_geom::solve_dae_index1(
                &mut f,
                &mut |t, y, z| eval_g(t, y, Some(z)),
                &y0,
                &z0,
                (0.0, (t1 - t0) / td),
                options,
            )
        };
        rhs.finish()?;
        constraint.finish()?;
        let result = result.map_err(error)?;
        let rows = result.t.len();
        Ok(Value::Record(RecordValue {
            fields: vec![
                (
                    "times".into(),
                    output_tensor(
                        vec![rows],
                        result.t.into_iter().map(|x| t0 + x * td).collect(),
                        d,
                    )?,
                ),
                (
                    "states".into(),
                    output_tensor(
                        vec![rows, n],
                        result.y.into_iter().flatten().map(|x| x * yd).collect(),
                        e,
                    )?,
                ),
                (
                    "algebraic".into(),
                    output_tensor(
                        vec![rows, m],
                        result.z.into_iter().flatten().map(|x| x * zd).collect(),
                        z,
                    )?,
                ),
                (
                    "rates".into(),
                    output_tensor(
                        vec![rows, n],
                        result
                            .dy
                            .into_iter()
                            .flatten()
                            .map(|x| x * yd / td)
                            .collect(),
                        rate,
                    )?,
                ),
                ("evaluations".into(), Value::int(result.nfev)),
                ("rejected".into(), Value::int(result.rejected)),
            ],
        }))
    }
}
