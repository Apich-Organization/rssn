//! Native constrained/global optimization in independently scaled physical units.
use super::CallCtx;
use super::Callbacks;
use super::DimExpr;
use super::EvalDyn;
use super::Installer;
use super::OpError;
use super::RuleError;
use super::Signature;
use super::TypeExpr;
use super::Value;
use super::bounded_size;
use super::count;
use super::equal;
use super::finite;
use super::function;
use super::positive;
use super::ratio;
use super::real;
use super::record;
use super::register;
use super::scalar;
use super::scale;
use super::tensor_type;
use super::typed;
use super::var;
use super::vec_value;
use super::vector;
use super::vector_dimension;
use crate::kernel::value::VariantValue;
use crate::kernels::optim;
use num_traits::ToPrimitive;
#[derive(Clone, Copy)]
enum Method {
    NelderMead,
    DifferentialEvolution,
    Annealing,
    Augmented,
    Sqp,
    SecondOrder,
}
impl Method {
    const fn constrained(self) -> bool {
        matches!(self, Self::Augmented | Self::Sqp | Self::SecondOrder)
    }

    const fn start(self) -> usize {
        match self {
            | Self::NelderMead => 7,
            | Self::Annealing => 9,
            | _ => 8,
        }
    }
}
fn optional_number(v: f64) -> Value {
    Value::Variant(VariantValue {
        tag: if v.is_finite() {
            "Some"
        } else {
            "None"
        }
        .into(),
        payload: if v.is_finite() {
            vec![Value::Float(v)]
        } else {
            vec![]
        },
    })
}
fn optional_float() -> TypeExpr {
    TypeExpr::Variant(vec![
        ("None".into(), vec![]),
        ("Some".into(), vec![TypeExpr::Float]),
    ])
}
impl EvalDyn for Method {
    // One match is the auditable ABI schema for this small native method family.
    #[allow(clippy::too_many_lines)]
    fn signature(&self) -> Signature {
        let d = var("D");
        let e = var("E");
        let c = var("C");
        let h = var("H");
        let x = tensor_type(1, d.clone());
        let mut args = match self {
            | Self::NelderMead => {
                vec![
                    x.clone(),
                    real(d.clone()),
                    real(e.clone()),
                    real(d),
                    TypeExpr::Float,
                    TypeExpr::Int,
                    TypeExpr::Int,
                ]
            },
            | Self::DifferentialEvolution => {
                vec![
                    x.clone(),
                    x.clone(),
                    real(d),
                    real(e.clone()),
                    TypeExpr::Int,
                    TypeExpr::Int,
                    TypeExpr::Int,
                    TypeExpr::Int,
                ]
            },
            | Self::Annealing => {
                vec![
                    x.clone(),
                    real(d.clone()),
                    real(e.clone()),
                    real(d),
                    real(e.clone()),
                    TypeExpr::Float,
                    TypeExpr::Int,
                    TypeExpr::Int,
                    TypeExpr::Int,
                ]
            },
            | Self::SecondOrder => {
                vec![
                    x.clone(),
                    tensor_type(
                        1,
                        DimExpr::Mul(vec![e.clone(), DimExpr::Inverse(Box::new(c.clone()))]),
                    ),
                    tensor_type(
                        1,
                        DimExpr::Mul(vec![e.clone(), DimExpr::Inverse(Box::new(h.clone()))]),
                    ),
                    real(d),
                    real(e.clone()),
                    real(c.clone()),
                    real(h.clone()),
                    TypeExpr::Int,
                ]
            },
            | _ => {
                vec![
                    x.clone(),
                    real(d),
                    real(e.clone()),
                    real(c.clone()),
                    real(h.clone()),
                    TypeExpr::Float,
                    TypeExpr::Int,
                    TypeExpr::Int,
                ]
            },
        };
        args.push(function(vec![x.clone()], real(e.clone())));
        if self.constrained() {
            args.push(function(vec![x.clone()], tensor_type(1, c.clone())));
            args.push(function(vec![x.clone()], tensor_type(1, h.clone())));
        }
        let mut fields = vec![
            ("x".into(), x),
            ("value".into(), real(e.clone())),
            ("iterations".into(), TypeExpr::Int),
            ("evaluations".into(), TypeExpr::Int),
            ("converged".into(), TypeExpr::Bool),
        ];
        if matches!(self, Self::Sqp) {
            fields.extend([
                (
                    "lambda".into(),
                    tensor_type(
                        1,
                        DimExpr::Mul(vec![e.clone(), DimExpr::Inverse(Box::new(c))]),
                    ),
                ),
                (
                    "mu".into(),
                    tensor_type(1, DimExpr::Mul(vec![e, DimExpr::Inverse(Box::new(h))])),
                ),
                ("status".into(), TypeExpr::String),
                ("second_order".into(), TypeExpr::String),
                ("normalized_violation".into(), optional_float()),
                ("normalized_kkt_residual".into(), optional_float()),
            ]);
        }
        typed(
            args,
            if matches!(self, Self::SecondOrder) {
                TypeExpr::String
            } else {
                TypeExpr::Record(fields)
            },
            if self.constrained() {
                &["D", "E", "C", "H"]
            } else {
                &["D", "E"]
            },
        )
    }

    // Keep each native method's validation, callback conversion and result units together.
    #[allow(clippy::too_many_lines)]
    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let expected = self.signature();
        super::arity(args, expected.params.len())?;
        let start = self.start();
        let budget = if matches!(self, Self::SecondOrder | Self::Annealing) {
            7
        } else if matches!(self, Self::DifferentialEvolution) {
            6
        } else {
            start - 1
        };
        let callbacks = Callbacks::new(ctx, args, &expected, start, count(&args[budget], budget)?)?;
        let sig = &callbacks.signatures[0];
        let d = vector_dimension(&sig.params[0], None)?;
        let e = super::super::scalar_dimension(&sig.ret)?;
        let (x_at, y_at) = if matches!(self, Self::DifferentialEvolution) {
            (2, 3)
        } else if matches!(self, Self::SecondOrder) {
            (3, 4)
        } else {
            (1, 2)
        };
        let xs = scale(&args[x_at], &d, x_at, ctx)?;
        let ys = scale(&args[y_at], &e, y_at, ctx)?;
        let mut initial = vector(&args[0], &d, None, 0, ctx)?;
        if initial.is_empty() {
            return Err(super::domain("optimization coordinates must be nonempty"));
        }
        let n = initial.len();
        bounded_size(n, n)?;
        equal(&vector_dimension(&sig.params[0], Some(n))?, &d)?;
        for v in &mut initial {
            *v /= xs;
        }
        let f = |x: &[f64]| {
            let x: Vec<_> = x.iter().map(|v| v * xs).collect();
            callbacks
                .call(0, &[vec_value(&x, &d)], |v| {
                    scalar(&v, &e, start, ctx).map(|v| v / ys)
                })
                .unwrap_or(f64::NAN)
        };
        let seed_at = if matches!(self, Self::DifferentialEvolution) {
            7
        } else {
            8
        };
        let seed = if matches!(self, Self::DifferentialEvolution | Self::Annealing) {
            let Value::Int(v) = &args[seed_at] else {
                return Err(super::domain(
                    "global optimization seed must be an exact integer",
                ));
            };
            v.to_u64()
                .ok_or_else(|| super::domain("global optimization seed must fit u64"))?
        } else {
            0
        };
        let mut constrained = None;
        let result = match self {
            | Self::NelderMead => {
                let step = scale(&args[3], &d, 3, ctx)? / xs;
                optim::nelder_mead(
                    f,
                    &initial,
                    step,
                    positive(&args[4], 4, ctx)?,
                    count(&args[5], 5)?,
                )
            },
            | Self::DifferentialEvolution => {
                let upper = vector(&args[1], &d, Some(n), 1, ctx)?;
                let mut bounds = Vec::with_capacity(n);
                for (lower, upper) in initial.iter().zip(upper) {
                    let upper = upper / xs;
                    if *lower > upper {
                        return Err(super::domain("lower bound exceeds upper bound"));
                    }
                    bounds.push((*lower, upper));
                }
                let population = count(&args[4], 4)?.max(4);
                bounded_size(population, n)?;
                optim::differential_evolution(f, &bounds, population, count(&args[5], 5)?, seed)
            },
            | Self::Annealing => {
                let step = scale(&args[3], &d, 3, ctx)? / xs;
                let temperature = scale(&args[4], &e, 4, ctx)? / ys;
                let cooling = positive(&args[5], 5, ctx)?;
                if cooling > 1.0 {
                    return Err(super::domain("annealing cooling must not exceed one"));
                }
                optim::simulated_annealing(
                    f,
                    &initial,
                    step,
                    temperature,
                    cooling,
                    count(&args[6], 6)?,
                    seed,
                )
            },
            | Self::Augmented | Self::Sqp | Self::SecondOrder => {
                let c = vector_dimension(&callbacks.signatures[1].ret, None)?;
                let h = vector_dimension(&callbacks.signatures[2].ret, None)?;
                for s in &callbacks.signatures[1..] {
                    equal(&vector_dimension(&s.params[0], Some(n))?, &d)?;
                }
                let (c_at, h_at) = if matches!(self, Self::SecondOrder) {
                    (5, 6)
                } else {
                    (3, 4)
                };
                let cs = scale(&args[c_at], &c, c_at, ctx)?;
                let hs = scale(&args[h_at], &h, h_at, ctx)?;
                let physical: Vec<_> = initial.iter().map(|v| v * xs).collect();
                let eq0 = callbacks
                    .call(1, &[vec_value(&physical, &d)], |v| {
                        vector(&v, &c, None, start + 1, ctx)
                    })
                    .unwrap_or_default();
                let ineq0 = callbacks
                    .call(2, &[vec_value(&physical, &d)], |v| {
                        vector(&v, &h, None, start + 2, ctx)
                    })
                    .unwrap_or_default();
                if callbacks
                    .runs
                    .iter()
                    .any(super::super::common::CallbackRun::failed)
                    || callbacks.exhausted.get()
                {
                    return Err(callbacks
                        .finish()
                        .err()
                        .unwrap_or_else(|| super::domain("constraint callback failed")));
                }
                let ne = eq0.len();
                let ni = ineq0.len();
                bounded_size(ne.saturating_add(ni), n)?;
                equal(
                    &vector_dimension(&callbacks.signatures[1].ret, Some(ne))?,
                    &c,
                )?;
                equal(
                    &vector_dimension(&callbacks.signatures[2].ret, Some(ni))?,
                    &h,
                )?;
                let constraint = |index: usize,
                                  x: &[f64],
                                  dimension: &DimExpr,
                                  length: usize,
                                  reference: f64| {
                    let x: Vec<_> = x.iter().map(|v| v * xs).collect();
                    callbacks
                        .call(index, &[vec_value(&x, &d)], |v| {
                            vector(&v, dimension, Some(length), start + index, ctx)
                                .map(|v| v.into_iter().map(|v| v / reference).collect())
                        })
                        .unwrap_or_else(|| vec![f64::NAN; length])
                };
                let eq = |x: &[f64]| constraint(1, x, &c, ne, cs);
                let ineq = |x: &[f64]| constraint(2, x, &h, ni, hs);
                if matches!(self, Self::SecondOrder) {
                    let lambda = vector(&args[1], &ratio(&e, &c)?, Some(ne), 1, ctx)?;
                    let mu = vector(&args[2], &ratio(&e, &h)?, Some(ni), 2, ctx)?;
                    let lambda: Vec<_> = lambda.into_iter().map(|v| v * cs / ys).collect();
                    let mu: Vec<_> = mu.into_iter().map(|v| v * hs / ys).collect();
                    let status = optim::second_order_check(f, eq, ineq, &initial, &lambda, &mu);
                    callbacks.finish()?;
                    return Ok(Value::String(format!("{status:?}")));
                }
                let tolerance = positive(&args[5], 5, ctx)?;
                let iterations = count(&args[6], 6)?;
                if matches!(self, Self::Augmented) {
                    optim::augmented_lagrangian(f, eq, ineq, &initial, tolerance, iterations)
                } else {
                    let out = optim::sqp(
                        f,
                        None::<&fn(&[f64]) -> Vec<f64>>,
                        eq,
                        ineq,
                        &initial,
                        tolerance,
                        iterations,
                    );
                    let result = optim::OptResult {
                        x: out.x.clone(),
                        fx: out.fx,
                        iterations: out.iterations,
                        evaluations: out.evaluations,
                        converged: matches!(out.status, optim::ConstrainedStatus::Optimal),
                    };
                    constrained = Some((out, c, h, cs, hs));
                    result
                }
            },
        };
        let used = callbacks.used.get();
        callbacks.finish()?;
        let mut fields = vec![
            (
                "x",
                vec_value(
                    &result.x.into_iter().map(|v| v * xs).collect::<Vec<_>>(),
                    &d,
                ),
            ),
            (
                "value",
                super::super::dimensional_float(finite(result.fx * ys)?, e.clone()),
            ),
            ("iterations", Value::int(result.iterations)),
            ("evaluations", Value::int(used)),
            ("converged", Value::Bool(result.converged)),
        ];
        if let Some((out, c, h, cs, hs)) = constrained {
            fields.extend([
                (
                    "lambda",
                    vec_value(
                        &out.lambda
                            .into_iter()
                            .map(|v| v * ys / cs)
                            .collect::<Vec<_>>(),
                        &ratio(&e, &c)?,
                    ),
                ),
                (
                    "mu",
                    vec_value(
                        &out.mu.into_iter().map(|v| v * ys / hs).collect::<Vec<_>>(),
                        &ratio(&e, &h)?,
                    ),
                ),
                ("status", Value::String(format!("{:?}", out.status))),
                (
                    "second_order",
                    Value::String(format!("{:?}", out.second_order)),
                ),
                ("normalized_violation", optional_number(out.violation)),
                ("normalized_kkt_residual", optional_number(out.kkt_residual)),
            ]);
        }
        Ok(record(fields))
    }
}
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    for (name, path, m) in [
        ("nelder_mead", "nelder_mead", Method::NelderMead),
        (
            "differential_evolution",
            "differential_evolution",
            Method::DifferentialEvolution,
        ),
        (
            "simulated_annealing",
            "simulated_annealing",
            Method::Annealing,
        ),
        (
            "augmented_lagrangian",
            "augmented_lagrangian",
            Method::Augmented,
        ),
        ("sqp", "sqp", Method::Sqp),
        (
            "second_order_check",
            "second_order_check",
            Method::SecondOrder,
        ),
    ] {
        register(
            i,
            &format!("num_native_{name}_scaled"),
            &format!("crate::kernels::optim::{path}"),
            "Native scaled global/constrained optimization with shared actual host evaluation limit. Nelder: x,xscale,yscale,step,tol,iterations,budget,f. Differential evolution: lo,hi,xscale,yscale,population,generations,budget,seed,f. Annealing: x,xscale,yscale,step,temperature,cooling,iterations,budget,seed,f. Constrained: x,xscale,yscale,equality_scale,inequality_scale,tol,iterations,budget,f,eq,ineq. Constraints independently carry C/H units; SQP returns physical multipliers E/C and E/H, native status and second-order verdict. Second-order: x,lambda,mu,xscale,yscale,eqscale,ineqscale,budget,f,eq,ineq. No global-optimality claim; native iteration/status semantics retained.",
            (m.start() + if m.constrained() { 3 } else { 1 }) as u8,
            m,
        )?;
    }
    Ok(())
}
