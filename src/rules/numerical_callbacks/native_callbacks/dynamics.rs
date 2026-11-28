//! Checked complex dynamics with a user-defined map and derivative.

use std::cell::Cell;

use crate::kernels::complex_dynamics as native;

use super::super::domain;
use super::super::record;
use super::CallCtx;
use super::Callbacks;
use super::Complex64;
use super::DimExpr;
use super::EvalDyn;
use super::Installer;
use super::OpError;
use super::RuleError;
use super::Signature;
use super::TypeExpr;
use super::Value;
use super::complex_type;
use super::count;
use super::decode;
use super::dimension;
use super::equal;
use super::finite;
use super::function;
use super::ratio;
use super::real;
use super::register;
use super::scale;
use super::typed;
use super::value;
use super::var;

#[derive(Clone, Copy)]
enum Method {
    Newton,
    Koenigs,
}

const MAX_STEPS: usize = 16_384;

fn finite_complex(z: Complex64) -> Result<Complex64, OpError> {
    finite(z.re)?;
    finite(z.im)?;
    Ok(z)
}

fn sample(
    callbacks: &Callbacks<'_, '_>,
    index: usize,
    point: Complex64,
    scales: (f64, f64),
    dimensions: (&DimExpr, &DimExpr),
    ctx: &CallCtx<'_>,
    invalid: &Cell<bool>,
) -> Complex64 {
    let (axis, output_scale) = scales;
    let (d, output) = dimensions;
    let point = point * axis;
    if invalid.get() || !point.is_finite() {
        invalid.set(true);
        return Complex64::new(f64::NAN, f64::NAN);
    }
    callbacks
        .call(index, &[value(point, d.clone())], |v| {
            finite_complex(decode(&v, output, index, ctx)? / output_scale)
        })
        .unwrap_or(Complex64::new(f64::NAN, f64::NAN))
}

fn checked_steps(
    steps: &Value,
    index: usize,
) -> Result<usize, OpError> {
    let steps = count(steps, index)?;
    if steps > MAX_STEPS {
        return Err(OpError::BudgetExceeded);
    }
    Ok(steps)
}

impl EvalDyn for Method {
    fn signature(&self) -> Signature {
        let d = var("D");
        let e = var("E");
        let point = complex_type(d.clone());
        match self {
            | Self::Newton => {
                let derivative =
                    DimExpr::Mul(vec![e.clone(), DimExpr::Inverse(Box::new(d.clone()))]);
                let root = TypeExpr::Record(vec![
                    ("point".into(), point.clone()),
                    ("iterations".into(), TypeExpr::Int),
                ]);
                let optional_root = TypeExpr::Variant(vec![
                    ("None".into(), Vec::new()),
                    ("Some".into(), vec![root]),
                ]);
                typed(
                    vec![
                        point.clone(),
                        real(d),
                        real(e.clone()),
                        real(e.clone()),
                        TypeExpr::Int,
                        TypeExpr::Int,
                        function(vec![point.clone()], complex_type(e)),
                        function(vec![point], complex_type(derivative)),
                    ],
                    TypeExpr::Record(vec![
                        ("root".into(), optional_root),
                        ("tested".into(), TypeExpr::Int),
                        ("evaluations".into(), TypeExpr::Int),
                    ]),
                    &["D", "E"],
                )
            },
            | Self::Koenigs => {
                typed(
                    vec![
                        point.clone(),
                        TypeExpr::Complex,
                        point.clone(),
                        real(d.clone()),
                        real(d.clone()),
                        TypeExpr::Int,
                        TypeExpr::Int,
                        function(vec![point.clone()], point.clone()),
                    ],
                    TypeExpr::Record(vec![
                        ("coordinate".into(), point),
                        ("fixed_point_error".into(), real(d)),
                        ("steps".into(), TypeExpr::Int),
                        ("evaluations".into(), TypeExpr::Int),
                    ]),
                    &["D"],
                )
            },
        }
    }

    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let expected = self.signature();
        super::super::arity(args, expected.params.len())?;
        match self {
            | Self::Newton => newton(args, ctx, &expected),
            | Self::Koenigs => koenigs(args, ctx, &expected),
        }
    }
}

fn newton(
    args: &[Value],
    ctx: &CallCtx<'_>,
    expected: &Signature,
) -> Result<Value, OpError> {
    let callbacks = Callbacks::new(ctx, args, expected, 6, count(&args[5], 5)?)?;
    let d = dimension(&callbacks.signatures[0].params[0])?;
    let e = dimension(&callbacks.signatures[0].ret)?;
    let derivative = ratio(&e, &d)?;
    equal(&dimension(&callbacks.signatures[1].params[0])?, &d)?;
    equal(&dimension(&callbacks.signatures[1].ret)?, &derivative)?;
    let axis = scale(&args[1], &d, 1, ctx)?;
    let amplitude = scale(&args[2], &e, 2, ctx)?;
    let derivative_scale = finite(amplitude / axis)?;
    let tolerance = finite(scale(&args[3], &e, 3, ctx)? / amplitude)?;
    if derivative_scale <= 0.0 || tolerance <= 0.0 {
        return Err(domain(
            "Newton derivative scale or normalized tolerance underflowed",
        ));
    }
    let max_steps = checked_steps(&args[4], 4)?;
    let initial = finite_complex(decode(&args[0], &d, 0, ctx)? / axis)?;
    let invalid = Cell::new(false);
    let tested = Cell::new(0);
    let f = |z| {
        let before = callbacks.used.get();
        let result = sample(&callbacks, 0, z, (axis, amplitude), (&d, &e), ctx, &invalid);
        tested.set(tested.get() + (callbacks.used.get() - before));
        result
    };
    let df = |z| {
        sample(
            &callbacks,
            1,
            z,
            (axis, derivative_scale),
            (&d, &derivative),
            ctx,
            &invalid,
        )
    };
    let candidate = native::newton_fractal_basin(f, df, initial, max_steps, tolerance);
    let evaluations = callbacks.used.get();
    callbacks.finish()?;
    crate::cancellation::checkpoint()?;
    if invalid.get() {
        return Err(domain("Newton orbit left the finite coordinate range"));
    }
    let root = candidate
        .map(|(root, iterations)| {
            Ok::<_, OpError>(record(vec![
                ("point", value(finite_complex(root * axis)?, d)),
                ("iterations", Value::int(iterations)),
            ]))
        })
        .transpose()?;
    Ok(record(vec![
        (
            "root",
            Value::Variant(crate::kernel::value::VariantValue {
                tag: if root.is_some() {
                    "Some".into()
                } else {
                    "None".into()
                },
                payload: root.into_iter().collect(),
            }),
        ),
        ("tested", Value::int(tested.get())),
        ("evaluations", Value::int(evaluations)),
    ]))
}

fn koenigs(
    args: &[Value],
    ctx: &CallCtx<'_>,
    expected: &Signature,
) -> Result<Value, OpError> {
    let callbacks = Callbacks::new(ctx, args, expected, 7, count(&args[6], 6)?)?;
    let d = dimension(&callbacks.signatures[0].params[0])?;
    equal(&dimension(&callbacks.signatures[0].ret)?, &d)?;
    let axis = scale(&args[3], &d, 3, ctx)?;
    let tolerance = finite(scale(&args[4], &d, 4, ctx)? / axis)?;
    if tolerance <= 0.0 {
        return Err(domain("normalized fixed-point tolerance underflowed"));
    }
    let max_steps = checked_steps(&args[5], 5)?;
    let fixed = finite_complex(decode(&args[0], &d, 0, ctx)? / axis)?;
    let lambda = decode(&args[1], &DimExpr::Dimless, 1, ctx)?;
    if !(0.0 < lambda.norm() && lambda.norm() < 1.0) {
        return Err(domain(
            "Koenigs coordinate needs a multiplier with magnitude strictly between zero and one",
        ));
    }
    let initial = finite_complex(decode(&args[2], &d, 2, ctx)? / axis)?;
    let invalid = Cell::new(false);
    let f = |z| sample(&callbacks, 0, z, (axis, axis), (&d, &d), ctx, &invalid);
    let residual = (f(fixed) - fixed).norm();
    if !residual.is_finite() || residual > tolerance {
        callbacks.finish()?;
        return Err(domain(
            "declared fixed point does not close within tolerance",
        ));
    }
    let before = callbacks.used.get();
    let coordinate = native::koenigs_coordinate(f, fixed, lambda, initial, max_steps);
    let evaluations = callbacks.used.get();
    callbacks.finish()?;
    crate::cancellation::checkpoint()?;
    if invalid.get() {
        return Err(domain("Koenigs orbit left the finite coordinate range"));
    }
    let coordinate = finite_complex(coordinate.map_err(domain)? * axis)?;
    Ok(record(vec![
        ("coordinate", value(coordinate, d.clone())),
        (
            "fixed_point_error",
            super::super::super::dimensional_float(finite(residual * axis)?, d),
        ),
        ("steps", Value::int(evaluations - before)),
        ("evaluations", Value::int(evaluations)),
    ]))
}

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    register(
        i,
        "num_complex_dynamics_newton_basin_checked",
        "crate::kernels::complex_dynamics::newton_fractal_basin",
        "Tests a Newton basin from an initial coordinate using positive coordinate and value scales, absolute residual tolerance, iteration and callback budgets, and function and derivative callbacks. Reports an optional measured root, the number of tested points and actual callbacks.",
        8,
        Method::Newton,
    )?;
    register(
        i,
        "num_complex_dynamics_koenigs_checked",
        "crate::kernels::complex_dynamics::koenigs_coordinate",
        "Returns a truncated measured Koenigs coordinate estimate with its fixed-point closure error and actual counts. The supplied multiplier must have magnitude strictly between zero and one. Uses positive coordinate scale and closure tolerance, and bounded iterations and callbacks.",
        8,
        Method::Koenigs,
    )
}
