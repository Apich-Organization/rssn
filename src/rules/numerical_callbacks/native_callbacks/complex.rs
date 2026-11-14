//! Native complex contour and analytic-root algorithms with independent units.
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
use super::count;
use super::equal;
use super::finite;
use super::function;
use super::positive;
use super::ratio;
use super::real;
use super::register;
use super::scale;
use super::typed;
use super::var;
use crate::kernel::signature::Exp;
use crate::kernel::signature::Rational;
use crate::kernel::signature::RepExpr;
use crate::kernel::value::ComplexValue;
use crate::kernel::value::DynValue;
use crate::kernels::complex as native;
use num_complex::Complex64;
use num_traits::ToPrimitive;
#[derive(Clone, Copy)]
enum Method {
    Circle,
    Contour,
    Residue,
    Derivative,
    Count,
    Newton,
    Cauchy,
}
fn complex_type(d: DimExpr) -> TypeExpr {
    if d == DimExpr::Dimless {
        TypeExpr::Complex
    } else {
        TypeExpr::Dim(Box::new(TypeExpr::Complex), d)
    }
}
fn dimension(t: &TypeExpr) -> Result<DimExpr, OpError> {
    match t {
        | TypeExpr::Complex => Ok(DimExpr::Dimless),
        | TypeExpr::Dim(inner, d) if matches!(inner.as_ref(), TypeExpr::Complex) => {
            super::super::canonical(d)?;
            Ok(d.clone())
        },
        | _ => Err(super::domain("callback requires a concrete complex scalar")),
    }
}
fn value(
    z: Complex64,
    d: DimExpr,
) -> Value {
    let v = Value::Complex(ComplexValue::Float { re: z.re, im: z.im });
    if d == DimExpr::Dimless {
        v
    } else {
        Value::Dyn(DynValue {
            value: Box::new(v),
            dim: d,
        })
    }
}
fn decode(
    v: &Value,
    d: &DimExpr,
    index: usize,
    ctx: &CallCtx<'_>,
) -> Result<Complex64, OpError> {
    fn payload(v: &Value) -> Result<Complex64, OpError> {
        match v {
            | Value::Dyn(v) => payload(&v.value),
            | Value::Complex(v) => {
                let (re, im) = v.to_f64();
                Ok(Complex64::new(super::finite(re)?, super::finite(im)?))
            },
            | _ => Err(super::domain("complex argument payload is malformed")),
        }
    }
    equal(&dimension(&ctx.value_type(v)?)?, d)?;
    let _ = index;
    payload(v)
}
impl Method {
    const fn start(self) -> usize {
        match self {
            | Self::Contour => 5,
            | Self::Derivative => 4,
            | Self::Cauchy => 7,
            | _ => 6,
        }
    }
}
impl EvalDyn for Method {
    fn signature(&self) -> Signature {
        let d = var("D");
        let e = var("E");
        let cd = complex_type(d.clone());
        let ce = complex_type(e.clone());
        let vec = TypeExpr::Tensor {
            rank: crate::kernel::signature::NatExpr::Lit(1),
            shape: None,
            elem: Box::new(cd.clone()),
        };
        let mut params = match self {
            | Self::Contour => {
                vec![
                    vec,
                    TypeExpr::Int,
                    real(d.clone()),
                    real(e.clone()),
                    TypeExpr::Int,
                ]
            },
            | Self::Derivative => vec![cd.clone(), real(d.clone()), real(e.clone()), TypeExpr::Int],
            | Self::Newton => {
                vec![
                    cd.clone(),
                    real(d.clone()),
                    real(e.clone()),
                    TypeExpr::Float,
                    TypeExpr::Int,
                    TypeExpr::Int,
                ]
            },
            | Self::Cauchy => {
                vec![
                    cd.clone(),
                    real(d.clone()),
                    TypeExpr::Int,
                    TypeExpr::Int,
                    real(d.clone()),
                    real(e.clone()),
                    TypeExpr::Int,
                ]
            },
            | _ => {
                vec![
                    cd.clone(),
                    real(d.clone()),
                    TypeExpr::Int,
                    real(d.clone()),
                    real(e.clone()),
                    TypeExpr::Int,
                ]
            },
        };
        params.push(function(vec![cd.clone()], ce));
        if matches!(self, Self::Newton) {
            params.push(function(
                vec![cd.clone()],
                complex_type(DimExpr::Mul(vec![
                    e.clone(),
                    DimExpr::Inverse(Box::new(d.clone())),
                ])),
            ));
        }
        let ret = match self {
            | Self::Count => TypeExpr::Int,
            | Self::Newton => cd,
            | Self::Derivative => {
                complex_type(DimExpr::Mul(vec![e, DimExpr::Inverse(Box::new(d))]))
            },
            | Self::Cauchy => {
                TypeExpr::Dyn(RepExpr {
                    repr: "complex".into(),
                    dim: None,
                })
            },
            | _ => complex_type(DimExpr::Mul(vec![d, e])),
        };
        typed(params, ret, &["D", "E"])
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
        let callbacks = Callbacks::new(
            ctx,
            args,
            &expected,
            start,
            count(&args[start - 1], start - 1)?,
        )?;
        let d = dimension(&callbacks.signatures[0].params[0])?;
        let e = dimension(&callbacks.signatures[0].ret)?;
        let (xs_at, ys_at) = match self {
            | Self::Contour => (2, 3),
            | Self::Derivative | Self::Newton => (1, 2),
            | Self::Cauchy => (4, 5),
            | _ => (3, 4),
        };
        let xs = scale(&args[xs_at], &d, xs_at, ctx)?;
        let ys = scale(&args[ys_at], &e, ys_at, ctx)?;
        let f = |z: Complex64| {
            callbacks
                .call(0, &[value(z * xs, d.clone())], |v| {
                    decode(&v, &e, start, ctx).map(|v| v / ys)
                })
                .unwrap_or(Complex64::new(f64::NAN, f64::NAN))
        };
        let center = if matches!(self, Self::Contour) {
            Complex64::new(0.0, 0.0)
        } else {
            decode(&args[0], &d, 0, ctx)? / xs
        };
        let radius = if matches!(
            self,
            Self::Circle | Self::Residue | Self::Count | Self::Cauchy
        ) {
            scale(&args[1], &d, 1, ctx)? / xs
        } else {
            0.0
        };
        let samples = if matches!(self, Self::Circle | Self::Residue | Self::Count) {
            count(&args[2], 2)?
        } else if matches!(self, Self::Cauchy) {
            count(&args[3], 3)?
        } else {
            0
        };
        let mut order = 0_u32;
        let z = match self {
            | Self::Circle => native::circle_integral(f, center, radius, samples),
            | Self::Residue => native::residue(f, center, radius, samples),
            | Self::Derivative => native::complex_derivative(&f, center),
            | Self::Count => {
                let out = native::count_zeros_poles(f, center, radius, samples);
                callbacks.finish()?;
                return out.map(Value::int).ok_or_else(|| {
                    super::domain("argument-principle contour is singular or undersampled")
                });
            },
            | Self::Cauchy => {
                let Value::Int(n) = &args[2] else {
                    return Err(super::domain("Cauchy derivative order must be an integer"));
                };
                order = n.to_u32().filter(|n| *n <= 170).ok_or_else(|| {
                    super::domain("Cauchy derivative order must be between zero and170")
                })?;
                native::cauchy_derivative(f, center, order, radius, samples)
            },
            | Self::Newton => {
                let derivative = ratio(&e, &d)?;
                equal(&dimension(&callbacks.signatures[1].params[0])?, &d)?;
                equal(&dimension(&callbacks.signatures[1].ret)?, &derivative)?;
                let df = |z| {
                    callbacks
                        .call(1, &[value(z * xs, d.clone())], |v| {
                            decode(&v, &derivative, start + 1, ctx).map(|v| v * xs / ys)
                        })
                        .unwrap_or(Complex64::new(f64::NAN, f64::NAN))
                };
                let out = native::newton(
                    f,
                    df,
                    center,
                    positive(&args[3], 3, ctx)?,
                    count(&args[4], 4)?,
                );
                callbacks.finish()?;
                return out
                    .map(|z| value(z * xs, d))
                    .ok_or_else(|| super::domain("complex Newton did not converge"));
            },
            | Self::Contour => {
                let Value::Tensor(t) = &args[0] else {
                    return Err(super::domain("contour path must be a complex vector"));
                };
                if t.shape != [t.data.len()] || t.data.len() < 2 {
                    return Err(super::domain("contour requires at least two path points"));
                }
                equal(&dimension(&t.elem)?, &d)?;
                let mut points = Vec::with_capacity(t.data.len());
                for v in &t.data {
                    let raw = match v {
                        | Value::Complex(c) => {
                            let (re, im) = c.to_f64();
                            Complex64::new(finite(re)?, finite(im)?)
                        },
                        | _ => decode(v, &d, 0, ctx)?,
                    };
                    points.push(raw / xs);
                }
                let panels = count(&args[1], 1)?;
                super::bounded_size(panels, points.len())?;
                native::contour_integral(f, &points, panels)
            },
        };
        callbacks.finish()?;
        let (z, dim) = match self {
            | Self::Derivative => (z * (ys / xs), ratio(&e, &d)?),
            | Self::Cauchy => {
                let coefficient = ys
                    / xs.powi(
                        i32::try_from(order)
                            .map_err(|_| super::domain("derivative order overflow"))?,
                    );
                let dim =
                    super::super::expression(&super::super::canonical(&DimExpr::Mul(vec![
                        e,
                        DimExpr::Pow(Box::new(d), Exp::Rat(Rational::new(-i64::from(order), 1))),
                    ]))?)?;
                (z * coefficient, dim)
            },
            | _ => {
                (
                    z * ys * xs,
                    super::super::expression(&super::super::canonical(&DimExpr::Mul(vec![d, e]))?)?,
                )
            },
        };
        finite(z.re)?;
        finite(z.im)?;
        Ok(value(z, dim))
    }
}
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    for (name, path, m) in [
        ("circle_integral", "circle_integral", Method::Circle),
        ("contour_integral", "contour_integral", Method::Contour),
        ("residue", "residue", Method::Residue),
        ("derivative", "complex_derivative", Method::Derivative),
        ("count_zeros_poles", "count_zeros_poles", Method::Count),
        ("newton", "newton", Method::Newton),
        ("cauchy_derivative", "cauchy_derivative", Method::Cauchy),
    ] {
        register(
            i,
            &format!("num_complex_{name}_scaled"),
            &format!("crate::kernels::complex::{path}"),
            "Native complex analysis with positive coordinate/value scales, independent D/E physical units and a shared actual callback budget. Circle/residue/count: center,radius,n,xscale,yscale,max_evaluations,f. Contour: path,panels,xscale,yscale,max_evaluations,f. Derivative: z,xscale,yscale,max_evaluations,f. Newton: z,xscale,yscale,tolerance,max_iterations,max_evaluations,f,df. Cauchy: z,radius,order,n,xscale,yscale,max_evaluations,f; dynamic result explicitly retains exact E/D^order units. Analyticity and contour isolation are caller assumptions; native singular/nonconvergence outcomes remain errors.",
            m.start() as u8
                + if matches!(m, Method::Newton) {
                    2
                } else {
                    1
                },
            m,
        )?;
    }
    Ok(())
}
