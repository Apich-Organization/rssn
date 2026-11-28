//! Unit-aware complex finite differences, shared across conformality tests.
use crate::kernels::complex::differentiation as native;

use super::super::CallCtx;
use super::super::Callbacks;
use super::super::DimExpr;
use super::super::EvalDyn;
use super::super::Installer;
use super::super::OpError;
use super::super::RuleError;
use super::super::Signature;
use super::super::TypeExpr;
use super::super::Value;
use super::super::count;
use super::super::domain;
use super::super::finite;
use super::super::function;
use super::super::ratio;
use super::super::real;
use super::super::record;
use super::super::register;
use super::super::scale;
use super::super::typed;
use super::super::var;
use super::complex_type;
use super::decode;
use super::dimension;
use super::value;

const D0: DimExpr = DimExpr::Dimless;

#[derive(Clone, Copy)]
enum Method {
    Dz,
    DzBar,
    CauchyRiemann,
    Laplacian,
    Polar,
    Conformal,
    ScaleRotation,
    ComplexStep,
}

impl Method {
    const fn count(self) -> u8 {
        match self {
            | Self::Polar | Self::CauchyRiemann | Self::Conformal => 7,
            | _ => 5,
        }
    }

    const fn uses_tolerance(self) -> bool {
        matches!(self, Self::CauchyRiemann | Self::Conformal)
    }

    fn derivative_dimension(
        self,
        input: DimExpr,
        output: DimExpr,
    ) -> DimExpr {
        let denominator = if matches!(self, Self::Laplacian) {
            DimExpr::Mul(vec![input.clone(), input])
        } else {
            input
        };
        DimExpr::Mul(vec![output, DimExpr::Inverse(Box::new(denominator))])
    }

    fn cartesian(
        self,
        args: &[Value],
        ctx: &CallCtx<'_>,
        callbacks: &Callbacks<'_, '_>,
        coordinate: &DimExpr,
        ordinate: &DimExpr,
    ) -> Result<Value, OpError> {
        let value_index = if self.uses_tolerance() {
            3
        } else {
            2
        };
        let y_scale = scale(&args[value_index], ordinate, value_index, ctx)?;

        let step = scale(&args[1], coordinate, 1, ctx)?;
        let input = |z| value(z, coordinate.clone());
        let callback = |z| {
            callbacks
                .call(0, &[input(z)], |arg| {
                    decode(&arg, ordinate, 0, ctx).map(|v| v / y_scale)
                })
                .unwrap_or(native::NAN)
        };

        let estimate = || {
            domain("complex derivative needs distinct representable coordinates and finite samples")
        };
        if matches!(self, Self::ComplexStep) {
            let point = super::super::scalar(&args[0], coordinate, 0, ctx)?;
            let result =
                y_scale * native::complex_step(callback, point, step).ok_or_else(estimate)?;
            let dim = ratio(ordinate, coordinate)?;
            return Ok(super::super::super::dimensional_float(finite(result)?, dim));
        }
        let point = decode(&args[0], coordinate, 0, ctx)?;
        if matches!(self, Self::Laplacian) {
            let result = y_scale * native::laplacian(callback, point, step).ok_or_else(estimate)?;
            let dim = ratio(
                ordinate,
                &DimExpr::Mul(vec![coordinate.clone(), coordinate.clone()]),
            )?;
            finite(result.re)?;
            finite(result.im)?;
            return Ok(value(result, dim));
        }
        let (x_scale, tolerance) = if self.uses_tolerance() {
            let x_scale = scale(&args[2], coordinate, 2, ctx)?;
            let tolerance = super::super::scalar(&args[4], &D0, 4, ctx)?;
            if tolerance < 0.0 {
                return Err(domain("conformality tolerance must be nonnegative"));
            }
            (x_scale, tolerance)
        } else {
            (1.0, 0.0)
        };
        let result = native::wirtinger(callback, point, step).ok_or_else(estimate)?;
        match self {
            | Self::CauchyRiemann => return Ok(Value::Bool(result.residual(x_scale) <= tolerance)),
            | Self::Conformal => {
                return Ok(Value::Bool(
                    result.residual(x_scale) <= tolerance
                        && (result.dz * x_scale).norm() > tolerance,
                ));
            },
            | _ => {},
        }
        let derivative = if matches!(self, Self::DzBar) {
            result.dzbar
        } else {
            result.dz
        };
        let rotation = derivative.arg();
        let derivative = derivative * y_scale;
        let dim = ratio(ordinate, coordinate)?;
        if matches!(self, Self::ScaleRotation) {
            Ok(record(vec![
                (
                    "scale",
                    super::super::super::dimensional_float(finite(derivative.norm())?, dim),
                ),
                ("rotation", Value::Float(finite(rotation)?)),
            ]))
        } else {
            finite(derivative.re)?;
            finite(derivative.im)?;
            Ok(value(derivative, dim))
        }
    }
}

fn polar(
    args: &[Value],
    ctx: &CallCtx<'_>,
    callbacks: &Callbacks<'_, '_>,
    coord: &DimExpr,
    ordinate: &DimExpr,
) -> Result<Value, OpError> {
    let y_scale = scale(&args[4], ordinate, 4, ctx)?;
    let radius = super::super::scalar(&args[0], coord, 0, ctx)?;
    let angle = super::super::scalar(&args[1], &D0, 1, ctx)?;
    let radial_step = scale(&args[2], coord, 2, ctx)?;
    let angular_step = scale(&args[3], &D0, 3, ctx)?;
    let callback = |r, theta| {
        let radial = super::super::super::dimensional_float(r, coord.clone());
        callbacks
            .call(0, &[radial, Value::Float(theta)], |arg| {
                decode(&arg, ordinate, 0, ctx).map(|v| v / y_scale)
            })
            .unwrap_or(native::NAN)
    };
    let (radial, angular) = native::polar(callback, radius, angle, radial_step, angular_step)
        .ok_or_else(|| {
            domain("polar stencil needs distinct finite samples and nonnegative radius")
        })?;
    Ok(record(vec![
        (
            "radial",
            super::super::super::dimensional_float(finite(radial * y_scale)?, ordinate.clone()),
        ),
        (
            "angular",
            super::super::super::dimensional_float(finite(angular * y_scale)?, ordinate.clone()),
        ),
    ]))
}

impl EvalDyn for Method {
    fn signature(&self) -> Signature {
        let d = var("D");
        let e = var("E");
        let (args, result) = if matches!(self, Self::Polar) {
            (
                vec![
                    real(d.clone()),
                    TypeExpr::Float,
                    real(d.clone()),
                    TypeExpr::Float,
                    real(e.clone()),
                    TypeExpr::Int,
                    function(vec![real(d), TypeExpr::Float], complex_type(e.clone())),
                ],
                TypeExpr::Record(vec![
                    ("radial".into(), real(e.clone())),
                    ("angular".into(), real(e)),
                ]),
            )
        } else {
            let derivative = self.derivative_dimension(d.clone(), e.clone());
            let point = if matches!(self, Self::ComplexStep) {
                real(d.clone())
            } else {
                complex_type(d.clone())
            };
            let mut args = vec![point, real(d.clone())];
            if self.uses_tolerance() {
                args.extend([real(d.clone()), real(e.clone()), TypeExpr::Float]);
            } else {
                args.push(real(e.clone()));
            }
            args.push(TypeExpr::Int);
            args.push(function(vec![complex_type(d)], complex_type(e)));
            let result = match self {
                | Self::CauchyRiemann | Self::Conformal => TypeExpr::Bool,
                | Self::ScaleRotation => {
                    TypeExpr::Record(vec![
                        ("scale".into(), real(derivative)),
                        ("rotation".into(), TypeExpr::Float),
                    ])
                },
                | Self::ComplexStep => real(derivative),
                | _ => complex_type(derivative),
            };
            (args, result)
        };
        typed(args, result, &["D", "E"])
    }

    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let signature = self.signature();
        super::super::arity(args, signature.params.len())?;
        let budget = usize::from(self.count()) - 2;
        let callbacks = Callbacks::new(
            ctx,
            args,
            &signature,
            usize::from(self.count()) - 1,
            count(&args[budget], budget)?,
        )?;
        let output = dimension(&callbacks.signatures[0].ret)?;
        let result = if matches!(self, Self::Polar) {
            let input = super::super::scalar_dimension(&callbacks.signatures[0].params[0])?;
            polar(args, ctx, &callbacks, &input, &output)
        } else {
            let input = dimension(&callbacks.signatures[0].params[0])?;
            self.cartesian(args, ctx, &callbacks, &input, &output)
        };
        callbacks.finish()?;
        result
    }
}

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    for (name, path, method, description) in [
        (
            "num_complex_wirtinger_dz_checked",
            "crate::kernels::complex::wirtinger_dz_numeric",
            Method::Dz,
            "Estimates the Wirtinger derivative with respect to z.",
        ),
        (
            "num_complex_wirtinger_dzbar_checked",
            "crate::kernels::complex::wirtinger_dzbar_numeric",
            Method::DzBar,
            "Estimates the Wirtinger derivative with respect to conjugate z.",
        ),
        (
            "num_complex_cauchy_riemann_checked",
            "crate::kernels::complex::check_cauchy_riemann_numeric",
            Method::CauchyRiemann,
            "Checks the relative Cauchy–Riemann residual in the normalized reference frame.",
        ),
        (
            "num_complex_laplacian_checked",
            "crate::kernels::complex::complex_laplacian_numeric",
            Method::Laplacian,
            "Estimates the Cartesian Laplacian in callback-value units over coordinate units squared.",
        ),
        (
            "num_complex_cauchy_riemann_polar_checked",
            "crate::kernels::complex::cauchy_riemann_polar_numeric",
            Method::Polar,
            "Returns radial and angular polar Cauchy–Riemann residuals with separate radial and angular steps.",
        ),
        (
            "num_conformal_numeric_checked",
            "crate::kernels::conformal::is_conformal_numeric",
            Method::Conformal,
            "Checks the relative Cauchy–Riemann residual and nonzero normalized magnification using one stencil.",
        ),
        (
            "num_conformal_scale_rotation_checked",
            "crate::kernels::conformal::conformal_scale_and_rotation",
            Method::ScaleRotation,
            "Returns measured local magnification and dimensionless rotation in radians.",
        ),
        (
            "num_complex_step_derivative_checked",
            "crate::kernels::complex::complex_step_derivative",
            Method::ComplexStep,
            "Estimates a real analytic derivative from the imaginary part of a single sample.",
        ),
    ] {
        let arguments = match method {
            | Method::Polar => {
                "radius, angle, radial_step, angular_step, value_scale, budget, callback"
            },
            | Method::ComplexStep => "x, step, value_scale, budget, callback",
            | Method::CauchyRiemann | Method::Conformal => {
                "point, step, coordinate_scale, value_scale, tolerance, budget, callback"
            },
            | _ => "point, step, value_scale, budget, callback",
        };
        let doc = format!(
            "`{name}({arguments})`: {description} Invalid or collapsed stencils and exhaustion of the total callback budget produce errors."
        );
        register(i, name, path, &doc, method.count(), method)?;
    }
    Ok(())
}
