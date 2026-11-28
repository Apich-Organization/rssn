//! Checked adaptive contour callbacks using portable oriented geometry.

use crate::kernel::signature::NatExpr;
use crate::kernels::contour_deform as native;
use crate::kernels::contour_deform::ContourQuadratureError;
use crate::kernels::quadrature::GkError;
use crate::rules::numerical_abi::codecs::checked_deformation;
use crate::rules::numerical_abi::codecs::contour_complex_positions;
use crate::rules::numerical_abi::codecs::contour_deformation_preflight;
use crate::rules::numerical_abi::codecs::contour_path_type;
use crate::rules::numerical_abi::codecs::decode_contour_path;
use num_complex::Complex64;

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
use super::complex::complex_type;
use super::complex::decode;
use super::complex::dimension;
use super::complex::value;
use super::count;
use super::domain;
use super::finite;
use super::function;
use super::real;
use super::record;
use super::register;
use super::scale;
use super::typed;
use super::var;
use super::vector;

#[derive(Clone, Copy)]
enum Method {
    Segment,
    Path,
    Residue,
    PrincipalValue,
    Deformed,
}

impl Method {
    const fn shape_count(self) -> usize {
        match self {
            | Self::PrincipalValue => 3,
            | Self::Deformed => 6,
            | _ => 1,
        }
    }
}

fn failure(error: ContourQuadratureError) -> OpError {
    match error {
        | ContourQuadratureError::Quadrature(GkError::PanelBudget) => OpError::BudgetExceeded,
        | ContourQuadratureError::Quadrature(GkError::Cancelled) => {
            crate::cancellation::Cancelled.into()
        },
        | error => domain(error.to_string()),
    }
}

impl EvalDyn for Method {
    fn signature(&self) -> Signature {
        let d = var("D");
        let e = var("E");
        let integral = DimExpr::Mul(vec![d.clone(), e.clone()]);
        let mut args = match self {
            | Self::Segment | Self::Path => vec![contour_path_type()],
            | Self::Residue => vec![real(d.clone())],
            | Self::Deformed => {
                vec![
                    complex_type(d.clone()),
                    complex_type(d.clone()),
                    contour_vector_type(complex_type(d.clone())),
                    contour_vector_type(complex_type(d.clone())),
                    real(d.clone()),
                    contour_vector_type(TypeExpr::Float),
                ]
            },
            | Self::PrincipalValue => {
                vec![
                    real(d.clone()),
                    real(d.clone()),
                    TypeExpr::Tensor {
                        rank: NatExpr::Lit(1),
                        shape: None,
                        elem: Box::new(real(d.clone())),
                    },
                ]
            },
        };
        args.extend([
            real(d.clone()),
            real(e.clone()),
            real(integral.clone()),
            TypeExpr::Int,
            TypeExpr::Int,
        ]);
        let coord = if matches!(self, Self::PrincipalValue) {
            real(d)
        } else {
            complex_type(d)
        };
        args.push(function(vec![coord], complex_type(e)));
        typed(
            args,
            TypeExpr::Record(vec![
                ("value".into(), complex_type(integral.clone())),
                ("error".into(), real(integral)),
                ("evaluations".into(), TypeExpr::Int),
                ("converged".into(), TypeExpr::Bool),
            ]),
            &["D", "E"],
        )
    }

    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let expected = self.signature();
        super::arity(args, expected.params.len())?;
        let start = self.shape_count();
        let callback = start + 5;
        let callbacks = Callbacks::new(
            ctx,
            args,
            &expected,
            callback,
            count(&args[start + 4], start + 4)?,
        )?;
        let d = if matches!(self, Self::PrincipalValue) {
            super::scalar_dimension(&callbacks.signatures[0].params[0])?
        } else {
            dimension(&callbacks.signatures[0].params[0])?
        };
        let e = dimension(&callbacks.signatures[0].ret)?;
        let integral = super::expression(&super::canonical(&DimExpr::Mul(vec![
            d.clone(),
            e.clone(),
        ]))?)?;
        let x_scale = scale(&args[start], &d, start, ctx)?;
        let y_scale = scale(&args[start + 1], &e, start + 1, ctx)?;
        let unit = finite(x_scale * y_scale)?;
        if unit <= 0.0 {
            return Err(domain("contour integral reference scale must be positive"));
        }
        let tol = scale(&args[start + 2], &integral, start + 2, ctx)? / unit;
        let panels = count(&args[start + 3], start + 3)?;
        let apply = |arg| {
            callbacks
                .call(0, &[arg], |input| {
                    decode(&input, &e, callback, ctx).map(|z| z / y_scale)
                })
                .unwrap_or(Complex64::new(f64::NAN, f64::NAN))
        };
        let complex = |z| apply(value(z * x_scale, d.clone()));
        let result = match self {
            | Self::Path | Self::Segment => {
                let path = decode_contour_path(&args[0], &d, x_scale)?;
                if matches!(self, Self::Segment) {
                    let [segment] = path.as_slice() else {
                        return Err(domain("one segment is required"));
                    };
                    native::integrate_segment_checked(complex, *segment, tol, panels)
                } else {
                    native::integrate_contour_checked(complex, &path, tol, panels)
                }
            },
            | Self::Deformed => {
                native::integrate_contour_checked(
                    complex,
                    &obstacle_path(args, &d, x_scale, ctx)?,
                    tol,
                    panels,
                )
            },
            | Self::Residue => {
                let radius = scale(&args[0], &d, 0, ctx)? / x_scale;
                native::residue_at_infinity_checked(complex, radius, tol, panels)
            },
            | Self::PrincipalValue => {
                super::check_dimension(ctx, &args[0], &d, 0)?;
                super::check_dimension(ctx, &args[1], &d, 1)?;
                let a = super::numeric(&args[0], 0)? / x_scale;
                let b = super::numeric(&args[1], 1)? / x_scale;
                let poles: Vec<_> = vector(&args[2], &d, None, 2, ctx)?
                    .into_iter()
                    .map(|p| p / x_scale)
                    .collect();
                native::cauchy_principal_value_checked(
                    |x| apply(super::super::dimensional_float(x * x_scale, d.clone())),
                    a,
                    b,
                    &poles,
                    tol,
                    panels,
                )
            },
        };
        callbacks.finish()?;
        let out = result.map_err(failure)?;
        let estimate = out.value * unit;
        finite(estimate.re)?;
        finite(estimate.im)?;
        Ok(record(vec![
            ("value", value(estimate, integral.clone())),
            (
                "error",
                super::super::dimensional_float(finite(out.error * unit)?, integral),
            ),
            ("evaluations", Value::int(out.evaluations)),
            ("converged", Value::Bool(out.converged)),
        ]))
    }
}

fn obstacle_path(
    args: &[Value],
    dim: &DimExpr,
    x_scale: f64,
    ctx: &CallCtx<'_>,
) -> Result<Vec<native::ContourSegment>, OpError> {
    contour_deformation_preflight(args)?;
    let [start, end, poles, origins, clearance, angles, ..] = args else {
        return Err(domain("incomplete deformed contour geometry"));
    };
    let start = decode(start, dim, 0, ctx)? / x_scale;
    let end = decode(end, dim, 1, ctx)? / x_scale;
    let poles = contour_complex_positions(poles, dim, 2, x_scale)?;
    let origins = contour_complex_positions(origins, dim, 3, x_scale)?;
    let clearance = scale(clearance, dim, 4, ctx)? / x_scale;
    let angles = vector(angles, &DimExpr::Dimless, None, 5, ctx)?;
    checked_deformation(start, end, &poles, &origins, clearance, &angles)
}

fn contour_vector_type(elem: TypeExpr) -> TypeExpr {
    TypeExpr::Tensor {
        rank: NatExpr::Lit(1),
        shape: None,
        elem: Box::new(elem),
    }
}

struct ValueOnly(Method);

impl EvalDyn for ValueOnly {
    fn signature(&self) -> Signature {
        let mut signature = self.0.signature();
        let start = self.0.shape_count();
        signature.params.remove(start + 3);
        signature.ret = complex_type(DimExpr::Mul(vec![var("D"), var("E")]));
        signature
    }

    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let signature = self.signature();
        super::arity(args, signature.params.len())?;
        let start = self.0.shape_count();
        let callback = start + 4;
        let callbacks = Callbacks::new(
            ctx,
            args,
            &signature,
            callback,
            count(&args[start + 3], start + 3)?,
        )?;
        let coord = &callbacks.signatures[0].params[0];
        let d = if matches!(self.0, Method::PrincipalValue) {
            super::scalar_dimension(coord)?
        } else {
            dimension(coord)?
        };
        let e = dimension(&callbacks.signatures[0].ret)?;
        let integral = super::expression(&super::canonical(&DimExpr::Mul(vec![
            d.clone(),
            e.clone(),
        ]))?)?;
        let x_scale = scale(&args[start], &d, start, ctx)?;
        let y_scale = scale(&args[start + 1], &e, start + 1, ctx)?;
        let unit = finite(x_scale * y_scale)?;
        if unit <= 0.0 {
            return Err(domain("contour integral reference scale must be positive"));
        }
        let tol = scale(&args[start + 2], &integral, start + 2, ctx)? / unit;
        let apply = |input| {
            callbacks
                .call(0, &[input], |output| {
                    decode(&output, &e, callback, ctx).map(|z| z / y_scale)
                })
                .unwrap_or(Complex64::new(f64::NAN, f64::NAN))
        };
        let complex = |z| apply(value(z * x_scale, d.clone()));
        let result = match self.0 {
            | Method::Segment => {
                let path = decode_contour_path(&args[0], &d, x_scale)?;
                let [segment] = path.as_slice() else {
                    return Err(domain("one segment is required"));
                };
                native::integrate_segment(&complex, *segment, tol)
            },
            | Method::Path => {
                let path = decode_contour_path(&args[0], &d, x_scale)?;
                native::integrate_contour(&complex, &path, tol)
            },
            | Method::PrincipalValue => {
                let a = super::scalar(&args[0], &d, 0, ctx)? / x_scale;
                let b = super::scalar(&args[1], &d, 1, ctx)? / x_scale;
                let poles: Vec<_> = vector(&args[2], &d, None, 2, ctx)?
                    .into_iter()
                    .map(|pole| pole / x_scale)
                    .collect();
                native::cauchy_principal_value(
                    |x| apply(super::super::dimensional_float(x * x_scale, d.clone())),
                    a,
                    b,
                    &poles,
                    tol,
                )
            },
            | _ => return Err(domain("invalid value-only contour entrypoint")),
        };
        callbacks.finish()?;
        crate::cancellation::checkpoint()?;
        let estimate = result * unit;
        finite(estimate.re)?;
        finite(estimate.im)?;
        Ok(value(estimate, integral))
    }
}

struct FixedResidue;

impl EvalDyn for FixedResidue {
    fn signature(&self) -> Signature {
        let d = var("D");
        let e = var("E");
        let integral = DimExpr::Mul(vec![d.clone(), e.clone()]);
        typed(
            vec![
                real(d.clone()),
                TypeExpr::Int,
                real(d.clone()),
                real(e.clone()),
                TypeExpr::Int,
                function(vec![complex_type(d)], complex_type(e)),
            ],
            complex_type(integral),
            &["D", "E"],
        )
    }

    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let signature = self.signature();
        super::arity(args, signature.params.len())?;
        let samples = count(&args[1], 1)?.max(32);
        let budget = count(&args[4], 4)?;
        if samples > 4096 || budget < samples {
            return Err(OpError::BudgetExceeded);
        }
        let callbacks = Callbacks::new(ctx, args, &signature, 5, budget)?;
        let d = dimension(&callbacks.signatures[0].params[0])?;
        let e = dimension(&callbacks.signatures[0].ret)?;
        let radius = scale(&args[0], &d, 0, ctx)?;
        let x_scale = scale(&args[2], &d, 2, ctx)?;
        let y_scale = scale(&args[3], &e, 3, ctx)?;
        let normalized_radius = radius / x_scale;
        crate::kernels::complex_roots::validate_circle(Complex64::ZERO, normalized_radius, samples)
            .map_err(domain)?;
        let scale = finite(x_scale * y_scale)?;
        if scale <= 0.0 {
            return Err(domain("contour integral reference scale must be positive"));
        }
        let apply = |z| {
            callbacks
                .call(0, &[value(z * x_scale, d.clone())], |input| {
                    decode(&input, &e, 5, ctx).map(|z| z / y_scale)
                })
                .unwrap_or(Complex64::new(f64::NAN, f64::NAN))
        };
        let result = native::residue_at_infinity_numeric(apply, normalized_radius, samples) * scale;
        callbacks.finish()?;
        crate::cancellation::checkpoint()?;
        finite(result.re)?;
        finite(result.im)?;
        Ok(value(
            result,
            super::expression(&super::canonical(&DimExpr::Mul(vec![d, e]))?)?,
        ))
    }
}

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    for (name, path, method, description) in [
        (
            "num_contour_integrate_segment_checked",
            "crate::kernels::contour_deform::integrate_segment_checked",
            Method::Segment,
            "Estimates one oriented line or arc from a single-segment portable contour record.",
        ),
        (
            "num_contour_integrate_checked",
            "crate::kernels::contour_deform::integrate_contour_checked",
            Method::Path,
            "Estimates an oriented line-and-arc integral along a portable contour record.",
        ),
        (
            "num_contour_residue_infinity_checked",
            "crate::kernels::contour_deform::residue_at_infinity_checked",
            Method::Residue,
            "Adaptively estimates the residue at infinity on a positive circle enclosing the finite singularities supplied by the caller.",
        ),
        (
            "num_contour_principal_value_checked",
            "crate::kernels::contour_deform::cauchy_principal_value_checked",
            Method::PrincipalValue,
            "Estimates a Cauchy principal value on two finite endpoints with a vector of strictly interior simple poles.",
        ),
        (
            "num_contour_integrate_deformed_checked",
            "crate::kernels::contour_deform::integrate_deformed",
            Method::Deformed,
            "Builds a verified deformation from measured endpoints, poles, ray origins, positive clearance and bare ray angles, then integrates without sampling unsafe geometry.",
        ),
    ] {
        let description = format!(
            "{description} Pass positive coordinate and callback-value scales, absolute tolerance, panel limit per component, total callback limit, then the callback. The result includes the value, heuristic absolute error, actual evaluations, and convergence status."
        );
        register(
            i,
            name,
            path,
            &description,
            (method.shape_count() + 6) as u8,
            method,
        )?;
    }
    for (name, path, method, doc) in [
        (
            "num_contour_integrate_segment",
            "crate::kernels::contour_deform::integrate_segment",
            Method::Segment,
            "Integrates one oriented measured segment, returning the value with the legacy tolerance and quadrature defaults. Pass path, coordinate scale, integrand scale, tolerance, callback budget and callback.",
        ),
        (
            "num_contour_integrate",
            "crate::kernels::contour_deform::integrate_contour",
            Method::Path,
            "Integrates an oriented measured path, returning the value with the legacy tolerance and quadrature defaults. Pass path, coordinate scale, integrand scale, tolerance, callback budget and callback.",
        ),
        (
            "num_contour_principal_value",
            "crate::kernels::contour_deform::cauchy_principal_value",
            Method::PrincipalValue,
            "Returns the legacy value-only principal value with interior pole filtering. Reversed intervals return zero. Pass endpoints, poles, coordinate scale, integrand scale, tolerance, callback budget and callback.",
        ),
    ] {
        register(
            i,
            name,
            path,
            doc,
            (method.shape_count() + 5) as u8,
            ValueOnly(method),
        )?;
    }
    register(
        i,
        "num_contour_residue_infinity_fixed",
        "crate::kernels::contour_deform::residue_at_infinity_numeric",
        "Estimates residue at infinity with exactly max(samples, 32) circular samples. Pass radius, requested samples, coordinate reference, callback-value reference, callback budget and callback.",
        6,
        FixedResidue,
    )?;
    Ok(())
}
