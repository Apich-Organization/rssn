//! Matrix-free native Newton-CG with exact user-supplied derivative actions.
use super::super::CallCtx;
use super::super::DimExpr;
use super::super::EvalDyn;
use super::super::Installer;
use super::super::OpError;
use super::super::RuleError;
use super::super::Signature;
use super::super::TypeExpr;
use super::super::Value;
use super::super::common::register;
use super::super::evolution::output_tensor;
use super::super::evolution::tensor_type;
use super::super::evolution::vector;
use super::super::native_callbacks::Callbacks;
use super::super::native_callbacks::count;
use super::super::native_callbacks::domain;
use super::super::native_callbacks::finite;
use super::super::native_callbacks::function;
use super::super::native_callbacks::real;
use super::super::native_callbacks::record;
use super::super::native_callbacks::scalar;
use super::super::native_callbacks::scale;
use super::super::native_callbacks::typed;
use super::super::native_callbacks::var;
use super::super::native_callbacks::vec_value;
use crate::kernels::newton_cg;
use crate::rules::numerical_abi::native_object_api as api;
#[derive(Clone, Copy)]
struct Newton;
impl EvalDyn for Newton {
    fn signature(&self) -> Signature {
        let d = var("D");
        let e = var("E");
        let v = tensor_type(1, DimExpr::Dimless);
        let config = TypeExpr::Record(vec![
            ("gradient_tolerance".into(), TypeExpr::Float),
            ("max_iterations".into(), TypeExpr::Int),
            ("max_cg_iterations".into(), TypeExpr::Int),
            ("cg_relative_tolerance".into(), TypeExpr::Float),
            ("initial_radius".into(), TypeExpr::Float),
            ("max_radius".into(), TypeExpr::Float),
            ("min_radius".into(), TypeExpr::Float),
            ("acceptance_threshold".into(), TypeExpr::Float),
        ]);
        let derivative = DimExpr::Mul(vec![e.clone(), DimExpr::Inverse(Box::new(d.clone()))]);
        let mut fields = vec![
            ("x".into(), tensor_type(1, d.clone())),
            ("value".into(), real(e.clone())),
            ("gradient_norm".into(), real(derivative)),
            ("normalized_gradient_norm".into(), TypeExpr::Float),
            ("status".into(), TypeExpr::String),
        ];
        for name in [
            "iterations",
            "objective_evaluations",
            "gradient_evaluations",
            "hessian_vector_evaluations",
            "negative_curvature_steps",
            "boundary_steps",
            "accepted_steps",
        ] {
            fields.push((name.into(), TypeExpr::Int));
        }
        typed(
            vec![
                tensor_type(1, d.clone()),
                real(d),
                real(e),
                config,
                TypeExpr::Int,
                function(vec![v.clone()], TypeExpr::Float),
                function(vec![v.clone()], v.clone()),
                function(vec![v.clone(), v.clone()], v),
            ],
            TypeExpr::Record(fields),
            &["D", "E"],
        )
    }

    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let callbacks = Callbacks::new(ctx, args, &self.signature(), 5, count(&args[4], 4)?)?;
        let TypeExpr::Tensor { elem, .. } = ctx.value_type(&args[0])? else {
            return Err(domain("Newton-CG point requires a vector"));
        };
        let d = super::super::scalar_dimension(&elem)?;
        let e = super::super::scalar_dimension(&ctx.value_type(&args[2])?)?;
        let xs = scale(&args[1], &d, 1, ctx)?;
        let ys = scale(&args[2], &e, 2, ctx)?;
        let x = vector(&args[0], &d, None, 0, ctx)?
            .into_iter()
            .map(|v| v / xs)
            .collect::<Vec<_>>();
        if x.is_empty() || x.len() > 4096 || x.iter().any(|v| !v.is_finite()) {
            return Err(domain("Newton-CG point dimension must be1..4096"));
        }
        let cfg = &args[3];
        let number = |name| scalar(api::field(cfg, name)?, &DimExpr::Dimless, 3, ctx);
        let options = newton_cg::NewtonCgOptions {
            gradient_tolerance: number("gradient_tolerance")?,
            max_iterations: count(api::field(cfg, "max_iterations")?, 3)?,
            max_cg_iterations: count(api::field(cfg, "max_cg_iterations")?, 3)?,
            cg_relative_tolerance: number("cg_relative_tolerance")?,
            initial_radius: number("initial_radius")?,
            max_radius: number("max_radius")?,
            min_radius: number("min_radius")?,
            acceptance_threshold: number("acceptance_threshold")?,
        };
        if options.max_cg_iterations > 4096 {
            return Err(domain("CG action limit exceeds4096"));
        }
        let f = |x: &[f64]| {
            callbacks
                .call(0, &[vec_value(x, &DimExpr::Dimless)], |v| {
                    scalar(&v, &DimExpr::Dimless, 5, ctx)
                })
                .unwrap_or(f64::NAN)
        };
        let g = |x: &[f64], out: &mut [f64]| {
            let v = callbacks.call(1, &[vec_value(x, &DimExpr::Dimless)], |v| {
                vector(&v, &DimExpr::Dimless, Some(out.len()), 6, ctx)
            });
            if let Some(v) = v {
                out.copy_from_slice(&v);
            } else {
                out.fill(f64::NAN);
            }
        };
        let h = |x: &[f64], direction: &[f64], out: &mut [f64]| {
            let v = callbacks.call(
                2,
                &[
                    vec_value(x, &DimExpr::Dimless),
                    vec_value(direction, &DimExpr::Dimless),
                ],
                |v| vector(&v, &DimExpr::Dimless, Some(out.len()), 7, ctx),
            );
            if let Some(v) = v {
                out.copy_from_slice(&v);
            } else {
                out.fill(f64::NAN);
            }
        };
        let result = newton_cg::minimize(f, g, h, &x, &options);
        callbacks.finish()?;
        let r = result.map_err(|e| domain(e.to_string()))?;
        let derivative = DimExpr::Mul(vec![e.clone(), DimExpr::Inverse(Box::new(d.clone()))]);
        Ok(record(vec![
            (
                "x",
                output_tensor(
                    vec![r.x.len()],
                    r.x.into_iter()
                        .map(|v| api::scale_product(v, &[xs], &[]))
                        .collect::<Result<_, _>>()?,
                    d,
                )?,
            ),
            (
                "value",
                super::super::dimensional_float(api::scale_product(r.value, &[ys], &[])?, e),
            ),
            (
                "gradient_norm",
                super::super::dimensional_float(
                    api::scale_product(r.gradient_norm, &[ys], &[xs])?,
                    derivative,
                ),
            ),
            (
                "normalized_gradient_norm",
                Value::Float(finite(r.gradient_norm)?),
            ),
            ("status", Value::String(format!("{:?}", r.status))),
            ("iterations", Value::int(r.iterations)),
            ("objective_evaluations", Value::int(r.objective_evaluations)),
            ("gradient_evaluations", Value::int(r.gradient_evaluations)),
            (
                "hessian_vector_evaluations",
                Value::int(r.hessian_vector_evaluations),
            ),
            (
                "negative_curvature_steps",
                Value::int(r.negative_curvature_steps),
            ),
            ("boundary_steps", Value::int(r.boundary_steps)),
            ("accepted_steps", Value::int(r.accepted_steps)),
        ]))
    }
}
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    register(
        i,
        "num_analysis_newton_cg_scaled",
        "crate::kernels::newton_cg::minimize",
        "Genuine matrix-free trust-region Newton-CG: physical point,positive coordinate/objective references,explicit normalized-coordinate metric options,shared callback evaluation limit,then normalized objective/gradient/Hessian-vector callbacks. Exact derivative actions come from the supplied callbacks; no finite-difference substitution. Returns physical x,value,gradient norm plus actual native status and work counters.",
        8,
        Newton,
    )
}
