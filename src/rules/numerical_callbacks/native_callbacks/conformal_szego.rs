//! Szegő preparation and one-shot mapping through bounded boundary callbacks.
use crate::kernels::conformal::SzegoKerzmanSteinSolver;
use crate::kernels::conformal::szego_kerzman_stein_map;
use crate::rules::numerical_abi::codecs::Output;
use crate::rules::numerical_abi::codecs::szego_model_type;
use num_complex::Complex64;

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
use super::super::arity;
use super::super::count;
use super::super::domain;
use super::super::function;
use super::super::record;
use super::super::register;
use super::super::typed;
use super::decode;

#[derive(Clone, Copy)]
struct Method(bool);

impl EvalDyn for Method {
    fn signature(&self) -> Signature {
        let boundary = function(vec![TypeExpr::Float], TypeExpr::Complex);
        let (mut args, result) = if self.0 {
            (
                vec![
                    TypeExpr::Complex,
                    TypeExpr::Complex,
                    TypeExpr::Int,
                    TypeExpr::Int,
                ],
                TypeExpr::Record(vec![
                    ("mapped".into(), TypeExpr::Complex),
                    ("szego".into(), TypeExpr::Complex),
                ]),
            )
        } else {
            (
                vec![TypeExpr::Complex, TypeExpr::Int, TypeExpr::Int],
                szego_model_type(),
            )
        };
        args.extend([boundary.clone(), boundary]);
        typed(args, result, &[])
    }

    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let signature = self.signature();
        arity(args, signature.params.len())?;
        let dimension = DimExpr::Dimless;
        let center = decode(&args[0], &dimension, 0, ctx)?;
        let start = if self.0 { 2 } else { 1 };
        let samples = count(&args[start], start)?;
        let budget = count(&args[start + 1], start + 1)?;
        if !(16..=192).contains(&samples) || budget < 2 * samples {
            return Err(OpError::BudgetExceeded);
        }
        let callbacks = Callbacks::new(ctx, args, &signature, start + 2, budget)?;
        let evaluate = |index, t| {
            callbacks
                .call(index, &[Value::Float(t)], |v| {
                    decode(&v, &dimension, index, ctx)
                })
                .unwrap_or(Complex64::new(f64::NAN, f64::NAN))
        };
        let gamma = |t| evaluate(0, t);
        let derivative = |t| evaluate(1, t);
        let result = if self.0 {
            let point = decode(&args[1], &dimension, 1, ctx)?;
            szego_kerzman_stein_map(gamma, derivative, samples, center, point)
                .map_err(domain)
                .and_then(|(mapped, szego)| {
                    Ok(record(vec![
                        ("mapped", mapped.encode()?),
                        ("szego", szego.encode()?),
                    ]))
                })
        } else {
            SzegoKerzmanSteinSolver::solve(gamma, derivative, samples, center)
                .map_err(domain)
                .and_then(Output::encode)
        };
        callbacks.finish()?;
        crate::cancellation::checkpoint()?;
        result
    }
}

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    register(
        i,
        "num_conformal_szego_prepare_checked",
        "crate::kernels::conformal::SzegoKerzmanSteinSolver::solve",
        "Prepares reusable Szegő quadrature from a counterclockwise boundary callback and its derivative on a normalized chart. Pass center, nodes, callback budget and callbacks.",
        5,
        Method(false),
    )?;
    register(
        i,
        "num_conformal_szego_map_checked",
        "crate::kernels::conformal::szego_kerzman_stein_map",
        "Returns a Riemann map and the Szegő kernel for one interior query. Pass center, query, nodes, callback budget, boundary and derivative callbacks.",
        6,
        Method(true),
    )
}
