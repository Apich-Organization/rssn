//! Measured directional jump across a branch cut, with checked callback work.
use crate::kernels::complex_coords::branch_jump_checked;
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
use super::super::real;
use super::super::register;
use super::super::scalar;
use super::super::scale;
use super::super::typed;
use super::super::var;
use super::complex_type;
use super::decode;
use super::dimension;
use super::value;

struct BranchJump;

impl EvalDyn for BranchJump {
    fn signature(&self) -> Signature {
        let d = var("D");
        let e = var("E");
        typed(
            vec![
                complex_type(d.clone()),
                real(d.clone()),
                TypeExpr::Float,
                real(e.clone()),
                TypeExpr::Int,
                function(vec![complex_type(d)], complex_type(e.clone())),
            ],
            complex_type(e),
            &["D", "E"],
        )
    }

    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let signature = self.signature();
        arity(args, signature.params.len())?;
        let budget = count(&args[4], 4)?;
        if budget < 4 {
            return Err(OpError::BudgetExceeded);
        }
        let callbacks = Callbacks::new(ctx, args, &signature, 5, budget)?;
        let input = dimension(&callbacks.signatures[0].params[0])?;
        let output = dimension(&callbacks.signatures[0].ret)?;
        let point = decode(&args[0], &input, 0, ctx)?;
        let step = scale(&args[1], &input, 1, ctx)?;
        let angle = scalar(&args[2], &DimExpr::Dimless, 2, ctx)?;
        let amplitude = scale(&args[3], &output, 3, ctx)?;
        let normal = Complex64::from_polar(step, angle);
        let sample = |z| {
            callbacks
                .call(0, &[value(z, input.clone())], |v| {
                    decode(&v, &output, 0, ctx).map(|v| v / amplitude)
                })
                .unwrap_or(Complex64::new(f64::NAN, f64::NAN))
        };
        let estimated = branch_jump_checked(sample, point, normal);
        callbacks.finish()?;
        crate::cancellation::checkpoint()?;
        let jump = amplitude
            * estimated.ok_or_else(|| {
                domain("branch jump needs four distinct finite rounded sample locations")
            })?;
        if jump.is_finite() {
            Ok(value(jump, output))
        } else {
            Err(domain("branch jump exceeds Float range"))
        }
    }
}

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    register(
        i,
        "num_complex_coords_branch_jump_checked",
        "crate::kernels::complex_coords::branch_cut_jump_numeric",
        "Estimates the positive-normal lip minus the negative-normal lip of a branch cut. Pass point, coordinate step, bare normal angle in radians, callback-value scale, callback budget and callback.",
        6,
        BranchJump,
    )
}
