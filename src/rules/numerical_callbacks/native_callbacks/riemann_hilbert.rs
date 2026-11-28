//! Checked matrix Riemann–Hilbert collocation with a dimensionless jump.

use super::CallCtx;
use super::Callbacks;
use super::EvalDyn;
use super::Installer;
use super::OpError;
use super::RuleError;
use super::Signature;
use super::TypeExpr;
use super::Value;
use super::arity;
use super::count;
use super::domain;
use super::function;
use super::record;
use super::register;
use super::typed;
use crate::kernels::riemann_hilbert::Matrix2;
use crate::kernels::riemann_hilbert::MatrixRhSolution;
use crate::kernels::riemann_hilbert::solve_matrix_rh_chebyshev;
use crate::rules::numerical_abi::codecs::Arg;
use crate::rules::numerical_abi::codecs::Output;
use num_complex::Complex64;

#[path = "riemann_hilbert_scalar.rs"]
mod scalar;

struct SolveMatrix;

impl EvalDyn for SolveMatrix {
    fn signature(&self) -> Signature {
        typed(
            vec![
                TypeExpr::Int,
                TypeExpr::Int,
                function(vec![TypeExpr::Float], <Matrix2 as Arg>::ty()),
            ],
            TypeExpr::Record(vec![
                ("solution".into(), <MatrixRhSolution as Output>::ty()),
                ("evaluations".into(), TypeExpr::Int),
            ]),
            &[],
        )
    }

    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let signature = self.signature();
        arity(args, signature.params.len())?;
        let n = count(&args[0], 0)?;
        let budget = count(&args[1], 1)?;
        if n < 4 {
            return Err(domain(
                "matrix Riemann–Hilbert collocation requires at least four nodes",
            ));
        }
        if n > 64 || budget < n + 2 {
            return Err(OpError::BudgetExceeded);
        }
        let callbacks = Callbacks::new(ctx, args, &signature, 2, budget)?;
        let sample = |point| callbacks.call(0, &[Value::Float(point)], |v| Matrix2::decode(&v, 2));
        let endpoint_identity = [-1.0, 1.0].into_iter().all(|point| {
            sample(point).is_some_and(|jump| jump.sub(Matrix2::identity()).norm() <= 1e-10)
        });
        if !endpoint_identity {
            callbacks.finish()?;
            return Err(domain(
                "matrix Riemann–Hilbert jump must approach the identity at both endpoints",
            ));
        }
        let invalid = Matrix2 {
            m: [Complex64::new(f64::NAN, f64::NAN); 4],
        };
        let solved = solve_matrix_rh_chebyshev(|point| sample(point).unwrap_or(invalid), n);
        let evaluations = callbacks.used.get();
        callbacks.finish()?;
        crate::cancellation::checkpoint()?;
        Ok(record(vec![
            ("solution", solved.map_err(domain)?.encode()?),
            ("evaluations", Value::int(evaluations)),
        ]))
    }
}

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    scalar::install(i)?;
    register(
        i,
        "num_riemann_hilbert_solve",
        "crate::kernels::riemann_hilbert::solve_matrix_rh_chebyshev",
        "Collocates a dimensionless 2x2 jump on the normalized interval (-1, 1). The callback must give the identity at both endpoints. Uses 4 to 64 Gauss–Legendre nodes and a host budget for both endpoint checks and each collocation sample. Returns the reusable collocation solution and actual callback count.",
        3,
        SolveMatrix,
    )
}
