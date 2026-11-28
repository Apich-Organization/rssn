//! Krawczyk and Rouché certification with callbacks enclosing whole disks.

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
use super::arity;
use super::complex::decode;
use super::count;
use super::domain;
use super::function;
use super::positive;
use super::register;
use super::typed;
use crate::kernels::complex_ball::ComplexBall;
use crate::kernels::complex_ball::certify_root_1d;
use crate::kernels::complex_ball::certify_root_nd;
use crate::kernels::complex_ball::certify_rouche_circle;
use crate::rules::numerical_abi::codecs::Arg;
use crate::rules::numerical_abi::codecs::CheckedComplexDisk as Disk;
use crate::rules::numerical_abi::codecs::CheckedComplexDiskJet as Jet;
use crate::rules::numerical_abi::codecs::CheckedComplexDiskMatrix as Matrix;
use crate::rules::numerical_abi::codecs::CheckedComplexDiskVector as Vector;
use crate::rules::numerical_abi::codecs::Output;

#[derive(Clone, Copy)]
enum Method {
    Scalar,
    System,
    Rouche,
}

const SYSTEM_LIMIT: usize = 64;

impl EvalDyn for Method {
    fn signature(&self) -> Signature {
        match self {
            | Self::Scalar => {
                typed(
                    vec![
                        <Disk as Arg>::ty(),
                        TypeExpr::Int,
                        function(vec![<Jet as Arg>::ty()], <Option<Jet> as Output>::ty()),
                    ],
                    <Option<Disk> as Output>::ty(),
                    &[],
                )
            },
            | Self::System => {
                typed(
                    vec![
                        <Vector as Arg>::ty(),
                        TypeExpr::Int,
                        function(
                            vec![<Vector as Arg>::ty()],
                            <Option<Vector> as Output>::ty(),
                        ),
                        function(
                            vec![<Vector as Arg>::ty()],
                            <Option<Matrix> as Output>::ty(),
                        ),
                    ],
                    <Option<Vector> as Output>::ty(),
                    &[],
                )
            },
            | Self::Rouche => {
                let callback = function(vec![<Disk as Arg>::ty()], <Option<Disk> as Output>::ty());
                typed(
                    vec![
                        TypeExpr::Complex,
                        TypeExpr::Float,
                        TypeExpr::Int,
                        TypeExpr::Int,
                        callback.clone(),
                        callback,
                    ],
                    TypeExpr::Bool,
                    &[],
                )
            },
        }
    }

    #[allow(
        clippy::indexing_slicing,
        clippy::arithmetic_side_effects,
        reason = "the signature checks argument arity and the bounded circle count fits twice in usize"
    )]
    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let signature = self.signature();
        arity(args, signature.params.len())?;
        match self {
            | Self::Scalar => {
                let ball = Disk::decode(&args[0], 0)?.0;
                let budget = count(&args[1], 1)?;
                if budget < 2 {
                    return Err(OpError::BudgetExceeded);
                }
                let callbacks = Callbacks::new(ctx, args, &signature, 2, budget)?;
                let certified = certify_root_1d(
                    |jet| {
                        let input = Jet(jet).encode().ok()?;
                        callbacks
                            .call(0, &[input], |value| Option::<Jet>::decode(&value, 2))
                            .flatten()
                            .map(|output| output.0)
                    },
                    ball,
                );
                callbacks.finish()?;
                crate::cancellation::checkpoint()?;
                certified.map(Disk).encode()
            },
            | Self::System => {
                let Value::Tensor(tensor) = &args[0] else {
                    return Err(domain("expected a rank-two complex disk box"));
                };
                let [n, 3] = tensor.shape.as_slice() else {
                    return Err(domain("complex disk boxes use shape [n, 3]"));
                };
                let budget = count(&args[1], 1)?;
                if *n > SYSTEM_LIMIT || budget < 3 {
                    return Err(OpError::BudgetExceeded);
                }
                let n = *n;
                let balls = Vector::decode(&args[0], 0)?.0;
                let callbacks = Callbacks::new(ctx, args, &signature, 2, budget)?;
                let points = |points: &[ComplexBall]| Vector(points.to_vec()).encode().ok();
                let certified = certify_root_nd(
                    |z| {
                        callbacks
                            .call(0, &[points(z)?], |value| {
                                let answer = Option::<Vector>::decode(&value, 2)?;
                                if answer.as_ref().is_some_and(|answer| answer.0.len() != n) {
                                    return Err(domain(
                                        "value callback returned a vector of the wrong length",
                                    ));
                                }
                                Ok(answer.map(|answer| answer.0))
                            })
                            .flatten()
                    },
                    |z| {
                        callbacks
                            .call(1, &[points(z)?], |value| {
                                let answer = Option::<Matrix>::decode(&value, 3)?;
                                if answer.as_ref().is_some_and(|answer| {
                                    answer.0.len() != n || answer.0.iter().any(|row| row.len() != n)
                                }) {
                                    return Err(domain(
                                        "Jacobian callback returned a matrix of the wrong shape",
                                    ));
                                }
                                Ok(answer.map(|answer| answer.0))
                            })
                            .flatten()
                    },
                    &balls,
                );
                callbacks.finish()?;
                crate::cancellation::checkpoint()?;
                certified.map(Vector).encode()
            },
            | Self::Rouche => {
                let center: Complex64 = decode(&args[0], &DimExpr::Dimless, 0, ctx)?;
                let radius = positive(&args[1], 1, ctx)?;
                let segments = count(&args[2], 2)?;
                let budget = count(&args[3], 3)?;
                if !(4..=4096).contains(&segments) {
                    return Err(domain(
                        "Rouché certification requires between 4 and 4096 arcs",
                    ));
                }
                if budget < 2 * segments {
                    return Err(OpError::BudgetExceeded);
                }
                let callbacks = Callbacks::new(ctx, args, &signature, 4, budget)?;
                let sample = |n, disk| {
                    callbacks
                        .call(n, &[Disk(disk).encode().ok()?], |value| {
                            Option::<Disk>::decode(&value, 4 + n)
                        })
                        .flatten()
                        .map(|output| output.0)
                };
                let certified = certify_rouche_circle(
                    |z| sample(0, z),
                    |z| sample(1, z),
                    center,
                    radius,
                    segments,
                );
                callbacks.finish()?;
                crate::cancellation::checkpoint()?;
                Ok(Value::Bool(certified))
            },
        }
    }
}

#[allow(
    clippy::single_call_fn,
    reason = "registration entrypoint for disk certification"
)]
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    register(
        i,
        "num_complex_disk_certify_root",
        "crate::kernels::complex_ball::certify_root_1d",
        "Krawczyk existence and uniqueness test. The callback must enclose the function and derivative on each entire normalized disk. None means the test could not certify a root.",
        3,
        Method::Scalar,
    )?;
    register(
        i,
        "num_complex_disk_certify_system",
        "crate::kernels::complex_ball::certify_root_nd",
        "Krawczyk existence and uniqueness test for a normalized complex box. The Jacobian callback must enclose every derivative on the entire box. None means no certificate was obtained.",
        4,
        Method::System,
    )?;
    register(
        i,
        "num_complex_disk_certify_rouche",
        "crate::kernels::complex_ball::certify_rouche_circle",
        "Checks the strict Rouché inequality on every covering disk of the normalized circle. Both callbacks must enclose their functions on the entire input disk and the functions must be holomorphic in the enclosed region.",
        6,
        Method::Rouche,
    )
}
