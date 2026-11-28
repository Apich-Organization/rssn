//! Reusable Sylvester and Lyapunov solvers with explicit matrix references.
use super::Array;
use super::ArrayD;
use super::Handle;
use super::Installer;
use super::IxDyn;
use super::NativeType;
use super::NewObject;
use super::OpError;
use super::RuleError;
use super::api;
use super::array;
use super::dim;
use super::err;
use super::lowrank::matrix;
use super::object_binding;
use super::product;
use super::sparse::inverse;
use crate::kernels::dense::Mat;
use crate::kernels::schur::SylvesterPlan;

struct PreparedSylvester {
    plan: SylvesterPlan,
    coefficient_reference: f64,
}

impl NativeType for PreparedSylvester {
    const KIND: &'static str = "rssn/native/sylvester-plan/v1";
    // Three square matrices per side, each at most 128 by 128.
    const MAX_BYTES: usize = 6 * 128 * 128 * size_of::<f64>() + 8192;

    fn validate(&self) -> Result<(), OpError> {
        if self.coefficient_reference.is_finite() && self.coefficient_reference > 0.0 {
            Ok(())
        } else {
            Err(err("invalid Sylvester coefficient reference"))
        }
    }

    fn retained_bytes(&self) -> usize {
        size_of::<Self>() + self.plan.heap_bytes()
    }
}

fn coefficient(
    a: &Array,
    reference: f64,
) -> Result<Mat, OpError> {
    let mut result = matrix(a)?;
    if a.scale.to_bits() != reference.to_bits() {
        for value in &mut result.data {
            *value = api::scale_without_underflow(*value, &[a.scale], &[reference])?;
        }
    }
    Ok(result)
}

object_binding!(
    sylvester_plan_new,
    "num_linear_sylvester_plan_new",
    "crate::kernels::schur::SylvesterPlan::new",
    "Prepares the Sylvester equation A X + X B = C from two square native arrays with one coefficient dimension. The larger of their positive reference scales is shared by both coefficients. A nonzero coefficient that underflows in that reference returns an error. Reuse the resulting plan for different right-hand sides.",
    true,
    (a: Handle<Array>, b: Handle<Array>) -> NewObject<PreparedSylvester>,
    |ctx, a, b| {
        let dimension = dim(ctx, &a.value)?;
        api::same_dimension(&dimension, &dim(ctx, &b.value)?)?;
        a.with(ctx, |a| b.with(ctx, |b| {
            let coefficient_reference = a.scale.max(b.scale);
            let plan = SylvesterPlan::new(
                &coefficient(a, coefficient_reference)?,
                &coefficient(b, coefficient_reference)?,
            ).map_err(err)?;
            Ok(NewObject::new(PreparedSylvester { plan, coefficient_reference }, vec![dimension]))
        }))
    }
);

object_binding!(
    lyapunov_plan_new,
    "num_linear_lyapunov_plan_new",
    "crate::kernels::schur::SylvesterPlan::for_lyapunov",
    "Prepares the continuous-time Lyapunov equation A X + X A^T = Q from one square native array. Reuse one real Schur decomposition for different right-hand sides.",
    true,
    (a: Handle<Array>) -> NewObject<PreparedSylvester>,
    |ctx, a| {
        let dimension = dim(ctx, &a.value)?;
        a.with(ctx, |a| {
            let plan = SylvesterPlan::for_lyapunov(&matrix(a)?).map_err(err)?;
            Ok(NewObject::new(PreparedSylvester {
                plan, coefficient_reference: a.scale,
            }, vec![dimension]))
        })
    }
);

object_binding!(
    sylvester_plan_solve,
    "num_linear_sylvester_plan_solve",
    "crate::kernels::schur::SylvesterPlan::solve",
    "Solves a prepared Sylvester or Lyapunov equation for a native right-hand side with the expected row and column counts. The solution uses the dimension and positive reference of the right-hand side divided by those of the coefficients.",
    true,
    (plan: Handle<PreparedSylvester>, rhs: Handle<Array>) -> NewObject<Array>,
    |ctx, plan, rhs| {
        let dimension = product(dim(ctx, &rhs.value)?, inverse(dim(ctx, &plan.value)?)?)?;
        plan.with(ctx, |plan| rhs.with(ctx, |rhs| {
            let solution = plan.plan.solve(&matrix(rhs)?).map_err(err)?;
            let data = ArrayD::from_shape_vec(IxDyn(&[solution.rows, solution.cols]), solution.data).map_err(err)?;
            let reference = api::scale_without_underflow(rhs.scale, &[], &[plan.coefficient_reference])?;
            Ok(array(data, reference, dimension))
        }))
    }
);

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &sylvester_plan_new::BINDING,
            &lyapunov_plan_new::BINDING,
            &sylvester_plan_solve::BINDING,
        ],
        "numerical_matrix",
    )
}
