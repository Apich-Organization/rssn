//! Special-function adapters with checked domains and reusable parameters.
use super::Installer;
use super::RuleError;
use super::native_object_api as api;
use crate::kernel::CallCtx;
use crate::kernel::OpError;
use crate::kernel::Value;
use crate::kernel::signature::TypeExpr;
use crate::kernels::fresnel::fresnel_sc_interval;
use crate::kernels::number_theory::PrimeSqrtPlan;
use crate::sim::elliptic::EllipticPiPlan;
use api::Handle;
use api::Many;
use api::Maybe;
use api::NativeArg;
use api::NativeType;
use api::NewObject;
use api::object_binding;

struct FresnelEndpoint(f64);

impl NativeArg for FresnelEndpoint {
    fn ty() -> TypeExpr {
        TypeExpr::Float
    }

    fn decode(
        value: &Value,
        _: usize,
        _: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        api::normalized(value)?;
        api::extended_numeric(value).map(Self)
    }
}

object_binding!(
    fresnel_interval,
    "num_fresnel_sc_interval",
    "crate::kernels::fresnel::fresnel_sc_interval",
    "Returns the dimensionless Fresnel sine and cosine integrals from a to b as item0 and item1. Narrow intervals and far-field tails use direct differences. Infinite bounds give improper integrals. Reversing the bounds negates both results.",
    true,
    (a: FresnelEndpoint, b: FresnelEndpoint) -> (f64, f64),
    |_, a, b| Ok(fresnel_sc_interval(a.0, b.0))
);

impl NativeType for EllipticPiPlan {
    const KIND: &'static str = "rssn/native/elliptic-pi-plan/v1";
    const MAX_BYTES: usize = 256;

    fn validate(&self) -> Result<(), OpError> {
        if self.complete().is_finite() && self.complete() > 0.0 {
            Ok(())
        } else {
            Err(api::error(
                "complete elliptic Pi must be finite and positive",
            ))
        }
    }
}

object_binding!(
    elliptic_pi_plan_new,
    "num_special_elliptic_pi_plan_new",
    "crate::sim::elliptic::EllipticPiPlan::new",
    "Prepares the third-kind elliptic integral at fixed characteristic n and parameter m. Use finite dimensionless values below one.",
    true,
    (n: f64, m: f64) -> NewObject<EllipticPiPlan>,
    |_, n, m| EllipticPiPlan::new(n, m).map(NewObject::normalized).map_err(api::error)
);

object_binding!(
    elliptic_pi_plan_complete,
    "num_special_elliptic_pi_plan_complete",
    "crate::sim::elliptic::EllipticPiPlan::complete",
    "Returns the dimensionless complete third-kind elliptic integral of the prepared parameters.",
    true,
    (plan: Handle<EllipticPiPlan>) -> f64,
    |ctx, plan| plan.with(ctx, |plan| Ok(plan.complete()))
);

object_binding!(
    elliptic_pi_plan_eval,
    "num_special_elliptic_pi_plan_eval",
    "crate::sim::elliptic::EllipticPiPlan::evaluate",
    "Evaluates the incomplete elliptic integral of the third kind at an amplitude in radians, tracking complete periods of the prepared parameters.",
    true,
    (plan: Handle<EllipticPiPlan>, phi: f64) -> f64,
    |ctx, plan, phi| plan.with(ctx, |plan| Ok(plan.evaluate(phi)))
);

object_binding!(
    elliptic_pi_plan_eval_many,
    "num_special_elliptic_pi_plan_eval_many",
    "crate::sim::elliptic::EllipticPiPlan::evaluate",
    "Evaluates the incomplete elliptic integral of the third kind in amplitude order. All amplitudes are in radians and share the same prepared characteristic, parameter and complete integral.",
    true,
    (plan: Handle<EllipticPiPlan>, amplitudes: Vec<f64>) -> Vec<f64>,
    |ctx, plan, amplitudes| {
        if amplitudes.len() > 65_536 {
            return Err(api::error("elliptic Pi amplitude batch exceeds the work budget"));
        }
        plan.with(ctx, |plan| Ok(crate::kernels::par::map_slice(
            &amplitudes, 256, |&phi| plan.evaluate(phi),
        )))
    }
);

object_binding!(
    bell_batch,
    "num_combinatorics_bell_batch",
    "crate::kernels::combinatorics::bell_batch",
    "Returns a Bell number for each unsigned integer index in input order. Finite values are Some, and counts that exceed Float are None. The batch reuses the same finite-count table as the scalar Bell function.",
    true,
    (indices: Vec<u64>) -> Many<Maybe<f64>>,
    |_, indices| {
        let counts = crate::kernels::combinatorics::bell_batch(&indices);
        Ok(Many(counts.into_iter().map(|value| Maybe(value.is_finite().then_some(value))).collect()))
    }
);

impl NativeType for PrimeSqrtPlan {
    const KIND: &'static str = "rssn/native/number-theory/prime-sqrt-plan/v1";
    const MAX_BYTES: usize = size_of::<Self>();

    fn validate(&self) -> Result<(), OpError> {
        Ok(())
    }
}

object_binding!(
    prime_sqrt_plan_new,
    "num_number_theory_prime_sqrt_plan_new",
    "crate::kernels::number_theory::PrimeSqrtPlan::new",
    "Checks a 64-bit modulus and prepares its square-root method. Returns None for a modulus other than a prime. Modulus two is supported.",
    true,
    (modulus: u64) -> Maybe<NewObject<PrimeSqrtPlan>>,
    |_, modulus| Ok(Maybe(PrimeSqrtPlan::new(modulus).map(NewObject::normalized)))
);

object_binding!(
    prime_sqrt_plan_modulus,
    "num_number_theory_prime_sqrt_plan_modulus",
    "crate::kernels::number_theory::PrimeSqrtPlan::modulus",
    "Returns the prime modulus retained by the prepared square-root plan.",
    true,
    (plan: Handle<PrimeSqrtPlan>) -> u64,
    |ctx, plan| plan.with(ctx, |plan| Ok(plan.modulus()))
);

object_binding!(
    prime_sqrt_plan_sqrt,
    "num_number_theory_prime_sqrt_plan_sqrt",
    "crate::kernels::number_theory::PrimeSqrtPlan::sqrt",
    "Reduces the unsigned input modulo the prepared prime. Returns the smaller square root, or None if the input is a quadratic nonresidue.",
    true,
    (plan: Handle<PrimeSqrtPlan>, value: u64) -> Option<u64>,
    |ctx, plan, value| plan.with(ctx, |plan| Ok(plan.sqrt(value)))
);

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &fresnel_interval::BINDING,
            &elliptic_pi_plan_new::BINDING,
            &elliptic_pi_plan_complete::BINDING,
            &elliptic_pi_plan_eval::BINDING,
            &elliptic_pi_plan_eval_many::BINDING,
            &bell_batch::BINDING,
            &prime_sqrt_plan_new::BINDING,
            &prime_sqrt_plan_modulus::BINDING,
            &prime_sqrt_plan_sqrt::BINDING,
        ],
        "numerical_special",
    )
}
