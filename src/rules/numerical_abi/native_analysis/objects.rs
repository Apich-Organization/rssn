//! Retained native quadrature generators and physically measured Gaussian models.
use crate::graph::RuleError;
use crate::graph::rule::Installer;
use crate::kernel::CallCtx;
use crate::kernel::OpError;
use crate::kernel::signature::DimExpr;
use crate::kernels::quadrature::Gaussian;
use crate::kernels::quadrature::Sobol;
use crate::rules::numerical_abi::native_object_api::Dynamic;
use crate::rules::numerical_abi::native_object_api::Handle;
use crate::rules::numerical_abi::native_object_api::Input;
use crate::rules::numerical_abi::native_object_api::NativeType;
use crate::rules::numerical_abi::native_object_api::NewObject;
use crate::rules::numerical_abi::native_object_api::object_binding;
use crate::rules::numerical_abi::native_object_api::{
    self as object,
};

impl NativeType for Gaussian {
    const KIND: &'static str = "rssn/native/analysis/gaussian/v1";
    const MAX_BYTES: usize = 64;

    fn validate(&self) -> Result<(), OpError> {
        if self.mean.is_finite() && self.sd.is_finite() && self.sd > 0.0 {
            Ok(())
        } else {
            Err(object::error(
                "Gaussian mean must be finite and standard deviation positive",
            ))
        }
    }
}
struct SobolObject {
    native: Sobol,
    dimension: usize,
}
impl NativeType for SobolObject {
    const KIND: &'static str = "rssn/native/analysis/sobol/v1";
    const MAX_BYTES: usize = 8192;

    fn validate(&self) -> Result<(), OpError> {
        if (1..=8).contains(&self.dimension) {
            Ok(())
        } else {
            Err(object::error("Sobol supports1..8 normalized dimensions"))
        }
    }
}
object_binding!(gaussian_new,"num_analysis_gaussian_new","crate::kernels::quadrature::Gaussian::new","Gaussian model with same-unit physical mean and positive standard deviation.",false,(mean:Input,sd:Input)->NewObject<Gaussian>,|_ctx:&CallCtx<'_>,mean:Input,sd:Input|{
 let d=object::dimension(&mean.0)?;object::same_dimension(&d,&object::dimension(&sd.0)?)?;
 let mean=super::super::numeric(&mean.0,0)?;let sd=super::super::numeric(&sd.0,1)?;
 Ok(NewObject::new(Gaussian::new(mean,sd),vec![d]))
});
object_binding!(gaussian_pdf,"num_analysis_gaussian_pdf","crate::kernels::quadrature::Gaussian::pdf","Physical Gaussian density has the reciprocal coordinate dimension.",false,(model:Handle<Gaussian>,x:Input)->Dynamic,|ctx:&CallCtx<'_>,model:Handle<Gaussian>,x:Input|{
 let d=model.dimensions(ctx)?.into_iter().next().ok_or_else(||object::error("Gaussian model lacks unit domain"))?;object::same_dimension(&d,&object::dimension(&x.0)?)?;
 let x=super::super::numeric(&x.0,1)?;let density=model.with(ctx,|m|Ok(m.pdf(x)))?;
 Ok(Dynamic(object::quantity(density,DimExpr::Inverse(Box::new(d)))?))
});
object_binding!(gaussian_weight,"num_analysis_gaussian_weight","crate::kernels::quadrature::Gaussian::weight","Dimensionless importance weight between measured Gaussian models in the same coordinate domain.",false,(target:Handle<Gaussian>,proposal:Handle<Gaussian>,x:Input)->f64,|ctx:&CallCtx<'_>,target:Handle<Gaussian>,proposal:Handle<Gaussian>,x:Input|{
 let d=target.dimensions(ctx)?.into_iter().next().ok_or_else(||object::error("Gaussian model lacks unit domain"))?;
 let pd=proposal.dimensions(ctx)?.into_iter().next().ok_or_else(||object::error("Gaussian proposal lacks unit domain"))?;
 object::same_dimension(&d,&pd)?;object::same_dimension(&d,&object::dimension(&x.0)?)?;let x=super::super::numeric(&x.0,2)?;
 target.with(ctx,|t|proposal.with(ctx,|p|Ok(t.weight(*p,x))))
});
object_binding!(sobol_new,"num_analysis_sobol_new","crate::kernels::quadrature::Sobol::new","Create a real native Sobol sequence in1..8 normalized dimensions; rejects unsupported dimensions instead of clamping.",true,(dimension:usize)->NewObject<SobolObject>,|_ctx:&CallCtx<'_>,dimension:usize|{
 if !(1..=8).contains(&dimension){return Err(object::error("Sobol dimension must be1..8"));}Ok(NewObject::normalized(SobolObject{native:Sobol::new(dimension),dimension}))
});
object_binding!(sobol_next,"num_analysis_sobol_next","crate::kernels::quadrature::Sobol::next_point","Advance the retained native Sobol sequence, returning normalized unit-cube coordinates.",true,(sequence:Handle<SobolObject>)->Vec<f64>,|ctx:&CallCtx<'_>,sequence:Handle<SobolObject>|sequence.with_mut(ctx,|s|Ok(s.native.next_point())));

object_binding!(gaussian_sample,"num_analysis_gaussian_sample","crate::kernels::quadrature::Gaussian::sample","Sample a measured Gaussian using the retained shared native RNG handle.",false,(model:Handle<Gaussian>,rng:Handle<crate::kernels::random::Rng>)->Dynamic,|ctx:&CallCtx<'_>,model:Handle<Gaussian>,rng:Handle<crate::kernels::random::Rng>|{
 let d=model.dimensions(ctx)?.into_iter().next().ok_or_else(||object::error("Gaussian model lacks unit domain"))?;
 let model=model.with(ctx,|m|Ok(*m))?;let value=rng.with_mut(ctx,|rng|Ok(model.sample(rng)))?;
 Ok(Dynamic(object::quantity(value,d)?))
});

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    object::install(
        i,
        &[
            &gaussian_new::BINDING,
            &gaussian_pdf::BINDING,
            &gaussian_weight::BINDING,
            &gaussian_sample::BINDING,
            &sobol_new::BINDING,
            &sobol_next::BINDING,
        ],
        "numerical_vector",
    )
}
