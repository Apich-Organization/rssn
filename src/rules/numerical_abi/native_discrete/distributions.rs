#![allow(clippy::wildcard_imports)] // Shared checked object method declarations.
use super::*;
use crate::kernels::stats::BinomialDist;
use crate::kernels::stats::ExponentialDist;
use crate::kernels::stats::GammaDist;
use crate::kernels::stats::NormalDist;
use crate::kernels::stats::PoissonDist;
use crate::kernels::stats::UniformDist;
impl NativeType for NormalDist {
    const KIND: &'static str = "rssn/native/discrete/legacy-normaldist/v1";
    const MAX_BYTES: usize = 512;

    fn validate(&self) -> Result<(), OpError> {
        Ok(())
    }
}
impl NativeType for UniformDist {
    const KIND: &'static str = "rssn/native/discrete/legacy-uniformdist/v1";
    const MAX_BYTES: usize = 512;

    fn validate(&self) -> Result<(), OpError> {
        Ok(())
    }
}
impl NativeType for ExponentialDist {
    const KIND: &'static str = "rssn/native/discrete/legacy-exponentialdist/v1";
    const MAX_BYTES: usize = 512;

    fn validate(&self) -> Result<(), OpError> {
        Ok(())
    }
}
impl NativeType for GammaDist {
    const KIND: &'static str = "rssn/native/discrete/legacy-gammadist/v1";
    const MAX_BYTES: usize = 512;

    fn validate(&self) -> Result<(), OpError> {
        Ok(())
    }
}
impl NativeType for PoissonDist {
    const KIND: &'static str = "rssn/native/discrete/legacy-poissondist/v1";
    const MAX_BYTES: usize = 512;

    fn validate(&self) -> Result<(), OpError> {
        Ok(())
    }
}
impl NativeType for BinomialDist {
    const KIND: &'static str = "rssn/native/discrete/legacy-binomialdist/v1";
    const MAX_BYTES: usize = 512;

    fn validate(&self) -> Result<(), OpError> {
        Ok(())
    }
}
object_binding!(normal_dist_new,"num_native_normal_dist_new","crate::kernels::stats::NormalDist::new","Checked native legacy distribution. Continuous parameters and queries retain coherent physical units; inverse CDF restores sample units and PDF has inverse sample units.",false,(mean:Input,sd:Input)->NewObject<NormalDist>,|_ctx:&CallCtx<'_>,mean:Input,sd:Input|{let(m,d)=scalar(&mean)?;let(s,e)=scalar(&sd)?;api::same_dimension(&d,&e)?;Ok(NewObject::new(NormalDist::new(m,s).map_err(error)?,vec![d]))});
object_binding!(normal_dist_cdf,"num_native_normal_dist_cdf","crate::kernels::stats::NormalDist::cdf","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<NormalDist>,x:Input)->f64,|ctx:&CallCtx<'_>,a:Handle<NormalDist>,x:Input|{let d=object_dimension(&a,ctx)?;api::same_dimension(&d,&api::dimension(&x.0)?)?;let x=api::extended_numeric(&x.0)?;a.with(ctx,|o|Ok(o.cdf(x)))});
object_binding!(normal_dist_pdf,"num_native_normal_dist_pdf","crate::kernels::stats::NormalDist::pdf","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<NormalDist>,x:Input)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<NormalDist>,x:Input|{let d=object_dimension(&a,ctx)?;api::same_dimension(&d,&api::dimension(&x.0)?)?;let x=api::extended_numeric(&x.0)?;a.with(ctx,|o|extended_quantity(o.pdf(x),canonical(&(DimExpr::Inverse(Box::new(d))))?))});
object_binding!(normal_dist_inverse_cdf,"num_native_normal_dist_inverse_cdf","crate::kernels::stats::NormalDist::inverse_cdf","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<NormalDist>,probability:f64)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<NormalDist>,probability:f64|{if !(0.0..=1.0).contains(&probability){return Err(error("CDF probability must lie in [0,1]"));}let d=object_dimension(&a,ctx)?;a.with(ctx,|o|extended_quantity(o.inverse_cdf(probability),d))});
object_binding!(uniform_dist_new,"num_native_uniform_dist_new","crate::kernels::stats::UniformDist::new","Checked native legacy distribution. Continuous parameters and queries retain coherent physical units; inverse CDF restores sample units and PDF has inverse sample units.",false,(lo:Input,hi:Input)->NewObject<UniformDist>,|_ctx:&CallCtx<'_>,lo:Input,hi:Input|{let(l,d)=scalar(&lo)?;let(h,e)=scalar(&hi)?;api::same_dimension(&d,&e)?;Ok(NewObject::new(UniformDist::new(l,h).map_err(error)?,vec![d]))});
object_binding!(uniform_dist_cdf,"num_native_uniform_dist_cdf","crate::kernels::stats::UniformDist::cdf","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<UniformDist>,x:Input)->f64,|ctx:&CallCtx<'_>,a:Handle<UniformDist>,x:Input|{let d=object_dimension(&a,ctx)?;api::same_dimension(&d,&api::dimension(&x.0)?)?;let x=api::extended_numeric(&x.0)?;a.with(ctx,|o|Ok(o.cdf(x)))});
object_binding!(uniform_dist_pdf,"num_native_uniform_dist_pdf","crate::kernels::stats::UniformDist::pdf","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<UniformDist>,x:Input)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<UniformDist>,x:Input|{let d=object_dimension(&a,ctx)?;api::same_dimension(&d,&api::dimension(&x.0)?)?;let x=api::extended_numeric(&x.0)?;a.with(ctx,|o|extended_quantity(o.pdf(x),canonical(&(DimExpr::Inverse(Box::new(d))))?))});
object_binding!(uniform_dist_inverse_cdf,"num_native_uniform_dist_inverse_cdf","crate::kernels::stats::UniformDist::inverse_cdf","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<UniformDist>,probability:f64)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<UniformDist>,probability:f64|{if !(0.0..=1.0).contains(&probability){return Err(error("CDF probability must lie in [0,1]"));}let d=object_dimension(&a,ctx)?;a.with(ctx,|o|extended_quantity(o.inverse_cdf(probability),d))});
object_binding!(exponential_dist_new,"num_native_exponential_dist_new","crate::kernels::stats::ExponentialDist::new","Checked native legacy distribution. Continuous parameters and queries retain coherent physical units; inverse CDF restores sample units and PDF has inverse sample units.",false,(rate:Input)->NewObject<ExponentialDist>,|_ctx:&CallCtx<'_>,rate:Input|{let(r,d)=scalar(&rate)?;Ok(NewObject::new(ExponentialDist::new(r).map_err(error)?,vec![canonical(&(DimExpr::Inverse(Box::new(d))))?]))});
object_binding!(exponential_dist_cdf,"num_native_exponential_dist_cdf","crate::kernels::stats::ExponentialDist::cdf","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<ExponentialDist>,x:Input)->f64,|ctx:&CallCtx<'_>,a:Handle<ExponentialDist>,x:Input|{let d=object_dimension(&a,ctx)?;api::same_dimension(&d,&api::dimension(&x.0)?)?;let x=api::extended_numeric(&x.0)?;a.with(ctx,|o|Ok(o.cdf(x)))});
object_binding!(exponential_dist_pdf,"num_native_exponential_dist_pdf","crate::kernels::stats::ExponentialDist::pdf","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<ExponentialDist>,x:Input)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<ExponentialDist>,x:Input|{let d=object_dimension(&a,ctx)?;api::same_dimension(&d,&api::dimension(&x.0)?)?;let x=api::extended_numeric(&x.0)?;a.with(ctx,|o|extended_quantity(o.pdf(x),canonical(&(DimExpr::Inverse(Box::new(d))))?))});
object_binding!(exponential_dist_inverse_cdf,"num_native_exponential_dist_inverse_cdf","crate::kernels::stats::ExponentialDist::inverse_cdf","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<ExponentialDist>,probability:f64)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<ExponentialDist>,probability:f64|{if !(0.0..=1.0).contains(&probability){return Err(error("CDF probability must lie in [0,1]"));}let d=object_dimension(&a,ctx)?;a.with(ctx,|o|extended_quantity(o.inverse_cdf(probability),d))});
object_binding!(gamma_dist_new,"num_native_gamma_dist_new","crate::kernels::stats::GammaDist::new","Checked native legacy distribution. Continuous parameters and queries retain coherent physical units; inverse CDF restores sample units and PDF has inverse sample units.",false,(shape:f64,rate:Input)->NewObject<GammaDist>,|_ctx:&CallCtx<'_>,shape:f64,rate:Input|{let(r,d)=scalar(&rate)?;Ok(NewObject::new(GammaDist::new(shape,r).map_err(error)?,vec![canonical(&(DimExpr::Inverse(Box::new(d))))?]))});
object_binding!(gamma_dist_cdf,"num_native_gamma_dist_cdf","crate::kernels::stats::GammaDist::cdf","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<GammaDist>,x:Input)->f64,|ctx:&CallCtx<'_>,a:Handle<GammaDist>,x:Input|{let d=object_dimension(&a,ctx)?;api::same_dimension(&d,&api::dimension(&x.0)?)?;let x=api::extended_numeric(&x.0)?;a.with(ctx,|o|Ok(o.cdf(x)))});
object_binding!(gamma_dist_pdf,"num_native_gamma_dist_pdf","crate::kernels::stats::GammaDist::pdf","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<GammaDist>,x:Input)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<GammaDist>,x:Input|{let d=object_dimension(&a,ctx)?;api::same_dimension(&d,&api::dimension(&x.0)?)?;let x=api::extended_numeric(&x.0)?;a.with(ctx,|o|extended_quantity(o.pdf(x),canonical(&(DimExpr::Inverse(Box::new(d))))?))});
object_binding!(gamma_dist_inverse_cdf,"num_native_gamma_dist_inverse_cdf","crate::kernels::stats::GammaDist::inverse_cdf","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<GammaDist>,probability:f64)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<GammaDist>,probability:f64|{if !(0.0..=1.0).contains(&probability){return Err(error("CDF probability must lie in [0,1]"));}let d=object_dimension(&a,ctx)?;a.with(ctx,|o|extended_quantity(o.inverse_cdf(probability),d))});
object_binding!(poisson_dist_new,"num_native_poisson_dist_new","crate::kernels::stats::PoissonDist::new","Checked native legacy distribution. Continuous parameters and queries retain coherent physical units; inverse CDF restores sample units and PDF has inverse sample units.",true,(rate:f64)->NewObject<PoissonDist>,|_ctx:&CallCtx<'_>,rate:f64|{Ok(NewObject::normalized(PoissonDist::new(rate).map_err(error)?))});
object_binding!(poisson_dist_cdf,"num_native_poisson_dist_cdf","crate::kernels::stats::PoissonDist::cdf","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<PoissonDist>,k:u64)->f64,|ctx:&CallCtx<'_>,a:Handle<PoissonDist>,k:u64|{a.with(ctx,|o|Ok(o.cdf(k)))});
object_binding!(poisson_dist_pmf,"num_native_poisson_dist_pmf","crate::kernels::stats::PoissonDist::pmf","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<PoissonDist>,k:u64)->f64,|ctx:&CallCtx<'_>,a:Handle<PoissonDist>,k:u64|{a.with(ctx,|o|Ok(o.pmf(k)))});
object_binding!(binomial_dist_new,"num_native_binomial_dist_new","crate::kernels::stats::BinomialDist::new","Checked native legacy distribution. Continuous parameters and queries retain coherent physical units; inverse CDF restores sample units and PDF has inverse sample units.",true,(n:u64,p:f64)->NewObject<BinomialDist>,|_ctx:&CallCtx<'_>,n:u64,p:f64|{if n>1_048_576{return Err(error("binomial trial count exceeds native work budget"));}Ok(NewObject::normalized(BinomialDist::new(n,p).map_err(error)?))});
object_binding!(binomial_dist_cdf,"num_native_binomial_dist_cdf","crate::kernels::stats::BinomialDist::cdf","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<BinomialDist>,k:u64)->f64,|ctx:&CallCtx<'_>,a:Handle<BinomialDist>,k:u64|{a.with(ctx,|o|Ok(o.cdf(k)))});
object_binding!(binomial_dist_pmf,"num_native_binomial_dist_pmf","crate::kernels::stats::BinomialDist::pmf","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<BinomialDist>,k:u64)->f64,|ctx:&CallCtx<'_>,a:Handle<BinomialDist>,k:u64|{a.with(ctx,|o|Ok(o.pmf(k)))});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &normal_dist_new::BINDING,
            &normal_dist_cdf::BINDING,
            &normal_dist_pdf::BINDING,
            &normal_dist_inverse_cdf::BINDING,
            &uniform_dist_new::BINDING,
            &uniform_dist_cdf::BINDING,
            &uniform_dist_pdf::BINDING,
            &uniform_dist_inverse_cdf::BINDING,
            &exponential_dist_new::BINDING,
            &exponential_dist_cdf::BINDING,
            &exponential_dist_pdf::BINDING,
            &exponential_dist_inverse_cdf::BINDING,
            &gamma_dist_new::BINDING,
            &gamma_dist_cdf::BINDING,
            &gamma_dist_pdf::BINDING,
            &gamma_dist_inverse_cdf::BINDING,
            &poisson_dist_new::BINDING,
            &poisson_dist_cdf::BINDING,
            &poisson_dist_pmf::BINDING,
            &binomial_dist_new::BINDING,
            &binomial_dist_cdf::BINDING,
            &binomial_dist_pmf::BINDING,
        ],
        "numerical_vector",
    )
}
