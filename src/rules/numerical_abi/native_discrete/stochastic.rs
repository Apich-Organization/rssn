#![allow(clippy::wildcard_imports)] // Shared checked method declarations.
use super::*;
use crate::kernels::online_stats::OnlineCovariance;
use crate::kernels::random::Rng;
impl NativeType for Rng {
    const KIND: &'static str = "rssn/native/rng/v1";
    const MAX_BYTES: usize = 256;

    fn validate(&self) -> Result<(), OpError> {
        Ok(())
    }
}
impl NativeType for OnlineCovariance {
    const KIND: &'static str = "rssn/native/discrete/online-covariance/v1";
    const MAX_BYTES: usize = 262144;

    fn validate(&self) -> Result<(), OpError> {
        if self.dimension() == 0 || self.dimension() > 128 {
            return Err(error("online covariance dimension must be 1..=128"));
        }
        Ok(())
    }
}
object_binding!(rng_new,"num_native_rng_new","crate::kernels::random::Rng::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(seed:u64)->NewObject<Rng>,|_ctx:&CallCtx<'_>,seed:u64|{Ok(NewObject::normalized(Rng::new(seed)))});
object_binding!(rng_uniform,"num_native_rng_uniform","crate::kernels::random::Rng::uniform","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(rng:Handle<Rng>)->f64,|ctx:&CallCtx<'_>,rng:Handle<Rng>|{rng.with_mut(ctx,|r|Ok(r.uniform()))});
object_binding!(rng_normal,"num_native_rng_normal","crate::kernels::random::Rng::normal","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(rng:Handle<Rng>)->f64,|ctx:&CallCtx<'_>,rng:Handle<Rng>|{rng.with_mut(ctx,|r|Ok(r.normal()))});
object_binding!(rng_next_u64,"num_native_rng_next_u64","crate::kernels::random::Rng::next_u64","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(rng:Handle<Rng>)->u64,|ctx:&CallCtx<'_>,rng:Handle<Rng>|{rng.with_mut(ctx,|r|Ok(r.next_u64()))});
object_binding!(rng_jump,"num_native_rng_jump","crate::kernels::random::Rng::jump","Advance the session-owned random generator by 2^128 raw steps and clear its cached normal sample. Returns true when the generator was updated.",true,(rng:Handle<Rng>)->bool,|ctx:&CallCtx<'_>,rng:Handle<Rng>|{rng.with_mut(ctx,|r|{r.jump();Ok(true)})});
object_binding!(rng_long_jump,"num_native_rng_long_jump","crate::kernels::random::Rng::long_jump","Advance the session-owned random generator by 2^192 raw steps and clear its cached normal sample. Returns true when the generator was updated.",true,(rng:Handle<Rng>)->bool,|ctx:&CallCtx<'_>,rng:Handle<Rng>|{rng.with_mut(ctx,|r|{r.long_jump();Ok(true)})});
object_binding!(rng_below,"num_native_rng_below","crate::kernels::random::Rng::below","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(rng:Handle<Rng>,n:u64)->u64,|ctx:&CallCtx<'_>,rng:Handle<Rng>,n:u64|{if n==0{return Err(error("uniform integer range must be nonempty"));}rng.with_mut(ctx,|r|Ok(r.below(n)))});
object_binding!(rng_beta,"num_native_rng_beta","crate::kernels::random::Rng::beta","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(rng:Handle<Rng>,a:f64,b:f64)->f64,|ctx:&CallCtx<'_>,rng:Handle<Rng>,a:f64,b:f64|{if a<=0.0||b<=0.0||a>1e8||b>1e8{return Err(error("beta shapes must be positive and within native sampling budget"));}rng.with_mut(ctx,|r|Ok(r.beta(a,b)))});
object_binding!(rng_binomial,"num_native_rng_binomial","crate::kernels::random::Rng::binomial","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(rng:Handle<Rng>,n:u64,p:f64)->u64,|ctx:&CallCtx<'_>,rng:Handle<Rng>,n:u64,p:f64|{if n>1_048_576||!(0.0..=1.0).contains(&p){return Err(error("binomial count/probability exceeds sampling domain or work budget"));}rng.with_mut(ctx,|r|Ok(r.binomial(n,p)))});
object_binding!(rng_poisson,"num_native_rng_poisson","crate::kernels::random::Rng::poisson","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(rng:Handle<Rng>,lambda:f64)->u64,|ctx:&CallCtx<'_>,rng:Handle<Rng>,lambda:f64|{if !(0.0..=1e6).contains(&lambda){return Err(error("Poisson mean must be nonnegative and within sampling budget"));}rng.with_mut(ctx,|r|Ok(r.poisson(lambda)))});
object_binding!(rng_uniform_range,"num_native_rng_uniform_range","crate::kernels::random::Rng::uniform_range","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(rng:Handle<Rng>,lo:Input,hi:Input)->Dynamic,|ctx:&CallCtx<'_>,rng:Handle<Rng>,lo:Input,hi:Input|{let(a,d)=scalar(&lo)?;let(b,e)=scalar(&hi)?;api::same_dimension(&d,&e)?;if a>=b{return Err(error("uniform range endpoints must increase"));}rng.with_mut(ctx,|r|Ok(Dynamic(api::quantity(r.uniform_range(a,b),d)?)))});
object_binding!(rng_normal_with,"num_native_rng_normal_with","crate::kernels::random::Rng::normal_with","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(rng:Handle<Rng>,mean:Input,sd:Input)->Dynamic,|ctx:&CallCtx<'_>,rng:Handle<Rng>,mean:Input,sd:Input|{let(m,d)=scalar(&mean)?;let(s,e)=scalar(&sd)?;api::same_dimension(&d,&e)?;if s<0.0{return Err(error("normal standard deviation must be nonnegative"));}rng.with_mut(ctx,|r|Ok(Dynamic(api::quantity(r.normal_with(m,s),d)?)))});
object_binding!(rng_exponential,"num_native_rng_exponential","crate::kernels::random::Rng::exponential","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(rng:Handle<Rng>,rate:Input)->Dynamic,|ctx:&CallCtx<'_>,rng:Handle<Rng>,rate:Input|{let(r,d)=scalar(&rate)?;if r<=0.0{return Err(error("exponential rate must be positive"));}let d=canonical(&(DimExpr::Inverse(Box::new(d))))?;rng.with_mut(ctx,|g|Ok(Dynamic(api::quantity(g.exponential(r),d)?)))});
object_binding!(rng_gamma,"num_native_rng_gamma","crate::kernels::random::Rng::gamma","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(rng:Handle<Rng>,shape:f64,scale:Input)->Dynamic,|ctx:&CallCtx<'_>,rng:Handle<Rng>,shape:f64,scale:Input|{let(s,d)=scalar(&scale)?;if shape<=0.0||shape>1e8||s<=0.0{return Err(error("gamma shape/scale must be positive and shape within sampling budget"));}rng.with_mut(ctx,|r|Ok(Dynamic(api::quantity(r.gamma(shape,s),d)?)))});
object_binding!(rng_fill_normal,"num_native_rng_fill_normal","crate::kernels::random::Rng::fill_normal","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(rng:Handle<Rng>,n:usize)->Vec<f64>,|ctx:&CallCtx<'_>,rng:Handle<Rng>,n:usize|{limit(n,262144)?;let mut values=vec![0.0;n];rng.with_mut(ctx,|r|{r.fill_normal(&mut values);Ok(values)})});
object_binding!(rng_shuffle,"num_native_rng_shuffle","crate::kernels::random::Rng::shuffle","Fisher-Yates specialization over arbitrary owned kernel Values. Element representations, physical unit witnesses and opaque handles are preserved; only the order is mutated in local owned storage.",false,(rng:Handle<Rng>,data:Input)->Dynamic,|ctx:&CallCtx<'_>,rng:Handle<Rng>,data:Input|{let Value::Tensor(mut tensor)=data.0 else{return Err(error("shuffle requires a rank-one tensor"));};
if tensor.shape.len()!=1||tensor.shape[0]!=tensor.data.len(){return Err(error("shuffle tensor shape is invalid"));}rng.with_mut(ctx,|r|{r.shuffle(&mut tensor.data);Ok(Dynamic(Value::Tensor(tensor)))})});
object_binding!(online_covariance_new,"num_native_online_covariance_new","crate::kernels::online_stats::OnlineCovariance::new","Create an online accumulator for homogeneous observation units specified by a positive reference quantity. Coherent means retain those units and covariance squares them.",false,(n:usize,unit:Input)->NewObject<OnlineCovariance>,|_ctx:&CallCtx<'_>,n:usize,unit:Input|{limit(n,128)?;let(_,d)=api::reference(&unit.0)?;Ok(NewObject::new(OnlineCovariance::new(n).map_err(error)?,vec![d]))});
object_binding!(online_covariance_count,"num_native_online_covariance_count","crate::kernels::online_stats::OnlineCovariance::count","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<OnlineCovariance>)->u64,|ctx:&CallCtx<'_>,a:Handle<OnlineCovariance>|{a.with(ctx,|o|Ok(o.count()))});
object_binding!(online_covariance_dimension,"num_native_online_covariance_dimension","crate::kernels::online_stats::OnlineCovariance::dimension","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<OnlineCovariance>)->usize,|ctx:&CallCtx<'_>,a:Handle<OnlineCovariance>|{a.with(ctx,|o|Ok(o.dimension()))});
object_binding!(online_covariance_push,"num_native_online_covariance_push","crate::kernels::online_stats::OnlineCovariance::push","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<OnlineCovariance>,data:Input)->bool,|ctx:&CallCtx<'_>,a:Handle<OnlineCovariance>,data:Input|{let d=object_dimension(&a,ctx)?;api::same_dimension(&d,&api::dimension(&data.0)?)?;let data=<Vec<f64> as Arg>::decode(&data.0,1)?;a.with_mut(ctx,|o|{o.push(&data).map_err(error)?;Ok(true)})});
object_binding!(online_covariance_merge,"num_native_online_covariance_merge","crate::kernels::online_stats::OnlineCovariance::merge","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<OnlineCovariance>,b:Handle<OnlineCovariance>)->bool,|ctx:&CallCtx<'_>,a:Handle<OnlineCovariance>,b:Handle<OnlineCovariance>|{same(&a,&b,ctx)?;let other=b.with(ctx,|o|Ok(o.clone()))?;a.with_mut(ctx,|o|{o.merge(&other).map_err(error)?;Ok(true)})});
object_binding!(online_covariance_mean,"num_native_online_covariance_mean","crate::kernels::online_stats::OnlineCovariance::mean","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<OnlineCovariance>)->Maybe<Dynamic>,|ctx:&CallCtx<'_>,a:Handle<OnlineCovariance>|{let d=object_dimension(&a,ctx)?;a.with(ctx,|o|Ok(Maybe(o.mean().map(|v|quantities(v.to_vec(),&(d))).transpose()?)))});
object_binding!(online_covariance_sample_covariance,"num_native_online_covariance_sample_covariance","crate::kernels::online_stats::OnlineCovariance::sample_covariance","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<OnlineCovariance>)->Maybe<Dynamic>,|ctx:&CallCtx<'_>,a:Handle<OnlineCovariance>|{let d=canonical(&(object_dimension(&a,ctx)?.pow(Exp::Rat(Rational::int(2)))))?;a.with(ctx,|o|{let n=o.dimension();Ok(Maybe(o.sample_covariance().map(|v|quantities_matrix(v.chunks_exact(n).map(<[f64]>::to_vec).collect(),&(d))).transpose()?))})});
object_binding!(online_covariance_population_covariance,"num_native_online_covariance_population_covariance","crate::kernels::online_stats::OnlineCovariance::population_covariance","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<OnlineCovariance>)->Maybe<Dynamic>,|ctx:&CallCtx<'_>,a:Handle<OnlineCovariance>|{let d=canonical(&(object_dimension(&a,ctx)?.pow(Exp::Rat(Rational::int(2)))))?;a.with(ctx,|o|{let n=o.dimension();Ok(Maybe(o.population_covariance().map(|v|quantities_matrix(v.chunks_exact(n).map(<[f64]>::to_vec).collect(),&(d))).transpose()?))})});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &rng_new::BINDING,
            &rng_uniform::BINDING,
            &rng_normal::BINDING,
            &rng_next_u64::BINDING,
            &rng_jump::BINDING,
            &rng_long_jump::BINDING,
            &rng_below::BINDING,
            &rng_beta::BINDING,
            &rng_binomial::BINDING,
            &rng_poisson::BINDING,
            &rng_uniform_range::BINDING,
            &rng_normal_with::BINDING,
            &rng_exponential::BINDING,
            &rng_gamma::BINDING,
            &rng_fill_normal::BINDING,
            &rng_shuffle::BINDING,
            &online_covariance_new::BINDING,
            &online_covariance_count::BINDING,
            &online_covariance_dimension::BINDING,
            &online_covariance_push::BINDING,
            &online_covariance_merge::BINDING,
            &online_covariance_mean::BINDING,
            &online_covariance_sample_covariance::BINDING,
            &online_covariance_population_covariance::BINDING,
        ],
        "numerical_vector",
    )
}
