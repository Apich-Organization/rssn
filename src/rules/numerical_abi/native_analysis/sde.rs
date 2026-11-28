//! Native SDE paths and explicit Brownian samples in independently scaled physical coordinates.
use crate::graph::RuleError;
use crate::graph::rule::Installer;
use crate::kernel::CallCtx;
use crate::kernel::OpError;
use crate::kernel::Value;
use crate::kernel::signature::DimExpr;
use crate::kernel::signature::Exp;
use crate::kernel::signature::Rational;
use crate::kernels::random::Rng;
use crate::kernels::sde as s;
use crate::rules::numerical_abi::codecs::Arg;
use crate::rules::numerical_abi::codecs::Output;
use crate::rules::numerical_abi::native_analysis_newton_cg::parameters;
use crate::rules::numerical_abi::native_analysis_ode_common::Budget;
use crate::rules::numerical_abi::native_analysis_ode_common::OdeFrame;
use crate::rules::numerical_abi::native_analysis_ode_common::OdeTape;
use crate::rules::numerical_abi::native_backend_extended::TapeSnapshot;
use crate::rules::numerical_abi::native_backend_extended::tape_snapshot;
use crate::rules::numerical_abi::native_object_api::Dynamic;
use crate::rules::numerical_abi::native_object_api::Handle;
use crate::rules::numerical_abi::native_object_api::Input;
use crate::rules::numerical_abi::native_object_api::Many;
use crate::rules::numerical_abi::native_object_api::NativeType;
use crate::rules::numerical_abi::native_object_api::NewObject;
use crate::rules::numerical_abi::native_object_api::object_binding;
use crate::rules::numerical_abi::native_object_api::{
    self as api,
};
use std::sync::Arc;
#[path = "sde_extended.rs"]
mod extended;
#[path = "sde_stochastic.rs"]
mod stochastic;
fn half(d: DimExpr) -> DimExpr {
    DimExpr::Pow(Box::new(d), Exp::Rat(Rational::new(1, 2)))
}
#[derive(Clone)]
struct Model {
    drift: OdeTape,
    diffusion: TapeSnapshot,
    params: Vec<f64>,
    m: usize,
    input_ratios: Vec<f64>,
    output_ratios: Vec<f64>,
}
impl NativeType for Model {
    const KIND: &'static str = "rssn/native/analysis/sde/model/v1";
    const MAX_BYTES: usize = 8 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        self.drift.validate()?;
        if self.m == 0
            || self.m > 16
            || self.output_ratios.len() != self.drift.frame.state_scales.len() * self.m
        {
            return Err(api::error("SDE noise dimension/shape budget mismatch"));
        }
        Ok(())
    }
}
impl Model {
    fn diffusion(
        &self,
        budget: &Budget,
        t: f64,
        y: &[f64],
        out: &mut [f64],
    ) {
        out.fill(f64::NAN);
        if out.len() != self.output_ratios.len() || y.len() != self.drift.frame.state_scales.len() {
            budget.fail(api::error("SDE diffusion shape mismatch"));
            return;
        }
        if !budget.charge() {
            return;
        }
        let x = std::iter::once(t)
            .chain(y.iter().copied())
            .zip(&self.input_ratios)
            .map(|(v, s)| api::scale_product(v, &[*s], &[]))
            .collect::<Result<Vec<_>, _>>();
        let x = match x {
            | Ok(x) => x,
            | Err(e) => {
                budget.fail(e);
                return;
            },
        };
        self.diffusion.tape.eval(&x, &self.params, out);
        for (v, s) in out.iter_mut().zip(&self.output_ratios) {
            match api::scale_product(*v, &[*s], &[]) {
                | Ok(x) => *v = x,
                | Err(e) => {
                    budget.fail(e);
                    *v = f64::NAN;
                },
            }
        }
    }
}
#[derive(Clone)]
struct Noise {
    sample: s::NoiseSample,
    m: usize,
    time_scale: f64,
    time_dimension: DimExpr,
    h: f64,
}
impl NativeType for Noise {
    const KIND: &'static str = "rssn/native/analysis/sde/noise/v1";
    const MAX_BYTES: usize = 8 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        noise_limit(self.m, self.sample.dws.len(), 1)?;
        if !self.h.is_finite()
            || self.h <= 0.0
            || self.sample.areas.len() != self.sample.dws.len()
            || self
                .sample
                .dws
                .iter()
                .any(|r| r.len() != self.m || r.iter().any(|x| !x.is_finite()))
            || self
                .sample
                .areas
                .iter()
                .any(|r| r.len() != self.m * self.m || r.iter().any(|x| !x.is_finite()))
        {
            return Err(api::error("invalid native Brownian/Levy sample shape"));
        }
        Ok(())
    }
}
struct Path {
    frame: OdeFrame,
    data: s::SdePath,
    evaluations: usize,
}
impl NativeType for Path {
    const KIND: &'static str = "rssn/native/analysis/sde/path/v1";
    const MAX_BYTES: usize = 8 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        if self.data.t.len() != self.data.y.len()
            || self
                .data
                .y
                .len()
                .checked_mul(self.frame.state_scales.len())
                .is_none_or(|n| n > 262144)
            || self.data.y.iter().any(|y| {
                y.len() != self.frame.state_scales.len() || y.iter().any(|v| !v.is_finite())
            })
        {
            return Err(api::error("SDE path storage/shape budget exceeded"));
        }
        Ok(())
    }
}
fn noise_limit(
    m: usize,
    steps: usize,
    terms: usize,
) -> Result<(), OpError> {
    if m == 0
        || m > 16
        || steps == 0
        || terms == 0
        || terms > 128
        || steps
            .checked_mul(m)
            .and_then(|n| n.checked_mul(m + terms))
            .is_none_or(|n| n > 262144)
    {
        Err(api::error("Brownian/Levy allocation/work budget exceeded"))
    } else {
        Ok(())
    }
}
object_binding!(model_new,"num_analysis_sde_model_new","","Owned native drift/diffusion tapes. Diffusion outputs are row-major n by m and have component dimensions state/sqrt(time); all input/output reference scales are explicit and may differ from the drift tape frame.",false,(drift:Handle<OdeTape>,diffusion:Input,params:Input,m:usize)->NewObject<Model>,|ctx,drift,diffusion,params,m|{let drift=drift.with(ctx,|d|Ok(d.clone()))?;let snapshot=tape_snapshot(ctx,&diffusion.0)?;let frame=&drift.frame;let n=frame.state_scales.len();if m==0||m>16||snapshot.input_scales.len()!=n+1||snapshot.output_scales.len()!=n*m{return Err(api::error("SDE diffusion tape shape mismatch"));}let mut dims=vec![frame.time_dimension.clone()];dims.extend(frame.state_dimensions.clone());let mut scales=vec![frame.time_scale];scales.extend(frame.state_scales.clone());let input_ratios=scales.iter().zip(&snapshot.input_scales).zip(dims.iter().zip(&snapshot.input_dimensions)).map(|((a,b),(da,db))|{api::same_dimension(da,db)?;api::scale_product(*a,&[],&[*b])}).collect::<Result<Vec<_>,OpError>>()?;let output_ratios=snapshot.output_scales.iter().enumerate().map(|(k,scale)|{let i=k/m;api::same_dimension(&snapshot.output_dimensions[k],&DimExpr::Mul(vec![frame.state_dimensions[i].clone(),DimExpr::Inverse(Box::new(half(frame.time_dimension.clone())))]))?;api::scale_product(*scale,&[frame.time_scale.sqrt()],&[frame.state_scales[i]])}).collect::<Result<Vec<_>,OpError>>()?;let params=parameters(&params,&snapshot.parameter_scales,&snapshot.parameter_dimensions)?;let dimensions=frame.dimensions();Ok(NewObject::new(Model{drift,diffusion:snapshot,params,m,input_ratios,output_ratios},dimensions))});
object_binding!(sample,"num_analysis_sde_sample_noise","crate::kernels::sde::sample_noise","Actual native conditional Brownian increments/Levy areas. Positive time reference defines normalized stochastic time; physical h has the same dimension. RNG state advances exactly through the native routine.",false,(rng:Handle<Rng>,m:usize,steps:usize,h:Input,time_reference:Input,terms:usize)->NewObject<Noise>,|ctx,rng,m,steps,h,time_reference,terms|{noise_limit(m,steps,terms)?;let(time_scale,time_dimension)=api::reference(&time_reference.0)?;api::same_dimension(&time_dimension,&api::dimension(&h.0)?)?;let h=api::scale_product(crate::rules::numerical_abi::numeric(&h.0,0)?,&[],&[time_scale])?;if h<=0.0{return Err(api::error("SDE step must be positive"));}let sample=rng.with_mut(ctx,|r|Ok(s::sample_noise(r,m,steps,h,terms)))?;Ok(NewObject::new(Noise{sample,m,time_scale,time_dimension:time_dimension.clone(),h},vec![time_dimension]))});
object_binding!(aggregate,"num_analysis_sde_aggregate_noise","crate::kernels::sde::aggregate_noise","Actual native Chen aggregation retains the same physical Brownian path and iterated integrals.",false,(noise:Handle<Noise>,group:usize)->NewObject<Noise>,|ctx,noise,group|noise.with(ctx,|n|{if group==0||group>n.sample.dws.len(){return Err(api::error("noise aggregation group out of range"));}
if !n.sample.dws.len().is_multiple_of(group){return Err(api::error("fixed-step noise aggregation requires complete equal groups"));}let sample=s::aggregate_noise(&n.sample,group,n.m);let h=api::scale_product(n.h,&[group as f64],&[])?;Ok(NewObject::new(Noise{sample,m:n.m,h,time_scale:n.time_scale,time_dimension:n.time_dimension.clone()},vec![n.time_dimension.clone()]))}));
#[allow(
    clippy::needless_pass_by_value,
    reason = "The typed native ABI transfers decoded arguments by value; this shared dispatcher mirrors its explicit source signature."
)]
fn fixed(
    ctx: &CallCtx<'_>,
    model: Handle<Model>,
    noise: Handle<Noise>,
    initial: Input,
    start: Input,
    maxeval: usize,
    kind: usize,
) -> Result<Many<Dynamic>, OpError> {
    let model = Arc::new(model.with(ctx, |m| Ok(m.clone()))?);
    let noise = noise.with(ctx, |n| Ok(n.clone()))?;
    let frame = &model.drift.frame;
    api::same_dimension(&frame.time_dimension, &noise.time_dimension)?;
    if noise.m != model.m {
        return Err(api::error("SDE noise dimension mismatch"));
    }
    let y = frame.state(&initial)?;
    let t = frame.time(&start)?;
    if noise
        .sample
        .dws
        .len()
        .checked_mul(y.len())
        .is_none_or(|n| n > 262000)
    {
        return Err(api::error("SDE path storage budget exceeded"));
    }
    let h = api::scale_product(noise.h, &[noise.time_scale], &[frame.time_scale])?;
    let ratio = api::scale_product(noise.time_scale.sqrt(), &[], &[frame.time_scale.sqrt()])?;
    let area_ratio = api::scale_product(noise.time_scale, &[], &[frame.time_scale])?;
    let mut sample = noise.sample;
    for row in &mut sample.dws {
        for v in row {
            *v = api::scale_product(*v, &[ratio], &[])?;
        }
    }
    for row in &mut sample.areas {
        for v in row {
            *v = api::scale_product(*v, &[area_ratio], &[])?;
        }
    }
    let budget = Budget::new(maxeval)?;
    let drift = |t, y: &[f64], out: &mut [f64]| model.drift.rate(&budget, t, y, out);
    let diff = |t, y: &[f64], out: &mut [f64]| model.diffusion(&budget, t, y, out);
    let m = model.m;
    let result = match kind {
        | 0 => s::euler_maruyama(drift, diff, m, t, h, &y, &sample.dws),
        | 1 => s::heun_stratonovich(drift, diff, m, t, h, &y, &sample.dws),
        | 2 => s::milstein(drift, diff, m, t, h, &y, &sample.dws),
        | 3 => s::milstein_stratonovich(drift, diff, m, t, h, &y, &sample.dws),
        | 4 => s::srk_platen(drift, diff, m, t, h, &y, &sample.dws),
        | 5 => s::milstein_noncommutative(drift, diff, m, t, h, &y, &sample),
        | 6 => s::milstein_noncommutative_stratonovich(drift, diff, m, t, h, &y, &sample),
        | 7 => s::srk_noncommutative(drift, diff, m, t, h, &y, &sample),
        | 8 => s::srk_noncommutative_stratonovich(drift, diff, m, t, h, &y, &sample),
        | _ => return Err(api::error("unknown fixed SDE method")),
    };
    budget.finish()?;
    Ok(Many(
        result
            .map_err(|e| api::error(format!("{e:?}")))?
            .iter()
            .map(|y| frame.restore(y).map(Dynamic))
            .collect::<Result<Vec<_>, _>>()?,
    ))
}
object_binding!(euler_maruyama,"num_analysis_sde_euler_maruyama","crate::kernels::sde::euler_maruyama","Actual native euler_maruyama over one retained physical Brownian sample. The caller explicitly selects Ito/Stratonovich and commutative/noncommutative methods; no scheme substitution. Independent state dimensions and a shared native evaluation budget are preserved.",false,(model:Handle<Model>,noise:Handle<Noise>,initial:Input,start:Input,maxeval:usize)->Many<Dynamic>,|ctx,model,noise,initial,start,maxeval|fixed(ctx,model,noise,initial,start,maxeval,0));
object_binding!(heun_stratonovich,"num_analysis_sde_heun_stratonovich","crate::kernels::sde::heun_stratonovich","Actual native heun_stratonovich over one retained physical Brownian sample. The caller explicitly selects Ito/Stratonovich and commutative/noncommutative methods; no scheme substitution. Independent state dimensions and a shared native evaluation budget are preserved.",false,(model:Handle<Model>,noise:Handle<Noise>,initial:Input,start:Input,maxeval:usize)->Many<Dynamic>,|ctx,model,noise,initial,start,maxeval|fixed(ctx,model,noise,initial,start,maxeval,1));
object_binding!(milstein,"num_analysis_sde_milstein","crate::kernels::sde::milstein","Actual native milstein over one retained physical Brownian sample. The caller explicitly selects Ito/Stratonovich and commutative/noncommutative methods; no scheme substitution. Independent state dimensions and a shared native evaluation budget are preserved.",false,(model:Handle<Model>,noise:Handle<Noise>,initial:Input,start:Input,maxeval:usize)->Many<Dynamic>,|ctx,model,noise,initial,start,maxeval|fixed(ctx,model,noise,initial,start,maxeval,2));
object_binding!(milstein_stratonovich,"num_analysis_sde_milstein_stratonovich","crate::kernels::sde::milstein_stratonovich","Actual native milstein_stratonovich over one retained physical Brownian sample. The caller explicitly selects Ito/Stratonovich and commutative/noncommutative methods; no scheme substitution. Independent state dimensions and a shared native evaluation budget are preserved.",false,(model:Handle<Model>,noise:Handle<Noise>,initial:Input,start:Input,maxeval:usize)->Many<Dynamic>,|ctx,model,noise,initial,start,maxeval|fixed(ctx,model,noise,initial,start,maxeval,3));
object_binding!(srk_platen,"num_analysis_sde_srk_platen","crate::kernels::sde::srk_platen","Actual native srk_platen over one retained physical Brownian sample. The caller explicitly selects Ito/Stratonovich and commutative/noncommutative methods; no scheme substitution. Independent state dimensions and a shared native evaluation budget are preserved.",false,(model:Handle<Model>,noise:Handle<Noise>,initial:Input,start:Input,maxeval:usize)->Many<Dynamic>,|ctx,model,noise,initial,start,maxeval|fixed(ctx,model,noise,initial,start,maxeval,4));
object_binding!(milstein_noncommutative,"num_analysis_sde_milstein_noncommutative","crate::kernels::sde::milstein_noncommutative","Actual native milstein_noncommutative over one retained physical Brownian sample. The caller explicitly selects Ito/Stratonovich and commutative/noncommutative methods; no scheme substitution. Independent state dimensions and a shared native evaluation budget are preserved.",false,(model:Handle<Model>,noise:Handle<Noise>,initial:Input,start:Input,maxeval:usize)->Many<Dynamic>,|ctx,model,noise,initial,start,maxeval|fixed(ctx,model,noise,initial,start,maxeval,5));
object_binding!(milstein_noncommutative_stratonovich,"num_analysis_sde_milstein_noncommutative_stratonovich","crate::kernels::sde::milstein_noncommutative_stratonovich","Actual native milstein_noncommutative_stratonovich over one retained physical Brownian sample. The caller explicitly selects Ito/Stratonovich and commutative/noncommutative methods; no scheme substitution. Independent state dimensions and a shared native evaluation budget are preserved.",false,(model:Handle<Model>,noise:Handle<Noise>,initial:Input,start:Input,maxeval:usize)->Many<Dynamic>,|ctx,model,noise,initial,start,maxeval|fixed(ctx,model,noise,initial,start,maxeval,6));
object_binding!(srk_noncommutative,"num_analysis_sde_srk_noncommutative","crate::kernels::sde::srk_noncommutative","Actual native srk_noncommutative over one retained physical Brownian sample. The caller explicitly selects Ito/Stratonovich and commutative/noncommutative methods; no scheme substitution. Independent state dimensions and a shared native evaluation budget are preserved.",false,(model:Handle<Model>,noise:Handle<Noise>,initial:Input,start:Input,maxeval:usize)->Many<Dynamic>,|ctx,model,noise,initial,start,maxeval|fixed(ctx,model,noise,initial,start,maxeval,7));
object_binding!(srk_noncommutative_stratonovich,"num_analysis_sde_srk_noncommutative_stratonovich","crate::kernels::sde::srk_noncommutative_stratonovich","Actual native srk_noncommutative_stratonovich over one retained physical Brownian sample. The caller explicitly selects Ito/Stratonovich and commutative/noncommutative methods; no scheme substitution. Independent state dimensions and a shared native evaluation budget are preserved.",false,(model:Handle<Model>,noise:Handle<Noise>,initial:Input,start:Input,maxeval:usize)->Many<Dynamic>,|ctx,model,noise,initial,start,maxeval|fixed(ctx,model,noise,initial,start,maxeval,8));
#[allow(
    clippy::needless_pass_by_value,
    clippy::too_many_arguments,
    reason = "The typed native ABI transfers decoded arguments by value; this shared dispatcher mirrors its explicit source signature."
)]
#[allow(
    deprecated,
    reason = "Expose the actual legacy native entry as requested, alongside its modern replacement."
)]
fn adaptive(
    ctx: &CallCtx<'_>,
    model: Handle<Model>,
    initial: Input,
    start: Input,
    end: Input,
    config: Input,
    maxeval: usize,
    kind: usize,
) -> Result<NewObject<Path>, OpError> {
    let model = model.with(ctx, |m| Ok(m.clone()))?;
    let frame = &model.drift.frame;
    let y = frame.state(&initial)?;
    let t0 = frame.time(&start)?;
    let t1 = frame.time(&end)?;
    let f = |name| {
        let v = api::field(&config.0, name)?;
        api::normalized(v)?;
        <f64 as Arg>::decode(v, 0)
    };
    let count = |name| <usize as Arg>::decode(api::field(&config.0, name)?, 0);
    let o = s::SdeAdaptiveOptions {
        rtol: f("rtol")?,
        atol: f("atol")?,
        seed: <u64 as Arg>::decode(api::field(&config.0, "seed")?, 0)?,
        max_steps: count("max_steps")?,
    };
    if o.max_steps == 0 || o.max_steps.checked_mul(y.len()).is_none_or(|n| n > 262000) {
        return Err(api::error("adaptive SDE path storage budget exceeded"));
    }
    let nc = if kind >= 3 {
        let cells = count("cells")?;
        let terms = count("terms")?;
        noise_limit(model.m, cells, terms)?;
        Some(s::SdeNoncommutativeOptions {
            base: o,
            cells,
            terms,
        })
    } else {
        None
    };
    let budget = Budget::new(maxeval)?;
    let drift = |t, y: &[f64], out: &mut [f64]| model.drift.rate(&budget, t, y, out);
    let diff = |t, y: &[f64], out: &mut [f64]| model.diffusion(&budget, t, y, out);
    let m = model.m;
    let result = match kind {
        | 0 => {
            s::sde_adaptive(
                drift,
                diff,
                m,
                t0,
                t1,
                &y,
                o.rtol,
                o.atol,
                o.seed,
                o.max_steps,
            )
        },
        | 1 => s::sde_adaptive_with(drift, diff, m, t0, t1, &y, &o),
        | 2 => s::sde_adaptive_stratonovich_with(drift, diff, m, t0, t1, &y, &o),
        | 3 => {
            s::sde_adaptive_noncommutative_with(drift, diff, m, t0, t1, &y, nc.as_ref().unwrap())
        },
        | 4 => {
            s::sde_adaptive_noncommutative_stratonovich_with(
                drift,
                diff,
                m,
                t0,
                t1,
                &y,
                nc.as_ref().unwrap(),
            )
        },
        | _ => return Err(api::error("unknown adaptive SDE entry")),
    };
    budget.finish()?;
    Ok(NewObject::new(
        Path {
            frame: frame.clone(),
            data: result.map_err(|e| api::error(format!("{e:?}")))?,
            evaluations: budget.used(),
        },
        frame.dimensions(),
    ))
}
object_binding!(sde_adaptive,"num_analysis_sde_adaptive","crate::kernels::sde::sde_adaptive","Actual native sde_adaptive with independently scaled physical state, normalized rtol/atol, seed, max_steps, and explicit evaluation budget. Noncommutative variants require bounded cells/terms and retain genuine conditional Levy areas.",false,(model:Handle<Model>,initial:Input,start:Input,end:Input,config:Input,maxeval:usize)->NewObject<Path>,|ctx,model,initial,start,end,config,maxeval|adaptive(ctx,model,initial,start,end,config,maxeval,0));
object_binding!(sde_adaptive_with,"num_analysis_sde_adaptive_with","crate::kernels::sde::sde_adaptive_with","Actual native sde_adaptive_with with independently scaled physical state, normalized rtol/atol, seed, max_steps, and explicit evaluation budget. Noncommutative variants require bounded cells/terms and retain genuine conditional Levy areas.",false,(model:Handle<Model>,initial:Input,start:Input,end:Input,config:Input,maxeval:usize)->NewObject<Path>,|ctx,model,initial,start,end,config,maxeval|adaptive(ctx,model,initial,start,end,config,maxeval,1));
object_binding!(sde_adaptive_stratonovich_with,"num_analysis_sde_adaptive_stratonovich_with","crate::kernels::sde::sde_adaptive_stratonovich_with","Actual native sde_adaptive_stratonovich_with with independently scaled physical state, normalized rtol/atol, seed, max_steps, and explicit evaluation budget. Noncommutative variants require bounded cells/terms and retain genuine conditional Levy areas.",false,(model:Handle<Model>,initial:Input,start:Input,end:Input,config:Input,maxeval:usize)->NewObject<Path>,|ctx,model,initial,start,end,config,maxeval|adaptive(ctx,model,initial,start,end,config,maxeval,2));
object_binding!(sde_adaptive_noncommutative_with,"num_analysis_sde_adaptive_noncommutative_with","crate::kernels::sde::sde_adaptive_noncommutative_with","Actual native sde_adaptive_noncommutative_with with independently scaled physical state, normalized rtol/atol, seed, max_steps, and explicit evaluation budget. Noncommutative variants require bounded cells/terms and retain genuine conditional Levy areas.",false,(model:Handle<Model>,initial:Input,start:Input,end:Input,config:Input,maxeval:usize)->NewObject<Path>,|ctx,model,initial,start,end,config,maxeval|adaptive(ctx,model,initial,start,end,config,maxeval,3));
object_binding!(sde_adaptive_noncommutative_stratonovich_with,"num_analysis_sde_adaptive_noncommutative_stratonovich_with","crate::kernels::sde::sde_adaptive_noncommutative_stratonovich_with","Actual native sde_adaptive_noncommutative_stratonovich_with with independently scaled physical state, normalized rtol/atol, seed, max_steps, and explicit evaluation budget. Noncommutative variants require bounded cells/terms and retain genuine conditional Levy areas.",false,(model:Handle<Model>,initial:Input,start:Input,end:Input,config:Input,maxeval:usize)->NewObject<Path>,|ctx,model,initial,start,end,config,maxeval|adaptive(ctx,model,initial,start,end,config,maxeval,4));
object_binding!(data,"num_analysis_sde_path_data","","Physical native accepted samples and original rejection/seed statistics; no interpolant is invented for stochastic paths.",false,(path:Handle<Path>)->Dynamic,|ctx,path|path.with(ctx,|p|{let times=Value::Tensor(crate::kernel::value::TensorValue{shape:vec![p.data.t.len()],elem:crate::kernel::signature::TypeExpr::Dim(Box::new(crate::kernel::signature::TypeExpr::Float),p.frame.time_dimension.clone()),data:p.data.t.iter().map(|v|api::scale_product(*v,&[p.frame.time_scale],&[]).map(Value::Float)).collect::<Result<_,_>>()?});let states=Value::Tensor(crate::kernel::value::TensorValue{shape:vec![p.data.y.len()],elem:api::dynamic_type(),data:p.data.y.iter().map(|y|p.frame.restore(y)).collect::<Result<_,_>>()?});Ok(Dynamic(api::record(vec![("times",times),("states",states),("accepted",Output::encode(p.data.accepted)?),("rejected",Output::encode(p.data.rejected)?),("seed",Output::encode(p.data.seed)?),("evaluations",Output::encode(p.evaluations)?)])))}));
object_binding!(ensemble_seed,"num_analysis_sde_ensemble_seed","crate::kernels::sde::ensemble_seed","Actual deterministic native ensemble member seed.",true,(seed:u64,index:usize)->u64,|_ctx,seed,index|Ok(s::ensemble_seed(seed,index)));

object_binding!(brownian,"num_analysis_sde_brownian_increments","crate::kernels::sde::brownian_increments","Actual native Brownian increments in physical sqrt(time) units. h must be strictly positive; shape and allocation budgets are checked.",false,(rng:Handle<Rng>,m:usize,steps:usize,h:Input)->Dynamic,|ctx,rng,m,steps,h|{noise_limit(m,steps,1)?;let time_dimension=api::dimension(&h.0)?;let h=crate::rules::numerical_abi::numeric(&h.0,0)?;if h<=0.0{return Err(api::error("Brownian step must be positive"));}let values=rng.with_mut(ctx,|r|Ok(s::brownian_increments(r,m,steps,h)))?;let data=values.into_iter().flatten().map(Value::Float).collect();Ok(Dynamic(Value::Tensor(crate::kernel::value::TensorValue{shape:vec![steps,m],elem:crate::kernel::signature::TypeExpr::Dim(Box::new(crate::kernel::signature::TypeExpr::Float),half(time_dimension)),data})))});
object_binding!(levy,"num_analysis_sde_levy_areas","crate::kernels::sde::levy_areas","Actual native conditional Levy areas. Supplied increments have sqrt(time) units and output flattened m by m Ito integrals have time units.",false,(rng:Handle<Rng>,m:usize,h:Input,dw:Input,terms:usize)->Dynamic,|ctx,rng,m,h,dw,terms|{noise_limit(m,1,terms)?;let dimension=api::dimension(&h.0)?;api::same_dimension(&half(dimension.clone()),&api::dimension(&dw.0)?)?;let dw=<Vec<f64> as Arg>::decode(&dw.0,0)?;let h=crate::rules::numerical_abi::numeric(&h.0,0)?;if h<=0.0||dw.len()!=m{return Err(api::error("Levy areas need positive h and m increments"));}let data=rng.with_mut(ctx,|r|Ok(s::levy_areas(r,m,h,&dw,terms)))?;Ok(Dynamic(Value::Tensor(crate::kernel::value::TensorValue{shape:vec![m,m],elem:crate::kernel::signature::TypeExpr::Dim(Box::new(crate::kernel::signature::TypeExpr::Float),dimension),data:data.into_iter().map(Value::Float).collect()})))});
object_binding!(noise_data,"num_analysis_sde_noise_data","","Physical Brownian/Levy arrays from a retained sample; increments have sqrt(time) units and areas have time units.",false,(noise:Handle<Noise>)->Dynamic,|ctx,noise|noise.with(ctx,|n|{let dw=n.sample.dws.iter().flatten().map(|v|api::scale_product(*v,&[n.time_scale.sqrt()],&[]).map(Value::Float)).collect::<Result<Vec<_>,_>>()?;let areas=n.sample.areas.iter().flatten().map(|v|api::scale_product(*v,&[n.time_scale],&[]).map(Value::Float)).collect::<Result<Vec<_>,_>>()?;Ok(Dynamic(api::record(vec![("increments",Value::Tensor(crate::kernel::value::TensorValue{shape:vec![n.sample.dws.len(),n.m],elem:crate::kernel::signature::TypeExpr::Dim(Box::new(crate::kernel::signature::TypeExpr::Float),half(n.time_dimension.clone())),data:dw})),("areas",Value::Tensor(crate::kernel::value::TensorValue{shape:vec![n.sample.areas.len(),n.m,n.m],elem:crate::kernel::signature::TypeExpr::Dim(Box::new(crate::kernel::signature::TypeExpr::Float),n.time_dimension.clone()),data:areas}))]))) }));
object_binding!(grid,"num_analysis_sde_noncommutative_noise_grid","crate::kernels::sde::noncommutative_noise_grid","Actual native reproducible noncommutative grid. The normalized [0,1] reference frame is explicitly a physical time reference equal to the requested positive span; state-independent noise retains genuine conditional areas.",false,(m:usize,start:Input,end:Input,cells:usize,terms:usize,seed:u64)->NewObject<Noise>,|_ctx,m,start,end,cells,terms,seed|{noise_limit(m,cells,terms)?;let dimension=api::dimension(&start.0)?;api::same_dimension(&dimension,&api::dimension(&end.0)?)?;let t0=crate::rules::numerical_abi::numeric(&start.0,0)?;let t1=crate::rules::numerical_abi::numeric(&end.0,0)?;let scale=t1-t0;if !scale.is_finite()||scale<=0.0{return Err(api::error("native noise grid needs a finite positive span"));}let options=s::SdeNoncommutativeOptions{base:s::SdeAdaptiveOptions{rtol:1e-6,atol:1e-9,seed,max_steps:cells},cells,terms};let sample=s::noncommutative_noise_grid(m,0.0,1.0,&options);Ok(NewObject::new(Noise{sample,m,time_scale:scale,time_dimension:dimension.clone(),h:1.0/cells as f64},vec![dimension]))});
fn ensemble_scheme(
    name: &str,
    terms: usize,
) -> Result<s::EnsembleScheme, OpError> {
    if terms == 0 || terms > 128 {
        return Err(api::error("ensemble Levy terms must be1..128"));
    }
    Ok(match name {
        | "EulerMaruyama" => s::EnsembleScheme::EulerMaruyama,
        | "Milstein" => s::EnsembleScheme::Milstein,
        | "SrkPlaten" => s::EnsembleScheme::SrkPlaten,
        | "HeunStratonovich" => s::EnsembleScheme::HeunStratonovich,
        | "MilsteinNoncommutative" => s::EnsembleScheme::MilsteinNoncommutative { terms },
        | "MilsteinNoncommutativeStratonovich" => {
            s::EnsembleScheme::MilsteinNoncommutativeStratonovich { terms }
        },
        | "SrkNoncommutative" => s::EnsembleScheme::SrkNoncommutative { terms },
        | "SrkNoncommutativeStratonovich" => {
            s::EnsembleScheme::SrkNoncommutativeStratonovich { terms }
        },
        | _ => return Err(api::error("unknown native EnsembleScheme variant")),
    })
}
object_binding!(ensemble,"num_analysis_sde_ensemble","crate::kernels::sde::sde_ensemble","Actual parallel native ensemble on immutable owned drift/diffusion tapes. The explicit global evaluation budget is shared atomically across real parallel members; original failures are preserved. Returned final states retain independent physical units.",false,(model:Handle<Model>,initial:Input,start:Input,h:Input,steps:usize,paths:usize,seed:u64,scheme:String,terms:usize,maxeval:usize)->Many<Dynamic>,|ctx,model,initial,start,h,steps,paths,seed,scheme,terms,maxeval|{let model=model.with(ctx,|m|Ok(m.clone()))?;let frame=&model.drift.frame;let y=frame.state(&initial)?;let t=frame.time(&start)?;let h=frame.time(&h)?;if h<=0.0||paths==0||paths>4096||paths.checked_mul(steps).and_then(|n|n.checked_mul(y.len().max(model.m*model.m))).is_none_or(|n|n>262144){return Err(api::error("SDE ensemble step/path/allocation work budget exceeded"));}let scheme=ensemble_scheme(&scheme,terms)?;let budget=Budget::new(maxeval)?;let drift=|t,y:&[f64],out:&mut[f64]|model.drift.rate(&budget,t,y,out);let diff=|t,y:&[f64],out:&mut[f64]|model.diffusion(&budget,t,y,out);let result=s::sde_ensemble(drift,diff,model.m,(t,h,steps),&y,paths,seed,scheme);budget.finish()?;Ok(Many(result.map_err(|e|api::error(format!("{e:?}")))?.iter().map(|y|frame.restore(y).map(Dynamic)).collect::<Result<Vec<_>,_>>()?))});
#[derive(Clone)]
struct Converted(Model);
impl NativeType for Converted {
    const KIND: &'static str = "rssn/native/analysis/sde/ito-drift/v1";
    const MAX_BYTES: usize = 8 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        self.0.validate()
    }
}
object_binding!(ito_new,"num_analysis_sde_stratonovich_to_ito_drift","crate::kernels::sde::stratonovich_to_ito_drift","Owned native Stratonovich-to-Ito drift recipe. Evaluation calls the actual native derivative-based conversion over the immutable physical reference frame, retaining no source callback lifetime.",false,(model:Handle<Model>)->NewObject<Converted>,|ctx,model|model.with(ctx,|m|{let budget=Budget::new(1)?;let drift=|t,y:&[f64],out:&mut[f64]|m.drift.rate(&budget,t,y,out);let diff=|t,y:&[f64],out:&mut[f64]|m.diffusion(&budget,t,y,out);let _native=s::stratonovich_to_ito_drift(&drift,&diff,m.m);Ok(NewObject::new(Converted(m.clone()),m.drift.frame.dimensions()))}));
object_binding!(ito_eval,"num_analysis_sde_ito_drift_eval","","Evaluate the actual native converted Ito drift. Each result component has its state unit divided by the time unit, with a checked shared evaluation budget.",false,(model:Handle<Converted>,time:Input,state:Input,maxeval:usize)->Dynamic,|ctx,model,time,state,maxeval|model.with(ctx,|c|{let m=&c.0;let frame=&m.drift.frame;let t=frame.time(&time)?;let y=frame.state(&state)?;let budget=Budget::new(maxeval)?;let drift=|t,y:&[f64],out:&mut[f64]|m.drift.rate(&budget,t,y,out);let diff=|t,y:&[f64],out:&mut[f64]|m.diffusion(&budget,t,y,out);let native=s::stratonovich_to_ito_drift(&drift,&diff,m.m);let mut out=vec![0.0;y.len()];native(t,&y,&mut out);budget.finish()?;let data=out.iter().zip(&frame.state_scales).zip(&frame.state_dimensions).map(|((v,s),d)|api::quantity(api::scale_product(*v,&[*s],&[frame.time_scale])?,DimExpr::Mul(vec![d.clone(),DimExpr::Inverse(Box::new(frame.time_dimension.clone()))]))).collect::<Result<Vec<_>,_>>()?;Ok(Dynamic(Value::Tensor(crate::kernel::value::TensorValue{shape:vec![data.len()],elem:api::dynamic_type(),data})))}));
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    extended::install(i)?;
    stochastic::install(i)?;
    api::install(
        i,
        &[
            &brownian::BINDING,
            &levy::BINDING,
            &noise_data::BINDING,
            &grid::BINDING,
            &ensemble::BINDING,
            &ito_new::BINDING,
            &ito_eval::BINDING,
            &model_new::BINDING,
            &sample::BINDING,
            &aggregate::BINDING,
            &euler_maruyama::BINDING,
            &heun_stratonovich::BINDING,
            &milstein::BINDING,
            &milstein_stratonovich::BINDING,
            &srk_platen::BINDING,
            &milstein_noncommutative::BINDING,
            &milstein_noncommutative_stratonovich::BINDING,
            &srk_noncommutative::BINDING,
            &srk_noncommutative_stratonovich::BINDING,
            &sde_adaptive::BINDING,
            &sde_adaptive_with::BINDING,
            &sde_adaptive_stratonovich_with::BINDING,
            &sde_adaptive_noncommutative_with::BINDING,
            &sde_adaptive_noncommutative_stratonovich_with::BINDING,
            &data::BINDING,
            &ensemble_seed::BINDING,
        ],
        "numerical_vector",
    )
}
