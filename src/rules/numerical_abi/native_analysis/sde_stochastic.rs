//! Native simulation stochastic options, streams, paths and method metadata.
use super::Budget;
use super::CallCtx;
use super::DimExpr;
use super::Dynamic;
use super::Handle;
use super::Input;
use super::Installer;
use super::Many;
use super::Model;
use super::NativeType;
use super::NewObject;
use super::OdeFrame;
use super::OpError;
use super::Output;
use super::RuleError;
use super::TapeSnapshot;
use super::Value;
use super::api;
use super::extended;
use super::half;
use super::noise_limit;
use super::object_binding;
use super::parameters;
use super::tape_snapshot;
use crate::rules::numerical_abi::native_object_api::Maybe;
use crate::sim::stochastic as st;
struct Stream(st::Stream);
impl NativeType for Stream {
    const KIND: &'static str = "rssn/native/analysis/stochastic/stream/v1";
    const MAX_BYTES: usize = 64;

    fn validate(&self) -> Result<(), OpError> {
        Ok(())
    }
}
#[derive(Clone)]
struct Options {
    data: st::Options,
    frame: OdeFrame,
}
impl NativeType for Options {
    const KIND: &'static str = "rssn/native/analysis/stochastic/options/v1";
    const MAX_BYTES: usize = 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        let n = self.frame.state_scales.len();
        if n == 0 || n > 128 || self.data.terms > 128 {
            return Err(api::error(
                "stochastic option shape/Levy terms budget exceeded",
            ));
        }
        let steps = match &self.data.schedule {
            | st::Schedule::Fixed { dt, steps } => {
                if *dt <= 0.0 || !dt.is_finite() {
                    return Err(api::error("stochastic step must be positive finite"));
                }
                *steps
            },
            | st::Schedule::Grid { times } => {
                if times.len() < 2
                    || times.iter().any(|t| !t.is_finite())
                    || times.windows(2).any(|r| r[1] <= r[0])
                {
                    return Err(api::error("native stochastic grid must increase"));
                }
                times.len()
            },
            | st::Schedule::Adaptive {
                t1,
                rtol,
                atol,
                max_steps,
                cells,
            } => {
                if *t1 <= self.data.t0
                    || !t1.is_finite()
                    || *rtol <= 0.0
                    || !rtol.is_finite()
                    || *atol < 0.0
                    || !atol.is_finite()
                    || *cells == 0
                    || *cells > 16384
                {
                    return Err(api::error(
                        "adaptive stochastic tolerance/span/cell bounds invalid",
                    ));
                }
                *max_steps
            },
        };
        if steps == 0
            || steps
                .checked_mul(n.max(self.data.structure.noise(n).pow(2)))
                .is_none_or(|k| k > 262000)
        {
            return Err(api::error(
                "stochastic path allocation/work budget exceeded",
            ));
        }
        Ok(())
    }
}
struct SimPath {
    data: st::Path,
    frame: OdeFrame,
}
impl NativeType for SimPath {
    const KIND: &'static str = "rssn/native/analysis/stochastic/path/v1";
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
            return Err(api::error("invalid native stochastic path shape/storage"));
        }
        Ok(())
    }
}
fn method(name: &str) -> Result<st::Method, OpError> {
    st::Method::from_name(name).ok_or_else(|| api::error("unknown native stochastic method name"))
}
fn structure(
    diagonal: bool,
    commutative: bool,
    m: usize,
    n: usize,
) -> Result<st::Structure, OpError> {
    if m == 0 || m > 16 || diagonal && m != n {
        return Err(api::error(
            "diagonal structure needs one noise channel per state; noise dimension1..16",
        ));
    }
    Ok(if diagonal {
        st::Structure::Diagonal
    } else {
        st::Structure::General {
            noise: m,
            commutative,
        }
    })
}
fn same_frame(
    a: &OdeFrame,
    b: &OdeFrame,
) -> Result<(), OpError> {
    if a.state_scales != b.state_scales
        || !crate::kernels::exact_eq(a.time_scale, b.time_scale)
        || a.state_dimensions.len() != b.state_dimensions.len()
    {
        return Err(api::error(
            "native options and model reference frames differ",
        ));
    }
    api::same_dimension(&a.time_dimension, &b.time_dimension)?;
    for (a, b) in a.state_dimensions.iter().zip(&b.state_dimensions) {
        api::same_dimension(a, b)?;
    }
    Ok(())
}
fn rates(
    model: &Model,
    budget: &Budget,
    t: f64,
    y: &[f64],
    diffusion: bool,
    structure: st::Structure,
) -> Result<Vec<f64>, st::SdeError> {
    let n = y.len();
    let mut out = vec![0.0; if diffusion { n * model.m } else { n }];
    if diffusion {
        model.diffusion(budget, t, y, &mut out);
    } else {
        model.drift.rate(budget, t, y, &mut out);
    }
    if let Err(e) = budget.finish() {
        return Err(st::SdeError::StepLimit(format!("{e:?}")));
    }
    if diffusion && structure == st::Structure::Diagonal {
        if out
            .iter()
            .enumerate()
            .any(|(k, v)| k / n != k % n && *v != 0.0)
        {
            budget.fail(api::error(
                "diffusion declared diagonal has nonzero off-diagonal coefficient",
            ));
            return Err(st::SdeError::Shape("nonzero off-diagonal diffusion".into()));
        }
        return Ok((0..n).map(|i| out[i * n + i]).collect());
    }
    Ok(out)
}
object_binding!(stream_new,"num_analysis_stochastic_stream_new","crate::sim::stochastic::Stream::new","Actual native stochastic method in explicit physical reference coordinates. Session/kind, storage, shape and evaluation limits are checked.",false,(seed:u64)->NewObject<Stream>,|_ctx,seed|Ok(NewObject::normalized(Stream(st::Stream::new(seed)))));
object_binding!(next_unit,"num_analysis_stochastic_next_unit","crate::sim::stochastic::Stream::next_unit","Actual native stochastic method in explicit physical reference coordinates. Session/kind, storage, shape and evaluation limits are checked.",false,(stream:Handle<Stream>)->f64,|ctx,stream|stream.with_mut(ctx,|s|Ok(s.0.next_unit())));
object_binding!(next_normal,"num_analysis_stochastic_next_normal","crate::sim::stochastic::Stream::next_normal","Actual native stochastic method in explicit physical reference coordinates. Session/kind, storage, shape and evaluation limits are checked.",false,(stream:Handle<Stream>)->f64,|ctx,stream|stream.with_mut(ctx,|s|Ok(s.0.next_normal())));
object_binding!(method_from_name,"num_analysis_stochastic_method_from_name","crate::sim::stochastic::Method::from_name","Actual native stochastic method in explicit physical reference coordinates. Session/kind, storage, shape and evaluation limits are checked.",false,(name:String)->Maybe<String>,|_ctx,name|Ok(Maybe(st::Method::from_name(&name).map(|m|m.name().into()))));
object_binding!(method_name,"num_analysis_stochastic_method_name","crate::sim::stochastic::Method::name","Actual native stochastic method in explicit physical reference coordinates. Session/kind, storage, shape and evaluation limits are checked.",false,(name:String)->String,|_ctx,name|Ok(method(&name)?.name().into()));
object_binding!(method_calculus,"num_analysis_stochastic_method_calculus","crate::sim::stochastic::Method::calculus","Actual native stochastic method in explicit physical reference coordinates. Session/kind, storage, shape and evaluation limits are checked.",false,(name:String)->String,|_ctx,name|Ok(format!("{:?}",method(&name)?.calculus())));
object_binding!(structure_noise,"num_analysis_stochastic_structure_noise","crate::sim::stochastic::Structure::noise","Actual native stochastic method in explicit physical reference coordinates. Session/kind, storage, shape and evaluation limits are checked.",false,(diagonal:bool,commutative:bool,m:usize,n:usize)->usize,|_ctx,diagonal,commutative,m,n|{if n==0||n>128{return Err(api::error("state size1..128required"));}Ok(structure(diagonal,commutative,m,n)?.noise(n))});
object_binding!(structure_coefficients,"num_analysis_stochastic_structure_coefficients","crate::sim::stochastic::Structure::coefficients","Actual native stochastic method in explicit physical reference coordinates. Session/kind, storage, shape and evaluation limits are checked.",false,(diagonal:bool,commutative:bool,m:usize,n:usize)->usize,|_ctx,diagonal,commutative,m,n|{if n==0||n>128{return Err(api::error("state size1..128required"));}Ok(structure(diagonal,commutative,m,n)?.coefficients(n))});
object_binding!(structure_is_commutative,"num_analysis_stochastic_structure_is_commutative","crate::sim::stochastic::Structure::is_commutative","Actual native stochastic method in explicit physical reference coordinates. Session/kind, storage, shape and evaluation limits are checked.",false,(diagonal:bool,commutative:bool,m:usize,n:usize)->bool,|_ctx,diagonal,commutative,m,n|{if n==0||n>128{return Err(api::error("state size1..128required"));}Ok(structure(diagonal,commutative,m,n)?.is_commutative())});
object_binding!(sde_order,"num_analysis_stochastic_sde_order","crate::sim::stochastic::sde_order","Actual native stochastic method in explicit physical reference coordinates. Session/kind, storage, shape and evaluation limits are checked.",false,(name:String)->Dynamic,|_ctx,name|{let o=st::sde_order(method(&name)?);Ok(Dynamic(api::record(vec![("strong",Output::encode(o.strong)?),("weak",Output::encode(o.weak)?)])))});
object_binding!(sde_order_for,"num_analysis_stochastic_sde_order_for","crate::sim::stochastic::sde_order_for","Actual native stochastic method in explicit physical reference coordinates. Session/kind, storage, shape and evaluation limits are checked.",false,(name:String,diagonal:bool,commutative:bool,m:usize,n:usize)->Dynamic,|_ctx,name,diagonal,commutative,m,n|{let o=st::sde_order_for(method(&name)?,structure(diagonal,commutative,m,n)?);Ok(Dynamic(api::record(vec![("strong",Output::encode(o.strong)?),("weak",Output::encode(o.weak)?)])))});
object_binding!(default_terms,"num_analysis_stochastic_default_terms","crate::sim::stochastic::default_terms","Actual native stochastic method in explicit physical reference coordinates. Session/kind, storage, shape and evaluation limits are checked.",false,(normalized_step:f64)->usize,|_ctx,normalized_step|{if normalized_step<=0.0||!normalized_step.is_finite(){return Err(api::error("default_terms needs positive normalized step"));}Ok(st::default_terms(normalized_step))});
object_binding!(ensemble_seed,"num_analysis_stochastic_ensemble_seed","crate::sim::stochastic::ensemble_seed","Actual native stochastic method in explicit physical reference coordinates. Session/kind, storage, shape and evaluation limits are checked.",false,(seed:u64,index:usize)->u64,|_ctx,seed,index|Ok(st::ensemble_seed(seed,index)));
object_binding!(options_fixed,"num_analysis_stochastic_options_fixed","crate::sim::stochastic::Options::fixed","Actual native stochastic method in explicit physical reference coordinates. Session/kind, storage, shape and evaluation limits are checked.",false,(model:Handle<Model>,name:String,diagonal:bool,commutative:bool,seed:u64,start:Input,h:Input,steps:usize)->NewObject<Options>,|ctx,model,name,diagonal,commutative,seed,start,h,steps|{let model=model.with(ctx,|m|Ok(m.clone()))?;let frame=model.drift.frame.clone();let method=method(&name)?;let structure=structure(diagonal,commutative,model.m,frame.state_scales.len())?;let t0=frame.time(&start)?;let data=st::Options::fixed(method,structure,seed,t0,frame.time(&h)?,steps);let dimensions=frame.dimensions();let options=Options{data,frame};options.validate()?;Ok(NewObject::new(options,dimensions))});
object_binding!(options_fixed_span,"num_analysis_stochastic_options_fixed_span","crate::sim::stochastic::Options::fixed_span","Actual native stochastic method in explicit physical reference coordinates. Session/kind, storage, shape and evaluation limits are checked.",false,(model:Handle<Model>,name:String,diagonal:bool,commutative:bool,seed:u64,start:Input,end:Input,h:Input)->NewObject<Options>,|ctx,model,name,diagonal,commutative,seed,start,end,h|{let model=model.with(ctx,|m|Ok(m.clone()))?;let frame=model.drift.frame.clone();let method=method(&name)?;let structure=structure(diagonal,commutative,model.m,frame.state_scales.len())?;let t0=frame.time(&start)?;let data={let t1=frame.time(&end)?;let h=frame.time(&h)?;if h<=0.0||t1<=t0||(t1-t0)/h>262000.0{return Err(api::error("fixed-span step budget exceeded"));}st::Options::fixed_span(method,structure,seed,t0,t1,h)};let dimensions=frame.dimensions();let options=Options{data,frame};options.validate()?;Ok(NewObject::new(options,dimensions))});
object_binding!(options_adaptive,"num_analysis_stochastic_options_adaptive","crate::sim::stochastic::Options::adaptive","Actual native stochastic method in explicit physical reference coordinates. Session/kind, storage, shape and evaluation limits are checked.",false,(model:Handle<Model>,name:String,diagonal:bool,commutative:bool,seed:u64,start:Input,end:Input,rtol:f64,atol:f64,max_steps:usize,cells:usize)->NewObject<Options>,|ctx,model,name,diagonal,commutative,seed,start,end,rtol,atol,max_steps,cells|{let model=model.with(ctx,|m|Ok(m.clone()))?;let frame=model.drift.frame.clone();let method=method(&name)?;let structure=structure(diagonal,commutative,model.m,frame.state_scales.len())?;let t0=frame.time(&start)?;let data=st::Options::adaptive(method,structure,seed,t0,frame.time(&end)?,rtol,atol,max_steps,cells);let dimensions=frame.dimensions();let options=Options{data,frame};options.validate()?;Ok(NewObject::new(options,dimensions))});
object_binding!(options_grid,"num_analysis_stochastic_options_grid","crate::sim::stochastic::Options::grid","Actual native stochastic method in explicit physical reference coordinates. Session/kind, storage, shape and evaluation limits are checked.",false,(model:Handle<Model>,name:String,diagonal:bool,commutative:bool,seed:u64,start:Input,times:Many<Input>)->NewObject<Options>,|ctx,model,name,diagonal,commutative,seed,start,times|{let model=model.with(ctx,|m|Ok(m.clone()))?;let frame=model.drift.frame.clone();let method=method(&name)?;let structure=structure(diagonal,commutative,model.m,frame.state_scales.len())?;let t0=frame.time(&start)?;let data={let times=times.0.iter().map(|t|frame.time(t)).collect::<Result<Vec<_>,_>>()?;if times.first()!=Some(&t0){return Err(api::error("native grid first time must equal start"));}st::Options::grid(method,structure,seed,times)};let dimensions=frame.dimensions();let options=Options{data,frame};options.validate()?;Ok(NewObject::new(options,dimensions))});
object_binding!(options_with_terms,"num_analysis_stochastic_options_with_terms","crate::sim::stochastic::Options::with_terms","Actual native stochastic method in explicit physical reference coordinates. Session/kind, storage, shape and evaluation limits are checked.",false,(options:Handle<Options>,terms:usize)->NewObject<Options>,|ctx,options,terms|options.with(ctx,|o|{if terms==0||terms>128{return Err(api::error("explicit Levy terms must be1..128"));}let data=o.data.clone().with_terms(terms);Ok(NewObject::new(Options{data,frame:o.frame.clone()},o.frame.dimensions()))}));
object_binding!(integrate,"num_analysis_stochastic_integrate","crate::sim::stochastic::integrate","Actual native stochastic method in explicit physical reference coordinates. Session/kind, storage, shape and evaluation limits are checked.",false,(model:Handle<Model>,options:Handle<Options>,initial:Input,maxeval:usize)->NewObject<SimPath>,|ctx,model,options,initial,maxeval|{let model=model.with(ctx,|m|Ok(m.clone()))?;let options=options.with(ctx,|o|Ok(o.clone()))?;same_frame(&model.drift.frame,&options.frame)?;let initial=options.frame.state(&initial)?;let budget=Budget::new(maxeval)?;let drift=|t,y:&[f64]|rates(&model,&budget,t,y,false,options.data.structure);let diff=|t,y:&[f64]|rates(&model,&budget,t,y,true,options.data.structure);let result=st::integrate(drift,diff,&initial,&options.data);budget.finish()?;Ok(NewObject::new(SimPath{data:result.map_err(api::error)?,frame:options.frame.clone()},options.frame.dimensions()))});
object_binding!(integrate_drift_implicit,"num_analysis_stochastic_integrate_drift_implicit","crate::sim::stochastic::integrate_drift_implicit","Actual native stochastic method in explicit physical reference coordinates. Session/kind, storage, shape and evaluation limits are checked.",false,(model:Handle<Model>,options:Handle<Options>,initial:Input,maxeval:usize,scheme:String,theta:f64)->NewObject<SimPath>,|ctx,model,options,initial,maxeval,scheme,theta|{let model=model.with(ctx,|m|Ok(m.clone()))?;let options=options.with(ctx,|o|Ok(o.clone()))?;same_frame(&model.drift.frame,&options.frame)?;let initial=options.frame.state(&initial)?;let budget=Budget::new(maxeval)?;let drift=|t,y:&[f64]|rates(&model,&budget,t,y,false,options.data.structure);let diff=|t,y:&[f64]|rates(&model,&budget,t,y,true,options.data.structure);let result=st::integrate_drift_implicit(drift,diff,&initial,super::extended::scheme(&scheme,theta)?,&options.data);budget.finish()?;Ok(NewObject::new(SimPath{data:result.map_err(api::error)?,frame:options.frame.clone()},options.frame.dimensions()))});
object_binding!(ensemble,"num_analysis_stochastic_ensemble","crate::sim::stochastic::ensemble","Actual native stochastic method in explicit physical reference coordinates. Session/kind, storage, shape and evaluation limits are checked.",false,(model:Handle<Model>,options:Handle<Options>,initial:Input,paths:usize,maxeval:usize)->Many<NewObject<SimPath>>,|ctx,model,options,initial,paths,maxeval|{let model=model.with(ctx,|m|Ok(m.clone()))?;let options=options.with(ctx,|o|Ok(o.clone()))?;same_frame(&model.drift.frame,&options.frame)?;let initial=options.frame.state(&initial)?;if paths==0||paths>7{return Err(api::error("retained ensemble path count must be1..7within64MiBbank"));}let budget=Budget::new(maxeval)?;let drift=|t,y:&[f64]|rates(&model,&budget,t,y,false,options.data.structure);let diff=|t,y:&[f64]|rates(&model,&budget,t,y,true,options.data.structure);let result=st::ensemble(drift,diff,&initial,&options.data,paths);budget.finish()?;Ok(Many(result.map_err(api::error)?.into_iter().map(|data|NewObject::new(SimPath{data,frame:options.frame.clone()},options.frame.dimensions())).collect()))});
object_binding!(path_last,"num_analysis_stochastic_path_last","crate::sim::stochastic::Path::last","Actual native stochastic method in explicit physical reference coordinates. Session/kind, storage, shape and evaluation limits are checked.",false,(path:Handle<SimPath>)->Maybe<Dynamic>,|ctx,path|path.with(ctx,|p|Ok(Maybe(p.data.last().map(|y|p.frame.restore(y).map(Dynamic)).transpose()?))));
object_binding!(path_eval,"num_analysis_stochastic_path_eval","crate::sim::stochastic::Path::eval","Actual native stochastic method in explicit physical reference coordinates. Session/kind, storage, shape and evaluation limits are checked.",false,(path:Handle<SimPath>,time:Input)->Maybe<Dynamic>,|ctx,path,time|path.with(ctx,|p|Ok(Maybe(p.data.eval(p.frame.time(&time)?).map(|y|p.frame.restore(&y).map(Dynamic)).transpose()?))));
object_binding!(path_resample,"num_analysis_stochastic_path_resample","crate::sim::stochastic::Path::resample","Actual native stochastic method in explicit physical reference coordinates. Session/kind, storage, shape and evaluation limits are checked.",false,(path:Handle<SimPath>,times:Many<Input>)->Many<Dynamic>,|ctx,path,times|path.with(ctx,|p|{if times.0.len().checked_mul(p.frame.state_scales.len()).is_none_or(|n|n>262144){return Err(api::error("stochastic resample storage budget exceeded"));}let ts=times.0.iter().map(|t|p.frame.time(t)).collect::<Result<Vec<_>,_>>()?;Ok(Many(p.data.resample(&ts).map_err(api::error)?.iter().map(|y|p.frame.restore(y).map(Dynamic)).collect::<Result<_,_>>()?))}));
object_binding!(wiener,"num_analysis_stochastic_wiener","crate::sim::stochastic::wiener","Actual native stochastic method in explicit physical reference coordinates. Session/kind, storage, shape and evaluation limits are checked.",false,(seed:u64,m:usize,steps:usize,h:Input)->Dynamic,|_ctx,seed,m,steps,h|{noise_limit(m,steps,1)?;let dimension=api::dimension(&h.0)?;let h=crate::rules::numerical_abi::numeric(&h.0,0)?;if h<=0.0{return Err(api::error("Wiener step must be positive"));}let values=st::wiener(seed,m,steps,h);let data=values.into_iter().flatten().map(|v|api::scale_product(v,&[1.0],&[]).map(Value::Float)).collect::<Result<_,_>>()?;Ok(Dynamic(Value::Tensor(crate::kernel::value::TensorValue{shape:vec![steps,m],elem:crate::kernel::signature::TypeExpr::Dim(Box::new(crate::kernel::signature::TypeExpr::Float),half(dimension)),data})))});
object_binding!(levy_increments,"num_analysis_stochastic_levy_increments","crate::sim::stochastic::levy_increments","Actual native stochastic method in explicit physical reference coordinates. Session/kind, storage, shape and evaluation limits are checked.",false,(seed:u64,m:usize,steps:usize,h:Input,alpha:f64,beta:f64,increment_reference:Input)->Dynamic,|_ctx,seed,m,steps,h,alpha,beta,increment_reference|{noise_limit(m,steps,1)?;let(reference,dimension)=api::reference(&increment_reference.0)?;api::normalized(&h.0)?;let h=crate::rules::numerical_abi::numeric(&h.0,0)?;if h<=0.0{return Err(api::error("Levy increments require positive explicitly normalized time step"));}let values=st::levy_increments(seed,m,steps,h,alpha,beta).map_err(api::error)?;let data=values.into_iter().flatten().map(|v|api::scale_product(v,&[reference],&[]).map(Value::Float)).collect::<Result<_,_>>()?;Ok(Dynamic(Value::Tensor(crate::kernel::value::TensorValue{shape:vec![steps,m],elem:crate::kernel::signature::TypeExpr::Dim(Box::new(crate::kernel::signature::TypeExpr::Float),dimension),data})))});
object_binding!(sample_poisson,"num_analysis_stochastic_sample_poisson","crate::sim::stochastic::sample_poisson","Actual native stochastic method in explicit physical reference coordinates. Session/kind, storage, shape and evaluation limits are checked.",false,(stream:Handle<Stream>,mean:f64)->usize,|ctx,stream,mean|{if !mean.is_finite()||!(0.0..=100000.0).contains(&mean){return Err(api::error("Poisson mean must be0..100000"));}stream.with_mut(ctx,|s|Ok(st::sample_poisson(&mut s.0,mean)))});
object_binding!(sample_alpha_stable,"num_analysis_stochastic_sample_alpha_stable","crate::sim::stochastic::sample_alpha_stable","Actual native stochastic method in explicit physical reference coordinates. Session/kind, storage, shape and evaluation limits are checked.",false,(stream:Handle<Stream>,alpha:f64,beta:f64,scale:Input,location:Input)->Dynamic,|ctx,stream,alpha,beta,scale,location|{let(scale,dim)=api::reference(&scale.0)?;api::same_dimension(&dim,&api::dimension(&location.0)?)?;let location=api::scale_product(crate::rules::numerical_abi::numeric(&location.0,0)?,&[],&[scale])?;let result=stream.with_mut(ctx,|s|st::sample_alpha_stable(&mut s.0,alpha,beta,1.0,location).map_err(api::error))?;Ok(Dynamic(api::quantity(api::scale_product(result,&[scale],&[])?,dim)?))});

struct Error(st::SdeError);
impl NativeType for Error {
    const KIND: &'static str = "rssn/native/analysis/stochastic/error/v1";
    const MAX_BYTES: usize = 4096;

    fn validate(&self) -> Result<(), OpError> {
        if self.0.to_string().len() > 4096 {
            return Err(api::error("stochastic error text exceeds4096bytes"));
        }
        Ok(())
    }
}
object_binding!(error_new,"num_analysis_stochastic_error_new","","Native stochastic error builder with exact Invalid/Shape/NonFinite/StepLimit/Unsupported variants; unsupported capability is an explicit closed known capability name.",false,(kind:String,message:String,capability:String)->NewObject<Error>,|_ctx,kind,message,capability|{if message.len()>2048{return Err(api::error("error text budget exceeded"));}let error=match kind.as_str(){"Invalid"=>st::SdeError::Invalid(message),"Shape"=>st::SdeError::Shape(message),"NonFinite"=>st::SdeError::NonFinite(message),"StepLimit"=>st::SdeError::StepLimit(message),"Unsupported"=>st::SdeError::Unsupported{capability:match capability.as_str(){"core.sde"=>"core.sde","sde.events"=>"sde.events","sde.noncommutative"=>"sde.noncommutative","sde.interpretation"=>"sde.interpretation",_=>return Err(api::error("unknown native capability name"))},detail:message},_=>return Err(api::error("unknown stochastic error variant"))};Ok(NewObject::normalized(Error(error)))});
object_binding!(error_code,"num_analysis_stochastic_error_code","crate::sim::stochastic::SdeError::code","Actual native stable diagnostic code.",false,(error:Handle<Error>)->String,|ctx,error|error.with(ctx,|e|Ok(e.0.code().into())));
object_binding!(error_capability,"num_analysis_stochastic_error_capability","crate::sim::stochastic::SdeError::capability","Actual native optional capability diagnostic.",false,(error:Handle<Error>)->Maybe<String>,|ctx,error|error.with(ctx,|e|Ok(Maybe(e.0.capability().map(str::to_owned)))));
object_binding!(events,"num_analysis_stochastic_events","crate::sim::stochastic::events","Actual native sign-change event detection on the native piecewise-linear stochastic path; this does not claim to detect unobserved between-node process crossings. Programs receive coherent physical time/state components.",false,(path:Handle<SimPath>,specs:Many<Input>,maxeval:usize)->Dynamic,|ctx,path,specs,maxeval|path.with(ctx,|p|{use crate::rules::numerical_abi::native_analysis_ode::{events,event_value};let programs=events(ctx,&specs,&p.frame,maxeval)?;let specs=programs.iter().map(|p|st::EventSpec{terminal:p.terminal,direction:p.direction}).collect::<Vec<_>>();let budget=Budget::new(maxeval)?;let mut event=|index,t,y:&[f64]|event_value(&programs[index],&p.frame,&budget,t,y);let (truncated,hits)=st::events(&p.data,&specs,&mut event);budget.finish()?;for p in &programs{p.program.finish()?;}let times=Value::Tensor(crate::kernel::value::TensorValue{shape:vec![truncated.t.len()],elem:crate::kernel::signature::TypeExpr::Dim(Box::new(crate::kernel::signature::TypeExpr::Float),p.frame.time_dimension.clone()),data:truncated.t.iter().map(|v|api::scale_product(*v,&[p.frame.time_scale],&[]).map(Value::Float)).collect::<Result<_,_>>()?});let hits=Value::Tensor(crate::kernel::value::TensorValue{shape:vec![hits.len()],elem:api::dynamic_type(),data:hits.iter().map(|h|Ok(api::record(vec![("index",Output::encode(h.index)?),("time",api::quantity(api::scale_product(h.t,&[p.frame.time_scale],&[])?,p.frame.time_dimension.clone())?),("state",p.frame.restore(&h.y)?)]))).collect::<Result<_,OpError>>()?});let states=Value::Tensor(crate::kernel::value::TensorValue{shape:vec![truncated.y.len()],elem:api::dynamic_type(),data:truncated.y.iter().map(|y|p.frame.restore(y)).collect::<Result<_,_>>()?});Ok(Dynamic(api::record(vec![("times",times),("states",states),("hits",hits)])))}));
#[derive(Clone)]
struct Distribution {
    data: st::JumpDistribution,
    scale: f64,
    dimension: DimExpr,
    n: usize,
}
impl NativeType for Distribution {
    const KIND: &'static str = "rssn/native/analysis/stochastic/jump-distribution/v1";
    const MAX_BYTES: usize = 4096;

    fn validate(&self) -> Result<(), OpError> {
        if self.n == 0 || self.n > 16 || !self.scale.is_finite() || self.scale <= 0.0 {
            return Err(api::error("invalid jump mark dimension/reference"));
        }
        Ok(())
    }
}
object_binding!(jump_distribution,"num_analysis_stochastic_jump_distribution","","All native JumpDistribution variants over an explicit positive physical mark reference. data1/data2 contain normalized vectors, controls contain normalized distribution parameters. Exact native validation is called before retention.",false,(kind:String,data1:Vec<f64>,data2:Vec<f64>,controls:Vec<f64>,mark_reference:Input,n:usize)->NewObject<Distribution>,|_ctx,kind,data1,data2,controls,mark_reference,n|{if n==0||n>16||data1.len()>16||data2.len()>16||controls.len()>4{return Err(api::error("jump distribution shape budget exceeded"));}let(scale,dimension)=api::reference(&mark_reference.0)?;let c=|i|controls.get(i).copied().ok_or_else(||api::error("missing jump distribution control"));let data=match kind.as_str(){"Constant"=>{if data1.len()!=n{return Err(api::error("constant mark shape mismatch"));}st::JumpDistribution::Constant(data1)},"Normal"=>{if data1.len()!=n||data2.len()!=n{return Err(api::error("normal mark shape mismatch"));}st::JumpDistribution::Normal{mean:data1,std:data2}},"DoubleExponential"=>st::JumpDistribution::DoubleExponential{p_up:c(0)?,eta_up:c(1)?,eta_down:c(2)?,dim:n},"Uniform"=>{if data1.len()!=n||data2.len()!=n{return Err(api::error("uniform mark shape mismatch"));}st::JumpDistribution::Uniform{low:data1,high:data2}},"AlphaStable"=>st::JumpDistribution::AlphaStable{alpha:c(0)?,beta:c(1)?,scale:c(2)?,location:c(3)?,dim:n},_=>return Err(api::error("unknown native jump distribution variant"))};let mut probe=st::Stream::new(0);let _=data.sample(&mut probe).map_err(api::error)?;Ok(NewObject::new(Distribution{data,scale,dimension:dimension.clone(),n},vec![dimension]))});
fn marks(
    d: &Distribution,
    values: Vec<f64>,
) -> Result<Dynamic, OpError> {
    let data = values
        .into_iter()
        .map(|v| api::scale_product(v, &[d.scale], &[]).map(Value::Float))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(Dynamic(Value::Tensor(crate::kernel::value::TensorValue {
        shape: vec![data.len()],
        elem: crate::kernel::signature::TypeExpr::Dim(
            Box::new(crate::kernel::signature::TypeExpr::Float),
            d.dimension.clone(),
        ),
        data,
    })))
}
object_binding!(jump_sample,"num_analysis_stochastic_jump_sample","crate::sim::stochastic::JumpDistribution::sample","Actual native jump mark sampling on the retained native stream, returning physical mark units.",false,(distribution:Handle<Distribution>,stream:Handle<Stream>)->Dynamic,|ctx,distribution,stream|distribution.with(ctx,|d|{let v=stream.with_mut(ctx,|s|d.data.sample(&mut s.0).map_err(api::error))?;marks(d,v)}));
object_binding!(jump_mean,"num_analysis_stochastic_jump_mean","crate::sim::stochastic::JumpDistribution::mean_mark","Actual native optional analytical mean mark; infinite or undefined stable means remain None.",false,(distribution:Handle<Distribution>)->Maybe<Dynamic>,|ctx,distribution|distribution.with(ctx,|d|Ok(Maybe(d.data.mean_mark().map(|v|marks(d,v)).transpose()?))));
object_binding!(mlmc,"num_analysis_stochastic_mlmc","crate::sim::stochastic::mlmc","Actual native MLMC with coupled fine/coarse Brownian paths. Payoff compiled program receives physical state components; target RMSE and estimate/standard-error have payoff units and variance has their square. Explicit hierarchy/sample/work caps bound all native levels.",false,(model:Handle<Model>,options:Handle<Options>,initial:Input,end:Input,payoff:Input,target_rmse:Input,hierarchy:Vec<usize>,maxeval:usize)->Dynamic,|ctx,model,options,initial,end,payoff,target_rmse,hierarchy,maxeval|{let model=model.with(ctx,|m|Ok(m.clone()))?;let options=options.with(ctx,|o|Ok(o.clone()))?;same_frame(&model.drift.frame,&options.frame)?;let frame=&options.frame;let initial=frame.state(&initial)?;if hierarchy.len()!=5{return Err(api::error("MLMC hierarchy requires base_steps,min_level,max_level,warmup,max_samples"));}let(base_steps,min_level,max_level,warmup_samples,max_samples_per_level)=(hierarchy[0],hierarchy[1],hierarchy[2],hierarchy[3],hierarchy[4]);if base_steps==0||max_level>8||min_level>max_level||warmup_samples<2||max_samples_per_level<warmup_samples||max_samples_per_level>4096||base_steps.checked_mul(1usize<<max_level).and_then(|n|n.checked_mul(max_samples_per_level)).and_then(|n|n.checked_mul(max_level+1)).and_then(|n|n.checked_mul(initial.len().max(model.m*model.m))).is_none_or(|n|n>2000000){return Err(api::error("MLMC hierarchy worst-case work budget exceeded"));}let payoff=crate::rules::numerical_abi::native_backend::compiled_snapshot(ctx,&payoff.0,maxeval)?;if payoff.input_dimensions.len()!=initial.len(){return Err(api::error("MLMC payoff input shape mismatch"));}for(a,b)in payoff.input_dimensions.iter().zip(&frame.state_dimensions){api::same_dimension(a,b)?;}api::same_dimension(&payoff.output_dimension,&api::dimension(&target_rmse.0)?)?;let target_rmse=crate::rules::numerical_abi::numeric(&target_rmse.0,0)?;let cfg=st::MlmcOptions{method:options.data.method,structure:options.data.structure,seed:options.data.seed,t0:options.data.t0,t1:frame.time(&end)?,base_steps,min_level,max_level,warmup_samples,max_samples_per_level,target_rmse};let budget=Budget::new(maxeval)?;let drift=|t,y:&[f64]|rates(&model,&budget,t,y,false,options.data.structure);let diff=|t,y:&[f64]|rates(&model,&budget,t,y,true,options.data.structure);let mut value=|y:&[f64]|{if !budget.charge(){return Err(st::SdeError::StepLimit("MLMC evaluation budget exhausted".into()));}let physical=y.iter().zip(&frame.state_scales).map(|(v,s)|api::scale_product(*v,&[*s],&[])).collect::<Result<Vec<_>,_>>();match physical.and_then(|x|payoff.evaluate_coherent(&x)){Ok(v)=>Ok(v),Err(e)=>{budget.fail(e);Err(st::SdeError::NonFinite("MLMC payoff failed".into()))}}};let result=st::mlmc(drift,diff,&mut value,&initial,&cfg);budget.finish()?;payoff.finish()?;let r=result.map_err(api::error)?;let dim=payoff.output_dimension.clone();let squared=DimExpr::Mul(vec![dim.clone(),dim.clone()]);let levels=Value::Tensor(crate::kernel::value::TensorValue{shape:vec![r.levels.len()],elem:api::dynamic_type(),data:r.levels.into_iter().map(|l|Ok(api::record(vec![("level",Output::encode(l.level)?),("steps",Output::encode(l.steps)?),("samples",Output::encode(l.samples)?),("mean_difference",api::quantity(l.mean_diff,dim.clone())?),("variance",api::quantity(l.variance,squared.clone())?)]))).collect::<Result<_,OpError>>()?});Ok(Dynamic(api::record(vec![("estimate",api::quantity(r.estimate,dim.clone())?),("variance",api::quantity(r.variance,squared)?),("standard_error",api::quantity(r.std_error,dim)?),("levels",levels),("total_cost",Output::encode(r.total_cost)?)])))});

struct JumpPath {
    data: st::JumpPath,
    frame: OdeFrame,
}
impl NativeType for JumpPath {
    const KIND: &'static str = "rssn/native/analysis/stochastic/jump-path/v1";
    const MAX_BYTES: usize = 8 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        if self.data.path.t.len() != self.data.path.y.len()
            || self
                .data
                .path
                .y
                .len()
                .checked_mul(self.frame.state_scales.len())
                .is_none_or(|n| n > 262144)
            || self.data.path.y.iter().any(|r| {
                r.len() != self.frame.state_scales.len() || r.iter().any(|v| !v.is_finite())
            })
        {
            return Err(api::error("jump trajectory shape/storage budget exceeded"));
        }
        Ok(())
    }
}
struct JumpTape {
    snapshot: TapeSnapshot,
    params: Vec<f64>,
    input_ratios: Vec<f64>,
    output_ratios: Vec<f64>,
}
impl JumpTape {
    fn new(
        ctx: &CallCtx<'_>,
        tape: &Input,
        params: &Input,
        frame: &OdeFrame,
        mark: Option<&Distribution>,
    ) -> Result<Self, OpError> {
        let snapshot = tape_snapshot(ctx, &tape.0)?;
        let mut scales = vec![frame.time_scale];
        scales.extend(frame.state_scales.clone());
        let mut dimensions = frame.dimensions();
        if let Some(mark) = mark {
            scales.extend(std::iter::repeat_n(mark.scale, mark.n));
            dimensions.extend(std::iter::repeat_n(mark.dimension.clone(), mark.n));
        }
        if snapshot.input_scales.len() != scales.len()
            || snapshot.output_scales.len() != frame.state_scales.len()
        {
            return Err(api::error("jump/expected-jump tape shape mismatch"));
        }
        let input_ratios = scales
            .iter()
            .zip(&snapshot.input_scales)
            .zip(dimensions.iter().zip(&snapshot.input_dimensions))
            .map(|((a, b), (da, db))| {
                api::same_dimension(da, db)?;
                api::scale_product(*a, &[], &[*b])
            })
            .collect::<Result<Vec<_>, OpError>>()?;
        let output_ratios = snapshot
            .output_scales
            .iter()
            .zip(&frame.state_scales)
            .zip(
                snapshot
                    .output_dimensions
                    .iter()
                    .zip(&frame.state_dimensions),
            )
            .map(|((a, b), (da, db))| {
                api::same_dimension(da, db)?;
                api::scale_product(*a, &[], &[*b])
            })
            .collect::<Result<Vec<_>, OpError>>()?;
        let params = parameters(
            params,
            &snapshot.parameter_scales,
            &snapshot.parameter_dimensions,
        )?;
        Ok(Self {
            snapshot,
            params,
            input_ratios,
            output_ratios,
        })
    }

    fn evaluate(
        &self,
        budget: &Budget,
        t: f64,
        y: &[f64],
        mark: &[f64],
    ) -> Result<Vec<f64>, st::SdeError> {
        let run = || {
            if !budget.charge() {
                return Err(OpError::BudgetExceeded);
            }
            let x = std::iter::once(t)
                .chain(y.iter().copied())
                .chain(mark.iter().copied())
                .zip(&self.input_ratios)
                .map(|(v, s)| api::scale_product(v, &[*s], &[]))
                .collect::<Result<Vec<_>, _>>()?;
            if x.len() != self.snapshot.input_scales.len() {
                return Err(api::error("jump tape argument size mismatch"));
            }
            let mut out = vec![0.0; self.output_ratios.len()];
            self.snapshot.tape.eval(&x, &self.params, &mut out);
            out.iter()
                .zip(&self.output_ratios)
                .map(|(v, s)| api::scale_product(*v, &[*s], &[]))
                .collect::<Result<Vec<_>, _>>()
        };
        match run() {
            | Ok(v) => Ok(v),
            | Err(e) => {
                budget.fail(e);
                Err(st::SdeError::NonFinite(
                    "native jump tape evaluation failed".into(),
                ))
            },
        }
    }
}
object_binding!(jump_compensated,"num_analysis_stochastic_jump_compensated","crate::sim::stochastic::integrate_jump_diffusion_compensated","Actual native compensated jump diffusion. Poisson intensity has inverse time units; retained marks have independent physical references. Jump tape receives(time,state,mark) and outputs state increments; expected-jump tape receives(time,state) and outputs the exact supplied expected state increment, rather than evaluating at the mean mark. Native arrivals/marks remain fixed before adaptive rejection.",false,(model:Handle<Model>,options:Handle<Options>,distribution:Handle<Distribution>,intensity:Input,initial:Input,jump_tape:Input,jump_params:Input,expected_tape:Input,expected_params:Input,maxeval:usize)->NewObject<JumpPath>,|ctx,model,options,distribution,intensity,initial,jump_tape,jump_params,expected_tape,expected_params,maxeval|{let model=model.with(ctx,|m|Ok(m.clone()))?;let options=options.with(ctx,|o|Ok(o.clone()))?;let distribution=distribution.with(ctx,|d|Ok(d.clone()))?;same_frame(&model.drift.frame,&options.frame)?;let frame=&options.frame;api::same_dimension(&DimExpr::Inverse(Box::new(frame.time_dimension.clone())),&api::dimension(&intensity.0)?)?;let intensity=api::scale_product(crate::rules::numerical_abi::numeric(&intensity.0,0)?,&[frame.time_scale],&[])?;let span=match &options.data.schedule{st::Schedule::Fixed{dt,steps}=>*dt * *steps as f64,st::Schedule::Grid{times}=>times.last().copied().unwrap_or(options.data.t0)-options.data.t0,st::Schedule::Adaptive{t1,..}=>*t1-options.data.t0};
if intensity<0.0||!span.is_finite()||span<=0.0||intensity*span>4096.0{return Err(api::error("jump expected arrival work budget exceeds4096or invalid span/intensity"));}let initial=frame.state(&initial)?;let jump=JumpTape::new(ctx,&jump_tape,&jump_params,frame,Some(&distribution))?;let expected=JumpTape::new(ctx,&expected_tape,&expected_params,frame,None)?;let budget=Budget::new(maxeval)?;let drift=|t,y:&[f64]|rates(&model,&budget,t,y,false,options.data.structure);let diff=|t,y:&[f64]|rates(&model,&budget,t,y,true,options.data.structure);let mut coefficient=|t,y:&[f64],z:&[f64]|jump.evaluate(&budget,t,y,z);let mut mean=|t,y:&[f64]|expected.evaluate(&budget,t,y,&[]);let spec=st::JumpSpec{intensity,distribution:distribution.data,compensated:true};let result=st::integrate_jump_diffusion_compensated(drift,diff,&mut coefficient,&mut mean,&initial,&spec,&options.data);budget.finish()?;Ok(NewObject::new(JumpPath{data:result.map_err(api::error)?,frame:frame.clone()},frame.dimensions()))});
object_binding!(jump_data,"num_analysis_stochastic_jump_data","","Physical native jump path samples and actual total/per-step jump counts.",false,(path:Handle<JumpPath>)->Dynamic,|ctx,path|path.with(ctx,|p|{let t=&p.data.path.t;let y=&p.data.path.y;let times=Value::Tensor(crate::kernel::value::TensorValue{shape:vec![t.len()],elem:crate::kernel::signature::TypeExpr::Dim(Box::new(crate::kernel::signature::TypeExpr::Float),p.frame.time_dimension.clone()),data:t.iter().map(|v|api::scale_product(*v,&[p.frame.time_scale],&[]).map(Value::Float)).collect::<Result<_,_>>()?});let states=Value::Tensor(crate::kernel::value::TensorValue{shape:vec![y.len()],elem:api::dynamic_type(),data:y.iter().map(|y|p.frame.restore(y)).collect::<Result<_,_>>()?});Ok(Dynamic(api::record(vec![("times",times),("states",states),("total_jumps",Output::encode(p.data.total_jumps)?),("jump_counts",Output::encode(p.data.jump_counts.clone())?),("rejected",Output::encode(p.data.path.rejected)?)])))}));
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &jump_compensated::BINDING,
            &jump_data::BINDING,
            &error_new::BINDING,
            &error_code::BINDING,
            &error_capability::BINDING,
            &events::BINDING,
            &jump_distribution::BINDING,
            &jump_sample::BINDING,
            &jump_mean::BINDING,
            &mlmc::BINDING,
            &stream_new::BINDING,
            &next_unit::BINDING,
            &next_normal::BINDING,
            &method_from_name::BINDING,
            &method_name::BINDING,
            &method_calculus::BINDING,
            &structure_noise::BINDING,
            &structure_coefficients::BINDING,
            &structure_is_commutative::BINDING,
            &sde_order::BINDING,
            &sde_order_for::BINDING,
            &default_terms::BINDING,
            &ensemble_seed::BINDING,
            &options_fixed::BINDING,
            &options_fixed_span::BINDING,
            &options_adaptive::BINDING,
            &options_grid::BINDING,
            &options_with_terms::BINDING,
            &integrate::BINDING,
            &integrate_drift_implicit::BINDING,
            &ensemble::BINDING,
            &path_last::BINDING,
            &path_eval::BINDING,
            &path_resample::BINDING,
            &wiener::BINDING,
            &levy_increments::BINDING,
            &sample_poisson::BINDING,
            &sample_alpha_stable::BINDING,
        ],
        "numerical_vector",
    )
}
