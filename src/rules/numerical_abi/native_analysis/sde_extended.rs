//! Native drift-implicit steps and stochastic samplers.
use std::num::NonZeroUsize;

use super::Budget;
use super::CallCtx;
use super::Dynamic;
use super::Handle;
use super::Input;
use super::Installer;
use super::Many;
use super::Model;
use super::NativeType;
use super::NewObject;
use super::Noise;
use super::OpError;
use super::Rng;
use super::RuleError;
use super::api;
use super::object_binding;
use crate::kernels::numeric_common::trajectory_sample_count;
use crate::kernels::sde_ext as e;
pub(super) fn scheme(
    name: &str,
    theta: f64,
) -> Result<e::DriftImplicitScheme, OpError> {
    Ok(match name {
        | "BackwardEuler" => e::DriftImplicitScheme::BackwardEuler,
        | "Trapezoidal" => e::DriftImplicitScheme::Trapezoidal,
        | "Theta" => {
            if !theta.is_finite() || !(0.0..=1.0).contains(&theta) {
                return Err(api::error("theta must be0..1"));
            }
            e::DriftImplicitScheme::Theta(theta)
        },
        | "DriftImplicitMilstein" => e::DriftImplicitScheme::DriftImplicitMilstein,
        | "SplitStepBackwardEuler" => e::DriftImplicitScheme::SplitStepBackwardEuler,
        | _ => return Err(api::error("unknown native drift implicit scheme")),
    })
}
#[derive(Clone, Copy)]
enum ImplicitSelection {
    FullPath,
    SampledPath(NonZeroUsize),
    Step(usize),
}

enum BrownianInput {
    Path(Vec<Vec<f64>>, Option<NonZeroUsize>),
    Step {
        increment: Vec<f64>,
        areas: Vec<f64>,
        area_scale: f64,
    },
}

fn scale_noise(
    values: &mut [f64],
    factor: f64,
) -> Result<(), OpError> {
    if factor.to_bits() != 1.0_f64.to_bits() {
        for value in values {
            *value = api::scale_product(*value, &[factor], &[])?;
        }
    }
    Ok(())
}

#[allow(
    clippy::too_many_arguments,
    reason = "the dispatcher mirrors the physical solver inputs"
)]
fn implicit(
    ctx: &CallCtx<'_>,
    model: &Handle<Model>,
    noise: &Handle<Noise>,
    initial: &Input,
    start: &Input,
    name: &str,
    theta: f64,
    maxeval: usize,
    selection: ImplicitSelection,
) -> Result<Many<Dynamic>, OpError> {
    let model = model.with(ctx, |m| Ok(m.clone()))?;
    let frame = &model.drift.frame;
    let y = frame.state(initial)?;
    let t = frame.time(start)?;
    let (h, increment_scale, input) = noise.with(ctx, |noise| {
        api::same_dimension(&frame.time_dimension, &noise.time_dimension)?;
        if noise.m != model.m {
            return Err(api::error("implicit SDE noise shape mismatch"));
        }
        let steps = noise.sample.dws.len();
        let retained = match selection {
            | ImplicitSelection::FullPath => Some(steps),
            | ImplicitSelection::SampledPath(stride) => trajectory_sample_count(steps, stride),
            | ImplicitSelection::Step(_) => Some(1),
        };
        if retained
            .and_then(|n| n.checked_mul(y.len()))
            .is_none_or(|n| n > 262000)
        {
            return Err(api::error("implicit SDE sample storage exceeded"));
        }
        let h = api::scale_product(noise.h, &[noise.time_scale], &[frame.time_scale])?;
        let increment_scale =
            api::scale_product(noise.time_scale.sqrt(), &[], &[frame.time_scale.sqrt()])?;
        let input = match selection {
            | ImplicitSelection::FullPath => BrownianInput::Path(noise.sample.dws.clone(), None),
            | ImplicitSelection::SampledPath(stride) => {
                BrownianInput::Path(noise.sample.dws.clone(), Some(stride))
            },
            | ImplicitSelection::Step(index) => {
                let (Some(increment), Some(areas)) =
                    (noise.sample.dws.get(index), noise.sample.areas.get(index))
                else {
                    return Err(api::error("noise step index out of range"));
                };
                BrownianInput::Step {
                    increment: increment.clone(),
                    areas: areas.clone(),
                    area_scale: api::scale_product(noise.time_scale, &[], &[frame.time_scale])?,
                }
            },
        };
        Ok((h, increment_scale, input))
    })?;
    let scheme = scheme(name, theta)?;
    let budget = Budget::new(maxeval)?;
    let drift = |t, y: &[f64], out: &mut [f64]| model.drift.rate(&budget, t, y, out);
    let diff = |t, y: &[f64], out: &mut [f64]| model.diffusion(&budget, t, y, out);
    let result = match input {
        | BrownianInput::Path(mut increments, stride) => {
            for row in &mut increments {
                scale_noise(row, increment_scale)?;
            }
            match stride {
                | Some(stride) => {
                    e::drift_implicit_sampled(
                        &drift,
                        &diff,
                        model.m,
                        t,
                        h,
                        &y,
                        &increments,
                        scheme,
                        stride.get(),
                    )
                },
                | None => {
                    e::drift_implicit_core(&drift, &diff, model.m, t, h, &y, &increments, scheme)
                },
            }
        },
        | BrownianInput::Step {
            mut increment,
            mut areas,
            area_scale,
        } => {
            scale_noise(&mut increment, increment_scale)?;
            scale_noise(&mut areas, area_scale)?;
            e::drift_implicit_step_with_areas(
                &drift, &diff, model.m, t, h, &y, &increment, &areas, scheme,
            )
            .map(|state| vec![state])
        },
    };
    budget.finish()?;
    Ok(Many(
        result
            .map_err(api::error)?
            .iter()
            .map(|y| frame.restore(y).map(Dynamic))
            .collect::<Result<_, _>>()?,
    ))
}
object_binding!(core,"num_analysis_sde_drift_implicit_core","crate::kernels::sde_ext::drift_implicit_core","Native drift-implicit path over retained Brownian increments. All five scheme variants use physical time and state references and share a checked evaluation budget.",false,(model:Handle<Model>,noise:Handle<Noise>,initial:Input,start:Input,scheme:String,theta:f64,maxeval:usize)->Many<Dynamic>,|ctx,model,noise,initial,start,scheme,theta,maxeval|implicit(ctx,&model,&noise,&initial,&start,&scheme,theta,maxeval,ImplicitSelection::FullPath));
object_binding!(step,"num_analysis_sde_drift_implicit_step_areas","crate::kernels::sde_ext::drift_implicit_step_with_areas","Native drift-implicit single step using the specified retained conditional Levy-area sample and its iterated integrals.",false,(model:Handle<Model>,noise:Handle<Noise>,initial:Input,start:Input,scheme:String,theta:f64,index:usize,maxeval:usize)->Many<Dynamic>,|ctx,model,noise,initial,start,scheme,theta,index,maxeval|implicit(ctx,&model,&noise,&initial,&start,&scheme,theta,maxeval,ImplicitSelection::Step(index)));
object_binding!(
    sampled,
    "num_analysis_sde_drift_implicit_sampled",
    "crate::kernels::sde_ext::drift_implicit_sampled",
    "Runs all drift-implicit steps over the retained Brownian increments and returns the initial state, every state at a positive multiple of the stride, and the final state once. State units and callback evaluation limits match the full path.",
    false,
    (model: Handle<Model>, noise: Handle<Noise>, initial: Input, start: Input, scheme: String, theta: f64, stride: usize, maxeval: usize) -> Many<Dynamic>,
    |ctx, model, noise, initial, start, scheme, theta, stride, maxeval| {
        let stride = NonZeroUsize::new(stride).ok_or_else(|| api::error("the snapshot stride must be positive"))?;
        implicit(ctx, &model, &noise, &initial, &start, &scheme, theta, maxeval, ImplicitSelection::SampledPath(stride))
    }
);
object_binding!(poisson,"num_analysis_sde_sample_poisson_with","crate::kernels::sde_ext::sample_poisson_with","Native callback-based Poisson sampler instantiated with the retained native RNG's actual unit/normal streams.",true,(rng:Handle<Rng>,mean:f64)->usize,|ctx,rng,mean|{if !mean.is_finite()||!(0.0..=100000.0).contains(&mean){return Err(api::error("Poisson mean must be finite0..100000"));}rng.with_mut(ctx,|rng|{let cell=std::cell::RefCell::new(rng);let mut unit=||cell.borrow_mut().uniform();let mut normal=||cell.borrow_mut().normal();Ok(e::sample_poisson_with(&mut unit,&mut normal,mean))})});
object_binding!(stable,"num_analysis_sde_sample_alpha_stable_with","crate::kernels::sde_ext::sample_alpha_stable_with","Native callback-based alpha-stable sampler using the retained RNG. Scale and location have the same physical domain; alpha/beta are dimensionless.",false,(rng:Handle<Rng>,alpha:f64,beta:f64,scale:Input,location:Input)->Dynamic,|ctx,rng,alpha,beta,scale,location|{let(scale,dim)=api::reference(&scale.0)?;api::same_dimension(&dim,&api::dimension(&location.0)?)?;let location=api::scale_product(crate::rules::numerical_abi::numeric(&location.0,0)?,&[],&[scale])?;let v=rng.with_mut(ctx,|rng|{let cell=std::cell::RefCell::new(rng);let mut unit=||cell.borrow_mut().uniform();let mut normal=||cell.borrow_mut().normal();e::sample_alpha_stable_with(&mut unit,&mut normal,alpha,beta,1.0,location).map_err(api::error)})?;Ok(Dynamic(api::quantity(api::scale_product(v,&[scale],&[])?,dim)?))});
struct Error(e::SdeExtError);
impl NativeType for Error {
    const KIND: &'static str = "rssn/native/analysis/sde/error/v1";
    const MAX_BYTES: usize = 4096;

    fn validate(&self) -> Result<(), OpError> {
        if self.0.to_string().len() > 4096 {
            Err(api::error("native error text budget exceeded"))
        } else {
            Ok(())
        }
    }
}
object_binding!(invalid,"num_analysis_sde_error_invalid","crate::kernels::sde_ext::SdeExtError::invalid","Construct the actual native invalid-parameter error object; text is bounded and preserved.",true,(message:String)->NewObject<Error>,|_ctx,message|{if message.len()>2048{return Err(api::error("error message budget exceeds2048bytes"));}Ok(NewObject::normalized(Error(e::SdeExtError::invalid(message))))});
object_binding!(message,"num_analysis_sde_error_message","","Inspect the retained native SDE extension error text.",true,(error:Handle<Error>)->String,|ctx,error|error.with(ctx,|e|Ok(e.0.to_string())));
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &core::BINDING,
            &step::BINDING,
            &sampled::BINDING,
            &poisson::BINDING,
            &stable::BINDING,
            &invalid::BINDING,
            &message::BINDING,
        ],
        "numerical_vector",
    )
}
