//! Native quaternion/rigid-body and Hamiltonian splitting interfaces.
use super::Arg;
use super::CallCtx;
use super::DimExpr;
use super::Dynamic;
use super::Exp;
use super::Frame;
use super::Handle;
use super::Input;
use super::Installer;
use super::NativeType;
use super::NewObject;
use super::OpError;
use super::Output;
use super::Rational;
use super::RuleError;
use super::TypeExpr;
use super::Value;
use super::angular;
use super::api;
use super::c;
use super::energy;
use super::inv;
use super::l;
use super::m;
use super::number;
use super::object_binding;
use super::out;
use super::power;
use super::t;
use super::vector;
use super::velocity;
use crate::rules::numerical_abi::native_analysis_newton_cg::parameters;
use crate::rules::numerical_abi::native_analysis_ode_common::Budget;
use crate::rules::numerical_abi::native_backend_extended::TapeSnapshot;
use crate::rules::numerical_abi::native_backend_extended::tape_snapshot;
use crate::sim::symplectic as s;
#[path = "symplectic_advanced.rs"]
mod advanced;
#[path = "symplectic_ks.rs"]
mod ks;
fn scheme(name: &str) -> Result<s::Scheme, OpError> {
    Ok(match name {
        | "Verlet" => s::Scheme::Verlet,
        | "ForestRuth" => s::Scheme::ForestRuth,
        | "Yoshida4" => s::Scheme::Yoshida4,
        | "Yoshida6" => s::Scheme::Yoshida6,
        | "Yoshida8" => s::Scheme::Yoshida8,
        | "BlanesMoan4" => s::Scheme::BlanesMoan4,
        | _ => return Err(api::error("unknown native splitting Scheme variant")),
    })
}
fn order(n: usize) -> Result<(), OpError> {
    if matches!(n, 2 | 4 | 6 | 8) {
        Ok(())
    } else {
        Err(api::error(
            "composition order must be2/4/6/8within bounded native stage count",
        ))
    }
}
struct Body {
    data: s::RigidBody,
    frame: Frame,
}
impl NativeType for Body {
    const KIND: &'static str = "rssn/native/analysis/symplectic/rigid-body/v1";
    const MAX_BYTES: usize = 1024;

    fn validate(&self) -> Result<(), OpError> {
        if self
            .data
            .inertia
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.0)
            || self
                .data
                .q
                .iter()
                .chain(&self.data.l_body)
                .any(|x| !x.is_finite())
        {
            return Err(api::error("invalid native rigid-body state/inertia"));
        }
        Ok(())
    }
}
fn inertia() -> DimExpr {
    DimExpr::Mul(vec![m(), power(l(), 2)])
}
fn inertia_scale(f: &Frame) -> Result<f64, OpError> {
    api::scale_product(1.0, &[f.mass, f.length, f.length], &[])
}
#[allow(
    clippy::needless_pass_by_value,
    reason = "The typed native ABI transfers decoded arguments by value; this shared dispatcher mirrors its explicit source signature."
)]
fn array<const N: usize>(
    v: &Input,
    scale: f64,
    dim: DimExpr,
) -> Result<[f64; N], OpError> {
    api::same_dimension(&api::dimension(&v.0)?, &dim)?;
    let mut a = <[f64; N] as Arg>::decode(&v.0, 0)?;
    for v in &mut a {
        *v = api::scale_product(*v, &[], &[scale])?;
    }
    Ok(a)
}
object_binding!(scheme_order,"num_analysis_symplectic_scheme_order","crate::sim::symplectic::Scheme::order","Actual native geometric operation; physical rigid-body moments/torques use explicit positive SI M/L/T references, while orientations and composition coefficients are dimensionless.",false,(name:String)->usize,|_ctx,name|Ok(scheme(&name)?.order()));
object_binding!(scheme_stages,"num_analysis_symplectic_scheme_stages","crate::sim::symplectic::Scheme::stages","Actual native geometric operation; physical rigid-body moments/torques use explicit positive SI M/L/T references, while orientations and composition coefficients are dimensionless.",false,(name:String)->Vec<[f64;2]>,|_ctx,name|Ok(scheme(&name)?.stages().into_iter().map(|(a,b)|[a,b]).collect()));
object_binding!(composition_weights,"num_analysis_symplectic_composition_weights","crate::sim::symplectic::composition_weights","Actual native geometric operation; physical rigid-body moments/torques use explicit positive SI M/L/T references, while orientations and composition coefficients are dimensionless.",false,(degree:usize)->Vec<f64>,|_ctx,degree|{order(degree)?;s::composition_weights(degree).map_err(api::error)});
object_binding!(triple_jump,"num_analysis_symplectic_triple_jump","crate::sim::symplectic::triple_jump","Actual native geometric operation; physical rigid-body moments/torques use explicit positive SI M/L/T references, while orientations and composition coefficients are dimensionless.",false,(base:Vec<[f64;2]>,degree:usize)->Vec<[f64;2]>,|_ctx,base,degree|{order(degree)?;if base.is_empty()||base.len()>128{return Err(api::error("triple-jump stage budget must be1..128"));}let base=base.into_iter().map(|[a,b]|(a,b)).collect::<Vec<_>>();Ok(s::triple_jump(&base,degree).into_iter().map(|(a,b)|[a,b]).collect())});
object_binding!(quat_mul,"num_analysis_symplectic_quat_mul","crate::sim::symplectic::quat_mul","Actual native geometric operation; physical rigid-body moments/torques use explicit positive SI M/L/T references, while orientations and composition coefficients are dimensionless.",false,(a:[f64;4],b:[f64;4])->[f64;4],|_ctx,a,b|Ok(s::quat_mul(&a,&b)));
object_binding!(quat_to_matrix,"num_analysis_symplectic_quat_to_matrix","crate::sim::symplectic::quat_to_matrix","Actual native geometric operation; physical rigid-body moments/torques use explicit positive SI M/L/T references, while orientations and composition coefficients are dimensionless.",false,(q:[f64;4])->[[f64;3];3],|_ctx,q|Ok(s::quat_to_matrix(&q)));
object_binding!(matrix_to_quat,"num_analysis_symplectic_matrix_to_quat","crate::sim::symplectic::matrix_to_quat","Actual native geometric operation; physical rigid-body moments/torques use explicit positive SI M/L/T references, while orientations and composition coefficients are dimensionless.",false,(matrix:[[f64;3];3])->[f64;4],|_ctx,matrix|Ok(s::matrix_to_quat(&matrix)));
object_binding!(body_new,"num_analysis_symplectic_body_new","","Actual native geometric operation; physical rigid-body moments/torques use explicit positive SI M/L/T references, while orientations and composition coefficients are dimensionless.",false,(moments:Input,orientation:[f64;4],angular_momentum:Input,length_ref:Input,mass_ref:Input,time_ref:Input)->NewObject<Body>,|_ctx,moments,orientation,angular_momentum,length_ref,mass_ref,time_ref|{let frame=Frame::new(&length_ref,&mass_ref,&time_ref)?;let inertia=array(&moments,inertia_scale(&frame)?,inertia())?;let l_body=array(&angular_momentum,frame.angular()?,angular())?;let norm=orientation.iter().map(|v|v*v).sum::<f64>();if (norm-1.0).abs()>1e-8{return Err(api::error("orientation must be a finite unit quaternion"));}Ok(NewObject::new(Body{data:s::RigidBody{inertia,q:orientation,l_body},frame},vec![l(),m(),t()]))});
object_binding!(kinetic_energy,"num_analysis_symplectic_kinetic_energy","crate::sim::symplectic::RigidBody::kinetic_energy","Actual native geometric operation; physical rigid-body moments/torques use explicit positive SI M/L/T references, while orientations and composition coefficients are dimensionless.",false,(body:Handle<Body>)->Dynamic,|ctx,body|body.with(ctx,|b|Ok(Dynamic(api::quantity(api::scale_product(b.data.kinetic_energy(),&[b.frame.energy()?],&[])?,energy())?))));
object_binding!(l_space,"num_analysis_symplectic_l_space","crate::sim::symplectic::RigidBody::l_space","Actual native geometric operation; physical rigid-body moments/torques use explicit positive SI M/L/T references, while orientations and composition coefficients are dimensionless.",false,(body:Handle<Body>)->Dynamic,|ctx,body|body.with(ctx,|b|Ok(Dynamic(out(&[b.data.l_space()],b.frame.angular()?,angular())?))));
object_binding!(free_flow_exact,"num_analysis_symplectic_free_flow_exact","crate::sim::symplectic::RigidBody::free_flow_exact","Actual native geometric operation; physical rigid-body moments/torques use explicit positive SI M/L/T references, while orientations and composition coefficients are dimensionless.",false,(body:Handle<Body>,h:Input)->(),|ctx,body,h|body.with_mut(ctx,|b|{let h=number(&h,b.frame.time,t())?;b.data.free_flow_exact(h);b.validate()}));
object_binding!(free_step,"num_analysis_symplectic_free_step","crate::sim::symplectic::RigidBody::free_step","Actual native geometric operation; physical rigid-body moments/torques use explicit positive SI M/L/T references, while orientations and composition coefficients are dimensionless.",false,(body:Handle<Body>,h:Input)->(),|ctx,body,h|body.with_mut(ctx,|b|{let h=number(&h,b.frame.time,t())?;b.data.free_step(h);b.validate()}));
object_binding!(free_step_top,"num_analysis_symplectic_free_step_top","crate::sim::symplectic::RigidBody::free_step_top","Actual native geometric operation; physical rigid-body moments/torques use explicit positive SI M/L/T references, while orientations and composition coefficients are dimensionless.",false,(body:Handle<Body>,h:Input)->(),|ctx,body,h|body.with_mut(ctx,|b|{let h=number(&h,b.frame.time,t())?;b.data.free_step_top(h);b.validate()}));
object_binding!(free_step_composed,"num_analysis_symplectic_free_step_composed","crate::sim::symplectic::RigidBody::free_step_composed","Actual native geometric operation; physical rigid-body moments/torques use explicit positive SI M/L/T references, while orientations and composition coefficients are dimensionless.",false,(body:Handle<Body>,h:Input,top:bool,degree:usize)->(),|ctx,body,h,top,degree|body.with_mut(ctx,|b|{order(degree)?;let h=number(&h,b.frame.time,t())?;b.data.free_step_composed(h,top,degree).map_err(api::error)?;b.validate()}));
object_binding!(symmetric_top_flow,"num_analysis_symplectic_symmetric_top_flow","crate::sim::symplectic::RigidBody::symmetric_top_flow","Actual native geometric operation; physical rigid-body moments/torques use explicit positive SI M/L/T references, while orientations and composition coefficients are dimensionless.",false,(body:Handle<Body>,axis:usize,h:Input)->(),|ctx,body,axis,h|body.with_mut(ctx,|b|{if axis>2{return Err(api::error("rigid-body axis must be0..2"));}let h=number(&h,b.frame.time,t())?;b.data.symmetric_top_flow(axis,h);b.validate()}));
object_binding!(kick_space,"num_analysis_symplectic_kick_space","crate::sim::symplectic::RigidBody::kick_space","Actual native geometric operation; physical rigid-body moments/torques use explicit positive SI M/L/T references, while orientations and composition coefficients are dimensionless.",false,(body:Handle<Body>,h:Input,torque:Input)->(),|ctx,body,h,torque|body.with_mut(ctx,|b|{let h=number(&h,b.frame.time,t())?;let torque=array(&torque,b.frame.energy()?,energy())?;b.data.kick_space(h,&torque);b.validate()}));
#[allow(
    clippy::needless_pass_by_value,
    reason = "The typed native ABI transfers decoded arguments by value; this shared dispatcher mirrors its explicit source signature."
)]
fn torque_step(
    ctx: &CallCtx<'_>,
    body: Handle<Body>,
    h: Input,
    torque: Input,
    params: Input,
    maxeval: usize,
    kind: usize,
) -> Result<(), OpError> {
    let snapshot = tape_snapshot(ctx, &torque.0)?;
    if snapshot.input_scales.len() != 4 || snapshot.output_scales.len() != 3 {
        return Err(api::error(
            "torque tape needs4quaternion inputs and3outputs",
        ));
    }
    for d in &snapshot.input_dimensions {
        api::same_dimension(d, &DimExpr::Dimless)?;
    }
    for d in &snapshot.output_dimensions {
        api::same_dimension(d, &energy())?;
    }
    let params = parameters(
        &params,
        &snapshot.parameter_scales,
        &snapshot.parameter_dimensions,
    )?;
    body.with_mut(ctx, |b| {
        let h = number(&h, b.frame.time, t())?;
        let energy = b.frame.energy()?;
        let budget = Budget::new(maxeval)?;
        let torque = |q: &[f64; 4]| {
            let run = || {
                if !budget.charge() {
                    return Err(OpError::BudgetExceeded);
                }
                let input = q
                    .iter()
                    .zip(&snapshot.input_scales)
                    .map(|(v, s)| api::scale_product(*v, &[], &[*s]))
                    .collect::<Result<Vec<_>, _>>()?;
                let mut out = vec![0.0; 3];
                snapshot.tape.eval(&input, &params, &mut out);
                let out = out
                    .iter()
                    .zip(&snapshot.output_scales)
                    .map(|(v, s)| api::scale_product(*v, &[*s], &[energy]))
                    .collect::<Result<Vec<_>, _>>()?;
                Ok::<[f64; 3], OpError>([out[0], out[1], out[2]])
            };
            match run() {
                | Ok(v) => v,
                | Err(e) => {
                    budget.fail(e);
                    [f64::NAN; 3]
                },
            }
        };
        let result = match kind {
            | 0 => {
                b.data.step(h, &torque);
                Ok(())
            },
            | 1 => {
                b.data.step_exact(h, &torque);
                Ok(())
            },
            | 2 => b.data.rattle_step(h, &torque),
            | _ => return Err(api::error("unknown native rigid-body step")),
        };
        budget.finish()?;
        result.map_err(api::error)?;
        b.validate()
    })
}
object_binding!(step,"num_analysis_symplectic_step","crate::sim::symplectic::RigidBody::step","Actual native geometric operation; physical rigid-body moments/torques use explicit positive SI M/L/T references, while orientations and composition coefficients are dimensionless.",false,(body:Handle<Body>,h:Input,torque:Input,params:Input,maxeval:usize)->(),|ctx,body,h,torque,params,maxeval|torque_step(ctx,body,h,torque,params,maxeval,0));
object_binding!(step_exact,"num_analysis_symplectic_step_exact","crate::sim::symplectic::RigidBody::step_exact","Actual native geometric operation; physical rigid-body moments/torques use explicit positive SI M/L/T references, while orientations and composition coefficients are dimensionless.",false,(body:Handle<Body>,h:Input,torque:Input,params:Input,maxeval:usize)->(),|ctx,body,h,torque,params,maxeval|torque_step(ctx,body,h,torque,params,maxeval,1));
object_binding!(rattle_step,"num_analysis_symplectic_rattle_step","crate::sim::symplectic::RigidBody::rattle_step","Actual native geometric operation; physical rigid-body moments/torques use explicit positive SI M/L/T references, while orientations and composition coefficients are dimensionless.",false,(body:Handle<Body>,h:Input,torque:Input,params:Input,maxeval:usize)->(),|ctx,body,h,torque,params,maxeval|torque_step(ctx,body,h,torque,params,maxeval,2));
object_binding!(body_data,"num_analysis_symplectic_body_data","","Actual native geometric operation; physical rigid-body moments/torques use explicit positive SI M/L/T references, while orientations and composition coefficients are dimensionless.",false,(body:Handle<Body>)->Dynamic,|ctx,body|body.with(ctx,|b|Ok(Dynamic(api::record(vec![("orientation",Output::encode(b.data.q)?),("angular_momentum",out(&[b.data.l_body],b.frame.angular()?,angular())?)])))));
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    advanced::install(i)?;
    ks::install(i)?;
    api::install(
        i,
        &[
            &scheme_order::BINDING,
            &scheme_stages::BINDING,
            &composition_weights::BINDING,
            &triple_jump::BINDING,
            &quat_mul::BINDING,
            &quat_to_matrix::BINDING,
            &matrix_to_quat::BINDING,
            &body_new::BINDING,
            &kinetic_energy::BINDING,
            &l_space::BINDING,
            &free_flow_exact::BINDING,
            &free_step::BINDING,
            &free_step_top::BINDING,
            &free_step_composed::BINDING,
            &symmetric_top_flow::BINDING,
            &kick_space::BINDING,
            &step::BINDING,
            &step_exact::BINDING,
            &rattle_step::BINDING,
            &body_data::BINDING,
        ],
        "numerical_vector",
    )
}
