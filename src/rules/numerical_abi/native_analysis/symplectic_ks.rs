//! Native KS regularization with explicit Cartesian and fictitious-time scales.
use super::Budget;
use super::CallCtx;
use super::DimExpr;
use super::Dynamic;
use super::Exp;
use super::Frame;
use super::Input;
use super::Installer;
use super::OpError;
use super::Rational;
use super::RuleError;
use super::TapeSnapshot;
use super::TypeExpr;
use super::Value;
use super::api;
use super::array;
use super::inv;
use super::l;
use super::number;
use super::object_binding;
use super::parameters;
use super::power;
use super::s;
use super::t;
use super::tape_snapshot;
use super::velocity;
fn frame(
    length_ref: &Input,
    time_ref: &Input,
) -> Result<Frame, OpError> {
    let (length, ld) = api::reference(&length_ref.0)?;
    let (time, td) = api::reference(&time_ref.0)?;
    api::same_dimension(&ld, &l())?;
    api::same_dimension(&td, &t())?;
    Ok(Frame {
        length,
        time,
        mass: 1.0,
    })
}
fn udim() -> DimExpr {
    l().pow(Exp::Rat(Rational::new(1, 2)))
}
fn updim() -> DimExpr {
    DimExpr::Mul(vec![l().pow(Exp::Rat(Rational::new(3, 2))), inv(t())])
}
fn up_scale(f: &Frame) -> Result<f64, OpError> {
    api::scale_product(f.length.sqrt(), &[f.length], &[f.time])
}
fn value<const N: usize>(
    a: &[f64; N],
    scale: f64,
    dim: DimExpr,
) -> Result<Value, OpError> {
    Ok(Value::Tensor(crate::kernel::value::TensorValue {
        shape: vec![N],
        elem: TypeExpr::Dim(Box::new(TypeExpr::Float), dim),
        data: a
            .iter()
            .map(|v| api::scale_product(*v, &[scale], &[]).map(Value::Float))
            .collect::<Result<_, _>>()?,
    }))
}
object_binding!(from_cartesian,"num_analysis_symplectic_ks_from_cartesian","crate::sim::symplectic::ks_from_cartesian","Actual native Cartesian-to-KS transform. Cartesian units are L,L/T; regular coordinates have sqrt(L) and du/ds has L^(3/2)/T, with fictitious time ds=T/L.",false,(position:Input,velocity:Input,length_ref:Input,time_ref:Input)->Dynamic,|_ctx,position,velocity,length_ref,time_ref|{let f=frame(&length_ref,&time_ref)?;let x=array(&position,f.length,l())?;let v=array(&velocity,f.vel()?,super::velocity())?;if x.iter().all(|v|*v==0.0){return Err(api::error("KS transform requires nonzero separation"));}let(u,up)=s::ks_from_cartesian(&x,&v);Ok(Dynamic(api::record(vec![("u",value(&u,f.length.sqrt(),udim())?),("up",value(&up,up_scale(&f)?,updim())?)])))});
object_binding!(to_cartesian,"num_analysis_symplectic_ks_to_cartesian","crate::sim::symplectic::ks_to_cartesian","Actual native KS-to-Cartesian transform with physical sqrt(L),L^(3/2)/T regular coordinates and L,L/T outputs.",false,(u:Input,up:Input,length_ref:Input,time_ref:Input)->Dynamic,|_ctx,u,up,length_ref,time_ref|{let f=frame(&length_ref,&time_ref)?;let u=array(&u,f.length.sqrt(),udim())?;let up=array(&up,up_scale(&f)?,updim())?;if u.iter().all(|v|*v==0.0){return Err(api::error("KS transform requires nonzero regular coordinate"));}let(x,v)=s::ks_to_cartesian(&u,&up);Ok(Dynamic(api::record(vec![("position",value(&x,f.length,l())?),("velocity",value(&v,f.vel()?,super::velocity())?)])))});
object_binding!(matrix,"num_analysis_symplectic_ks_matrix","crate::sim::symplectic::ks_matrix","Actual native KS matrix, linear in regular coordinates; input and matrix entries have sqrt(length) units.",false,(u:Input,length_ref:Input)->Dynamic,|_ctx,u,length_ref|{let(length,dimension)=api::reference(&length_ref.0)?;api::same_dimension(&dimension,&l())?;let u=array(&u,length.sqrt(),udim())?;let matrix=s::ks_matrix(&u);Ok(Dynamic(Value::Tensor(crate::kernel::value::TensorValue{shape:vec![4,4],elem:TypeExpr::Dim(Box::new(TypeExpr::Float),udim()),data:matrix.iter().flatten().map(|v|api::scale_product(*v,&[length.sqrt()],&[]).map(Value::Float)).collect::<Result<_,_>>()?})))});
fn perturbation(
    snapshot: &TapeSnapshot,
    params: &[f64],
    f: &Frame,
    budget: &Budget,
    x: &[f64; 3],
    v: &[f64; 3],
    time: f64,
) -> [f64; 3] {
    let run = || {
        if !budget.charge() {
            return Err(OpError::BudgetExceeded);
        }
        let scales = [
            f.length,
            f.length,
            f.length,
            f.vel()?,
            f.vel()?,
            f.vel()?,
            f.time,
        ];
        let input = x
            .iter()
            .chain(v)
            .copied()
            .chain(std::iter::once(time))
            .zip(scales)
            .zip(&snapshot.input_scales)
            .map(|((v, s), r)| api::scale_product(v, &[s], &[*r]))
            .collect::<Result<Vec<_>, _>>()?;
        let mut out = [0.0; 3];
        snapshot.tape.eval(&input, params, &mut out);
        for (v, r) in out.iter_mut().zip(&snapshot.output_scales) {
            *v = api::scale_product(*v, &[*r, f.time, f.time], &[f.length])?;
        }
        Ok::<_, OpError>(out)
    };
    match run() {
        | Ok(v) => v,
        | Err(e) => {
            budget.fail(e);
            [f64::NAN; 3]
        },
    }
}
fn tape(
    ctx: &CallCtx<'_>,
    input: &Input,
    params: &Input,
) -> Result<(TapeSnapshot, Vec<f64>), OpError> {
    let snapshot = tape_snapshot(ctx, &input.0)?;
    if snapshot.input_scales.len() != 7 || snapshot.output_scales.len() != 3 {
        return Err(api::error(
            "KS perturbation tape needs x3,v3,time inputs and acceleration3outputs",
        ));
    }
    for (i, d) in snapshot.input_dimensions.iter().enumerate() {
        api::same_dimension(
            d,
            &if i < 3 {
                l()
            } else if i < 6 {
                velocity()
            } else {
                t()
            },
        )?;
    }
    for d in &snapshot.output_dimensions {
        api::same_dimension(d, &DimExpr::Mul(vec![l(), power(t(), -2)]))?;
    }
    let params = parameters(
        params,
        &snapshot.parameter_scales,
        &snapshot.parameter_dimensions,
    )?;
    Ok((snapshot, params))
}
object_binding!(advance,"num_analysis_symplectic_ks_advance","crate::sim::symplectic::ks_advance","Actual native physical-time KS advance with explicit positive gravitational parameter mu(L³/T²), bounded eta and a genuine tape perturbation acceleration. Callback budget failures and incomplete native endpoint errors are preserved.",false,(position:Input,velocity:Input,mu:Input,perturbation_tape:Input,params:Input,dt:Input,eta:f64,length_ref:Input,time_ref:Input,maxeval:usize)->Dynamic,|ctx,position,velocity,mu,perturbation_tape,params,dt,eta,length_ref,time_ref,maxeval|{let f=frame(&length_ref,&time_ref)?;let x=array(&position,f.length,l())?;let v=array(&velocity,f.vel()?,super::velocity())?;let mu=number(&mu,api::scale_product(1.0,&[f.length,f.length,f.length],&[f.time,f.time])?,DimExpr::Mul(vec![power(l(),3),power(t(),-2)]))?;let dt=number(&dt,f.time,t())?;if !eta.is_finite()||!(1e-6..=1.0).contains(&eta){return Err(api::error("KS eta must be1e-6..1"));}let(snapshot,params)=tape(ctx,&perturbation_tape,&params)?;let budget=Budget::new(maxeval)?;let pert=|x:&[f64;3],v:&[f64;3],time:f64|perturbation(&snapshot,&params,&f,&budget,x,v,time);let result=s::ks_advance(&x,&v,mu,&pert,dt,eta);budget.finish()?;let(x,v)=result.map_err(api::error)?;Ok(Dynamic(api::record(vec![("position",value(&x,f.length,l())?),("velocity",value(&v,f.vel()?,super::velocity())?)])))});
object_binding!(propagate,"num_analysis_symplectic_ks_propagate","crate::sim::symplectic::ks_propagate","Actual native fixed fictitious-time KS propagation. Native central mu=1is normalized by an explicit positive physical mu reference: time reference=sqrt(Lref³/mu). Fictitious step has T/L units, outputs retain regular-coordinate/time/specific-energy units.",false,(position:Input,velocity:Input,mu_reference:Input,perturbation_tape:Input,params:Input,ds:Input,steps:usize,length_ref:Input,maxeval:usize)->Dynamic,|ctx,position,velocity,mu_reference,perturbation_tape,params,ds,steps,length_ref,maxeval|{let(length,ld)=api::reference(&length_ref.0)?;api::same_dimension(&ld,&l())?;let(mu,md)=api::reference(&mu_reference.0)?;api::same_dimension(&md,&DimExpr::Mul(vec![power(l(),3),power(t(),-2)]))?;let time=api::scale_product(length.sqrt(),&[length],&[mu.sqrt()])?;let f=Frame{length,time,mass:1.0};let x=array(&position,f.length,l())?;let v=array(&velocity,f.vel()?,super::velocity())?;if x.iter().all(|v|*v==0.0)||steps==0||steps>100000{return Err(api::error("KS initial separation/step budget invalid"));}let ds=number(&ds,api::scale_product(f.time,&[],&[f.length])?,DimExpr::Mul(vec![t(),inv(l())]))?;let(snapshot,params)=tape(ctx,&perturbation_tape,&params)?;let budget=Budget::new(maxeval)?;let pert=|x:&[f64;3],v:&[f64;3],time:f64|perturbation(&snapshot,&params,&f,&budget,x,v,time);let result=s::ks_propagate(&x,&v,&pert,ds,steps);budget.finish()?;Ok(Dynamic(api::record(vec![("u",value(&result.u,f.length.sqrt(),udim())?),("up",value(&result.up,up_scale(&f)?,updim())?),("time",api::quantity(api::scale_product(result.t,&[f.time],&[])?,t())?),("specific_energy",api::quantity(api::scale_product(result.energy,&[f.length,f.length],&[f.time,f.time])?,DimExpr::Mul(vec![power(l(),2),power(t(),-2)]))?)])))});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &from_cartesian::BINDING,
            &to_cartesian::BINDING,
            &matrix::BINDING,
            &advance::BINDING,
            &propagate::BINDING,
        ],
        "numerical_vector",
    )
}
