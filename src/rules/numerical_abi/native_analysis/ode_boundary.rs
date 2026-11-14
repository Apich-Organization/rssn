//! Actual native shooting and finite-difference boundary-value adapters.
use super::Budget;
use super::CompiledSnapshot;
use super::Dynamic;
use super::Handle;
use super::Input;
use super::Installer;
use super::Many;
use super::OdeTape;
use super::OpError;
use super::RuleError;
use super::Value;
use super::a;
use super::api;
use super::compiled_snapshot;
use super::object_binding;
use super::options;
use crate::kernel::signature::DimExpr;
use crate::kernels::ode_geom as g;
use crate::rules::numerical_abi::native_analysis_ode_common::VectorOracle;
use crate::rules::numerical_abi::native_analysis_ode_common::reference_frame;
object_binding!(shoot_initial,"num_analysis_ode_shoot_initial","crate::kernels::ode_geom::shoot_initial","Actual native damped shooting-parameter solve. Positive independent parameter/residual references define the normalized Newton metric; returned parameters retain their own physical units. Original tape failures take precedence over native nonconvergence.",false,(time_reference:Input,parameter_references:Many<Input>,residual_references:Many<Input>,residual:Input,params:Input,guess:Input,maxeval:usize)->Dynamic,|ctx,time_reference,parameter_references,residual_references,residual,params,guess,maxeval|{let p=reference_frame(&time_reference,&parameter_references)?;let r=reference_frame(&time_reference,&residual_references)?;if p.state_scales.len()!=r.state_scales.len(){return Err(api::error("shooting parameter/residual shape mismatch"));}let oracle=VectorOracle::new(ctx,&residual,&params,&p.state_scales,&p.state_dimensions,&r.state_scales,&r.state_dimensions)?;let initial=p.state(&guess)?;let budget=Budget::new(maxeval)?;let mut residual=|x:&[f64]|match oracle.evaluate(&budget,x){Ok(v)=>Ok(v),Err(e)=>{budget.fail(e.clone());Err(e)}};let result=g::shoot_initial(&mut residual,&initial);budget.finish()?;let x=result.map_err(|e|match e{g::ShootingError::Residual(e)=>e,g::ShootingError::Numerical(e)=>api::error(e)})?;Ok(Dynamic(p.restore(&x)?))});
object_binding!(shooting,"num_analysis_ode_shooting","crate::kernels::ode_adaptive::shooting","Actual native ODE shooting with true initial-state and terminal-residual tapes, independently measured parameters/residuals and original errors. RHS uses its mixed physical state frame; config controls inner native integration rather than an unrelated substitute.",false,(rhs:Handle<OdeTape>,parameter_references:Many<Input>,residual_references:Many<Input>,init:Input,init_params:Input,residual:Input,residual_params:Input,guess:Input,start:Input,end:Input,config:Input,maxeval:usize)->Dynamic,|ctx,rhs,parameter_references,residual_references,init,init_params,residual,residual_params,guess,start,end,config,maxeval|{let m=rhs.with(ctx,|m|Ok(m.clone()))?;let time_reference=Input(api::quantity(m.frame.time_scale,m.frame.time_dimension.clone())?);let p=reference_frame(&time_reference,&parameter_references)?;let r=reference_frame(&time_reference,&residual_references)?;if p.state_scales.len()!=r.state_scales.len(){return Err(api::error("shooting parameter/residual size mismatch"));}let initial=VectorOracle::new(ctx,&init,&init_params,&p.state_scales,&p.state_dimensions,&m.frame.state_scales,&m.frame.state_dimensions)?;let residual=VectorOracle::new(ctx,&residual,&residual_params,&m.frame.state_scales,&m.frame.state_dimensions,&r.state_scales,&r.state_dimensions)?;let guess=p.state(&guess)?;let o=options(&config,&m.frame)?;let budget=Budget::new(maxeval)?;let rate=|t,y:&[f64],out:&mut[f64]|m.rate(&budget,t,y,out);let init=|p:&[f64]|initial.sample(&budget,p);let residual=|y:&[f64]|residual.sample(&budget,y);let result=a::shooting(rate,m.frame.time(&start)?,m.frame.time(&end)?,init,residual,&guess,&o);budget.finish()?;Ok(Dynamic(p.restore(&result.map_err(api::error)?)?))});
object_binding!(bvp_fd,"num_analysis_ode_bvp_fd","crate::kernels::ode_adaptive::bvp_fd","Actual native finite-difference BVP y''=f(x,y,y') with physical coordinateD/stateE references. f inputs have domains D,E,E/D and output E/D²; guess is D->E. Boundary values, true native samples and explicit evaluation budget retain those domains.",false,(rhs:Input,guess:Input,start:Input,end:Input,y_start:Input,y_end:Input,coordinate_reference:Input,state_reference:Input,n:usize,maxeval:usize)->Dynamic,|ctx,rhs,guess,start,end,y_start,y_end,coordinate_reference,state_reference,n,maxeval|{let(x,xd)=api::reference(&coordinate_reference.0)?;let(y,yd)=api::reference(&state_reference.0)?;if n==0||n>256{return Err(api::error("BVP interior grid count must be1..256"));}let rhs=compiled_snapshot(ctx,&rhs.0,maxeval)?;let guess=compiled_snapshot(ctx,&guess.0,maxeval)?;let slope_dim=DimExpr::Mul(vec![yd.clone(),DimExpr::Inverse(Box::new(xd.clone()))]);let curvature_dim=DimExpr::Mul(vec![slope_dim.clone(),DimExpr::Inverse(Box::new(xd.clone()))]);if rhs.input_dimensions.len()!=3||guess.input_dimensions.len()!=1{return Err(api::error("BVP RHS needs3inputs and guess needs1"));}for(a,b)in rhs.input_dimensions.iter().zip([xd.clone(),yd.clone(),slope_dim]){api::same_dimension(a,&b)?;}api::same_dimension(&rhs.output_dimension,&curvature_dim)?;api::same_dimension(&guess.input_dimensions[0],&xd)?;api::same_dimension(&guess.output_dimension,&yd)?;let value=|v:&Input,scale:f64,dim:&DimExpr|{api::same_dimension(&api::dimension(&v.0)?,dim)?;api::scale_product(crate::rules::numerical_abi::numeric(&v.0,0)?,&[],&[scale])};let start=value(&start,x,&xd)?;let end=value(&end,x,&xd)?;let ya=value(&y_start,y,&yd)?;let yb=value(&y_end,y,&yd)?;if end<=start{return Err(api::error("BVP span must increase"));}let budget=Budget::new(maxeval)?;let sample=|f:&CompiledSnapshot,coords:Result<Vec<f64>,OpError>,numerator:&[f64],denominator:&[f64]|{let run=||{if !budget.charge(){return Err(OpError::BudgetExceeded);}let coords=coords?;api::scale_product(f.evaluate_coherent(&coords)?,numerator,denominator)};match run(){Ok(v)=>v,Err(e)=>{budget.fail(e);f64::NAN}}};let f=|xx,yy,dy|sample(&rhs,[api::scale_product(xx,&[x],&[]),api::scale_product(yy,&[y],&[]),api::scale_product(dy,&[y],&[x])].into_iter().collect(),&[x,x],&[y]);let g=|xx|sample(&guess,api::scale_product(xx,&[x],&[]).map(|v|vec![v]),&[],&[y]);let result=a::bvp_fd(f,start,end,ya,yb,n,g);budget.finish()?;rhs.finish()?;guess.finish()?;let(xs,ys)=result.map_err(api::error)?;let vector=|values:Vec<f64>,scale,dim|Ok::<_,OpError>(Value::Tensor(crate::kernel::value::TensorValue{shape:vec![values.len()],elem:crate::kernel::signature::TypeExpr::Dim(Box::new(crate::kernel::signature::TypeExpr::Float),dim),data:values.into_iter().map(|v|api::scale_product(v,&[scale],&[]).map(Value::Float)).collect::<Result<_,_>>()?}));Ok(Dynamic(api::record(vec![("coordinates",vector(xs,x,xd)?),("values",vector(ys,y,yd)?)])))});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &shoot_initial::BINDING,
            &shooting::BINDING,
            &bvp_fd::BINDING,
        ],
        "numerical_vector",
    )
}
