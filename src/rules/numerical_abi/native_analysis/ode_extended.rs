//! Native exponential, IMEX and variational ODE methods.
use super::Arg;
use super::Budget;
use super::CallCtx;
use super::Dynamic;
use super::Handle;
use super::Input;
use super::Installer;
use super::Many;
use super::NewObject;
use super::OdeFrame;
use super::OdeTape;
use super::OpError;
use super::RuleError;
use super::Solution;
use super::SolutionData;
use super::Value;
use super::api;
use super::object_binding;
use super::options;
use crate::kernels::dense::Mat;
use crate::kernels::ode_ext as e;
use crate::rules::numerical_abi::native_analysis_jacobian::MatrixData;
use crate::rules::numerical_abi::native_analysis_jacobian::MatrixObject;
use crate::rules::numerical_abi::native_analysis_ode_common::VectorOracle;
use crate::rules::numerical_abi::native_analysis_physics_ode::Acceleration;
fn matrix(
    ctx: &CallCtx<'_>,
    value: &Handle<MatrixObject>,
) -> Result<Mat, OpError> {
    value.with(ctx, |m| {
        let MatrixData::Dense(m) = &m.data else {
            return Err(api::error(
                "this native exponential method needs a dense normalized matrix",
            ));
        };
        Ok(m.clone())
    })
}
pub(super) fn same_frame(
    a: &OdeFrame,
    b: &OdeFrame,
) -> Result<(), OpError> {
    if a.state_scales != b.state_scales
        || !crate::kernels::exact_eq(a.time_scale, b.time_scale)
        || a.state_dimensions.len() != b.state_dimensions.len()
    {
        return Err(api::error("split RHS reference frames must agree"));
    }
    api::same_dimension(&a.time_dimension, &b.time_dimension)?;
    for (a, b) in a.state_dimensions.iter().zip(&b.state_dimensions) {
        api::same_dimension(a, b)?;
    }
    Ok(())
}
object_binding!(exp_matrix,"num_analysis_ode_exp_matrix","crate::kernels::ode_ext::exp_matrix","Actual native matrix exponential on explicitly dimensionless normalized matrix coefficients. Physical-rate Jacobians must first be multiplied by a matching physical time in the declared metric.",false,(a:Handle<MatrixObject>)->NewObject<MatrixObject>,|ctx,a|{let m=a.with(ctx,|m|{if m.frame.is_some(){return Err(api::error("exp_matrix requires an explicitly dimensionless matrix, not a rate Jacobian"));}Ok(m.clone())})?;let MatrixData::Dense(a)=m.data else{return Err(api::error("exp_matrix requires dense matrix"));};let out=e::exp_matrix(&a).map_err(api::error)?;Ok(NewObject::normalized(MatrixObject{data:MatrixData::Dense(out),frame:None,evaluations:0}))});
object_binding!(exprb32,"num_analysis_ode_exprb32","crate::kernels::ode_ext::exprb32","Actual native adaptive exponential Rosenbrock3/2 with exact native tape Jacobian, physical mixed state reconstruction and explicit attempts/evaluation budgets.",false,(rhs:Handle<OdeTape>,initial:Input,start:Input,end:Input,config:Input,maxeval:usize)->NewObject<Solution>,|ctx,rhs,initial,start,end,config,maxeval|{let m=rhs.with(ctx,|m|Ok(m.clone()))?;let frame=m.frame.clone();let y=frame.state(&initial)?;let o=options(&config,&frame)?;let budget=Budget::new(maxeval)?;let rate=|t,y:&[f64],out:&mut[f64]|m.rate(&budget,t,y,out);let jac=|t,y:&[f64]|m.jacobian(&budget,t,y);let result=e::exprb32(rate,Some(&jac),frame.time(&start)?,frame.time(&end)?,&y,&o);budget.finish()?;Ok(NewObject::new(Solution{frame:frame.clone(),data:SolutionData::Dense(result.map_err(api::error)?),evaluations:budget.used()},frame.dimensions()))});
#[allow(
    clippy::needless_pass_by_value,
    clippy::too_many_arguments,
    reason = "The typed native ABI transfers decoded arguments by value; this shared dispatcher mirrors its explicit source signature."
)]
fn exponential(
    ctx: &CallCtx<'_>,
    rhs: Handle<OdeTape>,
    linear: Option<Handle<MatrixObject>>,
    diagonal: Input,
    initial: Input,
    start: Input,
    end: Input,
    steps: usize,
    maxeval: usize,
    diag: bool,
) -> Result<Many<Dynamic>, OpError> {
    let m = rhs.with(ctx, |m| Ok(m.clone()))?;
    let frame = &m.frame;
    let y = frame.state(&initial)?;
    if steps == 0 || steps.checked_mul(y.len() + 2).is_none_or(|n| n > 200000) {
        return Err(api::error("ETDRK4 output/work bound exceeded"));
    }
    let budget = Budget::new(maxeval)?;
    let rate = |t, y: &[f64], out: &mut [f64]| m.rate(&budget, t, y, out);
    let t0 = frame.time(&start)?;
    let t1 = frame.time(&end)?;
    let result = if diag {
        api::same_dimension(
            &api::dimension(&diagonal.0)?,
            &crate::kernel::signature::DimExpr::Inverse(Box::new(frame.time_dimension.clone())),
        )?;
        let lambda = <Vec<f64> as Arg>::decode(&diagonal.0, 0)?
            .into_iter()
            .map(|v| api::scale_product(v, &[frame.time_scale], &[]))
            .collect::<Result<Vec<_>, _>>()?;
        if lambda.len() != y.len() {
            return Err(api::error("ETD diagonal/state size mismatch"));
        }
        e::etdrk4_diag(&lambda, rate, t0, t1, &y, steps)
    } else {
        let linear = linear
            .as_ref()
            .ok_or_else(|| api::error("ETD matrix missing"))?;
        let a = matrix(ctx, linear)?;
        if a.rows != y.len() || a.cols != y.len() {
            return Err(api::error("ETD matrix/state size mismatch"));
        }
        linear.with(ctx, |a| {
            if let Some(f) = &a.frame {
                same_frame(frame, f)?;
            }
            Ok(())
        })?;
        e::etdrk4(&a, rate, t0, t1, &y, steps)
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
object_binding!(etdrk4,"num_analysis_ode_etdrk4","crate::kernels::ode_ext::etdrk4","Actual native matrix ETDRK4 with independent physical state/time reference frame. Linear matrix coefficients are explicitly normalized in that same metric; RHS tape is the actual nonlinear term.",false,(rhs:Handle<OdeTape>,linear:Handle<MatrixObject>,initial:Input,start:Input,end:Input,steps:usize,maxeval:usize)->Many<Dynamic>,|ctx,rhs,linear,initial,start,end,steps,maxeval|exponential(ctx,rhs,Some(linear),Input(Value::Float(0.0)),initial,start,end,steps,maxeval,false));
object_binding!(etdrk4_diag,"num_analysis_ode_etdrk4_diag","crate::kernels::ode_ext::etdrk4_diag","Use real exponential Runge–Kutta weights for an independent linear rate per state component. The diagonal has inverse-time units. The retained tape encodes the nonlinear part of the RHS. The trajectory includes the initial and each accepted fixed-step state, with each component restored to its physical units. The evaluation budget bounds the four tape evaluations per step.",false,(rhs:Handle<OdeTape>,diagonal:Input,initial:Input,start:Input,end:Input,steps:usize,maxeval:usize)->Many<Dynamic>,|ctx,rhs,diagonal,initial,start,end,steps,maxeval|exponential(ctx,rhs,None,diagonal,initial,start,end,steps,maxeval,true));

#[allow(
    clippy::needless_pass_by_value,
    clippy::too_many_arguments,
    reason = "The typed native ABI transfers decoded arguments by value; this shared dispatcher mirrors its explicit source signature."
)]
#[allow(
    deprecated,
    reason = "Expose the actual legacy native entry as requested, alongside its modern replacement."
)]
fn imex(
    ctx: &CallCtx<'_>,
    explicit: Handle<OdeTape>,
    implicit: Handle<OdeTape>,
    initial: Input,
    start: Input,
    end: Input,
    steps: usize,
    name: String,
    maxeval: usize,
    with: bool,
) -> Result<Many<Dynamic>, OpError> {
    let ex = explicit.with(ctx, |m| Ok(m.clone()))?;
    let im = implicit.with(ctx, |m| Ok(m.clone()))?;
    same_frame(&ex.frame, &im.frame)?;
    let frame = &ex.frame;
    let y = frame.state(&initial)?;
    if steps == 0 || steps.checked_mul(y.len() + 2).is_none_or(|n| n > 200000) {
        return Err(api::error("IMEX output/work bound exceeded"));
    }
    let scheme = match name.as_str() {
        | "Euler" => e::ImexScheme::Euler,
        | "Ars222" => e::ImexScheme::Ars222,
        | "Ars343" => e::ImexScheme::Ars343,
        | _ => return Err(api::error("unknown native IMEX scheme")),
    };
    let budget = Budget::new(maxeval)?;
    let fex = |t, y: &[f64], out: &mut [f64]| ex.rate(&budget, t, y, out);
    let fim = |t, y: &[f64], out: &mut [f64]| im.rate(&budget, t, y, out);
    let jac = |t, y: &[f64]| im.jacobian(&budget, t, y);
    let t0 = frame.time(&start)?;
    let t1 = frame.time(&end)?;
    let result = if with {
        e::imex_rk_with(
            fex,
            fim,
            t0,
            t1,
            &y,
            &e::ImexOptions {
                jac_im: Some(&jac),
                steps,
                scheme,
            },
        )
    } else {
        e::imex_rk(fex, fim, Some(&jac), t0, t1, &y, steps, scheme)
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
object_binding!(imex_rk,"num_analysis_ode_imex_rk","crate::kernels::ode_ext::imex_rk","Actual native imex_rk on explicit/implicit physical RHS tapes sharing one declared metric, with all3native schemes, exact implicit tape Jacobian, checked steps and shared evaluation budget.",false,(explicit:Handle<OdeTape>,implicit:Handle<OdeTape>,initial:Input,start:Input,end:Input,steps:usize,name:String,maxeval:usize)->Many<Dynamic>,|ctx,explicit,implicit,initial,start,end,steps,name,maxeval|imex(ctx,explicit,implicit,initial,start,end,steps,name,maxeval,false));
object_binding!(imex_rk_with,"num_analysis_ode_imex_rk_with","crate::kernels::ode_ext::imex_rk_with","Actual native imex_rk_with on explicit/implicit physical RHS tapes sharing one declared metric, with all3native schemes, exact implicit tape Jacobian, checked steps and shared evaluation budget.",false,(explicit:Handle<OdeTape>,implicit:Handle<OdeTape>,initial:Input,start:Input,end:Input,steps:usize,name:String,maxeval:usize)->Many<Dynamic>,|ctx,explicit,implicit,initial,start,end,steps,name,maxeval|imex(ctx,explicit,implicit,initial,start,end,steps,name,maxeval,true));
object_binding!(symplectic_split,"num_analysis_ode_symplectic_split","crate::kernels::ode_ext::symplectic_split","Actual native symplectic split with pure-position acceleration; state contains physical positions then velocities. SplitMethod is explicit Verlet/TripleJump/Yoshida6/Yoshida8 and triple-jump order is bounded2..10.",false,(acceleration:Handle<Acceleration>,initial:Input,start:Input,end:Input,steps:usize,name:String,degree:usize,maxeval:usize)->Dynamic,|ctx,acceleration,initial,start,end,steps,name,degree,maxeval|acceleration.with(ctx,|a|{let frame=a.frame()?;let state=frame.state(&initial)?;let n=state.len()/2;if steps==0||steps>100000{return Err(api::error("split step budget must be1..100000"));}let method=match name.as_str(){"Verlet"=>e::SplitMethod::Verlet,"TripleJump"=>{if !matches!(degree,2|4|6|8|10){return Err(api::error("triple jump order must be2/4/6/8/10"));}e::SplitMethod::TripleJump(degree)},"Yoshida6"=>e::SplitMethod::Yoshida6,"Yoshida8"=>e::SplitMethod::Yoshida8,_=>return Err(api::error("unknown native SplitMethod"))};let budget=Budget::new(maxeval)?;let acceleration=|q:&[f64]|{let mut out=vec![0.0;n];a.evaluate(&budget,q,&mut out);out};let result=e::symplectic_split(acceleration,&state[..n],&state[n..],frame.time(&start)?,frame.time(&end)?,steps,method);budget.finish()?;let(mut q,p)=result.map_err(api::error)?;q.extend(p);Ok(Dynamic(frame.restore(&q)?))}));

#[allow(
    clippy::needless_pass_by_value,
    clippy::too_many_arguments,
    reason = "The typed native ABI transfers decoded arguments by value; this shared dispatcher mirrors its explicit source signature."
)]
#[allow(
    deprecated,
    reason = "Expose the actual legacy native entry as requested, alongside its modern replacement."
)]
fn variational(
    ctx: &CallCtx<'_>,
    mass: Input,
    mass_reference: Input,
    length_reference: Input,
    time_reference: Input,
    potential: Input,
    params: Input,
    positions: Input,
    momenta: Input,
    start: Input,
    end: Input,
    steps: usize,
    stages: usize,
    maxeval: usize,
    with: bool,
) -> Result<Dynamic, OpError> {
    use crate::kernel::signature::DimExpr as D;
    let (m, md) = api::reference(&mass_reference.0)?;
    let (l, ld) = api::reference(&length_reference.0)?;
    let (t, td) = api::reference(&time_reference.0)?;
    api::same_dimension(&md, &D::base("Mass"))?;
    api::same_dimension(&ld, &D::base("Length"))?;
    api::same_dimension(&td, &D::base("Time"))?;
    api::same_dimension(&api::dimension(&mass.0)?, &md)?;
    let mut matrix = <Mat as Arg>::decode(&mass.0, 0)?;
    let n = matrix.rows;
    if n == 0
        || n > 16
        || matrix.cols != n
        || steps == 0
        || steps.checked_mul(n * n).is_none_or(|n| n > 100000)
        || !(2..=8).contains(&stages)
    {
        return Err(api::error(
            "variational dimension/stage/step work bounds exceeded",
        ));
    }
    for v in &mut matrix.data {
        *v = api::scale_product(*v, &[], &[m])?;
    }
    let pscale = api::scale_product(m, &[l], &[t])?;
    let edim = D::Mul(vec![
        md.clone(),
        ld.clone(),
        ld.clone(),
        D::Inverse(Box::new(td.clone())),
        D::Inverse(Box::new(td.clone())),
    ]);
    let pdim = D::Mul(vec![md, ld.clone(), D::Inverse(Box::new(td.clone()))]);
    let escale = api::scale_product(m, &[l, l], &[t, t])?;
    let q = crate::rules::numerical_abi::native_analysis_newton_cg::parameters(
        &positions,
        &vec![l; n],
        &vec![ld.clone(); n],
    )?;
    let p = crate::rules::numerical_abi::native_analysis_newton_cg::parameters(
        &momenta,
        &vec![pscale; n],
        &vec![pdim.clone(); n],
    )?;
    let oracle = VectorOracle::new(
        ctx,
        &potential,
        &params,
        &vec![l; n],
        &vec![ld.clone(); n],
        &[escale],
        &[edim],
    )?;
    let frame = OdeFrame {
        time_scale: t,
        time_dimension: td,
        state_scales: [vec![l; n], vec![pscale; n]].concat(),
        state_dimensions: [vec![ld; n], vec![pdim; n]].concat(),
    };
    let t0 = frame.time(&start)?;
    let t1 = frame.time(&end)?;
    let budget = Budget::new(maxeval)?;
    let grad = |q: &[f64]| {
        match oracle.gradient(&budget, 0, q) {
            | Ok(v) => v,
            | Err(e) => {
                budget.fail(e);
                vec![f64::NAN; n]
            },
        }
    };
    let result = if with {
        e::variational_lobatto_with(
            &matrix,
            grad,
            &q,
            &p,
            t0,
            t1,
            &e::VariationalOptions { steps, stages },
        )
    } else {
        e::variational_lobatto(&matrix, grad, &q, &p, t0, t1, steps, stages)
    };
    budget.finish()?;
    let (mut q, p) = result.map_err(api::error)?;
    q.extend(p);
    Ok(Dynamic(frame.restore(&q)?))
}
object_binding!(variational_lobatto,"num_analysis_ode_variational_lobatto","crate::kernels::ode_ext::variational_lobatto","Actual native Galerkin-Lobatto variational integrator with physical mass matrix, position and momentum. A native tape reverse gradient supplies the potential force; explicit positive M/L/T references preserve energy/force units and stage/work budgets.",false,(mass:Input,mass_reference:Input,length_reference:Input,time_reference:Input,potential:Input,params:Input,positions:Input,momenta:Input,start:Input,end:Input,steps:usize,stages:usize,maxeval:usize)->Dynamic,|ctx,mass,mass_reference,length_reference,time_reference,potential,params,positions,momenta,start,end,steps,stages,maxeval|variational(ctx,mass,mass_reference,length_reference,time_reference,potential,params,positions,momenta,start,end,steps,stages,maxeval,false));
object_binding!(variational_lobatto_with,"num_analysis_ode_variational_lobatto_with","crate::kernels::ode_ext::variational_lobatto_with","Actual native Galerkin-Lobatto variational integrator with physical mass matrix, position and momentum. A native tape reverse gradient supplies the potential force; explicit positive M/L/T references preserve energy/force units and stage/work budgets.",false,(mass:Input,mass_reference:Input,length_reference:Input,time_reference:Input,potential:Input,params:Input,positions:Input,momenta:Input,start:Input,end:Input,steps:usize,stages:usize,maxeval:usize)->Dynamic,|ctx,mass,mass_reference,length_reference,time_reference,potential,params,positions,momenta,start,end,steps,stages,maxeval|variational(ctx,mass,mass_reference,length_reference,time_reference,potential,params,positions,momenta,start,end,steps,stages,maxeval,true));
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &variational_lobatto::BINDING,
            &variational_lobatto_with::BINDING,
            &exp_matrix::BINDING,
            &exprb32::BINDING,
            &etdrk4::BINDING,
            &etdrk4_diag::BINDING,
            &imex_rk::BINDING,
            &imex_rk_with::BINDING,
            &symplectic_split::BINDING,
        ],
        "numerical_vector",
    )
}
