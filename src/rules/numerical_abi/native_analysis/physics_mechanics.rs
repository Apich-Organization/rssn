//! Native mechanical ODE trait and symplectic integrator instantiations.
use super::Arc;
use super::Budget;
use super::CallCtx;
use super::Dynamic;
use super::Handle;
use super::Input;
use super::Installer;
use super::Many;
use super::NativeType;
use super::NewObject;
use super::OdeFrame;
use super::OpError;
use super::RuleError;
use super::Value;
use super::api;
use super::e;
use super::i;
use super::object_binding;
use super::output;
use crate::kernel::signature::DimExpr;
use crate::rules::numerical_abi::native_analysis_newton_cg::parameters;
use crate::rules::numerical_abi::native_backend_extended::TapeSnapshot;
use crate::rules::numerical_abi::native_backend_extended::tape_snapshot;
use crate::sim::physics_em::MechanicalSystem;
#[derive(Clone)]
pub(in crate::rules) struct Acceleration {
    snapshot: TapeSnapshot,
    params: Vec<f64>,
    time_scale: f64,
    time_dimension: DimExpr,
    ratios: Vec<f64>,
}
impl NativeType for Acceleration {
    const KIND: &'static str = "rssn/native/analysis/mechanics/acceleration/v1";
    const MAX_BYTES: usize = 4 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        if self.ratios.is_empty()
            || self.ratios.len() > 32
            || self.ratios.len() != self.snapshot.input_scales.len()
            || self.ratios.len() != self.snapshot.output_scales.len()
        {
            return Err(api::error(
                "mechanical acceleration needs1..32position inputs and matching outputs",
            ));
        }
        Ok(())
    }
}
impl Acceleration {
    pub(in crate::rules) fn frame(&self) -> Result<OdeFrame, OpError> {
        let mut scales = self.snapshot.input_scales.clone();
        scales.extend(
            self.snapshot
                .input_scales
                .iter()
                .map(|v| api::scale_product(*v, &[], &[self.time_scale]))
                .collect::<Result<Vec<_>, _>>()?,
        );
        let mut dims = self.snapshot.input_dimensions.clone();
        dims.extend(self.snapshot.input_dimensions.iter().map(|d| {
            DimExpr::Mul(vec![
                d.clone(),
                DimExpr::Inverse(Box::new(self.time_dimension.clone())),
            ])
        }));
        Ok(OdeFrame {
            time_scale: self.time_scale,
            time_dimension: self.time_dimension.clone(),
            state_scales: scales,
            state_dimensions: dims,
        })
    }

    pub(in crate::rules) fn evaluate(
        &self,
        budget: &Budget,
        x: &[f64],
        out: &mut [f64],
    ) {
        out.fill(f64::NAN);
        if x.len() != self.ratios.len() || out.len() != x.len() {
            budget.fail(api::error("mechanical acceleration shape mismatch"));
            return;
        }
        if !budget.charge() {
            return;
        }
        self.snapshot.tape.eval(x, &self.params, out);
        for (v, s) in out.iter_mut().zip(&self.ratios) {
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
struct Mechanical {
    data: Acceleration,
    budget: Arc<Budget>,
}
impl MechanicalSystem for Mechanical {
    fn spatial_dim(&self) -> usize {
        self.data.ratios.len()
    }

    fn eval_acceleration(
        &self,
        x: &[f64],
        a: &mut [f64],
    ) {
        self.data.evaluate(&self.budget, x, a);
    }
}
object_binding!(acceleration_new,"num_analysis_mechanics_acceleration_new","","Owned position-only acceleration tape. Positive time reference and each coordinate reference restore acceleration coordinate/time², preventing hidden dependence on velocity or time in the MechanicalSystem trait.",false,(tape:Input,params:Input,time_reference:Input)->NewObject<Acceleration>,|ctx,tape,params,time_reference|{let snapshot=tape_snapshot(ctx,&tape.0)?;let(time_scale,time_dimension)=api::reference(&time_reference.0)?;if snapshot.input_scales.is_empty()||snapshot.input_scales.len()!=snapshot.output_scales.len()||snapshot.input_scales.len()>32{return Err(api::error("acceleration tape shape mismatch"));}let ratios=snapshot.output_scales.iter().zip(&snapshot.input_scales).zip(snapshot.output_dimensions.iter().zip(&snapshot.input_dimensions)).map(|((o,i),(od,id))|{api::same_dimension(od,&DimExpr::Mul(vec![id.clone(),DimExpr::Inverse(Box::new(time_dimension.clone())),DimExpr::Inverse(Box::new(time_dimension.clone()))]))?;api::scale_product(*o,&[time_scale,time_scale],&[*i])}).collect::<Result<Vec<_>,OpError>>()?;let params=parameters(&params,&snapshot.parameter_scales,&snapshot.parameter_dimensions)?;let data=Acceleration{snapshot,params,time_scale,time_dimension,ratios};let dims=data.frame()?.dimensions();Ok(NewObject::new(data,dims))});
object_binding!(spatial_dim,"num_analysis_mechanics_spatial_dim","crate::sim::physics_em::MechanicalSystem::spatial_dim","Actual instantiated native MechanicalSystem spatial dimension.",false,(acceleration:Handle<Acceleration>)->usize,|ctx,acceleration|acceleration.with(ctx,|a|Ok(Mechanical{data:a.clone(),budget:Budget::new(1)?}.spatial_dim())));
object_binding!(eval,"num_analysis_mechanics_eval","crate::sim::physics_em::MechanicalSystem::eval_acceleration","Actual native MechanicalSystem acceleration evaluation with physical coordinate/time² result components.",false,(acceleration:Handle<Acceleration>,positions:Input,maxeval:usize)->Dynamic,|ctx,acceleration,positions,maxeval|acceleration.with(ctx,|a|{let x=parameters(&positions,&a.snapshot.input_scales,&a.snapshot.input_dimensions)?;let system=Mechanical{data:a.clone(),budget:Budget::new(maxeval)?};let mut out=vec![0.0;x.len()];system.eval_acceleration(&x,&mut out);system.budget.finish()?;let data=out.iter().zip(&a.snapshot.input_scales).zip(&a.snapshot.input_dimensions).map(|((v,s),d)|api::quantity(api::scale_product(*v,&[*s],&[a.time_scale,a.time_scale])?,DimExpr::Mul(vec![d.clone(),DimExpr::Inverse(Box::new(a.time_dimension.clone())),DimExpr::Inverse(Box::new(a.time_dimension.clone()))]))).collect::<Result<Vec<_>,_>>()?;Ok(Dynamic(Value::Tensor(crate::kernel::value::TensorValue{shape:vec![data.len()],elem:api::dynamic_type(),data})))}));
object_binding!(semi_implicit,"num_analysis_mechanics_semi_implicit_euler","crate::sim::physics_em::solve_semi_implicit_euler","Actual native mechanical semi-implicit Euler. State contains independently measured positions then corresponding velocities; source tape depends only on position.",false,(acceleration:Handle<Acceleration>,initial:Input,start:Input,end:Input,h:Input,maxeval:usize)->Many<Dynamic>,|ctx,acceleration,initial,start,end,h,maxeval|acceleration.with(ctx,|a|{let frame=a.frame()?;let y=frame.state(&initial)?;let t0=frame.time(&start)?;let t1=frame.time(&end)?;let h=frame.time(&h)?;if h<=0.0||t1<=t0||((t1-t0)/h).ceil()*(y.len()+4)as f64>262000.0{return Err(api::error("mechanical Euler output/work budget exceeded"));}let system=Mechanical{data:a.clone(),budget:Budget::new(maxeval)?};let result=e::solve_semi_implicit_euler(&system,&y,(t0,t1),h);system.budget.finish()?;output(&frame,result.map_err(api::error)?) }));
fn momentum_frame(
    a: &Acceleration,
    mass: &Input,
) -> Result<OdeFrame, OpError> {
    let (m, md) = api::reference(&mass.0)?;
    api::same_dimension(&md, &DimExpr::Base("Mass".into()))?;
    let mut f = a.frame()?;
    let n = a.ratios.len();
    for j in n..f.state_scales.len() {
        f.state_scales[j] = api::scale_product(f.state_scales[j], &[m], &[])?;
        f.state_dimensions[j] = DimExpr::Mul(vec![md.clone(), f.state_dimensions[j].clone()]);
    }
    Ok(f)
}
#[allow(
    clippy::needless_pass_by_value,
    clippy::too_many_arguments,
    reason = "The typed native ABI transfers decoded arguments by value; this shared dispatcher mirrors its explicit source signature."
)]
fn phase(
    ctx: &CallCtx<'_>,
    acceleration: Handle<Acceleration>,
    positions: Input,
    momenta: Input,
    mass: Input,
    h: Input,
    steps: usize,
    order: u32,
    maxeval: usize,
    kind: usize,
) -> Result<Many<Dynamic>, OpError> {
    let a = acceleration.with(ctx, |a| Ok(a.clone()))?;
    let frame = momentum_frame(&a, &mass)?;
    let n = a.ratios.len();
    let mut q = parameters(
        &positions,
        &frame.state_scales[..n],
        &frame.state_dimensions[..n],
    )?;
    let mut p = parameters(
        &momenta,
        &frame.state_scales[n..],
        &frame.state_dimensions[n..],
    )?;
    let h = frame.time(&h)?;
    if h <= 0.0
        || steps == 0
        || steps.checked_mul(2 * n + 4).is_none_or(|n| n > 262000)
        || !matches!(order, 2 | 4)
    {
        return Err(api::error(
            "symplectic trajectory order must be2or4and work bounded",
        ));
    }
    let budget = Budget::new(maxeval)?;
    let acceleration = |q: &[f64]| {
        let mut out = vec![0.0; n];
        a.evaluate(&budget, q, &mut out);
        out
    };
    let rows = match kind {
        | 0 => {
            i::verlet_step(&mut q, &mut p, 1.0, h, &acceleration);
            vec![(q, p)]
        },
        | 1 => {
            i::yoshida4_step(&mut q, &mut p, 1.0, h, &acceleration);
            vec![(q, p)]
        },
        | 2 => i::symplectic_trajectory(&q, &p, 1.0, h, steps, order, acceleration),
        | _ => return Err(api::error("unknown native symplectic operation")),
    };
    budget.finish()?;
    Ok(Many(
        rows.into_iter()
            .map(|(q, p)| {
                let mut state = q;
                state.extend(p);
                frame.restore(&state).map(Dynamic)
            })
            .collect::<Result<_, _>>()?,
    ))
}
object_binding!(verlet_step,"num_analysis_mechanics_verlet_step","crate::sim::integrators::verlet_step","Actual native verlet_step with physical position, momentum, mass and time. The mass itself defines the positive mass reference, so native mass1is a checked normalization and physical momentum is restored with m*coordinate/time. No numerical update is implemented in glue.",false,(acceleration:Handle<Acceleration>,positions:Input,momenta:Input,mass:Input,h:Input,steps:usize,order:u32,maxeval:usize)->Many<Dynamic>,|ctx,acceleration,positions,momenta,mass,h,steps,order,maxeval|phase(ctx,acceleration,positions,momenta,mass,h,steps,order,maxeval,0));
object_binding!(yoshida4_step,"num_analysis_mechanics_yoshida4_step","crate::sim::integrators::yoshida4_step","Actual native yoshida4_step with physical position, momentum, mass and time. The mass itself defines the positive mass reference, so native mass1is a checked normalization and physical momentum is restored with m*coordinate/time. No numerical update is implemented in glue.",false,(acceleration:Handle<Acceleration>,positions:Input,momenta:Input,mass:Input,h:Input,steps:usize,order:u32,maxeval:usize)->Many<Dynamic>,|ctx,acceleration,positions,momenta,mass,h,steps,order,maxeval|phase(ctx,acceleration,positions,momenta,mass,h,steps,order,maxeval,1));
object_binding!(symplectic_trajectory,"num_analysis_mechanics_symplectic_trajectory","crate::sim::integrators::symplectic_trajectory","Actual native symplectic_trajectory with physical position, momentum, mass and time. The mass itself defines the positive mass reference, so native mass1is a checked normalization and physical momentum is restored with m*coordinate/time. No numerical update is implemented in glue.",false,(acceleration:Handle<Acceleration>,positions:Input,momenta:Input,mass:Input,h:Input,steps:usize,order:u32,maxeval:usize)->Many<Dynamic>,|ctx,acceleration,positions,momenta,mass,h,steps,order,maxeval|phase(ctx,acceleration,positions,momenta,mass,h,steps,order,maxeval,2));

struct Linear {
    matrix: crate::kernels::matrix::Matrix<f64>,
    frame: OdeFrame,
}
impl NativeType for Linear {
    const KIND: &'static str = "rssn/native/analysis/mechanics/linear/v1";
    const MAX_BYTES: usize = 4 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        if self.frame.state_scales.len() > 128 {
            return Err(api::error("linear state budget exceeds128"));
        }
        Ok(())
    }
}
impl e::LinearOdeSystem for Linear {
    fn dim(&self) -> usize {
        self.frame.state_scales.len()
    }

    fn get_matrix(&self) -> crate::kernels::matrix::Matrix<f64> {
        self.matrix.clone()
    }
}
object_binding!(linear_new,"num_analysis_mechanics_linear_new","","Declare a constant linear ODE from a normalized native dense Jacobian carrying its physical RHS frame. The resulting coefficient matrix is the caller's explicit constant linear model, not a linearity inference from sampled values.",false,(matrix:Handle<crate::rules::numerical_abi::native_analysis_jacobian::MatrixObject>)->NewObject<Linear>,|ctx,matrix|matrix.with(ctx,|m|{let frame=m.frame.clone().ok_or_else(||api::error("physical linear model requires a Jacobian with RHS reference frame"))?;let crate::rules::numerical_abi::native_analysis_jacobian::MatrixData::Dense(m)=&m.data else{return Err(api::error("native linear Euler requires dense matrix"));};let data=Linear{matrix:crate::kernels::matrix::Matrix::new(m.rows,m.cols,m.data.clone()),frame:frame.clone()};Ok(NewObject::new(data,frame.dimensions()))}));
object_binding!(linear_dim,"num_analysis_mechanics_linear_dim","crate::sim::physics_em::LinearOdeSystem::dim","Actual native LinearOdeSystem dimension.",false,(linear:Handle<Linear>)->usize,|ctx,linear|linear.with(ctx,|l|Ok(e::LinearOdeSystem::dim(l))));
object_binding!(linear_matrix,"num_analysis_mechanics_linear_matrix","crate::sim::physics_em::LinearOdeSystem::get_matrix","Actual native LinearOdeSystem matrix in its declared normalized reference coordinate metric; mixed physical row/column domains remain on the retained model.",false,(linear:Handle<Linear>)->crate::rules::numerical_abi::native_object_api::NewObject<crate::rules::numerical_abi::native_analysis_jacobian::MatrixObject>,|ctx,linear|linear.with(ctx,|l|{let m=e::LinearOdeSystem::get_matrix(l);let n=e::LinearOdeSystem::dim(l);let data=(0..n).flat_map(|i|{let m=&m;(0..n).map(move|j|*m.get(i,j))}).collect();Ok(NewObject::new(crate::rules::numerical_abi::native_analysis_jacobian::MatrixObject{data:crate::rules::numerical_abi::native_analysis_jacobian::MatrixData::Dense(crate::kernels::dense::Mat{rows:n,cols:n,data}),frame:Some(l.frame.clone()),evaluations:0},l.frame.dimensions()))}));
object_binding!(backward_linear,"num_analysis_mechanics_backward_euler_linear","crate::sim::physics_em::solve_backward_euler_linear","Actual native generic backward Euler on the retained constant matrix with physical time/state reconstruction and bounded allocation.",false,(linear:Handle<Linear>,initial:Input,start:Input,end:Input,h:Input)->Many<Dynamic>,|ctx,linear,initial,start,end,h|linear.with(ctx,|l|{let frame=&l.frame;let y=frame.state(&initial)?;let t0=frame.time(&start)?;let t1=frame.time(&end)?;let h=frame.time(&h)?;if h<=0.0||t1<=t0||((t1-t0)/h).ceil()*(y.len()+4)as f64>262000.0{return Err(api::error("linear Euler step/output work budget exceeded"));}output(frame,e::solve_backward_euler_linear(l,&y,(t0,t1),h).map_err(api::error)?) }));
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &acceleration_new::BINDING,
            &spatial_dim::BINDING,
            &eval::BINDING,
            &semi_implicit::BINDING,
            &verlet_step::BINDING,
            &yoshida4_step::BINDING,
            &symplectic_trajectory::BINDING,
            &linear_new::BINDING,
            &linear_dim::BINDING,
            &linear_matrix::BINDING,
            &backward_linear::BINDING,
        ],
        "numerical_vector",
    )
}
