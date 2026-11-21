//! Concrete native argmin generic instantiations on measured immutable tapes.
use crate::backend::tape::Tape;
use crate::graph::RuleError;
use crate::graph::rule::Installer;
use crate::kernel::CallCtx;
use crate::kernel::OpError;
use crate::kernel::Value;
use crate::kernel::value::TensorValue;
use crate::kernels::optimize::EquationOptimizer;
use crate::kernels::optimize::OptimizationConfig;
use crate::kernels::optimize::ProblemType;
use crate::kernels::optimize::ResultAnalyzer;
use crate::kernels::optimize::{
    self,
};
use crate::rules::numerical_abi::codecs::Arg;
use crate::rules::numerical_abi::codecs::Output;
use crate::rules::numerical_abi::native_backend_extended::tape_snapshot;
use crate::rules::numerical_abi::native_object_api::Dynamic;
use crate::rules::numerical_abi::native_object_api::Handle;
use crate::rules::numerical_abi::native_object_api::Input;
use crate::rules::numerical_abi::native_object_api::NativeType;
use crate::rules::numerical_abi::native_object_api::NewObject;
use crate::rules::numerical_abi::native_object_api::object_binding;
use crate::rules::numerical_abi::native_object_api::{
    self as api,
};
use argmin::core::CostFunction;
use argmin::core::Gradient;
use argmin::core::IterState;
use argmin::core::State;
use ndarray::Array1;
use ndarray::Array2;
use std::sync::Arc;
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering;
#[derive(Clone)]
struct TapeCost {
    tape: Arc<Tape>,
    output: usize,
    params: Vec<f64>,
    remaining: Arc<AtomicUsize>,
}
impl TapeCost {
    fn charge(&self) -> Result<(), argmin::core::Error> {
        let mut current = self.remaining.load(Ordering::Relaxed);
        loop {
            let Some(next) = current.checked_sub(1) else {
                return Err(argmin::core::Error::msg(
                    "native optimization oracle budget exhausted",
                ));
            };
            match self.remaining.compare_exchange_weak(
                current,
                next,
                Ordering::Relaxed,
                Ordering::Relaxed,
            ) {
                | Ok(_) => return Ok(()),
                | Err(actual) => current = actual,
            }
        }
    }
}
impl CostFunction for TapeCost {
    type Output = f64;
    type Param = Array1<f64>;

    fn cost(
        &self,
        x: &Self::Param,
    ) -> Result<f64, argmin::core::Error> {
        self.charge()?;
        let x = x
            .as_slice()
            .ok_or_else(|| argmin::core::Error::msg("optimizer point must be contiguous"))?;
        let mut output = vec![f64::NAN; self.tape.outputs.len()];
        self.tape.eval(x, &self.params, &mut output);
        let value = output[self.output];
        if !value.is_finite() {
            return Err(argmin::core::Error::msg(
                "objective tape result is not finite",
            ));
        }
        Ok(value)
    }
}
impl Gradient for TapeCost {
    type Gradient = Array1<f64>;
    type Param = Array1<f64>;

    fn gradient(
        &self,
        x: &Self::Param,
    ) -> Result<Self::Gradient, argmin::core::Error> {
        self.charge()?;
        let x = x
            .as_slice()
            .ok_or_else(|| argmin::core::Error::msg("optimizer point must be contiguous"))?;
        let g = self.tape.gradient_reverse(self.output, x, &self.params);
        if g.len() != self.tape.inputs || g.iter().any(|v| !v.is_finite()) {
            return Err(argmin::core::Error::msg(
                "objective tape gradient is malformed or nonfinite",
            ));
        }
        Ok(Array1::from(g))
    }
}
type VectorState = IterState<Array1<f64>, Array1<f64>, (), (), (), f64>;
trait RetainedState {
    fn analyze(&self) -> String;
    fn print(&self);
    fn point(&self) -> Option<Vec<f64>>;
    fn cost(&self) -> f64;
    fn iterations(&self) -> u64;
    fn termination(&self) -> String;
}
struct Vector<S>(S);
impl<S: State<Param = Array1<f64>, Float = f64>> RetainedState for Vector<S> {
    fn analyze(&self) -> String {
        ResultAnalyzer::analyze_convergence(&self.0)
    }

    fn print(&self) {
        ResultAnalyzer::print_optimization_result(&self.0);
    }

    fn point(&self) -> Option<Vec<f64>> {
        self.0.get_best_param().map(ndarray::ArrayBase::to_vec)
    }

    fn cost(&self) -> f64 {
        self.0.get_best_cost()
    }

    fn iterations(&self) -> u64 {
        self.0.get_iter()
    }

    fn termination(&self) -> String {
        format!("{:?}", self.0.get_termination_status())
    }
}
struct Population(
    argmin::core::PopulationState<argmin::solver::particleswarm::Particle<Array1<f64>, f64>, f64>,
);
impl RetainedState for Population {
    fn analyze(&self) -> String {
        ResultAnalyzer::analyze_convergence(&self.0)
    }

    fn print(&self) {
        let mut projected = VectorState::new();
        projected.best_param = self.0.get_best_param().map(|p| p.position.clone());
        projected.best_cost = self.0.get_best_cost();
        projected.iter = self.0.get_iter();
        projected.counts.clone_from(self.0.get_func_counts());
        ResultAnalyzer::print_optimization_result(&projected);
    }

    fn point(&self) -> Option<Vec<f64>> {
        self.0.get_best_param().map(|p| p.position.to_vec())
    }

    fn cost(&self) -> f64 {
        self.0.get_best_cost()
    }

    fn iterations(&self) -> u64 {
        self.0.get_iter()
    }

    fn termination(&self) -> String {
        format!("{:?}", self.0.get_termination_status())
    }
}
struct ResultObject {
    state: Box<dyn RetainedState>,
    input_scales: Vec<f64>,
    input_dimensions: Vec<crate::kernel::signature::DimExpr>,
    output_scale: f64,
    output_dimension: crate::kernel::signature::DimExpr,
}
impl NativeType for ResultObject {
    const KIND: &'static str = "rssn/native/analysis/argmin/result/v1";
    const MAX_BYTES: usize = 4 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        if self.input_scales.len() > 128 || self.input_scales.len() != self.input_dimensions.len() {
            Err(api::error("optimizer state/reference shape mismatch"))
        } else {
            Ok(())
        }
    }
}
fn parameters(
    v: &Input,
    scales: &[f64],
    dimensions: &[crate::kernel::signature::DimExpr],
) -> Result<Vec<f64>, OpError> {
    super::native_analysis_newton_cg::parameters(v, scales, dimensions)
}
#[allow(
    clippy::needless_pass_by_value,
    clippy::too_many_arguments,
    reason = "The typed native ABI transfers decoded arguments by value; this shared dispatcher mirrors its explicit source signature."
)]
fn execute(
    ctx: &CallCtx<'_>,
    kind: usize,
    tape: Input,
    output: usize,
    initial: Input,
    upper: Input,
    params: Input,
    target: f64,
    iterations: u64,
    budget: usize,
) -> Result<NewObject<ResultObject>, OpError> {
    let frame = tape_snapshot(ctx, &tape.0)?;
    if frame.tape.inputs == 0
        || frame.tape.inputs > 128
        || output >= frame.output_scales.len()
        || iterations > 1_000_000
        || budget == 0
        || budget
            .checked_mul(
                frame
                    .tape
                    .insts
                    .len()
                    .saturating_add(frame.tape.outputs.len())
                    .saturating_add(frame.tape.args.len())
                    .max(1),
            )
            .is_none_or(|n| n > 100_000_000)
    {
        return Err(api::error(
            "native optimizer dimensions or execution budget exceeded",
        ));
    }
    let x = parameters(&initial, &frame.input_scales, &frame.input_dimensions)?;
    let params = parameters(
        &params,
        &frame.parameter_scales,
        &frame.parameter_dimensions,
    )?;
    let cfg = OptimizationConfig {
        max_iters: iterations,
        tolerance: target,
        problem_type: ProblemType::Custom,
        dimension: x.len(),
    };
    let cost = TapeCost {
        tape: Arc::new(frame.tape),
        output,
        params,
        remaining: Arc::new(AtomicUsize::new(budget)),
    };
    let state: Box<dyn RetainedState> = match kind {
        | 0 => {
            Box::new(Vector(
                EquationOptimizer::solve_with_gradient_descent(cost, Array1::from(x), &cfg)
                    .map_err(api::error)?
                    .state,
            ))
        },
        | 1 => {
            Box::new(Vector(
                EquationOptimizer::solve_with_bfgs(cost, Array1::from(x), &cfg)
                    .map_err(api::error)?
                    .state,
            ))
        },
        | 2 => {
            Box::new(Vector(
                EquationOptimizer::auto_solve_conjugate_gradient(cost, Array1::from(x), &cfg)
                    .map_err(api::error)?
                    .state,
            ))
        },
        | 3 => {
            let upper = parameters(&upper, &frame.input_scales, &frame.input_dimensions)?;
            if x.iter().zip(&upper).any(|(lo, hi)| lo >= hi) {
                return Err(api::error("PSO requires ordered physical bounds"));
            }
            Box::new(Population(
                EquationOptimizer::solve_with_pso(
                    cost,
                    (Array1::from(x), Array1::from(upper)),
                    &cfg,
                )
                .map_err(api::error)?
                .state,
            ))
        },
        | 4 => {
            let search = argmin::solver::linesearch::MoreThuenteLineSearch::<
                Array1<f64>,
                Array1<f64>,
                f64,
            >::new();
            let solver = argmin::solver::gradientdescent::SteepestDescent::new(search);
            let result = EquationOptimizer::auto_solve(cost, solver, |state: VectorState| {
                state
                    .param(Array1::from(x))
                    .max_iters(iterations)
                    .target_cost(target)
            })
            .map_err(api::error)?;
            Box::new(Vector(result.state))
        },
        | _ => return Err(api::error("unknown native optimizer instantiation")),
    };
    Ok(NewObject::new(
        ResultObject {
            state,
            input_scales: frame.input_scales,
            input_dimensions: frame.input_dimensions,
            output_scale: frame.output_scales[output],
            output_dimension: frame.output_dimensions[output].clone(),
        },
        vec![],
    ))
}

object_binding!(optimizer_solve_with_gradient_descent,"num_analysis_solve_with_gradient_descent_tape","crate::kernels::optimize::EquationOptimizer::solve_with_gradient_descent","Actual native argmin solve_with_gradient_descent instantiated with an owned native tape cost and exact reverse gradient. Inputs and parameters retain independent reference domains; target_cost is normalized objective target,not a gradient tolerance. Oracle budget covers both cost and gradient calls. PSO uses initial/upper as bounds and the native randomness semantics. Other methods ignore upper. Generic auto_solve explicitly instantiates the native steepest-descent solver.",false,(tape:Input,output:usize,initial:Input,upper:Input,params:Input,target_cost:f64,iterations:u64,budget:usize)->NewObject<ResultObject>,|ctx:&CallCtx<'_>,tape:Input,output:usize,initial:Input,upper:Input,params:Input,target_cost:f64,iterations:u64,budget:usize|execute(ctx,0,tape,output,initial,upper,params,target_cost,iterations,budget));

object_binding!(optimizer_solve_with_bfgs,"num_analysis_solve_with_bfgs_tape","crate::kernels::optimize::EquationOptimizer::solve_with_bfgs","Actual native argmin solve_with_bfgs instantiated with an owned native tape cost and exact reverse gradient. Inputs and parameters retain independent reference domains; target_cost is normalized objective target,not a gradient tolerance. Oracle budget covers both cost and gradient calls. PSO uses initial/upper as bounds and the native randomness semantics. Other methods ignore upper. Generic auto_solve explicitly instantiates the native steepest-descent solver.",false,(tape:Input,output:usize,initial:Input,upper:Input,params:Input,target_cost:f64,iterations:u64,budget:usize)->NewObject<ResultObject>,|ctx:&CallCtx<'_>,tape:Input,output:usize,initial:Input,upper:Input,params:Input,target_cost:f64,iterations:u64,budget:usize|execute(ctx,1,tape,output,initial,upper,params,target_cost,iterations,budget));

object_binding!(optimizer_auto_solve_conjugate_gradient,"num_analysis_auto_solve_conjugate_gradient_tape","crate::kernels::optimize::EquationOptimizer::auto_solve_conjugate_gradient","Actual native argmin auto_solve_conjugate_gradient instantiated with an owned native tape cost and exact reverse gradient. Inputs and parameters retain independent reference domains; target_cost is normalized objective target,not a gradient tolerance. Oracle budget covers both cost and gradient calls. PSO uses initial/upper as bounds and the native randomness semantics. Other methods ignore upper. Generic auto_solve explicitly instantiates the native steepest-descent solver.",false,(tape:Input,output:usize,initial:Input,upper:Input,params:Input,target_cost:f64,iterations:u64,budget:usize)->NewObject<ResultObject>,|ctx:&CallCtx<'_>,tape:Input,output:usize,initial:Input,upper:Input,params:Input,target_cost:f64,iterations:u64,budget:usize|execute(ctx,2,tape,output,initial,upper,params,target_cost,iterations,budget));

object_binding!(optimizer_solve_with_pso,"num_analysis_solve_with_pso_tape","crate::kernels::optimize::EquationOptimizer::solve_with_pso","Actual native argmin solve_with_pso instantiated with an owned native tape cost and exact reverse gradient. Inputs and parameters retain independent reference domains; target_cost is normalized objective target,not a gradient tolerance. Oracle budget covers both cost and gradient calls. PSO uses initial/upper as bounds and the native randomness semantics. Other methods ignore upper. Generic auto_solve explicitly instantiates the native steepest-descent solver.",false,(tape:Input,output:usize,initial:Input,upper:Input,params:Input,target_cost:f64,iterations:u64,budget:usize)->NewObject<ResultObject>,|ctx:&CallCtx<'_>,tape:Input,output:usize,initial:Input,upper:Input,params:Input,target_cost:f64,iterations:u64,budget:usize|execute(ctx,3,tape,output,initial,upper,params,target_cost,iterations,budget));

object_binding!(optimizer_auto_solve,"num_analysis_auto_solve_tape","crate::kernels::optimize::EquationOptimizer::auto_solve","Actual native argmin auto_solve instantiated with an owned native tape cost and exact reverse gradient. Inputs and parameters retain independent reference domains; target_cost is normalized objective target,not a gradient tolerance. Oracle budget covers both cost and gradient calls. PSO uses initial/upper as bounds and the native randomness semantics. Other methods ignore upper. Generic auto_solve explicitly instantiates the native steepest-descent solver.",false,(tape:Input,output:usize,initial:Input,upper:Input,params:Input,target_cost:f64,iterations:u64,budget:usize)->NewObject<ResultObject>,|ctx:&CallCtx<'_>,tape:Input,output:usize,initial:Input,upper:Input,params:Input,target_cost:f64,iterations:u64,budget:usize|execute(ctx,4,tape,output,initial,upper,params,target_cost,iterations,budget));

object_binding!(result_analyze,"num_analysis_optimization_analyze","crate::kernels::optimize::ResultAnalyzer::analyze_convergence","Native analyzer of retained optimizer state; native labels use normalized objective thresholds and are heuristics,not convergence certificates.",true,(result:Handle<ResultObject>)->String,|ctx:&CallCtx<'_>,result:Handle<ResultObject>|result.with(ctx,|r|Ok(r.state.analyze())));
object_binding!(result_print,"num_analysis_optimization_print","crate::kernels::optimize::ResultAnalyzer::print_optimization_result","Print native normalized optimizer-state diagnostics. Particle populations project their actual best position,cost,iteration and evaluation counts to the vector State required by the native printer.",true,(result:Handle<ResultObject>)->(),|ctx:&CallCtx<'_>,result:Handle<ResultObject>|result.with(ctx,|r|{r.state.print();Ok(())}));
object_binding!(result_data,"num_analysis_optimization_result","","Physical best point/objective and actual native termination from retained optimizer state. Missing native best values remain None; mixed coordinate domains are preserved per component.",false,(result:Handle<ResultObject>)->Dynamic,|ctx:&CallCtx<'_>,result:Handle<ResultObject>|result.with(ctx,|r|{
 let point=r.state.point().map(|x|{if x.len()!=r.input_scales.len(){return Err(api::error("optimizer point/reference count mismatch"));}let values=x.into_iter().zip(&r.input_scales).zip(&r.input_dimensions).map(|((v,s),d)|api::quantity(api::scale_product(v,&[*s],&[])?,d.clone())).collect::<Result<Vec<_>,_>>()?;Ok(Value::Tensor(TensorValue{shape:vec![values.len()],elem:api::dynamic_type(),data:values}))}).transpose()?;
 let cost=r.state.cost();let cost=if cost.is_finite(){Some(api::quantity(api::scale_product(cost,&[r.output_scale],&[])?,r.output_dimension.clone())?)}else{None};let option=|v:Option<Value>|Value::Variant(crate::kernel::value::VariantValue{tag:if v.is_some(){"Some"}else{"None"}.into(),payload:v.into_iter().collect()});
 Ok(Dynamic(api::record(vec![("x",option(point)),("value",option(cost)),("iterations",Output::encode(r.state.iterations())?),("termination",Value::String(r.state.termination()))])))
}));
struct RegressionObject {
    native: optimize::LinearRegression,
    xscale: f64,
    yscale: f64,
    dimensions: Vec<crate::kernel::signature::DimExpr>,
}
impl NativeType for RegressionObject {
    const KIND: &'static str = "rssn/native/analysis/argmin/regression/v1";
    const MAX_BYTES: usize = 4 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        if self.native.x.len() > 262144
            || self.native.x.ncols() > 128
            || self.native.x.nrows() != self.native.y.len()
            || self.dimensions.len() != 2
        {
            Err(api::error("native regression model shape/budget mismatch"))
        } else {
            Ok(())
        }
    }
}
object_binding!(regression_new,"num_analysis_linear_regression_new","crate::kernels::optimize::LinearRegression::new","Construct an actual native regression CostFunction with independently measured predictor/response domains and positive reference quantities; coefficients retain intercept E and slope E/D units.",false,(x:Input,y:Input,x_reference:Input,y_reference:Input)->NewObject<RegressionObject>,|_ctx:&CallCtx<'_>,x:Input,y:Input,x_reference:Input,y_reference:Input|{
 let(xs,dx)=api::reference(&x_reference.0)?;let(ys,dy)=api::reference(&y_reference.0)?;api::same_dimension(&dx,&api::dimension(&x.0)?)?;api::same_dimension(&dy,&api::dimension(&y.0)?)?;let x=<Vec<Vec<f64>> as Arg>::decode(&x.0,0)?;let y=<Vec<f64> as Arg>::decode(&y.0,1)?;let n=x.len();let m=x.first().map_or(0,Vec::len);if n==0||m==0||m>128||n.checked_mul(m).is_none_or(|n|n>262144)||n!=y.len()||x.iter().any(|r|r.len()!=m){return Err(api::error("native regression shape/allocation budget mismatch"));}
 let x=Array2::from_shape_vec((n,m),x.into_iter().flatten().map(|v|api::scale_product(v,&[],&[xs])).collect::<Result<_,_>>()?).map_err(api::error)?;let y=Array1::from(y.into_iter().map(|v|api::scale_product(v,&[],&[ys])).collect::<Result<Vec<_>,_>>()?);let native=optimize::LinearRegression::new(x,y).map_err(api::error)?;Ok(NewObject::new(RegressionObject{native,xscale:xs,yscale:ys,dimensions:vec![dx.clone(),dy.clone()]},vec![dx,dy]))
});
fn regression_parameters(
    v: &Input,
    r: &RegressionObject,
) -> Result<Array1<f64>, OpError> {
    let slope = crate::kernel::signature::DimExpr::Mul(vec![
        r.dimensions[1].clone(),
        crate::kernel::signature::DimExpr::Inverse(Box::new(r.dimensions[0].clone())),
    ]);
    let mut dimensions = vec![slope; r.native.x.ncols() + 1];
    dimensions[0] = r.dimensions[1].clone();
    let mut scales = vec![api::scale_product(r.yscale, &[], &[r.xscale])?; r.native.x.ncols() + 1];
    scales[0] = r.yscale;
    Ok(Array1::from(parameters(v, &scales, &dimensions)?))
}
fn regression_loss_value(
    value: f64,
    model: &RegressionObject,
) -> Result<Value, OpError> {
    let result_dimension = crate::kernel::signature::DimExpr::Mul(vec![
        model.dimensions[1].clone(),
        model.dimensions[1].clone(),
    ]);
    api::quantity(
        api::scale_product(value, &[model.yscale, model.yscale], &[])?,
        result_dimension,
    )
}

object_binding!(regression_cost,"num_analysis_regression_cost","","Evaluate the retained native regression CostFunction with measured intercept/slope coefficients; returns physical half mean-squared residual units.",false,(model:Handle<RegressionObject>,coefficients:Input)->Dynamic,|ctx:&CallCtx<'_>,model:Handle<RegressionObject>,coefficients:Input|model.with(ctx,|r|{let beta=regression_parameters(&coefficients,r)?;let value=r.native.cost(&beta).map_err(api::error)?;Ok(Dynamic(regression_loss_value(value,r)?))}));

object_binding!(
    regression_value_and_gradient, "num_analysis_regression_value_and_gradient",
    "crate::kernels::optimize::LinearRegression::value_and_gradient",
    "Return the half mean-squared residual and its coefficient gradient using the same residual pass. If predictor units are D and response units are E, coefficients use E then E/D, value uses E squared, and the gradient uses E then E times D. The gradient preserves the intercept-first coefficient order.",
    false,
    (model: Handle<RegressionObject>, coefficients: Input) -> Dynamic,
    |ctx, model, coefficients| model.with(ctx, |model| {
        let beta = regression_parameters(&coefficients, model)?;
        let (value, gradient) = model.native.value_and_gradient(&beta).map_err(api::error)?;
        let predicted_gradient_dimension = crate::kernel::signature::DimExpr::Mul(vec![
            model.dimensions[1].clone(), model.dimensions[0].clone(),
        ]);
        let data = gradient.iter().enumerate().map(|(i, &entry)| {
            if i == 0 {
                api::quantity(api::scale_product(entry, &[model.yscale], &[])?, model.dimensions[1].clone())
            } else {
                api::quantity(api::scale_product(entry, &[model.yscale, model.xscale], &[])?, predicted_gradient_dimension.clone())
            }
        }).collect::<Result<Vec<_>, _>>()?;
        Ok(Dynamic(api::record(vec![
            ("value", regression_loss_value(value, model)?),
            ("gradient", Value::Tensor(TensorValue {
                shape: vec![data.len()], elem: api::dynamic_type(), data,
            })),
        ])))
    })
);

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &optimizer_solve_with_gradient_descent::BINDING,
            &optimizer_solve_with_bfgs::BINDING,
            &optimizer_auto_solve_conjugate_gradient::BINDING,
            &optimizer_solve_with_pso::BINDING,
            &optimizer_auto_solve::BINDING,
            &result_analyze::BINDING,
            &result_print::BINDING,
            &result_data::BINDING,
            &regression_new::BINDING,
            &regression_cost::BINDING,
            &regression_value_and_gradient::BINDING,
        ],
        "numerical_vector",
    )
}
