//! Concrete native optimizer instantiations on measured immutable tapes.
//!
//! The solvers come from `basin` (see [`crate::kernels::optimize`]).
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
use basin::CostFunction;
use basin::Gradient;
use ndarray::Array1;
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
    fn charge(&self) -> Result<(), String> {
        let mut current = self.remaining.load(Ordering::Relaxed);
        loop {
            let Some(next) = current.checked_sub(1) else {
                return Err("native optimization oracle budget exhausted".to_string());
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

    fn cost_of(
        &self,
        x: &[f64],
    ) -> Result<f64, String> {
        self.charge()?;
        let mut output = vec![f64::NAN; self.tape.outputs.len()];
        self.tape.eval(x, &self.params, &mut output);
        let value = output[self.output];
        if !value.is_finite() {
            return Err("objective tape result is not finite".to_string());
        }
        Ok(value)
    }

    fn gradient_of(
        &self,
        x: &[f64],
    ) -> Result<Vec<f64>, String> {
        self.charge()?;
        let gradient = self.tape.gradient_reverse(self.output, x, &self.params);
        if gradient.len() != self.tape.inputs || gradient.iter().any(|value| !value.is_finite()) {
            return Err("objective tape gradient is malformed or nonfinite".to_string());
        }
        Ok(gradient)
    }
}
impl CostFunction for TapeCost {
    type Error = String;
    type Output = f64;
    type Param = Vec<f64>;

    fn cost(
        &self,
        x: &Vec<f64>,
    ) -> Result<f64, String> {
        self.cost_of(x)
    }
}
impl Gradient for TapeCost {
    type Gradient = Vec<f64>;

    fn gradient(
        &self,
        x: &Vec<f64>,
    ) -> Result<Vec<f64>, String> {
        self.gradient_of(x)
    }
}
/// The finished run a binding retains.
struct Finished {
    result: crate::kernels::optimize::MinimizeResult,
}
trait RetainedState {
    fn analyze(&self) -> String;
    fn print(&self);
    fn point(&self) -> Option<Vec<f64>>;
    fn cost(&self) -> f64;
    fn iterations(&self) -> u64;
    fn termination(&self) -> String;
}
impl RetainedState for Finished {
    fn analyze(&self) -> String {
        ResultAnalyzer::analyze_convergence(&self.result)
    }

    fn print(&self) {
        ResultAnalyzer::print_optimization_result(&self.result);
    }

    fn point(&self) -> Option<Vec<f64>> {
        Some(self.result.x.clone())
    }

    fn cost(&self) -> f64 {
        self.result.cost
    }

    fn iterations(&self) -> u64 {
        self.result.iterations
    }

    fn termination(&self) -> String {
        self.result.termination.clone()
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
    let tape = Arc::new(frame.tape);
    let bounds = match kind {
        | 3 => {
            let upper = parameters(&upper, &frame.input_scales, &frame.input_dimensions)?;
            if x.iter().zip(&upper).any(|(lo, hi)| lo >= hi) {
                return Err(api::error("PSO requires ordered physical bounds"));
            }
            Some((x.clone(), upper))
        },
        | _ => None,
    };
    let cost_for_value = TapeCost {
        tape: Arc::clone(&tape),
        output,
        params: params.clone(),
        remaining: Arc::new(AtomicUsize::new(budget)),
    };
    let cost_for_gradient = cost_for_value.clone();
    let problem = optimize::FnProblem::try_new(
        move |point: &[f64]| cost_for_value.cost_of(point),
        Some(move |point: &[f64]| cost_for_gradient.gradient_of(point)),
        x.len(),
        bounds.clone(),
    );
    let result = match kind {
        | 0 => EquationOptimizer::solve_with_gradient_descent(&problem, &x, &cfg),
        | 1 => EquationOptimizer::solve_with_bfgs(&problem, &x, &cfg),
        | 2 => EquationOptimizer::auto_solve_conjugate_gradient(&problem, &x, &cfg),
        | 3 => {
            EquationOptimizer::solve_with_pso(&problem, bounds.clone().unwrap_or_default(), &cfg)
        },
        | 4 => {
            EquationOptimizer::auto_solve(
                &problem,
                &x,
                optimize::SolverKind::GradientDescent { step: 0.1 },
                &cfg,
            )
        },
        | _ => return Err(api::error("unknown native optimizer instantiation")),
    }
    .map_err(api::error)?;
    Ok(NewObject::new(
        ResultObject {
            state: Box::new(Finished { result }),
            input_scales: frame.input_scales,
            input_dimensions: frame.input_dimensions,
            output_scale: frame.output_scales[output],
            output_dimension: frame.output_dimensions[output].clone(),
        },
        vec![],
    ))
}

object_binding!(optimizer_solve_with_gradient_descent,"num_analysis_solve_with_gradient_descent_tape","crate::kernels::optimize::EquationOptimizer::solve_with_gradient_descent","Actual native basin solve_with_gradient_descent instantiated with an owned native tape cost and exact reverse gradient. Inputs and parameters retain independent reference domains; target_cost is normalized objective target,not a gradient tolerance. Oracle budget covers both cost and gradient calls. PSO uses initial/upper as bounds and the native randomness semantics. Other methods ignore upper. Generic auto_solve explicitly instantiates the native steepest-descent solver.",false,(tape:Input,output:usize,initial:Input,upper:Input,params:Input,target_cost:f64,iterations:u64,budget:usize)->NewObject<ResultObject>,|ctx:&CallCtx<'_>,tape:Input,output:usize,initial:Input,upper:Input,params:Input,target_cost:f64,iterations:u64,budget:usize|execute(ctx,0,tape,output,initial,upper,params,target_cost,iterations,budget));

object_binding!(optimizer_solve_with_bfgs,"num_analysis_solve_with_bfgs_tape","crate::kernels::optimize::EquationOptimizer::solve_with_bfgs","Actual native basin solve_with_bfgs instantiated with an owned native tape cost and exact reverse gradient. Inputs and parameters retain independent reference domains; target_cost is normalized objective target,not a gradient tolerance. Oracle budget covers both cost and gradient calls. PSO uses initial/upper as bounds and the native randomness semantics. Other methods ignore upper. Generic auto_solve explicitly instantiates the native steepest-descent solver.",false,(tape:Input,output:usize,initial:Input,upper:Input,params:Input,target_cost:f64,iterations:u64,budget:usize)->NewObject<ResultObject>,|ctx:&CallCtx<'_>,tape:Input,output:usize,initial:Input,upper:Input,params:Input,target_cost:f64,iterations:u64,budget:usize|execute(ctx,1,tape,output,initial,upper,params,target_cost,iterations,budget));

object_binding!(optimizer_auto_solve_conjugate_gradient,"num_analysis_auto_solve_conjugate_gradient_tape","crate::kernels::optimize::EquationOptimizer::auto_solve_conjugate_gradient","Actual native basin auto_solve_conjugate_gradient instantiated with an owned native tape cost and exact reverse gradient. Inputs and parameters retain independent reference domains; target_cost is normalized objective target,not a gradient tolerance. Oracle budget covers both cost and gradient calls. PSO uses initial/upper as bounds and the native randomness semantics. Other methods ignore upper. Generic auto_solve explicitly instantiates the native steepest-descent solver.",false,(tape:Input,output:usize,initial:Input,upper:Input,params:Input,target_cost:f64,iterations:u64,budget:usize)->NewObject<ResultObject>,|ctx:&CallCtx<'_>,tape:Input,output:usize,initial:Input,upper:Input,params:Input,target_cost:f64,iterations:u64,budget:usize|execute(ctx,2,tape,output,initial,upper,params,target_cost,iterations,budget));

object_binding!(optimizer_solve_with_pso,"num_analysis_solve_with_pso_tape","crate::kernels::optimize::EquationOptimizer::solve_with_pso","Actual native basin solve_with_pso instantiated with an owned native tape cost and exact reverse gradient. Inputs and parameters retain independent reference domains; target_cost is normalized objective target,not a gradient tolerance. Oracle budget covers both cost and gradient calls. PSO uses initial/upper as bounds and the native randomness semantics. Other methods ignore upper. Generic auto_solve explicitly instantiates the native steepest-descent solver.",false,(tape:Input,output:usize,initial:Input,upper:Input,params:Input,target_cost:f64,iterations:u64,budget:usize)->NewObject<ResultObject>,|ctx:&CallCtx<'_>,tape:Input,output:usize,initial:Input,upper:Input,params:Input,target_cost:f64,iterations:u64,budget:usize|execute(ctx,3,tape,output,initial,upper,params,target_cost,iterations,budget));

object_binding!(optimizer_auto_solve,"num_analysis_auto_solve_tape","crate::kernels::optimize::EquationOptimizer::auto_solve","Actual native basin auto_solve instantiated with an owned native tape cost and exact reverse gradient. Inputs and parameters retain independent reference domains; target_cost is normalized objective target,not a gradient tolerance. Oracle budget covers both cost and gradient calls. PSO uses initial/upper as bounds and the native randomness semantics. Other methods ignore upper. Generic auto_solve explicitly instantiates the native steepest-descent solver.",false,(tape:Input,output:usize,initial:Input,upper:Input,params:Input,target_cost:f64,iterations:u64,budget:usize)->NewObject<ResultObject>,|ctx:&CallCtx<'_>,tape:Input,output:usize,initial:Input,upper:Input,params:Input,target_cost:f64,iterations:u64,budget:usize|execute(ctx,4,tape,output,initial,upper,params,target_cost,iterations,budget));

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
        let width = self.native.x.first().map_or(0, Vec::len);
        if self.native.x.len() > 262144
            || width > 128
            || self.native.x.len() != self.native.y.len()
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
 let x=x.into_iter().map(|row|row.into_iter().map(|v|api::scale_product(v,&[],&[xs])).collect::<Result<Vec<_>,_>>()).collect::<Result<Vec<_>,_>>()?;let y=y.into_iter().map(|v|api::scale_product(v,&[],&[ys])).collect::<Result<Vec<_>,_>>()?;let native=optimize::LinearRegression::new(x,y).map_err(api::error)?;Ok(NewObject::new(RegressionObject{native,xscale:xs,yscale:ys,dimensions:vec![dx.clone(),dy.clone()]},vec![dx,dy]))
});
fn regression_parameters(
    v: &Input,
    r: &RegressionObject,
) -> Result<Array1<f64>, OpError> {
    let slope = crate::kernel::signature::DimExpr::Mul(vec![
        r.dimensions[1].clone(),
        crate::kernel::signature::DimExpr::Inverse(Box::new(r.dimensions[0].clone())),
    ]);
    let width = r.native.x.first().map_or(0, Vec::len);
    let mut dimensions = vec![slope; width + 1];
    dimensions[0] = r.dimensions[1].clone();
    let mut scales = vec![api::scale_product(r.yscale, &[], &[r.xscale])?; width + 1];
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

object_binding!(regression_cost,"num_analysis_regression_cost","","Evaluate the retained native regression CostFunction with measured intercept/slope coefficients; returns physical half mean-squared residual units.",false,(model:Handle<RegressionObject>,coefficients:Input)->Dynamic,|ctx:&CallCtx<'_>,model:Handle<RegressionObject>,coefficients:Input|model.with(ctx,|r|{let beta=regression_parameters(&coefficients,r)?;let value=r.native.value_and_gradient(beta.as_slice().unwrap_or(&[])).map_err(api::error)?.0;Ok(Dynamic(regression_loss_value(value,r)?))}));

object_binding!(
    regression_value_and_gradient, "num_analysis_regression_value_and_gradient",
    "crate::kernels::optimize::LinearRegression::value_and_gradient",
    "Return the half mean-squared residual and its coefficient gradient using the same residual pass. If predictor units are D and response units are E, coefficients use E then E/D, value uses E squared, and the gradient uses E then E times D. The gradient preserves the intercept-first coefficient order.",
    false,
    (model: Handle<RegressionObject>, coefficients: Input) -> Dynamic,
    |ctx, model, coefficients| model.with(ctx, |model| {
        let beta = regression_parameters(&coefficients, model)?;
        let (value, gradient) = model.native.value_and_gradient(beta.as_slice().unwrap_or(&[])).map_err(api::error)?;
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
