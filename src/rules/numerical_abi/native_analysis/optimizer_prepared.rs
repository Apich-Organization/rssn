//! Prepared objective functions and scoped scalar root searches.

use super::RegressionObject;
use super::regression_parameters;
use crate::backend::tape::Inst;
use crate::backend::tape::Tape;
use crate::graph::RuleError;
use crate::graph::dim::DimRule;
use crate::graph::rule::Installer;
use crate::kernel::CallCtx;
use crate::kernel::OpError;
use crate::kernel::Value;
use crate::kernel::signature::DimExpr;
use crate::kernel::signature::TypeExpr;
use crate::kernel::value::CallbackRef;
use crate::kernel::value::TensorValue;
use crate::kernels::optimize::FnProblem;
use crate::kernels::optimize::MinimizeOptions;
use crate::kernels::optimize::MinimizeResult;
use crate::kernels::optimize::Rastrigin;
use crate::kernels::optimize::RootOptions;
use crate::kernels::optimize::RootReport;
use crate::kernels::optimize::RootSolver;
use crate::kernels::optimize::Rosenbrock;
use crate::kernels::optimize::SolverKind;
use crate::kernels::optimize::Sphere;
use crate::kernels::optimize::{
    self,
};
use crate::rules::numerical_abi::codecs::Arg;
use crate::rules::numerical_abi::codecs::Output;
use crate::rules::numerical_abi::native_backend_extended::TapeSnapshot;
use crate::rules::numerical_abi::native_backend_extended::tape_snapshot;
use crate::rules::numerical_abi::native_object_api::Handle;
use crate::rules::numerical_abi::native_object_api::Input;
use crate::rules::numerical_abi::native_object_api::NativeArg;
use crate::rules::numerical_abi::native_object_api::NativeOutput;
use crate::rules::numerical_abi::native_object_api::NativeType;
use crate::rules::numerical_abi::native_object_api::NewObject;
use crate::rules::numerical_abi::native_object_api::error;
use crate::rules::numerical_abi::native_object_api::object_binding;
use crate::rules::numerical_abi::native_object_api::{
    self as api,
};
use std::cell::Cell;
use std::cell::RefCell;
use std::mem::size_of;
use std::sync::Arc;
use std::sync::atomic::AtomicBool;
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering;

const MAX_DIMENSION: usize = 128;
const MAX_WORK: usize = 33_554_432;
const MAX_CALLS: usize = 65_536;

#[derive(Clone)]
struct Bounds {
    lower: Vec<f64>,
    upper: Vec<f64>,
}

impl Bounds {
    fn validate(
        &self,
        dimension: usize,
    ) -> Result<(), OpError> {
        if self.lower.len() != dimension
            || self.upper.len() != dimension
            || self
                .lower
                .iter()
                .zip(&self.upper)
                .any(|(lo, hi)| lo >= hi || !(hi - lo).is_finite())
        {
            return Err(error(
                "box bounds must be finite, increasing and match the objective dimension",
            ));
        }
        Ok(())
    }
}

#[derive(Clone)]
struct Problem {
    function: FnProblem,
    dimension: usize,
    work: usize,
    bytes: usize,
    bounds: Option<Bounds>,
}

impl Problem {
    const fn new(
        function: FnProblem,
        dimension: usize,
    ) -> Self {
        Self {
            function,
            dimension,
            work: dimension,
            bytes: size_of::<Self>(),
            bounds: None,
        }
    }

    fn check_point(
        &self,
        point: &[f64],
    ) -> Result<(), OpError> {
        if point.len() == self.dimension {
            Ok(())
        } else {
            Err(error("point length must match the objective dimension"))
        }
    }
}

impl NativeType for Problem {
    const KIND: &'static str = "rssn/native/analysis/optimizer/problem/v1";
    const MAX_BYTES: usize = 4 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        check_dimension(self.dimension, 1)?;
        if self.retained_bytes() > Self::MAX_BYTES || self.work > MAX_WORK {
            return Err(OpError::BudgetExceeded);
        }
        Ok(())
    }

    fn retained_bytes(&self) -> usize {
        let bounds = self.bounds.as_ref().map_or(0, |bounds| {
            bounds
                .lower
                .capacity()
                .saturating_add(bounds.upper.capacity())
                .saturating_mul(size_of::<f64>())
        });
        self.bytes.saturating_add(bounds)
    }
}

fn check_dimension(
    dimension: usize,
    minimum: usize,
) -> Result<(), OpError> {
    if (minimum..=MAX_DIMENSION).contains(&dimension) {
        Ok(())
    } else {
        Err(error("optimizer dimension is outside the supported range"))
    }
}

struct NormalizedTape(TapeSnapshot);
impl NativeArg for NormalizedTape {
    fn ty() -> TypeExpr {
        TypeExpr::Opaque("rssn/native/backend/tape/v1".into())
    }

    fn decode(
        value: &Value,
        _: usize,
        ctx: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        let snapshot = tape_snapshot(ctx, value)?;
        for dimension in snapshot
            .input_dimensions
            .iter()
            .chain(&snapshot.parameter_dimensions)
            .chain(&snapshot.output_dimensions)
        {
            api::same_dimension(dimension, &DimExpr::Dimless)?;
        }
        if snapshot
            .input_scales
            .iter()
            .chain(&snapshot.parameter_scales)
            .chain(&snapshot.output_scales)
            .any(|scale| scale.to_bits() != 1.0_f64.to_bits())
        {
            return Err(error(
                "prepared objectives require unit-scale normalized tapes",
            ));
        }
        Ok(Self(snapshot))
    }
}

#[derive(Clone)]
struct TapeOracle {
    tape: Arc<Tape>,
    parameters: Arc<[f64]>,
    output: usize,
}

impl TapeOracle {
    fn cost(
        &self,
        point: &[f64],
    ) -> f64 {
        let outputs = self.tape.outputs.len();
        if outputs == 1 {
            let mut output = [f64::NAN];
            self.tape.eval(point, &self.parameters, &mut output);
            let [value] = output;
            return value;
        }
        let mut output = crate::kernels::reuse::global_pool().scratch_f64(outputs);
        output.resize(outputs, f64::NAN);
        self.tape.eval(point, &self.parameters, &mut output);
        output.get(self.output).copied().unwrap_or(f64::NAN)
    }

    fn gradient(
        &self,
        point: &[f64],
    ) -> Vec<f64> {
        self.tape
            .gradient_reverse(self.output, point, &self.parameters)
    }
}

fn prepare_tape(
    tape: NormalizedTape,
    output: usize,
    parameters: Vec<f64>,
    analytic: bool,
    checked: bool,
) -> Result<NewObject<Problem>, OpError> {
    let tape = tape.0.tape;
    let dimension = tape.inputs;
    check_dimension(dimension, 1)?;
    if parameters.len() != tape.params || output >= tape.outputs.len() {
        return Err(error(
            "optimizer output or parameter list does not match the tape",
        ));
    }
    let work = tape
        .insts
        .len()
        .saturating_add(tape.args.len())
        .saturating_add(dimension);
    let bytes = size_of::<Problem>()
        .saturating_add(tape.insts.capacity().saturating_mul(size_of::<Inst>()))
        .saturating_add(
            tape.args
                .capacity()
                .saturating_add(tape.outputs.capacity())
                .saturating_mul(size_of::<usize>()),
        )
        .saturating_add(parameters.capacity().saturating_mul(size_of::<f64>()));
    let cost = TapeOracle {
        tape: Arc::new(tape),
        parameters: parameters.into(),
        output,
    };
    let gradient = cost.clone();
    let function = if checked {
        FnProblem::try_new(
            move |point| {
                let value = cost.cost(point);
                if value.is_finite() {
                    Ok(value)
                } else {
                    Err("tape objective is nonfinite".into())
                }
            },
            analytic.then_some(move |point: &[f64]| {
                let values = gradient.gradient(point);
                if values.len() == dimension && values.iter().all(|value| value.is_finite()) {
                    Ok(values)
                } else {
                    Err("tape gradient is malformed or nonfinite".into())
                }
            }),
            dimension,
            None,
        )
    } else {
        FnProblem::new(
            move |point| cost.cost(point),
            analytic.then_some(move |point: &[f64]| gradient.gradient(point)),
            dimension,
            None,
        )
    };
    Ok(NewObject::normalized(Problem {
        function,
        dimension,
        bytes,
        work,
        bounds: None,
    }))
}

object_binding!(
    problem_new, "num_optimizer_problem_new", "crate::kernels::optimize::FnProblem::new",
    "Prepares an objective from an immutable unit-scale tape, an output slot, dimensionless parameters, and an optional analytic tape gradient.",
    true, (tape: NormalizedTape, output: usize, parameters: Vec<f64>, analytic: bool) -> NewObject<Problem>,
    |_, tape, output, parameters, analytic| prepare_tape(tape, output, parameters, analytic, false)
);
object_binding!(
    problem_try_new, "num_optimizer_problem_try_new", "crate::kernels::optimize::FnProblem::try_new",
    "Prepares a fallible objective from an immutable unit-scale tape. Nonfinite objective and gradient evaluations return domain errors.",
    true, (tape: NormalizedTape, output: usize, parameters: Vec<f64>, analytic: bool) -> NewObject<Problem>,
    |_, tape, output, parameters, analytic| prepare_tape(tape, output, parameters, analytic, true)
);
object_binding!(
    has_gradient, "num_optimizer_has_gradient", "crate::kernels::optimize::FnProblem::has_gradient",
    "Reports whether a prepared problem carries an analytic gradient.",
    true, (problem: Handle<Problem>) -> bool,
    |ctx, problem| problem.with(ctx, |problem| Ok(problem.function.has_gradient()))
);
object_binding!(
    with_bounds, "num_optimizer_with_bounds", "crate::kernels::optimize::FnProblem::with_bounds",
    "Returns a new prepared problem with finite strictly ordered bounds, sharing its objective.",
    true, (problem: Handle<Problem>, lower: Vec<f64>, upper: Vec<f64>) -> NewObject<Problem>,
    |ctx, problem, lower, upper| problem.with(ctx, |problem| {
        let bounds = Bounds { lower, upper };
        bounds.validate(problem.dimension)?;
        let function = problem.function.with_bounds(bounds.lower.clone(), bounds.upper.clone());
        Ok(NewObject::normalized(Problem { function, bounds: Some(bounds), ..problem.clone() }))
    })
);

struct Evaluation(f64, Vec<f64>);
impl NativeOutput for Evaluation {
    fn ty() -> TypeExpr {
        TypeExpr::Record(vec![
            ("value".into(), TypeExpr::Float),
            ("gradient".into(), api::tensor_type(1, TypeExpr::Float)),
        ])
    }

    fn encode(
        self,
        _: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        Ok(api::record(vec![
            ("value", Output::encode(self.0)?),
            ("gradient", Output::encode(self.1)?),
        ]))
    }
}
object_binding!(
    value_gradient, "num_optimizer_value_and_gradient", "crate::kernels::optimize::FnProblem::value_and_gradient",
    "Evaluates a prepared objective and its analytic gradient at one normalized point.",
    true, (problem: Handle<Problem>, point: Vec<f64>) -> Evaluation,
    |ctx, problem, point| problem.with(ctx, |problem| {
        problem.check_point(&point)?;
        let (value, gradient) = problem.function.value_and_gradient(&point).map_err(error)?;
        Ok(Evaluation(value, gradient))
    })
);

macro_rules! analytic_model {
    ($model:ident, $factory:expr, $minimum:literal,
     $prepare:ident, $prepare_name:literal, $prepare_path:literal, $prepare_doc:literal,
     $cost:ident, $cost_name:literal, $cost_path:literal, $cost_doc:literal) => {
        object_binding!(
            $prepare, $prepare_name, $prepare_path, $prepare_doc,
            true, (dimension: usize) -> NewObject<Problem>,
            |_, dimension| {
                check_dimension(dimension, $minimum)?;
                let factory: fn(usize) -> FnProblem = $factory;
                Ok(NewObject::normalized(Problem::new(factory(dimension), dimension)))
            }
        );
        object_binding!(
            $cost, $cost_name, $cost_path, $cost_doc,
            true, (point: Vec<f64>) -> f64,
            |_, point| {
                check_dimension(point.len(), $minimum)?;
                Ok($model::cost(&point))
            }
        );
    };
}

analytic_model!(
    Rosenbrock,
    Rosenbrock::problem,
    2,
    rosenbrock,
    "num_optimizer_rosenbrock",
    "crate::kernels::optimize::Rosenbrock::problem",
    "Builds a normalized Rosenbrock objective with its analytic gradient.",
    rosenbrock_cost,
    "num_optimizer_rosenbrock_cost",
    "crate::kernels::optimize::Rosenbrock::cost",
    "Evaluates the Rosenbrock polynomial at a point with at least two coordinates."
);
analytic_model!(
    Sphere,
    Sphere::problem,
    1,
    sphere,
    "num_optimizer_sphere",
    "crate::kernels::optimize::Sphere::problem",
    "Builds a normalized sphere objective with its analytic gradient.",
    sphere_cost,
    "num_optimizer_sphere_cost",
    "crate::kernels::optimize::Sphere::cost",
    "Evaluates the sum of squared coordinates."
);
analytic_model!(
    Rastrigin,
    |dimension| Rastrigin { dimension }.problem(),
    1,
    rastrigin,
    "num_optimizer_rastrigin",
    "crate::kernels::optimize::Rastrigin::problem",
    "Builds a normalized Rastrigin objective with its analytic gradient.",
    rastrigin_cost,
    "num_optimizer_rastrigin_cost",
    "crate::kernels::optimize::Rastrigin::cost",
    "Evaluates the Rastrigin function at a normalized point."
);

object_binding!(
    regression_size, "num_optimizer_regression_coefficient_count", "crate::kernels::optimize::LinearRegression::coefficient_count",
    "Returns the coefficient count of the retained regression, including the intercept.",
    true, (model: Handle<RegressionObject>) -> usize,
    |ctx, model| model.with(ctx, |model| Ok(model.native.coefficient_count()))
);

struct ResponseReference(DimExpr);
impl NativeArg for ResponseReference {
    fn ty() -> TypeExpr {
        TypeExpr::Float.with_dim(DimExpr::Var("Response".into()))
    }

    fn decode(
        value: &Value,
        _: usize,
        _: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        let (_, dimension) = api::reference(value)?;
        Ok(Self(dimension))
    }
}

struct Residuals(TensorValue);
impl NativeOutput for Residuals {
    fn ty() -> TypeExpr {
        api::tensor_type(1, ResponseReference::ty())
    }

    fn encode(
        self,
        _: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        Ok(Value::Tensor(self.0))
    }
}

object_binding!(
    regression_residuals, "num_optimizer_regression_residuals", "crate::kernels::optimize::LinearRegression::residuals",
    "Returns fitted minus observed responses in coherent response units. The positive response reference attests the output dimension.",
    false, (model: Handle<RegressionObject>, coefficients: Input, response_reference: ResponseReference) -> Residuals,
    |ctx, model, coefficients, response_reference| model.with(ctx, |model| {
        let dimension = model.dimensions.get(1).ok_or_else(|| error("model lacks response units"))?;
        api::same_dimension(dimension, &response_reference.0)?;
        let beta = regression_parameters(&coefficients, model)?;
        let residuals = model.native.residuals(beta.as_slice().unwrap_or(&[]));
        let data = residuals.into_iter().map(|value| {
            api::quantity(api::scale_product(value, &[model.yscale], &[])?, dimension.clone())
        }).collect::<Result<Vec<_>, _>>()?;
        Ok(Residuals(TensorValue { shape: vec![data.len()], elem: TypeExpr::Float.with_dim(dimension.clone()), data }))
    })
);

object_binding!(
    regression_problem, "num_optimizer_regression_problem", "crate::kernels::optimize::LinearRegression::problem",
    "Prepares the regression objective in its stored normalized coordinates and response scale, with an intercept followed by slopes.",
    true, (model: Handle<RegressionObject>) -> NewObject<Problem>,
    |ctx, model| model.with(ctx, |model| {
        let dimension = model.native.coefficient_count();
        check_dimension(dimension, 2)?;
        let observations = model.native.y.len();
        let work = observations.saturating_mul(dimension);
        let bytes = size_of::<Problem>().saturating_add(work.saturating_mul(size_of::<f64>()))
            .saturating_add(observations.saturating_mul(size_of::<Vec<f64>>()));
        Ok(NewObject::normalized(Problem {
            function: model.native.problem(), dimension, work, bytes, bounds: None,
        }))
    })
);

struct Method(SolverKind);
impl Arg for Method {
    fn ty() -> TypeExpr {
        TypeExpr::String
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        const NAMES: [&str; SolverKind::ALL.len()] = [
            "gradient_descent",
            "bfgs",
            "lbfgs",
            "lbfgsb",
            "nonlinear_cg",
            "nelder_mead",
            "nelder_mead_projected",
            "bobyqa",
            "newuoa",
            "lincoa",
            "cobyla",
            "direct",
            "mads",
            "slsqp",
            "projected_gradient_descent",
            "tnc",
            "cma_es",
            "bounded_cma_es",
            "de",
            "global_best_pso",
            "random_search",
            "solis_wets",
            "ssga",
            "gbnm",
            "brent",
            "golden_section",
        ];
        let name = <String as Arg>::decode(value, index)?;
        NAMES
            .iter()
            .zip(SolverKind::ALL)
            .find_map(|(label, kind)| (*label == name).then_some(Self(kind)))
            .ok_or_else(|| error(format!("unknown optimizer solver {name:?}")))
    }
}

object_binding!(
    needs_bounds, "num_optimizer_needs_bounds", "crate::kernels::optimize::SolverKind::needs_bounds",
    "Reports whether the named solver requires an explicit finite box.",
    true, (method: Method) -> bool,
    |_, method| Ok(method.0.needs_bounds())
);
object_binding!(
    needs_gradient, "num_optimizer_needs_gradient", "crate::kernels::optimize::SolverKind::needs_gradient",
    "Reports whether the named solver uses the objective's analytic gradient.",
    true, (method: Method) -> bool,
    |_, method| Ok(method.0.needs_gradient())
);

struct Options {
    native: MinimizeOptions,
    budget: usize,
}
impl Arg for Options {
    fn ty() -> TypeExpr {
        TypeExpr::Record(vec![
            ("iterations".into(), TypeExpr::Int),
            ("evaluations".into(), TypeExpr::Int),
            ("tolerance".into(), TypeExpr::Float),
            ("gradient_tolerance".into(), TypeExpr::Float),
            ("seed".into(), TypeExpr::Int),
        ])
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        let iterations = <u64 as Arg>::decode(api::field(value, "iterations")?, index)?;
        let budget = <usize as Arg>::decode(api::field(value, "evaluations")?, index)?;
        let tolerance = <f64 as Arg>::decode(api::field(value, "tolerance")?, index)?;
        let gradient_tolerance =
            <f64 as Arg>::decode(api::field(value, "gradient_tolerance")?, index)?;
        let seed = <u64 as Arg>::decode(api::field(value, "seed")?, index)?;
        if !(1..=4_096).contains(&iterations) || !(1..=MAX_CALLS).contains(&budget) {
            return Err(OpError::BudgetExceeded);
        }
        if gradient_tolerance <= 0.0 {
            return Err(error("gradient tolerance must be positive"));
        }
        Ok(Self {
            native: MinimizeOptions {
                max_iter: iterations,
                tolerance,
                gradient_tolerance,
                seed,
                bounds: None,
            },
            budget,
        })
    }
}

struct Minimum(MinimizeResult);
impl NativeOutput for Minimum {
    fn ty() -> TypeExpr {
        TypeExpr::Record(vec![
            ("point".into(), api::tensor_type(1, TypeExpr::Float)),
            ("cost".into(), TypeExpr::Float),
            ("iterations".into(), TypeExpr::Int),
            ("cost_evals".into(), TypeExpr::Int),
            ("converged".into(), TypeExpr::Bool),
            ("termination".into(), TypeExpr::String),
        ])
    }

    fn encode(
        self,
        _: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let result = self.0;
        Ok(api::record(vec![
            ("point", Output::encode(result.x)?),
            ("cost", Output::encode(result.cost)?),
            ("iterations", Value::int(result.iterations)),
            ("cost_evals", Value::int(result.cost_evals)),
            ("converged", Value::Bool(result.converged)),
            ("termination", Value::String(result.termination)),
        ]))
    }
}

struct WorkBudget {
    remaining: AtomicUsize,
    exhausted: AtomicBool,
}
impl WorkBudget {
    fn charge(&self) -> Result<(), String> {
        crate::cancellation::checkpoint().map_err(|error| error.to_string())?;
        let reserved = self
            .remaining
            .try_update(Ordering::Relaxed, Ordering::Relaxed, |n| n.checked_sub(1));
        if reserved.is_err() {
            self.exhausted.store(true, Ordering::Relaxed);
            Err("optimizer oracle evaluation budget exhausted".into())
        } else {
            Ok(())
        }
    }
}

#[allow(
    clippy::single_call_fn,
    reason = "Keeps checked budget setup separate from adaptive solver dispatch"
)]
fn solve(
    problem: &Problem,
    start: &[f64],
    method: SolverKind,
    mut options: Options,
) -> Result<Minimum, OpError> {
    problem.check_point(start)?;
    if problem.work.saturating_mul(options.budget) > MAX_WORK {
        return Err(OpError::BudgetExceeded);
    }
    if let Some(bounds) = &problem.bounds {
        if start
            .iter()
            .zip(bounds.lower.iter().zip(&bounds.upper))
            .any(|(point, (lo, hi))| point < lo || point > hi)
        {
            return Err(error("initial point must lie in the prepared box"));
        }
        options.native.bounds = Some((bounds.lower.clone(), bounds.upper.clone()));
    } else if matches!(method, SolverKind::Brent | SolverKind::GoldenSection)
        && start.iter().any(|point| !(point.abs() * 11.0).is_finite())
    {
        return Err(error(
            "the scalar initial point cannot form a finite bracket",
        ));
    }
    let budget = Arc::new(WorkBudget {
        remaining: AtomicUsize::new(options.budget),
        exhausted: AtomicBool::new(false),
    });
    let guard = Arc::clone(&budget);
    let dimension = problem.dimension;
    let model = problem.function.with_oracle_guard(move |point| {
        guard.charge()?;
        if point.len() != dimension || point.iter().any(|value| !value.is_finite()) {
            return Err("optimizer produced an invalid objective point".into());
        }
        Ok(())
    });
    let result = optimize::minimize(&model, start, method, &options.native);
    crate::cancellation::checkpoint()?;
    if budget.exhausted.load(Ordering::Relaxed) {
        return Err(OpError::BudgetExceeded);
    }
    result.map(Minimum).map_err(error)
}

object_binding!(
    minimize, "num_optimizer_minimize", "crate::kernels::optimize::minimize",
    "Runs a named solver against a prepared normalized objective with explicit iteration, evaluation and seed controls.",
    true, (problem: Handle<Problem>, initial: Vec<f64>, method: Method, options: Options) -> Minimum,
    |ctx, problem, initial, method, options| problem.with(ctx, |problem| solve(problem, &initial, method.0, options))
);

struct RootMethod(RootSolver);
impl Arg for RootMethod {
    fn ty() -> TypeExpr {
        TypeExpr::String
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        let name = <String as Arg>::decode(value, index)?;
        match name.as_str() {
            | "bisection" => Ok(Self(RootSolver::Bisection)),
            | "brent" => Ok(Self(RootSolver::Brent)),
            | "secant" => Ok(Self(RootSolver::Secant)),
            | "toms748" => Ok(Self(RootSolver::Toms748)),
            | "newton" => Ok(Self(RootSolver::Newton)),
            | "halley" => Ok(Self(RootSolver::Halley)),
            | _ => Err(error(format!("unknown scalar root solver {name:?}"))),
        }
    }
}

struct ScalarFunction(CallbackRef);
impl NativeArg for ScalarFunction {
    fn ty() -> TypeExpr {
        TypeExpr::Function(vec![TypeExpr::Float], Box::new(TypeExpr::Float))
    }

    fn decode(
        value: &Value,
        _: usize,
        ctx: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        let Value::Function(callback) = value else {
            return Err(error("scalar root search requires a typed callback"));
        };
        let signature = ctx.callback_signature(*callback)?;
        if signature.params != [TypeExpr::Float] || signature.ret != TypeExpr::Float {
            return Err(error(
                "scalar root callbacks require dimensionless Float arguments and results",
            ));
        }
        Ok(Self(*callback))
    }
}

struct ScopedRoot<'a, 'b> {
    context: &'a CallCtx<'b>,
    budget: Cell<usize>,
    failure: RefCell<Option<OpError>>,
}
impl ScopedRoot<'_, '_> {
    fn call(
        &self,
        function: &ScalarFunction,
        point: f64,
    ) -> f64 {
        if self.failure.borrow().is_some() {
            return f64::NAN;
        }
        let invoke = || {
            crate::cancellation::checkpoint()?;
            let budget = self
                .budget
                .get()
                .checked_sub(1)
                .ok_or(OpError::BudgetExceeded)?;
            self.budget.set(budget);
            if !point.is_finite() {
                return Err(error("root solver attempted a nonfinite callback input"));
            }
            let result = self
                .context
                .invoke_callback(function.0, &[Value::Float(point)])?;
            api::same_dimension(&api::dimension(&result)?, &DimExpr::Dimless)?;
            <f64 as Arg>::decode(&result, 0)
        };
        match invoke() {
            | Ok(value) => value,
            | Err(error) => {
                *self.failure.borrow_mut() = Some(error);
                f64::NAN
            },
        }
    }
}

struct Root(RootReport);
impl NativeOutput for Root {
    fn ty() -> TypeExpr {
        TypeExpr::Record(vec![
            ("root".into(), TypeExpr::Float),
            ("residual".into(), TypeExpr::Float),
            ("iterations".into(), TypeExpr::Int),
            ("function_evals".into(), TypeExpr::Int),
            ("converged".into(), TypeExpr::Bool),
        ])
    }

    fn encode(
        self,
        _: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        Ok(api::record(vec![
            ("root", Output::encode(self.0.root)?),
            ("residual", Output::encode(self.0.residual)?),
            ("iterations", Value::int(self.0.iterations)),
            ("function_evals", Value::int(self.0.function_evals)),
            ("converged", Value::Bool(self.0.converged)),
        ]))
    }
}

struct RootConfig {
    solver: RootSolver,
    lower: f64,
    upper: f64,
    tolerance: f64,
    iterations: u64,
    evaluations: usize,
}
impl Arg for RootConfig {
    fn ty() -> TypeExpr {
        TypeExpr::Record(vec![
            ("solver".into(), TypeExpr::String),
            ("lower".into(), TypeExpr::Float),
            ("upper".into(), TypeExpr::Float),
            ("tolerance".into(), TypeExpr::Float),
            ("iterations".into(), TypeExpr::Int),
            ("evaluations".into(), TypeExpr::Int),
        ])
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        let solver = <RootMethod as Arg>::decode(api::field(value, "solver")?, index)?.0;
        let lower = <f64 as Arg>::decode(api::field(value, "lower")?, index)?;
        let upper = <f64 as Arg>::decode(api::field(value, "upper")?, index)?;
        let tolerance = <f64 as Arg>::decode(api::field(value, "tolerance")?, index)?;
        let iterations = <u64 as Arg>::decode(api::field(value, "iterations")?, index)?;
        let evaluations = <usize as Arg>::decode(api::field(value, "evaluations")?, index)?;
        if lower >= upper || !(upper - lower).is_finite() || tolerance <= 0.0 {
            return Err(error(
                "root search requires an ordered finite bracket and positive tolerance",
            ));
        }
        if !(1..=4_096).contains(&iterations) || !(1..=MAX_CALLS).contains(&evaluations) {
            return Err(OpError::BudgetExceeded);
        }
        Ok(Self {
            solver,
            lower,
            upper,
            tolerance,
            iterations,
            evaluations,
        })
    }
}

fn find_root(
    context: &CallCtx<'_>,
    function: &ScalarFunction,
    derivative: Option<&ScalarFunction>,
    config: &RootConfig,
) -> Result<Root, OpError> {
    let scoped = ScopedRoot {
        context,
        budget: Cell::new(config.evaluations),
        failure: RefCell::new(None),
    };
    let result = optimize::find_root(
        |point| scoped.call(function, point),
        derivative.map(|derivative| |point| scoped.call(derivative, point)),
        (config.lower, config.upper),
        config.solver,
        &RootOptions {
            max_iter: config.iterations,
            tolerance: config.tolerance,
        },
    );
    if let Some(error) = scoped.failure.into_inner() {
        return Err(error);
    }
    result.map(Root).map_err(error)
}

object_binding!(
    root, "num_optimizer_find_root", "crate::kernels::optimize::find_root",
    "Runs a bounded bracketed root search with a scoped normalized callback. Newton and Halley estimate a missing derivative.",
    true, (function: ScalarFunction, config: RootConfig) -> Root,
    |ctx, function, config| find_root(ctx, &function, None, &config)
);
object_binding!(
    root_derivative, "num_optimizer_find_root_with_derivative", "",
    "Runs a bounded bracketed root search with scoped normalized function and derivative callbacks.",
    true, (function: ScalarFunction, derivative: ScalarFunction, config: RootConfig) -> Root,
    |ctx, function, derivative, config| find_root(ctx, &function, Some(&derivative), &config)
);

#[allow(
    clippy::single_call_fn,
    reason = "Keeps objective and solver bindings with their checked adapters"
)]
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    const NAMESPACE: &str = "numerical_vector";
    api::install_with_dimension(
        i,
        &[
            &has_gradient::BINDING,
            &rosenbrock_cost::BINDING,
            &sphere_cost::BINDING,
            &rastrigin_cost::BINDING,
            &regression_size::BINDING,
            &needs_bounds::BINDING,
            &needs_gradient::BINDING,
        ],
        NAMESPACE,
        DimRule::Dimensionless,
    )?;
    api::install(
        i,
        &[
            &problem_new::BINDING,
            &problem_try_new::BINDING,
            &with_bounds::BINDING,
            &value_gradient::BINDING,
            &rosenbrock::BINDING,
            &sphere::BINDING,
            &rastrigin::BINDING,
            &regression_residuals::BINDING,
            &regression_problem::BINDING,
            &minimize::BINDING,
            &root::BINDING,
            &root_derivative::BINDING,
        ],
        NAMESPACE,
    )
}
