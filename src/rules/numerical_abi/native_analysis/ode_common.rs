//! Owned native tape right-hand sides in mixed physical state reference frames.
use crate::backend::tape::Tape;
use crate::graph::RuleError;
use crate::graph::rule::Installer;
use crate::kernel::CallCtx;
use crate::kernel::OpError;
use crate::kernel::Value;
use crate::kernel::signature::DimExpr;
use crate::kernels::dense::Mat;
use crate::rules::numerical_abi::native_analysis_newton_cg::parameters;
use crate::rules::numerical_abi::native_backend_extended::TapeSnapshot;
use crate::rules::numerical_abi::native_backend_extended::tape_snapshot;
use crate::rules::numerical_abi::native_object_api::Input;
use crate::rules::numerical_abi::native_object_api::NativeType;
use crate::rules::numerical_abi::native_object_api::NewObject;
use crate::rules::numerical_abi::native_object_api::object_binding;
use crate::rules::numerical_abi::native_object_api::{
    self as api,
};
use std::sync::Arc;
use std::sync::Mutex;
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering;
pub(in crate::rules) type RateFn = Box<dyn Fn(f64, &[f64], &mut [f64]) + Send + Sync>;
pub(in crate::rules) struct Budget {
    initial: usize,
    remaining: AtomicUsize,
    failure: Mutex<Option<OpError>>,
}
impl Budget {
    // Evaluation limits leave the high bit free. Using the same atomic for
    // failure and charging makes the hot path independent of the error lock.
    const FAILED: usize = 1usize << (usize::BITS - 1);

    pub fn new(n: usize) -> Result<Arc<Self>, OpError> {
        if n == 0 || n > 1_000_000 {
            return Err(api::error(
                "ODE/native callback evaluation budget must be1..1000000",
            ));
        }
        Ok(Arc::new(Self {
            initial: n,
            remaining: AtomicUsize::new(n),
            failure: Mutex::new(None),
        }))
    }

    pub fn fail(
        &self,
        error: OpError,
    ) {
        if let Ok(mut first) = self.failure.lock()
            && first.is_none()
        {
            *first = Some(error);
            self.remaining.fetch_or(Self::FAILED, Ordering::Relaxed);
        }
    }

    pub fn charge(&self) -> bool {
        let mut n = self.remaining.load(Ordering::Relaxed);
        loop {
            if n & Self::FAILED != 0 {
                return false;
            }
            let Some(next) = n.checked_sub(1) else {
                self.fail(OpError::BudgetExceeded);
                return false;
            };
            match self.remaining.compare_exchange_weak(
                n,
                next,
                Ordering::Relaxed,
                Ordering::Relaxed,
            ) {
                | Ok(_) => return true,
                | Err(actual) => n = actual,
            }
        }
    }

    pub fn used(&self) -> usize {
        self.initial
            .saturating_sub(self.remaining.load(Ordering::Relaxed) & !Self::FAILED)
    }

    pub fn finish(&self) -> Result<(), OpError> {
        let first = self
            .failure
            .lock()
            .map_err(|_| api::error("native callback failure state poisoned"))?;
        first.as_ref().map_or(Ok(()), |e| Err(e.clone()))
    }
}
#[derive(Clone)]
pub(in crate::rules) struct OdeFrame {
    pub time_scale: f64,
    pub time_dimension: DimExpr,
    pub state_scales: Vec<f64>,
    pub state_dimensions: Vec<DimExpr>,
}
impl OdeFrame {
    pub fn dimensions(&self) -> Vec<DimExpr> {
        let mut d = vec![self.time_dimension.clone()];
        d.extend(self.state_dimensions.clone());
        d
    }

    pub fn state(
        &self,
        input: &Input,
    ) -> Result<Vec<f64>, OpError> {
        parameters(input, &self.state_scales, &self.state_dimensions)
    }

    pub fn time(
        &self,
        input: &Input,
    ) -> Result<f64, OpError> {
        api::same_dimension(&self.time_dimension, &api::dimension(&input.0)?)?;
        api::scale_product(super::numeric(&input.0, 0)?, &[], &[self.time_scale])
    }

    pub fn restore(
        &self,
        values: &[f64],
    ) -> Result<Value, OpError> {
        if values.len() != self.state_scales.len() {
            return Err(api::error("ODE state/frame size mismatch"));
        }
        let data = values
            .iter()
            .zip(&self.state_scales)
            .zip(&self.state_dimensions)
            .map(|((v, s), d)| api::quantity(api::scale_product(*v, &[*s], &[])?, d.clone()))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Value::Tensor(crate::kernel::value::TensorValue {
            shape: vec![data.len()],
            elem: api::dynamic_type(),
            data,
        }))
    }
}
#[derive(Clone)]
pub(in crate::rules) struct OdeTape {
    pub tape: Arc<Tape>,
    pub params: Vec<f64>,
    pub frame: OdeFrame,
    pub rate_scales: Vec<f64>,
}
impl NativeType for OdeTape {
    const KIND: &'static str = "rssn/native/analysis/ode/rhs/v1";
    const MAX_BYTES: usize = 4 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        let n = self.frame.state_scales.len();
        if n == 0
            || n > 128
            || self.tape.inputs != n + 1
            || self.tape.outputs.len() != n
            || self.params.len() != self.tape.params
            || self.rate_scales.len() != n
            || self.tape.insts.len() > 4096
        {
            return Err(api::error(
                "ODE tape shape or128state/4096instruction budget mismatch",
            ));
        }
        Ok(())
    }
}
impl OdeTape {
    #[allow(
        clippy::needless_pass_by_value,
        reason = "The typed native ABI transfers decoded arguments by value; this shared dispatcher mirrors its explicit source signature."
    )]
    pub fn create(
        snapshot: TapeSnapshot,
        params: Input,
    ) -> Result<Self, OpError> {
        let n = snapshot
            .input_scales
            .len()
            .checked_sub(1)
            .ok_or_else(|| api::error("ODE tape needs time and state inputs"))?;
        if n == 0 || n > 128 || snapshot.output_scales.len() != n {
            return Err(api::error("ODE tape needs n+1 inputs and n outputs"));
        }
        let time = snapshot.input_dimensions[0].clone();
        let state = snapshot.input_dimensions[1..].to_vec();
        let mut rate_scales = Vec::with_capacity(n);
        for i in 0..n {
            let expected = DimExpr::Mul(vec![
                state[i].clone(),
                DimExpr::Inverse(Box::new(time.clone())),
            ]);
            api::same_dimension(&expected, &snapshot.output_dimensions[i])?;
            rate_scales.push(api::scale_product(
                snapshot.output_scales[i],
                &[snapshot.input_scales[0]],
                &[snapshot.input_scales[i + 1]],
            )?);
        }
        let params = parameters(
            &params,
            &snapshot.parameter_scales,
            &snapshot.parameter_dimensions,
        )?;
        Ok(Self {
            tape: Arc::new(snapshot.tape),
            params,
            frame: OdeFrame {
                time_scale: snapshot.input_scales[0],
                time_dimension: time,
                state_scales: snapshot.input_scales[1..].to_vec(),
                state_dimensions: state,
            },
            rate_scales,
        })
    }

    pub fn rate(
        &self,
        budget: &Budget,
        t: f64,
        y: &[f64],
        out: &mut [f64],
    ) {
        out.fill(f64::NAN);
        if y.len() != self.frame.state_scales.len() || out.len() != y.len() {
            budget.fail(api::error("native ODE rate shape mismatch"));
            return;
        }
        if !budget.charge() {
            return;
        }
        let mut x = Vec::with_capacity(y.len() + 1);
        x.push(t);
        x.extend_from_slice(y);
        self.tape.eval(&x, &self.params, out);
        for (v, s) in out.iter_mut().zip(&self.rate_scales) {
            match api::scale_product(*v, &[*s], &[]) {
                | Ok(x) => *v = x,
                | Err(e) => {
                    budget.fail(e);
                    *v = f64::NAN;
                },
            }
        }
    }

    pub fn jacobian(
        &self,
        budget: &Budget,
        t: f64,
        y: &[f64],
    ) -> Mat {
        let n = self.frame.state_scales.len();
        let mut matrix = Mat::zeros(n, n);
        if y.len() != n || !budget.charge() {
            matrix.data.fill(f64::NAN);
            return matrix;
        }
        let mut x = Vec::with_capacity(n + 1);
        x.push(t);
        x.extend_from_slice(y);
        for i in 0..n {
            let g = self.tape.gradient_reverse(i, &x, &self.params);
            if g.len() != n + 1 {
                budget.fail(api::error("native tape Jacobian shape mismatch"));
                matrix.data.fill(f64::NAN);
                return matrix;
            }
            for j in 0..n {
                match api::scale_product(g[j + 1], &[self.rate_scales[i]], &[]) {
                    | Ok(value) => matrix.set(i, j, value),
                    | Err(e) => {
                        budget.fail(e);
                        matrix.set(i, j, f64::NAN);
                    },
                }
            }
        }
        matrix
    }

    pub fn function(
        self: Arc<Self>,
        budget: Arc<Budget>,
    ) -> RateFn {
        Box::new(move |t, y, out| self.rate(&budget, t, y, out))
    }
}
object_binding!(rhs_new,"num_analysis_ode_rhs_new","","Build an owned native tape RHS: inputs(time,state components),outputs corresponding rates. Each state may have independent physical units; tape references define normalized solver coordinates. Parameter quantities retain their own domains.",false,(tape:Input,params:Input)->NewObject<OdeTape>,|ctx:&CallCtx<'_>,tape:Input,params:Input|{let model=OdeTape::create(tape_snapshot(ctx,&tape.0)?,params)?;let dimensions=model.frame.dimensions();Ok(NewObject::new(model,dimensions))});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(i, &[&rhs_new::BINDING], "numerical_vector")
}

/// An immutable native vector tape instantiated in an explicit caller reference frame.
/// Input/output conversions implement only changes of reference coordinates.
#[derive(Clone)]
pub(in crate::rules) struct VectorOracle {
    pub snapshot: TapeSnapshot,
    params: Vec<f64>,
    input_ratios: Vec<f64>,
    output_ratios: Vec<f64>,
}
impl VectorOracle {
    pub fn new(
        ctx: &CallCtx<'_>,
        tape: &Input,
        params: &Input,
        input_scales: &[f64],
        input_dimensions: &[DimExpr],
        output_scales: &[f64],
        output_dimensions: &[DimExpr],
    ) -> Result<Self, OpError> {
        let snapshot = tape_snapshot(ctx, &tape.0)?;
        if snapshot.input_scales.len() != input_scales.len()
            || snapshot.output_scales.len() != output_scales.len()
            || input_scales.len() != input_dimensions.len()
            || output_scales.len() != output_dimensions.len()
        {
            return Err(api::error("native vector tape reference arities disagree"));
        }
        let input_ratios = input_scales
            .iter()
            .zip(&snapshot.input_scales)
            .zip(input_dimensions.iter().zip(&snapshot.input_dimensions))
            .map(|((a, b), (da, db))| {
                api::same_dimension(da, db)?;
                api::scale_product(*a, &[], &[*b])
            })
            .collect::<Result<Vec<_>, OpError>>()?;
        let output_ratios = snapshot
            .output_scales
            .iter()
            .zip(output_scales)
            .zip(snapshot.output_dimensions.iter().zip(output_dimensions))
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

    pub fn evaluate(
        &self,
        budget: &Budget,
        input: &[f64],
    ) -> Result<Vec<f64>, OpError> {
        if input.len() != self.input_ratios.len() {
            return Err(api::error("native vector oracle input shape mismatch"));
        }
        if !budget.charge() {
            budget.finish()?;
            return Err(OpError::BudgetExceeded);
        }
        let input = input
            .iter()
            .zip(&self.input_ratios)
            .map(|(v, s)| api::scale_product(*v, &[*s], &[]))
            .collect::<Result<Vec<_>, _>>()?;
        let mut output = vec![0.0; self.output_ratios.len()];
        self.snapshot.tape.eval(&input, &self.params, &mut output);
        output
            .iter()
            .zip(&self.output_ratios)
            .map(|(v, s)| api::scale_product(*v, &[*s], &[]))
            .collect()
    }

    /// Exact native reverse gradient in the declared normalized caller coordinates.
    pub fn gradient(
        &self,
        budget: &Budget,
        output: usize,
        input: &[f64],
    ) -> Result<Vec<f64>, OpError> {
        if output >= self.output_ratios.len() || input.len() != self.input_ratios.len() {
            return Err(api::error("vector tape gradient shape/index mismatch"));
        }
        if !budget.charge() {
            budget.finish()?;
            return Err(OpError::BudgetExceeded);
        }
        let input = input
            .iter()
            .zip(&self.input_ratios)
            .map(|(v, r)| api::scale_product(*v, &[*r], &[]))
            .collect::<Result<Vec<_>, _>>()?;
        let gradient = self
            .snapshot
            .tape
            .gradient_reverse(output, &input, &self.params);
        gradient
            .iter()
            .zip(&self.input_ratios)
            .map(|(g, r)| api::scale_product(*g, &[self.output_ratios[output], *r], &[]))
            .collect()
    }

    pub fn sample(
        &self,
        budget: &Budget,
        input: &[f64],
    ) -> Vec<f64> {
        match self.evaluate(budget, input) {
            | Ok(v) => v,
            | Err(e) => {
                budget.fail(e);
                vec![f64::NAN; self.output_ratios.len()]
            },
        }
    }
}
pub(in crate::rules) fn reference_frame(
    time_reference: &Input,
    state_references: &crate::rules::numerical_abi::native_object_api::Many<Input>,
) -> Result<OdeFrame, OpError> {
    let (time_scale, time_dimension) = api::reference(&time_reference.0)?;
    if state_references.0.is_empty() || state_references.0.len() > 32 {
        return Err(api::error("ODE state references require1..32components"));
    }
    let refs = state_references
        .0
        .iter()
        .map(|r| api::reference(&r.0))
        .collect::<Result<Vec<_>, _>>()?;
    let (state_scales, state_dimensions) = refs.into_iter().unzip();
    Ok(OdeFrame {
        time_scale,
        time_dimension,
        state_scales,
        state_dimensions,
    })
}
