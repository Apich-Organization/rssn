// Measured linear state-space discretisation and simulation.
use crate::kernels::control::StateSpace;

struct DiscreteStateSpace {
    a: dense::Mat,
    b: dense::Mat,
    c: dense::Mat,
    d: dense::Mat,
    dt: f64,
}
record!(DiscreteStateSpace {
    a: dense::Mat,
    b: dense::Mat,
    c: dense::Mat,
    d: dense::Mat,
    dt: f64,
});

impl Arg for DiscreteStateSpace {
    fn ty() -> TypeExpr {
        <Self as Output>::ty()
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        use super::native_object_api::field;

        let Value::Record(record) = value else {
            return Err(OpError::Type {
                expected: Box::new(<Self as Arg>::ty()),
                got: Box::new(value.type_expr()),
                arg: index,
            });
        };
        if record.fields.len() != 5 {
            return Err(OpError::Domain(
                "a discrete state-space model needs fields a, b, c, d and dt".into(),
            ));
        }
        Ok(Self {
            a: dense::Mat::decode(field(value, "a")?, index)?,
            b: dense::Mat::decode(field(value, "b")?, index)?,
            c: dense::Mat::decode(field(value, "c")?, index)?,
            d: dense::Mat::decode(field(value, "d")?, index)?,
            dt: f64::decode(field(value, "dt")?, index)?,
        })
    }
}

struct StateSpaceShape {
    states: usize,
    inputs: usize,
    outputs: usize,
}

impl StateSpaceShape {
    fn check_zoh_workspace(&self) -> Result<(), OpError> {
        let size = self.states.saturating_add(self.inputs);
        if self.states != 0 && size.saturating_mul(size) >= super::native_object_api::VALUE_LIMIT {
            return Err(OpError::Domain(
                "augmented state-space matrix exceeds the native size limit".into(),
            ));
        }
        Ok(())
    }

    fn check_model_result(&self) -> Result<(), OpError> {
        let Self {
            states: n,
            inputs: m,
            outputs: p,
        } = *self;
        // One record, four tensor headers and two nodes for the measured dt.
        // A is dimensionless; each entry of B, C and D can need a unit witness.
        let measured = n
            .saturating_mul(m)
            .saturating_add(p.saturating_mul(n.saturating_add(m)));
        let nodes = n
            .saturating_mul(n)
            .saturating_add(measured.saturating_mul(2))
            .saturating_add(7);
        if nodes > super::native_object_api::VALUE_LIMIT {
            return Err(OpError::Domain(
                "discrete state-space model exceeds the native value limit".into(),
            ));
        }
        Ok(())
    }
}

fn state_space_rows(matrix: dense::Mat) -> Vec<Vec<f64>> {
    let mut data = matrix.data.into_iter();
    (0..matrix.rows)
        .map(|_| data.by_ref().take(matrix.cols).collect())
        .collect()
}

fn state_space_result(
    rows: Vec<Vec<f64>>,
    cols: usize,
) -> dense::Mat {
    dense::Mat {
        rows: rows.len(),
        cols,
        data: rows.into_iter().flatten().collect(),
    }
}

fn state_space_error(message: String) -> OpError {
    match crate::cancellation::checkpoint() {
        | Ok(()) => OpError::Domain(message),
        | Err(cancelled) => cancelled.into(),
    }
}

fn state_space_from_matrices(
    a: dense::Mat,
    b: dense::Mat,
    c: dense::Mat,
    d: dense::Mat,
    dt: Option<f64>,
) -> Result<(StateSpace, StateSpaceShape), OpError> {
    let n = a.rows;
    let m = b.cols;
    let p = c.rows;
    if a.cols != n || b.rows != n || c.cols != n || d.rows != p || d.cols != m {
        return Err(OpError::Domain(
            "state-space matrices need shapes [n,n], [n,m], [p,n] and [p,m]".into(),
        ));
    }
    if n == 0 && p == 0 && m != 0 {
        return Err(OpError::Domain(
            "a state-space model with inputs needs states or outputs".into(),
        ));
    }
    crate::cancellation::checkpoint()?;
    let [a, b, c, d] = [a, b, c, d].map(state_space_rows);
    let model = match dt {
        | Some(dt) => StateSpace::discrete(a, b, c, d, dt),
        | None => StateSpace::continuous(a, b, c, d),
    }
    .map_err(state_space_error)?;
    Ok((
        model,
        StateSpaceShape {
            states: n,
            inputs: m,
            outputs: p,
        },
    ))
}

fn state_space_trajectory(
    model: &StateSpace,
    shape: &StateSpaceShape,
    x0: &[f64],
    inputs: dense::Mat,
    times: Option<&[f64]>,
) -> Result<dense::Mat, OpError> {
    if inputs.cols != shape.inputs {
        return Err(OpError::Domain(
            "the input matrix width must match the state-space model".into(),
        ));
    }
    if inputs.rows.saturating_mul(shape.outputs.saturating_add(1))
        >= super::native_object_api::VALUE_LIMIT
    {
        return Err(OpError::Domain(
            "state-space trajectory exceeds the native tensor size limit".into(),
        ));
    }
    if times.is_some_and(|times| times.len() > 1) {
        shape.check_zoh_workspace()?;
    }
    let inputs = state_space_rows(inputs);
    let output = match times {
        | Some(times) => model.simulate_continuous(x0, &inputs, times),
        | None => model.simulate(x0, &inputs),
    }
    .map_err(state_space_error)?;
    crate::cancellation::checkpoint()?;
    Ok(state_space_result(output, shape.outputs))
}

binding!(
    control_c2d_zoh,
    "num_control_c2d_zoh",
    "numerical_matrix",
    "Discretise a measured continuous state-space model under a zero-order hold. Accepts positive sample time and returns a portable record with a, b, c, d and dt. Matrix a has inverse-time units before discretisation, and b gains one factor of Time.",
    "crate::kernels::control::state_space::StateSpace::c2d_zoh",
    DimensionPolicy::StateSpaceZoh,
    (a: dense::Mat, b: dense::Mat, c: dense::Mat, d: dense::Mat, dt: f64) -> Result<super::DiscreteStateSpace, OpError>,
    |a: dense::Mat, b: dense::Mat, c: dense::Mat, d: dense::Mat, dt: f64| {
        let (model, shape) = super::state_space_from_matrices(a, b, c, d, None)?;
        shape.check_zoh_workspace()?;
        shape.check_model_result()?;
        let crate::kernels::control::StateSpace { a, b, c, d, .. } = model.c2d_zoh(dt).map_err(super::state_space_error)?;
        crate::cancellation::checkpoint()?;
        Ok(super::DiscreteStateSpace {
            a: super::state_space_result(a, shape.states),
            b: super::state_space_result(b, shape.inputs),
            c: super::state_space_result(c, shape.states),
            d: super::state_space_result(d, shape.inputs),
            dt,
        })
    }
);

binding!(
    control_simulate_discrete,
    "num_control_simulate_discrete",
    "numerical_matrix",
    "Simulate a portable discrete state-space model from an initial state and one input row per sample. Return one measured output row per input row, including the initial output and direct feedthrough.",
    "crate::kernels::control::state_space::StateSpace::simulate",
    DimensionPolicy::StateSpaceDiscrete,
    (model: super::DiscreteStateSpace, x0: Vec<f64>, inputs: dense::Mat) -> Result<dense::Mat, OpError>,
    |model: super::DiscreteStateSpace, x0: Vec<f64>, inputs: dense::Mat| {
        let super::DiscreteStateSpace { a, b, c, d, dt } = model;
        let (model, shape) = super::state_space_from_matrices(a, b, c, d, Some(dt))?;
        super::state_space_trajectory(&model, &shape, &x0, inputs, None)
    }
);

binding!(
    control_simulate_continuous,
    "num_control_simulate_continuous",
    "numerical_matrix",
    "Simulate a measured continuous state-space model on a strictly increasing time grid. Hold each input row until the next time and return one measured output row per grid point.",
    "crate::kernels::control::state_space::StateSpace::simulate_continuous",
    DimensionPolicy::StateSpaceContinuous,
    (a: dense::Mat, b: dense::Mat, c: dense::Mat, d: dense::Mat, x0: Vec<f64>, inputs: dense::Mat, times: Vec<f64>) -> Result<dense::Mat, OpError>,
    |a: dense::Mat, b: dense::Mat, c: dense::Mat, d: dense::Mat, x0: Vec<f64>, inputs: dense::Mat, times: Vec<f64>| {
        let (model, shape) = super::state_space_from_matrices(a, b, c, d, None)?;
        super::state_space_trajectory(&model, &shape, &x0, inputs, Some(&times))
    }
);
