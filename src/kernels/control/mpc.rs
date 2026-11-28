//! Model-predictive control with a condensed quadratic program.
//!
//! The state trajectory over the prediction horizon is eliminated with the
//! prediction matrices, so the only decision variables are the control moves
//! over the control horizon. The resulting dense quadratic program is solved
//! by a primal active-set method with a Bland-style anti-cycling tie-break.
//!
//! Condensation builds state and input prediction maps from the discrete
//! model and horizon. A closed-loop rollout reuses those maps, the Hessian
//! and the constraint rows, solving a shifted QP at each state.
use super::linalg;
use super::state_space::StateSpace;

/// Options of an MPC step.
#[derive(Debug, Clone, PartialEq)]
pub struct MpcOptions {
    /// Prediction horizon in samples.
    pub horizon: usize,
    /// Control horizon in samples (at most [`MpcOptions::horizon`]).
    pub control_horizon: usize,
    /// State weights, one per state.
    pub q: Vec<f64>,
    /// Input weights, one per input.
    pub r: Vec<f64>,
    /// Lower input bounds.
    pub u_min: Option<Vec<f64>>,
    /// Upper input bounds.
    pub u_max: Option<Vec<f64>>,
    /// Lower state bounds.
    pub x_min: Option<Vec<f64>>,
    /// Upper state bounds.
    pub x_max: Option<Vec<f64>>,
}

/// Result of one MPC step.
#[derive(Debug, Clone, PartialEq)]
pub struct MpcStep {
    /// Optimal first control move.
    pub u: Vec<f64>,
    /// Predicted states over the horizon.
    pub predicted_states: Vec<Vec<f64>>,
    /// Optimal value of the stage cost.
    pub objective: f64,
    /// Number of active-set iterations.
    pub iterations: usize,
}

/// Dot product of two slices.
fn dot(
    a: &[f64],
    b: &[f64],
) -> f64 {
    crate::kernels::simd::dot_f64(a, b)
}

/// Solves `min 0.5 x^T H x + f^T x` subject to `A x <= b` by a primal
/// active-set method.
///
/// A feasible point is found by successive projection onto the violated
/// half-spaces. Equality-constrained subproblems are solved through their KKT
/// systems. Ties in the blocking-constraint test and in the multiplier test
/// are broken by the smallest constraint index, which prevents cycling.
fn solve_qp(
    h: &[Vec<f64>],
    f: &[f64],
    a: &[Vec<f64>],
    b: &[f64],
) -> Result<(Vec<f64>, usize), String> {
    let nv = h.len();
    if h.iter().any(|row| row.len() != nv) {
        return Err("the qp hessian must be square".to_string());
    }
    if f.len() != nv {
        return Err("the qp linear term has the wrong length".to_string());
    }
    if a.iter().any(|row| row.len() != nv) || b.len() != a.len() {
        return Err("the qp constraints have inconsistent dimensions".to_string());
    }
    if h.iter()
        .chain(a)
        .flatten()
        .chain(f)
        .chain(b)
        .any(|value| !value.is_finite())
    {
        return Err("the condensed QP contains non-finite data".into());
    }
    let mut x = vec![0.0; nv];
    for _ in 0..500 {
        check_cancelled()?;
        let mut worst = 0.0_f64;
        let mut worst_index = None;
        for (i, row) in a.iter().enumerate() {
            let violation = dot(row, &x) - b[i];
            if !violation.is_finite() {
                return Err("QP feasibility projection is outside Float range".into());
            }
            if violation > worst {
                worst = violation;
                worst_index = Some(i);
            }
        }
        if worst <= 1e-10 {
            break;
        }
        let i = match worst_index {
            | Some(i) => i,
            | None => break,
        };
        let norm = dot(&a[i], &a[i]);
        if !norm.is_finite() || norm <= 1e-300 {
            return Err("the qp constraint cannot be resolved at this scale".into());
        }
        let step = (dot(&a[i], &x) - b[i]) / norm;
        crate::kernels::simd::axpy_f64(-step, &a[i], &mut x);
    }
    for (i, row) in a.iter().enumerate() {
        if dot(row, &x) - b[i] > 1e-7 {
            return Err("the qp constraints are infeasible".to_string());
        }
    }
    let mut active: Vec<usize> = Vec::new();
    for (i, row) in a.iter().enumerate() {
        if (dot(row, &x) - b[i]).abs() <= 1e-8 {
            active.push(i);
        }
    }
    let mut iterations = 0;
    for _ in 0..1000 {
        check_cancelled()?;
        iterations += 1;
        let mut g = vec![0.0; nv];
        for i in 0..nv {
            g[i] = f[i];
            for j in 0..nv {
                g[i] += h[i][j] * x[j];
            }
        }
        if g.iter().any(|value| !value.is_finite()) {
            return Err("the QP gradient is outside Float range".into());
        }
        let k = active.len();
        let size = nv + k;
        let mut kkt = vec![vec![0.0; size]; size];
        for i in 0..nv {
            for j in 0..nv {
                kkt[i][j] = h[i][j];
            }
        }
        for (t, &ci) in active.iter().enumerate() {
            for j in 0..nv {
                kkt[nv + t][j] = a[ci][j];
                kkt[j][nv + t] = a[ci][j];
            }
        }
        let mut rhs = vec![0.0; size];
        for i in 0..nv {
            rhs[i] = -g[i];
        }
        let solution = match linalg::solve(&kkt, &rhs) {
            | Ok(s) => s,
            | Err(_) => {
                if active.is_empty() {
                    return Err("the qp kkt system is singular".to_string());
                }
                active.pop();
                continue;
            },
        };
        if solution.iter().any(|value| !value.is_finite()) {
            return Err("the QP KKT solution is outside Float range".into());
        }
        let p: Vec<f64> = solution[..nv].to_vec();
        let lambda: Vec<f64> = solution[nv..].to_vec();
        let p_norm = p.iter().fold(0.0_f64, |m, &v| m.max(v.abs()));
        if p_norm < 1e-10 {
            let mut worst = 0.0_f64;
            let mut worst_pos = None;
            for (t, &lam) in lambda.iter().enumerate() {
                if lam < worst {
                    worst = lam;
                    worst_pos = Some(t);
                }
            }
            if worst >= -1e-9 {
                if x.iter().any(|value| !value.is_finite())
                    || a.iter().zip(b).any(|(row, &bound)| {
                        let violation = dot(row, &x) - bound;
                        !violation.is_finite() || violation > 1e-7
                    })
                {
                    return Err("the QP candidate is non-finite or infeasible".into());
                }
                return Ok((x, iterations));
            }
            if let Some(t) = worst_pos {
                active.remove(t);
            }
        } else {
            let mut alpha = 1.0_f64;
            let mut blocking: Option<usize> = None;
            for i in 0..a.len() {
                if active.contains(&i) {
                    continue;
                }
                let slope = dot(&a[i], &p);
                if !slope.is_finite() {
                    return Err("the QP constraint direction is outside Float range".into());
                }
                if slope > 1e-12 {
                    let t = (b[i] - dot(&a[i], &x)) / slope;
                    if t < alpha - 1e-12 {
                        alpha = t;
                        blocking = Some(i);
                    } else if (t - alpha).abs() <= 1e-12 && blocking.is_none_or(|j| i < j) {
                        blocking = Some(i);
                    }
                }
            }
            if alpha < 0.0 {
                alpha = 0.0;
            }
            for j in 0..nv {
                x[j] += alpha * p[j];
            }
            if alpha < 1.0 - 1e-12 {
                if let Some(i) = blocking
                    && !active.contains(&i)
                {
                    active.push(i);
                }
            } else {
                for i in 0..a.len() {
                    if !active.contains(&i) && (dot(&a[i], &x) - b[i]).abs() <= 1e-8 {
                        active.push(i);
                    }
                }
            }
        }
    }
    Err("the qp active-set iteration did not converge".to_string())
}

#[derive(Clone, Copy)]
enum Bound {
    Control(f64),
    StateUpper { index: usize, limit: f64 },
    StateLower { index: usize, limit: f64 },
}

impl Bound {
    fn shifted(
        self,
        base: &[f64],
    ) -> f64 {
        match self {
            | Self::Control(limit) => limit,
            | Self::StateUpper { index, limit } => limit - base[index],
            | Self::StateLower { index, limit } => base[index] - limit,
        }
    }
}

fn check_limits(
    low: Option<&[f64]>,
    high: Option<&[f64]>,
    size: usize,
) -> Result<(), String> {
    if [low, high]
        .into_iter()
        .flatten()
        .any(|values| values.len() != size || values.iter().any(|value| !value.is_finite()))
        || low
            .zip(high)
            .is_some_and(|(low, high)| low.iter().zip(high).any(|(lo, hi)| lo > hi))
    {
        return Err("MPC bounds must be finite, ordered vectors of the model dimensions".into());
    }
    Ok(())
}

fn check_cancelled() -> Result<(), String> {
    crate::cancellation::checkpoint().map_err(|failure| failure.to_string())
}

/// State and input prediction maps and the invariant QP matrices.
pub(crate) struct CondensedMpc {
    a: Vec<Vec<f64>>,
    b: Vec<Vec<f64>>,
    horizon: usize,
    q: Vec<f64>,
    r: Vec<f64>,
    state_map: Vec<Vec<f64>>,
    input_map: Vec<Vec<f64>>,
    input_transpose: Vec<Vec<f64>>,
    hessian: Vec<Vec<f64>>,
    constraints: Vec<Vec<f64>>,
    bounds: Vec<Bound>,
}

impl CondensedMpc {
    pub(crate) fn new(
        ss: &StateSpace,
        opts: &MpcOptions,
    ) -> Result<Self, String> {
        if !ss.dt.is_some_and(|dt| dt.is_finite() && dt > 0.0) {
            return Err("MPC requires a finite positive discrete sample time".into());
        }
        let n = ss.a.len();
        let nu = ss.b.first().map_or(0, Vec::len);
        if n == 0
            || nu == 0
            || ss.b.len() != n
            || ss.a.iter().any(|row| row.len() != n)
            || ss.b.iter().any(|row| row.len() != nu)
            || ss
                .a
                .iter()
                .chain(&ss.b)
                .flatten()
                .any(|value| !value.is_finite())
        {
            return Err("MPC state and control matrix shapes or values are invalid".into());
        }
        let horizon = opts.horizon;
        let moves = opts.control_horizon;
        if !(1..=horizon).contains(&moves) {
            return Err("the control horizon must lie in 1..=horizon".into());
        }
        if opts.q.len() != n
            || opts.r.len() != nu
            || opts
                .q
                .iter()
                .chain(&opts.r)
                .any(|weight| !(weight.is_finite() && *weight >= 0.0))
        {
            return Err("MPC needs finite, nonnegative weights for each state and input".into());
        }
        check_limits(opts.u_min.as_deref(), opts.u_max.as_deref(), nu)?;
        check_limits(opts.x_min.as_deref(), opts.x_max.as_deref(), n)?;
        let rows = n
            .checked_mul(horizon)
            .ok_or_else(|| "MPC prediction dimensions overflow".to_string())?;
        let variables = nu
            .checked_mul(moves)
            .ok_or_else(|| "MPC control dimensions overflow".to_string())?;
        let mut powers = Vec::with_capacity(horizon);
        let mut blocks = Vec::with_capacity(horizon);
        let mut power = linalg::identity(n);
        for _ in 0..horizon {
            check_cancelled()?;
            blocks.push(linalg::matmul(&power, &ss.b)?);
            power = linalg::matmul(&ss.a, &power)?;
            powers.push(power.clone());
        }
        let state_map: Vec<_> = powers.into_iter().flatten().collect();
        let mut input_map = vec![vec![0.0; variables]; rows];
        for i in 1..=horizon {
            for j in 0..moves.min(i) {
                let block = &blocks[i - 1 - j];
                for (r, entries) in block.iter().enumerate() {
                    let column = j * nu;
                    input_map[(i - 1) * n + r][column..column + nu].copy_from_slice(entries);
                }
            }
        }
        let input_transpose = linalg::transpose(&input_map);
        let weighted: Vec<Vec<_>> = input_map
            .iter()
            .enumerate()
            .map(|(row, values)| values.iter().map(|value| value * opts.q[row % n]).collect())
            .collect();
        let mut hessian = linalg::matmul(&input_transpose, &weighted)?;
        for (i, row) in hessian.iter_mut().enumerate() {
            for value in row.iter_mut() {
                *value *= 2.0;
            }
            row[i] += 2.0 * opts.r[i % nu];
        }
        let mut constraints = Vec::new();
        let mut bounds = Vec::new();
        for (side, limits) in [(1.0, &opts.u_max), (-1.0, &opts.u_min)] {
            if let Some(limits) = limits {
                for step in 0..moves {
                    for (input, &limit) in limits.iter().enumerate() {
                        let mut row = vec![0.0; variables];
                        row[step * nu + input] = side;
                        constraints.push(row);
                        bounds.push(Bound::Control(side * limit));
                    }
                }
            }
        }
        for (upper, limits) in [(true, &opts.x_max), (false, &opts.x_min)] {
            if let Some(limits) = limits {
                for (index, row) in input_map.iter().enumerate() {
                    let limit = limits[index % n];
                    if upper {
                        constraints.push(row.clone());
                        bounds.push(Bound::StateUpper { index, limit });
                    } else {
                        constraints.push(row.iter().map(|&value| -value).collect());
                        bounds.push(Bound::StateLower { index, limit });
                    }
                }
            }
        }
        if state_map
            .iter()
            .chain(&input_map)
            .chain(&hessian)
            .flatten()
            .any(|value| !value.is_finite())
        {
            return Err("MPC prediction matrices are outside Float range".into());
        }
        Ok(Self {
            a: ss.a.clone(),
            b: ss.b.clone(),
            horizon,
            q: opts.q.clone(),
            r: opts.r.clone(),
            state_map,
            input_map,
            input_transpose,
            hessian,
            constraints,
            bounds,
        })
    }

    fn check_state(
        &self,
        state: &[f64],
        reference: &[f64],
    ) -> Result<(), String> {
        let n = self.a.len();
        if state.len() != n
            || reference.len() != n
            || state
                .iter()
                .chain(reference)
                .any(|value| !value.is_finite())
        {
            return Err("state and reference must be finite vectors of the state dimension".into());
        }
        Ok(())
    }

    fn solve(
        &self,
        state: &[f64],
        reference: &[f64],
    ) -> Result<(Vec<f64>, Vec<f64>, usize), String> {
        check_cancelled()?;
        let n = self.a.len();
        let base = linalg::matvec(&self.state_map, state);
        let weighted: Vec<_> = base
            .iter()
            .enumerate()
            .map(|(i, &value)| self.q[i % n] * (value - reference[i % n]))
            .collect();
        let mut gradient = linalg::matvec(&self.input_transpose, &weighted);
        for value in &mut gradient {
            *value *= 2.0;
        }
        let shifted: Vec<_> = self
            .bounds
            .iter()
            .map(|bound| bound.shifted(&base))
            .collect();
        let (moves, iterations) = solve_qp(&self.hessian, &gradient, &self.constraints, &shifted)?;
        Ok((base, moves, iterations))
    }

    fn step(
        &self,
        state: &[f64],
        reference: &[f64],
    ) -> Result<MpcStep, String> {
        self.check_state(state, reference)?;
        let (base, moves, iterations) = self.solve(state, reference)?;
        let n = self.a.len();
        let nu = self.r.len();
        let controlled = linalg::matvec(&self.input_map, &moves);
        let mut predicted_states = Vec::with_capacity(self.horizon);
        let mut objective = 0.0;
        for (open, actuated) in base.chunks_exact(n).zip(controlled.chunks_exact(n)) {
            let predicted: Vec<_> = open
                .iter()
                .zip(actuated)
                .map(|(unforced, forced)| unforced + forced)
                .collect();
            for (i, value) in predicted.iter().enumerate() {
                let delta = value - reference[i];
                objective += self.q[i] * delta * delta;
            }
            predicted_states.push(predicted);
        }
        for (i, &value) in moves.iter().enumerate() {
            objective += self.r[i % nu] * value * value;
        }
        if !objective.is_finite()
            || predicted_states
                .iter()
                .flatten()
                .any(|value| !value.is_finite())
        {
            return Err("MPC prediction is outside Float range".into());
        }
        Ok(MpcStep {
            u: moves[..nu].to_vec(),
            predicted_states,
            objective,
            iterations,
        })
    }

    pub(crate) fn closed_loop(
        &self,
        initial: &[f64],
        reference: &[f64],
        steps: usize,
    ) -> Result<Vec<Vec<f64>>, String> {
        self.check_state(initial, reference)?;
        let capacity = steps
            .checked_add(1)
            .ok_or_else(|| "MPC trajectory length overflow".to_string())?;
        let mut trajectory = Vec::with_capacity(capacity);
        let mut state = initial.to_vec();
        trajectory.push(state.clone());
        let nu = self.r.len();
        for _ in 0..steps {
            let (_, moves, _) = self.solve(&state, reference)?;
            let ax = linalg::matvec(&self.a, &state);
            let bu = linalg::matvec(&self.b, &moves[..nu]);
            for (value, (open, actuated)) in state.iter_mut().zip(ax.into_iter().zip(bu)) {
                *value = open + actuated;
            }
            if state.iter().any(|value| !value.is_finite()) {
                return Err("MPC trajectory is outside Float range".into());
            }
            trajectory.push(state.clone());
        }
        Ok(trajectory)
    }
}

/// Computes one MPC step.
///
/// The state trajectory is eliminated before solving the control moves.
/// Moves beyond the control horizon are held at zero.
///
/// # Errors
/// Returns an error for invalid data, an unsupported model, or a failed QP.
pub fn mpc_step(
    ss: &StateSpace,
    x0: &[f64],
    reference: &[f64],
    opts: &MpcOptions,
) -> Result<MpcStep, String> {
    CondensedMpc::new(ss, opts)?.step(x0, reference)
}

/// Runs the MPC controller in closed loop for `steps` samples. Condensation is
/// shared by all steps. The result starts with the initial state.
///
/// # Errors
/// As [`mpc_step`], or when the trajectory exceeds the finite Float range.
pub fn mpc_closed_loop(
    ss: &StateSpace,
    x0: &[f64],
    reference: &[f64],
    steps: usize,
    opts: &MpcOptions,
) -> Result<Vec<Vec<f64>>, String> {
    CondensedMpc::new(ss, opts)?.closed_loop(x0, reference, steps)
}

#[cfg(test)]
mod tests {
    use super::StateSpace;
    use super::*;

    fn double_integrator() -> StateSpace {
        StateSpace::discrete(
            vec![vec![1.0, 1.0], vec![0.0, 1.0]],
            vec![vec![0.5], vec![1.0]],
            vec![vec![1.0, 0.0], vec![0.0, 1.0]],
            vec![vec![0.0], vec![0.0]],
            1.0,
        )
        .unwrap()
    }

    fn options() -> MpcOptions {
        MpcOptions {
            horizon: 10,
            control_horizon: 10,
            q: vec![1.0, 1.0],
            r: vec![0.05],
            u_min: Some(vec![-1.0]),
            u_max: Some(vec![1.0]),
            x_min: None,
            x_max: None,
        }
    }

    #[test]
    fn step_respects_input_bounds() {
        let ss = double_integrator();
        let step = mpc_step(&ss, &[0.0, 0.0], &[1.0, 0.0], &options()).unwrap();
        assert!(step.u[0] >= -1.0 - 1e-9 && step.u[0] <= 1.0 + 1e-9);
        assert_eq!(step.predicted_states.len(), 10);
        assert!(step.objective.is_finite());
    }

    #[test]
    fn closed_loop_reaches_reference() {
        let ss = double_integrator();
        let trajectory = mpc_closed_loop(&ss, &[0.0, 0.0], &[1.0, 0.0], 20, &options()).unwrap();
        assert_eq!(trajectory.len(), 21);
        let final_state = trajectory.last().unwrap();
        assert!((final_state[0] - 1.0).abs() < 0.05);
        assert!(final_state[1].abs() < 0.05);
    }

    #[test]
    fn continuous_model_is_rejected() {
        let ss = StateSpace::continuous(
            vec![vec![0.0, 1.0], vec![0.0, 0.0]],
            vec![vec![0.0], vec![1.0]],
            vec![vec![1.0, 0.0], vec![0.0, 1.0]],
            vec![vec![0.0], vec![0.0]],
        )
        .unwrap();
        assert!(mpc_step(&ss, &[0.0, 0.0], &[1.0, 0.0], &options()).is_err());
    }

    #[test]
    fn state_constraints_are_respected() {
        let ss = double_integrator();
        let mut opts = options();
        opts.u_min = None;
        opts.u_max = None;
        opts.x_max = Some(vec![0.5, 10.0]);
        opts.x_min = Some(vec![-10.0, -10.0]);
        let step = mpc_step(&ss, &[0.0, 0.0], &[1.0, 0.0], &opts).unwrap();
        for state in &step.predicted_states {
            assert!(state[0] <= 0.5 + 1e-6);
        }
    }
}
