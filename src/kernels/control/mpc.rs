//! Model-predictive control with a condensed quadratic program.
//!
//! The state trajectory over the prediction horizon is eliminated with the
//! prediction matrices, so the only decision variables are the control moves
//! over the control horizon. The resulting dense quadratic program is solved
//! by a primal active-set method with a Bland-style anti-cycling tie-break.
//!
//! # Complexity
//!
//! The cost is documented as `O(horizon * n_u^3)` per step: linear in the
//! prediction horizon and independent of the state dimension in the
//! factorisation, because the states have been eliminated before the
//! factorisation.

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
    let mut x = vec![0.0; nv];
    for _ in 0..500 {
        let mut worst = 0.0_f64;
        let mut worst_index = None;
        for (i, row) in a.iter().enumerate() {
            let violation = dot(row, &x) - b[i];
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
        if norm <= 1e-300 {
            return Err("the qp constraints are infeasible".to_string());
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
        iterations += 1;
        let mut g = vec![0.0; nv];
        for i in 0..nv {
            g[i] = f[i];
            for j in 0..nv {
                g[i] += h[i][j] * x[j];
            }
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

/// Extracts the discrete `(A, B)` pair or reports that the model is
/// continuous.
fn discrete_pair(ss: &StateSpace) -> Result<(Vec<Vec<f64>>, Vec<Vec<f64>>), String> {
    if ss.dt.is_none() {
        return Err("mpc requires a discrete-time model".to_string());
    }
    Ok((ss.a.clone(), ss.b.clone()))
}

/// Computes one MPC step.
///
/// The quadratic program is condensed: the predicted states are eliminated
/// with the prediction matrices, leaving only the control moves over the
/// control horizon. Moves beyond the control horizon are held at zero.
///
/// # Errors
/// Returns an error for a continuous model, inconsistent option lengths, or a
/// failed quadratic program.
pub fn mpc_step(
    ss: &StateSpace,
    x0: &[f64],
    reference: &[f64],
    opts: &MpcOptions,
) -> Result<MpcStep, String> {
    let (a, b) = discrete_pair(ss)?;
    let n = a.len();
    let nu = b.first().map_or(0, Vec::len);
    if nu == 0 {
        return Err("the model must have at least one input".to_string());
    }
    if x0.len() != n {
        return Err("x0 must have n entries".to_string());
    }
    if reference.len() != n {
        return Err("reference must have n entries".to_string());
    }
    let big_n = opts.horizon;
    let nc = opts.control_horizon;
    if big_n == 0 {
        return Err("the horizon must be positive".to_string());
    }
    if nc == 0 || nc > big_n {
        return Err("the control horizon must lie in 1..=horizon".to_string());
    }
    if opts.q.len() != n {
        return Err("q must have n entries".to_string());
    }
    if opts.r.len() != nu {
        return Err("r must have nu entries".to_string());
    }
    let mut powers: Vec<Vec<Vec<f64>>> = Vec::with_capacity(big_n + 1);
    powers.push(linalg::identity(n));
    for k in 1..=big_n {
        let previous = &powers[k - 1];
        powers.push(linalg::matmul(&a, previous)?);
    }
    let rows = n * big_n;
    let cols = nu * nc;
    let mut sx = vec![vec![0.0; n]; rows];
    for i in 1..=big_n {
        for r in 0..n {
            for c in 0..n {
                sx[(i - 1) * n + r][c] = powers[i][r][c];
            }
        }
    }
    let mut su = vec![vec![0.0; cols]; rows];
    for i in 1..=big_n {
        for j in 0..nc.min(i) {
            let block = linalg::matmul(&powers[i - 1 - j], &b)?;
            for r in 0..n {
                for c in 0..nu {
                    su[(i - 1) * n + r][j * nu + c] = block[r][c];
                }
            }
        }
    }
    let base = linalg::matvec(&sx, x0);
    let mut error = vec![0.0; rows];
    for i in 0..big_n {
        for s in 0..n {
            error[i * n + s] = base[i * n + s] - reference[s];
        }
    }
    let mut weighted = vec![0.0; rows];
    for i in 0..rows {
        weighted[i] = opts.q[i % n] * error[i];
    }
    let su_t = linalg::transpose(&su);
    let mut f = linalg::matvec(&su_t, &weighted);
    for value in &mut f {
        *value *= 2.0;
    }
    let mut qsu = su.clone();
    for i in 0..rows {
        let w = opts.q[i % n];
        for j in 0..cols {
            qsu[i][j] *= w;
        }
    }
    let hq = linalg::matmul(&su_t, &qsu)?;
    let mut h = vec![vec![0.0; cols]; cols];
    for i in 0..cols {
        for j in 0..cols {
            h[i][j] = 2.0 * hq[i][j];
        }
        h[i][i] += 2.0 * opts.r[i % nu];
    }
    let mut a_rows: Vec<Vec<f64>> = Vec::new();
    let mut b_rhs: Vec<f64> = Vec::new();
    if let Some(umax) = &opts.u_max {
        if umax.len() != nu {
            return Err("u_max must have nu entries".to_string());
        }
        for j in 0..nc {
            for c in 0..nu {
                let mut row = vec![0.0; cols];
                row[j * nu + c] = 1.0;
                a_rows.push(row);
                b_rhs.push(umax[c]);
            }
        }
    }
    if let Some(umin) = &opts.u_min {
        if umin.len() != nu {
            return Err("u_min must have nu entries".to_string());
        }
        for j in 0..nc {
            for c in 0..nu {
                let mut row = vec![0.0; cols];
                row[j * nu + c] = -1.0;
                a_rows.push(row);
                b_rhs.push(-umin[c]);
            }
        }
    }
    if let Some(xmax) = &opts.x_max {
        if xmax.len() != n {
            return Err("x_max must have n entries".to_string());
        }
        for i in 0..big_n {
            for s in 0..n {
                a_rows.push(su[i * n + s].clone());
                b_rhs.push(xmax[s] - base[i * n + s]);
            }
        }
    }
    if let Some(xmin) = &opts.x_min {
        if xmin.len() != n {
            return Err("x_min must have n entries".to_string());
        }
        for i in 0..big_n {
            for s in 0..n {
                let row: Vec<f64> = su[i * n + s].iter().map(|&v| -v).collect();
                a_rows.push(row);
                b_rhs.push(-(xmin[s] - base[i * n + s]));
            }
        }
    }
    let (uvec, iterations) = solve_qp(&h, &f, &a_rows, &b_rhs)?;
    let su_u = linalg::matvec(&su, &uvec);
    let mut predicted = Vec::with_capacity(big_n);
    for i in 0..big_n {
        let mut xi = vec![0.0; n];
        for s in 0..n {
            xi[s] = base[i * n + s] + su_u[i * n + s];
        }
        predicted.push(xi);
    }
    let mut objective = 0.0;
    for xi in &predicted {
        for s in 0..n {
            let d = xi[s] - reference[s];
            objective += opts.q[s] * d * d;
        }
    }
    for j in 0..nc {
        for c in 0..nu {
            let u = uvec[j * nu + c];
            objective += opts.r[c] * u * u;
        }
    }
    let u = uvec[..nu].to_vec();
    Ok(MpcStep {
        u,
        predicted_states: predicted,
        objective,
        iterations,
    })
}

/// Runs the MPC controller in closed loop for `steps` samples.
///
/// Returns the state trajectory including the initial state, so the result has
/// `steps + 1` entries.
///
/// # Errors
/// As [`mpc_step`].
pub fn mpc_closed_loop(
    ss: &StateSpace,
    x0: &[f64],
    reference: &[f64],
    steps: usize,
    opts: &MpcOptions,
) -> Result<Vec<Vec<f64>>, String> {
    let (a, b) = discrete_pair(ss)?;
    let mut x = x0.to_vec();
    let mut trajectory = Vec::with_capacity(steps + 1);
    trajectory.push(x.clone());
    for _ in 0..steps {
        let step = mpc_step(ss, &x, reference, opts)?;
        let ax = linalg::matvec(&a, &x);
        let bu = linalg::matvec(&b, &step.u);
        for (i, value) in x.iter_mut().enumerate() {
            *value = ax[i] + bu.get(i).copied().unwrap_or(0.0);
        }
        trajectory.push(x.clone());
    }
    Ok(trajectory)
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
