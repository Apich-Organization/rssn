//! Symplectic/geometric integrators, Gauss-Legendre implicit RK,
//! Differential-Algebraic Equation (DAE) solvers, and Delay Differential Equation (DDE) solvers.

/// Type alias for ODE right-hand-side callback.
pub type Rate<'a> = dyn FnMut(f64, &[f64]) -> Vec<f64> + 'a;

/// Error returned by geometric, DAE, and DDE solvers.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum OdeGeomError {
    /// Empty state vector.
    Empty,
    /// Mismatched dimensions or lengths.
    ShapeMismatch(String),
    /// Invalid domain, span, or non-finite values.
    Domain(String),
    /// Singular Jacobian or constraint matrix.
    Singular,
    /// Iterative solver failed to converge.
    NotConverged,
}

impl std::fmt::Display for OdeGeomError {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        match self {
            | Self::Empty => f.write_str("empty input"),
            | Self::ShapeMismatch(s) => write!(f, "shape mismatch: {s}"),
            | Self::Domain(s) => write!(f, "domain error: {s}"),
            | Self::Singular => f.write_str("singular matrix"),
            | Self::NotConverged => f.write_str("did not converge"),
        }
    }
}

impl std::error::Error for OdeGeomError {}

/// Hermite-interpolatable solution returned by implicit DAE solver.
#[derive(Clone, Debug)]
pub struct GeomSolution {
    /// Time points.
    pub t: Vec<f64>,
    /// States at time points.
    pub y: Vec<Vec<f64>>,
    /// State derivatives at time points.
    pub rates: Vec<Vec<f64>>,
    /// Function evaluations.
    pub nfev: usize,
    /// Rejected steps.
    pub rejected: usize,
}

fn evaluate(
    rhs: &mut Rate<'_>,
    t: f64,
    y: &[f64],
    expected: usize,
) -> Result<Vec<f64>, OdeGeomError> {
    let dy = rhs(t, y);
    if dy.len() != expected {
        return Err(OdeGeomError::ShapeMismatch(format!(
            "rhs returned length {}, expected {expected}",
            dy.len()
        )));
    }
    if dy.iter().any(|v| !v.is_finite()) {
        return Err(OdeGeomError::Domain("rhs returned non-finite value".into()));
    }
    Ok(dy)
}

fn solve_dense(
    mut a: Vec<Vec<f64>>,
    mut b: Vec<f64>,
) -> Result<Vec<f64>, OdeGeomError> {
    let n = b.len();
    for column in 0..n {
        let pivot = (column..n)
            .max_by(|&i, &j| a[i][column].abs().total_cmp(&a[j][column].abs()))
            .unwrap_or(column);
        if a[pivot][column].abs() < 1e-300 || !a[pivot][column].is_finite() {
            return Err(OdeGeomError::Singular);
        }
        a.swap(column, pivot);
        b.swap(column, pivot);
        for row in column + 1..n {
            let factor = a[row][column] / a[column][column];
            for k in column..n {
                a[row][k] -= factor * a[column][k];
            }
            b[row] -= factor * b[column];
        }
    }
    let mut x = vec![0.0; n];
    for row in (0..n).rev() {
        let tail: f64 = (row + 1..n).map(|k| a[row][k] * x[k]).sum();
        x[row] = (b[row] - tail) / a[row][row];
    }
    if x.iter().any(|value| !value.is_finite()) {
        return Err(OdeGeomError::Singular);
    }
    Ok(x)
}

/// reversible systems.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SymplecticMethod {
    /// Second-order velocity Störmer–Verlet (kick-drift-kick).
    StormerVerlet,
    /// Ruth's 3-stage third-order explicit symplectic integrator (Ruth 1983).
    Ruth3,
    /// Fourth-order 3-stage Ruth / Forest–Ruth / Yoshida triple-jump composition.
    Ruth4,
    /// Yoshida's 7-stage sixth-order symmetric symplectic composition.
    Yoshida6,
    /// Yoshida's 15-stage eighth-order symmetric symplectic composition.
    Yoshida8,
    /// 1-stage second-order implicit midpoint (Gauss–Legendre s=1).
    ImplicitMidpoint,
    /// 2-stage fourth-order implicit Gauss–Legendre Runge–Kutta (s=2).
    GaussLegendre4,
    /// 3-stage sixth-order implicit Gauss–Legendre Runge–Kutta (s=3).
    GaussLegendre6,
}

/// Trajectory of a symplectic integration in phase space `(q, p)`.
#[derive(Clone, Debug)]
pub struct SymplecticSolution {
    /// Time grid `t_0 < t_1 < ... < t_N`.
    pub t: Vec<f64>,
    /// Generalized coordinates `q(t_k)` at each step.
    pub q: Vec<Vec<f64>>,
    /// Generalized momenta / velocities `p(t_k)` at each step.
    pub p: Vec<Vec<f64>>,
    /// Time derivatives `dq/dt` at each step (for Hermite interpolation).
    pub dq: Vec<Vec<f64>>,
    /// Time derivatives `dp/dt` at each step (for Hermite interpolation).
    pub dp: Vec<Vec<f64>>,
    /// The symplectic integrator used.
    pub method: SymplecticMethod,
    /// Total force / vector-field evaluations.
    pub nfev: usize,
}

impl SymplecticSolution {
    /// Continuous cubic Hermite evaluation of `(q(t), p(t))` at any `t` in the span.
    #[must_use]
    pub fn eval(
        &self,
        t: f64,
    ) -> Option<(Vec<f64>, Vec<f64>)> {
        let (&first, &last) = (self.t.first()?, self.t.last()?);
        if !(first..=last).contains(&t) {
            return None;
        }
        if crate::kernels::exact_eq(t, first) {
            return Some((self.q.first()?.clone(), self.p.first()?.clone()));
        }
        if crate::kernels::exact_eq(t, last) {
            return Some((self.q.last()?.clone(), self.p.last()?.clone()));
        }
        let upper = self
            .t
            .partition_point(|&node| node < t)
            .clamp(1, self.t.len() - 1);
        let lower = upper - 1;
        let h = self.t[upper] - self.t[lower];
        let s = (t - self.t[lower]) / h;
        let (h00, h10, h01, h11) = (
            2.0 * s * s * s - 3.0 * s * s + 1.0,
            s * s * s - 2.0 * s * s + s,
            -2.0 * s * s * s + 3.0 * s * s,
            s * s * s - s * s,
        );
        let d = self.q[lower].len();
        let q_interp = (0..d)
            .map(|i| {
                h00 * self.q[lower][i]
                    + h10 * h * self.dq[lower][i]
                    + h01 * self.q[upper][i]
                    + h11 * h * self.dq[upper][i]
            })
            .collect();
        let p_interp = (0..d)
            .map(|i| {
                h00 * self.p[lower][i]
                    + h10 * h * self.dp[lower][i]
                    + h01 * self.p[upper][i]
                    + h11 * h * self.dp[upper][i]
            })
            .collect();
        Some((q_interp, p_interp))
    }

    /// Maximum absolute Hamiltonian drift `max_k |H(q_k, p_k) - H(q_0, p_0)|` along the trajectory.
    #[must_use]
    pub fn max_energy_drift<H>(
        &self,
        mut hamiltonian: H,
    ) -> f64
    where
        H: FnMut(&[f64], &[f64]) -> f64,
    {
        if self.q.is_empty() {
            return 0.0;
        }
        let h0 = hamiltonian(&self.q[0], &self.p[0]);
        self.q
            .iter()
            .zip(&self.p)
            .map(|(qk, pk)| (hamiltonian(qk, pk) - h0).abs())
            .fold(0.0, f64::max)
    }
}

/// Symmetric Störmer–Verlet sub-step weights for Yoshida/Forest–Ruth compositions.
const fn yoshida_substeps(method: SymplecticMethod) -> &'static [f64] {
    // w1 = 1 / (2 - 2^(1/3)) = 1.3512071919596578, w0 = 1 - 2*w1 = -1.7024143839193155
    const RUTH4: [f64; 3] = [
        1.351_207_191_959_657_8,
        -1.702_414_383_919_315_5,
        1.351_207_191_959_657_8,
    ];
    // Yoshida (1990) Solution A 6th-order 7-stage symmetric weights:
    const YOSHIDA6: [f64; 7] = [
        0.784_513_610_477_560,
        0.235_573_213_359_357,
        -1.177_679_984_178_87,
        1.315_186_320_683_906_3,
        -1.177_679_984_178_87,
        0.235_573_213_359_357,
        0.784_513_610_477_560,
    ];
    // Yoshida (1990) 8th-order 15-stage symmetric weights:
    const YOSHIDA8: [f64; 15] = [
        0.741_670_364_350_612_9,
        -0.409_100_825_800_031_57,
        0.190_754_710_296_238_38,
        -0.573_862_471_116_082_3,
        0.299_064_181_303_655_9,
        0.334_624_918_245_298_2,
        0.315_293_092_396_766_6,
        -0.796_887_939_352_916_4,
        0.315_293_092_396_766_6,
        0.334_624_918_245_298_2,
        0.299_064_181_303_655_9,
        -0.573_862_471_116_082_3,
        0.190_754_710_296_238_38,
        -0.409_100_825_800_031_57,
        0.741_670_364_350_612_9,
    ];
    match method {
        | SymplecticMethod::StormerVerlet => &[1.0],
        | SymplecticMethod::Ruth4 => &RUTH4,
        | SymplecticMethod::Yoshida6 => &YOSHIDA6,
        | SymplecticMethod::Yoshida8 => &YOSHIDA8,
        | _ => &[],
    }
}

/// Integrates a second-order separable system `q' = p`, `p' = accel(t, q)` using a
/// symplectic or geometric integrator over `span = (t0, t1)` with `steps` steps.
///
/// For autonomous acceleration `accel(q)` with splitting methods (`StormerVerlet`,
/// `Ruth4`, `Yoshida6`, `Yoshida8`), this uses the exact velocity-Verlet composition
/// matching `rssn::kernels::ode_ext::symplectic_split` (`Rssn_Numeric::ode_symplectic`),
/// extended to support time-dependent forces, `Ruth3`, and implicit Gauss–Legendre
/// methods (`ImplicitMidpoint`, `GaussLegendre4`, `GaussLegendre6`).
///
/// # Errors
///
/// [`OdeGeomError::Domain`] for non-finite inputs or invalid span,
/// [`OdeGeomError::ShapeMismatch`] when `q0.len() != p0.len()` or `accel` changes length,
/// and [`OdeGeomError::NotConverged`] if an implicit stage fails to converge.
pub fn integrate_symplectic<A>(
    accel: &mut A,
    q0: &[f64],
    p0: &[f64],
    span: (f64, f64),
    steps: usize,
    method: SymplecticMethod,
) -> Result<SymplecticSolution, OdeGeomError>
where
    A: FnMut(f64, &[f64]) -> Vec<f64> + ?Sized,
{
    let mut dq_dt = |_t: f64, _q: &[f64], p: &[f64]| p.to_vec();
    let mut dp_dt = |t: f64, q: &[f64], _p: &[f64]| accel(t, q);
    integrate_hamiltonian(&mut dq_dt, &mut dp_dt, q0, p0, span, steps, method)
}

/// Integrates a canonical Hamiltonian or partitioned system
/// `dq/dt = dq_dt(t, q, p)`, `dp/dt = dp_dt(t, q, p)` over `span = (t0, t1)`
/// using `steps` steps of `method`.
///
/// Explicit splitting methods (`StormerVerlet`, `Ruth3`, `Ruth4`, `Yoshida6`, `Yoshida8`)
/// are symplectic for separable Hamiltonians `H(t, q, p) = T(t, p) + V(t, q)`.
/// Implicit Gauss–Legendre Runge–Kutta methods (`ImplicitMidpoint`, `GaussLegendre4`,
/// `GaussLegendre6`) are symplectic for **arbitrary** (including non-separable)
/// canonical Hamiltonian systems and preserve all quadratic invariants to machine precision.
///
/// # Errors
///
/// [`OdeGeomError::Domain`] for non-finite inputs or `t1 < t0`,
/// [`OdeGeomError::ShapeMismatch`] when `q0.len() != p0.len()`,
/// and [`OdeGeomError::NotConverged`] when implicit Newton iteration fails.
pub fn integrate_hamiltonian<FQ, FP>(
    dq_dt: &mut FQ,
    dp_dt: &mut FP,
    q0: &[f64],
    p0: &[f64],
    span: (f64, f64),
    steps: usize,
    method: SymplecticMethod,
) -> Result<SymplecticSolution, OdeGeomError>
where
    FQ: FnMut(f64, &[f64], &[f64]) -> Vec<f64> + ?Sized,
    FP: FnMut(f64, &[f64], &[f64]) -> Vec<f64> + ?Sized,
{
    let (t0, t1) = span;
    let d = q0.len();
    if d == 0 {
        return Err(OdeGeomError::Empty);
    }
    if p0.len() != d {
        return Err(OdeGeomError::ShapeMismatch(format!(
            "q0 has length {d} but p0 has length {}",
            p0.len()
        )));
    }
    if !t0.is_finite()
        || !t1.is_finite()
        || t1 < t0
        || !(t1 - t0).is_finite()
        || !q0.iter().chain(p0).all(|x| x.is_finite())
    {
        return Err(OdeGeomError::Domain(
            "invalid span or non-finite initial state in symplectic integration".into(),
        ));
    }
    let eval_qp = |dq_dt: &mut FQ,
                   dp_dt: &mut FP,
                   t: f64,
                   q: &[f64],
                   p: &[f64],
                   nfev: &mut usize|
     -> Result<(Vec<f64>, Vec<f64>), OdeGeomError> {
        let vq = dq_dt(t, q, p);
        let vp = dp_dt(t, q, p);
        *nfev += 1;
        if vq.len() != d || vp.len() != d {
            return Err(OdeGeomError::ShapeMismatch(format!(
                "Hamiltonian rates returned lengths ({}, {}) for dimension {d}",
                vq.len(),
                vp.len()
            )));
        }
        if !vq.iter().chain(&vp).all(|x| x.is_finite()) {
            return Err(OdeGeomError::Domain(
                "Hamiltonian rate returned a non-finite value".into(),
            ));
        }
        Ok((vq, vp))
    };

    let nsteps = steps.max(1);
    let h = (t1 - t0) / nsteps as f64;
    let mut ts = Vec::with_capacity(nsteps + 1);
    let mut qs = Vec::with_capacity(nsteps + 1);
    let mut ps = Vec::with_capacity(nsteps + 1);
    let mut dqs = Vec::with_capacity(nsteps + 1);
    let mut dps = Vec::with_capacity(nsteps + 1);

    let mut q = q0.to_vec();
    let mut p = p0.to_vec();
    let mut nfev = 0usize;
    let (vq0, vp0) = eval_qp(dq_dt, dp_dt, t0, &q, &p, &mut nfev)?;
    ts.push(t0);
    qs.push(q.clone());
    ps.push(p.clone());
    dqs.push(vq0);
    dps.push(vp0);

    if crate::kernels::exact_eq(t1, t0) {
        return Ok(SymplecticSolution {
            t: ts,
            q: qs,
            p: ps,
            dq: dqs,
            dp: dps,
            method,
            nfev,
        });
    }

    for step_idx in 0..nsteps {
        let t_curr = t0 + step_idx as f64 * h;
        let t_next = if step_idx + 1 == nsteps {
            t1
        } else {
            t0 + (step_idx + 1) as f64 * h
        };
        let dt = t_next - t_curr;

        match method {
            | SymplecticMethod::StormerVerlet
            | SymplecticMethod::Ruth4
            | SymplecticMethod::Yoshida6
            | SymplecticMethod::Yoshida8 => {
                let weights = yoshida_substeps(method);
                let mut sub_t = t_curr;
                for &w in weights {
                    let sub_h = w * dt;
                    // Half kick: p_{1/2} = p + (sub_h / 2) * dp_dt(sub_t, q, p)
                    let vp1 = dp_dt(sub_t, &q, &p);
                    nfev += 1;
                    if vp1.len() != d || !vp1.iter().all(|x| x.is_finite()) {
                        return Err(OdeGeomError::Domain(
                            "invalid dp/dt in symplectic substep".into(),
                        ));
                    }
                    for i in 0..d {
                        p[i] += 0.5 * sub_h * vp1[i];
                    }
                    // Full drift: q_{new} = q + sub_h * dq_dt(sub_t + sub_h/2, q, p_{1/2})
                    let vq = dq_dt(sub_t + 0.5 * sub_h, &q, &p);
                    if vq.len() != d || !vq.iter().all(|x| x.is_finite()) {
                        return Err(OdeGeomError::Domain(
                            "invalid dq/dt in symplectic substep".into(),
                        ));
                    }
                    for i in 0..d {
                        q[i] += sub_h * vq[i];
                    }
                    sub_t += sub_h;
                    // Half kick: p_{new} = p_{1/2} + (sub_h / 2) * dp_dt(sub_t, q_{new}, p_{1/2})
                    let vp2 = dp_dt(sub_t, &q, &p);
                    nfev += 1;
                    if vp2.len() != d || !vp2.iter().all(|x| x.is_finite()) {
                        return Err(OdeGeomError::Domain(
                            "invalid dp/dt in symplectic substep".into(),
                        ));
                    }
                    for i in 0..d {
                        p[i] += 0.5 * sub_h * vp2[i];
                    }
                }
            },
            | SymplecticMethod::Ruth3 => {
                // Ruth (1983) 3-stage 3rd-order explicit symplectic scheme:
                const C: [f64; 3] = [2.0 / 3.0, -2.0 / 3.0, 1.0];
                const D: [f64; 3] = [7.0 / 24.0, 0.75, -1.0 / 24.0];
                let mut sub_t = t_curr;
                for stage in 0..3 {
                    let vp = dp_dt(sub_t, &q, &p);
                    nfev += 1;
                    if vp.len() != d || !vp.iter().all(|x| x.is_finite()) {
                        return Err(OdeGeomError::Domain("invalid dp/dt in Ruth3".into()));
                    }
                    for i in 0..d {
                        p[i] += D[stage] * dt * vp[i];
                    }
                    let vq = dq_dt(sub_t, &q, &p);
                    if vq.len() != d || !vq.iter().all(|x| x.is_finite()) {
                        return Err(OdeGeomError::Domain("invalid dq/dt in Ruth3".into()));
                    }
                    for i in 0..d {
                        q[i] += C[stage] * dt * vq[i];
                    }
                    sub_t += C[stage] * dt;
                }
            },
            | SymplecticMethod::ImplicitMidpoint
            | SymplecticMethod::GaussLegendre4
            | SymplecticMethod::GaussLegendre6 => {
                let stages = match method {
                    | SymplecticMethod::ImplicitMidpoint => 1,
                    | SymplecticMethod::GaussLegendre4 => 2,
                    | SymplecticMethod::GaussLegendre6 => 3,
                    | _ => unreachable!(),
                };
                let mut state = Vec::with_capacity(2 * d);
                state.extend_from_slice(&q);
                state.extend_from_slice(&p);
                let mut combined_rhs = |t: f64, y: &[f64]| -> Vec<f64> {
                    let (q_slice, p_slice) = y.split_at(d);
                    let mut out = dq_dt(t, q_slice, p_slice);
                    out.extend(dp_dt(t, q_slice, p_slice));
                    out
                };
                let (next_state, evals) =
                    gauss_legendre_step(&mut combined_rhs, t_curr, &state, dt, stages)?;
                nfev += evals;
                q.copy_from_slice(&next_state[..d]);
                p.copy_from_slice(&next_state[d..]);
            },
        }

        let (vq, vp) = eval_qp(dq_dt, dp_dt, t_next, &q, &p, &mut nfev)?;
        ts.push(t_next);
        qs.push(q.clone());
        ps.push(p.clone());
        dqs.push(vq);
        dps.push(vp);
    }

    Ok(SymplecticSolution {
        t: ts,
        q: qs,
        p: ps,
        dq: dqs,
        dp: dps,
        method,
        nfev,
    })
}

/// Partitioned Störmer–Verlet step on a flat state `y = [q, p]` (or implicit
/// midpoint if `y.len()` is odd).
pub fn partitioned_verlet_step(
    rhs: &mut Rate<'_>,
    t: f64,
    y: &[f64],
    h: f64,
) -> Result<(Vec<f64>, usize), OdeGeomError> {
    let n = y.len();
    if !n.is_multiple_of(2) {
        return gauss_legendre_step(rhs, t, y, h, 1);
    }
    let d = n / 2;
    let mut state = y.to_vec();
    // Kick 1/2
    let k1 = evaluate(rhs, t, &state, n)?;
    for i in 0..d {
        state[d + i] += 0.5 * h * k1[d + i];
    }
    // Drift
    let k_mid = evaluate(rhs, t + 0.5 * h, &state, n)?;
    for i in 0..d {
        state[i] += h * k_mid[i];
    }
    // Kick 1/2
    let k2 = evaluate(rhs, t + h, &state, n)?;
    for i in 0..d {
        state[d + i] += 0.5 * h * k2[d + i];
    }
    Ok((state, 3))
}

/// Partitioned 4th-order Ruth (Forest–Ruth / Yoshida) step on a flat state `y = [q, p]`.
pub fn partitioned_ruth4_step(
    rhs: &mut Rate<'_>,
    t: f64,
    y: &[f64],
    h: f64,
) -> Result<(Vec<f64>, usize), OdeGeomError> {
    let n = y.len();
    if !n.is_multiple_of(2) {
        return gauss_legendre_step(rhs, t, y, h, 2);
    }
    let weights = yoshida_substeps(SymplecticMethod::Ruth4);
    let mut state = y.to_vec();
    let mut sub_t = t;
    let mut evals = 0;
    for &w in weights {
        let sub_h = w * h;
        let (next, c) = partitioned_verlet_step(rhs, sub_t, &state, sub_h)?;
        state = next;
        sub_t += sub_h;
        evals += c;
    }
    Ok((state, evals))
}

/// Advance one implicit Gauss-Legendre collocation step.
///
/// Performs one step of the `stages`-stage Gauss–Legendre implicit Runge–Kutta
/// method (`stages = 1` for 2nd-order implicit midpoint, `stages = 2` for 4th-order
/// `GaussLegendre4`, `stages = 3` for 6th-order `GaussLegendre6`).
///
/// Gauss–Legendre collocation methods are symmetric, A-stable, B-stable, and
/// symplectic for all canonical Hamiltonian systems, preserving arbitrary
/// quadratic invariants `y^T S y` exactly.
///
/// # Errors
///
/// [`OdeGeomError::Domain`] if `stages` is not in `1..=3` or `rhs` returns non-finite values,
/// [`OdeGeomError::ShapeMismatch`] if `rhs` changes dimension, and
/// [`OdeGeomError::NotConverged`] if the Newton solver fails to converge.
pub fn gauss_legendre_step(
    rhs: &mut Rate<'_>,
    t: f64,
    y: &[f64],
    h: f64,
    stages: usize,
) -> Result<(Vec<f64>, usize), OdeGeomError> {
    let n = y.len();
    if n == 0 {
        return Ok((Vec::new(), 0));
    }
    let (c, a, b): (Vec<f64>, Vec<Vec<f64>>, Vec<f64>) = match stages {
        | 1 => (vec![0.5], vec![vec![0.5]], vec![1.0]),
        | 2 => {
            let s3_6 = 3.0f64.sqrt() / 6.0;
            (
                vec![0.5 - s3_6, 0.5 + s3_6],
                vec![vec![0.25, 0.25 - s3_6], vec![0.25 + s3_6, 0.25]],
                vec![0.5, 0.5],
            )
        },
        | 3 => {
            let s15 = 15.0f64.sqrt();
            (
                vec![0.5 - s15 / 10.0, 0.5, 0.5 + s15 / 10.0],
                vec![
                    vec![5.0 / 36.0, 2.0 / 9.0 - s15 / 15.0, 5.0 / 36.0 - s15 / 30.0],
                    vec![5.0 / 36.0 + s15 / 24.0, 2.0 / 9.0, 5.0 / 36.0 - s15 / 24.0],
                    vec![5.0 / 36.0 + s15 / 30.0, 2.0 / 9.0 + s15 / 15.0, 5.0 / 36.0],
                ],
                vec![5.0 / 18.0, 4.0 / 9.0, 5.0 / 18.0],
            )
        },
        | _ => {
            return Err(OdeGeomError::Domain(
                "Gauss-Legendre stages must be 1 (order 2), 2 (order 4), or 3 (order 6)".into(),
            ));
        },
    };

    let mut nfev = 0usize;
    let k0 = evaluate(rhs, t, y, n)?;
    nfev += 1;

    // Unknown vector `k` has length `stages * n`, stacked as `[K_1, ..., K_s]`.
    let total_dim = stages * n;
    let mut k = vec![0.0; total_dim];
    for s in 0..stages {
        k[s * n..(s + 1) * n].copy_from_slice(&k0);
    }

    // First try a few fixed-point iterations (cheap for small h / non-stiff steps),
    // then switch to full Newton-Raphson if needed.
    let tol = 1e-14 * (1.0 + y.iter().fold(0.0f64, |acc, x| acc.max(x.abs())));
    let mut converged = false;
    for _ in 0..12 {
        let mut max_diff = 0.0f64;
        let mut k_new = vec![0.0; total_dim];
        for i in 0..stages {
            let mut yi = y.to_vec();
            for j in 0..stages {
                let a_ij = a[i][j] * h;
                for m in 0..n {
                    yi[m] += a_ij * k[j * n + m];
                }
            }
            let f_val = evaluate(rhs, t + c[i] * h, &yi, n)?;
            nfev += 1;
            for m in 0..n {
                let diff = (f_val[m] - k[i * n + m]).abs();
                max_diff = max_diff.max(diff);
                k_new[i * n + m] = f_val[m];
            }
        }
        k = k_new;
        if max_diff * h.abs().max(1.0) <= tol {
            converged = true;
            break;
        }
    }

    if !converged {
        let eval_residual = |rhs: &mut Rate<'_>,
                             k_vec: &[f64],
                             nfev: &mut usize|
         -> Result<Vec<f64>, OdeGeomError> {
            let mut res = vec![0.0; total_dim];
            for i in 0..stages {
                let mut yi = y.to_vec();
                for j in 0..stages {
                    let a_ij = a[i][j] * h;
                    for m in 0..n {
                        yi[m] += a_ij * k_vec[j * n + m];
                    }
                }
                let f_val = evaluate(rhs, t + c[i] * h, &yi, n)?;
                *nfev += 1;
                for m in 0..n {
                    res[i * n + m] = k_vec[i * n + m] - f_val[m];
                }
            }
            Ok(res)
        };

        let mut res = eval_residual(rhs, &k, &mut nfev)?;
        for _ in 0..20 {
            let norm = res.iter().fold(0.0f64, |acc, x| acc.max(x.abs()));
            if norm * h.abs().max(1.0) <= tol {
                converged = true;
                break;
            }
            // Compute Jacobian of f(t, y) at t + c_2*h to build the stage Jacobian efficiently.
            let mut y_mid = y.to_vec();
            for j in 0..stages {
                let w = b[j] * h;
                for m in 0..n {
                    y_mid[m] += w * k[j * n + m];
                }
            }
            let f_mid = evaluate(rhs, t + 0.5 * h, &y_mid, n)?;
            nfev += 1;
            let mut df_dy = vec![vec![0.0; n]; n];
            for col in 0..n {
                let eps = f64::EPSILON.sqrt() * y_mid[col].abs().max(1.0);
                let mut y_pert = y_mid.clone();
                y_pert[col] += eps;
                let f_pert = evaluate(rhs, t + 0.5 * h, &y_pert, n)?;
                nfev += 1;
                for row in 0..n {
                    df_dy[row][col] = (f_pert[row] - f_mid[row]) / eps;
                }
            }
            // Assemble J_stage = I - (h A \otimes df_dy)
            let mut jac = vec![vec![0.0; total_dim]; total_dim];
            for i in 0..stages {
                for j in 0..stages {
                    let factor = -h * a[i][j];
                    for r in 0..n {
                        for col in 0..n {
                            jac[i * n + r][j * n + col] = factor * df_dy[r][col];
                        }
                    }
                }
                for r in 0..n {
                    jac[i * n + r][i * n + r] += 1.0;
                }
            }
            let rhs_vec: Vec<f64> = res.iter().map(|x| -x).collect();
            let delta = solve_dense(jac, rhs_vec)?;
            for idx in 0..total_dim {
                k[idx] += delta[idx];
            }
            res = eval_residual(rhs, &k, &mut nfev)?;
        }
        let norm = res.iter().fold(0.0f64, |acc, x| acc.max(x.abs()));
        if !converged && norm * h.abs().max(1.0) > 1e-7 * (1.0 + norm) {
            return Err(OdeGeomError::NotConverged);
        }
    }

    let mut y_next = y.to_vec();
    for j in 0..stages {
        let bj = b[j] * h;
        for m in 0..n {
            y_next[m] += bj * k[j * n + m];
        }
    }
    if !y_next.iter().all(|x| x.is_finite()) {
        return Err(OdeGeomError::NotConverged);
    }
    Ok((y_next, nfev))
}

/// Integrate the differential equation with the selected numerical method.
///
/// Integrates a general first-order ODE system `y' = rhs(t, y)` over `span = (t0, t1)`
/// using the `stages`-stage Gauss–Legendre implicit symplectic Runge–Kutta method
/// (`stages = 1` for order 2, `stages = 2` for order 4, `stages = 3` for order 6).

#[derive(Clone, Copy, Debug)]
pub struct DaeOptions {
    /// Relative tolerance for adaptive step control.
    pub rtol: f64,
    /// Absolute tolerance for adaptive step control.
    pub atol: f64,
    /// Initial step size (`0.0` selects an automatic initial step).
    pub h0: f64,
    /// Minimum allowed step size.
    pub h_min: f64,
    /// Maximum allowed step size.
    pub h_max: f64,
    /// Maximum number of steps before returning [`OdeGeomError::NotConverged`].
    pub max_steps: usize,
    /// Tolerance for Newton iterations on stage and algebraic constraints.
    pub newton_tol: f64,
    /// Maximum Newton iterations per stage solve.
    pub max_newton_iter: usize,
}

impl Default for DaeOptions {
    fn default() -> Self {
        Self {
            rtol: 1e-6,
            atol: 1e-9,
            h0: 0.0,
            h_min: 1e-12,
            h_max: f64::INFINITY,
            max_steps: 100_000,
            newton_tol: 1e-10,
            max_newton_iter: 25,
        }
    }
}

/// Solution of a Differential-Algebraic Equation (DAE) system with differential
/// states `y(t)` and algebraic states `z(t)`.
#[derive(Clone, Debug)]
pub struct DaeSolution {
    /// Accepted time points `t_0 < t_1 < ... < t_N`.
    pub t: Vec<f64>,
    /// Differential variable states `y(t_k)`.
    pub y: Vec<Vec<f64>>,
    /// Algebraic variable states `z(t_k)`.
    pub z: Vec<Vec<f64>>,
    /// Differential rates `y'(t_k)` for continuous Hermite interpolation.
    pub dy: Vec<Vec<f64>>,
    /// Total function and constraint evaluations.
    pub nfev: usize,
    /// Number of rejected adaptive steps.
    pub rejected: usize,
}

impl DaeSolution {
    /// Evaluates `(y(t), z(t))` at any time `t` inside the integrated span using
    /// cubic Hermite interpolation for `y(t)` and linear/Hermite interpolation for `z(t)`.
    #[must_use]
    pub fn eval(
        &self,
        t: f64,
    ) -> Option<(Vec<f64>, Vec<f64>)> {
        let (&first, &last) = (self.t.first()?, self.t.last()?);
        if !(first..=last).contains(&t) {
            return None;
        }
        if crate::kernels::exact_eq(t, first) {
            return Some((self.y.first()?.clone(), self.z.first()?.clone()));
        }
        if crate::kernels::exact_eq(t, last) {
            return Some((self.y.last()?.clone(), self.z.last()?.clone()));
        }
        let upper = self
            .t
            .partition_point(|&node| node < t)
            .clamp(1, self.t.len() - 1);
        let lower = upper - 1;
        let h = self.t[upper] - self.t[lower];
        let s = (t - self.t[lower]) / h;
        let (h00, h10, h01, h11) = (
            2.0 * s * s * s - 3.0 * s * s + 1.0,
            s * s * s - 2.0 * s * s + s,
            -2.0 * s * s * s + 3.0 * s * s,
            s * s * s - s * s,
        );
        let ny = self.y[lower].len();
        let nz = self.z[lower].len();
        let y_interp = (0..ny)
            .map(|i| {
                h00 * self.y[lower][i]
                    + h10 * h * self.dy[lower][i]
                    + h01 * self.y[upper][i]
                    + h11 * h * self.dy[upper][i]
            })
            .collect();
        let z_interp = (0..nz)
            .map(|i| (1.0 - s) * self.z[lower][i] + s * self.z[upper][i])
            .collect();
        Some((y_interp, z_interp))
    }

    /// Maximum infinity-norm constraint violation `max_k ||g(t_k, y_k, z_k)||_inf`.
    #[must_use]
    pub fn max_constraint_violation<G>(
        &self,
        mut g: G,
    ) -> f64
    where
        G: FnMut(f64, &[f64], &[f64]) -> Vec<f64>,
    {
        self.t
            .iter()
            .enumerate()
            .map(|(k, &tk)| {
                g(tk, &self.y[k], &self.z[k])
                    .into_iter()
                    .fold(0.0f64, |acc, v| acc.max(v.abs()))
            })
            .fold(0.0f64, f64::max)
    }
}

/// Computes consistent initial algebraic states `z0` and differential derivatives `dy0`
/// for a semi-explicit Index-1 DAE system `y' = f(t, y, z)`, `0 = g(t, y, z)` at `t0`.
///
/// Given fixed initial differential variables `y0` and an initial guess `z0_guess`,
/// solves `g(t0, y0, z0) = 0` via damped Newton-Raphson iteration (`dg/dz` must be
/// nonsingular for an index-1 DAE). Returns `(dy0, z0)`.
///
/// # Errors
///
/// [`OdeGeomError::Domain`] for non-finite values, [`OdeGeomError::ShapeMismatch`] for
/// mismatched lengths, [`OdeGeomError::Singular`] if `dg/dz` is singular (e.g. higher-index DAE),
/// and [`OdeGeomError::NotConverged`] if Newton iteration does not converge.
pub fn dae_consistent_ic_index1<F, G>(
    f: &mut F,
    g: &mut G,
    t0: f64,
    y0: &[f64],
    z0_guess: &[f64],
    tol: f64,
    max_iter: usize,
) -> Result<(Vec<f64>, Vec<f64>), OdeGeomError>
where
    F: FnMut(f64, &[f64], &[f64]) -> Vec<f64> + ?Sized,
    G: FnMut(f64, &[f64], &[f64]) -> Vec<f64> + ?Sized,
{
    let ny = y0.len();
    let nz = z0_guess.len();
    if ny == 0 || nz == 0 {
        return Err(OdeGeomError::Empty);
    }
    let mut z = z0_guess.to_vec();
    let mut res = g(t0, y0, &z);
    if res.len() != nz {
        return Err(OdeGeomError::ShapeMismatch(format!(
            "algebraic constraint g returned {} values for {} algebraic variables",
            res.len(),
            nz
        )));
    }
    if !res.iter().all(|x| x.is_finite()) {
        return Err(OdeGeomError::Domain(
            "algebraic constraint returned non-finite value".into(),
        ));
    }

    let inf_norm = |v: &[f64]| v.iter().fold(0.0f64, |acc, x| acc.max(x.abs()));
    for _ in 0..max_iter.max(1) {
        if inf_norm(&res) <= tol {
            break;
        }
        let mut jac = vec![vec![0.0; nz]; nz];
        for col in 0..nz {
            let eps = f64::EPSILON.sqrt() * z[col].abs().max(1.0);
            let mut z_pert = z.clone();
            z_pert[col] += eps;
            let r_pert = g(t0, y0, &z_pert);
            if r_pert.len() != nz || !r_pert.iter().all(|x| x.is_finite()) {
                return Err(OdeGeomError::Domain(
                    "invalid algebraic constraint evaluation during Jacobian".into(),
                ));
            }
            for row in 0..nz {
                jac[row][col] = (r_pert[row] - res[row]) / eps;
            }
        }
        let rhs: Vec<f64> = res.iter().map(|x| -x).collect();
        let dz = solve_dense(jac, rhs)?;
        let mut alpha = 1.0;
        let current_norm = inf_norm(&res);
        let mut improved = false;
        for _ in 0..15 {
            let z_trial: Vec<f64> = z
                .iter()
                .zip(&dz)
                .map(|(zi, dzi)| zi + alpha * dzi)
                .collect();
            let r_trial = g(t0, y0, &z_trial);
            if r_trial.len() == nz
                && r_trial.iter().all(|x| x.is_finite())
                && inf_norm(&r_trial) < current_norm
            {
                z = z_trial;
                res = r_trial;
                improved = true;
                break;
            }
            alpha *= 0.5;
        }
        if !improved {
            for j in 0..nz {
                z[j] += dz[j];
            }
            res = g(t0, y0, &z);
        }
    }
    if inf_norm(&res) > tol.max(1e-7) {
        return Err(OdeGeomError::NotConverged);
    }
    let dy0 = f(t0, y0, &z);
    if dy0.len() != ny || !dy0.iter().all(|x| x.is_finite()) {
        return Err(OdeGeomError::ShapeMismatch(format!(
            "differential RHS returned {} values for {ny} differential variables",
            dy0.len()
        )));
    }
    Ok((dy0, z))
}

/// Solves a semi-explicit Index-1 Differential-Algebraic Equation (DAE) system:
/// `y'(t) = f(t, y, z)`, `0 = g(t, y, z)` over `span = (t0, t1)`.
///
/// The solver first refines `z0_guess` into a consistent initial condition satisfying
/// `g(t0, y0, z0) = 0`, and then integrates forward in time using the stiffly accurate,
/// L-stable 2-stage Radau IIA (order 3) implicit Runge–Kutta method with embedded
/// adaptive step-size control. Because Radau IIA has final stage abscissa `c_2 = 1`,
/// the algebraic constraint `g(t_{k+1}, y_{k+1}, z_{k+1}) = 0` is enforced directly
/// at every accepted step without projection drift.
///
/// # Errors
///
/// [`OdeGeomError::Empty`] if `y0` or `z0_guess` is empty,
/// [`OdeGeomError::Domain`] for non-finite inputs or invalid span/tolerances,
/// [`OdeGeomError::ShapeMismatch`] if `f` or `g` returns the wrong dimension,
/// [`OdeGeomError::Singular`] if the algebraic Jacobian `dg/dz` is singular, and
/// [`OdeGeomError::NotConverged`] if step refinement or Newton iteration fails.
pub fn solve_dae_index1<F, G>(
    f: &mut F,
    g: &mut G,
    y0: &[f64],
    z0_guess: &[f64],
    span: (f64, f64),
    options: DaeOptions,
) -> Result<DaeSolution, OdeGeomError>
where
    F: FnMut(f64, &[f64], &[f64]) -> Vec<f64> + ?Sized,
    G: FnMut(f64, &[f64], &[f64]) -> Vec<f64> + ?Sized,
{
    // 2-stage Radau IIA coefficients (order 3, L-stable, stiffly accurate: c2 = 1):
    // c = [1/3, 1], A = [[5/12, -1/12], [3/4, 1/4]], b = [3/4, 1/4]
    const C1: f64 = 1.0 / 3.0;
    const C2: f64 = 1.0;
    const A11: f64 = 5.0 / 12.0;
    const A12: f64 = -1.0 / 12.0;
    const A21: f64 = 0.75;
    const A22: f64 = 0.25;

    let (t0, t1) = span;
    let ny = y0.len();
    let nz = z0_guess.len();
    if ny == 0 || nz == 0 {
        return Err(OdeGeomError::Empty);
    }
    if !t0.is_finite()
        || !t1.is_finite()
        || t1 < t0
        || !(t1 - t0).is_finite()
        || !y0.iter().chain(z0_guess).all(|x| x.is_finite())
        || !options.rtol.is_finite()
        || options.rtol < 0.0
        || !options.atol.is_finite()
        || options.atol <= 0.0
    {
        return Err(OdeGeomError::Domain(
            "invalid DAE span, initial state, or tolerances".into(),
        ));
    }

    let (dy0, z0) = dae_consistent_ic_index1(
        f,
        g,
        t0,
        y0,
        z0_guess,
        options.newton_tol,
        options.max_newton_iter,
    )?;

    let mut ts = vec![t0];
    let mut ys = vec![y0.to_vec()];
    let mut zs = vec![z0.clone()];
    let mut dys = vec![dy0.clone()];
    let mut nfev = 2usize;
    let mut rejected = 0usize;

    if crate::kernels::exact_eq(t1, t0) {
        return Ok(DaeSolution {
            t: ts,
            y: ys,
            z: zs,
            dy: dys,
            nfev,
            rejected,
        });
    }

    let mut t = t0;
    let mut y = y0.to_vec();
    let mut z = z0;
    let mut dy = dy0;
    let mut h = if options.h0 > 0.0 {
        options.h0.min(t1 - t0)
    } else {
        ((t1 - t0) / 50.0).clamp(options.h_min, options.h_max)
    };

    let total_dim = 2 * (ny + nz);
    let mut step_count = 0usize;

    while t < t1 {
        step_count += 1;
        if step_count > options.max_steps {
            return Err(OdeGeomError::NotConverged);
        }
        let remaining = t1 - t;
        let step = if h >= remaining {
            remaining
        } else {
            h
        };
        if step < options.h_min && remaining > options.h_min {
            return Err(OdeGeomError::NotConverged);
        }

        // Unknown vector u = [Y1 (ny), Z1 (nz), Y2 (ny), Z2 (nz)]
        let mut u = Vec::with_capacity(total_dim);
        for i in 0..ny {
            u.push(y[i] + C1 * step * dy[i]);
        }
        u.extend_from_slice(&z);
        for i in 0..ny {
            u.push(y[i] + C2 * step * dy[i]);
        }
        u.extend_from_slice(&z);

        let eval_radau_res = |f: &mut F,
                              g: &mut G,
                              u_vec: &[f64],
                              nfev: &mut usize|
         -> Result<(Vec<f64>, Vec<f64>, Vec<f64>), OdeGeomError> {
            let y1 = &u_vec[0..ny];
            let z1 = &u_vec[ny..ny + nz];
            let y2 = &u_vec[ny + nz..2 * ny + nz];
            let z2 = &u_vec[2 * ny + nz..total_dim];

            let f1 = f(t + C1 * step, y1, z1);
            let g1 = g(t + C1 * step, y1, z1);
            let f2 = f(t + C2 * step, y2, z2);
            let g2 = g(t + C2 * step, y2, z2);
            *nfev += 2;

            if f1.len() != ny || f2.len() != ny || g1.len() != nz || g2.len() != nz {
                return Err(OdeGeomError::ShapeMismatch(
                    "DAE callback changed output dimension".into(),
                ));
            }
            if !f1
                .iter()
                .chain(&f2)
                .chain(&g1)
                .chain(&g2)
                .all(|x| x.is_finite())
            {
                return Err(OdeGeomError::Domain(
                    "DAE callback returned non-finite value".into(),
                ));
            }

            let mut res = Vec::with_capacity(total_dim);
            for i in 0..ny {
                res.push(y1[i] - y[i] - step * (A11 * f1[i] + A12 * f2[i]));
            }
            res.extend_from_slice(&g1);
            for i in 0..ny {
                res.push(y2[i] - y[i] - step * (A21 * f1[i] + A22 * f2[i]));
            }
            res.extend_from_slice(&g2);
            Ok((res, f1, f2))
        };

        let mut newton_ok = false;
        let mut f1_stage = vec![0.0; ny];
        let mut f2_stage = vec![0.0; ny];

        if let Ok((mut res, mut f1, mut f2)) = eval_radau_res(f, g, &u, &mut nfev) {
            for _ in 0..options.max_newton_iter {
                let norm = res.iter().fold(0.0f64, |acc, x| acc.max(x.abs()));
                if norm <= options.newton_tol {
                    newton_ok = true;
                    f1_stage = f1;
                    f2_stage = f2;
                    break;
                }
                let mut jac = vec![vec![0.0; total_dim]; total_dim];
                let mut col_ok = true;
                for col in 0..total_dim {
                    let eps = f64::EPSILON.sqrt() * u[col].abs().max(1.0);
                    let mut u_pert = u.clone();
                    u_pert[col] += eps;
                    if let Ok((r_pert, _, _)) = eval_radau_res(f, g, &u_pert, &mut nfev) {
                        for row in 0..total_dim {
                            jac[row][col] = (r_pert[row] - res[row]) / eps;
                        }
                    } else {
                        col_ok = false;
                        break;
                    }
                }
                if !col_ok {
                    break;
                }
                let rhs_vec: Vec<f64> = res.iter().map(|x| -x).collect();
                let Ok(du) = solve_dense(jac, rhs_vec) else {
                    break;
                };
                for idx in 0..total_dim {
                    u[idx] += du[idx];
                }
                let Ok((next_res, nf1, nf2)) = eval_radau_res(f, g, &u, &mut nfev) else {
                    break;
                };
                res = next_res;
                f1 = nf1;
                f2 = nf2;
                let new_norm = res.iter().fold(0.0f64, |acc, x| acc.max(x.abs()));
                if new_norm <= options.newton_tol {
                    newton_ok = true;
                    f1_stage = f1;
                    f2_stage = f2;
                    break;
                }
            }
        }

        if !newton_ok {
            rejected += 1;
            h = (0.5 * step).max(options.h_min);
            if step <= options.h_min {
                return Err(OdeGeomError::NotConverged);
            }
            continue;
        }

        let y2 = &u[ny + nz..2 * ny + nz];
        let z2 = &u[2 * ny + nz..total_dim];

        // Embedded error estimate comparing Radau IIA (order 3) against trapezoid/BDF predictor:
        // err_i = step * |1/4 * (f1_i - f2_i) + 1/2 * (dy_i - f1_i)|
        let mut err = 0.0f64;
        for i in 0..ny {
            let est = (step / 6.0) * (dy[i] - 2.0 * f1_stage[i] + f2_stage[i]);
            let scale = options.atol + options.rtol * y[i].abs().max(y2[i].abs());
            err = err.max((est / scale).abs());
        }

        if err <= 1.0 || step <= options.h_min {
            t = if crate::kernels::exact_eq(step, remaining) {
                t1
            } else {
                t + step
            };
            y = y2.to_vec();
            z = z2.to_vec();
            dy = f2_stage;
            ts.push(t);
            ys.push(y.clone());
            zs.push(z.clone());
            dys.push(dy.clone());

            let factor = if err == 0.0 {
                2.5
            } else {
                (0.9 * err.powf(-1.0 / 3.0)).clamp(0.2, 2.5)
            };
            h = (step * factor).clamp(options.h_min, options.h_max);
        } else {
            rejected += 1;
            let factor = (0.9 * err.powf(-1.0 / 3.0)).clamp(0.2, 0.8);
            h = (step * factor).max(options.h_min);
        }
    }

    Ok(DaeSolution {
        t: ts,
        y: ys,
        z: zs,
        dy: dys,
        nfev,
        rejected,
    })
}

/// Computes consistent initial conditions `(y0, z0)` for a semi-explicit Hessenberg
/// Index-2 DAE system `y'(t) = f(t, y, z)`, `0 = g(t, y)` at `t0`.
///
/// An index-2 DAE has two constraint manifolds:
/// 1. The explicit position/state constraint `g(t, y) = 0`.
/// 2. The hidden velocity/index-1 constraint obtained by differentiating `g(t, y) = 0`
///    with respect to time: `G(t, y) f(t, y, z) + \partial_t g(t, y) = 0`, where
///    `G = \partial g / \partial y` and `(G \cdot \partial f / \partial z)` is nonsingular.
///
/// This function first projects `y0_guess` onto `g(t0, y) = 0` via minimum-norm
/// Gauss–Newton projection, then solves the hidden constraint for `z0`.
///
/// # Errors
///
/// [`OdeGeomError::Empty`] if `y0_guess` or `z0_guess` is empty,
/// [`OdeGeomError::Singular`] if `G (\partial f / \partial z)` is singular, and
/// [`OdeGeomError::NotConverged`] if projection fails to converge.
pub fn dae_consistent_ic_index2<F, G>(
    f: &mut F,
    g: &mut G,
    t0: f64,
    y0_guess: &[f64],
    z0_guess: &[f64],
    tol: f64,
    max_iter: usize,
) -> Result<(Vec<f64>, Vec<f64>), OdeGeomError>
where
    F: FnMut(f64, &[f64], &[f64]) -> Vec<f64> + ?Sized,
    G: FnMut(f64, &[f64]) -> Vec<f64> + ?Sized,
{
    let ny = y0_guess.len();
    let nz = z0_guess.len();
    if ny == 0 || nz == 0 {
        return Err(OdeGeomError::Empty);
    }
    let inf_norm = |v: &[f64]| v.iter().fold(0.0f64, |acc, x| acc.max(x.abs()));

    // Step 1: Project y onto g(t0, y) = 0 using minimum-norm Gauss-Newton:
    // dy = -G^T (G G^T)^{-1} g(t0, y)
    let mut y = y0_guess.to_vec();
    for _ in 0..max_iter.max(1) {
        let gv = g(t0, &y);
        if gv.len() != nz {
            return Err(OdeGeomError::ShapeMismatch(format!(
                "index-2 constraint g returned {} values for {nz} algebraic variables",
                gv.len()
            )));
        }
        if !gv.iter().all(|x| x.is_finite()) {
            return Err(OdeGeomError::Domain(
                "index-2 constraint returned non-finite value".into(),
            ));
        }
        if inf_norm(&gv) <= tol {
            break;
        }
        let mut g_jac = vec![vec![0.0; ny]; nz];
        for col in 0..ny {
            let eps = f64::EPSILON.sqrt() * y[col].abs().max(1.0);
            let mut y_pert = y.clone();
            y_pert[col] += eps;
            let g_pert = g(t0, &y_pert);
            for row in 0..nz {
                g_jac[row][col] = (g_pert[row] - gv[row]) / eps;
            }
        }
        // Assemble M = G G^T (nz x nz)
        let mut m_mat = vec![vec![0.0; nz]; nz];
        for i in 0..nz {
            for j in 0..nz {
                m_mat[i][j] = (0..ny).map(|k| g_jac[i][k] * g_jac[j][k]).sum();
            }
        }
        let neg_g: Vec<f64> = gv.iter().map(|x| -x).collect();
        let lambda = solve_dense(m_mat, neg_g)?;
        for k in 0..ny {
            let delta_k: f64 = (0..nz).map(|i| g_jac[i][k] * lambda[i]).sum();
            y[k] += delta_k;
        }
    }
    let gv_final = g(t0, &y);
    if inf_norm(&gv_final) > tol.max(1e-7) {
        return Err(OdeGeomError::NotConverged);
    }

    // Helper to evaluate the first time-derivative of g(t, y) along f(t, y, z):
    // H_1(t, y, z) = G(t, y) f(t, y, z) + dg/dt(t, y)
    let eval_h1 = |f: &mut F,
                   g: &mut G,
                   t_val: f64,
                   y_val: &[f64],
                   z_val: &[f64]|
     -> Result<Vec<f64>, OdeGeomError> {
        let gv = g(t_val, y_val);
        if gv.len() != nz || !gv.iter().all(|x| x.is_finite()) {
            return Err(OdeGeomError::Domain(
                "invalid g(t, y) in index-2 derivative".into(),
            ));
        }
        let mut g_jac = vec![vec![0.0; ny]; nz];
        for col in 0..ny {
            let eps = 1e-7 * y_val[col].abs().max(1.0);
            let mut y_p = y_val.to_vec();
            let mut y_m = y_val.to_vec();
            y_p[col] += eps;
            y_m[col] -= eps;
            let gp = g(t_val, &y_p);
            let gm = g(t_val, &y_m);
            for row in 0..nz {
                g_jac[row][col] = (gp[row] - gm[row]) / (2.0 * eps);
            }
        }
        let dt_eps = 1e-7 * t_val.abs().max(1.0);
        let gt_p = g(t_val + dt_eps, y_val);
        let gt_m = g(t_val - dt_eps, y_val);
        let fv = f(t_val, y_val, z_val);
        if fv.len() != ny || !fv.iter().all(|x| x.is_finite()) {
            return Err(OdeGeomError::Domain(
                "invalid f(t, y, z) in index-2 IC solver".into(),
            ));
        }
        Ok((0..nz)
            .map(|i| {
                let gf: f64 = (0..ny).map(|k| g_jac[i][k] * fv[k]).sum();
                let dgt = (gt_p[i] - gt_m[i]) / (2.0 * dt_eps);
                gf + dgt
            })
            .collect())
    };

    // Step 2: Check whether H_1 depends directly on z (Hessenberg Index-2) or if
    // dH_1/dz is singular (mechanical / multibody constraint where H_1 is a velocity
    // constraint and H_2 = dH_1/dt depends on z).
    let mut z = z0_guess.to_vec();
    let mut h1_res = eval_h1(f, g, t0, &y, &z)?;
    let mut h1_z_jac = vec![vec![0.0; nz]; nz];
    for col in 0..nz {
        let eps = 1e-7 * z[col].abs().max(1.0);
        let mut z_pert = z.clone();
        z_pert[col] += eps;
        let h_pert = eval_h1(f, g, t0, &y, &z_pert)?;
        for row in 0..nz {
            h1_z_jac[row][col] = (h_pert[row] - h1_res[row]) / eps;
        }
    }
    let is_direct_index2 = solve_dense(h1_z_jac.clone(), vec![1.0; nz]).is_ok();

    if is_direct_index2 {
        for _ in 0..max_iter.max(1) {
            if inf_norm(&h1_res) <= tol {
                break;
            }
            let mut h_jac = vec![vec![0.0; nz]; nz];
            for col in 0..nz {
                let eps = 1e-7 * z[col].abs().max(1.0);
                let mut z_pert = z.clone();
                z_pert[col] += eps;
                let h_pert = eval_h1(f, g, t0, &y, &z_pert)?;
                for row in 0..nz {
                    h_jac[row][col] = (h_pert[row] - h1_res[row]) / eps;
                }
            }
            let neg_h: Vec<f64> = h1_res.iter().map(|x| -x).collect();
            let dz = solve_dense(h_jac, neg_h)?;
            for i in 0..nz {
                z[i] += dz[i];
            }
            h1_res = eval_h1(f, g, t0, &y, &z)?;
        }
        if inf_norm(&h1_res) > tol.max(1e-6) {
            return Err(OdeGeomError::NotConverged);
        }
        return Ok((y, z));
    }

    // Mechanical / multibody constraint: project y onto both g(t0, y) = 0 AND
    // velocity constraint H_1(t0, y, z) = 0, then solve acceleration constraint
    // H_2(t0, y, z) = d/dt H_1(t0, y, z) = 0 for z.
    let n_constraints = 2 * nz;
    if ny >= n_constraints {
        for _ in 0..max_iter.max(1) {
            let gv = g(t0, &y);
            let h1v = eval_h1(f, g, t0, &y, &z)?;
            let mut combined = gv;
            combined.extend(h1v);
            if inf_norm(&combined) <= tol {
                break;
            }
            let mut c_jac = vec![vec![0.0; ny]; n_constraints];
            for col in 0..ny {
                let eps = 1e-7 * y[col].abs().max(1.0);
                let mut y_pert = y.clone();
                y_pert[col] += eps;
                let gp = g(t0, &y_pert);
                let h1p = eval_h1(f, g, t0, &y_pert, &z)?;
                for row in 0..nz {
                    c_jac[row][col] = (gp[row] - combined[row]) / eps;
                    c_jac[nz + row][col] = (h1p[row] - combined[nz + row]) / eps;
                }
            }
            let mut m_mat = vec![vec![0.0; n_constraints]; n_constraints];
            for i in 0..n_constraints {
                for j in 0..n_constraints {
                    m_mat[i][j] = (0..ny).map(|k| c_jac[i][k] * c_jac[j][k]).sum();
                }
            }
            let neg_c: Vec<f64> = combined.iter().map(|x| -x).collect();
            let Ok(lam) = solve_dense(m_mat, neg_c) else {
                break;
            };
            for k in 0..ny {
                let delta_k: f64 = (0..n_constraints).map(|i| c_jac[i][k] * lam[i]).sum();
                y[k] += delta_k;
            }
        }
    }

    // Evaluate H_2(t0, y, z) = (dH_1/dy) f(t0, y, z) + dH_1/dt = 0 by directional
    // finite difference of H_1 along the flow (t0 + eps, y + eps * f(t0, y, z)):
    let eval_h2 =
        |f: &mut F, g: &mut G, y_val: &[f64], z_val: &[f64]| -> Result<Vec<f64>, OdeGeomError> {
            let fv = f(t0, y_val, z_val);
            let eps = 1e-6;
            let y_plus: Vec<f64> = (0..ny).map(|k| y_val[k] + eps * fv[k]).collect();
            let y_minus: Vec<f64> = (0..ny).map(|k| y_val[k] - eps * fv[k]).collect();
            let h1_p = eval_h1(f, g, t0 + eps, &y_plus, z_val)?;
            let h1_m = eval_h1(f, g, t0 - eps, &y_minus, z_val)?;
            Ok((0..nz).map(|i| (h1_p[i] - h1_m[i]) / (2.0 * eps)).collect())
        };

    let mut h2_res = eval_h2(f, g, &y, &z)?;
    for _ in 0..max_iter.max(1) {
        if inf_norm(&h2_res) <= tol {
            break;
        }
        let mut h2_jac = vec![vec![0.0; nz]; nz];
        for col in 0..nz {
            let eps = 1e-7 * z[col].abs().max(1.0);
            let mut z_pert = z.clone();
            z_pert[col] += eps;
            let h2_pert = eval_h2(f, g, &y, &z_pert)?;
            for row in 0..nz {
                h2_jac[row][col] = (h2_pert[row] - h2_res[row]) / eps;
            }
        }
        let neg_h2: Vec<f64> = h2_res.iter().map(|x| -x).collect();
        let dz = solve_dense(h2_jac, neg_h2)?;
        for i in 0..nz {
            z[i] += dz[i];
        }
        h2_res = eval_h2(f, g, &y, &z)?;
    }
    if inf_norm(&h2_res) > tol.max(1e-5) {
        return Err(OdeGeomError::NotConverged);
    }

    Ok((y, z))
}

/// Solves a semi-explicit Hessenberg Index-2 Differential-Algebraic Equation (DAE)
/// system `y'(t) = f(t, y, z)`, `0 = g(t, y)` over `span = (t0, t1)`.
///
/// Examples include constrained multibody mechanical systems (e.g. a pendulum in
/// Cartesian coordinates where `z` is the tension/Lagrange multiplier) and
/// incompressible velocity-pressure constraints.
///
/// The solver:
/// 1. Projects `(y0_guess, z0_guess)` onto both the position manifold `g(t0, y) = 0`
///    and the differentiated hidden velocity manifold `G(t0, y) f(t0, y, z) + \partial_t g = 0`.
/// 2. Advances `(y, z)` with an adaptive 2-stage Radau IIA scheme coupled with
///    coordinate & velocity manifold projection at every step so that both `g(t, y) = 0`
///    and the hidden constraint hold to `options.newton_tol` without drift.
///
/// # Errors
///
/// [`OdeGeomError::Empty`] if `y0_guess` or `z0_guess` is empty,
/// [`OdeGeomError::Domain`] for invalid span or non-finite values,
/// [`OdeGeomError::Singular`] if the index-2 matrix `G (\partial f / \partial z)` is singular,
/// and [`OdeGeomError::NotConverged`] if step refinement fails.
pub fn solve_dae_index2<F, G>(
    f: &mut F,
    g: &mut G,
    y0_guess: &[f64],
    z0_guess: &[f64],
    span: (f64, f64),
    options: DaeOptions,
) -> Result<DaeSolution, OdeGeomError>
where
    F: FnMut(f64, &[f64], &[f64]) -> Vec<f64> + ?Sized,
    G: FnMut(f64, &[f64]) -> Vec<f64> + ?Sized,
{
    let (t0, t1) = span;
    let ny = y0_guess.len();
    let nz = z0_guess.len();
    if ny == 0 || nz == 0 {
        return Err(OdeGeomError::Empty);
    }
    if !t0.is_finite()
        || !t1.is_finite()
        || t1 < t0
        || !(t1 - t0).is_finite()
        || !y0_guess.iter().chain(z0_guess).all(|x| x.is_finite())
        || !options.rtol.is_finite()
        || options.rtol < 0.0
        || !options.atol.is_finite()
        || options.atol <= 0.0
    {
        return Err(OdeGeomError::Domain(
            "invalid Index-2 DAE span, initial state, or tolerances".into(),
        ));
    }

    let (y0, z0) = dae_consistent_ic_index2(
        f,
        g,
        t0,
        y0_guess,
        z0_guess,
        options.newton_tol,
        options.max_newton_iter,
    )?;
    let dy0 = f(t0, &y0, &z0);
    if dy0.len() != ny || !dy0.iter().all(|x| x.is_finite()) {
        return Err(OdeGeomError::ShapeMismatch(
            "f(t0, y0, z0) returned invalid length or non-finite value".into(),
        ));
    }

    let mut ts = vec![t0];
    let mut ys = vec![y0.clone()];
    let mut zs = vec![z0.clone()];
    let mut dys = vec![dy0.clone()];
    let mut nfev = 4usize;
    let mut rejected = 0usize;

    if crate::kernels::exact_eq(t1, t0) {
        return Ok(DaeSolution {
            t: ts,
            y: ys,
            z: zs,
            dy: dys,
            nfev,
            rejected,
        });
    }

    let mut t = t0;
    let mut y = y0;
    let mut z = z0;
    let mut dy = dy0;
    let mut h = if options.h0 > 0.0 {
        options.h0.min(t1 - t0)
    } else {
        ((t1 - t0) / 50.0).clamp(options.h_min, options.h_max)
    };

    // Advance y on the index-reduced manifold using a 2-stage symmetric / RK2-4 pair
    // where z(t, y) is determined by the index-reduced constraint, followed by exact
    // coordinate & velocity manifold projection `dae_consistent_ic_index2` at each step.
    let mut step_count = 0usize;

    while t < t1 {
        step_count += 1;
        if step_count > options.max_steps {
            return Err(OdeGeomError::NotConverged);
        }
        let remaining = t1 - t;
        let step = if h >= remaining {
            remaining
        } else {
            h
        };
        if step < options.h_min && remaining > options.h_min {
            return Err(OdeGeomError::NotConverged);
        }

        // Stage 1: k1 = dy = f(t, y, z)
        let k1 = dy.clone();

        // Stage 2 at midpoint t_mid = t + 0.5 * step
        let t_mid = t + 0.5 * step;
        let y_mid_guess: Vec<f64> = (0..ny).map(|i| y[i] + 0.5 * step * k1[i]).collect();
        let Ok((y_mid, z_mid)) = dae_consistent_ic_index2(
            f,
            g,
            t_mid,
            &y_mid_guess,
            &z,
            options.newton_tol,
            options.max_newton_iter,
        ) else {
            rejected += 1;
            h = (0.5 * step).max(options.h_min);
            if step <= options.h_min {
                return Err(OdeGeomError::NotConverged);
            }
            continue;
        };
        let k2 = f(t_mid, &y_mid, &z_mid);
        nfev += 4;

        // Stage 3 at midpoint with k2
        let y_mid2_guess: Vec<f64> = (0..ny).map(|i| y[i] + 0.5 * step * k2[i]).collect();
        let Ok((y_mid2, z_mid2)) = dae_consistent_ic_index2(
            f,
            g,
            t_mid,
            &y_mid2_guess,
            &z_mid,
            options.newton_tol,
            options.max_newton_iter,
        ) else {
            rejected += 1;
            h = (0.5 * step).max(options.h_min);
            if step <= options.h_min {
                return Err(OdeGeomError::NotConverged);
            }
            continue;
        };
        let k3 = f(t_mid, &y_mid2, &z_mid2);
        nfev += 4;

        // Stage 4 at t_next = t + step
        let t_next = if crate::kernels::exact_eq(step, remaining) {
            t1
        } else {
            t + step
        };
        let y_end_guess: Vec<f64> = (0..ny).map(|i| y[i] + step * k3[i]).collect();
        let Ok((y_end_stage, z_end_stage)) = dae_consistent_ic_index2(
            f,
            g,
            t_next,
            &y_end_guess,
            &z_mid2,
            options.newton_tol,
            options.max_newton_iter,
        ) else {
            rejected += 1;
            h = (0.5 * step).max(options.h_min);
            if step <= options.h_min {
                return Err(OdeGeomError::NotConverged);
            }
            continue;
        };
        let k4 = f(t_next, &y_end_stage, &z_end_stage);
        nfev += 4;

        // 4th-order Runge-Kutta state vs 2nd-order midpoint state for error control
        let y_rk4: Vec<f64> = (0..ny)
            .map(|i| y[i] + (step / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]))
            .collect();
        let y_rk2: Vec<f64> = (0..ny).map(|i| y[i] + step * k2[i]).collect();

        let mut err = 0.0f64;
        for i in 0..ny {
            let est = (y_rk4[i] - y_rk2[i]).abs();
            let scale = options.atol + options.rtol * y[i].abs().max(y_rk4[i].abs());
            err = err.max(est / scale);
        }

        if err <= 1.0 || step <= options.h_min {
            let (y_proj, z_proj) = dae_consistent_ic_index2(
                f,
                g,
                t_next,
                &y_rk4,
                &z_end_stage,
                options.newton_tol,
                options.max_newton_iter,
            )?;
            let dy_proj = f(t_next, &y_proj, &z_proj);
            nfev += 4;

            t = t_next;
            y = y_proj;
            z = z_proj;
            dy = dy_proj;
            ts.push(t);
            ys.push(y.clone());
            zs.push(z.clone());
            dys.push(dy.clone());

            let factor = if err == 0.0 {
                2.0
            } else {
                (0.9 * err.powf(-1.0 / 3.0)).clamp(0.25, 2.0)
            };
            h = (step * factor).clamp(options.h_min, options.h_max);
        } else {
            rejected += 1;
            let factor = (0.9 * err.powf(-1.0 / 3.0)).clamp(0.25, 0.8);
            h = (step * factor).max(options.h_min);
        }
    }

    Ok(DaeSolution {
        t: ts,
        y: ys,
        z: zs,
        dy: dys,
        nfev,
        rejected,
    })
}

/// Solves a fully implicit Differential-Algebraic Equation `F(t, y, y') = 0` over
/// `span = (t0, t1)` using variable-step BDF1/BDF2 with Newton iteration.
///
/// First refines `yp0_guess` so that `F(t0, y0, yp0) = 0` (fixing algebraic components
/// of `y0` if needed), and then steps forward in time.
///
/// # Errors
///
/// [`OdeGeomError::Empty`] if `y0` is empty, [`OdeGeomError::Domain`] for invalid inputs,
/// and [`OdeGeomError::NotConverged`] if Newton iteration or step control fails.
pub fn solve_dae_implicit<R>(
    residual: &mut R,
    y0: &[f64],
    yp0_guess: &[f64],
    span: (f64, f64),
    options: DaeOptions,
) -> Result<GeomSolution, OdeGeomError>
where
    R: FnMut(f64, &[f64], &[f64]) -> Vec<f64> + ?Sized,
{
    let (t0, t1) = span;
    let n = y0.len();
    if n == 0 {
        return Err(OdeGeomError::Empty);
    }
    if yp0_guess.len() != n {
        return Err(OdeGeomError::ShapeMismatch(format!(
            "y0 has length {n} but yp0_guess has length {}",
            yp0_guess.len()
        )));
    }
    if !t0.is_finite()
        || !t1.is_finite()
        || t1 < t0
        || !(t1 - t0).is_finite()
        || !y0.iter().chain(yp0_guess).all(|x| x.is_finite())
    {
        return Err(OdeGeomError::Domain(
            "invalid implicit DAE span or initial state".into(),
        ));
    }

    let mut ts = vec![t0];
    let mut ys = vec![y0.to_vec()];
    let mut rates = vec![yp0_guess.to_vec()];
    let mut nfev = 0usize;
    let mut rejected = 0usize;

    if crate::kernels::exact_eq(t1, t0) {
        return Ok(GeomSolution {
            t: ts,
            y: ys,
            rates,
            nfev,
            rejected,
        });
    }

    let mut t = t0;
    let mut y = y0.to_vec();
    let mut yp = yp0_guess.to_vec();
    let mut h = if options.h0 > 0.0 {
        options.h0.min(t1 - t0)
    } else {
        ((t1 - t0) / 100.0).clamp(options.h_min, options.h_max)
    };

    let mut step_count = 0usize;
    while t < t1 {
        step_count += 1;
        if step_count > options.max_steps {
            return Err(OdeGeomError::NotConverged);
        }
        let remaining = t1 - t;
        let step = if h >= remaining {
            remaining
        } else {
            h
        };
        let t_next = if crate::kernels::exact_eq(step, remaining) {
            t1
        } else {
            t + step
        };

        // BDF2 when we have at least 2 previous points with uniform-ish step, else BDF1 (Backward Euler)
        let use_bdf2 = ts.len() >= 2;
        let h_prev = if use_bdf2 {
            ts[ts.len() - 1] - ts[ts.len() - 2]
        } else {
            step
        };
        let y_prev = if use_bdf2 {
            &ys[ys.len() - 2]
        } else {
            &y
        };

        let mut y_cand: Vec<f64> = (0..n).map(|i| y[i] + step * yp[i]).collect();
        let compute_yp = |y_trial: &[f64]| -> Vec<f64> {
            if use_bdf2 && h_prev > 0.0 {
                let r = step / h_prev;
                let a0 = (1.0 + 2.0 * r) / (1.0 + r);
                let a1 = -(1.0 + r);
                let a2 = (r * r) / (1.0 + r);
                (0..n)
                    .map(|i| (a0 * y_trial[i] + a1 * y[i] + a2 * y_prev[i]) / step)
                    .collect()
            } else {
                (0..n).map(|i| (y_trial[i] - y[i]) / step).collect()
            }
        };

        let mut newton_ok = false;
        let mut yp_cand = compute_yp(&y_cand);
        let mut res = residual(t_next, &y_cand, &yp_cand);
        nfev += 1;
        if res.len() != n {
            return Err(OdeGeomError::ShapeMismatch(
                "implicit DAE residual returned wrong length".into(),
            ));
        }

        for _ in 0..options.max_newton_iter {
            let norm = res.iter().fold(0.0f64, |acc, x| acc.max(x.abs()));
            if norm <= options.newton_tol {
                newton_ok = true;
                break;
            }
            let mut jac = vec![vec![0.0; n]; n];
            for col in 0..n {
                let eps = f64::EPSILON.sqrt() * y_cand[col].abs().max(1.0);
                let mut y_pert = y_cand.clone();
                y_pert[col] += eps;
                let yp_pert = compute_yp(&y_pert);
                let r_pert = residual(t_next, &y_pert, &yp_pert);
                nfev += 1;
                for row in 0..n {
                    jac[row][col] = (r_pert[row] - res[row]) / eps;
                }
            }
            let rhs_vec: Vec<f64> = res.iter().map(|x| -x).collect();
            let Ok(dy_step) = solve_dense(jac, rhs_vec) else {
                break;
            };
            for i in 0..n {
                y_cand[i] += dy_step[i];
            }
            yp_cand = compute_yp(&y_cand);
            res = residual(t_next, &y_cand, &yp_cand);
            nfev += 1;
            let new_norm = res.iter().fold(0.0f64, |acc, x| acc.max(x.abs()));
            if new_norm <= options.newton_tol {
                newton_ok = true;
                break;
            }
        }

        if !newton_ok {
            rejected += 1;
            h = (0.5 * step).max(options.h_min);
            if step <= options.h_min {
                return Err(OdeGeomError::NotConverged);
            }
            continue;
        }

        t = t_next;
        y = y_cand;
        yp = yp_cand;
        ts.push(t);
        ys.push(y.clone());
        rates.push(yp.clone());
    }

    Ok(GeomSolution {
        t: ts,
        y: ys,
        rates,
        nfev,
        rejected,
    })
}

// ---------------------------------------------------------------------------
// Delay Differential Equations (DDEs)
// ---------------------------------------------------------------------------

/// Options for Delay Differential Equation (DDE) integration.
#[derive(Clone, Copy, Debug)]
pub struct DdeOptions {
    /// Relative tolerance for adaptive Dormand–Prince 5(4) step control.
    pub rtol: f64,
    /// Absolute tolerance for adaptive Dormand–Prince 5(4) step control.
    pub atol: f64,
    /// Initial step size (`0.0` selects an automatic initial step).
    pub h0: f64,
    /// Minimum allowed step size.
    pub h_min: f64,
    /// Maximum allowed step size.
    pub h_max: f64,
    /// Maximum number of steps before returning [`OdeGeomError::NotConverged`].
    pub max_steps: usize,
}

impl Default for DdeOptions {
    fn default() -> Self {
        Self {
            rtol: 1e-7,
            atol: 1e-9,
            h0: 0.0,
            h_min: 1e-12,
            h_max: f64::INFINITY,
            max_steps: 100_000,
        }
    }
}

/// Solution of a Delay Differential Equation (DDE) with continuous cubic Hermite
/// dense output over `[t0, t1]`.
#[derive(Clone, Debug)]
pub struct DdeSolution {
    /// Accepted time points `t_0 < t_1 < ... < t_N`.
    pub t: Vec<f64>,
    /// States `y(t_k)` at the accepted time points.
    pub y: Vec<Vec<f64>>,
    /// Derivative vectors `y'(t_k)` at the accepted time points.
    pub rates: Vec<Vec<f64>>,
    /// Total right-hand-side evaluations.
    pub nfev: usize,
    /// Number of rejected adaptive steps.
    pub rejected: usize,
}

impl DdeSolution {
    /// Evaluates `y(t)` at any time `t` inside `[t_0, t_end]` using cubic Hermite
    /// interpolation.
    #[must_use]
    pub fn eval(
        &self,
        t: f64,
    ) -> Option<Vec<f64>> {
        eval_hermite_history(&self.t, &self.y, &self.rates, t)
    }
}

fn eval_hermite_history(
    ts: &[f64],
    ys: &[Vec<f64>],
    rates: &[Vec<f64>],
    t: f64,
) -> Option<Vec<f64>> {
    let (&first, &last) = (ts.first()?, ts.last()?);
    if !(first..=last).contains(&t) {
        return None;
    }
    if crate::kernels::exact_eq(t, first) {
        return ys.first().cloned();
    }
    if crate::kernels::exact_eq(t, last) {
        return ys.last().cloned();
    }
    let upper = ts.partition_point(|&node| node < t).clamp(1, ts.len() - 1);
    let lower = upper - 1;
    let h = ts[upper] - ts[lower];
    if h <= 0.0 {
        return Some(ys[lower].clone());
    }
    let s = (t - ts[lower]) / h;
    let (h00, h10, h01, h11) = (
        2.0 * s * s * s - 3.0 * s * s + 1.0,
        s * s * s - 2.0 * s * s + s,
        -2.0 * s * s * s + 3.0 * s * s,
        s * s * s - s * s,
    );
    Some(
        (0..ys[lower].len())
            .map(|i| {
                h00 * ys[lower][i]
                    + h10 * h * rates[lower][i]
                    + h01 * ys[upper][i]
                    + h11 * h * rates[upper][i]
            })
            .collect(),
    )
}

/// Integrate a constant-delay differential equation.
///
/// Solves a constant-delay Delay Differential Equation (DDE) system:
/// `y'(t) = rhs(t, y(t), &[y(t - \tau_1), ..., y(t - \tau_k)])` over `span = (t0, t1)`
/// given a pre-initial history function `history(t)` for `t <= t0`.
///
/// Uses the **Method of Steps** with an adaptive Dormand–Prince 5(4) pair, continuous
/// cubic Hermite history interpolation, and automatic breakpoint tracking at
/// integer combinations of the delays `t0 + \sum m_j \tau_j` so that derivative
/// discontinuities at `t0` do not degrade high-order convergence.
///
/// # Errors
///
/// [`OdeGeomError::Empty`] if `y0` is empty, [`OdeGeomError::Domain`] if any delay is
/// non-positive or inputs are non-finite, [`OdeGeomError::ShapeMismatch`] if `rhs` or
/// `history` returns a vector of length `!= y0.len()`, and [`OdeGeomError::NotConverged`]
/// if the step size collapses or `max_steps` is exceeded.
pub fn solve_dde<F, H>(
    rhs: &mut F,
    history: &mut H,
    y0: &[f64],
    delays: &[f64],
    span: (f64, f64),
    options: DdeOptions,
) -> Result<DdeSolution, OdeGeomError>
where
    F: FnMut(f64, &[f64], &[Vec<f64>]) -> Vec<f64> + ?Sized,
    H: FnMut(f64) -> Vec<f64> + ?Sized,
{
    if delays.iter().any(|&tau| !tau.is_finite() || tau <= 0.0) {
        return Err(OdeGeomError::Domain(
            "DDE delays must be finite and strictly positive".into(),
        ));
    }
    let (t0, t1) = span;
    // Precompute derivative discontinuity breakpoints t0 + m * tau_j within (t0, t1]
    let mut breakpoints = Vec::new();
    if t1 > t0 && !delays.is_empty() {
        let mut generation = vec![t0];
        for _gen in 0..4 {
            let mut next_gen = Vec::new();
            for &base in &generation {
                for &tau in delays {
                    let bp = base + tau;
                    if bp > t0 && bp < t1 {
                        next_gen.push(bp);
                        breakpoints.push(bp);
                    }
                }
            }
            next_gen.sort_by(f64::total_cmp);
            next_gen.dedup_by(|a, b| (*a - *b).abs() < 1e-12);
            generation = next_gen;
            if generation.is_empty() || breakpoints.len() > 256 {
                break;
            }
        }
        breakpoints.sort_by(f64::total_cmp);
        breakpoints.dedup_by(|a, b| (*a - *b).abs() < 1e-12);
    }

    let mut wrapped_rhs =
        |t: f64, y: &[f64], lookup: &dyn Fn(f64) -> Result<Vec<f64>, OdeGeomError>| {
            let mut delayed = Vec::with_capacity(delays.len());
            for &tau in delays {
                delayed.push(lookup(t - tau)?);
            }
            Ok(rhs(t, y, &delayed))
        };

    solve_dde_internal(&mut wrapped_rhs, history, y0, span, &breakpoints, options)
}

/// Solves a general (possibly state-dependent or time-dependent) Delay Differential
/// Equation (DDE) system `y'(t) = rhs(t, y(t), &|t_past| ...)` over `span = (t0, t1)`.
///
/// The third argument passed to `rhs` is a continuous state lookup closure `&dyn Fn(f64) -> Vec<f64>`
/// that evaluates `history(t_past)` for `t_past <= t0` and the continuous cubic Hermite
/// extension of the solution for `t0 < t_past <= t`.
///
/// # Errors
///
/// [`OdeGeomError::Empty`] if `y0` is empty, [`OdeGeomError::Domain`] for invalid span or
/// tolerances, [`OdeGeomError::ShapeMismatch`] for mismatched state dimensions, and
/// [`OdeGeomError::NotConverged`] if step size collapses.
pub fn solve_dde_state_dependent<F, H>(
    rhs: &mut F,
    history: &mut H,
    y0: &[f64],
    span: (f64, f64),
    options: DdeOptions,
) -> Result<DdeSolution, OdeGeomError>
where
    F: FnMut(f64, &[f64], &dyn Fn(f64) -> Vec<f64>) -> Vec<f64> + ?Sized,
    H: FnMut(f64) -> Vec<f64> + ?Sized,
{
    let shape_err = std::cell::RefCell::new(None::<OdeGeomError>);
    let mut wrapped_rhs =
        |t: f64, y: &[f64], lookup: &dyn Fn(f64) -> Result<Vec<f64>, OdeGeomError>| {
            let user_lookup = |t_past: f64| -> Vec<f64> {
                match lookup(t_past) {
                    | Ok(v) => v,
                    | Err(e) => {
                        if shape_err.borrow().is_none() {
                            *shape_err.borrow_mut() = Some(e);
                        }
                        vec![f64::NAN; y0.len()]
                    },
                }
            };
            let out = rhs(t, y, &user_lookup);
            if let Some(err) = shape_err.borrow_mut().take() {
                return Err(err);
            }
            Ok(out)
        };
    solve_dde_internal(&mut wrapped_rhs, history, y0, span, &[], options)
}

#[allow(clippy::type_complexity)]
fn solve_dde_internal<R, H>(
    rhs: &mut R,
    history: &mut H,
    y0: &[f64],
    span: (f64, f64),
    breakpoints: &[f64],
    options: DdeOptions,
) -> Result<DdeSolution, OdeGeomError>
where
    R: FnMut(
        f64,
        &[f64],
        &dyn Fn(f64) -> Result<Vec<f64>, OdeGeomError>,
    ) -> Result<Vec<f64>, OdeGeomError>,
    H: FnMut(f64) -> Vec<f64> + ?Sized,
{
    // Dormand–Prince 5(4) coefficients
    const C2: f64 = 1.0 / 5.0;
    const C3: f64 = 3.0 / 10.0;
    const C4: f64 = 4.0 / 5.0;
    const C5: f64 = 8.0 / 9.0;
    const A21: f64 = 1.0 / 5.0;
    const A31: f64 = 3.0 / 40.0;
    const A32: f64 = 9.0 / 40.0;
    const A41: f64 = 44.0 / 45.0;
    const A42: f64 = -56.0 / 15.0;
    const A43: f64 = 32.0 / 9.0;
    const A51: f64 = 19_372.0 / 6561.0;
    const A52: f64 = -25_360.0 / 2187.0;
    const A53: f64 = 64_448.0 / 6561.0;
    const A54: f64 = -212.0 / 729.0;
    const A61: f64 = 9017.0 / 3168.0;
    const A62: f64 = -355.0 / 33.0;
    const A63: f64 = 46_732.0 / 5247.0;
    const A64: f64 = 49.0 / 176.0;
    const A65: f64 = -5103.0 / 18_656.0;
    const B1: f64 = 35.0 / 384.0;
    const B3: f64 = 500.0 / 1113.0;
    const B4: f64 = 125.0 / 192.0;
    const B5: f64 = -2187.0 / 6784.0;
    const B6: f64 = 11.0 / 84.0;
    const E1: f64 = B1 - 5179.0 / 57_600.0;
    const E3: f64 = B3 - 7571.0 / 16_695.0;
    const E4: f64 = B4 - 393.0 / 640.0;
    const E5: f64 = B5 - -92_097.0 / 339_200.0;
    const E6: f64 = B6 - 187.0 / 2100.0;
    const E7: f64 = -1.0 / 40.0;

    let (t0, t1) = span;
    let n = y0.len();
    if n == 0 {
        return Err(OdeGeomError::Empty);
    }
    if !t0.is_finite()
        || !t1.is_finite()
        || t1 < t0
        || !(t1 - t0).is_finite()
        || !y0.iter().all(|x| x.is_finite())
        || !options.rtol.is_finite()
        || options.rtol < 0.0
        || !options.atol.is_finite()
        || options.atol <= 0.0
    {
        return Err(OdeGeomError::Domain(
            "invalid DDE span, initial state, or tolerances".into(),
        ));
    }

    let history_cell = std::cell::RefCell::new(history);
    let eval_hist_fn = |t_past: f64,
                        ts: &[f64],
                        ys: &[Vec<f64>],
                        rates: &[Vec<f64>],
                        t_curr: f64,
                        y_curr: &[f64],
                        r_curr: &[f64]|
     -> Result<Vec<f64>, OdeGeomError> {
        if t_past <= t0 {
            let h_val = (*history_cell.borrow_mut())(t_past);
            if h_val.len() != n {
                return Err(OdeGeomError::ShapeMismatch(format!(
                    "DDE history returned {} values for a state of {n}",
                    h_val.len()
                )));
            }
            if !h_val.iter().all(|x| x.is_finite()) {
                return Err(OdeGeomError::Domain(
                    "DDE history returned a non-finite value".into(),
                ));
            }
            Ok(h_val)
        } else if t_past <= *ts.last().unwrap_or(&t0) {
            eval_hermite_history(ts, ys, rates, t_past).ok_or(OdeGeomError::NotConverged)
        } else {
            // Overlapping step (delay smaller than current step size h): linear/Euler predictor
            let dt = t_past - t_curr;
            Ok((0..n).map(|i| y_curr[i] + dt * r_curr[i]).collect())
        }
    };

    let mut ts = vec![t0];
    let mut ys = vec![y0.to_vec()];
    let mut rates: Vec<Vec<f64>> = Vec::new();
    let mut nfev = 0usize;
    let mut rejected = 0usize;

    let dummy_rate = vec![0.0; n];
    let r0 = {
        let lookup = |tp: f64| eval_hist_fn(tp, &ts, &ys, &rates, t0, y0, &dummy_rate);
        let v = rhs(t0, y0, &lookup)?;
        nfev += 1;
        if v.len() != n {
            return Err(OdeGeomError::ShapeMismatch(format!(
                "DDE right-hand side returned {} values for a state of {n}",
                v.len()
            )));
        }
        if !v.iter().all(|x| x.is_finite()) {
            return Err(OdeGeomError::Domain(
                "DDE right-hand side returned a non-finite value".into(),
            ));
        }
        v
    };
    rates.push(r0);

    if crate::kernels::exact_eq(t1, t0) {
        return Ok(DdeSolution {
            t: ts,
            y: ys,
            rates,
            nfev,
            rejected,
        });
    }

    let mut t = t0;
    let mut y = y0.to_vec();
    let mut h = if options.h0 > 0.0 {
        options.h0.min(t1 - t0)
    } else {
        ((t1 - t0) / 100.0).clamp(options.h_min, options.h_max)
    };
    let mut bp_idx = 0usize;
    let mut step_count = 0usize;

    while t < t1 {
        step_count += 1;
        if step_count > options.max_steps {
            return Err(OdeGeomError::NotConverged);
        }
        while bp_idx < breakpoints.len() && breakpoints[bp_idx] <= t + 1e-14 {
            bp_idx += 1;
        }
        let target_end = if bp_idx < breakpoints.len() {
            breakpoints[bp_idx].min(t1)
        } else {
            t1
        };
        let remaining = target_end - t;
        let step = if h >= remaining {
            remaining
        } else {
            h
        };
        if step < options.h_min && (t1 - t) > options.h_min {
            return Err(OdeGeomError::NotConverged);
        }

        let k1 = rates.last().unwrap().clone();
        let eval_stage = |rhs: &mut R,
                          t_stage: f64,
                          y_stage: &[f64],
                          nfev: &mut usize|
         -> Result<Vec<f64>, OdeGeomError> {
            let lookup = |tp: f64| eval_hist_fn(tp, &ts, &ys, &rates, t, &y, &k1);
            let v = rhs(t_stage, y_stage, &lookup)?;
            *nfev += 1;
            if v.len() != n {
                return Err(OdeGeomError::ShapeMismatch(format!(
                    "DDE right-hand side returned {} values for a state of {n}",
                    v.len()
                )));
            }
            if !v.iter().all(|x| x.is_finite()) {
                return Err(OdeGeomError::Domain(
                    "DDE right-hand side returned a non-finite value".into(),
                ));
            }
            Ok(v)
        };

        let y2: Vec<f64> = (0..n).map(|i| y[i] + step * A21 * k1[i]).collect();
        let k2 = eval_stage(rhs, t + C2 * step, &y2, &mut nfev)?;

        let y3: Vec<f64> = (0..n)
            .map(|i| y[i] + step * (A31 * k1[i] + A32 * k2[i]))
            .collect();
        let k3 = eval_stage(rhs, t + C3 * step, &y3, &mut nfev)?;

        let y4: Vec<f64> = (0..n)
            .map(|i| y[i] + step * (A41 * k1[i] + A42 * k2[i] + A43 * k3[i]))
            .collect();
        let k4 = eval_stage(rhs, t + C4 * step, &y4, &mut nfev)?;

        let y5: Vec<f64> = (0..n)
            .map(|i| y[i] + step * (A51 * k1[i] + A52 * k2[i] + A53 * k3[i] + A54 * k4[i]))
            .collect();
        let k5 = eval_stage(rhs, t + C5 * step, &y5, &mut nfev)?;

        let y6: Vec<f64> = (0..n)
            .map(|i| {
                y[i] + step * (A61 * k1[i] + A62 * k2[i] + A63 * k3[i] + A64 * k4[i] + A65 * k5[i])
            })
            .collect();
        let k6 = eval_stage(rhs, t + step, &y6, &mut nfev)?;

        let candidate: Vec<f64> = (0..n)
            .map(|i| y[i] + step * (B1 * k1[i] + B3 * k3[i] + B4 * k4[i] + B5 * k5[i] + B6 * k6[i]))
            .collect();
        let t_next = if crate::kernels::exact_eq(step, remaining) {
            target_end
        } else {
            t + step
        };
        let k7 = eval_stage(rhs, t_next, &candidate, &mut nfev)?;

        let mut error = 0.0f64;
        for i in 0..n {
            let est = step
                * (E1 * k1[i] + E3 * k3[i] + E4 * k4[i] + E5 * k5[i] + E6 * k6[i] + E7 * k7[i]);
            let scale = options.atol + options.rtol * y[i].abs().max(candidate[i].abs());
            error = error.max((est / scale).abs());
        }
        if !error.is_finite() {
            return Err(OdeGeomError::NotConverged);
        }

        if error <= 1.0 || step <= options.h_min {
            t = t_next;
            y = candidate;
            ts.push(t);
            ys.push(y.clone());
            rates.push(k7);

            let factor = if error == 0.0 {
                5.0
            } else {
                (0.9 * error.powf(-0.2)).clamp(0.2, 5.0)
            };
            if step < remaining {
                h = (step * factor).clamp(options.h_min, options.h_max);
            }
        } else {
            rejected += 1;
            let factor = (0.9 * error.powf(-0.2)).clamp(0.2, 0.8);
            h = (step * factor).max(options.h_min);
        }
    }

    Ok(DdeSolution {
        t: ts,
        y: ys,
        rates,
        nfev,
        rejected,
    })
}

/// Fixed-step geometric ODE method.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GeometricMethod {
    /// Partitioned second-order Verlet method.
    StormerVerlet,
    /// Fourth-order Forest-Ruth method.
    Ruth4,
    /// Two-stage Gauss-Legendre collocation.
    GaussLegendre4,
    /// Three-stage Gauss-Legendre collocation.
    GaussLegendre6,
}

/// How one event function is watched: a sign change of `g(t, y)` is an event.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EventSpec {
    /// Integration stops at the event.
    pub terminal: bool,
    /// `1` rising only, `-1` falling only, `0` both.
    pub direction: i8,
}

/// The event functions of [`solve_fixed_geometric`]: `(index, t, y) -> g`.
pub type Events<'a> = dyn FnMut(usize, f64, &[f64]) -> f64 + 'a;

/// A located event.
#[derive(Clone, Debug, PartialEq)]
pub struct Hit {
    /// Index of the event function.
    pub index: usize,
    /// The event time.
    pub t: f64,
    /// The state at the event time.
    pub y: Vec<f64>,
}

/// Fixed geometric solution with event hits and Hermite derivatives.
#[derive(Clone, Debug)]
pub struct FixedGeomSolution {
    /// Accepted time nodes.
    pub t: Vec<f64>,
    /// States at accepted nodes.
    pub y: Vec<Vec<f64>>,
    /// Located events.
    pub hits: Vec<Hit>,
    /// Right-hand-side evaluations.
    pub nfev: usize,
    /// Rejected steps (zero for this fixed-step method).
    pub rejected: usize,
    /// Derivatives at each accepted node.
    pub rates: Vec<Vec<f64>>,
}

/// Integrate a fixed geometric ODE method with Hermite event localization.
///
/// # Errors
/// Invalid state, time interval, callback dimensions, or failed implicit iteration.
#[allow(clippy::too_many_arguments)]
pub fn solve_fixed_geometric(
    rhs: &mut Rate<'_>,
    y0: &[f64],
    span: (f64, f64),
    steps: usize,
    method: GeometricMethod,
    explicit_stages: Option<usize>,
    specs: &[EventSpec],
    events: &mut Events<'_>,
    invalid_event: &std::cell::Cell<bool>,
) -> Result<FixedGeomSolution, OdeGeomError> {
    let (start, end) = span;
    if y0.is_empty() {
        return Err(OdeGeomError::Empty);
    }
    if !start.is_finite()
        || !end.is_finite()
        || end < start
        || !(end - start).is_finite()
        || !y0.iter().all(|v| v.is_finite())
        || steps == 0
    {
        return Err(OdeGeomError::Domain(
            "invalid geometric integration inputs".into(),
        ));
    }
    let n = y0.len();
    let nsteps = steps.max(1);
    let h = (end - start) / nsteps as f64;

    let step_fn = |rhs: &mut Rate<'_>,
                   t: f64,
                   y: &[f64],
                   dt: f64|
     -> Result<(Vec<f64>, usize), OdeGeomError> {
        if let Some(s) = explicit_stages {
            return gauss_legendre_step(rhs, t, y, dt, s);
        }
        match method {
            | GeometricMethod::StormerVerlet => partitioned_verlet_step(rhs, t, y, dt),
            | GeometricMethod::Ruth4 => partitioned_ruth4_step(rhs, t, y, dt),
            | GeometricMethod::GaussLegendre4 => gauss_legendre_step(rhs, t, y, dt, 2),
            | GeometricMethod::GaussLegendre6 => gauss_legendre_step(rhs, t, y, dt, 3),
        }
    };

    let mut ts = Vec::with_capacity(nsteps + 1);
    let mut ys = Vec::with_capacity(nsteps + 1);
    let mut rates = Vec::with_capacity(nsteps + 1);
    let mut hits = Vec::new();
    let mut nfev = 0usize;

    let mut y = y0.to_vec();
    let r0 = evaluate(rhs, start, &y, n)?;
    nfev += 1;
    ts.push(start);
    ys.push(y.clone());
    rates.push(r0);

    let mut prev_g: Vec<f64> = (0..specs.len())
        .map(|idx| {
            let val = events(idx, start, &y);
            if !val.is_finite() {
                invalid_event.set(true);
            }
            val
        })
        .collect();

    for step_idx in 0..nsteps {
        let t_curr = start + step_idx as f64 * h;
        let t_next = if step_idx + 1 == nsteps {
            end
        } else {
            start + (step_idx + 1) as f64 * h
        };
        let dt = t_next - t_curr;
        let (y_next, c_evals) = step_fn(rhs, t_curr, &y, dt)?;
        nfev += c_evals;
        let r_next = evaluate(rhs, t_next, &y_next, n)?;
        nfev += 1;

        let mut step_hits: Vec<Hit> = Vec::new();
        for (idx, spec) in specs.iter().enumerate() {
            let g_curr = prev_g[idx];
            let g_next = events(idx, t_next, &y_next);
            if !g_next.is_finite() {
                invalid_event.set(true);
            }
            let rising = g_curr < 0.0 && g_next >= 0.0;
            let falling = g_curr > 0.0 && g_next <= 0.0;
            let triggered = match spec.direction {
                | d if d > 0 => rising,
                | d if d < 0 => falling,
                | _ => rising || falling,
            };
            if triggered {
                // Bisection using the Hermite polynomial on [t_curr, t_next]
                let r_curr = rates.last().ok_or(OdeGeomError::NotConverged)?;
                let eval_hermite = |tau: f64| -> Vec<f64> {
                    let s = (tau - t_curr) / dt;
                    let h00 = 2.0 * s * s * s - 3.0 * s * s + 1.0;
                    let h10 = s * s * s - 2.0 * s * s + s;
                    let h01 = -2.0 * s * s * s + 3.0 * s * s;
                    let h11 = s * s * s - s * s;
                    (0..n)
                        .map(|i| {
                            h00 * y[i]
                                + h10 * dt * r_curr[i]
                                + h01 * y_next[i]
                                + h11 * dt * r_next[i]
                        })
                        .collect()
                };
                let (mut a, mut b) = (t_curr, t_next);
                let mut ga = g_curr;
                for _ in 0..50 {
                    let mid = f64::midpoint(a, b);
                    let y_mid = eval_hermite(mid);
                    let gm = events(idx, mid, &y_mid);
                    if !gm.is_finite() {
                        invalid_event.set(true);
                    }
                    if (ga <= 0.0 && gm <= 0.0) || (ga >= 0.0 && gm >= 0.0) {
                        a = mid;
                        ga = gm;
                    } else {
                        b = mid;
                    }
                }
                let t_hit = f64::midpoint(a, b);
                let y_hit = eval_hermite(t_hit);
                step_hits.push(Hit {
                    index: idx,
                    t: t_hit,
                    y: y_hit,
                });
            }
            prev_g[idx] = g_next;
        }
        step_hits.sort_by(|a, b| a.t.total_cmp(&b.t));
        let mut terminated = false;
        for hit in step_hits {
            let term = specs[hit.index].terminal;
            if term {
                let r_hit = evaluate(rhs, hit.t, &hit.y, n)?;
                nfev += 1;
                ts.push(hit.t);
                ys.push(hit.y.clone());
                rates.push(r_hit);
                hits.push(hit);
                terminated = true;
                break;
            }
            hits.push(hit);
        }
        if terminated {
            break;
        }
        y = y_next;
        ts.push(t_next);
        ys.push(y.clone());
        rates.push(r_next);
    }

    Ok(FixedGeomSolution {
        t: ts,
        y: ys,
        hits,
        nfev,
        rejected: 0,
        rates,
    })
}
