//! State-space models: stability, controllability, observability,
//! discretisation and time-domain simulation.
//!
//! A model is continuous when [`StateSpace::dt`] is `None` and discrete when it
//! holds the sample time. The discretisation uses the zero-order-hold
//! exponential of the augmented matrix `[A B; 0 0]`, so the input matrix is
//! integrated exactly rather than approximated by a Euler step.

use num_complex::Complex64;

use super::linalg;

/// A linear time-invariant state-space model
///
/// ```text
/// x' = A x + B u     (or x[k+1] = A x[k] + B u[k])
/// y  = C x + D u
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct StateSpace {
    /// State matrix `A` (`n x n`).
    pub a: Vec<Vec<f64>>,
    /// Input matrix `B` (`n x m`).
    pub b: Vec<Vec<f64>>,
    /// Output matrix `C` (`p x n`).
    pub c: Vec<Vec<f64>>,
    /// Feedthrough matrix `D` (`p x m`).
    pub d: Vec<Vec<f64>>,
    /// Sample time; `None` marks a continuous-time model.
    pub dt: Option<f64>,
}

impl StateSpace {
    /// Builds a continuous-time model after checking the matrix dimensions.
    ///
    /// # Errors
    /// Returns an error when `a` is not square, when `b`/`c`/`d` are not
    /// consistent with `a`, or when an entry is not finite.
    pub fn new(
        a: Vec<Vec<f64>>,
        b: Vec<Vec<f64>>,
        c: Vec<Vec<f64>>,
        d: Vec<Vec<f64>>,
    ) -> Result<Self, String> {
        let n = a.len();
        if a.iter().any(|row| row.len() != n) {
            return Err("a must be square".to_string());
        }
        if c.iter().any(|row| row.len() != n) {
            return Err("c must have n columns".to_string());
        }
        let p = c.len();
        let m = if n > 0 {
            if b.len() != n {
                return Err("b must have n rows".to_string());
            }
            b.first().map_or(0, Vec::len)
        } else {
            if !b.is_empty() {
                return Err("a model with zero states must have an empty b".to_string());
            }
            d.first().map_or(0, Vec::len)
        };
        if n > 0 && b.iter().any(|row| row.len() != m) {
            return Err("b rows must have equal length".to_string());
        }
        if d.len() != p {
            return Err("d must have as many rows as c".to_string());
        }
        if d.iter().any(|row| row.len() != m) {
            return Err("d must have as many columns as b".to_string());
        }
        if a.iter()
            .chain(&b)
            .chain(&c)
            .chain(&d)
            .flat_map(|row| row.iter())
            .any(|v| !v.is_finite())
        {
            return Err("matrix entries must be finite".to_string());
        }
        Ok(Self { a, b, c, d, dt: None })
    }

    /// Builds a continuous-time model.
    ///
    /// # Errors
    /// As [`StateSpace::new`].
    pub fn continuous(
        a: Vec<Vec<f64>>,
        b: Vec<Vec<f64>>,
        c: Vec<Vec<f64>>,
        d: Vec<Vec<f64>>,
    ) -> Result<Self, String> {
        Self::new(a, b, c, d)
    }

    /// Builds a discrete-time model with sample time `dt`.
    ///
    /// # Errors
    /// As [`StateSpace::new`], plus an error when `dt` is not finite and
    /// positive.
    pub fn discrete(
        a: Vec<Vec<f64>>,
        b: Vec<Vec<f64>>,
        c: Vec<Vec<f64>>,
        d: Vec<Vec<f64>>,
        dt: f64,
    ) -> Result<Self, String> {
        if !(dt.is_finite() && dt > 0.0) {
            return Err("sample time must be positive".to_string());
        }
        let mut model = Self::new(a, b, c, d)?;
        model.dt = Some(dt);
        Ok(model)
    }

    /// Number of states `n`.
    pub fn n_states(&self) -> usize {
        self.a.len()
    }

    /// Number of inputs `m`.
    pub fn n_inputs(&self) -> usize {
        if self.b.is_empty() {
            self.d.first().map_or(0, Vec::len)
        } else {
            self.b.first().map_or(0, Vec::len)
        }
    }

    /// Number of outputs `p`.
    pub fn n_outputs(&self) -> usize {
        self.c.len()
    }

    /// Eigenvalues of `A`.
    ///
    /// Returns an empty vector for a static model or when the eigensolver
    /// fails.
    pub fn poles(&self) -> Vec<Complex64> {
        self.try_poles().unwrap_or_default()
    }

    /// Whether the model is asymptotically stable.
    ///
    /// A continuous model is stable when every pole has negative real part; a
    /// discrete model when every pole lies inside the unit circle. The
    /// predicate is false when a complete finite spectrum is unavailable.
    pub fn is_stable(&self) -> bool {
        self.try_is_stable().unwrap_or(false)
    }

    /// Checks asymptotic stability and reports failure to compute finite poles.
    ///
    /// # Errors
    /// Returns an error if any state lacks a finite pole from the eigensolver.
    pub fn try_is_stable(&self) -> Result<bool, String> {
        Ok(self.poles_are_stable(&self.try_poles()?))
    }

    /// Computes a finite eigenvalue for each state of the model.
    ///
    /// # Errors
    /// Returns an error if the eigensolver fails or a pole is not finite.
    pub fn try_poles(&self) -> Result<Vec<Complex64>, String> {
        if self.n_states() == 0 {
            return Ok(Vec::new());
        }
        let poles = linalg::eigenvalues(&self.a)?;
        if poles.len() != self.n_states() || poles.iter().any(|p| !p.is_finite()) {
            return Err("the state-space eigensolver did not return every finite pole".into());
        }
        Ok(poles)
    }

    fn poles_are_stable(
        &self,
        poles: &[Complex64],
    ) -> bool {
        match self.dt {
            | Some(_) => poles.iter().all(|p| p.norm() < 1.0),
            | None => poles.iter().all(|p| p.re < 0.0),
        }
    }

    /// Controllability matrix `[B, A B, ..., A^{n-1} B]` (`n x (n m)`).
    pub fn controllability_matrix(&self) -> Vec<Vec<f64>> {
        let n = self.n_states();
        let m = self.n_inputs();
        if n == 0 {
            return Vec::new();
        }
        let mut blocks: Vec<Vec<Vec<f64>>> = Vec::with_capacity(n);
        let mut current = self.b.clone();
        blocks.push(current.clone());
        for _ in 1..n {
            current = linalg::matmul_or_zero(&self.a, &current);
            blocks.push(current.clone());
        }
        let mut out = vec![vec![0.0; n * m]; n];
        for (bi, block) in blocks.iter().enumerate() {
            for i in 0..n {
                for j in 0..m {
                    out[i][bi * m + j] =
                        block.get(i).and_then(|r| r.get(j)).copied().unwrap_or(0.0);
                }
            }
        }
        out
    }

    /// Observability matrix `[C; C A; ...; C A^{n-1}]` (`(n p) x n`).
    pub fn observability_matrix(&self) -> Vec<Vec<f64>> {
        let n = self.n_states();
        let p = self.n_outputs();
        if n == 0 {
            return Vec::new();
        }
        let mut blocks: Vec<Vec<Vec<f64>>> = Vec::with_capacity(n);
        let mut current = self.c.clone();
        blocks.push(current.clone());
        for _ in 1..n {
            current = linalg::matmul_or_zero(&current, &self.a);
            blocks.push(current.clone());
        }
        let mut out = vec![vec![0.0; n]; n * p];
        for (bi, block) in blocks.iter().enumerate() {
            for i in 0..p {
                for j in 0..n {
                    out[bi * p + i][j] =
                        block.get(i).and_then(|r| r.get(j)).copied().unwrap_or(0.0);
                }
            }
        }
        out
    }

    /// Whether the pair `(A, B)` is controllable (numerical rank test with a
    /// relative tolerance of `1e-9`).
    pub fn is_controllable(&self) -> bool {
        let n = self.n_states();
        if n == 0 {
            return true;
        }
        linalg::rank_rel(&self.controllability_matrix()).is_ok_and(|r| r == n)
    }

    /// Whether the pair `(A, C)` is observable (numerical rank test with a
    /// relative tolerance of `1e-9`).
    pub fn is_observable(&self) -> bool {
        let n = self.n_states();
        if n == 0 {
            return true;
        }
        linalg::rank_rel(&self.observability_matrix()).is_ok_and(|r| r == n)
    }

    /// Zero-order-hold discretisation with sample time `dt`.
    ///
    /// Computes `exp([A B; 0 0] dt)` by scaling and squaring with a diagonal
    /// Pade approximant of order six. The upper-right block is the exactly
    /// integrated input matrix, not a Euler approximation.
    ///
    /// # Errors
    /// Returns an error when `dt` is not finite and positive, or when the
    /// matrix exponential fails.
    pub fn c2d_zoh(
        &self,
        dt: f64,
    ) -> Result<Self, String> {
        if !(dt.is_finite() && dt > 0.0) {
            return Err("sample time must be positive".to_string());
        }
        let n = self.n_states();
        let m = self.n_inputs();
        if n == 0 {
            return Ok(Self {
                a: Vec::new(),
                b: Vec::new(),
                c: self.c.clone(),
                d: self.d.clone(),
                dt: Some(dt),
            });
        }
        let size = n + m;
        let mut aug = linalg::zeros(size, size);
        for i in 0..n {
            for j in 0..n {
                aug[i][j] = self.a[i][j] * dt;
            }
            for j in 0..m {
                aug[i][n + j] = self.b[i][j] * dt;
            }
        }
        let exp = expm_pade(&aug)?;
        let mut ad = linalg::zeros(n, n);
        let mut bd = linalg::zeros(n, m);
        for i in 0..n {
            for j in 0..n {
                ad[i][j] = exp[i][j];
            }
            for j in 0..m {
                bd[i][j] = exp[i][n + j];
            }
        }
        Ok(Self {
            a: ad,
            b: bd,
            c: self.c.clone(),
            d: self.d.clone(),
            dt: Some(dt),
        })
    }

    /// Simulates a discrete-time model from `x0` under the input sequence `u`.
    ///
    /// The recursion `x[k+1] = A x[k] + B u[k]`, `y[k] = C x[k] + D u[k]` is
    /// applied directly. For continuous-time models, use
    /// [`StateSpace::simulate_continuous`] (with an explicit time grid) or
    /// [`StateSpace::simulate_with_dt`] (with a uniform sample time).
    ///
    /// Returns the output at each input sample.
    ///
    /// # Errors
    /// Returns an error for a wrong `x0` or input length, or when the model has
    /// no sample time.
    pub fn simulate(
        &self,
        x0: &[f64],
        u: &[Vec<f64>],
    ) -> Result<Vec<Vec<f64>>, String> {
        let n = self.n_states();
        let m = self.n_inputs();
        if x0.len() != n {
            return Err("x0 must have n entries".to_string());
        }
        if u.iter().any(|row| row.len() != m) {
            return Err("each input sample must have m entries".to_string());
        }
        if self.dt.is_none() {
            return Err("a sample time is required to simulate a continuous model".to_string());
        }
        let mut x = x0.to_vec();
        let mut next = Vec::with_capacity(n);
        let mut out = Vec::with_capacity(u.len());
        for uk in u {
            let mut y = Vec::with_capacity(self.c.len());
            linalg::sum_matvec_into(&self.c, &x, &self.d, uk, &mut y);
            out.push(y);
            linalg::sum_matvec_into(&self.a, &x, &self.b, uk, &mut next);
            std::mem::swap(&mut x, &mut next);
        }
        Ok(out)
    }

    /// Simulates a continuous-time model on the strictly increasing time grid
    /// `t` from initial state `x0` under piecewise-constant (zero-order-hold)
    /// inputs `u`, where `u[k]` is held across `[t[k], t[k + 1])`.
    ///
    /// Each step uses exact zero-order-hold discretisation via
    /// [`StateSpace::c2d_zoh`], caching `(Ad, Bd)` whenever consecutive step
    /// sizes match to a relative tolerance of `1e-12`. Uniform grids compute
    /// the matrix exponential only once.
    ///
    /// Returns the output `y[k] = C x(t[k]) + D u[k]` at each grid point.
    ///
    /// # Errors
    /// Returns an error when the model is discrete-time, when `x0` or `u`
    /// dimensions are wrong, when `u.len() != t.len()`, when any entry is
    /// non-finite, or when `t` is not strictly increasing.
    pub fn simulate_continuous(
        &self,
        x0: &[f64],
        u: &[Vec<f64>],
        t: &[f64],
    ) -> Result<Vec<Vec<f64>>, String> {
        if self.dt.is_some() {
            return Err("simulate_continuous requires a continuous-time model".to_string());
        }
        let n = self.n_states();
        let m = self.n_inputs();
        if x0.len() != n {
            return Err("x0 must have n entries".to_string());
        }
        if x0.iter().any(|v| !v.is_finite()) {
            return Err("x0 entries must be finite".to_string());
        }
        if u.len() != t.len() {
            return Err("input sequence and time grid must have the same length".to_string());
        }
        if u.iter().any(|row| row.len() != m) {
            return Err("each input sample must have m entries".to_string());
        }
        if u.iter().flatten().any(|v| !v.is_finite()) {
            return Err("input entries must be finite".to_string());
        }
        if t.iter().any(|v| !v.is_finite()) {
            return Err("time grid entries must be finite".to_string());
        }
        for k in 1..t.len() {
            if t[k].partial_cmp(&t[k - 1]) != Some(std::cmp::Ordering::Greater) {
                return Err("time grid must be strictly increasing".to_string());
            }
        }
        let mut x = x0.to_vec();
        let mut next = Vec::with_capacity(n);
        let mut out = Vec::with_capacity(t.len());
        let mut cached: Option<(f64, Vec<Vec<f64>>, Vec<Vec<f64>>)> = None;
        for k in 0..t.len() {
            let uk = &u[k];
            let mut y = Vec::with_capacity(self.c.len());
            linalg::sum_matvec_into(&self.c, &x, &self.d, uk, &mut y);
            out.push(y);
            if k + 1 < t.len() {
                let dt = t[k + 1] - t[k];
                let reuse = cached
                    .as_ref()
                    .is_some_and(|(prev_dt, _, _)| same_sample_time(dt, *prev_dt));
                if !reuse {
                    let disc = self.c2d_zoh(dt)?;
                    cached = Some((dt, disc.a, disc.b));
                }
                if let Some((_, ad, bd)) = cached.as_ref() {
                    linalg::sum_matvec_into(ad, &x, bd, uk, &mut next);
                    std::mem::swap(&mut x, &mut next);
                }
            }
        }
        Ok(out)
    }

    /// Simulates either a continuous-time model (by zero-order-hold
    /// discretisation at sample time `dt`) or a discrete-time model whose
    /// sample time matches `dt` to a relative tolerance of `1e-12`.
    ///
    /// # Errors
    /// Returns an error when `dt` is not finite and positive, when a discrete
    /// model's `dt` does not match the requested `dt`, or when [`StateSpace::c2d_zoh`]
    /// or [`StateSpace::simulate`] fails.
    pub fn simulate_with_dt(
        &self,
        x0: &[f64],
        u: &[Vec<f64>],
        dt: f64,
    ) -> Result<Vec<Vec<f64>>, String> {
        if !(dt.is_finite() && dt > 0.0) {
            return Err("sample time must be positive".to_string());
        }
        if let Some(model_dt) = self.dt {
            if !same_sample_time(model_dt, dt) {
                return Err("sample time does not match discrete model dt".to_string());
            }
            return self.simulate(x0, u);
        }
        let disc = self.c2d_zoh(dt)?;
        disc.simulate(x0, u)
    }
}

// Dividing by the larger duration avoids an absolute tolerance floor and
// keeps the comparison meaningful even for subnormal sample times.
fn same_sample_time(
    left: f64,
    right: f64,
) -> bool {
    left.is_finite()
        && right.is_finite()
        && left > 0.0
        && right > 0.0
        && (left - right).abs() / left.max(right) <= 1e-12
}

/// Coefficients of the order-six Pade approximant of the exponential.
const PADE6: [f64; 7] = [
    1.0,
    0.5,
    5.0 / 44.0,
    1.0 / 66.0,
    1.0 / 792.0,
    1.0 / 15840.0,
    1.0 / 665_280.0,
];

/// Matrix exponential by scaling and squaring with a Pade(6) approximant.
fn expm_pade(m: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, String> {
    let n = m.len();
    if n == 0 {
        return Ok(Vec::new());
    }
    let norm = m
        .iter()
        .map(|row| row.iter().map(|v| v.abs()).sum::<f64>())
        .fold(0.0_f64, f64::max);
    let mut squarings = 0_u32;
    let mut scale = 1.0_f64;
    while norm / scale > 0.5 && squarings < 60 {
        scale *= 2.0;
        squarings += 1;
    }
    let mut scaled = m.to_vec();
    for row in &mut scaled {
        for v in row.iter_mut() {
            *v /= scale;
        }
    }
    let mut u = linalg::identity(n);
    let mut v = linalg::identity(n);
    let mut power = linalg::identity(n);
    for k in 1..=6 {
        power = linalg::matmul(&power, &scaled)?;
        let coef = PADE6[k];
        let sign = if k % 2 == 0 { 1.0 } else { -1.0 };
        for i in 0..n {
            for j in 0..n {
                u[i][j] += coef * power[i][j];
                v[i][j] += sign * coef * power[i][j];
            }
        }
    }
    let mut x = linalg::solve_multi(&v, &u)?;
    for _ in 0..squarings {
        x = linalg::matmul(&x, &x)?;
    }
    Ok(x)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn double_integrator() -> StateSpace {
        StateSpace::continuous(
            vec![vec![0.0, 1.0], vec![0.0, 0.0]],
            vec![vec![0.0], vec![1.0]],
            vec![vec![1.0, 0.0]],
            vec![vec![0.0]],
        )
        .unwrap()
    }

    #[test]
    fn discrete_simulation_keeps_ragged_row_and_nan_semantics() {
        let model = StateSpace {
            a: vec![vec![0.0, 1.0, f64::NAN], vec![2.0]],
            b: vec![vec![2.0, 3.0]],
            c: vec![vec![1.0], vec![0.0, 1.0, f64::NAN], vec![]],
            d: vec![vec![4.0]],
            dt: Some(1.0),
        };
        let output = model
            .simulate(&[1.0, 2.0], &[vec![3.0, 5.0], vec![f64::NAN, 1.0]])
            .unwrap();
        assert_eq!(output[0], [13.0, 2.0, 0.0]);
        assert!(output[1][0].is_nan());
        assert_eq!(&output[1][1..], [2.0, 0.0]);

        let error = model.simulate(&[1.0, 2.0], &[vec![3.0]]).unwrap_err();
        assert_eq!(error, "each input sample must have m entries");
    }

    #[test]
    fn simulation_handles_static_and_autonomous_models() {
        let static_model =
            StateSpace::continuous(vec![], vec![], vec![vec![]], vec![vec![-2.0]]).unwrap();
        let samples = [vec![3.0], vec![-4.0]];
        let times = [0.0, 0.5];
        let expected = [vec![-6.0], vec![8.0]];
        assert_eq!(
            static_model
                .simulate_continuous(&[], &samples, &times)
                .unwrap(),
            expected
        );
        assert_eq!(
            static_model.simulate_with_dt(&[], &samples, 0.5).unwrap(),
            expected
        );

        let autonomous =
            StateSpace::continuous(vec![vec![-1.0]], vec![vec![]], vec![], vec![]).unwrap();
        let samples = [vec![], vec![]];
        let expected: [Vec<f64>; 2] = [vec![], vec![]];
        assert_eq!(
            autonomous
                .simulate_continuous(&[1.0], &samples, &times)
                .unwrap(),
            expected
        );
        assert_eq!(
            autonomous.simulate_with_dt(&[1.0], &samples, 0.5).unwrap(),
            expected
        );
    }

    #[test]
    fn dimensions_are_checked() {
        let bad = StateSpace::continuous(
            vec![vec![0.0, 1.0]],
            vec![vec![0.0]],
            vec![vec![1.0, 0.0]],
            vec![vec![0.0]],
        );
        assert!(bad.is_err());
    }

    #[test]
    fn double_integrator_is_controllable_but_not_stable() {
        let ss = double_integrator();
        assert!(ss.is_controllable());
        assert!(ss.is_observable());
        assert!(!ss.is_stable());
        assert_eq!(ss.n_states(), 2);
        assert_eq!(ss.n_inputs(), 1);
        assert_eq!(ss.n_outputs(), 1);
    }

    #[test]
    fn zoh_matches_analytic_double_integrator() {
        let ss = double_integrator();
        let dt = 0.5;
        let d = ss.c2d_zoh(dt).unwrap();
        // A = [[0,1],[0,0]] => Ad = [[1, dt],[0,1]], Bd = [dt^2/2, dt].
        assert!((d.a[0][0] - 1.0).abs() < 1e-12);
        assert!((d.a[0][1] - dt).abs() < 1e-12);
        assert!((d.a[1][1] - 1.0).abs() < 1e-12);
        assert!((d.b[0][0] - dt * dt / 2.0).abs() < 1e-12);
        assert!((d.b[1][0] - dt).abs() < 1e-12);
        assert_eq!(d.dt, Some(dt));
    }

    #[test]
    fn discrete_simulation_advances_state() {
        let ss = StateSpace::discrete(
            vec![vec![1.0]],
            vec![vec![1.0]],
            vec![vec![1.0]],
            vec![vec![0.0]],
            1.0,
        )
        .unwrap();
        let y = ss
            .simulate(&[0.0], &[vec![1.0], vec![1.0], vec![1.0]])
            .unwrap();
        assert_eq!(y.len(), 3);
        assert!((y[0][0] - 0.0).abs() < 1e-12);
        assert!((y[1][0] - 1.0).abs() < 1e-12);
        assert!((y[2][0] - 2.0).abs() < 1e-12);
    }

    #[test]
    fn continuous_simulation_needs_dt() {
        let ss = double_integrator();
        assert!(ss.simulate(&[0.0, 0.0], &[vec![1.0]]).is_err());
    }

    #[test]
    fn simulate_continuous_matches_rc_decay_and_harmonic_oscillator() {
        // First-order RC decay: x' = -2x, x(0) = 3 => y(t) = 3 exp(-2t).
        let rc = StateSpace::continuous(
            vec![vec![-2.0]],
            vec![vec![1.0]],
            vec![vec![1.0]],
            vec![vec![0.0]],
        )
        .unwrap();
        let t: Vec<f64> = (0..21).map(|i| f64::from(i) * 0.1).collect();
        let u = vec![vec![0.0]; t.len()];
        let y = rc.simulate_continuous(&[3.0], &u, &t).unwrap();
        for (tk, yk) in t.iter().zip(&y) {
            let expected = 3.0 * (-2.0 * tk).exp();
            assert!(
                (yk[0] - expected).abs() < 1e-12,
                "t={tk}: got {}, want {expected}",
                yk[0]
            );
        }

        // Harmonic oscillator: x1' = x2, x2' = -x1, x(0) = [1, 0] => y(t) = cos(t).
        let osc = StateSpace::continuous(
            vec![vec![0.0, 1.0], vec![-1.0, 0.0]],
            vec![vec![0.0], vec![1.0]],
            vec![vec![1.0, 0.0]],
            vec![vec![0.0]],
        )
        .unwrap();
        let y_osc = osc.simulate_continuous(&[1.0, 0.0], &u, &t).unwrap();
        for (tk, yk) in t.iter().zip(&y_osc) {
            assert!(
                (yk[0] - tk.cos()).abs() < 1e-11,
                "t={tk}: got {}, want {}",
                yk[0],
                tk.cos()
            );
        }
    }

    #[test]
    fn subpicosecond_sample_times_are_compared_relatively() {
        // An unforced RC circuit with time constant 1 ps satisfies
        // y(t) = exp(-t / 1 ps). The adjacent intervals here are 0.1 and 0.8 ps.
        let tau = 1e-12;
        let model = StateSpace::continuous(
            vec![vec![-1.0 / tau]],
            vec![vec![1.0]],
            vec![vec![1.0]],
            vec![vec![0.0]],
        )
        .unwrap();
        let times = [0.0, 1e-13, 9e-13];
        let zero_input = vec![vec![0.0]; times.len()];
        let output = model
            .simulate_continuous(&[1.0], &zero_input, &times)
            .unwrap();
        for (&time, sample) in times.iter().zip(&output) {
            let expected = (-time / tau).exp();
            assert!(
                (sample[0] - expected).abs() < 1e-12,
                "t={time}: got {}, want {expected}",
                sample[0]
            );
        }

        let discrete = model.c2d_zoh(times[1]).unwrap();
        let error = discrete
            .simulate_with_dt(&[1.0], &zero_input, times[2] - times[1])
            .unwrap_err();
        assert_eq!(error, "sample time does not match discrete model dt");
        assert_eq!(
            discrete
                .simulate_with_dt(&[1.0], &zero_input, times[1].next_up())
                .unwrap(),
            discrete.simulate(&[1.0], &zero_input).unwrap(),
        );
    }

    #[test]
    fn simulate_continuous_handles_non_uniform_grid_and_step_response() {
        // Double integrator on a non-uniform grid with piecewise-constant input.
        let ss = double_integrator();
        let t = [0.0, 0.5, 1.5, 2.0];
        let u = vec![vec![2.0], vec![-1.0], vec![0.0], vec![0.0]];
        let y = ss.simulate_continuous(&[0.0, 0.0], &u, &t).unwrap();
        assert!((y[0][0] - 0.0).abs() < 1e-12);
        assert!((y[1][0] - 0.25).abs() < 1e-12);
        assert!((y[2][0] - 0.75).abs() < 1e-12);
        assert!((y[3][0] - 0.75).abs() < 1e-12);

        // Compare continuous state-space step response against transfer::step_response.
        let g = crate::kernels::control::transfer::TransferFunction::new(
            vec![1.0],
            vec![1.0, 3.0, 2.0],
        )
        .unwrap();
        let g_ss = g.to_state_space().unwrap();
        let dt = 0.02;
        let grid: Vec<f64> = (0..101).map(|i| f64::from(i) * dt).collect();
        let ones = vec![vec![1.0]; grid.len()];
        let y_cont = g_ss.simulate_continuous(&[0.0, 0.0], &ones, &grid).unwrap();
        let y_dt = g_ss.simulate_with_dt(&[0.0, 0.0], &ones, dt).unwrap();
        let y_tf = crate::kernels::control::transfer::step_response(&g, &grid).unwrap();
        for k in 0..grid.len() {
            assert!((y_cont[k][0] - y_dt[k][0]).abs() < 1e-12);
            assert!(
                (y_cont[k][0] - y_tf[k]).abs() < 1e-6,
                "k={k}: ss={}, tf={}",
                y_cont[k][0],
                y_tf[k]
            );
        }
    }

    #[test]
    fn simulate_continuous_and_with_dt_validate_inputs() {
        let ss = double_integrator();
        // Empty time grid returns empty output.
        assert_eq!(
            ss.simulate_continuous(&[0.0, 0.0], &[], &[]).unwrap(),
            [] as [Vec<f64>; 0],
        );
        // Mismatched lengths, non-increasing times, non-finite inputs.
        assert!(
            ss.simulate_continuous(&[0.0, 0.0], &[vec![1.0]], &[0.0, 1.0])
                .is_err()
        );
        assert!(
            ss.simulate_continuous(&[0.0, 0.0], &[vec![1.0], vec![1.0]], &[1.0, 1.0])
                .is_err()
        );
        assert!(
            ss.simulate_continuous(&[0.0, 0.0], &[vec![1.0], vec![1.0]], &[1.0, 0.5])
                .is_err()
        );
        assert!(
            ss.simulate_continuous(&[f64::NAN, 0.0], &[vec![1.0]], &[0.0])
                .is_err()
        );
        assert!(ss.simulate_with_dt(&[0.0, 0.0], &[vec![1.0]], 0.0).is_err());
        assert!(
            ss.simulate_with_dt(&[0.0, 0.0], &[vec![1.0]], -0.1)
                .is_err()
        );

        // Discrete model rejects simulate_continuous and mismatched dt, accepts matching dt.
        let disc = ss.c2d_zoh(0.5).unwrap();
        assert!(
            disc.simulate_continuous(&[0.0, 0.0], &[vec![1.0]], &[0.0])
                .is_err()
        );
        assert!(
            disc.simulate_with_dt(&[0.0, 0.0], &[vec![1.0]], 0.25)
                .is_err()
        );
        let y_disc = disc
            .simulate_with_dt(&[0.0, 0.0], &[vec![1.0], vec![1.0]], 0.5)
            .unwrap();
        assert!((y_disc[1][0] - 0.125).abs() < 1e-12);
    }
}
