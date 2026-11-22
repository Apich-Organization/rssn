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
        if self.n_states() == 0 {
            return Vec::new();
        }
        linalg::eigenvalues(&self.a).unwrap_or_default()
    }

    /// Whether the model is asymptotically stable.
    ///
    /// A continuous model is stable when every pole has negative real part; a
    /// discrete model when every pole lies inside the unit circle.
    pub fn is_stable(&self) -> bool {
        let poles = self.poles();
        if poles.is_empty() {
            return true;
        }
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

    /// Simulates the model from `x0` under the input sequence `u`.
    ///
    /// The recursion `x[k+1] = A x[k] + B u[k]`, `y[k] = C x[k] + D u[k]` is
    /// applied directly. Continuous-time trajectories are produced by the
    /// time-grid integrators in [`crate::kernels::control::transfer`]; this
    /// routine requires a sample time because the input sequence has no time
    /// base of its own.
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
        let mut out = Vec::with_capacity(u.len());
        for uk in u {
            let mut y = linalg::matvec(&self.c, &x);
            let du = linalg::matvec(&self.d, uk);
            for (i, v) in y.iter_mut().enumerate() {
                *v += du.get(i).copied().unwrap_or(0.0);
            }
            out.push(y);
            let mut next = linalg::matvec(&self.a, &x);
            let bu = linalg::matvec(&self.b, uk);
            for (i, v) in next.iter_mut().enumerate() {
                *v += bu.get(i).copied().unwrap_or(0.0);
            }
            x = next;
        }
        Ok(out)
    }
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
}
