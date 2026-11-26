//! Laurent series via FFT, Padé analytic continuation, Wynn-epsilon sequence
//! acceleration, Borel-Padé resummation, and D-finite (holonomic) ODE analytic
//! continuation with monodromy matrices.

use std::f64::consts::TAU;

use num_complex::Complex64;

use crate::kernels::integrate::gauss_kronrod_any;

const ZERO: Complex64 = Complex64::new(0.0, 0.0);
const ONE: Complex64 = Complex64::new(1.0, 0.0);

/// Extracts Laurent coefficients `c_k` for `k in -neg_order..=pos_order` around
/// `center` on the circle `|z - center| = radius` using a single FFT.
#[must_use]
pub fn laurent_coeffs_fft<F>(
    f: F,
    center: Complex64,
    radius: f64,
    neg_order: usize,
    pos_order: usize,
    n_samples: usize,
) -> Vec<(i32, Complex64)>
where
    F: Fn(Complex64) -> Complex64,
{
    let min_n = (neg_order + pos_order + 2).max(8);
    let n = n_samples.max(min_n);
    let mut buffer = Vec::with_capacity(n);
    for j in 0..n {
        let theta = TAU * (j as f64) / (n as f64);
        let z = center + Complex64::from_polar(radius, theta);
        let val = f(z);
        buffer.push(rustfft::num_complex::Complex::new(val.re, val.im));
    }
    let mut planner = rustfft::FftPlanner::<f64>::new();
    let fft = planner.plan_fft_forward(n);
    fft.process(&mut buffer);

    let nf = n as f64;
    let mut out = Vec::with_capacity(neg_order + pos_order + 1);
    for k in -(neg_order as i32)..=(pos_order as i32) {
        let idx = k.rem_euclid(n as i32) as usize;
        let ck = buffer.get(idx).copied().unwrap_or_default();
        let coeff = Complex64::new(ck.re, ck.im) / (nf * radius.powi(k));
        out.push((k, coeff));
    }
    out
}

/// Solves a dense complex linear system `A x = b` via Gaussian elimination with partial pivoting.
fn solve_complex_linear(
    mut a: Vec<Vec<Complex64>>,
    mut b: Vec<Complex64>,
) -> Option<Vec<Complex64>> {
    let n = b.len();
    if n == 0 || a.len() != n || a.iter().any(|r| r.len() != n) {
        return None;
    }
    for col in 0..n {
        let mut piv = col;
        let mut best = 0.0_f64;
        for (r, row) in a.iter().enumerate().skip(col) {
            let v = row.get(col)?.norm();
            if v > best {
                best = v;
                piv = r;
            }
        }
        if best < 1e-30 || !best.is_finite() {
            return None;
        }
        if piv != col {
            a.swap(col, piv);
            b.swap(col, piv);
        }
        let diag = *a.get(col)?.get(col)?;
        for r in (col + 1)..n {
            let factor = *a.get(r)?.get(col)? / diag;
            if factor == ZERO {
                continue;
            }
            for j in col..n {
                let piv_val = *a.get(col)?.get(j)?;
                if let Some(cell) = a.get_mut(r).and_then(|row| row.get_mut(j)) {
                    *cell -= factor * piv_val;
                }
            }
            let b_col = *b.get(col)?;
            if let Some(slot) = b.get_mut(r) {
                *slot -= factor * b_col;
            }
        }
    }
    let mut x = vec![ZERO; n];
    for i in (0..n).rev() {
        let mut s = *b.get(i)?;
        for j in (i + 1)..n {
            s -= *a.get(i)?.get(j)? * *x.get(j)?;
        }
        let diag = *a.get(i)?.get(i)?;
        if let Some(slot) = x.get_mut(i) {
            *slot = s / diag;
        }
    }
    Some(x)
}

/// Constructs the `[m/n]` complex Padé approximant `P(z) / Q(z)` from Taylor
/// coefficients `c_0, c_1, ..., c_{m+n}` (with `Q(0) = 1`).
#[must_use]
pub fn pade_approximant_complex(
    taylor_coeffs: &[Complex64],
    m: usize,
    n: usize,
) -> Option<(Vec<Complex64>, Vec<Complex64>)> {
    if taylor_coeffs.len() <= m + n {
        return None;
    }
    let mut q = vec![ONE];
    if n > 0 {
        let mut mat = vec![vec![ZERO; n]; n];
        let mut rhs = vec![ZERO; n];
        for i in 0..n {
            let row_idx = m + 1 + i;
            if let Some(slot) = rhs.get_mut(i) {
                *slot = -*taylor_coeffs.get(row_idx)?;
            }
            for j in 0..n {
                let c_idx = row_idx.checked_sub(j + 1);
                if let Some(idx) = c_idx
                    && let Some(cell) = mat.get_mut(i).and_then(|r| r.get_mut(j))
                {
                    *cell = *taylor_coeffs.get(idx)?;
                }
            }
        }
        let Some(b_sol) = solve_complex_linear(mat, rhs) else {
            return pade_approximant_complex(taylor_coeffs, m, n - 1);
        };
        q.extend(b_sol);
    }
    let mut p = vec![ZERO; m + 1];
    for (k, p_slot) in p.iter_mut().enumerate() {
        let mut s = ZERO;
        for j in 0..=k.min(n) {
            s += *q.get(j)? * *taylor_coeffs.get(k - j)?;
        }
        *p_slot = s;
    }
    Some((p, q))
}

/// Evaluates rational function `P(z) / Q(z)` via Horner's method.
#[must_use]
pub fn eval_pade_complex(
    p: &[Complex64],
    q: &[Complex64],
    z: Complex64,
) -> Complex64 {
    let horner = |coeffs: &[Complex64]| -> Complex64 {
        let mut acc = ZERO;
        for &c in coeffs.iter().rev() {
            acc = acc * z + c;
        }
        acc
    };
    horner(p) / horner(q)
}

/// Wynn's epsilon algorithm for accelerating convergence of complex sequences.
#[must_use]
pub fn wynn_epsilon_complex(partial_sums: &[Complex64]) -> Complex64 {
    let n = partial_sums.len();
    if n == 0 {
        return ZERO;
    }
    if n == 1 {
        return partial_sums[0];
    }
    let mut prev = vec![ZERO; n + 1];
    let mut curr = partial_sums.to_vec();
    let mut best = *curr.last().unwrap_or(&ZERO);
    for k in 1..n {
        let len = n - k;
        let mut next = vec![ZERO; len];
        for i in 0..len {
            let diff =
                curr.get(i + 1).copied().unwrap_or(ZERO) - curr.get(i).copied().unwrap_or(ZERO);
            if diff.norm() < 1e-30 {
                return curr.get(i).copied().unwrap_or(best);
            }
            let val = prev.get(i + 1).copied().unwrap_or(ZERO) + ONE / diff;
            if let Some(slot) = next.get_mut(i) {
                *slot = val;
            }
        }
        if k % 2 == 0
            && let Some(&last) = next.last()
            && last.re.is_finite()
            && last.im.is_finite()
        {
            best = last;
        }
        prev = curr;
        curr = next;
    }
    best
}

/// Borel-Padé resummation of a divergent asymptotic series `sum_{k>=0} c_k z^k`
/// via Borel transform `B(t) = sum (c_k / k!) t^k`, Padé `[m/n]` approximant,
/// and Laplace integral `int_0^oo e^(-t) B(z * t) dt`.
#[must_use]
pub fn borel_pade_resum(
    coeffs: &[Complex64],
    z: Complex64,
    pade_m: usize,
    pade_n: usize,
) -> Option<Complex64> {
    if coeffs.len() <= pade_m + pade_n {
        return None;
    }
    let mut borel_coeffs = Vec::with_capacity(coeffs.len());
    let mut fact = 1.0_f64;
    for (k, &ck) in coeffs.iter().enumerate() {
        if k > 0 {
            fact *= k as f64;
        }
        borel_coeffs.push(ck / fact);
    }
    let (p, q) = pade_approximant_complex(&borel_coeffs, pade_m, pade_n)?;
    // Map t in [0, oo) along ray exp(i * ray_angle) to avoid poles on [0, oo).
    let eval_ray = |ray_angle: f64| -> Complex64 {
        let dir = Complex64::from_polar(1.0, ray_angle);
        let eval_laplace = |u: f64| -> Complex64 {
            if u >= 1.0 - 1e-12 {
                return ZERO;
            }
            let one_m_u = 1.0 - u;
            let s = u / one_m_u;
            let t = dir * s;
            let val = eval_pade_complex(&p, &q, z * t);
            if !val.re.is_finite() || !val.im.is_finite() {
                return ZERO;
            }
            let w = (-t).exp() * dir / (one_m_u * one_m_u);
            w * val
        };
        let re = gauss_kronrod_any(|u| eval_laplace(u).re, 0.0, 1.0 - 1e-12, 1e-11, 24).value;
        let im = gauss_kronrod_any(|u| eval_laplace(u).im, 0.0, 1.0 - 1e-12, 1e-11, 24).value;
        Complex64::new(re, im)
    };
    let on_axis = eval_ray(0.0);
    if on_axis.is_finite() {
        return Some(on_axis);
    }
    let lateral = eval_ray(0.2);
    lateral.is_finite().then_some(lateral)
}

/// Weierstrass power-series chain continuation: re-centers a Taylor series
/// `sum c_k (z - a_0)^k` along a sequence of overlapping disk centers `centers = [a_0, a_1, ...]`
/// and returns the Taylor coefficients around the final center `a_last`.
#[must_use]
pub fn weierstrass_chain_continue(
    initial_coeffs: &[Complex64],
    centers: &[Complex64],
) -> Option<Vec<Complex64>> {
    if initial_coeffs.is_empty() || centers.is_empty() {
        return None;
    }
    let mut coeffs = initial_coeffs.to_vec();
    let n = coeffs.len();
    for win in centers.windows(2) {
        let shift = win[1] - win[0];
        let mut next = vec![ZERO; n];
        // Re-center: c_new[j] = sum_{k=j}^{n-1} C(k, j) * c[k] * shift^{k - j}
        for (j, slot) in next.iter_mut().enumerate() {
            let mut acc = ZERO;
            let mut binom = 1.0_f64;
            let mut shift_pow = ONE;
            for k in j..n {
                acc += binom * *coeffs.get(k)? * shift_pow;
                shift_pow *= shift;
                binom *= ((k + 1) as f64) / ((k - j + 1) as f64);
            }
            *slot = acc;
        }
        coeffs = next;
    }
    Some(coeffs)
}

/// Linear homogeneous ODE `sum_{j=0}^r p_j(z) y^(j)(z) = 0` with polynomial
/// coefficients `coeffs[j]` in ascending powers of `z`.
#[derive(Clone, Debug)]
pub struct HolonomicOde {
    /// Polynomial coefficients `[p_0(z), p_1(z), ..., p_r(z)]`.
    pub coeffs: Vec<Vec<Complex64>>,
}

impl HolonomicOde {
    /// Order `r` of the differential equation.
    #[must_use]
    pub fn order(&self) -> usize {
        self.coeffs.len().saturating_sub(1)
    }

    /// Evaluates polynomial `p_j(z)` shifted to `w = z - z0` (Taylor coefficients around `z0`).
    fn shifted_poly(
        poly: &[Complex64],
        z0: Complex64,
    ) -> Vec<Complex64> {
        let deg = poly.len();
        let mut out = vec![ZERO; deg];
        for (k, &ck) in poly.iter().enumerate() {
            // ck * (z0 + w)^k
            let mut binom = 1.0_f64;
            for j in 0..=k {
                let term = ck * binom * z0.powi((k - j) as i32);
                if let Some(slot) = out.get_mut(j) {
                    *slot += term;
                }
                if j < k {
                    binom *= ((k - j) as f64) / ((j + 1) as f64);
                }
            }
        }
        out
    }

    /// Estimates the distance from `z0` to the nearest root of the leading coefficient `p_r(z)`.
    fn nearest_singularity_dist(
        &self,
        z0: Complex64,
    ) -> f64 {
        let Some(lead) = self.coeffs.last() else {
            return f64::INFINITY;
        };
        let mut deg = lead.len().saturating_sub(1);
        while deg > 0 && lead.get(deg).is_some_and(|c| c.norm() == 0.0) {
            deg -= 1;
        }
        if deg == 0 {
            return f64::INFINITY;
        }
        // Shift p_r(z0 + w) = sum b_k w^k; Cauchy root lower bound: |b_0| / (|b_0| + max |b_k|)
        let shifted = Self::shifted_poly(lead, z0);
        let b0 = shifted.first().map_or(0.0, |c| c.norm());
        if b0 == 0.0 {
            return 0.0;
        }
        let max_b = shifted
            .iter()
            .skip(1)
            .map(|c| c.norm())
            .fold(0.0_f64, f64::max);
        if max_b == 0.0 {
            return f64::INFINITY;
        }
        b0 / (b0 + max_b)
    }

    /// Expands `y(z0 + w) = sum_{n=0}^order a_n w^n` from `initial_state = [y(z0), y'(z0), ...]`.
    fn local_taylor(
        &self,
        z0: Complex64,
        initial_state: &[Complex64],
        taylor_order: usize,
    ) -> Option<Vec<Complex64>> {
        let r = self.order();
        if r == 0 || initial_state.len() != r {
            return None;
        }
        let shifted: Vec<Vec<Complex64>> = self
            .coeffs
            .iter()
            .map(|p| Self::shifted_poly(p, z0))
            .collect();
        let lead0 = shifted.last()?.first().copied()?;
        if lead0.norm() < 1e-30 {
            return None;
        }

        let mut a = vec![ZERO; taylor_order + 1];
        let mut fact = 1.0_f64;
        for k in 0..r {
            if k > 0 {
                fact *= k as f64;
            }
            if let Some(slot) = a.get_mut(k) {
                *slot = *initial_state.get(k)? / fact;
            }
        }

        // Solve for a_{n + r} from coefficient of w^n in sum_{j=0}^r p_j(z0 + w) y^(j)(z0 + w) = 0
        for n in 0..=taylor_order.saturating_sub(r) {
            let mut rhs = ZERO;
            for (j, pj) in shifted.iter().enumerate() {
                for (m, &pjm) in pj.iter().enumerate() {
                    if m > n {
                        break;
                    }
                    let k = n - m; // index in y^(j), corresponding to a_{k + j}
                    if j == r && m == 0 {
                        continue;
                    }
                    let idx = k + j;
                    if let Some(&akj) = a.get(idx) {
                        // falling factorial (k+j)! / k!
                        let mut fall = 1.0_f64;
                        for f in (k + 1)..=(k + j) {
                            fall *= f as f64;
                        }
                        rhs += pjm * fall * akj;
                    }
                }
            }
            let mut fall_r = 1.0_f64;
            for f in (n + 1)..=(n + r) {
                fall_r *= f as f64;
            }
            if let Some(slot) = a.get_mut(n + r) {
                *slot = -rhs / (lead0 * fall_r);
            }
        }
        Some(a)
    }
}

/// Analytically continues `[y, y', ..., y^(r-1)]` of a [`HolonomicOde`] along
/// the polygonal `path` using adaptive local Taylor series steps.
#[must_use]
pub fn continue_holonomic_along_path(
    ode: &HolonomicOde,
    path: &[Complex64],
    initial_state: &[Complex64],
    taylor_order: usize,
) -> Option<Vec<Complex64>> {
    let r = ode.order();
    if path.is_empty() || initial_state.len() != r {
        return None;
    }
    let mut state = initial_state.to_vec();
    let order = taylor_order.max(r + 8);

    for win in path.windows(2) {
        let mut cur = win[0];
        let target = win[1];
        let mut steps = 0_usize;
        while (target - cur).norm() > 1e-14 {
            steps += 1;
            if steps > 10_000 {
                return None;
            }
            let dist = ode.nearest_singularity_dist(cur);
            if dist <= 1e-12 {
                return None;
            }
            let rem = target - cur;
            let max_step = (0.45 * dist).min(rem.norm());
            let step = if rem.norm() <= max_step {
                rem
            } else {
                rem * (max_step / rem.norm())
            };
            let coeffs = ode.local_taylor(cur, &state, order)?;
            // Evaluate y^(j)(cur + step) for j = 0..r
            let mut next_state = vec![ZERO; r];
            for (j, slot) in next_state.iter_mut().enumerate() {
                let mut val = ZERO;
                let mut w_pow = ONE;
                for (n, &an) in coeffs.iter().enumerate().skip(j) {
                    let mut fall = 1.0_f64;
                    for f in (n - j + 1)..=n {
                        fall *= f as f64;
                    }
                    val += fall * an * w_pow;
                    w_pow *= step;
                }
                *slot = val;
            }
            state = next_state;
            cur += step;
        }
    }
    Some(state)
}

/// Computes the `r x r` monodromy matrix of a [`HolonomicOde`] around the
/// closed loop `loop_pts` (closing back to `loop_pts[0]`).
#[must_use]
pub fn holonomic_monodromy_matrix(
    ode: &HolonomicOde,
    loop_pts: &[Complex64],
    taylor_order: usize,
) -> Option<Vec<Vec<Complex64>>> {
    let r = ode.order();
    let &first = loop_pts.first()?;
    let mut closed = loop_pts.to_vec();
    if (closed.last().copied()? - first).norm() > 1e-12 {
        closed.push(first);
    }
    let mut cols = Vec::with_capacity(r);
    for j in 0..r {
        let mut e = vec![ZERO; r];
        if let Some(slot) = e.get_mut(j) {
            *slot = ONE;
        }
        cols.push(continue_holonomic_along_path(
            ode,
            &closed,
            &e,
            taylor_order,
        )?);
    }
    let mut mat = vec![vec![ZERO; r]; r];
    for (j, col) in cols.iter().enumerate() {
        for (i, &val) in col.iter().enumerate() {
            if let Some(cell) = mat.get_mut(i).and_then(|row| row.get_mut(j)) {
                *cell = val;
            }
        }
    }
    Some(mat)
}
