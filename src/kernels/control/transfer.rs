//! Transfer functions: roots, stability, state-space conversion, block
//! algebra and frequency response.
//!
//! Coefficients are stored in descending powers of `s`, the same convention as
//! [`crate::kernels::polynomial`]. Poles and zeros are the eigenvalues of the
//! companion matrix, so complex-conjugate roots are handled without a separate
//! complex root finder.

use crate::kernels::complex_scalar::Complex;
use num_complex::Complex64;
use std::sync::LazyLock;

use super::linalg;
use super::state_space::StateSpace;

/// A scalar rational transfer function `num(s) / den(s)`.
#[derive(Debug, Clone, PartialEq)]
pub struct TransferFunction {
    /// Numerator coefficients, descending powers of `s`.
    pub num: Vec<f64>,
    /// Denominator coefficients, descending powers of `s`.
    pub den: Vec<f64>,
}

/// Removes leading zero coefficients.
fn trim_leading(coeffs: &[f64]) -> &[f64] {
    let start = coeffs
        .iter()
        .position(|&c| c != 0.0)
        .unwrap_or(coeffs.len());
    &coeffs[start..]
}

fn trim_owned(mut coeffs: Vec<f64>) -> Vec<f64> {
    let leading = coeffs.len() - trim_leading(&coeffs).len();
    coeffs.drain(..leading);
    coeffs
}

/// Polynomial product in descending coefficients.
fn conv(
    a: &[f64],
    b: &[f64],
) -> Vec<f64> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut out = vec![0.0; a.len() + b.len() - 1];
    for (i, &ai) in a.iter().enumerate() {
        for (j, &bj) in b.iter().enumerate() {
            out[i + j] += ai * bj;
        }
    }
    out
}

/// Sum of two polynomials given in descending coefficients.
fn poly_add(
    a: &[f64],
    b: &[f64],
) -> Vec<f64> {
    let len = a.len().max(b.len());
    let mut out = vec![0.0; len];
    let off_a = len - a.len();
    let off_b = len - b.len();
    for i in 0..a.len() {
        out[off_a + i] += a[i];
    }
    for i in 0..b.len() {
        out[off_b + i] += b[i];
    }
    out
}

/// Roots of a polynomial given in descending coefficients.
fn poly_roots(coeffs: &[f64]) -> Result<Vec<Complex64>, String> {
    let c = trim_leading(coeffs);
    let [lead, rest @ ..] = c else {
        return Ok(Vec::new());
    };
    match rest {
        | [] => return Ok(Vec::new()),
        | [constant] => {
            let root = -constant / lead;
            if root.is_finite() {
                return Ok(vec![Complex64::new(root, 0.0)]);
            }
            return Err("polynomial root falls outside Float range".to_owned());
        },
        | _ => {},
    }
    let n = rest.len();
    let first_row = rest.iter().map(|&c| -c / lead).collect::<Vec<_>>();
    if first_row.iter().any(|c| !c.is_finite()) {
        return Err("normalized polynomial coefficients fall outside Float range".to_owned());
    }
    let mut companion = Vec::with_capacity(n);
    companion.push(first_row);
    for i in 1..n {
        let mut row = vec![0.0; n];
        row[i - 1] = 1.0;
        companion.push(row);
    }
    let roots = linalg::eigenvalues(&companion)?;
    if roots.len() != n
        || roots
            .iter()
            .any(|root| !root.re.is_finite() || !root.im.is_finite())
    {
        return Err("polynomial eigensolver did not return all finite roots".to_owned());
    }
    Ok(roots)
}

/// Evaluates a polynomial with complex argument by Horner's rule.
fn poly_eval(
    coeffs: &[f64],
    s: Complex64,
) -> Complex64 {
    coeffs.iter().fold(Complex64::new(0.0, 0.0), |acc, &c| {
        acc * s + Complex64::new(c, 0.0)
    })
}

impl TransferFunction {
    /// Builds a transfer function after trimming leading zero coefficients.
    ///
    /// Improper rational functions (`deg(num) > deg(den)`, such as a pure
    /// differentiator or PD controller) are permitted here so they can be used
    /// in frequency-response and block-algebra calculations; use
    /// [`TransferFunction::new_proper`] when a state-space realisation is
    /// required.
    ///
    /// # Errors
    /// Returns an error when the denominator is zero or when a coefficient is
    /// not finite.
    pub fn new(
        num: Vec<f64>,
        den: Vec<f64>,
    ) -> Result<Self, String> {
        if num.iter().chain(&den).any(|v| !v.is_finite()) {
            return Err("coefficients must be finite".to_string());
        }
        let den = trim_owned(den);
        if den.is_empty() {
            return Err("denominator must be nonzero".to_string());
        }
        let num = trim_owned(num);
        let num = if num.is_empty() {
            vec![0.0]
        } else {
            num
        };
        Ok(Self { num, den })
    }

    /// Builds a proper transfer function (`deg(num) <= deg(den)`) after
    /// trimming leading zero coefficients.
    ///
    /// # Errors
    /// Returns an error when [`TransferFunction::new`] fails or when the
    /// numerator has higher degree than the denominator.
    pub fn new_proper(
        num: Vec<f64>,
        den: Vec<f64>,
    ) -> Result<Self, String> {
        let tf = Self::new(num, den)?;
        if !tf.is_proper() {
            return Err("numerator must not have higher degree than denominator".to_string());
        }
        Ok(tf)
    }

    /// Whether the transfer function is proper (`deg(num) <= deg(den)`).
    ///
    /// A zero numerator is always proper.
    pub fn is_proper(&self) -> bool {
        let num = trim_leading(&self.num);
        let den = trim_leading(&self.den);
        !den.is_empty() && num.len() <= den.len()
    }

    /// Whether the transfer function is strictly proper (`deg(num) < deg(den)`).
    ///
    /// A zero numerator is always strictly proper.
    pub fn is_strictly_proper(&self) -> bool {
        let num = trim_leading(&self.num);
        let den = trim_leading(&self.den);
        !den.is_empty() && num.len() < den.len()
    }

    /// Order of the denominator (degree of `den`).
    pub fn order(&self) -> usize {
        self.den.len().saturating_sub(1)
    }

    /// Poles: roots of the denominator.
    pub fn poles(&self) -> Vec<Complex64> {
        self.poles_checked().unwrap_or_default()
    }

    pub(crate) fn poles_checked(&self) -> Result<Vec<Complex64>, String> {
        poly_roots(&self.den)
    }

    /// Zeros: roots of the numerator (empty for a zero numerator).
    pub fn zeros(&self) -> Vec<Complex64> {
        self.zeros_checked().unwrap_or_default()
    }

    pub(crate) fn zeros_checked(&self) -> Result<Vec<Complex64>, String> {
        poly_roots(&self.num)
    }

    /// Whether all poles can be computed and have negative real part.
    pub fn is_stable(&self) -> bool {
        self.poles_checked().is_ok_and(|poles| stable_poles(&poles))
    }

    /// Evaluates `G(s)` at a complex frequency.
    pub fn evaluate(
        &self,
        s: Complex64,
    ) -> Complex64 {
        let numerator = poly_eval(&self.num, s);
        let denominator = poly_eval(&self.den, s);
        let ratio = Complex::new(numerator.re, numerator.im)
            .div(Complex::new(denominator.re, denominator.im));
        Complex64::new(ratio.re, ratio.im)
    }

    pub(crate) fn evaluate_checked(
        &self,
        s: Complex64,
    ) -> Result<Complex64, String> {
        let response = self.evaluate(s);
        if response.re.is_finite() && response.im.is_finite() {
            Ok(response)
        } else {
            Err("transfer response is singular or outside Float range".to_owned())
        }
    }

    /// Converts to a continuous-time state-space model in controllable
    /// canonical form.
    ///
    /// # Errors
    /// Returns an error when the numerator is improper or the denominator has
    /// no nonzero leading coefficient.
    pub fn to_state_space(&self) -> Result<StateSpace, String> {
        if self.num.len() > self.den.len() {
            return Err("numerator must not have higher degree than denominator".to_string());
        }
        let Some(lead) = self.den.first().copied().filter(|&c| c != 0.0) else {
            return Err("denominator must be nonzero".to_string());
        };
        let n = self.den.len() - 1;
        let den: Vec<f64> = self.den.iter().map(|c| c / lead).collect();
        let mut num = vec![0.0; n + 1];
        let offset = n + 1 - self.num.len();
        for (i, &c) in self.num.iter().enumerate() {
            num[offset + i] = c / lead;
        }
        let direct = num[0];
        let mut a = vec![vec![0.0; n]; n];
        for j in 0..n {
            a[0][j] = -den[j + 1];
        }
        for i in 1..n {
            a[i][i - 1] = 1.0;
        }
        let mut b = vec![vec![0.0; 1]; n];
        if n > 0 {
            b[0][0] = 1.0;
        }
        let mut c = vec![vec![0.0; n]; 1];
        for j in 0..n {
            c[0][j] = num[j + 1] - direct * den[j + 1];
        }
        StateSpace::continuous(a, b, c, vec![vec![direct]])
    }
}

pub(crate) fn stable_poles(poles: &[Complex64]) -> bool {
    poles.iter().all(|pole| pole.re < 0.0)
}

/// Converts a state-space model to one transfer function per input/output
/// pair, ordered output-major (output `0` against every input, then output
/// `1`, and so on).
///
/// The denominator is obtained by the Faddeev-LeVerrier algorithm and the
/// numerator from the Markov parameters `C A^k B` plus the feedthrough term.
///
/// # Errors
/// Returns an error when the model has inconsistent dimensions or the
/// polynomial construction fails.
pub fn ss_to_tf(ss: &StateSpace) -> Result<Vec<TransferFunction>, String> {
    let p = ss.n_outputs();
    let m = ss.n_inputs();
    let n = ss.n_states();
    if n == 0 {
        let mut out = Vec::with_capacity(p * m);
        for i in 0..p {
            for j in 0..m {
                let d = ss.d.get(i).and_then(|r| r.get(j)).copied().unwrap_or(0.0);
                out.push(TransferFunction::new(vec![d], vec![1.0])?);
            }
        }
        return Ok(out);
    }
    let mut m_prev = linalg::identity(n);
    let mut den = vec![0.0; n + 1];
    den[0] = 1.0;
    let mut markov = Vec::with_capacity(n);
    markov.push(m_prev.clone());
    for k in 1..=n {
        let am = linalg::matmul(&ss.a, &m_prev)?;
        let mut trace = 0.0;
        for i in 0..n {
            trace += am[i][i];
        }
        let ck = -trace / k as f64;
        den[k] = ck;
        let mut mk = am;
        for i in 0..n {
            mk[i][i] += ck;
        }
        if k < n {
            markov.push(mk.clone());
        }
        m_prev = mk;
    }
    let mut out = Vec::with_capacity(p * m);
    for i in 0..p {
        for j in 0..m {
            let d = ss.d.get(i).and_then(|r| r.get(j)).copied().unwrap_or(0.0);
            let mut num = vec![0.0; n + 1];
            num[0] = d;
            for l in 1..=n {
                num[l] += d * den[l];
            }
            let b_col: Vec<f64> = (0..n)
                .map(|r| {
                    ss.b.get(r)
                        .and_then(|row| row.get(j))
                        .copied()
                        .unwrap_or(0.0)
                })
                .collect();
            for (k, mk) in markov.iter().enumerate() {
                let mb = linalg::matvec(mk, &b_col);
                let mut value = 0.0;
                for (c_idx, &bv) in mb.iter().enumerate() {
                    let ci =
                        ss.c.get(i)
                            .and_then(|row| row.get(c_idx))
                            .copied()
                            .unwrap_or(0.0);
                    value += ci * bv;
                }
                num[1 + k] += value;
            }
            out.push(TransferFunction::new(num, den.clone())?);
        }
    }
    Ok(out)
}

/// Series connection `a b`.
///
/// # Errors
/// As [`TransferFunction::new`].
pub fn series(
    a: &TransferFunction,
    b: &TransferFunction,
) -> Result<TransferFunction, String> {
    TransferFunction::new(conv(&a.num, &b.num), conv(&a.den, &b.den))
}

/// Parallel connection `a + b`.
///
/// # Errors
/// As [`TransferFunction::new`].
pub fn parallel(
    a: &TransferFunction,
    b: &TransferFunction,
) -> Result<TransferFunction, String> {
    let num = poly_add(&conv(&a.num, &b.den), &conv(&b.num, &a.den));
    TransferFunction::new(num, conv(&a.den, &b.den))
}

/// Negative feedback `a / (1 + a b)`.
///
/// # Errors
/// As [`TransferFunction::new`].
pub fn feedback(
    a: &TransferFunction,
    b: &TransferFunction,
) -> Result<TransferFunction, String> {
    let num = conv(&a.num, &b.den);
    let den = poly_add(&conv(&a.den, &b.den), &conv(&a.num, &b.num));
    TransferFunction::new(num, den)
}

/// Bode data: magnitude in decibels and unwrapped phase in degrees.
///
/// # Errors
/// Returns an error when a frequency is negative or not finite.
pub fn bode(
    tf: &TransferFunction,
    omega: &[f64],
) -> Result<(Vec<f64>, Vec<f64>), String> {
    let mut magnitude = Vec::with_capacity(omega.len());
    let mut phase = Vec::with_capacity(omega.len());
    for &w in omega {
        if !(w.is_finite() && w >= 0.0) {
            return Err("frequencies must be finite and non-negative".to_string());
        }
        let g = tf.evaluate(Complex64::new(0.0, w));
        magnitude.push(20.0 * g.norm().log10());
        phase.push(g.arg().to_degrees());
    }
    for i in 1..phase.len() {
        phase[i] = unwrap_phase(phase[i], phase[i - 1]);
    }
    Ok((magnitude, phase))
}

/// Nyquist plot points `(Re G(j w), Im G(j w))`.
///
/// # Errors
/// Returns an error when a frequency is negative or not finite.
pub fn nyquist(
    tf: &TransferFunction,
    omega: &[f64],
) -> Result<Vec<(f64, f64)>, String> {
    let mut out = Vec::with_capacity(omega.len());
    for &w in omega {
        if !(w.is_finite() && w >= 0.0) {
            return Err("frequencies must be finite and non-negative".to_string());
        }
        let g = tf.evaluate(Complex64::new(0.0, w));
        out.push((g.re, g.im));
    }
    Ok(out)
}

/// Unwraps a phase in degrees to lie within 180 degrees of `reference`.
fn unwrap_phase(
    mut phase: f64,
    reference: f64,
) -> f64 {
    while phase - reference > 180.0 {
        phase -= 360.0;
    }
    while phase - reference < -180.0 {
        phase += 360.0;
    }
    phase
}

fn phase_near(
    tf: &TransferFunction,
    w: f64,
    reference: f64,
) -> Result<f64, String> {
    let phase = tf
        .evaluate_checked(Complex64::new(0.0, w))?
        .arg()
        .to_degrees();
    Ok(unwrap_phase(phase, reference))
}

fn crossover_magnitude(response: Complex64) -> Result<f64, String> {
    let magnitude = response.norm();
    if magnitude.is_finite() {
        Ok(magnitude)
    } else {
        Err("loop magnitude is outside Float range".to_owned())
    }
}

/// Gain and phase margins of a loop transfer function.
///
/// Returns `(gain margin in dB, phase margin in degrees)`. Each margin is
/// `None` when the corresponding crossover does not occur on the internal
/// log-spaced grid from `1e-3` to `1e3` rad/s. Crossovers are refined by
/// bisection on `log10(omega)`.
///
/// # Errors
/// Returns an error when the transfer function cannot be evaluated.
pub fn gain_margin_phase_margin(
    tf: &TransferFunction
) -> Result<(Option<f64>, Option<f64>), String> {
    const POINTS: usize = 601;
    static GRID: LazyLock<[f64; POINTS]> = LazyLock::new(|| {
        std::array::from_fn(|i| 10.0_f64.powf(-3.0 + 6.0 * i as f64 / (POINTS - 1) as f64))
    });
    let grid = &*GRID;
    let mut phase = Vec::with_capacity(POINTS);
    let mut magnitude = Vec::with_capacity(POINTS);
    for &w in grid {
        let g = tf.evaluate_checked(Complex64::new(0.0, w))?;
        phase.push(g.arg().to_degrees());
        magnitude.push(crossover_magnitude(g)?);
    }
    for i in 1..phase.len() {
        phase[i] = unwrap_phase(phase[i], phase[i - 1]);
    }
    let mut gain_margin = None;
    for i in 0..POINTS - 1 {
        if phase[i] + 180.0 > 0.0 && phase[i + 1] + 180.0 <= 0.0 {
            let mut left = grid[i].log10();
            let mut right = grid[i + 1].log10();
            let mut p_left = phase[i];
            for _ in 0..100 {
                let mid = f64::midpoint(left, right);
                let p_mid = phase_near(tf, 10.0_f64.powf(mid), p_left)?;
                if p_mid + 180.0 > 0.0 {
                    left = mid;
                    p_left = p_mid;
                } else {
                    right = mid;
                }
            }
            let w = 10.0_f64.powf(f64::midpoint(left, right));
            let mag = crossover_magnitude(tf.evaluate_checked(Complex64::new(0.0, w))?)?;
            if mag > 0.0 {
                gain_margin = Some(-20.0 * mag.log10());
            }
            break;
        }
    }
    let mut phase_margin = None;
    for i in 0..POINTS - 1 {
        if magnitude[i] > 1.0 && magnitude[i + 1] <= 1.0 {
            let mut left = grid[i].log10();
            let mut right = grid[i + 1].log10();
            for _ in 0..100 {
                let mid = f64::midpoint(left, right);
                let mag = crossover_magnitude(
                    tf.evaluate_checked(Complex64::new(0.0, 10.0_f64.powf(mid)))?,
                )?;
                if mag > 1.0 {
                    left = mid;
                } else {
                    right = mid;
                }
            }
            let w = 10.0_f64.powf(f64::midpoint(left, right));
            let p = phase_near(tf, w, phase[i])?;
            phase_margin = Some(180.0 + p);
            break;
        }
    }
    Ok((gain_margin, phase_margin))
}

/// Closed-loop step response of `tf` on the time grid `t`.
///
/// The model is converted to state space and integrated with fourth-order
/// Runge-Kutta on each interval `t[k+1] - t[k]`.
///
/// # Errors
/// Returns an error when the time grid or its increments are not finite, the
/// grid is not strictly increasing, the conversion fails, or integration
/// exceeds the finite range.
pub fn step_response(
    tf: &TransferFunction,
    t: &[f64],
) -> Result<Vec<f64>, String> {
    if t.iter().any(|time| !time.is_finite()) {
        return Err("time vector must be finite".to_owned());
    }
    let ss = tf.to_state_space()?;
    let n = ss.n_states();
    let mut x = vec![0.0; n];
    let mut out = Vec::with_capacity(t.len());
    if t.is_empty() {
        return Ok(out);
    }
    out.push(output_scalar(&ss, &x, 1.0));
    for pair in t.windows(2) {
        let h = pair[1] - pair[0];
        if !(h.is_finite() && h > 0.0) {
            return Err("time vector must have finite, strictly positive increments".to_owned());
        }
        x = linalg::rk4_linear(&ss.a, &ss.b, &x, &[1.0], h);
        let sample = output_scalar(&ss, &x, 1.0);
        if !sample.is_finite() {
            return Err("step response is outside Float range".to_owned());
        }
        out.push(sample);
    }
    Ok(out)
}

/// Scalar output `C x + D u` for a single-input, single-output model.
fn output_scalar(
    ss: &StateSpace,
    x: &[f64],
    u: f64,
) -> f64 {
    let cx = linalg::matvec(&ss.c, x);
    let du = ss.d.first().and_then(|r| r.first()).copied().unwrap_or(0.0);
    cx.first().copied().unwrap_or(0.0) + du * u
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tf(
        num: Vec<f64>,
        den: Vec<f64>,
    ) -> TransferFunction {
        TransferFunction::new(num, den).unwrap()
    }

    #[test]
    fn leading_zeros_are_trimmed_and_zero_denominator_rejected() {
        let g = tf(vec![0.0, 1.0, 2.0], vec![0.0, 1.0, 3.0, 2.0]);
        assert_eq!(g.num, vec![1.0, 2.0]);
        assert_eq!(g.den, vec![1.0, 3.0, 2.0]);
        assert!(TransferFunction::new(vec![1.0], vec![0.0, 0.0]).is_err());
    }

    #[test]
    fn poles_and_zeros_of_second_order() {
        let g = tf(vec![1.0, 2.0], vec![1.0, 3.0, 2.0]);
        let mut poles: Vec<f64> = g.poles().iter().map(|p| p.re).collect();
        poles.sort_by(f64::total_cmp);
        assert!((poles[0] + 2.0).abs() < 1e-9);
        assert!((poles[1] + 1.0).abs() < 1e-9);
        assert!((g.zeros()[0].re + 2.0).abs() < 1e-9);
        assert!(g.is_stable());
    }

    #[test]
    fn round_trip_through_state_space() {
        let g = tf(vec![1.0, 2.0], vec![1.0, 3.0, 2.0]);
        let ss = g.to_state_space().unwrap();
        let back = ss_to_tf(&ss).unwrap();
        assert_eq!(back.len(), 1);
        for s in [
            Complex64::new(0.5, 0.0),
            Complex64::new(0.0, 2.0),
            Complex64::new(3.0, -1.0),
        ] {
            let a = g.evaluate(s);
            let b = back[0].evaluate(s);
            assert!((a - b).norm() < 1e-9, "s = {s:?}");
        }
    }

    #[test]
    fn bode_of_first_order_at_corner() {
        let g = tf(vec![1.0], vec![1.0, 1.0]);
        let (mag, phase) = bode(&g, &[1.0]).unwrap();
        assert!((mag[0] + 3.010_299_956_639_812).abs() < 1e-9);
        assert!((phase[0] + 45.0).abs() < 1e-9);
    }

    #[test]
    fn gain_margin_of_third_order_loop() {
        let g = tf(vec![1.0], vec![1.0, 3.0, 2.0, 0.0]);
        let (gm, pm) = gain_margin_phase_margin(&g).unwrap();
        let gm = gm.unwrap();
        // The gain margin as a ratio is 6, i.e. about 15.56 dB.
        assert!((10.0_f64.powf(gm / 20.0) - 6.0).abs() < 0.05);
        let pm = pm.unwrap();
        assert!(pm > 40.0 && pm < 70.0);
    }

    #[test]
    fn block_algebra_matches_pointwise_products() {
        let a = tf(vec![1.0], vec![1.0, 1.0]);
        let b = tf(vec![2.0], vec![1.0, 2.0]);
        let s = series(&a, &b).unwrap();
        let p = parallel(&a, &b).unwrap();
        let f = feedback(&a, &b).unwrap();
        for w in [0.0, 0.5, 2.0] {
            let z = Complex64::new(0.0, w);
            let ga = a.evaluate(z);
            let gb = b.evaluate(z);
            assert!((s.evaluate(z) - ga * gb).norm() < 1e-9);
            assert!((p.evaluate(z) - (ga + gb)).norm() < 1e-9);
            assert!((f.evaluate(z) - ga / (1.0 + ga * gb)).norm() < 1e-9);
        }
    }

    #[test]
    fn step_response_of_first_order() {
        let g = tf(vec![1.0], vec![1.0, 1.0]);
        let t: Vec<f64> = (0..101).map(|i| f64::from(i) * 0.01).collect();
        let y = step_response(&g, &t).unwrap();
        assert!((y[0]).abs() < 1e-12);
        assert!((y[100] - (1.0 - (-1.0_f64).exp())).abs() < 1e-6);
    }

    #[test]
    fn properness_helpers_and_new_proper() {
        let strictly_proper = tf(vec![1.0], vec![1.0, 1.0]);
        assert!(strictly_proper.is_proper());
        assert!(strictly_proper.is_strictly_proper());
        assert!(TransferFunction::new_proper(vec![1.0], vec![1.0, 1.0]).is_ok());

        let biproper = TransferFunction::new_proper(vec![2.0, 1.0], vec![1.0, 1.0]).unwrap();
        assert!(biproper.is_proper());
        assert!(!biproper.is_strictly_proper());
        let ss = biproper.to_state_space().unwrap();
        assert!((ss.d[0][0] - 2.0).abs() < 1e-12);
        let back = ss_to_tf(&ss).unwrap();
        let s = Complex64::new(0.0, 3.0);
        assert!((biproper.evaluate(s) - back[0].evaluate(s)).norm() < 1e-9);

        let zero_tf = TransferFunction::new_proper(vec![0.0], vec![1.0]).unwrap();
        assert!(zero_tf.is_proper());
        assert!(zero_tf.is_strictly_proper());

        // Improper transfer function is allowed by `new` (for Bode/Nyquist) but
        // rejected by `new_proper` and `to_state_space`.
        let improper = TransferFunction::new(vec![1.0, 0.0, 1.0], vec![1.0, 1.0]).unwrap();
        assert!(!improper.is_proper());
        assert!(!improper.is_strictly_proper());
        assert!(TransferFunction::new_proper(vec![1.0, 0.0, 1.0], vec![1.0, 1.0]).is_err());
        assert!(improper.to_state_space().is_err());
    }

    #[test]
    fn nyquist_and_step_response_edge_cases() {
        let g = tf(vec![1.0], vec![1.0, 1.0]);
        let omega = [0.0, 0.5, 1.0, 2.0];
        let pts = nyquist(&g, &omega).unwrap();
        for (&w, &(re, im)) in omega.iter().zip(&pts) {
            let denom = 1.0 + w * w;
            assert!((re - 1.0 / denom).abs() < 1e-12);
            assert!((im + w / denom).abs() < 1e-12);
        }
        assert!(nyquist(&g, &[-1.0]).is_err());
        assert!(bode(&g, &[f64::NAN]).is_err());

        // Empty and non-increasing time grids in step_response.
        assert_eq!(step_response(&g, &[]).unwrap(), [] as [f64; 0]);
        assert!(step_response(&g, &[0.0, 1.0, 1.0]).is_err());

        // Biproper (2s + 1) / (s + 1) = 2 - 1/(s + 1) has step response 1 + exp(-t).
        let biproper = tf(vec![2.0, 1.0], vec![1.0, 1.0]);
        let t: Vec<f64> = (0..101).map(|i| f64::from(i) * 0.01).collect();
        let y = step_response(&biproper, &t).unwrap();
        assert!((y[0] - 2.0).abs() < 1e-12);
        assert!((y[100] - (1.0 + (-1.0_f64).exp())).abs() < 1e-6);
    }

    #[test]
    fn ss_to_tf_handles_static_and_mimo_models() {
        let static_ss = StateSpace::continuous(
            Vec::new(),
            Vec::new(),
            vec![Vec::new()],
            vec![vec![3.0, -1.0]],
        )
        .unwrap();
        let tfs = ss_to_tf(&static_ss).unwrap();
        assert_eq!(tfs.len(), 2);
        assert!((tfs[0].evaluate(Complex64::new(0.0, 1.0)).re - 3.0).abs() < 1e-12);
        assert!((tfs[1].evaluate(Complex64::new(0.0, 1.0)).re + 1.0).abs() < 1e-12);

        // 1-state, 2-input, 2-output system: x' = -x + [1, 2] u, y = [1; 3] x.
        let mimo = StateSpace::continuous(
            vec![vec![-1.0]],
            vec![vec![1.0, 2.0]],
            vec![vec![1.0], vec![3.0]],
            vec![vec![0.0, 0.0], vec![0.0, 0.0]],
        )
        .unwrap();
        let mimo_tfs = ss_to_tf(&mimo).unwrap();
        assert_eq!(mimo_tfs.len(), 4);
        let s = Complex64::new(1.0, 0.0);
        // Output-major gains at s = 1: 1/2, 2/2, 3/2, 6/2.
        let expected = [0.5, 1.0, 1.5, 3.0];
        for (got, want) in mimo_tfs.iter().zip(expected) {
            assert!((got.evaluate(s).re - want).abs() < 1e-12);
        }
    }
}
