//! Complex root localization, bounds, polynomial solvers, argument principle,
//! Delves-Lyness, Rouché theorem checks, and 1D complex homotopy continuation.

use num_complex::Complex;
use std::f64::consts::PI;

fn eval_poly_and_deriv(
    coeffs: &[Complex<f64>],
    z: Complex<f64>,
) -> (Complex<f64>, Complex<f64>, Complex<f64>) {
    let mut p = Complex::new(0.0, 0.0);
    let mut dp = Complex::new(0.0, 0.0);
    let mut ddp = Complex::new(0.0, 0.0);
    for &c in coeffs.iter().rev() {
        ddp = ddp * z + dp * 2.0;
        dp = dp * z + p;
        p = p * z + c;
    }
    (p, dp, ddp)
}

fn trim_leading_zeros(coeffs: &[Complex<f64>]) -> &[Complex<f64>] {
    let mut end = coeffs.len();
    while end > 0 {
        if let Some(c) = coeffs.get(end.saturating_sub(1)) {
            if c.norm() > 0.0 {
                break;
            }
        }
        end = end.saturating_sub(1);
    }
    coeffs.get(..end).unwrap_or(&[])
}

/// Classic Cauchy upper bound on the modulus of all roots of
/// $p(z) = \sum_{k=0}^n c_k z^k$ (ascending coefficients):
/// $R = 1 + \max_{0 \le k < n} |c_k / c_n|$.
pub fn cauchy_root_bound(coeffs: &[Complex<f64>]) -> Option<f64> {
    let c = trim_leading_zeros(coeffs);
    if c.len() <= 1 {
        return None;
    }
    let cn = c.last()?.norm();
    if cn == 0.0 {
        return None;
    }
    let mut max_ratio = 0.0_f64;
    for ck in c.iter().take(c.len().saturating_sub(1)) {
        max_ratio = max_ratio.max(ck.norm() / cn);
    }
    Some(1.0 + max_ratio)
}

/// Sharp Cauchy root bound solving $|c_n| R^n = \sum_{k=0}^{n-1} |c_k| R^k$ via bisection.
pub fn cauchy_sharp_root_bound(coeffs: &[Complex<f64>]) -> Option<f64> {
    let c = trim_leading_zeros(coeffs);
    if c.len() <= 1 {
        return None;
    }
    let cn = c.last()?.norm();
    if cn == 0.0 {
        return None;
    }
    let mut hi = cauchy_root_bound(c)?;
    let mut lo = 0.0_f64;
    let eval_g = |r: f64| -> f64 {
        let mut sum = 0.0_f64;
        let mut r_pow = 1.0_f64;
        for ck in c.iter().take(c.len().saturating_sub(1)) {
            sum += ck.norm() * r_pow;
            r_pow *= r;
        }
        cn * r_pow - sum
    };
    if eval_g(0.0) >= 0.0 {
        return Some(0.0);
    }
    while eval_g(hi) < 0.0 {
        hi *= 2.0;
    }
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if eval_g(mid) >= 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    Some(hi)
}

/// Fujiwara upper bound on the modulus of all polynomial roots:
/// $R = 2 \max_{1 \le k \le n} \left|\frac{c_{n-k}}{c_n}\right|^{1/k}$ with factor $1/2$ on $k=n$.
pub fn fujiwara_root_bound(coeffs: &[Complex<f64>]) -> Option<f64> {
    let c = trim_leading_zeros(coeffs);
    let n = c.len().checked_sub(1)?;
    if n == 0 {
        return None;
    }
    let cn = c.last()?.norm();
    if cn == 0.0 {
        return None;
    }
    let mut m = 0.0_f64;
    for k in 1..=n {
        let idx = n.saturating_sub(k);
        let ck = c.get(idx)?.norm();
        let mut ratio = ck / cn;
        if k == n {
            ratio *= 0.5;
        }
        let term = ratio.powf(1.0 / (k as f64));
        m = m.max(term);
    }
    Some(2.0 * m)
}

/// Eneström–Kakeya annulus bounds $(r_{\min}, r_{\max})$ for a polynomial
/// with strictly positive real ascending coefficients $a_0, \dots, a_n > 0$.
pub fn enestrom_kakeya_bounds(real_positive_coeffs: &[f64]) -> Option<(f64, f64)> {
    if real_positive_coeffs.len() < 2 {
        return None;
    }
    if real_positive_coeffs
        .iter()
        .any(|&a| a <= 0.0 || !a.is_finite())
    {
        return None;
    }
    let mut r_min = f64::INFINITY;
    let mut r_max = 0.0_f64;
    for pair in real_positive_coeffs.windows(2) {
        let a_k = *pair.first()?;
        let a_k1 = *pair.get(1)?;
        let q = a_k / a_k1;
        r_min = r_min.min(q);
        r_max = r_max.max(q);
    }
    Some((r_min, r_max))
}

fn initial_circle_guesses(c: &[Complex<f64>]) -> Vec<Complex<f64>> {
    let n = c.len().saturating_sub(1);
    let r = cauchy_sharp_root_bound(c).unwrap_or(1.0).max(0.5);
    let center = if n >= 1 {
        let c_n = c.get(n).copied().unwrap_or(Complex::new(1.0, 0.0));
        let c_n1 = c.get(n.saturating_sub(1)).copied().unwrap_or_default();
        -c_n1 / (c_n * (n as f64))
    } else {
        Complex::new(0.0, 0.0)
    };
    (0..n)
        .map(|k| {
            let angle = 2.0 * PI * (k as f64) / (n as f64) + 0.137;
            center + Complex::from_polar(r * 0.6, angle)
        })
        .collect()
}

/// Finds all roots of a polynomial with ascending coefficients using the Durand–Kerner method.
pub fn durand_kerner(
    coeffs: &[Complex<f64>],
    max_iter: usize,
    tol: f64,
) -> Result<Vec<Complex<f64>>, String> {
    let c = trim_leading_zeros(coeffs);
    let n = c.len().saturating_sub(1);
    if n == 0 {
        return Ok(Vec::new());
    }
    let cn = c
        .last()
        .copied()
        .ok_or_else(|| "Empty coefficients".to_string())?;
    let monic: Vec<Complex<f64>> = c.iter().map(|&v| v / cn).collect();
    let mut roots = initial_circle_guesses(&monic);
    for _ in 0..max_iter {
        let mut max_step = 0.0_f64;
        for i in 0..n {
            let zi = *roots.get(i).ok_or_else(|| "Index error".to_string())?;
            let (p_val, _, _) = eval_poly_and_deriv(&monic, zi);
            let mut denom = Complex::new(1.0, 0.0);
            for (j, &zj) in roots.iter().enumerate() {
                if i != j {
                    denom *= zi - zj;
                }
            }
            if denom.norm() < 1e-30 {
                continue;
            }
            let delta = p_val / denom;
            max_step = max_step.max(delta.norm());
            if let Some(slot) = roots.get_mut(i) {
                *slot = zi - delta;
            }
        }
        if max_step < tol {
            return Ok(roots);
        }
    }
    Ok(roots)
}

/// Finds all roots of a polynomial with ascending coefficients using the Aberth–Ehrlich method.
pub fn aberth_ehrlich(
    coeffs: &[Complex<f64>],
    max_iter: usize,
    tol: f64,
) -> Result<Vec<Complex<f64>>, String> {
    let c = trim_leading_zeros(coeffs);
    let n = c.len().saturating_sub(1);
    if n == 0 {
        return Ok(Vec::new());
    }
    let cn = c
        .last()
        .copied()
        .ok_or_else(|| "Empty coefficients".to_string())?;
    let monic: Vec<Complex<f64>> = c.iter().map(|&v| v / cn).collect();
    let mut roots = initial_circle_guesses(&monic);
    for _ in 0..max_iter {
        let mut max_step = 0.0_f64;
        for i in 0..n {
            let zi = *roots.get(i).ok_or_else(|| "Index error".to_string())?;
            let (p_val, dp_val, _) = eval_poly_and_deriv(&monic, zi);
            if p_val.norm() < tol * 1e-3 {
                continue;
            }
            if dp_val.norm() < 1e-30 {
                if let Some(slot) = roots.get_mut(i) {
                    *slot = zi + Complex::new(1e-5, 1.3e-5);
                }
                continue;
            }
            let w = p_val / dp_val;
            let mut sum_inv = Complex::new(0.0, 0.0);
            for (j, &zj) in roots.iter().enumerate() {
                if i != j {
                    let diff = zi - zj;
                    if diff.norm() > 1e-30 {
                        sum_inv += Complex::new(1.0, 0.0) / diff;
                    }
                }
            }
            let denom = Complex::new(1.0, 0.0) - w * sum_inv;
            let delta = if denom.norm() > 1e-30 {
                w / denom
            } else {
                w
            };
            max_step = max_step.max(delta.norm());
            if let Some(slot) = roots.get_mut(i) {
                *slot = zi - delta;
            }
        }
        if max_step < tol {
            return Ok(roots);
        }
    }
    Ok(roots)
}

/// Finds a single root of a polynomial starting from `z0` via Laguerre's method.
pub fn laguerre_root(
    coeffs: &[Complex<f64>],
    z0: Complex<f64>,
    max_iter: usize,
    tol: f64,
) -> Result<Complex<f64>, String> {
    let c = trim_leading_zeros(coeffs);
    let n = c.len().saturating_sub(1);
    if n == 0 {
        return Err("Constant polynomial has no roots".to_string());
    }
    let n_f = n as f64;
    let mut z = z0;
    for _ in 0..max_iter {
        let (p, dp, ddp) = eval_poly_and_deriv(c, z);
        if p.norm() < tol {
            return Ok(z);
        }
        let g = dp / p;
        let h = g * g - ddp / p;
        let disc = ((n_f - 1.0) * (n_f * h - g * g)).sqrt();
        let d_plus = g + disc;
        let d_minus = g - disc;
        let denom = if d_plus.norm() >= d_minus.norm() {
            d_plus
        } else {
            d_minus
        };
        if denom.norm() < 1e-30 {
            return Err("Laguerre denominator vanished".to_string());
        }
        let a = Complex::new(n_f, 0.0) / denom;
        z -= a;
        if a.norm() < tol {
            return Ok(z);
        }
    }
    Ok(z)
}

/// Finds a complex root of an arbitrary analytic function `f` using Muller's quadratic interpolation.
pub fn muller_root<F>(
    f: F,
    mut z0: Complex<f64>,
    mut z1: Complex<f64>,
    mut z2: Complex<f64>,
    max_iter: usize,
    tol: f64,
) -> Result<Complex<f64>, String>
where
    F: Fn(Complex<f64>) -> Complex<f64>,
{
    for _ in 0..max_iter {
        let f0 = f(z0);
        let f1 = f(z1);
        let f2 = f(z2);
        if f2.norm() < tol {
            return Ok(z2);
        }
        let h1 = z1 - z0;
        let h2 = z2 - z1;
        if h1.norm() < 1e-30 || h2.norm() < 1e-30 {
            return Ok(z2);
        }
        let d1 = (f1 - f0) / h1;
        let d2 = (f2 - f1) / h2;
        let a = (d2 - d1) / (h2 + h1);
        let b = a * h2 + d2;
        let c = f2;
        let rad = (b * b - a * c * 4.0).sqrt();
        let den_plus = b + rad;
        let den_minus = b - rad;
        let den = if den_plus.norm() >= den_minus.norm() {
            den_plus
        } else {
            den_minus
        };
        if den.norm() < 1e-30 {
            return Err("Muller denominator vanished".to_string());
        }
        let dz = -c * 2.0 / den;
        let z3 = z2 + dz;
        if dz.norm() < tol {
            return Ok(z3);
        }
        z0 = z1;
        z1 = z2;
        z2 = z3;
    }
    Ok(z2)
}

/// Finds a complex root using Halley's cubically convergent method.
pub fn halley_root<F, DF, DDF>(
    f: F,
    df: DF,
    ddf: DDF,
    mut z: Complex<f64>,
    max_iter: usize,
    tol: f64,
) -> Result<Complex<f64>, String>
where
    F: Fn(Complex<f64>) -> Complex<f64>,
    DF: Fn(Complex<f64>) -> Complex<f64>,
    DDF: Fn(Complex<f64>) -> Complex<f64>,
{
    for _ in 0..max_iter {
        let fz = f(z);
        if fz.norm() < tol {
            return Ok(z);
        }
        let dfz = df(z);
        let ddfz = ddf(z);
        let denom = dfz * dfz * 2.0 - fz * ddfz;
        if denom.norm() < 1e-30 {
            return Err("Halley denominator vanished".to_string());
        }
        let dz = (fz * dfz * 2.0) / denom;
        z -= dz;
        if dz.norm() < tol {
            return Ok(z);
        }
    }
    Ok(z)
}

/// Computes power-sum moments $s_k = \frac{1}{2\pi i} \oint_{|z-c|=r} z^k \frac{f'(z)}{f(z)} dz$
/// for $k = 0, \dots, \text{max\_moment}$ using the periodic trapezoidal rule.
pub fn argument_principle_moments<F, DF>(
    f: F,
    df: DF,
    center: Complex<f64>,
    radius: f64,
    max_moment: usize,
    n_samples: usize,
) -> Result<Vec<Complex<f64>>, String>
where
    F: Fn(Complex<f64>) -> Complex<f64>,
    DF: Fn(Complex<f64>) -> Complex<f64>,
{
    if radius <= 0.0 || n_samples < 8 {
        return Err("Invalid circle radius or sample count".to_string());
    }
    let mut moments = vec![Complex::new(0.0, 0.0); max_moment.saturating_add(1)];
    let n_f = n_samples as f64;
    for j in 0..n_samples {
        let theta = 2.0 * PI * (j as f64) / n_f;
        let offset = Complex::from_polar(radius, theta);
        let z = center + offset;
        let fz = f(z);
        if fz.norm() < 1e-14 {
            return Err("Function has a zero/near-zero on the integration circle".to_string());
        }
        let log_deriv = df(z) / fz;
        // dz = i * offset * dtheta, divided by (2 * pi * i) -> offset / n_samples
        let weight = log_deriv * offset / n_f;
        let mut z_pow = Complex::new(1.0, 0.0);
        for slot in &mut moments {
            *slot += weight * z_pow;
            z_pow *= z;
        }
    }
    Ok(moments)
}

/// Computes the argument principle winding number $(Z - P) = \text{round}(\text{Re}(s_0))$.
pub fn argument_principle_winding<F, DF>(
    f: F,
    df: DF,
    center: Complex<f64>,
    radius: f64,
    n_samples: usize,
) -> Result<isize, String>
where
    F: Fn(Complex<f64>) -> Complex<f64>,
    DF: Fn(Complex<f64>) -> Complex<f64>,
{
    let m = argument_principle_moments(f, df, center, radius, 0, n_samples)?;
    let s0 = m
        .first()
        .copied()
        .ok_or_else(|| "Missing moment".to_string())?;
    if s0.im.abs() > 1e-2 {
        return Err(format!(
            "Non-integer winding number imaginary part: {}",
            s0.im
        ));
    }
    Ok(s0.re.round() as isize)
}

/// Delves–Lyness root localization inside $|z - \text{center}| < \text{radius}$ for an
/// analytic function without poles inside the disk. Uses argument principle moments,
/// Newton's identities to form the monic polynomial, and Aberth–Ehrlich to extract all zeros.
pub fn delves_lyness_roots<F, DF>(
    f: F,
    df: DF,
    center: Complex<f64>,
    radius: f64,
    n_samples: usize,
) -> Result<Vec<Complex<f64>>, String>
where
    F: Fn(Complex<f64>) -> Complex<f64>,
    DF: Fn(Complex<f64>) -> Complex<f64>,
{
    let s0_vec = argument_principle_moments(&f, &df, center, radius, 0, n_samples)?;
    let s0 = s0_vec
        .first()
        .copied()
        .ok_or_else(|| "Missing s0".to_string())?;
    let count_signed = s0.re.round() as isize;
    if count_signed < 0 {
        return Err("Negative winding number indicates poles inside contour".to_string());
    }
    let n = count_signed as usize;
    if n == 0 {
        return Ok(Vec::new());
    }
    // Shift origin to center for numerical stability of moments
    let shifted_moments = argument_principle_moments(
        |w| f(center + w),
        |w| df(center + w),
        Complex::new(0.0, 0.0),
        radius,
        n,
        n_samples,
    )?;
    // Newton identities: k * e_k = sum_{j=1}^k (-1)^{j-1} * e_{k-j} * s_j
    let mut e = vec![Complex::new(0.0, 0.0); n.saturating_add(1)];
    if let Some(e0) = e.first_mut() {
        *e0 = Complex::new(1.0, 0.0);
    }
    for k in 1..=n {
        let mut acc = Complex::new(0.0, 0.0);
        for j in 1..=k {
            let sign = if j % 2 == 1 { 1.0 } else { -1.0 };
            let e_prev = e.get(k.saturating_sub(j)).copied().unwrap_or_default();
            let s_j = shifted_moments.get(j).copied().unwrap_or_default();
            acc += e_prev * s_j * sign;
        }
        if let Some(ek) = e.get_mut(k) {
            *ek = acc / (k as f64);
        }
    }
    // Monic polynomial in ascending powers: sum_{m=0}^n (-1)^{n-m} e_{n-m} w^m
    let mut poly = Vec::with_capacity(n.saturating_add(1));
    for m in 0..=n {
        let idx = n.saturating_sub(m);
        let sign = if idx % 2 == 0 { 1.0 } else { -1.0 };
        let ek = e.get(idx).copied().unwrap_or_default();
        poly.push(ek * sign);
    }
    let shifted_roots = aberth_ehrlich(&poly, 200, 1e-13)?;
    Ok(shifted_roots.into_iter().map(|w| center + w).collect())
}

/// Checks Rouché's condition $|f(z)| > |g(z)|$ on $|z - \text{center}| = \text{radius}$.
pub fn rouche_check_circle<F, G>(
    f: F,
    g: G,
    center: Complex<f64>,
    radius: f64,
    n_samples: usize,
) -> Result<bool, String>
where
    F: Fn(Complex<f64>) -> Complex<f64>,
    G: Fn(Complex<f64>) -> Complex<f64>,
{
    if radius <= 0.0 || n_samples < 8 {
        return Err("Invalid circle radius or sample count".to_string());
    }
    let mut min_margin = f64::INFINITY;
    for j in 0..n_samples {
        let theta = 2.0 * PI * (j as f64) / (n_samples as f64);
        let z = center + Complex::from_polar(radius, theta);
        let margin = f(z).norm() - g(z).norm();
        min_margin = min_margin.min(margin);
    }
    Ok(min_margin > 0.0)
}

/// Tracks a root of a 1D complex homotopy $H(z(t), t) = 0$ from $t = 0$ to $t = 1$
/// starting from $z(0) = z_0$, using Euler predictor + Newton corrector.
pub fn homotopy_track_1d<H, DHZ, DHT>(
    h: H,
    dh_dz: DHZ,
    dh_dt: DHT,
    z0: Complex<f64>,
    n_steps: usize,
    newton_iters: usize,
    tol: f64,
) -> Result<Complex<f64>, String>
where
    H: Fn(Complex<f64>, f64) -> Complex<f64>,
    DHZ: Fn(Complex<f64>, f64) -> Complex<f64>,
    DHT: Fn(Complex<f64>, f64) -> Complex<f64>,
{
    if n_steps == 0 {
        return Err("n_steps must be positive".to_string());
    }
    let dt = 1.0 / (n_steps as f64);
    let mut z = z0;
    for step in 0..n_steps {
        let t = (step as f64) * dt;
        let t_next = ((step.saturating_add(1)) as f64) * dt;
        let hz = dh_dz(z, t);
        if hz.norm() < 1e-14 {
            return Err(format!("Singular Jacobian along homotopy at t = {t}"));
        }
        // Predictor: dz/dt = -(dH/dt)/(dH/dz)
        let dz_dt = -dh_dt(z, t) / hz;
        z += dz_dt * dt;
        // Corrector: Newton at t_next
        for _ in 0..newton_iters {
            let val = h(z, t_next);
            if val.norm() < tol {
                break;
            }
            let dval = dh_dz(z, t_next);
            if dval.norm() < 1e-14 {
                return Err(format!("Singular Jacobian in corrector at t = {t_next}"));
            }
            z -= val / dval;
        }
    }
    if h(z, 1.0).norm() > tol.max(1e-6) {
        return Err("Homotopy corrector failed to converge at t = 1".to_string());
    }
    Ok(z)
}
