//! Univalent function theory ($\mathcal{S}$ class), Bieberbach–de Branges conjecture checks,
//! Fekete–Szegő functional, Hankel determinants, Koebe growth/distortion/rotation theorems,
//! Schwarzian derivative, Nehari univalence criterion, Grunsky matrix, and Loewner slit ODE.

use crate::kernels::complex::cauchy_fft_derivatives;
use num_complex::Complex;
use std::f64::consts::PI;

/// Sharp radius of starlikeness for the class $\mathcal{S}$: $\tanh(\pi / 4) \approx 0.65579$.
pub const S_RADIUS_OF_STARLIKENESS: f64 = 0.655_794_202_632_672_4;

/// Sharp radius of convexity for the class $\mathcal{S}$: $2 - \sqrt{3} \approx 0.26795$.
pub const S_RADIUS_OF_CONVEXITY: f64 = 0.267_949_192_431_122_7;

/// Checks the Bieberbach–de Branges necessary condition $|a_n / a_1| \le n$ ($n \ge 2$)
/// for a normalized univalent function $f(z) = a_1 z + a_2 z^2 + \dots$ (with `coeffs[k]` = $a_k$).
pub fn check_bieberbach_de_branges(
    taylor_coeffs: &[Complex<f64>],
    tol: f64,
) -> bool {
    let a1 = match taylor_coeffs.get(1) {
        | Some(v) if v.norm() > 1e-15 => *v,
        | _ => return false,
    };
    let a1_norm = a1.norm();
    for (n, an) in taylor_coeffs.iter().enumerate().skip(2) {
        if (an.norm() / a1_norm) > (n as f64) + tol {
            return false;
        }
    }
    true
}

/// Computes the Fekete–Szegő functional $|a_3 - \mu a_2^2|$.
pub fn fekete_szego_functional(
    a2: Complex<f64>,
    a3: Complex<f64>,
    mu: Complex<f64>,
) -> f64 {
    (a3 - mu * a2 * a2).norm()
}

/// Computes the $q \times q$ Hankel determinant $H_q(n) = \det [c_{n + i + j}]_{i,j=0}^{q-1}$.
pub fn hankel_determinant(
    coeffs: &[Complex<f64>],
    q: usize,
    n: usize,
) -> Result<Complex<f64>, String> {
    if q == 0 {
        return Ok(Complex::new(1.0, 0.0));
    }
    let needed = n.saturating_add(2 * q.saturating_sub(1));
    if coeffs.len() <= needed {
        return Err("Insufficient coefficients for requested Hankel determinant".to_string());
    }
    let mut mat = vec![vec![Complex::new(0.0, 0.0); q]; q];
    for i in 0..q {
        for j in 0..q {
            let val = *coeffs.get(n + i + j).ok_or_else(|| "Index".to_string())?;
            if let Some(cell) = mat.get_mut(i).and_then(|r| r.get_mut(j)) {
                *cell = val;
            }
        }
    }
    // Gaussian elimination with partial pivoting for determinant
    let mut det = Complex::new(1.0, 0.0);
    for col in 0..q {
        let mut piv = col;
        let mut max_v = mat
            .get(col)
            .and_then(|r| r.get(col))
            .map(|c| c.norm())
            .unwrap_or(0.0);
        for row in (col + 1)..q {
            let v = mat
                .get(row)
                .and_then(|r| r.get(col))
                .map(|c| c.norm())
                .unwrap_or(0.0);
            if v > max_v {
                max_v = v;
                piv = row;
            }
        }
        if max_v < 1e-15 {
            return Ok(Complex::new(0.0, 0.0));
        }
        if piv != col {
            mat.swap(col, piv);
            det = -det;
        }
        let diag = mat
            .get(col)
            .and_then(|r| r.get(col))
            .copied()
            .unwrap_or_default();
        det *= diag;
        for row in (col + 1)..q {
            let factor = mat
                .get(row)
                .and_then(|r| r.get(col))
                .copied()
                .unwrap_or_default()
                / diag;
            for k in (col + 1)..q {
                let p_val = mat
                    .get(col)
                    .and_then(|r| r.get(k))
                    .copied()
                    .unwrap_or_default();
                if let Some(cell) = mat.get_mut(row).and_then(|r| r.get_mut(k)) {
                    *cell -= factor * p_val;
                }
            }
        }
    }
    Ok(det)
}

/// Koebe growth theorem bounds for $f \in \mathcal{S}$ at $|z| = r \in [0, 1)$:
/// $\frac{r}{(1+r)^2} \le |f(z)| \le \frac{r}{(1-r)^2}$.
pub fn koebe_growth_bounds(r: f64) -> Option<(f64, f64)> {
    if !(0.0..1.0).contains(&r) {
        return None;
    }
    let lo = r / ((1.0 + r) * (1.0 + r));
    let hi = r / ((1.0 - r) * (1.0 - r));
    Some((lo, hi))
}

/// Koebe distortion theorem bounds for $f \in \mathcal{S}$ at $|z| = r \in [0, 1)$:
/// $\frac{1-r}{(1+r)^3} \le |f'(z)| \le \frac{1+r}{(1-r)^3}$.
pub fn koebe_distortion_bounds(r: f64) -> Option<(f64, f64)> {
    if !(0.0..1.0).contains(&r) {
        return None;
    }
    let lo = (1.0 - r) / (1.0 + r).powi(3);
    let hi = (1.0 + r) / (1.0 - r).powi(3);
    Some((lo, hi))
}

/// Goluzin's sharp rotation theorem bound on $|\arg f'(z)|$ for $f \in \mathcal{S}$ at $|z| = r \in [0, 1)$:
/// $4 \arcsin(r)$ for $r \le 1/\sqrt{2}$, and $\pi + \ln\frac{r^2}{1 - r^2}$ for $1/\sqrt{2} < r < 1$.
pub fn koebe_rotation_bound(r: f64) -> Option<f64> {
    if !(0.0..1.0).contains(&r) {
        return None;
    }
    let thresh = std::f64::consts::FRAC_1_SQRT_2;
    if r <= thresh {
        Some(4.0 * r.asin())
    } else {
        Some(PI + (r * r / (1.0 - r * r)).ln())
    }
}

/// Computes the Schwarzian derivative $\{f, z\} = \frac{f'''(z)}{f'(z)} - \frac{3}{2} \left(\frac{f''(z)}{f'(z)}\right)^2$
/// using Cauchy–FFT derivatives on a circle of radius `h`.
pub fn schwarzian_derivative_numeric<F>(
    f: &F,
    z: Complex<f64>,
    h: f64,
) -> Result<Complex<f64>, String>
where
    F: Fn(Complex<f64>) -> Complex<f64>,
{
    let derivs = cauchy_fft_derivatives(f, z, h, 3, 64);
    let f1 = *derivs.get(1).ok_or_else(|| "Missing f'".to_string())?;
    let f2 = *derivs.get(2).ok_or_else(|| "Missing f''".to_string())?;
    let f3 = *derivs.get(3).ok_or_else(|| "Missing f'''".to_string())?;
    if f1.norm() < 1e-14 {
        return Err("Critical point f'(z) = 0 in Schwarzian derivative".to_string());
    }
    let ratio = f2 / f1;
    Ok(f3 / f1 - ratio * ratio * 1.5)
}

/// Checks Nehari's sufficient univalence criterion $(1 - |z|^2)^2 |\{f, z\}| \le 2$
/// across sample radii in $(0, 1)$ and `n_angles` angles per radius.
pub fn check_nehari_univalence<F>(
    f: &F,
    sample_radii: &[f64],
    n_angles: usize,
    h: f64,
) -> Result<bool, String>
where
    F: Fn(Complex<f64>) -> Complex<f64>,
{
    if n_angles == 0 {
        return Err("n_angles must be positive".to_string());
    }
    for &r in sample_radii {
        if !(0.0..1.0).contains(&r) {
            return Err("Sample radius must be in (0, 1)".to_string());
        }
        let weight = (1.0 - r * r).powi(2);
        for j in 0..n_angles {
            let theta = 2.0 * PI * (j as f64) / (n_angles as f64);
            let z = Complex::from_polar(r, theta);
            let s = schwarzian_derivative_numeric(f, z, h)?;
            if weight * s.norm() > 2.0 + 1e-9 {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

/// Numerically estimates the radius of starlikeness $\sup \{ r \in (0, r_{\max}] : \text{Re}(z f'(z)/f(z)) > 0 \text{ on } |z|=r \}$
/// via bisection.
pub fn radius_of_starlikeness_numeric<F, DF>(
    f: &F,
    df: &DF,
    r_max: f64,
    n_angles: usize,
) -> f64
where
    F: Fn(Complex<f64>) -> Complex<f64>,
    DF: Fn(Complex<f64>) -> Complex<f64>,
{
    let is_starlike_on_r = |r: f64| -> bool {
        if r <= 0.0 {
            return true;
        }
        for j in 0..n_angles {
            let theta = 2.0 * PI * (j as f64) / (n_angles as f64);
            let z = Complex::from_polar(r, theta);
            let fz = f(z);
            if fz.norm() < 1e-14 {
                return false;
            }
            let val = z * df(z) / fz;
            if val.re <= 0.0 {
                return false;
            }
        }
        true
    };
    let mut lo = 0.0_f64;
    let mut hi = r_max;
    for _ in 0..40 {
        let mid = 0.5 * (lo + hi);
        if is_starlike_on_r(mid) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    lo
}

/// Numerically estimates the radius of convexity $\sup \{ r \in (0, r_{\max}] : 1 + \text{Re}(z f''(z)/f'(z)) > 0 \text{ on } |z|=r \}$
/// via bisection.
pub fn radius_of_convexity_numeric<DF, DDF>(
    df: &DF,
    ddf: &DDF,
    r_max: f64,
    n_angles: usize,
) -> f64
where
    DF: Fn(Complex<f64>) -> Complex<f64>,
    DDF: Fn(Complex<f64>) -> Complex<f64>,
{
    let is_convex_on_r = |r: f64| -> bool {
        if r <= 0.0 {
            return true;
        }
        for j in 0..n_angles {
            let theta = 2.0 * PI * (j as f64) / (n_angles as f64);
            let z = Complex::from_polar(r, theta);
            let dfz = df(z);
            if dfz.norm() < 1e-14 {
                return false;
            }
            let val = 1.0 + (z * ddf(z) / dfz).re;
            if val <= 0.0 {
                return false;
            }
        }
        true
    };
    let mut lo = 0.0_f64;
    let mut hi = r_max;
    for _ in 0..40 {
        let mid = 0.5 * (lo + hi);
        if is_convex_on_r(mid) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    lo
}

/// Computes the $M \times M$ Grunsky matrix $G_{m,n} = \sqrt{mn} \, c_{m,n}$ for $1 \le m, n \le M$
/// from the Taylor coefficients $f(z) = z + a_2 z^2 + \dots$ (`taylor_coeffs[k]` = $a_k$, $a_1 = 1$),
/// defined via $\ln \frac{f(z) - f(\zeta)}{z - \zeta} = -\sum_{m,n \ge 1} c_{m,n} z^m \zeta^n$.
pub fn grunsky_matrix(
    taylor_coeffs: &[Complex<f64>],
    order: usize,
) -> Result<Vec<Vec<Complex<f64>>>, String> {
    if order == 0 {
        return Ok(Vec::new());
    }
    if taylor_coeffs.len() < 2 {
        return Err("Need at least a_0, a_1 coefficients".to_string());
    }
    let a1 = *taylor_coeffs
        .get(1)
        .ok_or_else(|| "Missing a_1".to_string())?;
    if a1.norm() < 1e-15 {
        return Err("a_1 must be non-zero".to_string());
    }
    // Normalized coefficients b_k = a_k / a_1 for k >= 1 (b_1 = 1)
    let max_deg = 2 * order + 1;
    if taylor_coeffs.len() <= max_deg {
        return Err(format!(
            "grunsky_matrix(order={order}) requires at least {} Taylor coefficients (got {})",
            max_deg + 1,
            taylor_coeffs.len()
        ));
    }
    let b = |k: usize| -> Complex<f64> {
        if k == 0 {
            Complex::new(0.0, 0.0)
        } else {
            taylor_coeffs.get(k).copied().unwrap_or_default() / a1
        }
    };

    // Difference quotient Q(z, w) = (f(z) - f(w)) / (z - w) = sum_{p,q >= 0} q_{p,q} z^p w^q
    // where q_{p,q} = b_{p+q+1}
    // Write Q(z, w) = 1 + R(z, w) where R has no constant term (q_{0,0} = b_1 = 1).
    let sz = order + 1;
    let mut r_mat = vec![vec![Complex::new(0.0, 0.0); sz]; sz];
    for p in 0..sz {
        for q in 0..sz {
            if p == 0 && q == 0 {
                continue;
            }
            if p + q + 1 <= max_deg {
                if let Some(cell) = r_mat.get_mut(p).and_then(|r| r.get_mut(q)) {
                    *cell = b(p + q + 1);
                }
            }
        }
    }

    // ln(1 + R) = sum_{k=1}^{2*order} (-1)^{k-1} R^k / k truncated to degree (order, order)
    let mut log_mat = vec![vec![Complex::new(0.0, 0.0); sz]; sz];
    let mut r_pow = r_mat.clone();
    for k in 1..=(2 * order) {
        let coeff = (if k % 2 == 1 { 1.0 } else { -1.0 }) / (k as f64);
        for p in 0..sz {
            for q in 0..sz {
                let v = r_pow
                    .get(p)
                    .and_then(|r| r.get(q))
                    .copied()
                    .unwrap_or_default();
                if let Some(cell) = log_mat.get_mut(p).and_then(|r| r.get_mut(q)) {
                    *cell += v * coeff;
                }
            }
        }
        // Multiply r_pow = r_pow * r_mat truncated to sz x sz
        let mut next = vec![vec![Complex::new(0.0, 0.0); sz]; sz];
        for p1 in 0..sz {
            for q1 in 0..sz {
                let v1 = r_pow
                    .get(p1)
                    .and_then(|r| r.get(q1))
                    .copied()
                    .unwrap_or_default();
                if v1.norm() == 0.0 {
                    continue;
                }
                for p2 in 0..(sz - p1) {
                    for q2 in 0..(sz - q1) {
                        let v2 = r_mat
                            .get(p2)
                            .and_then(|r| r.get(q2))
                            .copied()
                            .unwrap_or_default();
                        if let Some(cell) = next.get_mut(p1 + p2).and_then(|r| r.get_mut(q1 + q2)) {
                            *cell += v1 * v2;
                        }
                    }
                }
            }
        }
        r_pow = next;
    }

    // G_{m,n} = sqrt(m*n) * (-log_mat[m][n]) for 1 <= m, n <= order
    let mut g = vec![vec![Complex::new(0.0, 0.0); order]; order];
    for m in 1..=order {
        for n in 1..=order {
            let c_mn = -log_mat
                .get(m)
                .and_then(|r| r.get(n))
                .copied()
                .unwrap_or_default();
            let scale = ((m * n) as f64).sqrt();
            if let Some(cell) = g.get_mut(m - 1).and_then(|r| r.get_mut(n - 1)) {
                *cell = c_mn * scale;
            }
        }
    }
    Ok(g)
}

/// Verifies the Grunsky inequality $\left|\sum_{m,n=1}^M G_{m,n} x_m x_n\right| \le \sum_{m=1}^M |x_m|^2 + \text{tol}$.
pub fn check_grunsky_inequality(
    grunsky: &[Vec<Complex<f64>>],
    test_vec: &[Complex<f64>],
    tol: f64,
) -> Result<bool, String> {
    let m = grunsky.len();
    if test_vec.len() != m {
        return Err("Vector length must match Grunsky matrix dimension".to_string());
    }
    let mut quad = Complex::new(0.0, 0.0);
    let mut norm2 = 0.0_f64;
    for (i, &xi) in test_vec.iter().enumerate() {
        norm2 += xi.norm_sqr();
        let row = grunsky.get(i).ok_or_else(|| "Row".to_string())?;
        for (j, &xj) in test_vec.iter().enumerate() {
            let g_ij = *row.get(j).ok_or_else(|| "Col".to_string())?;
            quad += g_ij * xi * xj;
        }
    }
    if quad.norm() > norm2 + tol {
        return Ok(false);
    }
    // Also estimate spectral norm ||G||_2 via power iteration on G^* G
    if m > 0 {
        let mut v = vec![Complex::new(1.0 / (m as f64).sqrt(), 0.0); m];
        let mut sigma_sq = 0.0_f64;
        for _ in 0..24 {
            let mut gv = vec![Complex::new(0.0, 0.0); m];
            for i in 0..m {
                for j in 0..m {
                    gv[i] += grunsky[i][j] * v[j];
                }
            }
            let mut ghgv = vec![Complex::new(0.0, 0.0); m];
            for j in 0..m {
                for i in 0..m {
                    ghgv[j] += grunsky[i][j].conj() * gv[i];
                }
            }
            let norm = ghgv.iter().map(|c| c.norm_sqr()).sum::<f64>().sqrt();
            if norm < 1e-30 {
                sigma_sq = 0.0;
                break;
            }
            sigma_sq = norm;
            for j in 0..m {
                v[j] = ghgv[j] / norm;
            }
        }
        if sigma_sq.sqrt() > 1.0 + tol {
            return Ok(false);
        }
    }
    Ok(true)
}

/// Advances one RK4 step of the radial Loewner slit ODE
/// $\frac{dw}{dt} = w \frac{e^{i\kappa} + w}{e^{i\kappa} - w}$ with constant driving `driving_kappa` over `dt`.
pub fn loewner_slit_step(
    w: Complex<f64>,
    driving_kappa: f64,
    dt: f64,
) -> Complex<f64> {
    let u = Complex::from_polar(1.0, driving_kappa);
    if (u - w).norm() < 1e-12 {
        return w;
    }
    let rhs = |z: Complex<f64>| -> Complex<f64> { z * (u + z) / (u - z) };
    let k1 = rhs(w);
    let k2 = rhs(w + k1 * (0.5 * dt));
    let k3 = rhs(w + k2 * (0.5 * dt));
    let k4 = rhs(w + k3 * dt);
    w + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (dt / 6.0)
}

/// Checks Becker's sufficient univalence criterion `(1 - |z|^2) * |z * f''(z)/f'(z)| <= 1`.
#[must_use]
pub fn check_becker_univalence(
    z: Complex<f64>,
    df: Complex<f64>,
    ddf: Complex<f64>,
) -> bool {
    let r2 = z.norm_sqr();
    if r2 >= 1.0 || df.norm() < 1e-15 {
        return false;
    }
    (1.0 - r2) * (z * ddf / df).norm() <= 1.0
}
