//! Complex dynamics: orbits, periodic point classification & multipliers,
//! quadratic periodic points, Mandelbrot & Julia escape-time / Green's potential,
//! Böttcher & Koenigs linearization coordinates, Lyapunov exponents, and external rays.

use num_complex::Complex;
use std::f64::consts::PI;

/// Classification of a fixed or periodic point by its multiplier $|\lambda|$.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FixedPointKind {
    /// $\lambda = 0$ (critical periodic point).
    Superattracting,
    /// $0 < |\lambda| < 1$.
    Attracting,
    /// $|\lambda| = 1$.
    Neutral,
    /// $|\lambda| > 1$.
    Repelling,
}

/// Computes the forward orbit $[z_0, f(z_0), \dots, f^{\circ n}(z_0)]$.
pub fn orbit<F>(
    f: F,
    z0: Complex<f64>,
    n_steps: usize,
) -> Vec<Complex<f64>>
where
    F: Fn(Complex<f64>) -> Complex<f64>,
{
    let mut pts = Vec::with_capacity(n_steps.saturating_add(1));
    let mut z = z0;
    pts.push(z);
    for _ in 0..n_steps {
        z = f(z);
        pts.push(z);
    }
    pts
}

/// Computes the multiplier $\lambda = \prod_{k=0}^{p-1} f'(z_k)$ along a period-`p` orbit
/// starting at `z0`, and classifies the periodic point.
pub fn classify_periodic_point<F, DF>(
    f: F,
    df: DF,
    z0: Complex<f64>,
    period: usize,
) -> (Complex<f64>, FixedPointKind)
where
    F: Fn(Complex<f64>) -> Complex<f64>,
    DF: Fn(Complex<f64>) -> Complex<f64>,
{
    let mut z = z0;
    let mut mult = Complex::new(1.0, 0.0);
    for _ in 0..period.max(1) {
        mult *= df(z);
        z = f(z);
    }
    let m = mult.norm();
    let kind = if m < 1e-12 {
        FixedPointKind::Superattracting
    } else if m < 1.0 - 1e-12 {
        FixedPointKind::Attracting
    } else if m <= 1.0 + 1e-12 {
        FixedPointKind::Neutral
    } else {
        FixedPointKind::Repelling
    };
    (mult, kind)
}

/// Finds periodic points of $f_c(z) = z^2 + c$ of period `period` by running Newton's method
/// on $f_c^{\circ p}(z) - z = 0$ from the provided initial `seeds`, deduplicating within `tol * 1e3`.
pub fn find_periodic_points_quadratic(
    c: Complex<f64>,
    period: usize,
    seeds: &[Complex<f64>],
    max_iter: usize,
    tol: f64,
) -> Vec<Complex<f64>> {
    let p = period.max(1);
    let mut found: Vec<Complex<f64>> = Vec::new();
    let dedup_tol = (tol * 1e3).max(1e-7);
    for &seed in seeds {
        let mut z = seed;
        let mut converged = false;
        for _ in 0..max_iter {
            let mut w = z;
            let mut dw = Complex::new(1.0, 0.0);
            for _ in 0..p {
                dw = dw * w * 2.0;
                w = w * w + c;
            }
            let val = w - z;
            if val.norm() < tol {
                converged = true;
                break;
            }
            let dval = dw - Complex::new(1.0, 0.0);
            if dval.norm() < 1e-14 {
                break;
            }
            z -= val / dval;
        }
        if converged && found.iter().all(|&prev| (prev - z).norm() > dedup_tol) {
            found.push(z);
        }
    }
    found
}

fn quadratic_escape(
    z_start: Complex<f64>,
    c: Complex<f64>,
    max_iter: usize,
    escape_radius: f64,
) -> (usize, f64) {
    let r2 = escape_radius * escape_radius;
    let mut z = z_start;
    for n in 0..max_iter {
        let n2 = z.norm_sqr();
        if n2 > r2 {
            let log_z = 0.5 * n2.ln();
            let green = log_z * 2.0_f64.powi(-(n as i32));
            return (n, green);
        }
        z = z * z + c;
    }
    (max_iter, 0.0)
}

/// Mandelbrot escape time and Green's potential $G_c(c)$ for $z_0 = 0, z_{n+1} = z_n^2 + c$.
pub fn mandelbrot_escape_time(
    c: Complex<f64>,
    max_iter: usize,
    escape_radius: f64,
) -> (usize, f64) {
    quadratic_escape(Complex::new(0.0, 0.0), c, max_iter, escape_radius)
}

/// Julia escape time and Green's potential $G_c(z_0)$ for $z_{n+1} = z_n^2 + c$.
pub fn julia_escape_time(
    z0: Complex<f64>,
    c: Complex<f64>,
    max_iter: usize,
    escape_radius: f64,
) -> (usize, f64) {
    quadratic_escape(z0, c, max_iter, escape_radius)
}

/// Computes Böttcher's coordinate $\phi_c(z) = z \prod_{n=0}^{N-1} (1 + c / z_n^2)^{2^{-(n+1)}}$
/// in the basin of infinity of $f_c(z) = z^2 + c$, satisfying $\phi_c(z^2 + c) = \phi_c(z)^2$.
pub fn bottcher_coordinate(
    z: Complex<f64>,
    c: Complex<f64>,
    n_iters: usize,
) -> Result<Complex<f64>, String> {
    let mut zn = z;
    if zn.norm() <= 2.0_f64.max(c.norm() + 1.0) {
        // Still allow if orbit escapes; check after loop
    }
    let mut prod = Complex::new(1.0, 0.0);
    let mut weight = 0.5_f64;
    for _ in 0..n_iters {
        if zn.norm_sqr() < 1e-24 {
            return Err("Orbit hit 0 in Bottcher coordinate".to_string());
        }
        let factor = (Complex::new(1.0, 0.0) + c / (zn * zn)).powf(weight);
        prod *= factor;
        zn = zn * zn + c;
        weight *= 0.5;
        if zn.norm() > 1e150 {
            break;
        }
    }
    if zn.norm() <= 2.0 {
        return Err("Point is not clearly in the escaping basin of infinity".to_string());
    }
    Ok(z * prod)
}

/// Computes Koenigs's linearizing coordinate
/// $\psi(z) = \lim_{n\to\infty} \lambda^{-n} (f^{\circ n}(z) - z_0)$ near an attracting
/// fixed point $z_0$ ($0 < |\lambda| < 1$), satisfying $\psi(f(z)) = \lambda \psi(z)$.
pub fn koenigs_coordinate<F>(
    f: F,
    z0: Complex<f64>,
    lambda: Complex<f64>,
    z: Complex<f64>,
    n_iters: usize,
) -> Result<Complex<f64>, String>
where
    F: Fn(Complex<f64>) -> Complex<f64>,
{
    let lam_n = lambda.norm();
    if lam_n <= 0.0 || lam_n >= 1.0 {
        return Err("Koenigs coordinate requires 0 < |lambda| < 1".to_string());
    }
    let mut w = z;
    for _ in 0..n_iters {
        w = f(w);
    }
    let scale = lambda.powi(-(n_iters as i32));
    Ok(scale * (w - z0))
}

/// Estimates the Lyapunov exponent $\frac{1}{N} \sum_{n=\text{burn\_in}}^{\text{burn\_in}+N-1} \ln |f'(z_n)|$.
pub fn lyapunov_exponent<F, DF>(
    f: F,
    df: DF,
    z0: Complex<f64>,
    n_steps: usize,
    burn_in: usize,
) -> Result<f64, String>
where
    F: Fn(Complex<f64>) -> Complex<f64>,
    DF: Fn(Complex<f64>) -> Complex<f64>,
{
    if n_steps == 0 {
        return Err("n_steps must be positive".to_string());
    }
    let mut z = z0;
    for _ in 0..burn_in {
        z = f(z);
    }
    let mut sum = 0.0_f64;
    for _ in 0..n_steps {
        let d = df(z).norm();
        if d <= 0.0 || !d.is_finite() {
            return Err("Zero or non-finite derivative along orbit".to_string());
        }
        sum += d.ln();
        z = f(z);
    }
    Ok(sum / (n_steps as f64))
}

/// Traces a point on the external ray of $f_c(z) = z^2 + c$ at angle `theta_turns` $\in [0, 1)$
/// and Green potential `potential_r > 0`, satisfying $\phi_c(z) = e^{r + 2\pi i \theta}$,
/// by stepping inward from a large radius $R_0 = e^{2^k r}$ via $k$ square-root continuations.
pub fn external_ray_point(
    c: Complex<f64>,
    theta_turns: f64,
    potential_r: f64,
    n_steps: usize,
) -> Result<Complex<f64>, String> {
    if potential_r <= 0.0 || n_steps == 0 {
        return Err("potential_r and n_steps must be positive".to_string());
    }
    let k = n_steps.min(30);
    let r_top = potential_r * 2.0_f64.powi(k as i32);
    let theta_top = (theta_turns * 2.0_f64.powi(k as i32)).rem_euclid(1.0);
    // At large radius, phi_c(z) ≈ z
    let mut z = Complex::from_polar(r_top.exp(), 2.0 * PI * theta_top);
    // Pull back k times: z_{m-1}^2 + c = z_m, picking the root aligned with e^{2^{m-1}(r + 2pi i theta)}
    for m in (0..k).rev() {
        let r_m = potential_r * 2.0_f64.powi(m as i32);
        let th_m = (theta_turns * 2.0_f64.powi(m as i32)).rem_euclid(1.0);
        let approx = Complex::from_polar(r_m.exp(), 2.0 * PI * th_m);
        let root = (z - c).sqrt();
        z = if (root - approx).norm() <= (-root - approx).norm() {
            root
        } else {
            -root
        };
    }
    // Refine slightly with 4 Newton steps on bottcher_coordinate(z) = target
    let target = Complex::from_polar(potential_r.exp(), 2.0 * PI * theta_turns);
    for _ in 0..4 {
        let phi = bottcher_coordinate(z, c, 24)?;
        let err = phi - target;
        if err.norm() < 1e-12 {
            break;
        }
        let h = 1e-6;
        let dphi = (bottcher_coordinate(z + Complex::new(h, 0.0), c, 24)? - phi) / h;
        if dphi.norm() > 1e-14 {
            z -= err / dphi;
        }
    }
    Ok(z)
}
