//! Conformal mappings, Möbius geometry, Joukowsky/Kármán–Trefftz/Koebe maps,
//! Schwarz–Christoffel upper half-plane solver, Szegő/Kerzman–Stein disk mapping,
//! and Bergman kernel utilities.

use crate::kernels::complex::Mobius;
use crate::kernels::complex::check_cauchy_riemann_numeric;
use crate::kernels::complex::complex_derivative;
use crate::kernels::complex_scalar::Complex as ScalarComplex;
use crate::kernels::integrate::gauss_kronrod_any;
use num_complex::Complex;
use std::f64::consts::PI;

/// `cross_ratio(z1, z2, z3, z4)` computes `((z1 - z3) * (z2 - z4)) / ((z1 - z4) * (z2 - z3))`.
#[must_use]
pub fn cross_ratio(
    z1: Complex<f64>,
    z2: Complex<f64>,
    z3: Complex<f64>,
    z4: Complex<f64>,
) -> Complex<f64> {
    ((z1 - z3) * (z2 - z4)) / ((z1 - z4) * (z2 - z3))
}

fn mobius_to_01inf(z1: Complex<f64>, z2: Complex<f64>, z3: Complex<f64>) -> Mobius {
    // T(z) = ((z - z1)(z2 - z3)) / ((z - z3)(z2 - z1))
    let k1 = z2 - z3;
    let k2 = z2 - z1;
    Mobius {
        a: k1,
        b: -z1 * k1,
        c: k2,
        d: -z3 * k2,
    }
}

/// Constructs the unique Möbius transformation mapping `(z1, z2, z3) -> (w1, w2, w3)`.
///
/// # Errors
///
/// Returns an error if the source triple or the target triple contains repeated points.
pub fn mobius_from_three_points(
    z1: Complex<f64>,
    z2: Complex<f64>,
    z3: Complex<f64>,
    w1: Complex<f64>,
    w2: Complex<f64>,
    w3: Complex<f64>,
) -> Result<Mobius, String> {
    if (z1 - z2).norm() < 1e-14 || (z2 - z3).norm() < 1e-14 || (z1 - z3).norm() < 1e-14 {
        return Err("Source points must be distinct".to_string());
    }
    if (w1 - w2).norm() < 1e-14 || (w2 - w3).norm() < 1e-14 || (w1 - w3).norm() < 1e-14 {
        return Err("Target points must be distinct".to_string());
    }
    let tz = mobius_to_01inf(z1, z2, z3);
    let tw = mobius_to_01inf(w1, w2, w3);
    let tw_inv = tw.inverse();
    Ok(tw_inv.compose(&tz))
}

/// Classification of a non-singular Möbius transformation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MobiusClass {
    /// Identity map.
    Identity,
    /// Conjugated to a rotation ($\text{tr}^2 \in [0, 4)$).
    Elliptic,
    /// Single fixed point ($\text{tr}^2 = 4$).
    Parabolic,
    /// Conjugated to a pure scaling ($\text{tr}^2 \in (4, \infty)$).
    Hyperbolic,
    /// General spiral scaling ($\text{tr}^2 \notin [0, \infty)$).
    Loxodromic,
}

/// Classifies a Möbius transformation by its normalized trace squared $\tau = (a+d)^2 / (ad-bc)$.
///
/// # Errors
///
/// Returns an error for a degenerate transformation.
pub fn classify_mobius(m: &Mobius) -> Result<MobiusClass, String> {
    let det = m.a * m.d - m.b * m.c;
    if det.norm() < 1e-15 {
        return Err("Degenerate Möbius transformation".to_string());
    }
    if m.b.norm() < 1e-12 && m.c.norm() < 1e-12 && (m.a - m.d).norm() < 1e-12 {
        return Ok(MobiusClass::Identity);
    }
    let tr = m.a + m.d;
    let tau = (tr * tr) / det;
    if (tau - Complex::new(4.0, 0.0)).norm() < 1e-9 {
        return Ok(MobiusClass::Parabolic);
    }
    if tau.im.abs() < 1e-9 {
        if tau.re >= 0.0 && tau.re < 4.0 {
            return Ok(MobiusClass::Elliptic);
        }
        if tau.re > 4.0 {
            return Ok(MobiusClass::Hyperbolic);
        }
    }
    Ok(MobiusClass::Loxodromic)
}

/// Returns the finite fixed points of $m(z) = z$, i.e. roots of $c z^2 + (d - a) z - b = 0$.
#[must_use]
pub fn mobius_fixed_points(m: &Mobius) -> Vec<Complex<f64>> {
    if m.c.norm() < 1e-14 {
        let diff = m.d - m.a;
        if diff.norm() < 1e-14 {
            return Vec::new();
        }
        return vec![m.b / diff];
    }
    let da = m.d - m.a;
    let disc = (da * da + m.c * m.b * 4.0).sqrt();
    let two_c = m.c * 2.0;
    // Citardauq / stable quadratic roots avoiding cancellation
    let b_plus = -da + disc;
    let b_minus = -da - disc;
    let (z1, z2) = if b_plus.norm() >= b_minus.norm() {
        let r1 = b_plus / two_c;
        let r2 = if r1.norm() > 1e-30 {
            (-m.b / m.c) / r1
        } else {
            b_minus / two_c
        };
        (r1, r2)
    } else {
        let r1 = b_minus / two_c;
        let r2 = if r1.norm() > 1e-30 {
            (-m.b / m.c) / r1
        } else {
            b_plus / two_c
        };
        (r1, r2)
    };
    if (z1 - z2).norm() < 1e-12 {
        vec![z1]
    } else {
        vec![z1, z2]
    }
}

/// Unit disk automorphism $f(z) = e^{i\theta} \frac{z - \alpha}{1 - \bar{\alpha} z}$ for $|\alpha| < 1$.
///
/// # Errors
///
/// Returns an error if the translation parameter lies on or outside the unit circle.
pub fn disk_automorphism(alpha: Complex<f64>, theta: f64) -> Result<Mobius, String> {
    if alpha.norm() >= 1.0 {
        return Err("Disk automorphism requires |alpha| < 1".to_string());
    }
    let rot = Complex::from_polar(1.0, theta);
    Ok(Mobius {
        a: rot,
        b: -rot * alpha,
        c: -alpha.conj(),
        d: Complex::new(1.0, 0.0),
    })
}

/// `joukowsky(z)` returns the Joukowsky airfoil transformation `(z + 1/z)/2`.
#[must_use]
pub fn joukowsky(z: Complex<f64>) -> Complex<f64> {
    let reciprocal_half = ScalarComplex::new(0.5, 0.0).div(ScalarComplex::new(z.re, z.im));
    z * 0.5 + Complex::new(reciprocal_half.re, reciprocal_half.im)
}

/// Two branches of the inverse Joukowsky map $w \pm \sqrt{w^2 - 1}$.
#[must_use]
pub fn joukowsky_inverse(w: Complex<f64>) -> (Complex<f64>, Complex<f64>) {
    let scale = w.re.abs().max(w.im.abs());
    // Form the unscaled discriminant near the branch points at ±1.
    if !scale.is_finite() || scale <= 8.0 {
        let rad = (w * w - Complex::new(1.0, 0.0)).sqrt();
        return (w + rad, w - rad);
    }
    let normalized = w / scale;
    let reciprocal_scale = scale.recip();
    let rad = (normalized * normalized - reciprocal_scale * reciprocal_scale).sqrt();
    let plus = normalized + rad;
    let minus = normalized - rad;
    // The roots have product one. Take the reciprocal of the larger
    // root before undoing the scaling so the smaller one survives even
    // if the larger one overflows.
    let plus_is_larger = plus.re.abs().max(plus.im.abs()) >= minus.re.abs().max(minus.im.abs());
    let large = if plus_is_larger { plus } else { minus };
    let inverse = ScalarComplex::new(large.re, large.im).recip();
    let small = Complex::new(inverse.re, inverse.im) * reciprocal_scale;
    let large = large * scale;
    if plus_is_larger {
        (large, small)
    } else {
        (small, large)
    }
}

/// `karman_trefftz(z, exponent)` evaluates the Kármán–Trefftz map with principal complex powers.
/// For exponent `n`, the map is `n * ((z+1)^n + (z-1)^n) / ((z+1)^n - (z-1)^n)`.
#[must_use]
pub fn karman_trefftz(z: Complex<f64>, n: f64) -> Complex<f64> {
    let zp = (z + Complex::new(1.0, 0.0)).powf(n);
    let zm = (z - Complex::new(1.0, 0.0)).powf(n);
    Complex::new(n, 0.0) * (zp + zm) / (zp - zm)
}

/// `koebe_map(z)` returns the Koebe function `z / (1 - z)^2`.
#[must_use]
pub fn koebe_map(z: Complex<f64>) -> Complex<f64> {
    let one_m_z = Complex::new(1.0, 0.0) - z;
    z / (one_m_z * one_m_z)
}

/// Numerically verifies conformality at `z`: Cauchy–Riemann holds and $|f'(z)| > \text{tol}$.
pub fn is_conformal_numeric<F>(f: &F, z: Complex<f64>, h: f64, tol: f64) -> bool
where
    F: Fn(Complex<f64>) -> Complex<f64>,
{
    let cr_ok = check_cauchy_riemann_numeric(f, z, h, tol);
    let df = complex_derivative(f, z);
    cr_ok && df.norm() > tol
}

/// `cayley_transform(z)` returns the Möbius map `(z - i) / (z + i)`.
/// Its pole is at `-i`. On the upper half-plane its image is the open unit disk.
#[must_use]
pub fn cayley_transform(z: Complex<f64>) -> Complex<f64> {
    crate::kernels::complex_coords::uhp_to_disk(z).unwrap_or(Complex::new(f64::NAN, f64::NAN))
}

/// Returns the local conformal magnification $|f'(z)|$ and rotation angle $\arg f'(z)$.
pub fn conformal_scale_and_rotation<F>(f: &F, z: Complex<f64>, h: f64) -> (f64, f64)
where
    F: Fn(Complex<f64>) -> Complex<f64>,
{
    let df = if h.is_finite() && h > 0.0 {
        crate::kernels::complex::wirtinger_dz_numeric(f, z, h)
    } else {
        complex_derivative(f, z)
    };
    (df.norm(), df.arg())
}

fn integrate_complex_01<F>(f: F, tol: f64, max_evals: usize) -> Complex<f64>
where
    F: Fn(f64) -> Complex<f64>,
{
    let re = gauss_kronrod_any(|t| f(t).re, 0.0, 1.0, tol, max_evals).value;
    let im = gauss_kronrod_any(|t| f(t).im, 0.0, 1.0, tol, max_evals).value;
    Complex::new(re, im)
}

/// Schwarz–Christoffel map from the upper half-plane $\mathbb{H}^+$ to a polygon:
/// $f(z) = \text{shift} + \text{scale} \int_{`z_0}^z` \prod_{k=1}^m (\zeta - `x_k)^{\alpha_k` - 1} d\zeta$.
#[derive(Debug, Clone)]
pub struct SchwarzChristoffelMap {
    /// Prevertices $`x_1` < `x_2` < \dots < `x_m`$ on the real axis.
    pub prevertices: Vec<f64>,
    /// Interior angle factors `\alpha_k \in (0, 2)` (so exponent is `\alpha_k - 1`).
    pub interior_angles_pi: Vec<f64>,
    /// Complex scaling factor $C$.
    pub scale: Complex<f64>,
    /// Additive offset $`f(z_0)`$ at basepoint $`z_0` = i$.
    pub shift: Complex<f64>,
}

impl SchwarzChristoffelMap {
    fn integrand(&self, zeta: Complex<f64>) -> Complex<f64> {
        let mut prod = Complex::new(1.0, 0.0);
        for (&xk, &alpha) in self.prevertices.iter().zip(self.interior_angles_pi.iter()) {
            let diff = zeta - Complex::new(xk, 0.0);
            // Principal branch for Im(zeta) > 0 keeps diff in upper half-plane
            prod *= diff.powf(alpha - 1.0);
        }
        prod
    }

    /// Evaluates the Schwarz–Christoffel integral from $`z_0` = i$ to $z$ (with $\text{Im}(z) > 0$).
    ///
    /// # Errors
    ///
    /// This method currently returns `Ok` for every evaluation.
    pub fn eval(&self, z: Complex<f64>) -> Result<Complex<f64>, String> {
        let z0 = Complex::new(0.0, 1.0);
        let dir = z - z0;
        let integral = integrate_complex_01(
            |t| {
                let zeta = z0 + dir * t;
                self.integrand(zeta) * dir
            },
            1e-10,
            200,
        );
        Ok(self.shift + self.scale * integral)
    }

    fn segment_integral(&self, a: f64, b: f64) -> Complex<f64> {
        // Compound tanh-sinh-like / double-sin^2 substitution for acute endpoint singularities:
        // x(u) = a + span * sin^2(pi/2 * sin^2(pi/2 * u)) clusters doubly near endpoints.
        let span = b - a;
        integrate_complex_01(
            |u| {
                let s1 = (0.5 * PI * u).sin();
                let c1 = (0.5 * PI * u).cos();
                let v = s1 * s1;
                let dv = PI * s1 * c1;
                let s2 = (0.5 * PI * v).sin();
                let c2 = (0.5 * PI * v).cos();
                let x = a + span * s2 * s2;
                let dx = span * PI * s2 * c2 * dv;
                self.integrand(Complex::new(x, 1e-14)) * dx
            },
            1e-10,
            400,
        )
    }
}

fn polygon_interior_angles_pi(vertices: &[Complex<f64>]) -> Result<Vec<f64>, String> {
    let n = vertices.len();
    if n < 3 {
        return Err("Polygon must have at least 3 vertices".to_string());
    }
    let mut alphas = Vec::with_capacity(n);
    for i in 0..n {
        let prev = *vertices
            .get((i + n - 1) % n)
            .ok_or_else(|| "Index".to_string())?;
        let curr = *vertices.get(i).ok_or_else(|| "Index".to_string())?;
        let next = *vertices
            .get((i + 1) % n)
            .ok_or_else(|| "Index".to_string())?;
        let v_in = curr - prev;
        let v_out = next - curr;
        // Exterior turning angle for CCW polygon
        let turn = (v_out / v_in).arg();
        let alpha = 1.0 - turn / PI;
        alphas.push(alpha);
    }
    Ok(alphas)
}

/// Solves for the prevertices of a Schwarz–Christoffel map.
///
/// Maps the upper half-plane to the simple counterclockwise polygon `polygon_vertices`. Fixes
/// the first three prevertices at `(-1, 0, 1)` and finds the remaining `n - 3` by damped Newton
/// iteration.
///
/// # Errors
///
/// Returns an error for an invalid polygon or when the quadrature or Newton solve fails.
pub fn solve_sc_uhp(
    polygon_vertices: &[Complex<f64>],
    max_iter: usize,
    tol: f64,
) -> Result<SchwarzChristoffelMap, String> {
    let n = polygon_vertices.len();
    if n < 3 {
        return Err("Need at least 3 polygon vertices".to_string());
    }
    let alphas = polygon_interior_angles_pi(polygon_vertices)?;
    let mut prevertices = Vec::with_capacity(n);
    prevertices.push(-1.0);
    prevertices.push(0.0);
    for k in 2..n.saturating_sub(1) {
        prevertices.push(((k - 1) as f64) / ((n - 2) as f64));
    }
    prevertices.push(1.0);

    let mut target_lengths = Vec::with_capacity(n.saturating_sub(1));
    for i in 0..n.saturating_sub(1) {
        let a = *polygon_vertices.get(i).ok_or_else(|| "Index".to_string())?;
        let b = *polygon_vertices
            .get(i + 1)
            .ok_or_else(|| "Index".to_string())?;
        target_lengths.push((b - a).norm());
    }
    let len0 = *target_lengths
        .first()
        .ok_or_else(|| "Missing edge".to_string())?;
    let target_ratios: Vec<f64> = target_lengths.iter().skip(1).map(|&l| l / len0).collect();

    // Unconstrained coordinates u_k parameterizing 0 < x_2 < ... < x_{n-2} < 1 via softmax gaps
    let m_free = n.saturating_sub(3);
    let mut u = vec![0.0_f64; m_free];

    let build_prevertices = |u_vec: &[f64]| -> Vec<f64> {
        let mut xs = vec![-1.0, 0.0];
        let mut gaps = Vec::with_capacity(u_vec.len() + 1);
        let mut total = 1.0_f64; // last gap has exp(0) = 1
        for &uk in u_vec {
            let g = uk.clamp(-12.0, 12.0).exp();
            total += g;
            gaps.push(g);
        }
        let mut acc = 0.0_f64;
        for g in gaps {
            acc += g / total;
            xs.push(acc);
        }
        xs.push(1.0);
        xs
    };

    let eval_ratios = |xs: &[f64]| -> Result<(Vec<f64>, Complex<f64>), String> {
        let sc = SchwarzChristoffelMap {
            prevertices: xs.to_vec(),
            interior_angles_pi: alphas.clone(),
            scale: Complex::new(1.0, 0.0),
            shift: Complex::new(0.0, 0.0),
        };
        let mut segs = Vec::with_capacity(n.saturating_sub(1));
        for i in 0..n.saturating_sub(1) {
            let a = *xs.get(i).ok_or_else(|| "Index".to_string())?;
            let b = *xs.get(i + 1).ok_or_else(|| "Index".to_string())?;
            segs.push(sc.segment_integral(a, b));
        }
        let s0 = *segs.first().ok_or_else(|| "Missing seg".to_string())?;
        let s0_norm = s0.norm().max(1e-15);
        let ratios: Vec<f64> = segs.iter().skip(1).map(|s| s.norm() / s0_norm).collect();
        Ok((ratios, s0))
    };

    let mut converged = m_free == 0;
    if m_free > 0 {
        for _ in 0..max_iter {
            let xs = build_prevertices(&u);
            let (ratios, _) = eval_ratios(&xs)?;
            let mut res = Vec::with_capacity(m_free);
            let mut norm2 = 0.0_f64;
            for k in 0..m_free {
                let r_k = *ratios.get(k).ok_or_else(|| "Ratio index".to_string())?;
                let t_k = *target_ratios
                    .get(k)
                    .ok_or_else(|| "Target index".to_string())?;
                let diff = r_k - t_k;
                norm2 += diff * diff;
                res.push(diff);
            }
            if norm2.sqrt() < tol {
                converged = true;
                break;
            }
            // Finite-difference Jacobian
            let eps = 1e-5;
            let mut jac = vec![vec![0.0_f64; m_free]; m_free];
            for col in 0..m_free {
                let mut u_pert = u.clone();
                if let Some(slot) = u_pert.get_mut(col) {
                    *slot += eps;
                }
                let xs_p = build_prevertices(&u_pert);
                let (ratios_p, _) = eval_ratios(&xs_p)?;
                for row in 0..m_free {
                    let r_p = *ratios_p.get(row).ok_or_else(|| "Ratio p".to_string())?;
                    let r_0 = *ratios.get(row).ok_or_else(|| "Ratio 0".to_string())?;
                    if let Some(cell) = jac.get_mut(row).and_then(|r| r.get_mut(col)) {
                        *cell = (r_p - r_0) / eps;
                    }
                }
            }
            let delta = solve_real_linear(&jac, &res)?;
            let mut alpha_step = 1.0_f64;
            while alpha_step > 1e-4 {
                let u_try: Vec<f64> = u
                    .iter()
                    .zip(delta.iter())
                    .map(|(&uk, &dk)| uk - alpha_step * dk)
                    .collect();
                let xs_try = build_prevertices(&u_try);
                let (ratios_try, _) = eval_ratios(&xs_try)?;
                let mut try_norm2 = 0.0_f64;
                for k in 0..m_free {
                    let rk = *ratios_try.get(k).unwrap_or(&0.0);
                    let tk = *target_ratios.get(k).unwrap_or(&0.0);
                    try_norm2 += (rk - tk) * (rk - tk);
                }
                if try_norm2 < norm2 {
                    u = u_try;
                    break;
                }
                alpha_step *= 0.5;
            }
        }
        prevertices = build_prevertices(&u);
    }
    if !converged {
        return Err("Schwarz-Christoffel solver failed to converge within max_iter".to_string());
    }

    let (_, s0) = eval_ratios(&prevertices)?;
    let v0 = *polygon_vertices
        .first()
        .ok_or_else(|| "Index".to_string())?;
    let v1 = *polygon_vertices.get(1).ok_or_else(|| "Index".to_string())?;
    let scale = (v1 - v0) / s0;
    let mut sc = SchwarzChristoffelMap {
        prevertices,
        interior_angles_pi: alphas,
        scale,
        shift: Complex::new(0.0, 0.0),
    };
    // Integrate from z0 = i to a point near x0 = -1 along an interior approach
    let f_at_x0 = sc.eval(Complex::new(-1.0, 1e-6))?;
    sc.shift = v0 - f_at_x0;
    Ok(sc)
}

fn solve_real_linear(a: &[Vec<f64>], b: &[f64]) -> Result<Vec<f64>, String> {
    let n = b.len();
    let mut mat = a.to_vec();
    let mut rhs = b.to_vec();
    for col in 0..n {
        let mut piv = col;
        let mut max_v = mat
            .get(col)
            .and_then(|r| r.get(col))
            .copied()
            .unwrap_or(0.0)
            .abs();
        for row in (col + 1)..n {
            let v = mat
                .get(row)
                .and_then(|r| r.get(col))
                .copied()
                .unwrap_or(0.0)
                .abs();
            if v > max_v {
                max_v = v;
                piv = row;
            }
        }
        if max_v < 1e-14 {
            return Err("Singular Jacobian in SC solver".to_string());
        }
        if piv != col {
            mat.swap(col, piv);
            rhs.swap(col, piv);
        }
        let diag = mat
            .get(col)
            .and_then(|r| r.get(col))
            .copied()
            .unwrap_or(1.0);
        for row in (col + 1)..n {
            let factor = mat
                .get(row)
                .and_then(|r| r.get(col))
                .copied()
                .unwrap_or(0.0)
                / diag;
            for k in col..n {
                let p_val = mat.get(col).and_then(|r| r.get(k)).copied().unwrap_or(0.0);
                if let Some(cell) = mat.get_mut(row).and_then(|r| r.get_mut(k)) {
                    *cell -= factor * p_val;
                }
            }
            let rhs_col = *rhs.get(col).unwrap_or(&0.0);
            if let Some(r_slot) = rhs.get_mut(row) {
                *r_slot -= factor * rhs_col;
            }
        }
    }
    let mut x = vec![0.0_f64; n];
    for i in (0..n).rev() {
        let mut sum = *rhs.get(i).unwrap_or(&0.0);
        for j in (i + 1)..n {
            let a_ij = mat.get(i).and_then(|r| r.get(j)).copied().unwrap_or(0.0);
            let x_j = *x.get(j).unwrap_or(&0.0);
            sum -= a_ij * x_j;
        }
        let diag = mat.get(i).and_then(|r| r.get(i)).copied().unwrap_or(1.0);
        if let Some(slot) = x.get_mut(i) {
            *slot = sum / diag;
        }
    }
    Ok(x)
}

fn solve_complex_linear(
    a: &[Vec<Complex<f64>>],
    b: &[Complex<f64>],
) -> Result<Vec<Complex<f64>>, String> {
    let n = b.len();
    let mut mat = a.to_vec();
    let mut rhs = b.to_vec();
    for col in 0..n {
        let mut piv = col;
        let mut max_v = mat
            .get(col)
            .and_then(|r| r.get(col))
            .map_or(0.0, |c| c.norm());
        for row in (col + 1)..n {
            let v = mat
                .get(row)
                .and_then(|r| r.get(col))
                .map_or(0.0, |c| c.norm());
            if v > max_v {
                max_v = v;
                piv = row;
            }
        }
        if max_v < 1e-14 {
            return Err("Singular complex linear system".to_string());
        }
        if piv != col {
            mat.swap(col, piv);
            rhs.swap(col, piv);
        }
        let diag = mat
            .get(col)
            .and_then(|r| r.get(col))
            .copied()
            .unwrap_or_else(|| Complex::new(1.0, 0.0));
        for row in (col + 1)..n {
            let factor = mat
                .get(row)
                .and_then(|r| r.get(col))
                .copied()
                .unwrap_or_default()
                / diag;
            for k in col..n {
                let p_val = mat
                    .get(col)
                    .and_then(|r| r.get(k))
                    .copied()
                    .unwrap_or_default();
                if let Some(cell) = mat.get_mut(row).and_then(|r| r.get_mut(k)) {
                    *cell -= factor * p_val;
                }
            }
            let rhs_col = rhs.get(col).copied().unwrap_or_default();
            if let Some(r_slot) = rhs.get_mut(row) {
                *r_slot -= factor * rhs_col;
            }
        }
    }
    let mut x = vec![Complex::new(0.0, 0.0); n];
    for i in (0..n).rev() {
        let mut sum = rhs.get(i).copied().unwrap_or_default();
        for j in (i + 1)..n {
            let a_ij = mat
                .get(i)
                .and_then(|r| r.get(j))
                .copied()
                .unwrap_or_default();
            let x_j = x.get(j).copied().unwrap_or_default();
            sum -= a_ij * x_j;
        }
        let diag = mat
            .get(i)
            .and_then(|r| r.get(i))
            .copied()
            .unwrap_or_else(|| Complex::new(1.0, 0.0));
        if let Some(slot) = x.get_mut(i) {
            *slot = sum / diag;
        }
    }
    Ok(x)
}

/// Reusable Kerzman–Stein Szegő kernel and Riemann map solver for a smooth
/// counterclockwise Jordan domain $\Omega$ parameterized by `gamma(t)` for $t \in [0, 2\pi]$.
pub struct SzegoKerzmanSteinSolver {
    pts: Vec<Complex<f64>>,
    tangents: Vec<Complex<f64>>,
    speeds: Vec<f64>,
    dt: f64,
    a: Complex<f64>,
    s_boundary: Vec<Complex<f64>>,
    s_aa: Complex<f64>,
}

impl SzegoKerzmanSteinSolver {
    /// Solves $(I + A) `S_a` = \bar{H}_a$ once and caches the boundary Szegő kernel.
    ///
    /// # Errors
    ///
    /// Returns an error for fewer than 16 quadrature nodes, a vanishing boundary tangent, or a failed linear solve.
    pub fn solve<F, DF>(
        gamma: F,
        dgamma: DF,
        n_quad: usize,
        a: Complex<f64>,
    ) -> Result<Self, String>
    where
        F: Fn(f64) -> Complex<f64>,
        DF: Fn(f64) -> Complex<f64>,
    {
        if n_quad < 16 {
            return Err("n_quad must be at least 16".to_string());
        }
        let dt = 2.0 * PI / (n_quad as f64);
        let mut pts = Vec::with_capacity(n_quad);
        let mut tangents = Vec::with_capacity(n_quad);
        let mut speeds = Vec::with_capacity(n_quad);
        for j in 0..n_quad {
            let t = (j as f64) * dt;
            let p = gamma(t);
            let dp = dgamma(t);
            let sp = dp.norm();
            if sp < 1e-14 {
                return Err("Zero tangent speed on boundary".to_string());
            }
            pts.push(p);
            tangents.push(dp / sp);
            speeds.push(sp);
        }

        let two_pi_i = Complex::new(0.0, 2.0 * PI);
        let mut mat = vec![vec![Complex::new(0.0, 0.0); n_quad]; n_quad];
        let mut rhs = Vec::with_capacity(n_quad);

        for i in 0..n_quad {
            let zi = *pts.get(i).ok_or_else(|| "Index".to_string())?;
            let ti = *tangents.get(i).ok_or_else(|| "Index".to_string())?;
            let h_a_i = (ti / (zi - a)) / two_pi_i;
            rhs.push(h_a_i.conj());
            for j in 0..n_quad {
                let ds_j = *speeds.get(j).ok_or_else(|| "Index".to_string())? * dt;
                if i == j {
                    if let Some(cell) = mat.get_mut(i).and_then(|r| r.get_mut(j)) {
                        *cell = Complex::new(1.0, 0.0);
                    }
                } else {
                    let zj = *pts.get(j).ok_or_else(|| "Index".to_string())?;
                    let tj = *tangents.get(j).ok_or_else(|| "Index".to_string())?;
                    let h_ij = (tj / (zj - zi)) / two_pi_i;
                    let h_ji = (ti / (zi - zj)) / two_pi_i;
                    let k_ij = h_ji.conj() - h_ij;
                    if let Some(cell) = mat.get_mut(i).and_then(|r| r.get_mut(j)) {
                        *cell = k_ij * ds_j;
                    }
                }
            }
        }

        let s_boundary = solve_complex_linear(&mat, &rhs)?;
        let mut solver = Self {
            pts,
            tangents,
            speeds,
            dt,
            a,
            s_boundary,
            s_aa: Complex::new(0.0, 0.0),
        };
        solver.s_aa = solver.eval_szego(a);
        Ok(solver)
    }

    /// Evaluates the Szegő kernel $S(w, a)$ at an interior point `w`.
    #[must_use]
    pub fn eval_szego(&self, w: Complex<f64>) -> Complex<f64> {
        let two_pi_i = Complex::new(0.0, 2.0 * PI);
        let mut sum = Complex::new(0.0, 0.0);
        for j in 0..self.pts.len() {
            let zj = self.pts[j];
            let tj = self.tangents[j];
            let sj = self.s_boundary[j];
            let ds_j = self.speeds[j] * self.dt;
            sum += (sj * tj / (zj - w)) * ds_j;
        }
        sum / two_pi_i
    }

    /// Evaluates the Riemann map $f(z)$ and Szegő kernel $S(z, a)$ at `z`.
    #[must_use]
    pub fn eval_map(&self, z: Complex<f64>) -> (Complex<f64>, Complex<f64>) {
        let s_za = self.eval_szego(z);
        let dir = z - self.a;
        let integral_sq = integrate_complex_01(
            |u| {
                let w = self.a + dir * u;
                let sw = self.eval_szego(w);
                sw * sw * dir
            },
            1e-9,
            128,
        );
        let riemann_map_z = (Complex::new(2.0 * PI, 0.0) / self.s_aa) * integral_sq;
        (riemann_map_z, s_za)
    }

    /// Bergman kernel $K(z, a) = 4 \pi S(z, a)^2$ for simply-connected $\Omega$.
    #[must_use]
    pub fn eval_bergman(&self, z: Complex<f64>) -> Complex<f64> {
        let s = self.eval_szego(z);
        4.0 * PI * s * s
    }
}

/// Computes a Riemann map of a smooth Jordan domain.
///
/// Computes the Szegő kernel `S(\cdot, a)` on the boundary of a smooth Jordan domain
/// parameterized CCW by `gamma(t)` with derivative `dgamma(t)` for `t \in [0, 2\pi]`
/// via the Kerzman–Stein integral equation `(I + A) S_a = \bar{H}_a`, and evaluates
/// the Riemann map `f(z) = \frac{2\pi}{S(a,a)} \int_a^z S(\zeta, a)^2 d\zeta` at `z`.
///
/// # Errors
///
/// Returns an error if the boundary data or integral-equation solve is invalid.
pub fn szego_kerzman_stein_map<F, DF>(
    gamma: F,
    dgamma: DF,
    n_quad: usize,
    a: Complex<f64>,
    z: Complex<f64>,
) -> Result<(Complex<f64>, Complex<f64>), String>
where
    F: Fn(f64) -> Complex<f64>,
    DF: Fn(f64) -> Complex<f64>,
{
    let solver = SzegoKerzmanSteinSolver::solve(gamma, dgamma, n_quad, a)?;
    Ok(solver.eval_map(z))
}

/// Exact Bergman kernel of the unit disk: $K_{\mathbb{D}}(z, w) = \frac{1}{\pi (1 - z \bar{w})^2}$.
#[must_use]
pub fn bergman_kernel_disk(z: Complex<f64>, w: Complex<f64>) -> Complex<f64> {
    let denom = Complex::new(1.0, 0.0) - z * w.conj();
    Complex::new(1.0 / PI, 0.0) / (denom * denom)
}

/// Evaluates a Bergman kernel expansion `K(z, w) = \sum_{k=0}^M \phi_k(z) \overline{\phi_k(w)}`
/// from a list of orthonormal polynomial basis vectors (each in ascending power coefficients).
#[must_use]
pub fn bergman_kernel_polynomial_basis(
    orthonormal_polys: &[Vec<Complex<f64>>],
    z: Complex<f64>,
    w: Complex<f64>,
) -> Complex<f64> {
    let eval_p = |coeffs: &[Complex<f64>], x: Complex<f64>| -> Complex<f64> {
        let mut acc = Complex::new(0.0, 0.0);
        for &c in coeffs.iter().rev() {
            acc = acc * x + c;
        }
        acc
    };
    let mut sum = Complex::new(0.0, 0.0);
    for poly in orthonormal_polys {
        sum += eval_p(poly, z) * eval_p(poly, w).conj();
    }
    sum
}
