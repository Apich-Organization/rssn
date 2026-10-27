//! Numeric inverse scattering for the `KdV` equation from general initial data.
//!
//! `pde_kdv_scattering(u0(x), x, t)` evaluates the solution of
//! `u_t - 6 u u_x + u_xxx = 0`, `u(x, 0) = u0(x)`, at the point `(x, t)` by the
//! Gelfand–Levitan–Marchenko route:
//!
//! 1. the Jost solution `f(x, k) = e^{ikx}` (`x → +∞`) of `-ψ'' + u0 ψ = k² ψ`
//!    is integrated to `x = -L` by RK4, giving `f = A e^{ikx} + B e^{-ikx}`
//!    and the reflection coefficient `R(k) = -conj(B)/A`;
//! 2. the bound states `κ_n` are the roots of the growing coefficient of
//!    `f(x, iκ)`, with norming constants `c_n² = 1/∫ f(x, iκ_n)² dx`;
//! 3. the Marchenko kernel `F(s, t)` is the sum of `c_n² e^{8 κ_n³ t - κ_n s}`
//!    and `(1/π) ∫_0^∞ Re[R(k) e^{i(ks + 8k³t)}] dk`; it enters the integral
//!    equation `K(x, y) + F(x + y) + ∫_x^∞ K(x, z) F(z + y) dz = 0`, solved
//!    by Gauss–Legendre Nyström collocation, and `u = -2 d/dx K(x, x)`.
//!
//! The potential must decay (it is sampled on `[-12, 12]`); moderate
//! times `|t| ≲ 0.5` keep the dispersive phase `8 k³ t` resolved.

use super::curvilinear::gauss_legendre;
use crate::graph::Cx;
use crate::graph::Env;
use crate::graph::NodeId;
use crate::graph::Number;
use num_complex::Complex64;

/// Half-length of the scattering interval.
const HALF: f64 = 12.0;
/// RK4 step.
const STEP: f64 = 0.01;
/// Spectral grid step and cutoff.
const DK: f64 = 0.02;
const KMAX: f64 = 10.0;

/// RK4 for `y'' = q(x) y` backwards from `x = HALF` to `-HALF`, calling
/// `visit(x, y)` after every step.
fn integrate(
    potential: &[f64],
    shift: Complex64,
    start: (Complex64, Complex64),
    mut visit: impl FnMut(f64, Complex64),
) -> (Complex64, Complex64) {
    // `potential` samples V at x = HALF - j STEP / 2.
    let steps = potential.len() / 2;
    let (mut y, mut p) = start;
    let deriv = |v: f64, y: Complex64, p: Complex64| -> (Complex64, Complex64) { (p, (Complex64::new(v, 0.0) - shift) * y) };
    for j in 0..steps {
        let (v0, v1, v2) = (potential[2 * j], potential[2 * j + 1], potential[2 * j + 2]);
        let h = -STEP;
        let (k1y, k1p) = deriv(v0, y, p);
        let (k2y, k2p) = deriv(v1, y + k1y * (h / 2.0), p + k1p * (h / 2.0));
        let (k3y, k3p) = deriv(v1, y + k2y * (h / 2.0), p + k2p * (h / 2.0));
        let (k4y, k4p) = deriv(v2, y + k3y * h, p + k3p * h);
        y += (k1y + k2y * 2.0 + k3y * 2.0 + k4y) * (h / 6.0);
        p += (k1p + k2p * 2.0 + k3p * 2.0 + k4p) * (h / 6.0);
        visit(HALF - STEP * (j + 1) as f64, y);
    }
    (y, p)
}

/// The data for `GLM`.
struct Data {
    bound: Vec<(f64, f64)>,
    /// `(k, weight, R(k))` of a midpoint rule, finer near `k = 0` where weakly bound states shape `R`.
    reflection: Vec<(f64, f64, Complex64)>,
}

fn scattering_data(potential: &[f64]) -> Data {
    let kappa_max = potential.iter().fold(0.0_f64, |a, &v| a.max(-v)).sqrt() + 0.3;
    // Bound states: roots of g(κ) = κ f - f' at x = -L.
    let growing = |kappa: f64| -> f64 {
        let k = Complex64::new(0.0, kappa);
        let shift = k * k;
        let (f, p) = integrate(potential, shift, (Complex64::new(1.0, 0.0), Complex64::new(-kappa, 0.0)), |_, _| {});
        kappa * f.re - p.re
    };
    let mut bound = Vec::new();
    let steps = (kappa_max / 0.01).ceil() as usize;
    let mut previous = growing(0.005);
    let mut left = 0.005;
    for j in 1..=steps {
        let right = 0.005 + 0.01 * j as f64;
        let value = growing(right);
        if previous.signum() != value.signum() && previous != 0.0 {
            let (mut a, mut b, mut fa) = (left, right, previous);
            for _ in 0..60 {
                let mid = f64::midpoint(a, b);
                let fm = growing(mid);
                if fa.signum() == fm.signum() {
                    a = mid;
                    fa = fm;
                } else {
                    b = mid;
                }
            }
            let kappa = f64::midpoint(a, b);
            // Norming constant: ∫ f² by Simpson over the stored samples.
            let k = Complex64::new(0.0, kappa);
            let mut samples = vec![1.0_f64];
            integrate(potential, k * k, (Complex64::new(1.0, 0.0), Complex64::new(-kappa, 0.0)), |_, y| samples.push(y.re));
            let n = samples.len() - 1;
            let mut total = samples[0] * samples[0] + samples[n] * samples[n];
            for (i, v) in samples.iter().enumerate().take(n).skip(1) {
                total += if i % 2 == 1 { 4.0 } else { 2.0 } * v * v;
            }
            // The Simpson sum covers [-L, L]; the exponential tails of f beyond both ends are added.
            let integral = total * STEP / 3.0 + (1.0 + samples[n] * samples[n]) / (2.0 * kappa);
            let norming = (2.0 * kappa * HALF).exp() / integral;
            bound.push((kappa, norming));
        }
        previous = value;
        left = right;
    }
    // Reflection coefficient on the half line k > 0.
    let mut reflection = Vec::new();
    let mut quiet = 0;
    let fine = (0..160).map(|j| ((f64::from(j) + 0.5) * 0.0025, 0.0025));
    let coarse = (0..((KMAX - 0.4) / DK) as usize).map(|j| (0.4 + (j as f64 + 0.5) * DK, DK));
    for (k, weight) in fine.chain(coarse) {
        let kc = Complex64::new(k, 0.0);
        let (f, p) = integrate(potential, kc * kc, (Complex64::new(0.0, k * HALF).exp(), Complex64::new(0.0, k) * Complex64::new(0.0, k * HALF).exp()), |_, _| {});
        let i_k = Complex64::new(0.0, k);
        let x = -HALF;
        let a = (f + p / i_k) * Complex64::new(0.0, -k * x).exp() * 0.5;
        let b = (f - p / i_k) * Complex64::new(0.0, k * x).exp() * 0.5;
        let r = -b.conj() / a;
        quiet = if r.norm() < 1e-13 { quiet + 1 } else { 0 };
        reflection.push((k, weight, r));
        if quiet > 20 {
            break;
        }
    }
    Data { bound, reflection }
}

/// `F(s, t)`.
fn marchenko_kernel(
    data: &Data,
    s: f64,
    t: f64,
) -> f64 {
    let mut total: f64 = data.bound.iter().map(|&(kappa, c2)| c2 * (8.0 * kappa.powi(3) * t - kappa * s).exp()).sum();
    let mut integral = 0.0;
    for (k, weight, r) in &data.reflection {
        integral += weight * (r * Complex64::new(0.0, k * s + 8.0 * k.powi(3) * t).exp()).re;
    }
    total += integral / std::f64::consts::PI;
    total
}

/// Gaussian elimination with partial pivoting.
fn solve(
    mut a: Vec<Vec<f64>>,
    mut b: Vec<f64>,
) -> Option<Vec<f64>> {
    let n = b.len();
    for col in 0..n {
        let pivot = (col..n).max_by(|&i, &j| a[i][col].abs().total_cmp(&a[j][col].abs()))?;
        if a[pivot][col].abs() < 1e-300 {
            return None;
        }
        a.swap(col, pivot);
        b.swap(col, pivot);
        for row in (col + 1)..n {
            let factor = a[row][col] / a[col][col];
            for k in col..n {
                let v = a[col][k];
                a[row][k] -= factor * v;
            }
            let v = b[col];
            b[row] -= factor * v;
        }
    }
    let mut x = vec![0.0; n];
    for row in (0..n).rev() {
        let tail: f64 = ((row + 1)..n).map(|k| a[row][k] * x[k]).sum();
        x[row] = (b[row] - tail) / a[row][row];
    }
    Some(x)
}

/// The numeric solution `u(x, t)` of `KdV` for the initial datum `u0`.
pub(super) fn kdv_value(
    u0: &dyn Fn(f64) -> Option<f64>,
    x: f64,
    t: f64,
) -> Option<f64> {
    let samples = (2.0 * HALF / STEP) as usize;
    let mut potential = Vec::with_capacity(2 * samples + 1);
    for j in 0..=(2 * samples) {
        potential.push(u0(HALF - STEP * j as f64 / 2.0)?);
    }
    if potential[0].abs() > 1e-6 || potential[2 * samples].abs() > 1e-6 {
        return None;
    }
    let data = scattering_data(&potential);
    // The kernel table on a uniform grid with cubic interpolation.
    let h = 0.01;
    let kappa_min = data.bound.iter().map(|&(k, _)| k).fold(f64::INFINITY, f64::min);
    let span = if kappa_min.is_finite() { (30.0 / kappa_min).clamp(14.0, 60.0) } else { 14.0 };
    let lo = 2.0 * (x - 0.05) - 0.2;
    let count = ((2.0 * (x + 0.05) + 2.0 * span - lo) / h) as usize + 4;
    let table: Vec<f64> = (0..count).map(|i| marchenko_kernel(&data, lo + h * i as f64, t)).collect();
    let kernel = |s: f64| -> f64 {
        let position = (s - lo) / h;
        let i = (position.floor() as usize).clamp(1, count - 3);
        let u = position - i as f64;
        let (p0, p1, p2, p3) = (table[i - 1], table[i], table[i + 1], table[i + 2]);
        // Cubic Lagrange through four points.
        let a = -u * (u - 1.0) * (u - 2.0) / 6.0;
        let b = (u + 1.0) * (u - 1.0) * (u - 2.0) / 2.0;
        let c = -(u + 1.0) * u * (u - 2.0) / 2.0;
        let d = (u + 1.0) * u * (u - 1.0) / 6.0;
        a * p0 + b * p1 + c * p2 + d * p3
    };
    let panels = (span / 2.0).ceil() as usize;
    let nodes = gauss_legendre(12);
    let kxx = |x0: f64| -> Option<f64> {
        let width = span / panels as f64;
        let mut z = Vec::new();
        let mut w = Vec::new();
        for panel in 0..panels {
            let start = x0 + width * panel as f64;
            for &(node, weight) in &nodes {
                z.push(start + 0.5 * width * (node + 1.0));
                w.push(0.5 * width * weight);
            }
        }
        let m = z.len();
        let mut a = vec![vec![0.0; m]; m];
        let mut b = vec![0.0; m];
        for i in 0..m {
            for j in 0..m {
                a[i][j] = w[j] * kernel(z[j] + z[i]);
            }
            a[i][i] += 1.0;
            b[i] = -kernel(x0 + z[i]);
        }
        let k = solve(a, b)?;
        let tail: f64 = (0..m).map(|j| w[j] * k[j] * kernel(z[j] + x0)).sum();
        Some(-kernel(2.0 * x0) - tail)
    };
    let step = 0.01;
    let (m2, m1, p1, p2) = (kxx(x - 2.0 * step)?, kxx(x - step)?, kxx(x + step)?, kxx(x + 2.0 * step)?);
    let derivative = (-p2 + 8.0 * p1 - 8.0 * m1 + m2) / (12.0 * step);
    Some(-2.0 * derivative)
}

/// `pde_kdv_scattering(u0, x, t)` as a kernel body.
pub(super) fn kdv_scattering(
    cx: &mut Cx<'_>,
    u0: NodeId,
    x: NodeId,
    t: NodeId,
) -> Option<NodeId> {
    let (x_value, t_value) = (super::modes::exact_value(cx, x)?, super::modes::exact_value(cx, t)?);
    let x_node = cx.graph.sym("x");
    let symbol = cx.graph.symbol_of(x_node)?;
    let free = cx.graph.free_symbols(cx.graph.find(u0)).to_vec();
    if free.iter().any(|&s| s != symbol) {
        return None;
    }
    let datum = |v: f64| -> Option<f64> {
        let mut env = Env::numeric(0.0);
        env.bind(symbol, v);
        cx.graph.eval(u0, &env).filter(|r| r.is_finite())
    };
    let value = kdv_value(&datum, x_value, t_value)?;
    Some(cx.graph.num(Number::from(value)))
}
