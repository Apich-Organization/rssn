//! Numeric semantics of the eigenvalue operators used by the expansion
//! solvers: zeros of Bessel functions and roots of Robin eigenvalue
//! equations.
//!
//! * `bessel_zero(nu, m)`: the `m`-th positive zero of `J_ν`.
//! * `bessel_root(nu, h, m)`: the `m`-th positive root `z` of
//!   `z J_ν'(z) + h J_ν(z) = 0` (`m = 0` gives the trivial root 0). With
//!   `h = 0` these are the zeros of `J_ν'`; with `h = R·h_robin` they are
//!   the Robin eigenvalues of a disk of radius `R`.
//! * `sl_root(p0, q0, p1, q1, L, m)`: the `m`-th positive `k` for which
//!   `X'' = -k² X` has a nontrivial solution with `p0 X(0) + q0 X'(0) = 0`
//!   and `p1 X(L) + q1 X'(L) = 0`.
//! * `sl_neg_root(p0, q0, p1, q1, L)`: the `κ > 0` for which `X'' = κ² X`
//!   has such a solution (a negative eigenvalue `-κ²` of `-X''`, which
//!   occurs for Robin conditions of the wrong sign); NaN if there is none.
//! * `annulus_root(nu, pa, qa, a, pb, qb, b, m)`: the `m`-th positive `k`
//!   for which a cylinder function `F = A J_ν(k r) - B Y_ν(k r)` of order
//!   `ν` satisfies `pa F(a) + qa F'(a) = 0` and `pb F(b) + qb F'(b) = 0`
//!   (the eigenvalue condition of an annulus `a < r < b`: the zeros of the
//!   cross-product `J_ν(k a) Y_ν(k b) - J_ν(k b) Y_ν(k a)` for Dirichlet
//!   conditions).
//! * `sl_neg_root_n(p0, q0, p1, q1, L, m)`: the `m`-th `κ > 0` (ascending)
//!   for which `X'' = κ² X` has such a solution; Robin conditions of the
//!   wrong sign at both ends give two negative eigenvalues.
//! * `bessel_neg_root(nu, h)`: the `κ > 0` with `z I_ν'(z) + h I_ν(z) = 0`
//!   at `z = κ` (the radial mode `I_ν(κ r/R)`, eigenvalue `-κ²/R²`, of a
//!   disk or ball with a Robin condition of the wrong sign); NaN unless
//!   `h < -ν`.
//! * `annulus_neg_root(nu, pa, qa, a, pb, qb, b, m)`: the `m`-th `κ > 0`
//!   for which `F = A I_ν(κ r) + B K_ν(κ r)` meets both Robin conditions of
//!   an annulus (negative eigenvalue `-κ²`).
//! * `struve_h(nu, x)`, `struve_l(nu, x)`: the Struve functions `H_ν`,
//!   `L_ν` (needed by the Fourier–Bessel moments of odd powers).

use crate::graph::Arity;
use crate::graph::rule::Installer;
use crate::graph::OpDescriptor;
use crate::graph::RuleError;
use crate::kernels::special::bessel_i;
use crate::kernels::special::bessel_i0;
use crate::kernels::special::bessel_i1;
use crate::kernels::special::bessel_k0;
use crate::kernels::special::bessel_k1;
use crate::kernels::special::bessel_j;
use crate::kernels::special::bessel_k;
use crate::kernels::special::gamma_numerical;
use crate::kernels::special::bessel_y;
use crate::kernels::special::bessel_y1;

const BISECTIONS: u32 = 80;

thread_local! {
    static ROOTS: std::cell::RefCell<std::collections::HashMap<Vec<u64>, f64>> = std::cell::RefCell::new(std::collections::HashMap::new());
}

/// Memoises the root-finding operators (the same eigenvalue is requested
/// again at every quadrature node).
fn cached(
    tag: u64,
    args: &[f64],
    compute: impl FnOnce() -> f64,
) -> f64 {
    let mut key: Vec<u64> = args.iter().map(|v| v.to_bits()).collect();
    key.push(tag);
    if let Some(found) = ROOTS.with(|c| c.borrow().get(&key).copied()) {
        return found;
    }
    let value = compute();
    ROOTS.with(|c| {
        let mut map = c.borrow_mut();
        if map.len() > 8192 {
            map.clear();
        }
        map.insert(key, value);
    });
    value
}

/// The `m`-th sign change of `g` on `x ≥ start` (steps of `step`), refined
/// by bisection.
fn nth_sign_change(
    g: &dyn Fn(f64) -> f64,
    start: f64,
    step: f64,
    m: u32,
) -> f64 {
    let mut found = 0;
    let mut x = start;
    let mut gx = g(x);
    for _ in 0..200_000_u32 {
        let next = x + step;
        let gn = g(next);
        if gx.is_finite() && gn.is_finite() && (gx == 0.0 || gx * gn < 0.0) {
            found += 1;
            if found == m {
                let (mut lo, mut hi, mut glo) = (x, next, gx);
                for _ in 0..BISECTIONS {
                    let mid = f64::midpoint(lo, hi);
                    let gm = g(mid);
                    if glo * gm <= 0.0 {
                        hi = mid;
                    } else {
                        lo = mid;
                        glo = gm;
                    }
                }
                return f64::midpoint(lo, hi);
            }
        }
        x = next;
        gx = gn;
    }
    f64::NAN
}

#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss, clippy::float_cmp)]
fn index(m: f64) -> Option<u32> {
    if (0.0..=1e6).contains(&m) && m.fract() == 0.0 { Some(m as u32) } else { None }
}

/// `bessel_zero(nu, m)`.
pub(super) fn bessel_zero(args: &[f64]) -> f64 {
    let (Some(&nu), Some(m)) = (args.first(), args.get(1).and_then(|&m| index(m))) else {
        return f64::NAN;
    };
    if nu.is_nan() || nu < 0.0 || m == 0 {
        return f64::NAN;
    }
    let g = |x: f64| bessel_j(nu, x);
    cached(3, args, || nth_sign_change(&g, 0.5_f64.mul_add(nu, 0.01), 0.1, m))
}

/// `bessel_root(nu, h, m)`.
pub(super) fn bessel_root(args: &[f64]) -> f64 {
    let (Some(&nu), Some(&h), Some(m)) = (args.first(), args.get(1), args.get(2).and_then(|&m| index(m))) else {
        return f64::NAN;
    };
    if m == 0 {
        return 0.0;
    }
    let g = |x: f64| 0.5 * x * (bessel_j(nu - 1.0, x) - bessel_j(nu + 1.0, x)) + h * bessel_j(nu, x);
    cached(4, args, || nth_sign_change(&g, 0.05, 0.05, m))
}

/// `sl_root(p0, q0, p1, q1, L, m)`.
pub(super) fn sl_root(args: &[f64]) -> f64 {
    let (Some(&p0), Some(&q0), Some(&p1), Some(&q1), Some(&l), Some(m)) =
        (args.first(), args.get(1), args.get(2), args.get(3), args.get(4), args.get(5).and_then(|&m| index(m)))
    else {
        return f64::NAN;
    };
    if m == 0 || l <= 0.0 {
        return f64::NAN;
    }
    let g = |k: f64| (p1 * q0 - q1 * p0) * k * (k * l).cos() - (p1 * p0 + q1 * q0 * k * k) * (k * l).sin();
    let scale = std::f64::consts::PI / l;
    cached(5, args, || nth_sign_change(&g, 1e-3 * scale, scale / 40.0, m))
}

/// `(p0 q0 p1 q1 L)` scanned for the root `κ` of
/// `(p1 q0 - q1 p0) κ cosh κL - (p1 p0 - q1 q0 κ²) sinh κL`, scaled by
/// `2 e^{-κL}` to avoid overflow.
fn negative_condition(
    p0: f64,
    q0: f64,
    p1: f64,
    q1: f64,
    l: f64,
    kappa: f64,
) -> f64 {
    let decay = (-2.0 * kappa * l).exp();
    let a = p1.mul_add(q0, -q1 * p0);
    let b = p1.mul_add(p0, -(q1 * q0 * kappa * kappa));
    (a * kappa).mul_add(1.0 + decay, -b * (1.0 - decay))
}

/// `sl_neg_root(p0, q0, p1, q1, L)`.
pub(super) fn sl_neg_root(args: &[f64]) -> f64 {
    let (Some(&p0), Some(&q0), Some(&p1), Some(&q1), Some(&l)) = (args.first(), args.get(1), args.get(2), args.get(3), args.get(4)) else {
        return f64::NAN;
    };
    if l <= 0.0 {
        return f64::NAN;
    }
    let ratio = |p: f64, q: f64| if q == 0.0 { 0.0 } else { (p / q).abs() };
    let top = 8.0 / l + 4.0 * ratio(p0, q0).max(ratio(p1, q1));
    let g = |k: f64| negative_condition(p0, q0, p1, q1, l, k);
    let step = top / 4000.0;
    let found = nth_sign_change(&g, step * 0.05, step, 1);
    if found.is_finite() && found < top { found } else { f64::NAN }
}

/// The `m`-th root `κ > 0` of the negative-eigenvalue condition.
fn negative_root_index(
    p0: f64,
    q0: f64,
    p1: f64,
    q1: f64,
    l: f64,
    m: u32,
) -> f64 {
    if l <= 0.0 || m == 0 {
        return f64::NAN;
    }
    let ratio = |p: f64, q: f64| if q == 0.0 { 0.0 } else { (p / q).abs() };
    let top = 8.0 / l + 4.0 * ratio(p0, q0).max(ratio(p1, q1));
    let g = |k: f64| negative_condition(p0, q0, p1, q1, l, k);
    let step = top / 8000.0;
    let found = nth_sign_change(&g, step * 0.05, step, m);
    if found.is_finite() && found < top { found } else { f64::NAN }
}

/// The number of negative eigenvalues of the Robin problem (0, 1 or 2).
pub(super) fn negative_count(
    p0: f64,
    q0: f64,
    p1: f64,
    q1: f64,
    l: f64,
) -> u32 {
    let count = cached(8, &[p0, q0, p1, q1, l], || {
        let mut count = 0;
        while count < 4 && negative_root_index(p0, q0, p1, q1, l, count + 1).is_finite() {
            count += 1;
        }
        f64::from(count)
    });
    count as u32
}

/// `sl_neg_root_n(p0, q0, p1, q1, L, m)`.
pub(super) fn sl_neg_root_n(args: &[f64]) -> f64 {
    let (Some(&p0), Some(&q0), Some(&p1), Some(&q1), Some(&l), Some(m)) =
        (args.first(), args.get(1), args.get(2), args.get(3), args.get(4), args.get(5).and_then(|&m| index(m)))
    else {
        return f64::NAN;
    };
    negative_root_index(p0, q0, p1, q1, l, m)
}

/// `bessel_neg_root(nu, h)`.
pub(super) fn bessel_neg_root(args: &[f64]) -> f64 {
    let (Some(&nu), Some(&h)) = (args.first(), args.get(1)) else {
        return f64::NAN;
    };
    if nu.is_nan() || h.is_nan() || nu < 0.0 || nu + h >= 0.0 {
        return f64::NAN;
    }
    // z I_{ν+1}(z)/I_ν(z) + ν + h increases from ν + h < 0 to infinity.
    let g = |z: f64| z * bessel_i(nu + 1.0, z) / bessel_i(nu, z) + nu + h;
    let (mut lo, mut hi) = (1e-9, 2.0 * (-h) + 10.0);
    for _ in 0..200 {
        let mid = f64::midpoint(lo, hi);
        if g(mid) < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    f64::midpoint(lo, hi)
}

/// The power series of `H_ν` (`sign = -1`) or `L_ν` (`sign = 1`).
fn struve_series(
    nu: f64,
    x: f64,
    sign: f64,
) -> f64 {
    let half = 0.5 * x;
    let mut term = half.powf(nu + 1.0) / (gamma_numerical(1.5) * gamma_numerical(nu + 1.5));
    let mut total = term;
    for m in 0..2000_u32 {
        let k = f64::from(m);
        term *= sign * half * half / ((k + 1.5) * (k + nu + 1.5));
        total += term;
        if term.abs() <= 1e-17 * total.abs() {
            break;
        }
    }
    total
}

/// `H_ν(x)` by the integral representation, for `2ν` a non-negative integer.
fn struve_h_integral(
    nu: f64,
    x: f64,
) -> f64 {
    let steps = 2000 + (200.0 * x.abs()) as usize / 2 * 2;
    let width = std::f64::consts::FRAC_PI_2 / steps as f64;
    let order = (2.0 * nu).round() as i32;
    let f = |theta: f64| theta.cos().powi(order) * (x * theta.sin()).sin();
    let mut total = f(0.0) + f(std::f64::consts::FRAC_PI_2);
    for k in 1..steps {
        total += f(k as f64 * width) * if k % 2 == 1 { 4.0 } else { 2.0 };
    }
    let integral = total * width / 3.0;
    2.0 * (0.5 * x).powf(nu) / (std::f64::consts::PI.sqrt() * gamma_numerical(nu + 0.5)) * integral
}

/// `struve_h(nu, x)`.
pub(super) fn struve_h(args: &[f64]) -> f64 {
    let (Some(&nu), Some(&x)) = (args.first(), args.get(1)) else {
        return f64::NAN;
    };
    if nu.is_nan() || x.is_nan() {
        return f64::NAN;
    }
    if x < 0.0 {
        // H_ν(-x) = (-1)^{ν+1} H_ν(x) for integer ν.
        if nu.fract() != 0.0 {
            return f64::NAN;
        }
        let sign = if nu.rem_euclid(2.0) == 0.0 { -1.0 } else { 1.0 };
        return sign * struve_h(&[nu, -x]);
    }
    if x < 8.0 {
        return struve_series(nu, x, -1.0);
    }
    let twice = 2.0 * nu;
    if nu > -0.5 && (twice - twice.round()).abs() < 1e-12 {
        struve_h_integral(nu, x)
    } else {
        f64::NAN
    }
}

/// `struve_l(nu, x)`.
pub(super) fn struve_l(args: &[f64]) -> f64 {
    let (Some(&nu), Some(&x)) = (args.first(), args.get(1)) else {
        return f64::NAN;
    };
    if nu.is_nan() || x.is_nan() {
        return f64::NAN;
    }
    if x < 0.0 {
        if nu.fract() != 0.0 {
            return f64::NAN;
        }
        let sign = if nu.rem_euclid(2.0) == 0.0 { -1.0 } else { 1.0 };
        return sign * struve_l(&[nu, -x]);
    }
    struve_series(nu, x, 1.0)
}

/// `p F(r) + q F'(r)` for `F = J_ν(k r)` (`kind` false) or `Y_ν(k r)`.
fn cylinder_condition(
    nu: f64,
    p: f64,
    q: f64,
    k: f64,
    r: f64,
    second_kind: bool,
) -> f64 {
    let z = k * r;
    #[allow(clippy::float_cmp)]
    let f = |order: f64| {
        if !second_kind {
            bessel_j(order, z)
        } else if order == 1.0 {
            bessel_y1(z)
        } else {
            bessel_y(order, z)
        }
    };
    let value = f(nu);
    let slope = k * (-f(nu + 1.0) + nu / z * value);
    p.mul_add(value, q * slope)
}

/// `I_ν(x)` with the fast integer cases.
#[allow(clippy::float_cmp)]
fn fast_i(
    nu: f64,
    x: f64,
) -> f64 {
    if nu == 0.0 {
        bessel_i0(x)
    } else if nu == 1.0 {
        bessel_i1(x)
    } else {
        bessel_i(nu, x)
    }
}

/// `K_ν(x)` with the fast integer cases.
#[allow(clippy::float_cmp)]
fn fast_k(
    nu: f64,
    x: f64,
) -> f64 {
    if nu == 0.0 {
        bessel_k0(x)
    } else if nu == 1.0 {
        bessel_k1(x)
    } else {
        bessel_k(nu, x)
    }
}

/// `p F(r) + q F'(r)` for `F = I_ν(k r)` or `K_ν(k r)` (`second`).
fn modified_condition(
    nu: f64,
    p: f64,
    q: f64,
    k: f64,
    r: f64,
    second: bool,
) -> f64 {
    let z = k * r;
    let (value, slope) = if second {
        let v = fast_k(nu, z);
        (v, k * (-fast_k(nu + 1.0, z) + nu / z * v))
    } else {
        let v = fast_i(nu, z);
        (v, k * (fast_i(nu + 1.0, z) + nu / z * v))
    };
    p.mul_add(value, q * slope)
}

/// The `m`-th negative-eigenvalue root of the annulus problem.
fn annulus_negative_index(
    v: [f64; 7],
    m: u32,
) -> f64 {
    let [nu, pa, qa, a, pb, qb, b] = v;
    if m == 0 || a <= 0.0 || b <= a || nu < 0.0 {
        return f64::NAN;
    }
    let ratio = |p: f64, q: f64| if q == 0.0 { 0.0 } else { (p / q).abs() };
    let top = (8.0 / (b - a) + 4.0 * ratio(pa, qa).max(ratio(pb, qb)) + 2.0 * nu / a).min(600.0 / b);
    let g = |k: f64| {
        let (ia, ka) = (modified_condition(nu, pa, qa, k, a, false), modified_condition(nu, pa, qa, k, a, true));
        let (ib, kb) = (modified_condition(nu, pb, qb, k, b, false), modified_condition(nu, pb, qb, k, b, true));
        ka * ib - ia * kb
    };
    let step = top / 8000.0;
    let found = nth_sign_change(&g, step * 0.05, step, m);
    if found.is_finite() && found < top { found } else { f64::NAN }
}

/// The number of negative eigenvalues of the annulus problem.
pub(super) fn annulus_neg_count(v: [f64; 7]) -> u32 {
    let count = cached(7, &v, || {
        let mut count = 0;
        while count < 4 && annulus_negative_index(v, count + 1).is_finite() {
            count += 1;
        }
        f64::from(count)
    });
    count as u32
}

/// Whether the annulus problem has the eigenvalue zero (a mode
/// `A r^ν + B r^{-ν}`, or `A + B ln r` for `ν = 0`), decided by the
/// determinant of the two boundary conditions on that basis.
#[allow(clippy::float_cmp)]
pub(super) fn annulus_has_zero(v: [f64; 7]) -> bool {
    let [nu, pa, qa, a, pb, qb, b] = v;
    if a <= 0.0 || b <= a || nu < 0.0 {
        return false;
    }
    // The basis functions and their derivatives at r.
    let basis = |r: f64| -> [(f64, f64); 2] {
        if nu == 0.0 {
            [(1.0, 0.0), (r.ln(), 1.0 / r)]
        } else {
            [(r.powf(nu), nu * r.powf(nu - 1.0)), (r.powf(-nu), -nu * r.powf(-nu - 1.0))]
        }
    };
    let (fa, fb) = (basis(a), basis(b));
    let condition = |p: f64, q: f64, f: (f64, f64)| p * f.0 + q * f.1;
    let m = [[condition(pa, qa, fa[0]), condition(pa, qa, fa[1])], [condition(pb, qb, fb[0]), condition(pb, qb, fb[1])]];
    let determinant = m[0][0] * m[1][1] - m[0][1] * m[1][0];
    let scale = (m[0][0] * m[1][1]).abs() + (m[0][1] * m[1][0]).abs();
    scale > 0.0 && determinant.abs() <= 1e-9 * scale
}

/// `annulus_neg_root(nu, pa, qa, a, pb, qb, b, m)`.
pub(super) fn annulus_neg_root(args: &[f64]) -> f64 {
    let (Some(&nu), Some(&pa), Some(&qa), Some(&a), Some(&pb), Some(&qb), Some(&b), Some(m)) =
        (args.first(), args.get(1), args.get(2), args.get(3), args.get(4), args.get(5), args.get(6), args.get(7).and_then(|&m| index(m)))
    else {
        return f64::NAN;
    };
    cached(1, args, || annulus_negative_index([nu, pa, qa, a, pb, qb, b], m))
}

/// `annulus_root(nu, pa, qa, a, pb, qb, b, m)`.
pub(super) fn annulus_root(args: &[f64]) -> f64 {
    let (Some(&nu), Some(&pa), Some(&qa), Some(&a), Some(&pb), Some(&qb), Some(&b), Some(m)) =
        (args.first(), args.get(1), args.get(2), args.get(3), args.get(4), args.get(5), args.get(6), args.get(7).and_then(|&m| index(m)))
    else {
        return f64::NAN;
    };
    if m == 0 || a <= 0.0 || b <= a || nu < 0.0 {
        return f64::NAN;
    }
    let g = |k: f64| {
        let ya = cylinder_condition(nu, pa, qa, k, a, true);
        let ja = cylinder_condition(nu, pa, qa, k, a, false);
        let jb = cylinder_condition(nu, pb, qb, k, b, false);
        let yb = cylinder_condition(nu, pb, qb, k, b, true);
        ya.mul_add(jb, -(ja * yb))
    };
    let scale = std::f64::consts::PI / (b - a);
    cached(2, args, || nth_sign_change(&g, 1e-3 * scale, scale / 40.0, m))
}

/// Registers the operators.
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    i.op(OpDescriptor::new("bessel_zero", Arity::Fixed(2)).eval(bessel_zero))?;
    i.op(OpDescriptor::new("bessel_root", Arity::Fixed(3)).eval(bessel_root))?;
    i.op(OpDescriptor::new("sl_root", Arity::Fixed(6)).eval(sl_root))?;
    i.op(OpDescriptor::new("sl_neg_root", Arity::Fixed(5)).eval(sl_neg_root))?;
    i.op(OpDescriptor::new("sl_neg_root_n", Arity::Fixed(6)).eval(sl_neg_root_n))?;
    i.op(OpDescriptor::new("bessel_neg_root", Arity::Fixed(2)).eval(bessel_neg_root))?;
    i.op(OpDescriptor::new("struve_h", Arity::Fixed(2)).eval(struve_h))?;
    i.op(OpDescriptor::new("struve_l", Arity::Fixed(2)).eval(struve_l))?;
    i.op(OpDescriptor::new("annulus_neg_root", Arity::Fixed(8)).eval(annulus_neg_root))?;
    i.op(OpDescriptor::new("annulus_root", Arity::Fixed(8)).eval(annulus_root))?;
    Ok(())
}
