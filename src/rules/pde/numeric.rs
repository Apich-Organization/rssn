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

use crate::graph::Arity;
use crate::graph::rule::Installer;
use crate::graph::OpDescriptor;
use crate::graph::RuleError;
use crate::kernels::special::bessel_j;
use crate::kernels::special::bessel_y;
use crate::kernels::special::bessel_y1;

const BISECTIONS: u32 = 80;

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
    nth_sign_change(&g, 0.5_f64.mul_add(nu, 0.01), 0.1, m)
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
    nth_sign_change(&g, 0.05, 0.05, m)
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
    nth_sign_change(&g, 1e-3 * scale, scale / 40.0, m)
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
    nth_sign_change(&g, 1e-3 * scale, scale / 40.0, m)
}

/// Registers the operators.
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    i.op(OpDescriptor::new("bessel_zero", Arity::Fixed(2)).eval(bessel_zero))?;
    i.op(OpDescriptor::new("bessel_root", Arity::Fixed(3)).eval(bessel_root))?;
    i.op(OpDescriptor::new("sl_root", Arity::Fixed(6)).eval(sl_root))?;
    i.op(OpDescriptor::new("sl_neg_root", Arity::Fixed(5)).eval(sl_neg_root))?;
    i.op(OpDescriptor::new("annulus_root", Arity::Fixed(8)).eval(annulus_root))?;
    Ok(())
}
