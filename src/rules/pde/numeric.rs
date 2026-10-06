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

use crate::graph::Arity;
use crate::graph::rule::Installer;
use crate::graph::OpDescriptor;
use crate::graph::RuleError;
use crate::kernels::special::bessel_j;

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
                    let mid = 0.5 * (lo + hi);
                    let gm = g(mid);
                    if glo * gm <= 0.0 {
                        hi = mid;
                    } else {
                        lo = mid;
                        glo = gm;
                    }
                }
                return 0.5 * (lo + hi);
            }
        }
        x = next;
        gx = gn;
    }
    f64::NAN
}

fn index(m: f64) -> Option<u32> {
    (m >= 0.0 && m.fract() == 0.0 && m <= 1e6).then(|| m as u32)
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

/// Registers the operators.
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    i.op(OpDescriptor::new("bessel_zero", Arity::Fixed(2)).eval(bessel_zero))?;
    i.op(OpDescriptor::new("bessel_root", Arity::Fixed(3)).eval(bessel_root))?;
    i.op(OpDescriptor::new("sl_root", Arity::Fixed(6)).eval(sl_root))?;
    Ok(())
}
