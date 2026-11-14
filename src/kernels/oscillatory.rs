//! Highly oscillatory integrals by Levin collocation.
//!
//! `I = integral from a to b of f(x) exp(i w g(x)) dx` with a large frequency
//! `w` defeats every rule that samples the integrand: the number of points
//! must grow with `w`. Levin's observation is that if `p` solves the
//! *non-oscillatory* differential equation `p' + i w g' p = f`, then
//! `(p exp(i w g))' = f exp(i w g)` and
//!
//! ```text
//! I = p(b) exp(i w g(b)) - p(a) exp(i w g(a)).
//! ```
//!
//! [`levin`] finds `p` by collocation at Chebyshev points, where the slowly
//! varying particular solution is well approximated by a polynomial of low
//! degree *whatever* `w` is; the cost is one small dense solve, and the
//! accuracy improves as `w` grows. The method needs `g'` to keep one sign on
//! the interval (no stationary point); [`levin_adaptive`] reports an error
//! estimate, so a stationary point shows up as a failure to converge and not
//! as a wrong number.

#![allow(
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_precision_loss,
    clippy::many_single_char_names,
    clippy::module_name_repetitions,
    clippy::single_call_fn,
    clippy::missing_docs_in_private_items,
    clippy::needless_range_loop,
    reason = "small dense collocation matrices indexed within their own extents"
)]

use crate::kernels::dense;
use crate::kernels::dense::Mat;

/// A complex integral with an error estimate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Oscillatory {
    /// Real part.
    pub re: f64,
    /// Imaginary part.
    pub im: f64,
    /// An estimate of the absolute error (the difference to the next coarser
    /// collocation); `NAN` from [`levin`], which makes no estimate.
    pub error: f64,
}

/// The Chebyshev points `cos(j pi / (n - 1))` and their differentiation
/// matrix on `[-1, 1]` (Trefethen's construction, with the "negative sum
/// trick" for the diagonal).
fn chebyshev(n: usize) -> (Vec<f64>, Vec<Vec<f64>>) {
    let top = (n - 1) as f64;
    let x: Vec<f64> = (0..n)
        .map(|j| (std::f64::consts::PI * j as f64 / top).cos())
        .collect();
    let c = |j: usize| {
        if j == 0 || j == n - 1 {
            2.0
        } else {
            1.0
        }
    };
    let mut d = vec![vec![0.0; n]; n];
    for i in 0..n {
        let mut sum = 0.0;
        for j in 0..n {
            if i != j {
                let sign = if (i + j) % 2 == 0 {
                    1.0
                } else {
                    -1.0
                };
                d[i][j] = c(i) / c(j) * sign / (x[i] - x[j]);
                sum += d[i][j];
            }
        }
        d[i][i] = -sum;
    }
    (x, d)
}

/// `integral from a to b of f(x) exp(i w g(x)) dx` by Levin collocation at
/// `n` Chebyshev points (`n >= 2`); `dg` is the derivative of `g`.
///
/// # Errors
/// A message when `n < 2`, the data is not finite, or the collocation system
/// cannot be solved.
#[allow(clippy::float_cmp)] // `a == b` is the exact empty-interval sentinel
pub fn levin(
    f: &dyn Fn(f64) -> f64,
    g: &dyn Fn(f64) -> f64,
    dg: &dyn Fn(f64) -> f64,
    a: f64,
    b: f64,
    w: f64,
    n: usize,
) -> Result<Oscillatory, String> {
    if n < 2 {
        return Err("Levin collocation needs at least two points".into());
    }
    if !(a.is_finite() && b.is_finite() && w.is_finite()) {
        return Err("the interval and the frequency must be finite".into());
    }
    if a == b {
        return Ok(Oscillatory {
            re: 0.0,
            im: 0.0,
            error: f64::NAN,
        });
    }
    let (nodes, d) = chebyshev(n);
    let half = 0.5 * (b - a);
    let mid = f64::midpoint(a, b);
    let xs: Vec<f64> = nodes.iter().map(|t| mid + half * t).collect();
    // (D / half + i w diag(g')) (pr + i pi) = f, as a real 2n x 2n system.
    let mut system = Mat::zeros(2 * n, 2 * n);
    let mut rhs = vec![0.0; 2 * n];
    for i in 0..n {
        let slope = w * dg(xs[i]);
        for j in 0..n {
            let dij = d[i][j] / half;
            system.set(i, j, dij);
            system.set(n + i, n + j, dij);
        }
        system.set(i, n + i, system.at(i, n + i) - slope);
        system.set(n + i, i, system.at(n + i, i) + slope);
        rhs[i] = f(xs[i]);
    }
    if !system.is_finite() || rhs.iter().any(|v| !v.is_finite()) {
        return Err("the integrand or the phase is not finite on the interval".into());
    }
    // The system is close to singular when `w g'` is small (the constant is
    // then nearly in the kernel); the minimum-norm solution is the slowly
    // varying one, which is the one wanted.
    let p = match dense::solve(&system, &rhs) {
        | Ok(p) if p.iter().all(|v| v.is_finite()) => p,
        | _ => dense::lstsq_svd(&system, &rhs).map_err(|e| format!("collocation solve: {e}"))?,
    };
    // Node 0 is x = b, node n - 1 is x = a.
    let end = |index: usize, x: f64| {
        let (s, c) = (w * g(x)).sin_cos();
        let (pr, pi) = (p[index], p[n + index]);
        (pr * c - pi * s, pr * s + pi * c)
    };
    let (br, bi) = end(0, b);
    let (ar, ai) = end(n - 1, a);
    let (re, im) = (br - ar, bi - ai);
    if !(re.is_finite() && im.is_finite()) {
        return Err("the collocation solution is not finite".into());
    }
    Ok(Oscillatory {
        re,
        im,
        error: f64::NAN,
    })
}

/// [`levin`] with the number of points doubled until two successive results
/// agree to `tol` (absolute, or relative to the result when that is larger),
/// starting from 8 points and stopping at `max_points`.
///
/// # Errors
/// A message when no pair of successive collocations agrees: typically a
/// stationary point of the phase, or an integrand that is itself
/// oscillatory.
#[allow(clippy::too_many_arguments)] // flat signature mirroring `levin` plus the two stopping controls
pub fn levin_adaptive(
    f: &dyn Fn(f64) -> f64,
    g: &dyn Fn(f64) -> f64,
    dg: &dyn Fn(f64) -> f64,
    a: f64,
    b: f64,
    w: f64,
    tol: f64,
    max_points: usize,
) -> Result<Oscillatory, String> {
    let mut n = 8;
    let mut previous = levin(f, g, dg, a, b, w, n)?;
    let mut best = f64::INFINITY;
    while n < max_points {
        n = (2 * n).min(max_points.max(9));
        let next = levin(f, g, dg, a, b, w, n)?;
        let error = (next.re - previous.re).hypot(next.im - previous.im);
        let scale = next.re.hypot(next.im).max(1.0);
        best = best.min(error);
        if error <= tol * scale {
            return Ok(Oscillatory { error, ..next });
        }
        previous = next;
    }
    Err(format!(
        "Levin collocation did not converge with {max_points} points (best agreement {best:e}); \
         the phase may have a stationary point on the interval"
    ))
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

    use super::*;

    #[test]
    fn the_differentiation_matrix_is_exact_on_polynomials() {
        let (x, d) = chebyshev(6);
        // d/dx of x^3 at the nodes.
        for i in 0..6 {
            let got: f64 = (0..6).map(|j| d[i][j] * x[j].powi(3)).sum();
            assert!((got - 3.0 * x[i] * x[i]).abs() < 1e-12);
        }
    }

    #[test]
    fn linear_phase_matches_the_closed_forms() {
        for w in [10.0, 1.0e3, 1.0e6] {
            // integral of exp(i w x) over [0, 1] = (exp(i w) - 1) / (i w).
            let r = levin(&|_| 1.0, &|x| x, &|_| 1.0, 0.0, 1.0, w, 8).unwrap();
            let (s, c) = w.sin_cos();
            assert!(
                (r.re - s / w).abs() < 1e-12 && (r.im - (1.0 - c) / w).abs() < 1e-12,
                "w = {w}"
            );
            // integral of x exp(i w x) = exp(i w) (1/(i w) + 1/w^2) - 1/w^2.
            let r = levin(&|x| x, &|x| x, &|_| 1.0, 0.0, 1.0, w, 8).unwrap();
            let re = s / w + (c - 1.0) / (w * w);
            let im = -c / w + s / (w * w);
            assert!(
                (r.re - re).abs() < 1e-12 && (r.im - im).abs() < 1e-12,
                "w = {w}"
            );
        }
    }

    #[test]
    fn a_nonlinear_phase_costs_the_same_at_any_frequency() {
        // integral of 2 x exp(i w x^2) over [1, 2] = (exp(4 i w) - exp(i w)) / (i w).
        for w in [50.0, 5.0e4] {
            let r = levin_adaptive(
                &|x| 2.0 * x,
                &|x| x * x,
                &|x| 2.0 * x,
                1.0,
                2.0,
                w,
                1e-11,
                64,
            )
            .unwrap();
            let re = ((4.0 * w).sin() - w.sin()) / w;
            let im = (w.cos() - (4.0 * w).cos()) / w;
            assert!(
                (r.re - re).abs() < 1e-11 && (r.im - im).abs() < 1e-11,
                "w = {w}: {r:?}"
            );
        }
        // A smooth amplitude: compare with a brute-force Simpson rule.
        let w = 200.0;
        let f = |x: f64| 1.0 / (1.0 + x);
        let r = levin_adaptive(
            &f,
            &|x| x + 0.5 * x * x,
            &|x| 1.0 + x,
            0.0,
            1.0,
            w,
            1e-10,
            128,
        )
        .unwrap();
        let panels = 400_000;
        let h = 1.0 / f64::from(panels);
        let (mut re, mut im) = (0.0, 0.0);
        for k in 0..=panels {
            let x = f64::from(k) * h;
            let weight = if k == 0 || k == panels {
                1.0
            } else if k % 2 == 1 {
                4.0
            } else {
                2.0
            };
            let phase = w * (x + 0.5 * x * x);
            re += weight * f(x) * phase.cos();
            im += weight * f(x) * phase.sin();
        }
        assert!((r.re - re * h / 3.0).abs() < 1e-9 && (r.im - im * h / 3.0).abs() < 1e-9);
    }

    #[test]
    fn failures_are_reported() {
        // A stationary point at 0: collocation does not settle.
        assert!(
            levin_adaptive(
                &|_| 1.0,
                &|x| x * x,
                &|x| 2.0 * x,
                -1.0,
                1.0,
                400.0,
                1e-10,
                64
            )
            .is_err()
        );
        assert!(levin(&|_| 1.0, &|x| x, &|_| 1.0, 0.0, 1.0, 1.0, 1).is_err());
        assert!(levin(&|_| f64::NAN, &|x| x, &|_| 1.0, 0.0, 1.0, 1.0, 8).is_err());
        assert!(levin(&|_| 1.0, &|x| x, &|_| 1.0, 0.0, f64::INFINITY, 1.0, 8).is_err());
        let empty = levin(&|_| 1.0, &|x| x, &|_| 1.0, 2.0, 2.0, 5.0, 8).unwrap();
        assert_eq!((empty.re, empty.im), (0.0, 0.0));
    }
}
