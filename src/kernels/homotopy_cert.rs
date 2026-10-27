//! # Smale alpha-theory certificates for homotopy endpoints
//!
//! [`smale_alpha`] evaluates, for a square polynomial system `f` and a
//! point `x`, Smale's invariants
//!
//! * `beta = |Df(x)^{-1} f(x)|` (length of the Newton step),
//! * `gamma = max_{2 <= k <= d} (|Df(x)^{-1} D^k f(x) / k!|)^{1/(k-1)}`
//!   (the operator norms of the multilinear maps are bounded above by their
//!   Frobenius norms, so `gamma` is an upper bound; the derivative tensors of
//!   a polynomial are exact and finite),
//! * `alpha = beta gamma`.
//!
//! If `alpha < alpha_0 = (13 - 3 sqrt(17)) / 4 ~ 0.1577` then `x` is an
//! *approximate zero*: Newton's method started at `x` converges quadratically
//! to a zero `zeta` of `f` with `|x - zeta| <= 2 beta`
//! ([`Certificate::radius`]). Two approximate zeros whose balls
//! `B(x_i, 2 beta_i)` are disjoint belong to different zeros, so a set of
//! endpoints that are pairwise separated in this sense certifies that the
//! paths did not jump onto each other; when their number equals the Bezout
//! number the list of solutions is complete
//! ([`certify_endpoints`], [`homotopy_solve_certified`]).
//!
//! The invariants are evaluated in floating point (no directed rounding), so
//! the certificate is a numerical one up to the rounding level of `f64`; a
//! safety factor on `alpha` can be requested.

use num_complex::Complex64 as Cx;

use crate::kernels::homotopy::{HomotopyOptions, PathEnd, Precision, homotopy_solve};
use crate::kernels::rootfind::{PolyTerm, RootError};

/// `(13 - 3 sqrt(17)) / 4`: the threshold of Smale's alpha theorem.
pub const ALPHA_0: f64 = 0.157_670_780_786_754_6;

/// Smale invariants at a point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Certificate {
    /// Length of the Newton step.
    pub beta: f64,
    /// Upper bound of Smale's `gamma`.
    pub gamma: f64,
    /// `beta * gamma`.
    pub alpha: f64,
    /// `alpha < ALPHA_0`: the point is an approximate zero.
    pub certified: bool,
    /// Radius `2 beta` of a ball around the point that contains the
    /// associated zero (infinite when not certified).
    pub radius: f64,
}

/// Solves `A X = B` for a small complex system (`B` has several columns).
fn solve_multi(a: &[Vec<Cx>], b: &[Vec<Cx>]) -> Option<Vec<Vec<Cx>>> {
    let n = a.len();
    let mut a: Vec<Vec<Cx>> = a.to_vec();
    let mut b: Vec<Vec<Cx>> = b.to_vec();
    for k in 0..n {
        let p = (k..n).max_by(|&i, &j| a[i][k].norm().total_cmp(&a[j][k].norm()))?;
        if a[p][k].norm() < 1e-300 {
            return None;
        }
        a.swap(k, p);
        b.swap(k, p);
        for i in k + 1..n {
            let f = a[i][k] / a[k][k];
            for c in k..n {
                let t = a[k][c];
                a[i][c] -= f * t;
            }
            for c in 0..b[k].len() {
                let t = b[k][c];
                b[i][c] -= f * t;
            }
        }
    }
    let cols = b.first().map_or(0, Vec::len);
    let mut x = vec![vec![Cx::new(0.0, 0.0); cols]; n];
    for c in 0..cols {
        for i in (0..n).rev() {
            let mut s = b[i][c];
            for j in i + 1..n {
                s -= a[i][j] * x[j][c];
            }
            x[i][c] = s / a[i][i];
        }
    }
    Some(x)
}

/// Value and Jacobian of the system at `x`.
fn value_and_jacobian(sys: &[Vec<PolyTerm>], x: &[Cx]) -> (Vec<Cx>, Vec<Vec<Cx>>) {
    let n = x.len();
    let mut f = vec![Cx::new(0.0, 0.0); sys.len()];
    let mut j = vec![vec![Cx::new(0.0, 0.0); n]; sys.len()];
    for (i, eq) in sys.iter().enumerate() {
        for t in eq {
            f[i] += Cx::new(t.coef, 0.0) * monomial(x, &t.exps, &vec![0; n]);
            for k in 0..n {
                let mut m = vec![0u32; n];
                m[k] = 1;
                j[i][k] += Cx::new(t.coef, 0.0) * derivative_factor(&t.exps, &m) * monomial(x, &t.exps, &m);
            }
        }
    }
    (f, j)
}

/// `x^(e - m)` (zero when some `m_j > e_j`).
fn monomial(x: &[Cx], e: &[u32], m: &[u32]) -> Cx {
    let mut r = Cx::new(1.0, 0.0);
    for ((xj, &ej), &mj) in x.iter().zip(e).zip(m) {
        if mj > ej {
            return Cx::new(0.0, 0.0);
        }
        r *= xj.powu(ej - mj);
    }
    r
}

/// `prod_j e_j! / (e_j - m_j)!`.
fn derivative_factor(e: &[u32], m: &[u32]) -> f64 {
    let mut f = 1.0;
    for (&ej, &mj) in e.iter().zip(m) {
        if mj > ej {
            return 0.0;
        }
        for t in 0..mj {
            f *= f64::from(ej - t);
        }
    }
    f
}

/// Enumerates the multisets (multiplicity vectors) of size `k` over `n` indices.
fn multisets(n: usize, k: u32) -> Vec<Vec<u32>> {
    fn rec(n: usize, idx: usize, left: u32, cur: &mut Vec<u32>, out: &mut Vec<Vec<u32>>) {
        if idx + 1 == n {
            cur[idx] = left;
            out.push(cur.clone());
            return;
        }
        for c in 0..=left {
            cur[idx] = c;
            rec(n, idx + 1, left - c, cur, out);
        }
    }
    let mut out = Vec::new();
    let mut cur = vec![0u32; n];
    if n > 0 {
        rec(n, 0, k, &mut cur, &mut out);
    }
    out
}

fn factorial(k: u32) -> f64 {
    (1..=k).map(f64::from).product()
}

/// Smale's `beta`, `gamma` and `alpha` for the square polynomial system
/// `system` at the point `x` (one complex number per variable).
///
/// Returns `None` for a non-square or malformed system, or a singular
/// Jacobian at `x`.
#[must_use]
pub fn smale_alpha(system: &[Vec<PolyTerm>], x: &[(f64, f64)]) -> Option<Certificate> {
    let n = system.len();
    if n == 0 || x.len() != n || system.iter().any(|eq| eq.iter().any(|t| t.exps.len() != n)) {
        return None;
    }
    let pt: Vec<Cx> = x.iter().map(|&(re, im)| Cx::new(re, im)).collect();
    let (fv, jac) = value_and_jacobian(system, &pt);
    let rhs: Vec<Vec<Cx>> = fv.iter().map(|&v| vec![v]).collect();
    let step = solve_multi(&jac, &rhs)?;
    let beta = step.iter().map(|r| r[0].norm_sqr()).sum::<f64>().sqrt();
    let dmax = system
        .iter()
        .flat_map(|eq| eq.iter().map(|t| t.exps.iter().sum::<u32>()))
        .max()
        .unwrap_or(0);
    let mut gamma = 0.0_f64;
    for k in 2..=dmax {
        let kf = factorial(k);
        let mut fro2 = 0.0;
        for m in multisets(n, k) {
            // number of ordered index tuples with this multiplicity pattern
            let weight = kf / m.iter().map(|&c| factorial(c)).product::<f64>();
            // column of T^{(k)}: (1/k!) d^m f_i for every equation i
            let col: Vec<Vec<Cx>> = (0..n)
                .map(|i| {
                    let mut s = Cx::new(0.0, 0.0);
                    for t in &system[i] {
                        s += Cx::new(t.coef, 0.0) * derivative_factor(&t.exps, &m) * monomial(&pt, &t.exps, &m);
                    }
                    vec![s / kf]
                })
                .collect();
            let w = solve_multi(&jac, &col)?;
            fro2 += weight * w.iter().map(|r| r[0].norm_sqr()).sum::<f64>();
        }
        gamma = gamma.max(fro2.sqrt().powf(1.0 / f64::from(k - 1)));
    }
    let alpha = beta * gamma;
    let certified = alpha.is_finite() && alpha < ALPHA_0;
    Some(Certificate { beta, gamma, alpha, certified, radius: if certified { 2.0 * beta } else { f64::INFINITY } })
}

/// Outcome of certifying a list of path endpoints.
#[derive(Debug, Clone, PartialEq)]
pub struct CertificationReport {
    /// Certificate of every endpoint (`None` for endpoints at infinity,
    /// singular endpoints, or a singular Jacobian).
    pub certificates: Vec<Option<Certificate>>,
    /// Number of endpoints that are certified approximate zeros.
    pub certified: usize,
    /// Number of certified endpoints whose balls `B(x, 2 beta)` are disjoint
    /// from the balls of all other certified endpoints, i.e. certified to
    /// belong to pairwise different zeros.
    pub distinct: usize,
    /// Number of paths (the Bezout number for a total-degree homotopy).
    pub paths: usize,
    /// Every path ended in a certified, pairwise distinct zero: the solution
    /// list is complete and all solutions are simple.
    pub complete: bool,
}

/// Certifies the finite regular endpoints of a homotopy run against the
/// original `system` (Smale's alpha test plus pairwise ball disjointness).
#[must_use]
pub fn certify_endpoints(system: &[Vec<PolyTerm>], ends: &[PathEnd], safety: f64) -> CertificationReport {
    let certificates: Vec<Option<Certificate>> = ends
        .iter()
        .map(|p| {
            if p.at_infinity || p.x.iter().any(|c| !c.0.is_finite() || !c.1.is_finite()) {
                return None;
            }
            // use the double-double value when available
            let x: Vec<(f64, f64)> = p.x.iter().zip(&p.x_lo).map(|(h, l)| (h.0 + l.0, h.1 + l.1)).collect();
            smale_alpha(system, &x).map(|mut c| {
                c.certified = c.alpha * safety.max(1.0) < ALPHA_0;
                if !c.certified {
                    c.radius = f64::INFINITY;
                }
                c
            })
        })
        .collect();
    let certified = certificates.iter().flatten().filter(|c| c.certified).count();
    let mut distinct = 0;
    for (i, ci) in certificates.iter().enumerate() {
        let Some(ci) = ci.filter(|c| c.certified) else { continue };
        let separated = (0..ends.len()).all(|j| {
            if j == i {
                return true;
            }
            let Some(cj) = certificates[j].filter(|c| c.certified) else { return true };
            let d = ends[i]
                .x
                .iter()
                .zip(&ends[j].x)
                .zip(ends[i].x_lo.iter().zip(&ends[j].x_lo))
                .map(|((u, v), (ul, vl))| (u.0 - v.0 + ul.0 - vl.0).powi(2) + (u.1 - v.1 + ul.1 - vl.1).powi(2))
                .sum::<f64>()
                .sqrt();
            d > ci.radius + cj.radius
        });
        if separated {
            distinct += 1;
        }
    }
    CertificationReport { certificates, certified, distinct, paths: ends.len(), complete: distinct == ends.len() }
}

/// Result of [`homotopy_solve_certified`].
#[derive(Debug, Clone, PartialEq)]
pub struct CertifiedSolve {
    /// The path endpoints.
    pub ends: Vec<PathEnd>,
    /// Their certification.
    pub report: CertificationReport,
}

/// Runs [`homotopy_solve`] and certifies the endpoints.
///
/// Certification is [`certify_endpoints`]. If the run is in plain `f64` and not every path ends
/// in a certified, distinct zero, the whole solve is repeated with
/// [`crate::kernels::homotopy::Precision::Adaptive`] and the more complete
/// result is returned.
///
/// # Errors
/// As [`homotopy_solve`].
pub fn homotopy_solve_certified(system: &[Vec<PolyTerm>], opts: &HomotopyOptions) -> Result<CertifiedSolve, RootError> {
    let ends = homotopy_solve(system, opts)?;
    let report = certify_endpoints(system, &ends, 1.0);
    if report.complete || opts.precision != Precision::Double {
        return Ok(CertifiedSolve { ends, report });
    }
    let retry = HomotopyOptions { precision: Precision::Adaptive, ..*opts };
    let ends2 = homotopy_solve(system, &retry)?;
    let report2 = certify_endpoints(system, &ends2, 1.0);
    if report2.distinct > report.distinct {
        Ok(CertifiedSolve { ends: ends2, report: report2 })
    } else {
        Ok(CertifiedSolve { ends, report })
    }
}
