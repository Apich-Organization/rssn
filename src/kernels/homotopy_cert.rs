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
//! the alpha certificate is a numerical one up to the rounding level of
//! `f64`; a safety factor on `alpha` can be requested.
//!
//! # Rigorous certificates
//!
//! [`verify_zero`] proves a simple zero instead of estimating one: it applies
//! the Krawczyk operator
//!
//! `K(X) = y - Y f(y) + (I - Y J(X)) (X - y)`
//!
//! in the outward-rounded complex rectangle arithmetic of
//! [`crate::kernels::ball`], with `y` the centre of the box `X` and `Y` an
//! approximate inverse of the Jacobian at `y`. If `K(X)` lies in the interior
//! of `X`, then `f` has exactly one zero in `X` — a rigorous existence and
//! uniqueness statement, not a floating-point check. A negative result is
//! inconclusive (the box may be too wide, the zero may be multiple, or the
//! Jacobian singular), never a disproof.
//!
//! [`certify_endpoints_rigorous`] runs this test on the path endpoints,
//! inflated and contracted as far as the operator allows. Its report
//! (`RigorousReport`) counts the certified zeros and the pairwise disjoint
//! boxes; when the number of disjoint certified boxes equals the number of
//! paths (the Bezout number) the list of solutions is complete and every
//! solution is simple. [`homotopy_solve_verified`] is the end-to-end entry
//! point, with the same escalation to extended precision as
//! [`homotopy_solve_certified`].
//!
//! The tracking steps themselves can be certified too: with
//! [`HomotopyOptions::certify_steps`] every accepted step is checked by the
//! same operator on the homogenised system at the step's `t`, so each node of
//! the path is a certified zero of `H(·, t)` (see
//! [`crate::kernels::homotopy::PathEnd::verified_steps`]). That certifies the
//! nodes, not the infinitesimal continuation between them.
//!
//! Parallelism: the Smale and Krawczyk certificates of the endpoints, and the pairwise separation
//! tests, are computed concurrently (certificate lists in endpoint order; the reports are identical
//! to a one-thread run).

use crate::kernels::ball::CTerm;
use crate::kernels::ball::Cball;
use crate::kernels::ball::eval_jacobian_terms;
use crate::kernels::ball::eval_terms;
use crate::kernels::ball::jacobian_at_terms;
use crate::kernels::ball::real_terms;
use num_complex::Complex64 as Cx;

use crate::kernels::homotopy::HomotopyOptions;
use crate::kernels::homotopy::PathEnd;
use crate::kernels::homotopy::Precision;
use crate::kernels::homotopy::homotopy_solve;
use crate::kernels::rootfind::PolyTerm;
use crate::kernels::rootfind::RootError;

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
fn solve_multi(
    a: &[Vec<Cx>],
    b: &[Vec<Cx>],
) -> Option<Vec<Vec<Cx>>> {
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
fn value_and_jacobian(
    sys: &[Vec<PolyTerm>],
    x: &[Cx],
) -> (Vec<Cx>, Vec<Vec<Cx>>) {
    let n = x.len();
    let mut f = vec![Cx::new(0.0, 0.0); sys.len()];
    let mut j = vec![vec![Cx::new(0.0, 0.0); n]; sys.len()];
    for (i, eq) in sys.iter().enumerate() {
        for t in eq {
            f[i] += Cx::new(t.coef, 0.0) * monomial(x, &t.exps, &vec![0; n]);
            for k in 0..n {
                let mut m = vec![0u32; n];
                m[k] = 1;
                j[i][k] += Cx::new(t.coef, 0.0)
                    * derivative_factor(&t.exps, &m)
                    * monomial(x, &t.exps, &m);
            }
        }
    }
    (f, j)
}

/// `x^(e - m)` (zero when some `m_j > e_j`).
fn monomial(
    x: &[Cx],
    e: &[u32],
    m: &[u32],
) -> Cx {
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
fn derivative_factor(
    e: &[u32],
    m: &[u32],
) -> f64 {
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
fn multisets(
    n: usize,
    k: u32,
) -> Vec<Vec<u32>> {
    fn rec(
        n: usize,
        idx: usize,
        left: u32,
        cur: &mut Vec<u32>,
        out: &mut Vec<Vec<u32>>,
    ) {
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
pub fn smale_alpha(
    system: &[Vec<PolyTerm>],
    x: &[(f64, f64)],
) -> Option<Certificate> {
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
                        s += Cx::new(t.coef, 0.0)
                            * derivative_factor(&t.exps, &m)
                            * monomial(&pt, &t.exps, &m);
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
    Some(Certificate {
        beta,
        gamma,
        alpha,
        certified,
        radius: if certified {
            2.0 * beta
        } else {
            f64::INFINITY
        },
    })
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

/// Number of endpoints from which certification runs in parallel.
const CERT_GRAIN: usize = 16;

/// Certifies the finite regular endpoints of a homotopy run against the
/// original `system` (Smale's alpha test plus pairwise ball disjointness).
#[must_use]
pub fn certify_endpoints(
    system: &[Vec<PolyTerm>],
    ends: &[PathEnd],
    safety: f64,
) -> CertificationReport {
    let certificates: Vec<Option<Certificate>> = super::par::map_slice(ends, CERT_GRAIN, |p| {
        if p.at_infinity || p.x.iter().any(|c| !c.0.is_finite() || !c.1.is_finite()) {
            return None;
        }
        // use the double-double value when available
        let x: Vec<(f64, f64)> =
            p.x.iter()
                .zip(&p.x_lo)
                .map(|(h, l)| (h.0 + l.0, h.1 + l.1))
                .collect();
        smale_alpha(system, &x).map(|mut c| {
            c.certified = c.alpha * safety.max(1.0) < ALPHA_0;
            if !c.certified {
                c.radius = f64::INFINITY;
            }
            c
        })
    });
    let certified = certificates
        .iter()
        .flatten()
        .filter(|c| c.certified)
        .count();
    // pairwise separation, one row per endpoint (each row is a pure function)
    let sep = super::par::map_indexed(certificates.len(), CERT_GRAIN, |i| {
        let Some(ci) = certificates[i].filter(|c| c.certified) else {
            return false;
        };
        (0..ends.len()).all(|j| {
            if j == i {
                return true;
            }
            let Some(cj) = certificates[j].filter(|c| c.certified) else {
                return true;
            };
            let d = ends[i]
                .x
                .iter()
                .zip(&ends[j].x)
                .zip(ends[i].x_lo.iter().zip(&ends[j].x_lo))
                .map(|((u, v), (ul, vl))| {
                    (u.0 - v.0 + ul.0 - vl.0).powi(2) + (u.1 - v.1 + ul.1 - vl.1).powi(2)
                })
                .sum::<f64>()
                .sqrt();
            d > ci.radius + cj.radius
        })
    });
    let distinct = sep.iter().filter(|&&b| b).count();
    CertificationReport {
        certificates,
        certified,
        distinct,
        paths: ends.len(),
        complete: distinct == ends.len(),
    }
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
pub fn homotopy_solve_certified(
    system: &[Vec<PolyTerm>],
    opts: &HomotopyOptions,
) -> Result<CertifiedSolve, RootError> {
    let ends = homotopy_solve(system, opts)?;
    let report = certify_endpoints(system, &ends, 1.0);
    if report.complete || opts.precision != Precision::Double {
        return Ok(CertifiedSolve { ends, report });
    }
    let retry = HomotopyOptions {
        precision: Precision::Adaptive,
        ..*opts
    };
    let ends2 = homotopy_solve(system, &retry)?;
    let report2 = certify_endpoints(system, &ends2, 1.0);
    if report2.distinct > report.distinct {
        Ok(CertifiedSolve {
            ends: ends2,
            report: report2,
        })
    } else {
        Ok(CertifiedSolve { ends, report })
    }
}

// ----------------------------------------------------------------------
// Rigorous (interval) certificates
// ----------------------------------------------------------------------

/// A zero of a square polynomial system together with a box that provably
/// contains it and no other zero.
#[derive(Debug, Clone, PartialEq)]
pub struct RigorousBall {
    /// An enclosing rectangle per variable (real and imaginary parts).
    pub enclosure: Vec<Cball>,
    /// The centre of the box per variable.
    pub centre: Vec<(f64, f64)>,
    /// The largest half-width of the box (an upper bound of the distance
    /// from the centre in the infinity norm over the components).
    pub radius: f64,
}

impl RigorousBall {
    /// Whether the point `(re, im)` lies in every component of the box.
    #[must_use]
    pub fn contains(
        &self,
        re: f64,
        im: f64,
    ) -> bool {
        let point = Cx::new(re, im);
        self.enclosure.iter().all(|c| c.contains(point))
    }

    /// The largest distance from the centre to `(re, im)` over the
    /// components.
    #[must_use]
    pub fn distance(
        &self,
        re: f64,
        im: f64,
    ) -> f64 {
        let point = Cx::new(re, im);
        self.centre
            .iter()
            .map(|&(cr, ci)| (Cx::new(cr, ci) - point).norm())
            .fold(0.0, f64::max)
    }
}

/// The approximate inverse of the complex Jacobian at `point`.
fn inverse_at_terms(
    system: &[Vec<CTerm>],
    point: &[Cx],
) -> Option<Vec<Vec<Cx>>> {
    let jacobian = jacobian_at_terms(system, point);
    let n = jacobian.len();
    if n == 0 || point.len() != n || jacobian.iter().any(|row| row.len() != n) {
        return None;
    }
    let identity: Vec<Vec<Cx>> = (0..n)
        .map(|i| {
            (0..n)
                .map(|c| {
                    if i == c {
                        Cx::new(1.0, 0.0)
                    } else {
                        Cx::new(0.0, 0.0)
                    }
                })
                .collect()
        })
        .collect();
    solve_multi(&jacobian, &identity)
}

/// The Krawczyk image `K(X)` of the box, with `y` the centre of `X` and `Y`
/// its approximate Jacobian inverse. Every operation is outward rounded, so
/// the result encloses the true image.
fn krawczyk_terms(
    system: &[Vec<CTerm>],
    x: &[Cball],
    y: &[Cx],
    inverse: &[Vec<Cx>],
) -> Option<Vec<Cball>> {
    let n = x.len();
    if y.len() != n || inverse.len() != n {
        return None;
    }
    let at_y: Vec<Cball> = y.iter().map(|&c| Cball::of(c)).collect();
    let value = eval_terms(system, &at_y);
    let jacobian = eval_jacobian_terms(system, x);
    if value.len() != n || jacobian.len() != n {
        return None;
    }
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let row = inverse.get(i)?;
        if row.len() != n {
            return None;
        }
        let mut k = Cball::of(*y.get(i)?);
        for j in 0..n {
            k = k - Cball::of(*row.get(j)?) * *value.get(j)?;
        }
        for j in 0..n {
            let mut m = if i == j {
                Cball::one()
            } else {
                Cball::zero()
            };
            for l in 0..n {
                let entry = jacobian.get(l)?.get(j)?;
                m = m - Cball::of(*row.get(l)?) * *entry;
            }
            k = k + m * (*x.get(j)? - Cball::of(*y.get(j)?));
        }
        if !k.is_finite() {
            return None;
        }
        out.push(k);
    }
    Some(out)
}

/// The centres of a box.
fn centres(x: &[Cball]) -> Vec<Cx> {
    x.iter().map(Cball::centre).collect()
}

/// Verifies a zero near `guess` inside the box of relative half-width
/// `relative` (a complex-coefficient system).
///
/// Returns the contracted enclosure when the Krawczyk test proves that the
/// box contains exactly one zero; `None` is inconclusive.
#[must_use]
pub fn verify_zero_terms(
    system: &[Vec<CTerm>],
    guess: &[(f64, f64)],
    relative: f64,
) -> Option<RigorousBall> {
    let n = system.len();
    if n == 0 || guess.len() != n || !(relative.is_finite() && relative > 0.0) {
        return None;
    }
    let point: Vec<Cx> = guess.iter().map(|&(re, im)| Cx::new(re, im)).collect();
    let mut current: Vec<Cball> = point
        .iter()
        .map(|c| Cball::square(c.re, c.im, relative * (1.0 + c.norm())))
        .collect();
    let mut best: Option<Vec<Cball>> = None;
    for _ in 0..32 {
        let y = centres(&current);
        let Some(inverse) = inverse_at_terms(system, &y) else {
            break;
        };
        let Some(image) = krawczyk_terms(system, &current, &y, &inverse) else {
            break;
        };
        if !image
            .iter()
            .zip(&current)
            .all(|(a, b)| a.inside_strictly(b))
        {
            break;
        }
        let Some(next) = image
            .iter()
            .zip(&current)
            .map(|(a, b)| a.intersect(b))
            .collect::<Option<Vec<Cball>>>()
        else {
            break;
        };
        let scale = next
            .iter()
            .map(|c| 1.0 + c.centre().norm())
            .fold(1.0_f64, f64::max);
        let width = next.iter().map(Cball::half_width).fold(0.0_f64, f64::max);
        best = Some(next.clone());
        current = next;
        if width <= 1e-17 * scale {
            break;
        }
    }
    let enclosure = best?;
    let centre: Vec<(f64, f64)> = enclosure
        .iter()
        .map(|c| {
            let z = c.centre();
            (z.re, z.im)
        })
        .collect();
    let radius = enclosure
        .iter()
        .map(Cball::half_width)
        .fold(0.0_f64, f64::max);
    Some(RigorousBall {
        enclosure,
        centre,
        radius,
    })
}

/// The relative radii tried by [`certify_zero`], tightest first.
const RADII: [f64; 9] = [1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 1e-1];

/// Certifies a simple zero near `guess`, trying a ladder of box sizes.
///
/// Returns `None` when no box in the ladder passes the Krawczyk test (a
/// multiple zero, a singular Jacobian, or an endpoint that is not accurate
/// enough to inflate a certified box around).
#[must_use]
pub fn certify_zero_terms(
    system: &[Vec<CTerm>],
    guess: &[(f64, f64)],
) -> Option<RigorousBall> {
    RADII
        .iter()
        .find_map(|&relative| verify_zero_terms(system, guess, relative))
}

/// Verifies a zero of a real-coefficient system near `guess`.
#[must_use]
pub fn verify_zero(
    system: &[Vec<PolyTerm>],
    guess: &[(f64, f64)],
    relative: f64,
) -> Option<RigorousBall> {
    verify_zero_terms(&real_terms(system), guess, relative)
}

/// Certifies a simple zero of a real-coefficient system near `guess`.
#[must_use]
pub fn certify_zero(
    system: &[Vec<PolyTerm>],
    guess: &[(f64, f64)],
) -> Option<RigorousBall> {
    certify_zero_terms(&real_terms(system), guess)
}

/// Rigorous certification of a list of path endpoints.
#[derive(Debug, Clone, PartialEq)]
pub struct RigorousReport {
    /// The certified ball of every endpoint (`None` for endpoints at
    /// infinity, singular endpoints and endpoints that could not be
    /// certified).
    pub certificates: Vec<Option<RigorousBall>>,
    /// Number of endpoints with a proved simple zero.
    pub certified: usize,
    /// Number of certified endpoints whose boxes are pairwise disjoint from
    /// every other certified box, i.e. proved to be different zeros.
    pub distinct: usize,
    /// Number of paths (the Bezout number for a total-degree homotopy).
    pub paths: usize,
    /// Every path ended in a proved, pairwise different simple zero.
    pub complete: bool,
}

/// Whether two certified boxes are disjoint (some component separates them,
/// which proves that they cannot contain a common point).
#[must_use]
pub fn disjoint(
    a: &RigorousBall,
    b: &RigorousBall,
) -> bool {
    a.enclosure
        .iter()
        .zip(&b.enclosure)
        .any(|(p, q)| p.intersect(q).is_none())
}

/// Number of endpoints from which interval certification runs in parallel.
const RIGOROUS_GRAIN: usize = 4;

/// Certifies the finite endpoints of a homotopy run rigorously.
///
/// Endpoints at infinity, ill-conditioned endpoints and the centroids of
/// clusters of close roots are reported as uncertified; a certified endpoint
/// is a proved simple zero.
#[must_use]
pub fn certify_endpoints_rigorous(
    system: &[Vec<PolyTerm>],
    ends: &[PathEnd],
) -> RigorousReport {
    let certificates: Vec<Option<RigorousBall>> =
        super::par::map_slice(ends, RIGOROUS_GRAIN, |end| {
            if end.at_infinity || end.x.iter().any(|c| !c.0.is_finite() || !c.1.is_finite()) {
                return None;
            }
            let x: Vec<(f64, f64)> = end
                .x
                .iter()
                .zip(&end.x_lo)
                .map(|(h, l)| (h.0 + l.0, h.1 + l.1))
                .collect();
            certify_zero(system, &x)
        });
    let certified = certificates.iter().filter(|c| c.is_some()).count();
    let distinct = super::par::map_indexed(ends.len(), RIGOROUS_GRAIN * 8, |i| {
        let Some(bi) = certificates[i].as_ref() else {
            return false;
        };
        (0..ends.len())
            .all(|j| j == i || certificates[j].as_ref().is_none_or(|bj| disjoint(bi, bj)))
    })
    .into_iter()
    .filter(|&b| b)
    .count();
    RigorousReport {
        certificates,
        certified,
        distinct,
        paths: ends.len(),
        complete: distinct == ends.len(),
    }
}

/// Result of [`homotopy_solve_verified`].
#[derive(Debug, Clone, PartialEq)]
pub struct VerifiedSolve {
    /// The path endpoints.
    pub ends: Vec<PathEnd>,
    /// The rigorous certification.
    pub rigorous: RigorousReport,
    /// The floating-point certification (fast pre-filter).
    pub floating: CertificationReport,
}

/// Runs [`homotopy_solve`] and certifies the endpoints rigorously.
///
/// Tracking steps are certified as well
/// ([`HomotopyOptions::certify_steps`] is turned on for the run, which the
/// double-precision arithmetic supports). Runs with [`Precision::Double`]
/// escalate to [`Precision::Adaptive`] when the rigorous report is
/// incomplete and the escalated run certifies more distinct zeros.
///
/// # Errors
/// As [`homotopy_solve`].
pub fn homotopy_solve_verified(
    system: &[Vec<PolyTerm>],
    opts: &HomotopyOptions,
) -> Result<VerifiedSolve, RootError> {
    let certified_opts = HomotopyOptions {
        certify_steps: true,
        ..*opts
    };
    let ends = homotopy_solve(system, &certified_opts)?;
    let rigorous = certify_endpoints_rigorous(system, &ends);
    let floating = certify_endpoints(system, &ends, 1.0);
    if rigorous.complete || opts.precision != Precision::Double {
        return Ok(VerifiedSolve {
            ends,
            rigorous,
            floating,
        });
    }
    let retry = HomotopyOptions {
        precision: Precision::Adaptive,
        ..certified_opts
    };
    let ends2 = homotopy_solve(system, &retry)?;
    let rigorous2 = certify_endpoints_rigorous(system, &ends2);
    let floating2 = certify_endpoints(system, &ends2, 1.0);
    if rigorous2.distinct > rigorous.distinct {
        Ok(VerifiedSolve {
            ends: ends2,
            rigorous: rigorous2,
            floating: floating2,
        })
    } else {
        Ok(VerifiedSolve {
            ends,
            rigorous,
            floating,
        })
    }
}

/// The homogenised system `H(X, t) = 0` of the total-degree homotopy at a
/// complex `t`, as complex-coefficient terms in the `n + 1` homogeneous
/// variables.
///
/// `H_i = γ (1 - t)(X_{i+1}^{d_i} - X_0^{d_i}) + t f_i(X)` for
/// `i = 0..n`, plus the patch equation `-1 + Σ_k a_k X_k` (the same rows as
/// [`crate::kernels::homotopy`] evaluates). `None` for a malformed system.
#[must_use]
pub fn homotopy_node_terms(
    homogenised: &[Vec<PolyTerm>],
    degrees: &[u32],
    gamma: Cx,
    patch: &[Cx],
    t: Cx,
) -> Option<Vec<Vec<CTerm>>> {
    let n = homogenised.len();
    if n == 0 || degrees.len() != n || patch.len() != n + 1 || homogenised.iter().any(Vec::is_empty)
    {
        return None;
    }
    let one_minus_t = Cx::new(1.0, 0.0) - t;
    let mut out: Vec<Vec<CTerm>> = Vec::with_capacity(n + 1);
    for (i, equation) in homogenised.iter().enumerate() {
        let degree = *degrees.get(i)?;
        let mut row: Vec<CTerm> = Vec::with_capacity(equation.len() + 2);
        let mut high = vec![0_u32; n + 1];
        *high.get_mut(i + 1)? = degree;
        row.push(CTerm {
            coef: gamma * one_minus_t,
            exps: high,
        });
        let mut low = vec![0_u32; n + 1];
        *low.first_mut()? = degree;
        row.push(CTerm {
            coef: -(gamma * one_minus_t),
            exps: low,
        });
        for term in equation {
            if term.exps.len() != n + 1 {
                return None;
            }
            row.push(CTerm {
                coef: Cx::new(term.coef, 0.0) * t,
                exps: term.exps.clone(),
            });
        }
        out.push(row);
    }
    let mut patch_row: Vec<CTerm> = Vec::with_capacity(n + 2);
    patch_row.push(CTerm {
        coef: Cx::new(-1.0, 0.0),
        exps: vec![0_u32; n + 1],
    });
    for (k, &a) in patch.iter().enumerate() {
        let mut exps = vec![0_u32; n + 1];
        *exps.get_mut(k)? = 1;
        patch_row.push(CTerm { coef: a, exps });
    }
    out.push(patch_row);
    Some(out)
}

/// Certifies one node of a homotopy path: proves that `H(·, t)` has exactly
/// one zero in a small box around the homogeneous point `x`.
///
/// `homogenised` and `degrees` are the data of
/// [`crate::kernels::homotopy::homotopy_solve`]'s problem; `x` has `n + 1`
/// entries. Used by the tracker when
/// [`crate::kernels::homotopy::HomotopyOptions::certify_steps`] is set.
#[must_use]
pub fn certify_homotopy_node(
    homogenised: &[Vec<PolyTerm>],
    degrees: &[u32],
    gamma: Cx,
    patch: &[Cx],
    x: &[(f64, f64)],
    t: (f64, f64),
) -> Option<RigorousBall> {
    let system = homotopy_node_terms(homogenised, degrees, gamma, patch, Cx::new(t.0, t.1))?;
    certify_zero_terms(&system, x)
}
