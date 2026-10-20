//! Eigenfunction families beyond sines, cosines and plain Bessel modes, and
//! closed-form projection integrals.
//!
//! * Robin intervals with a zero or a negative eigenvalue: the Robin
//!   problem `X'' = -k² X`, `p₀ X + q₀ X' = 0`, `p₁ X + q₁ X' = 0` has, for
//!   conditions of the "wrong" sign, one mode with `k = 0` (`q₀ - p₀ ξ`) or
//!   `k = iκ` (`q₀ κ cosh κξ - p₀ sinh κξ`, eigenvalue `-κ²`), found
//!   numerically (`sl_neg_root`) and written out as the first mode.
//! * Annuli: the radial modes are cross-products
//!   `F = A J_ν(k r) - B Y_ν(k r)` with `A`, `B` fixed by the condition at
//!   the inner radius and `k = annulus_root(..)` by the one at the outer;
//!   shells with `l = 0` reduce to an interval for `r u`.
//! * Associated Legendre functions `P_l^m(cos θ) = sin^m θ P_l^{(m)}(cos θ)`
//!   as explicit polynomials, and projection integrals of polynomial data
//!   (in `cos θ` and `sin² θ`) in closed form.
//! * Fourier–Bessel projections of polynomial data `∫ r^k J_ν(κ r) dr` by
//!   the recurrences `∫ z^μ Z_ν = z^μ Z_{ν+1} - (μ-ν-1) ∫ z^{μ-1} Z_{ν+1}`
//!   and `∫ z^μ Z_ν = -z^μ Z_{ν-1} + (μ+ν-1) ∫ z^{μ-1} Z_{ν-1}`, which
//!   terminate when `μ - ν - 1` (or `μ + ν - 1`) is a non-negative even
//!   integer.
//! * Decaying and outgoing radial profiles outside a disk or a ball.

use super::spectral::End;
use super::spectral::Interval;
use super::spectral::Kind;
use super::spectral::interval_families;
use super::util::call;
use super::util::div;
use super::util::fraction;
use super::util::imaginary;
use super::util::is_zero_number;
use super::util::sqrt;
use crate::graph::Cx;
use crate::graph::Env;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::rules::calculus::derivative;
use crate::rules::complex::build::add;
use crate::rules::complex::build::mul;
use crate::rules::complex::build::neg;
use crate::rules::complex::build::pow;
use crate::rules::complex::build::powi;
use crate::rules::complex::build::sub;
use crate::rules::poly::repr::Gens;
use crate::rules::poly::repr::Limits;
use crate::rules::poly::repr::Poly;
use crate::rules::poly::repr::from_term;
use crate::rules::poly::repr::to_term;

/// A cylinder-function mode `r^(power2/2) Σ c_i Z_ν(κ r)` (`Z` the Bessel
/// function of the first kind, or of the second kind when the flag is
/// set): enough structure to project polynomial data in closed form.
#[derive(Clone, Debug)]
pub(super) struct Cylinder {
    pub nu: NodeId,
    pub kappa: NodeId,
    pub power2: i64,
    pub parts: Vec<(NodeId, bool)>,
}

/// One family of eigenfunctions along an axis, for a given index.
#[derive(Clone, Debug)]
pub(super) struct Family {
    pub mode: NodeId,
    /// The eigenvalue contributed to `-Δ`.
    pub eigen: NodeId,
    /// First index of the family.
    pub first: i64,
    /// How many leading indices are always written out.
    pub forced: i64,
    /// `∫ w φ² dx` when known in closed form for this index.
    pub norm: Option<NodeId>,
    pub cylinder: Option<Cylinder>,
}

impl Family {
    pub(super) const fn new(
        mode: NodeId,
        eigen: NodeId,
        first: i64,
        forced: i64,
        norm: Option<NodeId>,
    ) -> Self {
        Self { mode, eigen, first, forced, norm, cylinder: None }
    }
}

/// The value of a node without free symbols.
pub(super) fn exact_value(
    cx: &Cx<'_>,
    node: NodeId,
) -> Option<f64> {
    if !cx.graph.free_symbols(cx.graph.find(node)).is_empty() {
        return None;
    }
    cx.graph.eval(node, &Env::numeric(0.0)).filter(|v| v.is_finite())
}

/// The kind of condition `p u + q u'`.
pub(super) fn kind_of(
    cx: &mut Cx<'_>,
    p: NodeId,
    q: NodeId,
) -> Kind {
    if cx.is_zero(q) {
        Kind::Dirichlet
    } else if cx.is_zero(p) {
        Kind::Neumann
    } else {
        Kind::Robin
    }
}

// ----------------------------------------------------------------------
// Robin intervals
// ----------------------------------------------------------------------

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum Special {
    Zero,
    /// The `k`-th negative eigenvalue (`k ≥ 1`, ascending in `κ`).
    Negative(u32),
}

/// The modes of a Robin problem with eigenvalue zero or negative: whether
/// there is a zero mode, and how many negative ones (decided numerically
/// when the coefficients are numbers).
fn special_modes(
    cx: &Cx<'_>,
    ends: [NodeId; 5],
) -> Option<(bool, u32)> {
    let v: Vec<f64> = ends.iter().map(|&e| exact_value(cx, e)).collect::<Option<_>>()?;
    let (p0, q0, p1, q1, l) = (v[0], v[1], v[2], v[3], v[4]);
    let determinant = p0 * (p1 * l + q1) - q0 * p1;
    let scale = (p0.abs() + q0.abs()) * (p1.abs() * l + q1.abs()) + (q0 * p1).abs();
    let zero = determinant.abs() <= 1e-12 * scale.max(1e-300);
    let count = super::numeric::negative_count(p0, q0, p1, q1, l);
    (zero || count > 0).then_some((zero, count))
}

/// The families of a Robin interval (index `n`; `numeric` when it is a
/// number), including the zero and negative modes when there are some.
pub(super) fn robin_families(
    cx: &mut Cx<'_>,
    iv: &Interval,
    xi: NodeId,
    n: NodeId,
    numeric: bool,
) -> Option<Vec<Family>> {
    let (p0, q0, p1, q1, l) = (iv.left.p, iv.left.q, iv.right.p, iv.right.q, iv.length);
    let specials = special_modes(cx, [p0, q0, p1, q1, l]);
    let (zero_mode, negatives) = specials.unwrap_or((false, 0));
    let total = i64::from(zero_mode) + i64::from(negatives);
    let literal = if numeric { cx.graph.number_of(n).and_then(Number::to_i64) } else { None };
    let chosen = literal.filter(|&m| m < total).map(|m| {
        if zero_mode && m == 0 {
            Special::Zero
        } else {
            Special::Negative(u32::try_from(m - i64::from(zero_mode) + 1).unwrap_or(1))
        }
    });
    if let Some(kind) = chosen {
        let zero = cx.graph.int(0);
        return Some(vec![match kind {
            | Special::Zero => {
                // X = q₀ - p₀ ξ.
                let pxi = mul(cx.graph, &[p0, xi]);
                let mode = sub(cx.graph, q0, pxi);
                let (q2, p2, l2, l3) = (powi(cx.graph, q0, 2), powi(cx.graph, p0, 2), powi(cx.graph, l, 2), powi(cx.graph, l, 3));
                let third = fraction(cx, 1, 3)?;
                let a = mul(cx.graph, &[q2, l]);
                let b = mul(cx.graph, &[q0, p0, l2]);
                let c = mul(cx.graph, &[third, p2, l3]);
                let nb = neg(cx.graph, b);
                let norm = add(cx.graph, &[a, nb, c]);
                Family::new(cx.simplify(mode), zero, 0, total, Some(cx.simplify(norm)))
            },
            | Special::Negative(k) => {
                let placeholder = if k == 1 {
                    call(cx, "sl_neg_root", &[p0, q0, p1, q1, l])?
                } else {
                    let index = cx.graph.int(i64::from(k));
                    call(cx, "sl_neg_root_n", &[p0, q0, p1, q1, l, index])?
                };
                // A root that is a small rational (the drift gauge gives |α|)
                // is written exactly: its projection integrals then close.
                let exact = exact_value(cx, placeholder).and_then(|v| {
                    (1..=12_i64).find_map(|q| {
                        #[allow(clippy::cast_precision_loss, clippy::cast_possible_truncation)]
                        let scaled = (v * q as f64).round();
                        #[allow(clippy::cast_precision_loss, clippy::cast_possible_truncation)]
                        let close = (v * q as f64 - scaled).abs() < 1e-9;
                        #[allow(clippy::cast_possible_truncation)]
                        close.then_some((scaled as i64, q))
                    })
                });
                let kappa = match exact {
                    | Some((n, d)) => cx.graph.num(Number::fraction(n, d)?),
                    | None => placeholder,
                };
                let k_xi = mul(cx.graph, &[kappa, xi]);
                let (ch, sh) = (call(cx, "cosh", &[k_xi])?, call(cx, "sinh", &[k_xi])?);
                let a_coeff = mul(cx.graph, &[q0, kappa]);
                let b_coeff = neg(cx.graph, p0);
                let t1 = mul(cx.graph, &[a_coeff, ch]);
                let t2 = mul(cx.graph, &[b_coeff, sh]);
                let mode = add(cx.graph, &[t1, t2]);
                let k2 = powi(cx.graph, kappa, 2);
                let eigen = neg(cx.graph, k2);
                // ∫ X² with X = A cosh + B sinh.
                let two_kl = {
                    let two = cx.graph.int(2);
                    mul(cx.graph, &[two, kappa, l])
                };
                let (s2, c2) = (call(cx, "sinh", &[two_kl])?, call(cx, "cosh", &[two_kl])?);
                let four_k = {
                    let four = cx.graph.int(4);
                    mul(cx.graph, &[four, kappa])
                };
                let s_over = div(cx, s2, four_k);
                let half_l = {
                    let half = fraction(cx, 1, 2)?;
                    mul(cx.graph, &[half, l])
                };
                let minus_half_l = neg(cx.graph, half_l);
                let cosh_part = add(cx.graph, &[half_l, s_over]);
                let sinh_part = add(cx.graph, &[minus_half_l, s_over]);
                let one = cx.graph.int(1);
                let c2m = sub(cx.graph, c2, one);
                let cross = div(cx, c2m, four_k);
                let (a2, b2) = (powi(cx.graph, a_coeff, 2), powi(cx.graph, b_coeff, 2));
                let two = cx.graph.int(2);
                let n1 = mul(cx.graph, &[a2, cosh_part]);
                let n2 = mul(cx.graph, &[b2, sinh_part]);
                let n3 = mul(cx.graph, &[two, a_coeff, b_coeff, cross]);
                let norm = add(cx.graph, &[n1, n2, n3]);
                Family::new(cx.simplify(mode), cx.simplify(eigen), 0, total, Some(norm))
            },
        }]);
    }
    // A Robin end: X = q0 k cos(kξ) − p0 sin(kξ).
    // Positive roots are numbered from 1; the indices below `total` are the
    // zero and negative modes.
    let root_index = if total > 1 {
        let shift = cx.graph.int(1 - total);
        let shifted = add(cx.graph, &[n, shift]);
        cx.simplify(shifted)
    } else {
        n
    };
    let k = call(cx, "sl_root", &[p0, q0, p1, q1, l, root_index])?;
    let kx = mul(cx.graph, &[k, xi]);
    let (cos, sin) = (call(cx, "cos", &[kx])?, call(cx, "sin", &[kx])?);
    let a = mul(cx.graph, &[q0, k, cos]);
    let b = mul(cx.graph, &[p0, sin]);
    let mode = sub(cx.graph, a, b);
    let e = powi(cx.graph, k, 2);
    let e = cx.simplify(e);
    // ∫_0^L X² = (q₀²k² + p₀²) L/2 + (q₀²k² - p₀²) sin(2kL)/(4k) - q₀ p₀ sin²(kL).
    let norm = {
        let (q2, p2, k2) = (powi(cx.graph, q0, 2), powi(cx.graph, p0, 2), powi(cx.graph, k, 2));
        let q2k2 = mul(cx.graph, &[q2, k2]);
        let plus = add(cx.graph, &[q2k2, p2]);
        let minus = {
            let np2 = neg(cx.graph, p2);
            add(cx.graph, &[q2k2, np2])
        };
        let half_l = {
            let half = fraction(cx, 1, 2)?;
            mul(cx.graph, &[half, l])
        };
        let first_term = mul(cx.graph, &[plus, half_l]);
        let two_kl = {
            let two = cx.graph.int(2);
            mul(cx.graph, &[two, k, l])
        };
        let s2 = call(cx, "sin", &[two_kl])?;
        let four_k = {
            let four = cx.graph.int(4);
            mul(cx.graph, &[four, k])
        };
        let ratio = div(cx, s2, four_k);
        let second_term = mul(cx.graph, &[minus, ratio]);
        let kl = mul(cx.graph, &[k, l]);
        let s = call(cx, "sin", &[kl])?;
        let s_sq = powi(cx.graph, s, 2);
        let cross = mul(cx.graph, &[q0, p0, s_sq]);
        let third_term = neg(cx.graph, cross);
        add(cx.graph, &[first_term, second_term, third_term])
    };
    let (first, forced) = if total > 0 { (0, total) } else { (1, 0) };
    Some(vec![Family::new(mode, e, first, forced, Some(norm))])
}

/// The radial mode of a disk or ball with a non-positive eigenvalue: for
/// `sign = 0` (`h + ν = 0`) the power `r^ν`, for `sign = -1` (`h + ν < 0`)
/// the modified Bessel function `I_ν(κ r)` with `κ = bessel_neg_root(ν, h)/R`
/// and the eigenvalue `-κ²`. `shifted` is `h R` (with the ball's shift).
pub(super) fn radial_non_positive(
    cx: &mut Cx<'_>,
    r: NodeId,
    radius: NodeId,
    nu: NodeId,
    shifted: NodeId,
    sphere: bool,
    sign: i32,
) -> Option<Family> {
    let zero = cx.graph.int(0);
    let one = cx.graph.int(1);
    let two = cx.graph.int(2);
    let minus_half = fraction(cx, -1, 2)?;
    let with_shell = |cx: &mut Cx<'_>, f: NodeId| {
        if sphere {
            let root = pow(cx.graph, r, minus_half);
            mul(cx.graph, &[root, f])
        } else {
            f
        }
    };
    if sign == 0 {
        let power = pow(cx.graph, r, nu);
        let mode = with_shell(cx, power);
        // ∫ r^{2ν+1} dr = R^{2ν+2} / (2ν + 2).
        let two_nu = mul(cx.graph, &[two, nu]);
        let e = add(cx.graph, &[two_nu, two]);
        let top = pow(cx.graph, radius, e);
        let norm = div(cx, top, e);
        return Some(Family::new(cx.simplify(mode), zero, 0, 1, Some(cx.simplify(norm))));
    }
    let root = call(cx, "bessel_neg_root", &[nu, shifted])?;
    let inverse_radius = powi(cx.graph, radius, -1);
    let kappa = mul(cx.graph, &[root, inverse_radius]);
    let kr = mul(cx.graph, &[kappa, r]);
    let i_nu = call(cx, "besseli", &[nu, kr])?;
    let mode = with_shell(cx, i_nu);
    let k2 = powi(cx.graph, kappa, 2);
    let eigen = neg(cx.graph, k2);
    // ∫ r I_ν(κ r)² dr = R²/2 [(1 + ν²/z²) I_ν(z)² - I_ν'(z)²], z = κ R.
    let up = add(cx.graph, &[nu, one]);
    let (i_z, i_up) = (call(cx, "besseli", &[nu, root])?, call(cx, "besseli", &[up, root])?);
    let ratio = div(cx, nu, root);
    let slope = {
        let t = mul(cx.graph, &[ratio, i_z]);
        add(cx.graph, &[i_up, t])
    };
    let nu2_over_z2 = powi(cx.graph, ratio, 2);
    let factor = add(cx.graph, &[one, nu2_over_z2]);
    let lead = {
        let sq = powi(cx.graph, i_z, 2);
        mul(cx.graph, &[factor, sq])
    };
    let tail = {
        let sq = powi(cx.graph, slope, 2);
        neg(cx.graph, sq)
    };
    let bracket = add(cx.graph, &[lead, tail]);
    let half = fraction(cx, 1, 2)?;
    let r2 = powi(cx.graph, radius, 2);
    let norm = mul(cx.graph, &[half, r2, bracket]);
    Some(Family::new(cx.simplify(mode), cx.simplify(eigen), 0, 1, Some(norm)))
}

// ----------------------------------------------------------------------
// Annuli
// ----------------------------------------------------------------------

/// The `(p, q)` of a Robin end after the shift `u = r^(-1/2) F` of a ball.
fn shifted_p(
    cx: &mut Cx<'_>,
    end: &End,
    radius: NodeId,
    sphere: bool,
) -> NodeId {
    if !sphere {
        return end.p;
    }
    // p u + q u' = r^(-1/2) [(p - q/(2r)) F + q F'].
    let half = fraction(cx, -1, 2).unwrap_or(end.p);
    let q_over_r = div(cx, end.q, radius);
    let shift = mul(cx.graph, &[half, q_over_r]);
    let p = add(cx.graph, &[end.p, shift]);
    cx.simplify(p)
}

fn bessel(
    cx: &mut Cx<'_>,
    second: bool,
    order: NodeId,
    z: NodeId,
) -> Option<NodeId> {
    call(cx, if second { "bessely" } else { "besselj" }, &[order, z])
}

/// `Z_ν'(z) = -Z_{ν+1}(z) + (ν/z) Z_ν(z)`.
fn bessel_slope(
    cx: &mut Cx<'_>,
    second: bool,
    nu: NodeId,
    z: NodeId,
) -> Option<NodeId> {
    let one = cx.graph.int(1);
    let up = add(cx.graph, &[nu, one]);
    let z_up = bessel(cx, second, up, z)?;
    let z_nu = bessel(cx, second, nu, z)?;
    let ratio = div(cx, nu, z);
    let a = mul(cx.graph, &[ratio, z_nu]);
    let b = neg(cx.graph, z_up);
    Some(add(cx.graph, &[a, b]))
}

/// Whether the annulus has a zero mode and how many modes have a negative
/// eigenvalue (Robin conditions of the wrong sign). Symbolic radii are not
/// decided; a symbolic angular order is refused (`None`) when some integer
/// order has such a mode.
fn annulus_special_modes(
    cx: &Cx<'_>,
    nu: NodeId,
    ends: [NodeId; 6],
    robin: bool,
    sphere: bool,
) -> Option<(bool, u32)> {
    if !robin {
        return Some((false, 0));
    }
    let values: Option<Vec<f64>> = ends.iter().map(|&e| exact_value(cx, e)).collect();
    let Some(v) = values else {
        return Some((false, 0));
    };
    let [pa, qa, a, pb, qb, b] = <[f64; 6]>::try_from(v).ok()?;
    let probe = |order: f64| {
        let all = [order, pa, qa, a, pb, qb, b];
        (super::numeric::annulus_has_zero(all), super::numeric::annulus_neg_count(all))
    };
    match exact_value(cx, nu) {
        | Some(order) => Some(probe(order)),
        | None => {
            let shell = if sphere { 0.5 } else { 0.0 };
            for k in 0..=12_u32 {
                let (zero, negatives) = probe(f64::from(k) + shell);
                if zero || negatives > 0 {
                    return None;
                }
            }
            Some((false, 0))
        },
    }
}

/// The zero mode `A r^ν + B r^{-ν}` (`A + B ln r` for `ν = 0`, times
/// `r^{-1/2}` for a shell) meeting the inner condition by construction.
fn annulus_zero_mode(
    cx: &mut Cx<'_>,
    r: NodeId,
    nu: NodeId,
    ends: [NodeId; 6],
    sphere: bool,
) -> Option<Family> {
    let [pa, qa, a, _pb, _qb, b] = ends;
    let zero = cx.graph.int(0);
    let profile = if cx.graph.number_of(nu).is_some_and(Number::is_zero) {
        // A = -(p ln a + q/a), B = p.
        let ln_a = call(cx, "ln", &[a])?;
        let ln_r = call(cx, "ln", &[r])?;
        let t1 = mul(cx.graph, &[pa, ln_a]);
        let t2 = div(cx, qa, a);
        let sum = add(cx.graph, &[t1, t2]);
        let coefficient_a = neg(cx.graph, sum);
        let b_term = mul(cx.graph, &[pa, ln_r]);
        add(cx.graph, &[coefficient_a, b_term])
    } else {
        // A = -(p a^{-ν} - q ν a^{-ν-1}), B = p a^ν + q ν a^{ν-1}.
        let minus_nu = neg(cx.graph, nu);
        let a_minus = pow(cx.graph, a, minus_nu);
        let a_plus = pow(cx.graph, a, nu);
        let one = cx.graph.int(1);
        let neg_one = neg(cx.graph, one);
        let e_minus = add(cx.graph, &[minus_nu, neg_one]);
        let e_plus = sub(cx.graph, nu, one);
        let (a_m1, a_p1) = (pow(cx.graph, a, e_minus), pow(cx.graph, a, e_plus));
        let coefficient_a = {
            let first = mul(cx.graph, &[pa, a_minus]);
            let second = mul(cx.graph, &[qa, nu, a_m1]);
            let diff = sub(cx.graph, first, second);
            neg(cx.graph, diff)
        };
        let coefficient_b = {
            let first = mul(cx.graph, &[pa, a_plus]);
            let second = mul(cx.graph, &[qa, nu, a_p1]);
            add(cx.graph, &[first, second])
        };
        let (r_plus, r_minus) = (pow(cx.graph, r, nu), pow(cx.graph, r, minus_nu));
        let t1 = mul(cx.graph, &[coefficient_a, r_plus]);
        let t2 = mul(cx.graph, &[coefficient_b, r_minus]);
        add(cx.graph, &[t1, t2])
    };
    let mode = if sphere {
        let minus_half = fraction(cx, -1, 2)?;
        let root = pow(cx.graph, r, minus_half);
        mul(cx.graph, &[root, profile])
    } else {
        profile
    };
    // ∫ r F² dr over the annulus (the shell factor r^{-1} r² cancels).
    let square = powi(cx.graph, profile, 2);
    let body = mul(cx.graph, &[r, square]);
    let integral = super::util::defint(cx, body, r, a, b)?;
    let norm = cx.simplify(integral);
    Some(Family::new(cx.simplify(mode), zero, 0, 1, Some(norm)))
}

/// The `m`-th mode `F = C_K I_ν(κ r) - C_I K_ν(κ r)` (times `r^{-1/2}` for a
/// shell) of an annulus with eigenvalue `-κ²`, `κ = annulus_neg_root(..)`,
/// meeting the inner condition by construction.
fn annulus_negative_mode(
    cx: &mut Cx<'_>,
    r: NodeId,
    nu: NodeId,
    ends: [NodeId; 6],
    sphere: bool,
    m: i64,
) -> Option<Family> {
    let [pa, qa, a, pb, qb, b] = ends;
    let index = cx.graph.int(m);
    let kappa = call(cx, "annulus_neg_root", &[nu, pa, qa, a, pb, qb, b, index])?;
    let one = cx.graph.int(1);
    let nu_up = add(cx.graph, &[nu, one]);
    let modified = |cx: &mut Cx<'_>, second: bool, order: NodeId, z: NodeId| call(cx, if second { "besselk" } else { "besseli" }, &[order, z]);
    // `Z'(z)`: `I_ν' = I_{ν+1} + (ν/z) I_ν`, `K_ν' = -K_{ν+1} + (ν/z) K_ν`.
    let slope = |cx: &mut Cx<'_>, second: bool, z: NodeId| -> Option<NodeId> {
        let (z_nu, z_up) = (modified(cx, second, nu, z)?, modified(cx, second, nu_up, z)?);
        let ratio = div(cx, nu, z);
        let a_term = mul(cx.graph, &[ratio, z_nu]);
        let b_term = if second { neg(cx.graph, z_up) } else { z_up };
        Some(add(cx.graph, &[a_term, b_term]))
    };
    // `p Z(κ a) + q κ Z'(κ a)`.
    let ka = mul(cx.graph, &[kappa, a]);
    let condition = |cx: &mut Cx<'_>, second: bool| -> Option<NodeId> {
        let value = modified(cx, second, nu, ka)?;
        let derivative = slope(cx, second, ka)?;
        let t1 = mul(cx.graph, &[pa, value]);
        let t2 = mul(cx.graph, &[qa, kappa, derivative]);
        Some(add(cx.graph, &[t1, t2]))
    };
    let (cond_i, cond_k) = (condition(cx, false)?, condition(cx, true)?);
    let kr = mul(cx.graph, &[kappa, r]);
    let (i_r, k_r) = (modified(cx, false, nu, kr)?, modified(cx, true, nu, kr)?);
    let first = mul(cx.graph, &[cond_k, i_r]);
    let second = mul(cx.graph, &[cond_i, k_r]);
    let minus_second = neg(cx.graph, second);
    let profile = add(cx.graph, &[first, minus_second]);
    let mode = if sphere {
        let minus_half = fraction(cx, -1, 2)?;
        let root = pow(cx.graph, r, minus_half);
        mul(cx.graph, &[root, profile])
    } else {
        profile
    };
    // ∫ r F² dr = [r²/2 ((1 + ν²/z²) F² - F_z²)]_a^b.
    let nu2 = powi(cx.graph, nu, 2);
    let at = |cx: &mut Cx<'_>, point: NodeId| -> Option<NodeId> {
        let z = mul(cx.graph, &[kappa, point]);
        let (i_z, k_z) = (modified(cx, false, nu, z)?, modified(cx, true, nu, z)?);
        let (di, dk) = (slope(cx, false, z)?, slope(cx, true, z)?);
        let f = {
            let (x1, x2) = (mul(cx.graph, &[cond_k, i_z]), mul(cx.graph, &[cond_i, k_z]));
            sub(cx.graph, x1, x2)
        };
        let fz = {
            let (x1, x2) = (mul(cx.graph, &[cond_k, di]), mul(cx.graph, &[cond_i, dk]));
            sub(cx.graph, x1, x2)
        };
        let z2 = powi(cx.graph, z, 2);
        let ratio = div(cx, nu2, z2);
        let factor = add(cx.graph, &[one, ratio]);
        let f2 = powi(cx.graph, f, 2);
        let lead = mul(cx.graph, &[factor, f2]);
        let fz2 = powi(cx.graph, fz, 2);
        let tail = neg(cx.graph, fz2);
        let bracket = add(cx.graph, &[lead, tail]);
        let p2 = powi(cx.graph, point, 2);
        let half = fraction(cx, 1, 2)?;
        Some(mul(cx.graph, &[half, p2, bracket]))
    };
    let (at_b, at_a) = (at(cx, b)?, at(cx, a)?);
    let norm = sub(cx.graph, at_b, at_a);
    let k2 = powi(cx.graph, kappa, 2);
    let eigen = neg(cx.graph, k2);
    Some(Family::new(cx.simplify(mode), cx.simplify(eigen), 0, m, Some(norm)))
}

/// The modes of an annulus `a < r < b` for the Bessel order `nu` (the
/// disk: the angular index; the ball: `l + 1/2`), `sphere` for a shell.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(super) fn annulus_families(
    cx: &mut Cx<'_>,
    r: NodeId,
    inner: (NodeId, &End),
    outer: (NodeId, &End),
    sphere: bool,
    order: NodeId,
    n: NodeId,
    numeric: bool,
) -> Option<Vec<Family>> {
    let (a, end_a) = inner;
    let (b, end_b) = outer;
    let zero = cx.graph.int(0);
    let order_value = exact_value(cx, order);
    // Shells with l = 0: v = r u is a function on an interval.
    if sphere && order_value.is_some_and(|v| v.abs() < 1e-12) {
        // v = r u: p u + q u' = (p - q/r) v / r + q v' / r.
        let shift = |cx: &mut Cx<'_>, end: &End, point: NodeId| {
            let ratio = div(cx, end.q, point);
            let minus = neg(cx.graph, ratio);
            let p = add(cx.graph, &[end.p, minus]);
            cx.simplify(p)
        };
        let (pa, pb) = (shift(cx, end_a, a), shift(cx, end_b, b));
        let make = |cx: &mut Cx<'_>, end: &End, point: NodeId, p: NodeId| End { kind: kind_of(cx, p, end.q), point, p, q: end.q, value: zero };
        let left = make(cx, end_a, a, pa);
        let right = make(cx, end_b, b, pb);
        let length = sub(cx.graph, b, a);
        let length = cx.simplify(length);
        let iv = Interval { var: r, start: a, stop: b, length, left, right, periodic: false };
        let families = interval_families(cx, &iv, n, numeric)?;
        let mut out = Vec::new();
        for f in families {
            let mode = div(cx, f.mode, r);
            out.push(Family { mode, ..f });
        }
        return Some(out);
    }
    let half = fraction(cx, 1, 2)?;
    let nu = if sphere { add(cx.graph, &[order, half]) } else { order };
    let nu = cx.simplify(nu);
    let neumann_like = is_zero_number(cx.graph, end_a.p) && is_zero_number(cx.graph, end_b.p);
    let constant_mode = neumann_like && !sphere && order_value.is_some_and(|v| v.abs() < 1e-12);
    if numeric && constant_mode && cx.graph.number_of(n).is_some_and(Number::is_zero) {
        let one = cx.graph.int(1);
        let (a2, b2) = (powi(cx.graph, a, 2), powi(cx.graph, b, 2));
        let gap = sub(cx.graph, b2, a2);
        let norm = mul(cx.graph, &[half, gap]);
        return Some(vec![Family::new(one, zero, 0, 1, Some(cx.simplify(norm)))]);
    }
    let (pa, pb) = (shifted_p(cx, end_a, a, sphere), shifted_p(cx, end_b, b, sphere));
    let ends = [pa, end_a.q, a, pb, end_b.q, b];
    let robin = end_a.kind == Kind::Robin || end_b.kind == Kind::Robin;
    let (zero_mode, negatives) = annulus_special_modes(cx, nu, ends, robin, sphere)?;
    let specials = i64::from(negatives) + i64::from(zero_mode);
    let literal = if numeric { cx.graph.number_of(n).and_then(Number::to_i64) } else { None };
    if let Some(m) = literal.filter(|&m| m < specials) {
        if zero_mode && m == 0 {
            return annulus_zero_mode(cx, r, nu, ends, sphere).map(|f| vec![f]);
        }
        return annulus_negative_mode(cx, r, nu, ends, sphere, m + 1 - i64::from(zero_mode)).map(|f| vec![f]);
    }
    // Positive roots are numbered from 1; the indices below `specials` are
    // the modes with a zero or negative eigenvalue.
    let root_index = if specials > 1 {
        let shift = cx.graph.int(1 - specials);
        let shifted = add(cx.graph, &[n, shift]);
        cx.simplify(shifted)
    } else {
        n
    };
    let k = call(cx, "annulus_root", &[nu, pa, end_a.q, a, pb, end_b.q, b, root_index])?;
    let eigen = powi(cx.graph, k, 2);
    let ka = mul(cx.graph, &[k, a]);
    // A = (p + q ν/a) Y_ν(ka) - q k Y_{ν+1}(ka), B likewise with J.
    let one = cx.graph.int(1);
    let nu_up = add(cx.graph, &[nu, one]);
    let coefficient = |cx: &mut Cx<'_>, second: bool| -> Option<NodeId> {
        let z_nu = bessel(cx, second, nu, ka)?;
        let z_up = bessel(cx, second, nu_up, ka)?;
        let ratio = div(cx, nu, a);
        let qr = mul(cx.graph, &[end_a.q, ratio]);
        let lead = add(cx.graph, &[pa, qr]);
        let first = mul(cx.graph, &[lead, z_nu]);
        let second_term = mul(cx.graph, &[end_a.q, k, z_up]);
        Some(sub(cx.graph, first, second_term))
    };
    let (coef_a, coef_b) = (coefficient(cx, true)?, coefficient(cx, false)?);
    let kr = mul(cx.graph, &[k, r]);
    let (j, y) = (bessel(cx, false, nu, kr)?, bessel(cx, true, nu, kr)?);
    let (aj, by) = (mul(cx.graph, &[coef_a, j]), mul(cx.graph, &[coef_b, y]));
    let cross = sub(cx.graph, aj, by);
    let minus_b = neg(cx.graph, coef_b);
    let (mode, power2) = if sphere {
        let minus_half = fraction(cx, -1, 2)?;
        let root = pow(cx.graph, r, minus_half);
        (mul(cx.graph, &[root, cross]), -1)
    } else {
        (cross, 0)
    };
    // ∫ r F² dr = [r²/2 (F_z² + (1 - ν²/(k r)²) F²)]_a^b in the argument z = k r.
    let nu2 = powi(cx.graph, nu, 2);
    let k2 = powi(cx.graph, k, 2);
    let nu2_over_k2 = div(cx, nu2, k2);
    let at = |cx: &mut Cx<'_>, point: NodeId| -> Option<NodeId> {
        let z = mul(cx.graph, &[k, point]);
        let (j, y) = (bessel(cx, false, nu, z)?, bessel(cx, true, nu, z)?);
        let (jp, yp) = (bessel_slope(cx, false, nu, z)?, bessel_slope(cx, true, nu, z)?);
        let f = {
            let (x1, x2) = (mul(cx.graph, &[coef_a, j]), mul(cx.graph, &[coef_b, y]));
            sub(cx.graph, x1, x2)
        };
        let fp = {
            let (x1, x2) = (mul(cx.graph, &[coef_a, jp]), mul(cx.graph, &[coef_b, yp]));
            sub(cx.graph, x1, x2)
        };
        let p2 = powi(cx.graph, point, 2);
        let t1 = {
            let sq = powi(cx.graph, fp, 2);
            mul(cx.graph, &[p2, sq])
        };
        let t2 = {
            let gap = sub(cx.graph, p2, nu2_over_k2);
            let sq = powi(cx.graph, f, 2);
            mul(cx.graph, &[gap, sq])
        };
        Some(add(cx.graph, &[t1, t2]))
    };
    let (at_b, at_a) = (at(cx, b)?, at(cx, a)?);
    let gap = sub(cx.graph, at_b, at_a);
    let norm = mul(cx.graph, &[half, gap]);
    let first = i64::from(!constant_mode && specials == 0);
    let forced = if specials > 0 { specials } else { i64::from(constant_mode) };
    let cylinder = Cylinder { nu, kappa: k, power2, parts: vec![(coef_a, false), (minus_b, true)] };
    Some(vec![Family { mode, eigen, first, forced, norm: Some(norm), cylinder: Some(cylinder) }])
}

// ----------------------------------------------------------------------
// Legendre functions
// ----------------------------------------------------------------------

/// `P_l` as a polynomial in the generator `g`.
fn legendre_poly(
    g: u32,
    l: i64,
) -> Option<Poly> {
    let x = Poly::generator(g);
    let mut previous = Poly::constant(Number::from(1));
    if l == 0 {
        return Some(previous);
    }
    let mut current = x.clone();
    for k in 1..l {
        // (k + 1) P_{k+1} = (2k + 1) x P_k - k P_{k-1}.
        let a = x.mul(&current, 64)?.scale(&Number::from(2 * k + 1));
        let b = previous.scale(&Number::from(k));
        let next = a.sub(&b).scale(&Number::fraction(1, k + 1)?);
        previous = current;
        current = next;
    }
    Some(current)
}

/// `sin^m θ · P_l^{(m)}(cos θ)`: the associated Legendre function of `cos
/// θ` (without the Condon–Shortley sign).
pub(super) fn associated_mode(
    cx: &mut Cx<'_>,
    theta: NodeId,
    l: i64,
    m: i64,
) -> Option<NodeId> {
    let c = call(cx, "cos", &[theta])?;
    let mut gens = Gens::default();
    let g = gens.index(cx.graph, c);
    let mut poly = legendre_poly(g, l)?;
    for _ in 0..m {
        poly = poly.derivative(g);
    }
    let q = to_term(cx.graph, &gens, &poly);
    let s = call(cx, "sin", &[theta])?;
    let sm = powi(cx.graph, s, m);
    let mode = mul(cx.graph, &[sm, q]);
    Some(cx.simplify(mode))
}

/// `∫_{-1}^{1} P_l^m(x)² dx = 2/(2l+1) (l+m)!/(l-m)!`.
pub(super) fn polar_norm(
    cx: &mut Cx<'_>,
    l: i64,
    m: i64,
) -> Option<NodeId> {
    let mut ratio = Number::fraction(2, 2 * l + 1)?;
    for k in (l - m + 1)..=(l + m) {
        ratio = ratio.mul(&Number::from(k));
    }
    Some(cx.graph.num(ratio))
}

/// `C(n, k)` for small arguments.
fn binomial(
    n: u32,
    k: u32,
) -> i64 {
    let mut value: i64 = 1;
    for i in 0..k {
        value = value * i64::from(n - i) / i64::from(i + 1);
    }
    value
}

/// `∫_{-1}^{1} F dc` for `F` a polynomial in `c` (the node `c`) and, when
/// `s` is given, in an even power of `s` (with `s² = 1 - c²`, the substitution
/// `x = cos θ`); `None` if the integrand has another structure, or depends on
/// `var` otherwise.
pub(super) fn legendre_integral(
    cx: &mut Cx<'_>,
    integrand: NodeId,
    c: NodeId,
    s: Option<NodeId>,
    var: NodeId,
) -> Option<NodeId> {
    let mut gens = Gens::default();
    let poly = from_term(cx.graph, &mut gens, integrand, Limits { terms: 4096, exponent: 40 })?;
    let gc = gens.find(cx.graph, c);
    let gs = s.and_then(|s| gens.find(cx.graph, s));
    let mut result = Poly::zero();
    for (mono, coeff) in poly.terms() {
        let exponent_of = |g: Option<u32>| g.and_then(|g| mono.iter().find(|&&(h, _)| h == g)).map_or(0, |&(_, e)| e);
        let (ec, es) = (exponent_of(gc), exponent_of(gs));
        if es % 2 != 0 {
            return None;
        }
        let rest: Vec<(u32, u32)> = mono.iter().copied().filter(|&(h, _)| Some(h) != gc && Some(h) != gs).collect();
        let half = es / 2;
        let mut number = Number::from(0);
        for j in 0..=half {
            let power = ec + 2 * j;
            if power % 2 != 0 {
                continue;
            }
            let integral = Number::fraction(2, i64::from(power) + 1)?;
            let sign = if j % 2 == 0 { 1 } else { -1 };
            number = number.add(&integral.mul(&Number::from(sign * binomial(half, j))));
        }
        result = result.add(&Poly::monomial(rest, coeff.mul(&number)));
    }
    let out = to_term(cx.graph, &gens, &result);
    let symbol = cx.graph.symbol_of(var)?;
    if cx.graph.depends_on(cx.graph.find(out), symbol) {
        return None;
    }
    Some(cx.simplify(out))
}

// ----------------------------------------------------------------------
// Fourier–Bessel moments
// ----------------------------------------------------------------------

/// The node for the half-integer `v2 / 2`.
fn halves(
    cx: &mut Cx<'_>,
    v2: i64,
) -> Option<NodeId> {
    if v2 % 2 == 0 { Some(cx.graph.int(v2 / 2)) } else { fraction(cx, v2, 2) }
}

/// An antiderivative `G(z)` of `z^μ Z_ν(z)` (`μ = mu2/2`, `ν = nu2/2`) in
/// terms of Bessel functions, if a recurrence terminates.
fn bessel_antiderivative(
    cx: &mut Cx<'_>,
    second: bool,
    mu2: i64,
    nu2: i64,
    z: NodeId,
) -> Option<NodeId> {
    if let Some(found) = terminating_antiderivative(cx, second, mu2, nu2, z) {
        return Some(found);
    }
    if mu2 % 2 == 0 && nu2 % 2 == 0 && mu2 >= 0 && nu2 >= 0 && !(mu2 == 0 && (nu2 / 2) % 2 == 1) {
        return struve_antiderivative(cx, second, mu2 / 2, nu2 / 2, z);
    }
    None
}

/// An antiderivative of `z^μ Z_ν(z)` for integers `μ, ν ≥ 0` in terms of
/// Bessel and Struve functions: the order is lowered by
/// `∫ z^μ Z_ν = -z^μ Z_{ν-1} + (μ+ν-1) ∫ z^{μ-1} Z_{ν-1}` to zero, the
/// power by Lommel's `∫ z^μ Z_0 = z^μ Z_1 + (μ-1) z^{μ-1} Z_0 - (μ-1)² ∫ z^{μ-2} Z_0`,
/// and the base is `∫ Z_0 = z Z_0 + (π z/2)(Z_1 H_0 - Z_0 H_1)`.
fn struve_antiderivative(
    cx: &mut Cx<'_>,
    second: bool,
    mu: i64,
    nu: i64,
    z: NodeId,
) -> Option<NodeId> {
    let order = |cx: &mut Cx<'_>, k: i64| -> Option<NodeId> {
        let k = cx.graph.int(k);
        bessel(cx, second, k, z)
    };
    let power = |cx: &mut Cx<'_>, e: i64| powi(cx.graph, z, e);
    if nu == 0 {
        if mu == 0 {
            let (z0, z1) = (order(cx, 0)?, order(cx, 1)?);
            let (zero, one) = (cx.graph.int(0), cx.graph.int(1));
            let (h0, h1) = (call(cx, "struve_h", &[zero, z])?, call(cx, "struve_h", &[one, z])?);
            let pi = super::util::pi(cx)?;
            let half = fraction(cx, 1, 2)?;
            let a = mul(cx.graph, &[z, z0]);
            let b = mul(cx.graph, &[z1, h0]);
            let c = mul(cx.graph, &[z0, h1]);
            let nc = neg(cx.graph, c);
            let bracket = add(cx.graph, &[b, nc]);
            let d = mul(cx.graph, &[half, pi, z, bracket]);
            return Some(add(cx.graph, &[a, d]));
        }
        if mu == 1 {
            let z1 = order(cx, 1)?;
            return Some(mul(cx.graph, &[z, z1]));
        }
        let (z0, z1) = (order(cx, 0)?, order(cx, 1)?);
        let pm = power(cx, mu);
        let a = mul(cx.graph, &[pm, z1]);
        let c1 = cx.graph.int(mu - 1);
        let pm1 = power(cx, mu - 1);
        let b = mul(cx.graph, &[c1, pm1, z0]);
        let c2 = cx.graph.int(-(mu - 1) * (mu - 1));
        let rest = struve_antiderivative(cx, second, mu - 2, 0, z)?;
        let c = mul(cx.graph, &[c2, rest]);
        return Some(add(cx.graph, &[a, b, c]));
    }
    if mu == 0 {
        // Z_ν = Z_{ν-2} - 2 Z'_{ν-1}.
        if nu == 1 {
            let z0 = order(cx, 0)?;
            return Some(neg(cx.graph, z0));
        }
        let lower = struve_antiderivative(cx, second, 0, nu - 2, z)?;
        let top = order(cx, nu - 1)?;
        let m2 = cx.graph.int(-2);
        let b = mul(cx.graph, &[m2, top]);
        return Some(add(cx.graph, &[lower, b]));
    }
    let z_low = order(cx, nu - 1)?;
    let a = {
        let pm = power(cx, mu);
        let m = neg(cx.graph, pm);
        mul(cx.graph, &[m, z_low])
    };
    let factor = cx.graph.int(mu + nu - 1);
    let rest = struve_antiderivative(cx, second, mu - 1, nu - 1, z)?;
    let b = mul(cx.graph, &[factor, rest]);
    Some(add(cx.graph, &[a, b]))
}

fn terminating_antiderivative(
    cx: &mut Cx<'_>,
    second: bool,
    mu2: i64,
    nu2: i64,
    z: NodeId,
) -> Option<NodeId> {
    let up = |mu2: i64, nu2: i64| mu2 - nu2 - 2;
    let down = |mu2: i64, nu2: i64| mu2 + nu2 - 2;
    // d2 = 2 (μ ∓ ν - 1) must be a non-negative multiple of 4.
    let (use_up, mut d2) = if up(mu2, nu2) >= 0 && up(mu2, nu2) % 4 == 0 {
        (true, up(mu2, nu2))
    } else if down(mu2, nu2) >= 0 && down(mu2, nu2) % 4 == 0 {
        (false, down(mu2, nu2))
    } else {
        return None;
    };
    let (mut mu2, mut nu2) = (mu2, nu2);
    let mut factor = Number::from(1);
    let mut terms = Vec::new();
    loop {
        let order2 = if use_up { nu2 + 2 } else { nu2 - 2 };
        let order = halves(cx, order2)?;
        let z_order = bessel(cx, second, order, z)?;
        let power = {
            let mu = halves(cx, mu2)?;
            pow(cx.graph, z, mu)
        };
        let sign = Number::from(if use_up { 1 } else { -1 });
        let coefficient = cx.graph.num(factor.mul(&sign));
        terms.push(mul(cx.graph, &[coefficient, power, z_order]));
        // The next factor: -(μ - ν - 1) or +(μ + ν - 1), i.e. ∓ d2/2.
        let d = Number::fraction(d2, 2)?;
        if d.is_zero() {
            break;
        }
        factor = factor.mul(&if use_up { d.neg() } else { d });
        mu2 -= 2;
        nu2 += if use_up { 2 } else { -2 };
        d2 -= 4;
    }
    let total = add(cx.graph, &terms);
    Some(total)
}

/// `∫_{lo}^{hi} h(r) w(r) mode(r) dr` for a cylinder mode and polynomial
/// `h` in closed form (`weight2` the exponent of `w = r^{weight2/2}`
/// in halves).
pub(super) fn cylinder_integral(
    cx: &mut Cx<'_>,
    r: NodeId,
    lo: NodeId,
    hi: NodeId,
    weight2: i64,
    cylinder: &Cylinder,
    h: NodeId,
) -> Option<NodeId> {
    let nu_value = exact_value(cx, cylinder.nu)?;
    let nu2_float = 2.0 * nu_value;
    if (nu2_float - nu2_float.round()).abs() > 1e-9 || nu2_float.abs() > 1e6 {
        return None;
    }
    #[allow(clippy::cast_possible_truncation)]
    let nu2 = nu2_float.round() as i64;
    let r_symbol = cx.graph.symbol_of(r)?;
    // h as a polynomial in r.
    let coefficients: Vec<NodeId> = if cx.graph.depends_on(cx.graph.find(h), r_symbol) {
        let mut gens = Gens::default();
        let poly = from_term(cx.graph, &mut gens, h, Limits { terms: 4096, exponent: 40 })?;
        let g = gens.find(cx.graph, r)?;
        let parts = poly.coefficients_in(g);
        let mut out = Vec::new();
        for part in &parts {
            let term = to_term(cx.graph, &gens, part);
            if cx.graph.depends_on(cx.graph.find(term), r_symbol) {
                return None;
            }
            out.push(term);
        }
        out
    } else {
        vec![h]
    };
    let lo_is_zero = is_zero_number(cx.graph, lo);
    let mut total = Vec::new();
    for (j, &coefficient) in coefficients.iter().enumerate() {
        if cx.is_zero(coefficient) {
            continue;
        }
        let j = i64::try_from(j).ok()?;
        let mu2 = 2 * j + weight2 + cylinder.power2;
        let mut parts = Vec::new();
        for &(c, second) in &cylinder.parts {
            let z_hi = mul(cx.graph, &[cylinder.kappa, hi]);
            let upper = bessel_antiderivative(cx, second, mu2, nu2, z_hi)?;
            let value = if lo_is_zero {
                if second {
                    return None;
                }
                upper
            } else {
                let z_lo = mul(cx.graph, &[cylinder.kappa, lo]);
                let lower = bessel_antiderivative(cx, second, mu2, nu2, z_lo)?;
                sub(cx.graph, upper, lower)
            };
            parts.push(mul(cx.graph, &[c, value]));
        }
        let combined = add(cx.graph, &parts);
        let scale = {
            let exponent = halves(cx, -(mu2 + 2))?;
            pow(cx.graph, cylinder.kappa, exponent)
        };
        total.push(mul(cx.graph, &[coefficient, scale, combined]));
    }
    let sum = add(cx.graph, &total);
    Some(cx.simplify(sum))
}

// ----------------------------------------------------------------------
// Exterior profiles
// ----------------------------------------------------------------------

/// The radial profile (`r^n`, `r^l` inside when `interior`) of a decaying (`k² = 0`), outgoing (`k² > 0`:
/// Hankel functions of the first kind, the Sommerfeld condition) or
/// evanescent (`k² < 0`: Macdonald functions) solution of `Δu + k² u = 0`
/// outside a disk or ball for the angular index `index`, normalised so
/// that `p ρ(R) + q ρ'(R) = 1`.
#[allow(clippy::too_many_lines, clippy::too_many_arguments)]
pub(super) fn exterior_profile(
    cx: &mut Cx<'_>,
    r: NodeId,
    radius: NodeId,
    end: &End,
    sphere: bool,
    interior: bool,
    wavenumber2: NodeId,
    index: NodeId,
) -> Option<NodeId> {
    if interior && !cx.is_zero(wavenumber2) {
        return None;
    }
    let sign = if cx.is_zero(wavenumber2) {
        0
    } else {
        let v = super::util::sample(cx.graph, wavenumber2, 0)?;
        if v > 0.0 { 1 } else { -1 }
    };
    let literal = exact_value(cx, index).filter(|v| v.fract() == 0.0 && *v >= 0.0 && *v <= 40.0);
    let rho = match (sign, sphere) {
        | (0, false) => powi_node(cx, r, index, if interior { 1 } else { -1 }),
        | (0, true) => {
            let e = if interior {
                index
            } else {
                let one = cx.graph.int(1);
                let e = add(cx.graph, &[index, one]);
                neg(cx.graph, e)
            };
            pow(cx.graph, r, e)
        },
        | (_, false) => {
            let kk = if sign > 0 { wavenumber2 } else { neg(cx.graph, wavenumber2) };
            let k = sqrt(cx, kk)?;
            let kr = mul(cx.graph, &[k, r]);
            if sign > 0 {
                let (j, y) = (bessel(cx, false, index, kr)?, bessel(cx, true, index, kr)?);
                let i = imaginary(cx)?;
                let iy = mul(cx.graph, &[i, y]);
                add(cx.graph, &[j, iy])
            } else {
                call(cx, "besselk", &[index, kr])?
            }
        },
        | (_, true) => {
            let l = literal?;
            let l = i64::from(l as i32);
            let kk = if sign > 0 { wavenumber2 } else { neg(cx.graph, wavenumber2) };
            let k = sqrt(cx, kk)?;
            let kr = mul(cx.graph, &[k, r]);
            let two_kr = {
                let two = cx.graph.int(2);
                mul(cx.graph, &[two, kr])
            };
            // e^{±ikr}/r Σ_m (i or 1)^m (l+m)!/(m!(l-m)!) (2kr)^{-m}.
            let i = imaginary(cx)?;
            let mut series = Vec::new();
            for m in 0..=l {
                let mut ratio = Number::from(1);
                for q in 1..=m {
                    ratio = ratio.mul(&Number::fraction(l + q, q)?).mul(&Number::fraction(l - q + 1, 1)?);
                }
                let coefficient = cx.graph.num(ratio);
                let inverse = powi(cx.graph, two_kr, -m);
                let phase = if sign > 0 { powi(cx.graph, i, m) } else { cx.graph.int(1) };
                series.push(mul(cx.graph, &[coefficient, phase, inverse]));
            }
            let sum = add(cx.graph, &series);
            let exponent = if sign > 0 {
                mul(cx.graph, &[i, kr])
            } else {
                neg(cx.graph, kr)
            };
            let wave = call(cx, "exp", &[exponent])?;
            let inverse_r = powi(cx.graph, r, -1);
            mul(cx.graph, &[wave, inverse_r, sum])
        },
    };
    // Normalisation by p ρ(R) + q ρ'(R).
    let slope = derivative(cx.graph, rho, r)?;
    let at = |cx: &mut Cx<'_>, e: NodeId| {
        let s = cx.graph.substitute(e, r, radius);
        cx.simplify(s)
    };
    let (rho_r, slope_r) = (at(cx, rho), at(cx, slope));
    let a = mul(cx.graph, &[end.p, rho_r]);
    let b = mul(cx.graph, &[end.q, slope_r]);
    let boundary = add(cx.graph, &[a, b]);
    let boundary = cx.simplify(boundary);
    if cx.is_zero(boundary) {
        return None;
    }
    let profile = div(cx, rho, boundary);
    Some(cx.simplify(profile))
}

/// `r^(sign · index)`.
fn powi_node(
    cx: &mut Cx<'_>,
    r: NodeId,
    index: NodeId,
    sign: i64,
) -> NodeId {
    let s = cx.graph.int(sign);
    let e = mul(cx.graph, &[s, index]);
    pow(cx.graph, r, e)
}
