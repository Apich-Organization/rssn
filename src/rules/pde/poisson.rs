//! Poisson integrals: the harmonic function in a disk or ball with
//! arbitrary boundary values, in closed integral form.
//!
//! `pde_poisson_integral(g, R, 2)`: `g(phi)` on the circle `r = R`,
//! `u(r, theta) = (R² - r²)/(2π) ∫_0^{2π} g(φ) / (R² - 2 R r cos(θ - φ) + r²) dφ`.
//!
//! `pde_poisson_integral(g, R, 3)`: `g(theta, phi)` on the sphere,
//! `u(r, theta, phi) = R (R² - r²)/(4π) ∬ g sin θ' / (R² + r² - 2 R r cos γ)^{3/2} dθ' dφ'`
//! with `cos γ = cos θ cos θ' + sin θ sin θ' cos(φ - φ')`.

use super::util::call;
use super::util::defint;
use super::util::div;
use super::util::pi;
use super::util::sqrt;
use crate::graph::Cx;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::rules::complex::build::add;
use crate::rules::complex::build::mul;
use crate::rules::complex::build::neg;
use crate::rules::complex::build::pow;
use crate::rules::complex::build::powi;
use crate::rules::complex::build::sub;

/// The Poisson integral of the boundary datum `g` on the circle (`dimension`
/// 2) or the sphere (`dimension` 3) of radius `radius`.
pub(super) fn poisson_integral(
    cx: &mut Cx<'_>,
    g: NodeId,
    radius: NodeId,
    dimension: i64,
) -> Option<NodeId> {
    let (r, theta, phi) = (cx.graph.sym("r"), cx.graph.sym("theta"), cx.graph.sym("phi"));
    let (zero, two, pi) = (cx.graph.int(0), cx.graph.int(2), pi(cx)?);
    let two_pi = mul(cx.graph, &[two, pi]);
    let r2 = powi(cx.graph, radius, 2);
    let rho2 = powi(cx.graph, r, 2);
    let gap = sub(cx.graph, r2, rho2);
    match dimension {
        | 2 => {
            let alpha = cx.graph.sym("alpha");
            let datum = cx.graph.substitute(g, phi, alpha);
            let angle = sub(cx.graph, theta, alpha);
            let cosine = call(cx, "cos", &[angle])?;
            let cross = mul(cx.graph, &[two, radius, r, cosine]);
            let denominator = {
                let s = add(cx.graph, &[r2, rho2]);
                sub(cx.graph, s, cross)
            };
            let body = div(cx, datum, denominator);
            let integral = defint(cx, body, alpha, zero, two_pi)?;
            let prefactor = div(cx, gap, two_pi);
            let u = mul(cx.graph, &[prefactor, integral]);
            Some(cx.simplify(u))
        },
        | 3 => {
            let (alpha, beta) = (cx.graph.sym("alpha"), cx.graph.sym("beta"));
            let datum = {
                let first = cx.graph.substitute(g, theta, alpha);
                cx.graph.substitute(first, phi, beta)
            };
            let (sin_a, cos_a) = (call(cx, "sin", &[alpha])?, call(cx, "cos", &[alpha])?);
            let (sin_t, cos_t) = (call(cx, "sin", &[theta])?, call(cx, "cos", &[theta])?);
            let azimuth = {
                let d = sub(cx.graph, phi, beta);
                call(cx, "cos", &[d])?
            };
            let cos_gamma = {
                let polar = mul(cx.graph, &[cos_t, cos_a]);
                let equatorial = mul(cx.graph, &[sin_t, sin_a, azimuth]);
                add(cx.graph, &[polar, equatorial])
            };
            let cross = mul(cx.graph, &[two, radius, r, cos_gamma]);
            let base = {
                let s = add(cx.graph, &[r2, rho2]);
                sub(cx.graph, s, cross)
            };
            let root = sqrt(cx, base)?;
            let cube = powi(cx.graph, root, 3);
            let numerator = mul(cx.graph, &[datum, sin_a]);
            let body = div(cx, numerator, cube);
            let inner = defint(cx, body, beta, zero, two_pi)?;
            let outer = defint(cx, inner, alpha, zero, pi)?;
            let four_pi = {
                let four = cx.graph.num(Number::from(4));
                mul(cx.graph, &[four, pi])
            };
            let scaled = mul(cx.graph, &[gap, radius]);
            let prefactor = div(cx, scaled, four_pi);
            let u = mul(cx.graph, &[prefactor, outer]);
            Some(cx.simplify(u))
        },
        | _ => None,
    }
}

/// `pde_greens_ball(R, d)`: the Dirichlet Green's function of `-Δ` for the
/// disk (`d = 2`) or ball (`d = 3`) of radius `R` and, with the same formula,
/// for the exterior (`|x|, |y| > R`), in the variables `x, y[, z]` and the
/// source point `s_x, s_y[, s_z]`. With `ρ² = |x|²|s|²/R² - 2 x·s + R²`
/// (the distance from `x` to the Kelvin image of `s`, scaled by `|s|/R`):
///
/// * `d = 3`: `G = 1/(4π|x - s|) - 1/(4π ρ)`,
/// * `d = 2`: `G = ln(ρ²/|x - s|²)/(4π)`.
pub(super) fn greens_ball(
    cx: &mut Cx<'_>,
    radius: NodeId,
    dimension: i64,
) -> Option<NodeId> {
    if !(2..=3).contains(&dimension) {
        return None;
    }
    let names: &[(&str, &str)] = if dimension == 2 { &[("x", "s_x"), ("y", "s_y")] } else { &[("x", "s_x"), ("y", "s_y"), ("z", "s_z")] };
    let mut x_squared = Vec::new();
    let mut s_squared = Vec::new();
    let mut cross = Vec::new();
    let mut gaps = Vec::new();
    for &(a, b) in names {
        let (xa, sa) = (cx.graph.sym(a), cx.graph.sym(b));
        x_squared.push(powi(cx.graph, xa, 2));
        s_squared.push(powi(cx.graph, sa, 2));
        cross.push(mul(cx.graph, &[xa, sa]));
        let gap = sub(cx.graph, xa, sa);
        gaps.push(powi(cx.graph, gap, 2));
    }
    let (xx, ss, xs, gap2) = (add(cx.graph, &x_squared), add(cx.graph, &s_squared), add(cx.graph, &cross), add(cx.graph, &gaps));
    let r2 = powi(cx.graph, radius, 2);
    let rho2 = {
        let product = mul(cx.graph, &[xx, ss]);
        let scaled = div(cx, product, r2);
        let two = cx.graph.int(2);
        let double = mul(cx.graph, &[two, xs]);
        let t = sub(cx.graph, scaled, double);
        add(cx.graph, &[t, r2])
    };
    let pi = super::util::pi(cx)?;
    let four_pi = {
        let four = cx.graph.int(4);
        mul(cx.graph, &[four, pi])
    };
    let g = if dimension == 3 {
        let (d, rho) = (sqrt(cx, gap2)?, sqrt(cx, rho2)?);
        let one = cx.graph.int(1);
        let free = div(cx, one, d);
        let image = div(cx, one, rho);
        let difference = sub(cx.graph, free, image);
        div(cx, difference, four_pi)
    } else {
        let ratio = div(cx, rho2, gap2);
        let log = call(cx, "ln", &[ratio])?;
        div(cx, log, four_pi)
    };
    Some(cx.simplify(g))
}

/// The fundamental solution of `-Δ` in `d` dimensions at squared distance `q`.
fn fundamental_laplace(
    cx: &mut Cx<'_>,
    dimension: i64,
    q: NodeId,
) -> Option<NodeId> {
    let pi = super::util::pi(cx)?;
    if dimension == 3 {
        let four_pi = {
            let four = cx.graph.int(4);
            mul(cx.graph, &[four, pi])
        };
        let root = sqrt(cx, q)?;
        let denominator = mul(cx.graph, &[four_pi, root]);
        let one = cx.graph.int(1);
        Some(div(cx, one, denominator))
    } else {
        let four_pi = {
            let four = cx.graph.int(4);
            mul(cx.graph, &[four, pi])
        };
        let log = call(cx, "ln", &[q])?;
        let value = div(cx, log, four_pi);
        Some(neg(cx.graph, value))
    }
}

/// `pde_greens_robin(d, h)`: the Green's function of `-Δ` in the half-space
/// whose last coordinate (`y` for `d = 2`, `z` for `d = 3`) is positive, with
/// the Robin condition `∂_n u + h u = 0` (`∂_n = -∂_z`, `h ≥ 0`; `h = 0` is the
/// Neumann problem, the Green's function up to the usual constant):
///
/// `G = G₀(x - s) + G₀(x - s*) - 2h ∫_0^∞ e^{-hw} G₀(x_∥ - s_∥, z + s_z + w) dw`.
pub(super) fn greens_robin(
    cx: &mut Cx<'_>,
    dimension: i64,
    h: NodeId,
) -> Option<NodeId> {
    if !(2..=3).contains(&dimension) {
        return None;
    }
    let names: &[(&str, &str)] = if dimension == 2 { &[("x", "s_x"), ("y", "s_y")] } else { &[("x", "s_x"), ("y", "s_y"), ("z", "s_z")] };
    let last = names.len() - 1;
    let w = cx.graph.sym("w");
    let mut direct = Vec::new();
    let mut image = Vec::new();
    let mut integral = Vec::new();
    for (k, &(a, b)) in names.iter().enumerate() {
        let (xa, sa) = (cx.graph.sym(a), cx.graph.sym(b));
        if k == last {
            let minus = sub(cx.graph, xa, sa);
            direct.push(powi(cx.graph, minus, 2));
            let plus = add(cx.graph, &[xa, sa]);
            image.push(powi(cx.graph, plus, 2));
            let shifted = add(cx.graph, &[xa, sa, w]);
            integral.push(powi(cx.graph, shifted, 2));
        } else {
            let gap = sub(cx.graph, xa, sa);
            let square = powi(cx.graph, gap, 2);
            direct.push(square);
            image.push(square);
            integral.push(square);
        }
    }
    let (q_direct, q_image, q_integral) = (add(cx.graph, &direct), add(cx.graph, &image), add(cx.graph, &integral));
    let first = fundamental_laplace(cx, dimension, q_direct)?;
    let second = fundamental_laplace(cx, dimension, q_image)?;
    let mut terms = vec![first, second];
    if !cx.is_zero(h) {
        let g0 = fundamental_laplace(cx, dimension, q_integral)?;
        let decay = {
            let e = mul(cx.graph, &[h, w]);
            let e = neg(cx.graph, e);
            call(cx, "exp", &[e])?
        };
        let body = mul(cx.graph, &[decay, g0]);
        let zero = cx.graph.int(0);
        let oo = super::util::infinity(cx)?;
        let integral = defint(cx, body, w, zero, oo)?;
        let minus_two = cx.graph.int(-2);
        terms.push(mul(cx.graph, &[minus_two, h, integral]));
    }
    let g = add(cx.graph, &terms);
    Some(cx.simplify(g))
}

/// `pde_greens_ball_robin(R, h)`: the interior Green's function of `-Δ` in
/// the ball of radius `R` with `∂_r u + h u = 0` on `r = R` (`h > 0`):
///
/// `G = 1/(4π|x - s|) + Σ_l A_l (r s)^l P_l(cos γ)`,
/// `A_l = -(h - (l+1)/R) R^{-(2l+1)} / (4π (h + l/R))`.
pub(super) fn greens_ball_robin(
    cx: &mut Cx<'_>,
    radius: NodeId,
    h: NodeId,
) -> Option<NodeId> {
    let names = [("x", "s_x"), ("y", "s_y"), ("z", "s_z")];
    let mut gaps = Vec::new();
    let (mut rr, mut ss, mut dot) = (Vec::new(), Vec::new(), Vec::new());
    for &(a, b) in &names {
        let (xa, sa) = (cx.graph.sym(a), cx.graph.sym(b));
        let gap = sub(cx.graph, xa, sa);
        gaps.push(powi(cx.graph, gap, 2));
        rr.push(powi(cx.graph, xa, 2));
        ss.push(powi(cx.graph, sa, 2));
        dot.push(mul(cx.graph, &[xa, sa]));
    }
    let (q, xx, sss, xs) = (add(cx.graph, &gaps), add(cx.graph, &rr), add(cx.graph, &ss), add(cx.graph, &dot));
    let free = fundamental_laplace(cx, 3, q)?;
    let (r, s) = (sqrt(cx, xx)?, sqrt(cx, sss)?);
    let rs = mul(cx.graph, &[r, s]);
    let cosine = div(cx, xs, rs);
    let l = cx.graph.sym("l");
    let one = cx.graph.int(1);
    let l1 = add(cx.graph, &[l, one]);
    let pi = super::util::pi(cx)?;
    let four_pi = {
        let four = cx.graph.int(4);
        mul(cx.graph, &[four, pi])
    };
    // A_l (r s)^l P_l(cos γ).
    let numerator = {
        let ratio = div(cx, l1, radius);
        let bracket = sub(cx.graph, h, ratio);
        neg(cx.graph, bracket)
    };
    let power = {
        let two_l1 = {
            let two = cx.graph.int(2);
            let tl = mul(cx.graph, &[two, l]);
            add(cx.graph, &[tl, one])
        };
        let p = pow(cx.graph, radius, two_l1);
        div(cx, one, p)
    };
    let denominator = {
        let ratio = div(cx, l, radius);
        let bracket = add(cx.graph, &[h, ratio]);
        mul(cx.graph, &[four_pi, bracket])
    };
    let coefficient = {
        let top = mul(cx.graph, &[numerator, power]);
        div(cx, top, denominator)
    };
    let rs_l = pow(cx.graph, rs, l);
    let legendre = call(cx, "legendre", &[l, cosine])?;
    let body = mul(cx.graph, &[coefficient, rs_l, legendre]);
    let zero = cx.graph.int(0);
    let oo = super::util::infinity(cx)?;
    let series = call(cx, "sum", &[body, l, zero, oo])?;
    let g = add(cx.graph, &[free, series]);
    Some(g)
}
