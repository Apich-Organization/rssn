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
