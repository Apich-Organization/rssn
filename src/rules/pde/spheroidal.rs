//! Helmholtz and Laplace problems in prolate and oblate spheroidal
//! coordinates.
//!
//! The equation, in coordinates `(ξ, η, φ)` (or `(ζ, η, φ)` for the oblate
//! system, `(ξ, η)` without azimuthal dependence), multiplied by the common
//! scale factor,
//!
//! * prolate: `(ξ²-1) u_ξξ + 2ξ u_ξ + (1-η²) u_ηη - 2η u_η`, plus
//!   `(ξ²-η²)/((ξ²-1)(1-η²)) u_φφ`, plus `c² (ξ²-η²) u`, equals zero,
//! * oblate:  `(ζ²+1) u_ζζ + 2ζ u_ζ + (1-η²) u_ηη - 2η u_η`, plus
//!   `(ζ²+η²)/((ζ²+1)(1-η²)) u_φφ`, plus `c² (ζ²+η²) u`, equals zero,
//!
//! is recognised from its coefficients (any common factor, `c²` constant,
//! either sign). It separates into the angular equation
//! `(1-η²) S'' - 2η S' + (λ - γ η² - m²/(1-η²)) S = 0` with `γ = c²`
//! (prolate) or `γ = -c²` (oblate), solved by `spheroidal_eigen` /
//! `spheroidal_mode`, and the radial equation
//!
//! * prolate: `(ξ²-1) R'' + 2ξ R' - (λ - γ ξ² + m²/(ξ²-1)) R = 0`, regular at
//!   the focal segment `ξ = 1` (`R ~ (ξ²-1)^{m/2}`), integrated in
//!   `ξ = cosh τ`;
//! * oblate: `(ζ²+1) R'' + 2ζ R' - (λ + γ ζ² - m²/(ζ²+1)) R = 0` with the
//!   parity of `S` at the focal disk `ζ = 0`, integrated in `ζ = sinh σ`.
//!
//! The regular radial tables grow through the largest requested coordinate.
//! They share the fixed Runge–Kutta geometry and retain each finite prefix.
//!
//! The Dirichlet problem inside the spheroid `ξ < ξ₀` (`ζ < ζ₀`) with data
//! `u(ξ₀, η, φ) = g(η, φ)` is solved by the finite expansion
//! `u = Σ (a_mn cos mφ + b_mn sin mφ) S_mn(η) R_mn(ξ)/R_mn(ξ₀)` whose
//! coefficients are computed by quadrature; the expansion is returned when
//! it reproduces the data at sample points (so the data must be well
//! represented by the first 16 angular functions of every order `m ≤ 6`).
//! Exterior regions (radiation condition) are not covered.

use super::Conditions;
use super::Problem;
use super::coulomb::spheroidal_pair;
use super::util::call;
use super::util::div;
use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::DimRule;
use crate::graph::Env;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpDescriptor;
use crate::graph::RuleError;
use crate::graph::rule::Installer;
use crate::kernels::reuse::PlanCache;
use crate::rules::complex::build::add;
use crate::rules::complex::build::mul;
use std::cell::RefCell;
use std::f64::consts::PI;
use std::sync::Arc;
use std::sync::OnceLock;
use std::sync::RwLock;

const STEP: f64 = 0.002;
const T_MAX: f64 = 4.5;
const RADIAL_STEPS: usize = (T_MAX / STEP) as usize;
const GEOMETRY_BLOCK: usize = 64;
const GEOMETRY_BLOCKS: usize = RADIAL_STEPS.div_ceil(GEOMETRY_BLOCK);
const MODES: usize = 16;
const M_MAX: usize = 6;

#[derive(Clone, Copy)]
struct RadialStage {
    gamma_coordinate: f64,
    order_denominator: f64,
    damping: f64,
}

type RadialStep = [RadialStage; 3];
type GeometrySlot = OnceLock<[RadialStep; GEOMETRY_BLOCK]>;

#[allow(
    clippy::single_call_fn,
    reason = "shares the fixed RK4 coordinate geometry across radial plans"
)]
fn step_geometry(
    prolate: bool,
    index: usize,
) -> &'static RadialStep {
    static PROLATE: [GeometrySlot; GEOMETRY_BLOCKS] = [const { OnceLock::new() }; GEOMETRY_BLOCKS];
    static OBLATE: [GeometrySlot; GEOMETRY_BLOCKS] = [const { OnceLock::new() }; GEOMETRY_BLOCKS];
    let blocks = if prolate { &PROLATE } else { &OBLATE };
    let block = index / GEOMETRY_BLOCK;
    &blocks[block].get_or_init(|| {
        std::array::from_fn(|offset| {
            let t = STEP * (block * GEOMETRY_BLOCK + offset) as f64;
            [t, t + 0.5 * STEP, t + STEP].map(|t| {
                let (s, c) = (t.sinh(), t.cosh());
                if prolate {
                    RadialStage {
                        gamma_coordinate: c,
                        order_denominator: s * s,
                        damping: 1.0 / t.tanh(),
                    }
                } else {
                    RadialStage {
                        gamma_coordinate: s,
                        order_denominator: c * c,
                        damping: t.tanh(),
                    }
                }
            })
        })
    })[index % GEOMETRY_BLOCK]
}

struct RadialPath {
    prolate: bool,
    eigenvalue: f64,
    gamma: f64,
    order_squared: f64,
    values: Vec<f64>,
    derivative: f64,
    stopped: bool,
}

type RadialKey = (bool, usize, usize, u64);
type RadialTable = Arc<RwLock<RadialPath>>;
static RADIAL: PlanCache<RadialKey, RwLock<RadialPath>> = PlanCache::new(300);
const READ_WINDOW: usize = 32;

struct RadialWindow {
    start: usize,
    len: usize,
    values: [f64; READ_WINDOW],
}

type LocalRadial = (RadialKey, RadialTable, RadialWindow);

thread_local! {
    static LAST_RADIAL: RefCell<Option<LocalRadial>> = const { RefCell::new(None) };
}

/// The regular radial function on the uniform grid `τ_i = i STEP`
/// (`ξ = cosh τ` prolate, `ζ = sinh τ` oblate).
#[allow(
    clippy::single_call_fn,
    reason = "keeps radial plan initialization separate from progressive reads"
)]
fn radial_table(key: RadialKey) -> Option<RadialTable> {
    let (prolate, m, n, gamma_bits) = key;
    let gamma = f64::from_bits(gamma_bits);
    RADIAL
        .try_get_or_insert_with(key, || {
            let offset = n.checked_sub(m).ok_or(())?;
            let lambda = *spheroidal_pair(m, gamma).0.get(offset).ok_or(())?;
            if !lambda.is_finite() {
                return Err(());
            }
            let mf = m as f64;
            let (mut values, derivative) = if prolate {
                // Series at the focal segment: R = τ^m (1 + a τ²), ξ = cosh τ.
                let a = ((lambda - gamma) - mf * (mf + 1.0) / 3.0) / (4.0 * mf + 4.0);
                let start = 10usize;
                let tau0 = STEP * start as f64;
                let values = std::iter::once(if m == 0 { 1.0 } else { 0.0 })
                    .chain((1..=start).map(|i| {
                        let tau = STEP * i as f64;
                        tau.powi(m as i32) * (1.0 + a * tau * tau)
                    }))
                    .collect();
                let derivative = tau0.powi(m as i32 - 1)
                    * (mf * (1.0 + a * tau0 * tau0) + 2.0 * a * tau0 * tau0);
                (values, derivative)
            } else if offset.is_multiple_of(2) {
                (vec![1.0], 0.0)
            } else {
                (vec![0.0], 1.0)
            };
            let invalid = values.iter().position(|value| !value.is_finite());
            if let Some(invalid) = invalid {
                values.truncate(invalid);
            }
            Ok(Arc::new(RwLock::new(RadialPath {
                prolate,
                eigenvalue: lambda,
                gamma,
                order_squared: mf * mf,
                values,
                derivative,
                stopped: invalid.is_some() || !derivative.is_finite(),
            })))
        })
        .ok()
}

impl RadialPath {
    #[allow(
        clippy::single_call_fn,
        reason = "keeps radial RK4 progression and finite prefixes together"
    )]
    fn extend(
        &mut self,
        end: usize,
    ) {
        if self.stopped || self.values.len() > end {
            return;
        }
        let Some(&last) = self.values.last() else {
            return;
        };
        let (lambda, gamma, order, prolate) = (
            self.eigenvalue,
            self.gamma,
            self.order_squared,
            self.prolate,
        );
        let f = |stage: RadialStage, y: f64, p: f64| {
            let coordinate = stage.gamma_coordinate;
            let potential = if prolate {
                lambda - gamma * coordinate * coordinate + order / stage.order_denominator
            } else {
                lambda + gamma * coordinate * coordinate - order / stage.order_denominator
            };
            potential * y - stage.damping * p
        };
        let (mut y, mut p) = (last, self.derivative);
        let h = STEP;
        self.values.reserve(end + 1 - self.values.len());
        for i in self.values.len() - 1..end {
            let &[first, middle, last] = step_geometry(prolate, i);
            let (k1y, k1p) = (p, f(first, y, p));
            let (k2y, k2p) = (
                p + 0.5 * h * k1p,
                f(middle, y + 0.5 * h * k1y, p + 0.5 * h * k1p),
            );
            let (k3y, k3p) = (
                p + 0.5 * h * k2p,
                f(middle, y + 0.5 * h * k2y, p + 0.5 * h * k2p),
            );
            let (k4y, k4p) = (p + h * k3p, f(last, y + h * k3y, p + h * k3p));
            y += h / 6.0 * (k1y + 2.0 * k2y + 2.0 * k3y + k4y);
            p += h / 6.0 * (k1p + 2.0 * k2p + 2.0 * k3p + k4p);
            if !y.is_finite() {
                self.stopped = true;
                return;
            }
            self.values.push(y);
            self.derivative = p;
            if !p.is_finite() {
                self.stopped = true;
                return;
            }
        }
    }
}

impl RadialWindow {
    fn refill<const N: usize>(
        &mut self,
        path: &RadialPath,
        start: usize,
    ) -> Option<[f64; N]> {
        let values = path.values.get(start..)?;
        let samples = *values.first_chunk()?;
        self.start = start;
        self.len = READ_WINDOW.min(values.len());
        self.values[..self.len].copy_from_slice(&values[..self.len]);
        Some(samples)
    }
}

fn radial_samples<const N: usize>(
    key: RadialKey,
    start: usize,
) -> Option<[f64; N]> {
    LAST_RADIAL.with(|entry| {
        let mut entry = entry.borrow_mut();
        if entry.as_ref().is_none_or(|(cached, _, _)| *cached != key) {
            *entry = Some((
                key,
                radial_table(key)?,
                RadialWindow {
                    start: 0,
                    len: 0,
                    values: [0.0; READ_WINDOW],
                },
            ));
        }
        let (_, table, window) = entry.as_mut()?;
        if let Some(offset) = start.checked_sub(window.start)
            && let Some(samples) = window.values[..window.len]
                .get(offset..)
                .and_then(<[f64]>::first_chunk)
        {
            return Some(*samples);
        }
        let samples = window.refill(
            &table
                .read()
                .unwrap_or_else(std::sync::PoisonError::into_inner),
            start,
        );
        if samples.is_some() {
            return samples;
        }
        let mut path = table
            .write()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        path.extend(start + N - 1);
        let samples = window.refill(&path, start);
        drop(path);
        samples
    })
}

/// The regular radial function at `x` (`ξ ≥ 1` prolate, `ζ ≥ 0` oblate).
pub(super) fn radial_value(
    prolate: bool,
    m: usize,
    n: usize,
    gamma: f64,
    x: f64,
) -> f64 {
    if !x.is_finite() || !gamma.is_finite() {
        return f64::NAN;
    }
    let t = if prolate {
        if x < 1.0 {
            return f64::NAN;
        }
        x.acosh()
    } else {
        if x < 0.0 {
            return f64::NAN;
        }
        x.asinh()
    };
    let position = t / STEP;
    let limit = RADIAL_STEPS - 2;
    if position >= limit as f64 {
        return f64::NAN;
    }
    let key = (prolate, m, n, gamma.to_bits());
    let floor = position.floor() as usize;
    if crate::kernels::exact_eq(position, 0.0) {
        return radial_samples::<1>(key, floor).map_or(f64::NAN, |[sample]| sample);
    }
    let i = floor.max(1);
    let u = position - i as f64;
    let Some([p0, p1, p2, p3]) = radial_samples(key, i - 1) else {
        if crate::kernels::exact_eq(u, 0.0) {
            return radial_samples::<1>(key, i).map_or(f64::NAN, |[sample]| sample);
        }
        return f64::NAN;
    };
    let a = -u * (u - 1.0) * (u - 2.0) / 6.0;
    let b = (u + 1.0) * (u - 1.0) * (u - 2.0) / 2.0;
    let c = -(u + 1.0) * u * (u - 2.0) / 2.0;
    let d = (u + 1.0) * u * (u - 1.0) / 6.0;
    a * p0 + b * p1 + c * p2 + d * p3
}

/// `spheroidal_radial(kind, m, n, γ, x)`, `kind = 1` prolate, `-1` oblate.
pub(super) fn spheroidal_radial(args: &[f64]) -> f64 {
    let (Some(&kind), Some(&m), Some(&n), Some(&gamma), Some(&x)) = (
        args.first(),
        args.get(1),
        args.get(2),
        args.get(3),
        args.get(4),
    ) else {
        return f64::NAN;
    };
    if !kind.is_finite()
        || !(0.0..=1e4).contains(&m)
        || !(m..=1e4).contains(&n)
        || m.fract() != 0.0
        || n.fract() != 0.0
    {
        return f64::NAN;
    }
    radial_value(kind > 0.0, m as usize, n as usize, gamma, x)
}

/// Registers the operator.
///
/// The radial spheroidal function is a normalised numeric value in the
/// dimensionless spheroidal coordinates `(ξ, ζ)`, so it is dimensionless.
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let op = i.op(
        OpDescriptor::new("spheroidal_radial", Arity::Fixed(5))
            .eval(spheroidal_radial)
            .doc(
                "`spheroidal_radial(kind, m, n, gamma, xi)`: prolate (`kind > 0`) or oblate \
                 radial spheroidal wave function at `xi`. The kind, coordinate and `gamma` must be finite, \
                 and the degree and order must be integers with `0 <= m <= n <= 10000`. \
                 The radial coordinate is at least one for prolate and at least zero for oblate.",
            ),
    )?;
    i.dim(op, DimRule::Dimensionless);
    Ok(())
}

/// The recognised structure of a spheroidal Helmholtz equation.
struct Structure {
    prolate: bool,
    /// The angular parameter `γ`.
    gamma: f64,
    radial: usize,
    eta: usize,
    azimuth: Option<usize>,
}

fn evaluate(
    cx: &Cx<'_>,
    node: NodeId,
    bindings: &[(crate::graph::SymbolId, f64)],
) -> Option<f64> {
    let mut env = Env::numeric(0.0);
    for &(s, v) in bindings {
        env.bind(s, v);
    }
    cx.graph.eval(node, &env).filter(|v| v.is_finite())
}

/// Recognises the equation from the numeric values of its coefficients.
fn recognise(
    cx: &mut Cx<'_>,
    p: &Problem,
) -> Option<Structure> {
    let n = p.dimension();
    if !(2..=3).contains(&n) || p.nonlinear || !p.homogeneous(cx.graph) {
        return None;
    }
    let name = |cx: &Cx<'_>, j: usize| {
        cx.graph
            .symbol_of(p.vars[j])
            .map(|s| cx.graph.interner().symbol_name(s).to_owned())
    };
    let find = |cx: &Cx<'_>, wanted: &[&str]| {
        (0..n).find(|&j| name(cx, j).is_some_and(|nm| wanted.contains(&nm.as_str())))
    };
    let radial = find(cx, &["xi", "zeta"])?;
    let eta = find(cx, &["eta", "nu"])?;
    let azimuth = find(cx, &["phi", "varphi", "ph"]);
    if n == 3 && azimuth.is_none() {
        return None;
    }
    let symbols: Vec<crate::graph::SymbolId> = p
        .vars
        .iter()
        .map(|&v| cx.graph.symbol_of(v))
        .collect::<Option<_>>()?;
    let mut allowed = vec![
        p.unit(radial, 2),
        p.unit(radial, 1),
        p.unit(eta, 2),
        p.unit(eta, 1),
        vec![0; n],
    ];
    if let Some(a) = azimuth {
        allowed.push(p.unit(a, 2));
    }
    if !p.only(cx.graph, &allowed) {
        return None;
    }
    let coefficients: Vec<(Vec<u32>, NodeId)> = allowed
        .iter()
        .map(|i| (i.clone(), p.coefficient(cx.graph, i)))
        .collect();
    for prolate in [true, false] {
        let mut scale: Option<f64> = None;
        let mut gamma: Option<f64> = None;
        let mut ok = true;
        'samples: for &(xi, eta_v) in &[(1.7, 0.3), (2.4, -0.6), (3.1, 0.8), (1.3, 0.15)] {
            let mut bindings = vec![(symbols[radial], xi), (symbols[eta], eta_v)];
            if let Some(a) = azimuth {
                bindings.push((symbols[a], 0.4));
            }
            let (xx, ee) = (xi * xi, eta_v * eta_v);
            let radial_factor = if prolate {
                xx - 1.0
            } else {
                xx + 1.0
            };
            let mix = if prolate { xx - ee } else { xx + ee };
            for (index, node) in &coefficients {
                let Some(value) = evaluate(cx, *node, &bindings) else {
                    ok = false;
                    break 'samples;
                };
                let expected = if *index == p.unit(radial, 2) {
                    radial_factor
                } else if *index == p.unit(radial, 1) {
                    2.0 * xi
                } else if *index == p.unit(eta, 2) {
                    1.0 - ee
                } else if *index == p.unit(eta, 1) {
                    -2.0 * eta_v
                } else if azimuth.is_some_and(|a| *index == p.unit(a, 2)) {
                    mix / (radial_factor * (1.0 - ee))
                } else {
                    // The zero-order term: g × mix.
                    let a = scale.unwrap_or(1.0);
                    let g = value / (a * mix);
                    match gamma {
                        | None => gamma = Some(g),
                        | Some(old) if (old - g).abs() > 1e-9 * (1.0 + g.abs()) => {
                            ok = false;
                            break 'samples;
                        },
                        | Some(_) => {},
                    }
                    continue;
                };
                if *index == p.unit(radial, 2) && scale.is_none() {
                    scale = Some(value / expected);
                }
                let a = scale?;
                if (value - a * expected).abs() > 1e-9 * (1.0 + value.abs()) {
                    ok = false;
                    break 'samples;
                }
            }
        }
        if ok && let (Some(_), Some(g)) = (scale, gamma) {
            return Some(Structure {
                prolate,
                gamma: if prolate { g } else { -g },
                radial,
                eta,
                azimuth,
            });
        }
    }
    None
}

/// The Dirichlet problem inside a spheroid by the finite expansion.
#[allow(clippy::too_many_lines)]
pub(super) fn solve(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
) -> Option<NodeId> {
    let structure = recognise(cx, p)?;
    let [boundary] = conditions.0.as_slice() else {
        return None;
    };
    if boundary.on != structure.radial
        || boundary.robin.is_some()
        || boundary.derivative.iter().any(|&d| d != 0)
    {
        return None;
    }
    let x0 = super::modes::exact_value(cx, boundary.point)?;
    if (structure.prolate && x0 <= 1.0) || (!structure.prolate && x0 <= 0.0) || x0 > 20.0 {
        return None;
    }
    let (eta_node, radial_node) = (p.vars[structure.eta], p.vars[structure.radial]);
    let eta_symbol = cx.graph.symbol_of(eta_node)?;
    let phi_symbol = structure
        .azimuth
        .and_then(|a| cx.graph.symbol_of(p.vars[a]));
    let data = boundary.value;
    let free = cx.graph.free_symbols(cx.graph.find(data)).to_vec();
    if free
        .iter()
        .any(|&s| s != eta_symbol && Some(s) != phi_symbol)
    {
        return None;
    }
    let g = |cx: &Cx<'_>, eta_v: f64, phi_v: f64| -> Option<f64> {
        let mut bindings = vec![(eta_symbol, eta_v)];
        if let Some(s) = phi_symbol {
            bindings.push((s, phi_v));
        }
        evaluate(cx, data, &bindings)
    };
    // Quadrature nodes: Gauss–Legendre in η, trapezoid in φ.
    let eta_nodes = super::curvilinear::gauss_legendre(48);
    let phi_count = if phi_symbol.is_some() {
        64
    } else {
        1
    };
    let phi_values: Vec<f64> = (0..phi_count)
        .map(|k| 2.0 * PI * k as f64 / phi_count as f64)
        .collect();
    let mut grid = vec![vec![0.0; phi_count]; eta_nodes.len()];
    let mut scale = 1e-12_f64;
    for (i, &(e, _)) in eta_nodes.iter().enumerate() {
        for (j, &ph) in phi_values.iter().enumerate() {
            grid[i][j] = g(cx, e, ph)?;
            scale = scale.max(grid[i][j].abs());
        }
    }
    let gamma = structure.gamma;
    let m_top = if phi_symbol.is_some() {
        M_MAX
    } else {
        0
    };
    // (m, n, cos coefficient, sin coefficient)
    let mut coefficients: Vec<(usize, usize, f64, f64)> = Vec::new();
    for m in 0..=m_top {
        let pair = spheroidal_pair(m, gamma);
        let vectors = &pair.1;
        for j in 0..MODES {
            let n = m + j;
            let Some(d) = vectors.get(j) else {
                continue;
            };
            let sign = if d[j] < 0.0 { -1.0 } else { 1.0 };
            let mut cos_sum = 0.0;
            let mut sin_sum = 0.0;
            for (i, &(e, w)) in eta_nodes.iter().enumerate() {
                let s = super::coulomb::spheroidal_mode_value(m, gamma, e, d, sign);
                for (jj, &ph) in phi_values.iter().enumerate() {
                    cos_sum += w * s * grid[i][jj] * (m as f64 * ph).cos();
                    sin_sum += w * s * grid[i][jj] * (m as f64 * ph).sin();
                }
            }
            let phi_norm = if m == 0 {
                phi_count as f64
            } else {
                phi_count as f64 / 2.0
            };
            let (a, b) = (
                cos_sum / phi_norm,
                if m == 0 {
                    0.0
                } else {
                    sin_sum / phi_norm
                },
            );
            if a.abs() > 1e-10 * scale || b.abs() > 1e-10 * scale {
                coefficients.push((m, n, a, b));
            }
        }
    }
    // The expansion must reproduce the data.
    for &(e, ph) in &[(0.31, 0.7), (-0.62, 2.1), (0.83, 4.4), (-0.17, 5.5)] {
        let mut total = 0.0;
        for &(m, n, a, b) in &coefficients {
            let pair = spheroidal_pair(m, gamma);
            let j = n - m;
            let d = pair.1.get(j)?;
            let sign = if d[j] < 0.0 { -1.0 } else { 1.0 };
            let s = super::coulomb::spheroidal_mode_value(m, gamma, e, d, sign);
            let mp = m as f64 * ph;
            total += s * (a * mp.cos() + b * mp.sin());
        }
        if (total - g(cx, e, ph)?).abs() > 1e-4 * scale {
            return None;
        }
    }
    // Assemble the symbolic expansion.
    let kind = cx.graph.int(if structure.prolate { 1 } else { -1 });
    let gamma_node = cx.graph.num(Number::from(gamma));
    let x0_node = boundary.point;
    let mut terms = Vec::new();
    for &(m, n, a, b) in &coefficients {
        let (m_node, n_node) = (
            cx.graph.int(i64::try_from(m).ok()?),
            cx.graph.int(i64::try_from(n).ok()?),
        );
        let angular = call(
            cx,
            "spheroidal_mode",
            &[m_node, n_node, gamma_node, p.vars[structure.eta]],
        )?;
        let at = call(
            cx,
            "spheroidal_radial",
            &[kind, m_node, n_node, gamma_node, radial_node],
        )?;
        let at0 = call(
            cx,
            "spheroidal_radial",
            &[kind, m_node, n_node, gamma_node, x0_node],
        )?;
        let ratio = div(cx, at, at0);
        let mut azimuthal = Vec::new();
        if let Some(index) = structure.azimuth {
            let phi = p.vars[index];
            let m_phi = mul(cx.graph, &[m_node, phi]);
            if a.abs() > 0.0 {
                let c = call(cx, "cos", &[m_phi])?;
                let a_node = cx.graph.num(Number::from(a));
                azimuthal.push(mul(cx.graph, &[a_node, c]));
            }
            if b.abs() > 0.0 {
                let s = call(cx, "sin", &[m_phi])?;
                let b_node = cx.graph.num(Number::from(b));
                azimuthal.push(mul(cx.graph, &[b_node, s]));
            }
        } else {
            azimuthal.push(cx.graph.num(Number::from(a)));
        }
        let trig = add(cx.graph, &azimuthal);
        terms.push(mul(cx.graph, &[trig, angular, ratio]));
    }
    if terms.is_empty() {
        return None;
    }
    let u = add(cx.graph, &terms);
    Some(cx.simplify(u))
}
