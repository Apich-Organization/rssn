//! Periodic Sturm–Liouville problems with a smooth potential by Hill's
//! method (Mathieu, Hill and Lamé-type equations).
//!
//! `A u_xx + c(x) u` on a period `L` (`u(a) = u(a + L)`) is reduced to
//! `H y = -y'' + V y = μ y`, `V = -c/A`, `L y = -A μ y`. In the real
//! Fourier basis `{1, cos(2πkx/L), sin(2πkx/L)}` (normalised on the period)
//! the operator is a symmetric matrix — the Hill determinant is its
//! characteristic polynomial — whose eigenpairs are computed by cyclic
//! Jacobi rotations:
//!
//! * `hill_eigen(m, L, N, V₀…V_{M-1})`: the `m`-th eigenvalue (`m ≥ 1`,
//!   ascending, with multiplicity) for the potential sampled at
//!   `x_j = L j/M`; `N` is the number of harmonics (the matrix has
//!   `2N + 1` rows);
//! * `hill_mode(m, ξ, L, N, V…)`: the orthonormal eigenfunction at the
//!   offset `ξ` from the start of the period.
//!
//! For `V = 2 q cos 2x` these are the characteristic values of Mathieu's
//! equation `a_0, b_1, a_1, b_2, a_2, …` and the functions `ce_n`, `se_n`.

use super::Conditions;
use super::Problem;
use super::spectral::Axis;
use super::spectral::Engine;
use super::spectral::Item;
use super::spectral::Role;
use super::spectral::Shape;
use super::spectral::Symbol;
use super::spectral::Time;
use super::spectral::evolution_time;
use super::spectral::finish;
use super::spectral::interval;
use super::spectral::time_items;
use super::util::is_zero_number;
use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::Env;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpDescriptor;
use crate::graph::RuleError;
use crate::graph::rule::Installer;
use std::cell::RefCell;
use std::f64::consts::PI;
use std::rc::Rc;

const HARMONICS: usize = 24;
const SAMPLES: usize = 256;

/// Eigenpairs of a symmetric matrix by cyclic Jacobi rotations, sorted
/// ascending: `(values, vectors)` with `vectors[k]` the eigenvector of
/// `values[k]`.
pub(super) fn jacobi(mut a: Vec<Vec<f64>>) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = a.len();
    let mut v: Vec<Vec<f64>> = (0..n).map(|i| (0..n).map(|j| if i == j { 1.0 } else { 0.0 }).collect()).collect();
    for _ in 0..60 {
        let off: f64 = (0..n).flat_map(|i| (0..n).filter(move |&j| j != i).map(move |j| (i, j))).map(|(i, j)| a[i][j] * a[i][j]).sum();
        if off < 1e-26 {
            break;
        }
        for p in 0..n {
            for q in (p + 1)..n {
                if a[p][q].abs() < 1e-300 {
                    continue;
                }
                let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                let t = theta.signum() / (theta.abs() + theta.mul_add(theta, 1.0).sqrt());
                let t = if theta == 0.0 { 1.0 } else { t };
                let c = 1.0 / t.mul_add(t, 1.0).sqrt();
                let s = t * c;
                for k in 0..n {
                    let (akp, akq) = (a[k][p], a[k][q]);
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for k in 0..n {
                    let (apk, aqk) = (a[p][k], a[q][k]);
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
                for row in &mut v {
                    let (vkp, vkq) = (row[p], row[q]);
                    row[p] = c * vkp - s * vkq;
                    row[q] = s * vkp + c * vkq;
                }
            }
        }
    }
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| a[i][i].total_cmp(&a[j][j]));
    let values = order.iter().map(|&i| a[i][i]).collect();
    let vectors = order.iter().map(|&i| (0..n).map(|k| v[k][i]).collect()).collect();
    (values, vectors)
}

/// The value of the real Fourier basis function number `i` at `x`.
fn basis(
    i: usize,
    x: f64,
    l: f64,
) -> f64 {
    if i == 0 {
        return (1.0 / l).sqrt();
    }
    let k = i.div_ceil(2) as f64;
    let angle = 2.0 * PI * k * x / l;
    (2.0 / l).sqrt() * if i % 2 == 1 { angle.cos() } else { angle.sin() }
}

struct Solved {
    l: f64,
    values: Vec<f64>,
    vectors: Vec<Vec<f64>>,
}

type Cached = (Vec<u64>, Rc<Solved>);

thread_local! {
    static CACHE: RefCell<Vec<Cached>> = const { RefCell::new(Vec::new()) };
}

/// The decomposition for `L, N, V…`.
fn solved(args: &[f64]) -> Option<Rc<Solved>> {
    let key: Vec<u64> = args.iter().map(|v| v.to_bits()).collect();
    if let Some(found) = CACHE.with(|c| c.borrow().iter().find(|(k, _)| *k == key).map(|(_, v)| Rc::clone(v))) {
        return Some(found);
    }
    let (&l, &harmonics) = (args.first()?, args.get(1)?);
    let potential = args.get(2..)?;
    if l <= 0.0 || !(1.0..=64.0).contains(&harmonics) || potential.len() < 8 {
        return None;
    }
    let d = 2 * harmonics as usize + 1;
    let m = potential.len();
    let mut h = vec![vec![0.0; d]; d];
    for (i, row) in h.iter_mut().enumerate() {
        for (j, entry) in row.iter_mut().enumerate() {
            let mut total = 0.0;
            for (index, &v) in potential.iter().enumerate() {
                let x = l * index as f64 / m as f64;
                total += basis(i, x, l) * v * basis(j, x, l);
            }
            *entry = total * l / m as f64;
        }
        let k = i.div_ceil(2) as f64;
        row[i] += (2.0 * PI * k / l).powi(2);
    }
    let (values, vectors) = jacobi(h);
    let value = Rc::new(Solved { l, values, vectors });
    CACHE.with(|c| {
        let mut cache = c.borrow_mut();
        if cache.len() >= 32 {
            cache.remove(0);
        }
        cache.push((key, Rc::clone(&value)));
    });
    Some(value)
}

#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss, clippy::float_cmp)]
fn index(m: f64) -> Option<usize> {
    if (1.0..=1e4).contains(&m) && m.fract() == 0.0 { Some(m as usize) } else { None }
}

/// `hill_eigen(m, L, N, V…)`.
pub(super) fn hill_eigen(args: &[f64]) -> f64 {
    let (Some(m), Some(rest)) = (args.first().and_then(|&m| index(m)), args.get(1..)) else {
        return f64::NAN;
    };
    solved(rest).and_then(|s| s.values.get(m - 1).copied()).unwrap_or(f64::NAN)
}

/// `hill_mode(m, ξ, L, N, V…)`.
pub(super) fn hill_mode(args: &[f64]) -> f64 {
    let (Some(m), Some(&x), Some(rest)) = (args.first().and_then(|&m| index(m)), args.get(1), args.get(2..)) else {
        return f64::NAN;
    };
    let Some(s) = solved(rest) else {
        return f64::NAN;
    };
    let Some(vector) = s.vectors.get(m - 1) else {
        return f64::NAN;
    };
    vector.iter().enumerate().map(|(i, c)| c * basis(i, x, s.l)).sum()
}

/// Registers the operators.
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    i.op(OpDescriptor::new("hill_eigen", Arity::Variadic).eval(hill_eigen))?;
    i.op(OpDescriptor::new("hill_mode", Arity::Variadic).eval(hill_mode))?;
    Ok(())
}

/// Separation of variables for `A u_xx + c(x) u` on a periodic interval.
pub(super) fn solve(
    cx: &mut Cx<'_>,
    p: &Problem,
    conditions: &Conditions,
) -> Option<NodeId> {
    if p.nonlinear || p.dimension() != 2 {
        return None;
    }
    let time_index = evolution_time(cx, p)?;
    let j = 1 - time_index;
    let x = p.vars[j];
    let n = 2;
    let (ca, cb) = (p.coefficient(cx.graph, &p.unit(j, 2)), p.coefficient(cx.graph, &p.unit(j, 1)));
    let cc = p.coefficient(cx.graph, &vec![0; n]);
    if cx.is_zero(ca) || !p.constant(cx.graph, ca) || !cx.is_zero(cb) {
        return None;
    }
    let (a1, a2) = (p.coefficient(cx.graph, &p.unit(time_index, 1)), p.coefficient(cx.graph, &p.unit(time_index, 2)));
    if !p.constant(cx.graph, a1) || !p.constant(cx.graph, a2) {
        return None;
    }
    for (index, c) in &p.linear {
        let known = *index == p.unit(j, 2) || *index == vec![0; n] || *index == p.unit(time_index, 1) || *index == p.unit(time_index, 2);
        if !known && !is_zero_number(cx.graph, *c) {
            return None;
        }
    }
    let iv = interval(cx, p, conditions, j)?;
    if !iv.periodic {
        return None;
    }
    let value = |cx: &Cx<'_>, node: NodeId| super::modes::exact_value(cx, node);
    let (start, stop) = (value(cx, iv.start)?, value(cx, iv.stop)?);
    let a_value = value(cx, ca)?;
    let length = stop - start;
    if length <= 0.0 {
        return None;
    }
    let symbol = cx.graph.symbol_of(x)?;
    let mut potential = Vec::with_capacity(SAMPLES);
    for m in 0..SAMPLES {
        let mut env = Env::numeric(0.0);
        env.bind(symbol, start + length * m as f64 / SAMPLES as f64);
        let c = cx.graph.eval(cc, &env).filter(|v| v.is_finite())?;
        potential.push(-c / a_value);
    }
    let number = |cx: &mut Cx<'_>, v: f64| cx.graph.num(Number::from(v));
    let length_node = number(cx, length);
    let harmonics = cx.graph.int(i64::try_from(HARMONICS).ok()?);
    let mut params = vec![length_node, harmonics];
    for v in potential {
        params.push(number(cx, v));
    }
    let order = if is_zero_number(cx.graph, a2) { u32::from(!is_zero_number(cx.graph, a1)) } else { 2 };
    let time = Time { var: p.vars[time_index], order, a2, a1 };
    let axes = vec![Axis { var: x, shape: Shape::Hill { start: iv.start, stop: iv.stop, params } }];
    let mut items = time_items(cx, p, conditions, Some(time_index), order, n)?;
    if !p.homogeneous(cx.graph) {
        items.push(Item { role: Role::Source, h: p.source, special: None });
    }
    let engine = Engine { axes, time: Some(time), symbol: Symbol::Cartesian(vec![(vec![1], ca)]), face: None, exterior: None };
    finish(cx, p, conditions, &engine, &items, None, None)
}
