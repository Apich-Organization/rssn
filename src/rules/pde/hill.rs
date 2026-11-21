//! Periodic Sturm–Liouville problems with a smooth potential by Hill's
//! method (Mathieu, Hill and Lamé-type equations).
//!
//! `A u_xx + c(x) u` on a period `L` (`u(a) = u(a + L)`) is reduced to
//! `H y = -y'' + V y = μ y`, `V = -c/A`, `L y = -A μ y`. In the real
//! Fourier basis `{1, cos(2πkx/L), sin(2πkx/L)}` (normalised on the period)
//! the operator is a symmetric matrix. Fourier coefficients of the sampled
//! potential determine its entries through the product identities for sines
//! and cosines. Cyclic Jacobi rotations compute the eigenpairs once for use
//! by both the eigenvalue and mode queries.
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
use crate::graph::DimRule;
use crate::graph::Env;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpDescriptor;
use crate::graph::RuleError;
use crate::graph::rule::Installer;
use crate::kernels::reuse::PlanCache;
use crate::kernels::reuse::fft_plan;
use crate::kernels::reuse::fft_scratch_for;
use num_complex::Complex64;
use rustfft::FftDirection;
use std::cell::RefCell;
use std::f64::consts::SQRT_2;
use std::f64::consts::TAU;
use std::sync::Arc;

const HARMONICS: usize = 24;
const SAMPLES: usize = 256;
const MAX_HARMONICS: usize = 64;
const MAX_CACHED_SAMPLES: usize = 65_536;

/// Eigenpairs of a symmetric matrix by cyclic Jacobi rotations, sorted
/// ascending: `(values, vectors)` with `vectors[k]` the eigenvector of
/// `values[k]`.
pub(super) fn jacobi(mut a: Vec<Vec<f64>>) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = a.len();
    let mut v: Vec<Vec<f64>> = (0..n)
        .map(|i| (0..n).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
        .collect();
    for _ in 0..60 {
        let off: f64 = (0..n)
            .flat_map(|i| (0..n).filter(move |&j| j != i).map(move |j| (i, j)))
            .map(|(i, j)| a[i][j] * a[i][j])
            .sum();
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
    let vectors = order
        .iter()
        .map(|&i| (0..n).map(|k| v[k][i]).collect())
        .collect();
    (values, vectors)
}

/// Forms the real Galerkin matrix from the periodic coefficient table.
/// `RustFFT` uses a negative exponent, so the sine coefficient is `-im`.
#[allow(
    clippy::single_call_fn,
    reason = "keeps the Fourier product identities together"
)]
fn potential_entry(
    coefficients: &[Complex64],
    left: usize,
    right: usize,
) -> f64 {
    let (left, right) = (left.min(right), left.max(right));
    let frequency = right.div_ceil(2);
    if left == 0 {
        let entry = coefficients[frequency];
        return match right {
            | 0 => entry.re,
            | _ if right % 2 == 1 => SQRT_2 * entry.re,
            | _ => -SQRT_2 * entry.im,
        };
    }
    let left_frequency = left.div_ceil(2);
    let sum = coefficients[frequency + left_frequency];
    let difference = coefficients[frequency - left_frequency];
    match (left % 2, right % 2) {
        | (1, 1) => difference.re + sum.re,
        | (1, _) => -sum.im - difference.im,
        | (_, 1) => difference.im - sum.im,
        | _ => difference.re - sum.re,
    }
}

struct Solved {
    l: f64,
    values: Vec<f64>,
    vectors: Vec<Vec<f64>>,
}

type Cached = (Arc<[u64]>, Arc<Solved>);
static CACHE: PlanCache<Arc<[u64]>, Solved> = PlanCache::new(32);

thread_local! {
    // Repeated mode evaluations compare the key in place. The plan itself is
    // immutable and shared with other threads through the bounded pool.
    static LAST: RefCell<Option<Cached>> = const { RefCell::new(None) };
}

#[allow(
    clippy::single_call_fn,
    reason = "separates sampling and Fourier phase conventions from matrix assembly"
)]
fn coefficients(
    potential: &[f64],
    size: usize,
    scale: f64,
) -> Vec<Complex64> {
    let samples = potential.len();
    let mut coefficients = vec![Complex64::new(0.0, 0.0); size];
    if potential
        .iter()
        .all(|&sample| crate::kernels::exact_eq(sample, potential[0]))
    {
        // A constant sequence has only the zeroth coefficient. Frequencies
        // above the grid's Nyquist limit retain the original sampled aliasing.
        for entry in coefficients.iter_mut().step_by(samples) {
            entry.re = potential[0] / scale;
        }
        return coefficients;
    }
    let plan = fft_plan(samples, FftDirection::Forward);
    let mut spectrum = fft_scratch_for(samples);
    spectrum.extend(
        potential
            .iter()
            .map(|&value| Complex64::new(value / scale, 0.0)),
    );
    let mut scratch = fft_scratch_for(plan.get_inplace_scratch_len());
    scratch.resize(plan.get_inplace_scratch_len(), Complex64::new(0.0, 0.0));
    plan.process_with_scratch(&mut spectrum, &mut scratch);
    let weight = (samples as f64).recip();
    for (frequency, entry) in coefficients.iter_mut().enumerate() {
        let bin = frequency % samples;
        *entry = if bin > samples / 2 {
            spectrum[samples - bin].conj() * weight
        } else {
            spectrum[bin] * weight
        };
        if bin == 0 || 2 * bin == samples {
            entry.im = 0.0;
        }
    }
    coefficients
}

fn decompose<const PARALLEL: bool>(
    period: f64,
    harmonics: usize,
    potential: &[f64],
) -> Option<Arc<Solved>> {
    let kinetic = |frequency: usize| (TAU * frequency as f64 / period).powi(2);
    let largest_kinetic = kinetic(harmonics);
    if !largest_kinetic.is_finite() || potential.iter().any(|sample| !sample.is_finite()) {
        return None;
    }
    let largest = potential
        .iter()
        .fold(largest_kinetic, |max, value| max.max(value.abs()));
    // A binary scale retains normal mantissas and keeps Jacobi's absolute
    // stopping threshold relative to the assembled matrix.
    let exponent = f64::from_bits(largest.to_bits() & 0x7ff0_0000_0000_0000);
    let scale = if exponent > 0.0 {
        exponent
    } else if largest > 0.0 {
        largest
    } else {
        1.0
    };
    let d = 2 * harmonics + 1;
    let coefficients = coefficients(potential, d, scale);
    let mut h = vec![vec![0.0; d]; d];
    let row = |i: usize, row: &mut Vec<f64>| {
        for (j, entry) in row.iter_mut().enumerate() {
            *entry = potential_entry(&coefficients, i, j);
        }
        row[i] += kinetic(i.div_ceil(2)) / scale;
    };
    if PARALLEL {
        crate::kernels::par::for_each_mut(&mut h, 65, row);
    } else {
        for (i, entries) in h.iter_mut().enumerate() {
            row(i, entries);
        }
    }
    let (mut values, vectors) = jacobi(h);
    for value in &mut values {
        *value *= scale;
        if !value.is_finite() {
            return None;
        }
    }
    Some(Arc::new(Solved {
        l: period,
        values,
        vectors,
    }))
}

/// The decomposition for `L, N, V…`.
fn solved(args: &[f64]) -> Option<Arc<Solved>> {
    if let Some(found) = LAST.with(|entry| {
        let entry = entry.borrow();
        let (key, value) = entry.as_ref()?;
        key.iter()
            .copied()
            .eq(args.iter().map(|value| value.to_bits()))
            .then(|| Arc::clone(value))
    }) {
        return Some(found);
    }
    let &[period, harmonics, ref potential @ ..] = args else {
        return None;
    };
    let harmonics = index(harmonics).filter(|&value| value <= MAX_HARMONICS)?;
    if !period.is_finite() || period <= 0.0 || potential.len() < 8 {
        return None;
    }
    if potential.len() > MAX_CACHED_SAMPLES {
        return decompose::<true>(period, harmonics, potential);
    }
    let key: Arc<[u64]> = args.iter().map(|value| value.to_bits()).collect();
    // A pending plan can have callers waiting on Rayon workers. Its builder
    // completes on the calling thread, independently of those workers.
    let value = CACHE
        .try_get_or_insert_with(Arc::clone(&key), || {
            decompose::<false>(period, harmonics, potential).ok_or(())
        })
        .ok()?;
    LAST.with(|entry| *entry.borrow_mut() = Some((key, Arc::clone(&value))));
    Some(value)
}

#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::float_cmp
)]
fn index(m: f64) -> Option<usize> {
    if (1.0..=1e4).contains(&m) && m.fract() == 0.0 {
        Some(m as usize)
    } else {
        None
    }
}

/// `hill_eigen(m, L, N, V…)`.
pub(super) fn hill_eigen(args: &[f64]) -> f64 {
    let (Some(m), Some(rest)) = (args.first().and_then(|&m| index(m)), args.get(1..)) else {
        return f64::NAN;
    };
    solved(rest)
        .and_then(|s| s.values.get(m - 1).copied())
        .unwrap_or(f64::NAN)
}

/// `hill_mode(m, ξ, L, N, V…)`.
pub(super) fn hill_mode(args: &[f64]) -> f64 {
    let (Some(m), Some(&x), Some(rest)) = (
        args.first().and_then(|&m| index(m)),
        args.get(1),
        args.get(2..),
    ) else {
        return f64::NAN;
    };
    if !x.is_finite() {
        return f64::NAN;
    }
    let Some(s) = solved(rest) else {
        return f64::NAN;
    };
    let Some(vector) = s.vectors.get(m - 1) else {
        return f64::NAN;
    };
    let phase = TAU * ((x % s.l) / s.l);
    let mut sum = vector[0];
    for (frequency, &[cosine, sine]) in vector[1..].as_chunks::<2>().0.iter().enumerate() {
        let (sin, cos) = (phase * (frequency + 1) as f64).sin_cos();
        sum += SQRT_2 * cosine * cos;
        sum += SQRT_2 * sine * sin;
    }
    sum / s.l.sqrt()
}

/// Registers the operators.
///
/// Hill's method is a numeric Fourier-Galerkin solve: the eigenvalue `μ` and
/// the eigenfunction are pure numbers, so both operators are dimensionless.
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    for (name, eval, doc) in [
        (
            "hill_eigen",
            hill_eigen as fn(&[f64]) -> f64,
            "`hill_eigen(n, L, harmonics, potential_samples...)`: the `n`-th periodic Sturm-Liouville eigenvalue in ascending order. Uses 1–64 harmonics and at least eight finite, equally spaced potential samples on `[0, L)`.",
        ),
        (
            "hill_mode",
            hill_mode,
            "`hill_mode(n, x, L, harmonics, potential_samples...)`: the `n`-th normalized periodic eigenfunction at offset `x` from the start of the period. Shares the Fourier-Galerkin decomposition with `hill_eigen`.",
        ),
    ] {
        let op = i.op(OpDescriptor::new(name, Arity::Variadic).eval(eval).doc(doc))?;
        i.dim(op, DimRule::Dimensionless);
    }
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
    let (ca, cb) = (
        p.coefficient(cx.graph, &p.unit(j, 2)),
        p.coefficient(cx.graph, &p.unit(j, 1)),
    );
    let cc = p.coefficient(cx.graph, &vec![0; n]);
    if cx.is_zero(ca) || !p.constant(cx.graph, ca) || !cx.is_zero(cb) {
        return None;
    }
    let (a1, a2) = (
        p.coefficient(cx.graph, &p.unit(time_index, 1)),
        p.coefficient(cx.graph, &p.unit(time_index, 2)),
    );
    if !p.constant(cx.graph, a1) || !p.constant(cx.graph, a2) {
        return None;
    }
    for (index, c) in &p.linear {
        let known = *index == p.unit(j, 2)
            || *index == vec![0; n]
            || *index == p.unit(time_index, 1)
            || *index == p.unit(time_index, 2);
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
    let order = if is_zero_number(cx.graph, a2) {
        u32::from(!is_zero_number(cx.graph, a1))
    } else {
        2
    };
    let time = Time {
        var: p.vars[time_index],
        order,
        a2,
        a1,
    };
    let axes = vec![Axis {
        var: x,
        shape: Shape::Hill {
            start: iv.start,
            stop: iv.stop,
            params,
        },
    }];
    let mut items = time_items(cx, p, conditions, Some(time_index), order, n)?;
    if !p.homogeneous(cx.graph) {
        items.push(Item {
            role: Role::Source,
            h: p.source,
            special: None,
        });
    }
    let engine = Engine {
        axes,
        time: Some(time),
        symbol: Symbol::Cartesian(vec![(vec![1], ca)]),
        face: None,
        exterior: None,
    };
    finish(cx, p, conditions, &engine, &items, None, None)
}
