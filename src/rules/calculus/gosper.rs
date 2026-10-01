//! Gosper's algorithm: indefinite summation of hypergeometric terms.
//!
//! A term `t(k)` is hypergeometric when `t(k+1)/t(k)` is a rational
//! function of `k` — polynomials, powers `c^k`, factorials, binomial
//! coefficients and their products and quotients all are. Gosper's
//! algorithm decides whether such a term has a hypergeometric
//! antidifference `T(k)` with `T(k+1) - T(k) = t(k)` and finds it; then
//! `Σ_{k=a}^{b} t(k) = T(b+1) - T(a)` in closed form.
//!
//! The ratio is written `c A(k)/B(k)` with `A, B` over the rationals and a
//! constant `c` that may be symbolic. Gosper's normal form
//! `A/B = p(k+1)/p(k) · q(k)/r(k+1)` (with `q(k)` and `r(k+j)` coprime for
//! every `j ≥ 1`) reduces the problem to a polynomial solution `f` of
//! `c q(k) f(k+1) - r(k) f(k) = p(k)`; the linear system for its
//! coefficients is solved exactly over `Q(c)`. The result is checked
//! numerically before it is returned.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Zero;

use crate::graph::op::core;
use crate::graph::Cx;
use crate::graph::Env;
use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::rules::poly::best;
use crate::rules::poly::ratio;
use crate::rules::poly::repr::Gens;
use crate::rules::poly::repr::Limits;
use crate::rules::poly::repr::Poly;
use crate::rules::poly::univariate;
use crate::rules::poly::univariate::QPoly;

/// Largest shift tried when bringing the ratio into Gosper form.
const MAX_SHIFT: i64 = 48;
/// Largest degree tried for the polynomial `f`.
const MAX_DEGREE: usize = 24;

fn q(v: i64) -> BigRational {
    BigRational::from_integer(BigInt::from(v))
}

fn trim(mut p: QPoly) -> QPoly {
    while p.last().is_some_and(Zero::is_zero) {
        p.pop();
    }
    p
}

/// `P(k + s)`.
fn shift(
    p: &[BigRational],
    s: i64,
) -> QPoly {
    let mut out = vec![BigRational::zero(); p.len()];
    let s = q(s);
    for (i, a) in p.iter().enumerate() {
        // a (k + s)^i = a Σ_j C(i, j) s^(i-j) k^j
        let mut binomial = BigRational::one();
        for j in (0..=i).rev() {
            // C(i, j) s^(i - j)
            let power = num_traits::pow(s.clone(), i - j);
            out[j] += a * &binomial * power;
            if j > 0 {
                binomial = binomial * q(i64::try_from(j).unwrap_or(0)) / q(i64::try_from(i - j + 1).unwrap_or(1));
            }
        }
    }
    trim(out)
}

/// A rational function of one parameter γ over `Q`: numerator and
/// denominator in `Q[γ]`.
#[derive(Clone, Debug, PartialEq)]
pub(super) struct Frac {
    pub(super) n: QPoly,
    pub(super) d: QPoly,
}

impl Frac {
    pub(super) fn zero() -> Self {
        Self { n: Vec::new(), d: vec![BigRational::one()] }
    }

    pub(super) fn poly(p: QPoly) -> Self {
        Self { n: trim(p), d: vec![BigRational::one()] }
    }

    pub(super) const fn is_zero(&self) -> bool {
        self.n.is_empty()
    }

    pub(super) fn reduce(self) -> Self {
        if self.n.is_empty() {
            return Self::zero();
        }
        let g = univariate::gcd(&self.n, &self.d);
        let (Some((n, _)), Some((d, _))) = (univariate::divrem(&self.n, &g), univariate::divrem(&self.d, &g)) else {
            return self;
        };
        // Make the denominator monic.
        let lead = d.last().cloned().unwrap_or_else(BigRational::one);
        Self { n: n.iter().map(|c| c / &lead).collect(), d: d.iter().map(|c| c / &lead).collect() }
    }

    pub(super) fn sub(
        &self,
        o: &Self,
    ) -> Self {
        let n = univariate::sub(&univariate::mul(&self.n, &o.d), &univariate::mul(&o.n, &self.d));
        Self { n, d: univariate::mul(&self.d, &o.d) }.reduce()
    }

    pub(super) fn mul(
        &self,
        o: &Self,
    ) -> Self {
        Self { n: univariate::mul(&self.n, &o.n), d: univariate::mul(&self.d, &o.d) }.reduce()
    }

    pub(super) fn div(
        &self,
        o: &Self,
    ) -> Option<Self> {
        if o.is_zero() {
            return None;
        }
        Some(Self { n: univariate::mul(&self.n, &o.d), d: univariate::mul(&self.d, &o.n) }.reduce())
    }
}

/// A solution of `rows · x = rhs` over `Q(γ)` (free unknowns set to 0),
/// or `None` if the system is inconsistent.
pub(super) fn solve(
    mut rows: Vec<Vec<Frac>>,
    mut rhs: Vec<Frac>,
    n: usize,
) -> Option<Vec<Frac>> {
    let m = rows.len();
    let mut pivots = Vec::new();
    let mut row = 0;
    for col in 0..n {
        let Some(p) = (row..m).find(|&r| !rows[r][col].is_zero()) else {
            continue;
        };
        rows.swap(row, p);
        rhs.swap(row, p);
        for r in 0..m {
            if r == row || rows[r][col].is_zero() {
                continue;
            }
            let factor = rows[r][col].div(&rows[row][col])?;
            #[allow(clippy::needless_range_loop)] // index is used for more than one array / arithmetic; iterator form would not be clearer
            for c in col..n {
                let delta = factor.mul(&rows[row][c]);
                rows[r][c] = rows[r][c].sub(&delta);
            }
            let delta = factor.mul(&rhs[row]);
            rhs[r] = rhs[r].sub(&delta);
        }
        pivots.push((row, col));
        row += 1;
        if row == m {
            break;
        }
    }
    if rhs[row..].iter().any(|v| !v.is_zero()) {
        return None;
    }
    let mut x = vec![Frac::zero(); n];
    for &(r, c) in &pivots {
        x[c] = rhs[r].div(&rows[r][c])?;
    }
    Some(x)
}

/// The ratio `t(k+1)/t(k)` as `(c, A, B)`: `c` free of `k`, `A`, `B` over
/// `Q` with `B` monic.
fn hypergeometric_ratio(
    cx: &mut Cx<'_>,
    t: NodeId,
    k: NodeId,
) -> Option<(NodeId, QPoly, QPoly)> {
    let symbol = cx.graph.symbol_of(k)?;
    let one = cx.graph.int(1);
    let next_index = cx.graph.node(core::ADD, &[k, one]);
    let next = cx.graph.substitute(t, k, next_index);
    let minus_one = cx.graph.int(-1);
    let inverse = cx.graph.node(core::POW, &[t, minus_one]);
    let quotient = cx.graph.node(core::MUL, &[next, inverse]);
    let rho = cx.simplify(quotient);
    let graph = &mut *cx.graph;
    let mut gens = Gens::default();
    let gk = gens.index(graph, k);
    let fraction = ratio(graph, &mut gens, rho, Limits { terms: 256, exponent: 32 })?;
    for g in fraction.numer.support().into_iter().chain(fraction.denom.support()) {
        if g != gk && gens.node(g).is_some_and(|n| graph.depends_on(graph.find(n), symbol)) {
            return None;
        }
    }
    // Each side: a parameter-only leading coefficient times a polynomial
    // over Q in k.
    let split = |p: &Poly| -> Option<(Poly, QPoly)> {
        let coefficients = p.coefficients_in(gk);
        let lead = coefficients.iter().rev().find(|c| !c.is_zero())?.clone();
        let mut out = Vec::with_capacity(coefficients.len());
        for c in &coefficients {
            if c.is_zero() {
                out.push(BigRational::zero());
                continue;
            }
            let quotient = c.div_exact(&lead, 64)?;
            out.push(quotient.as_constant()?.to_rational()?);
        }
        Some((lead, trim(out)))
    };
    let (lead_n, a) = split(&fraction.numer)?;
    let (lead_d, b) = split(&fraction.denom)?;
    let top = crate::rules::poly::repr::to_term(graph, &gens, &lead_n);
    let bottom = crate::rules::poly::repr::to_term(graph, &gens, &lead_d);
    let inverse = graph.node(core::POW, &[bottom, minus_one]);
    let c = graph.node(core::MUL, &[top, inverse]);
    let c = cx.simplify(c);
    Some((c, a, b))
}

/// A hypergeometric antidifference of `t` in `k`: `T` with
/// `T(k+1) - T(k) = t(k)`, if one exists.
///
/// The index is an integer; with `nonnegative` it is also known to be
/// at least zero (sums from a non-negative lower limit), which lets
/// factorial and binomial ratios simplify.
pub fn antidifference(
    cx: &mut Cx<'_>,
    t: NodeId,
    k: NodeId,
    nonnegative: bool,
) -> Option<NodeId> {
    let t = best(cx.graph, t)?;
    // Work with a fresh index that carries what is known about it.
    let fresh = cx.graph.interner_mut().fresh_symbol("k");
    let facts = if nonnegative {
        crate::graph::Facts::INTEGER | crate::graph::Facts::NONNEGATIVE
    } else {
        crate::graph::Facts::INTEGER
    };
    cx.graph.assume(fresh, facts);
    let index = cx.graph.symbol_node(fresh);
    let local = cx.graph.substitute(t, k, index);
    let found = gosper(cx, local, index)?;
    Some(cx.graph.substitute(found, index, k))
}

fn gosper(
    cx: &mut Cx<'_>,
    t: NodeId,
    k: NodeId,
) -> Option<NodeId> {
    let (c, a, b) = hypergeometric_ratio(cx, t, k)?;
    if a.is_empty() || b.is_empty() {
        return None;
    }
    // A rational c is absorbed into q; a symbolic one becomes γ.
    let c_rational = cx.graph.number_of(c).and_then(Number::to_rational);
    let (mut qk, mut rk, mut pk): (QPoly, QPoly, QPoly) = (a, shift(&b, -1), vec![BigRational::one()]);
    // Gosper form: remove common factors of q(k) and r(k + j), j ≥ 1.
    let mut j = 1;
    while j <= MAX_SHIFT {
        let g = univariate::gcd(&qk, &shift(&rk, j));
        if g.len() > 1 {
            qk = univariate::divrem(&qk, &g)?.0;
            rk = univariate::divrem(&rk, &shift(&g, -j))?.0;
            for i in 1..j {
                pk = univariate::mul(&pk, &shift(&g, -i));
            }
            j = 1;
            continue;
        }
        j += 1;
    }
    // c q(k) f(k+1) - r(k) f(k) = p(k), f of degree d over Q(γ).
    let gamma: QPoly = match &c_rational {
        | Some(value) => vec![value.clone()],
        | None => vec![BigRational::zero(), BigRational::one()],
    };
    let degree_p = pk.len().saturating_sub(1);
    let mut solution = None;
    for d in 0..=(degree_p + MAX_DEGREE).min(MAX_DEGREE + 8) {
        let unknowns = d + 1;
        let size = (qk.len().max(rk.len()) + d).max(pk.len());
        let mut rows = vec![vec![Frac::zero(); unknowns]; size];
        let mut rhs = vec![Frac::zero(); size];
        for (m, value) in pk.iter().enumerate() {
            rhs[m] = Frac::poly(vec![value.clone()]);
        }
        for i in 0..unknowns {
            // k^i shifted: (k + 1)^i
            let mut monomial = vec![BigRational::zero(); i + 1];
            monomial[i] = BigRational::one();
            let shifted = shift(&monomial, 1);
            let left = univariate::mul(&qk, &shifted);
            let right = univariate::mul(&rk, &monomial);
            for (m, coefficient) in left.iter().enumerate() {
                if m < size {
                    let term = univariate::mul(&gamma, std::slice::from_ref(coefficient));
                    rows[m][i] = Frac::poly(univariate::add(&rows[m][i].n, &term));
                }
            }
            for (m, coefficient) in right.iter().enumerate() {
                if m < size {
                    rows[m][i] = Frac::poly(univariate::sub(&rows[m][i].n, std::slice::from_ref(coefficient)));
                }
            }
        }
        if let Some(x) = solve(rows, rhs, unknowns) {
            solution = Some(x);
            break;
        }
    }
    let f = solution?;
    // T(k) = r(k) f(k) / p(k) · t(k), with γ = c.
    let graph = &mut *cx.graph;
    let from_q = |graph: &mut Graph, p: &[BigRational], x: NodeId| -> NodeId {
        let terms: Vec<NodeId> = p
            .iter()
            .enumerate()
            .filter(|(_, c)| !c.is_zero())
            .map(|(i, c)| {
                let coefficient = graph.num(Number::rat(c.clone()));
                let e = graph.int(i64::try_from(i).unwrap_or(0));
                let power = graph.node(core::POW, &[x, e]);
                graph.node(core::MUL, &[coefficient, power])
            })
            .collect();
        match terms.as_slice() {
            | [] => graph.int(0),
            | [only] => *only,
            | _ => graph.node(core::ADD, &terms),
        }
    };
    let mut f_terms = Vec::with_capacity(f.len());
    for (i, coefficient) in f.iter().enumerate() {
        if coefficient.is_zero() {
            continue;
        }
        let numerator = from_q(graph, &coefficient.n, c);
        let denominator = from_q(graph, &coefficient.d, c);
        let minus = minus_one(graph);
        let inverse = graph.node(core::POW, &[denominator, minus]);
        let e = graph.int(i64::try_from(i).unwrap_or(0));
        let power = graph.node(core::POW, &[k, e]);
        f_terms.push(graph.node(core::MUL, &[numerator, inverse, power]));
    }
    let f_term = match f_terms.as_slice() {
        | [] => graph.int(0),
        | [only] => *only,
        | _ => graph.node(core::ADD, &f_terms),
    };
    let r_term = from_q(graph, &rk, k);
    let p_term = from_q(graph, &pk, k);
    let minus = minus_one(graph);
    let p_inverse = graph.node(core::POW, &[p_term, minus]);
    let big_t = graph.node(core::MUL, &[r_term, f_term, p_inverse, t]);
    let big_t = cx.simplify(big_t);
    verified(cx, big_t, t, k).then_some(big_t)
}

fn minus_one(graph: &mut Graph) -> NodeId {
    graph.int(-1)
}

/// `T(k+1) - T(k) = t(k)` at a few integer points with generic parameters.
fn verified(
    cx: &mut Cx<'_>,
    big_t: NodeId,
    t: NodeId,
    k: NodeId,
) -> bool {
    let Some(symbol) = cx.graph.symbol_of(k) else {
        return false;
    };
    let graph = &mut *cx.graph;
    let one = graph.int(1);
    let next_index = graph.node(core::ADD, &[k, one]);
    let next = graph.substitute(big_t, k, next_index);
    let minus = graph.int(-1);
    let previous = graph.node(core::MUL, &[minus, big_t]);
    let negated_t = graph.node(core::MUL, &[minus, t]);
    let residual = graph.node(core::ADD, &[next, previous, negated_t]);
    let mut checked = 0;
    for point in 3..9_i32 {
        let mut env = Env::numeric(0.0);
        for &s in graph.free_symbols(graph.find(residual)) {
            env.bind(s, 0.31 + 0.07 * f64::from(s.raw() % 11));
        }
        env.bind(symbol, f64::from(point));
        let (Some(value), Some(scale)) = (graph.eval(residual, &env), graph.eval(t, &env)) else {
            continue;
        };
        if !value.is_finite() || !scale.is_finite() {
            continue;
        }
        if value.abs() > 1e-8 * scale.abs().max(1.0) {
            return false;
        }
        checked += 1;
    }
    checked >= 3
}
