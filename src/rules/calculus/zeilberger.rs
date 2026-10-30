//! Zeilberger's creative telescoping for proper hypergeometric terms.
//!
//! A term `F(n, k)` is *proper hypergeometric* when it is a product of
//! rational functions, powers and Gamma functions (hence factorials and
//! binomial coefficients) of linear forms in `n` and `k`. For such a term
//! Zeilberger's algorithm finds a linear recurrence
//!
//! ```text
//! Σ_{j=0}^{J} p_j(n) F(n + j, k) = G(n, k + 1) - G(n, k)
//! ```
//!
//! with polynomial coefficients `p_j` and a rational *certificate* `G/F`;
//! summing over `k` makes the right-hand side telescope, so the definite sum
//! `S(n) = Σ_k F(n, k)` satisfies `Σ_j p_j(n) S(n + j) = 0`.
//!
//! The term is first factorised into a rational part and Gamma and
//! exponential factors. Its `k`-ratio `F(n, k+1)/F(n, k) = A/B` is then a
//! rational function, and the certificate is sought with denominator `B`.
//! Multiplying the telescoping identity by the common denominator turns it
//! into a polynomial identity in `n`, `k`, the parameters and the unknown
//! coefficients of the `p_j` and of the certificate numerator. Collecting
//! monomials gives a linear system over the rationals, solved exactly here;
//! every candidate is checked numerically before it is returned.
//!
//! The operators are
//!
//! | operator | value |
//! |---|---|
//! | `zeilberger(F, k, n)` | `list(list(p_0, …, p_J), G/F)`, the recurrence `Σ_j p_j(n) S(n+j) = 0` and the certificate |
//! | `definite_sum(F, k, n)` | the closed form of `Σ_{k=-oo}^{oo} F(n, k)` when the recurrence is of first order and solvable |
//!
//! `sum(F, k, -oo, oo)` also tries this route before the numeric kernel.
//!
//! Two searches are used.
//!
//! * A term with free parameters is searched with a bounded polynomial
//!   ansatz: the `p_j` and the certificate numerator are polynomials in `n`
//!   and the parameters of fixed degree (orders one and two).
//! * A parameter-free term is searched exactly over `Q(n)` for orders one to
//!   three, with no bound on the degrees in `n`: the system is solved for
//!   integer values of `n` and the solution is recovered as a rational
//!   function of `n` by Cauchy interpolation (see "The exact search over
//!   Q(n)" below). The certificate may have poles in `n`. This finds the
//!   second-order recurrences of Dixon's `Σ (-1)^k C(n,k)³` and Apéry's
//!   `Σ C(n,k)² C(n+k,k)²` and the recurrence of a Pfaff–Saalschütz 3F2.
//!
//! Limitations: the certificate is of the form `U/(B·Q)` with `U` a
//! polynomial in `k` (`B` the denominator of the `k`-ratio, `Q` a common
//! multiple of the denominators of the `n`-ratios); recurrences of order
//! above three, terms with symbolic parameters and orders above two, and
//! terms whose ratios are not rational (non-integer Gamma shifts) are not
//! covered. The recurrence is a recurrence of the summand: a sum without
//! natural boundaries (a 3F2 written with shifted factorials) also needs
//! the boundary values of the certificate before it becomes a recurrence of
//! the sum. `definite_sum` solves first-order recurrences only.

#![allow(
    clippy::single_call_fn,
    clippy::too_many_arguments,
    clippy::too_many_lines
)]

use std::collections::BTreeMap;

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::ToPrimitive;
use num_traits::Zero;

use crate::graph::Cx;
use crate::graph::Env;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::SymbolId;
use crate::graph::op::core;
use crate::rules::poly::best;
use crate::rules::poly::repr::Gens;
use crate::rules::poly::repr::Limits;
use crate::rules::poly::repr::from_term;
use crate::rules::poly::univariate;
use crate::rules::poly::univariate::QPoly;

use super::gosper::Frac;
use super::recurrence;

/// Largest shift applied to a Gamma factor while forming a ratio.
const MAX_GAMMA_SHIFT: i64 = 64;
/// Largest number of unknowns in one linear system.
const MAX_UNKNOWNS: usize = 72;
/// Largest degree of the parameter monomials in the ansatz.
const MAX_PARAMETER_DEGREE: usize = 2;
/// Number of further `n` values at which a closed form is checked.
const CHECKS: i64 = 5;

fn rat(value: i64) -> BigRational {
    BigRational::from_integer(BigInt::from(value))
}

fn add_of(
    graph: &mut Graph,
    terms: &[NodeId],
) -> NodeId {
    match terms {
        | [] => graph.int(0),
        | [only] => *only,
        | _ => graph.node(core::ADD, terms),
    }
}

fn mul_of(
    graph: &mut Graph,
    factors: &[NodeId],
) -> NodeId {
    match factors {
        | [] => graph.int(1),
        | [only] => *only,
        | _ => graph.node(core::MUL, factors),
    }
}

/// `base^exponent` for a non-negative literal exponent.
fn pow_of(
    graph: &mut Graph,
    base: NodeId,
    exponent: i64,
) -> NodeId {
    match exponent {
        | 0 => graph.int(1),
        | 1 => base,
        | e => {
            let n = graph.int(e);
            graph.node(core::POW, &[base, n])
        },
    }
}

/// `base^exponent` for a rational exponent.
fn pow_of_rat(
    graph: &mut Graph,
    base: NodeId,
    exponent: &BigRational,
) -> NodeId {
    if exponent.is_zero() {
        return graph.int(1);
    }
    if exponent.is_one() {
        return base;
    }
    let e = graph.num(Number::rat(exponent.clone()));
    graph.node(core::POW, &[base, e])
}

fn sub_of(
    graph: &mut Graph,
    a: NodeId,
    b: NodeId,
) -> NodeId {
    let minus_one = graph.int(-1);
    let negated = graph.node(core::MUL, &[minus_one, b]);
    graph.node(core::ADD, &[a, negated])
}

/// `r x` as a term.
fn scaled(
    graph: &mut Graph,
    x: NodeId,
    r: &BigRational,
) -> NodeId {
    if r.is_one() {
        return x;
    }
    let coefficient = graph.num(Number::rat(r.clone()));
    graph.node(core::MUL, &[coefficient, x])
}

/// `a n + b k + offset`.
fn linear_node(
    graph: &mut Graph,
    offset: NodeId,
    a: &BigRational,
    b: &BigRational,
    n: NodeId,
    k: NodeId,
) -> NodeId {
    let mut terms = vec![offset];
    if !a.is_zero() {
        terms.push(scaled(graph, n, a));
    }
    if !b.is_zero() {
        terms.push(scaled(graph, k, b));
    }
    add_of(graph, &terms)
}

/// The symbols of a telescoping problem.
#[derive(Copy, Clone)]
struct Symbols {
    n: NodeId,
    k: NodeId,
    ns: SymbolId,
    ks: SymbolId,
}

/// One Gamma factor `Γ(linear form)^exponent` of the term.
struct GammaFactor {
    n: i64,
    k: i64,
    offset: NodeId,
    exponent: i64,
}

/// One exponential factor `base^(n_coefficient·n + k_coefficient·k)`.
struct ExpFactor {
    base: NodeId,
    n: BigRational,
    k: BigRational,
}

/// A hypergeometric term split into a rational part and Gamma and
/// exponential factors. Constant factors cancel from every ratio and are
/// dropped by [`factorize`].
struct Factored {
    num: NodeId,
    den: NodeId,
    gammas: Vec<GammaFactor>,
    exps: Vec<ExpFactor>,
}

fn one(graph: &mut Graph) -> Factored {
    Factored {
        num: graph.int(1),
        den: graph.int(1),
        gammas: Vec::new(),
        exps: Vec::new(),
    }
}

fn constant(
    graph: &mut Graph,
    node: NodeId,
) -> Factored {
    Factored {
        num: node,
        den: graph.int(1),
        gammas: Vec::new(),
        exps: Vec::new(),
    }
}

fn combine(
    graph: &mut Graph,
    first: Factored,
    second: Factored,
) -> Factored {
    Factored {
        num: mul_of(graph, &[first.num, second.num]),
        den: mul_of(graph, &[first.den, second.den]),
        gammas: first.gammas.into_iter().chain(second.gammas).collect(),
        exps: first.exps.into_iter().chain(second.exps).collect(),
    }
}

fn integer_of(r: &BigRational) -> Option<i64> {
    r.is_integer().then(|| r.to_integer().to_i64()).flatten()
}

/// The rational coefficients `(a, b)` and constant term `c` of the linear
/// form `u = a n + b k + c`, or `None` when `u` is not linear.
fn linear_form(
    cx: &mut Cx<'_>,
    u: NodeId,
    symbols: Symbols,
) -> Option<(BigRational, BigRational, NodeId)> {
    fn at(
        cx: &mut Cx<'_>,
        u: NodeId,
        symbols: Symbols,
        nv: NodeId,
        kv: NodeId,
    ) -> Option<NodeId> {
        let substituted = cx.graph.substitute(u, symbols.n, nv);
        let substituted = cx.graph.substitute(substituted, symbols.k, kv);
        best(cx.graph, substituted)
    }
    let zero = cx.graph.int(0);
    let one = cx.graph.int(1);
    let constant = at(cx, u, symbols, zero, zero)?;
    let at_n = at(cx, u, symbols, one, zero)?;
    let at_k = at(cx, u, symbols, zero, one)?;
    let difference_n = sub_of(cx.graph, at_n, constant);
    let difference_k = sub_of(cx.graph, at_k, constant);
    let simplified_n = cx.simplify(difference_n);
    let simplified_k = cx.simplify(difference_k);
    let a = cx.graph.number_of(simplified_n)?.to_rational()?;
    let b = cx.graph.number_of(simplified_k)?.to_rational()?;
    let rebuilt = linear_node(cx.graph, constant, &a, &b, symbols.n, symbols.k);
    let residual = sub_of(cx.graph, u, rebuilt);
    cx.is_zero(residual).then_some((a, b, constant))
}

fn depends(
    graph: &Graph,
    node: NodeId,
    symbol: SymbolId,
) -> bool {
    graph.depends_on(graph.find(node), symbol)
}

fn free_of(
    graph: &Graph,
    node: NodeId,
    symbols: Symbols,
) -> bool {
    !depends(graph, node, symbols.ns) && !depends(graph, node, symbols.ks)
}

/// Factorises `term` for the ratios used by creative telescoping.
fn factorize(
    cx: &mut Cx<'_>,
    term: NodeId,
    symbols: Symbols,
) -> Option<Factored> {
    let term = best(cx.graph, term)?;
    factor_term(cx, term, symbols)
}

fn power_factored(
    graph: &mut Graph,
    factored: Factored,
    exponent: i64,
) -> Option<Factored> {
    let magnitude = exponent.checked_abs()?;
    let (num, den) = if exponent < 0 {
        (factored.den, factored.num)
    } else {
        (factored.num, factored.den)
    };
    let mut gammas = Vec::with_capacity(factored.gammas.len());
    for factor in factored.gammas {
        gammas.push(GammaFactor {
            n: factor.n,
            k: factor.k,
            offset: factor.offset,
            exponent: factor.exponent.checked_mul(exponent)?,
        });
    }
    let mut exps = Vec::with_capacity(factored.exps.len());
    let scale = rat(exponent);
    for factor in factored.exps {
        exps.push(ExpFactor {
            base: factor.base,
            n: &factor.n * &scale,
            k: &factor.k * &scale,
        });
    }
    Some(Factored {
        num: pow_of(graph, num, magnitude),
        den: pow_of(graph, den, magnitude),
        gammas,
        exps,
    })
}

fn gamma_factor(
    cx: &mut Cx<'_>,
    u: NodeId,
    symbols: Symbols,
    exponent: i64,
) -> Option<Factored> {
    let (a, b, offset) = linear_form(cx, u, symbols)?;
    let mut result = one(cx.graph);
    result.gammas.push(GammaFactor {
        n: integer_of(&a)?,
        k: integer_of(&b)?,
        offset,
        exponent,
    });
    Some(result)
}

fn factor_term(
    cx: &mut Cx<'_>,
    term: NodeId,
    symbols: Symbols,
) -> Option<Factored> {
    if free_of(cx.graph, term, symbols)
        || cx.graph.same(term, symbols.n)
        || cx.graph.same(term, symbols.k)
    {
        return Some(constant(cx.graph, term));
    }
    let op = cx.graph.op(term);
    let children = cx.graph.children(term).to_vec();
    if op == core::MUL {
        let mut accumulated = one(cx.graph);
        for child in children {
            let part = factor_term(cx, child, symbols)?;
            accumulated = combine(cx.graph, accumulated, part);
        }
        return Some(accumulated);
    }
    if op == core::POW {
        let &[base, exponent] = children.as_slice() else {
            return None;
        };
        if let Some(e) = cx.graph.number_of(exponent).and_then(Number::to_i64) {
            let inner = factor_term(cx, base, symbols)?;
            return power_factored(cx.graph, inner, e);
        }
        if free_of(cx.graph, base, symbols) {
            let (a, b, _) = linear_form(cx, exponent, symbols)?;
            let mut result = one(cx.graph);
            result.exps.push(ExpFactor { base, n: a, k: b });
            return Some(result);
        }
        return None;
    }
    let name = cx.graph.ops().get(op).name.to_string();
    match (name.as_str(), children.as_slice()) {
        | ("exp", &[u]) => {
            let (a, b, _) = linear_form(cx, u, symbols)?;
            let exp = cx.graph.ops().lookup("exp")?;
            let one_node = cx.graph.int(1);
            let base = cx.graph.node(exp, &[one_node]);
            let mut result = one(cx.graph);
            result.exps.push(ExpFactor { base, n: a, k: b });
            Some(result)
        },
        | ("factorial", &[u]) => {
            let (a, b, offset) = linear_form(cx, u, symbols)?;
            let one_node = cx.graph.int(1);
            let offset = cx.graph.node(core::ADD, &[offset, one_node]);
            let mut result = one(cx.graph);
            result.gammas.push(GammaFactor {
                n: integer_of(&a)?,
                k: integer_of(&b)?,
                offset,
                exponent: 1,
            });
            Some(result)
        },
        | ("gamma", &[u]) => gamma_factor(cx, u, symbols, 1),
        | ("binomial", &[u, v]) => {
            let (au, bu, cu) = linear_form(cx, u, symbols)?;
            let (av, bv, cv) = linear_form(cx, v, symbols)?;
            let one_node = cx.graph.int(1);
            let minus_one = cx.graph.int(-1);
            let top = cx.graph.node(core::ADD, &[cu, one_node]);
            let bottom = cx.graph.node(core::ADD, &[cv, one_node]);
            let negative = cx.graph.node(core::MUL, &[minus_one, cv]);
            let rest = cx.graph.node(core::ADD, &[cu, negative, one_node]);
            let mut result = one(cx.graph);
            result.gammas.push(GammaFactor {
                n: integer_of(&au)?,
                k: integer_of(&bu)?,
                offset: top,
                exponent: 1,
            });
            result.gammas.push(GammaFactor {
                n: integer_of(&av)?,
                k: integer_of(&bv)?,
                offset: bottom,
                exponent: -1,
            });
            result.gammas.push(GammaFactor {
                n: integer_of(&(&au - &av))?,
                k: integer_of(&(&bu - &bv))?,
                offset: rest,
                exponent: -1,
            });
            Some(result)
        },
        | _ => None,
    }
}

/// `Γ(l + d) / Γ(l)` with `l` the linear form of `factor`.
fn gamma_shift(
    cx: &mut Cx<'_>,
    factor: &GammaFactor,
    shift: i64,
    symbols: Symbols,
) -> Option<(NodeId, NodeId)> {
    if shift == 0 {
        return Some((cx.graph.int(1), cx.graph.int(1)));
    }
    if shift.unsigned_abs() > u64::try_from(MAX_GAMMA_SHIFT).ok()? {
        return None;
    }
    let line = linear_node(
        cx.graph,
        factor.offset,
        &rat(factor.n),
        &rat(factor.k),
        symbols.n,
        symbols.k,
    );
    let mut factors = Vec::new();
    if shift > 0 {
        for t in 0..shift {
            let add = cx.graph.int(t);
            factors.push(cx.graph.node(core::ADD, &[line, add]));
        }
    } else {
        for t in 1..=shift.unsigned_abs() {
            let minus = cx.graph.int(-i64::try_from(t).ok()?);
            factors.push(cx.graph.node(core::ADD, &[line, minus]));
        }
    }
    let product = mul_of(cx.graph, &factors);
    let magnitude = factor.exponent.checked_abs()?;
    let raised = pow_of(cx.graph, product, magnitude);
    if (shift > 0) == (factor.exponent > 0) {
        Some((raised, cx.graph.int(1)))
    } else {
        Some((cx.graph.int(1), raised))
    }
}

/// A quotient of two polynomial terms.
#[derive(Copy, Clone)]
struct Quotient {
    num: NodeId,
    den: NodeId,
}

/// `F(n + shift, k) / F(n, k)`.
fn n_ratio(
    cx: &mut Cx<'_>,
    factored: &Factored,
    shift: i64,
    symbols: Symbols,
) -> Option<Quotient> {
    let step = cx.graph.int(shift);
    let lifted = cx.graph.node(core::ADD, &[symbols.n, step]);
    let num_shift = cx.graph.substitute(factored.num, symbols.n, lifted);
    let den_shift = cx.graph.substitute(factored.den, symbols.n, lifted);
    let mut num = mul_of(cx.graph, &[num_shift, factored.den]);
    let mut den = mul_of(cx.graph, &[den_shift, factored.num]);
    for factor in &factored.gammas {
        let d = factor.n.checked_mul(shift)?;
        let (top, bottom) = gamma_shift(cx, factor, d, symbols)?;
        num = mul_of(cx.graph, &[num, top]);
        den = mul_of(cx.graph, &[den, bottom]);
    }
    for factor in &factored.exps {
        let exponent = &factor.n * &rat(shift);
        let top = pow_of_rat(cx.graph, factor.base, &exponent);
        num = mul_of(cx.graph, &[num, top]);
    }
    Some(Quotient { num, den })
}

/// `F(n, k + 1) / F(n, k)`.
fn k_ratio(
    cx: &mut Cx<'_>,
    factored: &Factored,
    symbols: Symbols,
) -> Option<Quotient> {
    let one_node = cx.graph.int(1);
    let lifted = cx.graph.node(core::ADD, &[symbols.k, one_node]);
    let num_shift = cx.graph.substitute(factored.num, symbols.k, lifted);
    let den_shift = cx.graph.substitute(factored.den, symbols.k, lifted);
    let mut num = mul_of(cx.graph, &[num_shift, factored.den]);
    let mut den = mul_of(cx.graph, &[den_shift, factored.num]);
    for factor in &factored.gammas {
        let (top, bottom) = gamma_shift(cx, factor, factor.k, symbols)?;
        num = mul_of(cx.graph, &[num, top]);
        den = mul_of(cx.graph, &[den, bottom]);
    }
    for factor in &factored.exps {
        let top = pow_of_rat(cx.graph, factor.base, &factor.k);
        num = mul_of(cx.graph, &[num, top]);
    }
    Some(Quotient { num, den })
}

/// One unknown of the linear system: a fresh symbol with a readable name.
fn unknown(
    graph: &mut Graph,
    prefix: &str,
) -> NodeId {
    let symbol = graph.interner_mut().fresh_symbol(prefix);
    graph.symbol_node(symbol)
}

/// The exponent vectors over the parameters, beginning with the constant
/// one.
fn parameter_monomials(
    count: usize,
    degree: usize,
) -> Vec<Vec<u32>> {
    if count == 0 {
        return vec![Vec::new()];
    }
    let limit = u32::try_from(degree).unwrap_or(0).saturating_add(1);
    let mut out = vec![vec![0_u32; count]];
    let mut current = vec![0_u32; count];
    'outer: loop {
        for position in (0..current.len()).rev() {
            let Some(slot) = current.get_mut(position) else {
                break 'outer;
            };
            *slot += 1;
            if *slot < limit {
                out.push(current.clone());
                continue 'outer;
            }
            *slot = 0;
        }
        break;
    }
    out
}

fn monomial_node(
    graph: &mut Graph,
    variables: &[NodeId],
    exponents: &[u32],
) -> NodeId {
    let factors: Vec<NodeId> = variables
        .iter()
        .zip(exponents)
        .filter(|&(_, &exponent)| exponent > 0)
        .map(|(&node, &exponent)| pow_of(graph, node, i64::from(exponent)))
        .collect();
    mul_of(graph, &factors)
}

/// The candidate unknowns of one ansatz.
struct Ansatz {
    order: usize,
    p_index: Vec<Vec<Vec<usize>>>,
    u_index: Vec<Vec<Vec<usize>>>,
    nodes: Vec<NodeId>,
    parameters: Vec<NodeId>,
    monomials: Vec<Vec<u32>>,
}

fn build_ansatz(
    graph: &mut Graph,
    parameters: &[NodeId],
    order: usize,
    n_degree: usize,
    u_n_degree: usize,
    u_k_degree: usize,
    parameter_degree: usize,
) -> Option<Ansatz> {
    let monomials = parameter_monomials(parameters.len(), parameter_degree);
    let count = (order + 1) * (n_degree + 1) * monomials.len()
        + (u_k_degree + 1) * (u_n_degree + 1) * monomials.len();
    if count > MAX_UNKNOWNS {
        return None;
    }
    let mut nodes = Vec::with_capacity(count);
    let mut p_index = Vec::with_capacity(order + 1);
    let mut u_index = Vec::with_capacity(u_k_degree + 1);
    for _ in 0..=order {
        let mut in_n = Vec::with_capacity(n_degree + 1);
        for _ in 0..=n_degree {
            let mut in_parameter = Vec::with_capacity(monomials.len());
            for _ in &monomials {
                in_parameter.push(nodes.len());
                nodes.push(unknown(graph, "c"));
            }
            in_n.push(in_parameter);
        }
        p_index.push(in_n);
    }
    for _ in 0..=u_k_degree {
        let mut in_k = Vec::with_capacity(u_n_degree + 1);
        for _ in 0..=u_n_degree {
            let mut in_parameter = Vec::with_capacity(monomials.len());
            for _ in &monomials {
                in_parameter.push(nodes.len());
                nodes.push(unknown(graph, "c"));
            }
            in_k.push(in_parameter);
        }
        u_index.push(in_k);
    }
    Some(Ansatz {
        order,
        p_index,
        u_index,
        nodes,
        parameters: parameters.to_vec(),
        monomials,
    })
}

impl Ansatz {
    /// `p_j(n)` as a term over the unknown symbols.
    fn coefficient(
        &self,
        graph: &mut Graph,
        j: usize,
        n: NodeId,
    ) -> NodeId {
        let Some(in_n) = self.p_index.get(j) else {
            return graph.int(0);
        };
        let mut terms = Vec::new();
        for (i, in_parameter) in in_n.iter().enumerate() {
            for (&node, exponents) in in_parameter.iter().zip(&self.monomials) {
                let mut factors = vec![
                    self.nodes
                        .get(node)
                        .copied()
                        .unwrap_or_else(|| graph.int(1)),
                ];
                if i > 0 {
                    factors.push(pow_of(graph, n, i64::try_from(i).unwrap_or(0)));
                }
                factors.push(monomial_node(graph, &self.parameters, exponents));
                terms.push(mul_of(graph, &factors));
            }
        }
        add_of(graph, &terms)
    }

    /// The certificate numerator as a term over the unknown symbols.
    fn numerator(
        &self,
        graph: &mut Graph,
        n: NodeId,
        k: NodeId,
    ) -> NodeId {
        let mut terms = Vec::new();
        for (m, in_k) in self.u_index.iter().enumerate() {
            for (s, in_parameter) in in_k.iter().enumerate() {
                for (&node, exponents) in in_parameter.iter().zip(&self.monomials) {
                    let mut factors = vec![
                        self.nodes
                            .get(node)
                            .copied()
                            .unwrap_or_else(|| graph.int(1)),
                    ];
                    if m > 0 {
                        factors.push(pow_of(graph, k, i64::try_from(m).unwrap_or(0)));
                    }
                    if s > 0 {
                        factors.push(pow_of(graph, n, i64::try_from(s).unwrap_or(0)));
                    }
                    factors.push(monomial_node(graph, &self.parameters, exponents));
                    terms.push(mul_of(graph, &factors));
                }
            }
        }
        add_of(graph, &terms)
    }

    /// The polynomial with numeric coefficients from a solution.
    fn numeric_polynomial(
        &self,
        graph: &mut Graph,
        values: &[BigRational],
        j: Option<usize>,
        n: NodeId,
        k: NodeId,
    ) -> NodeId {
        let mut terms = Vec::new();
        let value_at = |index: usize| values.get(index).cloned().unwrap_or_else(BigRational::zero);
        if let Some(j) = j {
            let Some(in_n) = self.p_index.get(j) else {
                return graph.int(0);
            };
            for (i, in_parameter) in in_n.iter().enumerate() {
                for (&index, exponents) in in_parameter.iter().zip(&self.monomials) {
                    let value = value_at(index);
                    if value.is_zero() {
                        continue;
                    }
                    let mut factors = vec![graph.num(Number::rat(value))];
                    if i > 0 {
                        factors.push(pow_of(graph, n, i64::try_from(i).unwrap_or(0)));
                    }
                    factors.push(monomial_node(graph, &self.parameters, exponents));
                    terms.push(mul_of(graph, &factors));
                }
            }
        } else {
            for (m, in_k) in self.u_index.iter().enumerate() {
                for (s, in_parameter) in in_k.iter().enumerate() {
                    for (&index, exponents) in in_parameter.iter().zip(&self.monomials) {
                        let value = value_at(index);
                        if value.is_zero() {
                            continue;
                        }
                        let mut factors = vec![graph.num(Number::rat(value))];
                        if m > 0 {
                            factors.push(pow_of(graph, k, i64::try_from(m).unwrap_or(0)));
                        }
                        if s > 0 {
                            factors.push(pow_of(graph, n, i64::try_from(s).unwrap_or(0)));
                        }
                        factors.push(monomial_node(graph, &self.parameters, exponents));
                        terms.push(mul_of(graph, &factors));
                    }
                }
            }
        }
        add_of(graph, &terms)
    }
}

/// The identity of one ansatz, as one unsimplified term.
fn identity_node(
    cx: &mut Cx<'_>,
    ansatz: &Ansatz,
    ratios: &[Quotient],
    kq: Quotient,
    symbols: Symbols,
) -> NodeId {
    let graph = &mut *cx.graph;
    let p: Vec<NodeId> = (0..=ansatz.order)
        .map(|j| ansatz.coefficient(graph, j, symbols.n))
        .collect();
    let u = ansatz.numerator(graph, symbols.n, symbols.k);
    let mut lhs_terms = Vec::with_capacity(ratios.len());
    for (j, ratio) in ratios.iter().enumerate() {
        let mut factors = vec![p.get(j).copied().unwrap_or_else(|| graph.int(0)), ratio.num];
        for (i, other) in ratios.iter().enumerate() {
            if i != j {
                factors.push(other.den);
            }
        }
        lhs_terms.push(mul_of(graph, &factors));
    }
    let lhs = add_of(graph, &lhs_terms);
    let common = mul_of(graph, &ratios.iter().map(|r| r.den).collect::<Vec<_>>());
    let one_node = graph.int(1);
    let next = graph.node(core::ADD, &[symbols.k, one_node]);
    let shifted_denominator = graph.substitute(kq.den, symbols.k, next);
    let shifted_u = graph.substitute(u, symbols.k, next);
    let shifted_common = graph.substitute(common, symbols.k, next);
    // The certificate is U/(B·Q), B the denominator of the k-ratio and Q
    // the common denominator of the n-ratios: the factors of Q are exactly
    // those that a shift of n would otherwise leave in the certificate.
    let upper_right = mul_of(graph, &[shifted_u, kq.num, common]);
    let lower_right = mul_of(graph, &[u, shifted_denominator, shifted_common]);
    let right = sub_of(graph, upper_right, lower_right);
    let scaled_lhs = mul_of(graph, &[lhs, kq.den, shifted_denominator, shifted_common]);
    sub_of(graph, scaled_lhs, right)
}

/// A linear system in reduced row echelon form.
struct System {
    rows: Vec<Vec<BigRational>>,
    pivots: Vec<(usize, usize)>,
    free: Vec<usize>,
    columns: usize,
}

fn solve_system(
    rows: Vec<Vec<BigRational>>,
    columns: usize,
) -> Option<System> {
    let mut rows = rows;
    let count = rows.len();
    let mut pivots = Vec::new();
    let mut row = 0;
    for column in 0..columns {
        if row >= count {
            break;
        }
        let Some(pivot) = (row..count).find(|&r| {
            rows.get(r)
                .and_then(|x| x.get(column))
                .is_some_and(|v| !v.is_zero())
        }) else {
            continue;
        };
        rows.swap(row, pivot);
        let lead = rows
            .get(row)
            .and_then(|x| x.get(column))
            .cloned()
            .unwrap_or_else(BigRational::zero);
        for r in 0..count {
            if r == row {
                continue;
            }
            let Some(factor) = rows
                .get(r)
                .and_then(|x| x.get(column))
                .filter(|v| !v.is_zero())
                .cloned()
            else {
                continue;
            };
            let factor = factor / &lead;
            for column_index in column..=columns {
                let offset = rows
                    .get(row)
                    .and_then(|x| x.get(column_index))
                    .cloned()
                    .unwrap_or_else(BigRational::zero);
                if let Some(slot) = rows.get_mut(r).and_then(|x| x.get_mut(column_index)) {
                    *slot -= offset * &factor;
                }
            }
        }
        pivots.push((row, column));
        row += 1;
    }
    for leftover in rows.iter().skip(row) {
        let inhomogeneous = leftover.iter().take(columns).all(Zero::is_zero)
            && !leftover.get(columns).is_some_and(Zero::is_zero);
        if inhomogeneous {
            return None;
        }
    }
    let free: Vec<usize> = (0..columns)
        .filter(|c| !pivots.iter().any(|&(_, p)| p == *c))
        .collect();
    Some(System {
        rows,
        pivots,
        free,
        columns,
    })
}

impl System {
    fn solution(
        &self,
        free_value: Option<(usize, BigRational)>,
    ) -> Vec<BigRational> {
        let mut x = vec![BigRational::zero(); self.columns];
        if let Some((position, value)) = free_value
            && let Some(slot) = x.get_mut(position)
        {
            *slot = value;
        }
        for &(row, column) in self.pivots.iter().rev() {
            let mut total = self
                .rows
                .get(row)
                .and_then(|r| r.get(self.columns))
                .cloned()
                .unwrap_or_else(BigRational::zero);
            for other in 0..self.columns {
                if other == column {
                    continue;
                }
                let coefficient = self
                    .rows
                    .get(row)
                    .and_then(|r| r.get(other))
                    .cloned()
                    .unwrap_or_else(BigRational::zero);
                let value = x.get(other).cloned().unwrap_or_else(BigRational::zero);
                total -= coefficient * value;
            }
            let lead = self
                .rows
                .get(row)
                .and_then(|r| r.get(column))
                .cloned()
                .unwrap_or_else(BigRational::zero);
            if lead.is_zero() {
                continue;
            }
            if let Some(slot) = x.get_mut(column) {
                *slot = total / lead;
            }
        }
        x
    }
}

/// The linear equations of the identity, grouped by monomial in `n`, `k`
/// and the parameters.
fn equations(
    cx: &mut Cx<'_>,
    identity: NodeId,
    ansatz: &Ansatz,
    symbols: Symbols,
) -> Option<Vec<Vec<BigRational>>> {
    let graph = &mut *cx.graph;
    let mut gens = Gens::default();
    let poly = from_term(
        graph,
        &mut gens,
        identity,
        Limits {
            terms: 400_000,
            exponent: 256,
        },
    )?;
    let mut unknown_index = BTreeMap::new();
    for (index, node) in ansatz.nodes.iter().enumerate() {
        let generator = gens.find(graph, *node)?;
        unknown_index.insert(generator, index);
    }
    let columns = ansatz.nodes.len();
    let mut grouped: BTreeMap<Vec<(u32, u32)>, Vec<BigRational>> = BTreeMap::new();
    for (monomial, coefficient) in poly.terms() {
        let mut key = Vec::with_capacity(monomial.len());
        let mut which = None;
        for &(generator, exponent) in monomial {
            if let Some(&index) = unknown_index.get(&generator) {
                if exponent != 1 || which.is_some() {
                    return None;
                }
                which = Some(index);
            } else {
                key.push((generator, exponent));
            }
        }
        for &(generator, _) in &key {
            let node = gens.node(generator)?;
            let known = graph.same(node, symbols.n) || graph.same(node, symbols.k);
            if !known && (depends(graph, node, symbols.ns) || depends(graph, node, symbols.ks)) {
                return None;
            }
        }
        let value = coefficient.to_rational()?;
        let row = grouped
            .entry(key)
            .or_insert_with(|| vec![BigRational::zero(); columns + 1]);
        if let Some(index) = which {
            if let Some(slot) = row.get_mut(index) {
                *slot += value;
            }
        } else if let Some(slot) = row.get_mut(columns) {
            *slot += value;
        }
    }
    Some(grouped.into_values().collect())
}

/// Scales a solution to primitive integer `p_j` coefficients.
fn normalize(
    values: &mut [BigRational],
    p_columns: &[usize],
    order_columns: &[usize],
) {
    let mut denominator = BigInt::one();
    for &column in p_columns {
        if let Some(value) = values.get(column) {
            let d = value.denom();
            if d > &BigInt::zero() {
                denominator = lcm(&denominator, d);
            }
        }
    }
    let scale = BigRational::from_integer(denominator);
    for value in values.iter_mut() {
        *value = value.clone() * &scale;
    }
    let mut divisor = BigInt::zero();
    for &column in p_columns {
        if let Some(value) = values.get(column) {
            divisor = gcd(&divisor, value.numer());
        }
    }
    if divisor > BigInt::one() {
        let divisor = BigRational::from_integer(divisor);
        for value in values.iter_mut() {
            *value = value.clone() / &divisor;
        }
    }
    let leading = order_columns
        .iter()
        .rev()
        .find_map(|&column| values.get(column).filter(|v| !v.is_zero()).cloned());
    if leading.is_some_and(|v| v.is_negative()) {
        for value in values.iter_mut() {
            *value = -value.clone();
        }
    }
}

fn gcd(
    a: &BigInt,
    b: &BigInt,
) -> BigInt {
    let (mut x, mut y) = (a.abs(), b.abs());
    while !y.is_zero() {
        let r = &x % &y;
        x = y;
        y = r;
    }
    x
}

fn lcm(
    a: &BigInt,
    b: &BigInt,
) -> BigInt {
    if a.is_zero() || b.is_zero() {
        return BigInt::zero();
    }
    (a / &gcd(a, b)) * b
}

/// The answer of a successful telescoping search.
pub(super) struct Telescoping {
    pub(super) coefficients: Vec<NodeId>,
    pub(super) certificate: NodeId,
}

/// The degree tuples tried for one order, richest first. A solution with
/// smaller degrees is contained in the larger ansatz, so a single generous
/// shape per order is enough; the shorter list keeps the exact linear
/// algebra affordable.
fn degrees(order: usize) -> Vec<(usize, usize, usize)> {
    if order == 1 {
        vec![(1, 1, 6), (2, 2, 6), (3, 2, 8), (3, 3, 8)]
    } else {
        vec![(2, 2, 8), (3, 2, 8), (3, 3, 8), (3, 3, 10)]
    }
}

/// The polynomial `Σ p_i n^i` as a term.
fn qpoly_node(
    graph: &mut Graph,
    p: &QPoly,
    n: NodeId,
) -> NodeId {
    let mut terms = Vec::new();
    for (i, coefficient) in p.iter().enumerate() {
        if coefficient.is_zero() {
            continue;
        }
        let c = graph.num(Number::rat(coefficient.clone()));
        let term = if i == 0 {
            c
        } else {
            let power = pow_of(graph, n, i64::try_from(i).unwrap_or(0));
            mul_of(graph, &[c, power])
        };
        terms.push(term);
    }
    add_of(graph, &terms)
}

/// A rational function of `n` as a term.
fn frac_node(
    graph: &mut Graph,
    f: &Frac,
    n: NodeId,
) -> NodeId {
    let numerator = qpoly_node(graph, &f.n, n);
    let denominator = qpoly_node(graph, &f.d, n);
    let minus_one = graph.int(-1);
    let inverse = graph.node(core::POW, &[denominator, minus_one]);
    let raw = graph.node(core::MUL, &[numerator, inverse]);
    match best(graph, raw) {
        | Some(simplified) => simplified,
        | None => raw,
    }
}

// ---------------------------------------------------------------------
// The exact search over Q(n)
// ---------------------------------------------------------------------
//
// The polynomial ansatz above fixes a degree in `n` for the recurrence and
// the certificate in advance. The coefficients of a second-order recurrence
// (Dixon, Apéry) and the numerator of the certificate are rational functions
// of `n` of degrees that are not known beforehand, and solving the linear
// system symbolically over `Q(n)` makes the coefficients explode.
//
// The search below works the other way round: `n` is specialised to an
// integer, which turns every ratio into a rational function of `k` alone, the
// linear system into one over `Q` that is small and cheap, and its (one
// dimensional) solution space into a vector that is normalised by
// `p_J = 1`. Repeating this for consecutive integers gives samples of the
// solution as a function of `n`, and each coordinate is recovered by
// rational function reconstruction (Cauchy interpolation through the
// extended Euclidean algorithm) with held-out points as a check. The result
// is the recurrence and a certificate numerator over `Q(n)`: poles in `n`
// are allowed. The caller verifies the final identity independently.

/// Highest recurrence order of the exact search.
const MAX_EXACT_ORDER: usize = 3;
/// The first integer `n` at which the term is specialised.
const FIRST_SAMPLE: i64 = 23;
/// Numbers of samples after which a reconstruction is attempted.
const SAMPLE_STAGES: [i64; 4] = [12, 20, 32, 48];
/// Largest degree in `k` of the certificate numerator.
const MAX_CERTIFICATE_DEGREE: usize = 48;

/// A polynomial in `n` and `k` over `Q`.
struct Bivariate {
    /// `(exponent of n, exponent of k, coefficient)`.
    terms: Vec<(u32, u32, BigRational)>,
}

impl Bivariate {
    fn from_node(
        graph: &mut Graph,
        node: NodeId,
        n: NodeId,
        k: NodeId,
    ) -> Option<Self> {
        let mut gens = Gens::default();
        let gn = gens.index(graph, n);
        let gk = gens.index(graph, k);
        let poly = from_term(
            graph,
            &mut gens,
            node,
            Limits {
                terms: 200_000,
                exponent: 256,
            },
        )?;
        let mut terms = Vec::with_capacity(poly.len());
        for (monomial, coefficient) in poly.terms() {
            let (mut n_power, mut k_power) = (0, 0);
            for &(generator, exponent) in monomial {
                if generator == gn {
                    n_power = exponent;
                } else if generator == gk {
                    k_power = exponent;
                } else {
                    return None;
                }
            }
            terms.push((n_power, k_power, coefficient.to_rational()?));
        }
        Some(Self { terms })
    }

    /// The polynomial in `k` at `n = n0`.
    fn at_n(
        &self,
        n0: &BigRational,
    ) -> QPoly {
        let mut powers = vec![BigRational::one()];
        let mut out: QPoly = Vec::new();
        for (n_power, k_power, coefficient) in &self.terms {
            let n_power = usize::try_from(*n_power).unwrap_or(0);
            let k_power = usize::try_from(*k_power).unwrap_or(0);
            while powers.len() <= n_power {
                let next = powers.last().cloned().unwrap_or_else(BigRational::one) * n0;
                powers.push(next);
            }
            if out.len() <= k_power {
                out.resize(k_power.saturating_add(1), BigRational::zero());
            }
            if let (Some(slot), Some(power)) = (out.get_mut(k_power), powers.get(n_power)) {
                *slot += coefficient * power;
            }
        }
        while out.last().is_some_and(Zero::is_zero) {
            out.pop();
        }
        out
    }

    fn k_degree(&self) -> usize {
        self.terms
            .iter()
            .map(|t| usize::try_from(t.1).unwrap_or(0))
            .max()
            .unwrap_or(0)
    }
}

/// The ratios of a term as polynomials in `n` and `k`:
/// `F(n, k+1)/F(n, k) = a/b` and `F(n+j, k)/F(n, k) = nums[j-1]/dens[j-1]`.
struct Problem {
    a: Bivariate,
    b: Bivariate,
    nums: Vec<Bivariate>,
    dens: Vec<Bivariate>,
}

/// `p(k + 1)`.
fn shift_k(p: &[BigRational]) -> QPoly {
    let step = vec![BigRational::one(), BigRational::one()];
    let mut out = QPoly::new();
    for c in p.iter().rev() {
        out = univariate::add(&univariate::mul(&out, &step), std::slice::from_ref(c));
    }
    out
}

fn product(polys: impl IntoIterator<Item = QPoly>) -> QPoly {
    polys
        .into_iter()
        .fold(vec![BigRational::one()], |acc, p| univariate::mul(&acc, &p))
}

/// The polynomials of the linear system at one `n`: column `j` of the
/// recurrence part is `P_j`, column `m` of the certificate part is `V_m`.
struct Slice {
    p: Vec<QPoly>,
    v: Vec<QPoly>,
}

/// The indices `i ≥ 1` of the denominators `D_i` whose product `Q` is a
/// common multiple of `D_1 .. D_order`: a `D_i` that divides the product of
/// the others (checked at two sample values of `n`) is dropped, which keeps
/// the degrees of the system low (for `C(n,k)^3` all of `D_1` hides in `D_2`).
fn denominator_set(
    problem: &Problem,
    order: usize,
) -> Vec<usize> {
    let mut kept: Vec<usize> = (1..=order).collect();
    for i in 1..=order {
        let others: Vec<usize> = kept.iter().copied().filter(|&l| l != i).collect();
        let divides = [FIRST_SAMPLE, FIRST_SAMPLE + 1].iter().all(|&n0| {
            let n0 = rat(n0);
            let rest = product(
                others
                    .iter()
                    .filter_map(|&l| problem.dens.get(l - 1))
                    .map(|d| d.at_n(&n0)),
            );
            problem
                .dens
                .get(i - 1)
                .map(|d| d.at_n(&n0))
                .is_some_and(|d| {
                    !d.is_empty()
                        && univariate::divrem(&rest, &d).is_some_and(|(_, r)| r.is_empty())
                })
        });
        if divides && kept.len() > 1 && kept.contains(&i) {
            kept = others;
        }
    }
    kept
}

/// The telescoping identity multiplied out for a fixed `n`:
/// `Σ_j p_j P_j - Σ_m u_m V_m = 0`, with
/// `P_j = N_j (Q/D_j) B B(k+1) Q(k+1)`,
/// `V_m = (k+1)^m A Q - k^m B(k+1) Q(k+1)` and `Q = Π_{i ∈ q} D_i`,
/// the certificate being `Σ_m u_m k^m / (B Q)`.
fn slice_at(
    problem: &Problem,
    order: usize,
    q_set: &[usize],
    certificate_degree: usize,
    n0: i64,
) -> Option<Slice> {
    let n0 = rat(n0);
    let a = problem.a.at_n(&n0);
    let b = problem.b.at_n(&n0);
    let mut nums = vec![vec![BigRational::one()]];
    let mut dens = vec![vec![BigRational::one()]];
    for j in 0..order {
        nums.push(problem.nums.get(j)?.at_n(&n0));
        dens.push(problem.dens.get(j)?.at_n(&n0));
    }
    if b.is_empty() || dens.iter().any(Vec::is_empty) {
        return None;
    }
    let q = product(q_set.iter().filter_map(|&i| dens.get(i)).cloned());
    let (b1, q1) = (shift_k(&b), shift_k(&q));
    let common = product([b, b1.clone(), q1.clone()]);
    let mut p = Vec::with_capacity(order.saturating_add(1));
    for (numerator, denominator) in nums.iter().zip(&dens) {
        let (cofactor, rest) = univariate::divrem(&q, denominator)?;
        if !rest.is_empty() {
            return None;
        }
        p.push(univariate::mul(
            &univariate::mul(numerator, &cofactor),
            &common,
        ));
    }
    let aq = univariate::mul(&a, &q);
    let bq = univariate::mul(&b1, &q1);
    let step = vec![BigRational::one(), BigRational::one()];
    let mut v = Vec::with_capacity(certificate_degree.saturating_add(1));
    let mut shifted_power = vec![BigRational::one()];
    for m in 0..=certificate_degree {
        let mut monomial = vec![BigRational::zero(); m];
        monomial.push(BigRational::one());
        let upper = univariate::mul(&shifted_power, &aq);
        let lower = univariate::mul(&monomial, &bq);
        v.push(univariate::sub(&upper, &lower));
        shifted_power = univariate::mul(&shifted_power, &step);
    }
    Some(Slice { p, v })
}

/// What the specialised system says about a recurrence of the given shape.
enum Null {
    /// No recurrence.
    Absent,
    /// The solution has `p_J = 0`.
    Degenerate,
    /// A single solution, normalised by `p_J = 1`: `(p_0..p_J, u_0..u_d)`.
    Unique(Vec<BigRational>),
    /// More than one solution.
    Many,
}

fn null_vector(
    slice: &Slice,
    order: usize,
    certificate_degree: usize,
) -> Null {
    let columns = order.saturating_add(certificate_degree).saturating_add(2);
    let height = slice
        .p
        .iter()
        .take(order.saturating_add(1))
        .chain(slice.v.iter().take(certificate_degree.saturating_add(1)))
        .map(Vec::len)
        .max()
        .unwrap_or(0);
    let mut rows = vec![vec![BigRational::zero(); columns.saturating_add(1)]; height];
    for (column, poly) in slice.p.iter().take(order.saturating_add(1)).enumerate() {
        for (s, c) in poly.iter().enumerate() {
            if let Some(slot) = rows.get_mut(s).and_then(|r| r.get_mut(column)) {
                *slot = c.clone();
            }
        }
    }
    for (m, poly) in slice
        .v
        .iter()
        .take(certificate_degree.saturating_add(1))
        .enumerate()
    {
        let column = m.saturating_add(order).saturating_add(1);
        for (s, c) in poly.iter().enumerate() {
            if let Some(slot) = rows.get_mut(s).and_then(|r| r.get_mut(column)) {
                *slot = -c.clone();
            }
        }
    }
    let Some(system) = solve_system(rows, columns) else {
        return Null::Absent;
    };
    match system.free.as_slice() {
        | [] => Null::Absent,
        | [free] => {
            let x = system.solution(Some((*free, BigRational::one())));
            let recurrence = x.get(..=order).unwrap_or_default();
            if recurrence.iter().all(Zero::is_zero) {
                return Null::Absent;
            }
            let Some(lead) = x.get(order).filter(|v| !v.is_zero()).cloned() else {
                return Null::Degenerate;
            };
            Null::Unique(x.iter().map(|v| v / &lead).collect())
        },
        | _ => Null::Many,
    }
}

/// The polynomial of degree `< xs.len()` through the points (Newton form).
fn interpolate(
    xs: &[BigRational],
    ys: &[BigRational],
) -> QPoly {
    let mut c = ys.to_vec();
    for j in 1..xs.len() {
        for i in (j..xs.len()).rev() {
            let step = &xs[i] - &xs[i - j];
            c[i] = (&c[i] - &c[i - 1]) / step;
        }
    }
    let Some(last) = c.last() else {
        return Vec::new();
    };
    let mut p = univariate::add(&[], std::slice::from_ref(last));
    for i in (0..xs.len().saturating_sub(1)).rev() {
        let factor = vec![-xs[i].clone(), BigRational::one()];
        p = univariate::add(&univariate::mul(&p, &factor), std::slice::from_ref(&c[i]));
    }
    p
}

/// The rational function `a/b` through the points, when its degrees leave
/// room for held-out points: the first points are fitted by the extended
/// Euclidean algorithm (every remainder is a candidate; the one of the
/// smallest total degree wins), the last three must agree.
fn rational_through(
    xs: &[BigRational],
    ys: &[BigRational],
) -> Option<(QPoly, QPoly)> {
    let m = xs.len();
    if m < 8 {
        return None;
    }
    let fit = m - 3;
    let (fit_x, fit_y) = (xs.get(..fit)?, ys.get(..fit)?);
    let verify = |num: &QPoly, den: &QPoly| {
        xs.iter().zip(ys).skip(fit).all(|(x, y)| {
            let d = univariate::eval(den, x);
            !d.is_zero() && univariate::eval(num, x) == y * d
        })
    };
    if fit_y.iter().all(Zero::is_zero) {
        let zero = (Vec::new(), vec![BigRational::one()]);
        return verify(&zero.0, &zero.1).then_some(zero);
    }
    let mut modulus = vec![BigRational::one()];
    for x in fit_x {
        modulus = univariate::mul(&modulus, &[-x.clone(), BigRational::one()]);
    }
    let (mut r_prev, mut r) = (modulus, interpolate(fit_x, fit_y));
    let (mut t_prev, mut t): (QPoly, QPoly) = (Vec::new(), vec![BigRational::one()]);
    let mut best: Option<(usize, QPoly, QPoly)> = None;
    while !r.is_empty() {
        // `score` is the total degree plus two.
        let score = r.len() + t.len();
        let admissible = !t.is_empty()
            && score < fit
            && fit_x.iter().all(|x| !univariate::eval(&t, x).is_zero());
        if admissible && best.as_ref().is_none_or(|b| score < b.0) {
            best = Some((score, r.clone(), t.clone()));
        }
        let (quotient, remainder) = univariate::divrem(&r_prev, &r)?;
        let t_next = univariate::sub(&t_prev, &univariate::mul(&quotient, &t));
        r_prev = r;
        r = remainder;
        t_prev = t;
        t = t_next;
    }
    let (_, num, den) = best?;
    let g = univariate::gcd(&num, &den);
    let (num, _) = univariate::divrem(&num, &g)?;
    let (den, _) = univariate::divrem(&den, &g)?;
    let lead = den.last()?.clone();
    let num: QPoly = num.iter().map(|c| c / &lead).collect();
    let den: QPoly = den.iter().map(|c| c / &lead).collect();
    verify(&num, &den).then_some((num, den))
}

/// The recurrence and certificate numerator found over `Q(n)`.
struct Exact {
    /// `p_0..p_J`, coprime polynomials in `n` with integer coefficients.
    p: Vec<QPoly>,
    /// The coefficients `u_m(n)` of `k^m` in the certificate numerator.
    u: Vec<Frac>,
    /// The denominators `D_i` whose product `Q` the certificate carries.
    q_set: Vec<usize>,
}

fn quotient_of(
    a: &[BigRational],
    b: &[BigRational],
) -> Option<QPoly> {
    univariate::divrem(a, b).map(|(q, _)| q)
}

/// The polynomial through the points when its degree leaves room for three
/// held-out points.
fn polynomial_through(
    xs: &[BigRational],
    ys: &[BigRational],
) -> Option<QPoly> {
    let fit = xs.len().checked_sub(3)?;
    let p = interpolate(xs.get(..fit)?, ys.get(..fit)?);
    let holds = xs
        .iter()
        .zip(ys)
        .skip(fit)
        .all(|(x, y)| univariate::eval(&p, x) == *y);
    (holds && p.len().saturating_add(2) <= fit).then_some(p)
}

/// Combines the samples into the normalised result: polynomial `p_j`
/// without common factor, integer and primitive, and the certificate
/// coefficients over their common denominator, scaled in step.
///
/// The recurrence coordinates are reconstructed one by one. The
/// certificate coordinates share one denominator, found from a fixed
/// integer combination of them, so that each of them is then a plain
/// polynomial interpolation.
fn assemble(
    xs: &[BigRational],
    vectors: &[Vec<BigRational>],
    order: usize,
    q_set: &[usize],
) -> Option<Exact> {
    let columns = vectors.first()?.len();
    let column = |c: usize| -> Vec<BigRational> {
        vectors
            .iter()
            .map(|v| v.get(c).cloned().unwrap_or_else(BigRational::zero))
            .collect()
    };
    let mut recurrence = Vec::with_capacity(order.saturating_add(1));
    for c in 0..=order {
        recurrence.push(rational_through(xs, &column(c))?);
    }
    let mut combined = vec![BigRational::zero(); xs.len()];
    for c in order.saturating_add(1)..columns {
        let weight =
            rat(i64::try_from(c.saturating_mul(c).saturating_mul(7).saturating_add(3)).ok()?);
        for (slot, y) in combined.iter_mut().zip(column(c)) {
            *slot += &weight * y;
        }
    }
    let (_, pole) = rational_through(xs, &combined)?;
    let pole_values: Vec<BigRational> = xs.iter().map(|x| univariate::eval(&pole, x)).collect();
    let mut u = Vec::with_capacity(columns.saturating_sub(order.saturating_add(1)));
    for c in order.saturating_add(1)..columns {
        let scaled: Vec<BigRational> = column(c)
            .iter()
            .zip(&pole_values)
            .map(|(y, d)| y * d)
            .collect();
        let numerator = polynomial_through(xs, &scaled)?;
        u.push(
            Frac {
                n: numerator,
                d: pole.clone(),
            }
            .reduce(),
        );
    }
    let mut common = vec![BigRational::one()];
    for (_, den) in &recurrence {
        let g = univariate::gcd(&common, den);
        common = univariate::mul(&common, &quotient_of(den, &g)?);
    }
    let mut scaled = Vec::with_capacity(recurrence.len());
    for (num, den) in &recurrence {
        scaled.push(univariate::mul(num, &quotient_of(&common, den)?));
    }
    let g = scaled
        .iter()
        .filter(|p| !p.is_empty())
        .fold(Vec::new(), |acc, p| univariate::gcd(&acc, p));
    let mut p = Vec::with_capacity(scaled.len());
    for poly in &scaled {
        p.push(quotient_of(poly, &g)?);
    }
    // A constant that makes the coefficients coprime integers.
    let mut denominator = BigInt::one();
    for c in p.iter().flatten() {
        denominator = num_integer::Integer::lcm(&denominator, c.denom());
    }
    let scale = BigRational::from_integer(denominator);
    let mut content = BigInt::zero();
    for c in p.iter().flatten() {
        content = gcd(&content, (c * &scale).numer());
    }
    if content.is_zero() {
        return None;
    }
    let mut constant = scale / BigRational::from_integer(content);
    if p.get(order)
        .and_then(|q| q.last())
        .is_some_and(Signed::is_negative)
    {
        constant = -constant;
    }
    let p: Vec<QPoly> = p
        .iter()
        .map(|q| q.iter().map(|c| c * &constant).collect())
        .collect();
    let factor = Frac { n: common, d: g }
        .reduce()
        .mul(&Frac::poly(vec![constant]));
    let u = u.iter().map(|f| f.mul(&factor)).collect();
    Some(Exact {
        p,
        u,
        q_set: q_set.to_vec(),
    })
}

/// Samples the specialised system at consecutive `n` and reconstructs.
fn reconstruct(
    problem: &Problem,
    order: usize,
    q_set: &[usize],
    degree: usize,
) -> Option<Exact> {
    let mut xs = Vec::new();
    let mut vectors = Vec::new();
    let mut next = 0_i64;
    for stage in SAMPLE_STAGES {
        while next < stage {
            let n0 = FIRST_SAMPLE.checked_add(next)?;
            next = next.checked_add(1)?;
            let Some(slice) = slice_at(problem, order, q_set, degree, n0) else {
                continue;
            };
            if let Null::Unique(x) = null_vector(&slice, order, degree) {
                xs.push(rat(n0));
                vectors.push(x);
            }
        }
        if let Some(found) = assemble(&xs, &vectors, order, q_set) {
            return Some(found);
        }
    }
    None
}

/// The smallest order, and for it the smallest certificate degree, at which
/// the specialised system has a solution; then the reconstruction over
/// `Q(n)`.
fn exact_search(problem: &Problem) -> Option<Exact> {
    for order in 1..=MAX_EXACT_ORDER.min(problem.nums.len()) {
        let q_set = denominator_set(problem, order);
        let k_degree = problem.b.k_degree()
            + problem.a.k_degree()
            + q_set
                .iter()
                .filter_map(|&i| problem.dens.get(i - 1))
                .map(Bivariate::k_degree)
                .sum::<usize>();
        let limit = k_degree.saturating_add(4).min(MAX_CERTIFICATE_DEGREE);
        for degree in 0..=limit {
            let slice = slice_at(problem, order, &q_set, degree, FIRST_SAMPLE)?;
            match null_vector(&slice, order, degree) {
                | Null::Absent | Null::Degenerate => {},
                | Null::Unique(_) | Null::Many => {
                    if let Some(found) = reconstruct(problem, order, &q_set, degree) {
                        return Some(found);
                    }
                    break;
                },
            }
        }
    }
    None
}

/// Creative telescoping over `Q(n)` for a parameter-free term.
fn telescoping_over_qn(
    cx: &mut Cx<'_>,
    f: NodeId,
    factors: &Factored,
    symbols: Symbols,
) -> Option<Telescoping> {
    let kq = k_ratio(cx, factors, symbols)?;
    let a = Bivariate::from_node(cx.graph, kq.num, symbols.n, symbols.k)?;
    let b = Bivariate::from_node(cx.graph, kq.den, symbols.n, symbols.k)?;
    let mut nums = Vec::new();
    let mut dens = Vec::new();
    let mut den_nodes = Vec::new();
    for j in 1..=MAX_EXACT_ORDER {
        let ratio = n_ratio(cx, factors, i64::try_from(j).ok()?, symbols)?;
        let (Some(num), Some(den)) = (
            Bivariate::from_node(cx.graph, ratio.num, symbols.n, symbols.k),
            Bivariate::from_node(cx.graph, ratio.den, symbols.n, symbols.k),
        ) else {
            break;
        };
        nums.push(num);
        dens.push(den);
        den_nodes.push(ratio.den);
    }
    let problem = Problem { a, b, nums, dens };
    let exact = exact_search(&problem)?;
    let coefficients: Vec<NodeId> = exact
        .p
        .iter()
        .map(|p| qpoly_node(cx.graph, p, symbols.n))
        .collect();
    let mut numerator_terms = Vec::new();
    for (m, u) in exact.u.iter().enumerate() {
        if u.is_zero() {
            continue;
        }
        let value = frac_node(cx.graph, u, symbols.n);
        let power = pow_of(cx.graph, symbols.k, i64::try_from(m).ok()?);
        numerator_terms.push(mul_of(cx.graph, &[value, power]));
    }
    let numerator = add_of(cx.graph, &numerator_terms);
    let mut denominator_factors = vec![kq.den];
    denominator_factors.extend(exact.q_set.iter().filter_map(|&i| den_nodes.get(i - 1)));
    let denominator = mul_of(cx.graph, &denominator_factors);
    let minus_one = cx.graph.int(-1);
    let inverse = cx.graph.node(core::POW, &[denominator, minus_one]);
    let raw = cx.graph.node(core::MUL, &[numerator, inverse]);
    let certificate = cx.simplify(raw);
    let found = Telescoping {
        coefficients,
        certificate,
    };
    verified(cx, f, &found, symbols, &[]).then_some(found)
}

/// Zeilberger's creative telescoping for `F(n, k)` in the index `k`.
pub(super) fn creative_telescoping(
    cx: &mut Cx<'_>,
    f: NodeId,
    k: NodeId,
    n: NodeId,
) -> Option<Telescoping> {
    let ns = cx.graph.symbol_of(n)?;
    let ks = cx.graph.symbol_of(k)?;
    let symbols = Symbols { n, k, ns, ks };
    let factors = factorize(cx, f, symbols)?;
    let free = cx.graph.free_symbols(cx.graph.find(f)).to_vec();
    let parameters: Vec<NodeId> = free
        .iter()
        .filter(|&&s| s != ns && s != ks)
        .map(|&s| cx.graph.symbol_node(s))
        .collect();
    if parameters.len() > 2 {
        return None;
    }
    if parameters.is_empty() {
        // Without parameters the exact search over `Q(n)` covers every
        // order, including certificates with poles in `n`.
        return telescoping_over_qn(cx, f, &factors, symbols);
    }
    let kq = k_ratio(cx, &factors, symbols)?;
    for order in 1..=2_usize {
        let ratios: Vec<Quotient> = (0..=order)
            .map(|j| n_ratio(cx, &factors, i64::try_from(j).unwrap_or(0), symbols))
            .collect::<Option<Vec<_>>>()?;
        for (n_degree, u_n_degree, u_k_degree) in degrees(order) {
            for parameter_degree in 0..=MAX_PARAMETER_DEGREE {
                let Some(ansatz) = build_ansatz(
                    cx.graph,
                    &parameters,
                    order,
                    n_degree,
                    u_n_degree,
                    u_k_degree,
                    parameter_degree,
                ) else {
                    continue;
                };
                let identity = identity_node(cx, &ansatz, &ratios, kq, symbols);
                let Some(rows) = equations(cx, identity, &ansatz, symbols) else {
                    continue;
                };
                let Some(system) = solve_system(rows, ansatz.nodes.len()) else {
                    continue;
                };
                let p_columns: Vec<usize> = ansatz
                    .p_index
                    .iter()
                    .flat_map(|in_n| in_n.iter().flatten().copied())
                    .collect();
                let order_columns: Vec<usize> = ansatz
                    .p_index
                    .get(order)
                    .into_iter()
                    .flatten()
                    .flatten()
                    .copied()
                    .collect();
                let mut candidates = vec![system.solution(None)];
                for &column in &system.free {
                    candidates.push(system.solution(Some((column, BigRational::one()))));
                }
                for mut values in candidates {
                    if !p_columns
                        .iter()
                        .any(|&c| values.get(c).is_some_and(|v| !v.is_zero()))
                    {
                        continue;
                    }
                    if order_columns
                        .iter()
                        .all(|&c| values.get(c).is_none_or(Zero::is_zero))
                    {
                        continue;
                    }
                    normalize(&mut values, &p_columns, &order_columns);
                    let coefficients: Vec<NodeId> = (0..=order)
                        .map(|j| ansatz.numeric_polynomial(cx.graph, &values, Some(j), n, k))
                        .collect();
                    let u = ansatz.numeric_polynomial(cx.graph, &values, None, n, k);
                    let minus_one = cx.graph.int(-1);
                    let common =
                        mul_of(cx.graph, &ratios.iter().map(|r| r.den).collect::<Vec<_>>());
                    let denominator = mul_of(cx.graph, &[kq.den, common]);
                    let inverse = cx.graph.node(core::POW, &[denominator, minus_one]);
                    let raw = cx.graph.node(core::MUL, &[u, inverse]);
                    let certificate = cx.simplify(raw);
                    let found = Telescoping {
                        coefficients,
                        certificate,
                    };
                    if verified(cx, f, &found, symbols, &parameters) {
                        return Some(found);
                    }
                }
            }
        }
    }
    None
}

/// Numeric check of `Σ_j p_j F(n+j,k) = G(n,k+1) - G(n,k)` at sample
/// points.
fn verified(
    cx: &mut Cx<'_>,
    f: NodeId,
    found: &Telescoping,
    symbols: Symbols,
    parameters: &[NodeId],
) -> bool {
    let graph = &mut *cx.graph;
    let one = graph.int(1);
    let mut lhs_terms = Vec::new();
    for (j, coefficient) in found.coefficients.iter().enumerate() {
        let step = graph.int(i64::try_from(j).unwrap_or(0));
        let lifted = graph.node(core::ADD, &[symbols.n, step]);
        let term = graph.substitute(f, symbols.n, lifted);
        lhs_terms.push(mul_of(graph, &[*coefficient, term]));
    }
    let lhs = add_of(graph, &lhs_terms);
    let next = graph.node(core::ADD, &[symbols.k, one]);
    let next_value = graph.substitute(f, symbols.k, next);
    let certificate_next = graph.substitute(found.certificate, symbols.k, next);
    let upper_right = mul_of(graph, &[certificate_next, next_value]);
    let lower_right = mul_of(graph, &[found.certificate, f]);
    let right = sub_of(graph, upper_right, lower_right);
    let residual = sub_of(graph, lhs, right);
    let mut checked = 0;
    for (n_value, k_value) in [(5_i64, 2_i64), (6, 3), (7, 1), (8, 4), (9, 2), (4, 1)] {
        let mut env = generic_env(graph, parameters);
        env.bind(symbols.ns, n_value as f64);
        env.bind(symbols.ks, k_value as f64);
        let (Some(value), Some(scale)) = (graph.eval(residual, &env), graph.eval(f, &env)) else {
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

fn generic_env(
    graph: &Graph,
    parameters: &[NodeId],
) -> Env {
    let mut env = Env::numeric(0.0);
    for (index, parameter) in parameters.iter().enumerate() {
        if let Some(symbol) = graph.symbol_of(*parameter) {
            let base = u32::try_from(index).unwrap_or(0) % 5;
            let extra = symbol.raw() % 3;
            env.bind(symbol, 2.0 + f64::from(base) + f64::from(extra));
        }
    }
    env
}

/// Kernel for `zeilberger(F, k, n)`.
pub(super) struct Zeilberger {
    pub(super) op: OpId,
}

impl Kernel for Zeilberger {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let &[f, k, n] = cx.graph.children(node) else {
            return Outcome::Pass;
        };
        let Some(found) = creative_telescoping(cx, f, k, n) else {
            return Outcome::Pass;
        };
        let recurrence = cx.graph.node(core::LIST, &found.coefficients);
        let answer = cx.graph.node(core::LIST, &[recurrence, found.certificate]);
        Outcome::Equal(answer)
    }
}

/// The closed form of `Σ_{k=-oo}^{oo} F(n, k)`, if one is found.
pub(super) fn sum_over_integers(
    cx: &mut Cx<'_>,
    f: NodeId,
    k: NodeId,
    n: NodeId,
) -> Option<NodeId> {
    let ns = cx.graph.symbol_of(n)?;
    let ks = cx.graph.symbol_of(k)?;
    let symbols = Symbols { n, k, ns, ks };
    let free = cx.graph.free_symbols(cx.graph.find(f)).to_vec();
    let parameters: Vec<NodeId> = free
        .iter()
        .filter(|&&s| s != ns && s != ks)
        .map(|&s| cx.graph.symbol_node(s))
        .collect();
    if parameters.len() > 2 {
        return None;
    }
    let found = creative_telescoping(cx, f, k, n)?;
    if found.coefficients.len() != 2 {
        return None;
    }
    let gamma = match parameters.as_slice() {
        | [] => None,
        | [only] => Some(*only),
        | _ => return None,
    };
    let coefficients: Vec<Vec<Frac>> = found
        .coefficients
        .iter()
        .map(|&p| polynomial_frac(cx, p, symbols.n, gamma))
        .collect::<Option<Vec<_>>>()?;
    let (start_n, start) = initial_value(cx, f, symbols, &parameters)?;
    let length = usize::try_from(start_n).ok()?.checked_add(1)?;
    let mut values = vec![Frac::zero(); length];
    if let Some(slot) = values.get_mut(usize::try_from(start_n).ok()?) {
        *slot = Frac::poly(vec![start]);
    }
    let closed = recurrence::solve_first_order(cx.graph, &coefficients, &values, symbols.n, gamma)?;
    let closed = cx.simplify(closed);
    verify_closed_form(cx, f, closed, symbols, &parameters, start_n)?;
    Some(closed)
}

/// The coefficients of the polynomial `p` in `n` as elements of `Q(γ)`.
fn polynomial_frac(
    cx: &mut Cx<'_>,
    p: NodeId,
    n: NodeId,
    gamma: Option<NodeId>,
) -> Option<Vec<Frac>> {
    let graph = &mut *cx.graph;
    let mut gens = Gens::default();
    let gn = gens.index(graph, n);
    let gg = gamma.map(|g| gens.index(graph, g));
    let poly = from_term(graph, &mut gens, p, Limits::default())?;
    for generator in poly.support() {
        if generator != gn && Some(generator) != gg {
            return None;
        }
    }
    let mut out = Vec::new();
    for coefficient in poly.coefficients_in(gn) {
        let q: QPoly = match gg {
            | Some(gg) => {
                coefficient
                    .univariate_in(gg)?
                    .iter()
                    .map(Number::to_rational)
                    .collect::<Option<QPoly>>()?
            },
            | None => vec![coefficient.as_constant()?.to_rational()?],
        };
        out.push(Frac::poly(q));
    }
    Some(out)
}

/// An initial value `S(n₀)` of the sum, computed directly.
fn initial_value(
    cx: &mut Cx<'_>,
    f: NodeId,
    symbols: Symbols,
    parameters: &[NodeId],
) -> Option<(i64, BigRational)> {
    let factors = factorize(cx, f, symbols)?;
    let env = generic_env(cx.graph, parameters);
    for n_value in 0..=4_i64 {
        let Some((low, high)) = support_bounds(cx.graph, &factors, n_value, &env) else {
            continue;
        };
        let (Some(low), Some(high)) = (low, high) else {
            continue;
        };
        let Some(value) = direct_sum(cx.graph, f, symbols, n_value, low, high, &env) else {
            continue;
        };
        if value.is_finite() && value.abs() < 1e12 {
            let Some(rational) = rationalize(value) else {
                continue;
            };
            if rational.is_zero() {
                continue;
            }
            return Some((n_value, rational));
        }
    }
    None
}

/// The finite support of the term in `k` at `n = n_value`, when the
/// denominator Gamma factors bound it.
fn support_bounds(
    graph: &Graph,
    factors: &Factored,
    n_value: i64,
    env: &Env,
) -> Option<(Option<i64>, Option<i64>)> {
    let mut low = None;
    let mut high = None;
    for factor in &factors.gammas {
        if factor.exponent >= 0 || factor.k == 0 {
            continue;
        }
        let offset = graph.eval(factor.offset, env)?;
        if !offset.is_finite() {
            return None;
        }
        #[allow(clippy::cast_precision_loss)]
        let line = factor.n as f64 * n_value as f64 + offset;
        #[allow(clippy::cast_precision_loss)]
        let b = factor.k as f64;
        let root = -line / b;
        if !root.is_finite() || root.abs() > 1e12 {
            return None;
        }
        if factor.k > 0 {
            let bound = root.floor();
            #[allow(clippy::cast_possible_truncation)]
            let bound = bound as i64;
            low = Some(low.map_or_else(
                || bound.saturating_add(1),
                |l: i64| l.max(bound.saturating_add(1)),
            ));
        } else {
            let bound = root.ceil();
            #[allow(clippy::cast_possible_truncation)]
            let bound = bound as i64;
            high = Some(high.map_or_else(
                || bound.saturating_sub(1),
                |h: i64| h.min(bound.saturating_sub(1)),
            ));
        }
    }
    Some((low, high))
}

/// `Σ_{k=low}^{high} F(n_value, k)` under the bindings of `env`.
fn direct_sum(
    graph: &mut Graph,
    f: NodeId,
    symbols: Symbols,
    n_value: i64,
    low: i64,
    high: i64,
    env: &Env,
) -> Option<f64> {
    if high < low || high.saturating_sub(low) > 400 {
        return None;
    }
    let n_node = graph.int(n_value);
    let term = graph.substitute(f, symbols.n, n_node);
    let mut total = 0.0;
    for k in low..=high {
        let k_node = graph.int(k);
        let at = graph.substitute(term, symbols.k, k_node);
        let value = graph.eval(at, env)?;
        if !value.is_finite() {
            return None;
        }
        total += value;
    }
    Some(total)
}

/// A small rational close to `value`.
fn rationalize(value: f64) -> Option<BigRational> {
    if !value.is_finite() {
        return None;
    }
    for denominator in 1..=1000_i64 {
        #[allow(clippy::cast_precision_loss)]
        let numerator = (value * denominator as f64).round();
        if !numerator.is_finite() || numerator.abs() > 1e12 {
            continue;
        }
        #[allow(clippy::cast_possible_truncation)]
        let numerator = numerator as i64;
        #[allow(clippy::cast_precision_loss)]
        let exact = numerator as f64 / denominator as f64;
        if (exact - value).abs() <= 1e-9 * value.abs().max(1.0) {
            return Some(BigRational::new(
                BigInt::from(numerator),
                BigInt::from(denominator),
            ));
        }
    }
    None
}

/// Checks the closed form against direct summation at further `n`.
fn verify_closed_form(
    cx: &mut Cx<'_>,
    f: NodeId,
    closed: NodeId,
    symbols: Symbols,
    parameters: &[NodeId],
    start_n: i64,
) -> Option<()> {
    let factors = factorize(cx, f, symbols)?;
    let env = generic_env(cx.graph, parameters);
    let mut checked = 0;
    for step in 1..=CHECKS {
        let n_value = start_n.checked_add(step)?;
        let Some((Some(low), Some(high))) = support_bounds(cx.graph, &factors, n_value, &env)
        else {
            continue;
        };
        let Some(want) = direct_sum(cx.graph, f, symbols, n_value, low, high, &env) else {
            continue;
        };
        let n_node = cx.graph.int(n_value);
        let got_node = cx.graph.substitute(closed, symbols.n, n_node);
        let mut bindings = env.clone();
        bindings.bind(symbols.ns, n_value as f64);
        let got = cx.graph.eval(got_node, &bindings)?;
        if !got.is_finite() || (got - want).abs() > 1e-7 * want.abs().max(1.0) {
            return None;
        }
        checked += 1;
    }
    (checked >= 2).then_some(())
}

/// Kernel for `definite_sum(F, k, n)`.
pub(super) struct DefiniteSum {
    pub(super) op: OpId,
}

impl Kernel for DefiniteSum {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let &[f, k, n] = cx.graph.children(node) else {
            return Outcome::Pass;
        };
        sum_over_integers(cx, f, k, n).map_or(Outcome::Pass, Outcome::Equal)
    }
}
