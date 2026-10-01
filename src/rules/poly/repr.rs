//! Sparse multivariate polynomials over arbitrary *generators*.
//!
//! A generator is any term the polynomial machinery does not look inside:
//! a symbol, `sin(x)`, `x^(1/2)`. Converting a term to a [`Poly`] expands
//! products and non-negative integer powers over sums; converting back
//! yields the expanded normal form. Everything polynomial in rssn —
//! expansion, collection, division, greatest common divisors,
//! factorisation, Gröbner bases — works on this representation.

use std::collections::BTreeMap;
use std::collections::HashMap;

use crate::graph::op::core;
use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::Number;

/// A power product: `(generator, exponent)` pairs sorted by generator, all
/// exponents positive. The empty monomial is `1`.
pub type Mono = Vec<(u32, u32)>;

/// The generators of a polynomial: index `i` stands for `nodes[i]`.
#[derive(Clone, Debug, Default)]
pub struct Gens {
    nodes: Vec<NodeId>,
}

impl Gens {
    /// Index of the generator equal to `node`, adding it if new. Nodes the
    /// graph knows to be equal share one generator.
    pub fn index(
        &mut self,
        graph: &Graph,
        node: NodeId,
    ) -> u32 {
        let position = self.nodes.iter().position(|&n| graph.same(n, node)).unwrap_or_else(|| {
            self.nodes.push(node);
            self.nodes.len().saturating_sub(1)
        });
        u32::try_from(position).unwrap_or(u32::MAX)
    }

    /// Index of an existing generator equal to `node`.
    #[must_use]
    pub fn find(
        &self,
        graph: &Graph,
        node: NodeId,
    ) -> Option<u32> {
        self.nodes.iter().position(|&n| graph.same(n, node)).and_then(|i| u32::try_from(i).ok())
    }

    /// The term generator `index` stands for.
    #[must_use]
    pub fn node(
        &self,
        index: u32,
    ) -> Option<NodeId> {
        self.nodes.get(index as usize).copied()
    }

    /// Number of generators.
    #[must_use]
    pub const fn len(&self) -> usize {
        self.nodes.len()
    }

    /// Whether there are no generators.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.nodes.is_empty()
    }
}

/// A polynomial: a map from monomials to non-zero coefficients.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Poly {
    terms: BTreeMap<Mono, Number>,
}

fn mono_mul(
    a: &Mono,
    b: &Mono,
) -> Mono {
    let mut out = Vec::with_capacity(a.len().saturating_add(b.len()));
    let (mut i, mut j) = (0, 0);
    loop {
        match (a.get(i), b.get(j)) {
            | (Some(&(ga, ea)), Some(&(gb, eb))) => {
                match ga.cmp(&gb) {
                    | std::cmp::Ordering::Equal => {
                        out.push((ga, ea.saturating_add(eb)));
                        i += 1;
                        j += 1;
                    },
                    | std::cmp::Ordering::Less => {
                        out.push((ga, ea));
                        i += 1;
                    },
                    | std::cmp::Ordering::Greater => {
                        out.push((gb, eb));
                        j += 1;
                    },
                }
            },
            | (Some(&x), None) => {
                out.push(x);
                i += 1;
            },
            | (None, Some(&y)) => {
                out.push(y);
                j += 1;
            },
            | (None, None) => return out,
        }
    }
}

/// Lexicographic monomial order: the exponent of the generator with the
/// smaller index decides first.
fn lex_cmp(
    a: &Mono,
    b: &Mono,
) -> std::cmp::Ordering {
    let (mut i, mut j) = (0, 0);
    loop {
        match (a.get(i), b.get(j)) {
            | (Some(&(ga, ea)), Some(&(gb, eb))) => {
                match ga.cmp(&gb) {
                    | std::cmp::Ordering::Equal => {
                        if ea != eb {
                            return ea.cmp(&eb);
                        }
                        i += 1;
                        j += 1;
                    },
                    | std::cmp::Ordering::Less => return std::cmp::Ordering::Greater,
                    | std::cmp::Ordering::Greater => return std::cmp::Ordering::Less,
                }
            },
            | (Some(_), None) => return std::cmp::Ordering::Greater,
            | (None, Some(_)) => return std::cmp::Ordering::Less,
            | (None, None) => return std::cmp::Ordering::Equal,
        }
    }
}

/// `a / b` if the monomial `b` divides `a`.
fn mono_div(
    a: &Mono,
    b: &Mono,
) -> Option<Mono> {
    let mut out = a.clone();
    for &(g, e) in b {
        let slot = out.iter_mut().find(|(h, _)| *h == g)?;
        slot.1 = slot.1.checked_sub(e)?;
    }
    out.retain(|&(_, e)| e > 0);
    Some(out)
}

impl Poly {
    /// The zero polynomial.
    #[must_use]
    pub fn zero() -> Self {
        Self::default()
    }

    /// A constant polynomial.
    #[must_use]
    pub fn constant(value: Number) -> Self {
        let mut p = Self::default();
        if !value.is_zero() {
            p.terms.insert(Vec::new(), value);
        }
        p
    }

    /// The polynomial consisting of generator `index` alone.
    #[must_use]
    pub fn generator(index: u32) -> Self {
        Self::monomial(vec![(index, 1)], Number::from(1))
    }

    /// A single term.
    #[must_use]
    pub fn monomial(
        mono: Mono,
        coeff: Number,
    ) -> Self {
        let mut p = Self::default();
        if !coeff.is_zero() {
            p.terms.insert(mono, coeff);
        }
        p
    }

    /// Whether this is the zero polynomial.
    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    /// Number of terms.
    #[must_use]
    pub fn len(&self) -> usize {
        self.terms.len()
    }

    /// Whether the polynomial has no terms (is zero).
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.terms.is_empty()
    }

    /// The terms, in an arbitrary but fixed order.
    pub fn terms(&self) -> impl Iterator<Item = (&Mono, &Number)> {
        self.terms.iter()
    }

    /// The value if the polynomial is constant.
    #[must_use]
    pub fn as_constant(&self) -> Option<Number> {
        match self.terms.len() {
            | 0 => Some(Number::from(0)),
            | 1 => self.terms.get(&Vec::new()).cloned(),
            | _ => None,
        }
    }

    /// Whether every coefficient is an exact number.
    #[must_use]
    pub fn is_exact(&self) -> bool {
        self.terms.values().all(Number::is_exact)
    }

    fn add_term(
        &mut self,
        mono: Mono,
        coeff: &Number,
    ) {
        match self.terms.get_mut(&mono) {
            | Some(slot) => {
                let sum = slot.add(coeff);
                if sum.is_zero() {
                    self.terms.remove(&mono);
                } else {
                    *slot = sum;
                }
            },
            | None => {
                if !coeff.is_zero() {
                    self.terms.insert(mono, coeff.clone());
                }
            },
        }
    }

    /// Sum.
    #[must_use]
    pub fn add(
        &self,
        other: &Self,
    ) -> Self {
        let mut out = self.clone();
        for (mono, coeff) in &other.terms {
            out.add_term(mono.clone(), coeff);
        }
        out
    }

    /// Additive inverse.
    #[must_use]
    pub fn neg(&self) -> Self {
        self.scale(&Number::from(-1))
    }

    /// Difference.
    #[must_use]
    pub fn sub(
        &self,
        other: &Self,
    ) -> Self {
        self.add(&other.neg())
    }

    /// Product with a number.
    #[must_use]
    pub fn scale(
        &self,
        factor: &Number,
    ) -> Self {
        if factor.is_zero() {
            return Self::zero();
        }
        Self { terms: self.terms.iter().map(|(m, c)| (m.clone(), c.mul(factor))).collect() }
    }

    /// Product. `None` when the result would exceed `cap` terms.
    #[must_use]
    pub fn mul(
        &self,
        other: &Self,
        cap: usize,
    ) -> Option<Self> {
        let mut out = Self::zero();
        for (ma, ca) in &self.terms {
            for (mb, cb) in &other.terms {
                out.add_term(mono_mul(ma, mb), &ca.mul(cb));
            }
            if out.terms.len() > cap {
                return None;
            }
        }
        Some(out)
    }

    /// Power by repeated squaring. `None` when a product exceeds `cap`.
    #[must_use]
    pub fn pow(
        &self,
        mut exponent: u32,
        cap: usize,
    ) -> Option<Self> {
        let mut result = Self::constant(Number::from(1));
        let mut base = self.clone();
        while exponent > 0 {
            if exponent & 1 == 1 {
                result = result.mul(&base, cap)?;
            }
            exponent >>= 1;
            if exponent > 0 {
                base = base.mul(&base, cap)?;
            }
        }
        Some(result)
    }

    /// The leading term in the lexicographic order of the generators.
    #[must_use]
    pub fn leading(&self) -> Option<(&Mono, &Number)> {
        self.terms.iter().max_by(|a, b| lex_cmp(a.0, b.0))
    }

    /// The largest monomial dividing every term (the empty monomial for
    /// zero).
    #[must_use]
    pub fn monomial_content(&self) -> Mono {
        let mut iter = self.terms.keys();
        let Some(first) = iter.next() else {
            return Vec::new();
        };
        let mut content = first.clone();
        for m in iter {
            content.retain_mut(|(g, e)| match m.iter().find(|(h, _)| h == g) {
                | Some(&(_, f)) => {
                    *e = (*e).min(f);
                    true
                },
                | None => false,
            });
        }
        content
    }

    /// Every term divided by the monomial `m`, which must divide them all.
    #[must_use]
    pub fn div_monomial(
        &self,
        m: &Mono,
    ) -> Self {
        Self { terms: self.terms.iter().filter_map(|(k, c)| Some((mono_div(k, m)?, c.clone()))).collect() }
    }

    /// The exact quotient `self / d`, or `None` if `d` does not divide
    /// `self` (or the division needs more than `cap` steps).
    #[must_use]
    pub fn div_exact(
        &self,
        d: &Self,
        cap: usize,
    ) -> Option<Self> {
        let (lead_mono, lead_coeff) = d.leading()?;
        let (lead_mono, inverse) = (lead_mono.clone(), lead_coeff.recip()?);
        let mut quotient = Self::zero();
        let mut rest = self.clone();
        for _ in 0..cap {
            let Some((mono, coeff)) = rest.leading() else {
                return Some(quotient);
            };
            let factor = mono_div(mono, &lead_mono)?;
            let c = coeff.mul(&inverse);
            let term = Self::monomial(factor, c);
            rest = rest.sub(&term.mul(d, cap.saturating_mul(4))?);
            quotient = quotient.add(&term);
        }
        rest.is_zero().then_some(quotient)
    }

    /// Degree in generator `gen`; zero for the zero polynomial.
    #[must_use]
    pub fn degree_in(
        &self,
        generator: u32,
    ) -> u32 {
        self.terms
            .keys()
            .filter_map(|m| m.iter().find(|&&(g, _)| g == generator).map(|&(_, e)| e))
            .max()
            .unwrap_or(0)
    }

    /// Total degree; zero for the zero polynomial.
    #[must_use]
    pub fn total_degree(&self) -> u32 {
        self.terms.keys().map(|m| m.iter().map(|&(_, e)| e).sum()).max().unwrap_or(0)
    }

    /// The generators that actually occur.
    #[must_use]
    pub fn support(&self) -> Vec<u32> {
        let mut gens: Vec<u32> = self.terms.keys().flat_map(|m| m.iter().map(|&(g, _)| g)).collect();
        gens.sort_unstable();
        gens.dedup();
        gens
    }

    /// Views the polynomial as univariate in `gen`: entry `k` is the
    /// coefficient of `gen^k`, itself a polynomial in the other generators.
    #[must_use]
    pub fn coefficients_in(
        &self,
        generator: u32,
    ) -> Vec<Self> {
        let mut out = vec![Self::zero(); self.degree_in(generator) as usize + 1];
        for (mono, coeff) in &self.terms {
            let power = mono.iter().find(|&&(g, _)| g == generator).map_or(0, |&(_, e)| e);
            let rest: Mono = mono.iter().copied().filter(|&(g, _)| g != generator).collect();
            if let Some(slot) = out.get_mut(power as usize) {
                slot.add_term(rest, coeff);
            }
        }
        out
    }

    /// The coefficients by degree if the polynomial involves no generator
    /// other than `gen`.
    #[must_use]
    pub fn univariate_in(
        &self,
        generator: u32,
    ) -> Option<Vec<Number>> {
        self.coefficients_in(generator).iter().map(Self::as_constant).collect()
    }

    /// Builds a univariate polynomial in `gen` from coefficients by degree.
    #[must_use]
    pub fn from_univariate(
        generator: u32,
        coeffs: &[Number],
    ) -> Self {
        let mut out = Self::zero();
        for (power, coeff) in coeffs.iter().enumerate() {
            let power = u32::try_from(power).unwrap_or(u32::MAX);
            let mono = if power == 0 { Vec::new() } else { vec![(generator, power)] };
            out.add_term(mono, coeff);
        }
        out
    }

    /// Partial derivative with respect to generator `gen`, treating the
    /// generators as independent.
    #[must_use]
    pub fn derivative(
        &self,
        generator: u32,
    ) -> Self {
        let mut out = Self::zero();
        for (mono, coeff) in &self.terms {
            let Some(&(_, power)) = mono.iter().find(|&&(g, _)| g == generator) else {
                continue;
            };
            let lowered: Mono = mono
                .iter()
                .filter_map(|&(g, e)| if g == generator { (e > 1).then_some((g, e - 1)) } else { Some((g, e)) })
                .collect();
            out.add_term(lowered, &coeff.mul(&Number::from(i64::from(power))));
        }
        out
    }
}

/// Limits for converting terms to polynomials.
#[derive(Copy, Clone, Debug)]
pub struct Limits {
    /// Maximum number of terms of any intermediate polynomial.
    pub terms: usize,
    /// Largest integer exponent that is expanded.
    pub exponent: u32,
}

impl Default for Limits {
    fn default() -> Self {
        Self { terms: 20_000, exponent: 256 }
    }
}

/// Converts the concrete term `node` to a polynomial over `gens`, expanding
/// products and non-negative integer powers.
///
/// Anything else becomes a generator — after its own arguments have been
/// expanded, so expansion is deep: `sin((x + 1)^2)` becomes the generator
/// `sin(x^2 + 2*x + 1)`.
///
/// Returns `None` when a limit is exceeded.
pub fn from_term(
    graph: &mut Graph,
    gens: &mut Gens,
    node: NodeId,
    limits: Limits,
) -> Option<Poly> {
    let mut memo: HashMap<NodeId, Poly> = HashMap::new();
    convert(graph, gens, node, limits, &mut memo, 0)
}

fn convert(
    graph: &mut Graph,
    gens: &mut Gens,
    node: NodeId,
    limits: Limits,
    memo: &mut HashMap<NodeId, Poly>,
    depth: usize,
) -> Option<Poly> {
    if let Some(done) = memo.get(&node) {
        return Some(done.clone());
    }
    if depth > 2_000 {
        return None;
    }
    let children = graph.children(node).to_vec();
    let poly = if let Some(n) = graph.number_of(node) {
        Poly::constant(n.clone())
    } else if graph.op(node) == core::ADD {
        let mut sum = Poly::zero();
        for child in children {
            sum = sum.add(&convert(graph, gens, child, limits, memo, depth + 1)?);
            if sum.len() > limits.terms {
                return None;
            }
        }
        sum
    } else if graph.op(node) == core::MUL {
        let mut product = Poly::constant(Number::from(1));
        for child in children {
            product = product.mul(&convert(graph, gens, child, limits, memo, depth + 1)?, limits.terms)?;
        }
        product
    } else {
        let power = match children.as_slice() {
            | [base, exp] if graph.op(node) == core::POW => graph
                .number_of(*exp)
                .and_then(Number::to_i64)
                .and_then(|e| u32::try_from(e).ok())
                .filter(|&e| e <= limits.exponent)
                .map(|e| (*base, e)),
            | _ => None,
        };
        match power {
            | Some((base, exponent)) => {
                convert(graph, gens, base, limits, memo, depth + 1)?.pow(exponent, limits.terms)?
            },
            | None => {
                // Not polynomial structure: a generator, with its arguments
                // expanded.
                let mut expanded = Vec::with_capacity(children.len());
                for child in &children {
                    let mut inner_gens = Gens::default();
                    let inner = from_term(graph, &mut inner_gens, *child, limits)?;
                    expanded.push(to_term(graph, &inner_gens, &inner));
                }
                let atom = if expanded == children { node } else { graph.try_node(graph.op(node), &expanded)? };
                Poly::generator(gens.index(graph, atom))
            },
        }
    };
    memo.insert(node, poly.clone());
    Some(poly)
}

/// Builds the term for one monomial times a coefficient.
fn term_of(
    graph: &mut Graph,
    gens: &Gens,
    mono: &Mono,
    coeff: &Number,
) -> NodeId {
    let mut factors = Vec::with_capacity(mono.len().saturating_add(1));
    if !coeff.is_one() || mono.is_empty() {
        factors.push(graph.num(coeff.clone()));
    }
    for &(generator, exponent) in mono {
        let base = gens.node(generator).unwrap_or(NodeId::NONE);
        if exponent == 1 {
            factors.push(base);
        } else {
            let e = graph.int(i64::from(exponent));
            factors.push(graph.node(core::POW, &[base, e]));
        }
    }
    match factors.as_slice() {
        | [only] => *only,
        | _ => graph.node(core::MUL, &factors),
    }
}

/// Converts a polynomial back to a term in expanded form.
pub fn to_term(
    graph: &mut Graph,
    gens: &Gens,
    poly: &Poly,
) -> NodeId {
    let terms: Vec<NodeId> = poly.terms().map(|(m, c)| (m.clone(), c.clone())).collect::<Vec<_>>()
        .iter()
        .map(|(mono, coeff)| term_of(graph, gens, mono, coeff))
        .collect();
    match terms.as_slice() {
        | [] => graph.int(0),
        | [only] => *only,
        | _ => graph.node(core::ADD, &terms),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn expand(src: &str) -> String {
        let mut g = Graph::new();
        let node = g.parse(src).unwrap_or_else(|e| panic!("{e}"));
        let mut gens = Gens::default();
        let poly = from_term(&mut g, &mut gens, node, Limits::default()).unwrap_or_else(|| panic!("limit"));
        let out = to_term(&mut g, &gens, &poly);
        g.display(out)
    }

    #[test]
    fn expansion() {
        assert_eq!(expand("(x + 1)^2"), "x^2 + 2*x + 1");
        assert_eq!(expand("(x + y) * (x - y)"), "x^2 - y^2");
        assert_eq!(expand("(a + b)^3"), "a^3 + 3*a^2*b + 3*a*b^2 + b^3");
        assert_eq!(expand("(x + 1) * (x + 2) * (x + 3)"), "x^3 + 6*x^2 + 11*x + 6");
        assert_eq!(expand("(x - 1)^2 - (x^2 - 2*x + 1)"), "0");
        assert_eq!(expand("(x/2 + 1/3)^2"), "1/4*x^2 + 1/3*x + 1/9");
    }

    #[test]
    fn non_polynomial_parts_are_generators() {
        assert_eq!(expand("(f(x) + 1)^2"), "f(x)^2 + 2*f(x) + 1");
        assert_eq!(expand("(x^(1/2) + 1)^2"), "(x^(1/2))^2 + 2*x^(1/2) + 1");
        assert_eq!(expand("(x + 1)^n * (x + 1)"), "x*(x + 1)^n + (x + 1)^n");
        assert_eq!(expand("f((x + 1)^2)"), "f(x^2 + 2*x + 1)", "expansion is deep");
    }

    #[test]
    fn limits_are_enforced() {
        let mut g = Graph::new();
        let node = g.parse("(a + b + c + d + e)^40").unwrap_or_else(|e| panic!("{e}"));
        let mut gens = Gens::default();
        assert!(from_term(&mut g, &mut gens, node, Limits { terms: 1_000, exponent: 256 }).is_none());
        let node = g.parse("(a + b)^1000").unwrap_or_else(|e| panic!("{e}"));
        let poly = from_term(&mut g, &mut gens, node, Limits::default());
        assert!(poly.is_some_and(|p| p.len() == 1), "beyond the exponent limit the power is a generator");
    }

    #[test]
    fn univariate_views() {
        let mut g = Graph::new();
        let node = g.parse("3*x^2*y + x*y^2 + 5*x + 7").unwrap_or_else(|e| panic!("{e}"));
        let mut gens = Gens::default();
        let poly = from_term(&mut g, &mut gens, node, Limits::default()).unwrap_or_default();
        let x_node = g.sym("x");
        let x = gens.find(&g, x_node).unwrap_or(u32::MAX);
        assert_eq!(poly.degree_in(x), 2);
        assert_eq!(poly.total_degree(), 3);
        let coeffs = poly.coefficients_in(x);
        let shown: Vec<String> = coeffs
            .iter()
            .map(|c| {
                let t = to_term(&mut g, &gens, c);
                g.display(t)
            })
            .collect();
        assert_eq!(shown, vec!["7", "y^2 + 5", "3*y"]);
        assert!(poly.univariate_in(x).is_none());
        let d = poly.derivative(x);
        let t = to_term(&mut g, &gens, &d);
        assert_eq!(g.display(t), "6*x*y + y^2 + 5");
    }

    #[test]
    fn arithmetic_laws() {
        let x = Poly::generator(0);
        let y = Poly::generator(1);
        let one = Poly::constant(Number::from(1));
        let a = x.add(&one);
        let b = y.sub(&x);
        let cap = 1_000;
        let left = a.mul(&b, cap).and_then(|p| p.mul(&a, cap));
        let right = a.pow(2, cap).and_then(|p| p.mul(&b, cap));
        assert_eq!(left, right);
        assert!(a.sub(&a).is_zero());
        assert_eq!(a.pow(0, cap), Some(one));
        assert_eq!(a.scale(&Number::from(0)), Poly::zero());
        assert_eq!(Poly::from_univariate(0, &[Number::from(1), Number::from(1)]), a);
    }
}
