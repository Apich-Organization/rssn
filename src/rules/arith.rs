//! Arithmetic: the numeric tower and the normal form of sums, products and
//! powers.
//!
//! Almost all of the work happens in [`Collect`], a window pass. It rewrites
//! a tree window in place into the usual canonical shape — like terms of a
//! sum merged, equal bases of a product merged, numbers folded — without
//! leaving intermediate terms in the graph. The e-graph side only carries
//! what must hold modulo known equalities: exact folding of all-numeric
//! nodes and the identities of 0 and 1.

use crate::graph::Cx;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::OpFlags;
use crate::graph::Number;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::RuleError;
use crate::graph::RuleSet;
use crate::graph::Tier;
use crate::graph::TreeWindow;
use crate::graph::WindowPass;
use crate::graph::op::core;
use crate::graph::rule::Installer;
use crate::graph::window::CellId;

/// The arithmetic rule set. Every other set depends on it.
#[must_use]
pub fn arith() -> RuleSet {
    RuleSet::new("arith", install)
}

fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    i.kernel("arith/fold", Tier::Normalize, Fold);
    i.kernel("arith/float", Tier::Normalize, FloatContagion);
    i.kernel("arith/radical", Tier::Normalize, Radical);
    i.pass(Collect);
    i.rewrites(
        Tier::Normalize,
        &[
            "arith/add-0: ?a + 0 => ?a",
            "arith/mul-1: ?a * 1 => ?a",
            "arith/mul-0: ?a * 0 => 0",
            "arith/pow-1: ?a ^ 1 => ?a",
            "arith/pow-0: ?a ^ 0 => 1",
            "arith/one-pow: 1 ^ ?a => 1",
            "arith/like-terms: ?c * ?a + ?d * ?a => (?c + ?d) * ?a if number(?c), number(?d)",
            "arith/like-terms-1: ?a + ?c * ?a => (1 + ?c) * ?a if number(?c)",
            "arith/double: ?a + ?a => 2 * ?a",
            "arith/square: ?a * ?a => ?a ^ 2",
            // The window pass merges powers it can see; these do the same
            // through equalities (`sec(x)` known to be `cos(x)^(-1)`).
            "arith/pow-merge-1: ?a * ?a ^ ?n => ?a ^ (?n + 1)",
            "arith/pow-merge: ?a ^ ?m * ?a ^ ?n => ?a ^ (?m + ?n)",
            "arith/pow-pow: (?a ^ ?m) ^ ?n => ?a ^ (?m * ?n) if integer(?n)",
            "arith/pow-pow-positive: (?a ^ ?m) ^ ?n => ?a ^ (?m * ?n) if positive(?a), real(?m)",
            "arith/pow-positive-factor: (?c * ?a) ^ ?e => ?c ^ ?e * ?a ^ ?e if number(?c), positive(?c)",
            "arith/pow-integer-factor: (?c * ?a) ^ ?n => ?c ^ ?n * ?a ^ ?n if number(?c), integer(?n)",
        ],
    )
}

/// Folds sums, products and powers whose operands are all literal numbers,
/// exactly where possible.
struct Fold;

impl Kernel for Fold {
    fn ops(&self) -> Vec<OpId> {
        vec![core::ADD, core::MUL, core::POW]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let values: Option<Vec<Number>> = cx
            .graph
            .children(node)
            .iter()
            .map(|&c| cx.graph.number_of(c).cloned())
            .collect();
        let Some(values) = values else {
            return Outcome::Pass;
        };
        let folded = match (cx.graph.op(node), values.as_slice()) {
            | (core::ADD, vs) => Some(vs.iter().fold(Number::from(0), |a, b| a.add(b))),
            | (core::MUL, vs) => Some(vs.iter().fold(Number::from(1), |a, b| a.mul(b))),
            | (core::POW, [b, e]) => b.pow(e),
            | _ => None,
        };
        folded.map_or(Outcome::Pass, |n| Outcome::Equal(cx.graph.num(n)))
    }

    fn revisit(&self) -> bool {
        true
    }
}

/// Evaluates any operator with scalar semantics once all its operands are
/// numbers and at least one of them is a float.
///
/// `sin(1)` is an exact value and stays symbolic; `sin(1.0)` has already
/// left the exact phase, so it is evaluated.
struct FloatContagion;

impl Kernel for FloatContagion {
    fn ops(&self) -> Vec<OpId> {
        vec![OpId::NONE]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let graph = &mut *cx.graph;
        let desc = graph.ops().get(graph.op(node));
        if desc.flags.has(OpFlags::PREDICATE) {
            return Outcome::Pass;
        }
        let Some(eval) = desc.eval else {
            return Outcome::Pass;
        };
        let children = graph.children(node);
        if children.is_empty() {
            return Outcome::Pass;
        }
        let mut any_float = false;
        let mut args = Vec::with_capacity(children.len());
        for &child in children {
            match graph.number_of(child) {
                | Some(n) => {
                    any_float |= !n.is_exact();
                    args.push(n.to_f64());
                },
                | None => return Outcome::Pass,
            }
        }
        if !any_float {
            return Outcome::Pass;
        }
        let value = eval(&args);
        if value.is_finite() {
            Outcome::Equal(graph.float(value))
        } else {
            Outcome::Pass
        }
    }

    fn revisit(&self) -> bool {
        true
    }
}

/// Window pass producing the canonical form of sums, products and powers.
///
/// * sum: numeric terms added up, terms with the same non-numeric part
///   merged (`2*x + 3*x` → `5*x`), zero terms dropped;
/// * product: a numeric coefficient distributed over a lone sum
///   (`2*(x + y)` → `2*x + 2*y`), numeric factors multiplied, equal bases merged by adding
///   numeric exponents (`x * x^2` → `x^3`), a zero factor annihilates;
/// * power: numeric powers folded exactly, `b^1`, `b^0`, `1^e`, nested
///   integer powers multiplied out, integer powers distributed over
///   products (`(2*x)^2` → `4*x^2`).
///
/// Like every computer algebra system this treats `x^0` as `1` and `x/x` as
/// `1`, i.e. it simplifies generically and ignores the measure-zero set
/// where the original expression was undefined.
pub struct Collect;

impl WindowPass for Collect {
    fn run(
        &self,
        graph: &mut Graph,
        window: &mut TreeWindow,
    ) -> bool {
        let mut changed = false;
        for cell in window.post_order() {
            if !window.is_live(cell) || window.atom(cell).is_some() {
                continue;
            }
            let rebuilt = match window.op(cell) {
                | core::ADD => Some(sum(graph, window, cell)),
                | core::MUL => Some(product(graph, window, cell)),
                | core::POW => power(graph, window, cell),
                | _ => None,
            };
            let Some(new) = rebuilt else {
                continue;
            };
            if window.fingerprint(graph, new) == window.fingerprint(graph, cell) {
                window.discard(new);
            } else {
                window.replace(cell, new);
                changed = true;
            }
        }
        changed
    }
}

fn number_atom(
    graph: &mut Graph,
    window: &mut TreeWindow,
    n: Number,
) -> CellId {
    let node = graph.num(n);
    window.new_atom(graph, node)
}

/// Finds the entry of `groups` whose key cells are structurally the same as
/// `key`, comparing element-wise after sorting by fingerprint.
fn find_group<T>(
    graph: &Graph,
    window: &TreeWindow,
    groups: &[(Vec<CellId>, T)],
    key: &[CellId],
) -> Option<usize> {
    groups.iter().position(|(other, _)| {
        other.len() == key.len()
            && other
                .iter()
                .zip(key)
                .all(|(&a, &b)| window.same(graph, a, b))
    })
}

fn sorted_by_fingerprint(
    graph: &Graph,
    window: &TreeWindow,
    mut cells: Vec<CellId>,
) -> Vec<CellId> {
    cells.sort_by_key(|&c| window.fingerprint(graph, c));
    cells
}

/// Builds a detached product cell of `coeff` and clones of `factors`.
fn build_product(
    graph: &mut Graph,
    window: &mut TreeWindow,
    coeff: &Number,
    factors: &[CellId],
) -> CellId {
    let mut pieces = Vec::with_capacity(factors.len().saturating_add(1));
    if !coeff.is_one() || factors.is_empty() {
        pieces.push(number_atom(graph, window, coeff.clone()));
    }
    for &f in factors {
        pieces.push(window.clone_subtree(f));
    }
    match pieces.as_slice() {
        | [only] => *only,
        | _ => window.new_node(core::MUL, &pieces),
    }
}

/// Splits a term of a sum into its numeric coefficient and the remaining
/// factors.
fn split_term(
    graph: &Graph,
    window: &TreeWindow,
    term: CellId,
) -> (Number, Vec<CellId>) {
    if let Some(n) = window.number(graph, term) {
        return (n.clone(), Vec::new());
    }
    if window.atom(term).is_some() || window.op(term) != core::MUL {
        return (Number::from(1), vec![term]);
    }
    let mut coeff = Number::from(1);
    let mut factors = Vec::new();
    for child in window.children(term) {
        match window.number(graph, child) {
            | Some(n) => coeff = coeff.mul(n),
            | None => factors.push(child),
        }
    }
    (coeff, factors)
}

fn sum(
    graph: &mut Graph,
    window: &mut TreeWindow,
    cell: CellId,
) -> CellId {
    let mut constant = Number::from(0);
    let mut groups: Vec<(Vec<CellId>, Number)> = Vec::new();
    for term in window.children(cell).collect::<Vec<_>>() {
        let (coeff, factors) = split_term(graph, window, term);
        if factors.is_empty() {
            constant = constant.add(&coeff);
            continue;
        }
        let key = sorted_by_fingerprint(graph, window, factors);
        match find_group(graph, window, &groups, &key) {
            | Some(i) => {
                if let Some(group) = groups.get_mut(i) {
                    group.1 = group.1.add(&coeff);
                }
            },
            | None => groups.push((key, coeff)),
        }
    }
    let mut pieces = Vec::with_capacity(groups.len().saturating_add(1));
    for (factors, coeff) in &groups {
        if !coeff.is_zero() {
            pieces.push(build_product(graph, window, coeff, factors));
        }
    }
    if !constant.is_zero() || pieces.is_empty() {
        pieces.push(number_atom(graph, window, constant));
    }
    match pieces.as_slice() {
        | [only] => *only,
        | _ => window.new_node(core::ADD, &pieces),
    }
}

/// Splits a factor of a product into base and numeric exponent.
fn split_factor(
    graph: &Graph,
    window: &TreeWindow,
    factor: CellId,
) -> (CellId, Number) {
    if window.atom(factor).is_none() && window.op(factor) == core::POW {
        let mut children = window.children(factor);
        if let (Some(base), Some(exp)) = (children.next(), children.next()) {
            if let Some(n) = window.number(graph, exp) {
                return (base, n.clone());
            }
        }
    }
    (factor, Number::from(1))
}

fn product(
    graph: &mut Graph,
    window: &mut TreeWindow,
    cell: CellId,
) -> CellId {
    let mut coeff = Number::from(1);
    let mut groups: Vec<(Vec<CellId>, Number)> = Vec::new();
    for factor in window.children(cell).collect::<Vec<_>>() {
        if let Some(n) = window.number(graph, factor) {
            coeff = coeff.mul(n);
            continue;
        }
        let (base, exp) = split_factor(graph, window, factor);
        if let Some(folded) = window.number(graph, base).and_then(|b| b.pow(&exp)) {
            coeff = coeff.mul(&folded);
            continue;
        }
        let key = [base];
        match find_group(graph, window, &groups, &key) {
            | Some(i) => {
                if let Some(group) = groups.get_mut(i) {
                    group.1 = group.1.add(&exp);
                }
            },
            | None => groups.push((key.to_vec(), exp)),
        }
    }
    if coeff.is_zero() {
        return number_atom(graph, window, coeff);
    }
    // A numeric coefficient is distributed over a sum so that like terms on
    // both sides of a parenthesis can meet: 2*(x + y) - 2*x.
    if let [(base, exp)] = groups.as_slice() {
        if let Some(&base) = base.first() {
            if exp.is_one()
                && !coeff.is_one()
                && window.atom(base).is_none()
                && window.op(base) == core::ADD
            {
                let terms: Vec<CellId> = window.children(base).collect();
                let pieces: Vec<CellId> = terms
                    .iter()
                    .map(|&t| build_product(graph, window, &coeff, &[t]))
                    .collect();
                return window.new_node(core::ADD, &pieces);
            }
        }
    }
    let mut pieces = Vec::with_capacity(groups.len().saturating_add(1));
    for (base, exp) in &groups {
        let Some(&base) = base.first() else {
            continue;
        };
        if exp.is_zero() {
            continue;
        }
        let base_copy = window.clone_subtree(base);
        if exp.is_one() {
            pieces.push(base_copy);
        } else {
            let exp_atom = number_atom(graph, window, exp.clone());
            pieces.push(window.new_node(core::POW, &[base_copy, exp_atom]));
        }
    }
    if !coeff.is_one() || pieces.is_empty() {
        pieces.insert(0, number_atom(graph, window, coeff));
    }
    match pieces.as_slice() {
        | [only] => *only,
        | _ => window.new_node(core::MUL, &pieces),
    }
}

fn power(
    graph: &mut Graph,
    window: &mut TreeWindow,
    cell: CellId,
) -> Option<CellId> {
    let mut children = window.children(cell);
    let (base, exp) = (children.next()?, children.next()?);
    drop(children);
    let base_num = window.number(graph, base).cloned();
    let exp_num = window.number(graph, exp).cloned();
    if let (Some(b), Some(e)) = (&base_num, &exp_num) {
        if let Some(folded) = b.pow(e) {
            return Some(number_atom(graph, window, folded));
        }
    }
    if let Some(e) = &exp_num {
        if e.is_one() {
            return Some(window.clone_subtree(base));
        }
        if e.is_zero() {
            return Some(number_atom(graph, window, Number::from(1)));
        }
    }
    if let Some(b) = &base_num {
        if b.is_one() {
            return Some(number_atom(graph, window, Number::from(1)));
        }
        if b.is_zero()
            && exp_num
                .as_ref()
                .is_some_and(|e| !e.is_negative() && !e.is_zero())
        {
            return Some(number_atom(graph, window, Number::from(0)));
        }
    }
    // The remaining rewrites need an integer exponent and a compound base.
    let n = exp_num.filter(Number::is_integer)?;
    if window.atom(base).is_some() {
        return None;
    }
    match window.op(base) {
        | core::POW => {
            // (b^e)^n = b^(e*n) for integer n.
            let mut inner = window.children(base);
            let (inner_base, inner_exp) = (inner.next()?, inner.next()?);
            drop(inner);
            let base_copy = window.clone_subtree(inner_base);
            let new_exp = match window.number(graph, inner_exp).cloned() {
                | Some(e) => number_atom(graph, window, e.mul(&n)),
                | None => {
                    let exp_copy = window.clone_subtree(inner_exp);
                    let n_atom = number_atom(graph, window, n);
                    window.new_node(core::MUL, &[n_atom, exp_copy])
                },
            };
            Some(window.new_node(core::POW, &[base_copy, new_exp]))
        },
        | core::MUL => {
            // (a*b)^n = a^n * b^n for integer n.
            let factors: Vec<CellId> = window.children(base).collect();
            let mut pieces = Vec::with_capacity(factors.len());
            for factor in factors {
                let copy = window.clone_subtree(factor);
                let n_atom = number_atom(graph, window, n.clone());
                pieces.push(window.new_node(core::POW, &[copy, n_atom]));
            }
            Some(window.new_node(core::MUL, &pieces))
        },
        | _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Budget;
    use crate::graph::ClosedForm;
    use crate::graph::Engine;
    use crate::graph::Env;
    use crate::graph::Extractor;
    use crate::graph::Saturate;

    fn simplify(src: &str) -> String {
        let mut g = Graph::new();
        let engine = Engine::install(&mut g, &[arith()]).unwrap_or_else(|e| panic!("{e}"));
        let root = g.parse(src).unwrap_or_else(|e| panic!("{e}"));
        engine.run(
            &mut g,
            &[root],
            &Env::symbolic(),
            &Saturate,
            &Budget::default(),
        );
        assert!(g.conflicts().is_empty(), "{:?}", g.conflicts());
        assert_eq!(g.validate(), Ok(()));
        let ex = Extractor::new(&g, &[root], &ClosedForm);
        ex.build(&mut g, root)
            .map_or_else(|| "<none>".to_owned(), |n| g.display(n))
    }

    #[test]
    fn numbers_fold_exactly() {
        assert_eq!(simplify("1 + 2 * 3"), "7");
        assert_eq!(simplify("1/3 + 1/6"), "1/2");
        assert_eq!(simplify("2^10"), "1024");
        assert_eq!(simplify("2^(-2)"), "1/4");
        assert_eq!(simplify("4^(1/2)"), "2");
        assert_eq!(simplify("8^(2/3)"), "4");
        assert_eq!(
            simplify("2^(1/2)"),
            "2^(1/2)",
            "algebraic numbers are not approximated"
        );
    }

    #[test]
    fn floats_are_contagious() {
        assert_eq!(simplify("0.5 + 1/2"), "1");
        assert_eq!(simplify("4.0^(1/2)"), "2");
        assert_eq!(simplify("x + 0.5 + 0.25"), "x + 0.75");
    }

    #[test]
    fn like_terms_are_collected() {
        assert_eq!(simplify("x + x"), "2*x");
        assert_eq!(simplify("2*x + 3*x"), "5*x");
        assert_eq!(simplify("x - x"), "0");
        assert_eq!(simplify("x*y + y*x"), "2*x*y");
        assert_eq!(simplify("3*a*b - a*b*3 + c"), "c");
        assert_eq!(simplify("x + 1 + x + 2"), "2*x + 3");
        assert_eq!(simplify("x/2 + x/2"), "x");
    }

    #[test]
    fn powers_are_merged() {
        assert_eq!(simplify("x * x"), "x^2");
        assert_eq!(simplify("x * x^2 * x^3"), "x^6");
        assert_eq!(simplify("x / x"), "1");
        assert_eq!(simplify("x^2 / x"), "x");
        assert_eq!(simplify("x^a * x^b"), "x^(a + b)");
        assert_eq!(simplify("2^(k + 1) / 2^k"), "2");
        assert_eq!(simplify("(x^2)^3"), "x^6");
        assert_eq!(simplify("(x^a)^2"), "x^(2*a)");
        assert_eq!(simplify("(2*x)^2"), "4*x^2");
        assert_eq!(simplify("(x*y)^2 / x"), "x*y^2");
    }

    #[test]
    fn identities_of_zero_and_one() {
        assert_eq!(simplify("0 * x + 1 * y + 0"), "y");
        assert_eq!(simplify("x^1"), "x");
        assert_eq!(simplify("x^0"), "1");
        assert_eq!(simplify("1^x"), "1");
        assert_eq!(simplify("0^2"), "0");
        assert_eq!(simplify("(x + 0)^(1 + 0)"), "x");
    }

    #[test]
    fn nested_structures_normalise_bottom_up() {
        assert_eq!(simplify("(a + a) * (b + b)"), "4*a*b");
        assert_eq!(simplify("2 * (x + x) + x"), "5*x");
        assert_eq!(simplify("(x + y) - (y + x)"), "0");
        assert_eq!(simplify("(x + y) * 2 - 2 * (y + x)"), "0");
        assert_eq!(simplify("3 * (x + 2*y) - x"), "2*x + 6*y");
        assert_eq!(simplify("-(a - b)"), "b - a");
    }

    #[test]
    fn shared_subterms_are_simplified_once_and_reused() {
        // `x + x` occurs under two different parents, so it is an atom of
        // the outer window and gets a window of its own.
        assert_eq!(
            simplify("f(x + x) + g(x + x)"),
            "f(2*x) + g(2*x)"
        );
    }
}

/// Exact powers of integers and fractions with rational exponents in
/// their normal form `c * m^(r/q)`: the integer part of the exponent split
/// off (`2^(3/2) = 2 * 2^(1/2)`), perfect powers taken out of the base
/// (`12^(1/2) = 2 * 3^(1/2)`) and fractions split (`(2/3)^(1/2) =
/// 2^(1/2) * 3^(-1/2)`).
struct Radical;

/// Largest `c` with `c^q | n` and the cofactor `n / c^q`, by trial
/// division (bases up to 10^12 only).
fn extract_power(
    n: u64,
    q: u32,
) -> Option<(u64, u64)> {
    if n > 1_000_000_000_000 || q < 2 {
        return None;
    }
    let (mut remaining, mut root, mut kept) = (n, 1_u64, 1_u64);
    let mut p = 2_u64;
    while p.checked_mul(p)? <= remaining {
        let mut count = 0_u32;
        while remaining % p == 0 {
            remaining /= p;
            count += 1;
        }
        root = root.checked_mul(p.checked_pow(count / q)?)?;
        kept = kept.checked_mul(p.checked_pow(count % q)?)?;
        p += 1;
    }
    Some((root, kept.checked_mul(remaining)?))
}

/// `n^e` for a positive integer `n` and rational non-integer `e` as
/// `(c, m, r/q)` with `n^e = c * m^(r/q)`, `0 < r < q`, `m` free of
/// `q`-th powers (`m = 1` when the power is exact).
fn radical_form(
    n: &num_bigint::BigInt,
    e: &num_rational::BigRational,
) -> Option<(num_rational::BigRational, u64, num_rational::BigRational)> {
    use num_traits::ToPrimitive;
    let n = n.to_u64()?;
    let q = e.denom().to_u32()?;
    let p = e.numer().clone();
    let k = num_integer::Integer::div_floor(&p, &num_bigint::BigInt::from(q));
    let r = &p - &k * num_bigint::BigInt::from(q);
    let (root, rest) = extract_power(n, q)?;
    let r_u32 = r.to_u32()?;
    let coefficient = num_rational::BigRational::from_integer(num_bigint::BigInt::from(n)).pow(k.to_i32()?)
        * num_rational::BigRational::from_integer(num_bigint::BigInt::from(root).pow(r_u32));
    let rest = if r_u32 == 0 { 1 } else { rest };
    Some((coefficient, rest, num_rational::BigRational::new(r, num_bigint::BigInt::from(q))))
}

impl Radical {
    /// `c * n1^e1 * n2^e2 * ... * rest` with every radical in normal form
    /// and all rational factors in `c`.
    fn product(
        graph: &mut Graph,
        node: NodeId,
    ) -> Outcome {
        use num_traits::One;
        use num_traits::Signed;
        let children = graph.children(node).to_vec();
        let mut coefficient = Number::from(1);
        let mut others = Vec::new();
        let mut radicals = Vec::new();
        let mut changed = false;
        for &c in &children {
            if let Some(n) = graph.as_number(c) {
                coefficient = coefficient.mul(n);
                continue;
            }
            let normal = (graph.op(c) == core::POW)
                .then(|| {
                    let &[b, e] = graph.children(c) else {
                        return None;
                    };
                    let b = graph.as_number(b)?.to_rational()?;
                    let e = graph.as_number(e)?.to_rational()?;
                    if !b.is_integer() || !b.is_positive() || e.is_integer() {
                        return None;
                    }
                    radical_form(&b.to_integer(), &e)
                })
                .flatten();
            match normal {
                | Some((c_part, rest, frac)) => {
                    if !c_part.is_one() {
                        changed = true;
                    }
                    coefficient = coefficient.mul(&Number::rat(c_part));
                    if rest != 1 {
                        radicals.push((rest, frac));
                    } else {
                        changed = true;
                    }
                },
                | None => others.push(c),
            }
        }
        if !changed || matches!(coefficient, Number::Float(_)) {
            return Outcome::Pass;
        }
        let mut factors = Vec::new();
        if !coefficient.is_one() {
            factors.push(graph.num(coefficient));
        }
        for (rest, frac) in radicals {
            let base = graph.int(i64::try_from(rest).unwrap_or(i64::MAX));
            let e = graph.num(Number::rat(frac));
            factors.push(graph.node(core::POW, &[base, e]));
        }
        factors.extend(others);
        Outcome::Equal(match factors.as_slice() {
            | [] => graph.int(1),
            | [only] => *only,
            | _ => graph.node(core::MUL, &factors),
        })
    }
}

impl Kernel for Radical {
    fn ops(&self) -> Vec<OpId> {
        vec![core::POW, core::MUL]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        use num_traits::Signed;
        use num_traits::ToPrimitive;
        use num_traits::Zero;
        let graph = &mut *cx.graph;
        if graph.op(node) == core::MUL {
            return Self::product(graph, node);
        }
        let &[base, exponent] = graph.children(node) else {
            return Outcome::Pass;
        };
        let (Some(b), Some(e)) = (graph.number_of(base).cloned(), graph.number_of(exponent).cloned()) else {
            return Outcome::Pass;
        };
        let (Some(b), Some(e)) = (b.to_rational(), e.to_rational()) else {
            return Outcome::Pass;
        };
        if e.is_integer() || !b.is_positive() {
            return Outcome::Pass;
        }
        // A fraction: numerator and denominator separately.
        if !b.is_integer() {
            let numer = graph.num(Number::rat(num_rational::BigRational::from_integer(b.numer().clone())));
            let denom = graph.num(Number::rat(num_rational::BigRational::from_integer(b.denom().clone())));
            let minus = graph.num(Number::rat(-e.clone()));
            let top = if b.numer() == &num_bigint::BigInt::from(1) { None } else { Some(graph.node(core::POW, &[numer, exponent])) };
            let bottom = graph.node(core::POW, &[denom, minus]);
            return Outcome::Equal(match top {
                | Some(t) => graph.node(core::MUL, &[t, bottom]),
                | None => bottom,
            });
        }
        let (Some(n), Some(q)) = (b.to_integer().to_u64(), e.denom().to_u32()) else {
            return Outcome::Pass;
        };
        let p = e.numer().clone();
        // e = k + r/q with 0 <= r < q
        let k = num_integer::Integer::div_floor(&p, &num_bigint::BigInt::from(q));
        let r = &p - &k * num_bigint::BigInt::from(q);
        let Some((root, rest)) = extract_power(n, q) else {
            return Outcome::Pass;
        };
        if k.is_zero() && root == 1 {
            return Outcome::Pass;
        }
        // n^(p/q) = n^k * root^r * rest^(r/q)
        let Some(r_u32) = r.to_u32() else {
            return Outcome::Pass;
        };
        let coefficient = num_rational::BigRational::from_integer(num_bigint::BigInt::from(n)).pow(k.to_i32().unwrap_or(0))
            * num_rational::BigRational::from_integer(num_bigint::BigInt::from(root).pow(r_u32));
        let c = graph.num(Number::rat(coefficient));
        if rest == 1 || r.is_zero() {
            return Outcome::Equal(c);
        }
        let rest_node = graph.int(i64::try_from(rest).unwrap_or(i64::MAX));
        let fraction = graph.num(Number::rat(num_rational::BigRational::new(r, num_bigint::BigInt::from(q))));
        let radical = graph.node(core::POW, &[rest_node, fraction]);
        Outcome::Equal(graph.node(core::MUL, &[c, radical]))
    }
}

#[cfg(test)]
mod radical_tests {
    use crate::rules::testing::simplify;

    #[test]
    fn radicals_in_normal_form() {
        let run = |s: &str| simplify(&[super::arith()], s);
        assert_eq!(run("-1/2*2^(3/2)"), "-2^(1/2)");
        assert_eq!(run("3*12^(1/2)"), "6*3^(1/2)");
        assert_eq!(run("8^(2/3)"), "4");
        assert_eq!(run("2*2^(-1/2)"), "2^(1/2)");
        assert_eq!(run("x*18^(1/2)/3"), "2^(1/2)*x");
    }
}
