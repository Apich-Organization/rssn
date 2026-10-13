//! Certified numeric roots by interval Newton / Krawczyk iteration.
//!
//! The equations and their symbolic Jacobian are lowered to tapes and
//! evaluated with the outward-rounded interval arithmetic of
//! [`crate::backend::tape`]. For a square system `F(x) = 0` and a box `X`
//! with midpoint `y` and `C ≈ J(y)^-1`, the Krawczyk operator
//!
//! `K(X) = y - C F(y) + (I - C J(X)) (X - y)`
//!
//! satisfies: a root in `X` lies in `K(X)`; if `K(X)` lies in the interior
//! of `X` then `X` contains exactly one root (and `K(X) ∩ X` is a tighter
//! enclosure of it).
//!
//! * `verify_near` turns an approximate root into a verified enclosure.
//! * `search` is an exhaustive bisection with exclusion (a component of
//!   `F(X)` that excludes 0, or `K(X) ∩ X` empty) and Krawczyk
//!   certification: it returns the boxes that provably contain exactly one
//!   root, the boxes it could not decide (roots of even multiplicity,
//!   singular Jacobians, regions where `F` is undefined, an exhausted
//!   budget) and whether the list of undecided boxes is empty, which is
//!   the proof that every root in the box was found.
//! * `nsolve_certified(eq, x, guess)` and `solve_box(eqs, vars, box)` expose
//!   them as operators.

use super::as_expression;
use crate::backend::tape::lower;
use crate::backend::tape::Interval;
use crate::backend::tape::Scalar;
use crate::backend::tape::Tape;
use crate::graph::op::core;
use crate::graph::Ball;
use crate::graph::Cx;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::SymbolId;
use crate::rules::calculus::derivative;
use crate::rules::poly::best;

/// A square or rectangular system with its Jacobian, compiled for interval
/// evaluation.
pub(super) struct System {
    tape: Tape,
    unknowns: usize,
    equations: usize,
    parameters: Vec<Interval>,
}

fn midpoint(i: Interval) -> f64 {
    if i.lo.is_infinite() || i.hi.is_infinite() {
        return 0.0;
    }
    0.5 * i.lo + 0.5 * i.hi
}

fn width(i: Interval) -> f64 {
    i.hi - i.lo
}

fn intersect(
    a: Interval,
    b: Interval,
) -> Option<Interval> {
    let (lo, hi) = (a.lo.max(b.lo), a.hi.min(b.hi));
    (lo <= hi).then_some(Interval { lo, hi })
}

fn interior(
    inner: Interval,
    outer: Interval,
) -> bool {
    outer.lo < inner.lo && inner.hi < outer.hi
}

impl System {
    /// Compiles `equations` (each `l = r` or an expression equal to zero)
    /// in `unknowns`; `bindings` give the parameters.
    pub(super) fn new(
        graph: &mut Graph,
        equations: &[NodeId],
        unknowns: &[SymbolId],
        bindings: &[(SymbolId, f64)],
    ) -> Option<Self> {
        let nodes: Vec<NodeId> = unknowns.iter().map(|&s| graph.symbol_node(s)).collect();
        let mut terms = Vec::new();
        for &e in equations {
            let expr = as_expression(graph, e);
            terms.push(best(graph, expr)?);
        }
        let mut roots = terms.clone();
        for &t in &terms {
            for &x in &nodes {
                roots.push(derivative(graph, t, x)?);
            }
        }
        let mut params = Vec::new();
        let mut values = Vec::new();
        for &(s, v) in bindings {
            if !unknowns.contains(&s) && !params.contains(&s) {
                params.push(s);
                values.push(Interval::point(v));
            }
        }
        let tape = lower(graph, &roots, unknowns, &params).ok()?;
        Some(Self { tape, unknowns: unknowns.len(), equations: equations.len(), parameters: values })
    }

    /// Enclosures of `F(X)` and of the Jacobian (row-major) over `x`.
    fn enclose(
        &self,
        x: &[Interval],
    ) -> (Vec<Interval>, Vec<Interval>) {
        let all = self.tape.enclose(x, &self.parameters);
        let (f, j) = all.split_at(self.equations.min(all.len()));
        (f.to_vec(), j.to_vec())
    }
}

fn invert(a: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    let n = a.len();
    let mut m: Vec<Vec<f64>> = a
        .iter()
        .enumerate()
        .map(|(i, row)| {
            let mut r = row.clone();
            r.extend((0..n).map(|j| if i == j { 1.0 } else { 0.0 }));
            r
        })
        .collect();
    for col in 0..n {
        let pivot = (col..n).max_by(|&i, &j| {
            let (x, y) = (m.get(i).and_then(|r| r.get(col)), m.get(j).and_then(|r| r.get(col)));
            x.copied().unwrap_or(0.0).abs().total_cmp(&y.copied().unwrap_or(0.0).abs())
        })?;
        m.swap(col, pivot);
        let p = *m.get(col)?.get(col)?;
        if !p.is_finite() || p.abs() < 1e-300 {
            return None;
        }
        for v in m.get_mut(col)?.iter_mut() {
            *v /= p;
        }
        let pivot_row = m.get(col)?.clone();
        for (r, row) in m.iter_mut().enumerate() {
            if r == col {
                continue;
            }
            let factor = *row.get(col)?;
            for (v, pv) in row.iter_mut().zip(&pivot_row) {
                *v -= factor * pv;
            }
        }
    }
    let out: Vec<Vec<f64>> = m.into_iter().map(|row| row.into_iter().skip(n).collect()).collect();
    out.iter().flatten().all(|v| v.is_finite()).then_some(out)
}

/// The Krawczyk image of `x`, or `None` when the system is not square or
/// the midpoint Jacobian is singular.
fn krawczyk(
    system: &System,
    x: &[Interval],
) -> Option<Vec<Interval>> {
    let n = system.unknowns;
    if system.equations != n {
        return None;
    }
    let y: Vec<f64> = x.iter().map(|&i| midpoint(i)).collect();
    let point: Vec<Interval> = y.iter().map(|&v| Interval::point(v)).collect();
    let (fy, _) = system.enclose(&point);
    let (_, jx) = system.enclose(x);
    let at = |r: usize, c: usize| jx.get(r * n + c).copied();
    let mid_j: Vec<Vec<f64>> =
        (0..n).map(|r| (0..n).map(|c| at(r, c).map_or(f64::NAN, midpoint)).collect()).collect();
    let c = invert(&mid_j)?;
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let mut k = Interval::point(*y.get(i)?);
        for (j, fj) in fy.iter().enumerate() {
            k = k.sub(Interval::point(*c.get(i)?.get(j)?).mul(*fj));
        }
        for kk in 0..n {
            let mut m = Interval::point(if i == kk { 1.0 } else { 0.0 });
            for j in 0..n {
                m = m.sub(Interval::point(*c.get(i)?.get(j)?).mul(at(j, kk)?));
            }
            let d = x.get(kk)?.sub(Interval::point(*y.get(kk)?));
            k = k.add(m.mul(d));
        }
        if k.lo.is_nan() || k.hi.is_nan() {
            return None;
        }
        out.push(k);
    }
    Some(out)
}

/// Contracts a box that is known to contain exactly one root as far as the
/// Krawczyk operator allows.
fn contract(
    system: &System,
    mut x: Vec<Interval>,
    relative: f64,
) -> Vec<Interval> {
    for _ in 0..60 {
        let Some(k) = krawczyk(system, &x) else { break };
        let next: Option<Vec<Interval>> = k.iter().zip(&x).map(|(&a, &b)| intersect(a, b)).collect();
        let Some(next) = next else { break };
        let progress = next.iter().zip(&x).any(|(a, b)| width(*a) < 0.9 * width(*b));
        x = next;
        if !progress || x.iter().all(|i| width(*i) <= relative * (1.0 + i.lo.abs().max(i.hi.abs()))) {
            break;
        }
    }
    x
}

/// A verified enclosure of a root near `guess`: a box that provably
/// contains exactly one root of the square system.
pub(super) fn verify_near(
    system: &System,
    guess: &[f64],
) -> Option<Vec<Interval>> {
    for scale in [1e-13_f64, 1e-11, 1e-9, 1e-7, 1e-5] {
        let x: Vec<Interval> = guess
            .iter()
            .map(|&g| {
                let r = scale * (1.0 + g.abs());
                Interval { lo: g - r, hi: g + r }
            })
            .collect();
        if let Some(k) = krawczyk(system, &x)
            && k.iter().zip(&x).all(|(&a, &b)| interior(a, b))
        {
            return Some(contract(system, x, 1e-15));
        }
    }
    None
}

/// The result of an exhaustive search.
pub(super) struct Search {
    /// Boxes each containing exactly one root.
    pub(super) verified: Vec<Vec<Interval>>,
    /// Boxes in which the search could not decide.
    pub(super) unresolved: Vec<Vec<Interval>>,
}

fn overlap(
    a: &[Interval],
    b: &[Interval],
) -> bool {
    a.iter().zip(b).all(|(p, q)| intersect(*p, *q).is_some())
}

/// Bisection with exclusion and Krawczyk certification over `region`.
pub(super) fn search(
    system: &System,
    region: &[Interval],
    tolerance: f64,
    max_boxes: usize,
) -> Search {
    let mut stack = vec![region.to_vec()];
    let mut verified: Vec<Vec<Interval>> = Vec::new();
    let mut unresolved: Vec<Vec<Interval>> = Vec::new();
    let mut processed = 0_usize;
    while let Some(x) = stack.pop() {
        processed += 1;
        if processed > max_boxes {
            unresolved.push(x);
            continue;
        }
        let (f, _) = system.enclose(&x);
        if f.iter().any(|i| i.lo > 0.0 || i.hi < 0.0) {
            continue;
        }
        let mut current = x;
        if let Some(k) = krawczyk(system, &current) {
            let clipped: Option<Vec<Interval>> = k.iter().zip(&current).map(|(&a, &b)| intersect(a, b)).collect();
            let Some(clipped) = clipped else { continue };
            if k.iter().zip(&current).all(|(&a, &b)| interior(a, b)) {
                let refined = contract(system, current, 1e-15);
                if !verified.iter().any(|v| overlap(v, &refined)) {
                    verified.push(refined);
                }
                continue;
            }
            current = clipped;
        }
        let Some((axis, widest)) = current
            .iter()
            .map(|i| width(*i))
            .enumerate()
            .max_by(|a, b| a.1.total_cmp(&b.1))
        else {
            continue;
        };
        let scale = current.iter().map(|i| 1.0 + i.lo.abs().max(i.hi.abs())).fold(1.0_f64, f64::max);
        if widest <= tolerance * scale || !widest.is_finite() {
            // A root on the boundary of the bisection can be certified by
            // an inflated box around the midpoint.
            let centre: Vec<f64> = current.iter().map(|&i| midpoint(i)).collect();
            if let Some(found) = verify_near(system, &centre) {
                if !verified.iter().any(|v| overlap(v, &found)) {
                    verified.push(found);
                }
            } else {
                unresolved.push(current);
            }
            continue;
        }
        let Some(interval) = current.get(axis).copied() else { continue };
        // Off-centre, so that symmetric problems do not put roots on cuts.
        let mid = interval.lo + 0.4142 * (interval.hi - interval.lo);
        let (mut left, mut right) = (current.clone(), current);
        if let Some(slot) = left.get_mut(axis) {
            *slot = Interval { lo: interval.lo, hi: mid };
        }
        if let Some(slot) = right.get_mut(axis) {
            *slot = Interval { lo: mid, hi: interval.hi };
        }
        stack.push(right);
        stack.push(left);
    }
    // Neighbouring undecided boxes that touch are one region.
    let mut kept: Vec<Vec<Interval>> = Vec::new();
    for u in unresolved {
        if verified.iter().any(|v| overlap(v, &u)) && u.iter().all(|i| width(*i) <= tolerance * (1.0 + i.lo.abs().max(i.hi.abs())) * 4.0)
        {
            // A sliver around a verified root: the root is already counted.
            continue;
        }
        kept.push(u);
    }
    Search { verified, unresolved: kept }
}

fn pair(
    graph: &mut Graph,
    i: Interval,
) -> NodeId {
    let (lo, hi) = (graph.float(i.lo), graph.float(i.hi));
    graph.node(core::LIST, &[lo, hi])
}

fn box_node(
    graph: &mut Graph,
    x: &[Interval],
) -> NodeId {
    let items: Vec<NodeId> = x.iter().map(|&i| pair(graph, i)).collect();
    graph.node(core::LIST, &items)
}

/// Reads `list(lo, hi)` pairs (or one pair) as a box.
fn read_box(
    graph: &mut Graph,
    cx_env: &crate::graph::Env,
    node: NodeId,
    dimension: usize,
) -> Option<Vec<Interval>> {
    let value = |graph: &mut Graph, n: NodeId| best(graph, n).and_then(|t| graph.eval(t, cx_env));
    let children = graph.children(node).to_vec();
    let pairs: Vec<NodeId> =
        if dimension == 1 && children.first().is_some_and(|&c| graph.op(c) != core::LIST) { vec![node] } else { children };
    if pairs.len() != dimension {
        return None;
    }
    pairs
        .iter()
        .map(|&p| {
            let ends = graph.children(p).to_vec();
            let &[a, b] = ends.as_slice() else { return None };
            let (lo, hi) = (value(graph, a)?, value(graph, b)?);
            (lo <= hi).then_some(Interval { lo, hi })
        })
        .collect()
}

/// `nsolve_certified(eq, x, guess)`: `list(lo, hi)`, an interval that
/// provably contains exactly one root; with lists of equations, unknowns
/// and guesses the list of such intervals.
pub(super) struct NsolveCertified {
    pub(super) op: OpId,
}

/// `solve_box(eqs, vars, list(list(lo, hi), ...))`: exhaustive certified
/// search; `list(verified, unresolved)`, each a list of boxes; `unresolved`
/// empty proves that the verified boxes are all the roots in the region.
pub(super) struct SolveBox {
    pub(super) op: OpId,
}

fn split_args(
    graph: &Graph,
    equations: NodeId,
    unknowns: NodeId,
) -> Option<(Vec<NodeId>, Vec<SymbolId>, bool)> {
    let scalar = graph.op(unknowns) != core::LIST;
    let eqs: Vec<NodeId> = if graph.op(equations) == core::LIST { graph.children(equations).to_vec() } else { vec![equations] };
    let names: Vec<NodeId> = if scalar { vec![unknowns] } else { graph.children(unknowns).to_vec() };
    let symbols = names.iter().map(|&u| graph.symbol_of(u)).collect::<Option<Vec<_>>>()?;
    Some((eqs, symbols, scalar))
}

impl Kernel for NsolveCertified {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let bindings = cx.env.bindings().to_vec();
        let graph = &mut *cx.graph;
        let &[equation, unknown, guess] = graph.children(node) else { return Outcome::Pass };
        let Some((equations, symbols, scalar)) = split_args(graph, equation, unknown) else { return Outcome::Pass };
        let Some(system) = System::new(graph, &equations, &symbols, &bindings) else { return Outcome::Pass };
        let guess_nodes: Vec<NodeId> = if scalar { vec![guess] } else { graph.children(guess).to_vec() };
        let start: Option<Vec<f64>> = guess_nodes.iter().map(|&g| best(graph, g).and_then(|t| graph.eval(t, cx.env))).collect();
        let Some(mut point) = start.filter(|s| s.len() == symbols.len()) else { return Outcome::Pass };
        // Polish with interval-free Newton steps through the Krawczyk
        // midpoint map: x <- mid(K(x)) converges quadratically.
        for _ in 0..40 {
            let radius = 1e-3;
            let x: Vec<Interval> = point.iter().map(|&g| Interval { lo: g - radius, hi: g + radius }).collect();
            let Some(k) = krawczyk(&system, &x) else { break };
            let next: Vec<f64> = k.iter().map(|i| midpoint(*i)).collect();
            if next.iter().any(|v| !v.is_finite()) {
                break;
            }
            let delta: f64 = next.iter().zip(&point).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
            point = next;
            if delta < 1e-15 {
                break;
            }
        }
        let Some(enclosure) = verify_near(&system, &point) else { return Outcome::Pass };
        if scalar {
            return match enclosure.first() {
                | Some(&i) => Outcome::Equal(pair(graph, i)),
                | None => Outcome::Pass,
            };
        }
        Outcome::Equal(box_node(graph, &enclosure))
    }

    fn revisit(&self) -> bool {
        true
    }
}

impl Kernel for SolveBox {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let bindings = cx.env.bindings().to_vec();
        let graph = &mut *cx.graph;
        let &[equation, unknown, region] = graph.children(node) else { return Outcome::Pass };
        let Some((equations, symbols, _)) = split_args(graph, equation, unknown) else { return Outcome::Pass };
        if equations.len() != symbols.len() {
            return Outcome::Pass;
        }
        let Some(bounds) = read_box(graph, cx.env, region, symbols.len()) else { return Outcome::Pass };
        let Some(system) = System::new(graph, &equations, &symbols, &bindings) else { return Outcome::Pass };
        let found = search(&system, &bounds, 1e-10, 200_000);
        let scalar = symbols.len() == 1;
        let render = |graph: &mut Graph, boxes: &[Vec<Interval>]| -> NodeId {
            let mut sorted: Vec<&Vec<Interval>> = boxes.iter().collect();
            sorted.sort_by(|a, b| {
                a.iter().zip(b.iter()).map(|(x, y)| x.lo.total_cmp(&y.lo)).find(|o| o.is_ne()).unwrap_or(std::cmp::Ordering::Equal)
            });
            let items: Vec<NodeId> = sorted
                .into_iter()
                .map(|b| match (scalar, b.first()) {
                    | (true, Some(&i)) => pair(graph, i),
                    | _ => box_node(graph, b),
                })
                .collect();
            graph.node(core::LIST, &items)
        };
        let verified = render(graph, &found.verified);
        let unresolved = render(graph, &found.unresolved);
        Outcome::Equal(graph.node(core::LIST, &[verified, unresolved]))
    }

    fn revisit(&self) -> bool {
        true
    }
}

/// A rigorous ball for a scalar root near `guess`, when one can be
/// verified.
pub(super) fn certified_ball(
    graph: &mut Graph,
    term_equation: NodeId,
    unknown: SymbolId,
    bindings: &[(SymbolId, f64)],
    root: f64,
) -> Option<Ball> {
    let system = System::new(graph, &[term_equation], &[unknown], bindings)?;
    let enclosure = verify_near(&system, &[root])?;
    let i = *enclosure.first()?;
    let mid = midpoint(i);
    Some(Ball { mid, rad: (i.hi - mid).max(mid - i.lo) })
}

#[cfg(test)]
mod tests {
    use crate::rules::testing::reduce_with;

    fn pairs(text: &str) -> Vec<(f64, f64)> {
        let numbers: Vec<f64> = text
            .replace("list", "")
            .split(|c: char| !(c.is_ascii_digit() || c == '.' || c == '-' || c == 'e' || c == 'E'))
            .filter(|s| !s.is_empty() && *s != "-" && *s != "e")
            .filter_map(|s| s.parse().ok())
            .collect();
        numbers.chunks(2).filter(|c| c.len() == 2).map(|c| (c[0], c[1])).collect()
    }

    #[test]
    fn nsolve_encloses_the_root() {
        let rules = crate::rules::standard();
        let (text, _) = reduce_with(&rules, "nsolve_certified(x^2 - 2, x, 1.4)", &[]);
        let p = pairs(&text);
        assert_eq!(p.len(), 1, "{text}");
        let (lo, hi) = p[0];
        assert!(lo <= 2.0_f64.sqrt() && 2.0_f64.sqrt() <= hi && hi - lo < 1e-12, "{text}");
    }

    #[test]
    fn exhaustive_box_search_finds_all_roots() {
        let rules = crate::rules::standard();
        // cos x = x/3 on [-6, 6]: three roots, certified complete.
        let (text, _) = reduce_with(&rules, "solve_box(cos(x) - x/3, x, list(-6, 6))", &[]);
        assert!(text.ends_with("list())"), "{text}");
        assert!(!pairs(&text).is_empty(), "{text}");
        let (text, _) = reduce_with(&rules, "solve_box(sin(x) - x/5, x, list(-8, 8))", &[]);
        assert_eq!(pairs(&text).len(), 3, "{text}");
        // A circle meets a line in two points.
        let (text, _) = reduce_with(
            &rules,
            "solve_box(list(x^2 + y^2 - 1, y - x), list(x, y), list(list(-2, 2), list(-2, 2)))",
            &[],
        );
        assert_eq!(pairs(&text).len(), 4, "{text}");
        assert!(text.ends_with("list())"), "{text}");
    }
}
