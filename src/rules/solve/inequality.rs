//! Inequalities and systems of conditions over the reals in one unknown.
//!
//! An atom `lt(f, g)` (likewise `le`, `gt`, `ge`) is solved by locating the
//! boundary points of `f - g` — its real roots, the zeros of its
//! denominator and the zeros of the arguments of logarithms and roots —
//! and testing the sign on each interval between them. The result is a
//! union of intervals; `list(...)` and `and(...)` intersect their members,
//! `or(...)` takes the union, and an equation `f = g` contributes its
//! roots as points. The answer is a formula in `x`: an `or` of
//! `and(lt(a, x), lt(x, b))` pieces (`le` where an endpoint is included),
//! an equation `x = a` for a single point, `true` or `false`.

use super::as_expression;
use super::difference;
use super::heuristics::nodes_with;
use super::solve_for;
use crate::graph::op::core;
use crate::graph::Env;
use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::OpId;
use crate::rules::poly::best;
use crate::rules::poly::ratio;
use crate::rules::poly::repr::to_term;
use crate::rules::poly::repr::Gens;
use crate::rules::poly::repr::Limits;

/// One end of an interval: its numeric value and the term that writes it.
#[derive(Copy, Clone)]
struct Bound {
    value: f64,
    node: NodeId,
    closed: bool,
}

/// An interval; `None` ends are infinite.
#[derive(Copy, Clone)]
struct Interval {
    lo: Option<Bound>,
    hi: Option<Bound>,
}

type Set = Vec<Interval>;

fn full() -> Set {
    vec![Interval { lo: None, hi: None }]
}

/// Whether `a` starts at or after `b` starts (lower ends).
fn lo_max(
    a: Option<Bound>,
    b: Option<Bound>,
) -> Option<Bound> {
    match (a, b) {
        | (None, x) | (x, None) => x,
        | (Some(x), Some(y)) => {
            if (x.value - y.value).abs() <= 1e-12 * (1.0 + x.value.abs()) {
                Some(Bound { closed: x.closed && y.closed, ..x })
            } else if x.value > y.value {
                Some(x)
            } else {
                Some(y)
            }
        },
    }
}

fn hi_min(
    a: Option<Bound>,
    b: Option<Bound>,
) -> Option<Bound> {
    match (a, b) {
        | (None, x) | (x, None) => x,
        | (Some(x), Some(y)) => {
            if (x.value - y.value).abs() <= 1e-12 * (1.0 + x.value.abs()) {
                Some(Bound { closed: x.closed && y.closed, ..x })
            } else if x.value < y.value {
                Some(x)
            } else {
                Some(y)
            }
        },
    }
}

fn nonempty(i: &Interval) -> bool {
    match (i.lo, i.hi) {
        | (Some(l), Some(h)) => {
            let tolerance = 1e-12 * (1.0 + l.value.abs());
            h.value > l.value + tolerance || ((h.value - l.value).abs() <= tolerance && l.closed && h.closed)
        },
        | _ => true,
    }
}

fn intersect(
    a: &Set,
    b: &Set,
) -> Set {
    let mut out = Vec::new();
    for x in a {
        for y in b {
            let i = Interval { lo: lo_max(x.lo, y.lo), hi: hi_min(x.hi, y.hi) };
            if nonempty(&i) {
                out.push(i);
            }
        }
    }
    sort(&mut out);
    out
}

fn sort(set: &mut Set) {
    set.sort_by(|a, b| {
        let key = |i: &Interval| i.lo.map_or(f64::NEG_INFINITY, |l| l.value);
        key(a).total_cmp(&key(b))
    });
}

fn union(
    a: &Set,
    b: &Set,
) -> Set {
    let mut all: Set = a.iter().chain(b).copied().collect();
    sort(&mut all);
    let mut out: Set = Vec::new();
    for i in all {
        match out.last_mut() {
            | Some(last) => {
                let touches = match (last.hi, i.lo) {
                    | (None, _) | (_, None) => true,
                    | (Some(h), Some(l)) => {
                        let tolerance = 1e-12 * (1.0 + h.value.abs());
                        l.value < h.value - tolerance || ((l.value - h.value).abs() <= tolerance && (h.closed || l.closed))
                    },
                };
                if touches {
                    let hi = match (last.hi, i.hi) {
                        | (None, _) | (_, None) => None,
                        | (Some(x), Some(y)) => {
                            if (x.value - y.value).abs() <= 1e-12 * (1.0 + x.value.abs()) {
                                Some(Bound { closed: x.closed || y.closed, ..x })
                            } else if x.value > y.value {
                                Some(x)
                            } else {
                                Some(y)
                            }
                        },
                    };
                    last.hi = hi;
                } else {
                    out.push(i);
                }
            },
            | None => out.push(i),
        }
    }
    out
}

/// The set where the comparison `kind` (0 `lt`, 1 `le`, 2 `gt`, 3 `ge`)
/// holds for `f` against zero.
fn comparison_set(
    graph: &mut Graph,
    f: NodeId,
    x: NodeId,
    kind: usize,
) -> Option<Set> {
    let symbol = graph.symbol_of(x)?;
    let holds = |v: f64| match kind {
        | 0 => v < 0.0,
        | 1 => v <= 1e-12,
        | 2 => v > 0.0,
        | _ => v >= -1e-12,
    };
    let strict = kind == 0 || kind == 2;
    let eval_at = |graph: &Graph, at: f64| -> Option<f64> {
        let mut env = Env::numeric(0.0);
        env.bind(symbol, at);
        graph.eval(f, &env).filter(|v| v.is_finite())
    };
    // Boundary points: (value, node).
    let mut boundary: Vec<(f64, NodeId)> = Vec::new();
    let mut add = |graph: &mut Graph, roots: Vec<NodeId>| {
        for r in roots {
            if let Some(v) = graph.eval(r, &Env::numeric(0.0)).filter(|v| v.is_finite()) {
                boundary.push((v, r));
            }
        }
    };
    let f_roots = solve_for(graph, f, x, 0)?;
    add(graph, f_roots);
    let mut gens = Gens::default();
    gens.index(graph, x);
    if let Some(fraction) = ratio(graph, &mut gens, f, Limits::default())
        && fraction.denom.as_constant().is_none() {
            let denominator = to_term(graph, &gens, &fraction.denom);
            let poles = solve_for(graph, denominator, x, 0)?;
            add(graph, poles);
        }
    // Domain boundaries of logarithms and roots.
    let ln = graph.ops().lookup("ln");
    let sqrt = graph.ops().lookup("sqrt");
    let domain = nodes_with(graph, f, |g, n| {
        let op = g.op(n);
        g.depends_on(g.find(n), symbol)
            && (Some(op) == ln
                || Some(op) == sqrt
                || (op == core::POW && g.children(n).get(1).and_then(|&e| g.number_of(e)).is_some_and(|e| !e.is_integer())))
    });
    for n in domain {
        if let Some(&arg) = graph.children(n).first()
            && let Some(roots) = solve_for(graph, arg, x, 0) {
                add(graph, roots);
            }
    }
    boundary.sort_by(|a, b| a.0.total_cmp(&b.0));
    boundary.dedup_by(|a, b| (a.0 - b.0).abs() <= 1e-12 * (1.0 + a.0.abs()));
    // Items alternate: region 0, point 0, region 1, ..., region n.
    let count = boundary.len();
    let mut included: Vec<bool> = Vec::new();
    for i in 0..=count {
        let sample = match (i.checked_sub(1).map(|j| boundary[j].0), boundary.get(i).map(|b| b.0)) {
            | (None, None) => 0.0,
            | (None, Some(b)) => b - 1.0 - b.abs(),
            | (Some(a), None) => a + 1.0 + a.abs(),
            | (Some(a), Some(b)) => f64::midpoint(a, b),
        };
        included.push(eval_at(graph, sample).is_some_and(holds));
        if let Some(&(value, _)) = boundary.get(i) {
            // The boundary point itself: defined and satisfying a
            // non-strict comparison.
            let ok = !strict && eval_at(graph, value).is_some_and(|v| v.abs() <= 1e-9 && holds(v));
            included.push(ok);
        }
    }
    let mut out: Set = Vec::new();
    let mut index = 0;
    while index < included.len() {
        if !included[index] {
            index += 1;
            continue;
        }
        let start = index;
        while index + 1 < included.len() && included[index + 1] {
            index += 1;
        }
        let end = index;
        index += 1;
        // Item k: even = region k/2, odd = point (k-1)/2.
        let lo = if start % 2 == 1 {
            let (value, node) = boundary[(start - 1) / 2];
            Some(Bound { value, node, closed: true })
        } else if start == 0 {
            None
        } else {
            let (value, node) = boundary[start / 2 - 1];
            Some(Bound { value, node, closed: false })
        };
        let hi = if end % 2 == 1 {
            let (value, node) = boundary[(end - 1) / 2];
            Some(Bound { value, node, closed: true })
        } else if end / 2 >= count {
            None
        } else {
            let (value, node) = boundary[end / 2];
            Some(Bound { value, node, closed: false })
        };
        out.push(Interval { lo, hi });
    }
    Some(out)
}

/// The outcome of rewriting a condition on `floor(u)` / `ceil(u)`.
enum Rounded {
    /// Never satisfied.
    Never,
    /// Satisfied exactly where this condition on `u` holds.
    When(NodeId),
}

/// A comparison or equation between `floor(u)` / `ceil(u)` and a number,
/// rewritten as a condition on `u` (an interval): `floor(u) < c` is
/// `u < ceil(c)`, `floor(u) = c` is `c <= u < c + 1` for an integer `c`,
/// and so on. `None` when `node` is not of that shape.
fn floor_condition(
    graph: &mut Graph,
    node: NodeId,
) -> Option<Rounded> {
    let (floor, ceil) = (graph.ops().lookup("floor")?, graph.ops().lookup("ceil")?);
    let op = graph.op(node);
    let names = ["lt", "le", "gt", "ge"];
    let kind = names.iter().position(|n| graph.ops().lookup(n) == Some(op));
    let is_eq = op == core::EQ;
    if kind.is_none() && !is_eq {
        return None;
    }
    let &[a, b] = graph.children(node) else {
        return None;
    };
    let rounding = |g: &Graph, n: NodeId| (g.op(n) == floor || g.op(n) == ceil) && g.children(n).len() == 1;
    // (rounding node, constant, kind seen from the rounding side)
    let (r, c, kind) = if rounding(graph, a) {
        (a, b, kind)
    } else if rounding(graph, b) {
        (b, a, kind.map(|k| [2, 3, 0, 1][k]))
    } else {
        return None;
    };
    let value = graph.eval(c, &Env::numeric(0.0)).filter(|v| v.is_finite())?;
    let u = *graph.children(r).first()?;
    let is_floor = graph.op(r) == floor;
    let (lt, le, gt, ge) = (
        graph.ops().lookup("lt")?,
        graph.ops().lookup("le")?,
        graph.ops().lookup("gt")?,
        graph.ops().lookup("ge")?,
    );
    let and = graph.ops().lookup("and")?;
    let int = |graph: &mut Graph, v: f64| graph.int(v as i64);
    let (lo, hi) = (value.floor(), value.ceil());
    match kind {
        | None => {
            if (value - lo).abs() > 0.0 {
                return Some(Rounded::Never);
            }
            let (low, high) = (int(graph, lo), int(graph, lo + 1.0));
            let (below, above) = if is_floor {
                (graph.node(ge, &[u, low]), graph.node(lt, &[u, high]))
            } else {
                let previous = int(graph, lo - 1.0);
                let top = int(graph, lo);
                (graph.node(gt, &[u, previous]), graph.node(le, &[u, top]))
            };
            Some(Rounded::When(graph.node(and, &[below, above])))
        },
        | Some(k) => {
            // The bound on `u` and the comparison it takes.
            let (bound, cmp) = match (is_floor, k) {
                | (true, 0) => (hi, lt),
                | (true, 1) => (lo + 1.0, lt),
                | (true, 2) => (lo + 1.0, ge),
                | (true, _) => (hi, ge),
                | (false, 0) => (hi - 1.0, le),
                | (false, 1) => (lo, le),
                | (false, 2) => (lo, gt),
                | (false, _) => (hi - 1.0, gt),
            };
            let _ = (le, gt);
            let b = int(graph, bound);
            Some(Rounded::When(graph.node(cmp, &[u, b])))
        },
    }
}

/// The set described by the condition `node`.
fn condition_set(
    graph: &mut Graph,
    node: NodeId,
    x: NodeId,
) -> Option<Set> {
    match floor_condition(graph, node) {
        | Some(Rounded::When(rewritten)) => return condition_set(graph, rewritten, x),
        | Some(Rounded::Never) => return Some(Vec::new()),
        | None => {},
    }
    let op = graph.op(node);
    let names = ["lt", "le", "gt", "ge"];
    let kinds: Vec<Option<OpId>> = names.iter().map(|n| graph.ops().lookup(n)).collect();
    if let Some(kind) = kinds.iter().position(|&o| o == Some(op)) {
        let &[lhs, rhs] = graph.children(node) else {
            return None;
        };
        let f = difference(graph, lhs, rhs);
        let f = best(graph, f)?;
        return comparison_set(graph, f, x, kind);
    }
    if graph.ops().lookup("ne") == Some(op) {
        let &[lhs, rhs] = graph.children(node) else {
            return None;
        };
        let f = difference(graph, lhs, rhs);
        let f = best(graph, f)?;
        let mut points: Vec<(f64, NodeId)> = Vec::new();
        for r in solve_for(graph, f, x, 0)? {
            if let Some(v) = graph.eval(r, &Env::numeric(0.0)).filter(|v| v.is_finite()) {
                points.push((v, r));
            }
        }
        points.sort_by(|a, b| a.0.total_cmp(&b.0));
        points.dedup_by(|a, b| (a.0 - b.0).abs() <= 1e-12 * (1.0 + a.0.abs()));
        let mut set: Set = Vec::new();
        let mut previous: Option<Bound> = None;
        for &(value, node) in &points {
            let b = Bound { value, node, closed: false };
            set.push(Interval { lo: previous, hi: Some(b) });
            previous = Some(b);
        }
        set.push(Interval { lo: previous, hi: None });
        return Some(set);
    }
    let and = graph.ops().lookup("and");
    let or = graph.ops().lookup("or");
    if op == core::LIST || Some(op) == and {
        let mut set = full();
        for c in graph.children(node).to_vec() {
            set = intersect(&set, &condition_set(graph, c, x)?);
        }
        return Some(set);
    }
    if Some(op) == or {
        let mut set: Set = Vec::new();
        for c in graph.children(node).to_vec() {
            set = union(&set, &condition_set(graph, c, x)?);
        }
        return Some(set);
    }
    // An equation: its roots.
    let expr = as_expression(graph, node);
    let roots = solve_for(graph, expr, x, 0)?;
    let mut set: Set = Vec::new();
    for r in roots {
        let value = graph.eval(r, &Env::numeric(0.0)).filter(|v| v.is_finite())?;
        let b = Bound { value, node: r, closed: true };
        set = union(&set, &vec![Interval { lo: Some(b), hi: Some(b) }]);
    }
    Some(set)
}

/// Whether `node` is a comparison, or a conjunction or disjunction of
/// conditions.
pub(super) fn is_condition(
    graph: &Graph,
    node: NodeId,
) -> bool {
    let op = graph.op(node);
    let rounded = |n: NodeId| {
        ["floor", "ceil"].iter().any(|name| graph.ops().lookup(name) == Some(graph.op(n)))
    };
    (op == core::EQ && graph.children(node).iter().any(|&c| rounded(c)))
        || ["lt", "le", "gt", "ge", "ne", "and", "or"].iter().any(|n| graph.ops().lookup(n) == Some(op))
        || (op == core::LIST
            && !graph.children(node).is_empty()
            && graph.children(node).iter().any(|&c| is_condition(graph, c)))
}

fn formula(
    graph: &mut Graph,
    set: &Set,
    x: NodeId,
) -> Option<NodeId> {
    let (lt, le) = (graph.ops().lookup("lt")?, graph.ops().lookup("le")?);
    let (and, or) = (graph.ops().lookup("and")?, graph.ops().lookup("or")?);
    let mut formulas = Vec::new();
    for i in set {
        if let (Some(a), Some(b)) = (i.lo, i.hi)
            && a.closed && b.closed && (a.value - b.value).abs() <= 1e-12 * (1.0 + a.value.abs()) {
                formulas.push(graph.node(core::EQ, &[x, a.node]));
                continue;
            }
        let mut conditions = Vec::new();
        if let Some(a) = i.lo {
            conditions.push(graph.node(if a.closed { le } else { lt }, &[a.node, x]));
        }
        if let Some(b) = i.hi {
            conditions.push(graph.node(if b.closed { le } else { lt }, &[x, b.node]));
        }
        formulas.push(match conditions.as_slice() {
            | [] => graph.node(graph.ops().lookup("true")?, &[]),
            | [only] => *only,
            | _ => graph.node(and, &conditions),
        });
    }
    Some(match formulas.as_slice() {
        | [] => graph.node(graph.ops().lookup("false")?, &[]),
        | [only] => *only,
        | _ => graph.node(or, &formulas),
    })
}

/// Solves the condition `node` for `x` as a formula.
pub(super) fn solve_condition(
    graph: &mut Graph,
    node: NodeId,
    x: NodeId,
) -> Option<NodeId> {
    let set = condition_set(graph, node, x)?;
    formula(graph, &set, x)
}
