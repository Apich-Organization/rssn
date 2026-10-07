//! Term rewriting systems and Knuth–Bendix completion.
//!
//! Terms are read syntactically: operators, undetermined functions
//! `f(x, y)` and constants are function symbols; the symbols listed as
//! *variables* are rewrite variables. Rules are written `lhs = rhs` and used
//! left to right.
//!
//! | operator | meaning |
//! |---|---|
//! | `knuth_bendix(list(l1 = r1, ...), list(x, y, ...))` | a convergent rewrite system for the equations, or no reduction when completion fails |
//! | `knuth_bendix(equations, variables, list(f1, f2, ...))` | with an explicit precedence `f1 < f2 < ...` |
//! | `rewrite_with(t, list(l = r, ...), list(x, ...))` | the normal form of `t` |
//! | `critical_pairs(list(l = r, ...), list(x, ...))` | the non-joinable critical pairs |
//!
//! Completion is Huet's procedure: orient each equation with the
//! lexicographic path order (LPO) induced by the precedence, add the rule,
//! inter-reduce, and add the normalised critical pairs of the new rule with
//! every rule as new equations, until none is left. An equation that cannot
//! be oriented makes completion fail. The default precedence ranks symbols
//! by arity, then by name.

use std::cmp::Ordering;
use std::collections::HashMap;

use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::OpDescriptor;
use crate::graph::OpFlags;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::RuleError;
use crate::graph::RuleSet;
use crate::graph::Tier;
use crate::graph::op::core;
use crate::graph::rule::Installer;

/// The rewriting rule set.
#[must_use]
pub fn rewriting() -> RuleSet {
    RuleSet::new("rewriting", install)
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum Request {
    Complete,
    Rewrite,
    CriticalPairs,
}

fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    for (name, arity, request) in [
        ("knuth_bendix", Arity::Variadic, Request::Complete),
        ("rewrite_with", Arity::Fixed(3), Request::Rewrite),
        ("critical_pairs", Arity::Fixed(2), Request::CriticalPairs),
    ] {
        let op = i.op(OpDescriptor::new(name, arity).flags(OpFlags::HEAVY).cost(100))?;
        i.kernel(&format!("rewriting/{name}"), Tier::Reduce, Rewriting { op, request });
    }
    Ok(())
}

/// A function symbol.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
enum Head {
    /// A registered operator.
    Op(OpId),
    /// An undetermined function, by its name symbol.
    Fun(NodeId),
    /// A constant leaf: a number, or a symbol that is not a variable.
    Leaf(NodeId),
}

/// A first-order term.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
enum T {
    Var(u32),
    App(Head, Vec<Self>),
}

type Subst = HashMap<u32, T>;
type Rule = (T, T);

/// Converts between graph terms and [`T`].
struct Codec {
    variables: Vec<NodeId>,
}

impl Codec {
    fn read(
        &self,
        graph: &Graph,
        node: NodeId,
    ) -> T {
        if let Some(i) = self.variables.iter().position(|&v| v == node) {
            return T::Var(u32::try_from(i).unwrap_or(u32::MAX));
        }
        let op = graph.op(node);
        let children = graph.children(node);
        if op == core::APPLY
            && let Some((&f, args)) = children.split_first()
                && graph.as_symbol(f).is_some() {
                    return T::App(Head::Fun(f), args.iter().map(|&a| self.read(graph, a)).collect());
                }
        if children.is_empty() && (op == core::LIT || op == core::SYM) {
            return T::App(Head::Leaf(node), Vec::new());
        }
        T::App(Head::Op(op), children.iter().map(|&a| self.read(graph, a)).collect())
    }

    fn write(
        &self,
        graph: &mut Graph,
        t: &T,
    ) -> Option<NodeId> {
        Some(match t {
            | T::Var(i) => *self.variables.get(*i as usize)?,
            | T::App(Head::Leaf(n), _) => *n,
            | T::App(Head::Fun(f), args) => {
                let mut children = vec![*f];
                for a in args {
                    children.push(self.write(graph, a)?);
                }
                graph.node(core::APPLY, &children)
            },
            | T::App(Head::Op(op), args) => {
                let children = args.iter().map(|a| self.write(graph, a)).collect::<Option<Vec<_>>>()?;
                graph.try_node(*op, &children)?
            },
        })
    }
}

/// The precedence on function symbols.
struct Precedence {
    explicit: HashMap<Head, usize>,
    names: HashMap<Head, String>,
}

impl Precedence {
    fn cmp(
        &self,
        f: Head,
        f_arity: usize,
        g: Head,
        g_arity: usize,
    ) -> Ordering {
        if f == g {
            return Ordering::Equal;
        }
        match (self.explicit.get(&f), self.explicit.get(&g)) {
            | (Some(a), Some(b)) => a.cmp(b),
            // Listed symbols rank above unlisted ones.
            | (Some(_), None) => Ordering::Greater,
            | (None, Some(_)) => Ordering::Less,
            | (None, None) => f_arity
                .cmp(&g_arity)
                .then_with(|| self.names.get(&f).cmp(&self.names.get(&g))),
        }
    }
}

fn occurs(
    v: u32,
    t: &T,
) -> bool {
    match t {
        | T::Var(w) => *w == v,
        | T::App(_, args) => args.iter().any(|a| occurs(v, a)),
    }
}

/// `s >_lpo t`.
fn lpo_greater(
    prec: &Precedence,
    s: &T,
    t: &T,
) -> bool {
    match (s, t) {
        | (T::Var(_), _) => false,
        | (_, T::Var(v)) => s != t && occurs(*v, s),
        | (T::App(f, ss), T::App(g, ts)) => {
            // Some argument of s is >= t.
            if ss.iter().any(|si| si == t || lpo_greater(prec, si, t)) {
                return true;
            }
            let dominates_all = || ts.iter().all(|tj| lpo_greater(prec, s, tj));
            match prec.cmp(*f, ss.len(), *g, ts.len()) {
                | Ordering::Greater => dominates_all(),
                | Ordering::Equal if ss.len() == ts.len() => {
                    // Lexicographic on arguments.
                    let first = ss.iter().zip(ts).find(|(a, b)| a != b);
                    first.is_some_and(|(a, b)| lpo_greater(prec, a, b)) && dominates_all()
                },
                | _ => false,
            }
        },
    }
}

fn apply(
    t: &T,
    sigma: &Subst,
) -> T {
    match t {
        | T::Var(v) => sigma.get(v).map_or_else(|| t.clone(), |s| apply(s, sigma)),
        | T::App(f, args) => T::App(*f, args.iter().map(|a| apply(a, sigma)).collect()),
    }
}

/// `t` with the variables bound by a matching substitution replaced, in
/// one pass (bound values are not substituted again).
fn instantiate(
    t: &T,
    sigma: &Subst,
) -> T {
    match t {
        | T::Var(v) => sigma.get(v).cloned().unwrap_or_else(|| t.clone()),
        | T::App(f, args) => T::App(*f, args.iter().map(|a| instantiate(a, sigma)).collect()),
    }
}

/// One-way matching of `pattern` against `t`.
fn matches(
    pattern: &T,
    t: &T,
    sigma: &mut Subst,
) -> bool {
    match pattern {
        | T::Var(v) => match sigma.get(v) {
            | Some(bound) => bound == t,
            | None => {
                sigma.insert(*v, t.clone());
                true
            },
        },
        | T::App(f, ps) => match t {
            | T::App(g, ts) if f == g && ps.len() == ts.len() => {
                ps.iter().zip(ts).all(|(p, u)| matches(p, u, sigma))
            },
            | _ => false,
        },
    }
}

/// The most general unifier of `a` and `b`.
fn unify(
    a: &T,
    b: &T,
) -> Option<Subst> {
    let mut sigma = Subst::new();
    let mut stack = vec![(a.clone(), b.clone())];
    while let Some((x, y)) = stack.pop() {
        let (x, y) = (apply(&x, &sigma), apply(&y, &sigma));
        if x == y {
            continue;
        }
        match (x, y) {
            | (T::Var(v), t) | (t, T::Var(v)) => {
                if occurs(v, &t) {
                    return None;
                }
                sigma.insert(v, t);
            },
            | (T::App(f, xs), T::App(g, ys)) => {
                if f != g || xs.len() != ys.len() {
                    return None;
                }
                stack.extend(xs.into_iter().zip(ys));
            },
        }
    }
    Some(sigma)
}

/// One rewrite step at the outermost-leftmost redex.
fn step(
    t: &T,
    rules: &[Rule],
) -> Option<T> {
    for (l, r) in rules {
        let mut sigma = Subst::new();
        if matches(l, t, &mut sigma) {
            return Some(instantiate(r, &sigma));
        }
    }
    if let T::App(f, args) = t {
        for (i, a) in args.iter().enumerate() {
            if let Some(b) = step(a, rules) {
                let mut new_args = args.clone();
                new_args[i] = b;
                return Some(T::App(*f, new_args));
            }
        }
    }
    None
}

/// The normal form of `t`, or `None` after too many steps.
fn normalize(
    t: &T,
    rules: &[Rule],
) -> Option<T> {
    let mut current = t.clone();
    for _ in 0..10_000 {
        match step(&current, rules) {
            | Some(next) => current = next,
            | None => return Some(current),
        }
    }
    None
}

fn size(t: &T) -> usize {
    match t {
        | T::Var(_) => 1,
        | T::App(_, args) => args.iter().map(size).sum::<usize>() + 1,
    }
}

fn shift(
    t: &T,
    by: u32,
) -> T {
    match t {
        | T::Var(v) => T::Var(v.saturating_add(by)),
        | T::App(f, args) => T::App(*f, args.iter().map(|a| shift(a, by)).collect()),
    }
}

fn max_var(t: &T) -> u32 {
    match t {
        | T::Var(v) => *v,
        | T::App(_, args) => args.iter().map(max_var).max().unwrap_or(0),
    }
}

/// Non-variable subterm positions of `t`, with the subterm.
fn positions(t: &T) -> Vec<(Vec<usize>, T)> {
    let mut out = Vec::new();
    let mut stack = vec![(Vec::new(), t.clone())];
    while let Some((path, s)) = stack.pop() {
        if let T::App(_, args) = &s {
            for (i, a) in args.iter().enumerate() {
                let mut p = path.clone();
                p.push(i);
                stack.push((p, a.clone()));
            }
            out.push((path, s));
        }
    }
    out
}

fn replace_at(
    t: &T,
    path: &[usize],
    with: &T,
) -> T {
    match (path.split_first(), t) {
        | (None, _) => with.clone(),
        | (Some((&i, rest)), T::App(f, args)) => {
            let mut new_args = args.clone();
            if let Some(slot) = new_args.get_mut(i) {
                *slot = replace_at(slot, rest, with);
            }
            T::App(*f, new_args)
        },
        | (Some(_), T::Var(_)) => t.clone(),
    }
}

/// Critical pairs of `outer` with `inner` overlapping into it.
fn overlaps(
    outer: &Rule,
    inner: &Rule,
    same: bool,
) -> Vec<(T, T)> {
    let offset = max_var(&outer.0).max(max_var(&outer.1)).saturating_add(1);
    let (l2, r2) = (shift(&inner.0, offset), shift(&inner.1, offset));
    let mut out = Vec::new();
    for (path, sub) in positions(&outer.0) {
        if same && path.is_empty() {
            continue;
        }
        if let Some(sigma) = unify(&sub, &l2) {
            let left = apply(&replace_at(&outer.0, &path, &r2), &sigma);
            let right = apply(&outer.1, &sigma);
            out.push((left, right));
        }
    }
    out
}

/// Huet's completion; `None` when an equation cannot be oriented or the
/// limits are exceeded.
fn complete(
    equations: Vec<(T, T)>,
    prec: &Precedence,
) -> Option<Vec<Rule>> {
    let mut pending = equations;
    let mut rules: Vec<Rule> = Vec::new();
    let mut work = 0_usize;
    // Fair selection: always the smallest pending equation.
    while let Some(pick) = (0..pending.len()).min_by_key(|&i| size(&pending[i].0) + size(&pending[i].1)) {
        let (s, t) = pending.swap_remove(pick);
        work = work.saturating_add(1);
        if work > 5_000 || rules.len() > 200 {
            return None;
        }
        let (s, t) = (normalize(&s, &rules)?, normalize(&t, &rules)?);
        if s == t {
            continue;
        }
        let rule = if lpo_greater(prec, &s, &t) {
            (s, t)
        } else if lpo_greater(prec, &t, &s) {
            (t, s)
        } else {
            return None;
        };
        let new = [rule.clone()];
        // Inter-reduce: rules whose left side the new rule rewrites become
        // equations again; right sides are normalised.
        let mut kept = Vec::with_capacity(rules.len());
        for (l, r) in rules {
            if step(&l, &new).is_some() {
                pending.push((l, r));
            } else {
                kept.push((l, r));
            }
        }
        rules = kept;
        rules.push(rule.clone());
        for i in 0..rules.len() {
            let r = normalize(&rules[i].1, &rules)?;
            rules[i].1 = r;
        }
        for other in rules.clone() {
            let same = other == rule;
            pending.extend(overlaps(&rule, &other, same));
            if !same {
                pending.extend(overlaps(&other, &rule, false));
            }
        }
    }
    Some(rules)
}

struct Rewriting {
    op: OpId,
    request: Request,
}

fn head_name(
    graph: &Graph,
    head: Head,
) -> String {
    match head {
        | Head::Op(op) => graph.ops().get(op).name.to_string(),
        | Head::Fun(n) | Head::Leaf(n) => graph.display(n),
    }
}

fn collect_heads(
    t: &T,
    out: &mut Vec<Head>,
) {
    if let T::App(f, args) = t {
        if !out.contains(f) {
            out.push(*f);
        }
        for a in args {
            collect_heads(a, out);
        }
    }
}

/// Reads `list(l = r, ...)` as pairs.
fn read_equations(
    graph: &Graph,
    codec: &Codec,
    list: NodeId,
) -> Option<Vec<(T, T)>> {
    if graph.op(list) != core::LIST {
        return None;
    }
    graph
        .children(list)
        .iter()
        .map(|&e| match *graph.children(e) {
            | [l, r] if graph.op(e) == core::EQ => Some((codec.read(graph, l), codec.read(graph, r))),
            | _ => None,
        })
        .collect()
}

/// Renames the variables of `t` to `0, 1, ...` in order of appearance.
fn canonical(
    t: &T,
    names: &mut Vec<u32>,
) -> T {
    match t {
        | T::Var(v) => {
            let index = names.iter().position(|w| w == v).unwrap_or_else(|| {
                names.push(*v);
                names.len() - 1
            });
            T::Var(u32::try_from(index).unwrap_or(u32::MAX))
        },
        | T::App(f, args) => T::App(*f, args.iter().map(|a| canonical(a, names)).collect()),
    }
}

fn write_rules(
    graph: &mut Graph,
    codec: &Codec,
    rules: &[Rule],
) -> Option<NodeId> {
    let mut items = Vec::with_capacity(rules.len());
    let mut codec = Codec { variables: codec.variables.clone() };
    for (l, r) in rules {
        let mut names = Vec::new();
        let (l, r) = (canonical(l, &mut names), canonical(r, &mut names));
        // More variables than were given: name the extra ones afresh.
        while codec.variables.len() < names.len() {
            let fresh = graph.interner_mut().fresh_symbol("v");
            codec.variables.push(graph.symbol_node(fresh));
        }
        let (l, r) = (codec.write(graph, &l)?, codec.write(graph, &r)?);
        items.push(graph.node(core::EQ, &[l, r]));
    }
    Some(graph.node(core::LIST, &items))
}

impl Kernel for Rewriting {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let graph = &mut *cx.graph;
        let args = graph.children(node).to_vec();
        let variables_of = |graph: &Graph, list: NodeId| -> Option<Vec<NodeId>> {
            (graph.op(list) == core::LIST).then(|| graph.children(list).to_vec())
        };
        let result = match (self.request, args.as_slice()) {
            | (Request::Complete, [equations, vars, rest @ ..]) if rest.len() <= 1 => (|| {
                let codec = Codec { variables: variables_of(graph, *vars)? };
                let equations = read_equations(graph, &codec, *equations)?;
                let mut heads = Vec::new();
                for (l, r) in &equations {
                    collect_heads(l, &mut heads);
                    collect_heads(r, &mut heads);
                }
                let mut explicit = HashMap::new();
                if let Some(&order) = rest.first() {
                    let symbols = variables_of(graph, order)?;
                    for (rank, s) in symbols.iter().enumerate() {
                        let name = graph.display(*s);
                        for &h in &heads {
                            if head_name(graph, h) == name {
                                explicit.insert(h, rank);
                            }
                        }
                    }
                }
                let names = heads.iter().map(|&h| (h, head_name(graph, h))).collect();
                let prec = Precedence { explicit, names };
                let rules = complete(equations, &prec)?;
                write_rules(graph, &codec, &rules)
            })(),
            | (Request::Rewrite, &[t, rules, vars]) => (|| {
                let codec = Codec { variables: variables_of(graph, vars)? };
                let rules = read_equations(graph, &codec, rules)?;
                // Variables of the term itself are ordinary constants.
                let plain = Codec { variables: Vec::new() };
                let term = plain.read(graph, t);
                let normal = normalize(&term, &rules)?;
                plain.write(graph, &normal)
            })(),
            | (Request::CriticalPairs, &[rules, vars]) => (|| {
                let codec = Codec { variables: variables_of(graph, vars)? };
                let rules = read_equations(graph, &codec, rules)?;
                let mut pairs = Vec::new();
                for (i, a) in rules.iter().enumerate() {
                    for (j, b) in rules.iter().enumerate() {
                        for (s, t) in overlaps(a, b, i == j) {
                            let (s, t) = (normalize(&s, &rules)?, normalize(&t, &rules)?);
                            if s != t && !pairs.contains(&(s.clone(), t.clone())) {
                                pairs.push((s, t));
                            }
                        }
                    }
                }
                write_rules(graph, &codec, &pairs)
            })(),
            | _ => None,
        };
        result.map_or(Outcome::Pass, Outcome::Pinned)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rules::testing::reduce_with;
    use crate::rules::testing::simplify;

    const GROUP: &str = "list(op(u, x) = x, op(inv(x), x) = u, op(op(x, y), z) = op(x, op(y, z)))";

    #[test]
    fn group_axioms_complete_to_ten_rules() {
        let sets = [rewriting()];
        let (rules, reduced) = reduce_with(&sets, &format!("knuth_bendix({GROUP}, list(x, y, z), list(u, op, inv))"), &[]);
        assert!(reduced, "{rules}");
        assert_eq!(rules.matches(" = ").count(), 10, "{rules}");
        let normal = |t: &str| simplify(&sets, &format!("rewrite_with({t}, {rules}, list(x, y, z))"));
        assert_eq!(normal("op(inv(op(a, b)), op(a, b))"), "u");
        assert_eq!(normal("inv(inv(a))"), "a");
        assert_eq!(normal("inv(op(a, b))"), "op(inv(b), inv(a))");
        assert_eq!(normal("op(a, op(inv(a), b))"), "b");
        assert_eq!(normal("inv(u)"), "u");
    }

    #[test]
    fn rewriting_and_critical_pairs() {
        let sets = [rewriting()];
        assert_eq!(simplify(&sets, "rewrite_with(f(f(a)), list(f(f(x)) = x), list(x))"), "a");
        assert_eq!(simplify(&sets, "rewrite_with(g(f(f(f(b)))), list(f(f(x)) = x), list(x))"), "g(f(b))");
        // f(f(f(x))) overlaps with itself: f(x) and f(x) — joinable, so no pairs.
        assert_eq!(simplify(&sets, "critical_pairs(list(f(f(x)) = x), list(x))"), "list()");
        let pairs = simplify(&sets, "critical_pairs(list(f(g(x)) = a, g(h(y)) = b), list(x, y))");
        assert!(pairs.contains("f(b)") && pairs.contains('a'), "{pairs}");
        // An equation that no path order can orient (commutativity) fails.
        let (_, reduced) = reduce_with(&sets, "knuth_bendix(list(m(x, y) = m(y, x)), list(x, y))", &[]);
        assert!(!reduced);
    }
}
