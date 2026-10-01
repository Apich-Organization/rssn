//! The heuristic scheduler: decides which rules run where and when.
//!
//! Plain equality saturation applies every rule everywhere until nothing
//! changes. That does not survive contact with a computer algebra workload,
//! so the [`Engine`] departs from it in five ways:
//!
//! 1. **Relevance.** Only classes reachable from the roots of the current
//!    request are considered; a long-lived graph does not slow a run down.
//! 2. **Tiers.** Heavy operators are reduced first and to a fixpoint
//!    ([`Tier::Reduce`]); simplification follows ([`Tier::Normalize`]);
//!    structure-changing identities ([`Tier::Explore`]) run last, under a
//!    node budget.
//! 3. **Windows.** Simplifying rewrites and window passes run destructively
//!    inside [`TreeWindow`]s, so the bulk of algebraic clean-up never
//!    creates e-nodes. An optional beam search uses the exploring rewrites
//!    as moves.
//! 4. **Back-off.** A rule that matches explosively is banned for a number
//!    of iterations that doubles each time.
//! 5. **Goals.** A run stops as soon as its [`Goal`] is reached — for
//!    example when the root has a closed form, or a numeric witness.
//!
//! When the environment asks for numerics, the engine also propagates
//! numeric witnesses bottom-up through every operator that declares scalar
//! semantics, so numeric kernels higher up find their inputs evaluated.

use std::collections::BTreeMap;
use std::collections::HashMap;
use std::collections::HashSet;
use std::sync::Arc;
use std::time::Duration;
use std::time::Instant;

use super::extract::ClosedForm;
use super::extract::Extractor;
use super::extract::SizeCost;
use super::extract::reachable;
use super::id::NodeId;
use super::id::OpId;
use super::op::OpFlags;
use super::rule::Action;
use super::rule::Cx;
use super::rule::Env;
use super::rule::Outcome;
use super::rule::Program;
use super::rule::Rewrite;
use super::rule::RuleError;
use super::rule::RuleSet;
use super::rule::Tier;
use super::store::Ball;
use super::store::Graph;
use super::window::TreeWindow;
use super::window::WindowPass;

/// Resource limits and heuristic knobs of a run.
#[derive(Clone, Debug)]
pub struct Budget {
    /// Maximum number of outer iterations.
    pub max_iterations: usize,
    /// Exploring rules stop firing once the graph holds this many nodes.
    pub max_nodes: usize,
    /// Matches a rewrite may produce per iteration before it is banned.
    pub match_limit: usize,
    /// Iterations a rule is banned for the first time it exceeds the limit.
    pub ban_length: usize,
    /// Inner fixpoint rounds for the reduce tier per iteration.
    pub reduce_rounds: usize,
    /// Maximum number of cells in a tree window.
    pub window_cells: usize,
    /// Maximum rewrites applied per window optimisation.
    pub window_steps: usize,
    /// Beam width for window search; `0` disables the beam.
    pub beam_width: usize,
    /// Beam depth for window search.
    pub beam_depth: usize,
    /// Search steps a rewrite may spend looking for matches per
    /// iteration before it is banned like a rule that matched too often.
    pub match_fuel: usize,
    /// Stop after this many consecutive iterations that changed the graph
    /// without making the best term of any root cheaper; `0` disables the
    /// check.
    pub patience: usize,
    /// Wall-clock limit.
    pub time_limit: Option<Duration>,
}

impl Default for Budget {
    fn default() -> Self {
        Self {
            max_iterations: 24,
            max_nodes: 20_000,
            match_limit: 1_000,
            ban_length: 2,
            reduce_rounds: 32,
            window_cells: 256,
            window_steps: 512,
            beam_width: 0,
            beam_depth: 3,
            match_fuel: 200_000,
            patience: 3,
            time_limit: None,
        }
    }
}

/// Decides when a run has produced what was asked for.
pub trait Goal {
    /// Whether the run can stop now.
    fn reached(
        &self,
        graph: &Graph,
        roots: &[NodeId],
        env: &Env,
    ) -> bool;
}

/// Never satisfied early: run until saturation or budget exhaustion.
#[derive(Copy, Clone, Debug, Default)]
pub struct Saturate;

impl Goal for Saturate {
    fn reached(
        &self,
        _graph: &Graph,
        _roots: &[NodeId],
        _env: &Env,
    ) -> bool {
        false
    }
}

/// Satisfied once every root has a literal value or a numeric witness at
/// least as tight as the environment's tolerance.
#[derive(Copy, Clone, Debug, Default)]
pub struct Evaluated;

impl Goal for Evaluated {
    fn reached(
        &self,
        graph: &Graph,
        roots: &[NodeId],
        env: &Env,
    ) -> bool {
        roots.iter().all(|&r| {
            graph
                .approx(graph.find(r))
                .is_some_and(|b| b.rad <= env.tolerance.max(0.0))
        })
    }
}

/// Why a run ended.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Stop {
    /// Nothing changed during a whole iteration.
    Saturated,
    /// The goal was reached.
    Goal,
    /// Exploration kept changing the graph without improving the answer.
    Plateau,
    /// The iteration budget ran out.
    Iterations,
    /// The time limit was hit.
    Time,
}

/// Summary of a run.
#[derive(Clone, Debug)]
pub struct Report {
    /// Why the run ended.
    pub stop: Stop,
    /// Outer iterations performed.
    pub iterations: usize,
    /// Classes merged during the run, congruence included.
    pub unions: u64,
    /// How often each rule changed the graph, in program order; rules that
    /// never fired are omitted.
    pub applied: Vec<(Arc<str>, usize)>,
    /// Whether every root has a term free of heavy operators.
    pub solved: bool,
}

#[derive(Default)]
struct RuleState {
    banned_until: usize,
    times_banned: u32,
    visited: HashSet<NodeId>,
    applied: usize,
}

/// A program together with the scheduling logic that runs it.
#[derive(Clone, Debug, Default)]
pub struct Engine {
    program: Program,
}

impl Engine {
    /// Installs `sets` (and their dependencies) into `graph`.
    ///
    /// # Errors
    /// Propagates the first [`RuleError`].
    pub fn install(
        graph: &mut Graph,
        sets: &[RuleSet],
    ) -> Result<Self, RuleError> {
        let mut engine = Self::default();
        engine.extend(graph, sets)?;
        Ok(engine)
    }

    /// Installs further rule sets; sets already present are skipped.
    ///
    /// # Errors
    /// Propagates the first [`RuleError`].
    pub fn extend(
        &mut self,
        graph: &mut Graph,
        sets: &[RuleSet],
    ) -> Result<(), RuleError> {
        for set in sets {
            set.install_into(graph, &mut self.program)?;
        }
        Ok(())
    }

    /// The installed program.
    #[must_use]
    pub const fn program(&self) -> &Program {
        &self.program
    }

    /// Runs the program on the classes reachable from `roots`.
    pub fn run(
        &self,
        graph: &mut Graph,
        roots: &[NodeId],
        env: &Env,
        goal: &dyn Goal,
        budget: &Budget,
    ) -> Report {
        let start = Instant::now();
        let unions_at_start = graph.union_count();
        let mut states: Vec<RuleState> = self
            .program
            .rules
            .iter()
            .map(|_| RuleState::default())
            .collect();
        let mut windowed: HashMap<NodeId, u64> = HashMap::new();
        let normalizing: Vec<&Rewrite> = self.rewrites_of(Tier::Normalize);
        let exploring: Vec<&Rewrite> = self.rewrites_of(Tier::Explore);
        let passes: Vec<&dyn WindowPass> = self.program.passes.iter().map(|p| &**p).collect();

        graph.rebuild();
        // A nested run works inside its parent's numeric context and must
        // leave the parent's witnesses alone.
        if env.depth == 0 {
            graph.clear_witnesses();
            for &(symbol, value) in env.bindings() {
                let node = graph.symbol_node(symbol);
                graph.set_approx(graph.find(node), Ball::exact(value));
            }
        }

        let mut stop = Stop::Iterations;
        let mut iterations = 0;
        let mut best_cost = u64::MAX;
        let mut stalled = 0_usize;
        for iteration in 0..budget.max_iterations {
            iterations = iteration.saturating_add(1);
            let mut progress = false;

            // Tier 0: reduce heavy operators to a local fixpoint.
            for _ in 0..budget.reduce_rounds {
                let changed = self.apply_tier(
                    graph,
                    roots,
                    env,
                    Tier::Reduce,
                    iteration,
                    budget,
                    &mut states,
                );
                graph.rebuild();
                progress |= changed;
                if !changed {
                    break;
                }
            }

            // Tier 1: destructive clean-up in windows, then the same
            // simplifications modulo the graph's equalities.
            progress |= Self::windows(
                graph,
                roots,
                budget,
                &normalizing,
                &exploring,
                &passes,
                &mut windowed,
            );
            graph.rebuild();
            progress |= self.apply_tier(
                graph,
                roots,
                env,
                Tier::Normalize,
                iteration,
                budget,
                &mut states,
            );
            graph.rebuild();

            if env.numeric {
                progress |= propagate_witnesses(graph, roots);
            }
            if goal.reached(graph, roots, env) {
                stop = Stop::Goal;
                break;
            }

            // Tier 2: exploration, while there is room.
            if graph.len() < budget.max_nodes {
                progress |= self.apply_tier(
                    graph,
                    roots,
                    env,
                    Tier::Explore,
                    iteration,
                    budget,
                    &mut states,
                );
                graph.rebuild();
            }

            if !progress {
                // Banned rules may still have work to do; a run is only
                // saturated when none is waiting.
                if states.iter().any(|s| s.banned_until > iteration) {
                    continue;
                }
                stop = Stop::Saturated;
                break;
            }
            if budget
                .time_limit
                .is_some_and(|limit| start.elapsed() >= limit)
            {
                stop = Stop::Time;
                break;
            }
            if budget.patience > 0 {
                // Progress is measured on the closed form when there is one
                // and on the best partially reduced term otherwise.
                let total = |extractor: &Extractor, graph: &Graph| {
                    roots.iter().try_fold(0_u64, |acc, &r| extractor.cost(graph, r).map(|c| acc.saturating_add(c)))
                };
                let cost = total(&Extractor::new(graph, roots, &ClosedForm), graph)
                    .or_else(|| total(&Extractor::new(graph, roots, &SizeCost), graph).map(|c| c.saturating_add(1 << 40)))
                    .unwrap_or(u64::MAX);
                if cost < best_cost {
                    best_cost = cost;
                    stalled = 0;
                } else {
                    stalled = stalled.saturating_add(1);
                    if stalled >= budget.patience {
                        stop = Stop::Plateau;
                        break;
                    }
                }
            }
        }

        let closed = Extractor::new(graph, roots, &ClosedForm);
        Report {
            stop,
            iterations,
            unions: graph.union_count().saturating_sub(unions_at_start),
            applied: self
                .program
                .rules
                .iter()
                .zip(&states)
                .filter(|(_, s)| s.applied > 0)
                .map(|(r, s)| (Arc::clone(&r.name), s.applied))
                .collect(),
            solved: roots.iter().all(|&r| closed.cost(graph, r).is_some()),
        }
    }

    fn rewrites_of(
        &self,
        tier: Tier,
    ) -> Vec<&Rewrite> {
        self.program
            .rules
            .iter()
            .filter(|r| r.tier == tier)
            .filter_map(|r| match &r.action {
                | Action::Rewrite(rw) => Some(rw),
                | Action::Kernel(_) => None,
            })
            .collect()
    }

    /// Applies every unbanned rule of `tier` once. Returns whether anything
    /// changed.
    #[allow(clippy::too_many_arguments)]
    fn apply_tier(
        &self,
        graph: &mut Graph,
        roots: &[NodeId],
        env: &Env,
        tier: Tier,
        iteration: usize,
        budget: &Budget,
        states: &mut [RuleState],
    ) -> bool {
        if !self.program.rules.iter().any(|r| r.tier == tier) {
            return false;
        }
        let heavy =
            |graph: &Graph, node: NodeId| graph.ops().get(graph.op(node)).flags.has(OpFlags::HEAVY);
        let mut candidates: BTreeMap<OpId, Vec<NodeId>> = BTreeMap::new();
        for class in reachable(graph, roots) {
            // A request that already has an answer is not reduced again
            // through its other spellings.
            if tier == Tier::Reduce && graph.enodes(class).any(|n| !heavy(graph, n)) {
                continue;
            }
            for enode in graph.enodes(class) {
                candidates.entry(graph.op(enode)).or_default().push(enode);
            }
        }
        let mut changed = false;
        for (rule, state) in self.program.rules.iter().zip(states.iter_mut()) {
            if rule.tier != tier || state.banned_until > iteration {
                continue;
            }
            match &rule.action {
                | Action::Rewrite(rewrite) => {
                    let Some(nodes) = candidates.get(&rewrite.trigger()) else {
                        continue;
                    };
                    let shift = state.times_banned.min(16);
                    let limit = budget.match_limit.saturating_mul(1_usize << shift);
                    let mut matches = Vec::new();
                    let mut fuel = budget.match_fuel;
                    for &node in nodes {
                        let room = limit.saturating_add(1).saturating_sub(matches.len());
                        if room == 0 {
                            break;
                        }
                        matches.extend(rewrite.lhs.matches_bounded(graph, node, rewrite.nvars, room, &mut fuel));
                        if fuel == 0 {
                            break;
                        }
                    }
                    if matches.len() > limit || fuel == 0 {
                        state.banned_until = iteration
                            .saturating_add(budget.ban_length.saturating_mul(1_usize << shift))
                            .saturating_add(1);
                        state.times_banned = state.times_banned.saturating_add(1);
                        continue;
                    }
                    for m in &matches {
                        if let Some(new) = rewrite.build(graph, m) {
                            if graph.union(m.root, new) {
                                state.applied = state.applied.saturating_add(1);
                                changed = true;
                            }
                        }
                    }
                },
                | Action::Kernel(kernel) => {
                    let revisit = kernel.revisit();
                    let wanted = kernel.ops();
                    let nodes: Vec<NodeId> = if wanted.contains(&OpId::NONE) {
                        candidates.values().flatten().copied().collect()
                    } else {
                        wanted
                            .iter()
                            .filter_map(|op| candidates.get(op))
                            .flatten()
                            .copied()
                            .collect()
                    };
                    {
                        for node in nodes {
                            if !revisit && !state.visited.insert(node) {
                                continue;
                            }
                            let outcome = kernel.reduce(&mut Cx { graph, env, engine: self }, node);
                            let fired = match outcome {
                                | Outcome::Pass => false,
                                | Outcome::Equal(other) => graph.union(node, other),
                                | Outcome::Pinned(other) => {
                                    let merged = graph.union(node, other);
                                    graph.pin(other);
                                    merged
                                },
                                | Outcome::Approx(ball) => graph.set_approx(graph.find(node), ball),
                            };
                            if fired {
                                state.applied = state.applied.saturating_add(1);
                                changed = true;
                            }
                        }
                    }
                },
            }
        }
        changed
    }

    /// Optimises one window per root and per shared class, children first.
    #[allow(clippy::too_many_arguments)]
    fn windows(
        graph: &mut Graph,
        roots: &[NodeId],
        budget: &Budget,
        normalizing: &[&Rewrite],
        exploring: &[&Rewrite],
        passes: &[&dyn WindowPass],
        done: &mut HashMap<NodeId, u64>,
    ) -> bool {
        if normalizing.is_empty() && passes.is_empty() {
            return false;
        }
        let classes = reachable(graph, roots);
        let extractor = Extractor::new(graph, roots, &SizeCost);
        let root_classes: Vec<_> = roots.iter().map(|&r| graph.find(r)).collect();
        // Build every window's term before merging anything: the extractor
        // describes the graph as it is now.
        let mut terms = Vec::new();
        for &class in classes.iter().rev() {
            if !root_classes.contains(&class) && graph.parents(class).len() <= 1 {
                // Unshared: some ancestor's window covers it.
                continue;
            }
            let Some(term) = extractor.build(graph, NodeId::from_raw(class.raw())) else {
                continue;
            };
            // A term is worth another look once the graph has learnt
            // something since its last window: what its atoms are known to
            // equal may have changed.
            let stamp = graph.union_count();
            if !graph.children(term).is_empty() && done.insert(term, stamp) != Some(stamp) {
                terms.push((NodeId::from_raw(class.raw()), term));
            }
        }
        let mut changed = false;
        for (origin, term) in terms {
            // Building the term may have flattened it into a shape the
            // graph does not yet know to be in the class it came from.
            changed |= graph.union(origin, term);
            let mut window = TreeWindow::carve(graph, term, budget.window_cells);
            let mut steps = window.optimize(graph, normalizing, passes, budget.window_steps);
            if budget.beam_width > 0 && !exploring.is_empty() {
                let before = window.cost(graph);
                let best = window.beam(
                    graph,
                    exploring,
                    normalizing,
                    passes,
                    budget.beam_width,
                    budget.beam_depth,
                    budget.window_steps,
                );
                if best.cost(graph) < before {
                    window = best;
                    steps = steps.saturating_add(1);
                }
            }
            if steps > 0 {
                let out = window.commit(graph);
                done.insert(out, graph.union_count());
                changed |= graph.union(term, out);
            }
        }
        changed
    }
}

/// Computes numeric witnesses bottom-up for every reachable class whose
/// e-nodes have scalar semantics and whose children are already evaluated.
/// Returns whether any witness was added or tightened.
fn propagate_witnesses(
    graph: &mut Graph,
    roots: &[NodeId],
) -> bool {
    let classes = reachable(graph, roots);
    let mut changed = false;
    let mut args: Vec<f64> = Vec::new();
    // Reverse discovery order approximates children-first; iterate until
    // stable to cover the rest (cycles cannot produce new witnesses).
    loop {
        let mut round = false;
        for &class in classes.iter().rev() {
            if graph.approx(class).is_some() {
                continue;
            }
            let mut found = None;
            for enode in graph.enodes(class) {
                let Some(eval) = graph.ops().get(graph.op(enode)).eval else {
                    continue;
                };
                args.clear();
                let mut rad = 0.0_f64;
                let complete = graph.children(enode).iter().all(|&c| {
                    graph.approx(graph.find(c)).is_some_and(|b| {
                        args.push(b.mid);
                        rad = rad.max(b.rad);
                        true
                    })
                });
                if !complete {
                    continue;
                }
                let value = eval(&args);
                if value.is_finite() {
                    found = Some(Ball { mid: value, rad });
                    break;
                }
            }
            if let Some(ball) = found {
                round |= graph.set_approx(class, ball);
            }
        }
        changed |= round;
        if !round {
            return changed;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::number::Number;
    use crate::graph::op::Arity;
    use crate::graph::op::OpDescriptor;
    use crate::graph::op::OpFlags;
    use crate::graph::op::core;
    use crate::graph::rule::Kernel;

    /// Folds sums, products and powers of literal numbers.
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

    fn calculus() -> RuleSet {
        RuleSet::new("mini-calculus", |i| {
            i.op(OpDescriptor::new("sin", Arity::Fixed(1))
                .eval(|a| a.first().map_or(f64::NAN, |x| x.sin())))?;
            i.op(OpDescriptor::new("cos", Arity::Fixed(1))
                .eval(|a| a.first().map_or(f64::NAN, |x| x.cos())))?;
            i.op(OpDescriptor::new("diff", Arity::Fixed(2))
                .flags(OpFlags::HEAVY)
                .cost(50))?;
            i.kernel("fold", Tier::Normalize, Fold);
            i.rewrites(
                Tier::Reduce,
                &[
                    "d-const: diff(?c, ?x) => 0 if free_of(?c, ?x)",
                    "d-var: diff(?x, ?x) => 1",
                    "d-sum: diff(?a + ?b, ?x) => diff(?a, ?x) + diff(?b, ?x)",
                    "d-prod: diff(?a * ?b, ?x) => diff(?a, ?x) * ?b + ?a * diff(?b, ?x)",
                    "d-pow: diff(?a ^ ?n, ?x) => ?n * ?a^(?n - 1) * diff(?a, ?x) if free_of(?n, ?x)",
                    "d-sin: diff(sin(?a), ?x) => cos(?a) * diff(?a, ?x)",
                    "d-cos: diff(cos(?a), ?x) => -sin(?a) * diff(?a, ?x)",
                ],
            )?;
            i.rewrites(
                Tier::Normalize,
                &[
                    "add-0: ?a + 0 => ?a",
                    "mul-1: ?a * 1 => ?a",
                    "mul-0: ?a * 0 => 0",
                    "pow-1: ?a ^ 1 => ?a",
                    "pow-0: ?a ^ 0 => 1 if nonzero(?a)",
                    "pythagoras: sin(?x)^2 + cos(?x)^2 => 1",
                ],
            )?;
            i.rewrites(
                Tier::Explore,
                &["distribute: ?a * (?b + ?c) <=> ?a*?b + ?a*?c"],
            )
        })
    }

    fn run(
        src: &str,
        env: &Env,
        goal: &dyn Goal,
    ) -> (Graph, NodeId, Report) {
        let mut g = Graph::new();
        let engine = Engine::install(&mut g, &[calculus()]).unwrap_or_else(|e| panic!("{e}"));
        let root = g.parse(src).unwrap_or_else(|e| panic!("{e}"));
        let report = engine.run(&mut g, &[root], env, goal, &Budget::default());
        assert_eq!(g.validate(), Ok(()));
        assert!(g.conflicts().is_empty(), "unsound: {:?}", g.conflicts());
        (g, root, report)
    }

    fn best(
        g: &mut Graph,
        root: NodeId,
    ) -> String {
        let ex = Extractor::new(g, &[root], &ClosedForm);
        ex.build(g, root)
            .map_or_else(|| "<unreduced>".to_owned(), |n| g.display(n))
    }

    #[test]
    fn derivative_is_reduced_and_cleaned_up() {
        let (mut g, root, report) = run("diff(x^2, x)", &Env::symbolic(), &Saturate);
        assert!(report.solved);
        assert_eq!(best(&mut g, root), "2*x");
    }

    #[test]
    fn chain_and_product_rules() {
        let (mut g, root, report) = run("diff(x * sin(x^2), x)", &Env::symbolic(), &Saturate);
        assert!(report.solved);
        let text = best(&mut g, root);
        // d/dx x·sin(x²) = sin(x²) + 2x²·cos(x²); check by evaluation rather
        // than by shape.
        let node = g.parse(&text).unwrap_or(NodeId::NONE);
        let x = g
            .interner()
            .find_symbol("x")
            .unwrap_or(crate::graph::SymbolId::NONE);
        let engine = Engine::install(&mut g, &[calculus()]).unwrap_or_else(|e| panic!("{e}"));
        let mut env = Env::numeric(1e-9);
        env.bind(x, 0.7);
        engine.run(&mut g, &[node], &env, &Evaluated, &Budget::default());
        let got = g.approx(g.find(node)).map_or(f64::NAN, |b| b.mid);
        let want = (0.49_f64).sin() + 2.0 * 0.49 * (0.49_f64).cos();
        assert!(
            (got - want).abs() < 1e-12,
            "{text} evaluated to {got}, expected {want}"
        );
    }

    #[test]
    fn unknown_dependencies_stay_unreduced() {
        // Nothing says how `apply(f, x)` depends on x: no rule may guess.
        let (mut g, root, report) = run("diff(f(x), x)", &Env::symbolic(), &Saturate);
        assert!(!report.solved);
        assert_eq!(best(&mut g, root), "<unreduced>");
    }

    #[test]
    fn numeric_goal_stops_early_with_a_witness() {
        let mut g = Graph::new();
        let engine = Engine::install(&mut g, &[calculus()]).unwrap_or_else(|e| panic!("{e}"));
        let root = g.parse("diff(sin(x)*cos(x), x)").unwrap_or(NodeId::NONE);
        let x = g
            .interner()
            .find_symbol("x")
            .unwrap_or(crate::graph::SymbolId::NONE);
        let mut env = Env::numeric(1e-9);
        env.bind(x, 0.0);
        let report = engine.run(&mut g, &[root], &env, &Evaluated, &Budget::default());
        assert_eq!(report.stop, Stop::Goal);
        let value = g.approx(g.find(root)).map(|b| b.mid);
        assert_eq!(value, Some(1.0), "cos(0)^2 - sin(0)^2");
        // The binding is an environment, not an equality.
        let zero = g.int(0);
        let xn = g.sym("x");
        assert!(!g.same(xn, zero));
    }

    #[test]
    fn witnesses_do_not_leak_between_runs() {
        let mut g = Graph::new();
        let engine = Engine::install(&mut g, &[calculus()]).unwrap_or_else(|e| panic!("{e}"));
        let root = g.parse("sin(x) + 1").unwrap_or(NodeId::NONE);
        let x = g
            .interner()
            .find_symbol("x")
            .unwrap_or(crate::graph::SymbolId::NONE);
        for value in [0.0, 1.0] {
            let mut env = Env::numeric(1e-9);
            env.bind(x, value);
            engine.run(&mut g, &[root], &env, &Evaluated, &Budget::default());
            let got = g.approx(g.find(root)).map_or(f64::NAN, |b| b.mid);
            assert!((got - (value.sin() + 1.0)).abs() < 1e-15);
        }
        engine.run(
            &mut g,
            &[root],
            &Env::symbolic(),
            &Saturate,
            &Budget::default(),
        );
        assert_eq!(g.approx(g.find(root)), None);
    }

    #[test]
    fn windows_do_the_cleanup_without_bloating_the_graph() {
        let mut g = Graph::new();
        let engine = Engine::install(&mut g, &[calculus()]).unwrap_or_else(|e| panic!("{e}"));
        let root = g
            .parse("(a + 0) * 1 * (sin(t)^2 + cos(t)^2 + b*0)")
            .unwrap_or(NodeId::NONE);
        let before = g.len();
        let budget = Budget {
            max_nodes: 0,
            ..Budget::default()
        };
        let report = engine.run(&mut g, &[root], &Env::symbolic(), &Saturate, &budget);
        assert_eq!(best(&mut g, root), "a");
        assert!(
            g.len() - before <= 12,
            "cleanup created {} nodes",
            g.len() - before
        );
        assert_eq!(report.stop, Stop::Saturated);
    }

    #[test]
    fn explosive_rules_are_banned_but_the_run_still_terminates() {
        let mut g = Graph::new();
        let engine = Engine::install(&mut g, &[calculus()]).unwrap_or_else(|e| panic!("{e}"));
        let root = g
            .parse("(a + b) * (c + d) * (e + f) * (g + h)")
            .unwrap_or(NodeId::NONE);
        let budget = Budget {
            match_limit: 8,
            max_iterations: 12,
            ..Budget::default()
        };
        let report = engine.run(&mut g, &[root], &Env::symbolic(), &Saturate, &budget);
        assert!(report.iterations <= 12);
        assert_eq!(g.validate(), Ok(()));
        assert!(g.conflicts().is_empty());
    }

    #[test]
    fn unrelated_terms_are_left_alone() {
        let mut g = Graph::new();
        let engine = Engine::install(&mut g, &[calculus()]).unwrap_or_else(|e| panic!("{e}"));
        let bystander = g.parse("diff(y^3, y)").unwrap_or(NodeId::NONE);
        let root = g.parse("diff(x^2, x)").unwrap_or(NodeId::NONE);
        engine.run(
            &mut g,
            &[root],
            &Env::symbolic(),
            &Saturate,
            &Budget::default(),
        );
        assert_eq!(
            g.class_size(g.find(bystander)),
            1,
            "only reachable classes are worked on"
        );
    }

    #[test]
    fn time_limit_is_honoured() {
        let mut g = Graph::new();
        let engine = Engine::install(&mut g, &[calculus()]).unwrap_or_else(|e| panic!("{e}"));
        let root = g
            .parse("(a + b) * (c + d) * (e + f)")
            .unwrap_or(NodeId::NONE);
        let budget = Budget {
            time_limit: Some(Duration::ZERO),
            ..Budget::default()
        };
        let report = engine.run(&mut g, &[root], &Env::symbolic(), &Saturate, &budget);
        assert!(matches!(report.stop, Stop::Time | Stop::Saturated));
        assert!(report.iterations <= 2);
    }

    #[test]
    fn plateau_stops_fruitless_exploration() {
        let mut g = Graph::new();
        let engine = Engine::install(&mut g, &[calculus()]).unwrap_or_else(|e| panic!("{e}"));
        let root = g
            .parse("(a + b) * (c + d) * (e + f)")
            .unwrap_or(NodeId::NONE);
        let budget = Budget {
            max_iterations: 50,
            patience: 2,
            ..Budget::default()
        };
        let report = engine.run(&mut g, &[root], &Env::symbolic(), &Saturate, &budget);
        assert_eq!(report.stop, Stop::Plateau);
        assert!(
            report.iterations <= 4,
            "ran {} iterations",
            report.iterations
        );
        assert_eq!(best(&mut g, root), "(a + b)*(c + d)*(e + f)");
    }
}
