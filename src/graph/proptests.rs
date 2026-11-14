//! Property tests for the kernel invariants.

use proptest::prelude::*;

use super::extract::Extractor;
use super::extract::SizeCost;
use super::id::NodeId;
use super::id::OpId;
use super::op::Arity;
use super::op::OpDescriptor;
use super::op::OpFlags;
use super::op::core;
use super::rule::Env;
use super::rule::Rewrite;
use super::rule::RuleSet;
use super::rule::Tier;
use super::schedule::Budget;
use super::schedule::Engine;
use super::schedule::Saturate;
use super::store::Graph;
use super::window::TreeWindow;

/// One step of a random term-building program: an opcode and two operand
/// picks (indices into the pool of nodes built so far).
type Step = (u8, usize, usize);

fn steps(max: usize) -> impl Strategy<Value = Vec<Step>> {
    prop::collection::vec((0_u8..6, 0_usize..1000, 0_usize..1000), 1..max)
}

/// Builds a pool of nodes over leaves `a..d`, a unary `f`, a binary `g`, a
/// commutative binary `h` and the AC `add`.
fn build_uninterpreted(
    g: &mut Graph,
    program: &[Step],
) -> Vec<NodeId> {
    let f = g
        .ops_mut()
        .register(OpDescriptor::new("f", Arity::Fixed(1)))
        .unwrap_or(OpId::NONE);
    let gg = g
        .ops_mut()
        .register(OpDescriptor::new("g", Arity::Fixed(2)))
        .unwrap_or(OpId::NONE);
    let h = g
        .ops_mut()
        .register(OpDescriptor::new("h", Arity::Fixed(2)).flags(OpFlags::COMMUTATIVE))
        .unwrap_or(OpId::NONE);
    let mut pool: Vec<NodeId> = ["a", "b", "c", "d"].iter().map(|n| g.sym(n)).collect();
    for &(op, i, j) in program {
        let x = pool[i % pool.len()];
        let y = pool[j % pool.len()];
        let node = match op {
            | 0 | 1 => g.node(f, &[x]),
            | 2 => g.node(gg, &[x, y]),
            | 3 => g.node(h, &[x, y]),
            | _ => g.node(core::ADD, &[x, y]),
        };
        pool.push(node);
    }
    pool
}

/// Textbook congruence closure over the concrete nodes of `g`, quadratic
/// and obviously correct.
fn naive_closure(
    g: &Graph,
    unions: &[(NodeId, NodeId)],
) -> Vec<usize> {
    fn find(
        parent: &mut [usize],
        mut x: usize,
    ) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }
    let n = g.len();
    let mut parent: Vec<usize> = (0..n).collect();
    for &(a, b) in unions {
        let (ra, rb) = (find(&mut parent, a.index()), find(&mut parent, b.index()));
        parent[ra] = rb;
    }
    loop {
        let mut changed = false;
        let keys: Vec<(OpId, Vec<usize>)> = (0..n)
            .map(|i| {
                let id = NodeId::from_raw(u32::try_from(i).unwrap_or(u32::MAX));
                let mut key: Vec<usize> = g
                    .children(id)
                    .iter()
                    .map(|c| find(&mut parent, c.index()))
                    .collect();
                if g.ops().get(g.op(id)).flags.has(OpFlags::COMMUTATIVE) {
                    key.sort_unstable();
                }
                (g.op(id), key)
            })
            .collect();
        for i in 0..n {
            for j in (i + 1)..n {
                let leaf = g
                    .children(NodeId::from_raw(u32::try_from(i).unwrap_or(u32::MAX)))
                    .is_empty();
                if !leaf && keys[i] == keys[j] {
                    let (ri, rj) = (find(&mut parent, i), find(&mut parent, j));
                    if ri != rj {
                        parent[ri] = rj;
                        changed = true;
                    }
                }
            }
        }
        if !changed {
            return (0..n).map(|i| find(&mut parent, i)).collect();
        }
    }
}

fn arithmetic_rules() -> RuleSet {
    RuleSet::new("prop-arith", |i| {
        i.rewrites(
            Tier::Normalize,
            &[
                "add-0: ?a + 0 => ?a",
                "mul-1: ?a * 1 => ?a",
                "mul-0: ?a * 0 => 0",
                "pow-1: ?a ^ 1 => ?a",
                "double: ?a + ?a => 2 * ?a",
                "square: ?a * ?a => ?a ^ 2",
            ],
        )?;
        i.rewrites(
            Tier::Explore,
            &[
                "factor: ?a*?b + ?a*?c => ?a*(?b + ?c)",
                "expand: ?a*(?b + ?c) => ?a*?b + ?a*?c",
            ],
        )
    })
}

/// Builds a pool of arithmetic terms over `x, y` and small integers.
fn build_arithmetic(
    g: &mut Graph,
    program: &[Step],
) -> Vec<NodeId> {
    let mut pool = vec![
        g.sym("x"),
        g.sym("y"),
        g.int(0),
        g.int(1),
        g.int(2),
        g.int(-1),
    ];
    for &(op, i, j) in program {
        let x = pool[i % pool.len()];
        let y = pool[j % pool.len()];
        let node = match op {
            | 0 | 1 => g.node(core::ADD, &[x, y]),
            | 2 | 3 => g.node(core::MUL, &[x, y]),
            | 4 => {
                let two = g.int(2);
                g.node(core::POW, &[x, two])
            },
            | _ => {
                let one = g.int(1);
                g.node(core::POW, &[x, one])
            },
        };
        pool.push(node);
    }
    pool
}

fn env_at(
    g: &mut Graph,
    x: f64,
    y: f64,
) -> Env {
    let mut env = Env::numeric(0.0);
    env.bind(g.interner_mut().symbol("x"), x);
    env.bind(g.interner_mut().symbol("y"), y);
    env
}

fn close(
    a: f64,
    b: f64,
) -> bool {
    (a - b).abs() <= 1e-6 * (1.0 + a.abs().max(b.abs())) || (!a.is_finite() && !b.is_finite())
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(200))]

    #[test]
    fn hash_consing_is_deterministic(program in steps(40)) {
        let mut g = Graph::new();
        let first = build_uninterpreted(&mut g, &program);
        let size = g.len();
        let second = build_uninterpreted(&mut g, &program);
        prop_assert_eq!(first, second);
        prop_assert_eq!(g.len(), size, "rebuilding the same terms must not allocate");
        prop_assert_eq!(g.validate(), Ok(()));
    }

    #[test]
    fn congruence_closure_agrees_with_the_naive_algorithm(
        program in steps(30),
        merges in prop::collection::vec((0_usize..1000, 0_usize..1000), 0..8),
    ) {
        let mut g = Graph::new();
        let pool = build_uninterpreted(&mut g, &program);
        let unions: Vec<(NodeId, NodeId)> =
            merges.iter().map(|&(i, j)| (pool[i % pool.len()], pool[j % pool.len()])).collect();
        for &(a, b) in &unions {
            g.union(a, b);
        }
        g.rebuild();
        prop_assert_eq!(g.validate(), Ok(()));
        let naive = naive_closure(&g, &unions);
        for i in 0..g.len() {
            for j in (i + 1)..g.len() {
                let (a, b) = (
                    NodeId::from_raw(u32::try_from(i).unwrap_or(u32::MAX)),
                    NodeId::from_raw(u32::try_from(j).unwrap_or(u32::MAX)),
                );
                prop_assert_eq!(g.same(a, b), naive[i] == naive[j], "nodes {} and {}", g.display(a), g.display(b));
            }
        }
    }

    #[test]
    fn rings_partition_the_nodes(
        program in steps(30),
        merges in prop::collection::vec((0_usize..1000, 0_usize..1000), 0..8),
    ) {
        let mut g = Graph::new();
        let pool = build_uninterpreted(&mut g, &program);
        for &(i, j) in &merges {
            g.union(pool[i % pool.len()], pool[j % pool.len()]);
        }
        g.rebuild();
        let mut seen = vec![0_u32; g.len()];
        let mut classes = 0;
        for i in 0..g.len() {
            let id = NodeId::from_raw(u32::try_from(i).unwrap_or(u32::MAX));
            if g.find(id).raw() == id.raw() {
                classes += 1;
                for member in g.members(g.find(id)) {
                    seen[member.index()] += 1;
                }
            }
        }
        prop_assert!(seen.iter().all(|&c| c == 1), "every node is in exactly one ring");
        prop_assert_eq!(classes, g.class_count());
    }

    #[test]
    fn window_commit_without_edits_is_the_identity(program in steps(30), cells in 1_usize..40) {
        let mut g = Graph::new();
        let pool = build_arithmetic(&mut g, &program);
        let root = *pool.last().unwrap_or(&NodeId::NONE);
        let window = TreeWindow::carve(&g, root, cells);
        prop_assert_eq!(window.commit(&mut g), root);
    }

    #[test]
    fn window_optimisation_preserves_value(program in steps(25), x in -3.0_f64..3.0, y in -3.0_f64..3.0) {
        let mut g = Graph::new();
        let engine = Engine::install(&mut g, &[arithmetic_rules()]).unwrap_or_else(|e| panic!("{e}"));
        let pool = build_arithmetic(&mut g, &program);
        let root = *pool.last().unwrap_or(&NodeId::NONE);
        let rewrites: Vec<&Rewrite> = engine
            .program()
            .rules
            .iter()
            .filter(|r| r.tier == Tier::Normalize)
            .filter_map(|r| match &r.action {
                | super::rule::Action::Rewrite(rw) => Some(rw),
                | super::rule::Action::Kernel(_) => None,
            })
            .collect();
        let mut window = TreeWindow::carve(&g, root, 128);
        window.optimize(&mut g, &rewrites, &[], 200);
        prop_assert!(window.len() == window.post_order().len());
        let out = window.commit(&mut g);
        let env = env_at(&mut g, x, y);
        let (want, got) = (g.eval(root, &env), g.eval(out, &env));
        prop_assert!(want.is_some() && got.is_some());
        prop_assert!(
            close(want.unwrap_or(f64::NAN), got.unwrap_or(f64::NAN)),
            "{} = {:?} but {} = {:?}", g.display(root), want, g.display(out), got
        );
    }

    #[test]
    fn engine_is_sound_and_extraction_never_gets_worse(
        program in steps(14), x in -2.0_f64..2.0, y in -2.0_f64..2.0,
    ) {
        let mut g = Graph::new();
        let engine = Engine::install(&mut g, &[arithmetic_rules()]).unwrap_or_else(|e| panic!("{e}"));
        let pool = build_arithmetic(&mut g, &program);
        let root = *pool.last().unwrap_or(&NodeId::NONE);
        let before = Extractor::new(&g, &[root], &SizeCost).cost(&g, root);
        let budget = Budget { max_iterations: 4, max_nodes: 2_000, ..Budget::default() };
        engine.run(&mut g, &[root], &Env::symbolic(), &Saturate, &budget);
        prop_assert_eq!(g.validate(), Ok(()));
        prop_assert!(g.conflicts().is_empty(), "two different numbers were merged: {:?}", g.conflicts());
        let extractor = Extractor::new(&g, &[root], &SizeCost);
        prop_assert!(extractor.cost(&g, root) <= before);
        let out = extractor.build(&mut g, root).unwrap_or(NodeId::NONE);
        let env = env_at(&mut g, x, y);
        let (want, got) = (g.eval(root, &env).unwrap_or(f64::NAN), g.eval(out, &env).unwrap_or(f64::NAN));
        prop_assert!(close(want, got), "{} = {} but {} = {}", g.display(root), want, g.display(out), got);
        // Every member of the root class must agree, not just the cheapest.
        for member in g.members(g.find(root)) {
            if let Some(value) = g.eval(member, &env) {
                prop_assert!(close(want, value), "{} = {} but member {} = {}", g.display(root), want, g.display(member), value);
            }
        }
    }
}
