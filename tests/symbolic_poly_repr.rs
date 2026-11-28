//! Shared opaque subterms and sparse polynomial boundary contracts.

use rssn::graph::Arity;
use rssn::graph::Graph;
use rssn::graph::Number;
use rssn::graph::OpDescriptor;
use rssn::graph::op::core;
use rssn::rules::poly::repr::Gens;
use rssn::rules::poly::repr::Limits;
use rssn::rules::poly::repr::Poly;
use rssn::rules::poly::repr::from_term;
use rssn::rules::poly::repr::to_term;

#[test]
fn shared_opaque_subterms_keep_their_expanded_graph_shape_across_callers() {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(4)
        .build()
        .unwrap();
    pool.install(|| {
        rssn::kernels::par::map_indexed(8, 2, |index| {
            let mut graph = Graph::new();
            let pair = graph
                .ops_mut()
                .register(OpDescriptor::new("opaque_pair", Arity::Fixed(2)))
                .unwrap();
            let mut raw = graph.parse("(x+1)^2").unwrap();
            let mut expected = graph.parse("x^2+2*x+1").unwrap();
            for _ in 0..30 + index {
                raw = graph.node(pair, &[raw, raw]);
                expected = graph.node(pair, &[expected, expected]);
            }
            let mut gens = Gens::default();
            let polynomial = from_term(&mut graph, &mut gens, raw, Limits::default()).unwrap();
            assert_eq!(polynomial, Poly::generator(0));
            let result = to_term(&mut graph, &gens, &polynomial);
            assert!(graph.same(result, expected));
        });
    });
}

#[test]
fn opaque_children_do_not_reset_the_expansion_depth_limit() {
    std::thread::Builder::new()
        .stack_size(2 * 1024 * 1024)
        .spawn(|| {
            let mut graph = Graph::new();
            let unary = graph
                .ops_mut()
                .register(OpDescriptor::new("opaque_unary", Arity::Fixed(1)))
                .unwrap();
            let one = graph.int(1);
            for mixed in [false, true] {
                let mut node = graph.sym("x");
                for index in 0..2100 {
                    node = if mixed && index % 2 == 1 {
                        graph.node(core::ADD, &[node, one])
                    } else {
                        graph.node(unary, &[node])
                    };
                    if index == 1999 {
                        assert!(
                            from_term(&mut graph, &mut Gens::default(), node, Limits::default())
                                .is_some()
                        );
                    }
                }
                assert!(
                    from_term(&mut graph, &mut Gens::default(), node, Limits::default()).is_none()
                );
            }
        })
        .unwrap()
        .join()
        .unwrap();
}

#[test]
fn preseeded_equal_generators_keep_the_earliest_index() {
    let mut graph = Graph::new();
    let mut gens = Gens::default();
    let nodes: Vec<_> = (0..24)
        .map(|index| {
            let node = graph.sym(&format!("g{index}"));
            assert_eq!(gens.index(&graph, node), index);
            node
        })
        .collect();
    assert!(graph.union(nodes[19], nodes[3]));
    let source = graph.node(core::ADD, &[nodes[3], nodes[19], nodes[22]]);
    let polynomial = from_term(&mut graph, &mut gens, source, Limits::default()).unwrap();
    assert_eq!(
        polynomial,
        Poly::generator(3)
            .scale(&Number::from(2))
            .add(&Poly::generator(22))
    );
    assert_eq!(gens.len(), 24);
    assert_eq!(gens.node(3), Some(nodes[3]));
}

#[test]
fn polynomial_views_keep_zero_and_exponent_boundaries_exact() {
    let degree_limit = Poly::monomial(vec![(2, u32::MAX)], Number::from(1));
    assert!(degree_limit.mul(&Poly::generator(2), 10).is_none());
    let mixed = degree_limit.mul(&Poly::generator(3), 10).unwrap();
    assert_eq!(mixed.total_degree(), u32::MAX);
    assert!(mixed.univariate_in(2).is_none());
    assert!(mixed.univariate_in(3).is_none());

    let coefficients = [Number::from(-3), Number::from(0), Number::from(2)];
    assert_eq!(
        Poly::from_univariate(7, &coefficients).univariate_in(7),
        Some(coefficients.to_vec())
    );
    assert_eq!(Poly::zero().univariate_in(7), Some(vec![Number::from(0)]));

    let tiny = Poly::monomial(vec![(3, 2)], Number::Float(f64::from_bits(1)));
    assert!(tiny.scale(&Number::Float(0.5)).is_zero());
}
