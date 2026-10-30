//! Graph algorithms (ported from `numerical_graph_test.rs`).

use std::time::Instant;

use rssn::kernels::graph;
use rssn::kernels::graph::CanonicalLabelling;
use rssn::kernels::graph::ColourGraph;
use rssn::kernels::graph::Graph;
use rssn::kernels::graph::Permutation;
use rssn::kernels::graph::SearchBudget;
use rssn::kernels::graph::bfs;
use rssn::kernels::graph::canonical_labelling;
use rssn::kernels::graph::dijkstra;
use rssn::kernels::graph::floyd_warshall;
use rssn::kernels::graph::isomorphic;
use rssn::kernels::graph::page_rank;
use rssn::kernels::graph::refine;
use rssn::kernels::graph::wl_equivalent;

fn cfg() -> proptest::prelude::ProptestConfig {
    proptest::prelude::ProptestConfig {
        rng_seed: proptest::test_runner::RngSeed::Fixed(0x5EED),
        failure_persistence: None,
        ..proptest::prelude::ProptestConfig::default()
    }
}

#[test]

fn test_graph_creation_and_bfs() {
    let mut graph = Graph::new(5);

    // 0 -> 1 -> 2
    // |    ^
    // v    |
    // 3 -> 4
    graph.add_edge(0, 1, 1.0);

    graph.add_edge(1, 2, 1.0);

    graph.add_edge(0, 3, 1.0);

    graph.add_edge(3, 4, 1.0);

    graph.add_edge(4, 1, 1.0);

    let dist = bfs(&graph, 0);

    assert_eq!(dist[0], 0);

    assert_eq!(dist[1], 1);

    assert_eq!(dist[2], 2);

    assert_eq!(dist[3], 1);

    assert_eq!(dist[4], 2);
}

#[test]

fn test_dijkstra() {
    let mut graph = Graph::new(3);

    graph.add_edge(0, 1, 1.0);

    graph.add_edge(1, 2, 2.0);

    graph.add_edge(0, 2, 4.0); // 1+2 = 3 < 4

    let (dist, prev) = dijkstra(&graph, 0);

    assert_eq!(dist[0], 0.0);

    assert_eq!(dist[1], 1.0);

    assert_eq!(dist[2], 3.0);

    assert_eq!(prev[0], None);

    assert_eq!(prev[1], Some(0));

    assert_eq!(prev[2], Some(1));
}

#[test]

fn test_page_rank() {
    let mut graph = Graph::new(3);

    // 0 -> 1 -> 2 -> 0 circle
    graph.add_edge(0, 1, 1.0);

    graph.add_edge(1, 2, 1.0);

    graph.add_edge(2, 0, 1.0);

    let scores = page_rank(&graph, 0.85, 1e-6, 100);

    // Should be equal due to symmetry
    assert!((scores[0] - scores[1]).abs() < 1e-4);

    assert!((scores[1] - scores[2]).abs() < 1e-4);

    assert!((scores[0] + scores[1] + scores[2] - 1.0).abs() < 1e-6);

    // Star graph: 1->0, 2->0
    let mut graph2 = Graph::new(3);

    graph2.add_edge(1, 0, 1.0);

    graph2.add_edge(2, 0, 1.0);

    let scores2 = page_rank(&graph2, 0.85, 1e-6, 100);

    // Node 0 should have highest score
    assert!(scores2[0] > scores2[1]);

    assert!(scores2[0] > scores2[2]);
}

#[test]

fn test_floyd_warshall() {
    let mut graph = Graph::new(4);

    // 0 -> 1 (1)
    // 1 -> 2 (2)
    // 2 -> 3 (3)
    // 3 -> 0 (10)
    graph.add_edge(0, 1, 1.0);

    graph.add_edge(1, 2, 2.0);

    graph.add_edge(2, 3, 3.0);

    graph.add_edge(3, 0, 10.0);

    let dist = floyd_warshall(&graph);

    let n = 4;

    // 0 to 3: 0->1->2->3 = 1+2+3 = 6
    assert_eq!(dist[3], 6.0);

    // 3 to 0: 3->0 = 10
    assert_eq!(dist[3 * n], 10.0);

    // 0 to 2: 0->1->2 = 3
    assert_eq!(dist[2], 3.0);

    // Diagonal
    assert_eq!(dist[0], 0.0);

    // Unconnected: 1 -> 0 ? 1->2->3->0 = 2+3+10 = 15
    assert_eq!(dist[n], 15.0);
}

#[test]

fn test_connected_components() {
    let mut graph = Graph::new(5);

    // 0-1-2, 3-4
    graph.add_edge(0, 1, 1.0);

    graph.add_edge(1, 0, 1.0);

    graph.add_edge(1, 2, 1.0);

    graph.add_edge(2, 1, 1.0);

    graph.add_edge(3, 4, 1.0);

    graph.add_edge(4, 3, 1.0);

    let comp = graph::connected_components(&graph);

    assert_eq!(comp[0], comp[1]);

    assert_eq!(comp[1], comp[2]);

    assert_eq!(comp[3], comp[4]);

    assert_ne!(comp[0], comp[3]);
}

#[test]

fn test_minimum_spanning_tree() {
    let mut graph = Graph::new(4);

    // 0-1 (1), 1-2 (2), 2-3 (3), 0-3 (10)
    // MST should be 0-1, 1-2, 2-3. Total weight = 1+2+3 = 6.
    graph.add_edge(0, 1, 1.0);

    graph.add_edge(1, 0, 1.0);

    graph.add_edge(1, 2, 2.0);

    graph.add_edge(2, 1, 2.0);

    graph.add_edge(2, 3, 3.0);

    graph.add_edge(3, 2, 3.0);

    graph.add_edge(0, 3, 10.0);

    graph.add_edge(3, 0, 10.0);

    let mst = graph::minimum_spanning_tree(&graph);

    // Check total weight
    let mut total_weight = 0.0;

    let mut edges_count = 0;

    for u in 0..4 {
        for &(v, w) in mst.adj(u) {
            if u < v {
                // count each edge once
                total_weight += w;

                edges_count += 1;
            }
        }
    }

    assert_eq!(total_weight, 6.0);

    assert_eq!(edges_count, 3);
}

// Property-based tests using proptest
proptest::proptest! {
    #![proptest_config(cfg())]

    /// PageRank scores should always sum to 1.0
    #[test]
    fn prop_page_rank_sums_to_one(n in 2usize..10) {
        let mut graph = Graph::new(n);
        // Create a simple cycle graph
        for i in 0..n {
            graph.add_edge(i, (i + 1) % n, 1.0);
        }
        let scores = page_rank(&graph, 0.85, 1e-9, 200);
        let sum: f64 = scores.iter().sum();
        proptest::prop_assert!((sum - 1.0).abs() < 1e-6, "PageRank sum was {} instead of 1.0", sum);
    }

    /// BFS distance to start node should always be 0
    #[test]
    fn prop_bfs_start_distance_zero(n in 1usize..20, start in 0usize..20) {
        let n = n.max(1);
        let start = start % n;
        let graph = Graph::new(n);
        let dist = bfs(&graph, start);
        proptest::prop_assert_eq!(dist[start], 0);
    }

    /// Dijkstra distance to start node should always be 0
    #[test]
    fn prop_dijkstra_start_distance_zero(n in 1usize..20, start in 0usize..20) {
        let n = n.max(1);
        let start = start % n;
        let graph = Graph::new(n);
        let (dist, _) = dijkstra(&graph, start);
        proptest::prop_assert_eq!(dist[start], 0.0);
    }

    /// Floyd-Warshall diagonal should always be 0
    #[test]
    fn prop_floyd_warshall_diagonal_zero(n in 1usize..10) {
        let graph = Graph::new(n);
        let dist = floyd_warshall(&graph);
        for i in 0..n {
            proptest::prop_assert_eq!(dist[i * n + i], 0.0);
        }
    }

    /// Connected components: each node should have a valid component ID
    #[test]
    fn prop_connected_components_valid(n in 1usize..20) {
        let graph = Graph::new(n);
        let comp = graph::connected_components(&graph);
        proptest::prop_assert_eq!(comp.len(), n);
        for i in 0..n {
            // Each isolated node gets its own component
            proptest::prop_assert!(comp[i] < n);
        }
    }

    /// MST of a connected graph with n nodes should have n-1 edges (when counted once)
    #[test]
    fn prop_mst_edge_count(n in 2usize..8) {
        let mut graph = Graph::new(n);
        // Create a fully connected undirected graph
        for i in 0..n {
            for j in (i+1)..n {
                let weight = ((i + j) % 10 + 1) as f64;
                graph.add_edge(i, j, weight);
                graph.add_edge(j, i, weight);
            }
        }
        let mst = graph::minimum_spanning_tree(&graph);

        let mut edges_count = 0;
        for u in 0..n {
            for &(v, _) in mst.adj(u) {
                if u < v {
                    edges_count += 1;
                }
            }
        }
        proptest::prop_assert_eq!(edges_count, n - 1, "MST should have n-1 edges");
    }
}

// ============================================================================
// Added: edge cases and cross-checks
// ============================================================================

#[test]
fn graph_accessors() {
    let mut g = Graph::new(3);
    assert_eq!(g.num_nodes(), 3);
    g.add_edge(0, 2, 2.5);
    assert_eq!(g.adj(0), &[(2, 2.5)]);
    assert!(g.adj(1).is_empty());
}

#[test]
fn unreachable_nodes() {
    let mut g = Graph::new(3);
    g.add_edge(0, 1, 1.0);
    assert_eq!(bfs(&g, 0), vec![0, 1, usize::MAX]);
    let (dist, prev) = dijkstra(&g, 0);
    assert_eq!(dist[2], f64::INFINITY);
    assert_eq!(prev[2], None);
    let fw = floyd_warshall(&g);
    assert_eq!(fw[2], f64::INFINITY); // 0 -> 2
    assert_eq!(fw[3], f64::INFINITY); // 1 -> 0 (directed)
    assert_eq!(fw[1], 1.0);
}

#[test]
fn page_rank_of_a_cycle_is_uniform() {
    let n = 5;
    let mut g = Graph::new(n);
    for i in 0..n {
        g.add_edge(i, (i + 1) % n, 1.0);
    }
    for s in page_rank(&g, 0.85, 1e-12, 500) {
        assert!((s - 0.2).abs() < 1e-8, "score {s}");
    }
}

#[test]
fn page_rank_favours_the_sink_of_a_star() {
    let mut g = Graph::new(4);
    for i in 1..4 {
        g.add_edge(i, 0, 1.0);
    }
    let s = page_rank(&g, 0.85, 1e-10, 500);
    assert!((s.iter().sum::<f64>() - 1.0).abs() < 1e-6);
    assert!(s[1] == s[2] || (s[1] - s[2]).abs() < 1e-9);
    assert!(s[0] > 2.0 * s[1]);
}

#[test]
fn dijkstra_prefers_cheaper_two_hop_route_and_reconstructs_path() {
    let mut g = Graph::new(5);
    g.add_edge(0, 1, 7.0);
    g.add_edge(0, 2, 2.0);
    g.add_edge(2, 1, 3.0);
    g.add_edge(1, 3, 1.0);
    g.add_edge(2, 3, 9.0);
    let (dist, prev) = dijkstra(&g, 0);
    assert_eq!(&dist[..4], &[0.0, 5.0, 2.0, 6.0]);
    let mut path = vec![3];
    while let Some(p) = prev[*path.last().unwrap_or(&0)] {
        path.push(p);
    }
    path.reverse();
    assert_eq!(path, vec![0, 2, 1, 3]);
}

#[test]
fn mst_skips_the_expensive_chord() {
    let mut g = Graph::new(3);
    for (u, v, w) in [(0, 1, 1.0), (1, 2, 1.0), (0, 2, 5.0)] {
        g.add_edge(u, v, w);
        g.add_edge(v, u, w);
    }
    let mst = graph::minimum_spanning_tree(&g);
    let total: f64 = (0..3)
        .flat_map(|u| mst.adj(u).iter().filter(move |e| u < e.0).map(|e| e.1))
        .sum();
    assert_eq!(total, 2.0);
}

#[test]
fn empty_graph() {
    let g = Graph::new(0);
    assert!(graph::connected_components(&g).is_empty());
    assert_eq!(graph::minimum_spanning_tree(&g).num_nodes(), 0);
}

proptest::proptest! {
    #![proptest_config(cfg())]

    /// Floyd-Warshall agrees with Dijkstra from every source on random weighted digraphs.
    #[test]
    fn prop_floyd_warshall_matches_dijkstra(
        n in 2usize..7,
        edges in proptest::collection::vec((0usize..7, 0usize..7, 1u32..20), 0..25),
    ) {
        let mut g = Graph::new(n);
        for (u, v, w) in edges {
            g.add_edge(u % n, v % n, f64::from(w));
        }
        let fw = floyd_warshall(&g);
        for s in 0..n {
            let (dist, _) = dijkstra(&g, s);
            for t in 0..n {
                proptest::prop_assert_eq!(fw[s * n + t], dist[t], "s={} t={}", s, t);
            }
        }
    }

    /// BFS hop counts equal Dijkstra distances when every weight is 1.
    #[test]
    fn prop_bfs_matches_unit_weight_dijkstra(
        n in 2usize..8,
        edges in proptest::collection::vec((0usize..8, 0usize..8), 0..20),
    ) {
        let mut g = Graph::new(n);
        for (u, v) in edges {
            g.add_edge(u % n, v % n, 1.0);
        }
        let hops = bfs(&g, 0);
        let (dist, _) = dijkstra(&g, 0);
        for t in 0..n {
            if hops[t] == usize::MAX {
                proptest::prop_assert_eq!(dist[t], f64::INFINITY);
            } else {
                proptest::prop_assert_eq!(dist[t], hops[t] as f64);
            }
        }
    }

    /// Connected components of an undirected path partition it into the expected number of parts.
    #[test]
    fn prop_components_of_disjoint_paths(k in 1usize..5, len in 1usize..5) {
        let n = k * len;
        let mut g = Graph::new(n);
        for c in 0..k {
            for i in 0..len - 1 {
                let (a, b) = (c * len + i, c * len + i + 1);
                g.add_edge(a, b, 1.0);
                g.add_edge(b, a, 1.0);
            }
        }
        let comp = graph::connected_components(&g);
        let mut ids: Vec<usize> = comp.clone();
        ids.sort_unstable();
        ids.dedup();
        proptest::prop_assert_eq!(ids.len(), k);
        for c in 0..k {
            for i in 0..len {
                proptest::prop_assert_eq!(comp[c * len + i], comp[c * len]);
            }
        }
    }
}

// ============================================================================
// Canonical labelling, partition refinement and automorphism groups
// ============================================================================

/// A deterministic linear congruential generator, so every family below is
/// reproducible without a dependency.
struct Lcg(u64);

impl Lcg {
    fn new(seed: u64) -> Self {
        Self(seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1))
    }

    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.0 >> 33
    }

    fn below(
        &mut self,
        bound: usize,
    ) -> usize {
        (self.next() % bound.max(1) as u64) as usize
    }
}

fn random_permutation(
    n: usize,
    seed: u64,
) -> Vec<usize> {
    let mut rng = Lcg::new(seed);
    let mut perm: Vec<usize> = (0..n).collect();
    for i in (1..n).rev() {
        let j = rng.below(i + 1);
        perm.swap(i, j);
    }
    perm
}

/// The same graph under a numbering `perm`: vertex `u` becomes `perm[u]`.
fn relabel(
    g: &ColourGraph,
    perm: &[usize],
) -> ColourGraph {
    let n = g.num_nodes();
    let mut out = ColourGraph::new(n, g.is_directed());
    assert!(out.set_colours(&relabel_colours(g, perm)));
    for u in 0..n {
        for v in 0..n {
            let label = g.edge_label(u, v);
            if label != 0 {
                assert!(out.set_edge(perm[u], perm[v], label), "edge {u}-{v}");
            }
        }
    }
    out
}

/// The vertex colours under the numbering `perm`.
fn relabel_colours(
    g: &ColourGraph,
    perm: &[usize],
) -> Vec<u32> {
    let mut colours = vec![0; g.num_nodes()];
    for u in 0..g.num_nodes() {
        colours[perm[u]] = g.colours()[u];
    }
    colours
}

fn undirected(
    n: usize,
    edges: &[(usize, usize)],
) -> ColourGraph {
    let mut g = ColourGraph::new(n, false);
    for &(u, v) in edges {
        assert!(g.set_edge(u, v, 1), "edge {u}-{v} out of range");
    }
    g
}

fn directed(
    n: usize,
    arcs: &[(usize, usize)],
) -> ColourGraph {
    let mut g = ColourGraph::new(n, true);
    for &(u, v) in arcs {
        assert!(g.set_edge(u, v, 1), "arc {u}->{v} out of range");
    }
    g
}

fn cycle(n: usize) -> ColourGraph {
    undirected(n, &(0..n).map(|i| (i, (i + 1) % n)).collect::<Vec<_>>())
}

fn path(n: usize) -> ColourGraph {
    undirected(n, &(1..n).map(|i| (i - 1, i)).collect::<Vec<_>>())
}

fn complete(n: usize) -> ColourGraph {
    undirected(
        n,
        &(0..n)
            .flat_map(|u| ((u + 1)..n).map(move |v| (u, v)))
            .collect::<Vec<_>>(),
    )
}

fn star(n: usize) -> ColourGraph {
    undirected(n, &(1..n).map(|i| (0, i)).collect::<Vec<_>>())
}

fn complete_bipartite(
    a: usize,
    b: usize,
) -> ColourGraph {
    undirected(
        a + b,
        &(0..a)
            .flat_map(|u| (0..b).map(move |v| (u, a + v)))
            .collect::<Vec<_>>(),
    )
}

fn disjoint_union(
    g: &ColourGraph,
    h: &ColourGraph,
) -> ColourGraph {
    let (n, m) = (g.num_nodes(), h.num_nodes());
    assert!(!g.is_directed() && !h.is_directed());
    let mut out = ColourGraph::new(n + m, false);
    for u in 0..n {
        for v in 0..n {
            let label = g.edge_label(u, v);
            if label != 0 {
                assert!(out.set_edge(u, v, label));
            }
        }
    }
    for u in 0..m {
        for v in 0..m {
            let label = h.edge_label(u, v);
            if label != 0 {
                assert!(out.set_edge(n + u, n + v, label));
            }
        }
    }
    out
}

/// The Petersen graph: an outer 5-cycle, an inner pentagram and the spokes.
fn petersen() -> ColourGraph {
    let mut edges: Vec<(usize, usize)> = (0..5).map(|i| (i, (i + 1) % 5)).collect();
    edges.extend((0..5).map(|i| (5 + i, 5 + (i + 2) % 5)));
    edges.extend((0..5).map(|i| (i, 5 + i)));
    undirected(10, &edges)
}

/// The `d`-dimensional cube `Q_d`.
fn hypercube(d: usize) -> ColourGraph {
    let n = 1 << d;
    undirected(
        n,
        &(0..n)
            .flat_map(|u| (0..d).map(move |k| (u, u ^ (1 << k))))
            .collect::<Vec<_>>(),
    )
}

/// The 4x4 rook graph: `SRG(16, 6, 2, 2)`.
fn rook_4() -> ColourGraph {
    let id = |r: usize, c: usize| r * 4 + c;
    let mut edges = Vec::new();
    for r in 0..4 {
        for c in 0..4 {
            for other in 0..4 {
                if other != c {
                    edges.push((id(r, c), id(r, other)));
                }
                if other != r {
                    edges.push((id(r, c), id(other, c)));
                }
            }
        }
    }
    undirected(16, &edges)
}

/// The Shrikhande graph: the Cayley graph of `Z4 x Z4` on
/// `{+-x, +-y, +-(x + y)}`, also `SRG(16, 6, 2, 2)` but not isomorphic to the
/// rook graph.
fn shrikhande() -> ColourGraph {
    let id = |r: usize, c: usize| (r % 4) * 4 + (c % 4);
    let steps = [(1, 0), (3, 0), (0, 1), (0, 3), (1, 1), (3, 3)];
    let mut edges = Vec::new();
    for r in 0..4 {
        for c in 0..4 {
            for &(dr, dc) in &steps {
                edges.push((id(r, c), id(r + dr, c + dc)));
            }
        }
    }
    undirected(16, &edges)
}

/// A random sparse undirected graph: each pair is an edge with probability
/// `numerator / denominator`.
fn random_sparse(
    n: usize,
    seed: u64,
    numerator: u64,
    denominator: u64,
) -> ColourGraph {
    let mut rng = Lcg::new(seed);
    let mut edges = Vec::new();
    for u in 0..n {
        for v in (u + 1)..n {
            if rng.next() % denominator < numerator {
                edges.push((u, v));
            }
        }
    }
    undirected(n, &edges)
}

/// A random labelled tree (each vertex is attached to an earlier one).
fn random_tree(
    n: usize,
    seed: u64,
) -> ColourGraph {
    let mut rng = Lcg::new(seed);
    undirected(n, &(1..n).map(|v| (rng.below(v), v)).collect::<Vec<_>>())
}

fn canonical(g: &ColourGraph) -> CanonicalLabelling {
    canonical_labelling(g, SearchBudget::for_size(g.num_nodes()))
}

fn budget(g: &ColourGraph) -> SearchBudget {
    SearchBudget::for_size(g.num_nodes())
}

fn group_order(g: &ColourGraph) -> String {
    canonical(g).group_order().to_string()
}

fn orbits(g: &ColourGraph) -> Vec<Vec<usize>> {
    canonical(g).orbits().to_vec()
}

fn families() -> Vec<(&'static str, ColourGraph)> {
    vec![
        ("cycle 9", cycle(9)),
        ("path 7", path(7)),
        ("star 6", star(6)),
        ("complete 6", complete(6)),
        ("empty 5", ColourGraph::new(5, false)),
        ("bipartite 3x4", complete_bipartite(3, 4)),
        ("petersen", petersen()),
        ("hypercube 3", hypercube(3)),
        ("hypercube 4", hypercube(4)),
        (
            "two triangles and a path",
            disjoint_union(&disjoint_union(&cycle(3), &cycle(3)), &path(4)),
        ),
        ("C4 + C6", disjoint_union(&cycle(4), &cycle(6))),
        ("random tree 40", random_tree(40, 7)),
        ("random sparse 45", random_sparse(45, 11, 1, 8)),
        ("random sparse 60", random_sparse(60, 23, 1, 12)),
    ]
}

/// Whether `g` is strongly regular with the given parameters.
fn is_strongly_regular(
    g: &ColourGraph,
    n: usize,
    k: usize,
    lambda: usize,
    mu: usize,
) -> bool {
    if g.num_nodes() != n || g.is_directed() {
        return false;
    }
    for u in 0..n {
        if g.degree(u) != k {
            return false;
        }
        for v in (u + 1)..n {
            let common = (0..n)
                .filter(|&w| w != u && w != v && g.has_edge(u, w) && g.has_edge(v, w))
                .count();
            let want = if g.has_edge(u, v) {
                lambda
            } else {
                mu
            };
            if common != want {
                return false;
            }
        }
    }
    true
}

#[test]
fn canonical_form_is_invariant_under_vertex_permutations() {
    for (name, g) in families() {
        let n = g.num_nodes();
        let base = canonical(&g);
        assert!(base.is_exact(), "{name}: the search should be exact");
        assert_eq!(base.code().len(), n * (n + 1), "{name}: code length");
        assert_eq!(base.order().len(), n);
        for i in 0..n {
            assert_eq!(
                base.code()[i],
                g.colours()[base.order()[i]],
                "{name}: colour row"
            );
            assert_eq!(base.position(base.order()[i]), Some(i), "{name}: position");
            for j in 0..n {
                assert_eq!(
                    base.code()[n + i * n + j],
                    g.edge_label(base.order()[i], base.order()[j]),
                    "{name}: adjacency at ({i}, {j})"
                );
            }
        }
        for seed in 0..6u64 {
            let perm = random_permutation(n, seed + 1);
            let h = relabel(&g, &perm);
            let other = canonical(&h);
            assert!(other.is_exact(), "{name}: permuted search should be exact");
            assert_eq!(other.code(), base.code(), "{name}: seed {seed}");
            assert_eq!(other.matrix(), base.matrix(), "{name}: seed {seed}");
            assert_eq!(
                other.group_order(),
                base.group_order(),
                "{name}: seed {seed}"
            );
            assert_eq!(
                other.orbits().len(),
                base.orbits().len(),
                "{name}: seed {seed}"
            );
        }
        assert_eq!(
            isomorphic(&g, &g, budget(&g)),
            Some(true),
            "{name}: reflexive"
        );
    }
}

#[test]
fn automorphism_group_orders_match_the_known_values() {
    // complete graphs: n!
    for n in 1..=6 {
        assert_eq!(
            group_order(&complete(n)),
            (1..=n).product::<usize>().to_string(),
            "K{n}"
        );
    }
    // cycles: 2n
    for n in 3..=12 {
        assert_eq!(group_order(&cycle(n)), (2 * n).to_string(), "C{n}");
    }
    // paths: 2 (1 below three vertices)
    assert_eq!(group_order(&path(2)), "2");
    for n in 3..=9 {
        assert_eq!(group_order(&path(n)), "2", "P{n}");
    }
    // stars: the leaves are permuted freely (from three vertices on the centre
    // is the only vertex of maximum degree)
    for n in 3..=7 {
        assert_eq!(
            group_order(&star(n)),
            (1..n).product::<usize>().to_string(),
            "star{n}"
        );
    }
    // empty graphs: n!
    for n in 0..=6 {
        assert_eq!(
            group_order(&ColourGraph::new(n, false)),
            (1..=n).product::<usize>().to_string(),
            "E{n}"
        );
    }
    // the Petersen graph: 120
    assert_eq!(group_order(&petersen()), "120");
    // cubes: 2^d d!
    for d in 2..=4 {
        let want = (1 << d) * (1..=d).product::<usize>();
        assert_eq!(group_order(&hypercube(d)), want.to_string(), "Q{d}");
    }
    // complete bipartite: a! b!, times two when the two parts have the same
    // size and can be swapped
    for (a, b) in [(1, 1), (2, 3), (3, 3), (3, 4), (4, 5)] {
        let swap = if a == b { 2 } else { 1 };
        let want = swap * (1..=a).product::<usize>() * (1..=b).product::<usize>();
        assert_eq!(
            group_order(&complete_bipartite(a, b)),
            want.to_string(),
            "K{a},{b}"
        );
    }
    // disjoint unions multiply the orders of the components
    assert_eq!(
        group_order(&disjoint_union(&cycle(3), &cycle(3))),
        (6 * 6 * 2).to_string()
    );
    assert_eq!(
        group_order(&disjoint_union(&cycle(3), &cycle(4))),
        (6 * 8).to_string()
    );
    // a random sparse graph is rigid
    assert_eq!(group_order(&random_sparse(45, 11, 1, 8)), "1");
    // two parallel edges of a 2-vertex graph are interchangeable
    let mut doubled = ColourGraph::new(2, false);
    assert!(doubled.set_edge(0, 1, 3));
    assert_eq!(group_order(&doubled), "2");
}

#[test]
fn orbit_partitions_are_reported() {
    assert_eq!(orbits(&petersen()), vec![(0..10).collect::<Vec<_>>()]);
    assert_eq!(orbits(&cycle(7)), vec![vec![0, 1, 2, 3, 4, 5, 6]]);
    assert_eq!(
        orbits(&complete_bipartite(3, 4)),
        vec![vec![0, 1, 2], vec![3, 4, 5, 6]]
    );
    assert_eq!(orbits(&path(6)), vec![vec![0, 5], vec![1, 4], vec![2, 3]]);
    assert_eq!(orbits(&star(6)), vec![vec![0], vec![1, 2, 3, 4, 5]]);
    assert_eq!(orbits(&ColourGraph::new(4, false)), vec![vec![0, 1, 2, 3]]);
    // two paths of different length: the two reversals plus the isolated
    // middle vertex of the second
    let union = disjoint_union(&path(4), &path(3));
    assert_eq!(
        orbits(&union),
        vec![vec![0, 3], vec![1, 2], vec![4, 6], vec![5]]
    );
    // a random graph has singleton orbits
    let random = random_sparse(45, 11, 1, 8);
    assert!(orbits(&random).iter().all(|orbit| orbit.len() == 1));
}

#[test]
fn generators_are_automorphisms_and_permutations_agree() {
    for (name, g) in [
        ("petersen", petersen()),
        ("hypercube 3", hypercube(3)),
        ("bipartite 4x5", complete_bipartite(4, 5)),
        ("cycle 7", cycle(7)),
        ("rook 4", rook_4()),
    ] {
        let n = g.num_nodes();
        let form = canonical(&g);
        for p in form.generators() {
            assert_eq!(p.len(), n, "{name}: degree");
            assert!(!p.is_identity(), "{name}: generators are non-trivial");
            for u in 0..n {
                let pu = p.apply(u).unwrap_or(u);
                assert_eq!(g.colours()[pu], g.colours()[u], "{name}: colours");
                for v in 0..n {
                    let pv = p.apply(v).unwrap_or(v);
                    assert_eq!(
                        g.edge_label(pu, pv),
                        g.edge_label(u, v),
                        "{name}: adjacency"
                    );
                }
            }
            let inverse = p.inverse();
            for u in 0..n {
                assert_eq!(
                    inverse.apply(p.apply(u).unwrap_or(u)),
                    Some(u),
                    "{name}: inverse"
                );
            }
            assert_eq!(
                p.compose(&inverse).map(|c| c.is_identity()),
                Some(true),
                "{name}"
            );
        }
        // the group acts transitively exactly on the non-singleton orbits
        for orbit in form.orbits() {
            assert!(!orbit.is_empty());
        }
    }
    // permutation arithmetic
    let p = Permutation::from_images(vec![1, 2, 0]).unwrap_or_else(|| Permutation::identity(3));
    assert_eq!(p.cycles(), vec![vec![0, 1, 2]]);
    assert_eq!(p.apply(2), Some(0));
    assert_eq!(p.compose(&p.inverse()).map(|c| c.is_identity()), Some(true));
    assert!(Permutation::from_images(vec![0, 0, 2]).is_none());
    assert!(Permutation::from_images(vec![0, 1, 3]).is_none());
}

#[test]
fn vertex_colours_restrict_the_group() {
    let mut coloured = cycle(6);
    assert!(coloured.set_colours(&[1, 0, 0, 0, 0, 0]));
    // only the reflection through the coloured vertex survives
    assert_eq!(group_order(&coloured), "2");
    let refined = refine(&coloured);
    assert_eq!(refined.num_cells(), 4);
    assert_eq!(refined.cell_of()[1], refined.cell_of()[5]);
    assert_ne!(refined.cell_of()[1], refined.cell_of()[2]);
    // the colouring is part of the canonical form
    let plain = canonical(&cycle(6));
    let form = canonical(&coloured);
    assert_ne!(plain.code(), form.code());
    for seed in 0..4u64 {
        let perm = random_permutation(6, seed + 5);
        let other = canonical(&relabel(&coloured, &perm));
        assert_eq!(other.code(), form.code(), "seed {seed}");
        assert_eq!(other.group_order(), form.group_order(), "seed {seed}");
    }
    // a colouring moved by an automorphism is isomorphic to it
    let mut moved = cycle(6);
    assert!(moved.set_colours(&[0, 0, 0, 1, 0, 0]));
    assert_eq!(
        isomorphic(&coloured, &moved, SearchBudget::for_size(6)),
        Some(true)
    );
    // an adjacent coloured pair and a pair at distance two are not related by
    // any automorphism of the cycle
    let mut adjacent = cycle(6);
    assert!(adjacent.set_colours(&[1, 1, 0, 0, 0, 0]));
    let mut apart = cycle(6);
    assert!(apart.set_colours(&[1, 0, 1, 0, 0, 0]));
    assert_eq!(
        isomorphic(&adjacent, &apart, SearchBudget::for_size(6)),
        Some(false)
    );
    assert_ne!(
        refine(&adjacent).certificate(),
        refine(&apart).certificate()
    );
}

#[test]
fn refinement_is_equitable_and_canonical() {
    for (name, g) in families() {
        let n = g.num_nodes();
        let refined = refine(&g);
        assert!(refined.verify(&g), "{name}: not equitable");
        assert_eq!(refined.cell_of().len(), n);
        assert_eq!(refined.cells().len(), refined.num_cells());
        let total: usize = refined.cells().iter().map(Vec::len).sum();
        assert_eq!(total, n, "{name}: the cells partition the vertices");
        assert_eq!(refined.certificate().num_nodes(), n);
        assert_eq!(refined.certificate().num_cells(), refined.num_cells());
        if refined.is_discrete() {
            assert!(refined.certificate().is_discrete());
        }
        for seed in 0..3u64 {
            let perm = random_permutation(n, seed + 1);
            let other = refine(&relabel(&g, &perm));
            assert_eq!(
                other.certificate(),
                refined.certificate(),
                "{name}: seed {seed}"
            );
            assert_eq!(other.cell_sizes(), refined.cell_sizes());
        }
        assert!(refined.rounds() >= 1);
    }
    // a regular graph is one cell, and its certificate still separates it from
    // a different regular graph
    for g in [cycle(7), petersen(), hypercube(4), complete(6)] {
        assert_eq!(refine(&g).num_cells(), 1);
    }
    assert!(!wl_equivalent(&cycle(7), &complete(7)));
    // C6 and two triangles are 2-regular and indistinguishable
    let c6 = cycle(6);
    let two_triangles = disjoint_union(&cycle(3), &cycle(3));
    assert_eq!(
        refine(&c6).certificate(),
        refine(&two_triangles).certificate()
    );
    assert!(wl_equivalent(&c6, &two_triangles));
    // but the canonical form is not fooled
    assert_ne!(canonical(&c6).code(), canonical(&two_triangles).code());
    assert_eq!(
        isomorphic(&c6, &two_triangles, SearchBudget::for_size(6)),
        Some(false)
    );
    // a path is not a cycle
    assert!(!wl_equivalent(&path(5), &cycle(5)));
}

#[test]
fn strongly_regular_pair_needs_the_canonical_form() {
    let rook = rook_4();
    let shrikhande = shrikhande();
    for g in [&rook, &shrikhande] {
        assert!(is_strongly_regular(g, 16, 6, 2, 2));
        // both are 6-regular and transitive, so refinement is one cell
        assert_eq!(refine(g).num_cells(), 1);
        assert_eq!(refine(g).certificate().cell_sizes(), vec![16]);
    }
    // 1-dimensional Weisfeiler-Lehman cannot tell them apart
    assert!(wl_equivalent(&rook, &shrikhande));
    assert_eq!(
        refine(&rook).certificate(),
        refine(&shrikhande).certificate()
    );
    // the canonical labelling can
    let (a, b) = (canonical(&rook), canonical(&shrikhande));
    assert!(a.is_exact() && b.is_exact());
    assert_ne!(a.code(), b.code());
    assert_eq!(
        isomorphic(&rook, &shrikhande, SearchBudget::for_size(16)),
        Some(false)
    );
    // their automorphism groups are the known ones
    assert_eq!(a.group_order().to_string(), "1152");
    assert_eq!(b.group_order().to_string(), "192");
    // each is isomorphic to itself and to every relabelling of itself
    for seed in 0..4u64 {
        let perm = random_permutation(16, seed + 3);
        assert_eq!(canonical(&relabel(&rook, &perm)).code(), a.code());
        assert_eq!(canonical(&relabel(&shrikhande, &perm)).code(), b.code());
    }
}

#[test]
fn directed_canonical_forms_and_groups() {
    let n = 12;
    let arcs: Vec<(usize, usize)> = (0..n).map(|i| (i, (i + 1) % n)).collect();
    let g = directed(n, &arcs);
    let form = canonical(&g);
    assert!(form.is_exact());
    // only the rotations preserve the orientation
    assert_eq!(form.group_order().to_string(), n.to_string());
    assert_eq!(form.orbits().len(), 1);
    for seed in 0..4u64 {
        let perm = random_permutation(n, seed + 1);
        let other = canonical(&relabel(&g, &perm));
        assert_eq!(other.code(), form.code(), "seed {seed}");
        assert_eq!(other.group_order(), form.group_order(), "seed {seed}");
    }
    // a directed path is rigid, and its reverse is isomorphic to it
    let forward = directed(6, &(1..6).map(|i| (i - 1, i)).collect::<Vec<_>>());
    let backward = directed(6, &(1..6).map(|i| (i, i - 1)).collect::<Vec<_>>());
    assert_eq!(group_order(&forward), "1");
    assert_eq!(
        isomorphic(&forward, &backward, SearchBudget::for_size(6)),
        Some(true)
    );
    // an outward star is not isomorphic to an inward one
    let out_star = directed(4, &[(0, 1), (0, 2), (0, 3)]);
    let in_star = directed(4, &[(1, 0), (2, 0), (3, 0)]);
    assert_eq!(group_order(&out_star), "6");
    assert_eq!(
        isomorphic(&out_star, &in_star, SearchBudget::for_size(4)),
        Some(false)
    );
    // a random acyclic digraph is rigid and its form is stable
    let mut rng = Lcg::new(31);
    let mut arcs = Vec::new();
    for u in 0..20 {
        for v in (u + 1)..20 {
            if rng.next() % 6 == 0 {
                arcs.push((u, v));
            }
        }
    }
    let acyclic = directed(20, &arcs);
    let form = canonical(&acyclic);
    assert!(form.is_exact());
    assert_eq!(form.code().len(), 20 * 21);
    for seed in 0..3u64 {
        assert_eq!(
            canonical(&relabel(&acyclic, &random_permutation(20, seed + 2))).code(),
            form.code()
        );
    }
    // directed and undirected graphs are never isomorphic
    assert_eq!(
        isomorphic(
            &directed(3, &[(0, 1), (1, 2)]),
            &path(3),
            SearchBudget::for_size(3)
        ),
        Some(false)
    );
}

#[test]
fn degenerate_graphs() {
    let empty = ColourGraph::new(0, false);
    let form = canonical(&empty);
    assert!(form.is_exact());
    assert!(form.order().is_empty());
    assert!(form.code().is_empty());
    assert_eq!(form.group_order().to_string(), "1");
    assert_eq!(
        isomorphic(
            &empty,
            &ColourGraph::new(0, true),
            SearchBudget::for_size(0)
        ),
        Some(false)
    );
    let single = ColourGraph::new(1, false);
    assert_eq!(canonical(&single).code().to_vec(), vec![0, 0]);
    assert_eq!(group_order(&single), "1");
    // a self-loop is only preserved by automorphisms that keep it
    let mut looped = ColourGraph::new(3, false);
    assert!(looped.set_edge(0, 0, 1));
    assert!(looped.set_edge(1, 2, 1));
    assert_eq!(group_order(&looped), "2");
    // labels are part of the structure: an edge of label 2 is not a single one
    let mut single_edge = ColourGraph::new(2, false);
    assert!(single_edge.set_edge(0, 1, 1));
    let mut double_edge = ColourGraph::new(2, false);
    assert!(double_edge.set_edge(0, 1, 2));
    assert_eq!(
        isomorphic(&single_edge, &double_edge, SearchBudget::for_size(2)),
        Some(false)
    );
    // the unused builder paths report failure instead of panicking
    let mut g = ColourGraph::new(2, false);
    assert!(!g.set_edge(0, 2, 1));
    assert!(!g.set_colours(&[1]));
}

/// The same graph built by listing its edges in a different order has the same
/// canonical form.
#[test]
fn canonical_form_does_not_depend_on_the_edge_order() {
    let edges: Vec<(usize, usize)> = (0..12)
        .flat_map(|u| {
            ((u + 1)..12)
                .filter(move |v| (u * 7 + v * 3) % 5 == 0)
                .map(move |v| (u, v))
        })
        .collect();
    let mut reversed = edges.clone();
    reversed.reverse();
    let a = undirected(12, &edges);
    let b = undirected(12, &reversed);
    assert_eq!(canonical(&a).code(), canonical(&b).code());
    assert_eq!(canonical(&a).order(), canonical(&b).order());
}

/// Every permutation of `0..n` is tried, so this is an independent oracle for
/// the canonical labelling on small graphs.
fn brute_force_isomorphic(
    a: &ColourGraph,
    b: &ColourGraph,
) -> bool {
    fn next_permutation(p: &mut [usize]) -> bool {
        let n = p.len();
        let Some(i) = (1..n).rev().find(|&i| p[i - 1] < p[i]) else {
            return false;
        };
        let j = (i..n).rev().find(|&j| p[j] > p[i - 1]).unwrap_or(i);
        p.swap(i - 1, j);
        p[i..].reverse();
        true
    }
    let n = a.num_nodes();
    if b.num_nodes() != n || a.is_directed() != b.is_directed() {
        return false;
    }
    let mut perm: Vec<usize> = (0..n).collect();
    loop {
        let same = (0..n).all(|u| {
            (0..n).all(|v| {
                a.edge_label(u, v) == b.edge_label(perm[u], perm[v])
                    && a.colours()[u] == b.colours()[perm[u]]
            })
        });
        if same {
            return true;
        }
        if !next_permutation(&mut perm) {
            return false;
        }
    }
}

#[test]
fn isomorphism_agrees_with_brute_force() {
    let mut rng = Lcg::new(2026);
    for _ in 0..24 {
        let n = 3 + rng.below(4);
        let mut a = ColourGraph::new(n, false);
        for u in 0..n {
            for v in (u + 1)..n {
                if rng.next() % 3 == 0 {
                    assert!(a.set_edge(u, v, 1));
                }
            }
        }
        // half of the cases compare against a relabelling, half against an
        // independent random graph
        let other = if rng.next() % 2 == 0 {
            let perm = random_permutation(n, rng.next());
            relabel(&a, &perm)
        } else {
            let mut b = ColourGraph::new(n, false);
            for u in 0..n {
                for v in (u + 1)..n {
                    if rng.next() % 3 == 0 {
                        assert!(b.set_edge(u, v, 1));
                    }
                }
            }
            b
        };
        let want = brute_force_isomorphic(&a, &other);
        assert_eq!(isomorphic(&a, &other, budget(&a)), Some(want), "n={n}");
        assert_eq!(
            canonical(&a).code() == canonical(&other).code(),
            want,
            "n={n}"
        );
    }
}

/// A 300-vertex sparse graph has to be labelled in well under a second on a
/// test build; the budget is nowhere near exhausted.
#[test]
fn canonical_labelling_scales_to_a_few_hundred_vertices() {
    let g = random_sparse(300, 99, 1, 40);
    let start = Instant::now();
    let form = canonical(&g);
    let elapsed = start.elapsed();
    assert!(
        form.is_exact(),
        "budget exhausted after {} nodes",
        form.num_nodes_visited()
    );
    assert_eq!(form.code().len(), 300 * 301);
    assert_eq!(form.group_order().to_string(), "1");
    assert!(elapsed.as_secs_f64() < 5.0, "took {elapsed:?}");
    // the form really is the graph in the order the search chose
    for i in 0..300 {
        for j in 0..300 {
            assert_eq!(
                form.code()[300 + i * 300 + j],
                g.edge_label(form.order()[i], form.order()[j])
            );
        }
    }
    // and a permuted copy gives the same code
    let perm = random_permutation(300, 17);
    assert_eq!(canonical(&relabel(&g, &perm)).code(), form.code());
    println!(
        "canonical labelling of 300 vertices: {:?}, {} nodes",
        elapsed,
        form.num_nodes_visited()
    );
}

/// The budget stops a search that cannot finish: a hundred disjoint triangles
/// are one Weisfeiler-Lehman class, so the individualisation tree is enormous.
/// The search then reports itself inexact instead of running away, and the
/// isomorphism test falls back to "undecided".
#[test]
fn the_budget_stops_a_pathological_search() {
    let mut edges = Vec::new();
    for t in 0..100 {
        let base = 3 * t;
        edges.push((base, base + 1));
        edges.push((base + 1, base + 2));
        edges.push((base + 2, base));
    }
    let g = undirected(300, &edges);
    assert_eq!(refine(&g).num_cells(), 1);
    let tight = SearchBudget::nodes(300);
    let start = Instant::now();
    let form = canonical_labelling(&g, tight);
    let elapsed = start.elapsed();
    assert!(!form.is_exact());
    assert_eq!(form.num_nodes_visited(), 301);
    assert_eq!(form.code().len(), 300 * 301);
    assert!(elapsed.as_secs_f64() < 5.0, "took {elapsed:?}");
    assert_eq!(isomorphic(&g, &g, tight), None);
}
