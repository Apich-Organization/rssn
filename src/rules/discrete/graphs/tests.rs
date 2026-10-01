use super::super::test_util::nested;
use super::super::test_util::nums;
use super::super::test_util::s;
use super::super::test_util::Lcg;

type Edge = (usize, usize, i64);

/// `graph(n, ...)` / `digraph(n, ...)` with explicit weights.
fn gt(
    n: usize,
    edges: &[Edge],
    directed: bool,
) -> String {
    let list: Vec<String> = edges.iter().map(|&(u, v, w)| format!("list({u}, {v}, {w})")).collect();
    format!("{}({n}, list({}))", if directed { "digraph" } else { "graph" }, list.join(", "))
}

/// The edges of a graph term, each with its weight (default 1).
#[allow(clippy::cast_sign_loss)] // operand is non-negative by construction (index/count)
fn parse_graph(text: &str) -> (usize, Vec<Edge>) {
    let n = nums(text).first().copied().unwrap_or(0) as usize;
    let start = text.find("list(").unwrap_or(0);
    let body = &text[start..text.len() - 1];
    let edges = nested(body)
        .into_iter()
        .map(|r| (r[0] as usize, r[1] as usize, r.get(2).copied().unwrap_or(1)))
        .collect();
    (n, edges)
}

/// Edge set normalised for comparison.
fn edge_set(
    edges: &[Edge],
    directed: bool,
) -> Vec<Edge> {
    let mut out: Vec<Edge> = edges.iter().map(|&(u, v, w)| if directed || u <= v { (u, v, w) } else { (v, u, w) }).collect();
    out.sort_unstable();
    out
}

#[allow(clippy::cast_sign_loss)] // operand is non-negative by construction (index/count)
fn random_graph(
    rng: &mut Lcg,
    n: usize,
    density: u64,
    directed: bool,
    max_w: i64,
) -> Vec<Edge> {
    let mut edges = Vec::new();
    for u in 0..n {
        for v in 0..n {
            if u == v || (!directed && u > v) {
                continue;
            }
            if rng.next(100) < density {
                edges.push((u, v, 1 + rng.next(max_w as u64) as i64));
            }
        }
    }
    edges
}

/// The top-level items of `list(a, b, ...)`.
fn top_level(text: &str) -> Vec<String> {
    let inner = text.strip_prefix("list(").and_then(|t| t.strip_suffix(')')).unwrap_or("");
    let mut parts = Vec::new();
    let (mut depth, mut start) = (0_i32, 0);
    for (i, c) in inner.char_indices() {
        match c {
            | '(' => depth += 1,
            | ')' => depth -= 1,
            | ',' if depth == 0 => {
                parts.push(inner[start..i].trim().to_string());
                start = i + 1;
            },
            | _ => {},
        }
    }
    parts.push(inner[start..].trim().to_string());
    parts
}

fn matrix_of(text: &str) -> Vec<Vec<i64>> {
    nested(&text.replace("oo", "9999999"))
}

fn adjacency(
    n: usize,
    edges: &[Edge],
    directed: bool,
) -> Vec<Vec<i64>> {
    let mut m = vec![vec![0; n]; n];
    for &(u, v, w) in edges {
        m[u][v] += w;
        if !directed {
            m[v][u] += w;
        }
    }
    m
}

fn union_find_components(
    n: usize,
    edges: &[Edge],
    skip_edge: Option<usize>,
    skip_vertex: Option<usize>,
) -> usize {
    fn find(
        p: &mut [usize],
        x: usize,
    ) -> usize {
        if p[x] != x {
            let r = find(p, p[x]);
            p[x] = r;
        }
        p[x]
    }
    let mut parent: Vec<usize> = (0..n).collect();
    for (i, &(u, v, _)) in edges.iter().enumerate() {
        if Some(i) == skip_edge || Some(u) == skip_vertex || Some(v) == skip_vertex {
            continue;
        }
        let (a, b) = (find(&mut parent, u), find(&mut parent, v));
        parent[a] = b;
    }
    (0..n).filter(|&x| Some(x) != skip_vertex && find(&mut parent, x) == x).count()
}

// ---------------------------------------------------------------- structure

#[test]
fn constructors_and_queries() {
    assert_eq!(s("graph_complete(4)"), "graph(4, list(list(0, 1), list(0, 2), list(0, 3), list(1, 2), list(1, 3), list(2, 3)))");
    assert_eq!(s("graph_path(3)"), "graph(3, list(list(0, 1), list(1, 2)))");
    assert_eq!(s("graph_cycle(3)"), "graph(3, list(list(0, 1), list(1, 2), list(2, 0)))");
    assert_eq!(s("graph_empty(2)"), "graph(2, list())");
    assert_eq!(s("graph_nodes(graph_path(3))"), "list(0, 1, 2)");
    assert_eq!(s("graph_node_count(graph_complete(5))"), "5");
    assert_eq!(s("graph_is_directed(graph_path(3))"), "false");
    assert_eq!(s("graph_is_directed(digraph(2, list(list(0, 1))))"), "true");
    assert_eq!(s("graph_node_id(graph_path(3), 2)"), "2");
    assert_eq!(s("graph_node_id(graph_path(3), 3)"), "graph_node_id(graph(3, list(list(0, 1), list(1, 2))), 3)");
    assert_eq!(s("graph_add_node(graph_path(2))"), "graph(3, list(list(0, 1)))");
    assert_eq!(s("graph_add_edge(graph_path(3), 0, 2)"), "graph(3, list(list(0, 1), list(1, 2), list(0, 2)))");
    assert_eq!(s("graph_add_edge(graph_path(3), 0, 2, w)"), "graph(3, list(list(0, 1), list(1, 2), list(0, 2, w)))");
    assert_eq!(s("graph_neighbors(graph(3, list(list(0, 1, a), list(0, 2))), 0)"), "list(list(1, a), list(2, 1))");
    assert_eq!(s("graph_out_degree(digraph(3, list(list(0, 1), list(0, 2), list(1, 2))), 0)"), "2");
    assert_eq!(s("graph_in_degree(digraph(3, list(list(0, 1), list(0, 2), list(1, 2))), 2)"), "2");
    assert_eq!(s("graph_in_degree(digraph(3, list(list(0, 1), list(0, 2), list(1, 2))), 0)"), "0");
    assert_eq!(s("graph_out_degree(graph_complete(4), 1)"), "3");
    assert_eq!(s("graph_edges(graph(3, list(list(2, 0, x), list(1, 1))))"), "list(list(0, 2, x), list(1, 1, 1))");
    assert_eq!(s("graph_edges(digraph(3, list(list(2, 0, x))))"), "list(list(2, 0, x))");
    // hyperedges are stored with the graph
    let h = s("graph_add_hyperedge(graph_path(3), list(0, 1, 2), w)");
    assert_eq!(h, "graph(3, list(list(0, 1), list(1, 2)), list(list(list(0, 1, 2), w)))");
    assert_eq!(s(&format!("graph_hyperedges({h})")), "list(list(list(0, 1, 2), w))");
    assert_eq!(s("graph_hyperedges(graph_path(3))"), "list()");
    // malformed graphs stay
    assert_eq!(s("graph_nodes(graph(2, list(list(0, 5))))"), "graph_nodes(graph(2, list(list(0, 5))))");
}

#[test]
fn matrices_with_symbolic_weights() {
    let tri = "graph(3, list(list(0, 1, a), list(1, 2, b), list(0, 2, c)))";
    assert_eq!(s(&format!("graph_adjacency({tri})")), "list(list(0, a, c), list(a, 0, b), list(c, b, 0))");
    assert_eq!(s("graph_adjacency(digraph(2, list(list(0, 1, 3))))"), "list(list(0, 3), list(0, 0))");
    // parallel edges add up
    assert_eq!(s("graph_adjacency(graph(2, list(list(0, 1, 2), list(0, 1, 3))))"), "list(list(0, 5), list(5, 0))");
    // composes with the linear algebra rules: det [[0,a,c],[a,0,b],[c,b,0]] = 2 a b c
    let det = s(&format!("det(graph_adjacency({tri}))"));
    assert!(det.contains('a') && det.contains('b') && det.contains('c'), "{det}");
    let numeric = s("det(graph_adjacency(graph(3, list(list(0, 1, 2), list(1, 2, 3), list(0, 2, 5)))))");
    assert_eq!(numeric, "60");
    // incidence: undirected 1/1, directed -1/+1
    assert_eq!(s("graph_incidence(graph_path(3))"), "list(list(1, 0), list(1, 1), list(0, 1))");
    assert_eq!(s("graph_incidence(digraph(3, list(list(0, 1), list(1, 2))))"), "list(list(-1, 0), list(1, -1), list(0, 1))");
    // Laplacian = D - A with weighted degrees
    assert_eq!(s("graph_laplacian(graph_path(3))"), "list(list(1, -1, 0), list(-1, 2, -1), list(0, -1, 1))");
    assert_eq!(
        s("graph_laplacian(graph(2, list(list(0, 1, w))))"),
        "list(list(w, -w), list(-w, w))"
    );
    // eigenvalues of the Laplacian of K3 are 0, 3, 3
    assert_eq!(s("eigenvals(graph_laplacian(graph_complete(3)))"), "list(0, 3, 3)");
    // matrix tree theorem: a cofactor of the Laplacian counts spanning trees (K4: 16)
    for (n, trees) in [(3_usize, 3_i64), (4, 16), (5, 125)] {
        let lap = matrix_of(&s(&format!("graph_laplacian(graph_complete({n}))")));
        let minor: Vec<String> = lap[1..]
            .iter()
            .map(|row| format!("list({})", row[1..].iter().map(ToString::to_string).collect::<Vec<_>>().join(", ")))
            .collect();
        assert_eq!(s(&format!("det(list({}))", minor.join(", "))), trees.to_string());
    }
    // weighted: a triangle with weights a, b, c has a*b + b*c + a*c weighted spanning trees
    let lap = s(&format!("graph_laplacian({tri})"));
    assert!(lap.contains("a + c") || lap.contains("c + a"), "{lap}");
}

// ---------------------------------------------------------------- traversal

#[test]
fn traversals_and_components() {
    let g = gt(6, &[(0, 1, 1), (0, 2, 1), (1, 3, 1), (2, 3, 1), (4, 5, 1)], false);
    assert_eq!(s(&format!("graph_bfs({g}, 0)")), "list(0, 1, 2, 3)");
    assert_eq!(s(&format!("graph_dfs({g}, 0)")), "list(0, 1, 3, 2)");
    assert_eq!(s(&format!("graph_components({g})")), "list(list(0, 1, 2, 3), list(4, 5))");
    assert_eq!(s(&format!("graph_is_connected({g})")), "false");
    assert_eq!(s("graph_is_connected(graph_cycle(5))"), "true");
    assert_eq!(s("graph_is_connected(graph_empty(0))"), "false");
    assert_eq!(s(&format!("graph_bfs({g}, 9)")), format!("graph_bfs({g}, 9)"));
    // components against union-find on random graphs
    let mut rng = Lcg(1);
    for _ in 0..25 {
        let n = 2 + rng.next(8) as usize;
        let edges = random_graph(&mut rng, n, 20, false, 1);
        let comps = nested(&s(&format!("graph_components({})", gt(n, &edges, false))));
        assert_eq!(comps.len(), union_find_components(n, &edges, None, None));
        let mut all: Vec<i64> = comps.into_iter().flatten().collect();
        all.sort_unstable();
        assert_eq!(all, (0..n as i64).collect::<Vec<_>>());
    }
}

#[test]
fn strongly_connected_components() {
    let g = gt(6, &[(0, 1, 1), (1, 2, 1), (2, 0, 1), (2, 3, 1), (3, 4, 1), (4, 3, 1)], true);
    assert_eq!(s(&format!("graph_scc({g})")), "list(list(0, 1, 2), list(3, 4), list(5))");
    let mut rng = Lcg(2);
    for _ in 0..25 {
        let n = 2 + rng.next(6) as usize;
        let edges = random_graph(&mut rng, n, 25, true, 1);
        // reference: mutual reachability by transitive closure
        let mut reach = vec![vec![false; n]; n];
        #[allow(clippy::needless_range_loop)] // index is used for more than one array / arithmetic; iterator form would not be clearer
        for i in 0..n {
            reach[i][i] = true;
        }
        for &(u, v, _) in &edges {
            reach[u][v] = true;
        }
        for k in 0..n {
            for i in 0..n {
                for j in 0..n {
                    reach[i][j] |= reach[i][k] && reach[k][j];
                }
            }
        }
        let mut want: Vec<Vec<i64>> = Vec::new();
        #[allow(clippy::needless_range_loop)] // index is used for more than one array / arithmetic; iterator form would not be clearer
        for i in 0..n {
            if want.iter().any(|c| c.contains(&(i as i64))) {
                continue;
            }
            want.push((0..n).filter(|&j| reach[i][j] && reach[j][i]).map(|j| j as i64).collect());
        }
        let got = nested(&s(&format!("graph_scc({})", gt(n, &edges, true))));
        assert_eq!(got, want, "{edges:?}");
    }
}

#[test]
#[allow(clippy::cast_sign_loss)] // operand is non-negative by construction (index/count)
fn cycles_bridges_and_articulation_points() {
    assert_eq!(s("graph_has_cycle(graph_path(5))"), "false");
    assert_eq!(s("graph_has_cycle(graph_cycle(5))"), "true");
    assert_eq!(s("graph_has_cycle(graph(2, list(list(0, 1), list(0, 1))))"), "true");
    assert_eq!(s("graph_has_cycle(graph(2, list(list(1, 1))))"), "true");
    assert_eq!(s("graph_has_cycle(digraph(3, list(list(0, 1), list(1, 2))))"), "false");
    assert_eq!(s("graph_has_cycle(digraph(3, list(list(0, 1), list(1, 2), list(2, 0))))"), "true");
    // two triangles joined by an edge: that edge is the only bridge
    let g = gt(6, &[(0, 1, 1), (1, 2, 1), (2, 0, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1), (5, 3, 1)], false);
    assert_eq!(s(&format!("graph_bridges({g})")), "list(list(list(2, 3)), list(2, 3))");
    assert_eq!(s("graph_bridges(graph_cycle(4))"), "list(list(), list())");
    assert_eq!(s("graph_bridges(graph_path(3))"), "list(list(list(1, 2), list(0, 1)), list(1))");
    let mut rng = Lcg(3);
    for _ in 0..40 {
        let n = 2 + rng.next(7) as usize;
        let edges = random_graph(&mut rng, n, 30, false, 1);
        let base = union_find_components(n, &edges, None, None);
        let want_bridges: Vec<(usize, usize)> = (0..edges.len())
            .filter(|&i| union_find_components(n, &edges, Some(i), None) > base)
            .map(|i| (edges[i].0.min(edges[i].1), edges[i].0.max(edges[i].1)))
            .collect();
        let want_ap: Vec<usize> = (0..n)
            .filter(|&v| union_find_components(n, &edges, None, Some(v)) > base - usize::from(!edges.iter().any(|e| e.0 == v || e.1 == v)))
            .collect();
        let out = s(&format!("graph_bridges({})", gt(n, &edges, false)));
        let parts = top_level(&out);
        let mut got_bridges: Vec<(usize, usize)> = nested(&parts[0])
            .iter()
            .map(|e| (e[0].min(e[1]) as usize, e[0].max(e[1]) as usize))
            .collect();
        got_bridges.sort_unstable();
        let mut want_bridges = want_bridges;
        want_bridges.sort_unstable();
        assert_eq!(got_bridges, want_bridges, "{edges:?}");
        let got_ap: Vec<usize> = nums(&parts[1]).into_iter().map(|x| x as usize).collect();
        assert_eq!(got_ap, want_ap, "{edges:?}");
        // cycle test against the edge count identity
        let has_cycle = edges.len() > n - base;
        assert_eq!(s(&format!("graph_has_cycle({})", gt(n, &edges, false))), has_cycle.to_string());
    }
}

// ---------------------------------------------------------------- spanning trees

fn brute_force_mst_weight(
    n: usize,
    edges: &[Edge],
) -> Option<i64> {
    let m = edges.len();
    let mut best: Option<i64> = None;
    for mask in 0u32..(1 << m) {
        if mask.count_ones() as usize != n - 1 {
            continue;
        }
        let chosen: Vec<Edge> = (0..m).filter(|&i| mask >> i & 1 == 1).map(|i| edges[i]).collect();
        if union_find_components(n, &chosen, None, None) == 1 {
            let w: i64 = chosen.iter().map(|e| e.2).sum();
            best = Some(best.map_or(w, |b: i64| b.min(w)));
        }
    }
    best
}

fn tree_weight(text: &str) -> i64 {
    nested(text).iter().map(|e| e[2]).sum()
}

#[test]
fn minimum_spanning_trees() {
    let g = gt(4, &[(0, 1, 1), (1, 2, 2), (2, 3, 3), (0, 3, 10), (0, 2, 5)], false);
    assert_eq!(s(&format!("graph_kruskal({g})")), "list(list(0, 1, 1), list(1, 2, 2), list(2, 3, 3))");
    assert_eq!(s(&format!("graph_prim({g}, 0)")), "list(list(0, 1, 1), list(1, 2, 2), list(2, 3, 3))");
    let mut rng = Lcg(4);
    for _ in 0..30 {
        let n = 3 + rng.next(4) as usize;
        let edges = random_graph(&mut rng, n, 60, false, 9);
        if edges.len() > 12 {
            continue;
        }
        let g = gt(n, &edges, false);
        let kruskal = s(&format!("graph_kruskal({g})"));
        let prim = s(&format!("graph_prim({g}, 0)"));
        match brute_force_mst_weight(n, &edges) {
            | Some(w) => {
                assert_eq!(nested(&kruskal).len(), n - 1);
                assert_eq!(tree_weight(&kruskal), w, "{edges:?}");
                assert_eq!(nested(&prim).len(), n - 1);
                assert_eq!(tree_weight(&prim), w, "{edges:?}");
            },
            | None => {
                // disconnected: a spanning forest has n - components edges
                let comps = union_find_components(n, &edges, None, None);
                assert_eq!(nested(&kruskal).len(), n - comps);
            },
        }
    }
    // symbolic weights are ordered last
    assert_eq!(s("graph_kruskal(graph(3, list(list(0, 1, x), list(1, 2, 4), list(0, 2, 3))))"), "list(list(0, 2, 3), list(1, 2, 4))");
}

// ---------------------------------------------------------------- flows

fn brute_force_min_cut(
    n: usize,
    edges: &[Edge],
    directed: bool,
    s: usize,
    t: usize,
) -> i64 {
    let cap = adjacency(n, edges, directed);
    let mut best = i64::MAX;
    for mask in 0u32..(1 << n) {
        if mask >> s & 1 == 0 || mask >> t & 1 == 1 {
            continue;
        }
        let mut c = 0;
        #[allow(clippy::needless_range_loop)] // index is used for more than one array / arithmetic; iterator form would not be clearer
        for u in 0..n {
            #[allow(clippy::needless_range_loop)] // index is used for more than one array / arithmetic; iterator form would not be clearer
            for v in 0..n {
                if mask >> u & 1 == 1 && mask >> v & 1 == 0 {
                    c += cap[u][v];
                }
            }
        }
        best = best.min(c);
    }
    best
}

#[test]
fn maximum_flow_equals_minimum_cut() {
    // the CLRS network
    let g = gt(
        6,
        &[(0, 1, 16), (0, 2, 13), (1, 2, 10), (2, 1, 4), (1, 3, 12), (3, 2, 9), (2, 4, 14), (4, 3, 7), (3, 5, 20), (4, 5, 4)],
        true,
    );
    assert_eq!(s(&format!("graph_edmonds_karp({g}, 0, 5)")), "23");
    assert_eq!(s(&format!("graph_dinic({g}, 0, 5)")), "23");
    let mut rng = Lcg(5);
    for round in 0..40 {
        let directed = round % 2 == 0;
        let n = 3 + rng.next(4) as usize;
        let edges = random_graph(&mut rng, n, 55, directed, 6);
        let want = brute_force_min_cut(n, &edges, directed, 0, n - 1);
        let g = gt(n, &edges, directed);
        assert_eq!(s(&format!("graph_edmonds_karp({g}, 0, {})", n - 1)), want.to_string(), "{edges:?} {directed}");
        assert_eq!(s(&format!("graph_dinic({g}, 0, {})", n - 1)), want.to_string(), "{edges:?} {directed}");
    }
    // rational and float capacities
    assert_eq!(s("graph_dinic(digraph(3, list(list(0, 1, 1/2), list(1, 2, 1/3))), 0, 2)"), "1/3");
    let f: f64 = s("graph_edmonds_karp(digraph(3, list(list(0, 1, 0.5), list(1, 2, 0.25))), 0, 2)").parse().unwrap_or(-1.0);
    assert!((f - 0.25).abs() < 1e-12);
    // source equal to sink, symbolic capacities: no answer
    assert_eq!(s("graph_dinic(digraph(2, list(list(0, 1, c))), 0, 1)"), "graph_dinic(digraph(2, list(list(0, 1, c))), 0, 1)");
    assert_eq!(s("graph_dinic(digraph(2, list(list(0, 1, 3))), 0, 0)"), "graph_dinic(digraph(2, list(list(0, 1, 3))), 0, 0)");
}

/// Maximum flow value and its minimum cost by enumerating every integral flow.
fn brute_force_min_cost_flow(
    n: usize,
    edges: &[(usize, usize, i64, i64)],
) -> (i64, i64) {
    fn go(
        i: usize,
        n: usize,
        edges: &[(usize, usize, i64, i64)],
        f: &mut Vec<i64>,
        best: &mut (i64, i64),
    ) {
        if i == edges.len() {
            let mut balance = vec![0_i64; n];
            for (k, &(u, v, _, _)) in edges.iter().enumerate() {
                balance[u] -= f[k];
                balance[v] += f[k];
            }
            if (1..n - 1).all(|x| balance[x] == 0) {
                let value = -balance[0];
                let cost: i64 = edges.iter().enumerate().map(|(k, e)| f[k] * e.3).sum();
                if value > best.0 || (value == best.0 && cost < best.1) {
                    *best = (value, cost);
                }
            }
            return;
        }
        for x in 0..=edges[i].2 {
            f[i] = x;
            go(i + 1, n, edges, f, best);
        }
    }
    let mut best = (0_i64, 0_i64);
    let m = edges.len();
    let mut f = vec![0_i64; m];
    go(0, n, edges, &mut f, &mut best);
    best
}

#[test]
fn minimum_cost_maximum_flow() {
    let g = "digraph(4, list(list(0, 1, list(2, 1)), list(0, 2, list(1, 2)), list(1, 3, list(1, 3)), list(2, 3, list(2, 1)), list(1, 2, list(1, 1))))";
    assert_eq!(s(&format!("graph_min_cost_flow({g}, 0, 3)")), "list(3, 10)");
    let mut rng = Lcg(6);
    for _ in 0..25 {
        let n = 4;
        let mut edges = Vec::new();
        for u in 0..n {
            for v in 0..n {
                if u != v && edges.len() < 6 && rng.next(100) < 40 {
                    edges.push((u, v, 1 + rng.next(2) as i64, rng.next(5) as i64));
                }
            }
        }
        let (flow, cost) = brute_force_min_cost_flow(n, &edges);
        let list: Vec<String> = edges.iter().map(|&(u, v, c, w)| format!("list({u}, {v}, list({c}, {w}))")).collect();
        let got = s(&format!("graph_min_cost_flow(digraph({n}, list({})), 0, {})", list.join(", "), n - 1));
        assert_eq!(got, format!("list({flow}, {cost})"), "{edges:?}");
    }
    // plain capacities are not a (capacity, cost) pair
    assert_eq!(s("graph_min_cost_flow(digraph(2, list(list(0, 1, 5))), 0, 1)"), "graph_min_cost_flow(digraph(2, list(list(0, 1, 5))), 0, 1)");
}

// ---------------------------------------------------------------- shortest paths

fn all_pairs_brute_force(
    n: usize,
    edges: &[Edge],
    directed: bool,
) -> Vec<Vec<Option<i64>>> {
    fn walk(
        u: usize,
        target: usize,
        adj: &[Vec<(usize, i64)>],
        seen: &mut Vec<bool>,
        cost: i64,
        best: &mut Option<i64>,
    ) {
        if u == target {
            *best = Some(best.map_or(cost, |b| b.min(cost)));
            return;
        }
        seen[u] = true;
        for &(v, w) in &adj[u] {
            if !seen[v] {
                walk(v, target, adj, seen, cost + w, best);
            }
        }
        seen[u] = false;
    }
    let mut adj = vec![Vec::new(); n];
    for &(u, v, w) in edges {
        adj[u].push((v, w));
        if !directed {
            adj[v].push((u, w));
        }
    }
    (0..n)
        .map(|a| {
            (0..n)
                .map(|b| {
                    let mut best = None;
                    walk(a, b, &adj, &mut vec![false; n], 0, &mut best);
                    best
                })
                .collect()
        })
        .collect()
}

fn dist_pairs(text: &str) -> Vec<Option<i64>> {
    nested(&text.replace("oo", "9999999")).iter().map(|r| (r[0] != 9_999_999).then_some(r[0])).collect()
}

#[test]
fn shortest_paths_against_brute_force() {
    let g = gt(5, &[(0, 1, 4), (0, 2, 1), (2, 1, 2), (1, 3, 1), (2, 3, 5)], true);
    assert_eq!(s(&format!("graph_dijkstra({g}, 0)")), "list(list(0, -1), list(3, 2), list(1, 0), list(4, 1), list(oo, -1))");
    assert_eq!(s(&format!("graph_bellman_ford({g}, 0)")), "list(list(0, -1), list(3, 2), list(1, 0), list(4, 1), list(oo, -1))");
    assert_eq!(
        s(&format!("graph_floyd_warshall({g})")),
        "list(list(0, 3, 1, 4, oo), list(oo, 0, oo, 1, oo), list(oo, 2, 0, 3, oo), list(oo, oo, oo, 0, oo), list(oo, oo, oo, oo, 0))"
    );
    let mut rng = Lcg(7);
    for round in 0..40 {
        let directed = round % 2 == 0;
        let n = 2 + rng.next(5) as usize;
        let edges = random_graph(&mut rng, n, 45, directed, 9);
        let want = all_pairs_brute_force(n, &edges, directed);
        let g = gt(n, &edges, directed);
        let fw = matrix_of(&s(&format!("graph_floyd_warshall({g})")));
        for a in 0..n {
            for b in 0..n {
                let got = (fw[a][b] != 9_999_999).then_some(fw[a][b]);
                assert_eq!(got, want[a][b], "FW {a}->{b} {edges:?} {directed}");
            }
            assert_eq!(dist_pairs(&s(&format!("graph_dijkstra({g}, {a})"))), want[a], "Dijkstra {a} {edges:?} {directed}");
            assert_eq!(dist_pairs(&s(&format!("graph_bellman_ford({g}, {a})"))), want[a], "BF {a} {edges:?} {directed}");
        }
    }
}

#[test]
fn negative_weights() {
    // a DAG with negative edges: Bellman-Ford and Floyd-Warshall agree with brute force
    let mut rng = Lcg(8);
    for _ in 0..30 {
        let n = 3 + rng.next(4) as usize;
        let mut edges = Vec::new();
        for u in 0..n {
            for v in (u + 1)..n {
                if rng.next(100) < 55 {
                    edges.push((u, v, rng.next(9) as i64 - 4));
                }
            }
        }
        let want = all_pairs_brute_force(n, &edges, true);
        let g = gt(n, &edges, true);
        let fw = matrix_of(&s(&format!("graph_floyd_warshall({g})")));
        for a in 0..n {
            assert_eq!(dist_pairs(&s(&format!("graph_bellman_ford({g}, {a})"))), want[a], "{edges:?}");
            for b in 0..n {
                assert_eq!((fw[a][b] != 9_999_999).then_some(fw[a][b]), want[a][b]);
            }
        }
    }
    // a negative cycle is reported
    assert_eq!(s("graph_bellman_ford(digraph(3, list(list(0, 1, 1), list(1, 2, -3), list(2, 1, 1))), 0)"), "false");
    // an undirected negative edge is a negative cycle
    assert_eq!(s("graph_bellman_ford(graph(2, list(list(0, 1, -1))), 0)"), "false");
}

#[test]
fn shortest_paths_with_symbolic_and_rational_weights() {
    assert_eq!(s("graph_dijkstra(digraph(3, list(list(0, 1, 1/2), list(1, 2, 1/3))), 0)"), "list(list(0, -1), list(1/2, 0), list(5/6, 1))");
    // a symbolic weight cannot be compared: the edge is not used (legacy behaviour)
    assert_eq!(s("graph_dijkstra(digraph(3, list(list(0, 1, a), list(0, 2, 2))), 0)"), "list(list(0, -1), list(oo, -1), list(2, 0))");
    // ... but a direct symbolic edge stays in the Floyd-Warshall matrix
    assert_eq!(s("graph_floyd_warshall(digraph(2, list(list(0, 1, a))))"), "list(list(0, a), list(oo, 0))");
    assert_eq!(s("graph_floyd_warshall(graph_empty(0))"), "list()");
}

#[test]
#[allow(clippy::cast_sign_loss)] // operand is non-negative by construction (index/count)
fn unweighted_shortest_paths() {
    let g = gt(5, &[(0, 1, 1), (1, 2, 1), (0, 3, 1), (3, 2, 1)], false);
    assert_eq!(s(&format!("graph_shortest_path_unweighted({g}, 0)")), "list(list(0, 0, -1), list(1, 1, 0), list(3, 1, 0), list(2, 2, 1))");
    let mut rng = Lcg(9);
    for _ in 0..25 {
        let n = 2 + rng.next(6) as usize;
        let edges = random_graph(&mut rng, n, 35, false, 1);
        let want = all_pairs_brute_force(n, &edges, false);
        let rows = nested(&s(&format!("graph_shortest_path_unweighted({}, 0)", gt(n, &edges, false))));
        let reachable = want[0].iter().filter(|d| d.is_some()).count();
        assert_eq!(rows.len(), reachable);
        for r in rows {
            assert_eq!(want[0][r[0] as usize], Some(r[1]));
        }
    }
}

// ---------------------------------------------------------------- matchings

fn brute_force_max_matching(
    n: usize,
    edges: &[Edge],
) -> usize {
    let m = edges.len();
    let mut best = 0;
    for mask in 0u32..(1 << m) {
        let mut used = vec![false; n];
        let mut ok = true;
        #[allow(clippy::needless_range_loop)] // index is used for more than one array / arithmetic; iterator form would not be clearer
        for i in 0..m {
            if mask >> i & 1 == 1 {
                let (u, v, _) = edges[i];
                if used[u] || used[v] {
                    ok = false;
                    break;
                }
                used[u] = true;
                used[v] = true;
            }
        }
        if ok {
            best = best.max(mask.count_ones() as usize);
        }
    }
    best
}

#[allow(clippy::cast_sign_loss)] // operand is non-negative by construction (index/count)
fn is_matching(
    pairs: &[Vec<i64>],
    edges: &[Edge],
) -> bool {
    let mut used = std::collections::HashSet::new();
    pairs.iter().all(|p| {
        let (u, v) = (p[0] as usize, p[1] as usize);
        used.insert(u) && used.insert(v) && edges.iter().any(|&(a, b, _)| (a, b) == (u, v) || (a, b) == (v, u))
    })
}

fn random_bipartite(
    rng: &mut Lcg,
    left: usize,
    right: usize,
) -> (Vec<Edge>, String) {
    let mut edges = Vec::new();
    for u in 0..left {
        for v in 0..right {
            if rng.next(100) < 45 {
                edges.push((u, left + v, 1));
            }
        }
    }
    let part: Vec<String> = (0..left + right).map(|i| usize::from(i >= left).to_string()).collect();
    (edges, format!("list({})", part.join(", ")))
}

#[test]
fn bipartite_matching_and_vertex_cover() {
    assert_eq!(s("graph_is_bipartite(graph_cycle(4))"), "list(0, 1, 0, 1)");
    assert_eq!(s("graph_is_bipartite(graph_cycle(5))"), "false");
    assert_eq!(s("graph_is_bipartite(graph(2, list(list(1, 1))))"), "false");
    assert_eq!(s("graph_is_bipartite(graph_empty(3))"), "list(0, 0, 0)");
    let mut rng = Lcg(10);
    for _ in 0..40 {
        let (left, right) = (1 + rng.next(4) as usize, 1 + rng.next(4) as usize);
        let (edges, part) = random_bipartite(&mut rng, left, right);
        if edges.len() > 12 {
            continue;
        }
        let n = left + right;
        let g = gt(n, &edges, false);
        let want = brute_force_max_matching(n, &edges);
        let kuhn = s(&format!("graph_bipartite_matching({g}, {part})"));
        let hk = s(&format!("graph_hopcroft_karp({g}, {part})"));
        let blossom = s(&format!("graph_blossom({g})"));
        for (name, text) in [("kuhn", &kuhn), ("hopcroft-karp", &hk), ("blossom", &blossom)] {
            let pairs = nested(text);
            assert_eq!(pairs.len(), want, "{name} {edges:?}");
            assert!(is_matching(&pairs, &edges), "{name} {edges:?}");
        }
        // Konig: a minimum vertex cover has the size of a maximum matching and covers every edge
        let cover = nums(&s(&format!("graph_vertex_cover({g}, {part}, {kuhn})")));
        assert_eq!(cover.len(), want, "{edges:?}");
        assert!(edges.iter().all(|&(u, v, _)| cover.contains(&(u as i64)) || cover.contains(&(v as i64))), "{edges:?}");
    }
    // a malformed partition stays
    assert_eq!(s("graph_hopcroft_karp(graph_path(3), list(0, 1))"), "graph_hopcroft_karp(graph(3, list(list(0, 1), list(1, 2))), list(0, 1))");
}

#[test]
fn blossom_on_general_graphs() {
    // an odd cycle with pendant edges needs a blossom contraction
    let g = gt(7, &[(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 0, 1), (2, 5, 1), (4, 6, 1)], false);
    assert_eq!(nested(&s(&format!("graph_blossom({g})"))).len(), 3);
    let petersen: Vec<Edge> = vec![
        (0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 0, 1),
        (0, 5, 1), (1, 6, 1), (2, 7, 1), (3, 8, 1), (4, 9, 1),
        (5, 7, 1), (7, 9, 1), (9, 6, 1), (6, 8, 1), (8, 5, 1),
    ];
    assert_eq!(nested(&s(&format!("graph_blossom({})", gt(10, &petersen, false)))).len(), 5);
    let mut rng = Lcg(11);
    for _ in 0..60 {
        let n = 3 + rng.next(5) as usize;
        let edges = random_graph(&mut rng, n, 40, false, 1);
        if edges.len() > 13 {
            continue;
        }
        let pairs = nested(&s(&format!("graph_blossom({})", gt(n, &edges, false))));
        assert_eq!(pairs.len(), brute_force_max_matching(n, &edges), "{edges:?}");
        assert!(is_matching(&pairs, &edges), "{edges:?}");
    }
}

// ---------------------------------------------------------------- ordering

#[test]
fn topological_sorts() {
    let g = "digraph(6, list(list(5, 2), list(5, 0), list(4, 0), list(4, 1), list(2, 3), list(3, 1)))";
    assert_eq!(s(&format!("graph_toposort_kahn({g})")), "list(4, 5, 2, 0, 3, 1)");
    assert_eq!(s(&format!("graph_toposort({g})")), "list(4, 5, 2, 0, 3, 1)");
    assert_eq!(s(&format!("graph_toposort_dfs({g})")), "list(5, 4, 2, 3, 1, 0)");
    assert_eq!(s("graph_toposort(digraph(2, list(list(0, 1), list(1, 0))))"), "false");
    assert_eq!(s("graph_toposort_dfs(digraph(2, list(list(0, 1), list(1, 0))))"), "false");
    assert_eq!(s("graph_toposort(graph_path(3))"), "false");
    let mut rng = Lcg(12);
    for _ in 0..30 {
        let n = 2 + rng.next(7) as usize;
        let mut edges = Vec::new();
        // a random DAG on a shuffled order
        let mut order: Vec<usize> = (0..n).collect();
        for i in (1..n).rev() {
            order.swap(i, rng.next(i as u64 + 1) as usize);
        }
        for i in 0..n {
            for j in (i + 1)..n {
                if rng.next(100) < 35 {
                    edges.push((order[i], order[j], 1));
                }
            }
        }
        for name in ["graph_toposort", "graph_toposort_dfs"] {
            let result = nums(&s(&format!("{name}({})", gt(n, &edges, true))));
            assert_eq!(result.len(), n);
            let pos: Vec<usize> = (0..n).map(|v| result.iter().position(|&x| x == v as i64).unwrap_or(0)).collect();
            assert!(edges.iter().all(|&(u, v, _)| pos[u] < pos[v]), "{name} {edges:?}");
        }
    }
}

// ---------------------------------------------------------------- spectra

#[test]
fn spectral_analysis_and_algebraic_connectivity() {
    let out = s("spectral_analysis(list(list(2, 0), list(0, 3)))");
    assert!(out.starts_with("list(list(2, 3)"), "{out}");
    let value = |src: &str| s(src).parse::<f64>().unwrap_or(f64::NAN);
    for n in 3..7_usize {
        let nf = n as f64;
        assert!((value(&format!("algebraic_connectivity(graph_complete({n}))")) - nf).abs() < 1e-9);
        let path = 2.0 * (1.0 - (std::f64::consts::PI / nf).cos());
        assert!((value(&format!("algebraic_connectivity(graph_path({n}))")) - path).abs() < 1e-9);
        let cycle = 2.0 * (1.0 - (2.0 * std::f64::consts::PI / nf).cos());
        assert!((value(&format!("algebraic_connectivity(graph_cycle({n}))")) - cycle).abs() < 1e-9);
    }
    // a disconnected graph has algebraic connectivity 0
    assert!(value("algebraic_connectivity(graph_empty(3))").abs() < 1e-12);
    assert_eq!(s("algebraic_connectivity(graph_empty(1))"), "algebraic_connectivity(graph(1, list()))");
    assert_eq!(s("algebraic_connectivity(graph(2, list(list(0, 1, w))))"), "algebraic_connectivity(graph(2, list(list(0, 1, w))))");
}

// ---------------------------------------------------------------- isomorphism and colouring

fn permute(
    edges: &[Edge],
    perm: &[usize],
) -> Vec<Edge> {
    edges.iter().map(|&(u, v, w)| (perm[u], perm[v], w)).collect()
}

fn brute_force_isomorphic(
    n: usize,
    a: &[Edge],
    b: &[Edge],
) -> bool {
    fn next_permutation(p: &mut [usize]) -> bool {
        let n = p.len();
        let Some(i) = (1..n).rev().find(|&i| p[i - 1] < p[i]) else { return false };
        let j = (i..n).rev().find(|&j| p[j] > p[i - 1]).unwrap_or(i);
        p.swap(i - 1, j);
        p[i..].reverse();
        true
    }
    let (ma, mb) = (adjacency(n, a, false), adjacency(n, b, false));
    let mut perm: Vec<usize> = (0..n).collect();
    loop {
        if (0..n).all(|u| (0..n).all(|v| ma[u][v] == mb[perm[u]][perm[v]])) {
            return true;
        }
        if !next_permutation(&mut perm) {
            return false;
        }
    }
}

#[test]
fn isomorphism() {
    let mut rng = Lcg(13);
    for _ in 0..30 {
        let n = 3 + rng.next(4) as usize;
        let a = random_graph(&mut rng, n, 45, false, 1);
        let b = random_graph(&mut rng, n, 45, false, 1);
        let mut perm: Vec<usize> = (0..n).collect();
        for i in (1..n).rev() {
            perm.swap(i, rng.next(i as u64 + 1) as usize);
        }
        let ga = gt(n, &a, false);
        let permuted = gt(n, &permute(&a, &perm), false);
        assert_eq!(s(&format!("graph_isomorphic({ga}, {permuted})")), "true", "{a:?}");
        assert_eq!(s(&format!("graph_isomorphic_heuristic({ga}, {permuted})")), "true", "{a:?}");
        let gb = gt(n, &b, false);
        let truth = brute_force_isomorphic(n, &a, &b);
        assert_eq!(s(&format!("graph_isomorphic({ga}, {gb})")), truth.to_string(), "{a:?} {b:?}");
        if truth {
            assert_eq!(s(&format!("graph_isomorphic_heuristic({ga}, {gb})")), "true");
        }
    }
    // the Weisfeiler-Lehman test cannot tell C6 from two triangles; the exact test can
    let c6 = "graph_cycle(6)";
    let two_triangles = "graph_disjoint_union(graph_cycle(3), graph_cycle(3))";
    assert_eq!(s(&format!("graph_isomorphic_heuristic({c6}, {two_triangles})")), "true");
    assert_eq!(s(&format!("graph_isomorphic({c6}, {two_triangles})")), "false");
    assert_eq!(s("graph_isomorphic(graph_path(4), graph_cycle(4))"), "false");
    assert_eq!(s("graph_isomorphic_heuristic(graph_path(4), graph_cycle(4))"), "false");
    assert_eq!(s("graph_isomorphic(graph_path(3), graph_path(4))"), "false");
    // directed graphs: a path and its reverse are isomorphic, a path and a star are not
    assert_eq!(s("graph_isomorphic(digraph(3, list(list(0, 1), list(1, 2))), digraph(3, list(list(2, 1), list(1, 0))))"), "true");
    assert_eq!(s("graph_isomorphic(digraph(3, list(list(0, 1), list(1, 2))), digraph(3, list(list(0, 1), list(0, 2))))"), "false");
    // a long cycle no longer blows up (the legacy colour strings did)
    assert_eq!(s("graph_isomorphic_heuristic(graph_cycle(30), graph_cycle(30))"), "true");
}

fn brute_force_chromatic(
    n: usize,
    edges: &[Edge],
) -> usize {
    fn go(
        v: usize,
        k: usize,
        n: usize,
        edges: &[Edge],
        colors: &mut Vec<usize>,
    ) -> bool {
        if v == n {
            return true;
        }
        for c in 0..k {
            if edges.iter().all(|&(a, b, _)| !((a == v && b < v && colors[b] == c) || (b == v && a < v && colors[a] == c))) {
                colors[v] = c;
                if go(v + 1, k, n, edges, colors) {
                    return true;
                }
            }
        }
        false
    }
    (1..=n.max(1)).find(|&k| go(0, k, n, edges, &mut vec![0; n])).unwrap_or(n)
}

#[test]
#[allow(clippy::cast_sign_loss)] // operand is non-negative by construction (index/count)
fn colouring() {
    for (g, chi) in [
        ("graph_complete(5)", "5"),
        ("graph_cycle(5)", "3"),
        ("graph_cycle(6)", "2"),
        ("graph_path(4)", "2"),
        ("graph_empty(3)", "1"),
        ("graph_empty(0)", "0"),
    ] {
        assert_eq!(s(&format!("graph_chromatic_number({g})")), chi, "{g}");
    }
    let petersen = gt(
        10,
        &[
            (0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 0, 1),
            (0, 5, 1), (1, 6, 1), (2, 7, 1), (3, 8, 1), (4, 9, 1),
            (5, 7, 1), (7, 9, 1), (9, 6, 1), (6, 8, 1), (8, 5, 1),
        ],
        false,
    );
    assert_eq!(s(&format!("graph_chromatic_number({petersen})")), "3");
    // a self-loop has no proper colouring
    assert_eq!(s("graph_chromatic_number(graph(1, list(list(0, 0))))"), "graph_chromatic_number(graph(1, list(list(0, 0))))");
    let mut rng = Lcg(14);
    for _ in 0..30 {
        let n = 2 + rng.next(6) as usize;
        let edges = random_graph(&mut rng, n, 45, false, 1);
        let g = gt(n, &edges, false);
        let chi = brute_force_chromatic(n, &edges);
        assert_eq!(s(&format!("graph_chromatic_number({g})")), chi.to_string(), "{edges:?}");
        // greedy colouring is proper and uses at most max degree + 1 colours
        let colors = nums(&s(&format!("graph_greedy_coloring({g})")));
        assert_eq!(colors.len(), n);
        assert!(edges.iter().all(|&(u, v, _)| colors[u] != colors[v]), "{edges:?} {colors:?}");
        let used = colors.iter().max().map_or(0, |m| m + 1) as usize;
        assert!(used >= chi);
        let max_degree = (0..n).map(|v| edges.iter().filter(|e| e.0 == v || e.1 == v).count()).max().unwrap_or(0);
        assert!(used <= max_degree + 1);
    }
}

// ---------------------------------------------------------------- graph operations

fn kron(
    a: &[Vec<i64>],
    b: &[Vec<i64>],
) -> Vec<Vec<i64>> {
    let (na, nb) = (a.len(), b.len());
    let mut m = vec![vec![0; na * nb]; na * nb];
    for i in 0..na {
        for j in 0..na {
            for k in 0..nb {
                for l in 0..nb {
                    m[i * nb + k][j * nb + l] = a[i][j] * b[k][l];
                }
            }
        }
    }
    m
}

fn identity(n: usize) -> Vec<Vec<i64>> {
    (0..n).map(|i| (0..n).map(|j| i64::from(i == j)).collect()).collect()
}

fn add(
    a: &[Vec<i64>],
    b: &[Vec<i64>],
) -> Vec<Vec<i64>> {
    a.iter().zip(b).map(|(r, q)| r.iter().zip(q).map(|(x, y)| x + y).collect()).collect()
}

#[test]
fn graph_operations() {
    assert_eq!(s("graph_complement(graph_path(3))"), "graph(3, list(list(0, 2)))");
    assert_eq!(s("graph_complement(graph_complete(4))"), "graph(4, list())");
    assert_eq!(s("graph_complement(digraph(2, list(list(0, 1))))"), "digraph(2, list(list(1, 0)))");
    assert_eq!(s("induced_subgraph(graph_complete(4), list(3, 1, 2))"), "graph(3, list(list(1, 2), list(1, 0), list(2, 0)))");
    assert_eq!(s("induced_subgraph(graph(3, list(list(0, 1, w), list(1, 2))), list(1, 0))"), "graph(2, list(list(1, 0, w)))");
    assert_eq!(s("graph_union(graph_path(3), graph(4, list(list(1, 2), list(2, 3))))"), "graph(4, list(list(0, 1), list(1, 2), list(2, 3)))");
    assert_eq!(s("graph_intersection(graph_complete(3), graph(4, list(list(0, 1), list(2, 3), list(1, 2, 7))))"), "graph(3, list(list(0, 1)))");
    assert_eq!(s("graph_disjoint_union(graph_path(2), graph_path(3))"), "graph(5, list(list(0, 1), list(2, 3), list(3, 4)))");
    assert_eq!(s("graph_join(graph_empty(2), graph_empty(1))"), "graph(3, list(list(0, 2), list(1, 2)))");
    assert_eq!(s("graph_cartesian(graph_path(2), graph_path(2))"), "graph(4, list(list(0, 2), list(1, 3), list(0, 1), list(2, 3)))");
    assert_eq!(s("graph_tensor(graph(2, list(list(0, 1, a))), graph(2, list(list(0, 1, b))))"), "graph(4, list(list(0, 3, a*b), list(1, 2, a*b)))");
    // the products satisfy their defining adjacency identities
    let mut rng = Lcg(15);
    for _ in 0..12 {
        let (n1, n2) = (2 + rng.next(3) as usize, 2 + rng.next(3) as usize);
        let (e1, e2) = (random_graph(&mut rng, n1, 50, false, 4), random_graph(&mut rng, n2, 50, false, 4));
        let (g1, g2) = (gt(n1, &e1, false), gt(n2, &e2, false));
        let (a1, a2) = (adjacency(n1, &e1, false), adjacency(n2, &e2, false));
        let adj = |g: &str| matrix_of(&s(&format!("graph_adjacency({g})")));
        let cart = s(&format!("graph_cartesian({g1}, {g2})"));
        assert_eq!(adj(&cart), add(&kron(&a1, &identity(n2)), &kron(&identity(n1), &a2)));
        let tensor = s(&format!("graph_tensor({g1}, {g2})"));
        assert_eq!(adj(&tensor), kron(&a1, &a2));
        // disjoint union: block diagonal; join: all cross edges
        let du = adj(&s(&format!("graph_disjoint_union({g1}, {g2})")));
        let join = adj(&s(&format!("graph_join({g1}, {g2})")));
        for u in 0..n1 + n2 {
            for v in 0..n1 + n2 {
                let (in1, in2) = (u < n1, v < n1);
                let expected_du = if in1 && in2 { a1[u][v] } else if !in1 && !in2 { a2[u - n1][v - n1] } else { 0 };
                assert_eq!(du[u][v], expected_du);
                let expected_join = if in1 == in2 { expected_du } else { 1 };
                assert_eq!(join[u][v], expected_join);
            }
        }
        // complement of the complement of a simple unweighted graph is itself
        let simple: Vec<Edge> = e1.iter().map(|&(u, v, _)| (u, v, 1)).collect();
        let g = gt(n1, &simple, false);
        let twice = s(&format!("graph_complement(graph_complement({g}))"));
        let (n, edges) = parse_graph(&twice);
        assert_eq!((n, edge_set(&edges, false)), (n1, edge_set(&simple, false)));
        // union and intersection of unweighted graphs are the set operations
        let other = random_graph(&mut rng, n1, 50, false, 1);
        let h = gt(n1, &other, false);
        let un: Vec<Edge> = {
            let mut all = edge_set(&simple, false);
            all.extend(edge_set(&other, false));
            all.sort_unstable();
            all.dedup();
            all
        };
        let (_, got) = parse_graph(&s(&format!("graph_union({g}, {h})")));
        assert_eq!(edge_set(&got, false), un);
        let inter: Vec<Edge> = edge_set(&simple, false).into_iter().filter(|e| edge_set(&other, false).contains(e)).collect();
        let (_, got) = parse_graph(&s(&format!("graph_intersection({g}, {h})")));
        assert_eq!(edge_set(&got, false), inter);
    }
    // directed operations
    assert_eq!(
        s("graph_join(digraph(1, list()), digraph(1, list()))"),
        "digraph(2, list(list(0, 1), list(1, 0)))"
    );
    assert_eq!(s("graph_cartesian(graph_path(2), digraph(2, list(list(0, 1))))"), "digraph(4, list(list(0, 2), list(1, 3), list(0, 1), list(2, 3)))");
}
