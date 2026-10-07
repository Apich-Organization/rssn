//! Graphs as terms, and the classical graph algorithms.
//!
//! # Representation
//!
//! A graph is `graph(n, edges)` (undirected) or `digraph(n, edges)`
//! (directed) with the vertices `0, ..., n - 1` and
//! `edges = list(list(u, v), list(u, v, w), ...)`: the weight `w` may be any
//! term (a number, a symbol, `list(capacity, cost)` for min-cost flow) and
//! defaults to `1`. An optional third argument holds hyperedges
//! `list(list(list(v, ...), w), ...)`. Parallel edges and self-loops are
//! allowed. Vertices are always integers; the products of two graphs number
//! the pair `(u, v)` as `u * n2 + v`.
//!
//! Results that are graphs come back in this form (a weight of `1` is
//! omitted); vertex sets, paths and matchings are lists of integers.
//! Matrices are `list(list(...))` terms that compose with the linear algebra
//! rules, with symbolic weights kept as they are:
//! `det(graph_adjacency(graph(3, list(list(0, 1, a), list(1, 2, b), list(0, 2, c)))))`.
//! "No result" answers are `false`.
//!
//! # Operators
//!
//! | operator | value |
//! |---|---|
//! | `graph_empty(n)`, `graph_complete(n)`, `graph_path(n)`, `graph_cycle(n)` | standard graphs |
//! | `graph_nodes(g)`, `graph_node_count(g)`, `graph_is_directed(g)`, `graph_node_id(g, v)` | vertices |
//! | `graph_add_node(g)`, `graph_add_edge(g, u, v[, w])`, `graph_add_hyperedge(g, vs, w)` | functional updates |
//! | `graph_neighbors(g, u)`, `graph_out_degree(g, u)`, `graph_in_degree(g, u)`, `graph_edges(g)`, `graph_hyperedges(g)` | `list(list(v, w), ...)`; degrees count adjacency entries; undirected edges are listed once (`u <= v`) as `list(u, v, w)` |
//! | `graph_adjacency(g)`, `graph_incidence(g)`, `graph_laplacian(g)` | matrices; the adjacency matrix adds up parallel edges; the incidence matrix has columns for `graph_edges` with `-1/+1` for directed and `1/1` for undirected edges; the Laplacian is `D - A` with the *weighted* degree matrix `D` |
//! | `graph_dfs(g, s)`, `graph_bfs(g, s)` | visit order |
//! | `graph_components(g)`, `graph_is_connected(g)`, `graph_scc(g)` | components in order of their first vertex (BFS order inside each of `graph_components`, sorted inside each of `graph_scc`) |
//! | `graph_has_cycle(g)`, `graph_bridges(g)` | cycle test; `list(bridges, articulation_points)` of an undirected graph |
//! | `graph_kruskal(g)`, `graph_prim(g, s)` | minimum spanning forest edges `list(u, v, w)` |
//! | `graph_edmonds_karp(g, s, t)`, `graph_dinic(g, s, t)` | maximum flow, edge weights are capacities (exact for exact capacities) |
//! | `graph_min_cost_flow(g, s, t)` | `list(flow, cost)` of a minimum-cost maximum flow, weights are `list(capacity, cost)` |
//! | `graph_dijkstra(g, s)`, `graph_bellman_ford(g, s)` | `list(list(dist, prev), ...)` per vertex (`oo`, `-1` when unreachable); Bellman-Ford gives `false` on a negative cycle |
//! | `graph_floyd_warshall(g)`, `graph_shortest_path_unweighted(g, s)` | distance matrix; `list(list(v, dist, prev), ...)` of the reachable vertices |
//! | `graph_is_bipartite(g)` | the 0/1 colouring or `false` |
//! | `graph_bipartite_matching(g, part)`, `graph_hopcroft_karp(g, part)`, `graph_vertex_cover(g, part, matching)`, `graph_blossom(g)` | maximum matchings as `list(u, v)`, minimum vertex cover (Konig), maximum matching of a general graph |
//! | `graph_toposort(g)`, `graph_toposort_kahn(g)`, `graph_toposort_dfs(g)` | a topological order of a digraph or `false` |
//! | `spectral_analysis(M)`, `algebraic_connectivity(g)` | `list(eigenvals(M), eigenvects(M))`; second smallest Laplacian eigenvalue (float, numeric weights) |
//! | `graph_isomorphic_heuristic(g, h)`, `graph_isomorphic(g, h)` | Weisfeiler-Lehman test (necessary condition); exact test by backtracking (up to 64 vertices) |
//! | `graph_greedy_coloring(g)`, `graph_chromatic_number(g)` | Welsh-Powell colours `0, 1, ...` per vertex; exact chromatic number (up to 40 vertices, no self-loops) |
//! | `induced_subgraph(g, vs)`, `graph_union`, `graph_intersection`, `graph_cartesian`, `graph_tensor`, `graph_complement`, `graph_disjoint_union`, `graph_join` | graph operations |
//!
//! Shortest-path algorithms order by the numeric value of the weights; a
//! symbolic weight counts as `oo` for the comparison (an edge with a symbolic
//! weight is never part of a computed path), as in the legacy code.

#[cfg(test)]
mod tests;

use std::collections::VecDeque;

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::Signed;
use num_traits::ToPrimitive;
use num_traits::Zero;

use super::apply;
use super::def;
use super::def_inert;
use super::def_request;
use super::idx;
use super::idxs;
use super::items;
use super::matrix;
use super::neg;
use super::prod;
use super::rational;
use super::sum;
use super::V;
use crate::graph::op::core;
use crate::graph::rule::Installer;
use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::RuleError;

/// One adjacency entry.
#[derive(Clone, Copy)]
struct Adj {
    to: usize,
    w: NodeId,
    id: usize,
}

/// A parsed graph term.
struct Gr {
    n: usize,
    directed: bool,
    edges: Vec<(usize, usize, NodeId)>,
    hyper: Vec<(Vec<usize>, NodeId)>,
    adj: Vec<Vec<Adj>>,
    radj: Vec<Vec<Adj>>,
}

impl Gr {
    fn new(
        n: usize,
        directed: bool,
    ) -> Self {
        Self {
            n,
            directed,
            edges: Vec::new(),
            hyper: Vec::new(),
            adj: vec![Vec::new(); n],
            radj: vec![Vec::new(); n],
        }
    }

    fn add_edge(
        &mut self,
        u: usize,
        v: usize,
        w: NodeId,
    ) {
        let id = self.edges.len();
        self.edges.push((u, v, w));
        self.adj[u].push(Adj { to: v, w, id });
        self.radj[v].push(Adj { to: u, w, id });
        if !self.directed {
            self.adj[v].push(Adj { to: u, w, id });
            self.radj[u].push(Adj { to: v, w, id });
        }
    }

    /// Each edge once: for an undirected graph `u <= v`.
    fn get_edges(&self) -> Vec<(usize, usize, NodeId)> {
        let mut out = Vec::new();
        let mut seen_loop = vec![false; self.edges.len()];
        for (u, list) in self.adj.iter().enumerate() {
            for a in list {
                if !self.directed && u > a.to {
                    continue;
                }
                if u == a.to && !self.directed {
                    if seen_loop[a.id] {
                        continue;
                    }
                    seen_loop[a.id] = true;
                }
                out.push((u, a.to, a.w));
            }
        }
        out
    }

    fn read(
        g: &Graph,
        node: NodeId,
        one: NodeId,
    ) -> Option<Self> {
        let graph_op = g.ops().lookup("graph")?;
        let digraph_op = g.ops().lookup("digraph")?;
        let is_graph = |n: NodeId| g.op(n) == graph_op || g.op(n) == digraph_op;
        let found = if is_graph(node) {
            node
        } else {
            g.enodes(g.find(node)).find(|&n| is_graph(n))?
        };
        let directed = g.op(found) == digraph_op;
        let children = g.children(found);
        if children.len() < 2 || children.len() > 3 {
            return None;
        }
        let n = idx(g, children[0])?;
        if n > 100_000 {
            return None;
        }
        let mut gr = Self::new(n, directed);
        for e in items(g, children[1])? {
            let parts = items(g, e)?;
            let (u, v) = (idx(g, *parts.first()?)?, idx(g, *parts.get(1)?)?);
            if u >= n || v >= n || parts.len() > 3 {
                return None;
            }
            gr.add_edge(u, v, parts.get(2).copied().unwrap_or(one));
        }
        if let Some(&h) = children.get(2) {
            for e in items(g, h)? {
                let parts = items(g, e)?;
                let [vs, w] = parts.as_slice() else { return None };
                let vs = idxs(g, *vs)?;
                if vs.iter().any(|&v| v >= n) {
                    return None;
                }
                gr.hyper.push((vs, *w));
            }
        }
        Some(gr)
    }

    fn build(
        &self,
        g: &mut Graph,
    ) -> NodeId {
        let int = |g: &mut Graph, v: usize| g.int(i64::try_from(v).unwrap_or(0));
        let edges: Vec<NodeId> = self
            .edges
            .iter()
            .map(|&(u, v, w)| {
                let (u, v) = (int(g, u), int(g, v));
                if g.number_of(w).is_some_and(Number::is_one) {
                    g.node(core::LIST, &[u, v])
                } else {
                    g.node(core::LIST, &[u, v, w])
                }
            })
            .collect();
        let edge_list = g.node(core::LIST, &edges);
        let n = int(g, self.n);
        let mut args = vec![n, edge_list];
        if !self.hyper.is_empty() {
            let hs: Vec<NodeId> = self
                .hyper
                .iter()
                .map(|(vs, w)| {
                    let vs: Vec<NodeId> = vs.iter().map(|&v| int(g, v)).collect();
                    let vs = g.node(core::LIST, &vs);
                    g.node(core::LIST, &[vs, *w])
                })
                .collect();
            args.push(g.node(core::LIST, &hs));
        }
        let name = if self.directed { "digraph" } else { "graph" };
        let op = g.ops().lookup(name).unwrap_or(core::LIST);
        g.try_node(op, &args).unwrap_or(edge_list)
    }
}

fn read(
    cx: &mut Cx<'_>,
    node: NodeId,
) -> Option<Gr> {
    let one = cx.graph.int(1);
    Gr::read(cx.graph, node, one)
}

fn graph_value(
    cx: &mut Cx<'_>,
    gr: &Gr,
) -> V {
    V::Node(gr.build(cx.graph))
}

fn indices(v: &[usize]) -> V {
    V::List(v.iter().map(|&x| V::uint(x)).collect())
}

fn pair(
    u: usize,
    v: usize,
) -> V {
    V::List(vec![V::uint(u), V::uint(v)])
}

fn vertex(
    cx: &Cx<'_>,
    node: NodeId,
    gr: &Gr,
) -> Option<usize> {
    idx(cx.graph, node).filter(|&v| v < gr.n)
}

/// The numeric value of a weight; `+oo` for anything symbolic.
fn key(
    g: &Graph,
    w: NodeId,
) -> f64 {
    g.number_of(w).map_or(f64::INFINITY, Number::to_f64)
}

/// `oo`.
fn infinity(g: &mut Graph) -> NodeId {
    match g.ops().lookup("oo") {
        | Some(op) => g.node(op, &[]),
        | None => g.sym("oo"),
    }
}

/// A distance: the term and its numeric key (`+oo` when unknown).
#[derive(Clone, Copy)]
struct Dist {
    term: NodeId,
    key: f64,
}

fn dist_add(
    g: &mut Graph,
    a: Dist,
    b: Dist,
) -> Dist {
    Dist {
        term: sum(g, &[a.term, b.term]),
        key: a.key + b.key,
    }
}

// ---------------- construction and queries ----------------

fn standard(
    cx: &mut Cx<'_>,
    a: &[NodeId],
    shape: fn(usize) -> Vec<(usize, usize)>,
) -> Option<V> {
    let n = idx(cx.graph, *a.first()?).filter(|&n| n <= 10_000)?;
    let mut gr = Gr::new(n, false);
    let one = cx.graph.int(1);
    for (u, v) in shape(n) {
        gr.add_edge(u, v, one);
    }
    Some(graph_value(cx, &gr))
}

fn graph_nodes(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let gr = read(cx, *a.first()?)?;
    Some(V::List((0..gr.n).map(V::uint).collect()))
}

fn graph_node_id(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [g, v] = a else { return None };
    let gr = read(cx, *g)?;
    vertex(cx, *v, &gr).map(V::uint)
}

fn graph_add_node(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let mut gr = read(cx, *a.first()?)?;
    gr.n += 1;
    gr.adj.push(Vec::new());
    gr.radj.push(Vec::new());
    Some(graph_value(cx, &gr))
}

fn graph_add_edge(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    if !(3..=4).contains(&a.len()) {
        return None;
    }
    let mut gr = read(cx, a[0])?;
    let (u, v) = (vertex(cx, a[1], &gr)?, vertex(cx, a[2], &gr)?);
    let w = match a.get(3) {
        | Some(&w) => w,
        | None => cx.graph.int(1),
    };
    gr.add_edge(u, v, w);
    Some(graph_value(cx, &gr))
}

fn graph_add_hyperedge(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [g, vs, w] = a else { return None };
    let mut gr = read(cx, *g)?;
    let vs = idxs(cx.graph, *vs)?;
    if vs.iter().any(|&v| v >= gr.n) {
        return None;
    }
    gr.hyper.push((vs, *w));
    Some(graph_value(cx, &gr))
}

fn graph_hyperedges(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let gr = read(cx, *a.first()?)?;
    Some(V::List(gr.hyper.iter().map(|(vs, w)| V::List(vec![indices(vs), V::Node(*w)])).collect()))
}

fn graph_neighbors(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [g, u] = a else { return None };
    let gr = read(cx, *g)?;
    let u = vertex(cx, *u, &gr)?;
    Some(V::List(gr.adj[u].iter().map(|e| V::List(vec![V::uint(e.to), V::Node(e.w)])).collect()))
}

fn graph_degree(
    cx: &mut Cx<'_>,
    a: &[NodeId],
    incoming: bool,
) -> Option<V> {
    let [g, u] = a else { return None };
    let gr = read(cx, *g)?;
    let u = vertex(cx, *u, &gr)?;
    Some(V::uint(if incoming { gr.radj[u].len() } else { gr.adj[u].len() }))
}

fn graph_edges(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let gr = read(cx, *a.first()?)?;
    Some(V::List(gr.get_edges().into_iter().map(edge_value).collect()))
}

/// For each ordered pair the weights of the parallel edges.
fn cells(gr: &Gr) -> Vec<Vec<Vec<NodeId>>> {
    let mut cells: Vec<Vec<Vec<NodeId>>> = vec![vec![Vec::new(); gr.n]; gr.n];
    for (u, list) in gr.adj.iter().enumerate() {
        for e in list {
            cells[u][e.to].push(e.w);
        }
    }
    cells
}

fn graph_adjacency(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let gr = read(cx, *a.first()?)?;
    let zero = cx.graph.int(0);
    let rows: Vec<Vec<NodeId>> = cells(&gr)
        .iter()
        .map(|row| row.iter().map(|c| if c.is_empty() { zero } else { sum(cx.graph, c) }).collect())
        .collect();
    Some(V::Node(matrix(cx.graph, &rows)))
}

fn graph_incidence(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let gr = read(cx, *a.first()?)?;
    let edges = gr.get_edges();
    let mut m = vec![vec![0_i64; edges.len()]; gr.n];
    for (j, &(u, v, _)) in edges.iter().enumerate() {
        m[u][j] = if gr.directed { -1 } else { 1 };
        m[v][j] = 1;
    }
    Some(V::List(m.into_iter().map(|r| V::List(r.into_iter().map(V::int).collect())).collect()))
}

fn graph_laplacian(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let gr = read(cx, *a.first()?)?;
    let zero = cx.graph.int(0);
    let cells = cells(&gr);
    let mut rows = Vec::with_capacity(gr.n);
    #[allow(clippy::needless_range_loop)] // index is used for more than one array / arithmetic; iterator form would not be clearer
    for u in 0..gr.n {
        let degree: Vec<NodeId> = gr.adj[u].iter().map(|e| e.w).collect();
        let degree = sum(cx.graph, &degree);
        let mut row = Vec::with_capacity(gr.n);
        #[allow(clippy::needless_range_loop)] // index is used for more than one array / arithmetic; iterator form would not be clearer
        for v in 0..gr.n {
            let cell = if cells[u][v].is_empty() { zero } else { sum(cx.graph, &cells[u][v]) };
            let minus = neg(cx.graph, cell);
            row.push(if u == v { sum(cx.graph, &[degree, minus]) } else { minus });
        }
        rows.push(row);
    }
    Some(V::Node(matrix(cx.graph, &rows)))
}

// ---------------- traversal ----------------

fn dfs_visit(
    gr: &Gr,
    u: usize,
    seen: &mut [bool],
    out: &mut Vec<usize>,
) {
    seen[u] = true;
    out.push(u);
    for e in &gr.adj[u] {
        if !seen[e.to] {
            dfs_visit(gr, e.to, seen, out);
        }
    }
}

fn bfs_order(
    gr: &Gr,
    s: usize,
) -> Vec<usize> {
    let mut seen = vec![false; gr.n];
    let mut out = Vec::new();
    let mut queue = VecDeque::from([s]);
    seen[s] = true;
    while let Some(u) = queue.pop_front() {
        out.push(u);
        for e in &gr.adj[u] {
            if !seen[e.to] {
                seen[e.to] = true;
                queue.push_back(e.to);
            }
        }
    }
    out
}

fn graph_dfs(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [g, s] = a else { return None };
    let gr = read(cx, *g)?;
    let s = vertex(cx, *s, &gr)?;
    let mut out = Vec::new();
    dfs_visit(&gr, s, &mut vec![false; gr.n], &mut out);
    Some(indices(&out))
}

fn graph_bfs(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [g, s] = a else { return None };
    let gr = read(cx, *g)?;
    let s = vertex(cx, *s, &gr)?;
    Some(indices(&bfs_order(&gr, s)))
}

fn components(gr: &Gr) -> Vec<Vec<usize>> {
    let mut seen = vec![false; gr.n];
    let mut out = Vec::new();
    for u in 0..gr.n {
        if !seen[u] {
            let comp = bfs_order(gr, u);
            for &v in &comp {
                seen[v] = true;
            }
            out.push(comp);
        }
    }
    out
}

fn graph_components(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let gr = read(cx, *a.first()?)?;
    Some(V::List(components(&gr).iter().map(|c| indices(c)).collect()))
}

fn graph_is_connected(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let gr = read(cx, *a.first()?)?;
    Some(V::Bool(components(&gr).len() == 1))
}

struct Tarjan {
    time: usize,
    disc: Vec<Option<usize>>,
    low: Vec<usize>,
    stack: Vec<usize>,
    on: Vec<bool>,
    out: Vec<Vec<usize>>,
}

fn tarjan(
    gr: &Gr,
    u: usize,
    s: &mut Tarjan,
) {
    s.disc[u] = Some(s.time);
    s.low[u] = s.time;
    s.time += 1;
    s.stack.push(u);
    s.on[u] = true;
    for e in &gr.adj[u] {
        let v = e.to;
        match s.disc[v] {
            | None => {
                tarjan(gr, v, s);
                s.low[u] = s.low[u].min(s.low[v]);
            },
            | Some(dv) if s.on[v] => s.low[u] = s.low[u].min(dv),
            | Some(_) => {},
        }
    }
    if Some(s.low[u]) == s.disc[u] {
        let mut comp = Vec::new();
        while let Some(v) = s.stack.pop() {
            s.on[v] = false;
            comp.push(v);
            if v == u {
                break;
            }
        }
        comp.sort_unstable();
        s.out.push(comp);
    }
}

fn graph_scc(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let gr = read(cx, *a.first()?)?;
    let mut s = Tarjan {
        time: 0,
        disc: vec![None; gr.n],
        low: vec![0; gr.n],
        stack: Vec::new(),
        on: vec![false; gr.n],
        out: Vec::new(),
    };
    for u in 0..gr.n {
        if s.disc[u].is_none() {
            tarjan(&gr, u, &mut s);
        }
    }
    let mut comps = s.out;
    comps.sort_by_key(|c| c.first().copied());
    Some(V::List(comps.iter().map(|c| indices(c)).collect()))
}

fn directed_cycle(
    gr: &Gr,
    u: usize,
    color: &mut [u8],
) -> bool {
    color[u] = 1;
    for e in &gr.adj[u] {
        if color[e.to] == 1 || (color[e.to] == 0 && directed_cycle(gr, e.to, color)) {
            return true;
        }
    }
    color[u] = 2;
    false
}

fn undirected_cycle(
    gr: &Gr,
    u: usize,
    parent_edge: Option<usize>,
    seen: &mut [bool],
) -> bool {
    seen[u] = true;
    for e in &gr.adj[u] {
        if Some(e.id) == parent_edge {
            continue;
        }
        if seen[e.to] || undirected_cycle(gr, e.to, Some(e.id), seen) {
            return true;
        }
    }
    false
}

fn has_cycle(gr: &Gr) -> bool {
    if gr.directed {
        let mut color = vec![0_u8; gr.n];
        (0..gr.n).any(|u| color[u] == 0 && directed_cycle(gr, u, &mut color))
    } else {
        let mut seen = vec![false; gr.n];
        (0..gr.n).any(|u| !seen[u] && undirected_cycle(gr, u, None, &mut seen))
    }
}

fn graph_has_cycle(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let gr = read(cx, *a.first()?)?;
    Some(V::Bool(has_cycle(&gr)))
}

struct Bridges {
    time: usize,
    disc: Vec<Option<usize>>,
    low: Vec<usize>,
    cuts: Vec<(usize, usize)>,
    ap: Vec<bool>,
}

fn bridge_visit(
    gr: &Gr,
    u: usize,
    parent_edge: Option<usize>,
    s: &mut Bridges,
) {
    s.disc[u] = Some(s.time);
    s.low[u] = s.time;
    s.time += 1;
    let mut children = 0;
    for e in &gr.adj[u] {
        if Some(e.id) == parent_edge {
            continue;
        }
        let v = e.to;
        if let Some(dv) = s.disc[v] {
            s.low[u] = s.low[u].min(dv);
        } else {
            children += 1;
            bridge_visit(gr, v, Some(e.id), s);
            s.low[u] = s.low[u].min(s.low[v]);
            if Some(s.low[v]) > s.disc[u] {
                s.cuts.push((u, v));
            }
            if parent_edge.is_some() && Some(s.low[v]) >= s.disc[u] {
                s.ap[u] = true;
            }
        }
    }
    if parent_edge.is_none() && children > 1 {
        s.ap[u] = true;
    }
}

fn graph_bridges(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let gr = read(cx, *a.first()?)?;
    let mut s = Bridges {
        time: 0,
        disc: vec![None; gr.n],
        low: vec![0; gr.n],
        cuts: Vec::new(),
        ap: vec![false; gr.n],
    };
    for u in 0..gr.n {
        if s.disc[u].is_none() {
            bridge_visit(&gr, u, None, &mut s);
        }
    }
    let aps: Vec<usize> = (0..gr.n).filter(|&u| s.ap[u]).collect();
    Some(V::List(vec![V::List(s.cuts.iter().map(|&(u, v)| pair(u, v)).collect()), indices(&aps)]))
}

// ---------------- spanning trees ----------------

struct Dsu(Vec<usize>);

impl Dsu {
    fn find(
        &mut self,
        x: usize,
    ) -> usize {
        if self.0[x] != x {
            let root = self.find(self.0[x]);
            self.0[x] = root;
        }
        self.0[x]
    }

    fn union(
        &mut self,
        a: usize,
        b: usize,
    ) -> bool {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra == rb {
            return false;
        }
        self.0[ra] = rb;
        true
    }
}

fn edge_value((u, v, w): (usize, usize, NodeId)) -> V {
    V::List(vec![V::uint(u), V::uint(v), V::Node(w)])
}

fn graph_kruskal(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let gr = read(cx, *a.first()?)?;
    let mut edges = gr.get_edges();
    edges.sort_by(|x, y| key(cx.graph, x.2).total_cmp(&key(cx.graph, y.2)));
    let mut dsu = Dsu((0..gr.n).collect());
    let tree: Vec<V> = edges.into_iter().filter(|&(u, v, _)| dsu.union(u, v)).map(edge_value).collect();
    Some(V::List(tree))
}

fn graph_prim(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [g, s] = a else { return None };
    let gr = read(cx, *g)?;
    let s = vertex(cx, *s, &gr)?;
    let mut inside = vec![false; gr.n];
    inside[s] = true;
    let mut tree = Vec::new();
    loop {
        let mut best: Option<(f64, usize, usize, NodeId)> = None;
        for u in (0..gr.n).filter(|&u| inside[u]) {
            for e in &gr.adj[u] {
                if !inside[e.to] {
                    let k = key(cx.graph, e.w);
                    if best.is_none_or(|b| k < b.0) {
                        best = Some((k, u, e.to, e.w));
                    }
                }
            }
        }
        let Some((_, u, v, w)) = best else { break };
        inside[v] = true;
        tree.push(edge_value((u, v, w)));
    }
    Some(V::List(tree))
}

// ---------------- flows ----------------

/// Capacities as exact rationals; `None` if any weight is not a number.
fn capacities(
    g: &Graph,
    gr: &Gr,
) -> Option<(Vec<Vec<BigRational>>, bool)> {
    let mut cap = vec![vec![BigRational::zero(); gr.n]; gr.n];
    let mut float = false;
    for (u, list) in gr.adj.iter().enumerate() {
        for e in list {
            let w = g.number_of(e.w)?;
            float |= !w.is_exact();
            cap[u][e.to] += rational(w)?;
        }
    }
    Some((cap, float))
}

fn number_value(
    cx: &mut Cx<'_>,
    r: &BigRational,
    float: bool,
) -> V {
    let number = if float { Number::Float(r.to_f64().unwrap_or(f64::NAN)) } else { Number::rat(r.clone()) };
    V::Node(cx.graph.num(number))
}

fn edmonds_karp(
    mut cap: Vec<Vec<BigRational>>,
    s: usize,
    t: usize,
) -> BigRational {
    let n = cap.len();
    let mut flow = BigRational::zero();
    loop {
        let mut parent = vec![usize::MAX; n];
        parent[s] = s;
        let mut queue = VecDeque::from([s]);
        while let Some(u) = queue.pop_front() {
            for v in 0..n {
                if parent[v] == usize::MAX && cap[u][v].is_positive() {
                    parent[v] = u;
                    queue.push_back(v);
                }
            }
        }
        if parent[t] == usize::MAX {
            return flow;
        }
        let mut bottleneck: Option<BigRational> = None;
        let mut v = t;
        while v != s {
            let c = cap[parent[v]][v].clone();
            bottleneck = Some(match bottleneck {
                | Some(b) => b.min(c),
                | None => c,
            });
            v = parent[v];
        }
        let bottleneck = bottleneck.unwrap_or_else(BigRational::zero);
        let mut v = t;
        while v != s {
            let u = parent[v];
            cap[u][v] -= &bottleneck;
            cap[v][u] += &bottleneck;
            v = u;
        }
        flow += bottleneck;
    }
}

fn dinic_push(
    cap: &mut [Vec<BigRational>],
    level: &[usize],
    it: &mut [usize],
    u: usize,
    t: usize,
    limit: &BigRational,
) -> BigRational {
    if u == t {
        return limit.clone();
    }
    while it[u] < cap.len() {
        let v = it[u];
        if level[v] == level[u].saturating_add(1) && cap[u][v].is_positive() {
            let room = limit.clone().min(cap[u][v].clone());
            let pushed = dinic_push(cap, level, it, v, t, &room);
            if pushed.is_positive() {
                cap[u][v] -= &pushed;
                cap[v][u] += &pushed;
                return pushed;
            }
        }
        it[u] += 1;
    }
    BigRational::zero()
}

fn dinic(
    mut cap: Vec<Vec<BigRational>>,
    s: usize,
    t: usize,
) -> BigRational {
    let n = cap.len();
    let mut flow = BigRational::zero();
    let mut total = BigRational::zero();
    for row in &cap {
        for c in row {
            total += c;
        }
    }
    loop {
        let mut level = vec![usize::MAX; n];
        level[s] = 0;
        let mut queue = VecDeque::from([s]);
        while let Some(u) = queue.pop_front() {
            for v in 0..n {
                if level[v] == usize::MAX && cap[u][v].is_positive() {
                    level[v] = level[u] + 1;
                    queue.push_back(v);
                }
            }
        }
        if level[t] == usize::MAX {
            return flow;
        }
        let mut it = vec![0; n];
        loop {
            let pushed = dinic_push(&mut cap, &level, &mut it, s, t, &total);
            if !pushed.is_positive() {
                break;
            }
            flow += pushed;
        }
    }
}

fn max_flow(
    cx: &mut Cx<'_>,
    a: &[NodeId],
    algorithm: fn(Vec<Vec<BigRational>>, usize, usize) -> BigRational,
) -> Option<V> {
    let [g, s, t] = a else { return None };
    let gr = read(cx, *g)?;
    let (s, t) = (vertex(cx, *s, &gr)?, vertex(cx, *t, &gr)?);
    if s == t {
        return None;
    }
    let (cap, float) = capacities(cx.graph, &gr)?;
    let flow = algorithm(cap, s, t);
    Some(number_value(cx, &flow, float))
}

struct Arc {
    to: usize,
    cap: BigRational,
    cost: BigRational,
}

fn graph_min_cost_flow(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [g, s, t] = a else { return None };
    let gr = read(cx, *g)?;
    let (s, t) = (vertex(cx, *s, &gr)?, vertex(cx, *t, &gr)?);
    if s == t {
        return None;
    }
    let mut arcs: Vec<Arc> = Vec::new();
    let mut out: Vec<Vec<usize>> = vec![Vec::new(); gr.n];
    let mut float = false;
    for (u, list) in gr.adj.iter().enumerate() {
        for e in list {
            let parts = items(cx.graph, e.w)?;
            let [cap, cost] = parts.as_slice() else { return None };
            let (cap, cost) = (cx.graph.number_of(*cap)?, cx.graph.number_of(*cost)?);
            float |= !cap.is_exact() || !cost.is_exact();
            let (cap, cost) = (rational(cap)?, rational(cost)?);
            out[u].push(arcs.len());
            arcs.push(Arc { to: e.to, cap, cost: cost.clone() });
            out[e.to].push(arcs.len());
            arcs.push(Arc { to: u, cap: BigRational::zero(), cost: -cost });
        }
    }
    let (mut flow, mut total) = (BigRational::zero(), BigRational::zero());
    loop {
        // Bellman-Ford over the residual arcs.
        let mut dist: Vec<Option<BigRational>> = vec![None; gr.n];
        let mut via = vec![usize::MAX; gr.n];
        dist[s] = Some(BigRational::zero());
        for _ in 0..gr.n {
            let mut changed = false;
            for u in 0..gr.n {
                let Some(du) = dist[u].clone() else { continue };
                for &id in &out[u] {
                    let arc = &arcs[id];
                    if arc.cap.is_positive() {
                        let candidate = &du + &arc.cost;
                        if dist[arc.to].as_ref().is_none_or(|d| candidate < *d) {
                            dist[arc.to] = Some(candidate);
                            via[arc.to] = id;
                            changed = true;
                        }
                    }
                }
            }
            if !changed {
                break;
            }
        }
        let Some(dt) = dist[t].clone() else { break };
        let mut bottleneck: Option<BigRational> = None;
        let mut v = t;
        while v != s {
            let id = via[v];
            let c = arcs[id].cap.clone();
            bottleneck = Some(match bottleneck {
                | Some(b) => b.min(c),
                | None => c,
            });
            v = arcs[id ^ 1].to;
        }
        let bottleneck = bottleneck?;
        let mut v = t;
        while v != s {
            let id = via[v];
            arcs[id].cap -= &bottleneck;
            arcs[id ^ 1].cap += &bottleneck;
            v = arcs[id ^ 1].to;
        }
        flow += &bottleneck;
        total += bottleneck * dt;
    }
    Some(V::List(vec![number_value(cx, &flow, float), number_value(cx, &total, float)]))
}

// ---------------- shortest paths ----------------

fn dist_list(
    dist: &[Dist],
    prev: &[Option<usize>],
) -> V {
    V::List(
        dist.iter()
            .zip(prev)
            .map(|(d, p)| V::List(vec![V::Node(d.term), V::Int(p.map_or(BigInt::from(-1), BigInt::from))]))
            .collect(),
    )
}

fn initial(
    cx: &mut Cx<'_>,
    n: usize,
    s: usize,
) -> Vec<Dist> {
    let inf = infinity(cx.graph);
    let zero = cx.graph.int(0);
    (0..n)
        .map(|v| if v == s { Dist { term: zero, key: 0.0 } } else { Dist { term: inf, key: f64::INFINITY } })
        .collect()
}

fn graph_dijkstra(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [g, s] = a else { return None };
    let gr = read(cx, *g)?;
    let s = vertex(cx, *s, &gr)?;
    let mut dist = initial(cx, gr.n, s);
    let mut prev: Vec<Option<usize>> = vec![None; gr.n];
    let mut done = vec![false; gr.n];
    loop {
        let mut pick: Option<usize> = None;
        for v in 0..gr.n {
            if !done[v] && dist[v].key.is_finite() && pick.is_none_or(|p| dist[v].key < dist[p].key) {
                pick = Some(v);
            }
        }
        let Some(u) = pick else { break };
        done[u] = true;
        for e in &gr.adj[u] {
            let w = Dist { term: e.w, key: key(cx.graph, e.w) };
            let candidate = dist_add(cx.graph, dist[u], w);
            if candidate.key.is_finite() && candidate.key < dist[e.to].key {
                dist[e.to] = candidate;
                prev[e.to] = Some(u);
            }
        }
    }
    Some(dist_list(&dist, &prev))
}

/// One round of relaxation over every edge; whether anything changed.
fn relax_all(
    cx: &mut Cx<'_>,
    gr: &Gr,
    dist: &mut [Dist],
    prev: &mut [Option<usize>],
) -> bool {
    let mut changed = false;
    for u in 0..gr.n {
        if !dist[u].key.is_finite() {
            continue;
        }
        for e in &gr.adj[u] {
            let w = Dist { term: e.w, key: key(cx.graph, e.w) };
            if !w.key.is_finite() {
                continue;
            }
            let candidate = dist_add(cx.graph, dist[u], w);
            if candidate.key < dist[e.to].key {
                dist[e.to] = candidate;
                prev[e.to] = Some(u);
                changed = true;
            }
        }
    }
    changed
}

fn graph_bellman_ford(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [g, s] = a else { return None };
    let gr = read(cx, *g)?;
    let s = vertex(cx, *s, &gr)?;
    let mut dist = initial(cx, gr.n, s);
    let mut prev: Vec<Option<usize>> = vec![None; gr.n];
    for _ in 1..gr.n {
        if !relax_all(cx, &gr, &mut dist, &mut prev) {
            break;
        }
    }
    let (mut probe, mut probe_prev) = (dist.clone(), prev.clone());
    if relax_all(cx, &gr, &mut probe, &mut probe_prev) {
        return Some(V::Bool(false));
    }
    Some(dist_list(&dist, &prev))
}

fn graph_floyd_warshall(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let gr = read(cx, *a.first()?)?;
    let inf = infinity(cx.graph);
    let zero = cx.graph.int(0);
    let n = gr.n;
    let mut d = vec![vec![Dist { term: inf, key: f64::INFINITY }; n]; n];
    // Edges with a symbolic weight are kept as they are (they cannot be
    // compared), a numeric weight replaces a larger one.
    let mut direct = vec![vec![false; n]; n];
    for (i, row) in d.iter_mut().enumerate() {
        row[i] = Dist { term: zero, key: 0.0 };
    }
    for (u, list) in gr.adj.iter().enumerate() {
        for e in list {
            let w = Dist { term: e.w, key: key(cx.graph, e.w) };
            if u == e.to {
                if w.key < 0.0 {
                    d[u][u] = w;
                }
            } else if !direct[u][e.to] || w.key < d[u][e.to].key {
                d[u][e.to] = w;
                direct[u][e.to] = true;
            }
        }
    }
    for k in 0..n {
        for i in 0..n {
            for j in 0..n {
                let candidate = dist_add(cx.graph, d[i][k], d[k][j]);
                if candidate.key < d[i][j].key {
                    d[i][j] = candidate;
                }
            }
        }
    }
    let rows: Vec<Vec<NodeId>> = d.iter().map(|r| r.iter().map(|x| x.term).collect()).collect();
    Some(V::Node(matrix(cx.graph, &rows)))
}

fn graph_shortest_path_unweighted(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [g, s] = a else { return None };
    let gr = read(cx, *g)?;
    let s = vertex(cx, *s, &gr)?;
    let mut dist: Vec<Option<usize>> = vec![None; gr.n];
    let mut prev: Vec<Option<usize>> = vec![None; gr.n];
    dist[s] = Some(0);
    let mut order = Vec::new();
    let mut queue = VecDeque::from([s]);
    while let Some(u) = queue.pop_front() {
        order.push(u);
        for e in &gr.adj[u] {
            if dist[e.to].is_none() {
                dist[e.to] = dist[u].map(|d| d + 1);
                prev[e.to] = Some(u);
                queue.push_back(e.to);
            }
        }
    }
    Some(V::List(
        order
            .into_iter()
            .map(|v| {
                V::List(vec![
                    V::uint(v),
                    V::uint(dist[v].unwrap_or(0)),
                    V::Int(prev[v].map_or_else(|| BigInt::from(-1), BigInt::from)),
                ])
            })
            .collect(),
    ))
}

// ---------------- bipartite graphs and matchings ----------------

fn bipartition(gr: &Gr) -> Option<Vec<u8>> {
    let mut color = vec![u8::MAX; gr.n];
    for s in 0..gr.n {
        if color[s] != u8::MAX {
            continue;
        }
        color[s] = 0;
        let mut queue = VecDeque::from([s]);
        while let Some(u) = queue.pop_front() {
            for e in &gr.adj[u] {
                if color[e.to] == u8::MAX {
                    color[e.to] = 1 - color[u];
                    queue.push_back(e.to);
                } else if color[e.to] == color[u] {
                    return None;
                }
            }
        }
    }
    Some(color)
}

fn graph_is_bipartite(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let gr = read(cx, *a.first()?)?;
    Some(match bipartition(&gr) {
        | Some(c) => V::List(c.into_iter().map(|x| V::int(i64::from(x))).collect()),
        | None => V::Bool(false),
    })
}

fn partition_of(
    cx: &Cx<'_>,
    node: NodeId,
    gr: &Gr,
) -> Option<Vec<u8>> {
    let p: Vec<u8> = idxs(cx.graph, node)?
        .into_iter()
        .map(|v| u8::try_from(v).ok().filter(|&b| b <= 1))
        .collect::<Option<_>>()?;
    (p.len() == gr.n).then_some(p)
}

fn kuhn(
    gr: &Gr,
    part: &[u8],
    u: usize,
    seen: &mut [bool],
    mate: &mut [Option<usize>],
) -> bool {
    for e in &gr.adj[u] {
        let v = e.to;
        if part[v] == 1 && !seen[v] {
            seen[v] = true;
            if mate[v].is_none_or(|w| kuhn(gr, part, w, seen, mate)) {
                mate[v] = Some(u);
                return true;
            }
        }
    }
    false
}

fn graph_bipartite_matching(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [g, p] = a else { return None };
    let gr = read(cx, *g)?;
    let part = partition_of(cx, *p, &gr)?;
    let mut mate = vec![None; gr.n];
    for u in (0..gr.n).filter(|&u| part[u] == 0) {
        kuhn(&gr, &part, u, &mut vec![false; gr.n], &mut mate);
    }
    Some(V::List(mate.iter().enumerate().filter_map(|(v, u)| u.map(|u| pair(u, v))).collect()))
}

struct HopcroftKarp<'a> {
    gr: &'a Gr,
    part: &'a [u8],
    mate: Vec<Option<usize>>,
    dist: Vec<usize>,
}

impl HopcroftKarp<'_> {
    fn bfs(&mut self) -> bool {
        let mut queue = VecDeque::new();
        for u in 0..self.gr.n {
            if self.part[u] == 0 && self.mate[u].is_none() {
                self.dist[u] = 0;
                queue.push_back(u);
            } else {
                self.dist[u] = usize::MAX;
            }
        }
        let mut found = false;
        while let Some(u) = queue.pop_front() {
            for e in &self.gr.adj[u] {
                if self.part[e.to] != 1 {
                    continue;
                }
                match self.mate[e.to] {
                    | None => found = true,
                    | Some(w) => {
                        if self.dist[w] == usize::MAX {
                            self.dist[w] = self.dist[u] + 1;
                            queue.push_back(w);
                        }
                    },
                }
            }
        }
        found
    }

    fn dfs(
        &mut self,
        u: usize,
    ) -> bool {
        let targets: Vec<usize> = self.gr.adj[u].iter().map(|e| e.to).collect();
        for v in targets {
            if self.part[v] != 1 {
                continue;
            }
            let ok = match self.mate[v] {
                | None => true,
                | Some(w) => self.dist[w] == self.dist[u].saturating_add(1) && self.dfs(w),
            };
            if ok {
                self.mate[v] = Some(u);
                self.mate[u] = Some(v);
                return true;
            }
        }
        self.dist[u] = usize::MAX;
        false
    }
}

fn graph_hopcroft_karp(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [g, p] = a else { return None };
    let gr = read(cx, *g)?;
    let part = partition_of(cx, *p, &gr)?;
    let mut hk = HopcroftKarp {
        gr: &gr,
        part: &part,
        mate: vec![None; gr.n],
        dist: vec![0; gr.n],
    };
    while hk.bfs() {
        #[allow(clippy::needless_range_loop)] // index is used for more than one array / arithmetic; iterator form would not be clearer
        for u in 0..gr.n {
            if part[u] == 0 && hk.mate[u].is_none() {
                hk.dfs(u);
            }
        }
    }
    Some(V::List(
        (0..gr.n).filter(|&u| part[u] == 0).filter_map(|u| hk.mate[u].map(|v| pair(u, v))).collect(),
    ))
}

fn graph_vertex_cover(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [g, p, m] = a else { return None };
    let gr = read(cx, *g)?;
    let part = partition_of(cx, *p, &gr)?;
    let mut mate: Vec<Option<usize>> = vec![None; gr.n];
    for e in items(cx.graph, *m)? {
        let uv = idxs(cx.graph, e)?;
        let [u, v] = uv.as_slice() else { return None };
        if *u >= gr.n || *v >= gr.n {
            return None;
        }
        mate[*u] = Some(*v);
        mate[*v] = Some(*u);
    }
    // Vertices reachable from unmatched left vertices by alternating paths.
    let mut reached = vec![false; gr.n];
    let mut queue = VecDeque::new();
    for u in 0..gr.n {
        if part[u] == 0 && mate[u].is_none() {
            reached[u] = true;
            queue.push_back(u);
        }
    }
    while let Some(u) = queue.pop_front() {
        if part[u] == 0 {
            for e in &gr.adj[u] {
                if part[e.to] == 1 && mate[u] != Some(e.to) && !reached[e.to] {
                    reached[e.to] = true;
                    queue.push_back(e.to);
                }
            }
        } else if let Some(w) = mate[u]
            && !reached[w] {
                reached[w] = true;
                queue.push_back(w);
            }
    }
    let cover: Vec<usize> =
        (0..gr.n).filter(|&v| (part[v] == 0 && !reached[v]) || (part[v] == 1 && reached[v])).collect();
    Some(indices(&cover))
}

/// Edmonds' blossom algorithm: a maximum matching of a general graph.
struct Blossom<'a> {
    adj: &'a [Vec<usize>],
    mate: Vec<usize>,
    parent: Vec<usize>,
    base: Vec<usize>,
}

const NONE: usize = usize::MAX;

impl Blossom<'_> {
    fn lca(
        &self,
        mut a: usize,
        mut b: usize,
    ) -> usize {
        let mut used = vec![false; self.base.len()];
        loop {
            a = self.base[a];
            used[a] = true;
            if self.mate[a] == NONE {
                break;
            }
            a = self.parent[self.mate[a]];
        }
        loop {
            b = self.base[b];
            if used[b] {
                return b;
            }
            b = self.parent[self.mate[b]];
        }
    }

    fn mark_path(
        &mut self,
        mut v: usize,
        b: usize,
        mut child: usize,
        blossom: &mut [bool],
    ) {
        while self.base[v] != b {
            blossom[self.base[v]] = true;
            blossom[self.base[self.mate[v]]] = true;
            self.parent[v] = child;
            child = self.mate[v];
            v = self.parent[self.mate[v]];
        }
    }

    fn find_path(
        &mut self,
        root: usize,
    ) -> usize {
        let n = self.base.len();
        let mut used = vec![false; n];
        self.parent = vec![NONE; n];
        for (i, b) in self.base.iter_mut().enumerate() {
            *b = i;
        }
        used[root] = true;
        let mut queue = vec![root];
        let mut head = 0;
        while head < queue.len() {
            let v = queue[head];
            head += 1;
            for k in 0..self.adj[v].len() {
                let to = self.adj[v][k];
                if to == v || self.base[v] == self.base[to] || self.mate[v] == to {
                    continue;
                }
                if to == root || (self.mate[to] != NONE && self.parent[self.mate[to]] != NONE) {
                    let current = self.lca(v, to);
                    let mut blossom = vec![false; n];
                    self.mark_path(v, current, to, &mut blossom);
                    self.mark_path(to, current, v, &mut blossom);
                    for i in 0..n {
                        if blossom[self.base[i]] {
                            self.base[i] = current;
                            if !used[i] {
                                used[i] = true;
                                queue.push(i);
                            }
                        }
                    }
                } else if self.parent[to] == NONE {
                    self.parent[to] = v;
                    if self.mate[to] == NONE {
                        return to;
                    }
                    let next = self.mate[to];
                    used[next] = true;
                    queue.push(next);
                }
            }
        }
        NONE
    }
}

fn graph_blossom(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let gr = read(cx, *a.first()?)?;
    let adj: Vec<Vec<usize>> = gr.adj.iter().map(|l| l.iter().map(|e| e.to).collect()).collect();
    let mut b = Blossom {
        adj: &adj,
        mate: vec![NONE; gr.n],
        parent: vec![NONE; gr.n],
        base: (0..gr.n).collect(),
    };
    for i in 0..gr.n {
        if b.mate[i] == NONE {
            let mut v = b.find_path(i);
            while v != NONE {
                let pv = b.parent[v];
                let ppv = b.mate[pv];
                b.mate[v] = pv;
                b.mate[pv] = v;
                v = ppv;
            }
        }
    }
    Some(V::List((0..gr.n).filter(|&u| b.mate[u] != NONE && u < b.mate[u]).map(|u| pair(u, b.mate[u])).collect()))
}

// ---------------- topological order ----------------

fn kahn(gr: &Gr) -> Option<Vec<usize>> {
    if !gr.directed {
        return None;
    }
    let mut indegree: Vec<usize> = gr.radj.iter().map(Vec::len).collect();
    let mut queue: VecDeque<usize> = (0..gr.n).filter(|&u| indegree[u] == 0).collect();
    let mut order = Vec::new();
    while let Some(u) = queue.pop_front() {
        order.push(u);
        for e in &gr.adj[u] {
            indegree[e.to] -= 1;
            if indegree[e.to] == 0 {
                queue.push_back(e.to);
            }
        }
    }
    (order.len() == gr.n).then_some(order)
}

fn topo_visit(
    gr: &Gr,
    u: usize,
    color: &mut [u8],
    post: &mut Vec<usize>,
) -> bool {
    color[u] = 1;
    for e in &gr.adj[u] {
        match color[e.to] {
            | 1 => return false,
            | 0 => {
                if !topo_visit(gr, e.to, color, post) {
                    return false;
                }
            },
            | _ => {},
        }
    }
    color[u] = 2;
    post.push(u);
    true
}

fn topological_dfs(gr: &Gr) -> Option<Vec<usize>> {
    if !gr.directed {
        return None;
    }
    let mut color = vec![0_u8; gr.n];
    let mut post = Vec::new();
    for u in 0..gr.n {
        if color[u] == 0 && !topo_visit(gr, u, &mut color, &mut post) {
            return None;
        }
    }
    post.reverse();
    Some(post)
}

fn toposort(
    cx: &mut Cx<'_>,
    a: &[NodeId],
    algorithm: fn(&Gr) -> Option<Vec<usize>>,
) -> Option<V> {
    let gr = read(cx, *a.first()?)?;
    Some(algorithm(&gr).map_or(V::Bool(false), |o| indices(&o)))
}

// ---------------- spectra ----------------

fn spectral_analysis(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let m = *a.first()?;
    let values = apply(cx.graph, "eigenvals", &[m])?;
    let vectors = apply(cx.graph, "eigenvects", &[m])?;
    Some(V::nodes(&[values, vectors]))
}

/// Eigenvalues of a symmetric matrix by cyclic Jacobi sweeps, ascending.
fn symmetric_eigenvalues(mut m: Vec<Vec<f64>>) -> Vec<f64> {
    let n = m.len();
    for _ in 0..100 {
        let mut off = 0.0;
        for (i, row) in m.iter().enumerate() {
            for (j, x) in row.iter().enumerate() {
                if i != j {
                    off += x * x;
                }
            }
        }
        if off < 1e-24 {
            break;
        }
        for p in 0..n {
            for q in (p + 1)..n {
                if m[p][q].abs() < 1e-300 {
                    continue;
                }
                let theta = (m[q][q] - m[p][p]) / (2.0 * m[p][q]);
                let t = if theta == 0.0 { 1.0 } else { theta.signum() / (theta.abs() + theta.hypot(1.0)) };
                let c = 1.0 / t.hypot(1.0);
                let s = t * c;
                for row in &mut m {
                    let (kp, kq) = (row[p], row[q]);
                    row[p] = c * kp - s * kq;
                    row[q] = s * kp + c * kq;
                }
                #[allow(clippy::needless_range_loop)] // index is used for more than one array / arithmetic; iterator form would not be clearer
                for k in 0..n {
                    let (pk, qk) = (m[p][k], m[q][k]);
                    m[p][k] = c * pk - s * qk;
                    m[q][k] = s * pk + c * qk;
                }
            }
        }
    }
    let mut values: Vec<f64> = (0..n).map(|i| m[i][i]).collect();
    values.sort_by(f64::total_cmp);
    values
}

fn algebraic_connectivity(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let gr = read(cx, *a.first()?)?;
    if gr.directed || gr.n < 2 {
        return None;
    }
    let mut l = vec![vec![0.0_f64; gr.n]; gr.n];
    for (u, list) in gr.adj.iter().enumerate() {
        for e in list {
            let w = cx.graph.number_of(e.w)?.to_f64();
            if e.to != u {
                l[u][e.to] -= w;
                l[u][u] += w;
            }
        }
    }
    let values = symmetric_eigenvalues(l);
    Some(V::Float(values[1]))
}

// ---------------- isomorphism and colouring ----------------

/// Weisfeiler-Lehman colour refinement of both graphs with shared colour
/// names; the graphs pass when the colour histograms agree.
fn wl_equivalent(
    g1: &Gr,
    g2: &Gr,
) -> bool {
    use std::collections::BTreeMap;
    if g1.n != g2.n || g1.get_edges().len() != g2.get_edges().len() {
        return false;
    }
    let n = g1.n;
    let mut c1: Vec<usize> = (0..n).map(|i| g1.radj[i].len()).collect();
    let mut c2: Vec<usize> = (0..n).map(|i| g2.radj[i].len()).collect();
    for _ in 0..n {
        let signature = |gr: &Gr, c: &[usize], i: usize| {
            let mut neighbours: Vec<usize> = gr.adj[i].iter().map(|e| c[e.to]).collect();
            neighbours.sort_unstable();
            (c[i], neighbours)
        };
        let sigs1: Vec<_> = (0..n).map(|i| signature(g1, &c1, i)).collect();
        let sigs2: Vec<_> = (0..n).map(|i| signature(g2, &c2, i)).collect();
        let mut names: BTreeMap<(usize, Vec<usize>), usize> = BTreeMap::new();
        for s in sigs1.iter().chain(&sigs2) {
            let next = names.len();
            names.entry(s.clone()).or_insert(next);
        }
        c1 = sigs1.iter().map(|s| names[s]).collect();
        c2 = sigs2.iter().map(|s| names[s]).collect();
    }
    c1.sort_unstable();
    c2.sort_unstable();
    c1 == c2
}

fn graph_isomorphic_heuristic(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [x, y] = a else { return None };
    let (g1, g2) = (read(cx, *x)?, read(cx, *y)?);
    Some(V::Bool(wl_equivalent(&g1, &g2)))
}

fn multiplicity(gr: &Gr) -> Vec<Vec<usize>> {
    let mut m = vec![vec![0_usize; gr.n]; gr.n];
    for (u, list) in gr.adj.iter().enumerate() {
        for e in list {
            m[u][e.to] += 1;
        }
    }
    m
}

fn extend(
    m1: &[Vec<usize>],
    m2: &[Vec<usize>],
    map: &mut Vec<usize>,
    used: &mut [bool],
) -> bool {
    let k = map.len();
    if k == m1.len() {
        return true;
    }
    for c in 0..m2.len() {
        if used[c] || m1[k][k] != m2[c][c] {
            continue;
        }
        if (0..k).all(|j| m1[k][j] == m2[c][map[j]] && m1[j][k] == m2[map[j]][c]) {
            map.push(c);
            used[c] = true;
            if extend(m1, m2, map, used) {
                return true;
            }
            used[c] = false;
            map.pop();
        }
    }
    false
}

fn graph_isomorphic(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [x, y] = a else { return None };
    let (g1, g2) = (read(cx, *x)?, read(cx, *y)?);
    if g1.n > 64 || g2.n > 64 {
        return None;
    }
    if g1.n != g2.n || g1.directed != g2.directed || !wl_equivalent(&g1, &g2) {
        return Some(V::Bool(false));
    }
    let (m1, m2) = (multiplicity(&g1), multiplicity(&g2));
    Some(V::Bool(extend(&m1, &m2, &mut Vec::new(), &mut vec![false; g1.n])))
}

fn has_loop(gr: &Gr) -> bool {
    gr.edges.iter().any(|&(u, v, _)| u == v)
}

fn graph_greedy_coloring(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let gr = read(cx, *a.first()?)?;
    if has_loop(&gr) {
        return None;
    }
    let mut order: Vec<usize> = (0..gr.n).collect();
    order.sort_by(|&x, &y| gr.adj[y].len().cmp(&gr.adj[x].len()));
    let mut color: Vec<Option<usize>> = vec![None; gr.n];
    let mut next = 0;
    for &u in &order {
        if color[u].is_some() {
            continue;
        }
        color[u] = Some(next);
        for &v in &order {
            if color[v].is_none() && gr.adj[v].iter().all(|e| color[e.to] != Some(next)) {
                color[v] = Some(next);
            }
        }
        next += 1;
    }
    Some(V::List(color.into_iter().map(|c| V::uint(c.unwrap_or(0))).collect()))
}

fn colorable(
    gr: &Gr,
    k: usize,
    colors: &mut [usize],
    u: usize,
) -> bool {
    if u == gr.n {
        return true;
    }
    for c in 1..=k {
        if gr.adj[u].iter().all(|e| colors[e.to] != c) {
            colors[u] = c;
            if colorable(gr, k, colors, u + 1) {
                return true;
            }
            colors[u] = 0;
        }
    }
    false
}

fn graph_chromatic_number(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let gr = read(cx, *a.first()?)?;
    if has_loop(&gr) || gr.n > 40 {
        return None;
    }
    if gr.n == 0 {
        return Some(V::int(0));
    }
    let k = (1..=gr.n).find(|&k| colorable(&gr, k, &mut vec![0; gr.n], 0))?;
    Some(V::uint(k))
}

// ---------------- graph operations ----------------

fn same_weight(
    g: &Graph,
    a: NodeId,
    b: NodeId,
) -> bool {
    g.find(a) == g.find(b)
}

fn two_graphs(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<(Gr, Gr)> {
    let [x, y] = a else { return None };
    Some((read(cx, *x)?, read(cx, *y)?))
}

fn induced_subgraph(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [g, vs] = a else { return None };
    let gr = read(cx, *g)?;
    let vs = idxs(cx.graph, *vs)?;
    let mut map = vec![None; gr.n];
    let mut count = 0;
    for &v in &vs {
        if v >= gr.n {
            return None;
        }
        if map[v].is_none() {
            map[v] = Some(count);
            count += 1;
        }
    }
    let mut sub = Gr::new(count, gr.directed);
    for &(u, v, w) in &gr.edges {
        if let (Some(a), Some(b)) = (map[u], map[v]) {
            sub.add_edge(a, b, w);
        }
    }
    Some(graph_value(cx, &sub))
}

fn graph_union(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (g1, g2) = two_graphs(cx, a)?;
    let mut out = Gr::new(g1.n.max(g2.n), g1.directed);
    for &(u, v, w) in &g1.edges {
        out.add_edge(u, v, w);
    }
    for (u, v, w) in g2.get_edges() {
        let present = out.edges.iter().any(|&(x, y, z)| {
            same_weight(cx.graph, z, w) && ((x, y) == (u, v) || (!out.directed && (x, y) == (v, u)))
        });
        if !present {
            out.add_edge(u, v, w);
        }
    }
    Some(graph_value(cx, &out))
}

fn graph_intersection(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (g1, g2) = two_graphs(cx, a)?;
    let n = g1.n.min(g2.n);
    let mut out = Gr::new(n, g1.directed && g2.directed);
    for (u, v, w) in g1.get_edges() {
        if u < n && v < n && g2.adj[u].iter().any(|e| e.to == v && same_weight(cx.graph, e.w, w)) {
            out.add_edge(u, v, w);
        }
    }
    Some(graph_value(cx, &out))
}

fn graph_complement(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let gr = read(cx, *a.first()?)?;
    let one = cx.graph.int(1);
    let mut out = Gr::new(gr.n, gr.directed);
    for i in 0..gr.n {
        for j in 0..gr.n {
            if i == j || (!gr.directed && i > j) {
                continue;
            }
            if !gr.adj[i].iter().any(|e| e.to == j) {
                out.add_edge(i, j, one);
            }
        }
    }
    Some(graph_value(cx, &out))
}

fn graph_cartesian(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (g1, g2) = two_graphs(cx, a)?;
    let mut out = Gr::new(g1.n * g2.n, g1.directed || g2.directed);
    let id = |u: usize, v: usize| u * g2.n + v;
    for v in 0..g2.n {
        for (u1, u2, w) in g1.get_edges() {
            out.add_edge(id(u1, v), id(u2, v), w);
        }
    }
    for u in 0..g1.n {
        for (v1, v2, w) in g2.get_edges() {
            out.add_edge(id(u, v1), id(u, v2), w);
        }
    }
    Some(graph_value(cx, &out))
}

fn graph_tensor(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (g1, g2) = two_graphs(cx, a)?;
    let mut out = Gr::new(g1.n * g2.n, g1.directed || g2.directed);
    let id = |u: usize, v: usize| u * g2.n + v;
    for (u1, u2, w1) in g1.get_edges() {
        for (v1, v2, w2) in g2.get_edges() {
            let w = prod(cx.graph, &[w1, w2]);
            out.add_edge(id(u1, v1), id(u2, v2), w);
            if !g1.directed && !g2.directed {
                out.add_edge(id(u1, v2), id(u2, v1), w);
            } else if !g1.directed {
                out.add_edge(id(u2, v1), id(u1, v2), w);
            } else if !g2.directed {
                out.add_edge(id(u1, v2), id(u2, v1), w);
            }
        }
    }
    Some(graph_value(cx, &out))
}

fn disjoint_union(
    g1: &Gr,
    g2: &Gr,
) -> Gr {
    let mut out = Gr::new(g1.n + g2.n, g1.directed && g2.directed);
    for (u, v, w) in g1.get_edges() {
        out.add_edge(u, v, w);
    }
    for (u, v, w) in g2.get_edges() {
        out.add_edge(g1.n + u, g1.n + v, w);
    }
    out
}

fn graph_disjoint_union(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (g1, g2) = two_graphs(cx, a)?;
    let out = disjoint_union(&g1, &g2);
    Some(graph_value(cx, &out))
}

fn graph_join(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (g1, g2) = two_graphs(cx, a)?;
    let one = cx.graph.int(1);
    let mut out = disjoint_union(&g1, &g2);
    for u in 0..g1.n {
        for v in 0..g2.n {
            out.add_edge(u, g1.n + v, one);
            if out.directed {
                out.add_edge(g1.n + v, u, one);
            }
        }
    }
    Some(graph_value(cx, &out))
}

pub(crate) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    def_inert(i, "graph", Arity::Variadic)?;
    def_inert(i, "digraph", Arity::Variadic)?;
    def(i, "graph_empty", Arity::Fixed(1), |cx, a| standard(cx, a, |_| Vec::new()))?;
    def(i, "graph_complete", Arity::Fixed(1), |cx, a| {
        standard(cx, a, |n| (0..n).flat_map(|u| ((u + 1)..n).map(move |v| (u, v))).collect())
    })?;
    def(i, "graph_path", Arity::Fixed(1), |cx, a| standard(cx, a, |n| (1..n).map(|v| (v - 1, v)).collect()))?;
    def(i, "graph_cycle", Arity::Fixed(1), |cx, a| {
        standard(cx, a, |n| (0..n).filter(|_| n >= 3).map(|u| (u, (u + 1) % n)).collect())
    })?;
    def(i, "graph_nodes", Arity::Fixed(1), graph_nodes)?;
    def(i, "graph_node_count", Arity::Fixed(1), |cx, a| {
        let gr = read(cx, *a.first()?)?;
        Some(V::uint(gr.n))
    })?;
    def(i, "graph_is_directed", Arity::Fixed(1), |cx, a| {
        let gr = read(cx, *a.first()?)?;
        Some(V::Bool(gr.directed))
    })?;
    def(i, "graph_node_id", Arity::Fixed(2), graph_node_id)?;
    def(i, "graph_add_node", Arity::Fixed(1), graph_add_node)?;
    def(i, "graph_add_edge", Arity::Variadic, graph_add_edge)?;
    def(i, "graph_add_hyperedge", Arity::Fixed(3), graph_add_hyperedge)?;
    def(i, "graph_hyperedges", Arity::Fixed(1), graph_hyperedges)?;
    def(i, "graph_neighbors", Arity::Fixed(2), graph_neighbors)?;
    def(i, "graph_out_degree", Arity::Fixed(2), |cx, a| graph_degree(cx, a, false))?;
    def(i, "graph_in_degree", Arity::Fixed(2), |cx, a| graph_degree(cx, a, true))?;
    def(i, "graph_edges", Arity::Fixed(1), graph_edges)?;
    def(i, "graph_adjacency", Arity::Fixed(1), graph_adjacency)?;
    def(i, "graph_incidence", Arity::Fixed(1), graph_incidence)?;
    def(i, "graph_laplacian", Arity::Fixed(1), graph_laplacian)?;

    def(i, "graph_dfs", Arity::Fixed(2), graph_dfs)?;
    def(i, "graph_bfs", Arity::Fixed(2), graph_bfs)?;
    def(i, "graph_components", Arity::Fixed(1), graph_components)?;
    def(i, "graph_is_connected", Arity::Fixed(1), graph_is_connected)?;
    def(i, "graph_scc", Arity::Fixed(1), graph_scc)?;
    def(i, "graph_has_cycle", Arity::Fixed(1), graph_has_cycle)?;
    def(i, "graph_bridges", Arity::Fixed(1), graph_bridges)?;
    def(i, "graph_kruskal", Arity::Fixed(1), graph_kruskal)?;
    def(i, "graph_prim", Arity::Fixed(2), graph_prim)?;
    def(i, "graph_edmonds_karp", Arity::Fixed(3), |cx, a| max_flow(cx, a, edmonds_karp))?;
    def(i, "graph_dinic", Arity::Fixed(3), |cx, a| max_flow(cx, a, dinic))?;
    def(i, "graph_min_cost_flow", Arity::Fixed(3), graph_min_cost_flow)?;
    def(i, "graph_dijkstra", Arity::Fixed(2), graph_dijkstra)?;
    def(i, "graph_bellman_ford", Arity::Fixed(2), graph_bellman_ford)?;
    def(i, "graph_floyd_warshall", Arity::Fixed(1), graph_floyd_warshall)?;
    def(i, "graph_shortest_path_unweighted", Arity::Fixed(2), graph_shortest_path_unweighted)?;
    def(i, "graph_is_bipartite", Arity::Fixed(1), graph_is_bipartite)?;
    def(i, "graph_bipartite_matching", Arity::Fixed(2), graph_bipartite_matching)?;
    def(i, "graph_hopcroft_karp", Arity::Fixed(2), graph_hopcroft_karp)?;
    def(i, "graph_vertex_cover", Arity::Fixed(3), graph_vertex_cover)?;
    def(i, "graph_blossom", Arity::Fixed(1), graph_blossom)?;
    def(i, "graph_toposort", Arity::Fixed(1), |cx, a| toposort(cx, a, kahn))?;
    def(i, "graph_toposort_kahn", Arity::Fixed(1), |cx, a| toposort(cx, a, kahn))?;
    def(i, "graph_toposort_dfs", Arity::Fixed(1), |cx, a| toposort(cx, a, topological_dfs))?;
    def_request(i, "spectral_analysis", Arity::Fixed(1), spectral_analysis)?;
    def(i, "algebraic_connectivity", Arity::Fixed(1), algebraic_connectivity)?;
    def(i, "graph_isomorphic_heuristic", Arity::Fixed(2), graph_isomorphic_heuristic)?;
    def(i, "graph_isomorphic", Arity::Fixed(2), graph_isomorphic)?;
    def(i, "graph_greedy_coloring", Arity::Fixed(1), graph_greedy_coloring)?;
    def(i, "graph_chromatic_number", Arity::Fixed(1), graph_chromatic_number)?;

    def(i, "induced_subgraph", Arity::Fixed(2), induced_subgraph)?;
    def(i, "graph_union", Arity::Fixed(2), graph_union)?;
    def(i, "graph_intersection", Arity::Fixed(2), graph_intersection)?;
    def(i, "graph_cartesian", Arity::Fixed(2), graph_cartesian)?;
    def(i, "graph_tensor", Arity::Fixed(2), graph_tensor)?;
    def(i, "graph_complement", Arity::Fixed(1), graph_complement)?;
    def(i, "graph_disjoint_union", Arity::Fixed(2), graph_disjoint_union)?;
    def(i, "graph_join", Arity::Fixed(2), graph_join)?;
    Ok(())
}
