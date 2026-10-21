//! # Multilevel nested dissection
//!
//! [`multilevel_nested_dissection_order`] orders a sparse symmetric pattern
//! by recursive vertex-separator bisection in the style of METIS:
//!
//! 1. **coarsening** by heavy-edge matching (vertex and edge weights are
//!    accumulated) until the graph has at most about 100 vertices;
//! 2. **initial bisection** of the coarsest graph by several greedy graph
//!    growing runs, each refined by Fiduccia-Mattheyses (FM);
//! 3. **uncoarsening**: the partition is projected to the next finer graph
//!    and refined by boundary FM (gain buckets with lazy invalidation, hill
//!    climbing with rollback to the best balanced state);
//! 4. the edge bisection is turned into a **vertex separator** by a minimum
//!    vertex cover of the bipartite boundary graph (maximum matching and
//!    Koenig's theorem), which is the smallest separator contained in the
//!    boundary;
//! 5. both halves are ordered recursively and the separator is numbered
//!    last; pieces of at most 100 vertices are ordered with
//!    [`crate::kernels::ordering::amd_order`].
//!
//! The top-level separators on grids are optimal (`k` for a `k x k` mesh,
//! `k^2` for a `k^3` cube). On meshes of 10^4 to 10^5 vertices the fill is
//! within a few percent of approximate minimum degree (and far below the
//! level-structure [`crate::kernels::ordering::nested_dissection_order`]);
//! [`hybrid_nested_dissection_order`] picks the best of both by exact
//! symbolic fill, so it is never worse than minimum degree. The elimination
//! tree of the nested-dissection orderings is balanced, which suits parallel
//! factorisation.

use crate::kernels::ordering::{amd_order, cholesky_fill};

/// Tuning parameters of [`multilevel_nested_dissection_with`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NdOptions {
    /// Pieces of at most this many vertices are ordered by minimum degree.
    pub leaf: usize,
    /// Allowed relative imbalance of a bisection (`0.1`: a side may hold up
    /// to 55 percent of the weight).
    pub imbalance: f64,
    /// Number of FM passes per level.
    pub fm_passes: usize,
}

impl Default for NdOptions {
    fn default() -> Self {
        Self { leaf: 2000, imbalance: 0.1, fm_passes: 8 }
    }
}

/// Vertex- and edge-weighted graph in compressed adjacency form.
struct WGraph {
    xadj: Vec<usize>,
    adj: Vec<usize>,
    ew: Vec<i64>,
    vw: Vec<i64>,
}

impl WGraph {
    fn n(&self) -> usize {
        self.vw.len()
    }

    fn from_lists(lists: &[Vec<usize>]) -> Self {
        let mut xadj = vec![0usize];
        let mut adj = Vec::new();
        for l in lists {
            adj.extend_from_slice(l);
            xadj.push(adj.len());
        }
        let ew = vec![1; adj.len()];
        Self { xadj, adj, ew, vw: vec![1; lists.len()] }
    }
}

struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
        self.0 >> 33
    }

    fn below(&mut self, n: usize) -> usize {
        (self.next() % n.max(1) as u64) as usize
    }
}

/// Heavy-edge matching; returns the coarse graph and `cmap` (fine -> coarse).
fn coarsen(g: &WGraph, rng: &mut Lcg) -> (WGraph, Vec<usize>) {
    let n = g.n();
    let mut order: Vec<usize> = (0..n).collect();
    for i in (1..n).rev() {
        order.swap(i, rng.below(i + 1));
    }
    let mut mate = vec![usize::MAX; n];
    for &v in &order {
        if mate[v] != usize::MAX {
            continue;
        }
        let mut best = usize::MAX;
        let mut bw = -1;
        for k in g.xadj[v]..g.xadj[v + 1] {
            let u = g.adj[k];
            if u != v && mate[u] == usize::MAX && g.ew[k] > bw {
                bw = g.ew[k];
                best = u;
            }
        }
        if best == usize::MAX {
            mate[v] = v;
        } else {
            mate[v] = best;
            mate[best] = v;
        }
    }
    let mut cmap = vec![usize::MAX; n];
    let mut nc = 0;
    for v in 0..n {
        if cmap[v] == usize::MAX {
            cmap[v] = nc;
            cmap[mate[v]] = nc;
            nc += 1;
        }
    }
    let mut vw = vec![0i64; nc];
    for v in 0..n {
        vw[cmap[v]] += g.vw[v];
    }
    let mut members: Vec<[usize; 2]> = vec![[usize::MAX; 2]; nc];
    for v in 0..n {
        let c = cmap[v];
        if members[c][0] == usize::MAX {
            members[c][0] = v;
        } else {
            members[c][1] = v;
        }
    }
    let mut xadj = vec![0usize];
    let mut adj: Vec<usize> = Vec::new();
    let mut ew: Vec<i64> = Vec::new();
    let mut pos = vec![usize::MAX; nc];
    for c in 0..nc {
        let start = adj.len();
        for &v in &members[c] {
            if v == usize::MAX {
                continue;
            }
            for k in g.xadj[v]..g.xadj[v + 1] {
                let cu = cmap[g.adj[k]];
                if cu == c {
                    continue;
                }
                if pos[cu] >= start && pos[cu] < adj.len() && adj[pos[cu]] == cu {
                    ew[pos[cu]] += g.ew[k];
                } else {
                    pos[cu] = adj.len();
                    adj.push(cu);
                    ew.push(g.ew[k]);
                }
            }
        }
        xadj.push(adj.len());
    }
    (WGraph { xadj, adj, ew, vw }, cmap)
}

fn edge_cut(g: &WGraph, part: &[u8]) -> i64 {
    let mut cut = 0;
    for v in 0..g.n() {
        for k in g.xadj[v]..g.xadj[v + 1] {
            if part[g.adj[k]] != part[v] {
                cut += g.ew[k];
            }
        }
    }
    cut / 2
}

/// Boundary Fiduccia-Mattheyses refinement of a bisection (in place).
fn fm_refine(g: &WGraph, part: &mut [u8], passes: usize, max_imbalance: f64) {
    let n = g.n();
    let total: i64 = g.vw.iter().sum();
    let limit = (0.5 * (1.0 + max_imbalance) * total as f64).ceil() as i64;
    let mut pw = [0i64; 2];
    for v in 0..n {
        pw[usize::from(part[v])] += g.vw[v];
    }
    for _ in 0..passes {
        let mut gain = vec![0i64; n];
        for v in 0..n {
            for k in g.xadj[v]..g.xadj[v + 1] {
                gain[v] += if part[g.adj[k]] == part[v] { -g.ew[k] } else { g.ew[k] };
            }
        }
        let mut heap = std::collections::BinaryHeap::new();
        for v in 0..n {
            // boundary vertices only
            if (g.xadj[v]..g.xadj[v + 1]).any(|k| part[g.adj[k]] != part[v]) {
                heap.push((gain[v], v));
            }
        }
        let mut locked = vec![false; n];
        let mut moves: Vec<usize> = Vec::new();
        let mut cut = edge_cut(g, part);
        let mut best_cut = cut;
        let mut best_len = 0usize;
        let mut best_imb = (pw[0] - pw[1]).abs();
        let mut since_best = 0usize;
        let mut deferred: Vec<(i64, usize)> = Vec::new();
        while let Some((gv, v)) = heap.pop() {
            if locked[v] || gv != gain[v] {
                continue;
            }
            let from = usize::from(part[v]);
            let to = 1 - from;
            let newto = pw[to] + g.vw[v];
            let balanced_after = newto <= limit;
            let helps_balance = pw[from] > limit;
            if !balanced_after && !helps_balance {
                deferred.push((gv, v));
                continue;
            }
            // apply the move
            locked[v] = true;
            part[v] = to as u8;
            pw[from] -= g.vw[v];
            pw[to] += g.vw[v];
            cut -= gv;
            moves.push(v);
            for k in g.xadj[v]..g.xadj[v + 1] {
                let u = g.adj[k];
                let w = g.ew[k];
                if usize::from(part[u]) == to {
                    gain[u] -= 2 * w;
                } else {
                    gain[u] += 2 * w;
                }
                if !locked[u] {
                    heap.push((gain[u], u));
                }
            }
            gain[v] = -gv;
            for d in deferred.drain(..) {
                if !locked[d.1] && d.0 == gain[d.1] {
                    heap.push(d);
                }
            }
            let imb = (pw[0] - pw[1]).abs();
            let ok = pw[0] <= limit && pw[1] <= limit;
            if ok && (cut < best_cut || (cut == best_cut && imb < best_imb)) {
                best_cut = cut;
                best_len = moves.len();
                best_imb = imb;
                since_best = 0;
            } else {
                since_best += 1;
                if since_best > 60.max(n / 50) {
                    break;
                }
            }
        }
        // roll back to the best state
        for &v in moves.iter().skip(best_len).rev() {
            let from = usize::from(part[v]);
            let to = 1 - from;
            part[v] = to as u8;
            pw[from] -= g.vw[v];
            pw[to] += g.vw[v];
        }
        if best_len == 0 {
            break;
        }
    }
}

/// Greedy graph growing from a random seed until half of the weight.
fn grow_partition(g: &WGraph, rng: &mut Lcg) -> Vec<u8> {
    let n = g.n();
    let total: i64 = g.vw.iter().sum();
    let mut part = vec![1u8; n];
    let mut w0 = 0i64;
    let mut visited = vec![false; n];
    let mut queue = std::collections::VecDeque::new();
    let seed = rng.below(n);
    queue.push_back(seed);
    visited[seed] = true;
    while w0 * 2 < total {
        let v = if let Some(v) = queue.pop_front() {
            v
        } else {
            // disconnected: jump to an unvisited vertex
            match (0..n).find(|&u| !visited[u]) {
                Some(u) => {
                    visited[u] = true;
                    u
                }
                None => break,
            }
        };
        part[v] = 0;
        w0 += g.vw[v];
        for k in g.xadj[v]..g.xadj[v + 1] {
            let u = g.adj[k];
            if !visited[u] {
                visited[u] = true;
                queue.push_back(u);
            }
        }
    }
    part
}

/// Multilevel edge bisection; returns the part (0/1) of every vertex.
fn bisect(g: &WGraph, rng: &mut Lcg, o: &NdOptions) -> Vec<u8> {
    let n = g.n();
    if n <= 60 {
        let mut best: Option<(i64, Vec<u8>)> = None;
        let tries = if n <= 8 { 4 } else { 10 };
        for _ in 0..tries {
            let mut p = grow_partition(g, rng);
            fm_refine(g, &mut p, o.fm_passes, o.imbalance);
            let c = edge_cut(g, &p);
            if best.as_ref().is_none_or(|b| c < b.0) {
                best = Some((c, p));
            }
        }
        return best.map_or_else(|| vec![0; n], |b| b.1);
    }
    let (cg, cmap) = coarsen(g, rng);
    if cg.n() * 10 > n * 9 {
        // no progress: bisect directly
        let mut p = grow_partition(g, rng);
        fm_refine(g, &mut p, o.fm_passes, o.imbalance);
        return p;
    }
    let cpart = bisect(&cg, rng, o);
    let mut part: Vec<u8> = (0..n).map(|v| cpart[cmap[v]]).collect();
    fm_refine(g, &mut part, o.fm_passes, o.imbalance);
    part
}

/// Minimum vertex cover of the bipartite boundary graph (Koenig).
fn vertex_separator(g: &WGraph, part: &[u8]) -> Vec<usize> {
    let n = g.n();
    let boundary: Vec<usize> = (0..n)
        .filter(|&v| (g.xadj[v]..g.xadj[v + 1]).any(|k| part[g.adj[k]] != part[v]))
        .collect();
    let mut lid = vec![usize::MAX; n];
    let (mut left, mut right) = (Vec::new(), Vec::new());
    for &v in &boundary {
        if part[v] == 0 {
            lid[v] = left.len();
            left.push(v);
        } else {
            lid[v] = right.len();
            right.push(v);
        }
    }
    let nbr = |v: usize| -> Vec<usize> {
        (g.xadj[v]..g.xadj[v + 1])
            .map(|k| g.adj[k])
            .filter(|&u| part[u] != part[v])
            .map(|u| lid[u])
            .collect()
    };
    let ladj: Vec<Vec<usize>> = left.iter().map(|&v| nbr(v)).collect();
    let mut match_l = vec![usize::MAX; left.len()];
    let mut match_r = vec![usize::MAX; right.len()];
    // augmenting-path matching (iterative DFS)
    for s in 0..left.len() {
        let mut visited = vec![false; right.len()];
        let mut stack: Vec<(usize, usize)> = vec![(s, 0)];
        let mut path: Vec<(usize, usize)> = Vec::new();
        let mut found = false;
        while let Some(&mut (u, ref mut idx)) = stack.last_mut() {
            if *idx < ladj[u].len() {
                let r = ladj[u][*idx];
                *idx += 1;
                if visited[r] {
                    continue;
                }
                visited[r] = true;
                path.push((u, r));
                if match_r[r] == usize::MAX {
                    found = true;
                    break;
                }
                stack.push((match_r[r], 0));
            } else {
                stack.pop();
                path.pop();
            }
        }
        if found {
            for &(u, r) in &path {
                match_l[u] = r;
                match_r[r] = u;
            }
        }
    }
    // Koenig: Z = vertices reachable from unmatched left vertices by alternating paths
    let mut zl = vec![false; left.len()];
    let mut zr = vec![false; right.len()];
    let mut stack: Vec<usize> = (0..left.len()).filter(|&u| match_l[u] == usize::MAX).collect();
    for &u in &stack {
        zl[u] = true;
    }
    while let Some(u) = stack.pop() {
        for &r in &ladj[u] {
            if !zr[r] && match_l[u] != r {
                zr[r] = true;
                let w = match_r[r];
                if w != usize::MAX && !zl[w] {
                    zl[w] = true;
                    stack.push(w);
                }
            }
        }
    }
    let mut sep = Vec::new();
    for (i, &v) in left.iter().enumerate() {
        if !zl[i] {
            sep.push(v);
        }
    }
    for (i, &v) in right.iter().enumerate() {
        if zr[i] {
            sep.push(v);
        }
    }
    sep
}

/// Vertex-separator FM refinement: `st[v]` is 0, 1 (the two sides) or 2 (the
/// separator). A separator vertex moves to side `p` at the price of pulling
/// its neighbours on the other side into the separator; moves are applied
/// greedily with hill climbing and rollback to the best state, and
/// separator vertices without neighbours on one side are released.
fn refine_separator(g: &WGraph, st: &mut [u8], passes: usize, imbalance: f64) {
    let n = g.n();
    let total: i64 = g.vw.iter().sum();
    let limit = (0.5 * (1.0 + imbalance) * total as f64).ceil() as i64;
    let mut w = [0i64; 3];
    for v in 0..n {
        w[usize::from(st[v])] += g.vw[v];
    }
    let gain = |st: &[u8], v: usize, p: u8| -> i64 {
        let mut s = g.vw[v];
        for k in g.xadj[v]..g.xadj[v + 1] {
            let u = g.adj[k];
            if st[u] == 1 - p {
                s -= g.vw[u];
            }
        }
        s
    };
    for _ in 0..passes {
        let mut heap = std::collections::BinaryHeap::new();
        for v in 0..n {
            if st[v] == 2 {
                for p in 0..2u8 {
                    heap.push((gain(st, v, p), v, p));
                }
            }
        }
        let mut locked = vec![false; n];
        let mut log: Vec<(usize, u8)> = Vec::new();
        let mut best_w2 = w[2];
        let mut best_len = 0usize;
        let mut stall = 0usize;
        while let Some((gv, v, p)) = heap.pop() {
            if st[v] != 2 || locked[v] {
                continue;
            }
            let g2 = gain(st, v, p);
            if g2 != gv {
                heap.push((g2, v, p));
                continue;
            }
            if w[usize::from(p)] + g.vw[v] > limit {
                continue;
            }
            log.push((v, 2));
            st[v] = p;
            w[usize::from(p)] += g.vw[v];
            w[2] -= g.vw[v];
            locked[v] = true;
            for k in g.xadj[v]..g.xadj[v + 1] {
                let u = g.adj[k];
                if st[u] == 1 - p {
                    log.push((u, 1 - p));
                    st[u] = 2;
                    w[usize::from(1 - p)] -= g.vw[u];
                    w[2] += g.vw[u];
                    for q in 0..2u8 {
                        heap.push((gain(st, u, q), u, q));
                    }
                }
            }
            if w[2] < best_w2 && w[0] <= limit && w[1] <= limit {
                best_w2 = w[2];
                best_len = log.len();
                stall = 0;
            } else {
                stall += 1;
                if stall > 80 {
                    break;
                }
            }
        }
        for &(x, old) in log.iter().skip(best_len).rev() {
            let cur = usize::from(st[x]);
            w[cur] -= g.vw[x];
            st[x] = old;
            w[usize::from(old)] += g.vw[x];
        }
        if best_len == 0 {
            break;
        }
    }
    // release separator vertices that touch only one side
    for v in 0..n {
        if st[v] != 2 {
            continue;
        }
        let (mut n0, mut n1) = (false, false);
        for k in g.xadj[v]..g.xadj[v + 1] {
            match st[g.adj[k]] {
                0 => n0 = true,
                1 => n1 = true,
                _ => {}
            }
        }
        if !n0 {
            st[v] = 0;
        } else if !n1 {
            st[v] = 1;
        }
    }
}

fn induced_lists(adj: &[Vec<usize>], verts: &[usize], local: &mut [usize]) -> Vec<Vec<usize>> {
    for (i, &v) in verts.iter().enumerate() {
        local[v] = i;
    }
    let lists = verts
        .iter()
        .map(|&v| adj[v].iter().filter(|&&u| local[u] != usize::MAX).map(|&u| local[u]).collect())
        .collect();
    for &v in verts {
        local[v] = usize::MAX;
    }
    lists
}

fn components(lists: &[Vec<usize>]) -> Vec<Vec<usize>> {
    let n = lists.len();
    let mut seen = vec![false; n];
    let mut out = Vec::new();
    for s in 0..n {
        if seen[s] {
            continue;
        }
        let mut comp = vec![s];
        seen[s] = true;
        let mut h = 0;
        while h < comp.len() {
            let v = comp[h];
            h += 1;
            for &u in &lists[v] {
                if !seen[u] {
                    seen[u] = true;
                    comp.push(u);
                }
            }
        }
        out.push(comp);
    }
    out
}

/// Multilevel nested dissection ordering of the pattern `adj` (see the
/// module documentation). Returns `perm` with `perm[new] = old`.
#[must_use]
pub fn multilevel_nested_dissection_order(adj: &[Vec<usize>]) -> Vec<usize> {
    multilevel_nested_dissection_with(adj, &NdOptions::default())
}

/// As [`multilevel_nested_dissection_order`] with explicit parameters.
#[must_use]
pub fn multilevel_nested_dissection_with(adj: &[Vec<usize>], opts: &NdOptions) -> Vec<usize> {
    let n = adj.len();
    let mut perm = Vec::with_capacity(n);
    let mut local = vec![usize::MAX; n];
    let mut rng = Lcg(0x9e37_79b9_7f4a_7c15);
    let all: Vec<usize> = (0..n).collect();
    order_piece(adj, &all, &mut local, &mut rng, &mut perm, opts);
    perm
}

fn order_piece(
    adj: &[Vec<usize>],
    verts: &[usize],
    local: &mut [usize],
    rng: &mut Lcg,
    perm: &mut Vec<usize>,
    o: &NdOptions,
) {
    if verts.is_empty() {
        return;
    }
    let lists = induced_lists(adj, verts, local);
    let comps = components(&lists);
    if comps.len() > 1 {
        for c in comps {
            let vs: Vec<usize> = c.iter().map(|&i| verts[i]).collect();
            order_piece(adj, &vs, local, rng, perm, o);
        }
        return;
    }
    if verts.len() <= o.leaf {
        order_leaf_with_halo(adj, verts, local, perm);
        return;
    }
    let g = WGraph::from_lists(&lists);
    let part = bisect(&g, rng, o);
    let sep0 = vertex_separator(&g, &part);
    let mut st = part.clone();
    for &x in &sep0 {
        st[x] = 2;
    }
    refine_separator(&g, &mut st, o.fm_passes, o.imbalance);
    let sep: Vec<usize> = (0..verts.len()).filter(|&i| st[i] == 2).collect();
    let side = |p: u8| -> Vec<usize> { (0..verts.len()).filter(|&i| st[i] == p).map(|i| verts[i]).collect() };
    let (a, b) = (side(0), side(1));
    if a.is_empty() || b.is_empty() || sep.len() * 2 > verts.len() {
        // degenerate split: fall back to minimum degree for this piece
        for i in amd_order(&lists) {
            perm.push(verts[i]);
        }
        return;
    }
    order_piece(adj, &a, local, rng, perm, o);
    order_piece(adj, &b, local, rng, perm, o);
    // separator vertices: order them by minimum degree within the separator graph
    let sv: Vec<usize> = sep.iter().map(|&i| verts[i]).collect();
    let sl = induced_lists(adj, &sv, local);
    for i in amd_order(&sl) {
        perm.push(sv[i]);
    }
}

/// Orders a leaf piece by minimum degree on the piece *together with its halo*
/// (the outside neighbours, which belong to separators eliminated later): the
/// halo vertices make the boundary vertices of the piece look as expensive as
/// they are, and are dropped from the resulting order (the usual substitute
/// for constrained minimum degree).
fn order_leaf_with_halo(adj: &[Vec<usize>], verts: &[usize], local: &mut [usize], perm: &mut Vec<usize>) {
    for (i, &v) in verts.iter().enumerate() {
        local[v] = i;
    }
    let mut all: Vec<usize> = verts.to_vec();
    for &v in verts {
        for &u in &adj[v] {
            if local[u] == usize::MAX {
                local[u] = all.len();
                all.push(u);
            }
        }
    }
    for &v in &all {
        local[v] = usize::MAX;
    }
    let lists = induced_lists(adj, &all, local);
    let nleaf = verts.len();
    for i in amd_order(&lists) {
        if i < nleaf {
            perm.push(verts[i]);
        }
    }
}

/// One multilevel vertex-separator bisection of a connected pattern:
/// returns `(left, right, separator)`; no edge joins `left` and `right`.
#[must_use]
pub fn vertex_separator_bisection(adj: &[Vec<usize>], opts: &NdOptions) -> (Vec<usize>, Vec<usize>, Vec<usize>) {
    let g = WGraph::from_lists(adj);
    let mut rng = Lcg(0x9e37_79b9_7f4a_7c15);
    let part = bisect(&g, &mut rng, opts);
    let sep0 = vertex_separator(&g, &part);
    let mut st = part;
    for &x in &sep0 {
        st[x] = 2;
    }
    refine_separator(&g, &mut st, opts.fm_passes, opts.imbalance);
    let pick = |p: u8| -> Vec<usize> { (0..adj.len()).filter(|&i| st[i] == p).collect() };
    (pick(0), pick(1), pick(2))
}

/// Best of approximate minimum degree and multilevel nested dissection (two
/// leaf sizes), judged by the exact symbolic Cholesky fill
/// ([`cholesky_fill`]); never worse than [`amd_order`] and, on large
/// two- and three-dimensional meshes, usually better.
#[must_use]
pub fn hybrid_nested_dissection_order(adj: &[Vec<usize>]) -> Vec<usize> {
    let mut best = amd_order(adj);
    let mut best_fill = cholesky_fill(adj, &best);
    for leaf in [2000usize, 400] {
        let cand = multilevel_nested_dissection_with(adj, &NdOptions { leaf, ..NdOptions::default() });
        let f = cholesky_fill(adj, &cand);
        if f < best_fill {
            best_fill = f;
            best = cand;
        }
    }
    best
}
