//! Fill-reducing orderings: approximate minimum degree and nested dissection.
//!
//! # Fill-reducing orderings for sparse symmetric patterns
//!
//! * [`Csr::rcm_ordering`](super::krylov::Csr::rcm_ordering): reverse
//!   Cuthill-McKee with a pseudo-peripheral start per component;
//! * [`amd_order`]: approximate minimum degree on the quotient graph
//!   (Amestoy, Davis and Duff). Eliminated vertices become *elements* that
//!   absorb the elements adjacent to them (including aggressive absorption),
//!   so no explicit clique is ever formed and the memory stays proportional
//!   to the number of nonzeros; degrees are the cheap external-degree upper
//!   bounds, indistinguishable variables are merged into supervariables
//!   (hash-based detection), vertices whose neighbourhood collapses to the
//!   pivot are mass-eliminated, and dense rows are moved to the end. Degree
//!   lists are bucket lists, so selecting the next pivot is `O(1)`
//!   amortised; the method scales to millions of vertices;
//! * [`nested_dissection_order`]: recursive graph bisection by BFS level
//!   structures from pseudo-peripheral vertices with separator refinement;
//!   small pieces are ordered with [`amd_order`]; separators are numbered
//!   last;
//! * [`cholesky_fill`]: symbolic Cholesky factor size for a given ordering.
//! * [`ordering_metrics`]: half-bandwidth and lower profile for a given
//!   ordering.
//!
//! The ordering routines take adjacency lists for a symmetric pattern and
//! return `perm` with `perm[new] = old`. `adj[i]` lists the neighbours of `i`.
//!
//! # Limitations and coverage
//!
//! * AMD, nested dissection and Cholesky fill assume every edge is listed
//!   in both directions with in-range vertex indices. AMD ignores duplicate
//!   neighbours and diagonal entries. For bandwidth and profile,
//!   [`ordering_metrics`] checks indices and treats each edge as undirected.
//! * [`amd_order`] is approximate minimum degree without a supernode
//!   post-processing pass; [`nested_dissection_order`] bisects by BFS level
//!   structures with separator refinement, not by a multilevel or spectral
//!   partition, so on very large irregular graphs the fill can be larger than
//!   a state-of-the-art ordering.
//! * [`cholesky_fill`] counts the symbolic factor of a symmetric pattern.
//!   The sparse LU in `krylov` applies a symmetric ordering before numerical
//!   row pivoting.
//! * Memory is proportional to the number of nonzeros for AMD; the nested
//!   dissection recursion stores each separator explicitly, and the routines
//!   are single-threaded and deterministic.

use super::numeric_common::NumericError;

pub(super) mod rcm;

const NIL: usize = usize::MAX;

/// Variable states of the quotient graph.
const VAR: u8 = 0;
const ELEM: u8 = 1;
const DEAD: u8 = 2;
const DENSE: u8 = 3;

/// Approximate minimum degree ordering of the pattern `adj`.
///
/// Returns `perm` with `perm[new] = old` (a permutation of `0..adj.len()`).
/// The pattern must be symmetric. Duplicate entries and self loops do not
/// affect the ordering.
///
/// # Panics
/// A neighbour index is outside `0..adj.len()`.
#[must_use]
pub fn amd_order(adj: &[Vec<usize>]) -> Vec<usize> {
    let n = adj.len();
    if adj.iter().all(Vec::is_empty) {
        return (0..n).rev().collect();
    }
    let mut mark = vec![0usize; n];
    let mut a: Vec<Vec<usize>> = adj
        .iter()
        .enumerate()
        .map(|(v, row)| {
            let tag = v + 1;
            mark[v] = tag;
            let mut neighbors = Vec::with_capacity(row.len().min(n));
            for &u in row {
                if mark[u] != tag {
                    mark[u] = tag;
                    neighbors.push(u);
                }
            }
            neighbors
        })
        .collect();
    let mut e: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut l: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut status = vec![VAR; n];
    let mut nv = vec![1usize; n];
    let mut degree: Vec<usize> = a.iter().map(Vec::len).collect();
    let mut elem_deg = vec![0usize; n];
    let mut parent = vec![NIL; n];
    let mut w = vec![0usize; n];
    // Initial row tags are 1..=n. Pivot tags start above them.
    let mut stamp = n;
    // dense rows are ordered last
    let dense_limit = ((10.0 * (n as f64).sqrt()) as usize).max(16);
    let mut dense: Vec<usize> = Vec::new();
    for i in 0..n {
        if degree[i] > dense_limit {
            status[i] = DENSE;
            dense.push(i);
        }
    }
    for i in 0..n {
        if status[i] == VAR {
            a[i].retain(|&j| status[j] == VAR);
            degree[i] = a[i].len();
        }
    }
    for &d in &dense {
        a[d] = Vec::new();
    }
    let n_active = n - dense.len();
    // bucket lists of variables by degree
    let mut head = vec![NIL; n + 1];
    let mut next = vec![NIL; n];
    let mut prev = vec![NIL; n];
    let insert =
        |v: usize, d: usize, head: &mut [usize], next: &mut [usize], prev: &mut [usize]| {
            let d = d.min(head.len() - 1);
            next[v] = head[d];
            prev[v] = NIL;
            if head[d] != NIL {
                prev[head[d]] = v;
            }
            head[d] = v;
        };
    let remove =
        |v: usize, d: usize, head: &mut [usize], next: &mut [usize], prev: &mut [usize]| {
            let d = d.min(head.len() - 1);
            if prev[v] == NIL {
                head[d] = next[v];
            } else {
                next[prev[v]] = next[v];
            }
            if next[v] != NIL {
                prev[next[v]] = prev[v];
            }
        };
    let mut elim = Vec::with_capacity(n);
    for v in 0..n {
        if status[v] == VAR {
            if degree[v] == 0 {
                status[v] = DEAD;
                elim.push(v);
            } else {
                insert(v, degree[v], &mut head, &mut next, &mut prev);
            }
        }
    }
    // A degree-zero vertex starts with an empty quotient neighborhood.
    // Preserve the descending tie order of the initial degree-zero bucket.
    elim.reverse();
    let mut mindeg = 0usize;
    let mut nel = elim.len();
    let mut lp = Vec::new();
    let mut hashes = Vec::new();
    let mut mass = Vec::new();
    while nel < n_active {
        while mindeg < n && head[mindeg] == NIL {
            mindeg += 1;
        }
        let p = head[mindeg];
        if p == NIL {
            break;
        }
        remove(p, degree[p], &mut head, &mut next, &mut prev);
        let mut nvpiv = nv[p];
        stamp += 1;
        mark[p] = stamp;
        // Lp = (A_p + union of L_e, e in E_p) minus p
        lp.clear();
        let mut element = std::mem::take(&mut a[p]);
        for &j in &element {
            if status[j] == VAR && mark[j] != stamp {
                mark[j] = stamp;
                lp.push(j);
            }
        }
        let ep = std::mem::take(&mut e[p]);
        for &el in &ep {
            if status[el] == ELEM {
                let le = std::mem::take(&mut l[el]);
                for &j in &le {
                    if status[j] == VAR && mark[j] != stamp {
                        mark[j] = stamp;
                        lp.push(j);
                    }
                }
                status[el] = DEAD;
            }
        }
        status[p] = ELEM;
        let mut lp_deg = 0usize;
        for &j in &lp {
            lp_deg += nv[j];
            remove(j, degree[j], &mut head, &mut next, &mut prev);
        }
        // |Le \ Lp| for the elements adjacent to Lp
        for &i in &lp {
            for &el in &e[i] {
                if status[el] == ELEM && el != p {
                    // Variables in Lp and adjacent elements have disjoint
                    // slots, so their marks share the current pivot tag.
                    if mark[el] != stamp {
                        mark[el] = stamp;
                        w[el] = elem_deg[el];
                    }
                    w[el] -= nv[i];
                }
            }
        }
        // update the variables of Lp
        hashes.clear();
        mass.clear();
        for &i in &lp {
            let mut ext = 0usize;
            let mut hsh = 0usize;
            e[i].retain(|&el| {
                if status[el] != ELEM || el == p {
                    return false;
                }
                let rem = w[el];
                if rem > 0 {
                    ext += rem;
                    hsh = hsh.wrapping_add(el);
                    true
                } else {
                    // Le is contained in Lp: aggressive absorption
                    status[el] = DEAD;
                    l[el] = Vec::new();
                    false
                }
            });
            a[i].retain(|&j| {
                let keep = status[j] == VAR && mark[j] != stamp;
                if keep {
                    ext += nv[j];
                    hsh = hsh.wrapping_add(j);
                }
                keep
            });
            if a[i].is_empty() && e[i].is_empty() {
                // only the pivot element remains: mass elimination
                mass.push(i);
                continue;
            }
            e[i].push(p);
            hsh = hsh.wrapping_add(p);
            let bound_n = n_active.saturating_sub(nel + nvpiv + nv[i]);
            let mut d = ext + (lp_deg - nv[i]);
            d = d.min(bound_n);
            d = d.min(degree[i] + lp_deg - nv[i]);
            degree[i] = d;
            hashes.push((hsh, i));
        }
        for &i in &mass {
            status[i] = DEAD;
            parent[i] = p;
            nvpiv += nv[i];
            a[i] = Vec::new();
            e[i] = Vec::new();
        }
        // supervariable detection among the surviving members of Lp
        hashes.sort_unstable();
        let mut s = 0;
        while s < hashes.len() {
            let mut t = s + 1;
            while t < hashes.len() && hashes[t].0 == hashes[s].0 {
                t += 1;
            }
            if t - s > 1 {
                for k in s..t {
                    let i = hashes[k].1;
                    a[i].sort_unstable();
                    e[i].sort_unstable();
                }
                for x in s..t {
                    let i = hashes[x].1;
                    if status[i] != VAR {
                        continue;
                    }
                    for y in x + 1..t {
                        let j = hashes[y].1;
                        if status[j] != VAR || a[i] != a[j] || e[i] != e[j] {
                            continue;
                        }
                        nv[i] += nv[j];
                        degree[i] = degree[i].saturating_sub(nv[j]);
                        status[j] = DEAD;
                        parent[j] = i;
                        a[j] = Vec::new();
                        e[j] = Vec::new();
                    }
                }
            }
            s = t;
        }
        nel += nvpiv;
        nv[p] = nvpiv;
        // finalise element p and reinsert the variables
        element.clear();
        elem_deg[p] = 0;
        for &i in lp.iter().filter(|&&i| status[i] == VAR) {
            element.push(i);
            elem_deg[p] += nv[i];
            insert(i, degree[i], &mut head, &mut next, &mut prev);
            mindeg = mindeg.min(degree[i].min(n));
        }
        l[p] = element;
        elim.push(p);
    }
    // numbering: pivot, then everything merged into it, depth first
    let mut first_child = head;
    let mut sibling = next;
    first_child.fill(NIL);
    sibling.fill(NIL);
    for (v, &p) in parent.iter().enumerate() {
        if p != NIL {
            sibling[v] = first_child[p];
            first_child[p] = v;
        }
    }
    let mut perm = Vec::with_capacity(n);
    let mut stack: Vec<usize> = Vec::new();
    for &p in &elim {
        stack.push(p);
        while let Some(v) = stack.pop() {
            perm.push(v);
            let mut child = first_child[v];
            while child != NIL {
                stack.push(child);
                child = sibling[child];
            }
        }
    }
    perm.extend(dense);
    perm
}

/// Work item of the iterative nested dissection: a piece still to split, or a separator to emit.
enum Task {
    Part(Vec<usize>),
    Emit(Vec<usize>),
}

/// Nested dissection ordering by recursive level-structure bisection.
///
/// Pieces of at most 300 vertices are ordered by [`amd_order`]; larger
/// connected pieces are split by a vertex separator taken from a BFS level
/// structure rooted at a pseudo-peripheral vertex (the level near the middle
/// with the best size/balance trade-off, pruned of vertices without
/// neighbours on the far side); separators are numbered after both halves.
#[must_use]
pub fn nested_dissection_order(adj: &[Vec<usize>]) -> Vec<usize> {
    let n = adj.len();
    let mut perm = Vec::with_capacity(n);
    let mut local = vec![NIL; n];
    let mut seen = vec![false; n];
    let mut tasks: Vec<Task> = Vec::new();
    // components of the whole graph
    let mut comps: Vec<Vec<usize>> = Vec::new();
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
            for &u in &adj[v] {
                if !seen[u] {
                    seen[u] = true;
                    comp.push(u);
                }
            }
        }
        comps.push(comp);
    }
    for c in comps.into_iter().rev() {
        tasks.push(Task::Part(c));
    }
    let mut level = vec![NIL; n];
    let mut inpart = vec![0usize; n];
    let mut tag = 0usize;
    while let Some(task) = tasks.pop() {
        let verts = match task {
            | Task::Emit(s) => {
                perm.extend(s);
                continue;
            },
            | Task::Part(v) => v,
        };
        if verts.len() <= 300 {
            let sub = induced(adj, &verts, &mut local);
            for i in amd_order(&sub) {
                perm.push(verts[i]);
            }
            continue;
        }
        tag += 1;
        for &v in &verts {
            inpart[v] = tag;
        }
        // pseudo-peripheral root by repeated BFS
        let bfs = |root: usize, level: &mut [usize]| -> Vec<usize> {
            for &v in &verts {
                level[v] = NIL;
            }
            let mut q = vec![root];
            level[root] = 0;
            let mut h = 0;
            while h < q.len() {
                let v = q[h];
                h += 1;
                for &u in &adj[v] {
                    if inpart[u] == tag && level[u] == NIL {
                        level[u] = level[v] + 1;
                        q.push(u);
                    }
                }
            }
            q
        };
        let mut root = verts[0];
        let mut q = bfs(root, &mut level);
        let mut ecc = level[*q.last().unwrap_or(&root)];
        for _ in 0..8 {
            let far = *q.last().unwrap_or(&root);
            let q2 = bfs(far, &mut level);
            let e2 = level[*q2.last().unwrap_or(&far)];
            if e2 > ecc {
                ecc = e2;
                root = far;
                q = q2;
            } else {
                q = bfs(root, &mut level);
                break;
            }
        }
        if q.len() < verts.len() {
            // disconnected piece (can only happen after separator removal)
            let rest: Vec<usize> = verts.iter().copied().filter(|&v| level[v] == NIL).collect();
            tasks.push(Task::Part(rest));
            tasks.push(Task::Part(q));
            continue;
        }
        let depth = ecc;
        if depth < 2 {
            // clique-like: no useful separator
            let sub = induced(adj, &verts, &mut local);
            for i in amd_order(&sub) {
                perm.push(verts[i]);
            }
            continue;
        }
        let mut counts = vec![0usize; depth + 1];
        for &v in &verts {
            counts[level[v]] += 1;
        }
        let total = verts.len();
        let mut best = 1usize;
        let mut best_cost = f64::INFINITY;
        let mut before = counts[0];
        for m in 1..depth {
            let after = total - before - counts[m];
            let imbalance = before.abs_diff(after) as f64 / total as f64;
            let cost = counts[m] as f64 * (1.0 + 4.0 * imbalance * imbalance) + 1e-3 * imbalance;
            if cost < best_cost {
                best_cost = cost;
                best = m;
            }
            before += counts[m];
        }
        let mut sep: Vec<usize> = Vec::new();
        let (mut pa, mut pb): (Vec<usize>, Vec<usize>) = (Vec::new(), Vec::new());
        for &v in &verts {
            match level[v].cmp(&best) {
                | std::cmp::Ordering::Less => pa.push(v),
                | std::cmp::Ordering::Greater => pb.push(v),
                | std::cmp::Ordering::Equal => {
                    if adj[v].iter().any(|&u| inpart[u] == tag && level[u] > best) {
                        sep.push(v);
                    } else {
                        pa.push(v);
                    }
                },
            }
        }
        tasks.push(Task::Emit(sep));
        tasks.push(Task::Part(pb));
        tasks.push(Task::Part(pa));
    }
    perm
}

/// Subgraph induced by `verts` with local numbering `0..verts.len()`.
fn induced(
    adj: &[Vec<usize>],
    verts: &[usize],
    local: &mut [usize],
) -> Vec<Vec<usize>> {
    for (i, &v) in verts.iter().enumerate() {
        local[v] = i;
    }
    let sub: Vec<Vec<usize>> = verts
        .iter()
        .map(|&v| {
            adj[v]
                .iter()
                .filter(|&&u| local[u] != NIL && verts.get(local[u]) == Some(&u))
                .map(|&u| local[u])
                .collect()
        })
        .collect();
    for &v in verts {
        local[v] = NIL;
    }
    sub
}

/// Number of nonzeros of the Cholesky factor `L` (including the diagonal)
/// of the pattern `adj` ordered by `perm` (`perm[new] = old`), by symbolic
/// factorisation along the elimination tree.
#[must_use]
pub fn cholesky_fill(
    adj: &[Vec<usize>],
    perm: &[usize],
) -> usize {
    let n = adj.len();
    let mut inv = vec![0usize; n];
    for (new, &old) in perm.iter().enumerate() {
        inv[old] = new;
    }
    let mut structs: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut children: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut total = 0usize;
    let mut mark = vec![NIL; n];
    for j in 0..n {
        let mut s: Vec<usize> = Vec::new();
        for &u in &adj[perm[j]] {
            let k = inv[u];
            if k > j && mark[k] != j {
                mark[k] = j;
                s.push(k);
            }
        }
        for &c in &children[j] {
            let sc = std::mem::take(&mut structs[c]);
            for &k in &sc {
                if k > j && mark[k] != j {
                    mark[k] = j;
                    s.push(k);
                }
            }
        }
        total += s.len() + 1;
        if let Some(&par) = s.iter().min() {
            children[par].push(j);
        }
        structs[j] = s;
    }
    total
}

/// Structural bandwidth and lower profile after a symmetric permutation.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct OrderingMetrics {
    /// Largest distance from the diagonal of any structural entry.
    pub half_bandwidth: usize,
    /// Sum of the distances to the first lower-triangular entry of each row.
    /// The diagonal contributes zero.
    pub profile: usize,
}

/// Measures the bandwidth and profile of a permuted sparse pattern.
///
/// `perm[new] = old` names the original vertex at each reordered position.
/// Edges may be listed in either or both directions. Duplicate and diagonal
/// entries leave both metrics unchanged. The work is linear in the number
/// of vertices and listed edges, with linear storage in the vertex count.
///
/// # Errors
/// The permutation is malformed, a neighbour is out of range, or the
/// accumulated profile does not fit a `usize`.
pub fn ordering_metrics(
    adj: &[Vec<usize>],
    perm: &[usize],
) -> Result<OrderingMetrics, NumericError> {
    measure_ordering(adj.iter().map(Vec::as_slice), perm)
}

pub(crate) fn measure_ordering<'a>(
    rows: impl ExactSizeIterator<Item = &'a [usize]>,
    perm: &[usize],
) -> Result<OrderingMetrics, NumericError> {
    let vertices = rows.len();
    if perm.len() != vertices {
        return Err(NumericError::ShapeMismatch(format!(
            "ordering has {} positions for {vertices} vertices",
            perm.len()
        )));
    }
    let mut positions = vec![NIL; vertices];
    for (new, &old) in perm.iter().enumerate() {
        match positions.get_mut(old) {
            | Some(position) if *position == NIL => *position = new,
            | Some(_) => {
                return Err(NumericError::Domain(format!(
                    "duplicate vertex {old} in ordering"
                )));
            },
            | None => {
                return Err(NumericError::Domain(format!(
                    "ordering vertex {old} is out of range"
                )));
            },
        }
    }
    let mut row_width = vec![0usize; vertices];
    for (row, neighbors) in rows.enumerate() {
        let i = positions[row];
        for &neighbor in neighbors {
            let Some(&j) = positions.get(neighbor) else {
                return Err(NumericError::Domain(format!(
                    "neighbour {neighbor} in row {row} is out of range",
                )));
            };
            let width = &mut row_width[i.max(j)];
            *width = (*width).max(i.abs_diff(j));
        }
    }
    let mut metrics = OrderingMetrics::default();
    for width in row_width {
        metrics.half_bandwidth = metrics.half_bandwidth.max(width);
        metrics.profile = metrics
            .profile
            .checked_add(width)
            .ok_or_else(|| NumericError::Domain("ordering profile exceeds usize".to_owned()))?;
    }
    Ok(metrics)
}
