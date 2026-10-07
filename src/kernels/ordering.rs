//! # Fill-reducing orderings for sparse symmetric patterns
//!
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
//!
//! All routines take the adjacency lists of a symmetric pattern without self
//! loops (`adj[i]` lists the neighbours of `i`, each edge in both lists) and
//! return `perm` with `perm[new] = old`.
#![allow(
    clippy::missing_const_for_fn,
    clippy::too_many_lines,
    clippy::cognitive_complexity,
    clippy::many_single_char_names,
    clippy::needless_range_loop,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::similar_names,
    clippy::while_float
)]

const NIL: usize = usize::MAX;

/// Variable states of the quotient graph.
const VAR: u8 = 0;
const ELEM: u8 = 1;
const DEAD: u8 = 2;
const DENSE: u8 = 3;

/// Approximate minimum degree ordering of the pattern `adj`.
///
/// Returns `perm` with `perm[new] = old` (a permutation of `0..adj.len()`).
#[must_use]
pub fn amd_order(adj: &[Vec<usize>]) -> Vec<usize> {
    let n = adj.len();
    if n == 0 {
        return Vec::new();
    }
    let mut a: Vec<Vec<usize>> = adj.to_vec();
    let mut e: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut l: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut status = vec![VAR; n];
    let mut nv = vec![1usize; n];
    let mut degree: Vec<usize> = a.iter().map(Vec::len).collect();
    let mut elem_deg = vec![0usize; n];
    let mut parent = vec![NIL; n];
    let mut mark = vec![0usize; n];
    let mut w = vec![0usize; n];
    let mut stamp = 0usize;
    let mut wbase = 0usize;
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
        a[d].clear();
    }
    let n_active = n - dense.len();
    // bucket lists of variables by degree
    let mut head = vec![NIL; n + 1];
    let mut next = vec![NIL; n];
    let mut prev = vec![NIL; n];
    let insert = |v: usize, d: usize, head: &mut [usize], next: &mut [usize], prev: &mut [usize]| {
        let d = d.min(head.len() - 1);
        next[v] = head[d];
        prev[v] = NIL;
        if head[d] != NIL {
            prev[head[d]] = v;
        }
        head[d] = v;
    };
    let remove = |v: usize, d: usize, head: &mut [usize], next: &mut [usize], prev: &mut [usize]| {
        let d = d.min(head.len() - 1);
        if prev[v] != NIL {
            next[prev[v]] = next[v];
        } else {
            head[d] = next[v];
        }
        if next[v] != NIL {
            prev[next[v]] = prev[v];
        }
    };
    for v in 0..n {
        if status[v] == VAR {
            insert(v, degree[v], &mut head, &mut next, &mut prev);
        }
    }
    let mut mindeg = 0usize;
    let mut nel = 0usize;
    let mut elim: Vec<usize> = Vec::with_capacity(n);
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
        let mut lp: Vec<usize> = Vec::new();
        let ap = std::mem::take(&mut a[p]);
        for &j in &ap {
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
        elem_deg[p] = lp_deg;
        // |Le \ Lp| for the elements adjacent to Lp
        wbase += n + 2;
        for &i in &lp {
            for &el in &e[i] {
                if status[el] == ELEM && el != p {
                    if w[el] < wbase {
                        w[el] = wbase + elem_deg[el];
                    }
                    w[el] -= nv[i];
                }
            }
        }
        // update the variables of Lp
        let mut hashes: Vec<(usize, usize)> = Vec::with_capacity(lp.len());
        let mut mass: Vec<usize> = Vec::new();
        for &i in &lp {
            let mut ext = 0usize;
            let mut hsh = 0usize;
            let old_e = std::mem::take(&mut e[i]);
            let mut new_e: Vec<usize> = Vec::with_capacity(old_e.len() + 1);
            for &el in &old_e {
                if status[el] != ELEM || el == p {
                    continue;
                }
                let rem = w[el] - wbase;
                if rem > 0 {
                    ext += rem;
                    hsh = hsh.wrapping_add(el);
                    new_e.push(el);
                } else {
                    // Le is contained in Lp: aggressive absorption
                    status[el] = DEAD;
                    l[el] = Vec::new();
                }
            }
            let old_a = std::mem::take(&mut a[i]);
            let mut new_a: Vec<usize> = Vec::with_capacity(old_a.len());
            for &j in &old_a {
                if status[j] == VAR && mark[j] != stamp {
                    ext += nv[j];
                    hsh = hsh.wrapping_add(j);
                    new_a.push(j);
                }
            }
            if new_a.is_empty() && new_e.is_empty() {
                // only the pivot element remains: mass elimination
                mass.push(i);
                a[i] = new_a;
                e[i] = new_e;
                continue;
            }
            new_e.push(p);
            hsh = hsh.wrapping_add(p);
            let bound_n = n_active.saturating_sub(nel + nvpiv + nv[i]);
            let mut d = ext + (lp_deg - nv[i]);
            d = d.min(bound_n);
            d = d.min(degree[i] + lp_deg - nv[i]);
            degree[i] = d;
            hashes.push((hsh, i));
            a[i] = new_a;
            e[i] = new_e;
        }
        for &i in &mass {
            status[i] = DEAD;
            parent[i] = p;
            nvpiv += nv[i];
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
        let live: Vec<usize> = lp.iter().copied().filter(|&i| status[i] == VAR).collect();
        elem_deg[p] = live.iter().map(|&i| nv[i]).sum();
        for &i in &live {
            insert(i, degree[i], &mut head, &mut next, &mut prev);
            mindeg = mindeg.min(degree[i].min(n));
        }
        l[p] = live;
        elim.push(p);
    }
    // numbering: pivot, then everything merged into it, depth first
    let mut children: Vec<Vec<usize>> = vec![Vec::new(); n];
    for v in 0..n {
        if parent[v] != NIL {
            children[parent[v]].push(v);
        }
    }
    let mut perm = Vec::with_capacity(n);
    let mut stack: Vec<usize> = Vec::new();
    for &p in &elim {
        stack.push(p);
        while let Some(v) = stack.pop() {
            perm.push(v);
            for &c in children[v].iter().rev() {
                stack.push(c);
            }
        }
    }
    perm.extend(dense);
    perm
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
    // iterative: (vertex set, emit-after separator)
    enum Task {
        Part(Vec<usize>),
        Emit(Vec<usize>),
    }
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
            Task::Emit(s) => {
                perm.extend(s);
                continue;
            }
            Task::Part(v) => v,
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
                std::cmp::Ordering::Less => pa.push(v),
                std::cmp::Ordering::Greater => pb.push(v),
                std::cmp::Ordering::Equal => {
                    if adj[v].iter().any(|&u| inpart[u] == tag && level[u] > best) {
                        sep.push(v);
                    } else {
                        pa.push(v);
                    }
                }
            }
        }
        tasks.push(Task::Emit(sep));
        tasks.push(Task::Part(pb));
        tasks.push(Task::Part(pa));
    }
    perm
}

/// Subgraph induced by `verts` with local numbering `0..verts.len()`.
fn induced(adj: &[Vec<usize>], verts: &[usize], local: &mut [usize]) -> Vec<Vec<usize>> {
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
pub fn cholesky_fill(adj: &[Vec<usize>], perm: &[usize]) -> usize {
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
