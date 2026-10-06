//! # Exact and general algebraic topology
//!
//! Integer homology with torsion from the Smith normal form of the boundary
//! matrices, cohomology (universal coefficients), homology over `Z/p`,
//! relative homology, standard triangulated spaces and the operations that
//! combine them, fundamental group presentations of 2-complexes, simplicial
//! maps and the maps they induce in homology, persistent homology by the
//! standard column reduction, the bottleneck distance of persistence
//! diagrams and cubical complexes of binary images.
//!
//! Conventions are those of [`super::topology`]: a simplex is a sorted
//! vector of vertex indices and a complex is a list of simplices closed
//! under taking faces.

use std::collections::BTreeSet;
use std::collections::HashMap;
use std::collections::HashSet;

use num_bigint::BigInt;
use num_traits::One;
use num_traits::Signed;
use num_traits::Zero;

use super::qlinalg as ql;
use super::qlinalg::Q;
use super::topology::boundary_matrix;
use super::topology::close_complex;
use super::topology::simplices_of_dim;
use super::topology::Simplex;

// ----------------------------------------------------------------------
// Smith normal form
// ----------------------------------------------------------------------

/// The nonzero invariant factors `d_1 | d_2 | ...` of an integer matrix
/// (its Smith normal form), exactly. Their number is the rank.
#[must_use]
pub fn smith_invariants(m: &[Vec<i64>]) -> Vec<BigInt> {
    let mut a: Vec<Vec<BigInt>> = m.iter().map(|r| r.iter().map(|&x| BigInt::from(x)).collect()).collect();
    let nr = a.len();
    let nc = a.first().map_or(0, Vec::len);
    let mut row_alive = vec![true; nr];
    let mut col_alive = vec![true; nc];
    let mut units = 0_usize;
    // Phase 1: eliminate entries equal to +-1 (they contribute invariant 1).
    for r in 0..nr {
        let Some(c) = (0..nc).find(|&c| col_alive[c] && a[r][c].abs().is_one()) else {
            continue;
        };
        let pivot_row = a[r].clone();
        let sign = pivot_row[c].clone();
        for (i, row) in a.iter_mut().enumerate() {
            if i == r || !row_alive[i] || row[c].is_zero() {
                continue;
            }
            let f = &row[c] * &sign;
            for (x, y) in row.iter_mut().zip(&pivot_row) {
                if !y.is_zero() {
                    *x -= &f * y;
                }
            }
        }
        row_alive[r] = false;
        col_alive[c] = false;
        units += 1;
    }
    let rest: Vec<Vec<BigInt>> = (0..nr)
        .filter(|&r| row_alive[r])
        .map(|r| (0..nc).filter(|&c| col_alive[c]).map(|c| a[r][c].clone()).collect())
        .collect();
    let mut out = vec![BigInt::one(); units];
    out.extend(smith_dense(rest));
    out
}

fn smith_dense(mut b: Vec<Vec<BigInt>>) -> Vec<BigInt> {
    let nr = b.len();
    let nc = b.first().map_or(0, Vec::len);
    let mut out = Vec::new();
    let mut t = 0;
    'pivots: while t < nr.min(nc) {
        let mut best: Option<(usize, usize)> = None;
        for i in t..nr {
            for j in t..nc {
                if !b[i][j].is_zero() && best.is_none_or(|(bi, bj)| b[i][j].abs() < b[bi][bj].abs()) {
                    best = Some((i, j));
                }
            }
        }
        let Some((pi, pj)) = best else { break };
        b.swap(t, pi);
        for row in &mut b {
            row.swap(t, pj);
        }
        for i in t + 1..nr {
            if b[i][t].is_zero() {
                continue;
            }
            let f = &b[i][t] / &b[t][t];
            let pivot_row = b[t].clone();
            for (x, y) in b[i].iter_mut().zip(&pivot_row).skip(t) {
                *x -= &f * y;
            }
            if !b[i][t].is_zero() {
                b.swap(t, i);
                continue 'pivots;
            }
        }
        for j in t + 1..nc {
            if b[t][j].is_zero() {
                continue;
            }
            let f = &b[t][j] / &b[t][t];
            for row in b.iter_mut().skip(t) {
                let d = &f * &row[t];
                row[j] -= d;
            }
            if !b[t][j].is_zero() {
                for row in &mut b {
                    row.swap(t, j);
                }
                continue 'pivots;
            }
        }
        let p = b[t][t].clone();
        let bad = (t + 1..nr).find(|&i| (t + 1..nc).any(|j| !(&b[i][j] % &p).is_zero()));
        if let Some(i) = bad {
            let row = b[i].clone();
            for (x, y) in b[t].iter_mut().zip(&row) {
                *x += y;
            }
            continue 'pivots;
        }
        out.push(b[t][t].abs());
        t += 1;
    }
    out
}

/// The rank of an integer matrix modulo the prime `p`.
#[must_use]
pub fn rank_mod_p(
    m: &[Vec<i64>],
    p: u64,
) -> usize {
    let pm = i128::from(p);
    let mut a: Vec<Vec<u64>> = m
        .iter()
        .map(|r| r.iter().map(|&x| u64::try_from(i128::from(x).rem_euclid(pm)).unwrap_or(0)).collect())
        .collect();
    let cols = a.first().map_or(0, Vec::len);
    let mut rank = 0;
    for c in 0..cols {
        let Some(piv) = (rank..a.len()).find(|&r| a[r][c] != 0) else { continue };
        a.swap(rank, piv);
        let inv = mod_pow(a[rank][c], p - 2, p);
        let pivot_row: Vec<u64> = a[rank].iter().map(|&x| mul_mod(x, inv, p)).collect();
        for (r, row) in a.iter_mut().enumerate().skip(rank + 1) {
            let f = row[c];
            if r == rank || f == 0 {
                continue;
            }
            for (x, &y) in row.iter_mut().zip(&pivot_row).skip(c) {
                *x = (*x + p - mul_mod(f, y, p)) % p;
            }
        }
        rank += 1;
        if rank == a.len() {
            break;
        }
    }
    rank
}

fn mul_mod(
    a: u64,
    b: u64,
    p: u64,
) -> u64 {
    u64::try_from(u128::from(a) * u128::from(b) % u128::from(p)).unwrap_or(0)
}

fn mod_pow(
    mut base: u64,
    mut exp: u64,
    p: u64,
) -> u64 {
    let mut acc = 1 % p;
    base %= p;
    while exp > 0 {
        if exp & 1 == 1 {
            acc = mul_mod(acc, base, p);
        }
        base = mul_mod(base, base, p);
        exp >>= 1;
    }
    acc
}

/// Whether `n` is prime (trial division; `n` is small).
#[must_use]
pub fn is_prime(n: u64) -> bool {
    n >= 2 && (2..=n.isqrt()).all(|d| !n.is_multiple_of(d))
}

// ----------------------------------------------------------------------
// Chain complexes
// ----------------------------------------------------------------------

/// A finite chain complex of free abelian groups: the rank of each `C_k`
/// and the integer matrices of the boundaries.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ChainComplex {
    /// `dims[k]` is the rank of `C_k`.
    pub dims: Vec<usize>,
    /// `d[k]` is the matrix of `C_k -> C_{k-1}` (`dims[k-1]` rows,
    /// `dims[k]` columns); `d[0]` is empty.
    pub d: Vec<Vec<Vec<i64>>>,
}

/// A homology group `Z^rank + Z/t_1 + ...`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Group {
    /// The free rank.
    pub rank: usize,
    /// The torsion coefficients (each greater than one, dividing the next).
    pub torsion: Vec<BigInt>,
}

/// The chain complex with basis `basis[k]` in degree `k` (lists of
/// `k`-simplices). Faces missing from `basis[k-1]` are dropped, which gives
/// the chain complex of a quotient.
#[must_use]
pub fn chain_complex(basis: &[Vec<Simplex>]) -> ChainComplex {
    let dims = basis.iter().map(Vec::len).collect();
    let d = basis.iter().enumerate().map(|(k, cols)| if k == 0 { Vec::new() } else { boundary_matrix(&basis[k - 1], cols) }).collect();
    ChainComplex { dims, d }
}

/// The chain complex of a closed simplicial complex.
#[must_use]
pub fn complex_chains(complex: &[Simplex]) -> ChainComplex {
    let top = complex.iter().map(Vec::len).max().unwrap_or(0);
    let basis: Vec<Vec<Simplex>> = (0..top).map(|k| simplices_of_dim(complex, k)).collect();
    chain_complex(&basis)
}

/// The chain complex of the pair `(K, L)`: simplices of `K` not in `L`.
#[must_use]
pub fn relative_chains(
    complex: &[Simplex],
    sub: &[Simplex],
) -> ChainComplex {
    let in_sub: HashSet<&Simplex> = sub.iter().collect();
    let top = complex.iter().map(Vec::len).max().unwrap_or(0);
    let basis: Vec<Vec<Simplex>> = (0..top).map(|k| simplices_of_dim(complex, k).into_iter().filter(|s| !in_sub.contains(s)).collect()).collect();
    chain_complex(&basis)
}

impl ChainComplex {
    fn dim_at(
        &self,
        k: usize,
    ) -> usize {
        self.dims.get(k).copied().unwrap_or(0)
    }

    fn invariants(
        &self,
        k: usize,
    ) -> Vec<BigInt> {
        if k == 0 || self.dim_at(k) == 0 || self.dim_at(k - 1) == 0 {
            return Vec::new();
        }
        smith_invariants(&self.d[k])
    }

    /// `H_k(C; Z)`.
    #[must_use]
    pub fn homology(
        &self,
        k: usize,
    ) -> Group {
        let below = self.invariants(k).len();
        let above = self.invariants(k + 1);
        let rank = self.dim_at(k).saturating_sub(below).saturating_sub(above.len());
        Group { rank, torsion: above.into_iter().filter(|t| !t.is_one()).collect() }
    }

    /// `H^k(C; Z)` by the universal coefficient theorem: the free part of
    /// `H_k` and the torsion of `H_{k-1}`.
    #[must_use]
    pub fn cohomology(
        &self,
        k: usize,
    ) -> Group {
        let rank = self.homology(k).rank;
        let torsion = if k == 0 { Vec::new() } else { self.homology(k - 1).torsion };
        Group { rank, torsion }
    }

    /// `dim H_k(C; Z/p)` for a prime `p`.
    #[must_use]
    pub fn homology_mod_p(
        &self,
        k: usize,
        p: u64,
    ) -> usize {
        let rk = |j: usize| if j == 0 || self.dim_at(j) == 0 || self.dim_at(j - 1) == 0 { 0 } else { rank_mod_p(&self.d[j], p) };
        self.dim_at(k).saturating_sub(rk(k)).saturating_sub(rk(k + 1))
    }

    /// The degree of the highest nonzero chain group, if any.
    #[must_use]
    pub fn top_degree(&self) -> Option<usize> {
        self.dims.iter().rposition(|&d| d > 0)
    }
}

/// The reduced homology: `H_0` loses one free summand when the complex is
/// nonempty.
#[must_use]
pub fn reduced_homology(
    cc: &ChainComplex,
    k: usize,
) -> Group {
    let mut g = cc.homology(k);
    if k == 0 && cc.dim_at(0) > 0 {
        g.rank = g.rank.saturating_sub(1);
    }
    g
}

/// Whether every simplex of `sub` belongs to `complex`.
#[must_use]
pub fn is_subcomplex(
    sub: &[Simplex],
    complex: &[Simplex],
) -> bool {
    let all: HashSet<&Simplex> = complex.iter().collect();
    sub.iter().all(|s| all.contains(s))
}

/// The union of two complexes.
#[must_use]
pub fn union(
    a: &[Simplex],
    b: &[Simplex],
) -> Vec<Simplex> {
    let mut all = a.to_vec();
    all.extend_from_slice(b);
    close_complex(&all)
}

/// The intersection of two complexes.
#[must_use]
pub fn intersection(
    a: &[Simplex],
    b: &[Simplex],
) -> Vec<Simplex> {
    let in_b: HashSet<&Simplex> = b.iter().collect();
    close_complex(&a.iter().filter(|s| in_b.contains(s)).cloned().collect::<Vec<_>>())
}

/// The Betti numbers `b_0 ..= b_top` of a complex.
#[must_use]
pub fn betti_list(complex: &[Simplex]) -> Vec<usize> {
    let cc = complex_chains(complex);
    (0..cc.dims.len()).map(|k| cc.homology(k).rank).collect()
}

/// Checks the Mayer-Vietoris sequence of `A ∪ B` at the level of ranks.
///
/// the Euler characteristics add up, and there are nonnegative ranks of the
/// connecting maps that make every stretch of the long exact sequence
/// exact. Returns the ranks of the connecting maps `H_k(A ∪ B) -> H_{k-1}(A ∩ B)`
/// for `k = 1 ..` when the check passes.
#[must_use]
pub fn mayer_vietoris(
    a: &[Simplex],
    b: &[Simplex],
) -> Option<Vec<usize>> {
    let u = union(a, b);
    let i = intersection(a, b);
    let chi = |c: &[Simplex]| super::topology::euler_characteristic(c);
    if chi(&u) != chi(a) + chi(b) - chi(&i) {
        return None;
    }
    let (bu, ba, bb, bi) = (betti_list(&u), betti_list(a), betti_list(b), betti_list(&i));
    let top = bu.len().max(ba.len()).max(bb.len()).max(bi.len()) + 1;
    let at = |v: &[usize], k: usize| i64::try_from(v.get(k).copied().unwrap_or(0)).unwrap_or(0);
    // d_k = rank of H_k(U) -> H_{k-1}(I); sigma_k = rank of H_k(A)+H_k(B) -> H_k(U);
    // rho_k = rank of H_k(I) -> H_k(A)+H_k(B).
    let mut connecting = vec![0_i64];
    for k in 0..top {
        let d_k = connecting[k];
        let sigma = at(&bu, k) - d_k;
        let rho = at(&ba, k) + at(&bb, k) - sigma;
        let d_next = at(&bi, k) - rho;
        if sigma < 0 || rho < 0 || d_next < 0 || sigma > at(&ba, k) + at(&bb, k) {
            return None;
        }
        connecting.push(d_next);
    }
    if connecting.last().copied() != Some(0) {
        return None;
    }
    connecting.into_iter().skip(1).map(|d| usize::try_from(d).ok()).collect()
}

// ----------------------------------------------------------------------
// Standard complexes
// ----------------------------------------------------------------------

/// The full `n`-simplex on vertices `0..=n`.
#[must_use]
pub fn full_simplex(n: usize) -> Vec<Simplex> {
    close_complex(&[(0..=n).collect()])
}

/// The `n`-sphere: the boundary of the `(n+1)`-simplex.
#[must_use]
pub fn sphere(n: usize) -> Vec<Simplex> {
    let facets: Vec<Simplex> = (0..=n + 1).map(|skip| (0..=n + 1).filter(|&v| v != skip).collect()).collect();
    close_complex(&facets)
}

fn next_vertex(complex: &[Simplex]) -> usize {
    complex.iter().flatten().max().map_or(0, |&m| m + 1)
}

/// The cone over a complex with a new apex.
#[must_use]
pub fn cone(complex: &[Simplex]) -> Vec<Simplex> {
    let apex = next_vertex(complex);
    let mut all = complex.to_vec();
    for s in complex {
        let mut t = s.clone();
        t.push(apex);
        all.push(t);
    }
    all.push(vec![apex]);
    close_complex(&all)
}

/// The suspension: two cones over the complex glued along it.
#[must_use]
pub fn suspension(complex: &[Simplex]) -> Vec<Simplex> {
    let (a, b) = (next_vertex(complex), next_vertex(complex) + 1);
    let mut all = complex.to_vec();
    for s in complex {
        for apex in [a, b] {
            let mut t = s.clone();
            t.push(apex);
            all.push(t);
        }
    }
    all.push(vec![a]);
    all.push(vec![b]);
    close_complex(&all)
}

/// The join `K * L` (the vertices of `L` are shifted past those of `K`).
#[must_use]
pub fn join(
    k: &[Simplex],
    l: &[Simplex],
) -> Vec<Simplex> {
    let shift = next_vertex(k);
    let l: Vec<Simplex> = l.iter().map(|s| s.iter().map(|v| v + shift).collect()).collect();
    let mut all = k.to_vec();
    all.extend(l.iter().cloned());
    for s in k {
        for t in &l {
            let mut u = s.clone();
            u.extend(t);
            all.push(u);
        }
    }
    close_complex(&all)
}

/// The wedge sum, gluing vertex `vk` of `K` to vertex `vl` of `L`.
#[must_use]
pub fn wedge(
    k: &[Simplex],
    l: &[Simplex],
    vk: usize,
    vl: usize,
) -> Vec<Simplex> {
    let shift = next_vertex(k);
    let map = |v: usize| if v == vl { vk } else { v + shift };
    let mut all = k.to_vec();
    all.extend(l.iter().map(|s| s.iter().map(|&v| map(v)).collect::<Simplex>()));
    close_complex(&all)
}

/// The triangulated product `K x L` (ordered product of simplices: the
/// vertex `(v, w)` is numbered `index(v) * |L| + index(w)`), or `None` when
/// it would exceed `limit` simplices.
#[must_use]
pub fn product(
    k: &[Simplex],
    l: &[Simplex],
    limit: usize,
) -> Option<Vec<Simplex>> {
    let vk: Vec<usize> = simplices_of_dim(k, 0).into_iter().map(|s| s[0]).collect();
    let vl: Vec<usize> = simplices_of_dim(l, 0).into_iter().map(|s| s[0]).collect();
    let in_k: HashSet<&Simplex> = k.iter().collect();
    let in_l: HashSet<&Simplex> = l.iter().collect();
    let mut found: BTreeSet<Simplex> = BTreeSet::new();
    let mut stack: Vec<Vec<(usize, usize)>> = Vec::new();
    for a in 0..vk.len() {
        for b in 0..vl.len() {
            stack.push(vec![(a, b)]);
        }
    }
    while let Some(chain) = stack.pop() {
        let mut label: Simplex = chain.iter().map(|&(a, b)| a * vl.len() + b).collect();
        label.sort_unstable();
        if !found.insert(label) {
            continue;
        }
        if found.len() > limit {
            return None;
        }
        let &(a0, b0) = chain.last()?;
        for a in a0..vk.len() {
            for b in b0..vl.len() {
                if (a, b) == (a0, b0) {
                    continue;
                }
                let mut next = chain.clone();
                next.push((a, b));
                let mut pk: Simplex = next.iter().map(|&(a, _)| vk[a]).collect();
                pk.sort_unstable();
                pk.dedup();
                let mut pl: Simplex = next.iter().map(|&(_, b)| vl[b]).collect();
                pl.sort_unstable();
                pl.dedup();
                if in_k.contains(&pk) && in_l.contains(&pl) {
                    stack.push(next);
                }
            }
        }
    }
    Some(found.into_iter().collect())
}

/// The triangulated Klein bottle from an `m x n` grid (`m, n >= 3`): the
/// columns are identified directly, the rows with a reflection.
#[must_use]
pub fn klein_bottle(
    m: usize,
    n: usize,
) -> Vec<Simplex> {
    let id = |i: usize, j: usize| -> usize {
        let (i, j) = if j >= n { ((m - i % m) % m, j - n) } else { (i % m, j) };
        i * n + j
    };
    let mut triangles = Vec::new();
    for i in 0..m {
        for j in 0..n {
            let (v0, v1, v2, v3) = (id(i, j), id(i, j + 1), id(i + 1, j), id(i + 1, j + 1));
            triangles.push(vec![v0, v1, v2]);
            triangles.push(vec![v1, v3, v2]);
        }
    }
    close_complex(&triangles)
}

/// A disc whose boundary has `n` edges, the boundary vertex at position `p`
/// being identified with the vertices of the same `class(p)`; an inner ring
/// and a cone point keep the result a simplicial complex.
fn glued_disc(
    n: usize,
    class: impl Fn(usize) -> usize,
) -> Vec<Simplex> {
    let inner = |p: usize| n + (p % n);
    let center = 2 * n;
    let mut triangles = Vec::new();
    for p in 0..n {
        let (b0, b1) = (class(p % n), class((p + 1) % n));
        let (r0, r1) = (inner(p), inner(p + 1));
        triangles.push(vec![b0, b1, r0]);
        triangles.push(vec![b1, r0, r1]);
        triangles.push(vec![center, r0, r1]);
    }
    close_complex(&triangles)
}

/// A triangulation of the dunce cap: a disc whose boundary, divided into
/// three arcs of `k >= 3` edges, is glued by the word `a a a^-1`.
#[must_use]
pub fn dunce_cap(k: usize) -> Vec<Simplex> {
    glued_disc(3 * k, |p| {
        if p < k {
            p
        } else if p < 2 * k {
            p - k
        } else {
            (3 * k - p) % k
        }
    })
}

/// A triangulation of the real projective plane: for `k = 0` the minimal
/// one with 6 vertices; for `k >= 3` a disc with `2k` boundary edges whose
/// antipodal boundary points are identified (word `a a`).
#[must_use]
pub fn projective_plane(k: usize) -> Vec<Simplex> {
    if k == 0 {
        return close_complex(&[
            vec![0, 1, 2],
            vec![0, 2, 3],
            vec![0, 3, 4],
            vec![0, 4, 5],
            vec![0, 1, 5],
            vec![1, 2, 4],
            vec![2, 3, 5],
            vec![1, 3, 4],
            vec![2, 4, 5],
            vec![1, 3, 5],
        ]);
    }
    glued_disc(2 * k, |p| p % k)
}

// ----------------------------------------------------------------------
// Fundamental group
// ----------------------------------------------------------------------

/// A group presentation: generators are numbered `1..=count` and a word is
/// a list of nonzero integers (negative for inverses).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Presentation {
    /// Number of generators.
    pub count: usize,
    /// The relators.
    pub relators: Vec<Vec<i64>>,
    /// For a presentation read from a complex: the edge of each generator.
    pub edges: Vec<(usize, usize)>,
}

/// The edge-path presentation of the fundamental group of the component of
/// the smallest vertex: the generators are the edges outside a spanning
/// tree, the relators the boundaries of the triangles.
#[must_use]
pub fn pi1_presentation(complex: &[Simplex]) -> Option<Presentation> {
    let vertices: Vec<usize> = simplices_of_dim(complex, 0).into_iter().map(|s| s[0]).collect();
    let base = *vertices.first()?;
    let edges = simplices_of_dim(complex, 1);
    let mut adj: HashMap<usize, Vec<usize>> = HashMap::new();
    for e in &edges {
        adj.entry(e[0]).or_default().push(e[1]);
        adj.entry(e[1]).or_default().push(e[0]);
    }
    let mut seen: HashSet<usize> = HashSet::from([base]);
    let mut tree: HashSet<(usize, usize)> = HashSet::new();
    let mut queue = std::collections::VecDeque::from([base]);
    while let Some(u) = queue.pop_front() {
        for &v in adj.get(&u).map_or(&[][..], Vec::as_slice) {
            if seen.insert(v) {
                tree.insert((u.min(v), u.max(v)));
                queue.push_back(v);
            }
        }
    }
    let mut generator: HashMap<(usize, usize), i64> = HashMap::new();
    let mut gens = Vec::new();
    for e in &edges {
        if seen.contains(&e[0]) && !tree.contains(&(e[0], e[1])) {
            gens.push((e[0], e[1]));
            generator.insert((e[0], e[1]), i64::try_from(gens.len()).ok()?);
        }
    }
    let letter = |u: usize, v: usize| -> Option<i64> {
        if u < v { generator.get(&(u, v)).copied() } else { generator.get(&(v, u)).map(|g| -g) }
    };
    let mut relators = Vec::new();
    for t in simplices_of_dim(complex, 2) {
        if !seen.contains(&t[0]) {
            continue;
        }
        let word: Vec<i64> = [letter(t[0], t[1]), letter(t[1], t[2]), letter(t[2], t[0])].into_iter().flatten().collect();
        let word = cyclic_reduce(word);
        if !word.is_empty() {
            relators.push(word);
        }
    }
    Some(Presentation { count: gens.len(), relators, edges: gens })
}

fn free_reduce(word: Vec<i64>) -> Vec<i64> {
    let mut out: Vec<i64> = Vec::with_capacity(word.len());
    for x in word {
        if out.last() == Some(&-x) {
            out.pop();
        } else {
            out.push(x);
        }
    }
    out
}

fn cyclic_reduce(word: Vec<i64>) -> Vec<i64> {
    let mut w = free_reduce(word);
    while w.len() >= 2 && w.first() == w.last().map(|x| -x).as_ref() {
        w.remove(0);
        w.pop();
    }
    w
}

fn invert_word(word: &[i64]) -> Vec<i64> {
    word.iter().rev().map(|x| -x).collect()
}

/// Tietze simplification: removes trivial relators and eliminates
/// generators that occur once in some relator or that a one-letter
/// relator kills. The result is presented with generators `1..=count`.
#[must_use]
pub fn simplify_presentation(p: &Presentation) -> Presentation {
    const LIMIT: usize = 4000;
    let mut alive: Vec<bool> = vec![true; p.count + 1];
    let mut rels: Vec<Vec<i64>> = p.relators.iter().cloned().map(cyclic_reduce).filter(|r| !r.is_empty()).collect();
    loop {
        rels.sort();
        rels.dedup();
        // a one-letter relator kills its generator
        if let Some(pos) = rels.iter().position(|r| r.len() == 1) {
            let g = rels[pos][0].abs();
            rels.remove(pos);
            alive[usize::try_from(g).unwrap_or(0)] = false;
            rels = rels.into_iter().map(|r| cyclic_reduce(r.into_iter().filter(|x| x.abs() != g).collect())).filter(|r| !r.is_empty()).collect();
            continue;
        }
        // a generator occurring once in a relator is expressed by the rest
        let mut choice: Option<(usize, usize, i64)> = None;
        for (ri, r) in rels.iter().enumerate() {
            for (pos, &x) in r.iter().enumerate() {
                let occurrences = r.iter().filter(|y| y.abs() == x.abs()).count();
                if occurrences == 1 && choice.is_none_or(|(c, _, _)| rels[c].len() > r.len()) {
                    choice = Some((ri, pos, x));
                }
            }
        }
        let Some((ri, pos, x)) = choice else { break };
        let r = &rels[ri];
        // r = u x v  =>  x = (u^-1 v^-1) for x > 0, x^-1 = u^-1 v^-1 otherwise
        let u = &r[..pos];
        let v = &r[pos + 1..];
        let mut rest: Vec<i64> = invert_word(u);
        rest.extend(invert_word(v));
        // `rest` is x (if x > 0) or x^-1.
        let image = if x > 0 { rest } else { invert_word(&rest) };
        let g = x.abs();
        let substituted: Vec<Vec<i64>> = rels
            .iter()
            .enumerate()
            .filter(|&(i, _)| i != ri)
            .map(|(_, w)| {
                let mut out = Vec::new();
                for &y in w {
                    if y == g {
                        out.extend(&image);
                    } else if y == -g {
                        out.extend(invert_word(&image));
                    } else {
                        out.push(y);
                    }
                }
                cyclic_reduce(out)
            })
            .filter(|w| !w.is_empty())
            .collect();
        if substituted.iter().map(Vec::len).sum::<usize>() > LIMIT {
            break;
        }
        rels = substituted;
        alive[usize::try_from(g).unwrap_or(0)] = false;
    }
    let mut renumber = vec![0_i64; p.count + 1];
    let mut count = 0;
    for g in 1..=p.count {
        if alive[g] {
            count += 1;
            renumber[g] = i64::try_from(count).unwrap_or(0);
        }
    }
    let relators = rels
        .into_iter()
        .map(|r| r.into_iter().map(|x| x.signum() * renumber[usize::try_from(x.abs()).unwrap_or(0)]).collect())
        .collect();
    Presentation { count, relators, edges: Vec::new() }
}

/// The abelianisation of a presented group: `Z^rank + torsion`.
#[must_use]
pub fn abelianization(p: &Presentation) -> Group {
    let rows: Vec<Vec<i64>> = p
        .relators
        .iter()
        .map(|r| {
            let mut row = vec![0_i64; p.count];
            for &x in r {
                if let Some(c) = row.get_mut(usize::try_from(x.abs()).unwrap_or(0).wrapping_sub(1)) {
                    *c += x.signum();
                }
            }
            row
        })
        .collect();
    let inv = if p.count == 0 || rows.is_empty() { Vec::new() } else { smith_invariants(&rows) };
    Group { rank: p.count - inv.len(), torsion: inv.into_iter().filter(|t| !t.is_one()).collect() }
}

// ----------------------------------------------------------------------
// Simplicial maps
// ----------------------------------------------------------------------

/// Whether the vertex map sends every simplex of `k` onto a simplex of `l`.
#[must_use]
pub fn is_simplicial_map<S: std::hash::BuildHasher>(
    k: &[Simplex],
    l: &[Simplex],
    f: &HashMap<usize, usize, S>,
) -> bool {
    let in_l: HashSet<&Simplex> = l.iter().collect();
    k.iter().all(|s| {
        let Some(mut image) = s.iter().map(|v| f.get(v).copied()).collect::<Option<Simplex>>() else { return false };
        image.sort_unstable();
        image.dedup();
        in_l.contains(&image)
    })
}

/// The matrix (over `L`'s `k`-simplices, columns `K`'s) of the chain map
/// induced by a simplicial map.
#[must_use]
pub fn induced_chain_map<S: std::hash::BuildHasher>(
    k: &[Simplex],
    l: &[Simplex],
    f: &HashMap<usize, usize, S>,
    dim: usize,
) -> Option<Vec<Vec<i64>>> {
    let (src, dst) = (simplices_of_dim(k, dim), simplices_of_dim(l, dim));
    let index: HashMap<&Simplex, usize> = dst.iter().enumerate().map(|(i, s)| (s, i)).collect();
    let mut m = vec![vec![0_i64; src.len()]; dst.len()];
    for (j, s) in src.iter().enumerate() {
        let image: Vec<usize> = s.iter().map(|v| f.get(v).copied()).collect::<Option<_>>()?;
        let mut sorted = image.clone();
        sorted.sort_unstable();
        if sorted.windows(2).any(|w| w[0] == w[1]) {
            continue;
        }
        let inversions = (0..image.len()).flat_map(|a| (a + 1..image.len()).map(move |b| (a, b))).filter(|&(a, b)| image[a] > image[b]).count();
        let i = *index.get(&sorted)?;
        m[i][j] = if inversions % 2 == 0 { 1 } else { -1 };
    }
    Some(m)
}

fn to_q(m: &[Vec<i64>]) -> ql::QMat {
    m.iter().map(|r| r.iter().map(|&x| ql::q(x)).collect()).collect()
}

/// Cycle representatives forming, together with a basis of the boundaries,
/// a basis of the cycles: `(boundary basis, homology representatives)` as
/// vectors of `C_k` over `Q`.
fn homology_basis(
    cc: &ChainComplex,
    k: usize,
) -> (ql::QMat, ql::QMat) {
    let n = cc.dim_at(k);
    let cycles = if k == 0 || cc.dim_at(k - 1) == 0 { ql::identity(n) } else { ql::nullspace(&to_q(&cc.d[k]), n) };
    let boundaries: ql::QMat = if cc.dim_at(k + 1) == 0 || n == 0 {
        Vec::new()
    } else {
        let d = to_q(&cc.d[k + 1]);
        ql::span_basis(&ql::transpose(&d, cc.dim_at(k + 1)), n)
    };
    let mut echelon = boundaries.clone();
    let mut pivots = ql::rref(&mut echelon, n);
    let mut reps = Vec::new();
    for z in cycles {
        let r = ql::reduce_mod(&z, &echelon, &pivots);
        if !ql::is_zero_vec(&r) {
            echelon.push(r.clone());
            pivots = ql::rref(&mut echelon, n);
            reps.push(z);
        }
    }
    (boundaries, reps)
}

/// The matrix of `f_*: H_k(K; Q) -> H_k(L; Q)` in the homology bases
/// chosen by this module (rows: `L`, columns: `K`).
#[must_use]
pub fn induced_homology<S: std::hash::BuildHasher>(
    k: &[Simplex],
    l: &[Simplex],
    f: &HashMap<usize, usize, S>,
    dim: usize,
) -> Option<ql::QMat> {
    if !is_simplicial_map(k, l, f) {
        return None;
    }
    let (ck, cl) = (complex_chains(k), complex_chains(l));
    let map = to_q(&induced_chain_map(k, l, f, dim)?);
    let (_, reps_k) = homology_basis(&ck, dim);
    let (bounds_l, reps_l) = homology_basis(&cl, dim);
    let n = cl.dim_at(dim);
    // columns: boundaries then representatives of L
    let cols: ql::QMat = bounds_l.iter().chain(&reps_l).cloned().collect();
    let system = ql::transpose(&cols, n);
    let nb = bounds_l.len();
    let mut out = vec![Vec::new(); reps_l.len()];
    for z in &reps_k {
        let image = ql::matvec(&map, z);
        let x = if cols.is_empty() { Vec::new() } else { ql::solve(&system, cols.len(), &image)? };
        for (row, value) in out.iter_mut().zip(x.into_iter().skip(nb)) {
            row.push(value);
        }
    }
    Some(out)
}

/// The rational matrix as `(rows, columns)` of entries.
#[must_use]
pub fn matrix_rank(m: &[Vec<Q>]) -> usize {
    ql::rank(m, m.first().map_or(0, Vec::len))
}

// ----------------------------------------------------------------------
// Persistent homology
// ----------------------------------------------------------------------

/// A persistence bar: the homological dimension, the birth and the death
/// (`f64::INFINITY` for a class that never dies).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bar {
    /// Homological dimension.
    pub dim: usize,
    /// Birth value.
    pub birth: f64,
    /// Death value.
    pub death: f64,
}

/// Persistent homology over `Z/2` of a filtered complex.
///
/// The input is
/// `(value, simplex)` pairs (any order; every face of a simplex must be
/// present with a value no larger). Bars of dimension up to `max_dim`
/// with positive length are returned, sorted by dimension, birth, death.
///
/// This is the standard algorithm: the columns of the boundary matrix,
/// ordered by (value, dimension, lexicographic), are reduced from left to
/// right until their lowest ones are distinct.
#[must_use]
pub fn persistence_bars(
    filtration: &[(f64, Simplex)],
    max_dim: usize,
) -> Option<Vec<Bar>> {
    let mut items: Vec<(f64, Simplex)> = filtration
        .iter()
        .map(|(v, s)| {
            let mut s = s.clone();
            s.sort_unstable();
            s.dedup();
            (*v, s)
        })
        .filter(|(_, s)| !s.is_empty())
        .collect();
    items.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.len().cmp(&b.1.len())).then_with(|| a.1.cmp(&b.1)));
    items.dedup_by(|a, b| a.1 == b.1);
    let index: HashMap<Simplex, usize> = items.iter().enumerate().map(|(i, (_, s))| (s.clone(), i)).collect();
    let m = items.len();
    let mut columns: Vec<Vec<usize>> = Vec::with_capacity(m);
    for (j, (_, s)) in items.iter().enumerate() {
        let mut col = Vec::new();
        if s.len() > 1 {
            for skip in 0..s.len() {
                let face: Simplex = s.iter().enumerate().filter(|&(i, _)| i != skip).map(|(_, &v)| v).collect();
                let i = *index.get(&face)?;
                if i >= j {
                    return None;
                }
                col.push(i);
            }
            col.sort_unstable();
        }
        columns.push(col);
    }
    let mut owner: HashMap<usize, usize> = HashMap::new();
    let mut paired = vec![false; m];
    let mut bars = Vec::new();
    for j in 0..m {
        while let Some(&low) = columns[j].last() {
            let Some(&other) = owner.get(&low) else {
                owner.insert(low, j);
                paired[low] = true;
                paired[j] = true;
                let (birth, death) = (items[low].0, items[j].0);
                if death > birth && items[low].1.len() <= max_dim + 1 {
                    bars.push(Bar { dim: items[low].1.len() - 1, birth, death });
                }
                break;
            };
            columns[j] = symmetric_difference(&columns[j], &columns[other]);
        }
    }
    for (i, (value, s)) in items.iter().enumerate() {
        if !paired[i] && s.len() <= max_dim + 1 {
            bars.push(Bar { dim: s.len() - 1, birth: *value, death: f64::INFINITY });
        }
    }
    bars.sort_by(|a, b| a.dim.cmp(&b.dim).then(a.birth.total_cmp(&b.birth)).then(a.death.total_cmp(&b.death)));
    Some(bars)
}

fn symmetric_difference(
    a: &[usize],
    b: &[usize],
) -> Vec<usize> {
    let mut out = Vec::with_capacity(a.len() + b.len());
    let (mut i, mut j) = (0, 0);
    while i < a.len() || j < b.len() {
        match (a.get(i), b.get(j)) {
            | (Some(&x), Some(&y)) if x == y => {
                i += 1;
                j += 1;
            },
            | (Some(&x), Some(&y)) if x < y => {
                out.push(x);
                i += 1;
            },
            | (Some(_) | None, Some(&y)) => {
                out.push(y);
                j += 1;
            },
            | (Some(&x), None) => {
                out.push(x);
                i += 1;
            },
            | (None, None) => break,
        }
    }
    out
}

/// The Vietoris-Rips filtration of a finite metric space.
///
/// The space is given by its
/// distance matrix: every clique of the graph of edges of length at most
/// `max_eps`, up to dimension `max_simplex_dim`, with value its diameter.
/// `None` when more than `limit` simplices would be needed.
#[must_use]
pub fn rips_filtration(
    dist: &[Vec<f64>],
    max_eps: f64,
    max_simplex_dim: usize,
    limit: usize,
) -> Option<Vec<(f64, Simplex)>> {
    let n = dist.len();
    let mut out: Vec<(f64, Simplex)> = (0..n).map(|v| (0.0, vec![v])).collect();
    let mut stack: Vec<(Simplex, f64, Vec<usize>)> = (0..n)
        .map(|v| (vec![v], 0.0, (v + 1..n).filter(|&w| dist[v][w] <= max_eps).collect()))
        .collect();
    while let Some((simplex, diameter, candidates)) = stack.pop() {
        if simplex.len() > max_simplex_dim {
            continue;
        }
        for (pos, &w) in candidates.iter().enumerate() {
            let d = simplex.iter().fold(diameter, |acc, &v| acc.max(dist[v][w]));
            let mut next = simplex.clone();
            next.push(w);
            out.push((d, next.clone()));
            if out.len() > limit {
                return None;
            }
            let rest: Vec<usize> = candidates[pos + 1..].iter().copied().filter(|&x| dist[w][x] <= max_eps).collect();
            stack.push((next, d, rest));
        }
    }
    Some(out)
}

/// The distance matrix of a point cloud.
#[must_use]
pub fn distance_matrix(points: &[Vec<f64>]) -> Vec<Vec<f64>> {
    points.iter().map(|p| points.iter().map(|q| super::topology::euclidean_distance(p, q)).collect()).collect()
}

/// The bottleneck distance of two persistence diagrams (lists of
/// `(birth, death)`). Infinite bars are matched among themselves; if their
/// numbers differ the distance is infinite.
#[must_use]
pub fn bottleneck_distance(
    a: &[(f64, f64)],
    b: &[(f64, f64)],
) -> f64 {
    let split = |d: &[(f64, f64)]| -> (Vec<(f64, f64)>, Vec<f64>) {
        (d.iter().copied().filter(|p| p.1.is_finite()).collect(), d.iter().filter(|p| !p.1.is_finite()).map(|p| p.0).collect())
    };
    let ((fa, mut ia), (fb, mut ib)) = (split(a), split(b));
    if ia.len() != ib.len() {
        return f64::INFINITY;
    }
    ia.sort_by(f64::total_cmp);
    ib.sort_by(f64::total_cmp);
    let infinite = ia.iter().zip(&ib).map(|(x, y)| (x - y).abs()).fold(0.0, f64::max);
    infinite.max(finite_bottleneck(&fa, &fb))
}

fn finite_bottleneck(
    a: &[(f64, f64)],
    b: &[(f64, f64)],
) -> f64 {
    let (na, nb) = (a.len(), b.len());
    if na + nb == 0 {
        return 0.0;
    }
    let linf = |p: (f64, f64), q: (f64, f64)| (p.0 - q.0).abs().max((p.1 - q.1).abs());
    let diag = |p: (f64, f64)| (p.1 - p.0).abs() / 2.0;
    let mut candidates = vec![0.0];
    for &p in a {
        candidates.push(diag(p));
        candidates.extend(b.iter().map(|&q| linf(p, q)));
    }
    candidates.extend(b.iter().map(|&q| diag(q)));
    candidates.sort_by(f64::total_cmp);
    candidates.dedup();
    let feasible = |delta: f64| -> bool {
        // left: a_0..a_{na-1}, then diagonal copies of b; right: b, then diagonal copies of a
        let n = na + nb;
        let mut adj = vec![Vec::new(); n];
        for i in 0..na {
            for (j, &bj) in b.iter().enumerate() {
                if linf(a[i], bj) <= delta {
                    adj[i].push(j);
                }
            }
            if diag(a[i]) <= delta {
                adj[i].push(nb + i);
            }
        }
        for j in 0..nb {
            if diag(b[j]) <= delta {
                adj[na + j].push(j);
            }
            for i in 0..na {
                adj[na + j].push(nb + i);
            }
        }
        let mut matched: Vec<Option<usize>> = vec![None; n];
        (0..n).all(|u| {
            let mut visited = vec![false; n];
            augment(u, &adj, &mut visited, &mut matched)
        })
    };
    let (mut lo, mut hi) = (0, candidates.len() - 1);
    while lo < hi {
        let mid = usize::midpoint(lo, hi);
        if feasible(candidates[mid]) {
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }
    candidates[lo]
}

fn augment(
    u: usize,
    adj: &[Vec<usize>],
    visited: &mut [bool],
    matched: &mut [Option<usize>],
) -> bool {
    for &v in &adj[u] {
        if visited[v] {
            continue;
        }
        visited[v] = true;
        if matched[v].is_none_or(|w| augment(w, adj, visited, matched)) {
            matched[v] = Some(u);
            return true;
        }
    }
    false
}

// ----------------------------------------------------------------------
// Cubical complexes of binary images
// ----------------------------------------------------------------------

/// Cells `(vertices, horizontal edges, vertical edges, squares)` of the
/// cubical complex of a binary image, each cell named by its lower-left
/// lattice point. With `eight`, the foreground is the union of the closed
/// pixel squares (8-connected foreground, 4-connected background);
/// otherwise vertices are pixels, edges join orthogonal neighbours and
/// squares are full 2x2 blocks (4-connected foreground).
type Cells = (BTreeSet<(i64, i64)>, BTreeSet<(i64, i64)>, BTreeSet<(i64, i64)>, BTreeSet<(i64, i64)>);

fn cubical_cells(
    image: &[Vec<bool>],
    eight: bool,
) -> Cells {
    let on = |x: i64, y: i64| -> bool {
        usize::try_from(y).ok().and_then(|y| image.get(y)).and_then(|row| usize::try_from(x).ok().and_then(|x| row.get(x))).copied().unwrap_or(false)
    };
    let (mut vs, mut hs, mut ws, mut sq) = (BTreeSet::new(), BTreeSet::new(), BTreeSet::new(), BTreeSet::new());
    for (y, row) in image.iter().enumerate() {
        for x in 0..row.len() {
            let (xi, yi) = (i64::try_from(x).unwrap_or(0), i64::try_from(y).unwrap_or(0));
            if !on(xi, yi) {
                continue;
            }
            if eight {
                sq.insert((xi, yi));
                vs.extend([(xi, yi), (xi + 1, yi), (xi, yi + 1), (xi + 1, yi + 1)]);
                hs.extend([(xi, yi), (xi, yi + 1)]);
                ws.extend([(xi, yi), (xi + 1, yi)]);
            } else {
                vs.insert((xi, yi));
                if on(xi + 1, yi) {
                    hs.insert((xi, yi));
                }
                if on(xi, yi + 1) {
                    ws.insert((xi, yi));
                }
                if on(xi + 1, yi) && on(xi, yi + 1) && on(xi + 1, yi + 1) {
                    sq.insert((xi, yi));
                }
            }
        }
    }
    (vs, hs, ws, sq)
}

const fn find(
    parent: &mut [usize],
    mut x: usize,
) -> usize {
    while parent[x] != x {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    x
}

/// The Betti numbers `(b_0, b_1)` of a binary image, from the Euler
/// characteristic of its cubical complex and the number of components (a
/// planar complex has no `H_2`).
#[must_use]
pub fn cubical_betti(
    image: &[Vec<bool>],
    eight: bool,
) -> (usize, usize) {
    let (vs, hs, ws, sq) = cubical_cells(image, eight);
    let vertices: Vec<(i64, i64)> = vs.iter().copied().collect();
    let position: HashMap<(i64, i64), usize> = vertices.iter().enumerate().map(|(i, &v)| (v, i)).collect();
    let mut parent: Vec<usize> = (0..vertices.len()).collect();
    let mut unite = |a: (i64, i64), b: (i64, i64)| {
        if let (Some(&i), Some(&j)) = (position.get(&a), position.get(&b)) {
            let (ri, rj) = (find(&mut parent, i), find(&mut parent, j));
            parent[ri] = rj;
        }
    };
    for &(x, y) in &hs {
        unite((x, y), (x + 1, y));
    }
    for &(x, y) in &ws {
        unite((x, y), (x, y + 1));
    }
    let components = (0..vertices.len()).filter(|&i| find(&mut parent, i) == i).count();
    let chi = i64::try_from(vs.len() + sq.len()).unwrap_or(0) - i64::try_from(hs.len() + ws.len()).unwrap_or(0);
    let b1 = i64::try_from(components).unwrap_or(0) - chi;
    (components, usize::try_from(b1).unwrap_or(0))
}

/// The chain complex of the cubical complex of a binary image (for
/// cross-checking [`cubical_betti`] with the Smith normal form).
#[must_use]
pub fn cubical_chains(
    image: &[Vec<bool>],
    eight: bool,
) -> ChainComplex {
    let (vs, hs, ws, sq) = cubical_cells(image, eight);
    let vpos: HashMap<(i64, i64), usize> = vs.iter().enumerate().map(|(i, &v)| (v, i)).collect();
    let edges: Vec<((i64, i64), bool)> = hs.iter().map(|&p| (p, true)).chain(ws.iter().map(|&p| (p, false))).collect();
    let epos: HashMap<((i64, i64), bool), usize> = edges.iter().enumerate().map(|(i, &e)| (e, i)).collect();
    let mut d1 = vec![vec![0_i64; edges.len()]; vs.len()];
    for (j, &((x, y), horizontal)) in edges.iter().enumerate() {
        let end = if horizontal { (x + 1, y) } else { (x, y + 1) };
        d1[vpos[&end]][j] += 1;
        d1[vpos[&(x, y)]][j] -= 1;
    }
    let squares: Vec<(i64, i64)> = sq.iter().copied().collect();
    let mut d2 = vec![vec![0_i64; squares.len()]; edges.len()];
    for (j, &(x, y)) in squares.iter().enumerate() {
        d2[epos[&((x, y), true)]][j] += 1;
        d2[epos[&((x + 1, y), false)]][j] += 1;
        d2[epos[&((x, y + 1), true)]][j] -= 1;
        d2[epos[&((x, y), false)]][j] -= 1;
    }
    ChainComplex { dims: vec![vs.len(), edges.len(), squares.len()], d: vec![Vec::new(), d1, d2] }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tors(g: &Group) -> Vec<i64> {
        use num_traits::ToPrimitive;
        g.torsion.iter().filter_map(ToPrimitive::to_i64).collect()
    }

    #[test]
    fn smith_form() {
        let m = vec![vec![2, 4, 4], vec![-6, 6, 12], vec![10, -4, -16]];
        let inv: Vec<i64> = smith_invariants(&m).iter().map(|x| i64::try_from(x.clone()).unwrap_or(0)).collect();
        assert_eq!(inv, vec![2, 6, 12]);
        assert_eq!(smith_invariants(&[vec![0, 0], vec![0, 0]]).len(), 0);
        assert_eq!(rank_mod_p(&[vec![2, 0], vec![0, 1]], 2), 1);
        assert_eq!(rank_mod_p(&[vec![2, 0], vec![0, 1]], 3), 2);
    }

    #[test]
    fn integer_homology_of_standard_spaces() {
        let rp2 = complex_chains(&projective_plane(0));
        assert_eq!(rp2.dims, vec![6, 15, 10]);
        assert_eq!((rp2.homology(0), rp2.homology(2).rank), (Group { rank: 1, torsion: vec![] }, 0));
        assert_eq!((rp2.homology(1).rank, tors(&rp2.homology(1))), (0, vec![2]));
        assert_eq!(tors(&rp2.cohomology(2)), vec![2]);
        assert_eq!(rp2.cohomology(1).torsion.len(), 0);
        // over Z/2 the projective plane has b_1 = b_2 = 1, over Z/3 it looks like a point
        assert_eq!((rp2.homology_mod_p(1, 2), rp2.homology_mod_p(2, 2)), (1, 1));
        assert_eq!((rp2.homology_mod_p(1, 3), rp2.homology_mod_p(2, 3)), (0, 0));
        let torus = complex_chains(&super::super::topology::torus_complex(3, 3));
        assert_eq!((torus.homology(1).rank, torus.homology(2).rank), (2, 1));
        assert!(torus.homology(1).torsion.is_empty() && torus.homology(2).torsion.is_empty());
        let klein = complex_chains(&klein_bottle(4, 4));
        assert_eq!((klein.homology(0).rank, klein.homology(1).rank, klein.homology(2).rank), (1, 1, 0));
        assert_eq!(tors(&klein.homology(1)), vec![2]);
        assert_eq!(super::super::topology::euler_characteristic(&klein_bottle(4, 4)), 0);
        let dunce = complex_chains(&dunce_cap(3));
        assert_eq!((dunce.homology(0).rank, dunce.homology(1), dunce.homology(2)), (1, Group { rank: 0, torsion: vec![] }, Group { rank: 0, torsion: vec![] }));
    }

    #[test]
    fn spheres_cones_joins_products() {
        for n in 0..4 {
            let cc = complex_chains(&sphere(n));
            assert_eq!(cc.homology(n).rank, if n == 0 { 2 } else { 1 });
            assert_eq!(reduced_homology(&cc, 0).rank, usize::from(n == 0));
        }
        let s1 = sphere(1);
        let c = complex_chains(&cone(&s1));
        assert_eq!((c.homology(0).rank, c.homology(1).rank, c.homology(2).rank), (1, 0, 0));
        assert_eq!(betti_list(&suspension(&s1)), vec![1, 0, 1]);
        assert_eq!(betti_list(&join(&sphere(1), &sphere(1))), vec![1, 0, 0, 1]);
        assert_eq!(betti_list(&wedge(&s1, &s1, 0, 0)), vec![1, 2]);
        let t = product(&s1, &s1, 100_000).expect("small");
        assert_eq!(betti_list(&t), vec![1, 2, 1]);
        let h = complex_chains(&t);
        assert!(h.homology(1).torsion.is_empty());
        assert_eq!(betti_list(&product(&full_simplex(2), &full_simplex(1), 100_000).expect("small")), vec![1, 0, 0, 0]);
        // Z/2 homology of RP^2 x S^1 has the right Euler characteristic
        let p = product(&projective_plane(0), &s1, 1_000_000).expect("fits");
        assert_eq!(super::super::topology::euler_characteristic(&p), 0);
    }

    #[test]
    fn relative_homology_and_mayer_vietoris() {
        let disc = full_simplex(2);
        let boundary = sphere(1);
        let rel = relative_chains(&disc, &boundary);
        assert_eq!((rel.homology(2).rank, rel.homology(1).rank, rel.homology(0).rank), (1, 0, 0));
        assert_eq!(reduced_homology(&complex_chains(&disc), 0).rank, 0);
        // S^2 as two discs glued along the equator
        let a = close_complex(&[vec![0, 1, 3], vec![1, 2, 3], vec![0, 2, 3]]);
        let b = close_complex(&[vec![0, 1, 4], vec![1, 2, 4], vec![0, 2, 4]]);
        assert_eq!(betti_list(&union(&a, &b)), vec![1, 0, 1]);
        assert_eq!(mayer_vietoris(&a, &b), Some(vec![0, 1, 0, 0]));
        // torus as the union of two annuli is not available directly; use two arcs of a circle
        let left = close_complex(&[vec![0, 1], vec![1, 2]]);
        let right = close_complex(&[vec![2, 3], vec![3, 0]]);
        assert_eq!(betti_list(&intersection(&left, &right)), vec![2]);
        assert_eq!(mayer_vietoris(&left, &right), Some(vec![1, 0, 0]));
        assert!(is_subcomplex(&boundary, &disc));
        assert!(!is_subcomplex(&disc, &boundary));
    }

    #[test]
    fn fundamental_groups() {
        let torus = super::super::topology::torus_complex(3, 3);
        let raw = pi1_presentation(&torus).expect("nonempty");
        assert_eq!(raw.count, 27 - 9 + 1);
        assert!(raw.relators.len() <= 18);
        let small = simplify_presentation(&raw);
        let ab = abelianization(&small);
        assert_eq!((ab.rank, ab.torsion.len()), (2, 0));
        assert_eq!(abelianization(&raw), ab);
        assert_eq!(small.count, 2);
        let rp2 = simplify_presentation(&pi1_presentation(&projective_plane(0)).expect("nonempty"));
        assert_eq!((rp2.count, rp2.relators.clone()), (1, vec![vec![1, 1]]));
        assert_eq!(tors(&abelianization(&rp2)), vec![2]);
        let klein = simplify_presentation(&pi1_presentation(&klein_bottle(4, 4)).expect("nonempty"));
        let ab = abelianization(&klein);
        assert_eq!((ab.rank, tors(&ab)), (1, vec![2]));
        let dunce = simplify_presentation(&pi1_presentation(&dunce_cap(3)).expect("nonempty"));
        assert_eq!(dunce.count, 0);
    }

    #[test]
    fn simplicial_maps_in_homology() {
        let circle = sphere(1);
        let disc = full_simplex(2);
        let id: HashMap<usize, usize> = (0..3).map(|v| (v, v)).collect();
        assert!(is_simplicial_map(&circle, &disc, &id));
        let m = induced_homology(&circle, &disc, &id, 1).expect("simplicial");
        assert!(m.is_empty());
        let m0 = induced_homology(&circle, &disc, &id, 0).expect("simplicial");
        assert_eq!(m0, vec![vec![ql::q(1)]]);
        let swap: HashMap<usize, usize> = [(0, 1), (1, 0), (2, 2)].into_iter().collect();
        let m1 = induced_homology(&circle, &circle, &swap, 1).expect("simplicial");
        assert_eq!(m1, vec![vec![ql::q(-1)]]);
        let bad: HashMap<usize, usize> = [(0, 0), (1, 1), (2, 5)].into_iter().collect();
        assert!(!is_simplicial_map(&circle, &circle, &bad));
        assert_eq!(matrix_rank(&m1), 1);
    }

    #[test]
    fn persistent_homology_of_a_circle() {
        let points: Vec<Vec<f64>> = (0..12).map(|i| {
            let t = std::f64::consts::TAU * f64::from(i) / 12.0;
            vec![t.cos(), t.sin()]
        }).collect();
        let f = rips_filtration(&distance_matrix(&points), 3.0, 2, 1_000_000).expect("small");
        let bars = persistence_bars(&f, 1).expect("valid");
        let h0: Vec<&Bar> = bars.iter().filter(|b| b.dim == 0).collect();
        assert_eq!(h0.len(), 12);
        assert_eq!(h0.iter().filter(|b| b.death.is_infinite()).count(), 1);
        let h1: Vec<&Bar> = bars.iter().filter(|b| b.dim == 1).collect();
        let long: Vec<&&Bar> = h1.iter().filter(|b| b.death - b.birth > 0.8).collect();
        assert_eq!(long.len(), 1, "{h1:?}");
        assert!((long[0].birth - 0.5176).abs() < 1e-3);
    }

    #[test]
    fn bottleneck() {
        let a = [(0.0, 1.0), (0.0, 3.0)];
        let b = [(0.0, 1.5), (0.0, 3.0)];
        assert!((bottleneck_distance(&a, &b) - 0.5).abs() < 1e-12);
        assert!((bottleneck_distance(&a, &[]) - 1.5).abs() < 1e-12);
        assert!(bottleneck_distance(&a, &a).abs() < 1e-12);
        assert!(bottleneck_distance(&[(0.0, f64::INFINITY)], &[]).is_infinite());
        assert!((bottleneck_distance(&[(0.0, f64::INFINITY), (1.0, 2.0)], &[(0.5, f64::INFINITY), (1.0, 2.0)]) - 0.5).abs() < 1e-12);
    }

    #[test]
    fn cubical() {
        let ring: Vec<Vec<bool>> = ["###", "#.#", "###"].iter().map(|r| r.chars().map(|c| c == '#').collect()).collect();
        assert_eq!(cubical_betti(&ring, true), (1, 1));
        assert_eq!(cubical_betti(&ring, false), (1, 1));
        let diag: Vec<Vec<bool>> = ["#.", ".#"].iter().map(|r| r.chars().map(|c| c == '#').collect()).collect();
        assert_eq!(cubical_betti(&diag, true), (1, 0));
        assert_eq!(cubical_betti(&diag, false), (2, 0));
        let two: Vec<Vec<bool>> = ["#.#", "..."].iter().map(|r| r.chars().map(|c| c == '#').collect()).collect();
        assert_eq!(cubical_betti(&two, true), (2, 0));
        // a checkerboard ring closed only diagonally
        let diag_ring: Vec<Vec<bool>> = [".#.", "#.#", ".#."].iter().map(|r| r.chars().map(|c| c == '#').collect()).collect();
        assert_eq!(cubical_betti(&diag_ring, true), (1, 1));
        assert_eq!(cubical_betti(&diag_ring, false), (4, 0));
        for (img, eight) in [(&ring, true), (&ring, false), (&diag_ring, true), (&two, true)] {
            let cc = cubical_chains(img, eight);
            let (b0, b1) = cubical_betti(img, eight);
            assert_eq!((cc.homology(0).rank, cc.homology(1).rank, cc.homology(2).rank), (b0, b1, 0));
        }
    }
}
