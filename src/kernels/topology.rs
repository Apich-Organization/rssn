//! # Numerical Computational Topology
//!
//! Pure numeric algorithms for computational topology: Euclidean distances,
//! Vietoris-Rips complexes, exact Betti numbers of simplicial complexes
//! (ranks of integer boundary matrices are computed exactly over the
//! rationals), exact persistence diagrams (boundary-matrix reduction over
//! `Z/2` with clearing) and connected components of a
//! graph given by adjacency lists.
//!
//! A simplex is a sorted vector of distinct vertex indices; a complex is a
//! list of simplices closed under taking faces, sorted by dimension and then
//! lexicographically (see [`close_complex`]).

use std::collections::BTreeSet;
use std::collections::HashMap;
use std::collections::VecDeque;

use num_rational::BigRational;
use num_traits::Zero;

/// A simplex: a sorted vector of distinct vertex indices.
pub type Simplex = Vec<usize>;

/// A persistence interval (birth, death).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PersistenceInterval {
    /// Birth radius.
    pub birth: f64,
    /// Death radius.
    pub death: f64,
}

/// A persistence diagram for one homological dimension.
#[derive(Debug, Clone, PartialEq)]
pub struct PersistenceDiagram {
    /// Homological dimension.
    pub dimension: usize,
    /// The intervals of the diagram.
    pub intervals: Vec<PersistenceInterval>,
}

/// The connected components of a graph given by adjacency lists (BFS).
///
/// Neighbours outside `0..adj.len()` are ignored.
#[must_use]
pub fn find_connected_components(adj: &[Vec<usize>]) -> Vec<Vec<usize>> {
    let mut visited = vec![false; adj.len()];
    let mut components = Vec::new();
    for start in 0..adj.len() {
        if visited[start] {
            continue;
        }
        let mut component = Vec::new();
        let mut queue = VecDeque::new();
        visited[start] = true;
        queue.push_back(start);
        while let Some(u) = queue.pop_front() {
            component.push(u);
            for &v in &adj[u] {
                if v < adj.len() && !visited[v] {
                    visited[v] = true;
                    queue.push_back(v);
                }
            }
        }
        components.push(component);
    }
    components
}

/// The Euclidean distance between two points.
#[must_use]
pub fn euclidean_distance(
    p1: &[f64],
    p2: &[f64],
) -> f64 {
    p1.iter().zip(p2).map(|(a, b)| (a - b).powi(2)).sum::<f64>().sqrt()
}

/// The Vietoris-Rips complex of a point cloud.
///
/// It holds every set of points that are
/// pairwise within `epsilon`, up to dimension `max_dim`. Simplices are
/// returned by increasing dimension, lexicographically within a dimension.
#[must_use]
pub fn vietoris_rips_complex(
    points: &[Vec<f64>],
    epsilon: f64,
    max_dim: usize,
) -> Vec<Simplex> {
    let n = points.len();
    let mut simplices: Vec<Simplex> = (0..n).map(|i| vec![i]).collect();
    if max_dim == 0 {
        return simplices;
    }
    let mut current = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            if euclidean_distance(&points[i], &points[j]) <= epsilon {
                current.push(vec![i, j]);
            }
        }
    }
    simplices.extend(current.iter().cloned());
    for _ in 2..=max_dim {
        let mut next = Vec::new();
        for simplex in &current {
            let Some(&last) = simplex.last() else { continue };
            for i in (last + 1)..n {
                if simplex.iter().all(|&v| euclidean_distance(&points[v], &points[i]) <= epsilon) {
                    let mut bigger = simplex.clone();
                    bigger.push(i);
                    next.push(bigger);
                }
            }
        }
        if next.is_empty() {
            break;
        }
        simplices.extend(next.iter().cloned());
        current = next;
    }
    simplices
}

/// The simplices of `simplices` closed under faces, without duplicates,
/// sorted by dimension and then lexicographically. Each simplex is sorted
/// and deduplicated first; empty simplices are dropped.
#[must_use]
pub fn close_complex(simplices: &[Simplex]) -> Vec<Simplex> {
    let mut seen: BTreeSet<Simplex> = BTreeSet::new();
    let mut stack: Vec<Simplex> = simplices
        .iter()
        .map(|s| {
            let mut s = s.clone();
            s.sort_unstable();
            s.dedup();
            s
        })
        .filter(|s| !s.is_empty())
        .collect();
    while let Some(s) = stack.pop() {
        if !seen.insert(s.clone()) {
            continue;
        }
        if s.len() > 1 {
            for skip in 0..s.len() {
                let face: Simplex = s.iter().enumerate().filter(|&(i, _)| i != skip).map(|(_, &v)| v).collect();
                if !seen.contains(&face) {
                    stack.push(face);
                }
            }
        }
    }
    let mut out: Vec<Simplex> = seen.into_iter().collect();
    out.sort_by(|a, b| a.len().cmp(&b.len()).then_with(|| a.cmp(b)));
    out
}

/// The dimension of a closed complex (`None` if empty).
#[must_use]
pub fn complex_dimension(complex: &[Simplex]) -> Option<usize> {
    complex.iter().map(|s| s.len().saturating_sub(1)).max()
}

/// The `k`-simplices of a complex, in order.
#[must_use]
pub fn simplices_of_dim(
    complex: &[Simplex],
    k: usize,
) -> Vec<Simplex> {
    complex.iter().filter(|s| s.len() == k + 1).cloned().collect()
}

/// The Euler characteristic `sum (-1)^k #k-simplices`.
#[must_use]
pub fn euler_characteristic(complex: &[Simplex]) -> i64 {
    complex.iter().map(|s| if s.len() % 2 == 1 { 1 } else { -1 }).sum()
}

/// The boundary of an oriented simplex: `(sign, face)` pairs, with sign
/// `(-1)^i` for the face omitting vertex `i`. Empty for a vertex.
#[must_use]
pub fn simplex_boundary(simplex: &[usize]) -> Vec<(i64, Simplex)> {
    if simplex.len() < 2 {
        return Vec::new();
    }
    (0..simplex.len())
        .map(|i| {
            let face = simplex.iter().enumerate().filter(|&(j, _)| j != i).map(|(_, &v)| v).collect();
            (if i % 2 == 0 { 1 } else { -1 }, face)
        })
        .collect()
}

/// The boundary matrix of `cols` (k-simplices) onto `rows` ((k-1)-simplices);
/// `rows.len()` rows and `cols.len()` columns.
#[must_use]
pub fn boundary_matrix(
    rows: &[Simplex],
    cols: &[Simplex],
) -> Vec<Vec<i64>> {
    let index: HashMap<&Simplex, usize> = rows.iter().enumerate().map(|(i, s)| (s, i)).collect();
    let mut m = vec![vec![0; cols.len()]; rows.len()];
    for (j, s) in cols.iter().enumerate() {
        for (sign, face) in simplex_boundary(s) {
            if let Some(&i) = index.get(&face) {
                m[i][j] = sign;
            }
        }
    }
    m
}

/// The product of two integer matrices (`a` is `r x m`, `b` is `m x c`).
#[must_use]
pub fn int_matmul(
    a: &[Vec<i64>],
    b: &[Vec<i64>],
) -> Vec<Vec<i64>> {
    let cols = b.first().map_or(0, Vec::len);
    a.iter()
        .map(|row| (0..cols).map(|j| row.iter().zip(b).map(|(x, brow)| x * brow[j]).sum()).collect())
        .collect()
}

/// The transpose of an integer matrix with `cols` columns.
#[must_use]
pub fn transpose(
    m: &[Vec<i64>],
    cols: usize,
) -> Vec<Vec<i64>> {
    (0..cols).map(|j| m.iter().map(|row| row[j]).collect()).collect()
}

/// The rank of an integer matrix, computed exactly over the rationals by
/// Gaussian elimination.
#[must_use]
pub fn integer_rank(m: &[Vec<i64>]) -> usize {
    let mut a: Vec<Vec<BigRational>> =
        m.iter().map(|r| r.iter().map(|&x| BigRational::from_integer(x.into())).collect()).collect();
    let cols = a.first().map_or(0, Vec::len);
    let mut rank = 0;
    for c in 0..cols {
        let Some(p) = (rank..a.len()).find(|&r| !a[r][c].is_zero()) else { continue };
        a.swap(rank, p);
        let pivot = a[rank][c].clone();
        for r in (rank + 1)..a.len() {
            if a[r][c].is_zero() {
                continue;
            }
            let f = &a[r][c] / &pivot;
            #[allow(clippy::needless_range_loop)] // index is used for more than one array / arithmetic; iterator form would not be clearer
            for k in c..cols {
                let d = &f * &a[rank][k];
                a[r][k] -= d;
            }
        }
        rank += 1;
        if rank == a.len() {
            break;
        }
    }
    rank
}

/// The homology Betti number `b_k` of a closed complex:
/// `dim ker d_k - rank d_{k+1}`.
#[must_use]
pub fn betti_number(
    complex: &[Simplex],
    k: usize,
) -> usize {
    let sk = simplices_of_dim(complex, k);
    if sk.is_empty() {
        return 0;
    }
    let rank_k = if k == 0 {
        0
    } else {
        integer_rank(&boundary_matrix(&simplices_of_dim(complex, k - 1), &sk))
    };
    let above = simplices_of_dim(complex, k + 1);
    let rank_up = if above.is_empty() { 0 } else { integer_rank(&boundary_matrix(&sk, &above)) };
    sk.len().saturating_sub(rank_k).saturating_sub(rank_up)
}

/// The cohomology Betti number `b^k`, from the ranks of the coboundary
/// matrices (the transposes of the boundary matrices).
#[must_use]
pub fn cohomology_betti_number(
    complex: &[Simplex],
    k: usize,
) -> usize {
    let sk = simplices_of_dim(complex, k);
    if sk.is_empty() {
        return 0;
    }
    let above = simplices_of_dim(complex, k + 1);
    let rank_d = if above.is_empty() {
        0
    } else {
        integer_rank(&transpose(&boundary_matrix(&sk, &above), above.len()))
    };
    let rank_prev = if k == 0 {
        0
    } else {
        let below = simplices_of_dim(complex, k - 1);
        integer_rank(&transpose(&boundary_matrix(&below, &sk), sk.len()))
    };
    sk.len().saturating_sub(rank_d).saturating_sub(rank_prev)
}

/// Whether `d_k d_{k+1} = 0` for every `k`.
#[must_use]
pub fn verify_boundary_property(complex: &[Simplex]) -> bool {
    check_squares(complex, false)
}

/// Whether the coboundaries satisfy `d^{k+1} d^k = 0` for every `k`.
#[must_use]
pub fn verify_coboundary_property(complex: &[Simplex]) -> bool {
    check_squares(complex, true)
}

fn check_squares(
    complex: &[Simplex],
    dual: bool,
) -> bool {
    let Some(top) = complex_dimension(complex) else { return true };
    (1..top).all(|k| {
        let (a, b, c) = (
            simplices_of_dim(complex, k - 1),
            simplices_of_dim(complex, k),
            simplices_of_dim(complex, k + 1),
        );
        let dk = boundary_matrix(&a, &b);
        let dk1 = boundary_matrix(&b, &c);
        let product = if dual {
            int_matmul(&transpose(&dk1, c.len()), &transpose(&dk, b.len()))
        } else {
            int_matmul(&dk, &dk1)
        };
        product.iter().all(|r| r.iter().all(|&x| x == 0))
    })
}

/// The triangulated `width x height` grid: two triangles per cell.
#[must_use]
pub fn grid_complex(
    width: usize,
    height: usize,
) -> Vec<Simplex> {
    let mut triangles = Vec::new();
    for i in 0..height {
        for j in 0..width {
            let v0 = i * (width + 1) + j;
            let v1 = v0 + 1;
            let v2 = (i + 1) * (width + 1) + j;
            let v3 = v2 + 1;
            triangles.push(vec![v0, v1, v2]);
            triangles.push(vec![v1, v3, v2]);
        }
    }
    close_complex(&triangles)
}

/// The triangulated torus obtained from an `m x n` grid by identifying
/// opposite edges (a valid simplicial complex for `m, n >= 3`).
#[must_use]
pub fn torus_complex(
    m: usize,
    n: usize,
) -> Vec<Simplex> {
    let mut triangles = Vec::new();
    for i in 0..m {
        for j in 0..n {
            let v0 = i * n + j;
            let v1 = i * n + (j + 1) % n;
            let v2 = ((i + 1) % m) * n + j;
            let v3 = ((i + 1) % m) * n + (j + 1) % n;
            triangles.push(vec![v0, v1, v2]);
            triangles.push(vec![v1, v3, v2]);
        }
    }
    close_complex(&triangles)
}

fn radius(
    max_epsilon: f64,
    step: usize,
    steps: usize,
) -> f64 {
    if steps == 0 { max_epsilon } else { max_epsilon * (step as f64 / steps as f64) }
}

/// The Vietoris-Rips filtration restricted to its 1-skeleton: for each of
/// `steps + 1` radii `max_epsilon * s / steps` the vertices and the edges
/// no longer than the radius.
#[must_use]
pub fn vietoris_rips_filtration(
    points: &[Vec<f64>],
    max_epsilon: f64,
    steps: usize,
) -> Vec<(f64, Vec<Simplex>)> {
    (0..=steps)
        .map(|step| {
            let eps = radius(max_epsilon, step, steps);
            (eps, close_complex(&vietoris_rips_complex(points, eps, 1)))
        })
        .collect()
}

/// The Betti numbers `b_0 ..= b_max_dim` of the Vietoris-Rips complex at
/// radius `epsilon`.
#[must_use]
pub fn betti_numbers_at_radius(
    points: &[Vec<f64>],
    epsilon: f64,
    max_dim: usize,
) -> Vec<usize> {
    let complex = close_complex(&vietoris_rips_complex(points, epsilon, max_dim));
    (0..=max_dim).map(|k| betti_number(&complex, k)).collect()
}

/// Errors of the exact persistent-homology routines.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PersistenceError {
    /// A simplex has a face missing from the filtration, or a face that
    /// enters later than the simplex itself.
    InvalidFiltration,
    /// The Vietoris-Rips filtration would need more simplices than allowed.
    TooManySimplices,
}

/// Default cap on the number of simplices of a Vietoris-Rips filtration.
pub const DEFAULT_SIMPLEX_LIMIT: usize = 20_000_000;

fn xor_sorted(
    a: &[usize],
    b: &[usize],
) -> Vec<usize> {
    let mut out = Vec::with_capacity(a.len() + b.len());
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => {
                out.push(a[i]);
                i += 1;
            }
            std::cmp::Ordering::Greater => {
                out.push(b[j]);
                j += 1;
            }
            std::cmp::Ordering::Equal => {
                i += 1;
                j += 1;
            }
        }
    }
    out.extend_from_slice(&a[i..]);
    out.extend_from_slice(&b[j..]);
    out
}

/// Exact persistent homology over `Z/2` of a filtered simplicial complex.
///
/// `filtration` lists `(value, simplex)` pairs in any order; every proper
/// face of a simplex must be present with a value no larger than the
/// simplex's. The simplices are ordered by value, then dimension, then
/// lexicographically, and the boundary matrix is reduced by the standard
/// left-to-right column algorithm with *clearing* (the twist
/// optimisation): dimensions are processed from the top down, and a column
/// that is the pivot of a reduced column of the next dimension is known to
/// reduce to zero and is skipped.
///
/// Returns one [`PersistenceDiagram`] per dimension `0..=max_dim`; each
/// interval is an exact `(birth, death)` pair of filtration values and a
/// class that never dies has `death = f64::INFINITY`. Zero-length
/// intervals are omitted. Intervals are sorted by birth then death.
///
/// # Errors
/// [`PersistenceError::InvalidFiltration`] if the input is not a filtration.
pub fn persistent_homology(
    filtration: &[(f64, Simplex)],
    max_dim: usize,
) -> Result<Vec<PersistenceDiagram>, PersistenceError> {
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
    let m = items.len();
    let index: HashMap<&Simplex, usize> = items.iter().enumerate().map(|(i, (_, s))| (s, i)).collect();
    let dim_of = |i: usize| items[i].1.len() - 1;

    // Boundary columns (sorted row indices), grouped by dimension.
    let mut columns: Vec<Vec<usize>> = vec![Vec::new(); m];
    let mut by_dim: Vec<Vec<usize>> = vec![Vec::new(); max_dim + 2];
    for j in 0..m {
        let s = &items[j].1;
        let d = s.len() - 1;
        if d <= max_dim + 1 {
            by_dim[d].push(j);
        }
        if s.len() > 1 && d <= max_dim + 1 {
            let mut col = Vec::with_capacity(s.len());
            for skip in 0..s.len() {
                let face: Simplex = s.iter().enumerate().filter(|&(i, _)| i != skip).map(|(_, &v)| v).collect();
                let i = *index.get(&face).ok_or(PersistenceError::InvalidFiltration)?;
                if i >= j {
                    return Err(PersistenceError::InvalidFiltration);
                }
                col.push(i);
            }
            col.sort_unstable();
            columns[j] = col;
        }
    }

    let mut pivot_owner: Vec<Option<usize>> = vec![None; m];
    let mut paired = vec![false; m];
    let mut cleared = vec![false; m];
    let mut diagrams: Vec<PersistenceDiagram> =
        (0..=max_dim).map(|d| PersistenceDiagram { dimension: d, intervals: Vec::new() }).collect();

    for d in (1..=max_dim + 1).rev() {
        for &j in &by_dim[d] {
            if cleared[j] {
                continue;
            }
            while let Some(&low) = columns[j].last() {
                if let Some(other) = pivot_owner[low] {
                    let reduced = xor_sorted(&columns[j], &columns[other]);
                    columns[j] = reduced;
                    continue;
                }
                pivot_owner[low] = Some(j);
                paired[low] = true;
                paired[j] = true;
                cleared[low] = true;
                let (birth, death) = (items[low].0, items[j].0);
                if death > birth {
                    diagrams[dim_of(low)].intervals.push(PersistenceInterval { birth, death });
                }
                break;
            }
        }
    }
    for i in 0..m {
        let d = dim_of(i);
        if !paired[i] && d <= max_dim {
            diagrams[d].intervals.push(PersistenceInterval { birth: items[i].0, death: f64::INFINITY });
        }
    }
    for dg in &mut diagrams {
        dg.intervals.sort_by(|a, b| a.birth.total_cmp(&b.birth).then(a.death.total_cmp(&b.death)));
    }
    Ok(diagrams)
}

/// Exact Vietoris-Rips persistent homology of a point cloud.
///
/// The filtration contains every simplex up to dimension `max_dim + 1`
/// whose edges are all at most `max_epsilon` long, entering at its diameter
/// (the longest edge); it is processed by [`persistent_homology`], so every
/// `(birth, death)` pair is exact (a distance between two input points).
/// Classes alive at `max_epsilon` get `death = f64::INFINITY`.
///
/// # Errors
/// [`PersistenceError::TooManySimplices`] if more than `simplex_limit`
/// simplices are needed (use [`DEFAULT_SIMPLEX_LIMIT`] when in doubt).
pub fn persistent_homology_rips(
    points: &[Vec<f64>],
    max_epsilon: f64,
    max_dim: usize,
    simplex_limit: usize,
) -> Result<Vec<PersistenceDiagram>, PersistenceError> {
    let dist = super::homology::distance_matrix(points);
    let filtration = super::homology::rips_filtration(&dist, max_epsilon, max_dim + 1, simplex_limit)
        .ok_or(PersistenceError::TooManySimplices)?;
    persistent_homology(&filtration, max_dim)
}

/// The persistence diagrams of a point cloud, `max_dim` and below.
///
/// This is a compatibility wrapper around [`persistent_homology_rips`]: the
/// diagrams are now exact, not sampled on a grid of radii, so `steps` is
/// ignored. Intervals of classes still alive at `max_epsilon` are closed at
/// `max_epsilon` (the convention of the former grid-based routine; use
/// [`persistent_homology_rips`] to get `f64::INFINITY` instead). If the
/// filtration is too large to build, empty diagrams are returned.
#[must_use]
pub fn compute_persistence(
    points: &[Vec<f64>],
    max_epsilon: f64,
    _steps: usize,
    max_dim: usize,
) -> Vec<PersistenceDiagram> {
    let mut diagrams = persistent_homology_rips(points, max_epsilon, max_dim, DEFAULT_SIMPLEX_LIMIT)
        .unwrap_or_else(|_| (0..=max_dim).map(|d| PersistenceDiagram { dimension: d, intervals: Vec::new() }).collect());
    for dg in &mut diagrams {
        for iv in &mut dg.intervals {
            if iv.death.is_infinite() {
                iv.death = max_epsilon;
            }
        }
    }
    diagrams
}
