//! # Numerical Computational Topology
//!
//! Pure numeric algorithms for computational topology: Euclidean distances,
//! Vietoris-Rips complexes, exact Betti numbers of simplicial complexes
//! (ranks of integer boundary matrices are computed exactly over the
//! rationals), naive persistence diagrams and connected components of a
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

/// The (naive) persistence diagrams of a point cloud.
///
/// Betti numbers are
/// sampled at `steps + 1` radii and a rise of `b_k` opens an interval, a
/// fall closes the most recently opened one; intervals still open at
/// `max_epsilon` end there.
#[must_use]
pub fn compute_persistence(
    points: &[Vec<f64>],
    max_epsilon: f64,
    steps: usize,
    max_dim: usize,
) -> Vec<PersistenceDiagram> {
    let mut diagrams: Vec<PersistenceDiagram> =
        (0..=max_dim).map(|d| PersistenceDiagram { dimension: d, intervals: Vec::new() }).collect();
    let mut prev = vec![0; max_dim + 1];
    let mut open: Vec<Vec<f64>> = vec![Vec::new(); max_dim + 1];
    for step in 0..=steps {
        let eps = radius(max_epsilon, step, steps);
        let current = betti_numbers_at_radius(points, eps, max_dim);
        for d in 0..=max_dim {
            if current[d] > prev[d] {
                for _ in 0..(current[d] - prev[d]) {
                    open[d].push(eps);
                }
            } else {
                for _ in 0..(prev[d] - current[d]) {
                    if let Some(birth) = open[d].pop() {
                        diagrams[d].intervals.push(PersistenceInterval { birth, death: eps });
                    }
                }
            }
        }
        prev = current;
    }
    for d in 0..=max_dim {
        while let Some(birth) = open[d].pop() {
            diagrams[d].intervals.push(PersistenceInterval { birth, death: max_epsilon });
        }
    }
    diagrams
}
