//! # Helmholtz boundary elements: exterior acoustic scattering (Burton–Miller)
//!
//! Scattering of a plane wave `p_inc = e^{i k d·x}` (time dependence `e^{-iωt}`)
//! by a closed triangulated body, with piecewise-constant elements and
//! collocation at the triangle centroids. The free-space Green function is
//! `G = e^{ikr} / (4π r)`.
//!
//! The classical conventional boundary integral equation fails at the
//! irregular frequencies (the eigenfrequencies of the interior Dirichlet or
//! Neumann problem). The default formulation is the **Burton–Miller** one,
//! the conventional equation plus `α` times its normal derivative with
//! `α = i / k` (unique for every `k`):
//!
//! * sound-hard (`∂p/∂n = 0`): `(½ - K + α N) p = p_inc - α ∂_n p_inc`, unknown
//!   total surface pressure `p`;
//! * sound-soft (`p = 0`): `(S + α (½ + K')) q = p_inc + α ∂_n p_inc`, unknown
//!   normal derivative `q = ∂_n p`,
//!
//! with the single layer `S`, double layer `K`, its adjoint `K'` and the
//! hypersingular operator `N`. Off-diagonal entries use a 7-point rule with
//! adaptive subdivision near the collocation point; on the (flat) own
//! element `K` and `K'` vanish and `S` and the Hadamard finite part of `N`
//! are integrated in polar coordinates about the centroid. `α = 0` gives the
//! conventional equation (to demonstrate the irregular frequencies).
//!
//! The dense complex system is solved by LU; the accuracy for a sphere is
//! first order in the element size (flat facets), compared with the Mie
//! series in the tests.
//!
//! [`solve_helmholtz_exterior_aca`] solves the same system without forming the matrix: the
//! collocation points form a cluster tree (recursive bisection), the admissible blocks
//! `min(diam_s, diam_t) <= η dist(s, t)` are compressed by **complex adaptive cross approximation**
//! with partial pivoting, the near field is stored densely, and the system is solved by
//! restarted complex GMRES (Jacobi preconditioned) with hierarchical products.
//!
//! Limitations: constant elements only, closed orientable surfaces with outward normals, no
//! impedance boundary conditions; the compressed solver pays off from a few thousand elements
//! and only while the body is not much larger than a few wavelengths (the ranks of the
//! oscillatory kernel grow with `k·diam`, there is no directional/FMM-type compression);
//! `solve_helmholtz_exterior` is dense `O(n²)`.

use crate::sim::physics_bem::SurfaceMesh3D;
use num_complex::Complex64;
use rayon::prelude::*;
use std::f64::consts::PI;

type V3 = [f64; 3];

/// Boundary condition of the scatterer.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Scatterer {
    /// Rigid, `∂p/∂n = 0`.
    SoundHard,
    /// Pressure release, `p = 0`.
    SoundSoft,
}

/// Plane incident wave `e^{i k d·x}`.
#[derive(Clone, Copy, Debug)]
pub struct PlaneWave {
    /// Propagation direction (normalised internally).
    pub direction: V3,
}

/// Solution of the exterior scattering problem.
#[derive(Clone, Debug)]
pub struct HelmholtzSolution {
    /// Wave number.
    pub k: f64,
    /// Boundary condition.
    pub kind: Scatterer,
    /// Surface unknown per triangle: the total pressure (hard) or the normal derivative (soft).
    pub density: Vec<Complex64>,
    dir: V3,
}

fn sub(
    a: V3,
    b: V3,
) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(
    a: V3,
    b: V3,
) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(
    a: V3,
    b: V3,
) -> V3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn norm(a: V3) -> f64 {
    dot(a, a).sqrt()
}

#[derive(Clone, Copy)]
enum Kernel {
    Single,
    Double,
    Adjoint,
    Hyper,
}

const RADON_A: f64 = 0.101_286_507_323_456_3;
const RADON_B: f64 = 0.470_142_064_105_115_1;

fn radon7() -> [([f64; 3], f64); 7] {
    let (wa, wb) = (2.0 * 0.062_969_590_272_413_6, 2.0 * 0.066_197_076_394_253_1);
    [
        ([1.0 / 3.0; 3], 0.225),
        ([RADON_A, RADON_A, 1.0 - 2.0 * RADON_A], wa),
        ([RADON_A, 1.0 - 2.0 * RADON_A, RADON_A], wa),
        ([1.0 - 2.0 * RADON_A, RADON_A, RADON_A], wa),
        ([RADON_B, RADON_B, 1.0 - 2.0 * RADON_B], wb),
        ([RADON_B, 1.0 - 2.0 * RADON_B, RADON_B], wb),
        ([1.0 - 2.0 * RADON_B, RADON_B, RADON_B], wb),
    ]
}

/// Value of a kernel for source point `y` (normal `ny`) and field point `x` (normal `nx`).
fn kernel_value(
    kind: Kernel,
    k: f64,
    x: V3,
    nx: V3,
    y: V3,
    ny: V3,
) -> Complex64 {
    let d = sub(y, x);
    let r = norm(d);
    let ikr = Complex64::new(0.0, k * r);
    let g = ikr.exp() / (4.0 * PI * r);
    match kind {
        | Kernel::Single => g,
        | Kernel::Double | Kernel::Adjoint | Kernel::Hyper => {
            let gp = g * (Complex64::new(0.0, k) - 1.0 / r);
            let dny = dot(d, ny) / r;
            let dnx = -dot(d, nx) / r;
            match kind {
                | Kernel::Double => gp * dny,
                | Kernel::Adjoint => gp * dnx,
                | _ => {
                    let gpp = g
                        * (Complex64::new(-k * k, 0.0) - Complex64::new(0.0, 2.0 * k) / r
                            + 2.0 / (r * r));
                    let d2r = -dot(nx, ny) / r + dot(d, ny) * dot(d, nx) / (r * r * r);
                    gpp * dnx * dny + gp * d2r
                },
            }
        },
    }
}

/// Integral of a kernel over triangle `tri` by 7-point rules with recursive subdivision.
fn integrate(
    kind: Kernel,
    k: f64,
    x: V3,
    nx: V3,
    tri: [V3; 3],
    ny: V3,
    depth: usize,
) -> Complex64 {
    let centre = [
        (tri[0][0] + tri[1][0] + tri[2][0]) / 3.0,
        (tri[0][1] + tri[1][1] + tri[2][1]) / 3.0,
        (tri[0][2] + tri[1][2] + tri[2][2]) / 3.0,
    ];
    let diam = norm(sub(tri[0], tri[1]))
        .max(norm(sub(tri[1], tri[2])))
        .max(norm(sub(tri[2], tri[0])));
    let dist = norm(sub(centre, x));
    let limit = match kind {
        | Kernel::Hyper => 3.0,
        | _ => 2.0,
    };
    let max_depth = match kind {
        | Kernel::Hyper => 5,
        | _ => 4,
    };
    if dist < limit * diam && depth < max_depth {
        let mid = |a: V3, b: V3| {
            [
                f64::midpoint(a[0], b[0]),
                f64::midpoint(a[1], b[1]),
                f64::midpoint(a[2], b[2]),
            ]
        };
        let (m01, m12, m20) = (
            mid(tri[0], tri[1]),
            mid(tri[1], tri[2]),
            mid(tri[2], tri[0]),
        );
        return [
            [tri[0], m01, m20],
            [m01, tri[1], m12],
            [m20, m12, tri[2]],
            [m12, m20, m01],
        ]
        .iter()
        .map(|t| integrate(kind, k, x, nx, *t, ny, depth + 1))
        .sum();
    }
    let area = 0.5 * norm(cross(sub(tri[1], tri[0]), sub(tri[2], tri[0])));
    let mut s = Complex64::new(0.0, 0.0);
    for (b, w) in radon7() {
        let y: V3 = std::array::from_fn(|a| b[0] * tri[0][a] + b[1] * tri[1][a] + b[2] * tri[2][a]);
        s += kernel_value(kind, k, x, nx, y, ny) * (area * w);
    }
    s
}

/// Gauss–Legendre nodes and weights on `[-1, 1]` (8 points).
const GL8_X: [f64; 8] = [
    -0.960_289_856_497_536_3,
    -0.796_666_477_413_626_7,
    -0.525_532_409_916_329,
    -0.183_434_642_495_649_8,
    0.183_434_642_495_649_8,
    0.525_532_409_916_329,
    0.796_666_477_413_626_7,
    0.960_289_856_497_536_3,
];
const GL8_W: [f64; 8] = [
    0.101_228_536_290_376_26,
    0.222_381_034_453_374_5,
    0.313_706_645_877_887_3,
    0.362_683_783_378_362,
    0.362_683_783_378_362,
    0.313_706_645_877_887_3,
    0.222_381_034_453_374_5,
    0.101_228_536_290_376_26,
];

/// Self terms for the flat triangle `tri` with the collocation point at its centroid:
/// `(∫ G dS, FP ∫ ∂²G/∂n∂n' dS)`.
fn self_terms(
    k: f64,
    tri: [V3; 3],
) -> (Complex64, Complex64) {
    let c: V3 = std::array::from_fn(|a| (tri[0][a] + tri[1][a] + tri[2][a]) / 3.0);
    let mut s_int = Complex64::new(0.0, 0.0);
    let mut n_static = 0.0;
    let mut n_diff = Complex64::new(0.0, 0.0);
    let ik = Complex64::new(0.0, k);
    for e in 0..3 {
        let (a, b) = (tri[e], tri[(e + 1) % 3]);
        let edge = sub(b, a);
        let len = norm(edge);
        let t = [edge[0] / len, edge[1] / len, edge[2] / len];
        // foot of the perpendicular from the centroid
        let s0 = dot(sub(c, a), t);
        let foot: V3 = std::array::from_fn(|i| a[i] + s0 * t[i]);
        let dperp = norm(sub(foot, c));
        let (s1, s2) = (-s0, len - s0);
        let (th1, th2) = ((s1 / dperp).atan(), (s2 / dperp).atan());
        let (half, mid) = (0.5 * (th2 - th1), f64::midpoint(th2, th1));
        for q in 0..8 {
            let th = mid + half * GL8_X[q];
            let wth = half * GL8_W[q];
            let rmax = dperp / th.cos();
            // ∫_0^R e^{ikr}/(4π) dr
            let int_s = if k * rmax < 1e-8 {
                Complex64::new(rmax, 0.0)
            } else {
                ((ik * rmax).exp() - 1.0) / ik
            };
            s_int += int_s * (wth / (4.0 * PI));
            n_static += wth / rmax;
            // ∫_0^R [e^{ikr}(1 - ikr) - 1] / (4π r²) dr by Gauss-Legendre in r
            let mut acc = Complex64::new(0.0, 0.0);
            for m in 0..8 {
                let r = 0.5 * rmax * (1.0 + GL8_X[m]);
                let f = if k * r < 0.1 {
                    let z = ik * r;
                    let mut sum = Complex64::new(0.0, 0.0);
                    let mut zn = z * z;
                    let mut fact = 2.0;
                    for n in 2..9 {
                        sum += zn * ((1.0 - f64::from(n)) / fact);
                        zn *= z;
                        fact *= f64::from(n + 1);
                    }
                    sum / (r * r)
                } else {
                    ((ik * r).exp() * (1.0 - ik * r) - 1.0) / (r * r)
                };
                acc += f * (0.5 * rmax * GL8_W[m]);
            }
            n_diff += acc * (wth / (4.0 * PI));
        }
    }
    (s_int, n_diff - Complex64::new(n_static / (4.0 * PI), 0.0))
}

fn plane_wave(
    k: f64,
    dir: V3,
    x: V3,
) -> Complex64 {
    Complex64::new(0.0, k * dot(dir, x)).exp()
}

/// Dense complex linear solve (LU with partial pivoting).
fn solve_dense(
    mut a: Vec<Complex64>,
    mut b: Vec<Complex64>,
    n: usize,
) -> Result<Vec<Complex64>, String> {
    for c in 0..n {
        let mut piv = c;
        let mut best = a[c * n + c].norm();
        for r in c + 1..n {
            if a[r * n + c].norm() > best {
                best = a[r * n + c].norm();
                piv = r;
            }
        }
        if best < 1e-300 {
            return Err("singular BEM system".into());
        }
        if piv != c {
            for j in 0..n {
                a.swap(c * n + j, piv * n + j);
            }
            b.swap(c, piv);
        }
        let d = a[c * n + c];
        let (head, tail) = a.split_at_mut((c + 1) * n);
        let row_c = &head[c * n..(c + 1) * n];
        let bc = b[c];
        let updates: Vec<(usize, Complex64)> = tail
            .par_chunks_mut(n)
            .enumerate()
            .map(|(ri, row)| {
                let f = row[c] / d;
                if f.norm() > 0.0 {
                    for j in c..n {
                        row[j] -= f * row_c[j];
                    }
                }
                (c + 1 + ri, f)
            })
            .collect();
        for (r, f) in updates {
            b[r] -= f * bc;
        }
    }
    for c in (0..n).rev() {
        let mut s = b[c];
        for j in c + 1..n {
            s -= a[c * n + j] * b[j];
        }
        b[c] = s / a[c * n + c];
    }
    Ok(b)
}

/// Geometry and formulation of one scattering problem: the entries of the collocation matrix and
/// the right-hand side, shared by the dense and the compressed solver.
struct Setup {
    k: f64,
    kind: Scatterer,
    alpha: Complex64,
    tris: Vec<[V3; 3]>,
    normals: Vec<V3>,
    cents: Vec<V3>,
}

impl Setup {
    /// Entry `(i, j)` of the collocation matrix.
    fn entry(
        &self,
        i: usize,
        j: usize,
    ) -> Complex64 {
        let (k, alpha) = (self.k, self.alpha);
        let (x, nx) = (self.cents[i], self.normals[i]);
        let own = i == j;
        let (s_ij, n_ij) = if own {
            self_terms(k, self.tris[j])
        } else {
            (Complex64::new(0.0, 0.0), Complex64::new(0.0, 0.0))
        };
        match self.kind {
            | Scatterer::SoundHard => {
                let kk = if own {
                    Complex64::new(0.0, 0.0)
                } else {
                    integrate(Kernel::Double, k, x, nx, self.tris[j], self.normals[j], 0)
                };
                let nn = if own {
                    n_ij
                } else if alpha.norm() > 0.0 {
                    integrate(Kernel::Hyper, k, x, nx, self.tris[j], self.normals[j], 0)
                } else {
                    Complex64::new(0.0, 0.0)
                };
                let half = if own { 0.5 } else { 0.0 };
                Complex64::new(half, 0.0) - kk + alpha * nn
            },
            | Scatterer::SoundSoft => {
                let ss = if own {
                    s_ij
                } else {
                    integrate(Kernel::Single, k, x, nx, self.tris[j], self.normals[j], 0)
                };
                let kp = if own {
                    Complex64::new(0.0, 0.0)
                } else {
                    integrate(Kernel::Adjoint, k, x, nx, self.tris[j], self.normals[j], 0)
                };
                let half = if own { 0.5 } else { 0.0 };
                ss + alpha * (Complex64::new(half, 0.0) + kp)
            },
        }
    }

    /// Right-hand side for the incident direction `dir`.
    fn rhs(
        &self,
        dir: V3,
    ) -> Vec<Complex64> {
        (0..self.cents.len())
            .map(|i| {
                let pinc = plane_wave(self.k, dir, self.cents[i]);
                let dpinc = pinc * Complex64::new(0.0, self.k * dot(dir, self.normals[i]));
                match self.kind {
                    | Scatterer::SoundHard => pinc - self.alpha * dpinc,
                    | Scatterer::SoundSoft => pinc + self.alpha * dpinc,
                }
            })
            .collect()
    }
}

/// Options of [`solve_helmholtz_exterior_aca`].
#[derive(Clone, Copy, Debug)]
pub struct HelmholtzAcaOptions {
    /// Maximum number of collocation points in a leaf cluster.
    pub leaf_size: usize,
    /// Admissibility parameter `η`.
    pub eta: f64,
    /// Relative ACA tolerance of the low-rank blocks.
    pub aca_tol: f64,
    /// GMRES relative tolerance.
    pub gmres_tol: f64,
    /// GMRES restart length.
    pub restart: usize,
    /// Maximum number of GMRES iterations.
    pub max_iter: usize,
}

impl Default for HelmholtzAcaOptions {
    fn default() -> Self {
        Self {
            leaf_size: 32,
            eta: 1.0,
            aca_tol: 1e-6,
            gmres_tol: 1e-8,
            restart: 80,
            max_iter: 800,
        }
    }
}

/// Storage statistics of the compressed matrix.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HelmholtzAcaStats {
    /// Number of low-rank blocks.
    pub low_rank_blocks: usize,
    /// Number of dense blocks.
    pub dense_blocks: usize,
    /// Stored complex numbers.
    pub stored: usize,
    /// `stored / n²`.
    pub compression: f64,
    /// GMRES iterations used.
    pub iterations: usize,
}

enum HBlock {
    Dense {
        rows: Vec<usize>,
        cols: Vec<usize>,
        data: Vec<Complex64>,
    },
    LowRank {
        rows: Vec<usize>,
        cols: Vec<usize>,
        /// `u[l]` has `rows.len()` entries, `v[l]` has `cols.len()`; the block is `sum_l u_l v_l^T`.
        u: Vec<Vec<Complex64>>,
        v: Vec<Vec<Complex64>>,
    },
}

/// Cluster of collocation indices with its bounding box and children.
struct HCluster {
    idx: Vec<usize>,
    lo: V3,
    hi: V3,
    kids: Option<Box<(Self, Self)>>,
}

fn build_cluster(
    idx: Vec<usize>,
    pts: &[V3],
    leaf: usize,
) -> HCluster {
    let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for &i in &idx {
        for a in 0..3 {
            lo[a] = lo[a].min(pts[i][a]);
            hi[a] = hi[a].max(pts[i][a]);
        }
    }
    if idx.len() <= leaf {
        return HCluster {
            idx,
            lo,
            hi,
            kids: None,
        };
    }
    let axis = (0..3)
        .max_by(|&a, &b| (hi[a] - lo[a]).total_cmp(&(hi[b] - lo[b])))
        .unwrap_or(0);
    let mut sorted = idx.clone();
    sorted.sort_by(|&a, &b| pts[a][axis].total_cmp(&pts[b][axis]).then(a.cmp(&b)));
    let right = sorted.split_off(sorted.len() / 2);
    let kids = Box::new((
        build_cluster(sorted, pts, leaf),
        build_cluster(right, pts, leaf),
    ));
    HCluster {
        idx,
        lo,
        hi,
        kids: Some(kids),
    }
}

fn cluster_diam(c: &HCluster) -> f64 {
    norm(sub(c.hi, c.lo))
}

fn cluster_dist(
    a: &HCluster,
    b: &HCluster,
) -> f64 {
    let mut d2 = 0.0;
    for k in 0..3 {
        let gap = (a.lo[k] - b.hi[k]).max(b.lo[k] - a.hi[k]).max(0.0);
        d2 += gap * gap;
    }
    d2.sqrt()
}

/// Block pairs of the hierarchical partition: `(rows, cols, admissible)`.
fn partition_pairs<'a>(
    s: &'a HCluster,
    t: &'a HCluster,
    eta: f64,
    out: &mut Vec<(&'a [usize], &'a [usize], bool)>,
) {
    if cluster_diam(s).min(cluster_diam(t)) <= eta * cluster_dist(s, t) {
        out.push((&s.idx, &t.idx, true));
        return;
    }
    match (&s.kids, &t.kids) {
        | (None, None) => out.push((&s.idx, &t.idx, false)),
        | (Some(ks), None) => {
            partition_pairs(&ks.0, t, eta, out);
            partition_pairs(&ks.1, t, eta, out);
        },
        | (None, Some(kt)) => {
            partition_pairs(s, &kt.0, eta, out);
            partition_pairs(s, &kt.1, eta, out);
        },
        | (Some(ks), Some(kt)) => {
            for a in [&ks.0, &ks.1] {
                for b in [&kt.0, &kt.1] {
                    partition_pairs(a, b, eta, out);
                }
            }
        },
    }
}

fn herm_dot(
    a: &[Complex64],
    b: &[Complex64],
) -> Complex64 {
    a.iter().zip(b).map(|(x, y)| x.conj() * y).sum()
}

fn vec_norm(a: &[Complex64]) -> f64 {
    a.iter().map(Complex64::norm_sqr).sum::<f64>().sqrt()
}

/// Factors `(u, v)` of a low-rank block `sum_l u_l v_l^T`.
type CrossFactors = (Vec<Vec<Complex64>>, Vec<Vec<Complex64>>);

/// Adaptive cross approximation with partial pivoting of the `m x n` block with entries
/// `entry(i, j)`; returns `(u, v)` with block `~ sum_l u_l v_l^T`, or `None` when the rank reaches
/// `min(m, n) / 2 + 1` (the block is then stored densely).
fn aca_complex(
    m: usize,
    n: usize,
    tol: f64,
    entry: &dyn Fn(usize, usize) -> Complex64,
) -> Option<CrossFactors> {
    let max_rank = m.min(n) / 2 + 1;
    let (mut us, mut vs): (Vec<Vec<Complex64>>, Vec<Vec<Complex64>>) = (Vec::new(), Vec::new());
    let (mut row_used, mut col_used) = (vec![false; m], vec![false; n]);
    let mut s2 = 0.0_f64;
    let mut i = 0;
    let mut misses = 0;
    while us.len() < max_rank {
        row_used[i] = true;
        let mut row: Vec<Complex64> = (0..n).map(|j| entry(i, j)).collect();
        for (u, v) in us.iter().zip(&vs) {
            let c = u[i];
            for (r, x) in row.iter_mut().zip(v) {
                *r -= c * x;
            }
        }
        let jstar = (0..n)
            .filter(|&j| !col_used[j])
            .max_by(|&a, &b| row[a].norm().total_cmp(&row[b].norm()));
        let Some(j) = jstar else { break };
        let pivot = row[j];
        if pivot.norm() < 1e-300 {
            // a vanishing residual row: try another row, give up after a few
            misses += 1;
            match (0..m).find(|&r| !row_used[r]) {
                | Some(next) if misses < 4 => {
                    i = next;
                    continue;
                },
                | _ => break,
            }
        }
        col_used[j] = true;
        let v: Vec<Complex64> = row.iter().map(|x| x / pivot).collect();
        let mut col: Vec<Complex64> = (0..m).map(|r| entry(r, j)).collect();
        for (u, vv) in us.iter().zip(&vs) {
            let c = vv[j];
            for (x, y) in col.iter_mut().zip(u) {
                *x -= c * y;
            }
        }
        let (nu, nv) = (vec_norm(&col), vec_norm(&v));
        let mut cross = Complex64::new(0.0, 0.0);
        for (u, vv) in us.iter().zip(&vs) {
            cross += herm_dot(u, &col) * herm_dot(vv, &v);
        }
        s2 += nu * nu * nv * nv + 2.0 * cross.re;
        let next = (0..m)
            .filter(|&r| !row_used[r])
            .max_by(|&a, &b| col[a].norm().total_cmp(&col[b].norm()));
        us.push(col);
        vs.push(v);
        if nu * nv <= tol * s2.max(0.0).sqrt() {
            return Some((us, vs));
        }
        match next {
            | Some(r) => i = r,
            | None => return Some((us, vs)),
        }
    }
    if us.len() >= max_rank {
        None
    } else {
        Some((us, vs))
    }
}

struct HMat {
    n: usize,
    blocks: Vec<HBlock>,
}

impl HMat {
    fn matvec(
        &self,
        x: &[Complex64],
    ) -> Vec<Complex64> {
        // each block adds into a private vector; the sum is taken in block order
        let parts: Vec<Vec<(usize, Complex64)>> =
            crate::kernels::par::map_slice(&self.blocks, 16, |b| match b {
                | HBlock::Dense { rows, cols, data } => {
                    let nc = cols.len();
                    rows.iter()
                        .enumerate()
                        .map(|(r, &gi)| {
                            let mut acc = Complex64::new(0.0, 0.0);
                            for (c, &gj) in cols.iter().enumerate() {
                                acc += data[r * nc + c] * x[gj];
                            }
                            (gi, acc)
                        })
                        .collect()
                },
                | HBlock::LowRank { rows, cols, u, v } => {
                    let mut acc = vec![Complex64::new(0.0, 0.0); rows.len()];
                    for (ul, vl) in u.iter().zip(v) {
                        let mut c = Complex64::new(0.0, 0.0);
                        for (a, &gj) in vl.iter().zip(cols) {
                            c += a * x[gj];
                        }
                        for (t, w) in acc.iter_mut().zip(ul) {
                            *t += c * w;
                        }
                    }
                    rows.iter().copied().zip(acc).collect()
                },
            });
        let mut y = vec![Complex64::new(0.0, 0.0); self.n];
        for part in parts {
            for (i, v) in part {
                y[i] += v;
            }
        }
        y
    }

    fn stored(&self) -> usize {
        self.blocks
            .iter()
            .map(|b| match b {
                | HBlock::Dense { data, .. } => data.len(),
                | HBlock::LowRank { rows, cols, u, .. } => u.len() * (rows.len() + cols.len()),
            })
            .sum()
    }
}

/// Restarted GMRES for a complex operator with a diagonal (Jacobi) left preconditioner.
/// Returns the solution and the iteration count; `None` when it did not converge.
fn gmres_complex(
    apply: &dyn Fn(&[Complex64]) -> Vec<Complex64>,
    diag: &[Complex64],
    b: &[Complex64],
    tol: f64,
    restart: usize,
    max_iter: usize,
) -> Option<(Vec<Complex64>, usize)> {
    let n = b.len();
    let prec = |v: &mut Vec<Complex64>| {
        for (x, d) in v.iter_mut().zip(diag) {
            if d.norm() > 0.0 {
                *x /= d;
            }
        }
    };
    let mut pb = b.to_vec();
    prec(&mut pb);
    let bnorm = vec_norm(&pb).max(1e-300);
    let mut x = vec![Complex64::new(0.0, 0.0); n];
    let mut total = 0;
    while total < max_iter {
        let mut r = apply(&x);
        for (ri, bi) in r.iter_mut().zip(b) {
            *ri = bi - *ri;
        }
        prec(&mut r);
        let beta = vec_norm(&r);
        if beta <= tol * bnorm {
            return Some((x, total));
        }
        let m = restart.min(n);
        let mut basis: Vec<Vec<Complex64>> = vec![r.iter().map(|v| v / beta).collect()];
        let mut h = vec![vec![Complex64::new(0.0, 0.0); m]; m + 1];
        let (mut cs, mut sn) = (vec![0.0_f64; m], vec![Complex64::new(0.0, 0.0); m]);
        let mut g = vec![Complex64::new(0.0, 0.0); m + 1];
        g[0] = Complex64::new(beta, 0.0);
        let mut used = 0;
        for j in 0..m {
            let mut w = apply(&basis[j]);
            prec(&mut w);
            for (i, q) in basis.iter().enumerate() {
                let hij = herm_dot(q, &w);
                h[i][j] = hij;
                for (wk, qk) in w.iter_mut().zip(q) {
                    *wk -= hij * qk;
                }
            }
            let hn = vec_norm(&w);
            h[j + 1][j] = Complex64::new(hn, 0.0);
            // apply the previous rotations to the new column
            for i in 0..j {
                let t = cs[i] * h[i][j] + sn[i] * h[i + 1][j];
                h[i + 1][j] = -sn[i].conj() * h[i][j] + cs[i] * h[i + 1][j];
                h[i][j] = t;
            }
            let (a, bb) = (h[j][j], h[j + 1][j]);
            let denom = (a.norm_sqr() + bb.norm_sqr()).sqrt();
            if denom > 0.0 {
                cs[j] = a.norm() / denom;
                sn[j] = if a.norm() > 0.0 {
                    (a / a.norm()) * bb.conj() / denom
                } else {
                    bb.conj() / denom
                };
            } else {
                cs[j] = 1.0;
            }
            h[j][j] = cs[j] * a + sn[j] * bb;
            h[j + 1][j] = Complex64::new(0.0, 0.0);
            let gj = g[j];
            g[j + 1] = -sn[j].conj() * gj;
            g[j] = cs[j] * gj;
            used = j + 1;
            total += 1;
            if hn > 1e-300 {
                basis.push(w.iter().map(|v| v / hn).collect());
            }
            if g[j + 1].norm() <= tol * bnorm || hn <= 1e-300 || total >= max_iter {
                break;
            }
        }
        // back substitution
        let mut y = vec![Complex64::new(0.0, 0.0); used];
        for i in (0..used).rev() {
            let mut sum = g[i];
            for k in i + 1..used {
                sum -= h[i][k] * y[k];
            }
            y[i] = sum / h[i][i];
        }
        for (i, yi) in y.iter().enumerate() {
            for (xk, qk) in x.iter_mut().zip(&basis[i]) {
                *xk += yi * qk;
            }
        }
    }
    let mut r = apply(&x);
    for (ri, bi) in r.iter_mut().zip(b) {
        *ri = bi - *ri;
    }
    prec(&mut r);
    (vec_norm(&r) <= 10.0 * tol * bnorm).then_some((x, total))
}

/// [`solve_helmholtz_exterior`] with a hierarchical matrix: complex ACA for the admissible blocks
/// and complex GMRES; no dense `n x n` matrix is formed. Returns the solution and storage
/// statistics.
///
/// # Errors
/// As [`solve_helmholtz_exterior`]; also when GMRES does not converge.
pub fn solve_helmholtz_exterior_aca(
    mesh: &SurfaceMesh3D,
    k: f64,
    kind: Scatterer,
    wave: PlaneWave,
    coupling: Option<Complex64>,
    opts: &HelmholtzAcaOptions,
) -> Result<(HelmholtzSolution, HelmholtzAcaStats), String> {
    mesh.validate()?;
    if !(k > 0.0 && k.is_finite()) {
        return Err("the wave number must be positive".into());
    }
    let dn = norm(wave.direction);
    if dn == 0.0 {
        return Err("the incident direction must be non-zero".into());
    }
    let dir = [
        wave.direction[0] / dn,
        wave.direction[1] / dn,
        wave.direction[2] / dn,
    ];
    let alpha = coupling.unwrap_or_else(|| Complex64::new(0.0, 1.0 / k));
    let nt = mesh.triangles.len();
    let tris: Vec<[V3; 3]> = (0..nt).filter_map(|t| mesh.corners(t)).collect();
    let normals: Vec<V3> = (0..nt).filter_map(|t| mesh.normal(t)).collect();
    let cents: Vec<V3> = (0..nt).filter_map(|t| mesh.centroid(t)).collect();
    if tris.len() != nt || normals.len() != nt || cents.len() != nt {
        return Err("degenerate triangle in the mesh".into());
    }
    let setup = Setup {
        k,
        kind,
        alpha,
        tris,
        normals,
        cents,
    };
    let root = build_cluster((0..nt).collect(), &setup.cents, opts.leaf_size.max(2));
    let mut pairs = Vec::new();
    partition_pairs(&root, &root, opts.eta, &mut pairs);
    let entry = |i: usize, j: usize| setup.entry(i, j);
    // blocks are independent: built in parallel, kept in pair order
    let blocks: Vec<HBlock> = crate::kernels::par::map_slice(&pairs, 4, |&(rows, cols, far)| {
        let dense = || {
            let data = rows
                .iter()
                .flat_map(|&i| cols.iter().map(move |&j| (i, j)))
                .map(|(i, j)| entry(i, j))
                .collect();
            HBlock::Dense {
                rows: rows.to_vec(),
                cols: cols.to_vec(),
                data,
            }
        };
        if !far {
            return dense();
        }
        let local = |a: usize, b: usize| entry(rows[a], cols[b]);
        match aca_complex(rows.len(), cols.len(), opts.aca_tol, &local) {
            | Some((u, v)) => HBlock::LowRank {
                rows: rows.to_vec(),
                cols: cols.to_vec(),
                u,
                v,
            },
            | None => dense(),
        }
    });
    let h = HMat { n: nt, blocks };
    let diag: Vec<Complex64> = (0..nt).map(|i| entry(i, i)).collect();
    let rhs = setup.rhs(dir);
    let (density, iterations) = gmres_complex(
        &|v| h.matvec(v),
        &diag,
        &rhs,
        opts.gmres_tol,
        opts.restart,
        opts.max_iter,
    )
    .ok_or("GMRES did not converge")?;
    let stats = HelmholtzAcaStats {
        low_rank_blocks: h
            .blocks
            .iter()
            .filter(|b| matches!(b, HBlock::LowRank { .. }))
            .count(),
        dense_blocks: h
            .blocks
            .iter()
            .filter(|b| matches!(b, HBlock::Dense { .. }))
            .count(),
        stored: h.stored(),
        compression: h.stored() as f64 / (nt as f64 * nt as f64),
        iterations,
    };
    Ok((
        HelmholtzSolution {
            k,
            kind,
            density,
            dir,
        },
        stats,
    ))
}

/// Solves the exterior scattering problem.
///
/// `coupling` is the Burton–Miller parameter `α` (`None` for the
/// recommended `i / k`, `Some(0)` for the conventional equation, which is
/// singular at the irregular frequencies).
///
/// # Errors
/// For an invalid mesh, `k <= 0` or a singular system.
pub fn solve_helmholtz_exterior(
    mesh: &SurfaceMesh3D,
    k: f64,
    kind: Scatterer,
    wave: PlaneWave,
    coupling: Option<Complex64>,
) -> Result<HelmholtzSolution, String> {
    mesh.validate()?;
    if !(k > 0.0 && k.is_finite()) {
        return Err("the wave number must be positive".into());
    }
    let dn = norm(wave.direction);
    if dn == 0.0 {
        return Err("the incident direction must be non-zero".into());
    }
    let dir = [
        wave.direction[0] / dn,
        wave.direction[1] / dn,
        wave.direction[2] / dn,
    ];
    let alpha = coupling.unwrap_or_else(|| Complex64::new(0.0, 1.0 / k));
    let nt = mesh.triangles.len();
    let tris: Vec<[V3; 3]> = (0..nt).filter_map(|t| mesh.corners(t)).collect();
    let normals: Vec<V3> = (0..nt).filter_map(|t| mesh.normal(t)).collect();
    let cents: Vec<V3> = (0..nt).filter_map(|t| mesh.centroid(t)).collect();
    if tris.len() != nt || normals.len() != nt || cents.len() != nt {
        return Err("degenerate triangle in the mesh".into());
    }
    let setup = Setup {
        k,
        kind,
        alpha,
        tris,
        normals,
        cents,
    };
    let rows: Vec<Vec<Complex64>> = (0..nt)
        .into_par_iter()
        .map(|i| (0..nt).map(|j| setup.entry(i, j)).collect())
        .collect();
    let a: Vec<Complex64> = rows.into_iter().flatten().collect();
    let rhs = setup.rhs(dir);
    let density = solve_dense(a, rhs, nt)?;
    Ok(HelmholtzSolution {
        k,
        kind,
        density,
        dir,
    })
}

impl HelmholtzSolution {
    /// Scattered pressure at a point away from the surface.
    #[must_use]
    pub fn scattered_field(
        &self,
        mesh: &SurfaceMesh3D,
        x: V3,
    ) -> Complex64 {
        let mut s = Complex64::new(0.0, 0.0);
        for (j, &dens) in self.density.iter().enumerate() {
            let (Some(tri), Some(ny)) = (mesh.corners(j), mesh.normal(j)) else {
                continue;
            };
            let nx = [0.0; 3];
            match self.kind {
                | Scatterer::SoundHard => {
                    s += integrate(Kernel::Double, self.k, x, nx, tri, ny, 0) * dens;
                },
                | Scatterer::SoundSoft => {
                    s -= integrate(Kernel::Single, self.k, x, nx, tri, ny, 0) * dens;
                },
            }
        }
        s
    }

    /// Incident field at a point.
    #[must_use]
    pub fn incident_field(
        &self,
        x: V3,
    ) -> Complex64 {
        plane_wave(self.k, self.dir, x)
    }
}
