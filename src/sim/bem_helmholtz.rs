//! Helmholtz BEM: exterior acoustic scattering with Burton-Miller.
//!
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

use crate::kernels::complex_scalar::Complex;
use crate::sim::physics_bem::SurfaceMesh3D;
use num_complex::Complex64;
use rayon::prelude::*;
use std::alloc::Layout;
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
    /// Tolerance for the actual left-preconditioned relative residual.
    pub gmres_tol: f64,
    /// Positive GMRES restart length.
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

impl HBlock {
    fn rows(&self) -> &[usize] {
        match self {
            | Self::Dense { rows, .. } | Self::LowRank { rows, .. } => rows,
        }
    }

    fn product(
        &self,
        x: &[Complex64],
        out: &mut [Complex64],
    ) {
        match self {
            | Self::Dense { cols, data, .. } => {
                for (result, row) in out.iter_mut().zip(data.chunks_exact(cols.len())) {
                    let mut acc = Complex64::ZERO;
                    for (a, &j) in row.iter().zip(cols) {
                        acc += a * x[j];
                    }
                    *result = acc;
                }
            },
            | Self::LowRank { cols, u, v, .. } => {
                out.fill(Complex64::ZERO);
                for (ul, vl) in u.iter().zip(v) {
                    let mut scale = Complex64::ZERO;
                    for (a, &j) in vl.iter().zip(cols) {
                        scale += a * x[j];
                    }
                    for (value, weight) in out.iter_mut().zip(ul) {
                        *value += scale * weight;
                    }
                }
            },
        }
    }
}

const H_PRODUCT_BLOCK_GRAIN: usize = 16;
const H_PRODUCT_ROW_GRAIN: usize = 2048;

enum HProductWorkspace {
    Serial(Vec<Complex64>),
    Parallel {
        parts: Vec<Vec<Complex64>>,
        contributions: Option<Vec<Vec<(usize, usize)>>>,
    },
}

impl HProductWorkspace {
    fn new(matrix: &HMat) -> Self {
        if !crate::execution::should_split(matrix.blocks.len(), H_PRODUCT_BLOCK_GRAIN) {
            let rows = matrix
                .blocks
                .iter()
                .map(|block| block.rows().len())
                .max()
                .unwrap_or(0);
            return Self::Serial(vec![Complex64::ZERO; rows]);
        }
        let contributions =
            crate::execution::should_split(matrix.n, H_PRODUCT_ROW_GRAIN).then(|| {
                let mut contributions = vec![Vec::new(); matrix.n];
                for (block, data) in matrix.blocks.iter().enumerate() {
                    for (local, &row) in data.rows().iter().enumerate() {
                        contributions[row].push((block, local));
                    }
                }
                contributions
            });
        let parts = matrix
            .blocks
            .iter()
            .map(|block| vec![Complex64::ZERO; block.rows().len()])
            .collect();
        Self::Parallel { parts, contributions }
    }
}

impl HMat {
    fn matvec(
        &self,
        x: &[Complex64],
        y: &mut [Complex64],
        scratch: &mut HProductWorkspace,
    ) {
        debug_assert_eq!(y.len(), self.n);
        match scratch {
            | HProductWorkspace::Serial(part) => {
                y.fill(Complex64::ZERO);
                for block in &self.blocks {
                    let rows = block.rows();
                    let part = &mut part[..rows.len()];
                    block.product(x, part);
                    for (&row, &value) in rows.iter().zip(part.iter()) {
                        y[row] += value;
                    }
                }
            },
            | HProductWorkspace::Parallel { parts, contributions } => {
                let blocks = &self.blocks;
                crate::kernels::par::for_each_mut(parts, H_PRODUCT_BLOCK_GRAIN, |index, part| {
                    blocks[index].product(x, part);
                });
                if let Some(contributions) = contributions {
                    // Each row consumes contributions in the original block order.
                    crate::kernels::par::for_each_mut(y, H_PRODUCT_ROW_GRAIN, |row, value| {
                        *value = contributions[row]
                            .iter()
                            .fold(Complex64::ZERO, |sum, &(block, local)| {
                                sum + parts[block][local]
                            });
                    });
                } else {
                    y.fill(Complex64::ZERO);
                    for (block, part) in blocks.iter().zip(parts) {
                        for (&row, &value) in block.rows().iter().zip(part.iter()) {
                            y[row] += value;
                        }
                    }
                }
            },
        }
    }

    fn stored(&self) -> usize {
        self.blocks
            .iter()
            .map(|b| {
                match b {
                    | HBlock::Dense { data, .. } => data.len(),
                    | HBlock::LowRank { rows, cols, u, .. } => u.len() * (rows.len() + cols.len()),
                }
            })
            .sum()
    }
}

const fn gmres_finite(z: Complex64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}

fn gmres_div(
    numerator: Complex64,
    denominator: Complex64,
) -> Complex64 {
    let q =
        Complex::new(numerator.re, numerator.im).div(Complex::new(denominator.re, denominator.im));
    Complex64::new(q.re, q.im)
}

fn gmres_peak(values: &[Complex64]) -> f64 {
    values
        .iter()
        .fold(0.0_f64, |peak, z| peak.max(z.re.abs()).max(z.im.abs()))
}

fn gmres_norm(values: &[Complex64]) -> f64 {
    let ordinary = vec_norm(values);
    if ordinary.is_nan() || (1e-150..=1e150).contains(&ordinary) {
        return ordinary;
    }
    let peak = gmres_peak(values);
    if peak == 0.0 || !peak.is_finite() {
        return peak;
    }
    peak * values
        .iter()
        .map(|z| (z / peak).norm_sqr())
        .sum::<f64>()
        .sqrt()
}

const fn gmres_controls_valid(
    tol: f64,
    restart: usize,
) -> bool {
    restart > 0 && tol.is_finite() && tol >= 0.0
}

/// Restarted complex GMRES using the actual left-preconditioned residual.
fn gmres_complex(
    apply: &mut dyn FnMut(&[Complex64], &mut [Complex64]),
    diag: &[Complex64],
    b: &[Complex64],
    tol: f64,
    restart: usize,
    max_iter: usize,
) -> Option<(Vec<Complex64>, usize)> {
    let n = b.len();
    if !gmres_controls_valid(tol, restart)
        || diag.len() != n
        || b.iter().chain(diag).any(|&z| !gmres_finite(z))
    {
        return None;
    }
    let prec = |v: &mut [Complex64]| {
        for (x, &d) in v.iter_mut().zip(diag) {
            let value = if d == Complex64::ZERO {
                *x
            } else {
                gmres_div(*x, d)
            };
            if !gmres_finite(value) || (value == Complex64::ZERO && *x != Complex64::ZERO) {
                return false;
            }
            *x = value;
        }
        true
    };
    let mut pb = b.to_vec();
    if !prec(&mut pb) {
        return None;
    }
    let mut x = vec![Complex64::ZERO; n];
    let estimated_norm = vec_norm(&pb);
    let scale = if (1e-100..=1e100).contains(&estimated_norm) {
        None
    } else {
        let peak = gmres_peak(&pb);
        if peak == 0.0 {
            return b.iter().all(|&z| z == Complex64::ZERO).then_some((x, 0));
        }
        for value in &mut pb {
            *value /= peak;
        }
        Some(peak)
    };
    let bnorm = if scale.is_some() {
        gmres_norm(&pb)
    } else {
        estimated_norm
    };
    let mut apply_into = |input: &[Complex64], output: &mut [Complex64]| {
        output.fill(Complex64::new(f64::NAN, f64::NAN));
        apply(input, output);
    };
    let residual = |apply_into: &mut dyn FnMut(&[Complex64], &mut [Complex64]),
                    input: &[Complex64],
                    output: &mut [Complex64]| {
        apply_into(input, output);
        for (value, &rhs) in output.iter_mut().zip(b) {
            *value = rhs - *value;
        }
        prec(output)
    };
    let m = restart.min(n).min(max_iter);
    let stride = m.checked_add(1)?;
    let small_len = m.checked_mul(stride)?;
    Layout::array::<Complex64>(small_len).ok()?;
    let mut basis = vec![vec![Complex64::ZERO; n]];
    let mut w = vec![Complex64::ZERO; n];
    let mut h = vec![Complex64::ZERO; small_len];
    let mut cs = vec![0.0; m];
    let mut sn = vec![Complex64::ZERO; m];
    let mut g = vec![Complex64::ZERO; stride];
    let mut total = 0;
    loop {
        if scale.is_some() {
            // Iterate on the scaled preconditioned equation. Applying A to
            // the scaled unknown avoids losing a uniformly tiny right side.
            apply_into(&x, &mut w);
            if !prec(&mut w) {
                return None;
            }
            for (value, &rhs) in w.iter_mut().zip(&pb) {
                *value = rhs - *value;
            }
        } else if !residual(&mut apply_into, &x, &mut w) {
            return None;
        }
        let beta = gmres_norm(&w);
        if !beta.is_finite() {
            return None;
        }
        let exhausted = total >= max_iter;
        let threshold = if exhausted {
            10.0 * tol * bnorm
        } else {
            tol * bnorm
        };
        if beta <= threshold {
            if let Some(scale) = scale {
                for value in &mut x {
                    *value *= scale;
                }
                if !x.iter().copied().all(gmres_finite) || !residual(&mut apply_into, &x, &mut w) {
                    return None;
                }
                for value in &mut w {
                    *value /= scale;
                }
                let restored = gmres_norm(&w);
                if !restored.is_finite() || restored > threshold {
                    return None;
                }
            }
            return Some((x, total));
        }
        if exhausted {
            return None;
        }
        for (value, direction) in w.iter().zip(&mut basis[0]) {
            *direction = value / beta;
        }
        g[0] = Complex64::new(beta, 0.0);
        let mut used = 0;
        for j in 0..m {
            apply_into(&basis[j], &mut w);
            total += 1;
            if !prec(&mut w) {
                return None;
            }
            let col = &mut h[j * stride..(j + 1) * stride];
            for (i, q) in basis[..=j].iter().enumerate() {
                let hij = herm_dot(q, &w);
                col[i] = hij;
                for (wk, qk) in w.iter_mut().zip(q) {
                    *wk -= hij * qk;
                }
            }
            let hn = gmres_norm(&w);
            col[j + 1] = Complex64::new(hn, 0.0);
            // apply the previous rotations to the new column
            for i in 0..j {
                let t = cs[i] * col[i] + sn[i] * col[i + 1];
                col[i + 1] = -sn[i].conj() * col[i] + cs[i] * col[i + 1];
                col[i] = t;
            }
            let (a, bb) = (col[j], col[j + 1]);
            let a_norm = a.norm();
            let square = a.norm_sqr() + bb.norm_sqr();
            let denom = if square.is_normal() {
                square.sqrt()
            } else {
                a_norm.hypot(bb.norm())
            };
            if !denom.is_finite() {
                return None;
            }
            if denom == 0.0 {
                break;
            }
            cs[j] = a_norm / denom;
            sn[j] = if a_norm > 0.0 {
                (a / a_norm) * bb.conj() / denom
            } else {
                bb.conj() / denom
            };
            col[j] = cs[j] * a + sn[j] * bb;
            col[j + 1] = Complex64::ZERO;
            let gj = g[j];
            g[j + 1] = -sn[j].conj() * gj;
            g[j] = cs[j] * gj;
            used = j + 1;
            if g[j + 1].norm() <= tol * bnorm || hn == 0.0 || total >= max_iter {
                break;
            }
            if used < m {
                if used == basis.len() {
                    basis.push(vec![Complex64::ZERO; n]);
                }
                for (value, direction) in w.iter().zip(&mut basis[used]) {
                    *direction = value / hn;
                }
            }
        }
        if used == 0 {
            return None;
        }
        // Back substitute in the Givens right-hand side, reusing its slots.
        for i in (0..used).rev() {
            let mut sum = g[i];
            for k in i + 1..used {
                sum -= h[k * stride + i] * g[k];
            }
            g[i] = gmres_div(sum, h[i * stride + i]);
            if !gmres_finite(g[i]) {
                return None;
            }
        }
        for (i, yi) in g[..used].iter().enumerate() {
            for (xk, qk) in x.iter_mut().zip(&basis[i]) {
                *xk += yi * qk;
            }
        }
        if !x.iter().copied().all(gmres_finite) {
            return None;
        }
    }
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
    if !gmres_controls_valid(opts.gmres_tol, opts.restart) {
        return Err(
            "GMRES requires a positive restart length and finite nonnegative tolerance".into(),
        );
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
            | Some((u, v)) => {
                HBlock::LowRank {
                    rows: rows.to_vec(),
                    cols: cols.to_vec(),
                    u,
                    v,
                }
            },
            | None => dense(),
        }
    });
    let h = HMat { n: nt, blocks };
    let diag: Vec<Complex64> = (0..nt).map(|i| entry(i, i)).collect();
    let rhs = setup.rhs(dir);
    let mut scratch = HProductWorkspace::new(&h);
    let (density, iterations) = gmres_complex(
        &mut |v, out| h.matvec(v, out, &mut scratch),
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

#[cfg(test)]
mod gmres_tests {
    use super::gmres_complex;
    use num_complex::Complex64 as C;

    #[test]
    fn hierarchical_product_reuses_buffers_and_preserves_block_order() {
        use super::H_PRODUCT_ROW_GRAIN;
        use super::HBlock;
        use super::HMat;
        use super::HProductWorkspace;

        for n in [3, H_PRODUCT_ROW_GRAIN] {
            let blocks = (0..18)
                .map(|index| {
                    let coefficient = match index % 3 {
                        | 0 => (1e16, -1e16),
                        | 1 => (1.0, 1.0),
                        | _ => (-1e16, 1e16),
                    };
                    let value = C::new(coefficient.0, coefficient.1);
                    if index % 2 == 0 {
                        HBlock::Dense {
                            rows: vec![n - 1, 0],
                            cols: vec![0],
                            data: vec![value, value],
                        }
                    } else {
                        HBlock::LowRank {
                            rows: vec![0, n - 1],
                            cols: vec![0],
                            u: vec![vec![C::ONE, C::ONE]],
                            v: vec![vec![value]],
                        }
                    }
                })
                .collect();
            let matrix = HMat { n, blocks };
            for workers in [1, 2] {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(workers)
                    .build()
                    .unwrap();
                pool.install(|| {
                    let mut scratch = HProductWorkspace::new(&matrix);
                    let mut serial = HProductWorkspace::Serial(vec![C::ZERO; 2]);
                    let mut actual = vec![C::I; n];
                    let mut expected = vec![C::ZERO; n];
                    for input in [C::ONE, C::ZERO, C::I, C::new(2.0, -3.0)] {
                        matrix.matvec(&[input], &mut expected, &mut serial);
                        matrix.matvec(&[input], &mut actual, &mut scratch);
                        for (a, b) in actual.iter().zip(&expected) {
                            assert_eq!(
                                (a.re.to_bits(), a.im.to_bits()),
                                (b.re.to_bits(), b.im.to_bits())
                            );
                        }
                    }
                });
            }
        }
    }

    fn identity(
        x: &[C],
        y: &mut [C],
    ) {
        y.copy_from_slice(x);
    }

    #[test]
    fn complex_gmres_recovers_both_extreme_equation_scales() {
        let expected = [C::new(1.0, 0.5), C::new(-1.0, -0.25)];
        for scale in [2.0_f64.powi(-600), 2.0_f64.powi(600)] {
            let rhs = expected.map(|z| z * scale);
            let (x, steps) =
                gmres_complex(&mut identity, &[C::ONE; 2], &rhs, 1e-12, 2, 12).unwrap();
            assert!(steps > 0);
            for (&actual, want) in x.iter().zip(expected) {
                assert!((actual / scale - want).norm() < 1e-14);
            }

            let phase = C::new(0.75, 0.5);
            let diagonal = phase * scale;
            let mut apply = |x: &[C], y: &mut [C]| {
                for (v, &z) in y.iter_mut().zip(x) {
                    *v = diagonal * z;
                }
            };
            let rhs = expected.map(|z| diagonal * z);
            let (x, _) = gmres_complex(&mut apply, &[diagonal; 2], &rhs, 1e-12, 2, 12).unwrap();
            for (&actual, want) in x.iter().zip(expected) {
                assert!((actual - want).norm() < 1e-14);
            }

            let rhs = expected.map(|z| phase * z);
            let (x, _) = gmres_complex(&mut apply, &[C::ZERO; 2], &rhs, 1e-12, 2, 12).unwrap();
            for (&actual, want) in x.iter().zip(expected) {
                assert!((actual * scale - want).norm() < 1e-14);
            }
        }
    }

    #[test]
    fn complex_gmres_checks_controls_and_the_restored_iterate() {
        assert!(gmres_complex(&mut identity, &[C::ONE], &[C::ONE], 1e-12, 0, 3).is_none());
        assert!(gmres_complex(&mut identity, &[], &[C::ONE], 1e-12, 1, 3).is_none());
        assert!(gmres_complex(&mut identity, &[C::ONE], &[C::ONE], f64::NAN, 1, 3).is_none());
        assert!(gmres_complex(&mut identity, &[C::ONE], &[C::ONE], 1e-12, 1, 0).is_none());
        assert!(
            gmres_complex(
                &mut |_, y| y.fill(C::ZERO),
                &[C::ONE],
                &[C::ONE],
                1e-12,
                1,
                3
            )
            .is_none()
        );
        assert!(gmres_complex(&mut |_, _| (), &[C::ONE], &[C::ONE], 1e-12, 1, 3).is_none());

        let tiny = C::new(f64::from_bits(1), 0.0);
        assert!(
            gmres_complex(
                &mut |x, y| y[0] = 2.0 * x[0],
                &[C::ONE],
                &[tiny],
                1e-12,
                1,
                3,
            )
            .is_none()
        );
        let (x, _) = gmres_complex(&mut identity, &[C::ONE], &[tiny], 1e-12, 1, 3).unwrap();
        assert_eq!(x, [tiny]);
    }

    #[test]
    fn reused_restart_storage_solves_a_nonnormal_complex_system() {
        let diagonal = [C::new(2.0, 0.5), C::new(1.5, -0.25), C::new(3.0, 0.25)];
        let mut apply = |x: &[C], y: &mut [C]| {
            y.copy_from_slice(&[
                diagonal[0] * x[0] + C::new(0.25, -0.125) * x[1],
                diagonal[1] * x[1] + C::new(-0.375, 0.25) * x[2],
                diagonal[2] * x[2] + C::new(0.25, 0.125) * x[0],
            ]);
        };
        let expected = [C::new(0.75, 0.5), C::new(-1.0, -0.25), C::new(0.25, 1.5)];
        let mut rhs = [C::ZERO; 3];
        apply(&expected, &mut rhs);
        for restart in [1, 2, 3] {
            let (x, steps) =
                gmres_complex(&mut apply, &diagonal, &rhs, 1e-12, restart, 100).unwrap();
            assert!(steps > 1);
            for (&actual, want) in x.iter().zip(expected) {
                assert!((actual - want).norm() < 5e-12);
            }
        }
    }
}
