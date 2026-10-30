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
//! Limitations: constant elements only, dense `O(n²)` storage, closed
//! orientable surfaces with outward normals, no impedance boundary
//! conditions.

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
    let rows: Vec<Vec<Complex64>> = (0..nt)
        .into_par_iter()
        .map(|i| {
            let (x, nx) = (cents[i], normals[i]);
            (0..nt)
                .map(|j| {
                    let own = i == j;
                    let (s_ij, n_ij) = if own {
                        self_terms(k, tris[j])
                    } else {
                        (Complex64::new(0.0, 0.0), Complex64::new(0.0, 0.0))
                    };
                    match kind {
                        | Scatterer::SoundHard => {
                            let kk = if own {
                                Complex64::new(0.0, 0.0)
                            } else {
                                integrate(Kernel::Double, k, x, nx, tris[j], normals[j], 0)
                            };
                            let nn = if own {
                                n_ij
                            } else if alpha.norm() > 0.0 {
                                integrate(Kernel::Hyper, k, x, nx, tris[j], normals[j], 0)
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
                                integrate(Kernel::Single, k, x, nx, tris[j], normals[j], 0)
                            };
                            let kp = if own {
                                Complex64::new(0.0, 0.0)
                            } else {
                                integrate(Kernel::Adjoint, k, x, nx, tris[j], normals[j], 0)
                            };
                            let half = if own { 0.5 } else { 0.0 };
                            ss + alpha * (Complex64::new(half, 0.0) + kp)
                        },
                    }
                })
                .collect()
        })
        .collect();
    let a: Vec<Complex64> = rows.into_iter().flatten().collect();
    let rhs: Vec<Complex64> = (0..nt)
        .map(|i| {
            let pinc = plane_wave(k, dir, cents[i]);
            let dpinc = pinc * Complex64::new(0.0, k * dot(dir, normals[i]));
            match kind {
                | Scatterer::SoundHard => pinc - alpha * dpinc,
                | Scatterer::SoundSoft => pinc + alpha * dpinc,
            }
        })
        .collect();
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
                    s += integrate(Kernel::Double, self.k, x, nx, tri, ny, 0) * dens
                },
                | Scatterer::SoundSoft => {
                    s -= integrate(Kernel::Single, self.k, x, nx, tri, ny, 0) * dens
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
