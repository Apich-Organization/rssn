//! # Ewald summation and smooth particle-mesh Ewald (PME)
//!
//! Long-range Coulomb interactions in a periodic 3-D box (orthorhombic or general triclinic cell).
//! The total energy splits into a screened real-space part (use
//! [`PairTerm::Coulomb`] with `alpha > 0` in the MD force field), a
//! reciprocal-space part and the self/neutralising terms provided here.
//!
//! | Method | Cost | Notes |
//! |---|---|---|
//! | [`ewald_reciprocal`] | `O(N K^3)` | direct reciprocal sum, energy, forces and virial |
//! | [`pme_reciprocal`] | `O(N p^3 + K^3 log K)` | cardinal B-spline (even order) charge spreading, 3-D FFT from [`crate::kernels::fftconv`], analytic spline-derivative forces |
//! | [`ewald_energy_forces`] | `O(N^2 + N K^3)` | full reference (real space with minimum image, needs `r_cut <= L/2`) |
//! | [`ewald_reciprocal_cell`], [`pme_reciprocal_cell`], [`ewald_energy_forces_cell`] | as above | **triclinic cells** (rows of `h` are the lattice vectors): reciprocal basis, fractional-coordinate B-splines with the `h^-1` chain rule, real-space sum over all images within `r_cut` |
//!
//! Dimensions: 3-D only (the Ewald kernel `1/r` is the 3-D Coulomb law).
//! Conventions: Gaussian-like units with coupling constant `k_c`; the
//! neutralising background correction is included for charged cells.

use std::f64::consts::PI;

use crate::kernels::fftconv::Complex;
use crate::kernels::fftconv::fft;
use crate::kernels::special::erfc_numerical;
use crate::sim::md_engine::ExtraForce;
use crate::sim::md_engine::PairTerm;

/// Result of a reciprocal-space evaluation.
#[derive(Clone, Debug)]
pub struct Recip {
    /// Energy (reciprocal + self + background).
    pub energy: f64,
    /// Forces.
    pub forces: Vec<[f64; 3]>,
    /// Virial tensor `-dE/d(strain)` of the reciprocal part.
    pub virial: [[f64; 3]; 3],
}

fn charge_terms(
    q: &[f64],
    vol: f64,
    alpha: f64,
    kc: f64,
) -> f64 {
    let q2: f64 = q.iter().map(|x| x * x).sum();
    let qt: f64 = q.iter().sum();
    -kc * alpha / PI.sqrt() * q2 - kc * PI * qt * qt / (2.0 * vol * alpha * alpha)
}

/// Diagonal cell matrix of an orthorhombic box.
#[must_use]
pub fn diag_cell(len: [f64; 3]) -> [[f64; 3]; 3] {
    [[len[0], 0.0, 0.0], [0.0, len[1], 0.0], [0.0, 0.0, len[2]]]
}

/// Determinant of a 3x3 matrix.
fn det3(h: &[[f64; 3]; 3]) -> f64 {
    h[0][0] * (h[1][1] * h[2][2] - h[1][2] * h[2][1]) - h[0][1] * (h[1][0] * h[2][2] - h[1][2] * h[2][0])
        + h[0][2] * (h[1][0] * h[2][1] - h[1][1] * h[2][0])
}

/// Inverse of a 3x3 matrix (`None` if singular).
fn inv3(h: &[[f64; 3]; 3]) -> Option<[[f64; 3]; 3]> {
    let d = det3(h);
    if d.abs() < 1e-300 {
        return None;
    }
    let c = |a: usize, b: usize| {
        let (a1, a2) = ((a + 1) % 3, (a + 2) % 3);
        let (b1, b2) = ((b + 1) % 3, (b + 2) % 3);
        h[a1][b1] * h[a2][b2] - h[a1][b2] * h[a2][b1]
    };
    let mut inv = [[0.0; 3]; 3];
    for a in 0..3 {
        for b in 0..3 {
            inv[b][a] = c(a, b) / d;
        }
    }
    Some(inv)
}

/// Reciprocal lattice vectors (rows `b_i` with `b_i . a_j = 2 pi delta_ij`) of the cell `h`
/// whose rows are the lattice vectors.
///
/// # Errors
/// Fails for a singular cell.
pub fn reciprocal_basis(h: &[[f64; 3]; 3]) -> Result<[[f64; 3]; 3], String> {
    let inv = inv3(h).ok_or("singular cell matrix")?;
    let mut b = [[0.0; 3]; 3];
    for i in 0..3 {
        for c in 0..3 {
            b[i][c] = 2.0 * PI * inv[c][i];
        }
    }
    Ok(b)
}

/// Direct reciprocal-space Ewald sum of an orthorhombic box with all wave vectors `|n_a| <= kmax`,
/// plus self and background terms. See [`ewald_reciprocal_cell`] for general cells.
#[must_use]
pub fn ewald_reciprocal(
    pos: &[[f64; 3]],
    q: &[f64],
    len: [f64; 3],
    alpha: f64,
    kmax: i32,
    kc: f64,
) -> Recip {
    ewald_reciprocal_cell(pos, q, &diag_cell(len), alpha, kmax, kc).unwrap_or(Recip {
        energy: f64::NAN,
        forces: vec![[f64::NAN; 3]; pos.len()],
        virial: [[f64::NAN; 3]; 3],
    })
}

/// Direct reciprocal-space Ewald sum for a **general triclinic cell** `h` (rows are the lattice
/// vectors): wave vectors `k = n_1 b_1 + n_2 b_2 + n_3 b_3` with `|n_i| <= kmax`, energy, forces,
/// the reciprocal virial `-dE/d(strain)` and the self and neutralising-background terms.
///
/// # Errors
/// Fails for a singular cell.
pub fn ewald_reciprocal_cell(
    pos: &[[f64; 3]],
    q: &[f64],
    h: &[[f64; 3]; 3],
    alpha: f64,
    kmax: i32,
    kc: f64,
) -> Result<Recip, String> {
    let vol = det3(h).abs();
    let b = reciprocal_basis(h)?;
    let n = pos.len();
    let mut energy = charge_terms(q, vol, alpha, kc);
    let mut forces = vec![[0.0; 3]; n];
    let mut virial = [[0.0; 3]; 3];
    let pref = 2.0 * PI * kc / vol;
    for nx in -kmax..=kmax {
        for ny in -kmax..=kmax {
            for nz in -kmax..=kmax {
                if nx == 0 && ny == 0 && nz == 0 {
                    continue;
                }
                let (fx, fy, fz) = (f64::from(nx), f64::from(ny), f64::from(nz));
                let k = [
                    fx * b[0][0] + fy * b[1][0] + fz * b[2][0],
                    fx * b[0][1] + fy * b[1][1] + fz * b[2][1],
                    fx * b[0][2] + fy * b[1][2] + fz * b[2][2],
                ];
                let k2 = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
                let a = (-k2 / (4.0 * alpha * alpha)).exp() / k2;
                let (mut sc, mut ss) = (0.0, 0.0);
                let mut cs = Vec::with_capacity(n);
                for i in 0..n {
                    let ph = k[0] * pos[i][0] + k[1] * pos[i][1] + k[2] * pos[i][2];
                    let (s, c) = ph.sin_cos();
                    sc += q[i] * c;
                    ss += q[i] * s;
                    cs.push((c, s));
                }
                let ek = pref * a * (sc * sc + ss * ss);
                energy += ek;
                let w = -2.0 * (1.0 / k2 + 1.0 / (4.0 * alpha * alpha));
                for x in 0..3 {
                    for y in 0..3 {
                        virial[x][y] += ek * (if x == y { 1.0 } else { 0.0 } + w * k[x] * k[y]);
                    }
                }
                for i in 0..n {
                    let (c, s) = cs[i];
                    let im = sc * s - ss * c;
                    let f = 2.0 * pref * a * q[i] * im;
                    for x in 0..3 {
                        forces[i][x] += f * k[x];
                    }
                }
            }
        }
    }
    Ok(Recip { energy, forces, virial })
}

/// Complete Ewald energy and forces for a **triclinic cell**: the real-space sum runs over all
/// periodic images within `rcut` (no minimum-image restriction, self-images included), the
/// reciprocal part is [`ewald_reciprocal_cell`].
///
/// # Errors
/// Fails for a singular cell.
pub fn ewald_energy_forces_cell(
    pos: &[[f64; 3]],
    q: &[f64],
    h: &[[f64; 3]; 3],
    alpha: f64,
    kmax: i32,
    rcut: f64,
    kc: f64,
) -> Result<(f64, Vec<[f64; 3]>), String> {
    let rec = ewald_reciprocal_cell(pos, q, h, alpha, kmax, kc)?;
    let b = reciprocal_basis(h)?;
    let mut e = rec.energy;
    let mut f = rec.forces;
    // number of images needed along each lattice direction
    let nmax: Vec<i32> = (0..3)
        .map(|i| {
            let bn = (b[i][0] * b[i][0] + b[i][1] * b[i][1] + b[i][2] * b[i][2]).sqrt();
            (rcut * bn / (2.0 * PI)).ceil() as i32 + 1
        })
        .collect();
    let term = PairTerm::Coulomb { k: kc, alpha };
    let n = pos.len();
    for i in 0..n {
        for j in i..n {
            for nx in -nmax[0]..=nmax[0] {
                for ny in -nmax[1]..=nmax[1] {
                    for nz in -nmax[2]..=nmax[2] {
                        if i == j && (nx, ny, nz) <= (0, 0, 0) {
                            continue;
                        }
                        let (fx, fy, fz) = (f64::from(nx), f64::from(ny), f64::from(nz));
                        let mut d = [0.0; 3];
                        for a in 0..3 {
                            d[a] = pos[i][a] - pos[j][a] + fx * h[0][a] + fy * h[1][a] + fz * h[2][a];
                        }
                        let r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
                        if r2 < rcut * rcut && r2 > 0.0 {
                            let (u, fr) = term.eval(r2, q[i] * q[j]);
                            e += u;
                            for a in 0..3 {
                                f[i][a] += fr * d[a];
                                f[j][a] -= fr * d[a];
                            }
                        }
                    }
                }
            }
        }
    }
    Ok((e, f))
}

/// Complete Ewald energy and forces (real space with minimum image requires `rcut <= min(L)/2`).
#[must_use]
pub fn ewald_energy_forces(
    pos: &[[f64; 3]],
    q: &[f64],
    len: [f64; 3],
    alpha: f64,
    kmax: i32,
    rcut: f64,
    kc: f64,
) -> (f64, Vec<[f64; 3]>) {
    let rec = ewald_reciprocal(pos, q, len, alpha, kmax, kc);
    let mut e = rec.energy;
    let mut f = rec.forces;
    let term = PairTerm::Coulomb { k: kc, alpha };
    for i in 0..pos.len() {
        for j in i + 1..pos.len() {
            let mut d = [0.0; 3];
            for a in 0..3 {
                d[a] = pos[i][a] - pos[j][a];
                d[a] -= len[a] * (d[a] / len[a]).round();
            }
            let r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
            if r2 < rcut * rcut {
                let (u, fr) = term.eval(r2, q[i] * q[j]);
                e += u;
                for a in 0..3 {
                    f[i][a] += fr * d[a];
                    f[j][a] -= fr * d[a];
                }
            }
        }
    }
    (e, f)
}

/// Cardinal B-spline `M_n(u)` (support `(0, n)`).
#[must_use]
pub fn bspline(
    n: usize,
    u: f64,
) -> f64 {
    if u <= 0.0 || u >= n as f64 {
        return 0.0;
    }
    if n == 2 {
        return 1.0 - (u - 1.0).abs();
    }
    let nm = n as f64;
    u / (nm - 1.0) * bspline(n - 1, u) + (nm - u) / (nm - 1.0) * bspline(n - 1, nm - u)
}

fn bspline_deriv(
    n: usize,
    u: f64,
) -> f64 {
    bspline(n - 1, u) - bspline(n - 1, u - 1.0)
}

fn fft3(
    data: &mut [Complex],
    dims: [usize; 3],
    inverse: bool,
) {
    let [n0, n1, n2] = dims;
    let idx = |i: usize, j: usize, k: usize| (i * n1 + j) * n2 + k;
    for i in 0..n0 {
        for j in 0..n1 {
            let line: Vec<Complex> = (0..n2).map(|k| data[idx(i, j, k)]).collect();
            for (k, v) in fft(&line, inverse).into_iter().enumerate() {
                data[idx(i, j, k)] = v;
            }
        }
    }
    for i in 0..n0 {
        for k in 0..n2 {
            let line: Vec<Complex> = (0..n1).map(|j| data[idx(i, j, k)]).collect();
            for (j, v) in fft(&line, inverse).into_iter().enumerate() {
                data[idx(i, j, k)] = v;
            }
        }
    }
    for j in 0..n1 {
        for k in 0..n2 {
            let line: Vec<Complex> = (0..n0).map(|i| data[idx(i, j, k)]).collect();
            for (i, v) in fft(&line, inverse).into_iter().enumerate() {
                data[idx(i, j, k)] = v;
            }
        }
    }
}

/// `|b(m)|^2` of the PME structure-factor correction for grid size `k` and spline order `n`.
fn b2(
    n: usize,
    k: usize,
    m: usize,
) -> f64 {
    let (mut re, mut im) = (0.0, 0.0);
    for j in 0..n - 1 {
        let ang = 2.0 * PI * (m * j) as f64 / k as f64;
        let w = bspline(n, (j + 1) as f64);
        re += w * ang.cos();
        im += w * ang.sin();
    }
    1.0 / (re * re + im * im).max(1e-300)
}

/// B-spline weights and derivatives of one particle.
struct Spl {
    base: [i64; 3],
    w: [Vec<f64>; 3],
    dw: [Vec<f64>; 3],
}

/// Smooth PME reciprocal energy, forces and virial (plus self/background terms) of an
/// orthorhombic box. See [`pme_reciprocal_cell`] for general cells.
///
/// `grid` are the mesh sizes (any size, powers of two are fastest), `order` the (even) B-spline order.
#[must_use]
pub fn pme_reciprocal(
    pos: &[[f64; 3]],
    q: &[f64],
    len: [f64; 3],
    alpha: f64,
    grid: [usize; 3],
    order: usize,
    kc: f64,
) -> Recip {
    pme_reciprocal_cell(pos, q, &diag_cell(len), alpha, grid, order, kc).unwrap_or(Recip {
        energy: f64::NAN,
        forces: vec![[f64::NAN; 3]; pos.len()],
        virial: [[f64::NAN; 3]; 3],
    })
}

/// Smooth particle-mesh Ewald for a **general triclinic cell** `h` (rows are the lattice
/// vectors): charges are spread with cardinal B-splines on a mesh in the fractional coordinates
/// `s = r h^-1`, the convolution is done by FFT, and the forces use the analytic spline
/// derivatives (chain rule through `h^-1`). The structure-factor correction and the
/// reciprocal-space influence function are those of the sheared lattice (`k = m_1 b_1 + m_2 b_2 +
/// m_3 b_3`).
///
/// # Errors
/// Fails for a singular cell, an odd or too small spline order, or a mesh smaller than the order.
pub fn pme_reciprocal_cell(
    pos: &[[f64; 3]],
    q: &[f64],
    h: &[[f64; 3]; 3],
    alpha: f64,
    grid: [usize; 3],
    order: usize,
    kc: f64,
) -> Result<Recip, String> {
    if order < 2 || order % 2 != 0 || grid.iter().any(|&g| g < order) {
        return Err("PME needs an even spline order >= 2 and meshes at least as large".into());
    }
    let vol = det3(h).abs();
    let hinv = inv3(h).ok_or("singular cell matrix")?;
    let bvec = reciprocal_basis(h)?;
    let n = pos.len();
    let total = grid[0] * grid[1] * grid[2];
    let idx = |i: usize, j: usize, k: usize| (i * grid[1] + j) * grid[2] + k;
    let mut qg = vec![(0.0, 0.0); total];
    let mut spl = Vec::with_capacity(n);
    for i in 0..n {
        let mut base = [0i64; 3];
        let mut w: [Vec<f64>; 3] = [vec![], vec![], vec![]];
        let mut dw: [Vec<f64>; 3] = [vec![], vec![], vec![]];
        for a in 0..3 {
            let s = pos[i][0] * hinv[0][a] + pos[i][1] * hinv[1][a] + pos[i][2] * hinv[2][a];
            let u = (s - s.floor()) * grid[a] as f64;
            let fl = u.floor();
            let fr = u - fl;
            base[a] = fl as i64;
            for j in 0..order {
                w[a].push(bspline(order, fr + j as f64));
                // derivative with respect to the fractional coordinate s
                dw[a].push(bspline_deriv(order, fr + j as f64) * grid[a] as f64);
            }
        }
        for jx in 0..order {
            for jy in 0..order {
                for jz in 0..order {
                    let ix = (base[0] - jx as i64).rem_euclid(grid[0] as i64) as usize;
                    let iy = (base[1] - jy as i64).rem_euclid(grid[1] as i64) as usize;
                    let iz = (base[2] - jz as i64).rem_euclid(grid[2] as i64) as usize;
                    qg[idx(ix, iy, iz)].0 += q[i] * w[0][jx] * w[1][jy] * w[2][jz];
                }
            }
        }
        spl.push(Spl { base, w, dw });
    }
    fft3(&mut qg, grid, false);
    let b2x: Vec<f64> = (0..grid[0]).map(|m| b2(order, grid[0], m)).collect();
    let b2y: Vec<f64> = (0..grid[1]).map(|m| b2(order, grid[1], m)).collect();
    let b2z: Vec<f64> = (0..grid[2]).map(|m| b2(order, grid[2], m)).collect();
    let signed = |m: usize, k: usize| if m <= k / 2 { m as f64 } else { m as f64 - k as f64 };
    let pref = 2.0 * PI * kc / vol;
    let mut energy = charge_terms(q, vol, alpha, kc);
    let mut virial = [[0.0; 3]; 3];
    let mut conv = vec![(0.0, 0.0); total];
    for i in 0..grid[0] {
        for j in 0..grid[1] {
            for k in 0..grid[2] {
                if i == 0 && j == 0 && k == 0 {
                    continue;
                }
                let (mi, mj, mk) = (signed(i, grid[0]), signed(j, grid[1]), signed(k, grid[2]));
                let kv = [
                    mi * bvec[0][0] + mj * bvec[1][0] + mk * bvec[2][0],
                    mi * bvec[0][1] + mj * bvec[1][1] + mk * bvec[2][1],
                    mi * bvec[0][2] + mj * bvec[1][2] + mk * bvec[2][2],
                ];
                let k2 = kv[0] * kv[0] + kv[1] * kv[1] + kv[2] * kv[2];
                let c = (-k2 / (4.0 * alpha * alpha)).exp() / k2 * b2x[i] * b2y[j] * b2z[k];
                let f = qg[idx(i, j, k)];
                let mag = f.0 * f.0 + f.1 * f.1;
                let ek = pref * c * mag;
                energy += ek;
                let w = -2.0 * (1.0 / k2 + 1.0 / (4.0 * alpha * alpha));
                for x in 0..3 {
                    for y in 0..3 {
                        virial[x][y] += ek * (if x == y { 1.0 } else { 0.0 } + w * kv[x] * kv[y]);
                    }
                }
                conv[idx(i, j, k)] = (c * f.0, c * f.1);
            }
        }
    }
    fft3(&mut conv, grid, true);
    let ntot = total as f64;
    let mut forces = vec![[0.0; 3]; n];
    for i in 0..n {
        let s = &spl[i];
        for jx in 0..order {
            for jy in 0..order {
                for jz in 0..order {
                    let ix = (s.base[0] - jx as i64).rem_euclid(grid[0] as i64) as usize;
                    let iy = (s.base[1] - jy as i64).rem_euclid(grid[1] as i64) as usize;
                    let iz = (s.base[2] - jz as i64).rem_euclid(grid[2] as i64) as usize;
                    let de_dq = 2.0 * pref * conv[idx(ix, iy, iz)].0 * ntot;
                    // gradient with respect to the fractional coordinates, then chain rule via h^-1
                    let gs = [s.dw[0][jx] * s.w[1][jy] * s.w[2][jz], s.w[0][jx] * s.dw[1][jy] * s.w[2][jz], s.w[0][jx] * s.w[1][jy] * s.dw[2][jz]];
                    for c in 0..3 {
                        let gc = gs[0] * hinv[c][0] + gs[1] * hinv[c][1] + gs[2] * hinv[c][2];
                        forces[i][c] -= q[i] * de_dq * gc;
                    }
                }
            }
        }
    }
    Ok(Recip { energy, forces, virial })
}

/// Force callback adding the direct-Ewald reciprocal, self and background terms to an MD force evaluation.
#[must_use]
pub fn ewald_extra(
    alpha: f64,
    kmax: i32,
    kc: f64,
) -> ExtraForce<3> {
    Box::new(move |sys, out| {
        let h = sys.bx.cell.unwrap_or_else(|| diag_cell(sys.bx.len));
        let r = ewald_reciprocal_cell(&sys.pos, &sys.charge, &h, alpha, kmax, kc).unwrap_or(Recip {
            energy: f64::NAN,
            forces: vec![[f64::NAN; 3]; sys.pos.len()],
            virial: [[f64::NAN; 3]; 3],
        });
        out.epot += r.energy;
        for (f, g) in out.f.iter_mut().zip(&r.forces) {
            for a in 0..3 {
                f[a] += g[a];
            }
        }
        for a in 0..3 {
            for b in 0..3 {
                out.virial[a][b] += r.virial[a][b];
            }
        }
    })
}

/// Force callback adding the PME reciprocal, self and background terms.
#[must_use]
pub fn pme_extra(
    alpha: f64,
    grid: [usize; 3],
    order: usize,
    kc: f64,
) -> ExtraForce<3> {
    Box::new(move |sys, out| {
        let h = sys.bx.cell.unwrap_or_else(|| diag_cell(sys.bx.len));
        let r = pme_reciprocal_cell(&sys.pos, &sys.charge, &h, alpha, grid, order, kc).unwrap_or(Recip {
            energy: f64::NAN,
            forces: vec![[f64::NAN; 3]; sys.pos.len()],
            virial: [[f64::NAN; 3]; 3],
        });
        out.epot += r.energy;
        for (f, g) in out.f.iter_mut().zip(&r.forces) {
            for a in 0..3 {
                f[a] += g[a];
            }
        }
        for a in 0..3 {
            for b in 0..3 {
                out.virial[a][b] += r.virial[a][b];
            }
        }
    })
}

/// Real-space screened Coulomb energy of a pair (convenience for checks).
#[must_use]
pub fn real_space_pair(
    r: f64,
    qq: f64,
    alpha: f64,
    kc: f64,
) -> f64 {
    kc * qq * erfc_numerical(alpha * r) / r
}
