//! # Ewald summation and smooth particle-mesh Ewald (PME)
//!
//! Long-range Coulomb interactions in a periodic orthorhombic 3-D box.
//! The total energy splits into a screened real-space part (use
//! [`PairTerm::Coulomb`] with `alpha > 0` in the MD force field), a
//! reciprocal-space part and the self/neutralising terms provided here.
//!
//! | Method | Cost | Notes |
//! |---|---|---|
//! | [`ewald_reciprocal`] | `O(N K^3)` | direct reciprocal sum, energy, forces and virial |
//! | [`pme_reciprocal`] | `O(N p^3 + K^3 log K)` | cardinal B-spline (even order) charge spreading, 3-D FFT from [`crate::kernels::fftconv`], analytic spline-derivative forces |
//! | [`ewald_energy_forces`] | `O(N^2 + N K^3)` | full reference (real space with minimum image, needs `r_cut <= L/2`) |
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

/// Direct reciprocal-space Ewald sum with all wave vectors `|n_a| <= kmax`, plus self and background terms.
#[must_use]
pub fn ewald_reciprocal(
    pos: &[[f64; 3]],
    q: &[f64],
    len: [f64; 3],
    alpha: f64,
    kmax: i32,
    kc: f64,
) -> Recip {
    let vol = len[0] * len[1] * len[2];
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
                let k = [2.0 * PI * f64::from(nx) / len[0], 2.0 * PI * f64::from(ny) / len[1], 2.0 * PI * f64::from(nz) / len[2]];
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
    Recip { energy, forces, virial }
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

/// Smooth PME reciprocal energy, forces and virial (plus self/background terms).
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
    let vol = len[0] * len[1] * len[2];
    let n = pos.len();
    let total = grid[0] * grid[1] * grid[2];
    let idx = |i: usize, j: usize, k: usize| (i * grid[1] + j) * grid[2] + k;
    // spread charges
    let mut qg = vec![(0.0, 0.0); total];
    struct Spl {
        base: [i64; 3],
        w: [Vec<f64>; 3],
        dw: [Vec<f64>; 3],
    }
    let mut spl = Vec::with_capacity(n);
    for i in 0..n {
        let mut base = [0i64; 3];
        let mut w: [Vec<f64>; 3] = [vec![], vec![], vec![]];
        let mut dw: [Vec<f64>; 3] = [vec![], vec![], vec![]];
        for a in 0..3 {
            let u = (pos[i][a] / len[a] - (pos[i][a] / len[a]).floor()) * grid[a] as f64;
            let fl = u.floor();
            let fr = u - fl;
            base[a] = fl as i64;
            for j in 0..order {
                w[a].push(bspline(order, fr + j as f64));
                dw[a].push(bspline_deriv(order, fr + j as f64) * grid[a] as f64 / len[a]);
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
                let kv = [2.0 * PI * signed(i, grid[0]) / len[0], 2.0 * PI * signed(j, grid[1]) / len[1], 2.0 * PI * signed(k, grid[2]) / len[2]];
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
                    // dE/dQ = 2 * pref * Re Psi, Psi = Ntot * IFFT(C FQ)
                    let de_dq = 2.0 * pref * conv[idx(ix, iy, iz)].0 * ntot;
                    let g = [s.dw[0][jx] * s.w[1][jy] * s.w[2][jz], s.w[0][jx] * s.dw[1][jy] * s.w[2][jz], s.w[0][jx] * s.w[1][jy] * s.dw[2][jz]];
                    for a in 0..3 {
                        forces[i][a] -= q[i] * de_dq * g[a];
                    }
                }
            }
        }
    }
    Recip { energy, forces, virial }
}

/// Force callback adding the direct-Ewald reciprocal, self and background terms to an MD force evaluation.
#[must_use]
pub fn ewald_extra(
    alpha: f64,
    kmax: i32,
    kc: f64,
) -> ExtraForce<3> {
    Box::new(move |sys, out| {
        let r = ewald_reciprocal(&sys.pos, &sys.charge, sys.bx.len, alpha, kmax, kc);
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
        let r = pme_reciprocal(&sys.pos, &sys.charge, sys.bx.len, alpha, grid, order, kc);
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
