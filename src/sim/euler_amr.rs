//! # Block-structured adaptive mesh refinement for Euler and shallow water (1D and 2D)
//!
//! Berger–Colella style AMR on top of the finite-volume kernels of
//! [`euler_nd`](crate::sim::euler_nd):
//!
//! * a hierarchy of up to `max_level` refinement levels with ratio 2 in each
//!   direction; level 0 covers the whole domain, finer levels are lists of
//!   non-overlapping rectangular **patches** that are re-created from error
//!   indicators (density jump) every `regrid_interval` coarse steps by
//!   recursive bisection of the tagged region (efficiency-driven clustering,
//!   buffer layers, proper nesting);
//! * **time sub-cycling**: every level takes steps of half the coarser time
//!   step (same CFL number), integrated with SSP-RK3;
//! * ghost cells of fine patches come from neighbouring patches of the same
//!   level or by conservative slope-limited linear interpolation in space and
//!   linear interpolation in time from the coarser level (recursively, so
//!   patches need not be strictly nested);
//! * conservation: after a fine level has caught up, fine data are **averaged
//!   down** into the covered coarse cells, and the **flux registers** correct the coarse
//!   cells next to a coarse–fine interface with the difference between the
//!   coarse flux and the time-integrated fine fluxes (reflux), so mass,
//!   momentum and energy are conserved to round-off with periodic or
//!   reflective boundaries;
//! * reconstruction (first order, MUSCL, WENO5), Riemann solvers and boundary
//!   conditions (transmissive, reflective, periodic, fixed) are those of
//!   [`FvConfig`].
//!
//! Limitations: 1D and 2D only, refinement ratio 2, characteristic
//! non-reflecting boundaries and 3D are not supported, patches are
//! aligned with the coarse grid and the regridding re-creates the whole fine
//! hierarchy (data are copied where the old patches overlap the new ones).

use crate::sim::euler_nd::FvConfig;
use crate::sim::euler_nd::FvGrid;
use crate::sim::euler_nd::GasBc;
use crate::sim::euler_nd::Law;
use crate::sim::euler_nd::State;
use crate::sim::euler_nd::admissible;
use crate::sim::euler_nd::face_states;
use crate::sim::euler_nd::riemann_flux;
use crate::sim::euler_nd::sound_speed;
use crate::sim::euler_nd::swap_dir;
use crate::sim::euler_nd::to_cons;
use crate::sim::euler_nd::to_prim;

const GH: usize = 3;

/// AMR parameters.
#[derive(Clone, Copy, Debug)]
pub struct AmrConfig {
    /// Finite-volume scheme and boundary conditions.
    pub fv: FvConfig,
    /// Number of refinement levels above the base grid.
    pub max_level: usize,
    /// Coarse steps between regridding operations.
    pub regrid_interval: usize,
    /// A cell is tagged where `|ρ_{i+1} - ρ_{i-1}| / (2 ρ_i)` exceeds this value.
    pub threshold: f64,
    /// Layers of cells added around tagged cells.
    pub buffer: usize,
    /// Minimum fraction of tagged cells in a patch before it is bisected.
    pub efficiency: f64,
}

impl AmrConfig {
    /// Default AMR parameters around a scheme: two levels, regrid every 4 steps, threshold 0.05.
    #[must_use]
    pub const fn new(fv: FvConfig) -> Self {
        Self {
            fv,
            max_level: 2,
            regrid_interval: 4,
            threshold: 0.05,
            buffer: 2,
            efficiency: 0.6,
        }
    }
}

struct Patch {
    lo: [usize; 2],
    n: [usize; 2],
    q: Vec<State>,
    old: Vec<State>,
    new: Vec<State>,
    /// Time-averaged (RK-weighted) numerical flux per axis of the last step.
    flux: [Vec<State>; 2],
    /// Flux registers `[axis][side]` over the coarse transverse faces.
    reg: [[Vec<State>; 2]; 2],
}

impl Patch {
    fn new(
        lo: [usize; 2],
        n: [usize; 2],
        q: Vec<State>,
    ) -> Self {
        Self {
            lo,
            n,
            old: q.clone(),
            new: q.clone(),
            q,
            flux: [Vec::new(), Vec::new()],
            reg: [[Vec::new(), Vec::new()], [Vec::new(), Vec::new()]],
        }
    }

    const fn flat(
        &self,
        i: [usize; 2],
    ) -> usize {
        (i[0] - self.lo[0]) + self.n[0] * (i[1] - self.lo[1])
    }

    const fn contains(
        &self,
        i: [usize; 2],
    ) -> bool {
        i[0] >= self.lo[0]
            && i[0] < self.lo[0] + self.n[0]
            && i[1] >= self.lo[1]
            && i[1] < self.lo[1] + self.n[1]
    }
}

struct Level {
    patches: Vec<Patch>,
    t_old: f64,
    t_new: f64,
}

enum Mapped {
    Cell([usize; 2], [bool; 2]),
    Fixed(State),
}

struct Amr {
    cfg: AmrConfig,
    dim: usize,
    n0: [usize; 2],
    h0: [f64; 2],
    nv: usize,
    levels: Vec<Level>,
}

fn minmod(
    a: f64,
    b: f64,
) -> f64 {
    if a * b <= 0.0 {
        0.0
    } else if a.abs() < b.abs() {
        a
    } else {
        b
    }
}

impl Amr {
    const fn dom(
        &self,
        l: usize,
    ) -> [isize; 2] {
        let f = 1isize << l;
        [
            self.n0[0] as isize * f,
            if self.dim == 2 {
                self.n0[1] as isize * f
            } else {
                1
            },
        ]
    }

    fn h(
        &self,
        l: usize,
        a: usize,
    ) -> f64 {
        self.h0[a] / (1u64 << l) as f64
    }

    fn find(
        &self,
        l: usize,
        i: [usize; 2],
    ) -> Option<(usize, usize)> {
        self.levels[l]
            .patches
            .iter()
            .enumerate()
            .find(|(_, p)| p.contains(i))
            .map(|(k, p)| (k, p.flat(i)))
    }

    /// Maps an index outside the domain through the boundary conditions.
    fn bc_map(
        &self,
        l: usize,
        idx: [isize; 2],
    ) -> Mapped {
        let dom = self.dom(l);
        let law = self.cfg.fv.law;
        let mut j = idx;
        let mut flip = [false; 2];
        for a in 0..self.dim {
            if j[a] >= 0 && j[a] < dom[a] {
                continue;
            }
            let hi = j[a] >= dom[a];
            match self.cfg.fv.bcs[2 * a + usize::from(hi)] {
                | GasBc::Periodic => j[a] = j[a].rem_euclid(dom[a]),
                | GasBc::Reflective => {
                    j[a] = if hi {
                        2 * dom[a] - 1 - j[a]
                    } else {
                        -1 - j[a]
                    };
                    flip[a] = !flip[a];
                },
                | GasBc::Fixed(w) => return Mapped::Fixed(to_cons(law, w)),
                | _ => j[a] = j[a].clamp(0, dom[a] - 1),
            }
        }
        Mapped::Cell([j[0] as usize, j[1].max(0) as usize], flip)
    }

    /// Conserved state of cell `idx` of level `l` at time `tau`.
    ///
    /// `stage` holds the stage arrays of the level currently being stepped (same-level
    /// neighbours use them instead of the time interpolation between `old` and `new`).
    fn sample(
        &self,
        l: usize,
        idx: [isize; 2],
        tau: f64,
        stage: Option<(usize, &[Vec<State>])>,
    ) -> State {
        let dom = self.dom(l);
        let inside = (0..2).all(|a| idx[a] >= 0 && idx[a] < dom[a]);
        if !inside {
            return match self.bc_map(l, idx) {
                | Mapped::Fixed(s) => s,
                | Mapped::Cell(j, flip) => {
                    let mut s = self.sample(l, [j[0] as isize, j[1] as isize], tau, stage);
                    for a in 0..2 {
                        if flip[a] {
                            s[1 + a] = -s[1 + a];
                        }
                    }
                    s
                },
            };
        }
        let u = [idx[0] as usize, idx[1] as usize];
        if let Some((k, flat)) = self.find(l, u) {
            if let Some((sl, arrays)) = stage
                && sl == l
            {
                return arrays[k][flat];
            }
            let lev = &self.levels[l];
            let p = &lev.patches[k];
            let w = if lev.t_new > lev.t_old {
                ((tau - lev.t_old) / (lev.t_new - lev.t_old)).clamp(0.0, 1.0)
            } else {
                1.0
            };
            let mut s = [0.0; 5];
            for c in 0..5 {
                s[c] = (1.0 - w) * p.old[flat][c] + w * p.new[flat][c];
            }
            return s;
        }
        // not covered at this level: conservative slope-limited prolongation from the parent
        let parent = [idx[0].div_euclid(2), idx[1].div_euclid(2)];
        let q0 = self.sample(l - 1, parent, tau, stage);
        let mut out = q0;
        for a in 0..self.dim {
            let mut ep = parent;
            ep[a] += 1;
            let mut em = parent;
            em[a] -= 1;
            let (qp, qm) = (
                self.sample(l - 1, ep, tau, stage),
                self.sample(l - 1, em, tau, stage),
            );
            let off = if idx[a].rem_euclid(2) == 0 {
                -0.25
            } else {
                0.25
            };
            for c in 0..5 {
                out[c] += off * minmod(qp[c] - q0[c], q0[c] - qm[c]);
            }
        }
        out
    }

    /// Padded copy (with ghost cells) of a patch's stage array.
    fn padded(
        &self,
        l: usize,
        k: usize,
        stage: &[Vec<State>],
        tau: f64,
    ) -> Vec<State> {
        let p = &self.levels[l].patches[k];
        let g = [GH, if self.dim == 2 { GH } else { 0 }];
        let np = [p.n[0] + 2 * g[0], p.n[1] + 2 * g[1]];
        let mut out = vec![[0.0; 5]; np[0] * np[1]];
        for j in 0..np[1] {
            for i in 0..np[0] {
                let gi = [
                    p.lo[0] as isize + i as isize - g[0] as isize,
                    p.lo[1] as isize + j as isize - g[1] as isize,
                ];
                let interior = i >= g[0] && i < g[0] + p.n[0] && j >= g[1] && j < g[1] + p.n[1];
                out[i + np[0] * j] = if interior {
                    stage[k][(i - g[0]) + p.n[0] * (j - g[1])]
                } else {
                    self.sample(l, gi, tau, Some((l, stage)))
                };
            }
        }
        out
    }

    /// Numerical fluxes of a patch from its padded state.
    fn fluxes(
        &self,
        padded: &[State],
        n: [usize; 2],
    ) -> [Vec<State>; 2] {
        let cfg = &self.cfg.fv;
        let law = cfg.law;
        let g = [GH, if self.dim == 2 { GH } else { 0 }];
        let np = [n[0] + 2 * g[0], n[1] + 2 * g[1]];
        let wp: Vec<State> = padded.iter().map(|&c| to_prim(law, c)).collect();
        let mut out: [Vec<State>; 2] = [Vec::new(), Vec::new()];
        for d in 0..self.dim {
            let o = 1 - d;
            let nlines = if self.dim == 2 { n[o] } else { 1 };
            let nf = n[d] + 1;
            let mut flux = vec![[0.0; 5]; nf * nlines];
            let mut line = vec![[0.0; 5]; np[d]];
            for t in 0..nlines {
                for k in 0..np[d] {
                    let mut ij = [0usize; 2];
                    ij[d] = k;
                    ij[o] = t + g[o];
                    line[k] = wp[ij[0] + np[0] * ij[1]];
                }
                for f in 0..nf {
                    let (mut l, mut r) = face_states(cfg.recon, &line, GH + f, self.nv);
                    if !admissible(law, &l) || !admissible(law, &r) {
                        l = line[GH + f - 1];
                        r = line[GH + f];
                    }
                    swap_dir(&mut l, d);
                    swap_dir(&mut r, d);
                    let mut fl = riemann_flux(law, cfg.riemann, &l, &r);
                    swap_dir(&mut fl, d);
                    // flux layout: axis 0 -> f + nf t ; axis 1 -> t + n0 f
                    let at = if d == 0 {
                        f + nf * t
                    } else {
                        t + n[0] * f
                    };
                    flux[at] = fl;
                }
            }
            out[d] = flux;
        }
        out
    }

    fn residual(
        &self,
        n: [usize; 2],
        flux: &[Vec<State>; 2],
        l: usize,
    ) -> Vec<State> {
        let mut res = vec![[0.0; 5]; n[0] * n[1]];
        for j in 0..n[1] {
            for i in 0..n[0] {
                for c in 0..5 {
                    let mut r = -(flux[0][i + 1 + (n[0] + 1) * j][c]
                        - flux[0][i + (n[0] + 1) * j][c])
                        / self.h(l, 0);
                    if self.dim == 2 {
                        r -= (flux[1][i + n[0] * (j + 1)][c] - flux[1][i + n[0] * j][c])
                            / self.h(l, 1);
                    }
                    res[i + n[0] * j][c] = r;
                }
            }
        }
        res
    }

    /// One SSP-RK3 step of level `l` from `t` over `dt`.
    fn step_level(
        &mut self,
        l: usize,
        t: f64,
        dt: f64,
    ) -> Result<(), String> {
        let np = self.levels[l].patches.len();
        for p in &mut self.levels[l].patches {
            p.old.clone_from(&p.q);
        }
        self.levels[l].t_old = t;
        self.levels[l].t_new = t + dt;
        let q0: Vec<Vec<State>> = self.levels[l].patches.iter().map(|p| p.q.clone()).collect();
        let mut eff: Vec<[Vec<State>; 2]> = Vec::new();
        let weights = [1.0 / 6.0, 1.0 / 6.0, 4.0 / 6.0];
        let times = [t, t + dt, t + 0.5 * dt];
        let mut qcur = q0.clone();
        for s in 0..3 {
            let mut resid = Vec::with_capacity(np);
            let mut fl_all = Vec::with_capacity(np);
            for k in 0..np {
                let pad = self.padded(l, k, &qcur, times[s]);
                let n = self.levels[l].patches[k].n;
                let fl = self.fluxes(&pad, n);
                resid.push(self.residual(n, &fl, l));
                fl_all.push(fl);
            }
            if s == 0 {
                eff = fl_all
                    .iter()
                    .map(|f| {
                        [
                            f[0].iter().map(|v| v.map(|x| x * weights[0])).collect(),
                            f[1].iter().map(|v| v.map(|x| x * weights[0])).collect(),
                        ]
                    })
                    .collect();
            } else {
                for k in 0..np {
                    for a in 0..2 {
                        for (e, f) in eff[k][a].iter_mut().zip(&fl_all[k][a]) {
                            for c in 0..5 {
                                e[c] += weights[s] * f[c];
                            }
                        }
                    }
                }
            }
            let mut next = Vec::with_capacity(np);
            for k in 0..np {
                let v: Vec<State> = (0..qcur[k].len())
                    .map(|c| {
                        let mut o = [0.0; 5];
                        for m in 0..5 {
                            o[m] = match s {
                                | 0 => q0[k][c][m] + dt * resid[k][c][m],
                                | 1 => {
                                    0.75 * q0[k][c][m]
                                        + 0.25 * (qcur[k][c][m] + dt * resid[k][c][m])
                                },
                                | _ => {
                                    q0[k][c][m] / 3.0
                                        + 2.0 / 3.0 * (qcur[k][c][m] + dt * resid[k][c][m])
                                },
                            };
                        }
                        o
                    })
                    .collect();
                next.push(v);
            }
            qcur = next;
        }
        for (k, p) in self.levels[l].patches.iter_mut().enumerate() {
            p.q.clone_from(&qcur[k]);
            p.new.clone_from(&qcur[k]);
            p.flux = std::mem::take(&mut eff[k]);
            if p.q.iter().any(|c| c.iter().any(|v| !v.is_finite())) {
                return Err("non-finite state in the AMR solver".into());
            }
        }
        // accumulate the flux registers of fine patches
        if l >= 1 {
            let dim = self.dim;
            for p in &mut self.levels[l].patches {
                for a in 0..dim {
                    let o = 1 - a;
                    let ntc = if dim == 2 { p.n[o] / 2 } else { 1 };
                    for s in 0..2 {
                        if p.reg[a][s].len() != ntc {
                            p.reg[a][s] = vec![[0.0; 5]; ntc];
                        }
                        let face = if s == 0 { 0 } else { p.n[a] };
                        for tc in 0..ntc {
                            let sub = if dim == 2 { 2 } else { 1 };
                            let mut sum = [0.0; 5];
                            for u in 0..sub {
                                let t = tc * sub + u;
                                let at = if a == 0 {
                                    face + (p.n[0] + 1) * t
                                } else {
                                    t + p.n[0] * face
                                };
                                for c in 0..5 {
                                    sum[c] += p.flux[a][at][c];
                                }
                            }
                            for c in 0..5 {
                                p.reg[a][s][tc][c] += dt * sum[c] / sub as f64;
                            }
                        }
                    }
                }
            }
        }
        Ok(())
    }

    fn advance(
        &mut self,
        l: usize,
        t: f64,
        dt: f64,
    ) -> Result<(), String> {
        self.step_level(l, t, dt)?;
        if l < self.cfg.max_level
            && l + 1 < self.levels.len()
            && !self.levels[l + 1].patches.is_empty()
        {
            for p in &mut self.levels[l + 1].patches {
                p.reg = [[Vec::new(), Vec::new()], [Vec::new(), Vec::new()]];
            }
            self.advance(l + 1, t, 0.5 * dt)?;
            self.advance(l + 1, t + 0.5 * dt, 0.5 * dt)?;
            self.average_down(l + 1);
            self.reflux(l + 1, dt);
        }
        Ok(())
    }

    /// Replaces the coarse cells under the patches of level `lf` by the mean of their children.
    fn average_down(
        &mut self,
        lf: usize,
    ) {
        let sub = if self.dim == 2 { 4.0 } else { 2.0 };
        let fine = std::mem::take(&mut self.levels[lf].patches);
        for fp in &fine {
            let (cn, clo) = (
                [
                    fp.n[0] / 2,
                    if self.dim == 2 {
                        fp.n[1] / 2
                    } else {
                        1
                    },
                ],
                [
                    fp.lo[0] / 2,
                    if self.dim == 2 {
                        fp.lo[1] / 2
                    } else {
                        0
                    },
                ],
            );
            for cj in 0..cn[1] {
                for ci in 0..cn[0] {
                    let cidx = [clo[0] + ci, clo[1] + cj];
                    let mut sum = [0.0; 5];
                    for dj in 0..if self.dim == 2 { 2 } else { 1 } {
                        for di in 0..2 {
                            let v = fp.q[(2 * ci + di) + fp.n[0] * (2 * cj + dj)];
                            for c in 0..5 {
                                sum[c] += v[c];
                            }
                        }
                    }
                    if let Some((k, flat)) = self.find(lf - 1, cidx) {
                        let cp = &mut self.levels[lf - 1].patches[k];
                        for c in 0..5 {
                            cp.q[flat][c] = sum[c] / sub;
                            cp.new[flat][c] = sum[c] / sub;
                        }
                    }
                }
            }
        }
        self.levels[lf].patches = fine;
    }

    /// Flux-register correction of the coarse cells next to the patches of level `lf`.
    fn reflux(
        &mut self,
        lf: usize,
        dt_c: f64,
    ) {
        let lc = lf - 1;
        let dom = self.dom(lc);
        let fine = std::mem::take(&mut self.levels[lf].patches);
        for fp in &fine {
            for a in 0..self.dim {
                let o = 1 - a;
                let ntc = if self.dim == 2 {
                    fp.n[o] / 2
                } else {
                    1
                };
                let clo = [
                    fp.lo[0] / 2,
                    if self.dim == 2 {
                        fp.lo[1] / 2
                    } else {
                        0
                    },
                ];
                for s in 0..2 {
                    if fp.reg[a][s].len() != ntc {
                        continue;
                    }
                    let mut ca = if s == 0 {
                        clo[a] as isize - 1
                    } else {
                        (clo[a] + fp.n[a] / 2) as isize
                    };
                    if ca < 0 || ca >= dom[a] {
                        if self.cfg.fv.bcs[2 * a] != GasBc::Periodic {
                            continue;
                        }
                        ca = ca.rem_euclid(dom[a]);
                    }
                    for tc in 0..ntc {
                        let mut cidx = [0usize; 2];
                        cidx[a] = ca as usize;
                        cidx[o] = if self.dim == 2 {
                            clo[o] + tc
                        } else {
                            0
                        };
                        // skip when the neighbour is itself refined
                        let fidx = [cidx[0] * 2, cidx[1] * 2];
                        if fine.iter().any(|p| p.contains(fidx)) {
                            continue;
                        }
                        let Some((k, flat)) = self.find(lc, cidx) else {
                            continue;
                        };
                        let hc = self.h(lc, a);
                        let cp = &mut self.levels[lc].patches[k];
                        let local = [cidx[0] - cp.lo[0], cidx[1] - cp.lo[1]];
                        let face = if s == 0 {
                            local[a] + 1
                        } else {
                            local[a]
                        };
                        let at = if a == 0 {
                            face + (cp.n[0] + 1) * local[1]
                        } else {
                            local[0] + cp.n[0] * face
                        };
                        let sign = if s == 0 { -1.0 } else { 1.0 };
                        for c in 0..5 {
                            let corr = sign * (fp.reg[a][s][tc][c] - dt_c * cp.flux[a][at][c]) / hc;
                            cp.q[flat][c] += corr;
                            cp.new[flat][c] += corr;
                        }
                    }
                }
            }
        }
        self.levels[lf].patches = fine;
    }

    fn max_rate(
        &self,
        l: usize,
    ) -> f64 {
        let law = self.cfg.fv.law;
        let mut best: f64 = 0.0;
        for p in &self.levels[l].patches {
            for &c in &p.q {
                let w = to_prim(law, c);
                let cs = sound_speed(law, &w);
                let r: f64 = (0..self.dim)
                    .map(|a| (w[1 + a].abs() + cs) / self.h(l, a))
                    .sum();
                best = best.max(r);
            }
        }
        best
    }

    /// Conserved state of a cell of level `l` without any time interpolation (current data).
    fn current(
        &self,
        l: usize,
        idx: [isize; 2],
    ) -> State {
        self.sample_now(l, idx)
    }

    fn sample_now(
        &self,
        l: usize,
        idx: [isize; 2],
    ) -> State {
        // the current data of every level equal `new` (= `q`); set the interpolation weight to one
        let dom = self.dom(l);
        let inside = (0..2).all(|a| idx[a] >= 0 && idx[a] < dom[a]);
        if !inside {
            return match self.bc_map(l, idx) {
                | Mapped::Fixed(s) => s,
                | Mapped::Cell(j, flip) => {
                    let mut s = self.sample_now(l, [j[0] as isize, j[1] as isize]);
                    for a in 0..2 {
                        if flip[a] {
                            s[1 + a] = -s[1 + a];
                        }
                    }
                    s
                },
            };
        }
        let u = [idx[0] as usize, idx[1] as usize];
        if let Some((k, flat)) = self.find(l, u) {
            return self.levels[l].patches[k].q[flat];
        }
        let parent = [idx[0].div_euclid(2), idx[1].div_euclid(2)];
        let q0 = self.sample_now(l - 1, parent);
        let mut out = q0;
        for a in 0..self.dim {
            let mut ep = parent;
            ep[a] += 1;
            let mut em = parent;
            em[a] -= 1;
            let (qp, qm) = (self.sample_now(l - 1, ep), self.sample_now(l - 1, em));
            let off = if idx[a].rem_euclid(2) == 0 {
                -0.25
            } else {
                0.25
            };
            for c in 0..5 {
                out[c] += off * minmod(qp[c] - q0[c], q0[c] - qm[c]);
            }
        }
        out
    }

    /// Creates the fine hierarchy from error indicators.
    fn regrid(&mut self) {
        let t = self.levels[0].t_new;
        let dim = self.dim;
        let old_levels: Vec<Vec<Patch>> = self
            .levels
            .iter_mut()
            .skip(1)
            .map(|lv| std::mem::take(&mut lv.patches))
            .collect();
        for l in 0..self.cfg.max_level {
            let dom = self.dom(l);
            let nd = [dom[0] as usize, dom[1] as usize];
            // coverage of level l and tags
            let mut covered = vec![false; nd[0] * nd[1]];
            let mut tags = vec![false; nd[0] * nd[1]];
            for p in &self.levels[l].patches {
                for j in 0..p.n[1] {
                    for i in 0..p.n[0] {
                        covered[(p.lo[0] + i) + nd[0] * (p.lo[1] + j)] = true;
                    }
                }
            }
            let law = self.cfg.fv.law;
            let rho = |idx: [isize; 2]| to_prim(law, self.current(l, idx))[0];
            for j in 0..nd[1] {
                for i in 0..nd[0] {
                    if !covered[i + nd[0] * j] {
                        continue;
                    }
                    let c = [i as isize, j as isize];
                    let mut ind: f64 = 0.0;
                    for a in 0..dim {
                        let (mut p, mut m) = (c, c);
                        p[a] += 1;
                        m[a] -= 1;
                        ind = ind.max((rho(p) - rho(m)).abs() / (2.0 * rho(c).abs().max(1e-12)));
                    }
                    if ind > self.cfg.threshold {
                        tags[i + nd[0] * j] = true;
                    }
                }
            }
            // buffer
            let b = self.cfg.buffer as isize;
            let mut buffered = tags.clone();
            for j in 0..nd[1] {
                for i in 0..nd[0] {
                    if !tags[i + nd[0] * j] {
                        continue;
                    }
                    for dj in if dim == 2 { -b..=b } else { 0..=0 } {
                        for di in -b..=b {
                            let (ii, jj) = (i as isize + di, j as isize + dj);
                            if ii >= 0 && ii < dom[0] && jj >= 0 && jj < dom[1] {
                                buffered[ii as usize + nd[0] * jj as usize] = true;
                            }
                        }
                    }
                }
            }
            // proper nesting: only cells whose +-2 neighbourhood is covered by level l may be refined
            let mut nested = vec![false; nd[0] * nd[1]];
            for j in 0..nd[1] {
                for i in 0..nd[0] {
                    if !covered[i + nd[0] * j] {
                        continue;
                    }
                    let mut ok = true;
                    for dj in if dim == 2 { -2..=2 } else { 0..=0 } {
                        for di in -2..=2 {
                            let (ii, jj) = (i as isize + di, j as isize + dj);
                            if ii >= 0
                                && ii < dom[0]
                                && jj >= 0
                                && jj < dom[1]
                                && !covered[ii as usize + nd[0] * jj as usize]
                            {
                                ok = false;
                            }
                        }
                    }
                    nested[i + nd[0] * j] = ok;
                }
            }
            let allowed: Vec<bool> = buffered
                .iter()
                .zip(&nested)
                .map(|(a, b)| *a && *b)
                .collect();
            let boxes = cluster(&allowed, &nested, nd, self.cfg.efficiency);
            let mut patches = Vec::new();
            for (clo, cn) in boxes {
                let flo = [clo[0] * 2, if dim == 2 { clo[1] * 2 } else { 0 }];
                let fnn = [cn[0] * 2, if dim == 2 { cn[1] * 2 } else { 1 }];
                let mut q = Vec::with_capacity(fnn[0] * fnn[1]);
                for j in 0..fnn[1] {
                    for i in 0..fnn[0] {
                        let gi = [flo[0] + i, flo[1] + j];
                        // copy from an old patch of the same level when possible
                        let from_old = old_levels.get(l).and_then(|ps| {
                            ps.iter().find(|p| p.contains(gi)).map(|p| p.q[p.flat(gi)])
                        });
                        q.push(from_old.unwrap_or_else(|| {
                            self.sample_now(l + 1, [gi[0] as isize, gi[1] as isize])
                        }));
                    }
                }
                patches.push(Patch::new(flo, fnn, q));
            }
            self.levels[l + 1].patches = patches;
            self.levels[l + 1].t_old = t;
            self.levels[l + 1].t_new = t;
        }
    }
}

/// Recursive bisection of the tagged cells into boxes `(lo, n)`.
fn cluster(
    tags: &[bool],
    valid: &[bool],
    nd: [usize; 2],
    eff: f64,
) -> Vec<([usize; 2], [usize; 2])> {
    let mut out = Vec::new();
    let tagged = |i: usize, j: usize| tags[i + nd[0] * j];
    // bounding box of a region
    let bbox = |lo: [usize; 2], hi: [usize; 2]| -> Option<([usize; 2], [usize; 2], usize)> {
        let (mut bl, mut bh, mut cnt) = ([usize::MAX; 2], [0usize; 2], 0usize);
        for j in lo[1]..hi[1] {
            for i in lo[0]..hi[0] {
                if tagged(i, j) {
                    bl[0] = bl[0].min(i);
                    bl[1] = bl[1].min(j);
                    bh[0] = bh[0].max(i + 1);
                    bh[1] = bh[1].max(j + 1);
                    cnt += 1;
                }
            }
        }
        (cnt > 0).then_some((bl, bh, cnt))
    };
    let mut stack = vec![([0usize, 0usize], nd)];
    while let Some((lo, hi)) = stack.pop() {
        let Some((bl, bh, cnt)) = bbox(lo, hi) else {
            continue;
        };
        let size = [bh[0] - bl[0], bh[1] - bl[1]];
        let fraction = cnt as f64 / (size[0] * size[1]) as f64;
        let axis = usize::from(size[1] > size[0]);
        let all_valid = (bl[1]..bh[1]).all(|j| (bl[0]..bh[0]).all(|i| valid[i + nd[0] * j]));
        if all_valid && (fraction >= eff || size[axis] < 4) {
            out.push((bl, size));
            continue;
        }
        if size[axis] < 2 {
            continue;
        }
        // signature along the long axis: split at the emptiest slice in the middle half
        let (a0, a1) = (bl[axis], bh[axis]);
        let mid_cut = usize::midpoint(a0, a1);
        let mut best = (usize::MAX, mid_cut);
        for cut in a0 + size[axis] / 4..=a1 - size[axis] / 4 {
            if cut <= a0 || cut >= a1 {
                continue;
            }
            let mut count = 0;
            let slice = if cut == 0 { 0 } else { cut - 1 };
            for t in bl[1 - axis]..bh[1 - axis] {
                let (i, j) = if axis == 0 {
                    (slice, t)
                } else {
                    (t, slice)
                };
                if tagged(i, j) {
                    count += 1;
                }
            }
            if count < best.0
                || (count == best.0 && cut.abs_diff(mid_cut) < best.1.abs_diff(mid_cut))
            {
                best = (count, cut);
            }
        }
        let cut = best.1;
        let (mut hi1, mut lo2) = (bh, bl);
        hi1[axis] = cut;
        lo2[axis] = cut;
        stack.push((bl, hi1));
        stack.push((lo2, bh));
    }
    out
}

/// One refinement level of an AMR solution.
#[derive(Clone, Debug)]
pub struct AmrLevelData {
    /// Level number (0 is the base grid).
    pub level: usize,
    /// Patches as `(lower index, size, conserved states)` in the index space of the level.
    pub patches: Vec<([usize; 2], [usize; 2], Vec<State>)>,
}

/// Final state of an AMR run.
#[derive(Clone, Debug)]
pub struct AmrSolution {
    /// Base grid.
    pub grid: FvGrid,
    /// Hierarchy.
    pub levels: Vec<AmrLevelData>,
    /// Final time.
    pub t: f64,
    /// Number of coarse steps.
    pub steps: usize,
    law: Law,
}

impl AmrSolution {
    /// Number of cells per level.
    #[must_use]
    pub fn cells_per_level(&self) -> Vec<usize> {
        self.levels
            .iter()
            .map(|l| l.patches.iter().map(|p| p.1[0] * p.1[1]).sum())
            .collect()
    }

    fn covered(
        &self,
        level: usize,
        idx: [usize; 2],
    ) -> bool {
        self.levels.get(level).is_some_and(|l| {
            l.patches.iter().any(|(lo, n, _)| {
                idx[0] >= lo[0] && idx[0] < lo[0] + n[0] && idx[1] >= lo[1] && idx[1] < lo[1] + n[1]
            })
        })
    }

    /// Integrals of the conserved variables over the composite grid (every cell counted at its finest level).
    #[must_use]
    pub fn totals(&self) -> State {
        let mut s = [0.0; 5];
        for lv in &self.levels {
            let vol: f64 = (0..self.grid.dim)
                .map(|a| self.grid.h(a) / (1u64 << lv.level) as f64)
                .product();
            for (lo, n, q) in &lv.patches {
                for j in 0..n[1] {
                    for i in 0..n[0] {
                        let idx = [2 * (lo[0] + i), 2 * (lo[1] + j)];
                        if self.covered(lv.level + 1, idx) {
                            continue;
                        }
                        for c in 0..5 {
                            s[c] += q[i + n[0] * j][c] * vol;
                        }
                    }
                }
            }
        }
        s
    }

    /// `L1` norm of `component(prim) - exact(x)` over the composite grid (every cell counted at its finest level).
    pub fn composite_l1_error(
        &self,
        component: usize,
        exact: impl Fn([f64; 3]) -> f64,
    ) -> f64 {
        let dim = self.grid.dim;
        let mut sum = 0.0;
        for lv in &self.levels {
            let f = (1u64 << lv.level) as f64;
            let h = [
                self.grid.h(0) / f,
                if dim == 2 {
                    self.grid.h(1) / f
                } else {
                    1.0
                },
            ];
            let vol = h[0] * if dim == 2 { h[1] } else { 1.0 };
            for (lo, n, q) in &lv.patches {
                for j in 0..n[1] {
                    for i in 0..n[0] {
                        let idx = [2 * (lo[0] + i), 2 * (lo[1] + j)];
                        if self.covered(lv.level + 1, idx) {
                            continue;
                        }
                        let x = [
                            self.grid.lo[0] + ((lo[0] + i) as f64 + 0.5) * h[0],
                            if dim == 2 {
                                self.grid.lo[1] + ((lo[1] + j) as f64 + 0.5) * h[1]
                            } else {
                                0.0
                            },
                            0.0,
                        ];
                        let w = to_prim(self.law, q[i + n[0] * j]);
                        sum += (w[component] - exact(x)).abs() * vol;
                    }
                }
            }
        }
        sum
    }

    /// Primitive variables on the uniform grid of level `level` (coarser data are injected piecewise constantly).
    #[must_use]
    pub fn primitive_on_level(
        &self,
        level: usize,
    ) -> Vec<State> {
        let f = 1usize << level;
        let n = [
            self.grid.n[0] * f,
            if self.grid.dim == 2 {
                self.grid.n[1] * f
            } else {
                1
            },
        ];
        let mut out = vec![[0.0; 5]; n[0] * n[1]];
        for j in 0..n[1] {
            for i in 0..n[0] {
                // finest covering level not above `level`
                let mut found = None;
                for l in (0..=level.min(self.levels.len() - 1)).rev() {
                    let sh = level - l;
                    let idx = [i >> sh, j >> sh];
                    if let Some((lo, nn, q)) = self.levels[l].patches.iter().find(|(lo, nn, _)| {
                        idx[0] >= lo[0]
                            && idx[0] < lo[0] + nn[0]
                            && idx[1] >= lo[1]
                            && idx[1] < lo[1] + nn[1]
                    }) {
                        found = Some(q[(idx[0] - lo[0]) + nn[0] * (idx[1] - lo[1])]);
                        break;
                    }
                }
                out[i + n[0] * j] = to_prim(self.law, found.unwrap_or([1.0, 0.0, 0.0, 0.0, 1.0]));
            }
        }
        out
    }
}

/// Runs the AMR solver from primitive initial data (evaluated at the cell centres of every level) to `t_end`.
///
/// # Errors
/// For 3D grids, unsupported boundary conditions, a CFL number outside
/// `(0, 1)`, non-admissible initial data or a blown-up solution.
pub fn solve_amr(
    grid: &FvGrid,
    cfg: &AmrConfig,
    init: &dyn Fn([f64; 3]) -> State,
    t_end: f64,
) -> Result<AmrSolution, String> {
    if grid.dim > 2 {
        return Err("AMR is available in 1D and 2D".into());
    }
    if !(cfg.fv.cfl > 0.0 && cfg.fv.cfl < 1.0) {
        return Err("cfl must be in (0, 1)".into());
    }
    if cfg
        .fv
        .bcs
        .iter()
        .any(|b| matches!(b, GasBc::NonReflecting(_)))
    {
        return Err("non-reflecting boundaries are not supported by the AMR solver".into());
    }
    for a in 0..grid.dim {
        if (cfg.fv.bcs[2 * a] == GasBc::Periodic) != (cfg.fv.bcs[2 * a + 1] == GasBc::Periodic) {
            return Err("periodic must be set on both faces of an axis".into());
        }
    }
    if !grid.n[0].is_multiple_of(2) || (grid.dim == 2 && !grid.n[1].is_multiple_of(2)) {
        return Err("the base grid needs an even number of cells per axis".into());
    }
    let nv = match cfg.fv.law {
        | Law::Euler { .. } => 5,
        | Law::ShallowWater { .. } => 1 + grid.dim,
    };
    let n0 = [
        grid.n[0],
        if grid.dim == 2 {
            grid.n[1]
        } else {
            1
        },
    ];
    let h0 = [
        grid.h(0),
        if grid.dim == 2 {
            grid.h(1)
        } else {
            1.0
        },
    ];
    let mut q = Vec::with_capacity(n0[0] * n0[1]);
    for j in 0..n0[1] {
        for i in 0..n0[0] {
            let x = [
                grid.lo[0] + (i as f64 + 0.5) * h0[0],
                if grid.dim == 2 {
                    grid.lo[1] + (j as f64 + 0.5) * h0[1]
                } else {
                    0.0
                },
                0.0,
            ];
            let w = init(x);
            if !admissible(cfg.fv.law, &w) {
                return Err("initial state is not admissible".into());
            }
            q.push(to_cons(cfg.fv.law, w));
        }
    }
    let mut levels = vec![Level {
        patches: vec![Patch::new([0, 0], n0, q)],
        t_old: 0.0,
        t_new: 0.0,
    }];
    for _ in 0..cfg.max_level {
        levels.push(Level {
            patches: Vec::new(),
            t_old: 0.0,
            t_new: 0.0,
        });
    }
    let mut amr = Amr {
        cfg: *cfg,
        dim: grid.dim,
        n0,
        h0,
        nv,
        levels,
    };
    // re-evaluate the initial data on the refined levels from the exact function
    for _ in 0..=cfg.max_level {
        amr.regrid();
        amr.init_fine(init, grid);
    }
    let (mut t, mut steps) = (0.0, 0usize);
    while t < t_end * (1.0 - 1e-14) {
        if steps > 0 && steps % cfg.regrid_interval.max(1) == 0 {
            amr.regrid();
        }
        let mut dt = f64::INFINITY;
        for l in 0..=cfg.max_level {
            if amr.levels[l].patches.is_empty() {
                continue;
            }
            let rate = amr.max_rate(l);
            if !rate.is_finite() || rate <= 0.0 {
                return Err("non-finite wave speed".into());
            }
            dt = dt.min(cfg.fv.cfl / rate * (1u64 << l) as f64);
        }
        let dt = dt.min(t_end - t);
        amr.advance(0, t, dt)?;
        t += dt;
        steps += 1;
        if steps > 2_000_000 {
            return Err("too many steps".into());
        }
    }
    let levels = amr
        .levels
        .iter()
        .enumerate()
        .map(|(level, lv)| {
            AmrLevelData {
                level,
                patches: lv
                    .patches
                    .iter()
                    .map(|p| (p.lo, p.n, p.q.clone()))
                    .collect(),
            }
        })
        .collect();
    Ok(AmrSolution {
        grid: *grid,
        levels,
        t,
        steps,
        law: cfg.fv.law,
    })
}

impl Amr {
    /// Overwrites fine patches by the initial function sampled at their own cell centres.
    fn init_fine(
        &mut self,
        init: &dyn Fn([f64; 3]) -> State,
        grid: &FvGrid,
    ) {
        let law = self.cfg.fv.law;
        let dim = self.dim;
        for l in 1..self.levels.len() {
            let hs = [self.h(l, 0), self.h(l, 1)];
            for p in &mut self.levels[l].patches {
                for j in 0..p.n[1] {
                    for i in 0..p.n[0] {
                        let x = [
                            grid.lo[0] + ((p.lo[0] + i) as f64 + 0.5) * hs[0],
                            if dim == 2 {
                                grid.lo[1] + ((p.lo[1] + j) as f64 + 0.5) * hs[1]
                            } else {
                                0.0
                            },
                            0.0,
                        ];
                        let w = init(x);
                        if admissible(law, &w) {
                            p.q[i + p.n[0] * j] = to_cons(law, w);
                        }
                    }
                }
                p.old.clone_from(&p.q);
                p.new.clone_from(&p.q);
            }
        }
        // consistent coarse data below the fine patches
        for l in (1..self.levels.len()).rev() {
            self.average_down(l);
        }
    }
}
