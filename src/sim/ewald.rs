//! Ewald summation and particle-mesh Ewald.
//!
//! # Ewald summation and smooth particle-mesh Ewald (PME)
//!
//! Long-range Coulomb interactions in a periodic 3-D box (orthorhombic or general triclinic cell).
//! The total energy splits into a screened real-space part (use
//! [`PairTerm::Coulomb`] with `alpha > 0` in the MD force field), a
//! reciprocal-space part and the self/neutralising terms provided here.
//!
//! | Method | Cost | Notes |
//! |---|---|---|
//! | [`ewald_reciprocal`], [`EwaldReciprocalPlan`] | `O(N K^3)` | direct reciprocal sum, energy, forces and virial |
//! | [`pme_reciprocal`] | `O(N p^3 + K^3 log K)` | cardinal B-spline (even order) charge spreading, 3-D FFT from [`crate::kernels::fftconv`], analytic spline-derivative forces |
//! | [`p3m_reciprocal_cell`] | `O(N p^3 + K^3 log K)` | **P3M** with the optimal Hockney-Eastwood influence function (alias sums), `ik` differentiation, triclinic cells; clearly smaller force errors than smooth PME on the same mesh |
//! | [`ewald_energy_forces`] | `O(N^2 + N K^3)` | full reference (real space with minimum image, needs `r_cut <= L/2`) |
//! | [`ewald_reciprocal_cell`], [`pme_reciprocal_cell`], [`ewald_energy_forces_cell`] | as above | **triclinic cells** (rows of `h` are the lattice vectors): reciprocal basis, fractional-coordinate B-splines with the `h^-1` chain rule, real-space sum over all images within `r_cut` |
//!
//! Dimensions: 3-D only (the Ewald kernel `1/r` is the 3-D Coulomb law).
//! Conventions: Gaussian-like units with coupling constant `k_c`; the
//! neutralising background correction is included for charged cells.
//!
//! Parallelism: [`ewald_reciprocal_cell`] cuts the wave vectors into fixed blocks of 32 whose
//! partial sums are computed concurrently and added in block order, and the orthorhombic real-space
//! sum of [`ewald_energy_forces`] evaluates one row of pairs per task, both above a work threshold;
//! the result depends only on the problem size, never on the thread count.

use std::f64::consts::PI;
use std::f64::consts::TAU;
use std::sync::Arc;

use crate::kernels::dd::Dd;
use crate::kernels::fftconv::Complex;
use crate::kernels::fftconv::fft;
use crate::kernels::numeric_common::scale_product;
use crate::kernels::reuse::PlanCache;
use crate::kernels::reuse::ScratchPool;
use crate::kernels::reuse::f64_scratch;
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

/// Result used when a cell-based evaluation fails (all entries NaN).
fn failed_recip(n: usize) -> Recip {
    Recip {
        energy: f64::NAN,
        forces: vec![[f64::NAN; 3]; n],
        virial: [[f64::NAN; 3]; 3],
    }
}

/// Diagonal cell matrix of an orthorhombic box.
#[must_use]
pub const fn diag_cell(len: [f64; 3]) -> [[f64; 3]; 3] {
    [[len[0], 0.0, 0.0], [0.0, len[1], 0.0], [0.0, 0.0, len[2]]]
}

/// Determinant of a 3x3 matrix.
fn det3(h: &[[f64; 3]; 3]) -> f64 {
    h[0][0] * (h[1][1] * h[2][2] - h[1][2] * h[2][1])
        - h[0][1] * (h[1][0] * h[2][2] - h[1][2] * h[2][0])
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

/// `wave pairs x particles` from which the reciprocal Ewald sum runs in parallel.
const RECIPROCAL_PAR_WORK: usize = 200_000;
/// Wave vectors per block of the parallel reciprocal sum.
const RECIPROCAL_BLOCK: usize = 32;
/// `n^2` from which the orthorhombic real-space sum runs in parallel.
const REAL_SPACE_PAR_PAIRS: usize = 4_000_000;

static RECIPROCAL_VECTORS: ScratchPool<[f64; 3]> = ScratchPool::new(16, 4 * 1024 * 1024);
static RECIPROCAL_PHASES: ScratchPool<(f64, f64)> = ScratchPool::new(16, 4 * 1024 * 1024);

#[derive(Debug)]
struct ReciprocalWave {
    index: [f64; 3],
    direction: [f64; 3],
    energy: f64,
    force: f64,
    strain: f64,
}

#[derive(Default)]
struct ReciprocalTotals {
    energy: f64,
    virial: [[f64; 3]; 3],
}

/// Prepared direct reciprocal Ewald sum for a fixed periodic cell and cutoff.
///
/// The cell stores the three Cartesian lattice vectors as rows. The immutable
/// plan shares reciprocal geometry across particle configurations and calls.
/// Each retained wave represents a conjugate pair, with periodic phases
/// evaluated in fractional coordinates.
#[derive(Debug)]
pub struct EwaldReciprocalPlan {
    cell: [[f64; 3]; 3],
    inverse: [[f64; 3]; 3],
    diagonal: Option<[f64; 3]>,
    length: f64,
    volume: f64,
    alpha: f64,
    waves: Vec<ReciprocalWave>,
}

impl EwaldReciprocalPlan {
    pub(crate) const fn heap_bytes(&self) -> usize {
        self.waves.capacity() * size_of::<ReciprocalWave>()
    }

    /// Prepares all nonzero reciprocal wave pairs with `|n_i| <= kmax`.
    ///
    /// # Errors
    ///
    /// Invalid cell, screening or cutoff, unrepresentable reciprocal
    /// geometry, or inability to allocate the wave list.
    pub fn new(
        h: &[[f64; 3]; 3],
        alpha: f64,
        kmax: i32,
    ) -> Result<Self, String> {
        let length = h.iter().flatten().map(|x| x.abs()).fold(0.0, f64::max);
        if !h.iter().flatten().all(|x| x.is_finite())
            || length == 0.0
            || !alpha.is_finite()
            || alpha <= 0.0
            || kmax < 0
        {
            return Err("invalid direct Ewald cell, screening or cutoff".into());
        }
        // The raw cell volume need not fit in a float.
        let cell = h.map(|row| row.map(|value| value / length));
        let volume = det3(&cell).abs();
        let inverse = inv3(&cell).ok_or("singular cell matrix")?;
        let alpha = alpha * length;
        if !inverse.iter().flatten().all(|x| x.is_finite()) || !alpha.is_finite() || alpha == 0.0 {
            return Err("direct Ewald cell or reduced screening is not representable".into());
        }
        let b: [[f64; 3]; 3] =
            std::array::from_fn(|i| std::array::from_fn(|c| TAU * inverse[c][i]));
        let n = 2 * kmax as usize + 1;
        let pairs = n.checked_pow(3).ok_or("direct Ewald wave count overflow")? / 2;
        let mut waves = Vec::new();
        waves
            .try_reserve_exact(pairs.min(131_072))
            .map_err(|_| "direct Ewald wave list allocation failed")?;
        'modes: for nx in -kmax..=0 {
            for ny in -kmax..=kmax {
                for nz in -kmax..=kmax {
                    if (nx, ny, nz) == (0, 0, 0) {
                        break 'modes;
                    }
                    let index = [nx, ny, nz].map(f64::from);
                    let k = std::array::from_fn(|c| {
                        index[0] * b[0][c] + index[1] * b[1][c] + index[2] * b[2][c]
                    });
                    if let Some(wave) = reciprocal_wave(index, k, volume, alpha)? {
                        if waves.len() == waves.capacity() {
                            waves
                                .try_reserve(1)
                                .map_err(|_| "direct Ewald wave list allocation failed")?;
                        }
                        waves.push(wave);
                    }
                }
            }
        }
        let diagonal = (0..3)
            .all(|i| (0..3).all(|j| i == j || h[i][j] == 0.0))
            .then(|| [h[0][0], h[1][1], h[2][2]]);
        Ok(Self {
            cell: *h,
            inverse,
            diagonal,
            length,
            volume,
            alpha,
            waves,
        })
    }

    fn angles(
        &self,
        point: [f64; 3],
    ) -> Result<[f64; 3], String> {
        if let Some(periods) = self.diagonal {
            return Ok(std::array::from_fn(|a| {
                let period = periods[a];
                let value = point[a];
                let wrapped = if value.abs() > period.abs() {
                    value % period
                } else {
                    value
                };
                let fraction = wrapped / period;
                TAU * (fraction - fraction.round())
            }));
        }
        let fraction = |point: [f64; 3]| -> [f64; 3] {
            let [x, y, z] = point.map(|value| value / self.length);
            std::array::from_fn(|a| {
                x * self.inverse[0][a] + y * self.inverse[1][a] + z * self.inverse[2][a]
            })
        };
        let mut f = fraction(point);
        if f.iter().any(|value| !value.is_finite()) {
            return Err("direct Ewald fractional coordinates are not representable".into());
        }
        if f.iter().any(|x| x.abs() > 4.0) {
            let images = f.map(|value| Dd::from_f64(value.round()));
            let reduced = std::array::from_fn(|c| {
                let mut coordinate = Dd::from_f64(point[c]);
                for (i, shift) in images.iter().enumerate() {
                    coordinate -= *shift * Dd::from_f64(self.cell[i][c]);
                }
                coordinate.to_f64()
            });
            let reduced = fraction(reduced);
            if reduced.iter().all(|x| x.is_finite()) {
                f = reduced;
            }
        }
        Ok(f.map(|value| TAU * (value - value.round())))
    }

    /// Computes reciprocal forces and virial, with self and background energy.
    ///
    /// The charges correspond one-to-one to the Cartesian positions. The
    /// coupling constant scales energies, forces and virial. Fixed-order
    /// reduction keeps repeated evaluations bitwise independent of worker
    /// count.
    ///
    /// # Errors
    ///
    /// Wrong charge length, non-finite input, or unrepresentable phases or
    /// results.
    pub fn evaluate(
        &self,
        pos: &[[f64; 3]],
        q: &[f64],
        kc: f64,
    ) -> Result<Recip, String> {
        if q.len() != pos.len()
            || !kc.is_finite()
            || q.iter().any(|x| !x.is_finite())
            || pos.iter().flatten().any(|x| !x.is_finite())
        {
            return Err("invalid direct Ewald positions, charges or coupling".into());
        }
        let n = pos.len();
        let mut forces = vec![[0.0; 3]; n];
        let mut totals = ReciprocalTotals::default();
        let maximum_charge = q.iter().map(|value| value.abs()).fold(0.0, f64::max);
        if kc == 0.0 || maximum_charge == 0.0 {
            return Ok(Recip {
                energy: totals.energy,
                virial: totals.virial,
                forces,
            });
        }
        // Exact binary rescaling preserves cancellation in neutral cells.
        let exponent = maximum_charge.to_bits() & 0x7ff0_0000_0000_0000;
        let charge = f64::from_bits(exponent).max(f64::MIN_POSITIVE);
        let mut normalized = f64_scratch(n);
        normalized.extend(q.iter().map(|value| *value / charge));
        if !self.waves.is_empty() {
            let mut angles = RECIPROCAL_VECTORS.checkout_for(n);
            for &point in pos {
                angles.push(self.angles(point)?);
            }
            self.sum_waves(&angles, &normalized, &mut totals, &mut forces);
        }
        let energy_scale = [kc.abs(), charge, charge];
        let (charge_squared, net_charge) = total_charge(&normalized);
        let self_energy = -scaled_ewald(
            charge_squared / PI.sqrt(),
            &[energy_scale[0], charge, charge, self.alpha],
            &[self.length],
        )?;
        let background = if net_charge == 0.0 {
            0.0
        } else {
            scaled_ewald(
                0.5 * PI,
                &[
                    energy_scale[0],
                    charge,
                    charge,
                    net_charge.abs(),
                    net_charge.abs(),
                ],
                &[self.length, self.volume, self.alpha, self.alpha],
            )?
        };
        let reciprocal = scaled_ewald(totals.energy, &energy_scale, &[self.length])?;
        let energy = (self_energy - background + reciprocal) * kc.signum();
        if !energy.is_finite() {
            return Err("direct Ewald energy is not representable".into());
        }
        for x in 0..3 {
            for y in x..3 {
                let value =
                    scaled_ewald(totals.virial[x][y], &energy_scale, &[self.length])? * kc.signum();
                totals.virial[x][y] = value;
                totals.virial[y][x] = value;
            }
        }
        for (force, &q) in forces.iter_mut().zip(q) {
            if q == 0.0 {
                force.fill(0.0);
            } else {
                let sign = kc.signum() * q.signum();
                for component in force {
                    *component = scaled_ewald(
                        *component,
                        &[kc.abs(), charge, q.abs()],
                        &[self.length, self.length],
                    )? * sign;
                }
            }
        }
        Ok(Recip {
            energy,
            virial: totals.virial,
            forces,
        })
    }

    fn sum_waves(
        &self,
        angles: &[[f64; 3]],
        q: &[f64],
        totals: &mut ReciprocalTotals,
        forces: &mut [[f64; 3]],
    ) {
        let n = q.len();
        if self.waves.len().saturating_mul(n) < RECIPROCAL_PAR_WORK {
            let mut phases = RECIPROCAL_PHASES.checkout_for(n);
            phases.resize(n, (0.0, 0.0));
            reciprocal_block(&self.waves, angles, q, totals, forces, &mut phases);
            return;
        }

        // Group already ordered blocks into memory-bounded batches.
        // Group size affects scheduling and scratch use only.
        let per_block = n.saturating_mul(size_of::<[f64; 3]>() + size_of::<(f64, f64)>());
        let max_blocks = (16 * 1024 * 1024 / per_block).clamp(1, 32);
        let blocks = self
            .waves
            .len()
            .div_ceil(RECIPROCAL_BLOCK)
            .min(max_blocks)
            .min(crate::execution::cpu_threads().saturating_mul(4));
        let mut workspace: Vec<_> = (0..blocks)
            .map(|_| {
                let mut forces = RECIPROCAL_VECTORS.checkout_for(n);
                forces.resize(n, [0.0; 3]);
                let mut phases = RECIPROCAL_PHASES.checkout_for(n);
                phases.resize(n, (0.0, 0.0));
                (ReciprocalTotals::default(), forces, phases)
            })
            .collect();
        for group in self.waves.chunks(blocks * RECIPROCAL_BLOCK) {
            let active = &mut workspace[..group.len().div_ceil(RECIPROCAL_BLOCK)];
            crate::kernels::par::for_each_mut(active, 2, |i, (totals, forces, phases)| {
                let start = i * RECIPROCAL_BLOCK;
                let end = (start + RECIPROCAL_BLOCK).min(group.len());
                *totals = ReciprocalTotals::default();
                forces.fill([0.0; 3]);
                reciprocal_block(&group[start..end], angles, q, totals, forces, phases);
            });
            for (part, block_forces, _) in active {
                totals.energy += part.energy;
                for a in 0..3 {
                    for b in a..3 {
                        totals.virial[a][b] += part.virial[a][b];
                    }
                }
                for (force, add) in forces.iter_mut().zip(block_forces.iter()) {
                    for a in 0..3 {
                        force[a] += add[a];
                    }
                }
            }
        }
    }
}

fn reciprocal_wave(
    index: [f64; 3],
    k: [f64; 3],
    volume: f64,
    alpha: f64,
) -> Result<Option<ReciprocalWave>, String> {
    let [x, y, z] = k;
    let square = x * x + y * y + z * z;
    let norm = if square.is_normal() {
        square.sqrt()
    } else {
        x.hypot(y).hypot(z)
    };
    if !norm.is_finite() || norm == 0.0 {
        return Err("direct Ewald wave vector is not representable".into());
    }
    let screen = (0.5 * norm / alpha).powi(2);
    let exponential = (-screen).exp();
    if exponential == 0.0 {
        return Ok(None);
    }
    let direct = (4.0 * PI / volume) * exponential / square;
    let weight = if direct.is_normal() && square.is_normal() {
        direct
    } else {
        scaled_ewald(4.0 * PI * exponential, &[], &[volume, norm, norm])?
    };
    let force = (2.0 * weight) * norm;
    let strain = -2.0 * (1.0 + screen) * weight;
    if !force.is_finite() || !strain.is_finite() {
        return Err("direct Ewald wave weights are not representable".into());
    }
    Ok(Some(ReciprocalWave {
        index,
        direction: k.map(|value| value / norm),
        energy: weight,
        force,
        strain,
    }))
}

fn reciprocal_block(
    waves: &[ReciprocalWave],
    angles: &[[f64; 3]],
    charges: &[f64],
    totals: &mut ReciprocalTotals,
    forces: &mut [[f64; 3]],
    phases: &mut [(f64, f64)],
) {
    for wave in waves {
        let [x, y, z] = wave.index;
        let (mut sc, mut ss) = (0.0, 0.0);
        for ((&[a, b, c], &q), phase) in angles.iter().zip(charges).zip(&mut *phases) {
            let (s, c) = (x * a + y * b + z * c).sin_cos();
            sc += q * c;
            ss += q * s;
            *phase = (s, c);
        }
        let s2 = sc * sc + ss * ss;
        let energy = wave.energy * s2;
        totals.energy += energy;
        let strain = wave.strain * s2;
        for a in 0..3 {
            for b in a..3 {
                let diagonal = if a == b { energy } else { 0.0 };
                totals.virial[a][b] += diagonal + strain * wave.direction[a] * wave.direction[b];
            }
        }
        for ([fx, fy, fz], &(s, c)) in forces.iter_mut().zip(&*phases) {
            let f = wave.force * (sc * s - ss * c);
            *fx += f * wave.direction[0];
            *fy += f * wave.direction[1];
            *fz += f * wave.direction[2];
        }
    }
}

fn total_charge(q: &[f64]) -> (f64, f64) {
    let (mut q2, mut total, mut correction) = (0.0, 0.0_f64, 0.0);
    for &q in q {
        q2 += q * q;
        let next = total + q;
        correction += if total.abs() >= q.abs() {
            (total - next) + q
        } else {
            (q - next) + total
        };
        total = next;
    }
    (q2, total + correction)
}

fn scaled_ewald(
    value: f64,
    numerator: &[f64],
    denominator: &[f64],
) -> Result<f64, String> {
    scale_product(value, numerator, denominator).map_err(|error| error.to_string())
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
    ewald_reciprocal_cell(pos, q, &diag_cell(len), alpha, kmax, kc)
        .unwrap_or_else(|_| failed_recip(pos.len()))
}

/// Direct reciprocal-space Ewald sum for a general triclinic cell.
///
/// The rows of `h` are lattice vectors. Returns energy, forces and the
/// reciprocal virial, including self and neutralising terms in the energy.
/// Small fixed-cell plans are retained for repeated evaluations.
///
/// # Errors
///
/// Invalid cell or particle data, or an unrepresentable result.
pub fn ewald_reciprocal_cell(
    pos: &[[f64; 3]],
    q: &[f64],
    h: &[[f64; 3]; 3],
    alpha: f64,
    kmax: i32,
    kc: f64,
) -> Result<Recip, String> {
    static PLANS: PlanCache<([[u64; 3]; 3], u64, i32), EwaldReciprocalPlan> = PlanCache::new(2);
    // A negative legacy cutoff has no reciprocal modes.
    let kmax = kmax.max(0);
    if kmax > 24 {
        return EwaldReciprocalPlan::new(h, alpha, kmax)?.evaluate(pos, q, kc);
    }
    let key = (h.map(|row| row.map(f64::to_bits)), alpha.to_bits(), kmax);
    let plan = PLANS.try_get_or_insert_with(key, || {
        EwaldReciprocalPlan::new(h, alpha, kmax).map(Arc::new)
    })?;
    plan.evaluate(pos, q, kc)
}

/// Complete Ewald energy and forces for a triclinic cell.
///
/// The real-space sum runs over all
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
                            d[a] =
                                pos[i][a] - pos[j][a] + fx * h[0][a] + fy * h[1][a] + fz * h[2][a];
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
    let n = pos.len();
    if n * n >= REAL_SPACE_PAR_PAIRS {
        // One task per particle: row i visits every other particle j and
        // accumulates only the force on i (each pair is evaluated twice,
        // the energy counted with weight 1/2), in a fixed j order; the row
        // energies are added in particle order. Independent of the thread
        // count.
        let rows: Vec<(f64, [f64; 3])> = crate::kernels::par::map_indexed(n, 8, |i| {
            let (mut e_i, mut f_i) = (0.0, [0.0; 3]);
            for j in 0..n {
                if j == i {
                    continue;
                }
                let mut d = [0.0; 3];
                for a in 0..3 {
                    d[a] = pos[i][a] - pos[j][a];
                    d[a] -= len[a] * (d[a] / len[a]).round();
                }
                let r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
                if r2 < rcut * rcut {
                    let (u, fr) = term.eval(r2, q[i] * q[j]);
                    e_i += 0.5 * u;
                    for a in 0..3 {
                        f_i[a] += fr * d[a];
                    }
                }
            }
            (e_i, f_i)
        });
        for (i, (e_i, f_i)) in rows.into_iter().enumerate() {
            e += e_i;
            for a in 0..3 {
                f[i][a] += f_i[a];
            }
        }
        return (e, f);
    }
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
    pme_reciprocal_cell(pos, q, &diag_cell(len), alpha, grid, order, kc)
        .unwrap_or_else(|_| failed_recip(pos.len()))
}

/// Smooth particle-mesh Ewald for a general triclinic cell.
///
/// The cell `h` has the lattice vectors as rows. Charges are spread with cardinal B-splines on a mesh in the fractional coordinates
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
    if order < 2 || !order.is_multiple_of(2) || grid.iter().any(|&g| g < order) {
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
    let signed = |m: usize, k: usize| {
        if m <= k / 2 {
            m as f64
        } else {
            m as f64 - k as f64
        }
    };
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
                    let gs = [
                        s.dw[0][jx] * s.w[1][jy] * s.w[2][jz],
                        s.w[0][jx] * s.dw[1][jy] * s.w[2][jz],
                        s.w[0][jx] * s.w[1][jy] * s.dw[2][jz],
                    ];
                    for c in 0..3 {
                        let gc = gs[0] * hinv[c][0] + gs[1] * hinv[c][1] + gs[2] * hinv[c][2];
                        forces[i][c] -= q[i] * de_dq * gc;
                    }
                }
            }
        }
    }
    Ok(Recip {
        energy,
        forces,
        virial,
    })
}

/// P3M reciprocal-space energy, forces and virial with the **optimal (Hockney-Eastwood)
/// influence function** for a general triclinic cell.
///
/// Charges are assigned with cardinal B-splines of even `order` on a mesh in fractional
/// coordinates; the influence function minimises the rms force error over all charge
/// configurations for `ik` differentiation,
/// `G(k) = sum_m U^2(k_m) (k . k_m) R(k_m) / (k^2 [sum_m U^2(k_m)]^2)` with
/// `R(k) = 4 pi exp(-k^2/4 alpha^2)/k^2`, `U` the transform of the assignment function and the sum
/// running over the `(2 n_alias + 1)^3` aliases `k_m = k + m K b`. The field is obtained by
/// `ik` differentiation of the mesh potential (three inverse FFTs) and interpolated with the
/// assignment function; the virial uses the analytic Ewald reciprocal form with the optimal `G`
/// (exact up to the aliasing error). The self and neutralising-background terms are included.
///
/// # Errors
/// Fails for a singular cell, an odd or too small spline order, or a mesh smaller than the order.
pub fn p3m_reciprocal_cell(
    pos: &[[f64; 3]],
    q: &[f64],
    h: &[[f64; 3]; 3],
    alpha: f64,
    grid: [usize; 3],
    order: usize,
    kc: f64,
) -> Result<Recip, String> {
    const N_ALIAS: i64 = 2;
    if order < 2 || !order.is_multiple_of(2) || grid.iter().any(|&g| g < order) {
        return Err("P3M needs an even spline order >= 2 and meshes at least as large".into());
    }
    let vol = det3(h).abs();
    let hinv = inv3(h).ok_or("singular cell matrix")?;
    let bvec = reciprocal_basis(h)?;
    let n = pos.len();
    let total = grid[0] * grid[1] * grid[2];
    let idx = |i: usize, j: usize, k: usize| (i * grid[1] + j) * grid[2] + k;
    // charge assignment
    let mut qg = vec![(0.0, 0.0); total];
    let mut spl: Vec<([i64; 3], [Vec<f64>; 3])> = Vec::with_capacity(n);
    for i in 0..n {
        let mut base = [0i64; 3];
        let mut w: [Vec<f64>; 3] = [vec![], vec![], vec![]];
        for a in 0..3 {
            let s = pos[i][0] * hinv[0][a] + pos[i][1] * hinv[1][a] + pos[i][2] * hinv[2][a];
            let u = (s - s.floor()) * grid[a] as f64;
            let fl = u.floor();
            base[a] = fl as i64;
            for j in 0..order {
                w[a].push(bspline(order, u - fl + j as f64));
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
        spl.push((base, w));
    }
    fft3(&mut qg, grid, false);
    let signed = |m: usize, k: usize| {
        if m <= k / 2 {
            m as i64
        } else {
            m as i64 - k as i64
        }
    };
    let sinc_pow = |x: f64| -> f64 {
        if x.abs() < 1e-12 {
            1.0
        } else {
            (x.sin() / x).powi(2 * order as i32)
        }
    };
    let kvec = |m: [i64; 3]| -> [f64; 3] {
        std::array::from_fn(|c| {
            m[0] as f64 * bvec[0][c] + m[1] as f64 * bvec[1][c] + m[2] as f64 * bvec[2][c]
        })
    };
    let four_alpha2 = 4.0 * alpha * alpha;
    let pref = kc / (2.0 * vol);
    let mut energy = charge_terms(q, vol, alpha, kc);
    let mut virial = [[0.0; 3]; 3];
    let mut phi = vec![(0.0, 0.0); total];
    for i in 0..grid[0] {
        for j in 0..grid[1] {
            for k in 0..grid[2] {
                if i == 0 && j == 0 && k == 0 {
                    continue;
                }
                let nn = [signed(i, grid[0]), signed(j, grid[1]), signed(k, grid[2])];
                let kv = kvec(nn);
                let k2 = kv[0] * kv[0] + kv[1] * kv[1] + kv[2] * kv[2];
                let (mut num, mut den) = (0.0, 0.0);
                for ax in -N_ALIAS..=N_ALIAS {
                    for ay in -N_ALIAS..=N_ALIAS {
                        for az in -N_ALIAS..=N_ALIAS {
                            let m = [
                                nn[0] + ax * grid[0] as i64,
                                nn[1] + ay * grid[1] as i64,
                                nn[2] + az * grid[2] as i64,
                            ];
                            let km = kvec(m);
                            let km2 = km[0] * km[0] + km[1] * km[1] + km[2] * km[2];
                            let u2 = sinc_pow(std::f64::consts::PI * m[0] as f64 / grid[0] as f64)
                                * sinc_pow(std::f64::consts::PI * m[1] as f64 / grid[1] as f64)
                                * sinc_pow(std::f64::consts::PI * m[2] as f64 / grid[2] as f64);
                            den += u2;
                            num += u2
                                * (kv[0] * km[0] + kv[1] * km[1] + kv[2] * km[2])
                                * 4.0
                                * PI
                                * (-km2 / four_alpha2).exp()
                                / km2;
                        }
                    }
                }
                let g = num / (k2 * den * den);
                let f = qg[idx(i, j, k)];
                let mag = f.0 * f.0 + f.1 * f.1;
                let ek = pref * g * mag;
                energy += ek;
                let w = -2.0 * (1.0 / k2 + 1.0 / four_alpha2);
                for x in 0..3 {
                    for y in 0..3 {
                        virial[x][y] += ek * (if x == y { 1.0 } else { 0.0 } + w * kv[x] * kv[y]);
                    }
                }
                phi[idx(i, j, k)] = (g * f.0, g * f.1);
            }
        }
    }
    // ik differentiation: E_c(k) = -i k_c phi(k); Nyquist components of the derivative vanish
    let ntot = total as f64;
    let mut fields: Vec<Vec<(f64, f64)>> = Vec::with_capacity(3);
    for c in 0..3 {
        let mut ec = vec![(0.0, 0.0); total];
        for i in 0..grid[0] {
            for j in 0..grid[1] {
                for k in 0..grid[2] {
                    let nn = [signed(i, grid[0]), signed(j, grid[1]), signed(k, grid[2])];
                    let nyq = (grid[0].is_multiple_of(2) && i == grid[0] / 2)
                        || (grid[1].is_multiple_of(2) && j == grid[1] / 2)
                        || (grid[2].is_multiple_of(2) && k == grid[2] / 2);
                    if nyq {
                        continue;
                    }
                    let kv = kvec(nn);
                    let p = phi[idx(i, j, k)];
                    // -i k_c (a + i b) = k_c b - i k_c a
                    ec[idx(i, j, k)] = (kv[c] * p.1, -kv[c] * p.0);
                }
            }
        }
        fft3(&mut ec, grid, true);
        fields.push(ec);
    }
    let mut forces = vec![[0.0; 3]; n];
    for i in 0..n {
        let (base, w) = &spl[i];
        for jx in 0..order {
            for jy in 0..order {
                for jz in 0..order {
                    let ix = (base[0] - jx as i64).rem_euclid(grid[0] as i64) as usize;
                    let iy = (base[1] - jy as i64).rem_euclid(grid[1] as i64) as usize;
                    let iz = (base[2] - jz as i64).rem_euclid(grid[2] as i64) as usize;
                    let wt = w[0][jx] * w[1][jy] * w[2][jz];
                    for c in 0..3 {
                        forces[i][c] += q[i] * wt * fields[c][idx(ix, iy, iz)].0 * ntot * kc / vol;
                    }
                }
            }
        }
    }
    Ok(Recip {
        energy,
        forces,
        virial,
    })
}

/// Force callback adding the direct-Ewald reciprocal, self and background terms to an MD force evaluation.
pub fn ewald_extra(
    alpha: f64,
    kmax: i32,
    kc: f64,
) -> ExtraForce<3> {
    Box::new(move |sys, out| {
        let h = sys.bx.cell.unwrap_or_else(|| diag_cell(sys.bx.len));
        let r = ewald_reciprocal_cell(&sys.pos, &sys.charge, &h, alpha, kmax, kc)
            .unwrap_or_else(|_| failed_recip(sys.pos.len()));
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
pub fn pme_extra(
    alpha: f64,
    grid: [usize; 3],
    order: usize,
    kc: f64,
) -> ExtraForce<3> {
    Box::new(move |sys, out| {
        let h = sys.bx.cell.unwrap_or_else(|| diag_cell(sys.bx.len));
        let r = pme_reciprocal_cell(&sys.pos, &sys.charge, &h, alpha, grid, order, kc)
            .unwrap_or_else(|_| failed_recip(sys.pos.len()));
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

/// Force callback adding the P3M reciprocal (optimal influence function), self and background terms.
pub fn p3m_extra(
    alpha: f64,
    grid: [usize; 3],
    order: usize,
    kc: f64,
) -> ExtraForce<3> {
    Box::new(move |sys, out| {
        let h = sys.bx.cell.unwrap_or_else(|| diag_cell(sys.bx.len));
        let r = p3m_reciprocal_cell(&sys.pos, &sys.charge, &h, alpha, grid, order, kc)
            .unwrap_or_else(|_| failed_recip(sys.pos.len()));
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
