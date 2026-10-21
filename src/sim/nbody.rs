//! # Tree and multipole N-body solvers for gravity and Coulomb
//!
//! Long-range interactions with the Gauss-law kernel of the dimension: the
//! field of a source `s_j` at `x_j` on a point `x` is
//! `s_j (x_j - x) / (|x_j - x|^2 + eps^2)^(D/2)`, i.e. `1/r^(D-1)` forces
//! (constant sheets in 1-D, `1/r` in 2-D, `1/r^2` in 3-D).
//!
//! | Solver | Dimensions | Accuracy |
//! |---|---|---|
//! | [`BhTree`] Barnes-Hut (monopole + dipole about the `|s|`-weighted centre, so signed Coulomb charges work) | 1, 2, 3 | `O(theta^2)` |
//! | [`direct_field`] | any | exact reference, `O(N^2)` |
//! | [`fmm2d_field`] adaptive fast multipole method (quadtree, dual-tree traversal, M2M/M2L/L2L, order `p`) | 2 | geometric in `p` (about `1e-10` for `p = 24`) |
//! | [`fmm3d_field`] / [`crate::sim::fmm::fmm3d`] adaptive octree FMM (Cartesian Taylor expansions, M2M/M2L/L2L, potential and field) | 3 | geometric in `p` and the acceptance `theta` |
//!
//! Helpers [`gravity_accelerations`] and [`coulomb_forces`] wrap the fields
//! for use with the integrators of [`crate::sim::symplectic`].
//! The FMM algorithms live in [`crate::sim::fmm`]: the far field uses the same Plummer-softened
//! kernel as the near field, and [`crate::sim::fmm::periodic_fmm3d`] provides the periodic 3-D FMM
//! (triclinic cells). Hermite integrators with block time steps, 1PN terms and KS binaries are in
//! [`crate::sim::hermite`]. Limitations: Barnes-Hut and the 2-D FMM are open-boundary only;
//! the periodic FMM treats the lattice beyond the near images with the unsoftened kernel.

const NONE: u32 = u32::MAX;
const MAX_DEPTH: usize = 48;

#[derive(Clone, Debug)]
struct Node<const D: usize> {
    center: [f64; D],
    half: f64,
    s: f64,
    sabs: f64,
    com: [f64; D],
    dip: [f64; D],
    children: [u32; 8],
    bodies: Vec<u32>,
}

/// Barnes-Hut tree (`D <= 3`).
#[derive(Clone, Debug)]
pub struct BhTree<const D: usize> {
    nodes: Vec<Node<D>>,
    pos: Vec<[f64; D]>,
    src: Vec<f64>,
    /// Opening angle.
    pub theta: f64,
    /// Plummer softening.
    pub eps: f64,
}

fn kernel<const D: usize>(
    d: &[f64; D],
    eps: f64,
) -> [f64; D] {
    let r2: f64 = d.iter().map(|x| x * x).sum::<f64>() + eps * eps;
    let inv = r2.powf(-(D as f64) / 2.0);
    let mut o = [0.0; D];
    for a in 0..D {
        o[a] = d[a] * inv;
    }
    o
}

impl<const D: usize> BhTree<D> {
    /// Builds the tree for sources `src` at `pos`.
    ///
    /// # Errors
    /// Fails for `D > 3`, empty input or length mismatch.
    pub fn new(
        pos: &[[f64; D]],
        src: &[f64],
        theta: f64,
        eps: f64,
    ) -> Result<Self, String> {
        if D == 0 || D > 3 {
            return Err("Barnes-Hut supports dimensions 1 to 3".into());
        }
        if pos.is_empty() || pos.len() != src.len() {
            return Err("positions and sources must be non-empty and of equal length".into());
        }
        let mut lo = [f64::INFINITY; D];
        let mut hi = [f64::NEG_INFINITY; D];
        for p in pos {
            for a in 0..D {
                lo[a] = lo[a].min(p[a]);
                hi[a] = hi[a].max(p[a]);
            }
        }
        let mut center = [0.0; D];
        let mut half = 0.0_f64;
        for a in 0..D {
            center[a] = f64::midpoint(lo[a], hi[a]);
            half = half.max(0.5 * (hi[a] - lo[a]));
        }
        half = half * 1.000_001 + 1e-12;
        let root = Node { center, half, s: 0.0, sabs: 0.0, com: center, dip: [0.0; D], children: [NONE; 8], bodies: Vec::new() };
        let mut t = Self { nodes: vec![root], pos: pos.to_vec(), src: src.to_vec(), theta, eps };
        for i in 0..pos.len() {
            t.insert(0, i as u32, 0);
        }
        t.moments(0);
        Ok(t)
    }

    fn octant(
        &self,
        node: usize,
        p: &[f64; D],
    ) -> usize {
        let mut o = 0;
        for a in 0..D {
            if p[a] >= self.nodes[node].center[a] {
                o |= 1 << a;
            }
        }
        o
    }

    fn insert(
        &mut self,
        node: usize,
        body: u32,
        depth: usize,
    ) {
        let is_leaf = self.nodes[node].children.iter().all(|&c| c == NONE);
        if is_leaf {
            if self.nodes[node].bodies.is_empty() || depth >= MAX_DEPTH {
                self.nodes[node].bodies.push(body);
                return;
            }
            // split: push the existing bodies down
            let old = std::mem::take(&mut self.nodes[node].bodies);
            for b in old {
                self.insert_child(node, b, depth);
            }
        }
        self.insert_child(node, body, depth);
    }

    fn insert_child(
        &mut self,
        node: usize,
        body: u32,
        depth: usize,
    ) {
        let p = self.pos[body as usize];
        let o = self.octant(node, &p);
        if self.nodes[node].children[o] == NONE {
            let h = self.nodes[node].half * 0.5;
            let mut c = self.nodes[node].center;
            for a in 0..D {
                c[a] += if (o >> a) & 1 == 1 { h } else { -h };
            }
            self.nodes.push(Node { center: c, half: h, s: 0.0, sabs: 0.0, com: c, dip: [0.0; D], children: [NONE; 8], bodies: Vec::new() });
            let id = (self.nodes.len() - 1) as u32;
            self.nodes[node].children[o] = id;
        }
        let child = self.nodes[node].children[o] as usize;
        self.insert(child, body, depth + 1);
    }

    fn moments(
        &mut self,
        node: usize,
    ) {
        let children = self.nodes[node].children;
        let bodies = self.nodes[node].bodies.clone();
        for &c in &children {
            if c != NONE {
                self.moments(c as usize);
            }
        }
        // gather the contributing point sources: bodies, or child (s, centre) pairs
        let mut s = 0.0;
        let mut sabs = 0.0;
        let mut wc = [0.0; D];
        for &b in &bodies {
            let (sj, p) = (self.src[b as usize], self.pos[b as usize]);
            s += sj;
            sabs += sj.abs();
            for a in 0..D {
                wc[a] += sj.abs() * p[a];
            }
        }
        for &c in &children {
            if c != NONE {
                let ch = &self.nodes[c as usize];
                s += ch.s;
                sabs += ch.sabs;
                for a in 0..D {
                    wc[a] += ch.sabs * ch.com[a];
                }
            }
        }
        let com = if sabs > 0.0 {
            let mut c = [0.0; D];
            for a in 0..D {
                c[a] = wc[a] / sabs;
            }
            c
        } else {
            self.nodes[node].center
        };
        let mut dip = [0.0; D];
        for &b in &bodies {
            let (sj, p) = (self.src[b as usize], self.pos[b as usize]);
            for a in 0..D {
                dip[a] += sj * (p[a] - com[a]);
            }
        }
        for &c in &children {
            if c != NONE {
                let ch = &self.nodes[c as usize];
                for a in 0..D {
                    dip[a] += ch.dip[a] + ch.s * (ch.com[a] - com[a]);
                }
            }
        }
        let n = &mut self.nodes[node];
        n.s = s;
        n.sabs = sabs;
        n.com = com;
        n.dip = dip;
    }

    /// Field `sum_j s_j K(x_j - x)` at `x`, skipping the body `skip` (pass `usize::MAX` for none).
    #[must_use]
    pub fn field(
        &self,
        x: &[f64; D],
        skip: usize,
    ) -> [f64; D] {
        let mut out = [0.0; D];
        let mut stack = vec![0usize];
        while let Some(ni) = stack.pop() {
            let n = &self.nodes[ni];
            let is_leaf = n.children.iter().all(|&c| c == NONE);
            if is_leaf {
                for &b in &n.bodies {
                    if b as usize == skip {
                        continue;
                    }
                    let mut d = [0.0; D];
                    for a in 0..D {
                        d[a] = self.pos[b as usize][a] - x[a];
                    }
                    let k = kernel(&d, self.eps);
                    for a in 0..D {
                        out[a] += self.src[b as usize] * k[a];
                    }
                }
                continue;
            }
            let mut d = [0.0; D];
            for a in 0..D {
                d[a] = n.com[a] - x[a];
            }
            let dist = d.iter().map(|v| v * v).sum::<f64>().sqrt();
            // the target must lie outside the node's bounding box for the expansion to be used
            let mut inside = true;
            for a in 0..D {
                if (x[a] - n.center[a]).abs() > n.half {
                    inside = false;
                }
            }
            if !inside && 2.0 * n.half < self.theta * dist {
                let k = kernel(&d, self.eps);
                let dn = dist.powi(D as i32 + 2);
                let d0 = dist.powi(D as i32);
                for b in 0..D {
                    let mut dp = 0.0;
                    for a in 0..D {
                        let dk = if a == b { 1.0 / d0 } else { 0.0 } - D as f64 * d[a] * d[b] / dn;
                        dp += n.dip[a] * dk;
                    }
                    out[b] += n.s * k[b] + dp;
                }
            } else {
                for &c in &n.children {
                    if c != NONE {
                        stack.push(c as usize);
                    }
                }
            }
        }
        out
    }
}

/// Exact `O(N^2)` field of all sources at every body (self excluded).
#[must_use]
pub fn direct_field<const D: usize>(
    pos: &[[f64; D]],
    src: &[f64],
    eps: f64,
) -> Vec<[f64; D]> {
    let n = pos.len();
    let mut out = vec![[0.0; D]; n];
    for i in 0..n {
        for j in 0..n {
            if i == j {
                continue;
            }
            let mut d = [0.0; D];
            for a in 0..D {
                d[a] = pos[j][a] - pos[i][a];
            }
            let k = kernel(&d, eps);
            for a in 0..D {
                out[i][a] += src[j] * k[a];
            }
        }
    }
    out
}

/// Barnes-Hut gravitational accelerations `a_i = G sum_j m_j K(x_j - x_i)`.
///
/// # Errors
/// Propagates tree construction errors.
pub fn gravity_accelerations<const D: usize>(
    pos: &[[f64; D]],
    mass: &[f64],
    g: f64,
    eps: f64,
    theta: f64,
) -> Result<Vec<[f64; D]>, String> {
    let t = BhTree::new(pos, mass, theta, eps)?;
    Ok((0..pos.len())
        .map(|i| {
            let mut f = t.field(&pos[i], i);
            for c in &mut f {
                *c *= g;
            }
            f
        })
        .collect())
}

/// Barnes-Hut Coulomb forces `F_i = -k q_i sum_j q_j K(x_j - x_i)` (like charges repel).
///
/// # Errors
/// Propagates tree construction errors.
pub fn coulomb_forces<const D: usize>(
    pos: &[[f64; D]],
    charge: &[f64],
    k: f64,
    eps: f64,
    theta: f64,
) -> Result<Vec<[f64; D]>, String> {
    let t = BhTree::new(pos, charge, theta, eps)?;
    Ok((0..pos.len())
        .map(|i| {
            let mut f = t.field(&pos[i], i);
            for c in &mut f {
                *c *= -k * charge[i];
            }
            f
        })
        .collect())
}


// ---------------------------------------------------------------------------
// FMM front ends (the algorithms live in `crate::sim::fmm`)
// ---------------------------------------------------------------------------

/// Field `sum_j s_j (x_j - x_i) / |x_j - x_i|^2` at every point (self excluded) by an
/// adaptive 2-D FMM of order `p` with about `ppb` points per leaf box.
///
/// # Errors
/// Fails for empty input or length mismatch.
pub fn fmm2d_field(
    pos: &[[f64; 2]],
    src: &[f64],
    p: usize,
    ppb: usize,
) -> Result<Vec<[f64; 2]>, String> {
    crate::sim::fmm::fmm2d(pos, src, &crate::sim::fmm::FmmOptions { order: p, leaf_size: ppb, theta: 0.4, eps: 0.0 })
}

/// Field `sum_j s_j (x_j - x_i) / |x_j - x_i|^3` at every point (self excluded) by the
/// adaptive 3-D FMM of order `p` with about `ppb` points per leaf box.
///
/// # Errors
/// Fails for empty input or length mismatch.
pub fn fmm3d_field(
    pos: &[[f64; 3]],
    src: &[f64],
    p: usize,
    ppb: usize,
) -> Result<Vec<[f64; 3]>, String> {
    let r = crate::sim::fmm::fmm3d(pos, src, &crate::sim::fmm::FmmOptions { order: p, leaf_size: ppb, theta: 0.5, eps: 0.0 })?;
    Ok(r.field)
}
