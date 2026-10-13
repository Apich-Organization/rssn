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
//! | [`fmm2d_field`] fast multipole method (uniform quadtree, M2M/M2L/L2L, order `p`) | 2 | geometric in `p` (about `1e-10` for `p = 24`) |
//!
//! Helpers [`gravity_accelerations`] and [`coulomb_forces`] wrap the fields
//! for use with the integrators of [`crate::sim::symplectic`].
//! Limitations: the FMM is two-dimensional (log potential) and non-adaptive;
//! the 3-D FMM is not implemented (use Barnes-Hut).
#![allow(
    clippy::arithmetic_side_effects,
    clippy::assigning_clones,
    clippy::cast_lossless,
    clippy::cast_possible_truncation,
    clippy::chunks_exact_to_as_chunks,
    clippy::cognitive_complexity,
    clippy::derive_partial_eq_without_eq,
    clippy::doc_markdown,
    clippy::excessive_precision,
    clippy::float_cmp,
    clippy::if_not_else,
    clippy::imprecise_flops,
    clippy::indexing_slicing,
    clippy::int_plus_one,
    clippy::items_after_statements,
    clippy::manual_let_else,
    clippy::manual_memcpy,
    clippy::manual_midpoint,
    clippy::many_single_char_names,
    clippy::missing_const_for_fn,
    clippy::module_name_repetitions,
    clippy::needless_pass_by_value,
    clippy::needless_range_loop,
    clippy::option_if_let_else,
    clippy::similar_names,
    clippy::single_call_fn,
    clippy::single_match_else,
    clippy::suboptimal_flops,
    clippy::suspicious_operation_groupings,
    clippy::too_long_first_doc_paragraph,
    clippy::too_many_arguments,
    clippy::too_many_lines,
    clippy::type_complexity,
    clippy::unreadable_literal,
    clippy::useless_let_if_seq,
    clippy::while_float
)]

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
            center[a] = 0.5 * (lo[a] + hi[a]);
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
// 2-D fast multipole method
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Debug, Default)]
struct Cx {
    re: f64,
    im: f64,
}

impl Cx {
    const fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }

    fn add(
        self,
        o: Self,
    ) -> Self {
        Self::new(self.re + o.re, self.im + o.im)
    }

    fn sub(
        self,
        o: Self,
    ) -> Self {
        Self::new(self.re - o.re, self.im - o.im)
    }

    fn mul(
        self,
        o: Self,
    ) -> Self {
        Self::new(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)
    }

    fn scale(
        self,
        s: f64,
    ) -> Self {
        Self::new(self.re * s, self.im * s)
    }

    fn inv(self) -> Self {
        let n = self.re * self.re + self.im * self.im;
        Self::new(self.re / n, -self.im / n)
    }
}

fn binom_table(p: usize) -> Vec<Vec<f64>> {
    let n = 2 * p + 2;
    let mut c = vec![vec![0.0; n + 1]; n + 1];
    for i in 0..=n {
        c[i][0] = 1.0;
        for j in 1..=i {
            c[i][j] = c[i - 1][j - 1] + if j <= i - 1 { c[i - 1][j] } else { 0.0 };
        }
    }
    c
}

/// Field `sum_j s_j (x_j - x_i) / |x_j - x_i|^2` at every point (self excluded) by a 2-D FMM of order `p`.
///
/// Uses a uniform quadtree with about `ppb` points per finest box.
///
/// # Errors
/// Fails for empty input or length mismatch.
pub fn fmm2d_field(
    pos: &[[f64; 2]],
    src: &[f64],
    p: usize,
    ppb: usize,
) -> Result<Vec<[f64; 2]>, String> {
    let n = pos.len();
    if n == 0 || src.len() != n {
        return Err("positions and sources must be non-empty and of equal length".into());
    }
    let (mut lo, mut hi) = ([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]);
    for q in pos {
        for a in 0..2 {
            lo[a] = lo[a].min(q[a]);
            hi[a] = hi[a].max(q[a]);
        }
    }
    let size = (hi[0] - lo[0]).max(hi[1] - lo[1]) * 1.000_001 + 1e-12;
    let mut levels = 2usize;
    while (n as f64) / (4.0_f64.powi(levels as i32)) > ppb as f64 && levels < 9 {
        levels += 1;
    }
    let nb = 1usize << levels;
    let bsize = size / nb as f64;
    let binom = binom_table(p);
    let centre = |lvl: usize, ix: usize, iy: usize| -> Cx {
        let bs = size / (1usize << lvl) as f64;
        Cx::new(lo[0] + (ix as f64 + 0.5) * bs, lo[1] + (iy as f64 + 0.5) * bs)
    };
    // assign points to leaf boxes
    let mut leaf: Vec<Vec<usize>> = vec![Vec::new(); nb * nb];
    for (i, q) in pos.iter().enumerate() {
        let ix = (((q[0] - lo[0]) / bsize) as usize).min(nb - 1);
        let iy = (((q[1] - lo[1]) / bsize) as usize).min(nb - 1);
        leaf[iy * nb + ix].push(i);
    }
    // multipoles per level: mp[lvl][box][k]
    let mut mp: Vec<Vec<Vec<Cx>>> = Vec::new();
    for lvl in 0..=levels {
        let m = 1usize << lvl;
        mp.push(vec![vec![Cx::default(); p + 1]; m * m]);
    }
    for iy in 0..nb {
        for ix in 0..nb {
            let c = centre(levels, ix, iy);
            let b = &mut mp[levels][iy * nb + ix];
            for &i in &leaf[iy * nb + ix] {
                let dz = Cx::new(pos[i][0], pos[i][1]).sub(c);
                let mut pw = Cx::new(1.0, 0.0);
                for k in 0..=p {
                    b[k] = b[k].add(pw.scale(src[i]));
                    pw = pw.mul(dz);
                }
            }
        }
    }
    // upward pass (M2M)
    for lvl in (0..levels).rev() {
        let m = 1usize << lvl;
        for iy in 0..m {
            for ix in 0..m {
                let cp = centre(lvl, ix, iy);
                let mut acc = vec![Cx::default(); p + 1];
                for (dx, dy) in [(0, 0), (1, 0), (0, 1), (1, 1)] {
                    let (cx_, cy_) = (2 * ix + dx, 2 * iy + dy);
                    let t = centre(lvl + 1, cx_, cy_).sub(cp);
                    let child = &mp[lvl + 1][cy_ * (2 * m) + cx_];
                    let mut tp = vec![Cx::new(1.0, 0.0); p + 1];
                    for k in 1..=p {
                        tp[k] = tp[k - 1].mul(t);
                    }
                    for k in 0..=p {
                        for l in 0..=k {
                            acc[k] = acc[k].add(tp[k - l].mul(child[l]).scale(binom[k][l]));
                        }
                    }
                }
                mp[lvl][iy * m + ix] = acc;
            }
        }
    }
    // downward pass: locals
    let mut local: Vec<Vec<Vec<Cx>>> = Vec::new();
    for lvl in 0..=levels {
        let m = 1usize << lvl;
        local.push(vec![vec![Cx::default(); p + 1]; m * m]);
    }
    for lvl in 2..=levels {
        let m = 1usize << lvl;
        for iy in 0..m {
            for ix in 0..m {
                let ct = centre(lvl, ix, iy);
                let mut d = vec![Cx::default(); p + 1];
                // L2L from parent
                {
                    let (px, py) = (ix / 2, iy / 2);
                    let cp = centre(lvl - 1, px, py);
                    let t = ct.sub(cp);
                    let par = &local[lvl - 1][py * (m / 2) + px];
                    let mut tp = vec![Cx::new(1.0, 0.0); p + 1];
                    for k in 1..=p {
                        tp[k] = tp[k - 1].mul(t);
                    }
                    for mm in 0..=p {
                        for l in mm..=p {
                            d[mm] = d[mm].add(par[l].mul(tp[l - mm]).scale(binom[l][mm]));
                        }
                    }
                }
                // M2L from the interaction list
                let (px, py) = ((ix / 2) as i64, (iy / 2) as i64);
                let mh = (m / 2) as i64;
                for ny in (py - 1)..=(py + 1) {
                    for nx in (px - 1)..=(px + 1) {
                        if nx < 0 || ny < 0 || nx >= mh || ny >= mh {
                            continue;
                        }
                        for (cdx, cdy) in [(0, 0), (1, 0), (0, 1), (1, 1)] {
                            let (sx, sy) = (2 * nx + cdx, 2 * ny + cdy);
                            if (sx - ix as i64).abs() <= 1 && (sy - iy as i64).abs() <= 1 {
                                continue;
                            }
                            let cs = centre(lvl, sx as usize, sy as usize);
                            let u = ct.sub(cs);
                            let ui = u.inv();
                            let src_mp = &mp[lvl][sy as usize * m + sx as usize];
                            // d_l = sum_k b_k (-1)^l C(k+l, l) / u^(k+l+1)
                            let mut upow = vec![Cx::new(1.0, 0.0); 2 * p + 2];
                            for k in 1..2 * p + 2 {
                                upow[k] = upow[k - 1].mul(ui);
                            }
                            for l in 0..=p {
                                let sign = if l % 2 == 0 { 1.0 } else { -1.0 };
                                for k in 0..=p {
                                    d[l] = d[l].add(src_mp[k].mul(upow[k + l + 1]).scale(sign * binom[k + l][l]));
                                }
                            }
                        }
                    }
                }
                local[lvl][iy * m + ix] = d;
            }
        }
    }
    // evaluation
    let mut out = vec![[0.0; 2]; n];
    for iy in 0..nb {
        for ix in 0..nb {
            let ct = centre(levels, ix, iy);
            let d = &local[levels][iy * nb + ix];
            for &i in &leaf[iy * nb + ix] {
                let z = Cx::new(pos[i][0], pos[i][1]);
                let w = z.sub(ct);
                // Horner for the local expansion
                let mut acc = Cx::default();
                for l in (0..=p).rev() {
                    acc = acc.mul(w).add(d[l]);
                }
                // direct near field
                for ny in iy.saturating_sub(1)..=(iy + 1).min(nb - 1) {
                    for nx in ix.saturating_sub(1)..=(ix + 1).min(nb - 1) {
                        for &j in &leaf[ny * nb + nx] {
                            if j == i {
                                continue;
                            }
                            let dz = z.sub(Cx::new(pos[j][0], pos[j][1]));
                            acc = acc.add(dz.inv().scale(src[j]));
                        }
                    }
                }
                // field f = -conj(W)
                out[i] = [-acc.re, acc.im];
            }
        }
    }
    Ok(out)
}
