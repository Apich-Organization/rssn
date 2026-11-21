//! Cut-cell geometry (open fractions of segments, faces and cells) from a level set.
//!
//! # Cut-cell geometry from a level set
//!
//! Open (fluid) fractions of grid segments, faces and cells that are cut by a
//! body given as a level set `φ` (`φ > 0` is the open region, `φ <= 0` the
//! body). All routines are dimension agnostic: a sub-box is described by its
//! corner `lo`, its edge lengths `len` and the list of axes it spans
//! (1 axis: a segment, 2 axes: a face or a 2D cell, 3 axes: a 3D cell).
//!
//! * segments: every sign change of `φ` on a fine sampling of the segment is
//!   located by bisection on the true level set (exact for any smooth body up
//!   to the root tolerance);
//! * faces: the square is subdivided adaptively where the corner signs
//!   differ; on the finest sub-squares `φ` is taken linear on each of the four
//!   triangles that join the corners with the centre and the positive area is
//!   clipped **exactly** (planar interfaces are reproduced to round-off);
//! * cells: the same with the six Kuhn tetrahedra of each sub-cube.
//!
//! The clipped measure is assembled from the edge crossings between positive
//! and nonpositive vertices. Repeated and zero vertex values are allowed.
//!
//! These are the apertures and volume fractions needed by cut-cell finite
//! volume solvers such as [`ns_mac`](crate::sim::ns_mac).

/// Level-set function: positive in the open region.
pub type LevelSet<'a> = &'a dyn Fn([f64; 3]) -> f64;

/// Fraction of a point, segment, triangle or tetrahedron where the linearly
/// interpolated vertex values are positive.
///
/// Returns NaN for any other vertex count or a nonfinite value.
#[must_use]
pub fn simplex_positive_fraction(f: &[f64]) -> f64 {
    if !(1..=4).contains(&f.len()) {
        return f64::NAN;
    }
    let (mut pos, mut neg) = ([0.0; 4], [0.0; 4]);
    let (mut npos, mut nneg) = (0, 0);
    for &value in f {
        if !value.is_finite() {
            return f64::NAN;
        }
        if value > 0.0 {
            pos[npos] = value;
            npos += 1;
        } else {
            neg[nneg] = value.abs();
            nneg += 1;
        }
    }
    let positive = &pos[..npos];
    let negative = &neg[..nneg];
    match (positive, negative) {
        | ([], _) => 0.0,
        | (_, []) => 1.0,
        | ([p], others) => others.iter().map(|&q| edge_fractions(*p, q).0).product(),
        | (others, [q]) => {
            // Expand the complement of the closed corner as a sum of
            // nonnegative terms. Subtracting its measure from one can erase
            // a thin open region even when the result is representable.
            let (mut open, mut closed) = (0.0, 1.0);
            for &p in others {
                let (open_edge, closed_edge) = edge_fractions(p, *q);
                open += closed * open_edge;
                closed *= closed_edge;
            }
            open.min(1.0)
        },
        | ([a, b], [c, d]) => {
            let (ac, ca) = edge_fractions(*a, *c);
            let (ad, da) = edge_fractions(*a, *d);
            let (bc, _) = edge_fractions(*b, *c);
            let (bd, _) = edge_fractions(*b, *d);
            // The open portion splits into three tetrahedra. Expanding their
            // barycentric determinants gives these positive products.
            (ac * ad + da * ac * bd + ca * bc * bd).min(1.0)
        },
        | _ => unreachable!("at most four vertices were classified"),
    }
}

fn edge_fractions(
    p: f64,
    q: f64,
) -> (f64, f64) {
    let sum = p + q;
    if sum.is_finite() {
        (p / sum, q / sum)
    } else {
        let (p, q) = (p * 0.5, q * 0.5);
        let sum = p + q;
        (p / sum, q / sum)
    }
}

/// Open fraction of a segment `x0 + t e_axis len`, `t in [0, 1]`.
///
/// Returns NaN for invalid geometry or nonfinite level-set samples.
#[must_use]
pub fn open_length_fraction(
    phi: LevelSet<'_>,
    x0: [f64; 3],
    axis: usize,
    len: f64,
) -> f64 {
    const N: usize = 16;
    if axis >= 3 || x0.iter().any(|x| !x.is_finite()) || !(x0[axis] + len).is_finite() {
        return f64::NAN;
    }
    let at = |t: f64| {
        let mut x = x0;
        x[axis] += t * len;
        phi(x)
    };
    let vals: [f64; N + 1] = std::array::from_fn(|i| at(i as f64 / N as f64));
    if vals.iter().any(|value| !value.is_finite()) {
        return f64::NAN;
    }
    let mut open = 0.0;
    for i in 0..N {
        let (t0, t1) = (i as f64 / N as f64, (i + 1) as f64 / N as f64);
        let (f0, f1) = (vals[i], vals[i + 1]);
        let (p0, p1) = (f0 > 0.0, f1 > 0.0);
        if p0 && p1 {
            open += t1 - t0;
        } else if p0 != p1 {
            // bisection on the true level set
            let (mut a, mut b) = (t0, t1);
            for _ in 0..60 {
                let m = f64::midpoint(a, b);
                if crate::kernels::exact_eq(m, a) || crate::kernels::exact_eq(m, b) {
                    break;
                }
                let fm = at(m);
                if !fm.is_finite() {
                    return f64::NAN;
                }
                if (fm > 0.0) == p0 {
                    a = m;
                } else {
                    b = m;
                }
            }
            let root = f64::midpoint(a, b);
            open += if p0 { root - t0 } else { t1 - root };
        }
    }
    open
}

fn recurse(
    phi: LevelSet<'_>,
    lo: [f64; 3],
    len: [f64; 3],
    axes: &[usize],
    depth: usize,
) -> f64 {
    let nc = 1usize << axes.len();
    let mut values = [0.0; 8];
    let f = &mut values[..nc];
    for (m, value) in f.iter_mut().enumerate() {
        let mut x = lo;
        for (k, &a) in axes.iter().enumerate() {
            if m >> k & 1 == 1 {
                x[a] += len[a];
            }
        }
        *value = phi(x);
    }
    // a smooth body can sit entirely inside a box (all corners open): probe the centre
    let mut centre = lo;
    for &a in axes {
        centre[a] += 0.5 * len[a];
    }
    let fc = phi(centre);
    if !fc.is_finite() || f.iter().any(|x| !x.is_finite()) {
        return f64::NAN;
    }
    if f.iter().all(|&v| v > 0.0) && fc > 0.0 {
        return 1.0;
    }
    if f.iter().all(|&v| v <= 0.0) && fc <= 0.0 {
        return 0.0;
    }
    if depth > 0 {
        let mut half = len;
        for &a in axes {
            half[a] *= 0.5;
        }
        let mut s = 0.0;
        for m in 0..nc {
            let mut l = lo;
            for (k, &a) in axes.iter().enumerate() {
                if m >> k & 1 == 1 {
                    l[a] += half[a];
                }
            }
            s += recurse(phi, l, half, axes, depth - 1);
        }
        return s / nc as f64;
    }
    // leaf: exact clipping of the piecewise-linear interpolant
    match axes.len() {
        | 1 => simplex_positive_fraction(&[f[0], f[1]]),
        | 2 => {
            // corners: 0=(0,0) 1=(1,0) 2=(0,1) 3=(1,1); fan around the centre
            let tris = [[0, 1], [1, 3], [3, 2], [2, 0]];
            tris.iter()
                .map(|t| simplex_positive_fraction(&[f[t[0]], f[t[1]], fc]))
                .sum::<f64>()
                / 4.0
        },
        | _ => {
            // Kuhn decomposition along the diagonal 0 -> 7
            const PATHS: [[usize; 4]; 6] = [
                [0, 1, 3, 7],
                [0, 1, 5, 7],
                [0, 2, 3, 7],
                [0, 2, 6, 7],
                [0, 4, 5, 7],
                [0, 4, 6, 7],
            ];
            PATHS
                .iter()
                .map(|p| simplex_positive_fraction(&[f[p[0]], f[p[1]], f[p[2]], f[p[3]]]))
                .sum::<f64>()
                / 6.0
        },
    }
}

/// Open fraction of the axis-aligned box that starts at `lo`, has edge lengths
/// `len` along the listed `axes` (one to three of them) and is sampled with
/// `depth` levels of adaptive refinement.
///
/// One axis gives the length fraction of a segment, two axes the area
/// fraction of a face (or of a 2D cell) and three axes the volume fraction of
/// a 3D cell. Boxes whose corners and centre are all open (or all closed) are
/// returned immediately.
///
/// Returns NaN for repeated or invalid axes, nonfinite geometry, or a
/// nonfinite level-set sample.
#[must_use]
pub fn open_fraction(
    phi: LevelSet<'_>,
    lo: [f64; 3],
    len: [f64; 3],
    axes: &[usize],
    depth: usize,
) -> f64 {
    let mut used = 0_u8;
    for &axis in axes {
        if axis >= 3 || used & (1 << axis) != 0 || !(lo[axis] + len[axis]).is_finite() {
            return f64::NAN;
        }
        used |= 1 << axis;
    }
    if lo.iter().any(|x| !x.is_finite()) {
        return f64::NAN;
    }
    match axes.len() {
        | 0 => simplex_positive_fraction(&[phi(lo)]),
        | 1 => open_length_fraction(phi, lo, axes[0], len[axes[0]]),
        | _ => recurse(phi, lo, len, axes, depth),
    }
}

/// Default refinement depth for faces and cells.
#[must_use]
pub const fn default_depth(naxes: usize) -> usize {
    if naxes >= 3 { 3 } else { 4 }
}
