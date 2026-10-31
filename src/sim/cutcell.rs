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
//! The clipped measure of a simplex with vertex values `f_i` is
//! `Σ_i max(f_i, 0)^d / Π_{j≠i} (f_i - f_j)` (fraction of the simplex), which
//! is used for triangles and tetrahedra alike.
//!
//! These are the apertures and volume fractions needed by cut-cell finite
//! volume solvers such as [`ns_mac`](crate::sim::ns_mac).

/// Level-set function: positive in the open region.
pub type LevelSet<'a> = &'a dyn Fn([f64; 3]) -> f64;

/// Fraction of a simplex with vertex values `f` (length `d + 1`) where the
/// linearly interpolated function is positive.
#[must_use]
pub fn simplex_positive_fraction(f: &[f64]) -> f64 {
    let d = f.len() - 1;
    let pos: Vec<usize> = (0..=d).filter(|&i| f[i] > 0.0).collect();
    let neg: Vec<usize> = (0..=d).filter(|&i| f[i] <= 0.0).collect();
    if neg.is_empty() {
        return 1.0;
    }
    if pos.is_empty() {
        return 0.0;
    }
    // fraction of the simplex cut off at the vertex `v` against the opposite-sign vertices `others`
    let corner = |v: usize, others: &[usize]| -> f64 {
        others.iter().map(|&o| f[v] / (f[v] - f[o])).product()
    };
    if pos.len() == 1 {
        return corner(pos[0], &neg);
    }
    if neg.len() == 1 {
        return 1.0 - corner(neg[0], &pos);
    }
    // tetrahedron with two positive vertices a, b and two non-positive vertices c, e:
    // the positive part is the union of three tetrahedra on the cut points
    let (a, b, c, e) = (pos[0], pos[1], neg[0], neg[1]);
    let cut = |p: usize, q: usize| -> [f64; 4] {
        // barycentric coordinates of the zero crossing on the edge p-q
        let t = f[p] / (f[p] - f[q]);
        let mut w = [0.0; 4];
        w[p] = 1.0 - t;
        w[q] = t;
        w
    };
    let unit = |i: usize| -> [f64; 4] {
        let mut w = [0.0; 4];
        w[i] = 1.0;
        w
    };
    let (va, vb) = (unit(a), unit(b));
    let (pac, pae, pbc, pbe) = (cut(a, c), cut(a, e), cut(b, c), cut(b, e));
    let det4 = |r: [[f64; 4]; 4]| -> f64 {
        // determinant by cofactor expansion
        let det3 = |m: [[f64; 3]; 3]| {
            m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
                - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
        };
        let mut s = 0.0;
        for col in 0..4 {
            let mut m = [[0.0; 3]; 3];
            for i in 1..4 {
                let mut cc = 0;
                for j in 0..4 {
                    if j != col {
                        m[i - 1][cc] = r[i][j];
                        cc += 1;
                    }
                }
            }
            s += if col % 2 == 0 { 1.0 } else { -1.0 } * r[0][col] * det3(m);
        }
        s.abs()
    };
    det4([va, vb, pac, pae]) + det4([vb, pae, pac, pbe]) + det4([vb, pac, pbc, pbe])
}

/// Open fraction of a segment `x0 + t e_axis len`, `t in [0, 1]`.
#[must_use]
pub fn open_length_fraction(
    phi: LevelSet<'_>,
    x0: [f64; 3],
    axis: usize,
    len: f64,
) -> f64 {
    const N: usize = 16;
    let at = |t: f64| {
        let mut x = x0;
        x[axis] += t * len;
        phi(x)
    };
    let vals: Vec<f64> = (0..=N).map(|i| at(i as f64 / N as f64)).collect();
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
                if (at(m) > 0.0) == p0 {
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

fn corners(
    lo: [f64; 3],
    len: [f64; 3],
    axes: &[usize],
) -> Vec<[f64; 3]> {
    let nc = 1usize << axes.len();
    (0..nc)
        .map(|m| {
            let mut x = lo;
            for (k, &a) in axes.iter().enumerate() {
                if m >> k & 1 == 1 {
                    x[a] += len[a];
                }
            }
            x
        })
        .collect()
}

fn recurse(
    phi: LevelSet<'_>,
    lo: [f64; 3],
    len: [f64; 3],
    axes: &[usize],
    depth: usize,
) -> f64 {
    let cs = corners(lo, len, axes);
    let f: Vec<f64> = cs.iter().map(|&x| phi(x)).collect();
    // a smooth body can sit entirely inside a box (all corners open): probe the centre
    let mut centre = lo;
    for &a in axes {
        centre[a] += 0.5 * len[a];
    }
    let fc = phi(centre);
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
        let nc = 1usize << axes.len();
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
#[must_use]
pub fn open_fraction(
    phi: LevelSet<'_>,
    lo: [f64; 3],
    len: [f64; 3],
    axes: &[usize],
    depth: usize,
) -> f64 {
    match axes.len() {
        | 0 => f64::from(u8::from(phi(lo) > 0.0)),
        | 1 => open_length_fraction(phi, lo, axes[0], len[axes[0]]),
        | _ => recurse(phi, lo, len, axes, depth),
    }
}

/// Default refinement depth for faces and cells.
#[must_use]
pub const fn default_depth(naxes: usize) -> usize {
    if naxes >= 3 { 3 } else { 4 }
}
