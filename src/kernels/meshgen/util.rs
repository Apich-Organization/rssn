//! Shared geometry for the mesh-generation kernels.
//!
//! Adaptive-precision geometric predicates (Shewchuk's exact floating-point
//! expansion arithmetic with fast `f64` error-bound filters), deterministic
//! Simulation of Simplicity (SoS) tie-breakers for cocircular and cospherical
//! inputs, circumcircle geometry, triangle quality measures, polygon helpers
//! and the [`TriMesh2`] container used by [`super::delaunay`],
//! [`super::refine`], [`super::advancing_front`] and [`super::tetra`].

/// Half the IEEE-754 binary64 unit in the last place (`2^-53`), matching
/// Shewchuk's `epsilon`.
const EPS: f64 = f64::EPSILON * 0.5;

/// Error bound for the fast `f64` filter of [`orient2d`].
const CCW_ERRBOUND_A: f64 = (3.0 + 16.0 * EPS) * EPS;

/// Error bound for the fast `f64` filter of [`orient3d`].
const O3D_ERRBOUND_A: f64 = (7.0 + 56.0 * EPS) * EPS;

/// Error bound for the fast `f64` filter of [`incircle`].
const ICC_ERRBOUND_A: f64 = (10.0 + 96.0 * EPS) * EPS;

/// Error bound for the fast `f64` filter of [`insphere`].
const ISP_ERRBOUND_A: f64 = (16.0 + 224.0 * EPS) * EPS;

/// An unstructured triangle mesh in the plane.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct TriMesh2 {
    /// Vertex coordinates, one `[x, y]` per vertex.
    pub points: Vec<[f64; 2]>,
    /// Counter-clockwise triangles, given as triples of vertex indices.
    pub triangles: Vec<[usize; 3]>,
}

// ---------------------------------------------------------------------------
// Shewchuk's exact adaptive-precision floating-point expansion arithmetic
// ---------------------------------------------------------------------------

/// Error-free transformation of `a + b` assuming `|a| >= |b|`.
///
/// Returns `(sum, err)` such that `sum + err == a + b` exactly.
const fn fast_two_sum(
    a: f64,
    b: f64,
) -> (f64, f64) {
    let x = a + b;
    let b_virtual = x - a;
    let y = b - b_virtual;
    (x, y)
}

/// Error-free transformation of `a + b` for arbitrary `a` and `b`.
///
/// Returns `(sum, err)` such that `sum + err == a + b` exactly.
const fn two_sum(
    a: f64,
    b: f64,
) -> (f64, f64) {
    let x = a + b;
    let b_virtual = x - a;
    let a_virtual = x - b_virtual;
    let b_roundoff = b - b_virtual;
    let a_roundoff = a - a_virtual;
    let y = a_roundoff + b_roundoff;
    (x, y)
}

/// Error-free transformation of `a - b` for arbitrary `a` and `b`.
///
/// Returns `(diff, err)` such that `diff + err == a - b` exactly.
const fn two_diff(
    a: f64,
    b: f64,
) -> (f64, f64) {
    let x = a - b;
    let b_virtual = a - x;
    let a_virtual = x + b_virtual;
    let b_roundoff = b_virtual - b;
    let a_roundoff = a - a_virtual;
    let y = a_roundoff + b_roundoff;
    (x, y)
}

/// Exact two-component expansion `[lo, hi]` representing `a - b`.
const fn two_diff_exp(
    a: f64,
    b: f64,
) -> [f64; 2] {
    let (hi, lo) = two_diff(a, b);
    [lo, hi]
}

/// Error-free transformation of `a * b` using IEEE-754 fused multiply-add.
///
/// Returns `(prod, err)` such that `prod + err == a * b` exactly.
fn two_product(
    a: f64,
    b: f64,
) -> (f64, f64) {
    let x = a * b;
    let y = a.mul_add(b, -x);
    (x, y)
}

/// Adds a scalar `b` to the non-overlapping expansion `e` (least-significant
/// component first), eliminating zero components.
fn grow_expansion(
    e: &[f64],
    b: f64,
) -> Vec<f64> {
    let mut h = Vec::with_capacity(e.len() + 1);
    let mut q = b;
    for &enow in e {
        let (q_new, hh) = two_sum(q, enow);
        q = q_new;
        if hh != 0.0 {
            h.push(hh);
        }
    }
    if q != 0.0 || h.is_empty() {
        h.push(q);
    }
    h
}

/// Multiplies the non-overlapping expansion `e` by a scalar `b`, eliminating
/// zero components.
fn scale_expansion(
    e: &[f64],
    b: f64,
) -> Vec<f64> {
    if e.is_empty() || b == 0.0 {
        return vec![0.0];
    }
    let mut h = Vec::with_capacity(2 * e.len());
    let (mut q, hh) = two_product(e[0], b);
    if hh != 0.0 {
        h.push(hh);
    }
    for &enow in &e[1..] {
        let (p1, p0) = two_product(enow, b);
        let (sum, hh1) = two_sum(q, p0);
        if hh1 != 0.0 {
            h.push(hh1);
        }
        let (q_new, hh2) = fast_two_sum(p1, sum);
        q = q_new;
        if hh2 != 0.0 {
            h.push(hh2);
        }
    }
    if q != 0.0 || h.is_empty() {
        h.push(q);
    }
    h
}

/// Sums two non-overlapping expansions `e` and `f` in linear time (Shewchuk's
/// `fast_expansion_sum_zeroelim`), eliminating zero components.
fn expansion_sum(
    e: &[f64],
    f: &[f64],
) -> Vec<f64> {
    if e.is_empty() {
        return f.to_vec();
    }
    if f.is_empty() {
        return e.to_vec();
    }
    if f.len() == 1 {
        return grow_expansion(e, f[0]);
    }
    if e.len() == 1 {
        return grow_expansion(f, e[0]);
    }

    let mut h = Vec::with_capacity(e.len() + f.len());
    let mut ei = 0_usize;
    let mut fi = 0_usize;
    let mut enow = e[0];
    let mut fnow = f[0];

    let mut q = if fnow.abs() > enow.abs() {
        ei += 1;
        if ei < e.len() {
            enow = e[ei];
        }
        e[0]
    } else {
        fi += 1;
        if fi < f.len() {
            fnow = f[fi];
        }
        f[0]
    };

    if ei < e.len() && fi < f.len() {
        let (q_new, hh) = if fnow.abs() > enow.abs() {
            let res = fast_two_sum(enow, q);
            ei += 1;
            if ei < e.len() {
                enow = e[ei];
            }
            res
        } else {
            let res = fast_two_sum(fnow, q);
            fi += 1;
            if fi < f.len() {
                fnow = f[fi];
            }
            res
        };
        q = q_new;
        if hh != 0.0 {
            h.push(hh);
        }

        while ei < e.len() && fi < f.len() {
            let (q_new, hh) = if fnow.abs() > enow.abs() {
                let res = two_sum(q, enow);
                ei += 1;
                if ei < e.len() {
                    enow = e[ei];
                }
                res
            } else {
                let res = two_sum(q, fnow);
                fi += 1;
                if fi < f.len() {
                    fnow = f[fi];
                }
                res
            };
            q = q_new;
            if hh != 0.0 {
                h.push(hh);
            }
        }
    }

    while ei < e.len() {
        let (q_new, hh) = two_sum(q, enow);
        q = q_new;
        if hh != 0.0 {
            h.push(hh);
        }
        ei += 1;
        if ei < e.len() {
            enow = e[ei];
        }
    }

    while fi < f.len() {
        let (q_new, hh) = two_sum(q, fnow);
        q = q_new;
        if hh != 0.0 {
            h.push(hh);
        }
        fi += 1;
        if fi < f.len() {
            fnow = f[fi];
        }
    }

    if q != 0.0 || h.is_empty() {
        h.push(q);
    }
    h
}

/// Computes the exact difference `e - f` of two non-overlapping expansions.
fn expansion_diff(
    e: &[f64],
    f: &[f64],
) -> Vec<f64> {
    let neg_f: Vec<f64> = f.iter().map(|&x| -x).collect();
    expansion_sum(e, &neg_f)
}

/// Computes the exact product `e * f` of two non-overlapping expansions.
fn expansion_product(
    e: &[f64],
    f: &[f64],
) -> Vec<f64> {
    let mut acc = vec![0.0];
    for &fi in f {
        if fi != 0.0 {
            let term = scale_expansion(e, fi);
            acc = expansion_sum(&acc, &term);
        }
    }
    acc
}

/// Computes a one-ulp estimate of the expansion `e` (preserving the exact sign
/// for any zero-eliminated non-overlapping expansion).
fn estimate(e: &[f64]) -> f64 {
    let mut q = 0.0;
    for &v in e {
        q += v;
    }
    q
}

// ---------------------------------------------------------------------------
// Adaptive-precision geometric predicates
// ---------------------------------------------------------------------------

/// Twice the signed area of the triangle `a`, `b`, `c`.
///
/// The value is positive when `a`, `b`, `c` wind counter-clockwise, negative
/// when they wind clockwise, and `0.0` when the three points are collinear.
/// Uses a fast `f64` filter with Shewchuk's exact expansion fallback.
pub(super) fn orient2d(
    a: [f64; 2],
    b: [f64; 2],
    c: [f64; 2],
) -> f64 {
    let detleft = (a[0] - c[0]) * (b[1] - c[1]);
    let detright = (a[1] - c[1]) * (b[0] - c[0]);
    let det = detleft - detright;

    let detsum = if detleft > 0.0 {
        if detright <= 0.0 {
            return det;
        }
        detleft + detright
    } else if detleft < 0.0 {
        if detright >= 0.0 {
            return det;
        }
        -detleft - detright
    } else {
        return det;
    };

    if det.abs() >= CCW_ERRBOUND_A * detsum {
        return det;
    }

    orient2d_exact(a, b, c)
}

/// Exact expansion evaluation of [`orient2d`].
fn orient2d_exact(
    a: [f64; 2],
    b: [f64; 2],
    c: [f64; 2],
) -> f64 {
    let acx = two_diff_exp(a[0], c[0]);
    let acy = two_diff_exp(a[1], c[1]);
    let bcx = two_diff_exp(b[0], c[0]);
    let bcy = two_diff_exp(b[1], c[1]);

    let left = expansion_product(&acx, &bcy);
    let right = expansion_product(&acy, &bcx);
    let diff = expansion_diff(&left, &right);
    estimate(&diff)
}

/// Six times the signed volume of the tetrahedron `a`, `b`, `c`, `d`.
///
/// Equals `det([b - a; c - a; d - a])`, which is positive for the standard
/// right-handed simplex `(0,0,0), (1,0,0), (0,1,0), (0,0,1)`. Uses a fast
/// `f64` filter with Shewchuk's exact expansion fallback.
pub(super) fn orient3d(
    a: [f64; 3],
    b: [f64; 3],
    c: [f64; 3],
    d: [f64; 3],
) -> f64 {
    let ux = b[0] - a[0];
    let uy = b[1] - a[1];
    let uz = b[2] - a[2];
    let vx = c[0] - a[0];
    let vy = c[1] - a[1];
    let vz = c[2] - a[2];
    let wx = d[0] - a[0];
    let wy = d[1] - a[1];
    let wz = d[2] - a[2];

    let vywz = vy * wz;
    let vzwy = vz * wy;
    let vzwx = vz * wx;
    let vxwz = vx * wz;
    let vxwy = vx * wy;
    let vywx = vy * wx;

    let det = ux * (vywz - vzwy) + uy * (vzwx - vxwz) + uz * (vxwy - vywx);
    let permanent = (vywz.abs() + vzwy.abs()) * ux.abs()
        + (vzwx.abs() + vxwz.abs()) * uy.abs()
        + (vxwy.abs() + vywx.abs()) * uz.abs();

    if det.abs() > O3D_ERRBOUND_A * permanent {
        return det;
    }

    orient3d_exact(a, b, c, d)
}

/// Exact expansion evaluation of [`orient3d`].
fn orient3d_exact(
    a: [f64; 3],
    b: [f64; 3],
    c: [f64; 3],
    d: [f64; 3],
) -> f64 {
    let ux = two_diff_exp(b[0], a[0]);
    let uy = two_diff_exp(b[1], a[1]);
    let uz = two_diff_exp(b[2], a[2]);
    let vx = two_diff_exp(c[0], a[0]);
    let vy = two_diff_exp(c[1], a[1]);
    let vz = two_diff_exp(c[2], a[2]);
    let wx = two_diff_exp(d[0], a[0]);
    let wy = two_diff_exp(d[1], a[1]);
    let wz = two_diff_exp(d[2], a[2]);

    let m0 = expansion_diff(&expansion_product(&vy, &wz), &expansion_product(&vz, &wy));
    let m1 = expansion_diff(&expansion_product(&vz, &wx), &expansion_product(&vx, &wz));
    let m2 = expansion_diff(&expansion_product(&vx, &wy), &expansion_product(&vy, &wx));

    let t0 = expansion_product(&ux, &m0);
    let t1 = expansion_product(&uy, &m1);
    let t2 = expansion_product(&uz, &m2);
    let sum = expansion_sum(&expansion_sum(&t0, &t1), &t2);
    estimate(&sum)
}

/// The in-circle determinant of `d` against the counter-clockwise triangle
/// `a`, `b`, `c`.
///
/// The value is positive when `d` lies strictly inside the circumcircle of the
/// triangle, zero when it lies on it and negative when it lies outside.
/// Uses a fast `f64` filter with Shewchuk's exact expansion fallback.
pub(super) fn incircle(
    a: [f64; 2],
    b: [f64; 2],
    c: [f64; 2],
    d: [f64; 2],
) -> f64 {
    let adx = a[0] - d[0];
    let ady = a[1] - d[1];
    let bdx = b[0] - d[0];
    let bdy = b[1] - d[1];
    let cdx = c[0] - d[0];
    let cdy = c[1] - d[1];

    let bdxcdy = bdx * cdy;
    let cdxbdy = cdx * bdy;
    let alift = adx * adx + ady * ady;

    let cdxady = cdx * ady;
    let adxcdy = adx * cdy;
    let blift = bdx * bdx + bdy * bdy;

    let adxbdy = adx * bdy;
    let bdxady = bdx * ady;
    let clift = cdx * cdx + cdy * cdy;

    let det = alift * (bdxcdy - cdxbdy) + blift * (cdxady - adxcdy) + clift * (adxbdy - bdxady);

    let permanent = (bdxcdy.abs() + cdxbdy.abs()) * alift
        + (cdxady.abs() + adxcdy.abs()) * blift
        + (adxbdy.abs() + bdxady.abs()) * clift;

    if det.abs() > ICC_ERRBOUND_A * permanent {
        return det;
    }

    incircle_exact(a, b, c, d)
}

/// Exact expansion evaluation of [`incircle`].
fn incircle_exact(
    a: [f64; 2],
    b: [f64; 2],
    c: [f64; 2],
    d: [f64; 2],
) -> f64 {
    let adx = two_diff_exp(a[0], d[0]);
    let ady = two_diff_exp(a[1], d[1]);
    let bdx = two_diff_exp(b[0], d[0]);
    let bdy = two_diff_exp(b[1], d[1]);
    let cdx = two_diff_exp(c[0], d[0]);
    let cdy = two_diff_exp(c[1], d[1]);

    let alift = expansion_sum(
        &expansion_product(&adx, &adx),
        &expansion_product(&ady, &ady),
    );
    let blift = expansion_sum(
        &expansion_product(&bdx, &bdx),
        &expansion_product(&bdy, &bdy),
    );
    let clift = expansion_sum(
        &expansion_product(&cdx, &cdx),
        &expansion_product(&cdy, &cdy),
    );

    let bc = expansion_diff(
        &expansion_product(&bdx, &cdy),
        &expansion_product(&cdx, &bdy),
    );
    let ca = expansion_diff(
        &expansion_product(&cdx, &ady),
        &expansion_product(&adx, &cdy),
    );
    let ab = expansion_diff(
        &expansion_product(&adx, &bdy),
        &expansion_product(&bdx, &ady),
    );

    let t1 = expansion_product(&alift, &bc);
    let t2 = expansion_product(&blift, &ca);
    let t3 = expansion_product(&clift, &ab);
    let total = expansion_sum(&expansion_sum(&t1, &t2), &t3);
    estimate(&total)
}

/// The in-sphere determinant of `e` against the positively oriented
/// tetrahedron `a`, `b`, `c`, `d`.
///
/// The value is positive when `e` lies strictly inside the circumsphere of
/// `a`, `b`, `c`, `d`, zero when it lies on it, and negative when it lies
/// outside. Uses a fast `f64` filter with Shewchuk's exact expansion fallback.
pub(super) fn insphere(
    a: [f64; 3],
    b: [f64; 3],
    c: [f64; 3],
    d: [f64; 3],
    e: [f64; 3],
) -> f64 {
    let aex = a[0] - e[0];
    let aey = a[1] - e[1];
    let aez = a[2] - e[2];
    let bex = b[0] - e[0];
    let bey = b[1] - e[1];
    let bez = b[2] - e[2];
    let cex = c[0] - e[0];
    let cey = c[1] - e[1];
    let cez = c[2] - e[2];
    let dex = d[0] - e[0];
    let dey = d[1] - e[1];
    let dez = d[2] - e[2];

    let aexbey = aex * bey;
    let bexaey = bex * aey;
    let ab = aexbey - bexaey;
    let bexcey = bex * cey;
    let cexbey = cex * bey;
    let bc = bexcey - cexbey;
    let cexdey = cex * dey;
    let dexcey = dex * cey;
    let cd = cexdey - dexcey;
    let dexaey = dex * aey;
    let aexdey = aex * dey;
    let da = dexaey - aexdey;
    let aexcey = aex * cey;
    let cexaey = cex * aey;
    let ac = aexcey - cexaey;
    let bexdey = bex * dey;
    let dexbey = dex * bey;
    let bd = bexdey - dexbey;

    let abc = aez * bc - bez * ac + cez * ab;
    let bcd = bez * cd - cez * bd + dez * bc;
    let cda = cez * da + dez * ac + aez * cd;
    let dab = dez * ab + aez * bd + bez * da;

    let alift = aex * aex + aey * aey + aez * aez;
    let blift = bex * bex + bey * bey + bez * bez;
    let clift = cex * cex + cey * cey + cez * cez;
    let dlift = dex * dex + dey * dey + dez * dez;

    let det = alift * bcd - blift * cda + clift * dab - dlift * abc;

    let permanent = ((cexdey.abs() + dexcey.abs()) * bez.abs()
        + (dexbey.abs() + bexdey.abs()) * cez.abs()
        + (bexcey.abs() + cexbey.abs()) * dez.abs())
        * alift
        + ((dexaey.abs() + aexdey.abs()) * cez.abs()
            + (aexcey.abs() + cexaey.abs()) * dez.abs()
            + (cexdey.abs() + dexcey.abs()) * aez.abs())
            * blift
        + ((aexbey.abs() + bexaey.abs()) * dez.abs()
            + (bexdey.abs() + dexbey.abs()) * aez.abs()
            + (dexaey.abs() + aexdey.abs()) * bez.abs())
            * clift
        + ((bexcey.abs() + cexbey.abs()) * aez.abs()
            + (cexaey.abs() + aexcey.abs()) * bez.abs()
            + (aexbey.abs() + bexaey.abs()) * cez.abs())
            * dlift;

    if det.abs() > ISP_ERRBOUND_A * permanent {
        return det;
    }

    insphere_exact(a, b, c, d, e)
}

/// Exact expansion evaluation of [`insphere`].
fn insphere_exact(
    a: [f64; 3],
    b: [f64; 3],
    c: [f64; 3],
    d: [f64; 3],
    e: [f64; 3],
) -> f64 {
    let aex = two_diff_exp(a[0], e[0]);
    let aey = two_diff_exp(a[1], e[1]);
    let aez = two_diff_exp(a[2], e[2]);
    let bex = two_diff_exp(b[0], e[0]);
    let bey = two_diff_exp(b[1], e[1]);
    let bez = two_diff_exp(b[2], e[2]);
    let cex = two_diff_exp(c[0], e[0]);
    let cey = two_diff_exp(c[1], e[1]);
    let cez = two_diff_exp(c[2], e[2]);
    let dex = two_diff_exp(d[0], e[0]);
    let dey = two_diff_exp(d[1], e[1]);
    let dez = two_diff_exp(d[2], e[2]);

    let ab = expansion_diff(
        &expansion_product(&aex, &bey),
        &expansion_product(&bex, &aey),
    );
    let bc = expansion_diff(
        &expansion_product(&bex, &cey),
        &expansion_product(&cex, &bey),
    );
    let cd = expansion_diff(
        &expansion_product(&cex, &dey),
        &expansion_product(&dex, &cey),
    );
    let da = expansion_diff(
        &expansion_product(&dex, &aey),
        &expansion_product(&aex, &dey),
    );
    let ac = expansion_diff(
        &expansion_product(&aex, &cey),
        &expansion_product(&cex, &aey),
    );
    let bd = expansion_diff(
        &expansion_product(&bex, &dey),
        &expansion_product(&dex, &bey),
    );

    let abc = expansion_sum(
        &expansion_diff(&expansion_product(&aez, &bc), &expansion_product(&bez, &ac)),
        &expansion_product(&cez, &ab),
    );
    let bcd = expansion_sum(
        &expansion_diff(&expansion_product(&bez, &cd), &expansion_product(&cez, &bd)),
        &expansion_product(&dez, &bc),
    );
    let cda = expansion_sum(
        &expansion_sum(&expansion_product(&cez, &da), &expansion_product(&dez, &ac)),
        &expansion_product(&aez, &cd),
    );
    let dab = expansion_sum(
        &expansion_sum(&expansion_product(&dez, &ab), &expansion_product(&aez, &bd)),
        &expansion_product(&bez, &da),
    );

    let alift = expansion_sum(
        &expansion_sum(
            &expansion_product(&aex, &aex),
            &expansion_product(&aey, &aey),
        ),
        &expansion_product(&aez, &aez),
    );
    let blift = expansion_sum(
        &expansion_sum(
            &expansion_product(&bex, &bex),
            &expansion_product(&bey, &bey),
        ),
        &expansion_product(&bez, &bez),
    );
    let clift = expansion_sum(
        &expansion_sum(
            &expansion_product(&cex, &cex),
            &expansion_product(&cey, &cey),
        ),
        &expansion_product(&cez, &cez),
    );
    let dlift = expansion_sum(
        &expansion_sum(
            &expansion_product(&dex, &dex),
            &expansion_product(&dey, &dey),
        ),
        &expansion_product(&dez, &dez),
    );

    let t1 = expansion_product(&alift, &bcd);
    let t2 = expansion_product(&blift, &cda);
    let t3 = expansion_product(&clift, &dab);
    let t4 = expansion_product(&dlift, &abc);
    let total = expansion_diff(&expansion_sum(&expansion_diff(&t1, &t2), &t3), &t4);
    estimate(&total)
}

/// Sorts a small slice in place and returns `+1` for an even permutation or
/// `-1` for an odd permutation.
fn sort_with_parity(slice: &mut [usize]) -> i8 {
    let mut sign = 1_i8;
    for i in 1..slice.len() {
        let mut j = i;
        while j > 0 && slice[j - 1] > slice[j] {
            slice.swap(j - 1, j);
            sign = -sign;
            j -= 1;
        }
    }
    sign
}

/// Exact in-circle predicate with deterministic Simulation of Simplicity (SoS)
/// tie-breaking for cocircular points.
///
/// Perturbs each vertex's paraboloid lifting coordinate by `+eps^(2^index)`
/// (smallest index lifted most). Returns `true` when `points[d]` lies inside
/// the circumcircle of the counter-clockwise triangle `(a, b, c)`.
pub(super) fn incircle_sos(
    points: &[[f64; 2]],
    a: usize,
    b: usize,
    c: usize,
    d: usize,
) -> bool {
    if a == b || a == c || a == d || b == c || b == d || c == d {
        return false;
    }
    let val = incircle(points[a], points[b], points[c], points[d]);
    if val > 0.0 {
        return true;
    }
    if val < 0.0 {
        return false;
    }

    let mut idx = [a, b, c, d];
    let parity = sort_with_parity(&mut idx);

    let minors = [
        (
            1_i8,
            orient2d(points[idx[1]], points[idx[2]], points[idx[3]]),
        ),
        (
            -1_i8,
            orient2d(points[idx[0]], points[idx[2]], points[idx[3]]),
        ),
        (
            1_i8,
            orient2d(points[idx[0]], points[idx[1]], points[idx[3]]),
        ),
        (
            -1_i8,
            orient2d(points[idx[0]], points[idx[1]], points[idx[2]]),
        ),
    ];

    for (cofactor_sign, det) in minors {
        if det != 0.0 {
            let det_sign = if det > 0.0 { 1_i8 } else { -1_i8 };
            return (parity * cofactor_sign * det_sign) > 0;
        }
    }
    false
}

/// Exact in-sphere predicate with deterministic Simulation of Simplicity (SoS)
/// tie-breaking for cospherical points.
///
/// Perturbs each vertex's paraboloid lifting coordinate by `+eps^(2^index)`
/// (smallest index lifted most). Returns `true` when `points[e]` lies inside
/// the circumsphere of the positively oriented tetrahedron `(a, b, c, d)`.
pub(super) fn insphere_sos(
    points: &[[f64; 3]],
    a: usize,
    b: usize,
    c: usize,
    d: usize,
    e: usize,
) -> bool {
    let mut idx = [a, b, c, d, e];
    let parity = sort_with_parity(&mut idx);
    for w in idx.windows(2) {
        if w[0] == w[1] {
            return false;
        }
    }

    let val = insphere(points[a], points[b], points[c], points[d], points[e]);
    if val > 0.0 {
        return true;
    }
    if val < 0.0 {
        return false;
    }

    let minors = [
        (
            -1_i8,
            orient3d(
                points[idx[1]],
                points[idx[2]],
                points[idx[3]],
                points[idx[4]],
            ),
        ),
        (
            1_i8,
            orient3d(
                points[idx[0]],
                points[idx[2]],
                points[idx[3]],
                points[idx[4]],
            ),
        ),
        (
            -1_i8,
            orient3d(
                points[idx[0]],
                points[idx[1]],
                points[idx[3]],
                points[idx[4]],
            ),
        ),
        (
            1_i8,
            orient3d(
                points[idx[0]],
                points[idx[1]],
                points[idx[2]],
                points[idx[4]],
            ),
        ),
        (
            -1_i8,
            orient3d(
                points[idx[0]],
                points[idx[1]],
                points[idx[2]],
                points[idx[3]],
            ),
        ),
    ];

    for (cofactor_sign, det) in minors {
        if det != 0.0 {
            let det_sign = if det > 0.0 { 1_i8 } else { -1_i8 };
            return (parity * cofactor_sign * det_sign) > 0;
        }
    }
    false
}

/// The circumcentre and squared circumradius of `a`, `b`, `c`.
///
/// Returns [`None`] when the three points are collinear to working precision,
/// in which case the circumcircle is not defined.
pub(super) fn circumcircle(
    a: [f64; 2],
    b: [f64; 2],
    c: [f64; 2],
) -> Option<([f64; 2], f64)> {
    let bx = b[0] - a[0];
    let by = b[1] - a[1];
    let cx = c[0] - a[0];
    let cy = c[1] - a[1];

    let det = 2.0 * (bx * cy - by * cx);
    let scale = (bx * bx + by * by)
        .max(cx * cx + cy * cy)
        .max(f64::MIN_POSITIVE);
    if det.abs() <= 1e-14 * scale {
        return None;
    }

    let b2 = bx * bx + by * by;
    let c2 = cx * cx + cy * cy;
    let ux = (cy * b2 - by * c2) / det;
    let uy = (bx * c2 - cx * b2) / det;

    Some(([a[0] + ux, a[1] + uy], ux * ux + uy * uy))
}

/// The normalised key of the undirected edge between `a` and `b`.
pub(super) fn edge_key(
    a: usize,
    b: usize,
) -> (usize, usize) {
    if a <= b { (a, b) } else { (b, a) }
}

/// The area of the triangle `a`, `b`, `c`, always non-negative.
pub(super) fn triangle_area(
    a: [f64; 2],
    b: [f64; 2],
    c: [f64; 2],
) -> f64 {
    0.5 * orient2d(a, b, c).abs()
}

/// The angle in radians at `p` in the corner `q`-`p`-`r`.
pub(super) fn angle_at(
    p: [f64; 2],
    q: [f64; 2],
    r: [f64; 2],
) -> f64 {
    let ux = q[0] - p[0];
    let uy = q[1] - p[1];
    let vx = r[0] - p[0];
    let vy = r[1] - p[1];

    let du = ux.hypot(uy);
    let dv = vx.hypot(vy);
    if du <= 0.0 || dv <= 0.0 {
        return 0.0;
    }

    let cosine = ((ux * vx + uy * vy) / (du * dv)).clamp(-1.0, 1.0);
    cosine.acos()
}

/// The smallest angle in radians of the triangle `a`, `b`, `c`.
pub(super) fn triangle_min_angle(
    a: [f64; 2],
    b: [f64; 2],
    c: [f64; 2],
) -> f64 {
    angle_at(a, b, c)
        .min(angle_at(b, c, a))
        .min(angle_at(c, a, b))
}

/// The centroid of the triangle `a`, `b`, `c`.
pub(super) fn centroid(
    a: [f64; 2],
    b: [f64; 2],
    c: [f64; 2],
) -> [f64; 2] {
    [(a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0]
}

/// The midpoint of the segment `a`-`b`.
pub(super) fn midpoint(
    a: [f64; 2],
    b: [f64; 2],
) -> [f64; 2] {
    [f64::midpoint(a[0], b[0]), f64::midpoint(a[1], b[1])]
}

/// The squared distance between `a` and `b`.
pub(super) fn distance2(
    a: [f64; 2],
    b: [f64; 2],
) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    dx * dx + dy * dy
}

/// Returns `t` with its vertices reordered so that it winds counter-clockwise.
pub(super) fn as_ccw(
    points: &[[f64; 2]],
    mut t: [usize; 3],
) -> [usize; 3] {
    if orient2d(points[t[0]], points[t[1]], points[t[2]]) < 0.0 {
        t.swap(1, 2);
    }
    t
}

/// Whether `p` lies strictly inside the counter-clockwise triangle `a`, `b`, `c`.
pub(super) fn point_in_triangle(
    p: [f64; 2],
    a: [f64; 2],
    b: [f64; 2],
    c: [f64; 2],
    eps: f64,
) -> bool {
    orient2d(a, b, p) > eps && orient2d(b, c, p) > eps && orient2d(c, a, p) > eps
}

/// Whether the two open segments `a`-`b` and `c`-`d` cross properly.
///
/// Endpoint contact is not reported; only a transversal crossing of the two
/// interiors counts.
pub(super) fn segments_cross(
    a: [f64; 2],
    b: [f64; 2],
    c: [f64; 2],
    d: [f64; 2],
) -> bool {
    let o1 = orient2d(a, b, c);
    let o2 = orient2d(a, b, d);
    let o3 = orient2d(c, d, a);
    let o4 = orient2d(c, d, b);

    let opposite = |u: f64, v: f64| (u > 0.0 && v < 0.0) || (u < 0.0 && v > 0.0);
    opposite(o1, o2) && opposite(o3, o4)
}

/// The distance from `p` to the segment `a`-`b`.
pub(super) fn point_segment_distance(
    p: [f64; 2],
    a: [f64; 2],
    b: [f64; 2],
) -> f64 {
    let abx = b[0] - a[0];
    let aby = b[1] - a[1];
    let len2 = abx * abx + aby * aby;
    if len2 <= 0.0 {
        return distance2(p, a).sqrt();
    }

    let t = (((p[0] - a[0]) * abx + (p[1] - a[1]) * aby) / len2).clamp(0.0, 1.0);
    let q = [a[0] + t * abx, a[1] + t * aby];
    distance2(p, q).sqrt()
}

/// Whether `p` lies inside the closed polygon `poly` by the even-odd rule.
pub(super) fn point_in_polygon(
    p: [f64; 2],
    poly: &[[f64; 2]],
) -> bool {
    if poly.len() < 3 {
        return false;
    }

    let mut inside = false;
    let n = poly.len();
    let mut j = n - 1;
    for i in 0..n {
        let pi = poly[i];
        let pj = poly[j];
        if (pi[1] > p[1]) != (pj[1] > p[1]) {
            let t = (p[1] - pi[1]) / (pj[1] - pi[1]);
            if p[0] < pi[0] + t * (pj[0] - pi[0]) {
                inside = !inside;
            }
        }
        j = i;
    }
    inside
}

/// The signed area of the polygon `poly`, positive when it winds counter-clockwise.
pub(super) fn polygon_signed_area(poly: &[[f64; 2]]) -> f64 {
    if poly.len() < 3 {
        return 0.0;
    }

    let mut acc = 0.0;
    let n = poly.len();
    let mut j = n - 1;
    for i in 0..n {
        acc += poly[j][0] * poly[i][1] - poly[i][0] * poly[j][1];
        j = i;
    }
    acc / 2.0
}

/// Whether `p` lies on the convex-hull boundary of `points` up to `tol`.
///
/// A point is on the hull when every other point sits on one side of the line
/// through it and at least two other points are not collinear with it.
pub(super) fn is_hull_edge(
    p: [f64; 2],
    q: [f64; 2],
    points: &[[f64; 2]],
    tol: f64,
) -> bool {
    let mut has_left = false;
    let mut has_right = false;
    for &r in points {
        let o = orient2d(p, q, r);
        if o > tol {
            has_left = true;
        } else if o < -tol {
            has_right = true;
        }
        if has_left && has_right {
            return false;
        }
    }
    has_left || has_right
}

/// The area of every triangle in `triangles`.
///
/// Triangles whose indices do not address `points` are reported as `0.0`.
#[must_use]
pub fn triangle_areas(
    points: &[[f64; 2]],
    triangles: &[[usize; 3]],
) -> Vec<f64> {
    let n = points.len();
    triangles
        .iter()
        .map(|t| {
            if t[0] < n && t[1] < n && t[2] < n {
                triangle_area(points[t[0]], points[t[1]], points[t[2]])
            } else {
                0.0
            }
        })
        .collect()
}

/// The smallest angle in degrees over all `triangles`.
///
/// Returns `0.0` for an empty mesh or a degenerate triangle.
#[must_use]
pub fn min_angle_deg(
    points: &[[f64; 2]],
    triangles: &[[usize; 3]],
) -> f64 {
    let n = points.len();
    let mut best = f64::INFINITY;
    for t in triangles {
        if t[0] >= n || t[1] >= n || t[2] >= n {
            continue;
        }
        let angle = triangle_min_angle(points[t[0]], points[t[1]], points[t[2]]);
        best = best.min(angle);
    }
    if best.is_finite() {
        best.to_degrees()
    } else {
        0.0
    }
}

/// The number of unique undirected edges spanned by `triangles`.
#[cfg(test)]
pub(super) fn count_edges(triangles: &[[usize; 3]]) -> usize {
    let mut seen: std::collections::HashSet<(usize, usize)> = std::collections::HashSet::new();
    for t in triangles {
        for k in 0..3 {
            seen.insert(edge_key(t[k], t[(k + 1) % 3]));
        }
    }
    seen.len()
}

/// Triangulates the simple counter-clockwise polygon `poly` by ear clipping.
///
/// The result is a list of counter-clockwise triangles that cover the polygon.
/// A polygon that cannot be reduced to ears (because it is not simple) yields
/// whatever triangles could be clipped before the guard ran out.
pub(super) fn ear_clip(
    points: &[[f64; 2]],
    poly: &[usize],
) -> Vec<[usize; 3]> {
    let mut ring: Vec<usize> = poly.to_vec();
    if ring.len() >= 2 && ring[0] == ring[ring.len() - 1] {
        ring.pop();
    }
    if ring.len() < 3 {
        return Vec::new();
    }

    let mut out: Vec<[usize; 3]> = Vec::new();
    let mut guard = 0;
    let limit = 4 * ring.len() * ring.len() + 16;

    while ring.len() > 3 {
        guard += 1;
        if guard > limit {
            return out;
        }

        let n = ring.len();
        let mut clipped = false;
        for i in 0..n {
            let prev = ring[(i + n - 1) % n];
            let cur = ring[i];
            let next = ring[(i + 1) % n];

            if orient2d(points[prev], points[cur], points[next]) <= 1e-14 {
                continue;
            }

            let mut blocked = false;
            for &q in &ring {
                if q == prev || q == cur || q == next {
                    continue;
                }
                if point_in_triangle(points[q], points[prev], points[cur], points[next], 1e-14) {
                    blocked = true;
                    break;
                }
            }
            if blocked {
                continue;
            }

            out.push([prev, cur, next]);
            ring.remove(i);
            clipped = true;
            break;
        }

        if !clipped {
            // Drop a collinear vertex so the loop can still make progress.
            let n = ring.len();
            let mut removed = false;
            for i in 0..n {
                let prev = ring[(i + n - 1) % n];
                let cur = ring[i];
                let next = ring[(i + 1) % n];
                if orient2d(points[prev], points[cur], points[next]).abs() <= 1e-14 {
                    ring.remove(i);
                    removed = true;
                    break;
                }
            }
            if !removed {
                return out;
            }
        }
    }

    if ring.len() == 3 && orient2d(points[ring[0]], points[ring[1]], points[ring[2]]) > 0.0 {
        out.push([ring[0], ring[1], ring[2]]);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::TriMesh2;
    use super::count_edges;
    use super::incircle;
    use super::incircle_sos;
    use super::insphere;
    use super::insphere_sos;
    use super::min_angle_deg;
    use super::orient2d;
    use super::orient3d;
    use super::point_in_polygon;
    use super::triangle_areas;
    use super::triangle_min_angle;

    const UNIT_SQUARE: [[f64; 2]; 4] = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];

    #[test]
    fn orient2d_signs_follow_the_winding() {
        assert!(orient2d(UNIT_SQUARE[0], UNIT_SQUARE[1], UNIT_SQUARE[2]) > 0.0);
        assert!(orient2d(UNIT_SQUARE[0], UNIT_SQUARE[2], UNIT_SQUARE[1]) < 0.0);
        assert_eq!(orient2d([0.0, 0.0], [1.0, 1.0], [2.0, 2.0]), 0.0);
    }

    #[test]
    fn orient2d_exact_on_nearly_collinear_large_coordinates() {
        // Exact determinant is (1e8 + 1)*(1e8 - 1) - 1e8*1e8 = -1.0,
        // whereas plain f64 rounds 1e16 - 1 to 1e16 and yields 0.0.
        let a = [0.0, 0.0];
        let b = [1e8 + 1.0, 1e8];
        let c = [1e8, 1e8 - 1.0];
        assert!(orient2d(a, b, c) < 0.0);
        assert!(orient2d(a, c, b) > 0.0);
        assert_eq!(orient2d(a, [1e8, 1e8], [2e8, 2e8]), 0.0);
    }

    #[test]
    fn incircle_detects_the_interior_of_the_circumcircle() {
        let a = [0.0, 0.0];
        let b = [1.0, 0.0];
        let c = [0.0, 1.0];
        // The circumcentre of the unit right triangle is (0.5, 0.5) with
        // squared radius 0.5, so (0.25, 0.25) is inside and (2, 2) outside.
        assert!(incircle(a, b, c, [0.25, 0.25]) > 0.0);
        assert_eq!(incircle(a, b, c, [1.0, 1.0]), 0.0);
        assert!(incircle(a, b, c, [2.0, 2.0]) < 0.0);
    }

    #[test]
    fn incircle_exact_on_nearly_cocircular_large_coordinates() {
        let r = 1e5;
        let a = [r, 0.0];
        let b = [0.0, r];
        let c = [-r, 0.0];
        assert_eq!(incircle(a, b, c, [0.0, -r]), 0.0);
        assert!(incircle(a, b, c, [0.0, -r + 1e-10]) > 0.0);
        assert!(incircle(a, b, c, [0.0, -r - 1e-10]) < 0.0);
    }

    #[test]
    fn orient3d_exact_on_nearly_coplanar_large_coordinates() {
        let a = [0.0, 0.0, 0.0];
        let b = [1e8 + 1.0, 1e8, 0.0];
        let c = [1e8, 1e8 - 1.0, 0.0];
        let d = [0.0, 0.0, 1.0];
        assert!(orient3d(a, b, c, d) < 0.0);
        assert!(orient3d(a, c, b, d) > 0.0);
        assert_eq!(
            orient3d(a, [1e8, 1e8, 0.0], [2e8, 2e8, 0.0], [0.0, 0.0, 1.0]),
            0.0
        );
    }

    #[test]
    fn insphere_exact_on_nearly_cospherical_points() {
        let r = 1e4;
        let a = [r, 0.0, 0.0];
        let b = [0.0, r, 0.0];
        let c = [-r, 0.0, 0.0];
        let d = [0.0, 0.0, r];
        assert!(orient3d(a, b, c, d) > 0.0);
        assert_eq!(insphere(a, b, c, d, [0.0, 0.0, -r]), 0.0);
        assert!(insphere(a, b, c, d, [0.0, 0.0, -r + 1e-10]) > 0.0);
        assert!(insphere(a, b, c, d, [0.0, 0.0, -r - 1e-10]) < 0.0);
    }

    #[test]
    fn sos_tie_breakers_are_consistent_on_cocircular_and_cospherical_grids() {
        // Cocircular unit square: cyclic shift (0,1,2,3) -> (1,2,3,0) is an odd
        // permutation, so exactly one of the two diagonals is strictly preferred.
        let s0 = incircle_sos(&UNIT_SQUARE, 0, 1, 2, 3);
        let s1 = incircle_sos(&UNIT_SQUARE, 1, 2, 3, 0);
        assert_ne!(s0, s1);
        // Even permutation of the triangle vertices preserves the result.
        assert_eq!(incircle_sos(&UNIT_SQUARE, 1, 2, 0, 3), s0);
        assert_eq!(incircle_sos(&UNIT_SQUARE, 2, 0, 1, 3), s0);

        // Cospherical unit cube vertices.
        let cube = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
        ];
        assert!(orient3d(cube[0], cube[1], cube[2], cube[3]) > 0.0);
        assert_eq!(insphere(cube[0], cube[1], cube[2], cube[3], cube[4]), 0.0);
        let in0 = insphere_sos(&cube, 0, 1, 2, 3, 4);
        // Even permutation of (0,1,2,3) preserves the SoS result.
        assert_eq!(insphere_sos(&cube, 1, 2, 0, 3, 4), in0);
        // Odd permutation of (0,1,2,3) flips the SoS result.
        assert_ne!(insphere_sos(&cube, 1, 0, 2, 3, 4), in0);
    }

    #[test]
    fn triangle_areas_and_min_angle_match_the_square() {
        let points = UNIT_SQUARE.to_vec();
        let triangles = [[0, 1, 2], [0, 2, 3]];
        let areas = triangle_areas(&points, &triangles);
        assert!((areas[0] - 0.5).abs() < 1e-12);
        assert!((areas[1] - 0.5).abs() < 1e-12);
        assert!((min_angle_deg(&points, &triangles) - 45.0).abs() < 1e-9);
        assert_eq!(count_edges(&triangles), 5);
    }

    #[test]
    fn degenerate_triangle_has_zero_minimum_angle() {
        let angle = triangle_min_angle([0.0, 0.0], [1.0, 0.0], [2.0, 0.0]);
        assert!(angle < 1e-12);
        assert_eq!(min_angle_deg(&UNIT_SQUARE, &[]), 0.0);
    }

    #[test]
    fn point_in_polygon_accepts_the_interior_only() {
        let poly = UNIT_SQUARE.to_vec();
        assert!(point_in_polygon([0.5, 0.5], &poly));
        assert!(!point_in_polygon([1.5, 0.5], &poly));
        assert!(!point_in_polygon([-0.1, -0.1], &poly));
    }

    #[test]
    fn empty_mesh_is_the_default() {
        let mesh = TriMesh2::default();
        assert_eq!(mesh.points.len(), 0);
        assert_eq!(mesh.triangles.len(), 0);
    }
}
