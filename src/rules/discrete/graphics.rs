//! Computer graphics: homogeneous transformation matrices, curves,
//! quaternions, meshes and ray queries.
//!
//! Points and vectors are `list(x, y, z)`, matrices `list(list(..), ..)`
//! (row major, acting on column vectors, so `M * p` is `matmul(M, p)`), and
//! quaternions `list(w, x, y, z)`. Everything that can stay symbolic does:
//! `rotation_2d(theta)` has `cos(theta)` entries, `bezier(.., t)` accepts a
//! symbolic parameter, and numbers fold to numbers
//! (`rotation_2d(pi/2)` is an exact permutation matrix). The ray, slerp,
//! refraction and barycentric queries are numeric and wrap
//! [`crate::kernels::computer_graphics`]; they stay unreduced for symbolic
//! arguments.
//!
//! | operator | value |
//! |---|---|
//! | `translation_2d(tx, ty)`, `translation_3d(tx, ty, tz)` | 3x3 / 4x4 translation |
//! | `scaling_2d(sx, sy)`, `scaling_3d(sx, sy, sz)`, `shear_2d(shx, shy)` | scaling and shear |
//! | `rotation_2d(theta)`, `rotation_3d_x/y/z(theta)` | rotation about the origin / a coordinate axis |
//! | `rotation_axis_angle(list(ax, ay, az), theta)` | Rodrigues rotation about the (normalised) axis, 4x4 |
//! | `reflection_2d(theta)` | reflection in the line through the origin at angle `theta`, 3x3 |
//! | `reflection_3d(list(a, b, c))` | reflection in the plane through the origin with normal `(a, b, c)`, 4x4 |
//! | `perspective(fovy, aspect, near, far)`, `orthographic(l, r, b, t, n, f)` | OpenGL style projections |
//! | `look_at(eye, center, up)` | view matrix (entries may be symbolic) |
//! | `apply_transform(M, p)`, `apply_transform_vector(M, v)` | `p` of the matrix size is multiplied as is (homogeneous); one shorter, a 1 (point, followed by the perspective divide) or 0 (vector) is appended |
//! | `bezier(P, t)`, `bezier_derivative(P, t)`, `bezier_split(P, t0)` | de Casteljau evaluation, derivative, split into `list(left, right)` control polygons |
//! | `bspline(P, degree, knots, t)` | de Boor evaluation for a numeric `t` in `[knots[degree], knots[n]]` |
//! | `catmull_rom(p0, p1, p2, p3, t)` | Catmull-Rom segment between `p1` and `p2` |
//! | `quat_mul(q, r)`, `quat_conj(q)`, `quat_inverse(q)`, `quat_norm(q)`, `quat_normalize(q)` | quaternion algebra |
//! | `quat_from_axis_angle(axis, theta)`, `quat_rotate(q, v)`, `quat_to_matrix(q)` | rotations (`quat_rotate` computes `q v q*`; `q` should be a unit quaternion) |
//! | `quat_slerp(q1, q2, t)` | numeric spherical interpolation |
//! | `mesh_transform(vertices, M)`, `mesh_normals(vertices, faces)`, `mesh_triangulate(faces)` | meshes are a vertex list and a face list of vertex indices; transformed vertices get the perspective divide, normals are the unit normals of the first three vertices of each face with at least three, triangulation is a fan |
//! | `ray_sphere(o, d, center, r)`, `ray_plane(o, d, point, n)`, `ray_triangle(o, d, a, b, c)` | numeric: `list(t, point, normal)` of the nearest hit, `false` for a miss |
//! | `reflect(d, n)` | `d - 2 (d.n) n` (symbolic) |
//! | `refract(d, n, eta)` | numeric Snell refraction, `false` on total internal reflection |
//! | `barycentric(p, a, b, c)` | numeric barycentric coordinates `list(u, v, w)` |

use super::apply;
use super::def;
use super::def_request;
use super::div;
use super::float;
use super::idx;
use super::items;
use super::matmul;
use super::matrix;
use super::neg;
use super::prod;
use super::rows;
use super::sum;
use super::V;
use crate::graph::rule::Installer;
use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::RuleError;
use crate::kernels::computer_graphics as cg;

type Row = Vec<NodeId>;

// ----------------------------------------------------------------------
// Symbolic scalar helpers (every result is simplified)
// ----------------------------------------------------------------------

fn int(
    cx: &mut Cx<'_>,
    v: i64,
) -> NodeId {
    cx.graph.int(v)
}

fn add(
    cx: &mut Cx<'_>,
    a: NodeId,
    b: NodeId,
) -> NodeId {
    let s = sum(cx.graph, &[a, b]);
    cx.simplify(s)
}

fn sub(
    cx: &mut Cx<'_>,
    a: NodeId,
    b: NodeId,
) -> NodeId {
    let nb = neg(cx.graph, b);
    add(cx, a, nb)
}

fn mul(
    cx: &mut Cx<'_>,
    a: NodeId,
    b: NodeId,
) -> NodeId {
    let p = prod(cx.graph, &[a, b]);
    cx.simplify(p)
}

fn quot(
    cx: &mut Cx<'_>,
    a: NodeId,
    b: NodeId,
) -> NodeId {
    let d = div(cx.graph, a, b);
    cx.simplify(d)
}

fn negate(
    cx: &mut Cx<'_>,
    a: NodeId,
) -> NodeId {
    let n = neg(cx.graph, a);
    cx.simplify(n)
}

/// `name(x)`, simplified.
fn call1(
    cx: &mut Cx<'_>,
    name: &str,
    x: NodeId,
) -> Option<NodeId> {
    let n = apply(cx.graph, name, &[x])?;
    Some(cx.simplify(n))
}

fn sqrt(
    cx: &mut Cx<'_>,
    x: NodeId,
) -> Option<NodeId> {
    call1(cx, "sqrt", x)
}

/// `sum(coef_i * a_i * b_i)`.
fn combine(
    cx: &mut Cx<'_>,
    terms: &[(i64, NodeId, NodeId)],
) -> NodeId {
    let nodes: Vec<NodeId> = terms
        .iter()
        .map(|&(k, a, b)| {
            let k = cx.graph.int(k);
            prod(cx.graph, &[k, a, b])
        })
        .collect();
    let s = sum(cx.graph, &nodes);
    cx.simplify(s)
}

fn dot(
    cx: &mut Cx<'_>,
    a: &[NodeId],
    b: &[NodeId],
) -> NodeId {
    let terms: Vec<(i64, NodeId, NodeId)> = a.iter().zip(b).map(|(&x, &y)| (1, x, y)).collect();
    combine(cx, &terms)
}

fn cross(
    cx: &mut Cx<'_>,
    a: &[NodeId],
    b: &[NodeId],
) -> Option<Row> {
    let ([a0, a1, a2], [b0, b1, b2]) = (<[NodeId; 3]>::try_from(a).ok()?, <[NodeId; 3]>::try_from(b).ok()?);
    Some(vec![
        combine(cx, &[(1, a1, b2), (-1, a2, b1)]),
        combine(cx, &[(1, a2, b0), (-1, a0, b2)]),
        combine(cx, &[(1, a0, b1), (-1, a1, b0)]),
    ])
}

fn normalize(
    cx: &mut Cx<'_>,
    v: &[NodeId],
) -> Option<Row> {
    let sq = dot(cx, v, v);
    let len = sqrt(cx, sq)?;
    Some(v.iter().map(|&x| quot(cx, x, len)).collect())
}

fn minus(
    cx: &mut Cx<'_>,
    a: &[NodeId],
    b: &[NodeId],
) -> Row {
    a.iter().zip(b).map(|(&x, &y)| sub(cx, x, y)).collect()
}

fn identity(
    cx: &mut Cx<'_>,
    n: usize,
) -> Vec<Row> {
    let zero = int(cx, 0);
    let one = int(cx, 1);
    (0..n).map(|i| (0..n).map(|j| if i == j { one } else { zero }).collect()).collect()
}

fn done(
    cx: &mut Cx<'_>,
    rows: &[Row],
) -> V {
    V::Node(matrix(cx.graph, rows))
}

fn point(
    g: &Graph,
    n: NodeId,
) -> Option<Row> {
    items(g, n)
}

fn arg<const N: usize>(a: &[NodeId]) -> Option<[NodeId; N]> {
    <[NodeId; N]>::try_from(a).ok()
}

fn fixed_vec<const N: usize>(
    g: &Graph,
    n: NodeId,
) -> Option<[NodeId; N]> {
    arg(&items(g, n)?)
}

// ----------------------------------------------------------------------
// Affine transformations
// ----------------------------------------------------------------------

/// `I` with the last column's first `t.len()` entries replaced by `t`.
fn translation(
    cx: &mut Cx<'_>,
    t: &[NodeId],
) -> Option<V> {
    let mut m = identity(cx, t.len() + 1);
    for (i, &x) in t.iter().enumerate() {
        m[i][t.len()] = x;
    }
    Some(done(cx, &m))
}

/// A diagonal matrix with `d` followed by 1.
fn diagonal(
    cx: &mut Cx<'_>,
    d: &[NodeId],
) -> Option<V> {
    let mut m = identity(cx, d.len() + 1);
    for (i, &x) in d.iter().enumerate() {
        m[i][i] = x;
    }
    Some(done(cx, &m))
}

fn shear_2d(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [shx, shy] = arg(a)?;
    let mut m = identity(cx, 3);
    m[0][1] = shx;
    m[1][0] = shy;
    Some(done(cx, &m))
}

/// `(cos a, sin a, -sin a)`.
fn trig(
    cx: &mut Cx<'_>,
    angle: NodeId,
) -> Option<(NodeId, NodeId, NodeId)> {
    let c = call1(cx, "cos", angle)?;
    let s = call1(cx, "sin", angle)?;
    let ns = negate(cx, s);
    Some((c, s, ns))
}

/// Rotation in the plane of coordinates `(i, j)` of an `n x n` identity.
fn plane_rotation(
    cx: &mut Cx<'_>,
    n: usize,
    angle: NodeId,
    (i, j): (usize, usize),
    flip: bool,
) -> Option<V> {
    let (c, s, ns) = trig(cx, angle)?;
    let mut m = identity(cx, n);
    m[i][i] = c;
    m[j][j] = c;
    m[i][j] = if flip { s } else { ns };
    m[j][i] = if flip { ns } else { s };
    Some(done(cx, &m))
}

fn rotation_axis_angle(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [axis, angle] = arg(a)?;
    let u = normalize(cx, &fixed_vec::<3>(cx.graph, axis)?)?;
    let (c, s, _) = trig(cx, angle)?;
    let one = int(cx, 1);
    let omc = sub(cx, one, c);
    let mut m = identity(cx, 4);
    for i in 0..3 {
        for j in 0..3 {
            let uu = mul(cx, u[i], u[j]);
            let first = mul(cx, uu, omc);
            m[i][j] = if i == j {
                add(cx, c, first)
            } else {
                // The skew part: +-u_k sin(angle) with (i, j, k) cyclic.
                let k = 3 - i - j;
                let positive = (j + 1) % 3 == i;
                let us = mul(cx, u[k], s);
                if positive { add(cx, first, us) } else { sub(cx, first, us) }
            };
        }
    }
    Some(done(cx, &m))
}

fn reflection_2d(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [theta] = arg(a)?;
    let two = int(cx, 2);
    let double = mul(cx, two, theta);
    let (c, s, _) = trig(cx, double)?;
    let nc = negate(cx, c);
    let mut m = identity(cx, 3);
    m[0][0] = c;
    m[0][1] = s;
    m[1][0] = s;
    m[1][1] = nc;
    Some(done(cx, &m))
}

fn reflection_3d(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [normal] = arg(a)?;
    let n = fixed_vec::<3>(cx.graph, normal)?;
    let len2 = dot(cx, &n, &n);
    let one = int(cx, 1);
    let mut m = identity(cx, 4);
    for i in 0..3 {
        for j in 0..3 {
            let nn = mul(cx, n[i], n[j]);
            let ratio = quot(cx, nn, len2);
            let twice = combine(cx, &[(2, ratio, one)]);
            m[i][j] = if i == j {
                sub(cx, one, twice)
            } else {
                negate(cx, twice)
            };
        }
    }
    Some(done(cx, &m))
}

// ----------------------------------------------------------------------
// Projections and the camera
// ----------------------------------------------------------------------

fn perspective(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [fovy, aspect, near, far] = arg(a)?;
    let two = int(cx, 2);
    let one = int(cx, 1);
    let half = quot(cx, fovy, two);
    let f = call1(cx, "tan", half)?;
    let fa = mul(cx, f, aspect);
    let range = sub(cx, near, far);
    let sum_nf = add(cx, near, far);
    let nf = mul(cx, near, far);
    let two_nf = mul(cx, two, nf);
    let mut m = identity(cx, 4);
    m[0][0] = quot(cx, one, fa);
    m[1][1] = quot(cx, one, f);
    m[2][2] = quot(cx, sum_nf, range);
    m[2][3] = quot(cx, two_nf, range);
    m[3][2] = int(cx, -1);
    m[3][3] = int(cx, 0);
    Some(done(cx, &m))
}

fn orthographic(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [l, r, b, t, n, f] = arg(a)?;
    let two = int(cx, 2);
    let mut m = identity(cx, 4);
    for (k, (lo, hi)) in [(l, r), (b, t), (n, f)].into_iter().enumerate() {
        let span = sub(cx, hi, lo);
        let total = add(cx, hi, lo);
        let scale = quot(cx, two, span);
        let shift = quot(cx, total, span);
        // The z axis looks down -z: scale -2/(f-n), shift -(f+n)/(f-n).
        m[k][k] = if k == 2 { negate(cx, scale) } else { scale };
        m[k][3] = negate(cx, shift);
    }
    Some(done(cx, &m))
}

fn look_at(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [eye, center, up] = arg(a)?;
    let (eye, center, up) = (
        fixed_vec::<3>(cx.graph, eye)?,
        fixed_vec::<3>(cx.graph, center)?,
        fixed_vec::<3>(cx.graph, up)?,
    );
    let dir = minus(cx, &center, &eye);
    let f = normalize(cx, &dir)?;
    let side = cross(cx, &f, &up)?;
    let s = normalize(cx, &side)?;
    let u = cross(cx, &s, &f)?;
    let mut m = identity(cx, 4);
    for k in 0..3 {
        m[0][k] = s[k];
        m[1][k] = u[k];
        m[2][k] = negate(cx, f[k]);
    }
    let se = dot(cx, &s, &eye);
    let ue = dot(cx, &u, &eye);
    m[0][3] = negate(cx, se);
    m[1][3] = negate(cx, ue);
    m[2][3] = dot(cx, &f, &eye);
    Some(done(cx, &m))
}

// ----------------------------------------------------------------------
// Transforming points and meshes
// ----------------------------------------------------------------------

/// `M * p`: homogeneous when `p` has the matrix size, otherwise `p` gets a
/// trailing `w` (1 for points, which are divided by the result's `w`, and 0
/// for vectors).
fn transformed(
    cx: &mut Cx<'_>,
    m: &[Row],
    p: &[NodeId],
    is_point: bool,
) -> Option<Row> {
    let n = m.len();
    if m.iter().any(|r| r.len() != n) {
        return None;
    }
    let mut col = p.to_vec();
    let homogeneous = p.len() == n;
    if !homogeneous {
        if p.len() + 1 != n {
            return None;
        }
        col.push(int(cx, i64::from(is_point)));
    }
    let column: Vec<Row> = col.iter().map(|&x| vec![x]).collect();
    let out: Row = matmul(cx, m, &column)?.into_iter().map(|r| r[0]).collect();
    if homogeneous {
        return Some(out);
    }
    let (coords, w) = out.split_at(n - 1);
    if !is_point || cx.graph.number_of(w[0]).is_some_and(crate::graph::Number::is_one) {
        return Some(coords.to_vec());
    }
    if cx.graph.number_of(w[0]).is_some_and(crate::graph::Number::is_zero) {
        return None;
    }
    Some(coords.iter().map(|&x| quot(cx, x, w[0])).collect())
}

fn transform(
    cx: &mut Cx<'_>,
    a: &[NodeId],
    is_point: bool,
) -> Option<V> {
    let [m, p] = arg(a)?;
    let m = rows(cx.graph, m)?;
    let p = point(cx.graph, p)?;
    Some(V::nodes(&transformed(cx, &m, &p, is_point)?))
}

fn mesh_transform(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [vertices, m] = arg(a)?;
    let m = rows(cx.graph, m)?;
    let mut out = Vec::new();
    for v in rows(cx.graph, vertices)? {
        out.push(V::nodes(&transformed(cx, &m, &v, true)?));
    }
    Some(V::List(out))
}

fn faces(
    g: &Graph,
    n: NodeId,
) -> Option<Vec<Vec<usize>>> {
    items(g, n)?
        .into_iter()
        .map(|f| items(g, f)?.into_iter().map(|i| idx(g, i)).collect())
        .collect()
}

fn mesh_normals(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [vertices, fs] = arg(a)?;
    let vs = rows(cx.graph, vertices)?;
    let mut out = Vec::new();
    for f in faces(cx.graph, fs)? {
        if f.len() < 3 {
            continue;
        }
        let (v0, v1, v2) = (vs.get(f[0])?, vs.get(f[1])?, vs.get(f[2])?);
        let e1 = minus(cx, v1, v0);
        let e2 = minus(cx, v2, v0);
        let n = cross(cx, &e1, &e2)?;
        out.push(V::nodes(&normalize(cx, &n)?));
    }
    Some(V::List(out))
}

fn mesh_triangulate(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [fs] = arg(a)?;
    let mut out = Vec::new();
    for f in faces(cx.graph, fs)? {
        if f.len() <= 3 {
            out.push(V::List(f.into_iter().map(V::uint).collect()));
        } else {
            for i in 1..f.len() - 1 {
                out.push(V::List(vec![V::uint(f[0]), V::uint(f[i]), V::uint(f[i + 1])]));
            }
        }
    }
    Some(V::List(out))
}

// ----------------------------------------------------------------------
// Curves
// ----------------------------------------------------------------------

/// `(1 - t) a + t b`, coordinatewise.
fn lerp(
    cx: &mut Cx<'_>,
    a: &[NodeId],
    b: &[NodeId],
    t: NodeId,
) -> Row {
    let one = int(cx, 1);
    let omt = sub(cx, one, t);
    a.iter()
        .zip(b)
        .map(|(&x, &y)| {
            let left = mul(cx, omt, x);
            let right = mul(cx, t, y);
            add(cx, left, right)
        })
        .collect()
}

/// A control polygon: at least one point, all of one dimension.
fn polygon(
    g: &Graph,
    n: NodeId,
) -> Option<Vec<Row>> {
    let pts = rows(g, n)?;
    let dim = pts.first()?.len();
    pts.iter().all(|p| p.len() == dim).then_some(pts)
}

/// The de Casteljau pyramid; level `k` has `n - k` points.
fn pyramid(
    cx: &mut Cx<'_>,
    pts: &[Row],
    t: NodeId,
) -> Vec<Vec<Row>> {
    let mut levels = vec![pts.to_vec()];
    for k in 1..pts.len() {
        let prev = &levels[k - 1];
        let next: Vec<Row> = (0..pts.len() - k).map(|i| (i, prev[i].clone(), prev[i + 1].clone())).collect::<Vec<_>>().into_iter().map(|(_, a, b)| lerp(cx, &a, &b, t)).collect();
        levels.push(next);
    }
    levels
}

fn bezier(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [cps, t] = arg(a)?;
    let pts = polygon(cx.graph, cps)?;
    let levels = pyramid(cx, &pts, t);
    Some(V::nodes(&levels.last()?[0]))
}

fn bezier_derivative(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [cps, t] = arg(a)?;
    let pts = polygon(cx.graph, cps)?;
    let degree = i64::try_from(pts.len() - 1).ok()?;
    if degree == 0 {
        let zero = int(cx, 0);
        return Some(V::nodes(&vec![zero; pts[0].len()]));
    }
    let n = int(cx, degree);
    let diffs: Vec<Row> = pts
        .windows(2)
        .map(|w| {
            let d = minus(cx, &w[1], &w[0]);
            d.into_iter().map(|x| mul(cx, n, x)).collect()
        })
        .collect();
    let levels = pyramid(cx, &diffs, t);
    Some(V::nodes(&levels.last()?[0]))
}

fn bezier_split(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [cps, t] = arg(a)?;
    let pts = polygon(cx.graph, cps)?;
    let n = pts.len();
    let levels = pyramid(cx, &pts, t);
    let left: Vec<V> = (0..n).map(|k| V::nodes(&levels[k][0])).collect();
    let right: Vec<V> = (0..n).map(|k| V::nodes(&levels[n - 1 - k][k])).collect();
    Some(V::List(vec![V::List(left), V::List(right)]))
}

fn bspline(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [cps, degree, knots, t] = arg(a)?;
    let pts = polygon(cx.graph, cps)?;
    let p = idx(cx.graph, degree)?;
    let knot_nodes = items(cx.graph, knots)?;
    let ks: Vec<f64> = knot_nodes.iter().map(|&k| float(cx.graph, k)).collect::<Option<_>>()?;
    let tv = float(cx.graph, t)?;
    let n = pts.len();
    if ks.len() != n + p + 1 || n <= p || ks.windows(2).any(|w| w[0] > w[1]) {
        return None;
    }
    if tv < ks[p] || tv > ks[n] {
        return None;
    }
    // The span k with knots[k] <= t < knots[k + 1]; the end of the domain
    // belongs to the last non-empty span.
    let k = (p..n).rev().find(|&k| ks[k] <= tv && ks[k] < ks[k + 1])?;
    let mut d: Vec<Row> = (0..=p).map(|j| pts[j + k - p].clone()).collect();
    for r in 1..=p {
        for j in (r..=p).rev() {
            let i = j + k - p;
            let lo = knot_nodes[i];
            let hi = knot_nodes[i + p + 1 - r];
            let num = sub(cx, t, lo);
            let den = sub(cx, hi, lo);
            let alpha = quot(cx, num, den);
            let (prev, cur) = (d[j - 1].clone(), d[j].clone());
            d[j] = lerp(cx, &prev, &cur, alpha);
        }
    }
    Some(V::nodes(&d[p]))
}

fn catmull_rom(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [p0, p1, p2, p3, t] = arg(a)?;
    let (p0, p1, p2, p3) = (point(cx.graph, p0)?, point(cx.graph, p1)?, point(cx.graph, p2)?, point(cx.graph, p3)?);
    if [p1.len(), p2.len(), p3.len()].iter().any(|&l| l != p0.len()) {
        return None;
    }
    let t2 = mul(cx, t, t);
    let t3 = mul(cx, t2, t);
    let half = {
        let (one, two) = (int(cx, 1), int(cx, 2));
        quot(cx, one, two)
    };
    let mut out = Vec::new();
    for k in 0..p0.len() {
        let (a0, a1, a2, a3) = (p0[k], p1[k], p2[k], p3[k]);
        let one = int(cx, 1);
        let c1 = combine(cx, &[(-1, a0, one), (1, a2, one)]);
        let c2 = combine(cx, &[(2, a0, one), (-5, a1, one), (4, a2, one), (-1, a3, one)]);
        let c3 = combine(cx, &[(-1, a0, one), (3, a1, one), (-3, a2, one), (1, a3, one)]);
        let two_a1 = combine(cx, &[(2, a1, one)]);
        let t1 = mul(cx, c1, t);
        let tt2 = mul(cx, c2, t2);
        let tt3 = mul(cx, c3, t3);
        let total = sum(cx.graph, &[two_a1, t1, tt2, tt3]);
        let total = cx.simplify(total);
        out.push(mul(cx, half, total));
    }
    Some(V::nodes(&out))
}

// ----------------------------------------------------------------------
// Quaternions
// ----------------------------------------------------------------------

fn quat(
    g: &Graph,
    n: NodeId,
) -> Option<[NodeId; 4]> {
    fixed_vec(g, n)
}

fn qmul(
    cx: &mut Cx<'_>,
    [w1, x1, y1, z1]: [NodeId; 4],
    [w2, x2, y2, z2]: [NodeId; 4],
) -> [NodeId; 4] {
    [
        combine(cx, &[(1, w1, w2), (-1, x1, x2), (-1, y1, y2), (-1, z1, z2)]),
        combine(cx, &[(1, w1, x2), (1, x1, w2), (1, y1, z2), (-1, z1, y2)]),
        combine(cx, &[(1, w1, y2), (-1, x1, z2), (1, y1, w2), (1, z1, x2)]),
        combine(cx, &[(1, w1, z2), (1, x1, y2), (-1, y1, x2), (1, z1, w2)]),
    ]
}

fn qconj(
    cx: &mut Cx<'_>,
    [w, x, y, z]: [NodeId; 4],
) -> [NodeId; 4] {
    [w, negate(cx, x), negate(cx, y), negate(cx, z)]
}

fn quat_binary(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [p, q] = arg(a)?;
    let (p, q) = (quat(cx.graph, p)?, quat(cx.graph, q)?);
    Some(V::nodes(&qmul(cx, p, q)))
}

fn quat_conj(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let q = quat(cx.graph, *a.first()?)?;
    Some(V::nodes(&qconj(cx, q)))
}

fn quat_inverse(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let q = quat(cx.graph, *a.first()?)?;
    let n2 = dot(cx, &q, &q);
    let c = qconj(cx, q);
    Some(V::nodes(&c.iter().map(|&x| quot(cx, x, n2)).collect::<Vec<_>>()))
}

fn quat_norm(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let q = quat(cx.graph, *a.first()?)?;
    let n2 = dot(cx, &q, &q);
    Some(V::Node(sqrt(cx, n2)?))
}

fn quat_normalize(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let q = quat(cx.graph, *a.first()?)?;
    Some(V::nodes(&normalize(cx, &q)?))
}

fn quat_from_axis_angle(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [axis, angle] = arg(a)?;
    let n = normalize(cx, &fixed_vec::<3>(cx.graph, axis)?)?;
    let two = int(cx, 2);
    let half = quot(cx, angle, two);
    let c = call1(cx, "cos", half)?;
    let s = call1(cx, "sin", half)?;
    let mut out = vec![c];
    out.extend(n.iter().map(|&x| mul(cx, x, s)));
    Some(V::nodes(&out))
}

fn quat_rotate(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [q, v] = arg(a)?;
    let (q, v) = (quat(cx.graph, q)?, fixed_vec::<3>(cx.graph, v)?);
    let zero = int(cx, 0);
    let qv = qmul(cx, q, [zero, v[0], v[1], v[2]]);
    let c = qconj(cx, q);
    let r = qmul(cx, qv, c);
    Some(V::nodes(&r[1..]))
}

fn quat_to_matrix(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [w, x, y, z] = quat(cx.graph, *a.first()?)?;
    let n2 = dot(cx, &[w, x, y, z], &[w, x, y, z]);
    // 2 a b / |q|^2.
    let pair = |cx: &mut Cx<'_>, a: NodeId, b: NodeId| {
        let ab = combine(cx, &[(2, a, b)]);
        quot(cx, ab, n2)
    };
    let (xx, yy, zz) = (pair(cx, x, x), pair(cx, y, y), pair(cx, z, z));
    let (xy, xz, yz) = (pair(cx, x, y), pair(cx, x, z), pair(cx, y, z));
    let (wx, wy, wz) = (pair(cx, w, x), pair(cx, w, y), pair(cx, w, z));
    let one = int(cx, 1);
    let diag = |cx: &mut Cx<'_>, p: NodeId, q: NodeId| {
        let s = add(cx, p, q);
        sub(cx, one, s)
    };
    let mut m = identity(cx, 4);
    m[0][0] = diag(cx, yy, zz);
    m[1][1] = diag(cx, xx, zz);
    m[2][2] = diag(cx, xx, yy);
    m[0][1] = sub(cx, xy, wz);
    m[0][2] = add(cx, xz, wy);
    m[1][0] = add(cx, xy, wz);
    m[1][2] = sub(cx, yz, wx);
    m[2][0] = sub(cx, xz, wy);
    m[2][1] = add(cx, yz, wx);
    Some(done(cx, &m))
}

// ----------------------------------------------------------------------
// Numeric queries (kernels)
// ----------------------------------------------------------------------

/// A 2D or 3D literal vector (a missing `z` is 0).
fn vec3(
    g: &Graph,
    n: NodeId,
) -> Option<(f64, f64, f64)> {
    let v: Vec<f64> = items(g, n)?.into_iter().map(|c| float(g, c)).collect::<Option<_>>()?;
    match v.as_slice() {
        | [x, y] => Some((*x, *y, 0.0)),
        | [x, y, z] => Some((*x, *y, *z)),
        | _ => None,
    }
}

fn floats3(v: (f64, f64, f64)) -> V {
    V::List(vec![V::Float(v.0), V::Float(v.1), V::Float(v.2)])
}

const fn vector3(v: (f64, f64, f64)) -> cg::Vector3D {
    cg::Vector3D::new(v.0, v.1, v.2)
}

const fn point3(v: (f64, f64, f64)) -> cg::Point3D {
    cg::Point3D::new(v.0, v.1, v.2)
}

fn hit(h: Option<cg::Intersection>) -> V {
    h.map_or(V::Bool(false), |h| {
        V::List(vec![
            V::Float(h.t),
            floats3((h.point.x, h.point.y, h.point.z)),
            floats3((h.normal.x, h.normal.y, h.normal.z)),
        ])
    })
}

fn ray(
    g: &Graph,
    origin: NodeId,
    dir: NodeId,
) -> Option<cg::Ray> {
    Some(cg::Ray::new(point3(vec3(g, origin)?), vector3(vec3(g, dir)?)))
}

fn ray_sphere(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [o, d, c, r] = arg(a)?;
    let g = &*cx.graph;
    let sphere = cg::Sphere::new(point3(vec3(g, c)?), float(g, r)?);
    Some(hit(cg::ray_sphere_intersection(&ray(g, o, d)?, &sphere)))
}

fn ray_plane(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [o, d, p, n] = arg(a)?;
    let g = &*cx.graph;
    let plane = cg::Plane::new(point3(vec3(g, p)?), vector3(vec3(g, n)?));
    Some(hit(cg::ray_plane_intersection(&ray(g, o, d)?, &plane)))
}

fn ray_triangle(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [o, d, v0, v1, v2] = arg(a)?;
    let g = &*cx.graph;
    let (v0, v1, v2) = (point3(vec3(g, v0)?), point3(vec3(g, v1)?), point3(vec3(g, v2)?));
    Some(hit(cg::ray_triangle_intersection(&ray(g, o, d)?, &v0, &v1, &v2)))
}

fn refract(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [d, n, eta] = arg(a)?;
    let g = &*cx.graph;
    let out = cg::refract(&vector3(vec3(g, d)?), &vector3(vec3(g, n)?), float(g, eta)?);
    Some(out.map_or(V::Bool(false), |v| floats3((v.x, v.y, v.z))))
}

fn barycentric(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [p, v0, v1, v2] = arg(a)?;
    let g = &*cx.graph;
    let (u, v, w) = cg::barycentric_coordinates(&point3(vec3(g, p)?), &point3(vec3(g, v0)?), &point3(vec3(g, v1)?), &point3(vec3(g, v2)?));
    Some(V::List(vec![V::Float(u), V::Float(v), V::Float(w)]))
}

fn quat_slerp(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [p, q, t] = arg(a)?;
    let g = &*cx.graph;
    let read = |n: NodeId| -> Option<cg::Quaternion> {
        let v: Vec<f64> = items(g, n)?.into_iter().map(|c| float(g, c)).collect::<Option<_>>()?;
        let [w, x, y, z] = <[f64; 4]>::try_from(v).ok()?;
        Some(cg::Quaternion::new(w, x, y, z))
    };
    let r = read(p)?.slerp(&read(q)?, float(g, t)?);
    Some(V::List(vec![V::Float(r.w), V::Float(r.x), V::Float(r.y), V::Float(r.z)]))
}

/// `d - 2 (d . n) n`.
fn reflect(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [d, n] = arg(a)?;
    let (d, n) = (point(cx.graph, d)?, point(cx.graph, n)?);
    if d.len() != n.len() {
        return None;
    }
    let k = dot(cx, &d, &n);
    let one = int(cx, 1);
    let out: Row = d
        .iter()
        .zip(&n)
        .map(|(&di, &ni)| combine(cx, &[(1, di, one), (-2, k, ni)]))
        .collect();
    Some(V::nodes(&out))
}

/// Registers the graphics operators.
pub(crate) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    def(i, "translation_2d", Arity::Fixed(2), translation)?;
    def(i, "translation_3d", Arity::Fixed(3), translation)?;
    def(i, "scaling_2d", Arity::Fixed(2), diagonal)?;
    def(i, "scaling_3d", Arity::Fixed(3), diagonal)?;
    def(i, "shear_2d", Arity::Fixed(2), shear_2d)?;
    def(i, "rotation_2d", Arity::Fixed(1), |cx, a| {
        let [t] = arg(a)?;
        plane_rotation(cx, 3, t, (0, 1), false)
    })?;
    def(i, "rotation_3d_x", Arity::Fixed(1), |cx, a| {
        let [t] = arg(a)?;
        plane_rotation(cx, 4, t, (1, 2), false)
    })?;
    def(i, "rotation_3d_y", Arity::Fixed(1), |cx, a| {
        let [t] = arg(a)?;
        plane_rotation(cx, 4, t, (0, 2), true)
    })?;
    def(i, "rotation_3d_z", Arity::Fixed(1), |cx, a| {
        let [t] = arg(a)?;
        plane_rotation(cx, 4, t, (0, 1), false)
    })?;
    def(i, "rotation_axis_angle", Arity::Fixed(2), rotation_axis_angle)?;
    def(i, "reflection_2d", Arity::Fixed(1), reflection_2d)?;
    def(i, "reflection_3d", Arity::Fixed(1), reflection_3d)?;
    def(i, "perspective", Arity::Fixed(4), perspective)?;
    def(i, "orthographic", Arity::Fixed(6), orthographic)?;
    def(i, "look_at", Arity::Fixed(3), look_at)?;
    def(i, "apply_transform", Arity::Fixed(2), |cx, a| transform(cx, a, true))?;
    def(i, "apply_transform_vector", Arity::Fixed(2), |cx, a| transform(cx, a, false))?;
    def(i, "bezier", Arity::Fixed(2), bezier)?;
    def(i, "bezier_derivative", Arity::Fixed(2), bezier_derivative)?;
    def(i, "bezier_split", Arity::Fixed(2), bezier_split)?;
    def(i, "bspline", Arity::Fixed(4), bspline)?;
    def(i, "catmull_rom", Arity::Fixed(5), catmull_rom)?;
    def(i, "quat_mul", Arity::Fixed(2), quat_binary)?;
    def(i, "quat_conj", Arity::Fixed(1), quat_conj)?;
    def(i, "quat_inverse", Arity::Fixed(1), quat_inverse)?;
    def_request(i, "quat_norm", Arity::Fixed(1), quat_norm)?;
    def(i, "quat_normalize", Arity::Fixed(1), quat_normalize)?;
    def(i, "quat_from_axis_angle", Arity::Fixed(2), quat_from_axis_angle)?;
    def(i, "quat_rotate", Arity::Fixed(2), quat_rotate)?;
    def(i, "quat_to_matrix", Arity::Fixed(1), quat_to_matrix)?;
    def(i, "quat_slerp", Arity::Fixed(3), quat_slerp)?;
    def(i, "mesh_transform", Arity::Fixed(2), mesh_transform)?;
    def(i, "mesh_normals", Arity::Fixed(2), mesh_normals)?;
    def(i, "mesh_triangulate", Arity::Fixed(1), mesh_triangulate)?;
    def(i, "ray_sphere", Arity::Fixed(4), ray_sphere)?;
    def(i, "ray_plane", Arity::Fixed(4), ray_plane)?;
    def(i, "ray_triangle", Arity::Fixed(5), ray_triangle)?;
    def(i, "reflect", Arity::Fixed(2), reflect)?;
    def(i, "refract", Arity::Fixed(3), refract)?;
    def(i, "barycentric", Arity::Fixed(4), barycentric)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::super::test_util::s;

    /// Every number in `text`, in order.
    fn floats(text: &str) -> Vec<f64> {
        text.split(|c: char| !(c.is_ascii_digit() || matches!(c, '.' | '-' | 'e')))
            .filter_map(|t| t.parse().ok())
            .collect()
    }

    fn close(
        text: &str,
        expected: &[f64],
    ) {
        let got = floats(text);
        assert_eq!(got.len(), expected.len(), "{text}");
        for (g, e) in got.iter().zip(expected) {
            assert!((g - e).abs() < 1e-9, "{text}: {g} vs {e}");
        }
    }

    #[test]
    fn translation_scaling_and_shear() {
        assert_eq!(s("translation_2d(a, b)"), "list(list(1, 0, a), list(0, 1, b), list(0, 0, 1))");
        assert_eq!(s("translation_3d(1, 2, 3)"), "list(list(1, 0, 0, 1), list(0, 1, 0, 2), list(0, 0, 1, 3), list(0, 0, 0, 1))");
        assert_eq!(s("scaling_2d(2, 3)"), "list(list(2, 0, 0), list(0, 3, 0), list(0, 0, 1))");
        assert_eq!(s("scaling_3d(2, 3, 4)"), "list(list(2, 0, 0, 0), list(0, 3, 0, 0), list(0, 0, 4, 0), list(0, 0, 0, 1))");
        assert_eq!(s("shear_2d(a, b)"), "list(list(1, a, 0), list(b, 1, 0), list(0, 0, 1))");
        assert_eq!(s("apply_transform(translation_2d(1, 2), list(3, 4))"), "list(4, 6)");
        assert_eq!(s("apply_transform_vector(translation_2d(1, 2), list(3, 4))"), "list(3, 4)");
    }

    #[test]
    fn rotations_stay_symbolic_and_fold_for_numbers() {
        assert_eq!(s("rotation_2d(t)"), "list(list(cos(t), -sin(t), 0), list(sin(t), cos(t), 0), list(0, 0, 1))");
        assert_eq!(s("rotation_2d(pi/2)"), "list(list(0, -1, 0), list(1, 0, 0), list(0, 0, 1))");
        assert_eq!(s("apply_transform(rotation_2d(pi/2), list(1, 0, 1))"), "list(0, 1, 1)");
        assert_eq!(s("apply_transform(rotation_2d(pi/2), list(1, 0))"), "list(0, 1)");
        assert_eq!(s("rotation_3d_x(t)"), "list(list(1, 0, 0, 0), list(0, cos(t), -sin(t), 0), list(0, sin(t), cos(t), 0), list(0, 0, 0, 1))");
        assert_eq!(s("rotation_3d_y(pi/2)"), "list(list(0, 0, 1, 0), list(0, 1, 0, 0), list(-1, 0, 0, 0), list(0, 0, 0, 1))");
        assert_eq!(s("rotation_3d_z(pi/2)"), "list(list(0, -1, 0, 0), list(1, 0, 0, 0), list(0, 0, 1, 0), list(0, 0, 0, 1))");
        assert_eq!(s("rotation_3d_z(0)"), "list(list(1, 0, 0, 0), list(0, 1, 0, 0), list(0, 0, 1, 0), list(0, 0, 0, 1))");
    }

    #[test]
    fn rodrigues_matches_the_axis_rotations() {
        assert_eq!(s("rotation_axis_angle(list(0, 0, 1), pi/2)"), s("rotation_3d_z(pi/2)"));
        assert_eq!(s("rotation_axis_angle(list(0, 1, 0), t)"), s("rotation_3d_y(t)"));
        assert_eq!(s("rotation_axis_angle(list(2, 0, 0), t)"), s("rotation_3d_x(t)"), "the axis is normalised");
        assert_eq!(s("apply_transform(rotation_axis_angle(list(0, 0, 1), pi), list(1, 2, 3))"), "list(-1, -2, 3)");
    }

    #[test]
    fn reflections() {
        assert_eq!(s("reflection_2d(pi/4)"), "list(list(0, 1, 0), list(1, 0, 0), list(0, 0, 1))");
        assert_eq!(s("reflection_2d(0)"), "list(list(1, 0, 0), list(0, -1, 0), list(0, 0, 1))");
        assert_eq!(s("reflection_3d(list(0, 0, 1))"), "list(list(1, 0, 0, 0), list(0, 1, 0, 0), list(0, 0, -1, 0), list(0, 0, 0, 1))");
        assert_eq!(s("reflection_3d(list(0, 0, 2))"), s("reflection_3d(list(0, 0, 1))"), "non-unit normals are normalised");
        assert_eq!(s("reflection_3d(list(1, 1, 0))"), "list(list(0, -1, 0, 0), list(-1, 0, 0, 0), list(0, 0, 1, 0), list(0, 0, 0, 1))");
        assert_eq!(s("apply_transform(reflection_3d(list(1, 0, 0)), list(2, 3, 4))"), "list(-2, 3, 4)");
        assert_eq!(s("reflect(list(1, -1, 0), list(0, 1, 0))"), "list(1, 1, 0)");
        assert_eq!(s("reflect(list(a, b), list(0, 1))"), "list(a, -b)");
    }

    #[test]
    fn projections_and_the_camera() {
        assert_eq!(s("perspective(pi/2, 1, 1, 3)"), "list(list(1, 0, 0, 0), list(0, 1, 0, 0), list(0, 0, -2, -3), list(0, 0, -1, 0))");
        let wide = s("perspective(pi/2, a, 1, 3)");
        assert!(wide.starts_with("list(list(1/a, 0, 0, 0)"), "{wide}");
        assert_eq!(s("orthographic(-1, 1, -1, 1, 1, 3)"), "list(list(1, 0, 0, 0), list(0, 1, 0, 0), list(0, 0, -1, -2), list(0, 0, 0, 1))");
        assert_eq!(s("look_at(list(0, 0, 5), list(0, 0, 0), list(0, 1, 0))"), "list(list(1, 0, 0, 0), list(0, 1, 0, 0), list(0, 0, 1, -5), list(0, 0, 0, 1))");
        // The eye maps to the origin, the centre to the negative z axis.
        assert_eq!(s("apply_transform(look_at(list(1, 2, 3), list(1, 2, 0), list(0, 1, 0)), list(1, 2, 3))"), "list(0, 0, 0)");
        assert_eq!(s("apply_transform(look_at(list(1, 2, 3), list(1, 2, 0), list(0, 1, 0)), list(1, 2, 0))"), "list(0, 0, -3)");
        // Perspective divide.
        assert_eq!(s("apply_transform(perspective(pi/2, 1, 1, 3), list(0, 0, -2))"), "list(0, 0, 1/2)");
    }

    #[test]
    fn matrices_compose_with_matmul() {
        assert_eq!(s("matmul(translation_2d(1, 0), rotation_2d(pi/2))"), "list(list(0, -1, 1), list(1, 0, 0), list(0, 0, 1))");
    }

    #[test]
    fn bezier_curves() {
        assert_eq!(s("bezier(list(list(0, 0), list(1, 2), list(2, 0)), 1/2)"), "list(1, 1)");
        assert_eq!(s("bezier(list(list(0, 0), list(1, 2), list(2, 0)), 0)"), "list(0, 0)");
        assert_eq!(s("bezier(list(list(0, 0), list(1, 2), list(2, 0)), 1)"), "list(2, 0)");
        assert_eq!(s("bezier(list(list(0), list(1)), t)"), "list(t)");
        assert_eq!(s("bezier(list(list(0, 0), list(1, 2), list(2, 0)), t)"), "list(2*t, 4*t - 4*t^2)");
        assert_eq!(s("bezier(list(list(5, 5)), t)"), "list(5, 5)");
        assert_eq!(s("bezier_derivative(list(list(0, 0), list(1, 2), list(2, 0)), 0)"), "list(2, 4)");
        assert_eq!(s("bezier_derivative(list(list(0, 0), list(1, 2), list(2, 0)), t)"), "list(2, 4 - 8*t)");
        assert_eq!(s("bezier_derivative(list(list(1, 1)), t)"), "list(0, 0)");
        let split = "bezier_split(list(list(0, 0), list(2, 2), list(4, 0)), 1/2)";
        assert_eq!(s(split), "list(list(list(0, 0), list(1, 1), list(2, 1)), list(list(2, 1), list(3, 1), list(4, 0)))");
        // The halves reproduce the curve.
        assert_eq!(s("bezier(list(list(0, 0), list(1, 1), list(2, 1)), 1)"), s("bezier(list(list(0, 0), list(2, 2), list(4, 0)), 1/2)"));
    }

    #[test]
    fn bsplines_and_catmull_rom() {
        let cps = "list(list(0, 0), list(1, 1), list(2, 0), list(3, 1))";
        let knots = "list(0, 0, 0, 1, 2, 2, 2)";
        assert_eq!(s(&format!("bspline({cps}, 2, {knots}, 0)")), "list(0, 0)");
        assert_eq!(s(&format!("bspline({cps}, 2, {knots}, 2)")), "list(3, 1)");
        assert_eq!(s(&format!("bspline({cps}, 2, {knots}, 1)")), "list(3/2, 1/2)");
        // Degree 1 interpolates linearly.
        assert_eq!(s("bspline(list(list(0), list(2), list(2)), 1, list(0, 0, 1, 2, 2), 1/2)"), "list(1)");
        // Outside the domain, or with the wrong number of knots: unreduced.
        assert!(s(&format!("bspline({cps}, 2, {knots}, 3)")).starts_with("bspline("));
        assert!(s(&format!("bspline({cps}, 2, list(0, 1), 1)")).starts_with("bspline("));
        assert_eq!(s("catmull_rom(list(0, 0), list(1, 0), list(2, 1), list(3, 1), 0)"), "list(1, 0)");
        assert_eq!(s("catmull_rom(list(0, 0), list(1, 0), list(2, 1), list(3, 1), 1)"), "list(2, 1)");
        assert_eq!(s("catmull_rom(list(0, 0), list(1, 0), list(2, 1), list(3, 1), 1/2)"), "list(3/2, 1/2)");
        assert_eq!(s("catmull_rom(list(0), list(1), list(2), list(3), t)"), "list(t + 1)");
    }

    #[test]
    fn quaternion_algebra() {
        assert_eq!(s("quat_mul(list(1, 2, 3, 4), list(5, 6, 7, 8))"), "list(-60, 12, 30, 24)");
        assert_eq!(s("quat_mul(list(0, 1, 0, 0), list(0, 0, 1, 0))"), "list(0, 0, 0, 1)", "i j = k");
        assert_eq!(s("quat_conj(list(1, 2, 3, 4))"), "list(1, -2, -3, -4)");
        assert_eq!(s("quat_norm(list(1, 2, 2, 4))"), "5");
        assert_eq!(s("quat_norm(list(a, 0, 0, 0))"), "sqrt(a^2)");
        assert_eq!(s("quat_normalize(list(0, 3, 0, 4))"), "list(0, 3/5, 0, 4/5)");
        assert_eq!(s("quat_inverse(list(0, 0, 0, 2))"), "list(0, 0, 0, -1/2)");
        assert_eq!(s("quat_mul(list(1, 2, 2, 4), quat_inverse(list(1, 2, 2, 4)))"), "list(1, 0, 0, 0)");
        assert_eq!(s("quat_mul(list(w, x, y, z), quat_conj(list(w, x, y, z)))"), "list(w^2 + x^2 + y^2 + z^2, 0, 0, 0)");
    }

    #[test]
    fn quaternion_rotations() {
        assert_eq!(s("quat_from_axis_angle(list(0, 0, 1), pi)"), "list(0, 0, 0, 1)");
        assert_eq!(s("quat_from_axis_angle(list(0, 0, 3), 0)"), "list(1, 0, 0, 0)");
        assert_eq!(s("quat_rotate(quat_from_axis_angle(list(0, 0, 1), pi/2), list(1, 0, 0))"), "list(0, 1, 0)");
        assert_eq!(s("quat_to_matrix(list(1, 0, 0, 0))"), "list(list(1, 0, 0, 0), list(0, 1, 0, 0), list(0, 0, 1, 0), list(0, 0, 0, 1))");
        assert_eq!(s("quat_to_matrix(quat_from_axis_angle(list(0, 0, 1), pi))"), "list(list(-1, 0, 0, 0), list(0, -1, 0, 0), list(0, 0, 1, 0), list(0, 0, 0, 1))");
        // Non-unit quaternions are normalised by the matrix.
        assert_eq!(s("quat_to_matrix(list(0, 0, 0, 5))"), s("quat_to_matrix(list(0, 0, 0, 1))"));
        close(&s("quat_slerp(list(1, 0, 0, 0), list(0, 1, 0, 0), 0.5)"), &[0.5_f64.sqrt(), 0.5_f64.sqrt(), 0.0, 0.0]);
        close(&s("quat_slerp(list(1, 0, 0, 0), list(0, 1, 0, 0), 0)"), &[1.0, 0.0, 0.0, 0.0]);
        assert!(s("quat_slerp(list(1, 0, 0, 0), list(0, 1, 0, 0), t)").starts_with("quat_slerp("));
    }

    #[test]
    fn meshes() {
        assert_eq!(s("mesh_transform(list(list(1, 0, 0), list(0, 1, 0)), translation_3d(1, 2, 3))"), "list(list(2, 2, 3), list(1, 3, 3))");
        assert_eq!(s("mesh_transform(list(list(1, 0, 0)), rotation_3d_z(pi/2))"), "list(list(0, 1, 0))");
        let tri = "list(list(0, 0, 0), list(1, 0, 0), list(0, 1, 0))";
        assert_eq!(s(&format!("mesh_normals({tri}, list(list(0, 1, 2)))")), "list(list(0, 0, 1))");
        assert_eq!(s(&format!("mesh_normals({tri}, list(list(0, 2, 1), list(0, 1)))")), "list(list(0, 0, -1))", "short faces are skipped");
        assert_eq!(s("mesh_triangulate(list(list(0, 1, 2, 3), list(4, 5, 6)))"), "list(list(0, 1, 2), list(0, 2, 3), list(4, 5, 6))");
        assert_eq!(s("mesh_triangulate(list(list(0, 1, 2, 3, 4)))"), "list(list(0, 1, 2), list(0, 2, 3), list(0, 3, 4))");
    }

    #[test]
    fn ray_queries() {
        assert_eq!(s("ray_sphere(list(0, 0, -5), list(0, 0, 1), list(0, 0, 0), 1)"), "list(4, list(0, 0, -1), list(0, 0, -1))");
        assert_eq!(s("ray_sphere(list(0, 5, -5), list(0, 0, 1), list(0, 0, 0), 1)"), "false");
        assert_eq!(s("ray_sphere(list(0, 0, 5), list(0, 0, 1), list(0, 0, 0), 1)"), "false", "the sphere is behind the ray");
        assert_eq!(s("ray_plane(list(0, 0, 5), list(0, 0, -1), list(0, 0, 0), list(0, 0, 1))"), "list(5, list(0, 0, 0), list(0, 0, 1))");
        assert_eq!(s("ray_plane(list(0, 0, 5), list(1, 0, 0), list(0, 0, 0), list(0, 0, 1))"), "false", "parallel");
        let tri = "list(0, 0, 0), list(1, 0, 0), list(0, 1, 0)";
        assert_eq!(s(&format!("ray_triangle(list(0.25, 0.25, 1), list(0, 0, -1), {tri})")), "list(1, list(0.25, 0.25, 0), list(0, 0, 1))");
        assert_eq!(s(&format!("ray_triangle(list(2, 2, 1), list(0, 0, -1), {tri})")), "false");
        assert!(s("ray_sphere(list(a, 0, 0), list(0, 0, 1), list(0, 0, 0), 1)").starts_with("ray_sphere("));
    }

    #[test]
    fn refraction_and_barycentric_coordinates() {
        close(&s("refract(list(0, -1, 0), list(0, 1, 0), 0.5)"), &[0.0, -1.0, 0.0]);
        close(&s("refract(list(1, -1, 0), list(0, 1, 0), 1)"), &[0.5_f64.sqrt(), -0.5_f64.sqrt(), 0.0]);
        assert_eq!(s("refract(list(1, -0.1, 0), list(0, 1, 0), 1.5)"), "false", "total internal reflection");
        close(&s("barycentric(list(0.25, 0.25), list(0, 0), list(1, 0), list(0, 1))"), &[0.5, 0.25, 0.25]);
        close(&s("barycentric(list(1, 0, 0), list(0, 0, 0), list(1, 0, 0), list(0, 1, 0))"), &[0.0, 1.0, 0.0]);
    }

    #[test]
    fn malformed_arguments_stay_unreduced() {
        assert!(s("rotation_2d(t)").starts_with("list("));
        assert!(s("apply_transform(rotation_2d(t), list(1, 2, 3, 4))").starts_with("apply_transform("));
        assert!(s("quat_inverse(list(1, 2, 3))").starts_with("quat_inverse("));
        assert!(s("bezier(list(list(0, 0), list(1)), t)").starts_with("bezier("));
        assert!(s("look_at(list(0, 0), list(0, 0, 0), list(0, 1, 0))").starts_with("look_at("));
    }
}
