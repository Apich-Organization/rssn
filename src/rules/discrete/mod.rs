//! Discrete mathematics: cryptography, coding theory, finite fields, graphs,
//! algebraic topology, fractals and chaos, computer graphics and finite
//! groups.
//!
//! These peripheral domains are utility wrappers around the graph engine:
//! their inputs and outputs are ordinary terms (numbers, symbols, lists), so
//! they compose with the rest of the system. A graph's adjacency matrix is a
//! `list(list(...))` that `det`, `charpoly` or `eigenvals` of the linear
//! algebra rules accept directly; a transformation matrix of the graphics
//! module multiplies with `matmul`.
//!
//! All operators are *pure functions of their arguments*: nothing is
//! random. Where the legacy code drew a random number (key pairs, ECDSA
//! nonces) the value is now an explicit argument.
//!
//! Every operator is reduced by a kernel that fires when its arguments have
//! the expected literal shape; anything else (a symbolic modulus, a
//! malformed graph, a failing precondition) stays unreduced.
//!
//! # Conventions
//!
//! * Polynomials over a finite field are coefficient lists, **highest degree
//!   first**, without leading zeros (the zero polynomial is `list()`).
//! * "No result" answers are the boolean `false` (no bipartition, no
//!   topological order, negative cycle).
//!
//! # Operators
//!
//! | module | operators |
//! |---|---|
//! | [`crypto`] | `ec_curve`, `ec_on_curve`, `ec_neg`, `ec_double`, `ec_add`, `ec_mul`, `ec_order`, `ec_is_infinity`, `ec_x`, `ec_y`, `ecdh_public`, `ecdh_shared`, `ec_compress`, `ec_decompress`, `ecdsa_sign`, `ecdsa_verify`, `rsa_keygen`, `rsa_encrypt`, `rsa_decrypt` |
//! | [`coding`] | `hamming_distance`, `hamming_weight`, `hamming_encode`, `hamming_check`, `hamming_decode`, `rs_encode`, `rs_check`, `rs_decode`, `rs_error_count`, `bch_encode`, `bch_decode`, `crc32`, `crc32_verify`, `crc32_update`, `crc32_finalize`, `crc16`, `crc8`, `gf256_*`, `gf256_poly_*` |
//! | [`finite_field`] | `gf_*` (prime field), `gfp_*` (polynomials over GF(p)), `gfx_*` (extension fields GF(p)\[x\]/(m)) |
//! | [`graphs`] | `graph`/`digraph` terms, `graph_*` queries, traversals, shortest paths, spanning trees, flows, matchings, colouring, isomorphism, products |
//! | [`topology`] | `sc_*` simplicial complexes |
//! | [`fractal`] | `mandelbrot_escape`, `julia_escape`, `burning_ship_escape`, `multibrot_escape`, `newton_fractal_root`, `mandelbrot_iterate`, `mandelbrot_orbit`, `mandelbrot_fixed_points`, `mandelbrot_stability`, `complex_map_fixed_points`, `complex_map_stability`, `map_fixed_points`, `map_stability`, `lyapunov_exponent`, `logistic_iterate`, `logistic_bifurcation`, `logistic_lyapunov`, `lorenz`, `lorenz_orbit`, `lorenz_lyapunov`, `rossler_orbit`, `henon_orbit`, `tinkerbell_orbit`, `ifs_apply`, `ifs_generate`, `similarity_dimension`, `moran_dimension`, `box_counting`, `correlation_dimension`, `orbit_density`, `orbit_entropy` |
//! | [`graphics`] | `translation_2d/3d`, `scaling_2d/3d`, `shear_2d`, `rotation_2d`, `rotation_3d_x/y/z`, `rotation_axis_angle`, `reflection_2d/3d`, `perspective`, `orthographic`, `look_at`, `transform_point`, `transform_vector`, `bezier`, `bezier_derivative`, `bezier_split`, `bspline`, `catmull_rom`, `quat_mul`, `quat_conj`, `quat_inverse`, `quat_norm`, `quat_normalize`, `quat_from_axis_angle`, `quat_rotate`, `quat_to_matrix`, `quat_slerp`, `mesh_transform`, `mesh_normals`, `mesh_triangulate`, `ray_sphere`, `ray_plane`, `ray_triangle`, `reflect`, `refract`, `barycentric` |
//! | [`groups`] | `cyclic_group`, `dihedral_group`, `symmetric_group`, `klein_four_group`, `group_*`, `perm_*` |
//!
//! See the documentation of each module for the exact term formats.

#![allow(dead_code)] // TEMP

pub mod coding;
pub mod crypto;
pub mod finite_field;
pub mod fractal;
pub mod graphics;
pub mod graphs;
pub mod groups;
pub mod topology;

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::ToPrimitive;
use num_traits::Zero;

use crate::graph::op::core;
use crate::graph::rule::Installer;
use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpDescriptor;
use crate::graph::OpFlags;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::Payload;
use crate::graph::RuleError;
use crate::graph::RuleSet;
use crate::graph::Tier;

use super::elementary::elementary;
use super::linalg::linalg;
use super::number_theory::number_theory;

/// The discrete mathematics rule set.
#[must_use]
pub fn discrete() -> RuleSet {
    RuleSet::new("discrete", install)
        .needs(elementary())
        .needs(number_theory())
        .needs(linalg())
}

fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    crypto::install(i)?;
    finite_field::install(i)?;
    coding::install(i)?;
    graphs::install(i)?;
    topology::install(i)?;
    fractal::install(i)?;
    graphics::install(i)?;
    groups::install(i)?;
    Ok(())
}

// ----------------------------------------------------------------------
// Values and kernels
// ----------------------------------------------------------------------

/// What a discrete kernel computed.
pub(crate) enum V {
    /// An already built term.
    Node(NodeId),
    /// An exact integer.
    Int(BigInt),
    /// A truth value.
    Bool(bool),
    /// A float.
    Float(f64),
    /// `list(...)`.
    List(Vec<V>),
}

impl V {
    pub(crate) fn int(v: i64) -> Self {
        Self::Int(BigInt::from(v))
    }

    pub(crate) fn uint(v: usize) -> Self {
        Self::Int(BigInt::from(v))
    }

    pub(crate) fn ints<T: Into<BigInt>>(items: impl IntoIterator<Item = T>) -> Self {
        Self::List(items.into_iter().map(|v| Self::Int(v.into())).collect())
    }

    pub(crate) fn nodes(items: &[NodeId]) -> Self {
        Self::List(items.iter().map(|&n| Self::Node(n)).collect())
    }

    pub(crate) fn build(
        self,
        graph: &mut Graph,
    ) -> NodeId {
        match self {
            | Self::Node(n) => n,
            | Self::Int(v) => graph.num(Number::Int(v)),
            | Self::Bool(b) => graph.lit(Payload::Bool(b)),
            | Self::Float(x) => graph.float(x),
            | Self::List(items) => {
                let nodes: Vec<NodeId> = items.into_iter().map(|v| v.build(graph)).collect();
                graph.node(core::LIST, &nodes)
            },
        }
    }
}

/// A computation on the argument nodes of an operator.
pub(crate) type Run = fn(&mut Cx<'_>, &[NodeId]) -> Option<V>;

/// Reduces `op` with a [`Run`] function.
struct Fun {
    op: OpId,
    run: Run,
    /// Lists are larger than the request; pin them so that extraction
    /// returns them.
    pin: bool,
}

impl Kernel for Fun {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let args = cx.graph.children(node).to_vec();
        let Some(value) = (self.run)(cx, &args) else {
            return Outcome::Pass;
        };
        let result = value.build(cx.graph);
        let head = cx.graph.op(result);
        let structured = head == core::LIST || matches!(&*cx.graph.ops().get(head).name, "graph" | "digraph" | "group");
        if self.pin && structured {
            Outcome::Pinned(result)
        } else {
            Outcome::Equal(result)
        }
    }

    fn revisit(&self) -> bool {
        true
    }
}

/// Registers an operator with a kernel; list results are pinned.
pub(crate) fn def(
    i: &mut Installer<'_>,
    name: &str,
    arity: Arity,
    run: Run,
) -> Result<OpId, RuleError> {
    let op = i.op(OpDescriptor::new(name, arity))?;
    i.kernel(&format!("discrete/{name}"), Tier::Normalize, Fun { op, run, pin: true });
    Ok(op)
}

/// Registers a *request* operator: it stands for a computation whose answer
/// is a term still to be simplified (a formula with `sin`, `ln`, ...), so
/// it is flagged heavy and the answer is not pinned.
pub(crate) fn def_request(
    i: &mut Installer<'_>,
    name: &str,
    arity: Arity,
    run: Run,
) -> Result<OpId, RuleError> {
    let op = i.op(OpDescriptor::new(name, arity).flags(OpFlags::HEAVY).cost(100))?;
    i.kernel(&format!("discrete/{name}"), Tier::Reduce, Fun { op, run, pin: false });
    Ok(op)
}

/// Registers a data constructor with no kernel (an inert term).
pub(crate) fn def_inert(
    i: &mut Installer<'_>,
    name: &str,
    arity: Arity,
) -> Result<OpId, RuleError> {
    i.op(OpDescriptor::new(name, arity))
}

// ----------------------------------------------------------------------
// Reading arguments
// ----------------------------------------------------------------------

/// A literal integer.
pub(crate) fn big(
    g: &Graph,
    n: NodeId,
) -> Option<BigInt> {
    match g.number_of(n)? {
        | Number::Int(v) => Some(v.clone()),
        | _ => None,
    }
}

/// A literal integer that fits `i64`.
pub(crate) fn small(
    g: &Graph,
    n: NodeId,
) -> Option<i64> {
    big(g, n)?.to_i64()
}

/// A literal non-negative integer.
pub(crate) fn idx(
    g: &Graph,
    n: NodeId,
) -> Option<usize> {
    big(g, n)?.to_usize()
}

/// A literal number as a float.
pub(crate) fn float(
    g: &Graph,
    n: NodeId,
) -> Option<f64> {
    g.number_of(n).map(Number::to_f64)
}

/// The children of a `list`.
pub(crate) fn items(
    g: &Graph,
    n: NodeId,
) -> Option<Vec<NodeId>> {
    if g.op(n) == core::LIST {
        return Some(g.children(n).to_vec());
    }
    // The argument may be a request that has since been reduced to a list.
    g.enodes(g.find(n))
        .find(|&e| g.op(e) == core::LIST)
        .map(|e| g.children(e).to_vec())
}

/// A `list` of literal integers.
pub(crate) fn bigs(
    g: &Graph,
    n: NodeId,
) -> Option<Vec<BigInt>> {
    items(g, n)?.into_iter().map(|c| big(g, c)).collect()
}

/// A `list` of literal non-negative integers.
pub(crate) fn idxs(
    g: &Graph,
    n: NodeId,
) -> Option<Vec<usize>> {
    items(g, n)?.into_iter().map(|c| idx(g, c)).collect()
}

/// A `list` of bytes (integers in `0..=255`).
pub(crate) fn bytes(
    g: &Graph,
    n: NodeId,
) -> Option<Vec<u8>> {
    items(g, n)?
        .into_iter()
        .map(|c| small(g, c).and_then(|v| u8::try_from(v).ok()))
        .collect()
}

/// A `list` of floats.
pub(crate) fn floats(
    g: &Graph,
    n: NodeId,
) -> Option<Vec<f64>> {
    items(g, n)?.into_iter().map(|c| float(g, c)).collect()
}

/// A list of lists.
pub(crate) fn rows(
    g: &Graph,
    n: NodeId,
) -> Option<Vec<Vec<NodeId>>> {
    items(g, n)?.into_iter().map(|r| items(g, r)).collect()
}

/// Bytes as a list value.
pub(crate) fn byte_list(data: &[u8]) -> V {
    V::ints(data.iter().map(|&b| i64::from(b)))
}

// ----------------------------------------------------------------------
// Building terms
// ----------------------------------------------------------------------

/// Whether `n` is a literal zero.
pub(crate) fn is_zero(
    g: &Graph,
    n: NodeId,
) -> bool {
    g.number_of(n).is_some_and(Number::is_zero)
}

/// Whether `n` is a literal one.
pub(crate) fn is_one(
    g: &Graph,
    n: NodeId,
) -> bool {
    g.number_of(n).is_some_and(Number::is_one)
}

/// The sum of `terms`; literal numbers are folded.
pub(crate) fn sum(
    g: &mut Graph,
    terms: &[NodeId],
) -> NodeId {
    let mut acc: Option<Number> = None;
    let mut rest = Vec::new();
    for &t in terms {
        match g.number_of(t).cloned() {
            | Some(v) => {
                acc = Some(match acc {
                    | Some(a) => a.add(&v),
                    | None => v,
                });
            },
            | None => rest.push(t),
        }
    }
    let mut all = Vec::new();
    if let Some(a) = acc {
        if !a.is_zero() || rest.is_empty() {
            all.push(g.num(a));
        }
    }
    all.extend(rest);
    match all.as_slice() {
        | [] => g.int(0),
        | [only] => *only,
        | _ => g.node(core::ADD, &all),
    }
}

/// The product of `factors`; literal numbers are folded.
pub(crate) fn prod(
    g: &mut Graph,
    factors: &[NodeId],
) -> NodeId {
    let mut acc = Number::from(1_i64);
    let mut rest = Vec::new();
    for &f in factors {
        match g.number_of(f).cloned() {
            | Some(v) => acc = acc.mul(&v),
            | None => rest.push(f),
        }
    }
    if acc.is_zero() {
        return g.int(0);
    }
    let mut all = Vec::new();
    if !acc.is_one() || rest.is_empty() {
        all.push(g.num(acc));
    }
    all.extend(rest);
    match all.as_slice() {
        | [only] => *only,
        | _ => g.node(core::MUL, &all),
    }
}

/// `-n`.
pub(crate) fn neg(
    g: &mut Graph,
    n: NodeId,
) -> NodeId {
    let minus_one = g.int(-1);
    prod(g, &[minus_one, n])
}

/// `n^e`.
pub(crate) fn pow(
    g: &mut Graph,
    n: NodeId,
    e: NodeId,
) -> NodeId {
    g.node(core::POW, &[n, e])
}

/// `1 / n`.
pub(crate) fn inv(
    g: &mut Graph,
    n: NodeId,
) -> NodeId {
    let minus_one = g.int(-1);
    pow(g, n, minus_one)
}

/// `a / b`.
pub(crate) fn div(
    g: &mut Graph,
    a: NodeId,
    b: NodeId,
) -> NodeId {
    let b = inv(g, b);
    prod(g, &[a, b])
}

/// `name(args...)` for an operator registered by another rule set.
pub(crate) fn apply(
    g: &mut Graph,
    name: &str,
    args: &[NodeId],
) -> Option<NodeId> {
    let op = g.ops().lookup(name)?;
    g.try_node(op, args)
}

/// `list(items...)`.
pub(crate) fn list(
    g: &mut Graph,
    items: &[NodeId],
) -> NodeId {
    g.node(core::LIST, items)
}

/// A matrix term from rows of nodes.
pub(crate) fn matrix(
    g: &mut Graph,
    rows: &[Vec<NodeId>],
) -> NodeId {
    let rows: Vec<NodeId> = rows.iter().map(|r| list(g, r)).collect();
    list(g, &rows)
}

/// Simplifies every node with the running program.
pub(crate) fn simplified(
    cx: &mut Cx<'_>,
    nodes: &[NodeId],
) -> Vec<NodeId> {
    nodes.iter().map(|&n| cx.simplify(n)).collect()
}

/// The product of two matrices of nodes; entries are simplified.
pub(crate) fn matmul(
    cx: &mut Cx<'_>,
    a: &[Vec<NodeId>],
    b: &[Vec<NodeId>],
) -> Option<Vec<Vec<NodeId>>> {
    let inner = b.len();
    let cols = b.first().map_or(0, Vec::len);
    let mut out = Vec::with_capacity(a.len());
    for row in a {
        if row.len() != inner || b.iter().any(|r| r.len() != cols) {
            return None;
        }
        let mut new_row = Vec::with_capacity(cols);
        for j in 0..cols {
            let terms: Vec<NodeId> = (0..inner)
                .map(|k| prod(cx.graph, &[row[k], b[k][j]]))
                .collect();
            let s = sum(cx.graph, &terms);
            new_row.push(cx.simplify(s));
        }
        out.push(new_row);
    }
    Some(out)
}

// ----------------------------------------------------------------------
// Integer helpers shared by several modules
// ----------------------------------------------------------------------

/// Remainder in `[0, |m|)`.
pub(crate) fn modulo(
    a: &BigInt,
    m: &BigInt,
) -> BigInt {
    let m = m.abs();
    let r = a % &m;
    if r.is_negative() { r + m } else { r }
}

/// The inverse of `a` modulo `m`, if `gcd(a, m) = 1`.
pub(crate) fn mod_inverse(
    a: &BigInt,
    m: &BigInt,
) -> Option<BigInt> {
    if !m.is_positive() {
        return None;
    }
    let (mut r0, mut r1) = (modulo(a, m), m.clone());
    let (mut s0, mut s1) = (BigInt::one(), BigInt::zero());
    while !r1.is_zero() {
        let q = &r0 / &r1;
        let r2 = &r0 - &q * &r1;
        let s2 = &s0 - &q * &s1;
        (r0, r1) = (r1, r2);
        (s0, s1) = (s1, s2);
    }
    if r0.is_one() {
        Some(modulo(&s0, m))
    } else if m.is_one() {
        Some(BigInt::zero())
    } else {
        None
    }
}

/// An exact rational from a literal number; floats are converted exactly.
pub(crate) fn rational(n: &Number) -> Option<BigRational> {
    match n {
        | Number::Float(x) => BigRational::from_float(*x),
        | other => other.to_rational(),
    }
}

/// Helpers for the tests of the submodules.
#[cfg(test)]
pub(crate) mod test_util {
    use super::discrete;
    use crate::rules::testing::simplify;

    /// The closed form of `src`.
    pub(crate) fn s(src: &str) -> String {
        simplify(&[discrete()], src)
    }

    /// A small deterministic generator.
    pub(crate) struct Lcg(pub u64);

    impl Lcg {
        pub(crate) fn next(
            &mut self,
            bound: u64,
        ) -> u64 {
            self.0 = self.0.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
            (self.0 >> 33) % bound
        }
    }

    /// Every integer in `text`, in order.
    pub(crate) fn nums(text: &str) -> Vec<i64> {
        text.split(|c: char| !(c.is_ascii_digit() || c == '-'))
            .filter(|t| !t.is_empty() && *t != "-")
            .filter_map(|t| t.parse().ok())
            .collect()
    }

    /// The rows of a `list(list(..), list(..))`.
    pub(crate) fn nested(text: &str) -> Vec<Vec<i64>> {
        let inner = text.strip_prefix("list(").and_then(|t| t.strip_suffix(')')).unwrap_or("");
        inner.split("list(").skip(1).map(nums).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rules::testing::reduce_with;
    use crate::rules::testing::simplify;

    #[test]
    fn discrete_installs_with_the_standard_sets() {
        let s = simplify(&[discrete()], "ec_add(list(2, 3, 97), list(3, 6), list(3, 6))");
        assert_eq!(s, "list(80, 10)");
        let (text, reduced) = reduce_with(&[discrete()], "hamming_weight(list(1, 0, 1, 1))", &[]);
        assert!(reduced);
        assert_eq!(text, "3");
    }
}
