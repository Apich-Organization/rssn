//! Lie theory made general.
//!
//! Root systems and representations of the semisimple Lie algebras, the
//! structure theory of Lie algebras given by a basis, the classical matrix
//! algebras, Casimir operators, one-parameter subgroups and the matrix
//! logarithm.
//!
//! The exact algorithms live in [`crate::kernels::lie_structure`]; this
//! module reads and writes terms. It complements [`super::lie`] (brackets,
//! structure constants, the exponential map), which it needs.
//!
//! # Types, weights
//!
//! A *type* `T` is a name such as `A2`, `E8`, `G2`, `B3` (a bare symbol), or a Cartan matrix `list(list(2, -1), list(-1, 2))` (any
//! numbering of the Dynkin diagram, reducible matrices allowed where noted).
//! Cartan matrices use Bourbaki's numbering and the convention
//! `a_ij = <alpha_i, alpha_j^vee>`: `B_n` has `-2` at `(n-1, n)`, `C_n` at
//! `(n, n-1)`, `G2 = list(list(2, -1), list(-3, 2))`. Roots are written in
//! *simple-root coordinates*; weights by their *Dynkin labels*
//! `<lambda, alpha_i^vee>`, so `list(1, 1)` is the adjoint of `A2` and
//! `list(1, 0)` its fundamental representation `3`.
//!
//! | operator | value |
//! |---|---|
//! | `lie_cartan_matrix(T)` | the Cartan matrix |
//! | `lie_cartan_type(M)` | the type of a Cartan matrix as a string, `"A2+A1"` for a reducible one; unevaluated if `M` is not the Cartan matrix of a semisimple algebra |
//! | `lie_rank(T)`, `lie_dimension(T)`, `lie_num_positive_roots(T)` | rank, dimension `n + 2N`, number of positive roots |
//! | `lie_simple_roots(T)` | the simple roots in Euclidean coordinates (Bourbaki; standard Cartan matrices of the types `A..G`) |
//! | `lie_positive_roots(T)` | all positive roots, by height |
//! | `lie_highest_root(T)` | the highest root (irreducible types) |
//! | `lie_weyl_group_order(T)`, `lie_exponents(T)` | `\|W\|` and the exponents, from the heights of the roots |
//! | `lie_coxeter_number(T)`, `lie_dual_coxeter_number(T)` | `h` and `h^vee` (irreducible) |
//! | `lie_dynkin_diagram(T)` | edges `list(i, j, multiplicity)` (1-based; for a multiple edge `alpha_i` is the longer root: the arrow points to `j`) |
//! | `lie_fundamental_weights(T)` | the fundamental weights in simple-root coordinates |
//! | `lie_gram(T)` | the Gram matrix of the simple roots, long roots of squared length 2 |
//! | `dim_irrep(T, hw)` | the Weyl dimension formula: `dim_irrep(A2, list(1, 1)) = 8` |
//! | `weight_multiplicity(T, hw, w)` | the multiplicity of the weight `w` (Freudenthal) |
//! | `irrep_weights(T, hw)` | the dominant weights with multiplicities, `list(list(w, m), ...)` |
//! | `tensor_decomposition(T, hw1, hw2)` | `V(hw1) ⊗ V(hw2) = ⊕ m V(hw)` by Brauer-Klimyk, `list(list(hw, m), ...)` |
//! | `casimir_eigenvalue(T, hw)` | `(L, L + 2 rho)` with long roots of squared length 2 |
//! | `dynkin_index(T, hw)` | `dim(V) C(V) / dim(g)`: 1 for the fundamental of `su(n)` |
//!
//! # Structure theory
//!
//! A Lie algebra `L` is a list of independent matrices with rational
//! entries closed under the commutator, or `lie_sc(c)` for structure
//! constants `c[i][j][k]` (`[e_i, e_j] = sum_k c_ijk e_k`; for a basis with
//! complex entries use `lie_sc(structure_constants(basis))`). Subspaces and
//! elements are returned as coordinate vectors in the basis of `L`.
//!
//! | operator | value |
//! |---|---|
//! | `lie_sc_of(L)` | `lie_sc(...)` of a matrix basis |
//! | `lie_derived_series(L)`, `lie_lower_central_series(L)`, `lie_upper_central_series(L)` | the series as lists of bases (the upper series without the leading 0) |
//! | `lie_derived_dims(L)`, `lie_lower_central_dims(L)`, `lie_upper_central_dims(L)` | their dimensions |
//! | `lie_center(L)`, `lie_radical(L)` | bases of the center and of the radical (`[g,g]` orthogonal for the Killing form) |
//! | `lie_is_solvable(L)`, `lie_is_nilpotent(L)`, `lie_is_semisimple(L)` | the tests (series; Killing form nondegenerate) |
//! | `lie_cartan_solvable_test(L)` | Cartan's criterion `K(g, [g,g]) = 0` |
//! | `lie_levi(L)` | `list(levi factor, radical)` |
//! | `lie_cartan_subalgebra(L)` | a Cartan subalgebra (generalised zero eigenspace of a regular element) |
//! | `lie_root_decomposition(L, h)` | `list(list(list(root, basis)...), missing)`: the root spaces for rational roots and the dimension left over (non-split algebras) |
//! | `lie_is_subalgebra(L, vectors)` | closure of a span under the bracket |
//! | `lie_elements(L, vectors)` | the matrices with the given coordinates (matrix bases) |
//! | `lie_casimir_operator(L[, rho])` | `sum K^ab rho(e_a) rho(e_b)` for the Killing form `K`; `rho` is a list of matrices, one per basis element (the basis matrices themselves by default) |
//!
//! # Standard bases, one-parameter subgroups
//!
//! | operator | value |
//! |---|---|
//! | `gl_basis(n)`, `sl_basis(n)`, `so_basis(n)`, `sp_basis(n)`, `su_basis(n)` | bases of `gl(n)`, `sl(n)` (off-diagonal units, then `H_k = E_kk - E_{k+1,k+1}`), `so(n)` (`E_ij - E_ji`), `sp(2n)` (`[[A, B], [C, -A^T]]`, `B`, `C` symmetric; `n` is the half size), `su(n)` (`E_ij - E_ji`, `i(E_ij + E_ji)`, `i H_k`) |
//! | `one_parameter_subgroup(A, t)` | `exp_map(t A)`, reduced by [`super::lie`] |
//! | `matrix_log(g)` | the principal logarithm for the cases `exp_map` covers: identity, unipotent, diagonal with positive entries, `2 x 2` matrices with positive determinant and real eigenvalues or a rotation, rotations of `R^3` |

use std::collections::BTreeMap;

use num_bigint::BigInt;
use num_traits::One;
use num_traits::Signed;
use num_traits::Zero;

use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::RuleError;
use crate::graph::RuleSet;
use crate::graph::rule::Installer;
use crate::kernels::lie_structure as ks;
use crate::kernels::lie_structure::Family;
use crate::kernels::lie_structure::Lie;
use crate::kernels::lie_structure::RootSystem;
use crate::kernels::qlinalg as ql;
use crate::kernels::qlinalg::Q;
use crate::kernels::qlinalg::QMat;

use super::discrete::apply;
use super::discrete::def;
use super::discrete::def_inert;
use super::discrete::def_request;
use super::discrete::idx;
use super::discrete::items;
use super::discrete::prod;
use super::discrete::rational;
use super::discrete::rows;
use super::discrete::small;
use super::discrete::sum;
use super::discrete::V;
use super::lie::lie;

/// The Lie structure rule set.
#[must_use]
pub fn lie_structure() -> RuleSet {
    RuleSet::new("lie_structure", install).needs(lie())
}

/// Largest weight cache / orbit sizes are handled in the kernel; these
/// bound the size of what is read from terms.
const MAX_BASIS: usize = 80;
const MAX_MATRIX: usize = 40;

// ----------------------------------------------------------------------
// Reading
// ----------------------------------------------------------------------

fn read_type(
    g: &Graph,
    n: NodeId,
) -> Option<ks::Cartan> {
    if let Some(matrix) = rows(g, n) {
        return matrix.into_iter().map(|r| r.into_iter().map(|x| small(g, x)).collect::<Option<Vec<_>>>()).collect();
    }
    ks::parse_type(g.display(n).trim_matches('"'))
}

fn root_system(
    g: &Graph,
    n: NodeId,
) -> Option<RootSystem> {
    RootSystem::new(read_type(g, n)?)
}

fn labels(
    g: &Graph,
    n: NodeId,
    rank: usize,
) -> Option<Vec<i64>> {
    let l: Vec<i64> = items(g, n)?.into_iter().map(|x| small(g, x)).collect::<Option<_>>()?;
    (l.len() == rank).then_some(l)
}

fn rat(
    g: &Graph,
    n: NodeId,
) -> Option<Q> {
    g.number_of(n).and_then(rational)
}

fn rat_vec(
    g: &Graph,
    n: NodeId,
) -> Option<Vec<Q>> {
    items(g, n)?.into_iter().map(|x| rat(g, x)).collect()
}

fn rat_matrix(
    g: &Graph,
    n: NodeId,
) -> Option<QMat> {
    let m: QMat = rows(g, n)?.into_iter().map(|r| r.into_iter().map(|x| rat(g, x)).collect::<Option<Vec<_>>>()).collect::<Option<_>>()?;
    let width = m.first()?.len();
    (width > 0 && m.iter().all(|r| r.len() == width)).then_some(m)
}

fn read_constants(
    g: &Graph,
    n: NodeId,
) -> Option<Vec<Vec<Vec<Q>>>> {
    items(g, n)?
        .into_iter()
        .map(|plane| items(g, plane)?.into_iter().map(|r| rat_vec(g, r)).collect::<Option<Vec<_>>>())
        .collect()
}

/// The algebra of a term and, for a matrix basis, the matrices.
fn read_lie(
    g: &Graph,
    n: NodeId,
) -> Option<(Lie, Option<Vec<QMat>>)> {
    let sc = g.ops().lookup("lie_sc")?;
    let wrapped = if g.op(n) == sc { Some(n) } else { g.enodes(g.find(n)).find(|&e| g.op(e) == sc) };
    if let Some(node) = wrapped {
        let [c] = g.children(node) else { return None };
        let c = read_constants(g, *c)?;
        if c.len() > MAX_BASIS {
            return None;
        }
        return Some((Lie::from_constants(c)?, None));
    }
    let basis: Vec<QMat> = items(g, n)?.into_iter().map(|m| rat_matrix(g, m)).collect::<Option<_>>()?;
    if basis.is_empty() || basis.len() > MAX_BASIS || basis.iter().any(|m| m.len() != m[0].len() || m.len() > MAX_MATRIX) {
        return None;
    }
    let lie = Lie::from_matrices(&basis)?;
    Some((lie, Some(basis)))
}

fn lie_arg(
    cx: &Cx<'_>,
    a: &[NodeId],
) -> Option<Lie> {
    read_lie(cx.graph, *a.first()?).map(|(l, _)| l)
}

// ----------------------------------------------------------------------
// Writing
// ----------------------------------------------------------------------

fn q_v(x: &Q) -> V {
    V::Rat(x.clone())
}

fn vec_v(v: &[Q]) -> V {
    V::List(v.iter().map(q_v).collect())
}

fn mat_v(m: &[Vec<Q>]) -> V {
    V::List(m.iter().map(|r| vec_v(r)).collect())
}

fn int_vec_v(v: &[i64]) -> V {
    V::ints(v.iter().copied())
}

fn pairs_v(pairs: impl IntoIterator<Item = (Vec<i64>, BigInt)>) -> V {
    V::List(pairs.into_iter().map(|(w, m)| V::List(vec![int_vec_v(&w), V::Int(m)])).collect())
}

// ----------------------------------------------------------------------
// Root systems
// ----------------------------------------------------------------------

fn lie_cartan_matrix(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let c = read_type(cx.graph, *a.first()?)?;
    ks::recognize(&c)?;
    Some(V::List(c.iter().map(|r| int_vec_v(r)).collect()))
}

fn lie_cartan_type(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let parts = ks::recognize(&read_type(cx.graph, *a.first()?)?)?;
    Some(V::Str(parts.iter().map(|(l, n)| format!("{l}{n}")).collect::<Vec<_>>().join("+")))
}

fn lie_rank(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    Some(V::uint(root_system(cx.graph, *a.first()?)?.n))
}

fn lie_dimension(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    Some(V::uint(root_system(cx.graph, *a.first()?)?.dimension()))
}

fn lie_num_positive_roots(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    Some(V::uint(root_system(cx.graph, *a.first()?)?.num_positive()))
}

fn lie_simple_roots(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let c = read_type(cx.graph, *a.first()?)?;
    let n = c.len();
    let (letter, _) = *ks::recognize(&c)?.first()?;
    (ks::recognize(&c)?.len() == 1 && ks::cartan_matrix(letter, n)? == c).then_some(())?;
    Some(mat_v(&ks::euclidean_simple_roots(letter, n)?))
}

fn lie_positive_roots(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let r = root_system(cx.graph, *a.first()?)?;
    Some(V::List(r.positive.iter().map(|x| int_vec_v(x)).collect()))
}

fn lie_highest_root(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    Some(int_vec_v(&root_system(cx.graph, *a.first()?)?.highest_root()?))
}

fn lie_weyl_group_order(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    Some(V::Int(root_system(cx.graph, *a.first()?)?.weyl_order()))
}

fn lie_exponents(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    Some(V::ints(root_system(cx.graph, *a.first()?)?.exponents().into_iter().map(|e| e as u64)))
}

fn lie_coxeter_number(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    Some(V::uint(root_system(cx.graph, *a.first()?)?.coxeter_number()?))
}

fn lie_dual_coxeter_number(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    Some(q_v(&root_system(cx.graph, *a.first()?)?.dual_coxeter_number()?))
}

fn lie_dynkin_diagram(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let r = root_system(cx.graph, *a.first()?)?;
    let mut edges = Vec::new();
    for i in 0..r.n {
        for j in 0..r.n {
            let (aij, aji) = (-r.cartan[i][j], -r.cartan[j][i]);
            if i == j || aij <= 0 || aji <= 0 {
                continue;
            }
            let multiplicity = aij.max(aji);
            // the longer root first; simple edges once
            if aij > aji || (aij == aji && i < j) {
                edges.push(V::ints([i as i64 + 1, j as i64 + 1, multiplicity]));
            }
        }
    }
    Some(V::List(edges))
}

fn lie_fundamental_weights(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    Some(mat_v(root_system(cx.graph, *a.first()?)?.fundamental_weights()))
}

fn lie_gram(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    Some(mat_v(&root_system(cx.graph, *a.first()?)?.gram()))
}

fn dim_irrep(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [t, hw] = a else { return None };
    let r = root_system(cx.graph, *t)?;
    Some(V::Int(r.weyl_dim(&labels(cx.graph, *hw, r.n)?)?))
}

fn weight_multiplicity(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [t, hw, w] = a else { return None };
    let r = root_system(cx.graph, *t)?;
    let (top, w) = (labels(cx.graph, *hw, r.n)?, labels(cx.graph, *w, r.n)?);
    Some(V::Int(r.weight_multiplicity(&top, &w)?))
}

fn irrep_weights(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [t, hw] = a else { return None };
    let r = root_system(cx.graph, *t)?;
    let top = labels(cx.graph, *hw, r.n)?;
    let table: BTreeMap<Vec<i64>, BigInt> = r.dominant_multiplicities(&top, 20_000)?;
    let mut ordered: Vec<(Vec<i64>, BigInt)> = table.into_iter().collect();
    ordered.reverse();
    Some(pairs_v(ordered))
}

fn tensor_decomposition(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [t, x, y] = a else { return None };
    let r = root_system(cx.graph, *t)?;
    let (x, y) = (labels(cx.graph, *x, r.n)?, labels(cx.graph, *y, r.n)?);
    Some(pairs_v(r.tensor(&x, &y)?))
}

fn casimir_eigenvalue(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [t, hw] = a else { return None };
    let r = root_system(cx.graph, *t)?;
    Some(q_v(&r.casimir(&labels(cx.graph, *hw, r.n)?)?))
}

fn dynkin_index(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [t, hw] = a else { return None };
    let r = root_system(cx.graph, *t)?;
    Some(q_v(&r.dynkin_index(&labels(cx.graph, *hw, r.n)?)?))
}

// ----------------------------------------------------------------------
// Structure theory
// ----------------------------------------------------------------------

fn lie_sc_of(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let lie = lie_arg(cx, a)?;
    let planes = V::List(lie.c.iter().map(|plane| V::List(plane.iter().map(|r| vec_v(r)).collect())).collect());
    let node = planes.build(cx.graph);
    apply(cx.graph, "lie_sc", &[node]).map(V::Node)
}

fn bases_v(series: &[QMat]) -> V {
    V::List(series.iter().map(|b| mat_v(b)).collect())
}

fn dims_v(series: &[QMat]) -> V {
    V::ints(series.iter().map(|b| b.len() as u64))
}

fn lie_derived_series(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    Some(bases_v(&lie_arg(cx, a)?.derived_series()))
}

fn lie_derived_dims(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    Some(dims_v(&lie_arg(cx, a)?.derived_series()))
}

fn lie_lower_central_series(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    Some(bases_v(&lie_arg(cx, a)?.lower_central_series()))
}

fn lie_lower_central_dims(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    Some(dims_v(&lie_arg(cx, a)?.lower_central_series()))
}

fn lie_upper_central_series(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    Some(bases_v(&lie_arg(cx, a)?.upper_central_series()))
}

fn lie_upper_central_dims(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    Some(dims_v(&lie_arg(cx, a)?.upper_central_series()))
}

fn lie_center(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    Some(mat_v(&lie_arg(cx, a)?.center()))
}

fn lie_radical(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    Some(mat_v(&lie_arg(cx, a)?.radical()))
}

fn lie_is_solvable(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    Some(V::Bool(lie_arg(cx, a)?.is_solvable()))
}

fn lie_is_nilpotent(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    Some(V::Bool(lie_arg(cx, a)?.is_nilpotent()))
}

fn lie_is_semisimple(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    Some(V::Bool(lie_arg(cx, a)?.is_semisimple()))
}

fn lie_cartan_solvable_test(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    Some(V::Bool(lie_arg(cx, a)?.cartan_solvable_test()))
}

fn lie_levi(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (levi, radical) = lie_arg(cx, a)?.levi()?;
    Some(V::List(vec![mat_v(&levi), mat_v(&radical)]))
}

fn lie_cartan_subalgebra(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    Some(mat_v(&lie_arg(cx, a)?.cartan_subalgebra()))
}

fn lie_root_decomposition(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [l, h] = a else { return None };
    let (lie, _) = read_lie(cx.graph, *l)?;
    let h: QMat = items(cx.graph, *h)?.into_iter().map(|v| rat_vec(cx.graph, v)).collect::<Option<_>>()?;
    if h.iter().any(|v| v.len() != lie.n) {
        return None;
    }
    let (spaces, missing) = lie.root_decomposition(&h)?;
    let spaces = V::List(spaces.iter().map(|s| V::List(vec![vec_v(&s.root), mat_v(&s.space)])).collect());
    Some(V::List(vec![spaces, V::uint(missing)]))
}

fn lie_is_subalgebra(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [l, v] = a else { return None };
    let (lie, _) = read_lie(cx.graph, *l)?;
    let vectors: QMat = items(cx.graph, *v)?.into_iter().map(|x| rat_vec(cx.graph, x)).collect::<Option<_>>()?;
    (vectors.iter().all(|x| x.len() == lie.n)).then(|| V::Bool(lie.is_subalgebra(&vectors)))
}

fn lie_elements(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [l, v] = a else { return None };
    let (lie, basis) = read_lie(cx.graph, *l)?;
    let basis = basis?;
    let size = basis.first()?.len();
    let mut out = Vec::new();
    for x in items(cx.graph, *v)? {
        let coords = rat_vec(cx.graph, x)?;
        if coords.len() != lie.n {
            return None;
        }
        let mut m = vec![ql::zeros(size); size];
        for (c, b) in coords.iter().zip(&basis) {
            for (mrow, brow) in m.iter_mut().zip(b) {
                for (x, y) in mrow.iter_mut().zip(brow) {
                    *x += c * y;
                }
            }
        }
        out.push(mat_v(&m));
    }
    Some(V::List(out))
}

fn lie_casimir_operator(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (lie, basis) = read_lie(cx.graph, *a.first()?)?;
    let rho: Vec<QMat> = match a {
        | [_] => basis?,
        | [_, r] => items(cx.graph, *r)?.into_iter().map(|m| rat_matrix(cx.graph, m)).collect::<Option<_>>()?,
        | _ => return None,
    };
    Some(mat_v(&lie.casimir(&rho)?))
}

// ----------------------------------------------------------------------
// Standard bases, subgroups, logarithm
// ----------------------------------------------------------------------

fn gauss(
    cx: &mut Cx<'_>,
    (re, im): (i64, i64),
) -> Option<NodeId> {
    let re_node = cx.graph.int(re);
    if im == 0 {
        return Some(re_node);
    }
    let unit = apply(cx.graph, "I", &[])?;
    let im_node = cx.graph.int(im);
    let imaginary = prod(cx.graph, &[im_node, unit]);
    let total = if re == 0 { imaginary } else { sum(cx.graph, &[re_node, imaginary]) };
    Some(cx.simplify(total))
}

fn standard(
    cx: &mut Cx<'_>,
    a: &[NodeId],
    family: Family,
) -> Option<V> {
    let n = idx(cx.graph, *a.first()?)?;
    let basis = ks::standard_basis(family, n)?;
    let mut out = Vec::with_capacity(basis.len());
    for m in &basis {
        let mut rows_v = Vec::with_capacity(m.len());
        for row in m {
            let mut cells = Vec::with_capacity(row.len());
            for &entry in row {
                cells.push(V::Node(gauss(cx, entry)?));
            }
            rows_v.push(V::List(cells));
        }
        out.push(V::List(rows_v));
    }
    Some(V::List(out))
}

fn gl_basis(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    standard(cx, a, Family::Gl)
}

fn sl_basis(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    standard(cx, a, Family::Sl)
}

fn so_basis(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    standard(cx, a, Family::So)
}

fn sp_basis(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    standard(cx, a, Family::Sp)
}

fn su_basis(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    standard(cx, a, Family::Su)
}

fn one_parameter_subgroup(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [m, t] = a else { return None };
    let matrix = rows(cx.graph, *m)?;
    let n = matrix.len();
    if n == 0 || n > MAX_MATRIX || matrix.iter().any(|r| r.len() != n) {
        return None;
    }
    let mut scaled = Vec::with_capacity(n);
    for row in matrix {
        let cells: Vec<V> = row
            .into_iter()
            .map(|x| {
                let p = prod(cx.graph, &[*t, x]);
                V::Node(cx.simplify(p))
            })
            .collect();
        scaled.push(V::List(cells));
    }
    let term = V::List(scaled).build(cx.graph);
    apply(cx.graph, "exp_map", &[term]).map(V::Node)
}

fn literal(
    g: &Graph,
    n: NodeId,
    v: i64,
) -> bool {
    g.number_of(n).is_some_and(|x| *x == Number::from(v))
}

fn call1(
    cx: &mut Cx<'_>,
    name: &str,
    arg: NodeId,
) -> Option<NodeId> {
    let t = apply(cx.graph, name, &[arg])?;
    Some(cx.simplify(t))
}

fn scaled(
    cx: &mut Cx<'_>,
    factor: NodeId,
    entry: NodeId,
) -> NodeId {
    let p = prod(cx.graph, &[factor, entry]);
    cx.simplify(p)
}

fn power(
    cx: &mut Cx<'_>,
    base: NodeId,
    exponent: Number,
) -> NodeId {
    let e = cx.graph.num(exponent);
    let p = cx.graph.node(crate::graph::op::core::POW, &[base, e]);
    cx.simplify(p)
}

/// `x I + y M` for a square matrix `M` of nodes.
fn affine(
    cx: &mut Cx<'_>,
    x: NodeId,
    y: NodeId,
    m: &[Vec<NodeId>],
) -> V {
    let n = m.len();
    let mut out = Vec::with_capacity(n);
    for (i, row) in m.iter().enumerate() {
        let mut cells = Vec::with_capacity(n);
        for (j, &entry) in row.iter().enumerate() {
            let mut terms = vec![scaled(cx, y, entry)];
            if i == j {
                terms.push(x);
            }
            let total = sum(cx.graph, &terms);
            cells.push(V::Node(cx.simplify(total)));
        }
        out.push(V::List(cells));
    }
    V::List(out)
}

fn matrix_log(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let g = rows(cx.graph, *a.first()?)?;
    let n = g.len();
    if n == 0 || n > MAX_MATRIX || g.iter().any(|r| r.len() != n) {
        return None;
    }
    let off_diagonal_zero = g.iter().enumerate().all(|(i, r)| r.iter().enumerate().all(|(j, &x)| i == j || literal(cx.graph, x, 0)));
    let zero = cx.graph.int(0);
    if off_diagonal_zero {
        let mut out = Vec::with_capacity(n);
        for (i, row) in g.iter().enumerate() {
            let mut cells = Vec::with_capacity(n);
            for (j, &x) in row.iter().enumerate() {
                if i != j || literal(cx.graph, x, 1) {
                    cells.push(V::Node(zero));
                } else if cx.graph.number_of(x).is_some_and(|v| v.to_f64() <= 0.0) {
                    return None;
                } else {
                    cells.push(V::Node(call1(cx, "ln", x)?));
                }
            }
            out.push(V::List(cells));
        }
        return Some(V::List(out));
    }
    // unipotent: N = g - I is nilpotent, log g = sum (-1)^(k+1) N^k / k
    let minus_one = cx.graph.int(-1);
    let nil: Vec<Vec<NodeId>> = g
        .iter()
        .enumerate()
        .map(|(i, row)| {
            row.iter()
                .enumerate()
                .map(|(j, &x)| {
                    let shifted = if i == j { sum(cx.graph, &[x, minus_one]) } else { x };
                    cx.simplify(shifted)
                })
                .collect()
        })
        .collect();
    let mut power_k = nil.clone();
    let mut total = vec![vec![Vec::new(); n]; n];
    let mut nilpotent = false;
    for k in 1..=n {
        if power_k.iter().flatten().all(|&x| literal(cx.graph, x, 0)) {
            nilpotent = true;
            break;
        }
        let coefficient = cx.graph.num(Number::fraction(if k % 2 == 1 { 1 } else { -1 }, i64::try_from(k).ok()?)?);
        for (i, row) in power_k.iter().enumerate() {
            for (j, &x) in row.iter().enumerate() {
                total[i][j].push(prod(cx.graph, &[coefficient, x]));
            }
        }
        power_k = super::discrete::matmul(cx, &power_k, &nil)?;
    }
    if nilpotent || power_k.iter().flatten().all(|&x| literal(cx.graph, x, 0)) {
        let cells: Vec<V> = total
            .into_iter()
            .map(|row| {
                V::List(
                    row.into_iter()
                        .map(|terms| {
                            let s = sum(cx.graph, &terms);
                            V::Node(cx.simplify(s))
                        })
                        .collect(),
                )
            })
            .collect();
        return Some(V::List(cells));
    }
    // rational matrices with a closed-form logarithm
    let m = rat_matrix(cx.graph, *a.first()?)?;
    if n == 2 {
        return log_two_by_two(cx, &g, &m);
    }
    if n == 3 {
        return log_rotation(cx, &g, &m);
    }
    None
}

fn log_two_by_two(
    cx: &mut Cx<'_>,
    g: &[Vec<NodeId>],
    m: &QMat,
) -> Option<V> {
    let trace = &m[0][0] + &m[1][1];
    let det = &m[0][0] * &m[1][1] - &m[0][1] * &m[1][0];
    if !det.is_positive() {
        return None;
    }
    let delta = cx.graph.num(Number::rat(det.clone()));
    let log_delta = call1(cx, "ln", delta)?;
    let half = cx.graph.num(Number::fraction(1, 2)?);
    let mean = scaled(cx, half, log_delta);
    let inv_sqrt = power(cx, delta, Number::fraction(-1, 2)?);
    let c_exact = &trace * &trace - Q::from_integer(BigInt::from(4)) * &det;
    // g' = g / sqrt(det); c = tr g' / 2
    let tr = cx.graph.num(Number::rat(trace));
    let c = {
        let t = scaled(cx, inv_sqrt, tr);
        scaled(cx, half, t)
    };
    let shifted: Vec<Vec<NodeId>> = g.iter().map(|r| r.iter().map(|&x| scaled(cx, inv_sqrt, x)).collect()).collect();
    let minus_c = {
        let m1 = cx.graph.int(-1);
        scaled(cx, m1, c)
    };
    let centred: Vec<Vec<NodeId>> = shifted
        .iter()
        .enumerate()
        .map(|(i, r)| {
            r.iter()
                .enumerate()
                .map(|(j, &x)| {
                    let s = if i == j { sum(cx.graph, &[x, minus_c]) } else { x };
                    cx.simplify(s)
                })
                .collect()
        })
        .collect();
    let one = cx.graph.int(1);
    if c_exact.is_zero() {
        // parabolic: X = mean I + (g' - I)
        let m1 = cx.graph.int(-1);
        let plain: Vec<Vec<NodeId>> = shifted
            .iter()
            .enumerate()
            .map(|(i, r)| {
                r.iter()
                    .enumerate()
                    .map(|(j, &x)| {
                        let s = if i == j { sum(cx.graph, &[x, m1]) } else { x };
                        cx.simplify(s)
                    })
                    .collect()
            })
            .collect();
        return Some(affine(cx, mean, one, &plain));
    }
    let c_squared = power(cx, c, Number::from(2));
    let m1 = cx.graph.int(-1);
    let one_minus = {
        let neg = scaled(cx, m1, c_squared);
        let s = sum(cx.graph, &[one, neg]);
        cx.simplify(s)
    };
    if c_exact.is_negative() {
        // elliptic: theta = acos(c), X = mean I + theta / sin(theta) (g' - c I)
        let theta = call1(cx, "acos", c)?;
        let sin = power(cx, one_minus, Number::fraction(1, 2)?);
        let ratio = {
            let inv = power(cx, sin, Number::from(-1));
            scaled(cx, theta, inv)
        };
        return Some(affine(cx, mean, ratio, &centred));
    }
    if (&m[0][0] + &m[1][1]).is_positive() {
        // hyperbolic with positive trace: s = acosh(c)
        let s = call1(cx, "acosh", c)?;
        let root = {
            let neg = scaled(cx, m1, one_minus);
            power(cx, neg, Number::fraction(-1, 2)?)
        };
        let ratio = scaled(cx, s, root);
        return Some(affine(cx, mean, ratio, &centred));
    }
    None
}

fn log_rotation(
    cx: &mut Cx<'_>,
    g: &[Vec<NodeId>],
    m: &QMat,
) -> Option<V> {
    let transpose = ql::transpose(m, 3);
    let product = ql::matmul(m, &transpose);
    if product != ql::identity(3) || !ql::det(m).is_one() {
        return None;
    }
    let trace = &m[0][0] + &m[1][1] + &m[2][2];
    let cos = (trace - Q::one()) / ql::q(2);
    if cos == -Q::one() {
        return None;
    }
    let c = cx.graph.num(Number::rat(cos));
    let theta = call1(cx, "acos", c)?;
    let c2 = power(cx, c, Number::from(2));
    let m1 = cx.graph.int(-1);
    let one = cx.graph.int(1);
    let sin = {
        let neg = scaled(cx, m1, c2);
        let s = sum(cx.graph, &[one, neg]);
        let s = cx.simplify(s);
        power(cx, s, Number::fraction(1, 2)?)
    };
    let factor = {
        let inv = power(cx, sin, Number::from(-1));
        let t = scaled(cx, theta, inv);
        let half = cx.graph.num(Number::fraction(1, 2)?);
        scaled(cx, half, t)
    };
    let mut out = Vec::new();
    for (i, row) in g.iter().enumerate() {
        let mut cells = Vec::new();
        for (j, &entry) in row.iter().enumerate() {
            let neg = scaled(cx, m1, g[j][i]);
            let diff = sum(cx.graph, &[entry, neg]);
            let diff = cx.simplify(diff);
            cells.push(V::Node(scaled(cx, factor, diff)));
        }
        out.push(V::List(cells));
    }
    Some(V::List(out))
}

fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    def(i, "lie_cartan_matrix", Arity::Fixed(1), lie_cartan_matrix)?;
    def(i, "lie_cartan_type", Arity::Fixed(1), lie_cartan_type)?;
    def(i, "lie_rank", Arity::Fixed(1), lie_rank)?;
    def(i, "lie_dimension", Arity::Fixed(1), lie_dimension)?;
    def(i, "lie_num_positive_roots", Arity::Fixed(1), lie_num_positive_roots)?;
    def(i, "lie_simple_roots", Arity::Fixed(1), lie_simple_roots)?;
    def(i, "lie_positive_roots", Arity::Fixed(1), lie_positive_roots)?;
    def(i, "lie_highest_root", Arity::Fixed(1), lie_highest_root)?;
    def(i, "lie_weyl_group_order", Arity::Fixed(1), lie_weyl_group_order)?;
    def(i, "lie_exponents", Arity::Fixed(1), lie_exponents)?;
    def(i, "lie_coxeter_number", Arity::Fixed(1), lie_coxeter_number)?;
    def(i, "lie_dual_coxeter_number", Arity::Fixed(1), lie_dual_coxeter_number)?;
    def(i, "lie_dynkin_diagram", Arity::Fixed(1), lie_dynkin_diagram)?;
    def(i, "lie_fundamental_weights", Arity::Fixed(1), lie_fundamental_weights)?;
    def(i, "lie_gram", Arity::Fixed(1), lie_gram)?;
    def(i, "dim_irrep", Arity::Fixed(2), dim_irrep)?;
    def(i, "weight_multiplicity", Arity::Fixed(3), weight_multiplicity)?;
    def(i, "irrep_weights", Arity::Fixed(2), irrep_weights)?;
    def(i, "tensor_decomposition", Arity::Fixed(3), tensor_decomposition)?;
    def(i, "casimir_eigenvalue", Arity::Fixed(2), casimir_eigenvalue)?;
    def(i, "dynkin_index", Arity::Fixed(2), dynkin_index)?;
    def_inert(i, "lie_sc", Arity::Fixed(1))?;
    def(i, "lie_sc_of", Arity::Fixed(1), lie_sc_of)?;
    def(i, "lie_derived_series", Arity::Fixed(1), lie_derived_series)?;
    def(i, "lie_derived_dims", Arity::Fixed(1), lie_derived_dims)?;
    def(i, "lie_lower_central_series", Arity::Fixed(1), lie_lower_central_series)?;
    def(i, "lie_lower_central_dims", Arity::Fixed(1), lie_lower_central_dims)?;
    def(i, "lie_upper_central_series", Arity::Fixed(1), lie_upper_central_series)?;
    def(i, "lie_upper_central_dims", Arity::Fixed(1), lie_upper_central_dims)?;
    def(i, "lie_center", Arity::Fixed(1), lie_center)?;
    def(i, "lie_radical", Arity::Fixed(1), lie_radical)?;
    def(i, "lie_is_solvable", Arity::Fixed(1), lie_is_solvable)?;
    def(i, "lie_is_nilpotent", Arity::Fixed(1), lie_is_nilpotent)?;
    def(i, "lie_is_semisimple", Arity::Fixed(1), lie_is_semisimple)?;
    def(i, "lie_cartan_solvable_test", Arity::Fixed(1), lie_cartan_solvable_test)?;
    def(i, "lie_levi", Arity::Fixed(1), lie_levi)?;
    def(i, "lie_cartan_subalgebra", Arity::Fixed(1), lie_cartan_subalgebra)?;
    def(i, "lie_root_decomposition", Arity::Fixed(2), lie_root_decomposition)?;
    def(i, "lie_is_subalgebra", Arity::Fixed(2), lie_is_subalgebra)?;
    def(i, "lie_elements", Arity::Fixed(2), lie_elements)?;
    def(i, "lie_casimir_operator", Arity::Variadic, lie_casimir_operator)?;
    def(i, "gl_basis", Arity::Fixed(1), gl_basis)?;
    def(i, "sl_basis", Arity::Fixed(1), sl_basis)?;
    def(i, "so_basis", Arity::Fixed(1), so_basis)?;
    def(i, "sp_basis", Arity::Fixed(1), sp_basis)?;
    def(i, "su_basis", Arity::Fixed(1), su_basis)?;
    def_request(i, "one_parameter_subgroup", Arity::Fixed(2), one_parameter_subgroup)?;
    def_request(i, "matrix_log", Arity::Fixed(1), matrix_log)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rules::testing::simplify;

    fn run(src: &str) -> String {
        simplify(&[lie_structure()], src)
    }

    #[test]
    fn weyl_dimensions_of_su3() {
        assert_eq!(run("dim_irrep(A2, list(1, 0))"), "3");
        assert_eq!(run("dim_irrep(A2, list(0, 1))"), "3");
        assert_eq!(run("dim_irrep(A2, list(2, 0))"), "6");
        assert_eq!(run("dim_irrep(A2, list(1, 1))"), "8");
        assert_eq!(run("dim_irrep(A2, list(3, 0))"), "10");
        assert_eq!(run("dim_irrep(A2, list(2, 2))"), "27");
        assert_eq!(run("dim_irrep(list(list(2, -1), list(-1, 2)), list(1, 1))"), "8");
        assert_eq!(run("dim_irrep(E8, list(0, 0, 0, 0, 0, 0, 0, 1))"), "248");
        assert_eq!(run("dim_irrep(G2, list(0, 1))"), "14");
        assert_eq!(run("dim_irrep(A2, list(1, 1, 1))"), "dim_irrep(A2, list(1, 1, 1))");
        assert_eq!(run("dim_irrep(A2, list(-1, 1))"), "dim_irrep(A2, list(-1, 1))");
    }

    #[test]
    fn classification_data() {
        assert_eq!(run("lie_weyl_group_order(G2)"), "12");
        assert_eq!(run("lie_weyl_group_order(E8)"), "696729600");
        assert_eq!(run("lie_weyl_group_order(A4)"), "120");
        assert_eq!(run("lie_dimension(E8)"), "248");
        assert_eq!(run("lie_dimension(E7)"), "133");
        assert_eq!(run("lie_dimension(F4)"), "52");
        assert_eq!(run("lie_coxeter_number(E6)"), "12");
        assert_eq!(run("lie_dual_coxeter_number(B3)"), "5");
        assert_eq!(run("lie_highest_root(G2)"), "list(3, 2)");
        assert_eq!(run("lie_highest_root(A3)"), "list(1, 1, 1)");
        assert_eq!(run("lie_positive_roots(A2)"), "list(list(1, 0), list(0, 1), list(1, 1))");
        assert_eq!(run("lie_exponents(G2)"), "list(1, 5)");
        assert_eq!(run("lie_cartan_matrix(G2)"), "list(list(2, -1), list(-3, 2))");
        assert_eq!(run("lie_cartan_matrix(B3)"), "list(list(2, -1, 0), list(-1, 2, -2), list(0, -1, 2))");
        assert_eq!(run("lie_dynkin_diagram(G2)"), "list(list(2, 1, 3))");
        assert_eq!(run("lie_dynkin_diagram(B3)"), "list(list(1, 2, 1), list(2, 3, 2))");
        assert_eq!(run("lie_fundamental_weights(A2)"), "list(list(2/3, 1/3), list(1/3, 2/3))");
        assert_eq!(run("lie_gram(B2)"), "list(list(2, -1), list(-1, 1))");
        assert_eq!(run("lie_simple_roots(B2)"), "list(list(1, -1), list(0, 1))");
        assert_eq!(run("lie_simple_roots(G2)"), "list(list(1, -1, 0), list(-2, 1, 1))");
    }

    #[test]
    fn cartan_matrix_types() {
        assert_eq!(run("lie_cartan_type(list(list(2, -1), list(-1, 2)))"), "\"A2\"");
        assert_eq!(run("lie_cartan_type(list(list(2, -1), list(-2, 2)))"), "\"B2\"");
        assert_eq!(run("lie_cartan_type(list(list(2, -3), list(-1, 2)))"), "\"G2\"");
        assert_eq!(run("lie_cartan_type(lie_cartan_matrix(D5))"), "\"D5\"");
        assert_eq!(run("lie_cartan_type(lie_cartan_matrix(E7))"), "\"E7\"");
        assert_eq!(run("lie_cartan_type(lie_cartan_matrix(F4))"), "\"F4\"");
        assert_eq!(run("lie_cartan_type(lie_cartan_matrix(C4))"), "\"C4\"");
        assert_eq!(run("lie_cartan_type(list(list(2, 0, -1), list(0, 2, 0), list(-1, 0, 2)))"), "\"A2+A1\"");
        // not the Cartan matrix of a semisimple algebra
        let bad = "lie_cartan_type(list(list(2, -2), list(-2, 2)))";
        assert_eq!(crate::rules::testing::reduce_with(&[lie_structure()], bad, &[]).0, bad);
    }

    #[test]
    fn multiplicities_and_tensor_products() {
        assert_eq!(run("weight_multiplicity(A2, list(1, 1), list(0, 0))"), "2");
        assert_eq!(run("weight_multiplicity(A2, list(1, 1), list(2, -1))"), "1");
        assert_eq!(run("weight_multiplicity(A2, list(2, 2), list(0, 0))"), "3");
        assert_eq!(run("irrep_weights(A2, list(1, 1))"), "list(list(list(1, 1), 1), list(list(0, 0), 2))");
        // 3 x 3bar = 8 + 1
        assert_eq!(run("tensor_decomposition(A2, list(1, 0), list(0, 1))"), "list(list(list(0, 0), 1), list(list(1, 1), 1))");
        assert_eq!(run("tensor_decomposition(A2, list(1, 0), list(1, 0))"), "list(list(list(0, 1), 1), list(list(2, 0), 1))");
        assert_eq!(run("tensor_decomposition(A1, list(2), list(3))"), "list(list(list(1), 1), list(list(3), 1), list(list(5), 1))");
        assert_eq!(run("casimir_eigenvalue(A2, list(1, 0))"), "8/3");
        assert_eq!(run("casimir_eigenvalue(A2, list(1, 1))"), "6");
        assert_eq!(run("dynkin_index(A2, list(1, 0))"), "1");
        assert_eq!(run("dynkin_index(A2, list(1, 1))"), "6");
    }

    const UPPER: &str = "list(list(list(1, 0, 0), list(0, 0, 0), list(0, 0, 0)), list(list(0, 1, 0), list(0, 0, 0), list(0, 0, 0)), list(list(0, 0, 1), list(0, 0, 0), list(0, 0, 0)), list(list(0, 0, 0), list(0, 1, 0), list(0, 0, 0)), list(list(0, 0, 0), list(0, 0, 1), list(0, 0, 0)), list(list(0, 0, 0), list(0, 0, 0), list(0, 0, 1)))";
    const HEISENBERG: &str = "list(list(list(0, 1, 0), list(0, 0, 0), list(0, 0, 0)), list(list(0, 0, 0), list(0, 0, 1), list(0, 0, 0)), list(list(0, 0, 1), list(0, 0, 0), list(0, 0, 0)))";

    #[test]
    fn solvable_nilpotent_semisimple() {
        assert_eq!(run(&format!("lie_is_solvable({UPPER})")), "true");
        assert_eq!(run(&format!("lie_is_nilpotent({UPPER})")), "false");
        assert_eq!(run(&format!("lie_is_semisimple({UPPER})")), "false");
        assert_eq!(run(&format!("lie_cartan_solvable_test({UPPER})")), "true");
        assert_eq!(run(&format!("lie_derived_dims({UPPER})")), "list(6, 3, 1, 0)");
        assert_eq!(run(&format!("lie_is_nilpotent({HEISENBERG})")), "true");
        assert_eq!(run(&format!("lie_lower_central_dims({HEISENBERG})")), "list(3, 1, 0)");
        assert_eq!(run(&format!("lie_upper_central_dims({HEISENBERG})")), "list(1, 3)");
        assert_eq!(run(&format!("lie_center({HEISENBERG})")), "list(list(0, 0, 1))");
        for basis in ["sl_basis(3)", "so_basis(4)", "sp_basis(2)"] {
            assert_eq!(run(&format!("lie_is_semisimple({basis})")), "true", "{basis}");
            assert_eq!(run(&format!("lie_is_solvable({basis})")), "false", "{basis}");
            assert_eq!(run(&format!("lie_cartan_solvable_test({basis})")), "false", "{basis}");
            assert_eq!(run(&format!("lie_radical({basis})")), "list()", "{basis}");
        }
        assert_eq!(run("lie_is_semisimple(gl_basis(2))"), "false");
        assert_eq!(run("lie_radical(gl_basis(2))"), "list(list(1, 0, 0, 1))");
        // a list that is not closed is not a Lie algebra
        let open = "lie_is_solvable(list(list(list(0, 1), list(0, 0)), list(list(0, 0), list(1, 0))))";
        assert_eq!(crate::rules::testing::reduce_with(&[lie_structure()], open, &[]).0, open);
        // structure constants of so(3)
        assert_eq!(run("lie_is_semisimple(lie_sc(structure_constants(so3_basis())))"), "true");
        assert_eq!(run("lie_is_semisimple(lie_sc(structure_constants(su2_basis())))"), "true");
        assert_eq!(run("lie_derived_dims(lie_sc_of(so_basis(3)))"), "list(3)");
    }

    #[test]
    fn levi_cartan_roots() {
        // sl2 acting on C^2 inside 3x3 matrices
        let basis = "list(list(list(1, 0, 0), list(0, -1, 0), list(0, 0, 0)), list(list(0, 1, 0), list(0, 0, 0), list(0, 0, 0)), list(list(0, 0, 0), list(1, 0, 0), list(0, 0, 0)), list(list(0, 0, 1), list(0, 0, 0), list(0, 0, 0)), list(list(0, 0, 0), list(0, 0, 1), list(0, 0, 0)))";
        let levi = run(&format!("lie_levi({basis})"));
        assert!(levi.starts_with("list(list(list("), "{levi}");
        assert_eq!(run(&format!("lie_radical({basis})")), "list(list(0, 0, 0, 1, 0), list(0, 0, 0, 0, 1))");
        let h = run("lie_cartan_subalgebra(sl_basis(3))");
        assert_eq!(h.matches("list(").count(), 3, "{h}");
        let roots = run("lie_root_decomposition(sl_basis(3), lie_cartan_subalgebra(sl_basis(3)))");
        assert!(roots.ends_with(", 0)"), "{roots}");
        assert!(roots.matches("list(list(").count() >= 7);
        let so3 = run("lie_root_decomposition(sl_basis(2), lie_cartan_subalgebra(sl_basis(2)))");
        assert!(so3.ends_with(", 0)"), "{so3}");
        let elements = run("lie_elements(sl_basis(2), list(list(1, 0, 0), list(0, 0, 1)))");
        assert_eq!(elements, "list(list(list(0, 1), list(0, 0)), list(list(1, 0), list(0, -1)))");
    }

    #[test]
    fn standard_bases_and_casimir() {
        assert_eq!(run("sl_basis(2)"), "list(list(list(0, 1), list(0, 0)), list(list(0, 0), list(1, 0)), list(list(1, 0), list(0, -1)))");
        assert_eq!(run("so_basis(3)"), "list(list(list(0, 1, 0), list(-1, 0, 0), list(0, 0, 0)), list(list(0, 0, 1), list(0, 0, 0), list(-1, 0, 0)), list(list(0, 0, 0), list(0, 0, 1), list(0, -1, 0)))");
        assert_eq!(run("su_basis(2)"), "list(list(list(0, 1), list(-1, 0)), list(list(0, I), list(I, 0)), list(list(I, 0), list(0, -I)))");
        assert_eq!(run("lie_dimension(A2)"), "8");
        for (basis, dim) in [("sl_basis(4)", 15), ("so_basis(5)", 10), ("sp_basis(3)", 21), ("gl_basis(3)", 9), ("su_basis(3)", 8)] {
            let text = run(&format!("is_lie_algebra({basis})"));
            assert_eq!(text, "true", "{basis}");
            let n = run(&format!("lie_derived_dims(lie_sc(structure_constants({basis})))"));
            assert!(n.starts_with(&format!("list({dim}")), "{basis}: {n}");
        }
        assert_eq!(run("lie_casimir_operator(sl_basis(2))"), "list(list(3/8, 0), list(0, 3/8))");
        assert_eq!(run("lie_casimir_operator(so_basis(3))"), "list(list(1, 0, 0), list(0, 1, 0), list(0, 0, 1))");
    }

    #[test]
    fn subgroups_and_logarithms() {
        assert_eq!(run("one_parameter_subgroup(list(list(0, -1), list(1, 0)), t)"), "list(list(cos(t), -sin(t)), list(sin(t), cos(t)))");
        assert_eq!(run("one_parameter_subgroup(list(list(1, 0), list(0, -1)), t)"), "list(list(exp(t), 0), list(0, exp(-t)))");
        assert_eq!(run("matrix_log(list(list(1, 0), list(0, 1)))"), "list(list(0, 0), list(0, 0))");
        assert_eq!(run("matrix_log(list(list(1, 1), list(0, 1)))"), "list(list(0, 1), list(0, 0))");
        assert_eq!(run("matrix_log(list(list(1, 2, 3), list(0, 1, 4), list(0, 0, 1)))"), "list(list(0, 2, -1), list(0, 0, 4), list(0, 0, 0))");
        assert_eq!(run("matrix_log(list(list(2, 0), list(0, 1)))"), "list(list(ln(2), 0), list(0, 0))");
        assert_eq!(run("exp_map(matrix_log(list(list(1, 2, 3), list(0, 1, 4), list(0, 0, 1))))"), "list(list(1, 2, 3), list(0, 1, 4), list(0, 0, 1))");
        assert_eq!(run("matrix_log(list(list(0, -1), list(1, 0)))"), "list(list(0, -1/2*pi), list(1/2*pi, 0))");
        assert_eq!(run("matrix_log(list(list(0, -1, 0), list(1, 0, 0), list(0, 0, 1)))"), "list(list(0, -1/2*pi, 0), list(1/2*pi, 0, 0), list(0, 0, 0))");
        // not covered: unevaluated
        let open = "matrix_log(list(list(-1, 0), list(0, 1)))";
        assert_eq!(crate::rules::testing::reduce_with(&[lie_structure()], open, &[]).0, open);
    }
}
