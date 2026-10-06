//! Tensors, Riemannian geometry, curvilinear coordinates and differential
//! forms.
//!
//! A tensor is a nested list: a vector `list(a, b)`, a rank-2 tensor
//! `list(list(..), ..)`, and so on; which indices are upper or lower is the
//! caller's convention, made explicit by the metric-using operators. A
//! metric is the rank-2 tensor `g_ij` in the coordinates `vars`.
//!
//! | operator | value |
//! |---|---|
//! | `tensor_rank(T)`, `component(T, list(i, ...))` | rank; one component (indices from 1) |
//! | `tensor_add(A, B)`, `tensor_smul(c, T)`, `tensor_outer(A, B)` | componentwise sum, scalar multiple, outer product |
//! | `contract(T, i, j)` | trace over indices `i` and `j` (from 1) |
//! | `lower_index(V, g)`, `raise_index(W, g)` | `g_ij V^j`, `g^ij W_j` |
//! | `christoffel1(g, vars)`, `christoffel(g, vars)` | `Γ_kij`, `Γ^k_ij` as `[k][i][j]` |
//! | `riemann(g, vars)` | `R^r_smn` as `[r][s][m][n]` |
//! | `ricci(g, vars)`, `ricci_scalar(g, vars)`, `einstein(g, vars)` | `R_sn = R^r_srn`, `R = g^sn R_sn`, `G = Ric - R g / 2` |
//! | `covariant_derivative(V, g, vars)` | `∇_k V^i = ∂_k V^i + Γ^i_kj V^j` as `[i][k]` |
//! | `to_cartesian(sys, vars)`, `from_cartesian(sys, list(x, y[, z]))` | coordinate maps; `sys` is `cartesian`, `polar`, `cylindrical` or `spherical` (`r, theta, phi` with `theta` from the z axis) |
//! | `coordinate_metric(sys, vars)`, `scale_factors(sys, vars)` | `g = Jᵀ J`; `h_i = sqrt(g_ii)` |
//! | `transform_point(from, to, point)` | a point's coordinates in another system |
//! | `transform_expression(f, from, from_vars, to, to_vars)` | `f` rewritten in the other system's variables |
//! | `transform_vector(V, ...)`, `transform_covector(W, ...)`, `transform_tensor2(T, ...)` | contravariant, covariant and contravariant rank-2 components, same trailing arguments as `transform_expression` |
//! | `grad_in(f, sys, vars)`, `div_in(F, ...)`, `curl_in(F, ...)`, `laplacian_in(f, ...)` | vector calculus in orthogonal coordinates (physical components) |
//! | `form(list(list(c, list(i, ...)), ...))` | a differential form `Σ c dx_i ∧ ...` (indices from 1) |
//! | `exterior_d(ω, vars)`, `wedge(ω, η)` | exterior derivative and product |
//! | `greens_theorem(P, Q, list(x, y), list(list(x0, x1), list(y0, y1)))` | the circulation of `(P, Q)` around a rectangle, as `∬ (Q_x - P_y)` |
//! | `gauss_theorem(F, vars, bounds)` | the flux out of a box, as `∭ div F` |
//! | `stokes_theorem(F, surface, u, v, u0, u1, v0, v1)` | the circulation around the boundary of a parametrised surface, as the flux of `curl F` |
//! | `cell(φ, params, bounds)` | an oriented `k`-cell: the image of the box `bounds = list(list(a1, b1), ...)` of the parameters `params` under the chart `φ = list(x1(u), ..., xn(u))` (inert) |
//! | `chain(list(list(s, cell), ...))` | a formal sum of cells with integer multiplicities (inert) |
//! | `boundary(M)` | `∂M` of a cell (constant bounds) or chain: `Σ_i Σ_α (-1)^(i+α) M|_{u_i = a_i / b_i}` as a chain of `(k-1)`-cells; `boundary(boundary(M))` integrates to zero |
//! | `pullback(ω, vars, φ, params)` | `φ*ω`: the form in the parameters (`dx_i = Σ_j ∂φ_i/∂u_j du_j`) |
//! | `integrate_form(ω, vars, M)` | `∫_M ω` over a cell or chain: the pull-back's top-degree coefficient integrated over the box (iterated `defint`, first parameter innermost); a 0-form over a 0-cell is its value there |
//! | `geodesic_acceleration(g, vars, v)`, `kretschmann(g, vars)`, `gaussian_curvature(g, vars)`, `volume_element(g, vars)`, `laplace_beltrami(f, g, vars)`, `covariant_divergence(V, g, vars)`, `lie_bracket(V, W, vars)`, `killing_tensor(xi, g, vars)`, `is_killing(xi, g, vars)` | geodesics, curvature invariants, Beltrami operators, Lie brackets and Killing fields (see the `ext` submodule's table) |
//! | `generalized_stokes(ω, vars, M)` | `integrate_form(exterior_d(ω), vars, M) = integrate_form(ω, vars, boundary(M))`, both sides evaluated |

mod ext;

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
use crate::graph::RuleError;
use crate::graph::RuleSet;
use crate::graph::Tier;
use crate::graph::op::core;
use crate::graph::rule::Installer;
use crate::rules::calculus::derivative;
use crate::rules::complex::build::add;
use crate::rules::complex::build::mul;
use crate::rules::complex::build::neg;
use crate::rules::complex::build::pow;
use crate::rules::complex::build::powi;
use crate::rules::complex::build::sub;
use crate::rules::linalg::invert;
use crate::rules::linalg::linalg;
use crate::rules::linalg::normal_form;
use crate::rules::poly::best;

/// The geometry rule set.
#[must_use]
pub fn geometry() -> RuleSet {
    RuleSet::new("geometry", install).needs(linalg())
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum Request {
    Rank,
    Component,
    Add,
    Smul,
    Outer,
    Contract,
    Lower,
    Raise,
    Christoffel1,
    Christoffel2,
    Riemann,
    Ricci,
    RicciScalar,
    Einstein,
    Covariant,
    ToCartesian,
    FromCartesian,
    Metric,
    ScaleFactors,
    TransformPoint,
    TransformExpression,
    TransformVector,
    TransformCovector,
    TransformTensor2,
    GradIn,
    DivIn,
    CurlIn,
    LaplacianIn,
    ExteriorD,
    Wedge,
    Greens,
    Gauss,
    Stokes,
    Boundary,
    Pullback,
    IntegrateForm,
    GeneralizedStokes,
}

fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let table: [(&str, u8, Request); 37] = [
        ("tensor_rank", 1, Request::Rank),
        ("component", 2, Request::Component),
        ("tensor_add", 2, Request::Add),
        ("tensor_smul", 2, Request::Smul),
        ("tensor_outer", 2, Request::Outer),
        ("contract", 3, Request::Contract),
        ("lower_index", 2, Request::Lower),
        ("raise_index", 2, Request::Raise),
        ("christoffel1", 2, Request::Christoffel1),
        ("christoffel", 2, Request::Christoffel2),
        ("riemann", 2, Request::Riemann),
        ("ricci", 2, Request::Ricci),
        ("ricci_scalar", 2, Request::RicciScalar),
        ("einstein", 2, Request::Einstein),
        ("covariant_derivative", 3, Request::Covariant),
        ("to_cartesian", 2, Request::ToCartesian),
        ("from_cartesian", 2, Request::FromCartesian),
        ("coordinate_metric", 2, Request::Metric),
        ("scale_factors", 2, Request::ScaleFactors),
        ("transform_point", 3, Request::TransformPoint),
        ("transform_expression", 5, Request::TransformExpression),
        ("transform_vector", 5, Request::TransformVector),
        ("transform_covector", 5, Request::TransformCovector),
        ("transform_tensor2", 5, Request::TransformTensor2),
        ("grad_in", 3, Request::GradIn),
        ("div_in", 3, Request::DivIn),
        ("curl_in", 3, Request::CurlIn),
        ("laplacian_in", 3, Request::LaplacianIn),
        ("exterior_d", 2, Request::ExteriorD),
        ("wedge", 2, Request::Wedge),
        ("greens_theorem", 4, Request::Greens),
        ("gauss_theorem", 3, Request::Gauss),
        ("stokes_theorem", 8, Request::Stokes),
        ("boundary", 1, Request::Boundary),
        ("pullback", 4, Request::Pullback),
        ("integrate_form", 3, Request::IntegrateForm),
        ("generalized_stokes", 3, Request::GeneralizedStokes),
    ];
    i.op(OpDescriptor::new("form", Arity::Fixed(1)))?;
    i.op(OpDescriptor::new("cell", Arity::Fixed(3)))?;
    i.op(OpDescriptor::new("chain", Arity::Fixed(1)))?;
    for (name, arity, request) in table {
        let op = i.op(OpDescriptor::new(name, Arity::Fixed(arity)).flags(OpFlags::HEAVY).cost(100))?;
        i.kernel(&format!("geometry/{name}"), Tier::Reduce, Geometry { op, request });
    }
    ext::install(i)?;
    Ok(())
}

// ----------------------------------------------------------------------
// Nested lists
// ----------------------------------------------------------------------

/// A dense tensor of terms in row-major order.
#[derive(Clone, Debug)]
struct Tensor {
    shape: Vec<usize>,
    data: Vec<NodeId>,
}

impl Tensor {
    /// Reads a nested list; a non-list is a scalar (rank 0).
    fn read(
        graph: &mut Graph,
        node: NodeId,
    ) -> Option<Self> {
        let term = best(graph, node)?;
        if graph.op(term) != core::LIST {
            return Some(Self { shape: Vec::new(), data: vec![term] });
        }
        let items = graph.children(term).to_vec();
        let mut parts = Vec::with_capacity(items.len());
        for item in items {
            parts.push(Self::read(graph, item)?);
        }
        let inner = parts.first().map_or_else(Vec::new, |p| p.shape.clone());
        if parts.iter().any(|p| p.shape != inner) {
            return None;
        }
        let mut shape = vec![parts.len()];
        shape.extend(inner);
        Some(Self { shape, data: parts.into_iter().flat_map(|p| p.data).collect() })
    }

    fn write(
        &self,
        graph: &mut Graph,
    ) -> NodeId {
        fn build(
            graph: &mut Graph,
            shape: &[usize],
            data: &[NodeId],
        ) -> NodeId {
            match shape.split_first() {
                | None => data.first().copied().unwrap_or(NodeId::NONE),
                | Some((&n, rest)) => {
                    let stride: usize = rest.iter().product();
                    let items: Vec<NodeId> = (0..n).map(|k| build(graph, rest, &data[k * stride..(k + 1) * stride])).collect();
                    graph.node(core::LIST, &items)
                },
            }
        }
        build(graph, &self.shape, &self.data)
    }

    fn zeros(
        graph: &mut Graph,
        shape: Vec<usize>,
    ) -> Self {
        let zero = graph.int(0);
        let len = shape.iter().product();
        Self { shape, data: vec![zero; len] }
    }

    fn offset(
        &self,
        index: &[usize],
    ) -> usize {
        index.iter().zip(&self.shape).fold(0, |acc, (&i, &n)| acc * n + i)
    }

    fn get(
        &self,
        index: &[usize],
    ) -> NodeId {
        self.data[self.offset(index)]
    }

    fn set(
        &mut self,
        index: &[usize],
        value: NodeId,
    ) {
        let k = self.offset(index);
        self.data[k] = value;
    }
}

/// All multi-indices of a shape.
fn indices(shape: &[usize]) -> Vec<Vec<usize>> {
    let mut out = vec![Vec::new()];
    for &n in shape {
        out = out.into_iter().flat_map(|prefix| (0..n).map(move |k| [prefix.clone(), vec![k]].concat())).collect();
    }
    out
}

fn vector(
    graph: &mut Graph,
    node: NodeId,
) -> Option<Vec<NodeId>> {
    let t = Tensor::read(graph, node)?;
    (t.shape.len() == 1).then_some(t.data)
}

fn square(
    graph: &mut Graph,
    node: NodeId,
) -> Option<Vec<Vec<NodeId>>> {
    let t = Tensor::read(graph, node)?;
    let [n, m] = t.shape.as_slice() else {
        return None;
    };
    (n == m).then(|| t.data.chunks(*n).map(<[NodeId]>::to_vec).collect())
}

fn list(
    graph: &mut Graph,
    items: &[NodeId],
) -> NodeId {
    graph.node(core::LIST, items)
}

fn matrix_term(
    graph: &mut Graph,
    rows: &[Vec<NodeId>],
) -> NodeId {
    let rows: Vec<NodeId> = rows.iter().map(|r| list(graph, r)).collect();
    list(graph, &rows)
}

fn small_index(
    graph: &Graph,
    node: NodeId,
) -> Option<usize> {
    usize::try_from(graph.number_of(node)?.to_i64()?).ok()?.checked_sub(1)
}

/// `∂f/∂x`, simplified.
fn d(
    cx: &mut Cx<'_>,
    f: NodeId,
    x: NodeId,
) -> Option<NodeId> {
    let term = best(cx.graph, f)?;
    let derivative = derivative(cx.graph, term, x)?;
    Some(cx.simplify(derivative))
}

fn half(graph: &mut Graph) -> NodeId {
    graph.num(Number::fraction(1, 2).unwrap_or_else(|| Number::from(0)))
}

fn sqrt(
    graph: &mut Graph,
    x: NodeId,
) -> NodeId {
    let h = half(graph);
    pow(graph, x, h)
}

// ----------------------------------------------------------------------
// Curvature
// ----------------------------------------------------------------------

/// `Γ_kij = (∂_i g_jk + ∂_j g_ik - ∂_k g_ij) / 2`.
fn christoffel_first(
    cx: &mut Cx<'_>,
    g: &[Vec<NodeId>],
    vars: &[NodeId],
) -> Option<Tensor> {
    let n = vars.len();
    if g.len() != n {
        return None;
    }
    // dg[a][b][c] = ∂_c g_ab
    let mut dg = vec![vec![vec![NodeId::NONE; n]; n]; n];
    for a in 0..n {
        for b in 0..n {
            for c in 0..n {
                dg[a][b][c] = d(cx, g[a][b], vars[c])?;
            }
        }
    }
    let mut out = Tensor::zeros(cx.graph, vec![n, n, n]);
    let h = half(cx.graph);
    #[allow(clippy::needless_range_loop)] // index is used for more than one array / arithmetic; iterator form would not be clearer
    for k in 0..n {
        #[allow(clippy::needless_range_loop)] // index is used for more than one array / arithmetic; iterator form would not be clearer
        for i in 0..n {
            #[allow(clippy::needless_range_loop)] // index is used for more than one array / arithmetic; iterator form would not be clearer
            for j in 0..n {
                let first = add(cx.graph, &[dg[j][k][i], dg[i][k][j]]);
                let difference = sub(cx.graph, first, dg[i][j][k]);
                let value = mul(cx.graph, &[h, difference]);
                out.set(&[k, i, j], cx.simplify(value));
            }
        }
    }
    Some(out)
}

/// `Γ^k_ij = g^kl Γ_lij`.
fn christoffel_second(
    cx: &mut Cx<'_>,
    g: &[Vec<NodeId>],
    vars: &[NodeId],
) -> Option<Tensor> {
    let n = vars.len();
    let first = christoffel_first(cx, g, vars)?;
    let inverse = invert(cx, g)?;
    let mut out = Tensor::zeros(cx.graph, vec![n, n, n]);
    #[allow(clippy::needless_range_loop)] // index is used for more than one array / arithmetic; iterator form would not be clearer
    for k in 0..n {
        for i in 0..n {
            for j in 0..n {
                let terms: Vec<NodeId> = (0..n).map(|l| mul(cx.graph, &[inverse[k][l], first.get(&[l, i, j])])).collect();
                let sum = add(cx.graph, &terms);
                out.set(&[k, i, j], normal_form(cx, sum));
            }
        }
    }
    Some(out)
}

/// `R^r_smn = ∂_m Γ^r_ns - ∂_n Γ^r_ms + Γ^r_ml Γ^l_ns - Γ^r_nl Γ^l_ms`.
fn riemann(
    cx: &mut Cx<'_>,
    g: &[Vec<NodeId>],
    vars: &[NodeId],
) -> Option<Tensor> {
    let n = vars.len();
    let gamma = christoffel_second(cx, g, vars)?;
    let mut out = Tensor::zeros(cx.graph, vec![n, n, n, n]);
    for r in 0..n {
        for s in 0..n {
            for m in 0..n {
                for v in 0..n {
                    if v <= m {
                        continue;
                    }
                    let first = d(cx, gamma.get(&[r, v, s]), vars[m])?;
                    let second = d(cx, gamma.get(&[r, m, s]), vars[v])?;
                    let mut terms = vec![first, neg(cx.graph, second)];
                    for l in 0..n {
                        let p = mul(cx.graph, &[gamma.get(&[r, m, l]), gamma.get(&[l, v, s])]);
                        let q = mul(cx.graph, &[gamma.get(&[r, v, l]), gamma.get(&[l, m, s])]);
                        terms.push(p);
                        terms.push(neg(cx.graph, q));
                    }
                    let sum = add(cx.graph, &terms);
                    let value = normal_form(cx, sum);
                    out.set(&[r, s, m, v], value);
                    // Antisymmetric in the last pair.
                    let negated = neg(cx.graph, value);
                    let negated = cx.simplify(negated);
                    out.set(&[r, s, v, m], negated);
                }
            }
        }
    }
    Some(out)
}

fn ricci(
    cx: &mut Cx<'_>,
    g: &[Vec<NodeId>],
    vars: &[NodeId],
) -> Option<Vec<Vec<NodeId>>> {
    let n = vars.len();
    let r = riemann(cx, g, vars)?;
    let mut out = vec![vec![NodeId::NONE; n]; n];
    #[allow(clippy::needless_range_loop)] // index is used for more than one array / arithmetic; iterator form would not be clearer
    for s in 0..n {
        #[allow(clippy::needless_range_loop)] // index is used for more than one array / arithmetic; iterator form would not be clearer
        for v in 0..n {
            let terms: Vec<NodeId> = (0..n).map(|k| r.get(&[k, s, k, v])).collect();
            let sum = add(cx.graph, &terms);
            out[s][v] = normal_form(cx, sum);
        }
    }
    Some(out)
}

fn ricci_scalar(
    cx: &mut Cx<'_>,
    g: &[Vec<NodeId>],
    vars: &[NodeId],
) -> Option<NodeId> {
    let n = vars.len();
    let ric = ricci(cx, g, vars)?;
    let inverse = invert(cx, g)?;
    let mut terms = Vec::new();
    for s in 0..n {
        for v in 0..n {
            terms.push(mul(cx.graph, &[inverse[s][v], ric[s][v]]));
        }
    }
    let sum = add(cx.graph, &terms);
    Some(normal_form(cx, sum))
}

// ----------------------------------------------------------------------
// Coordinate systems
// ----------------------------------------------------------------------

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum System {
    Cartesian,
    Polar,
    Cylindrical,
    Spherical,
}

fn system(
    graph: &Graph,
    node: NodeId,
) -> Option<System> {
    let symbol = graph.symbol_of(node)?;
    match graph.interner().symbol_name(symbol) {
        | "cartesian" => Some(System::Cartesian),
        | "polar" => Some(System::Polar),
        | "cylindrical" => Some(System::Cylindrical),
        | "spherical" => Some(System::Spherical),
        | _ => None,
    }
}

fn named(
    graph: &mut Graph,
    name: &str,
    args: &[NodeId],
) -> Option<NodeId> {
    let op = graph.ops().lookup(name)?;
    graph.try_node(op, args)
}

/// Cartesian coordinates as functions of the system's variables.
fn to_cartesian(
    graph: &mut Graph,
    sys: System,
    vars: &[NodeId],
) -> Option<Vec<NodeId>> {
    Some(match (sys, vars) {
        | (System::Cartesian, _) => vars.to_vec(),
        | (System::Polar, &[r, t]) => {
            let (c, s) = (named(graph, "cos", &[t])?, named(graph, "sin", &[t])?);
            vec![mul(graph, &[r, c]), mul(graph, &[r, s])]
        },
        | (System::Cylindrical, &[r, t, z]) => {
            let (c, s) = (named(graph, "cos", &[t])?, named(graph, "sin", &[t])?);
            vec![mul(graph, &[r, c]), mul(graph, &[r, s]), z]
        },
        | (System::Spherical, &[r, t, p]) => {
            let (ct, st) = (named(graph, "cos", &[t])?, named(graph, "sin", &[t])?);
            let (cp, sp) = (named(graph, "cos", &[p])?, named(graph, "sin", &[p])?);
            vec![mul(graph, &[r, st, cp]), mul(graph, &[r, st, sp]), mul(graph, &[r, ct])]
        },
        | _ => return None,
    })
}

/// The system's variables as functions of Cartesian coordinates.
fn from_cartesian(
    graph: &mut Graph,
    sys: System,
    xyz: &[NodeId],
) -> Option<Vec<NodeId>> {
    let radius = |graph: &mut Graph, parts: &[NodeId]| {
        let squares: Vec<NodeId> = parts.iter().map(|&p| powi(graph, p, 2)).collect();
        let sum = add(graph, &squares);
        sqrt(graph, sum)
    };
    Some(match (sys, xyz) {
        | (System::Cartesian, _) => xyz.to_vec(),
        | (System::Polar, &[x, y]) => vec![radius(graph, &[x, y]), named(graph, "atan2", &[y, x])?],
        | (System::Cylindrical, &[x, y, z]) => vec![radius(graph, &[x, y]), named(graph, "atan2", &[y, x])?, z],
        | (System::Spherical, &[x, y, z]) => {
            let r = radius(graph, &[x, y, z]);
            let inverse = powi(graph, r, -1);
            let ratio = mul(graph, &[z, inverse]);
            vec![r, named(graph, "acos", &[ratio])?, named(graph, "atan2", &[y, x])?]
        },
        | _ => return None,
    })
}

/// `J[a][i] = ∂x_a / ∂q_i` of the map to Cartesian coordinates.
fn jacobian(
    cx: &mut Cx<'_>,
    sys: System,
    vars: &[NodeId],
) -> Option<Vec<Vec<NodeId>>> {
    let cartesian = to_cartesian(cx.graph, sys, vars)?;
    let mut out = Vec::with_capacity(cartesian.len());
    for &x in &cartesian {
        let mut row = Vec::with_capacity(vars.len());
        for &q in vars {
            row.push(d(cx, x, q)?);
        }
        out.push(row);
    }
    Some(out)
}

fn metric(
    cx: &mut Cx<'_>,
    sys: System,
    vars: &[NodeId],
) -> Option<Vec<Vec<NodeId>>> {
    let j = jacobian(cx, sys, vars)?;
    let n = vars.len();
    let mut g = vec![vec![NodeId::NONE; n]; n];
    for a in 0..n {
        for b in 0..n {
            let terms: Vec<NodeId> = j.iter().map(|row| mul(cx.graph, &[row[a], row[b]])).collect();
            let sum = add(cx.graph, &terms);
            g[a][b] = cx.simplify(sum);
        }
    }
    Some(g)
}

/// Scale factors of an orthogonal system; `None` if the metric is not
/// diagonal.
fn scale_factors(
    cx: &mut Cx<'_>,
    sys: System,
    vars: &[NodeId],
) -> Option<Vec<NodeId>> {
    // The textbook scale factors, valid on the usual coordinate ranges
    // (r ≥ 0, 0 ≤ θ ≤ π), where sqrt(g_ii) has a sign-free closed form.
    let one = cx.graph.int(1);
    match (sys, vars) {
        | (System::Cartesian, _) => return Some(vec![one; vars.len()]),
        | (System::Polar, &[r, _]) | (System::Cylindrical, &[r, _, _]) => {
            let mut h = vec![one, r];
            if vars.len() == 3 {
                h.push(one);
            }
            return Some(h);
        },
        | (System::Spherical, &[r, t, _]) => {
            let sin = cx.graph.ops().lookup("sin")?;
            let s = cx.graph.node(sin, &[t]);
            let rs = mul(cx.graph, &[r, s]);
            return Some(vec![one, r, rs]);
        },
        | _ => {},
    }
    let g = metric(cx, sys, vars)?;
    for (a, row) in g.iter().enumerate() {
        for (b, &entry) in row.iter().enumerate() {
            if a != b && !cx.graph.number_of(entry).is_some_and(Number::is_zero) {
                return None;
            }
        }
    }
    let mut out = Vec::with_capacity(g.len());
    for (k, row) in g.iter().enumerate() {
        let root = sqrt(cx.graph, row[k]);
        out.push(cx.simplify(root));
    }
    Some(out)
}

/// Arguments `(from, from_vars, to, to_vars)` of the transform operators.
struct Change {
    from: System,
    from_vars: Vec<NodeId>,
    to: System,
    to_vars: Vec<NodeId>,
}

impl Change {
    fn read(
        graph: &mut Graph,
        args: &[NodeId],
    ) -> Option<Self> {
        let (&from, &from_vars, &to, &to_vars) = (args.get(1)?, args.get(2)?, args.get(3)?, args.get(4)?);
        Some(Self {
            from: system(graph, from)?,
            from_vars: vector(graph, from_vars)?,
            to: system(graph, to)?,
            to_vars: vector(graph, to_vars)?,
        })
    }

    /// The `from` variables as functions of the `to` variables.
    fn old_in_new(
        &self,
        graph: &mut Graph,
    ) -> Option<Vec<NodeId>> {
        let cartesian = to_cartesian(graph, self.to, &self.to_vars)?;
        from_cartesian(graph, self.from, &cartesian)
    }

    /// `f(from_vars)` rewritten in the `to` variables, simplified.
    fn rewrite(
        &self,
        cx: &mut Cx<'_>,
        f: NodeId,
    ) -> Option<NodeId> {
        let olds = self.old_in_new(cx.graph)?;
        let mut term = best(cx.graph, f)?;
        // Substitute through fresh symbols so that variables shared by both
        // systems (cylindrical z) are not substituted twice.
        let mut fresh = Vec::with_capacity(self.from_vars.len());
        for &v in &self.from_vars {
            let symbol = cx.graph.interner_mut().fresh_symbol("q");
            let node = cx.graph.symbol_node(symbol);
            term = cx.graph.substitute(term, v, node);
            fresh.push(node);
        }
        for (&node, &value) in fresh.iter().zip(&olds) {
            term = cx.graph.substitute(term, node, value);
        }
        Some(cx.simplify(term))
    }

    /// `Λ[a][i] = ∂q'_a / ∂q_i` in the `to` variables: `J_to^-1 J_from`.
    fn lambda(
        &self,
        cx: &mut Cx<'_>,
    ) -> Option<Vec<Vec<NodeId>>> {
        let j_to = jacobian(cx, self.to, &self.to_vars)?;
        let j_from = jacobian(cx, self.from, &self.from_vars)?;
        let j_from: Vec<Vec<NodeId>> =
            j_from.iter().map(|row| row.iter().map(|&e| self.rewrite(cx, e)).collect::<Option<_>>()).collect::<Option<_>>()?;
        let inverse = invert(cx, &j_to)?;
        product(cx, &inverse, &j_from)
    }
}

fn product(
    cx: &mut Cx<'_>,
    a: &[Vec<NodeId>],
    b: &[Vec<NodeId>],
) -> Option<Vec<Vec<NodeId>>> {
    let inner = b.len();
    if a.iter().any(|r| r.len() != inner) {
        return None;
    }
    let cols = b.first().map_or(0, Vec::len);
    let mut out = Vec::with_capacity(a.len());
    for row in a {
        let mut line = Vec::with_capacity(cols);
        #[allow(clippy::needless_range_loop)] // index is used for more than one array / arithmetic; iterator form would not be clearer
        for c in 0..cols {
            let terms: Vec<NodeId> = (0..inner).map(|k| mul(cx.graph, &[row[k], b[k][c]])).collect();
            let sum = add(cx.graph, &terms);
            line.push(normal_form(cx, sum));
        }
        out.push(line);
    }
    Some(out)
}

fn transpose(m: &[Vec<NodeId>]) -> Vec<Vec<NodeId>> {
    let cols = m.first().map_or(0, Vec::len);
    (0..cols).map(|c| m.iter().map(|row| row[c]).collect()).collect()
}

// ----------------------------------------------------------------------
// Differential forms
// ----------------------------------------------------------------------

/// Terms `(coefficient, sorted indices)` of a form.
fn read_form(
    graph: &mut Graph,
    node: NodeId,
) -> Option<Vec<(NodeId, Vec<usize>)>> {
    let term = best(graph, node)?;
    if graph.ops().get(graph.op(term)).name.as_ref() != "form" {
        return None;
    }
    let &[terms] = graph.children(term) else {
        return None;
    };
    let mut out = Vec::new();
    for entry in vector(graph, terms).or_else(|| graph.children(terms).first().map(|_| graph.children(terms).to_vec()))? {
        let &[c, idx] = graph.children(entry) else {
            return None;
        };
        let indices: Vec<usize> = graph.children(idx).iter().map(|&i| small_index(graph, i)).collect::<Option<_>>()?;
        out.push((c, indices));
    }
    Some(out)
}

/// Sorts an index list, returning the permutation sign, or `None` when an
/// index repeats (the wedge vanishes).
fn sort_indices(indices: &mut [usize]) -> Option<bool> {
    let mut even = true;
    for i in 0..indices.len() {
        for j in 0..indices.len() - 1 - i {
            if indices[j] == indices[j + 1] {
                return None;
            }
            if indices[j] > indices[j + 1] {
                indices.swap(j, j + 1);
                even = !even;
            }
        }
    }
    if indices.windows(2).any(|w| w[0] == w[1]) {
        return None;
    }
    Some(even)
}

fn write_form(
    cx: &mut Cx<'_>,
    terms: Vec<(NodeId, Vec<usize>)>,
) -> Option<NodeId> {
    // Combine equal index sets, drop zeros.
    let mut combined: Vec<(Vec<usize>, Vec<NodeId>)> = Vec::new();
    for (c, idx) in terms {
        match combined.iter_mut().find(|(i, _)| *i == idx) {
            | Some((_, cs)) => cs.push(c),
            | None => combined.push((idx, vec![c])),
        }
    }
    combined.sort_by(|a, b| a.0.cmp(&b.0));
    let mut entries = Vec::new();
    for (idx, cs) in combined {
        let sum = add(cx.graph, &cs);
        let c = cx.simplify(sum);
        if cx.graph.number_of(c).is_some_and(Number::is_zero) {
            continue;
        }
        let index_nodes: Vec<NodeId> = idx.iter().map(|&i| cx.graph.int(i64::try_from(i + 1).unwrap_or(0))).collect();
        let index_list = list(cx.graph, &index_nodes);
        entries.push(list(cx.graph, &[c, index_list]));
    }
    let terms = list(cx.graph, &entries);
    named(cx.graph, "form", &[terms])
}

// ----------------------------------------------------------------------
// The kernel
// ----------------------------------------------------------------------

struct Geometry {
    op: OpId,
    request: Request,
}

impl Kernel for Geometry {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let args = cx.graph.children(node).to_vec();
        self.compute(cx, &args).map_or(Outcome::Pass, Outcome::Equal)
    }

    fn revisit(&self) -> bool {
        // Arguments are often requests that reduce later.
        true
    }
}

impl Geometry {
    #[allow(clippy::too_many_lines)]
    fn compute(
        &self,
        cx: &mut Cx<'_>,
        args: &[NodeId],
    ) -> Option<NodeId> {
        let arg = |k: usize| args.get(k).copied();
        match self.request {
            | Request::Rank => {
                let t = Tensor::read(cx.graph, arg(0)?)?;
                Some(cx.graph.int(i64::try_from(t.shape.len()).ok()?))
            },
            | Request::Component => {
                let t = Tensor::read(cx.graph, arg(0)?)?;
                let index: Vec<usize> =
                    vector(cx.graph, arg(1)?)?.iter().map(|&i| small_index(cx.graph, i)).collect::<Option<_>>()?;
                (index.len() == t.shape.len() && index.iter().zip(&t.shape).all(|(i, n)| i < n)).then(|| t.get(&index))
            },
            | Request::Add => {
                let (a, b) = (Tensor::read(cx.graph, arg(0)?)?, Tensor::read(cx.graph, arg(1)?)?);
                if a.shape != b.shape {
                    return None;
                }
                let data = a.data.iter().zip(&b.data).map(|(&x, &y)| add(cx.graph, &[x, y])).collect();
                Some(Tensor { shape: a.shape, data }.write(cx.graph))
            },
            | Request::Smul => {
                let c = arg(0)?;
                let t = Tensor::read(cx.graph, arg(1)?)?;
                let data = t.data.iter().map(|&x| mul(cx.graph, &[c, x])).collect();
                Some(Tensor { shape: t.shape, data }.write(cx.graph))
            },
            | Request::Outer => {
                let (a, b) = (Tensor::read(cx.graph, arg(0)?)?, Tensor::read(cx.graph, arg(1)?)?);
                let mut data = Vec::with_capacity(a.data.len() * b.data.len());
                for &x in &a.data {
                    for &y in &b.data {
                        data.push(mul(cx.graph, &[x, y]));
                    }
                }
                let shape = [a.shape, b.shape].concat();
                Some(Tensor { shape, data }.write(cx.graph))
            },
            | Request::Contract => {
                let t = Tensor::read(cx.graph, arg(0)?)?;
                let (i, j) = (small_index(cx.graph, arg(1)?)?, small_index(cx.graph, arg(2)?)?);
                if i == j || i >= t.shape.len() || j >= t.shape.len() || t.shape[i] != t.shape[j] {
                    return None;
                }
                let rest: Vec<usize> =
                    t.shape.iter().enumerate().filter(|&(k, _)| k != i && k != j).map(|(_, &n)| n).collect();
                let mut out = Tensor::zeros(cx.graph, rest.clone());
                for index in indices(&rest) {
                    let mut terms = Vec::with_capacity(t.shape[i]);
                    for k in 0..t.shape[i] {
                        let mut full = Vec::with_capacity(t.shape.len());
                        let mut rest_iter = index.iter();
                        for position in 0..t.shape.len() {
                            full.push(if position == i || position == j { k } else { *rest_iter.next()? });
                        }
                        terms.push(t.get(&full));
                    }
                    let sum = add(cx.graph, &terms);
                    out.set(&index, sum);
                }
                Some(out.write(cx.graph))
            },
            | Request::Lower | Request::Raise => {
                let v = vector(cx.graph, arg(0)?)?;
                let g = square(cx.graph, arg(1)?)?;
                let m = if self.request == Request::Lower { g } else { invert(cx, &g)? };
                let column: Vec<Vec<NodeId>> = v.iter().map(|&x| vec![x]).collect();
                let result = product(cx, &m, &column)?;
                let flat: Vec<NodeId> = result.into_iter().map(|r| r[0]).collect();
                Some(list(cx.graph, &flat))
            },
            | Request::Christoffel1 | Request::Christoffel2 | Request::Riemann => {
                let g = square(cx.graph, arg(0)?)?;
                let vars = vector(cx.graph, arg(1)?)?;
                let t = match self.request {
                    | Request::Christoffel1 => christoffel_first(cx, &g, &vars)?,
                    | Request::Christoffel2 => christoffel_second(cx, &g, &vars)?,
                    | _ => riemann(cx, &g, &vars)?,
                };
                Some(t.write(cx.graph))
            },
            | Request::Ricci => {
                let g = square(cx.graph, arg(0)?)?;
                let vars = vector(cx.graph, arg(1)?)?;
                let ric = ricci(cx, &g, &vars)?;
                Some(matrix_term(cx.graph, &ric))
            },
            | Request::RicciScalar => {
                let g = square(cx.graph, arg(0)?)?;
                let vars = vector(cx.graph, arg(1)?)?;
                ricci_scalar(cx, &g, &vars)
            },
            | Request::Einstein => {
                let g = square(cx.graph, arg(0)?)?;
                let vars = vector(cx.graph, arg(1)?)?;
                let ric = ricci(cx, &g, &vars)?;
                let scalar = ricci_scalar(cx, &g, &vars)?;
                let h = half(cx.graph);
                let mut out = ric.clone();
                for (a, row) in out.iter_mut().enumerate() {
                    for (b, entry) in row.iter_mut().enumerate() {
                        let scaled = mul(cx.graph, &[h, scalar, g[a][b]]);
                        let difference = sub(cx.graph, ric[a][b], scaled);
                        *entry = normal_form(cx, difference);
                    }
                }
                Some(matrix_term(cx.graph, &out))
            },
            | Request::Covariant => {
                let v = vector(cx.graph, arg(0)?)?;
                let g = square(cx.graph, arg(1)?)?;
                let vars = vector(cx.graph, arg(2)?)?;
                let n = vars.len();
                if v.len() != n {
                    return None;
                }
                let gamma = christoffel_second(cx, &g, &vars)?;
                let mut out = vec![vec![NodeId::NONE; n]; n];
                for i in 0..n {
                    for k in 0..n {
                        let mut terms = vec![d(cx, v[i], vars[k])?];
                        #[allow(clippy::needless_range_loop)] // index is used for more than one array / arithmetic; iterator form would not be clearer
                        for j in 0..n {
                            terms.push(mul(cx.graph, &[gamma.get(&[i, k, j]), v[j]]));
                        }
                        let sum = add(cx.graph, &terms);
                        out[i][k] = normal_form(cx, sum);
                    }
                }
                Some(matrix_term(cx.graph, &out))
            },
            | Request::ToCartesian | Request::FromCartesian => {
                let sys = system(cx.graph, arg(0)?)?;
                let vars = vector(cx.graph, arg(1)?)?;
                let out = if self.request == Request::ToCartesian {
                    to_cartesian(cx.graph, sys, &vars)?
                } else {
                    from_cartesian(cx.graph, sys, &vars)?
                };
                Some(list(cx.graph, &out))
            },
            | Request::Metric => {
                let sys = system(cx.graph, arg(0)?)?;
                let vars = vector(cx.graph, arg(1)?)?;
                let g = metric(cx, sys, &vars)?;
                Some(matrix_term(cx.graph, &g))
            },
            | Request::ScaleFactors => {
                let sys = system(cx.graph, arg(0)?)?;
                let vars = vector(cx.graph, arg(1)?)?;
                let h = scale_factors(cx, sys, &vars)?;
                Some(list(cx.graph, &h))
            },
            | Request::TransformPoint => {
                let (from, to) = (system(cx.graph, arg(0)?)?, system(cx.graph, arg(1)?)?);
                let point = vector(cx.graph, arg(2)?)?;
                let cartesian = to_cartesian(cx.graph, from, &point)?;
                let out = from_cartesian(cx.graph, to, &cartesian)?;
                Some(list(cx.graph, &out))
            },
            | Request::TransformExpression => {
                let change = Change::read(cx.graph, args)?;
                change.rewrite(cx, arg(0)?)
            },
            | Request::TransformVector | Request::TransformCovector | Request::TransformTensor2 => {
                let change = Change::read(cx.graph, args)?;
                let lambda = change.lambda(cx)?;
                match self.request {
                    | Request::TransformVector => {
                        // V'^a = Λ^a_i V^i
                        let v = vector(cx.graph, arg(0)?)?;
                        let column: Vec<Vec<NodeId>> =
                            v.iter().map(|&x| change.rewrite(cx, x).map(|y| vec![y])).collect::<Option<_>>()?;
                        let result = product(cx, &lambda, &column)?;
                        let flat: Vec<NodeId> = result.into_iter().map(|r| r[0]).collect();
                        Some(list(cx.graph, &flat))
                    },
                    | Request::TransformCovector => {
                        // W'_a = (Λ^-1)^i_a W_i
                        let w = vector(cx.graph, arg(0)?)?;
                        let inverse = invert(cx, &lambda)?;
                        let row: Vec<NodeId> = w.iter().map(|&x| change.rewrite(cx, x)).collect::<Option<_>>()?;
                        let result = product(cx, &[row], &inverse)?;
                        Some(list(cx.graph, &result[0]))
                    },
                    | _ => {
                        // T'^ab = Λ^a_i Λ^b_j T^ij
                        let t = square(cx.graph, arg(0)?)?;
                        let t: Vec<Vec<NodeId>> = t
                            .iter()
                            .map(|row| row.iter().map(|&x| change.rewrite(cx, x)).collect::<Option<_>>())
                            .collect::<Option<_>>()?;
                        let left = product(cx, &lambda, &t)?;
                        let result = product(cx, &left, &transpose(&lambda))?;
                        Some(matrix_term(cx.graph, &result))
                    },
                }
            },
            | Request::GradIn | Request::DivIn | Request::CurlIn | Request::LaplacianIn => {
                let sys = system(cx.graph, arg(1)?)?;
                let vars = vector(cx.graph, arg(2)?)?;
                let h = scale_factors(cx, sys, &vars)?;
                let n = vars.len();
                let all = mul(cx.graph, &h);
                match self.request {
                    | Request::GradIn => {
                        let f = arg(0)?;
                        let mut out = Vec::with_capacity(n);
                        for k in 0..n {
                            let df = d(cx, f, vars[k])?;
                            let inverse = powi(cx.graph, h[k], -1);
                            let component = mul(cx.graph, &[inverse, df]);
                            out.push(normal_form(cx, component));
                        }
                        Some(list(cx.graph, &out))
                    },
                    | Request::DivIn => {
                        // 1/(h1..hn) Σ ∂_k (F_k h1..hn / h_k)
                        let field = vector(cx.graph, arg(0)?)?;
                        if field.len() != n {
                            return None;
                        }
                        let mut terms = Vec::with_capacity(n);
                        for k in 0..n {
                            let inverse = powi(cx.graph, h[k], -1);
                            let flux = mul(cx.graph, &[field[k], all, inverse]);
                            terms.push(d(cx, flux, vars[k])?);
                        }
                        let sum = add(cx.graph, &terms);
                        let inverse = powi(cx.graph, all, -1);
                        let value = mul(cx.graph, &[inverse, sum]);
                        Some(normal_form(cx, value))
                    },
                    | Request::LaplacianIn => {
                        let f = arg(0)?;
                        let mut terms = Vec::with_capacity(n);
                        for k in 0..n {
                            let df = d(cx, f, vars[k])?;
                            let h2 = powi(cx.graph, h[k], -2);
                            let flux = mul(cx.graph, &[all, h2, df]);
                            terms.push(d(cx, flux, vars[k])?);
                        }
                        let sum = add(cx.graph, &terms);
                        let inverse = powi(cx.graph, all, -1);
                        let value = mul(cx.graph, &[inverse, sum]);
                        Some(normal_form(cx, value))
                    },
                    | _ => {
                        // (curl F)_1 = 1/(h2 h3) [∂_2 (h3 F3) - ∂_3 (h2 F2)], cyclic
                        let field = vector(cx.graph, arg(0)?)?;
                        if n != 3 || field.len() != 3 {
                            return None;
                        }
                        let mut out = Vec::with_capacity(3);
                        for a in 0..3 {
                            let (b, c) = ((a + 1) % 3, (a + 2) % 3);
                            let hc_fc = mul(cx.graph, &[h[c], field[c]]);
                            let hb_fb = mul(cx.graph, &[h[b], field[b]]);
                            let first = d(cx, hc_fc, vars[b])?;
                            let second = d(cx, hb_fb, vars[c])?;
                            let difference = sub(cx.graph, first, second);
                            let scale = mul(cx.graph, &[h[b], h[c]]);
                            let inverse = powi(cx.graph, scale, -1);
                            let value = mul(cx.graph, &[inverse, difference]);
                            out.push(normal_form(cx, value));
                        }
                        Some(list(cx.graph, &out))
                    },
                }
            },
            | Request::ExteriorD => {
                let terms = read_form(cx.graph, arg(0)?)?;
                let vars = vector(cx.graph, arg(1)?)?;
                let mut out = Vec::new();
                for (c, idx) in terms {
                    for (j, &x) in vars.iter().enumerate() {
                        let dc = d(cx, c, x)?;
                        if cx.graph.number_of(dc).is_some_and(Number::is_zero) {
                            continue;
                        }
                        let mut full = vec![j];
                        full.extend(&idx);
                        let Some(even) = sort_indices(&mut full) else {
                            continue;
                        };
                        out.push((if even { dc } else { neg(cx.graph, dc) }, full));
                    }
                }
                write_form(cx, out)
            },
            | Request::Wedge => {
                let (a, b) = (read_form(cx.graph, arg(0)?)?, read_form(cx.graph, arg(1)?)?);
                let mut out = Vec::new();
                for (c1, i1) in &a {
                    for (c2, i2) in &b {
                        let mut full = [i1.clone(), i2.clone()].concat();
                        let Some(even) = sort_indices(&mut full) else {
                            continue;
                        };
                        let product = mul(cx.graph, &[*c1, *c2]);
                        out.push((if even { product } else { neg(cx.graph, product) }, full));
                    }
                }
                write_form(cx, out)
            },
            | Request::Greens => {
                let (p, q) = (arg(0)?, arg(1)?);
                let vars = vector(cx.graph, arg(2)?)?;
                let &[x, y] = vars.as_slice() else {
                    return None;
                };
                let qx = d(cx, q, x)?;
                let py = d(cx, p, y)?;
                let integrand = sub(cx.graph, qx, py);
                named(cx.graph, "volume_integral", &[integrand, arg(2)?, arg(3)?])
            },
            | Request::Gauss => {
                let divergence = named(cx.graph, "div", &[arg(0)?, arg(1)?])?;
                let divergence = cx.simplify(divergence);
                named(cx.graph, "volume_integral", &[divergence, arg(1)?, arg(2)?])
            },
            | Request::Stokes => {
                // ∬ curl F(r(u, v)) · (r_u × r_v) du dv
                let field = vector(cx.graph, arg(0)?)?;
                let surface = vector(cx.graph, arg(1)?)?;
                let (u, v) = (arg(2)?, arg(3)?);
                if field.len() != 3 || surface.len() != 3 {
                    return None;
                }
                let xyz: Vec<NodeId> = ["x", "y", "z"].iter().map(|name| cx.graph.sym(name)).collect();
                let field_node = list(cx.graph, &field);
                let xyz_node = list(cx.graph, &xyz);
                let curl = named(cx.graph, "curl", &[field_node, xyz_node])?;
                let curl = cx.simplify(curl);
                let curl = vector(cx.graph, curl)?;
                let mut ru = Vec::with_capacity(3);
                let mut rv = Vec::with_capacity(3);
                for &c in &surface {
                    ru.push(d(cx, c, u)?);
                    rv.push(d(cx, c, v)?);
                }
                let normal = [
                    sub_products(cx.graph, ru[1], rv[2], ru[2], rv[1]),
                    sub_products(cx.graph, ru[2], rv[0], ru[0], rv[2]),
                    sub_products(cx.graph, ru[0], rv[1], ru[1], rv[0]),
                ];
                let mut terms = Vec::with_capacity(3);
                for (k, &component) in curl.iter().enumerate() {
                    let mut along = component;
                    for (&symbol, &value) in xyz.iter().zip(&surface) {
                        along = cx.graph.substitute(along, symbol, value);
                    }
                    terms.push(mul(cx.graph, &[along, normal[k]]));
                }
                let integrand = add(cx.graph, &terms);
                let integrand = cx.simplify(integrand);
                let inner = named(cx.graph, "defint", &[integrand, u, arg(4)?, arg(5)?])?;
                named(cx.graph, "defint", &[inner, v, arg(6)?, arg(7)?])
            },
            | Request::Boundary => {
                let mut faces = Vec::new();
                for (sign, cell) in read_chain(cx.graph, arg(0)?)? {
                    for (s, face) in cell_boundary(cx, &cell)? {
                        faces.push((sign * s, face));
                    }
                }
                Some(write_chain(cx.graph, &faces))
            },
            | Request::Pullback => {
                let terms = read_form(cx.graph, arg(0)?)?;
                let vars = items(cx.graph, arg(1)?)?;
                let chart = items(cx.graph, arg(2)?)?;
                let params = items(cx.graph, arg(3)?)?;
                let pulled = pullback(cx, &terms, &vars, &chart, &params)?;
                write_form(cx, pulled)
            },
            | Request::IntegrateForm => {
                let terms = read_form(cx.graph, arg(0)?)?;
                let vars = items(cx.graph, arg(1)?)?;
                let mut total = Vec::new();
                for (sign, cell) in read_chain(cx.graph, arg(2)?)? {
                    let value = integrate_over_cell(cx, &terms, &vars, &cell)?;
                    let m = cx.graph.int(sign);
                    total.push(mul(cx.graph, &[m, value]));
                }
                let sum = add(cx.graph, &total);
                Some(cx.simplify(sum))
            },
            | Request::GeneralizedStokes => {
                let (omega, vars, manifold) = (arg(0)?, arg(1)?, arg(2)?);
                let d_omega = named(cx.graph, "exterior_d", &[omega, vars])?;
                let d_omega = cx.simplify(d_omega);
                let boundary = named(cx.graph, "boundary", &[manifold])?;
                let boundary = cx.simplify(boundary);
                let lhs = named(cx.graph, "integrate_form", &[d_omega, vars, manifold])?;
                let rhs = named(cx.graph, "integrate_form", &[omega, vars, boundary])?;
                let (lhs, rhs) = (cx.simplify(lhs), cx.simplify(rhs));
                Some(cx.graph.node(core::EQ, &[lhs, rhs]))
            },
        }
    }
}

// ----------------------------------------------------------------------
// Chains, boundaries and integration of forms
// ----------------------------------------------------------------------

/// The items of a `list(...)` (possibly empty).
fn items(
    graph: &mut Graph,
    node: NodeId,
) -> Option<Vec<NodeId>> {
    let node = best(graph, node)?;
    (graph.op(node) == core::LIST).then(|| graph.children(node).to_vec())
}

/// A cell: chart, parameters and `(a, b)` bounds.
struct Cell {
    chart: Vec<NodeId>,
    params: Vec<NodeId>,
    bounds: Vec<(NodeId, NodeId)>,
}

fn read_cell(
    graph: &mut Graph,
    node: NodeId,
) -> Option<Cell> {
    let node = best(graph, node)?;
    if graph.ops().get(graph.op(node)).name.as_ref() != "cell" {
        return None;
    }
    let &[chart, params, bounds] = graph.children(node) else {
        return None;
    };
    let (chart, params) = (items(graph, chart)?, items(graph, params)?);
    let mut pairs = Vec::new();
    for b in items(graph, bounds)? {
        let &[lo, hi] = graph.children(b) else {
            return None;
        };
        pairs.push((lo, hi));
    }
    (pairs.len() == params.len()).then_some(Cell { chart, params, bounds: pairs })
}

fn write_cell(
    graph: &mut Graph,
    cell: &Cell,
) -> NodeId {
    let chart = list(graph, &cell.chart);
    let params = list(graph, &cell.params);
    let bounds: Vec<NodeId> = cell.bounds.iter().map(|&(a, b)| list(graph, &[a, b])).collect();
    let bounds = list(graph, &bounds);
    named(graph, "cell", &[chart, params, bounds]).unwrap_or(chart)
}

/// A cell (multiplicity 1) or a chain, as `(multiplicity, cell)` pairs.
fn read_chain(
    graph: &mut Graph,
    node: NodeId,
) -> Option<Vec<(i64, Cell)>> {
    let node = best(graph, node)?;
    if graph.ops().get(graph.op(node)).name.as_ref() != "chain" {
        return Some(vec![(1, read_cell(graph, node)?)]);
    }
    let &[terms] = graph.children(node) else {
        return None;
    };
    let mut out = Vec::new();
    for t in items(graph, terms)? {
        let &[m, c] = graph.children(t) else {
            return None;
        };
        let m = graph.number_of(m)?.to_i64()?;
        out.push((m, read_cell(graph, c)?));
    }
    Some(out)
}

fn write_chain(
    graph: &mut Graph,
    cells: &[(i64, Cell)],
) -> NodeId {
    let mut terms = Vec::new();
    for (m, c) in cells {
        if *m == 0 {
            continue;
        }
        let m = graph.int(*m);
        let c = write_cell(graph, c);
        terms.push(list(graph, &[m, c]));
    }
    let terms = list(graph, &terms);
    named(graph, "chain", &[terms]).unwrap_or(terms)
}

/// `∂` of a cell with constant bounds: the faces `u_i = a_i` (sign
/// `(-1)^i`) and `u_i = b_i` (sign `(-1)^(i+1)`), `i` from 1.
fn cell_boundary(
    cx: &mut Cx<'_>,
    cell: &Cell,
) -> Option<Vec<(i64, Cell)>> {
    let constant = cell.bounds.iter().all(|&(a, b)| {
        cell.params.iter().all(|&u| {
            let Some(sym) = cx.graph.as_symbol(u) else {
                return false;
            };
            !cx.graph.depends_on(cx.graph.find(a), sym) && !cx.graph.depends_on(cx.graph.find(b), sym)
        })
    });
    if !constant || cell.params.is_empty() {
        return None;
    }
    let mut faces = Vec::new();
    for (i, &u) in cell.params.iter().enumerate() {
        let (a, b) = cell.bounds[i];
        let parity = if i % 2 == 0 { 1 } else { -1 }; // (-1)^(i+1) with i from 1
        for (value, sign) in [(a, -parity), (b, parity)] {
            let chart: Vec<NodeId> = cell
                .chart
                .iter()
                .map(|&c| {
                    let s = cx.graph.substitute(c, u, value);
                    cx.simplify(s)
                })
                .collect();
            let mut params = cell.params.clone();
            params.remove(i);
            let mut bounds = cell.bounds.clone();
            bounds.remove(i);
            faces.push((sign, Cell { chart, params, bounds }));
        }
    }
    Some(faces)
}

/// `φ*ω` in the parameters `params`.
fn pullback(
    cx: &mut Cx<'_>,
    terms: &[(NodeId, Vec<usize>)],
    vars: &[NodeId],
    chart: &[NodeId],
    params: &[NodeId],
) -> Option<Vec<(NodeId, Vec<usize>)>> {
    if chart.len() != vars.len() {
        return None;
    }
    // Simultaneous substitution through fresh symbols.
    let fresh: Vec<NodeId> = (0..vars.len())
        .map(|_| {
            let s = cx.graph.interner_mut().fresh_symbol("p");
            cx.graph.symbol_node(s)
        })
        .collect();
    let mut jacobian = Vec::with_capacity(chart.len());
    for &c in chart {
        let row: Vec<NodeId> = params.iter().map(|&u| d(cx, c, u)).collect::<Option<_>>()?;
        jacobian.push(row);
    }
    let mut out = Vec::new();
    for (c, idx) in terms {
        let mut coefficient = *c;
        for (&x, &f) in vars.iter().zip(&fresh) {
            coefficient = cx.graph.substitute(coefficient, x, f);
        }
        for (&f, &value) in fresh.iter().zip(chart) {
            coefficient = cx.graph.substitute(coefficient, f, value);
        }
        // Wedge the pulled-back differentials dx_i one at a time.
        let mut partial: Vec<(NodeId, Vec<usize>)> = vec![(coefficient, Vec::new())];
        for &i in idx {
            let row = jacobian.get(i)?;
            let mut next = Vec::new();
            for (coef, js) in &partial {
                for (j, &entry) in row.iter().enumerate() {
                    if cx.graph.number_of(entry).is_some_and(Number::is_zero) {
                        continue;
                    }
                    let mut full = js.clone();
                    full.push(j);
                    let Some(even) = sort_indices(&mut full) else {
                        continue;
                    };
                    let product = mul(cx.graph, &[*coef, entry]);
                    next.push((if even { product } else { neg(cx.graph, product) }, full));
                }
            }
            partial = next;
        }
        out.extend(partial);
    }
    Some(out)
}

/// `∫_cell ω`.
fn integrate_over_cell(
    cx: &mut Cx<'_>,
    terms: &[(NodeId, Vec<usize>)],
    vars: &[NodeId],
    cell: &Cell,
) -> Option<NodeId> {
    let pulled = pullback(cx, terms, vars, &cell.chart, &cell.params)?;
    let top: Vec<usize> = (0..cell.params.len()).collect();
    let coefficients: Vec<NodeId> = pulled.into_iter().filter(|(_, idx)| *idx == top).map(|(c, _)| c).collect();
    let mut integrand = add(cx.graph, &coefficients);
    integrand = cx.simplify(integrand);
    for (&u, &(a, b)) in cell.params.iter().zip(&cell.bounds) {
        let request = named(cx.graph, "defint", &[integrand, u, a, b])?;
        integrand = cx.simplify(request);
    }
    Some(integrand)
}

/// `a b - c d`.
fn sub_products(
    graph: &mut Graph,
    a: NodeId,
    b: NodeId,
    c: NodeId,
    d: NodeId,
) -> NodeId {
    let first = mul(graph, &[a, b]);
    let second = mul(graph, &[c, d]);
    sub(graph, first, second)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Budget;
    use crate::graph::Engine;
    use crate::graph::Env;
    use crate::graph::Extractor;
    use crate::graph::Facts;
    use crate::graph::Saturate;
    use crate::graph::SizeCost;
    use crate::rules::testing::numeric;
    use crate::rules::testing::simplify;

    fn run(src: &str) -> String {
        simplify(&[geometry()], src)
    }

    /// Runs with `r` assumed positive.
    fn run_positive(src: &str) -> String {
        let mut g = Graph::new();
        let engine = Engine::install(&mut g, &[geometry()]).unwrap_or_else(|e| panic!("{e}"));
        for name in ["r", "a"] {
            let s = g.interner_mut().symbol(name);
            g.assume(s, Facts::POSITIVE);
        }
        let root = g.parse(src).unwrap_or_else(|e| panic!("{e}"));
        engine.run(&mut g, &[root], &Env::symbolic(), &Saturate, &Budget::default());
        let term = Extractor::new(&g, &[root], &SizeCost).build(&mut g, root).unwrap_or(root);
        g.display(term)
    }

    #[test]
    fn tensor_algebra() {
        assert_eq!(run("tensor_rank(list(list(list(1, 2)), list(list(3, 4))))"), "3");
        assert_eq!(run("component(list(list(a, b), list(c, d)), list(2, 1))"), "c");
        assert_eq!(run("tensor_add(list(1, x), list(2, y))"), "list(3, x + y)");
        assert_eq!(run("tensor_smul(2, list(list(1, x)))"), "list(list(2, 2*x))");
        assert_eq!(run("tensor_outer(list(1, 2), list(x, y))"), "list(list(x, y), list(2*x, 2*y))");
        assert_eq!(run("contract(list(list(a, b), list(c, d)), 1, 2)"), "a + d");
        assert_eq!(run("contract(tensor_outer(list(1, 2), list(3, 4)), 1, 2)"), "11");
        assert_eq!(run("lower_index(list(1, 2), list(list(1, 0), list(0, r^2)))"), "list(1, 2*r^2)");
        assert_eq!(run("raise_index(list(1, 2), list(list(1, 0), list(0, r^2)))"), "list(1, 2/r^2)");
    }

    #[test]
    fn christoffel_symbols_of_the_plane_in_polar_coordinates() {
        let g = "list(list(1, 0), list(0, r^2))";
        // Γ^r_θθ = -r, Γ^θ_rθ = Γ^θ_θr = 1/r
        assert_eq!(
            run(&format!("christoffel({g}, list(r, theta))")),
            "list(list(list(0, 0), list(0, -r)), list(list(0, 1/r), list(1/r, 0)))"
        );
        assert_eq!(run(&format!("ricci_scalar({g}, list(r, theta))")), "0");
    }

    #[test]
    fn curvature_of_the_sphere() {
        let g = "list(list(a^2, 0), list(0, a^2*sin(theta)^2))";
        let vars = "list(theta, phi)";
        assert_eq!(run_positive(&format!("ricci_scalar({g}, {vars})")), "2/a^2");
        assert_eq!(run_positive(&format!("component(riemann({g}, {vars}), list(1, 2, 1, 2))")), "sin(theta)^2");
        assert_eq!(run_positive(&format!("ricci({g}, {vars})")), "list(list(1, 0), list(0, sin(theta)^2))");
        // In two dimensions the Einstein tensor vanishes identically.
        assert_eq!(run_positive(&format!("einstein({g}, {vars})")), "list(list(0, 0), list(0, 0))");
    }

    #[test]
    fn covariant_derivative_of_a_constant_cartesian_field() {
        // The field e_x in polar components: (cos θ, -sin θ / r) is
        // parallel, so its covariant derivative vanishes.
        let text = run_positive("covariant_derivative(list(cos(theta), -sin(theta)/r), list(list(1, 0), list(0, r^2)), list(r, theta))");
        assert_eq!(text, "list(list(0, 0), list(0, 0))");
    }

    #[test]
    fn coordinate_systems() {
        assert_eq!(run("to_cartesian(polar, list(r, t))"), "list(r*cos(t), r*sin(t))");
        assert_eq!(run_positive("coordinate_metric(spherical, list(r, t, p))"), "list(list(1, 0, 0), list(0, r^2, 0), list(0, 0, r^2*sin(t)^2))");
        assert_eq!(run_positive("scale_factors(cylindrical, list(r, t, z))"), "list(1, r, 1)");
        assert_eq!(run("transform_point(cartesian, polar, list(0, 2))"), "list(2, 1/2*pi)");
        assert_eq!(run("transform_point(polar, cartesian, list(2, pi))"), "list(-2, 0)");
        let (value, _) = numeric(&[geometry()], "component(transform_point(spherical, cylindrical, list(2, 1, 3)), list(1))", &[], 1e-12);
        assert!((value - 2.0 * 1.0_f64.sin()).abs() < 1e-12);
    }

    #[test]
    fn transformations_of_fields() {
        assert_eq!(run_positive("transform_expression(x^2 + y^2, cartesian, list(x, y), polar, list(r, t))"), "r^2");
        // The Cartesian field (x, y) is r e_r: polar components (r, 0).
        assert_eq!(run_positive("transform_vector(list(x, y), cartesian, list(x, y), polar, list(r, t))"), "list(r, 0)");
        // The gradient of x is a covector (1, 0): polar (cos t, -r sin t).
        assert_eq!(run_positive("transform_covector(list(1, 0), cartesian, list(x, y), polar, list(r, t))"), "list(cos(t), -r*sin(t))");
        // The Euclidean metric's inverse transforms to diag(1, 1/r^2).
        assert_eq!(
            run_positive("transform_tensor2(list(list(1, 0), list(0, 1)), cartesian, list(x, y), polar, list(r, t))"),
            "list(list(1, 0), list(0, 1/r^2))"
        );
    }

    #[test]
    fn vector_calculus_in_curvilinear_coordinates() {
        assert_eq!(run_positive("grad_in(r^2*cos(t), polar, list(r, t))"), "list(2*r*cos(t), -r*sin(t))");
        assert_eq!(run_positive("div_in(list(r, 0, 0), spherical, list(r, t, p))"), "3");
        assert_eq!(run_positive("laplacian_in(1/r, spherical, list(r, t, p))"), "0");
        assert_eq!(run_positive("laplacian_in(r^2, cylindrical, list(r, t, z))"), "4");
        // The field r e_phi in cylindrical coordinates is a rigid rotation:
        // curl = 2 e_z.
        assert_eq!(run_positive("curl_in(list(0, r, 0), cylindrical, list(r, t, z))"), "list(0, 0, 2)");
    }

    #[test]
    fn differential_forms() {
        // d(x dy) = dx ∧ dy
        assert_eq!(run("exterior_d(form(list(list(x, list(2)))), list(x, y))"), "form(list(list(1, list(1, 2))))");
        // d(df) = 0
        let df = run("exterior_d(form(list(list(x^2*y*z, list()))), list(x, y, z))");
        assert_eq!(run(&format!("exterior_d({df}, list(x, y, z))")), "form(list())");
        // dx ∧ dy = -dy ∧ dx, dx ∧ dx = 0
        assert_eq!(run("wedge(form(list(list(1, list(2)))), form(list(list(1, list(1)))))"), "form(list(list(-1, list(1, 2))))");
        assert_eq!(run("wedge(form(list(list(1, list(1)))), form(list(list(1, list(1)))))"), "form(list())");
    }

    #[test]
    fn integral_theorems() {
        // Circulation of (-y, x) around the unit square: 2 * area = 2.
        assert_eq!(run("greens_theorem(-y, x, list(x, y), list(list(0, 1), list(0, 1)))"), "2");
        // Flux of (x, y, z) out of the unit cube: 3.
        assert_eq!(run("gauss_theorem(list(x, y, z), list(x, y, z), list(list(0, 1), list(0, 1), list(0, 1)))"), "3");
        // Stokes on the unit disc (z = 0): circulation of (-y, x, 0) = 2 pi.
        assert_eq!(run("stokes_theorem(list(-y, x, 0), list(u*cos(v), u*sin(v), 0), u, v, 0, 1, 0, 2*pi)"), "2*pi");
    }

    #[test]
    fn chains_boundaries_and_the_generalized_stokes_theorem() {
        let square = "cell(list(x, y), list(x, y), list(list(0, 1), list(0, 1)))";
        // ∂ of the unit square: four oriented edges.
        let edges = run(&format!("boundary({square})"));
        assert_eq!(edges.matches("cell(").count(), 4, "{edges}");
        // ∫_∂I² x dy = ∫_I² dx∧dy = 1.
        let x_dy = "form(list(list(x, list(2))))";
        assert_eq!(run(&format!("integrate_form({x_dy}, list(x, y), boundary({square}))")), "1");
        assert_eq!(run(&format!("generalized_stokes({x_dy}, list(x, y), {square})")), "1 = 1");
        // ∂∂ = 0: a 0-form integrates to zero over the boundary of a boundary.
        assert_eq!(run(&format!("integrate_form(form(list(list(x^2*y + 3, list()))), list(x, y), boundary(boundary({square})))")), "0");
        // The unit disc in polar parameters, ω = -y dx + x dy: both sides 2 pi.
        let disc = "cell(list(r*cos(t), r*sin(t)), list(r, t), list(list(0, 1), list(0, 2*pi)))";
        let omega = "form(list(list(-y, list(1)), list(x, list(2))))";
        assert_eq!(run(&format!("generalized_stokes({omega}, list(x, y), {disc})")), "2*pi = 2*pi");
        // A 3-cell: the unit cube, ω = x dy∧dz (flux of (x, 0, 0)) = 1.
        let cube = "cell(list(x, y, z), list(x, y, z), list(list(0, 1), list(0, 1), list(0, 1)))";
        assert_eq!(run(&format!("generalized_stokes(form(list(list(x, list(2, 3)))), list(x, y, z), {cube})")), "1 = 1");
        // A curved surface in R³: the upper hemisphere, ω = x dy (circulation of (0, x, 0)) = pi.
        let hemisphere = "cell(list(sin(p)*cos(t), sin(p)*sin(t), cos(p)), list(p, t), list(list(0, pi/2), list(0, 2*pi)))";
        assert_eq!(run(&format!("generalized_stokes(form(list(list(x, list(2)))), list(x, y, z), {hemisphere})")), "pi = pi");
        // Pull-back of dx∧dy under polar coordinates: r dr∧dt.
        assert_eq!(run("pullback(form(list(list(1, list(1, 2)))), list(x, y), list(r*cos(t), r*sin(t)), list(r, t))"), "form(list(list(r, list(1, 2))))");
    }

    #[test]
    fn christoffel_symbols_of_the_first_kind() {
        // The plane in polar coordinates: Γ_{rθθ} = -r, Γ_{θrθ} = Γ_{θθr} = r.
        let g = "list(list(1, 0), list(0, r^2))";
        assert_eq!(
            run(&format!("christoffel1({g}, list(r, theta))")),
            "list(list(list(0, 0), list(0, -r)), list(list(0, r), list(r, 0)))"
        );
        // A constant metric has none.
        assert_eq!(
            run("christoffel1(list(list(1, 2), list(2, 5)), list(x, y))"),
            "list(list(list(0, 0), list(0, 0)), list(list(0, 0), list(0, 0)))"
        );
        // The sphere: Γ_{θφφ} = -a² sinθ cosθ, Γ_{φθφ} = a² sinθ cosθ.
        let sphere = "list(list(a^2, 0), list(0, a^2*sin(theta)^2))";
        let vars = "list(theta, phi)";
        let value = |i: u8, j: u8, k: u8| {
            numeric(&[geometry()], &format!("component(christoffel1({sphere}, {vars}), list({i}, {j}, {k}))"), &[("a", 2.0), ("theta", 0.7)], 1e-12).0
        };
        let expected = 4.0 * 0.7_f64.sin() * 0.7_f64.cos();
        assert!((value(1, 2, 2) + expected).abs() < 1e-12);
        assert!((value(2, 1, 2) - expected).abs() < 1e-12);
        assert!((value(2, 2, 1) - expected).abs() < 1e-12);
        assert!(value(1, 1, 1).abs() < 1e-12);
        // Γ^φ_θφ = g^{φφ} Γ_{φθφ} = cot θ.
        let second = numeric(&[geometry()], &format!("component(christoffel({sphere}, {vars}), list(2, 1, 2))"), &[("a", 2.0), ("theta", 0.7)], 1e-12).0;
        assert!((second - 1.0 / 0.7_f64.tan()).abs() < 1e-12);
        assert!((second - value(2, 1, 2) / (4.0 * 0.7_f64.sin().powi(2))).abs() < 1e-12);
    }
}
