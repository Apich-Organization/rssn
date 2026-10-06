//! Further Riemannian geometry: geodesics, curvature invariants, Lie
//! brackets and Killing fields, the Laplace–Beltrami operator.
//!
//! | operator | value |
//! |---|---|
//! | `metric_geodesic_acceleration(g, vars, v)` | `a^k = -Γ^k_ij v^i v^j`: the coordinate acceleration of a geodesic with velocity `v` |
//! | `kretschmann(g, vars)` | the curvature invariant `R_abcd R^abcd` |
//! | `gaussian_curvature(g, vars)` | `R / 2` of a two-dimensional metric |
//! | `volume_element(g, vars)` | `sqrt(det g)` |
//! | `laplace_beltrami(f, g, vars)` | `(1/√g) ∂_i (√g g^ij ∂_j f)` |
//! | `covariant_divergence(V, g, vars)` | `(1/√g) ∂_i (√g V^i)` |
//! | `killing_tensor(xi, g, vars)` | the Lie derivative `L_ξ g_ab = ξ^c ∂_c g_ab + g_cb ∂_a ξ^c + g_ac ∂_b ξ^c` (a matrix) |
//! | `is_killing(xi, g, vars)` | whether `L_ξ g` vanishes: `ξ` generates an isometry |

use super::add;
use super::christoffel_second;
use super::d;
use super::mul;
use super::neg;
use super::riemann;
use super::ricci_scalar;
use super::square;
use super::vector;
use super::Tensor;
use crate::graph::op::core;
use crate::graph::rule::Installer;
use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::OpDescriptor;
use crate::graph::OpFlags;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::Payload;
use crate::graph::RuleError;
use crate::graph::Tier;
use crate::rules::linalg::invert;
use crate::rules::linalg::normal_form;

#[derive(Copy, Clone)]
enum Kind {
    Geodesic,
    Kretschmann,
    Gaussian,
    Volume,
    Beltrami,
    Divergence,
    KillingTensor,
    IsKilling,
}

struct Extension {
    op: OpId,
    kind: Kind,
}

fn simplified_sum(
    cx: &mut Cx<'_>,
    terms: &[NodeId],
) -> NodeId {
    let total = add(cx.graph, terms);
    normal_form(cx, total)
}

fn volume(
    cx: &mut Cx<'_>,
    g: &[Vec<NodeId>],
) -> Option<NodeId> {
    let det = cx.graph.ops().lookup("det")?;
    let rows: Vec<NodeId> = g.iter().map(|r| cx.graph.node(core::LIST, r)).collect();
    let matrix = cx.graph.node(core::LIST, &rows);
    let request = cx.graph.try_node(det, &[matrix])?;
    let determinant = cx.simplify(request);
    let half = cx.graph.num(crate::graph::Number::fraction(1, 2)?);
    let root = cx.graph.node(core::POW, &[determinant, half]);
    Some(cx.simplify(root))
}

fn kretschmann(
    cx: &mut Cx<'_>,
    g: &[Vec<NodeId>],
    vars: &[NodeId],
) -> Option<NodeId> {
    let n = vars.len();
    let r = riemann(cx, g, vars)?;
    let inverse = invert(cx, g)?;
    // R_abcd = g_ae R^e_bcd, kept only where non-zero.
    let mut lower = Tensor::zeros(cx.graph, vec![n, n, n, n]);
    for (a, g_row) in g.iter().enumerate() {
        for b in 0..n {
            for c in 0..n {
                for e in 0..n {
                    let terms: Vec<NodeId> = (0..n).map(|f| mul(cx.graph, &[g_row[f], r.get(&[f, b, c, e])])).collect();
                    let value = simplified_sum(cx, &terms);
                    lower.set(&[a, b, c, e], value);
                }
            }
        }
    }
    // R^abcd = g^bf g^cg g^dh R^a_fgh: raise the last three indices in turn.
    let mut up = r;
    for slot in 1..4 {
        let mut next = Tensor::zeros(cx.graph, vec![n, n, n, n]);
        for idx in super::indices(&[n, n, n, n]) {
            let terms: Vec<NodeId> = (0..n)
                .map(|f| {
                    let mut at = idx.clone();
                    at[slot] = f;
                    mul(cx.graph, &[inverse[idx[slot]][f], up.get(&at)])
                })
                .collect();
            let value = simplified_sum(cx, &terms);
            next.set(&idx, value);
        }
        up = next;
    }
    let mut terms = Vec::new();
    for idx in super::indices(&[n, n, n, n]) {
        let (a, b, c, e) = (idx[0], idx[1], idx[2], idx[3]);
        terms.push(mul(cx.graph, &[lower.get(&[a, b, c, e]), up.get(&[a, b, c, e])]));
    }
    Some(simplified_sum(cx, &terms))
}

fn derivative_of(
    cx: &mut Cx<'_>,
    components: &[NodeId],
    vars: &[NodeId],
) -> Option<Vec<Vec<NodeId>>> {
    // out[i][j] = ∂_j components[i]
    components.iter().map(|&c| vars.iter().map(|&v| d(cx, c, v)).collect()).collect()
}

impl Extension {
    #[allow(clippy::too_many_lines)]
    fn compute(
        &self,
        cx: &mut Cx<'_>,
        args: &[NodeId],
    ) -> Option<NodeId> {
        let arg = |k: usize| args.get(k).copied();
        match self.kind {
            | Kind::Geodesic => {
                let g = square(cx.graph, arg(0)?)?;
                let vars = vector(cx.graph, arg(1)?)?;
                let v = vector(cx.graph, arg(2)?)?;
                let n = vars.len();
                if v.len() != n {
                    return None;
                }
                let gamma = christoffel_second(cx, &g, &vars)?;
                let mut out = Vec::with_capacity(n);
                for k in 0..n {
                    let mut terms = Vec::new();
                    for i in 0..n {
                        for j in 0..n {
                            terms.push(mul(cx.graph, &[gamma.get(&[k, i, j]), v[i], v[j]]));
                        }
                    }
                    let sum = simplified_sum(cx, &terms);
                    let negated = neg(cx.graph, sum);
                    out.push(cx.simplify(negated));
                }
                Some(cx.graph.node(core::LIST, &out))
            },
            | Kind::Kretschmann => {
                let g = square(cx.graph, arg(0)?)?;
                let vars = vector(cx.graph, arg(1)?)?;
                if vars.len() > 4 {
                    return None;
                }
                kretschmann(cx, &g, &vars)
            },
            | Kind::Gaussian => {
                let g = square(cx.graph, arg(0)?)?;
                let vars = vector(cx.graph, arg(1)?)?;
                if vars.len() != 2 {
                    return None;
                }
                let scalar = ricci_scalar(cx, &g, &vars)?;
                let half = cx.graph.num(crate::graph::Number::fraction(1, 2)?);
                let value = mul(cx.graph, &[half, scalar]);
                Some(normal_form(cx, value))
            },
            | Kind::Volume => {
                let g = square(cx.graph, arg(0)?)?;
                volume(cx, &g)
            },
            | Kind::Beltrami => {
                let f = arg(0)?;
                let g = square(cx.graph, arg(1)?)?;
                let vars = vector(cx.graph, arg(2)?)?;
                let inverse = invert(cx, &g)?;
                let root = volume(cx, &g)?;
                let mut gradient = Vec::with_capacity(vars.len());
                for &v in &vars {
                    gradient.push(d(cx, f, v)?);
                }
                let mut terms = Vec::new();
                for (i, &vi) in vars.iter().enumerate() {
                    let flux: Vec<NodeId> = (0..vars.len())
                        .map(|j| mul(cx.graph, &[root, inverse[i][j], gradient[j]]))
                        .collect();
                    let flux = simplified_sum(cx, &flux);
                    terms.push(d(cx, flux, vi)?);
                }
                let total = simplified_sum(cx, &terms);
                let inverse_root = {
                    let minus_one = cx.graph.int(-1);
                    cx.graph.node(core::POW, &[root, minus_one])
                };
                let value = mul(cx.graph, &[inverse_root, total]);
                Some(normal_form(cx, value))
            },
            | Kind::Divergence => {
                let field = vector(cx.graph, arg(0)?)?;
                let g = square(cx.graph, arg(1)?)?;
                let vars = vector(cx.graph, arg(2)?)?;
                if field.len() != vars.len() {
                    return None;
                }
                let root = volume(cx, &g)?;
                let mut terms = Vec::new();
                for (&component, &v) in field.iter().zip(&vars) {
                    let flux = mul(cx.graph, &[root, component]);
                    let flux = normal_form(cx, flux);
                    terms.push(d(cx, flux, v)?);
                }
                let total = simplified_sum(cx, &terms);
                let minus_one = cx.graph.int(-1);
                let inverse_root = cx.graph.node(core::POW, &[root, minus_one]);
                let value = mul(cx.graph, &[inverse_root, total]);
                Some(normal_form(cx, value))
            },
            | Kind::KillingTensor | Kind::IsKilling => {
                let xi = vector(cx.graph, arg(0)?)?;
                let g = square(cx.graph, arg(1)?)?;
                let vars = vector(cx.graph, arg(2)?)?;
                let n = vars.len();
                if xi.len() != n {
                    return None;
                }
                let dxi = derivative_of(cx, &xi, &vars)?;
                let mut rows = Vec::with_capacity(n);
                for a in 0..n {
                    let mut row = Vec::with_capacity(n);
                    for b in 0..n {
                        let mut terms = Vec::new();
                        for c in 0..n {
                            let dg = d(cx, g[a][b], vars[c])?;
                            terms.push(mul(cx.graph, &[xi[c], dg]));
                            // ∂_a ξ^c = dxi[c][a]
                            terms.push(mul(cx.graph, &[g[c][b], dxi[c][a]]));
                            terms.push(mul(cx.graph, &[g[a][c], dxi[c][b]]));
                        }
                        row.push(simplified_sum(cx, &terms));
                    }
                    rows.push(row);
                }
                if matches!(self.kind, Kind::IsKilling) {
                    let zero = rows.iter().flatten().all(|&e| cx.is_zero(e));
                    return Some(cx.graph.lit(Payload::Bool(zero)));
                }
                let items: Vec<NodeId> = rows.iter().map(|r| cx.graph.node(core::LIST, r)).collect();
                Some(cx.graph.node(core::LIST, &items))
            },
        }
    }
}

impl Kernel for Extension {
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
        true
    }
}

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    for (name, arity, kind) in [
        ("metric_geodesic_acceleration", 3, Kind::Geodesic),
        ("kretschmann", 2, Kind::Kretschmann),
        ("gaussian_curvature", 2, Kind::Gaussian),
        ("volume_element", 2, Kind::Volume),
        ("laplace_beltrami", 3, Kind::Beltrami),
        ("covariant_divergence", 3, Kind::Divergence),
        ("killing_tensor", 3, Kind::KillingTensor),
        ("is_killing", 3, Kind::IsKilling),
    ] {
        let op = i.op(OpDescriptor::new(name, Arity::Fixed(arity)).flags(OpFlags::HEAVY).cost(100))?;
        i.kernel(&format!("geometry/{name}"), Tier::Reduce, Extension { op, kind });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::rules::geometry::geometry;
    use crate::rules::testing::simplify;

    fn run(src: &str) -> String {
        simplify(&[geometry()], src)
    }

    /// With `r` positive.
    fn run_positive(src: &str) -> String {
        let (text, reduced) = crate::rules::testing::reduce_with(&[geometry()], src, &[("r", crate::graph::Facts::POSITIVE)]);
        assert!(reduced, "{text}");
        text
    }

    const SPHERE: &str = "list(list(1, 0), list(0, sin(t)^2))";

    #[test]
    fn curvature_invariants() {
        assert_eq!(run(&format!("gaussian_curvature({SPHERE}, list(t, p))")), "1");
        assert_eq!(run(&format!("kretschmann({SPHERE}, list(t, p))")), "4");
        assert_eq!(run("kretschmann(list(list(1, 0), list(0, 1)), list(x, y))"), "0");
        // A sphere of radius r has Gaussian curvature 1/r^2 and Kretschmann 4/r^4.
        assert_eq!(run_positive("gaussian_curvature(list(list(r^2, 0), list(0, r^2*sin(t)^2)), list(t, p))"), "1/r^2");
    }

    #[test]
    fn geodesics_and_operators() {
        // On the sphere a meridian is a geodesic: no acceleration for v = (1, 0).
        assert_eq!(run(&format!("metric_geodesic_acceleration({SPHERE}, list(t, p), list(1, 0))")), "list(0, 0)");
        // Moving along a latitude needs a centripetal term.
        assert_eq!(
            run(&format!("metric_geodesic_acceleration({SPHERE}, list(t, p), list(0, 1))")),
            "list(cos(t)*sin(t), 0)"
        );
        assert_eq!(run_positive("volume_element(list(list(1, 0), list(0, r^2)), list(r, t))"), "r");
        // The Laplace–Beltrami operator in polar coordinates.
        assert_eq!(run_positive("laplace_beltrami(r^2, list(list(1, 0), list(0, r^2)), list(r, t))"), "4");
        assert_eq!(run_positive("covariant_divergence(list(r, 0), list(list(1, 0), list(0, r^2)), list(r, t))"), "2");
    }

    #[test]
    fn killing_fields() {
        // Rotation is Killing for the flat plane; scaling is not.
        assert_eq!(run("is_killing(list(-y, x), list(list(1, 0), list(0, 1)), list(x, y))"), "true");
        assert_eq!(run("is_killing(list(x, y), list(list(1, 0), list(0, 1)), list(x, y))"), "false");
        assert_eq!(run("killing_tensor(list(x, y), list(list(1, 0), list(0, 1)), list(x, y))"), "list(list(2, 0), list(0, 2))");
        // d/dphi is Killing on the sphere.
        assert_eq!(run(&format!("is_killing(list(0, 1), {SPHERE}, list(t, p))")), "true");
    }
}
