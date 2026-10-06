//! Hilbert-space approximation and operator algebra: Gram matrices, best
//! approximation in a subspace, expansion coefficients in an orthogonal
//! basis, eigenvalues of degenerate integral operators and formal adjoints
//! of differential operators.
//!
//! | operator | value |
//! |---|---|
//! | `gram_matrix(list(f1, ...), x, a, b)` | the matrix of inner products `<f_i, f_j>` on `[a, b]` |
//! | `orthogonal_coefficients(f, list(e1, ...), x, a, b)` | `<e_i, f> / <e_i, e_i>`: the coefficients of `f` in an orthogonal family |
//! | `best_approximation(f, list(e1, ...), x, a, b)` | the orthogonal projection of `f` onto the span of `e_i` (normal equations `G c = <e, f>`) |
//! | `l2_distance(f, g, x, a, b)` | `‖f - g‖₂` |
//! | `sobolev_h1(f, x, a, b)` | `(‖f‖² + ‖f'‖²)^(1/2)` |
//! | `fredholm_eigenvalues(list(a1, ...), list(b1, ...), x, t, a, b)` | the non-zero eigenvalues `μ` of the integral operator with the degenerate kernel `Σ a_i(x) b_i(t)`: the eigenvalues of the matrix `∫ b_i a_j`; the characteristic values are their reciprocals |
//! | `operator_adjoint(A)` | the formal adjoint of an operator term (`op_d(x)* = -op_d(x)`, `op_mul(g)* = op_mul(conj g)`, `(AB)* = B* A*`, the Laplacian and the identity are self-adjoint); `op_int` has no adjoint in this algebra |
//! | `op_commutator(A, B)` | the operator `A B - B A` |

use super::op_algebra;
use super::vector;
use super::OpAlgebra;
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
use crate::graph::RuleError;
use crate::graph::Tier;
use crate::rules::poly::best;

#[derive(Copy, Clone)]
enum Kind {
    GramMatrix,
    Coefficients,
    BestApproximation,
    FredholmEigenvalues,
    Adjoint,
}

struct Extension {
    op: OpId,
    kind: Kind,
}

fn request(
    cx: &mut Cx<'_>,
    name: &str,
    args: &[NodeId],
) -> Option<NodeId> {
    let op = cx.graph.ops().lookup(name)?;
    cx.graph.try_node(op, args)
}

fn inner(
    cx: &mut Cx<'_>,
    f: NodeId,
    g: NodeId,
    interval: &[NodeId; 3],
) -> Option<NodeId> {
    let [x, a, b] = *interval;
    let node = request(cx, "inner_product", &[f, g, x, a, b])?;
    Some(cx.simplify(node))
}

fn adjoint(
    cx: &mut Cx<'_>,
    algebra: OpAlgebra,
    a: NodeId,
    depth: usize,
) -> Option<NodeId> {
    if depth > 32 {
        return None;
    }
    let a = best(cx.graph, a)?;
    let op = cx.graph.op(a);
    let args = cx.graph.children(a).to_vec();
    if op == algebra.mul {
        let conjugate = request(cx, "conj", &[args[0]])?;
        let conjugate = cx.simplify(conjugate);
        return Some(cx.graph.node(algebra.mul, &[conjugate]));
    }
    if op == algebra.d {
        let minus_one = cx.graph.int(-1);
        return Some(cx.graph.node(algebra.scale, &[minus_one, a]));
    }
    if op == algebra.laplacian || op == algebra.identity {
        return Some(a);
    }
    if op == algebra.add {
        let parts: Option<Vec<NodeId>> = args.iter().map(|&p| adjoint(cx, algebra, p, depth + 1)).collect();
        return Some(cx.graph.node(algebra.add, &parts?));
    }
    if op == algebra.scale {
        let conjugate = request(cx, "conj", &[args[0]])?;
        let conjugate = cx.simplify(conjugate);
        let inner = adjoint(cx, algebra, args[1], depth + 1)?;
        return Some(cx.graph.node(algebra.scale, &[conjugate, inner]));
    }
    if op == algebra.compose {
        let parts: Option<Vec<NodeId>> = args.iter().rev().map(|&p| adjoint(cx, algebra, p, depth + 1)).collect();
        return Some(cx.graph.node(algebra.compose, &parts?));
    }
    if op == algebra.power {
        let inner = adjoint(cx, algebra, args[0], depth + 1)?;
        return Some(cx.graph.node(algebra.power, &[inner, args[1]]));
    }
    None
}

impl Extension {
    fn compute(
        &self,
        cx: &mut Cx<'_>,
        args: &[NodeId],
    ) -> Option<NodeId> {
        match self.kind {
            | Kind::GramMatrix => {
                let &[fs, x, a, b] = args else { return None };
                let fs = vector(cx.graph, fs)?;
                let interval = [x, a, b];
                let mut rows = Vec::with_capacity(fs.len());
                for &fi in &fs {
                    let mut row = Vec::with_capacity(fs.len());
                    for &fj in &fs {
                        row.push(inner(cx, fi, fj, &interval)?);
                    }
                    rows.push(cx.graph.node(core::LIST, &row));
                }
                Some(cx.graph.node(core::LIST, &rows))
            },
            | Kind::Coefficients => {
                let &[f, es, x, a, b] = args else { return None };
                let es = vector(cx.graph, es)?;
                let interval = [x, a, b];
                let mut out = Vec::with_capacity(es.len());
                for &e in &es {
                    let numerator = inner(cx, e, f, &interval)?;
                    let denominator = inner(cx, e, e, &interval)?;
                    let minus_one = cx.graph.int(-1);
                    let reciprocal = cx.graph.node(core::POW, &[denominator, minus_one]);
                    let value = cx.graph.node(core::MUL, &[numerator, reciprocal]);
                    out.push(cx.simplify(value));
                }
                Some(cx.graph.node(core::LIST, &out))
            },
            | Kind::BestApproximation => {
                let &[f, es, x, a, b] = args else { return None };
                let es = vector(cx.graph, es)?;
                let interval = [x, a, b];
                let mut gram = Vec::with_capacity(es.len());
                let mut rhs = Vec::with_capacity(es.len());
                for &ei in &es {
                    let mut row = Vec::with_capacity(es.len());
                    for &ej in &es {
                        row.push(inner(cx, ei, ej, &interval)?);
                    }
                    gram.push(cx.graph.node(core::LIST, &row));
                    rhs.push(inner(cx, ei, f, &interval)?);
                }
                let gram = cx.graph.node(core::LIST, &gram);
                let rhs = cx.graph.node(core::LIST, &rhs);
                let system = request(cx, "linsolve", &[gram, rhs])?;
                let solved = cx.simplify(system);
                if cx.graph.op(solved) != core::LIST || cx.graph.children(solved).len() != es.len() {
                    return None;
                }
                let coefficients = cx.graph.children(solved).to_vec();
                let terms: Vec<NodeId> =
                    coefficients.iter().zip(&es).map(|(&c, &e)| cx.graph.node(core::MUL, &[c, e])).collect();
                let total = cx.graph.node(core::ADD, &terms);
                Some(cx.simplify(total))
            },
            | Kind::FredholmEigenvalues => {
                let &[us, vs, x, t, a, b] = args else { return None };
                let (us, vs) = (vector(cx.graph, us)?, vector(cx.graph, vs)?);
                if us.len() != vs.len() {
                    return None;
                }
                let mut rows = Vec::with_capacity(us.len());
                for &v in &vs {
                    let mut row = Vec::with_capacity(us.len());
                    for &u in &us {
                        let u_of_t = cx.graph.substitute(u, x, t);
                        let body = cx.graph.node(core::MUL, &[v, u_of_t]);
                        let integral = request(cx, "defint", &[body, t, a, b])?;
                        row.push(cx.simplify(integral));
                    }
                    rows.push(cx.graph.node(core::LIST, &row));
                }
                let matrix = cx.graph.node(core::LIST, &rows);
                let eigenvalues = request(cx, "eigenvals", &[matrix])?;
                let solved = cx.simplify(eigenvalues);
                (cx.graph.op(solved) == core::LIST).then_some(solved)
            },
            | Kind::Adjoint => {
                let algebra = op_algebra(cx.graph)?;
                adjoint(cx, algebra, *args.first()?, 0)
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
        ("gram_matrix", 4, Kind::GramMatrix),
        ("orthogonal_coefficients", 5, Kind::Coefficients),
        ("best_approximation", 5, Kind::BestApproximation),
        ("fredholm_eigenvalues", 6, Kind::FredholmEigenvalues),
        ("operator_adjoint", 1, Kind::Adjoint),
    ] {
        let op = i.op(OpDescriptor::new(name, Arity::Fixed(arity)).flags(OpFlags::HEAVY).cost(100))?;
        i.kernel(&format!("functional/{name}"), Tier::Reduce, Extension { op, kind });
    }
    i.define(&[
        "l2_distance(f, g, x, a, b) := l2_norm(f - g, x, a, b)",
        "sobolev_h1(f, x, a, b) := (l2_norm(f, x, a, b)^2 + l2_norm(diff(f, x), x, a, b)^2)^(1/2)",
        "op_commutator(A, B) := op_add(op_compose(A, B), op_scale(-1, op_compose(B, A)))",
    ])
}

#[cfg(test)]
mod tests {
    use crate::rules::functional::functional;
    use crate::rules::testing::simplify;

    fn run(src: &str) -> String {
        simplify(&[functional()], src)
    }

    #[test]
    fn gram_matrices_and_projections() {
        assert_eq!(run("gram_matrix(list(1, x), x, 0, 1)"), "list(list(1, 1/2), list(1/2, 1/3))");
        assert_eq!(run("orthogonal_coefficients(x, list(1, 2*x - 1), x, 0, 1)"), "list(1/2, 1/2)");
        // The best linear approximation of x^2 on [0, 1] is x - 1/6.
        assert_eq!(run("best_approximation(x^2, list(1, x), x, 0, 1)"), "x - 1/6");
        // sin on [0, pi] is its own best approximation in span(sin).
        assert_eq!(run("best_approximation(sin(x), list(sin(x)), x, 0, pi)"), "sin(x)");
        assert_eq!(run("l2_distance(x, 0, x, 0, 1)"), "1/3^(1/2)");
    }

    #[test]
    fn degenerate_kernels_and_adjoints() {
        assert_eq!(run("fredholm_eigenvalues(list(x), list(t), x, t, 0, 1)"), "list(1/3)");
        assert_eq!(run("fredholm_eigenvalues(list(sin(x)), list(sin(t)), x, t, 0, pi)"), "list(1/2*pi)");
        assert_eq!(run("operator_adjoint(op_d(x))"), "op_scale(-1, op_d(x))");
        assert_eq!(run("operator_adjoint(op_compose(op_d(x), op_mul(2)))"), "op_compose(op_mul(2), op_scale(-1, op_d(x)))");
        assert_eq!(run("operator_adjoint(op_laplacian(list(x, y)))"), "op_laplacian(list(x, y))");
        // [d/dx, x] f = f.
        assert_eq!(run("op_apply(op_commutator(op_d(x), op_mul(x)), f(x))"), "f(x)");
    }
}
