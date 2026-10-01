//! Functional analysis and integral equations.
//!
//! # Linear operators
//!
//! An operator on functions is a term built from `op_mul(g)`
//! (multiplication by `g`), `op_d(x)` (`d/dx`), `op_int(a, x)`
//! (`f ↦ ∫_a^x f`), `op_laplacian(list(x, ...))`, `op_identity`,
//! `op_add(A, B, ...)`, `op_scale(c, A)`, `op_compose(A, B, ...)` (`B`
//! acts first) and `op_power(A, n)`; a matrix acts on a vector by
//! `matmul`. `op_apply(A, f)` applies it. The quantum-mechanical
//! operators of [`physics`](super::physics) are built from these.
//!
//! # Function spaces
//!
//! `L²` and `L^p` on an interval `[a, b]` of the variable `x`; functions
//! may be complex (the inner product conjugates its first argument).
//!
//! | operator | value |
//! |---|---|
//! | `inner_product(f, g, x, a, b)` | `∫_a^b conj(f) g dx` |
//! | `l2_norm(f, x, a, b)`, `lp_norm(f, p, x, a, b)` | `‖f‖₂`, `(∫ |f|^p)^(1/p)` |
//! | `are_orthogonal(f, g, x, a, b)` | `true` when the inner product vanishes |
//! | `project_onto(f, g, x, a, b)` | the component of `f` along `g` |
//! | `gram_schmidt(list(f1, ...), x, a, b)`, `gram_schmidt_orthonormal(...)` | an orthogonal (orthonormal) basis of the span, dependent members dropped |
//!
//! # Integral equations
//!
//! For `y(x) = f(x) + λ ∫ K(x, t) y(t) dt` (Fredholm on `[a, b]`, Volterra
//! on `[a, x]`):
//!
//! | operator | value |
//! |---|---|
//! | `fredholm_neumann(f, lambda, K, x, t, a, b, n)` | the `n`-th iterate of the Neumann series |
//! | `fredholm_separable(f, lambda, list(a1(x), ...), list(b1(t), ...), x, t, a, b)` | the exact solution for the degenerate kernel `K = Σ a_i(x) b_i(t)` |
//! | `fredholm_solve(f, lambda, K, x, t, a, b)` | exact, when `K` separates into products (found by expanding it) |
//! | `volterra_successive(f, lambda, K, x, t, a, n)` | the `n`-th Picard iterate |
//! | `volterra_to_ode(f, lambda, K, x, t, a, y(x))` | the equivalent differential equation `y' = f' + λ K(x, x) y + λ ∫ ∂K/∂x y dt` (second kind, differentiated once) |
//! | `volterra_solve(f, lambda, K, x, t, a)` | exact for kernels `K = k(x - t)` with a Laplace-transformable `k`, by the convolution theorem; for kernels free of `t` and `x`, by the ODE |
//! | `airfoil_equation(f, x, t)` | the inversion of the finite Hilbert transform, `y(x) = -1/(π sqrt(1 - x²)) ∫_{-1}^{1} sqrt(1 - t²) f(t)/(t - x) dt + C/sqrt(1 - x²)` |

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
use crate::graph::RuleError;
use crate::graph::RuleSet;
use crate::graph::Tier;
use crate::rules::complex::build::add;
use crate::rules::complex::build::mul;
use crate::rules::complex::build::neg;
use crate::rules::complex::build::powi;
use crate::rules::complex::build::sub;
use crate::rules::poly::best;
use crate::rules::poly::repr::from_term;
use crate::rules::poly::repr::Gens;
use crate::rules::poly::repr::Limits;

use super::calculus::derivative;
use super::complex::complex;
use super::linalg::linalg;
use super::logic::logic;
use super::solve::solve_linear;
use super::transforms::transforms;

/// The functional-analysis rule set.
#[must_use]
pub fn functional() -> RuleSet {
    RuleSet::new("functional", install).needs(linalg()).needs(complex()).needs(logic()).needs(transforms())
}

/// The operator-algebra constructors of one graph.
#[derive(Copy, Clone, Debug)]
pub(crate) struct OpAlgebra {
    mul: OpId,
    d: OpId,
    int: OpId,
    laplacian: OpId,
    identity: OpId,
    add: OpId,
    scale: OpId,
    compose: OpId,
    power: OpId,
    matmul: OpId,
    defint: OpId,
}

/// The operator algebra of `graph`, once this rule set is installed.
pub(crate) fn op_algebra(graph: &Graph) -> Option<OpAlgebra> {
    let ops = graph.ops();
    Some(OpAlgebra {
        mul: ops.lookup("op_mul")?,
        d: ops.lookup("op_d")?,
        int: ops.lookup("op_int")?,
        laplacian: ops.lookup("op_laplacian")?,
        identity: ops.lookup("op_identity")?,
        add: ops.lookup("op_add")?,
        scale: ops.lookup("op_scale")?,
        compose: ops.lookup("op_compose")?,
        power: ops.lookup("op_power")?,
        matmul: ops.lookup("matmul")?,
        defint: ops.lookup("defint")?,
    })
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum Request {
    Apply,
    Orthogonal,
    GramSchmidt { normalise: bool },
    NeumannSeries,
    Separable,
    FredholmSolve,
    Successive,
    VolterraToOde,
    VolterraSolve,
}

fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let plain = |name: &str, arity: Arity| OpDescriptor::new(name, arity);
    for (name, arity) in [
        ("op_mul", Arity::Fixed(1)),
        ("op_d", Arity::Fixed(1)),
        ("op_int", Arity::Fixed(2)),
        ("op_laplacian", Arity::Fixed(1)),
        ("op_identity", Arity::Fixed(0)),
        ("op_add", Arity::Variadic),
        ("op_scale", Arity::Fixed(2)),
        ("op_compose", Arity::Variadic),
        ("op_power", Arity::Fixed(2)),
    ] {
        i.op(plain(name, arity))?;
    }
    let algebra = op_algebra(i.graph()).ok_or(RuleError::Invalid {
        rule: "functional".to_owned(),
        reason: "needs the linear algebra and calculus rule sets",
    })?;
    let heavy = |name: &str, arity: u8| OpDescriptor::new(name, Arity::Fixed(arity)).flags(OpFlags::HEAVY).cost(100);
    for (name, arity, request) in [
        ("op_apply", 2, Request::Apply),
        ("are_orthogonal", 5, Request::Orthogonal),
        ("gram_schmidt", 4, Request::GramSchmidt { normalise: false }),
        ("gram_schmidt_orthonormal", 4, Request::GramSchmidt { normalise: true }),
        ("fredholm_neumann", 8, Request::NeumannSeries),
        ("fredholm_separable", 8, Request::Separable),
        ("fredholm_solve", 7, Request::FredholmSolve),
        ("volterra_successive", 7, Request::Successive),
        ("volterra_to_ode", 7, Request::VolterraToOde),
        ("volterra_solve", 6, Request::VolterraSolve),
    ] {
        let op = i.op(heavy(name, arity))?;
        i.kernel(&format!("functional/{name}"), Tier::Reduce, Functional { op, request, algebra });
    }
    i.define(&[
        "inner_product(f, g, x, a, b) := defint(conj(f) * g, x, a, b)",
        "l2_norm(f, x, a, b) := inner_product(f, f, x, a, b)^(1/2)",
        "lp_norm(f, p, x, a, b) := defint(abs(f)^p, x, a, b)^(1/p)",
        "project_onto(f, g, x, a, b) := inner_product(g, f, x, a, b) / inner_product(g, g, x, a, b) * g",
        "airfoil_equation(f, x, t) := -1/(pi*(1 - x^2)^(1/2)) * defint((1 - t^2)^(1/2) * f/(t - x), t, -1, 1) + C/(1 - x^2)^(1/2)",
    ])
}

struct Functional {
    op: OpId,
    request: Request,
    algebra: OpAlgebra,
}

impl Kernel for Functional {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let args = cx.graph.children(node).to_vec();
        let result = match self.request {
            | Request::Apply => apply(cx, self.algebra, args[0], args[1], 0),
            | Request::Orthogonal => orthogonal(cx, self.algebra, &args),
            | Request::GramSchmidt { normalise } => gram_schmidt(cx, self.algebra, &args, normalise),
            | Request::NeumannSeries => neumann(cx, self.algebra, &args),
            | Request::Separable => separable(cx, self.algebra, &args),
            | Request::FredholmSolve => fredholm_solve(cx, self.algebra, &args),
            | Request::Successive => successive(cx, self.algebra, &args),
            | Request::VolterraToOde => volterra_to_ode(cx, self.algebra, &args),
            | Request::VolterraSolve => volterra_solve(cx, self.algebra, &args),
        };
        result.map_or(Outcome::Pass, Outcome::Equal)
    }
}

// ----------------------------------------------------------------------
// Operators
// ----------------------------------------------------------------------

/// Whether a heavy operator occurs anywhere in `node`.
fn is_heavy(
    graph: &Graph,
    node: NodeId,
) -> bool {
    let mut stack = vec![node];
    let mut seen = std::collections::HashSet::new();
    while let Some(n) = stack.pop() {
        if !seen.insert(n) {
            continue;
        }
        if graph.ops().get(graph.op(n)).flags.has(OpFlags::HEAVY) {
            return true;
        }
        stack.extend_from_slice(graph.children(n));
    }
    false
}

fn vector(
    graph: &mut Graph,
    node: NodeId,
) -> Option<Vec<NodeId>> {
    let term = best(graph, node)?;
    (graph.op(term) == core::LIST).then(|| graph.children(term).to_vec())
}

/// `A f` for an operator term `A`.
pub(crate) fn apply(
    cx: &mut Cx<'_>,
    algebra: OpAlgebra,
    a: NodeId,
    f: NodeId,
    depth: usize,
) -> Option<NodeId> {
    if depth > 32 {
        return None;
    }
    let a = best(cx.graph, a)?;
    let f = best(cx.graph, f)?;
    let op = cx.graph.op(a);
    let args = cx.graph.children(a).to_vec();
    let graph = &mut *cx.graph;
    if op == algebra.mul {
        return Some(mul(graph, &[args[0], f]));
    }
    if op == algebra.d {
        graph.symbol_of(args[0])?;
        return derivative(graph, f, args[0]);
    }
    if op == algebra.int {
        // f ↦ ∫_a^x f(t) dt
        let (lower, x) = (args[0], args[1]);
        graph.symbol_of(x)?;
        let fresh = graph.interner_mut().fresh_symbol("t");
        let t = graph.symbol_node(fresh);
        let body = graph.substitute(f, x, t);
        return Some(graph.node(algebra.defint, &[body, t, lower, x]));
    }
    if op == algebra.identity {
        return Some(f);
    }
    if op == algebra.laplacian {
        let vars = vector(graph, args[0])?;
        let mut terms = Vec::with_capacity(vars.len());
        for x in vars {
            graph.symbol_of(x)?;
            let first = derivative(graph, f, x)?;
            terms.push(derivative(graph, first, x)?);
        }
        return Some(add(graph, &terms));
    }
    if op == algebra.add {
        let mut terms = Vec::with_capacity(args.len());
        for part in args {
            terms.push(apply(cx, algebra, part, f, depth + 1)?);
        }
        return Some(add(cx.graph, &terms));
    }
    if op == algebra.scale {
        let inner = apply(cx, algebra, args[1], f, depth + 1)?;
        return Some(mul(cx.graph, &[args[0], inner]));
    }
    if op == algebra.compose {
        let mut state = f;
        for &part in args.iter().rev() {
            state = apply(cx, algebra, part, state, depth + 1)?;
        }
        return Some(state);
    }
    if op == algebra.power {
        let n = graph.number_of(args[1]).and_then(Number::to_i64).filter(|n| (0..=16).contains(n))?;
        let mut state = f;
        for _ in 0..n {
            state = apply(cx, algebra, args[0], state, depth + 1)?;
        }
        return Some(state);
    }
    if op == core::LIST {
        return Some(graph.node(algebra.matmul, &[a, f]));
    }
    // Any other term multiplies, once it is a closed form (an operator
    // definition that has not been expanded yet is not a multiplier).
    if is_heavy(graph, a) {
        return None;
    }
    Some(mul(graph, &[a, f]))
}

// ----------------------------------------------------------------------
// Function spaces
// ----------------------------------------------------------------------

/// `∫_a^b conj(f) g dx`, simplified (computed by a nested run).
fn inner(
    cx: &mut Cx<'_>,
    algebra: OpAlgebra,
    f: NodeId,
    g: NodeId,
    x: NodeId,
    a: NodeId,
    b: NodeId,
) -> Option<NodeId> {
    let conj = cx.graph.ops().lookup("conj")?;
    let cf = cx.graph.node(conj, &[f]);
    let product = mul(cx.graph, &[cf, g]);
    let integral = cx.graph.node(algebra.defint, &[product, x, a, b]);
    let value = cx.simplify(integral);
    (!is_heavy(cx.graph, value)).then_some(value)
}

fn truth(
    graph: &mut Graph,
    value: bool,
) -> Option<NodeId> {
    let op = graph.ops().lookup(if value { "true" } else { "false" })?;
    Some(graph.node(op, &[]))
}

fn orthogonal(
    cx: &mut Cx<'_>,
    algebra: OpAlgebra,
    args: &[NodeId],
) -> Option<NodeId> {
    let &[f, g, x, a, b] = args else {
        return None;
    };
    let product = inner(cx, algebra, f, g, x, a, b)?;
    let zero = cx.is_zero(product);
    if !zero && cx.graph.number_of(product).is_none() && !cx.graph.free_symbols(cx.graph.find(product)).is_empty() {
        // A parameter-dependent product that does not simplify to zero.
        return truth(cx.graph, false);
    }
    truth(cx.graph, zero)
}

fn gram_schmidt(
    cx: &mut Cx<'_>,
    algebra: OpAlgebra,
    args: &[NodeId],
    normalise: bool,
) -> Option<NodeId> {
    let &[functions, x, a, b] = args else {
        return None;
    };
    let functions = vector(cx.graph, functions)?;
    let mut basis: Vec<(NodeId, NodeId)> = Vec::new(); // (u, <u, u>)
    for f in functions {
        let mut u = f;
        for &(e, norm2) in &basis {
            let c = inner(cx, algebra, e, f, x, a, b)?;
            let norm2_inv = powi(cx.graph, norm2, -1);
            let ratio = mul(cx.graph, &[c, norm2_inv]);
            let component = mul(cx.graph, &[ratio, e]);
            u = sub(cx.graph, u, component);
        }
        let u = cx.simplify(u);
        let norm2 = inner(cx, algebra, u, u, x, a, b)?;
        if cx.is_zero(norm2) {
            continue;
        }
        basis.push((u, norm2));
    }
    let half = cx.graph.num(Number::fraction(-1, 2)?);
    let out: Vec<NodeId> = basis
        .into_iter()
        .map(|(u, norm2)| {
            if normalise {
                let scale = cx.graph.node(core::POW, &[norm2, half]);
                let value = mul(cx.graph, &[scale, u]);
                cx.simplify(value)
            } else {
                u
            }
        })
        .collect();
    Some(cx.graph.node(core::LIST, &out))
}

// ----------------------------------------------------------------------
// Integral equations
// ----------------------------------------------------------------------

/// `f + λ ∫_lower^upper K(x, t) y(t) dt` for a concrete `y(x)`.
#[allow(clippy::too_many_arguments)]
fn iterate(
    cx: &mut Cx<'_>,
    algebra: OpAlgebra,
    f: NodeId,
    lambda: NodeId,
    kernel: NodeId,
    x: NodeId,
    t: NodeId,
    lower: NodeId,
    upper: NodeId,
    y: NodeId,
) -> Option<NodeId> {
    let y_t = cx.graph.substitute(y, x, t);
    let integrand = mul(cx.graph, &[kernel, y_t]);
    let integral = cx.graph.node(algebra.defint, &[integrand, t, lower, upper]);
    let scaled = mul(cx.graph, &[lambda, integral]);
    let next = add(cx.graph, &[f, scaled]);
    let value = cx.simplify(next);
    (!is_heavy(cx.graph, value)).then_some(value)
}

fn neumann(
    cx: &mut Cx<'_>,
    algebra: OpAlgebra,
    args: &[NodeId],
) -> Option<NodeId> {
    let &[f, lambda, kernel, x, t, a, b, n] = args else {
        return None;
    };
    let n = cx.graph.number_of(n).and_then(Number::to_i64).filter(|n| (0..=12).contains(n))?;
    let mut y = f;
    for _ in 0..n {
        y = iterate(cx, algebra, f, lambda, kernel, x, t, a, b, y)?;
    }
    Some(y)
}

fn successive(
    cx: &mut Cx<'_>,
    algebra: OpAlgebra,
    args: &[NodeId],
) -> Option<NodeId> {
    let &[f, lambda, kernel, x, t, a, n] = args else {
        return None;
    };
    let n = cx.graph.number_of(n).and_then(Number::to_i64).filter(|n| (0..=12).contains(n))?;
    let mut y = f;
    for _ in 0..n {
        y = iterate(cx, algebra, f, lambda, kernel, x, t, a, x, y)?;
    }
    Some(y)
}

/// Degenerate kernel `Σ a_i(x) b_i(t)`: with `c_k = ∫ b_k(t) y(t) dt`,
/// `c_k - λ Σ_i c_i ∫ b_k a_i = ∫ b_k f`, a linear system.
#[allow(clippy::too_many_arguments)]
fn solve_degenerate(
    cx: &mut Cx<'_>,
    algebra: OpAlgebra,
    f: NodeId,
    lambda: NodeId,
    a_funcs: &[NodeId],
    b_funcs: &[NodeId],
    x: NodeId,
    t: NodeId,
    lower: NodeId,
    upper: NodeId,
) -> Option<NodeId> {
    let m = a_funcs.len();
    if m != b_funcs.len() || m == 0 {
        return None;
    }
    let unknowns: Vec<NodeId> = (0..m)
        .map(|k| {
            let s = cx.graph.interner_mut().fresh_symbol(&format!("c{k}"));
            cx.graph.symbol_node(s)
        })
        .collect();
    let integrate = |cx: &mut Cx<'_>, body: NodeId| -> Option<NodeId> {
        let integral = cx.graph.node(algebra.defint, &[body, t, lower, upper]);
        let value = cx.simplify(integral);
        (!is_heavy(cx.graph, value)).then_some(value)
    };
    let mut equations = Vec::with_capacity(m);
    for k in 0..m {
        let b_k = cx.graph.substitute(b_funcs[k], x, t);
        let f_t = cx.graph.substitute(f, x, t);
        let beta_body = mul(cx.graph, &[b_k, f_t]);
        let beta = integrate(cx, beta_body)?;
        let mut terms = vec![unknowns[k]];
        for i in 0..m {
            let a_i = cx.graph.substitute(a_funcs[i], x, t);
            let body = mul(cx.graph, &[b_k, a_i]);
            let alpha = integrate(cx, body)?;
            let product = mul(cx.graph, &[lambda, alpha, unknowns[i]]);
            terms.push(neg(cx.graph, product));
        }
        terms.push(neg(cx.graph, beta));
        let equation = add(cx.graph, &terms);
        equations.push(cx.simplify(equation));
    }
    let values = solve_linear(cx.graph, &equations, &unknowns)?;
    let mut sum = Vec::with_capacity(m);
    for (value, &a_i) in values.into_iter().zip(a_funcs) {
        sum.push(mul(cx.graph, &[value, a_i]));
    }
    let combination = add(cx.graph, &sum);
    let scaled = mul(cx.graph, &[lambda, combination]);
    let solution = add(cx.graph, &[f, scaled]);
    Some(cx.simplify(solution))
}

fn separable(
    cx: &mut Cx<'_>,
    algebra: OpAlgebra,
    args: &[NodeId],
) -> Option<NodeId> {
    let &[f, lambda, a_list, b_list, x, t, a, b] = args else {
        return None;
    };
    let a_funcs = vector(cx.graph, a_list)?;
    let b_funcs = vector(cx.graph, b_list)?;
    solve_degenerate(cx, algebra, f, lambda, &a_funcs, &b_funcs, x, t, a, b)
}

/// Splits `K(x, t)` into `Σ a_i(x) b_i(t)` by expanding it: every term of
/// the expansion is a product of a factor in `x` and a factor in `t`.
fn split_kernel(
    graph: &mut Graph,
    kernel: NodeId,
    x: NodeId,
    t: NodeId,
) -> Option<(Vec<NodeId>, Vec<NodeId>)> {
    let (xs, ts) = (graph.symbol_of(x)?, graph.symbol_of(t)?);
    let kernel = best(graph, kernel)?;
    let mut gens = Gens::default();
    let poly = from_term(graph, &mut gens, kernel, Limits { terms: 64, exponent: 16 })?;
    let (mut a_funcs, mut b_funcs) = (Vec::new(), Vec::new());
    for (mono, coeff) in poly.terms() {
        let (mut in_x, mut in_t) = (vec![graph.num(coeff.clone())], Vec::new());
        for &(g, e) in mono {
            let node = gens.node(g)?;
            let factor = powi(graph, node, i64::from(e));
            let (dx, dt) = (graph.depends_on(graph.find(node), xs), graph.depends_on(graph.find(node), ts));
            match (dx, dt) {
                | (true, true) => return None,
                | (false, true) => in_t.push(factor),
                | _ => in_x.push(factor),
            }
        }
        a_funcs.push(mul(graph, &in_x));
        b_funcs.push(mul(graph, &in_t));
    }
    Some((a_funcs, b_funcs))
}

fn fredholm_solve(
    cx: &mut Cx<'_>,
    algebra: OpAlgebra,
    args: &[NodeId],
) -> Option<NodeId> {
    let &[f, lambda, kernel, x, t, a, b] = args else {
        return None;
    };
    let (a_funcs, b_funcs) = split_kernel(cx.graph, kernel, x, t)?;
    solve_degenerate(cx, algebra, f, lambda, &a_funcs, &b_funcs, x, t, a, b)
}

/// Differentiating `y = f + λ ∫_a^x K y dt` once.
fn volterra_to_ode(
    cx: &mut Cx<'_>,
    algebra: OpAlgebra,
    args: &[NodeId],
) -> Option<NodeId> {
    let &[f, lambda, kernel, x, t, a, y] = args else {
        return None;
    };
    let y = best(cx.graph, y)?;
    let y_prime = derivative(cx.graph, y, x)?;
    let f_prime = derivative(cx.graph, f, x)?;
    let k_xx = cx.graph.substitute(kernel, t, x);
    let direct = mul(cx.graph, &[lambda, k_xx, y]);
    let dk = derivative(cx.graph, kernel, x)?;
    let dk = cx.simplify(dk);
    let mut terms = vec![f_prime, direct];
    if !cx.graph.number_of(dk).is_some_and(Number::is_zero) {
        let y_t = cx.graph.replace_subterm(y, x, t);
        let body = mul(cx.graph, &[dk, y_t]);
        let integral = cx.graph.node(algebra.defint, &[body, t, a, x]);
        terms.push(mul(cx.graph, &[lambda, integral]));
    }
    let rhs = add(cx.graph, &terms);
    let rhs = cx.simplify(rhs);
    Some(cx.graph.node(core::EQ, &[y_prime, rhs]))
}

/// Volterra equations of convolution type, `K = k(x - t)`, on `[0, x]`:
/// `Y = F / (1 - λ k̂)` in Laplace space.
fn volterra_solve(
    cx: &mut Cx<'_>,
    _algebra: OpAlgebra,
    args: &[NodeId],
) -> Option<NodeId> {
    let &[f, lambda, kernel, x, t, a] = args else {
        return None;
    };
    if !cx.graph.number_of(a).is_some_and(Number::is_zero) {
        return None;
    }
    // k(u) with u = x - t: substitute t = x - u and check x drops out.
    let u_symbol = cx.graph.interner_mut().fresh_symbol("u");
    let u = cx.graph.symbol_node(u_symbol);
    let shifted = sub(cx.graph, x, u);
    let k_u = cx.graph.substitute(kernel, t, shifted);
    let k_u = cx.simplify(k_u);
    let x_symbol = cx.graph.symbol_of(x)?;
    if cx.graph.depends_on(cx.graph.find(k_u), x_symbol) {
        return None;
    }
    let (laplace, inverse) = (cx.graph.ops().lookup("laplace")?, cx.graph.ops().lookup("inverse_laplace")?);
    let s_symbol = cx.graph.interner_mut().fresh_symbol("s");
    let s = cx.graph.symbol_node(s_symbol);
    let k_hat = cx.graph.node(laplace, &[k_u, u, s]);
    let f_hat = cx.graph.node(laplace, &[f, x, s]);
    let lk = mul(cx.graph, &[lambda, k_hat]);
    let one = cx.graph.int(1);
    let denominator = sub(cx.graph, one, lk);
    let denominator_inv = powi(cx.graph, denominator, -1);
    let y_hat = mul(cx.graph, &[f_hat, denominator_inv]);
    let y_hat = cx.simplify(y_hat);
    let request = cx.graph.node(inverse, &[y_hat, s, x]);
    let solution = cx.simplify(request);
    (!is_heavy(cx.graph, solution)).then_some(solution)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Facts;
    use crate::rules::testing::reduce_with;

    fn run(src: &str) -> String {
        let (text, reduced) = reduce_with(&[functional()], src, &[("x", Facts::REAL), ("t", Facts::REAL)]);
        assert!(reduced, "`{src}` was not fully reduced: {text}");
        text
    }

    #[test]
    fn operators() {
        assert_eq!(run("op_apply(op_compose(op_d(x), op_mul(x)), x^2)"), "3*x^2");
        assert_eq!(run("op_apply(op_int(0, x), cos(x))"), "sin(x)");
        assert_eq!(run("op_apply(op_add(op_identity, op_scale(2, op_d(x))), exp(x))"), "3*exp(x)");
        assert_eq!(run("op_apply(op_power(op_d(x), 3), x^4)"), "24*x");
    }

    #[test]
    fn function_spaces() {
        assert_eq!(run("inner_product(sin(x), sin(x), x, 0, pi)"), "1/2*pi");
        assert_eq!(run("are_orthogonal(sin(x), cos(x), x, -pi, pi)"), "true");
        assert_eq!(run("are_orthogonal(x, x^2 + 1, x, 0, 1)"), "false");
        assert_eq!(run("l2_norm(1, x, 0, 4)"), "2");
        assert_eq!(run("lp_norm(x, 1, x, 0, 2)"), "2");
        assert_eq!(run("project_onto(x, 1, x, 0, 1)"), "1/2");
        // Gram–Schmidt on 1, x, x² over [-1, 1]: Legendre polynomials.
        assert_eq!(run("gram_schmidt(list(1, x, x^2), x, -1, 1)"), run("list(1, x, x^2 - 1/3)"));
        let orthonormal = run("gram_schmidt_orthonormal(list(1, x), x, -1, 1)");
        let difference = run(&format!("{orthonormal} - list(2^(-1/2), x*(3/2)^(1/2))"));
        assert!(difference == "0" || difference == "list(0, 0)", "{orthonormal}");
        // A dependent member is dropped.
        assert_eq!(run("gram_schmidt(list(x, 2*x), x, 0, 1)"), "list(x)");
    }

    #[test]
    fn integral_equations() {
        // y = x + ∫_0^1 x t y dt  has the solution y = 3x/2.
        assert_eq!(run("fredholm_separable(x, 1, list(x), list(t), x, t, 0, 1)"), "3/2*x");
        assert_eq!(run("fredholm_solve(x, 1, x*t, x, t, 0, 1)"), "3/2*x");
        // Two-term kernel: y = 1 + λ ∫_0^1 (x + t) y dt with λ = 1.
        let two = run("fredholm_solve(1, 1, x + t, x, t, 0, 1)");
        assert_eq!(two, run("-12*x - 6"));
        // Neumann iterates converge to the same answer for a small λ.
        assert_eq!(run("fredholm_neumann(x, 0, x*t, x, t, 0, 1, 3)"), "x");
        assert_eq!(run("fredholm_neumann(x, 1, x*t, x, t, 0, 1, 2)"), run("x + 1/3*x + 1/9*x"));
        // y = 1 + ∫_0^x y dt  is solved by exp(x).
        assert_eq!(run("volterra_successive(1, 1, 1, x, t, 0, 3)"), run("1 + x + x^2/2 + x^3/6"));
        assert_eq!(run("volterra_solve(1, 1, 1, x, t, 0)"), "exp(x)");
        // y = sin(x) - ∫_0^x (x - t) y dt  (Laplace: Y = s/(s² + 1)² · ...)
        let conv = run("volterra_solve(x, -1, x - t, x, t, 0)");
        assert_eq!(conv, "sin(x)");
        assert_eq!(run("volterra_to_ode(1, 1, 1, x, t, 0, y(x))"), "diff(y(x), x) = y(x)");
    }
}
