//! The Meijer G-function: operator, exact reductions, differentiation,
//! integration and series.
//!
//! The operator is `meijerg(m, n, list(a_1, ..., a_p), list(b_1, ..., b_q), z)`
//! with `m`, `n` literal integers. Reductions to elementary and special
//! functions, the derivative and antiderivative formulas, the series at zero
//! and asymptotics at infinity, and the numeric value are installed here.
//!
//! # Operators
//!
//! | operator | meaning |
//! |---|---|
//! | `meijerg(m, n, list(a), list(b), z)` | `G^{m,n}_{p,q}(z | a; b)`, `p = len(a)`, `q = len(b)`; the first `n` entries of `a` and the first `m` entries of `b` are the contour-defining sets |
//! | `meijer_mellin(m, n, list(a), list(b), s)` | the Mellin transform `int_0^inf x^{s-1} G(x) dx`, as the exact gamma product |
//! | `meijer_series(m, n, list(a), list(b), z, k[, direction])` | the first `k` terms of the residue series at `z = 0` (`direction = 0`, the default) or of the expansion at infinity (`direction = oo`) |
//!
//! `meijer_mellin` is the exact formula
//! ```text
//! int_0^inf x^{s-1} G^{m,n}_{p,q}(x | a; b) dx
//!   = prod_{j=1}^m Gamma(b_j+s) prod_{j=1}^n Gamma(1-a_j-s)
//!     / (prod_{j=m+1}^q Gamma(1-b_j-s) prod_{j=n+1}^p Gamma(a_j+s)),
//! ```
//! valid in a vertical strip of `s` (between the pole families of the two
//! numerator products); the operator returns the formula, not a promise that
//! the strip is non-empty for the given parameters.
//!
//! # Exact reductions
//!
//! Every reduction is a kernel, not a text rewrite, so list arithmetic is
//! exact. Before a reduction is returned it is checked numerically against
//! the numeric kernel at several arguments and parameter values; a reduction
//! whose check fails is skipped (the request stays unevaluated rather than
//! becoming wrong). The table:
//!
//! | parameters | value |
//! |---|---|
//! | `G^{1,0}_{0,1}(z | b)` | `z^b e^{-z}` |
//! | `G^{1,1}_{1,1}(z | a; b)` | `Gamma(1+b-a) z^b (1+z)^{a-b-1}` |
//! | `G^{1,0}_{0,2}(z | b, b+1/2)` | `z^b cos(2 sqrt z)/sqrt(pi)` |
//! | `G^{1,0}_{0,2}(z | b+1/2, b)` | `z^b sin(2 sqrt z)/sqrt(pi)` |
//! | `G^{1,1}_{1,2}(z | 1; 1/2, 0)` | `sqrt(pi) erf(sqrt z)` |
//! | `G^{2,0}_{1,2}(z | 1; 1/2, 0)` | `sqrt(pi) erfc(sqrt z)` |
//! | `G^{1,0}_{0,2}(z | nu/2, -nu/2)` | `J_nu(2 sqrt z)` |
//! | `G^{1,0}_{0,2}(-w | nu/2, -nu/2)` | `i^nu I_nu(2 sqrt w)` |
//! | `G^{2,0}_{0,2}(z | nu/2, -nu/2)` | `2 K_nu(2 sqrt z)` |
//! | `G^{1,0}_{0,2}(z | b, b)` | `z^b J_0(2 sqrt z)` |
//! | `G^{2,0}_{1,3}(z | (-nu-1)/2; nu/2, -nu/2, (-nu-1)/2)` | `Y_nu(2 sqrt z)` |
//! | `G^{1,2}_{2,2}(z | 1, 1; 1, 0)` | `ln(1+z)` |
//! | `G^{1,2}_{2,2}(z | 1/2, 1; 1/2, 0)` | `2 atan(sqrt z)` |
//! | any `p < q` with a terminating numerator parameter | the finite Slater `pF_{q-1}` series as an explicit polynomial |
//!
//! The general terminating case builds the full Slater series
//! `sum_h z^{b_h} C_h pF_{q-1}(A_h; B_h; (-1)^{p-m-n} z)` and truncates it
//! exactly where a numerator parameter `A_i = 1 + b_h - a_i` is a
//! non-positive integer. For `p = q` with `n = 0` the G-function is `0` for
//! `|z| > 1` (the defining contour encloses no pole there), so its finite
//! Slater series is *not* a global identity; such reductions are
//! deliberately not installed. `p = q` reductions with `n > 0` are
//! installed, because their value is the analytic continuation of the same
//! closed form.
//!
//! # Differentiation and integration
//!
//! The kernel on `diff` nodes uses
//! `z d/dz G^{m,n}_{p,q}(z|a;b) = G^{m,n+1}_{p+1,q+1}(z | 0, a; b, 1)`
//! (DLMF 16.18.4). The antiderivative follows from the same formula with
//! `h = -1`:
//! `int G^{m,n}_{p,q}(z|a;b) dz = z G^{m,n+1}_{p+1,q+1}(z | 0, a; b, -1)`.
//! A linear argument `c x` is divided by `c`. Both kernels verify themselves
//! against finite differences when the parameters are numeric.
//! `Partials` cannot express the parameter shift (it has to prepend to a
//! list), so these kernels, not the attribute, are the source of truth.
//!
//! # Numeric phase
//!
//! A kernel on `meijerg` returns `Outcome::Approx` when a numeric answer is
//! wanted, every free symbol is bound and all parameters are real numbers:
//! the value comes from [`crate::kernels::meijer::meijer_g`]. Otherwise it
//! returns `Pass`, so a symbolic request is never replaced by a float.
//! The symbolic reduction kernel is inactive in the numeric phase, where
//! the direct evaluation is both faster and free of branch surprises.
//!
//! # Truncation
//!
//! `meijer_series` returns a finite sum and, when the `big_o` operator is
//! installed and the leading exponents are zero, appends `big_o(z^k)` (at
//! zero) or `big_o(z^(-k))` (at infinity). With non-zero exponents the
//! remainder is `O(z^{b_min+k})` (respectively `O(z^{a_min-1-k})`) and is
//! documented rather than written, because `big_o` would misstate the
//! power. Numeric parameters are checked against the kernel at a small (or,
//! at infinity, large) argument before the series is returned. The
//! `direction = oo` form exists only where the `1/z` residue series applies
//! (`p > q`, or `p = q` with `|z| > 1`); for the entire cases (`p < q`) the
//! `1/z` series is divergent, the check fails, and the request stays
//! unevaluated instead of receiving an asymptotic answer.

use std::collections::HashMap;

use num_complex::Complex64;

use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::Env;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpDescriptor;
use crate::graph::OpFlags;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::RuleError;
use crate::graph::SymbolId;
use crate::graph::Tier;
use crate::graph::op::core;
use crate::graph::rule::Installer;
use crate::kernels::meijer as numeric;
use crate::rules::poly::repr::Gens;
use crate::rules::poly::repr::Limits;
use crate::rules::poly::repr::Poly;
use crate::rules::poly::repr::from_term;
use crate::rules::poly::repr::to_term;

type C = Complex64;

/// The operators the kernels build their answers from.
#[derive(Copy, Clone)]
struct Ops {
    meijerg: OpId,
    exp: OpId,
    ln: OpId,
    gamma: Option<OpId>,
    sin: Option<OpId>,
    cos: Option<OpId>,
    atan: Option<OpId>,
    erf: Option<OpId>,
    erfc: Option<OpId>,
    besselj: Option<OpId>,
    bessely: Option<OpId>,
    besselk: Option<OpId>,
    besseli: Option<OpId>,
    pi: Option<OpId>,
    big_o: Option<OpId>,
    oo: Option<OpId>,
}

/// Installs the Meijer G operators and kernels.
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let meijerg = i.op(OpDescriptor::new("meijerg", Arity::Fixed(5))
        .flags(OpFlags::HEAVY)
        .cost(100))?;
    let mellin = i.op(OpDescriptor::new("meijer_mellin", Arity::Fixed(5))
        .flags(OpFlags::HEAVY)
        .cost(100))?;
    let series = i.op(OpDescriptor::new("meijer_series", Arity::Variadic)
        .flags(OpFlags::HEAVY)
        .cost(100))?;
    let invalid = |reason: &'static str| RuleError::Invalid {
        rule: "special/meijer".to_owned(),
        reason,
    };
    let exp = i
        .graph()
        .ops()
        .lookup("exp")
        .ok_or_else(|| invalid("`exp` is not registered"))?;
    let ln = i
        .graph()
        .ops()
        .lookup("ln")
        .ok_or_else(|| invalid("`ln` is not registered"))?;
    let ops = Ops {
        meijerg,
        exp,
        ln,
        gamma: i.graph().ops().lookup("gamma"),
        sin: i.graph().ops().lookup("sin"),
        cos: i.graph().ops().lookup("cos"),
        atan: i.graph().ops().lookup("atan"),
        erf: i.graph().ops().lookup("erf"),
        erfc: i.graph().ops().lookup("erfc"),
        besselj: i.graph().ops().lookup("besselj"),
        bessely: i.graph().ops().lookup("bessely"),
        besselk: i.graph().ops().lookup("besselk"),
        besseli: i.graph().ops().lookup("besseli"),
        pi: i.graph().ops().lookup("pi"),
        big_o: i.graph().ops().lookup("big_o"),
        oo: i.graph().ops().lookup("oo"),
    };
    i.kernel("special/meijerg-numeric", Tier::Reduce, Numeric { meijerg });
    i.kernel("special/meijerg-reduce", Tier::Reduce, Reduce { ops });
    i.kernel(
        "special/meijer-mellin",
        Tier::Reduce,
        Mellin { mellin, ops },
    );
    i.kernel(
        "special/meijer-series",
        Tier::Reduce,
        Series { series, ops },
    );
    if let Some(diff) = i.graph().ops().lookup("diff") {
        i.kernel(
            "special/meijerg-diff",
            Tier::Reduce,
            Differentiate { meijerg, diff, ops },
        );
    }
    if let Some(integral) = i.graph().ops().lookup("integral") {
        i.kernel(
            "special/meijerg-integral",
            Tier::Reduce,
            Integrate {
                meijerg,
                integral,
                ops,
            },
        );
    }
    Ok(())
}

// ------------------------------------------------------------ term building

fn sum_of(
    graph: &mut Graph,
    terms: &[NodeId],
) -> NodeId {
    match terms {
        | [] => graph.int(0),
        | [only] => *only,
        | _ => graph.node(core::ADD, terms),
    }
}

fn product_of(
    graph: &mut Graph,
    factors: &[NodeId],
) -> NodeId {
    match factors {
        | [] => graph.int(1),
        | [only] => *only,
        | _ => graph.node(core::MUL, factors),
    }
}

fn one(graph: &mut Graph) -> NodeId {
    graph.int(1)
}

fn int_node(
    graph: &mut Graph,
    value: i64,
) -> NodeId {
    graph.int(value)
}

fn half(graph: &mut Graph) -> NodeId {
    graph.num(Number::fraction(1, 2).unwrap_or_else(|| Number::from(1)))
}

fn add(
    graph: &mut Graph,
    a: NodeId,
    b: NodeId,
) -> NodeId {
    sum_of(graph, &[a, b])
}

fn mul(
    graph: &mut Graph,
    a: NodeId,
    b: NodeId,
) -> NodeId {
    product_of(graph, &[a, b])
}

fn mul_many(
    graph: &mut Graph,
    factors: &[NodeId],
) -> NodeId {
    product_of(graph, factors)
}

fn neg(
    graph: &mut Graph,
    a: NodeId,
) -> NodeId {
    let minus = int_node(graph, -1);
    product_of(graph, &[minus, a])
}

fn sub(
    graph: &mut Graph,
    a: NodeId,
    b: NodeId,
) -> NodeId {
    let minus = neg(graph, b);
    sum_of(graph, &[a, minus])
}

fn pow(
    graph: &mut Graph,
    base: NodeId,
    exponent: NodeId,
) -> NodeId {
    graph.node(core::POW, &[base, exponent])
}

fn recip(
    graph: &mut Graph,
    a: NodeId,
) -> NodeId {
    let minus = int_node(graph, -1);
    pow(graph, a, minus)
}

fn div(
    graph: &mut Graph,
    a: NodeId,
    b: NodeId,
) -> NodeId {
    let inverse = recip(graph, b);
    mul(graph, a, inverse)
}

fn sqrt_of(
    graph: &mut Graph,
    a: NodeId,
) -> NodeId {
    let h = half(graph);
    pow(graph, a, h)
}

fn plus_int(
    graph: &mut Graph,
    a: NodeId,
    k: i64,
) -> NodeId {
    let kk = int_node(graph, k);
    add(graph, a, kk)
}

fn sub_int(
    graph: &mut Graph,
    a: NodeId,
    k: i64,
) -> NodeId {
    let kk = int_node(graph, k);
    sub(graph, a, kk)
}

fn one_plus(
    graph: &mut Graph,
    a: NodeId,
) -> NodeId {
    plus_int(graph, a, 1)
}

fn one_minus(
    graph: &mut Graph,
    a: NodeId,
) -> NodeId {
    let n = neg(graph, a);
    plus_int(graph, n, 1)
}

fn exp_of(
    graph: &mut Graph,
    ops: &Ops,
    a: NodeId,
) -> NodeId {
    graph.node(ops.exp, &[a])
}

fn ln_of(
    graph: &mut Graph,
    ops: &Ops,
    a: NodeId,
) -> NodeId {
    graph.node(ops.ln, &[a])
}

fn gamma_of(
    graph: &mut Graph,
    ops: &Ops,
    a: NodeId,
) -> Option<NodeId> {
    Some(graph.node(ops.gamma?, &[a]))
}

fn pi_node(
    graph: &mut Graph,
    ops: &Ops,
) -> Option<NodeId> {
    Some(graph.node(ops.pi?, &[]))
}

fn sqrt_pi(
    graph: &mut Graph,
    ops: &Ops,
) -> Option<NodeId> {
    let pi = pi_node(graph, ops)?;
    Some(sqrt_of(graph, pi))
}

/// The exact real value of a parameter node.
fn value_of(
    graph: &Graph,
    node: NodeId,
) -> Option<f64> {
    if let Some(n) = graph.number_of(node) {
        let v = n.to_f64();
        return v.is_finite().then_some(v);
    }
    graph
        .eval(node, &Env::numeric(0.0))
        .filter(|v| v.is_finite())
}

fn near(
    x: f64,
    want: f64,
) -> bool {
    (x - want).abs() < 1e-9
}

/// Whether a node is exactly zero after simplification.
fn is_zero(
    cx: &mut Cx<'_>,
    node: NodeId,
) -> bool {
    let simplified = cx.simplify(node);
    cx.graph.number_of(simplified).is_some_and(Number::is_zero)
}

fn whole_index(
    graph: &Graph,
    node: NodeId,
) -> Option<usize> {
    graph
        .number_of(node)
        .and_then(Number::to_i64)
        .and_then(|v| usize::try_from(v).ok())
}

fn numeric_parameters(
    graph: &Graph,
    args: &GArgs,
) -> Option<(Vec<f64>, Vec<f64>)> {
    let a: Option<Vec<f64>> = args.a.iter().map(|&x| value_of(graph, x)).collect();
    let b: Option<Vec<f64>> = args.b.iter().map(|&x| value_of(graph, x)).collect();
    Some((a?, b?))
}

// ------------------------------------------------------- argument parsing

/// One `meijerg` request.
#[derive(Clone)]
struct GArgs {
    m: usize,
    n: usize,
    a: Vec<NodeId>,
    b: Vec<NodeId>,
    z: NodeId,
}

fn list_items(
    graph: &Graph,
    node: NodeId,
) -> Option<Vec<NodeId>> {
    if graph.op(node) != core::LIST {
        return None;
    }
    Some(graph.children(node).to_vec())
}

fn parse_g(
    graph: &Graph,
    node: NodeId,
) -> Option<GArgs> {
    let children = graph.children(node).to_vec();
    let m = whole_index(graph, *children.first()?)?;
    let n = whole_index(graph, *children.get(1)?)?;
    let a = list_items(graph, *children.get(2)?)?;
    let b = list_items(graph, *children.get(3)?)?;
    let z = *children.get(4)?;
    if m > b.len() || n > a.len() {
        return None;
    }
    Some(GArgs { m, n, a, b, z })
}

fn list_of(
    graph: &mut Graph,
    items: &[NodeId],
) -> NodeId {
    graph.node(core::LIST, items)
}

fn g_node(
    graph: &mut Graph,
    ops: &Ops,
    m: usize,
    n: usize,
    a: &[NodeId],
    b: &[NodeId],
    z: NodeId,
) -> Option<NodeId> {
    let m = int_node(graph, i64::try_from(m).ok()?);
    let n = int_node(graph, i64::try_from(n).ok()?);
    let alist = list_of(graph, a);
    let blist = list_of(graph, b);
    graph.try_node(ops.meijerg, &[m, n, alist, blist, z])
}

// -------------------------------------------------------------- numeric phase

/// Numeric witness of a fully bound `meijerg` request.
struct Numeric {
    meijerg: OpId,
}

impl Kernel for Numeric {
    fn ops(&self) -> Vec<OpId> {
        vec![self.meijerg]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        if !cx.env.numeric {
            return Outcome::Pass;
        }
        let graph = &mut *cx.graph;
        let Some(args) = parse_g(graph, node) else {
            return Outcome::Pass;
        };
        let free = graph.free_symbols(graph.find(node)).to_vec();
        if free.iter().any(|&s| cx.env.value(s).is_none()) {
            return Outcome::Pass;
        }
        let mut a = Vec::with_capacity(args.a.len());
        for &item in &args.a {
            let Some(v) = graph.eval(item, cx.env) else {
                return Outcome::Pass;
            };
            if !v.is_finite() {
                return Outcome::Pass;
            }
            a.push(v);
        }
        let mut b = Vec::with_capacity(args.b.len());
        for &item in &args.b {
            let Some(v) = graph.eval(item, cx.env) else {
                return Outcome::Pass;
            };
            if !v.is_finite() {
                return Outcome::Pass;
            }
            b.push(v);
        }
        let Some(at) = graph.eval(args.z, cx.env) else {
            return Outcome::Pass;
        };
        if !at.is_finite() {
            return Outcome::Pass;
        }
        let Some(value) = numeric::meijer_g_real(args.m, args.n, &a, &b, at) else {
            return Outcome::Pass;
        };
        if !value.is_finite() {
            return Outcome::Pass;
        }
        let rad = 1e-11 * value.abs().max(1.0);
        Outcome::Approx(crate::graph::Ball { mid: value, rad })
    }

    fn revisit(&self) -> bool {
        true
    }
}

// ------------------------------------------------------------ reductions

/// A constructor for a closed form, used both for the answer and for the
/// numeric self-check with concrete parameters.
type Build = fn(&mut Graph, &Ops, &[NodeId], &[NodeId], NodeId) -> Option<NodeId>;

/// A numeric sample: `a` values, `b` values and one or more `z` values.
type Sample = (&'static [f64], &'static [f64], &'static [f64]);

const EXP_SAMPLES: &[Sample] = &[
    (&[], &[0.0], &[0.31, -0.7, 2.4]),
    (&[], &[1.0], &[0.31, 1.9]),
    (&[], &[-1.3], &[0.31, 1.9]),
];

const FIRST_ORDER_SAMPLES: &[Sample] = &[
    (&[0.5], &[0.2], &[0.35, 1.4]),
    (&[2.3], &[1.0], &[0.35, 2.2]),
    (&[-0.4], &[1.4], &[0.35, 1.1]),
];

const COS_SAMPLES: &[Sample] = &[
    (&[], &[0.0, 0.5], &[0.3, 1.4, 3.2]),
    (&[], &[1.1, 1.6], &[0.6, 2.1]),
    (&[], &[-0.4, 0.1], &[0.4, 1.7]),
];

const SIN_SAMPLES: &[Sample] = &[
    (&[], &[0.5, 0.0], &[0.3, 1.4, 3.2]),
    (&[], &[1.6, 1.1], &[0.6, 2.1]),
    (&[], &[0.1, -0.4], &[0.4, 1.7]),
];

const ERF_SAMPLES: &[Sample] = &[(&[1.0], &[0.5, 0.0], &[0.2, 0.8, 2.1])];

const ERFC_SAMPLES: &[Sample] = &[(&[1.0], &[0.5, 0.0], &[0.2, 0.8, 2.1])];

const J_SAMPLES: &[Sample] = &[
    (&[], &[0.0, 0.0], &[0.4, 1.5]),
    (&[], &[0.25, -0.25], &[0.4, 1.5]),
    (&[], &[0.65, -0.65], &[0.4, 2.2]),
    (&[], &[1.2, -1.2], &[0.4, 1.5]),
];

const J0_SAMPLES: &[Sample] = &[
    (&[], &[0.7, 0.7], &[0.4, 1.5]),
    (&[], &[-0.3, -0.3], &[0.5, 2.4]),
];

const I_SAMPLES: &[Sample] = &[
    (&[], &[0.0, 0.0], &[-0.5, -1.2]),
    (&[], &[0.25, -0.25], &[-0.4, -1.5]),
    (&[], &[0.65, -0.65], &[-0.7, -2.0]),
];

const K_SAMPLES: &[Sample] = &[
    (&[], &[0.0, 0.0], &[0.4, 1.5]),
    (&[], &[0.35, -0.35], &[0.4, 1.5]),
    (&[], &[1.0, -1.0], &[0.4, 2.0]),
];

const Y_SAMPLES: &[Sample] = &[
    (&[-0.5], &[0.0, 0.0, -0.5], &[0.4, 1.5]),
    (&[-1.0], &[0.5, -0.5, -1.0], &[0.4, 1.5]),
    (&[-1.5], &[1.0, -1.0, -1.5], &[0.4, 1.5]),
    (&[-0.85], &[0.35, -0.35, -0.85], &[0.4, 1.5]),
];

const LOG_SAMPLES: &[Sample] = &[(&[1.0, 1.0], &[1.0, 0.0], &[-0.5, 0.5, 2.4])];

const ATAN_SAMPLES: &[Sample] = &[(&[0.5, 1.0], &[0.5, 0.0], &[0.2, 1.1, 4.0])];

/// Runs a reduction's numeric self-check at the given sample parameters.
fn check_build(
    graph: &mut Graph,
    ops: &Ops,
    m: usize,
    n: usize,
    build: Build,
    samples: &[Sample],
    tolerance: f64,
) -> bool {
    for &(sa, sb, sz) in samples {
        let ac: Vec<C> = sa.iter().map(|&v| C::new(v, 0.0)).collect();
        let bc: Vec<C> = sb.iter().map(|&v| C::new(v, 0.0)).collect();
        let an: Vec<NodeId> = sa.iter().map(|&v| graph.float(v)).collect();
        let bn: Vec<NodeId> = sb.iter().map(|&v| graph.float(v)).collect();
        for &zv in sz {
            let zn = graph.float(zv);
            let Some(candidate) = build(graph, ops, &an, &bn, zn) else {
                return false;
            };
            let Some(got) = graph.eval_complex(candidate, &HashMap::new()) else {
                return false;
            };
            let Some(want) = numeric::meijer_g(m, n, &ac, &bc, C::new(zv, 0.0)) else {
                return false;
            };
            let scale = want.norm().max(1.0);
            if (got - want).norm() > tolerance * scale {
                return false;
            }
        }
    }
    true
}

// ---- individual closed forms

fn build_exp(
    graph: &mut Graph,
    ops: &Ops,
    a: &[NodeId],
    b: &[NodeId],
    z: NodeId,
) -> Option<NodeId> {
    if !a.is_empty() || b.len() != 1 {
        return None;
    }
    let exponent = *b.first()?;
    let power = pow(graph, z, exponent);
    let minus_z = neg(graph, z);
    let decay = exp_of(graph, ops, minus_z);
    Some(mul(graph, power, decay))
}

fn build_first_order(
    graph: &mut Graph,
    ops: &Ops,
    a: &[NodeId],
    b: &[NodeId],
    z: NodeId,
) -> Option<NodeId> {
    if a.len() != 1 || b.len() != 1 {
        return None;
    }
    let a1 = *a.first()?;
    let b1 = *b.first()?;
    let difference = sub(graph, b1, a1);
    let gamma_arg = one_plus(graph, difference);
    let coefficient = gamma_of(graph, ops, gamma_arg)?;
    let z_power = pow(graph, z, b1);
    let base = one_plus(graph, z);
    let difference = sub(graph, a1, b1);
    let exponent = sub_int(graph, difference, 1);
    let power = pow(graph, base, exponent);
    Some(mul_many(graph, &[coefficient, z_power, power]))
}

/// `hi - lo - 1/2` close to zero.
fn half_step_numeric(
    graph: &Graph,
    hi: NodeId,
    lo: NodeId,
) -> bool {
    let (Some(h), Some(l)) = (value_of(graph, hi), value_of(graph, lo)) else {
        return false;
    };
    (h - l - 0.5).abs() < 1e-9
}

fn build_cos(
    graph: &mut Graph,
    ops: &Ops,
    a: &[NodeId],
    b: &[NodeId],
    z: NodeId,
) -> Option<NodeId> {
    if !a.is_empty() || b.len() != 2 {
        return None;
    }
    let b1 = *b.first()?;
    let b2 = *b.get(1)?;
    if !half_step_numeric(graph, b2, b1) {
        return None;
    }
    let cos = ops.cos?;
    let root = sqrt_of(graph, z);
    let two = int_node(graph, 2);
    let argument = mul(graph, two, root);
    let cosine = graph.node(cos, &[argument]);
    let scale_base = sqrt_pi(graph, ops)?;
    let scale = recip(graph, scale_base);
    let z_power = pow(graph, z, b1);
    Some(mul_many(graph, &[z_power, cosine, scale]))
}

fn build_sin(
    graph: &mut Graph,
    ops: &Ops,
    a: &[NodeId],
    b: &[NodeId],
    z: NodeId,
) -> Option<NodeId> {
    if !a.is_empty() || b.len() != 2 {
        return None;
    }
    let b1 = *b.first()?;
    let b2 = *b.get(1)?;
    if !half_step_numeric(graph, b1, b2) {
        return None;
    }
    let sin = ops.sin?;
    let root = sqrt_of(graph, z);
    let two = int_node(graph, 2);
    let argument = mul(graph, two, root);
    let sine = graph.node(sin, &[argument]);
    let hv = half(graph);
    let exponent = sub(graph, b1, hv);
    let scale_base = sqrt_pi(graph, ops)?;
    let scale = recip(graph, scale_base);
    let z_power = pow(graph, z, exponent);
    Some(mul_many(graph, &[z_power, sine, scale]))
}

fn erf_form(
    graph: &mut Graph,
    ops: &Ops,
    a: &[NodeId],
    b: &[NodeId],
    z: NodeId,
) -> Option<NodeId> {
    if a.len() != 1 || b.len() != 2 {
        return None;
    }
    let (a1, b1, b2) = (*a.first()?, *b.first()?, *b.get(1)?);
    if !near(value_of(graph, a1)?, 1.0)
        || !near(value_of(graph, b1)?, 0.5)
        || !near(value_of(graph, b2)?, 0.0)
    {
        return None;
    }
    let erf = ops.erf?;
    let root = sqrt_of(graph, z);
    let call = graph.node(erf, &[root]);
    let scale = sqrt_pi(graph, ops)?;
    Some(mul(graph, scale, call))
}

fn erfc_form(
    graph: &mut Graph,
    ops: &Ops,
    a: &[NodeId],
    b: &[NodeId],
    z: NodeId,
) -> Option<NodeId> {
    if a.len() != 1 || b.len() != 2 {
        return None;
    }
    let (a1, b1, b2) = (*a.first()?, *b.first()?, *b.get(1)?);
    if !near(value_of(graph, a1)?, 1.0)
        || !near(value_of(graph, b1)?, 0.5)
        || !near(value_of(graph, b2)?, 0.0)
    {
        return None;
    }
    let erfc = ops.erfc?;
    let root = sqrt_of(graph, z);
    let call = graph.node(erfc, &[root]);
    let scale = sqrt_pi(graph, ops)?;
    Some(mul(graph, scale, call))
}

fn bessel_argument(
    graph: &mut Graph,
    z: NodeId,
) -> NodeId {
    let root = sqrt_of(graph, z);
    let two = int_node(graph, 2);
    mul(graph, two, root)
}

/// Writes a node that is plainly negative as `-w` and returns `w`: a
/// negative literal, or a product with a negative numeric coefficient.
fn decompose_negative(
    graph: &mut Graph,
    node: NodeId,
) -> Option<NodeId> {
    if let Some(n) = graph.number_of(node) {
        let value = n.to_f64();
        return (value < 0.0).then(|| graph.float(-value));
    }
    if graph.op(node) != core::MUL {
        return None;
    }
    let children = graph.children(node).to_vec();
    let position = children
        .iter()
        .position(|&child| graph.number_of(child).is_some_and(|n| n.to_f64() < 0.0))?;
    let coefficient = graph.number_of(*children.get(position)?)?.to_f64();
    let mut rest: Vec<NodeId> = children
        .iter()
        .enumerate()
        .filter(|&(i, _)| i != position)
        .map(|(_, &child)| child)
        .collect();
    rest.push(graph.float(-coefficient));
    Some(product_of(graph, &rest))
}

/// `G^{1,0}_{0,2}(-w | nu/2, -nu/2) = i^nu I_nu(2 sqrt(w))`, the modified
/// Bessel form (`w = -z`).
fn build_bessel_i(
    graph: &mut Graph,
    ops: &Ops,
    a: &[NodeId],
    b: &[NodeId],
    z: NodeId,
) -> Option<NodeId> {
    if !a.is_empty() || b.len() != 2 {
        return None;
    }
    let w = decompose_negative(graph, z)?;
    let b1 = *b.first()?;
    let b2 = *b.get(1)?;
    let order = sub(graph, b1, b2);
    let minus_one = int_node(graph, -1);
    let imaginary = sqrt_of(graph, minus_one);
    let phase = pow(graph, imaginary, order);
    let argument = bessel_argument(graph, w);
    let call = graph.node(ops.besseli?, &[order, argument]);
    Some(mul(graph, phase, call))
}

fn build_bessel_j(
    graph: &mut Graph,
    ops: &Ops,
    a: &[NodeId],
    b: &[NodeId],
    z: NodeId,
) -> Option<NodeId> {
    if !a.is_empty() || b.len() != 2 {
        return None;
    }
    let b1 = *b.first()?;
    let b2 = *b.get(1)?;
    let j = ops.besselj?;
    if graph.same(b1, b2) {
        let argument = bessel_argument(graph, z);
        let order = int_node(graph, 0);
        let call = graph.node(j, &[order, argument]);
        let power = pow(graph, z, b1);
        return Some(mul(graph, power, call));
    }
    let order = sub(graph, b1, b2);
    let argument = bessel_argument(graph, z);
    Some(graph.node(j, &[order, argument]))
}

fn build_bessel_k(
    graph: &mut Graph,
    ops: &Ops,
    a: &[NodeId],
    b: &[NodeId],
    z: NodeId,
) -> Option<NodeId> {
    if !a.is_empty() || b.len() != 2 {
        return None;
    }
    let b1 = *b.first()?;
    let b2 = *b.get(1)?;
    let order = sub(graph, b1, b2);
    let argument = bessel_argument(graph, z);
    let call = graph.node(ops.besselk?, &[order, argument]);
    let scale = int_node(graph, 2);
    Some(mul(graph, scale, call))
}

fn build_bessel_y(
    graph: &mut Graph,
    ops: &Ops,
    a: &[NodeId],
    b: &[NodeId],
    z: NodeId,
) -> Option<NodeId> {
    if a.len() != 1 || b.len() != 3 {
        return None;
    }
    let b1 = *b.first()?;
    let b2 = *b.get(1)?;
    let order = sub(graph, b1, b2);
    let argument = bessel_argument(graph, z);
    Some(graph.node(ops.bessely?, &[order, argument]))
}

fn literal_pair(
    graph: &Graph,
    a: &[NodeId],
    b: &[NodeId],
    av: &[f64],
    bv: &[f64],
) -> bool {
    a.len() == av.len()
        && b.len() == bv.len()
        && a.iter()
            .zip(av)
            .all(|(&node, &want)| value_of(graph, node).is_some_and(|v| near(v, want)))
        && b.iter()
            .zip(bv)
            .all(|(&node, &want)| value_of(graph, node).is_some_and(|v| near(v, want)))
}

fn build_log(
    graph: &mut Graph,
    ops: &Ops,
    a: &[NodeId],
    b: &[NodeId],
    z: NodeId,
) -> Option<NodeId> {
    if !literal_pair(graph, a, b, &[1.0, 1.0], &[1.0, 0.0]) {
        return None;
    }
    let inside = one_plus(graph, z);
    Some(ln_of(graph, ops, inside))
}

fn build_atan(
    graph: &mut Graph,
    ops: &Ops,
    a: &[NodeId],
    b: &[NodeId],
    z: NodeId,
) -> Option<NodeId> {
    if !literal_pair(graph, a, b, &[0.5, 1.0], &[0.5, 0.0]) {
        return None;
    }
    let atan = ops.atan?;
    let root = sqrt_of(graph, z);
    let call = graph.node(atan, &[root]);
    let two = int_node(graph, 2);
    Some(mul(graph, two, call))
}

// ---- general terminating Slater series

fn factorial_i64(k: usize) -> Option<i64> {
    let mut acc: i64 = 1;
    for j in 2..=k {
        acc = acc.checked_mul(i64::try_from(j).ok()?)?;
    }
    Some(acc)
}

/// The full left Slater series when a numerator parameter terminates it.
/// Only installed for `p < q`, where the series is a global identity.
fn terminating_series(
    graph: &mut Graph,
    ops: &Ops,
    m: usize,
    n: usize,
    a: &[NodeId],
    b: &[NodeId],
    z: NodeId,
) -> Option<NodeId> {
    let p = a.len();
    let q = b.len();
    if p >= q || m == 0 || m > q || n > p {
        return None;
    }
    let gamma = ops.gamma?;
    let av: Vec<f64> = a
        .iter()
        .map(|&x| value_of(graph, x))
        .collect::<Option<Vec<f64>>>()?;
    let bv: Vec<f64> = b
        .iter()
        .map(|&x| value_of(graph, x))
        .collect::<Option<Vec<f64>>>()?;
    let parity = (p as i64 - m as i64 - n as i64).rem_euclid(2);
    let weight = if parity == 1 {
        neg(graph, z)
    } else {
        z
    };
    let mut families = Vec::with_capacity(m);
    for h in 0..m {
        let bh = *b.get(h)?;
        let mut avec = Vec::with_capacity(p);
        let mut aval = Vec::with_capacity(p);
        for (i, &ai) in a.iter().enumerate() {
            let one_plus_bh = one_plus(graph, bh);
            avec.push(sub(graph, one_plus_bh, ai));
            let bhv = *bv.get(h)?;
            let aiv = *av.get(i)?;
            aval.push(1.0 + bhv - aiv);
        }
        let mut degree: Option<usize> = None;
        for &v in &aval {
            if v <= 0.0 && near(v, v.round()) {
                let d = (-v).round() as usize;
                degree = Some(degree.map_or(d, |x| x.min(d)));
            }
        }
        let Some(limit) = degree else {
            // A non-terminating family makes the whole Slater series
            // infinite: only a terminating one can be written down.
            return None;
        };
        if limit > 24 {
            return None;
        }
        let mut bvec = Vec::with_capacity(q.saturating_sub(1));
        for (j, &bj) in b.iter().enumerate() {
            if j != h {
                let one_plus_bh = one_plus(graph, bh);
                bvec.push(sub(graph, one_plus_bh, bj));
            }
        }
        let mut numerator = Vec::with_capacity(m + n);
        for (j, &bj) in b.iter().enumerate().take(m) {
            if j != h {
                let difference = sub(graph, bj, bh);
                numerator.push(graph.node(gamma, &[difference]));
            }
        }
        for &ai in a.iter().take(n) {
            let one_plus_bh = one_plus(graph, bh);
            let difference = sub(graph, one_plus_bh, ai);
            numerator.push(graph.node(gamma, &[difference]));
        }
        let mut denominator = Vec::with_capacity(p + q);
        for &bj in b.iter().take(q).skip(m) {
            let one_plus_bh = one_plus(graph, bh);
            let difference = sub(graph, one_plus_bh, bj);
            denominator.push(graph.node(gamma, &[difference]));
        }
        for &ai in a.iter().take(p).skip(n) {
            let difference = sub(graph, ai, bh);
            denominator.push(graph.node(gamma, &[difference]));
        }
        let top = mul_many(graph, &numerator);
        let bottom = mul_many(graph, &denominator);
        let coefficient = div(graph, top, bottom);
        let mut terms = Vec::with_capacity(limit.saturating_add(1));
        let mut numerator_poch: Vec<NodeId> = vec![one(graph); avec.len()];
        let mut denominator_poch: Vec<NodeId> = vec![one(graph); bvec.len()];
        for k in 0..=limit {
            let kk = int_node(graph, i64::try_from(k).ok()?);
            let wk = pow(graph, weight, kk);
            let factorial = int_node(graph, factorial_i64(k)?);
            let np = mul_many(graph, &numerator_poch);
            let dp = mul_many(graph, &denominator_poch);
            let denominator_value = mul(graph, dp, factorial);
            let ratio = div(graph, np, denominator_value);
            terms.push(mul_many(graph, &[coefficient, wk, ratio]));
            let shift = int_node(graph, i64::try_from(k).ok()?);
            for (slot, &base) in numerator_poch.iter_mut().zip(&avec) {
                let next = add(graph, base, shift);
                *slot = mul(graph, *slot, next);
            }
            for (slot, &base) in denominator_poch.iter_mut().zip(&bvec) {
                let next = add(graph, base, shift);
                *slot = mul(graph, *slot, next);
            }
        }
        let z_power = pow(graph, z, bh);
        let inner = sum_of(graph, &terms);
        families.push(mul(graph, z_power, inner));
    }
    if families.is_empty() {
        return None;
    }
    Some(sum_of(graph, &families))
}

fn try_terminating(
    graph: &mut Graph,
    ops: &Ops,
    args: &GArgs,
) -> Option<NodeId> {
    let candidate = terminating_series(graph, ops, args.m, args.n, &args.a, &args.b, args.z)?;
    let Some((av, bv)) = numeric_parameters(graph, args) else {
        return Some(candidate);
    };
    let ac: Vec<C> = av.iter().map(|&v| C::new(v, 0.0)).collect();
    let bc: Vec<C> = bv.iter().map(|&v| C::new(v, 0.0)).collect();
    for z_test in [0.4, 1.3] {
        let expected = numeric::meijer_g(args.m, args.n, &ac, &bc, C::new(z_test, 0.0))?;
        let z_node = graph.float(z_test);
        let form = terminating_series(graph, ops, args.m, args.n, &args.a, &args.b, z_node)?;
        let got = graph.eval_complex(form, &HashMap::new())?;
        if (got - expected).norm() > 1e-7 * expected.norm().max(1.0) {
            return None;
        }
    }
    Some(candidate)
}

// ---- the reduction kernel

/// Exact reductions of `meijerg`, each numerically self-checked.
struct Reduce {
    ops: Ops,
}

impl Kernel for Reduce {
    fn ops(&self) -> Vec<OpId> {
        vec![self.ops.meijerg]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        if cx.env.numeric {
            return Outcome::Pass;
        }
        let Some(args) = parse_g(cx.graph, node) else {
            return Outcome::Pass;
        };
        reduce_request(cx, &self.ops, &args).map_or(Outcome::Pass, Outcome::Equal)
    }
}

fn reduce_request(
    cx: &mut Cx<'_>,
    ops: &Ops,
    args: &GArgs,
) -> Option<NodeId> {
    let (m, n, p, q) = (args.m, args.n, args.a.len(), args.b.len());
    // Symbolic zero tests need the engine, so they run before the builder
    // borrows the graph for the rest of the function.
    let mut order_sum_zero = [false, false, false];
    let mut third_equals_first = false;
    if p == 0
        && q == 2
        && let (Some(&b1), Some(&b2)) = (args.b.first(), args.b.get(1))
    {
        let sum = add(&mut *cx.graph, b1, b2);
        order_sum_zero[0] = is_zero(cx, sum);
        order_sum_zero[1] = order_sum_zero[0];
    }
    if p == 1
        && q == 3
        && let (&[a1], &[b1, b2, _]) = (args.a.as_slice(), args.b.as_slice())
    {
        third_equals_first = cx.graph.same(a1, args.b.get(2).copied()?);
        let sum = add(&mut *cx.graph, b1, b2);
        order_sum_zero[2] = is_zero(cx, sum);
    }
    let graph = &mut *cx.graph;
    // Each table entry is tried in turn and accepted only when its numeric
    // self-check passes.
    macro_rules! attempt {
        ($condition:expr, $build:ident, $samples:expr, $tolerance:expr) => {
            if $condition
                && let Some(candidate) = $build(graph, ops, &args.a, &args.b, args.z)
                && check_build(graph, ops, m, n, $build, $samples, $tolerance)
            {
                return Some(candidate);
            }
        };
    }
    attempt!(
        m == 1 && n == 0 && p == 0 && q == 1,
        build_exp,
        EXP_SAMPLES,
        1e-9
    );
    attempt!(
        m == 1 && n == 1 && p == 1 && q == 1,
        build_first_order,
        FIRST_ORDER_SAMPLES,
        1e-8
    );
    attempt!(
        m == 1 && n == 0 && p == 0 && q == 2,
        build_cos,
        COS_SAMPLES,
        1e-9
    );
    attempt!(
        m == 1 && n == 0 && p == 0 && q == 2,
        build_sin,
        SIN_SAMPLES,
        1e-9
    );
    attempt!(
        m == 1 && n == 0 && p == 0 && q == 2 && order_sum_zero[0],
        build_bessel_i,
        I_SAMPLES,
        1e-8
    );
    if m == 1
        && n == 0
        && p == 0
        && q == 2
        && let (Some(&b0), Some(&b1)) = (args.b.first(), args.b.get(1))
    {
        let same_pair = graph.same(b0, b1);
        if order_sum_zero[0] || same_pair {
            let samples = if same_pair {
                J0_SAMPLES
            } else {
                J_SAMPLES
            };
            if let Some(candidate) = build_bessel_j(graph, ops, &args.a, &args.b, args.z)
                && check_build(graph, ops, m, n, build_bessel_j, samples, 1e-8)
            {
                return Some(candidate);
            }
        }
    }
    attempt!(
        m == 2 && n == 0 && p == 0 && q == 2 && order_sum_zero[1],
        build_bessel_k,
        K_SAMPLES,
        1e-7
    );
    attempt!(
        m == 2 && n == 0 && p == 1 && q == 3 && order_sum_zero[2] && third_equals_first,
        build_bessel_y,
        Y_SAMPLES,
        1e-6
    );
    attempt!(
        m == 1 && n == 1 && p == 1 && q == 2,
        erf_form,
        ERF_SAMPLES,
        1e-8
    );
    attempt!(
        m == 2 && n == 0 && p == 1 && q == 2,
        erfc_form,
        ERFC_SAMPLES,
        1e-8
    );
    attempt!(
        m == 1 && n == 2 && p == 2 && q == 2,
        build_log,
        LOG_SAMPLES,
        1e-8
    );
    attempt!(
        m == 1 && n == 2 && p == 2 && q == 2,
        build_atan,
        ATAN_SAMPLES,
        1e-8
    );
    try_terminating(graph, ops, args)
}

// -------------------------------------------------------------- derivative

/// `diff(meijerg(...), x)` by the parameter shift of DLMF 16.18.4.
struct Differentiate {
    meijerg: OpId,
    diff: OpId,
    ops: Ops,
}

impl Kernel for Differentiate {
    fn ops(&self) -> Vec<OpId> {
        vec![self.diff]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        if cx.env.numeric {
            return Outcome::Pass;
        }
        let &[body, var] = cx.graph.children(node) else {
            return Outcome::Pass;
        };
        if cx.graph.op(body) != self.meijerg {
            return Outcome::Pass;
        }
        if has_reduced_form(cx.graph, body) {
            return Outcome::Pass;
        }
        let Some(args) = parse_g(cx.graph, body) else {
            return Outcome::Pass;
        };
        let Some(x) = cx.graph.symbol_of(var) else {
            return Outcome::Pass;
        };
        if !cx.graph.depends_on(cx.graph.find(args.z), x) {
            return Outcome::Equal(cx.graph.int(0));
        }
        let slope = if cx.graph.same(args.z, var) {
            Some(cx.graph.int(1))
        } else {
            linear_slope(cx, args.z, x)
        };
        let Some(slope) = slope else {
            return Outcome::Pass;
        };
        if cx.graph.number_of(args.z).is_some_and(Number::is_zero) {
            return Outcome::Pass;
        }
        let Some(shifted) = self.shifted(cx, &args) else {
            return Outcome::Pass;
        };
        if !verify_derivative(cx.graph, &args) {
            return Outcome::Pass;
        }
        let graph = &mut *cx.graph;
        let one = graph.int(1);
        let derivative = if graph.same(slope, one) {
            div(graph, shifted, args.z)
        } else {
            let scaled = mul(graph, slope, shifted);
            div(graph, scaled, args.z)
        };
        Outcome::Equal(derivative)
    }
}

impl Differentiate {
    /// `G^{m,n+1}_{p+1,q+1}(z | 0, a; b, 1)`.
    fn shifted(
        &self,
        cx: &mut Cx<'_>,
        args: &GArgs,
    ) -> Option<NodeId> {
        let graph = &mut *cx.graph;
        let zero = graph.int(0);
        let mut a = Vec::with_capacity(args.a.len().saturating_add(1));
        a.push(zero);
        a.extend_from_slice(&args.a);
        let mut b = Vec::with_capacity(args.b.len().saturating_add(1));
        b.extend_from_slice(&args.b);
        let one = graph.int(1);
        b.push(one);
        g_node(
            graph,
            &self.ops,
            args.m,
            args.n.saturating_add(1),
            &a,
            &b,
            args.z,
        )
    }
}

/// Finite-difference check of `z G' = G^{m,n+1}_{p+1,q+1}(z|0,a;b,1)`.
fn verify_derivative(
    graph: &Graph,
    args: &GArgs,
) -> bool {
    let Some((a, b)) = numeric_parameters(graph, args) else {
        return true;
    };
    let ac: Vec<C> = a.iter().map(|&v| C::new(v, 0.0)).collect();
    let bc: Vec<C> = b.iter().map(|&v| C::new(v, 0.0)).collect();
    let mut a2: Vec<C> = vec![C::new(0.0, 0.0)];
    a2.extend_from_slice(&ac);
    let mut b2 = bc.clone();
    b2.push(C::new(1.0, 0.0));
    for z in [0.41, 1.3, 2.2] {
        let h = 1e-5;
        let up = numeric::meijer_g(args.m, args.n, &ac, &bc, C::new(z + h, 0.0));
        let down = numeric::meijer_g(args.m, args.n, &ac, &bc, C::new(z - h, 0.0));
        let (Some(up), Some(down)) = (up, down) else {
            continue;
        };
        let lhs = (up - down) * (z / (2.0 * h));
        let Some(rhs) = numeric::meijer_g(args.m, args.n + 1, &a2, &b2, C::new(z, 0.0)) else {
            // The shifted G is not numerically available at this point; the
            // formula itself is the verified DLMF identity.
            continue;
        };
        if (lhs - rhs).norm() > 1e-6 * rhs.norm().max(1.0) {
            return false;
        }
    }
    true
}

/// Whether a class has a member whose root operator is not a heavy request
/// (i.e. a reduced form).
fn has_reduced_form(
    graph: &Graph,
    node: NodeId,
) -> bool {
    graph
        .enodes(graph.find(node))
        .any(|n| !graph.ops().get(graph.op(n)).flags.has(OpFlags::HEAVY))
}

/// The slope of `node = c x + d` in `x`, with `c`, `d` free of `x`.
fn linear_slope(
    cx: &mut Cx<'_>,
    node: NodeId,
    x: SymbolId,
) -> Option<NodeId> {
    if !cx.graph.depends_on(cx.graph.find(node), x) {
        return None;
    }
    let graph = &mut *cx.graph;
    let x_node = graph.symbol_node(x);
    let mut gens = Gens::default();
    let generator = gens.index(graph, x_node);
    let poly: Poly = from_term(graph, &mut gens, node, Limits::default())?;
    if poly.degree_in(generator) != 1 {
        return None;
    }
    for g in poly.support() {
        if g != generator
            && gens
                .node(g)
                .is_some_and(|n| graph.depends_on(graph.find(n), x))
        {
            return None;
        }
    }
    let coefficients = poly.coefficients_in(generator);
    let slope = coefficients.get(1)?.clone();
    Some(to_term(graph, &gens, &slope))
}

// ------------------------------------------------------------ antiderivative

/// `integral(meijerg(...), x)` by the `h = -1` parameter shift.
struct Integrate {
    meijerg: OpId,
    integral: OpId,
    ops: Ops,
}

impl Kernel for Integrate {
    fn ops(&self) -> Vec<OpId> {
        vec![self.integral]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        if cx.env.numeric {
            return Outcome::Pass;
        }
        let &[body, var] = cx.graph.children(node) else {
            return Outcome::Pass;
        };
        if cx.graph.op(body) != self.meijerg {
            return Outcome::Pass;
        }
        if has_reduced_form(cx.graph, body) {
            return Outcome::Pass;
        }
        let Some(args) = parse_g(cx.graph, body) else {
            return Outcome::Pass;
        };
        let Some(x) = cx.graph.symbol_of(var) else {
            return Outcome::Pass;
        };
        if !cx.graph.depends_on(cx.graph.find(args.z), x) {
            let graph = &mut *cx.graph;
            return Outcome::Equal(mul(graph, body, var));
        }
        let slope = if cx.graph.same(args.z, var) {
            Some(cx.graph.int(1))
        } else {
            linear_slope(cx, args.z, x)
        };
        let Some(slope) = slope else {
            return Outcome::Pass;
        };
        let Some(shifted) = self.shifted(cx, &args) else {
            return Outcome::Pass;
        };
        if !verify_primitive(cx.graph, &args) {
            return Outcome::Pass;
        }
        let graph = &mut *cx.graph;
        let one = graph.int(1);
        let primitive = mul(graph, args.z, shifted);
        let candidate = if graph.same(slope, one) {
            primitive
        } else {
            div(graph, primitive, slope)
        };
        Outcome::Equal(candidate)
    }
}

impl Integrate {
    /// `z G^{m,n+1}_{p+1,q+1}(z | 0, a; b, -1)`.
    fn shifted(
        &self,
        cx: &mut Cx<'_>,
        args: &GArgs,
    ) -> Option<NodeId> {
        let graph = &mut *cx.graph;
        let zero = graph.int(0);
        let mut a = Vec::with_capacity(args.a.len().saturating_add(1));
        a.push(zero);
        a.extend_from_slice(&args.a);
        let mut b = Vec::with_capacity(args.b.len().saturating_add(1));
        b.extend_from_slice(&args.b);
        let minus_one = graph.int(-1);
        b.push(minus_one);
        g_node(
            graph,
            &self.ops,
            args.m,
            args.n.saturating_add(1),
            &a,
            &b,
            args.z,
        )
    }
}

/// Checks `(t H(t))' = G(t)` for the shifted `H` at a few `t`.
fn verify_primitive(
    graph: &Graph,
    args: &GArgs,
) -> bool {
    let Some((a, b)) = numeric_parameters(graph, args) else {
        return true;
    };
    let ac: Vec<C> = a.iter().map(|&v| C::new(v, 0.0)).collect();
    let bc: Vec<C> = b.iter().map(|&v| C::new(v, 0.0)).collect();
    let mut a2: Vec<C> = vec![C::new(0.0, 0.0)];
    a2.extend_from_slice(&ac);
    let mut b2 = bc.clone();
    b2.push(C::new(-1.0, 0.0));
    let anti = |t: f64| -> Option<C> {
        let h = numeric::meijer_g(args.m, args.n + 1, &a2, &b2, C::new(t, 0.0))?;
        Some(h * t)
    };
    for t in [0.42, 1.3, 2.1] {
        let step = 1e-5;
        let (Some(up), Some(down)) = (anti(t + step), anti(t - step)) else {
            continue;
        };
        let derivative = (up - down) / (2.0 * step);
        let Some(want) = numeric::meijer_g(args.m, args.n, &ac, &bc, C::new(t, 0.0)) else {
            continue;
        };
        if (derivative - want).norm() > 1e-6 * want.norm().max(1.0) {
            return false;
        }
    }
    true
}

// ------------------------------------------------------------------ Mellin

/// `meijer_mellin(m, n, list(a), list(b), s)` to its gamma product.
struct Mellin {
    mellin: OpId,
    ops: Ops,
}

impl Kernel for Mellin {
    fn ops(&self) -> Vec<OpId> {
        vec![self.mellin]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let graph = &mut *cx.graph;
        let Some(args) = parse_g(graph, node) else {
            return Outcome::Pass;
        };
        let Some(gamma) = self.ops.gamma else {
            return Outcome::Pass;
        };
        let Some(&s) = graph.children(node).get(4) else {
            return Outcome::Pass;
        };
        let mut numerator = Vec::new();
        for &bj in args.b.iter().take(args.m) {
            let argument = add(graph, bj, s);
            numerator.push(graph.node(gamma, &[argument]));
        }
        for &ai in args.a.iter().take(args.n) {
            let one_minus_ai = one_minus(graph, ai);
            let argument = sub(graph, one_minus_ai, s);
            numerator.push(graph.node(gamma, &[argument]));
        }
        let mut denominator = Vec::new();
        for &bj in args.b.iter().skip(args.m) {
            let one_minus_bj = one_minus(graph, bj);
            let argument = sub(graph, one_minus_bj, s);
            denominator.push(graph.node(gamma, &[argument]));
        }
        for &ai in args.a.iter().skip(args.n) {
            let argument = add(graph, ai, s);
            denominator.push(graph.node(gamma, &[argument]));
        }
        let top = mul_many(graph, &numerator);
        let bottom = mul_many(graph, &denominator);
        Outcome::Equal(div(graph, top, bottom))
    }
}

// ------------------------------------------------------------------ series

/// `meijer_series(m, n, list(a), list(b), z, k[, direction])`.
struct Series {
    series: OpId,
    ops: Ops,
}

impl Kernel for Series {
    fn ops(&self) -> Vec<OpId> {
        vec![self.series]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let graph = &mut *cx.graph;
        let children = graph.children(node).to_vec();
        if !(children.len() == 6 || children.len() == 7) {
            return Outcome::Pass;
        }
        let Some(args) = parse_g(graph, node) else {
            return Outcome::Pass;
        };
        let Some(k_node) = children.get(5).copied() else {
            return Outcome::Pass;
        };
        let Some(count) = whole_index(graph, k_node) else {
            return Outcome::Pass;
        };
        if !(1..=12).contains(&count) {
            return Outcome::Pass;
        }
        let at_infinity = match children.get(6).copied() {
            | None => false,
            | Some(direction) => {
                if graph.number_of(direction).is_some_and(Number::is_zero) {
                    false
                } else if graph.op(direction) == self.ops.oo.unwrap_or(OpId::NONE) {
                    true
                } else {
                    return Outcome::Pass;
                }
            },
        };
        // The `1/z` series exists only when the contour can be closed to the
        // left; for `p < q` the G-function is entire and that series
        // diverges, so the request must stay unevaluated.
        if at_infinity && args.a.len() < args.b.len() {
            return Outcome::Pass;
        }
        let expression = if at_infinity {
            series_at_infinity(graph, &self.ops, &args, count)
        } else {
            series_at_zero(graph, &self.ops, &args, count)
        };
        let Some(expression) = expression else {
            return Outcome::Pass;
        };
        if !series_checks(graph, &self.ops, &args, count, at_infinity) {
            return Outcome::Pass;
        }
        if let Some(big_o) = self.ops.big_o
            && let Some(order) = series_order(graph, &args, count, at_infinity)
        {
            let term = graph.node(big_o, &[order]);
            let combined = sum_of(graph, &[expression, term]);
            return Outcome::Equal(combined);
        }
        Outcome::Equal(expression)
    }
}

/// Checks the truncated series against the kernel at a small (or large)
/// numeric argument, when the parameters are numeric.
fn series_checks(
    graph: &mut Graph,
    ops: &Ops,
    args: &GArgs,
    count: usize,
    at_infinity: bool,
) -> bool {
    let Some((av, bv)) = numeric_parameters(graph, args) else {
        return true;
    };
    let ac: Vec<C> = av.iter().map(|&v| C::new(v, 0.0)).collect();
    let bc: Vec<C> = bv.iter().map(|&v| C::new(v, 0.0)).collect();
    let z_test = if at_infinity { 1e4 } else { 1e-3 };
    let Some(want) = numeric::meijer_g(args.m, args.n, &ac, &bc, C::new(z_test, 0.0)) else {
        return false;
    };
    let zn = graph.float(z_test);
    let test_args = GArgs {
        m: args.m,
        n: args.n,
        a: args.a.clone(),
        b: args.b.clone(),
        z: zn,
    };
    let Some(form) = (if at_infinity {
        series_at_infinity(graph, ops, &test_args, count)
    } else {
        series_at_zero(graph, ops, &test_args, count)
    }) else {
        return false;
    };
    let Some(got) = graph.eval_complex(form, &HashMap::new()) else {
        return false;
    };
    (got - want).norm() <= 1e-4 * want.norm().max(1.0)
}

fn series_order(
    graph: &mut Graph,
    args: &GArgs,
    count: usize,
    at_infinity: bool,
) -> Option<NodeId> {
    let exponents: Vec<NodeId> = if at_infinity {
        args.a.iter().take(args.n).copied().collect()
    } else {
        args.b.iter().take(args.m).copied().collect()
    };
    if exponents
        .iter()
        .all(|&e| graph.number_of(e).is_some_and(Number::is_zero))
    {
        let k = graph.int(i64::try_from(count).ok()?);
        Some(if at_infinity {
            let power = pow(graph, args.z, k);
            recip(graph, power)
        } else {
            pow(graph, args.z, k)
        })
    } else {
        None
    }
}

/// The residue series at zero: `sum_h z^{b_h} C_h sum_j c_j z^j`.
fn series_at_zero(
    graph: &mut Graph,
    ops: &Ops,
    args: &GArgs,
    count: usize,
) -> Option<NodeId> {
    let gamma = ops.gamma?;
    let m = args.m;
    let n = args.n;
    let p = args.a.len();
    let q = args.b.len();
    let parity = (p as i64 - m as i64 - n as i64).rem_euclid(2);
    let weight = if parity == 1 {
        neg(graph, args.z)
    } else {
        args.z
    };
    let mut families = Vec::with_capacity(m);
    for h in 0..m {
        let bh = *args.b.get(h)?;
        let mut numerator = Vec::new();
        let mut denominator = Vec::new();
        for (j, &bj) in args.b.iter().enumerate().take(m) {
            if j != h {
                let difference = sub(graph, bj, bh);
                numerator.push(graph.node(gamma, &[difference]));
            }
        }
        for &ai in args.a.iter().take(n) {
            let one_plus_bh = one_plus(graph, bh);
            let difference = sub(graph, one_plus_bh, ai);
            numerator.push(graph.node(gamma, &[difference]));
        }
        for &bj in args.b.iter().skip(m) {
            let one_plus_bh = one_plus(graph, bh);
            let difference = sub(graph, one_plus_bh, bj);
            denominator.push(graph.node(gamma, &[difference]));
        }
        for &ai in args.a.iter().skip(n) {
            let difference = sub(graph, ai, bh);
            denominator.push(graph.node(gamma, &[difference]));
        }
        let top = mul_many(graph, &numerator);
        let bottom = mul_many(graph, &denominator);
        let coefficient = div(graph, top, bottom);
        let mut avec = Vec::with_capacity(p);
        for &ai in &args.a {
            let one_plus_bh = one_plus(graph, bh);
            avec.push(sub(graph, one_plus_bh, ai));
        }
        let mut bvec = Vec::with_capacity(q.saturating_sub(1));
        for (j, &bj) in args.b.iter().enumerate() {
            if j != h {
                let one_plus_bh = one_plus(graph, bh);
                bvec.push(sub(graph, one_plus_bh, bj));
            }
        }
        let mut terms = Vec::with_capacity(count);
        let mut numerator_poch: Vec<NodeId> = vec![one(graph); avec.len()];
        let mut denominator_poch: Vec<NodeId> = vec![one(graph); bvec.len()];
        for j in 0..count {
            let jj = int_node(graph, i64::try_from(j).ok()?);
            let power = pow(graph, weight, jj);
            let factorial = int_node(graph, factorial_i64(j)?);
            let np = mul_many(graph, &numerator_poch);
            let dp = mul_many(graph, &denominator_poch);
            let denominator_value = mul(graph, dp, factorial);
            let ratio = div(graph, np, denominator_value);
            terms.push(mul_many(graph, &[power, ratio]));
            let shift = int_node(graph, i64::try_from(j).ok()?);
            for (slot, &base) in numerator_poch.iter_mut().zip(&avec) {
                let next = add(graph, base, shift);
                *slot = mul(graph, *slot, next);
            }
            for (slot, &base) in denominator_poch.iter_mut().zip(&bvec) {
                let next = add(graph, base, shift);
                *slot = mul(graph, *slot, next);
            }
        }
        let z_power = pow(graph, args.z, bh);
        let inner = sum_of(graph, &terms);
        let family = mul(graph, z_power, coefficient);
        families.push(mul(graph, family, inner));
    }
    Some(sum_of(graph, &families))
}

/// The expansion at infinity: `sum_h z^{a_h-1} D_h sum_j d_j z^{-j}`.
fn series_at_infinity(
    graph: &mut Graph,
    ops: &Ops,
    args: &GArgs,
    count: usize,
) -> Option<NodeId> {
    let gamma = ops.gamma?;
    let m = args.m;
    let n = args.n;
    let q = args.b.len();
    let parity = (q as i64 - m as i64 - n as i64).rem_euclid(2);
    let inverse = recip(graph, args.z);
    let weight = if parity == 1 {
        neg(graph, inverse)
    } else {
        inverse
    };
    let mut families = Vec::with_capacity(n);
    for h in 0..n {
        let ah = *args.a.get(h)?;
        let mut numerator = Vec::new();
        let mut denominator = Vec::new();
        for (j, &aj) in args.a.iter().enumerate().take(n) {
            if j != h {
                let difference = sub(graph, ah, aj);
                numerator.push(graph.node(gamma, &[difference]));
            }
        }
        for &bj in args.b.iter().take(m) {
            let one_minus_ah = one_minus(graph, ah);
            let argument = add(graph, one_minus_ah, bj);
            numerator.push(graph.node(gamma, &[argument]));
        }
        for &aj in args.a.iter().skip(n) {
            let one_minus_ah = one_minus(graph, ah);
            let argument = add(graph, one_minus_ah, aj);
            denominator.push(graph.node(gamma, &[argument]));
        }
        for &bj in args.b.iter().skip(m) {
            let difference = sub(graph, ah, bj);
            denominator.push(graph.node(gamma, &[difference]));
        }
        let top = mul_many(graph, &numerator);
        let bottom = mul_many(graph, &denominator);
        let coefficient = div(graph, top, bottom);
        let mut avec = Vec::with_capacity(q);
        for &bj in &args.b {
            let one_minus_ah = one_minus(graph, ah);
            avec.push(add(graph, one_minus_ah, bj));
        }
        let mut bvec = Vec::with_capacity(args.a.len().saturating_sub(1));
        for (j, &aj) in args.a.iter().enumerate() {
            if j != h {
                let one_minus_ah = one_minus(graph, ah);
                bvec.push(add(graph, one_minus_ah, aj));
            }
        }
        let mut terms = Vec::with_capacity(count);
        let mut numerator_poch: Vec<NodeId> = vec![one(graph); avec.len()];
        let mut denominator_poch: Vec<NodeId> = vec![one(graph); bvec.len()];
        for j in 0..count {
            let jj = int_node(graph, i64::try_from(j).ok()?);
            let power = pow(graph, weight, jj);
            let factorial = int_node(graph, factorial_i64(j)?);
            let np = mul_many(graph, &numerator_poch);
            let dp = mul_many(graph, &denominator_poch);
            let denominator_value = mul(graph, dp, factorial);
            let ratio = div(graph, np, denominator_value);
            terms.push(mul_many(graph, &[power, ratio]));
            let shift = int_node(graph, i64::try_from(j).ok()?);
            for (slot, &base) in numerator_poch.iter_mut().zip(&avec) {
                let next = add(graph, base, shift);
                *slot = mul(graph, *slot, next);
            }
            for (slot, &base) in denominator_poch.iter_mut().zip(&bvec) {
                let next = add(graph, base, shift);
                *slot = mul(graph, *slot, next);
            }
        }
        let exponent = sub_int(graph, ah, 1);
        let z_power = pow(graph, args.z, exponent);
        let inner = sum_of(graph, &terms);
        let family = mul(graph, z_power, coefficient);
        families.push(mul(graph, family, inner));
    }
    Some(sum_of(graph, &families))
}

// ------------------------------------------------------------------- tests

#[cfg(test)]
mod tests {
    use crate::graph::Budget;
    use crate::graph::ClosedForm;
    use crate::graph::Engine;
    use crate::graph::Env;
    use crate::graph::Extractor;
    use crate::graph::Graph;
    use crate::graph::Saturate;
    use crate::graph::SizeCost;
    use crate::rules::standard;

    fn run_with(
        src: &str,
        numeric: bool,
    ) -> (String, bool, Option<f64>) {
        let mut graph = Graph::new();
        let sets = standard();
        let engine = Engine::install(&mut graph, &sets).unwrap_or_else(|e| panic!("{e}"));
        let root = graph
            .parse(src)
            .unwrap_or_else(|e| panic!("cannot parse `{src}`: {e}"));
        let env = if numeric {
            Env::numeric(1e-12)
        } else {
            Env::symbolic()
        };
        engine.run(&mut graph, &[root], &env, &Saturate, &Budget::default());
        assert!(graph.conflicts().is_empty(), "unsound: {src}");
        let closed = Extractor::new(&graph, &[root], &ClosedForm).build(&mut graph, root);
        let reduced = closed.is_some();
        let best = closed.unwrap_or_else(|| {
            Extractor::new(&graph, &[root], &SizeCost)
                .build(&mut graph, root)
                .unwrap_or(root)
        });
        let value = graph.approx(graph.find(root)).map(|b| b.mid);
        (graph.display(best), reduced, value)
    }

    fn closed(src: &str) -> String {
        let (text, reduced, _) = run_with(src, false);
        assert!(reduced, "`{src}` was not reduced: {text}");
        text
    }

    /// Numeric run with bindings, returning the root's witness.
    fn witness(
        src: &str,
        bindings: &[(&str, f64)],
    ) -> Option<f64> {
        let mut graph = Graph::new();
        let sets = standard();
        let engine = Engine::install(&mut graph, &sets).unwrap_or_else(|e| panic!("{e}"));
        let root = graph.parse(src).unwrap_or_else(|e| panic!("{e}"));
        let mut env = Env::numeric(1e-12);
        for &(name, value) in bindings {
            env.bind(graph.interner_mut().symbol(name), value);
        }
        engine.run(&mut graph, &[root], &env, &Saturate, &Budget::default());
        graph.approx(graph.find(root)).map(|b| b.mid)
    }

    fn eval_at(
        src: &str,
        bindings: &[(&str, f64)],
    ) -> f64 {
        let mut graph = Graph::new();
        let sets = standard();
        assert!(Engine::install(&mut graph, &sets).is_ok());
        let root = graph.parse(src).unwrap_or_else(|e| panic!("{e}"));
        let mut env = Env::numeric(0.0);
        for &(name, value) in bindings {
            env.bind(graph.interner_mut().symbol(name), value);
        }
        graph.eval(root, &env).unwrap_or(f64::NAN)
    }

    #[test]
    fn reductions_to_elementary_functions() {
        let cases: &[(&str, &[(&str, f64)])] = &[
            ("meijerg(1, 0, list(), list(0), z)", &[("z", 0.7)]),
            ("meijerg(1, 0, list(), list(3/2), z)", &[("z", 0.7)]),
            ("meijerg(1, 1, list(1/2), list(1/5), z)", &[("z", 0.7)]),
            ("meijerg(1, 0, list(), list(0, 1/2), z)", &[("z", 0.7)]),
            ("meijerg(1, 0, list(), list(1/2, 0), z)", &[("z", 0.7)]),
            ("meijerg(1, 1, list(1), list(1/2, 0), z)", &[("z", 0.7)]),
            ("meijerg(2, 0, list(1), list(1/2, 0), z)", &[("z", 0.7)]),
            ("meijerg(1, 2, list(1, 1), list(1, 0), z)", &[("z", 0.7)]),
            (
                "meijerg(1, 2, list(1/2, 1), list(1/2, 0), z)",
                &[("z", 0.7)],
            ),
        ];
        for &(src, bindings) in cases {
            let text = closed(src);
            assert!(!text.contains("meijerg"), "{src} -> {text}");
            let from_closed = eval_at(&text, bindings);
            let want = witness(src, bindings).unwrap_or(f64::NAN);
            assert!(
                (from_closed - want).abs() < 1e-7 * (1.0 + want.abs()),
                "{src} -> {text} = {from_closed}, kernel {want}"
            );
        }
    }

    #[test]
    fn bessel_reductions() {
        for nu in [0.0, 0.7, 2.0] {
            let src = format!("meijerg(1, 0, list(), list({nu}/2, -{nu}/2), z)");
            let text = closed(&src);
            assert!(text.contains("besselj"), "{src} -> {text}");
        }
        let text = closed("meijerg(2, 0, list(), list(nu/2, -nu/2), z)");
        assert!(text.contains("besselk"), "{text}");
        let text = closed("meijerg(2, 0, list((-nu - 1)/2), list(nu/2, -nu/2, (-nu - 1)/2), z)");
        assert!(text.contains("bessely"), "{text}");
        // A negative argument is the modified Bessel form: for nu = 0 it has
        // no complex phase and can be evaluated directly.
        let text = closed("meijerg(1, 0, list(), list(0, 0), -z)");
        assert!(text.contains("besseli"), "{text}");
        assert!(
            (eval_at(&text, &[("z", 0.5)])
                - crate::kernels::special::bessel_i(0.0, 2.0 * 0.5_f64.sqrt()))
            .abs()
                < 1e-9,
            "{text}"
        );
        let text = closed("meijerg(1, 0, list(), list(nu/2, -nu/2), -z)");
        assert!(text.contains("besseli"), "{text}");
    }

    #[test]
    fn derivative_and_antiderivative() {
        let derivative = closed("diff(meijerg(1, 0, list(), list(0), x), x)");
        assert!(!derivative.contains("diff"), "{derivative}");
        assert!(
            (eval_at(&derivative, &[("x", 0.6)]) + (-0.6_f64).exp()).abs() < 1e-9,
            "{derivative}"
        );
        let primitive = closed("integral(meijerg(1, 0, list(), list(0), x), x)");
        assert!(!primitive.contains("integral"), "{primitive}");
        let h = 1e-6;
        let up = eval_at(&primitive, &[("x", 0.6 + h)]);
        let down = eval_at(&primitive, &[("x", 0.6 - h)]);
        let numerical = (up - down) / (2.0 * h);
        assert!(
            (numerical - (-0.6_f64).exp()).abs() < 1e-5,
            "{primitive}: {numerical}"
        );
        let derivative = closed("diff(meijerg(1, 0, list(), list(1/3), x), x)");
        let up = eval_at("x^(1/3)*exp(-x)", &[("x", 1.2 + h)]);
        let down = eval_at("x^(1/3)*exp(-x)", &[("x", 1.2 - h)]);
        let numerical = (up - down) / (2.0 * h);
        assert!(
            (eval_at(&derivative, &[("x", 1.2)]) - numerical).abs() < 1e-5,
            "{derivative}: {numerical}"
        );
        let primitive = closed("integral(meijerg(1, 0, list(), list(0), 2*x), x)");
        assert!(!primitive.contains("integral"), "{primitive}");
        let up = eval_at(&primitive, &[("x", 0.3 + h)]);
        let down = eval_at(&primitive, &[("x", 0.3 - h)]);
        let numerical = (up - down) / (2.0 * h);
        assert!((numerical - (-0.6_f64).exp()).abs() < 1e-5, "{primitive}");
    }

    #[test]
    fn mellin_formula() {
        let text = closed("meijer_mellin(1, 0, list(), list(0), s)");
        assert!(text.contains("gamma"), "{text}");
        for s in [0.7, 1.3, 2.4] {
            let value = eval_at(&text, &[("s", s)]);
            let want = crate::kernels::special::gamma_numerical(s);
            assert!((value - want).abs() < 1e-10, "{s}: {value} vs {want}");
        }
    }

    #[test]
    fn residue_series_at_zero() {
        let text = closed("meijer_series(1, 0, list(), list(0), z, 4)");
        assert!(
            text.contains("big_o") || !text.contains("meijerg"),
            "{text}"
        );
        let value = eval_at(&text, &[("z", 1e-3)]);
        assert!((value - (-1e-3_f64).exp()).abs() < 1e-8, "{text}: {value}");
        let text = closed("meijer_series(1, 0, list(), list(1/2), z, 3)");
        let value = eval_at(&text, &[("z", 1e-4)]);
        let want = (1e-4_f64).sqrt() * (-1e-4_f64).exp();
        assert!((value - want).abs() < 1e-9, "{text}: {value} vs {want}");
        // The 1/z expansion does not exist for an entire case (p < q): the
        // request stays unevaluated rather than returning a divergent series.
        let (text, reduced, _) = run_with("meijer_series(1, 0, list(), list(0), z, 3, oo)", false);
        assert!(!reduced, "a divergent 1/z series was installed: {text}");
        assert!(text.contains("meijer_series"), "{text}");
    }

    #[test]
    fn invalid_orders_stay_unevaluated() {
        let (text, reduced, _) = run_with("meijerg(3, 0, list(0), list(0), z)", false);
        assert!(!reduced, "invalid order was reduced to {text}");
        assert!(text.contains("meijerg"), "{text}");
    }

    #[test]
    fn numeric_phase_uses_the_kernel() {
        let (text, _, value) = run_with("meijerg(1, 0, list(), list(0.5), 1.7)", true);
        let value = value.unwrap_or(f64::NAN);
        let want = 1.7_f64.sqrt() * (-1.7_f64).exp();
        assert!((value - want).abs() < 1e-10, "{text}: {value} vs {want}");
    }
}
