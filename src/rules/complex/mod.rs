//! Complex numbers.
//!
//! The imaginary unit is the nullary operator `I`; complex values are
//! ordinary terms such as `2 + 3*I` or `exp(I*x)`. Nothing is assumed about
//! a symbol unless the session says so: `re(z)` of a symbol `z` stays
//! `re(z)`, while for a symbol assumed real it is `z` itself.
//!
//! | operator | value |
//! |---|---|
//! | `I` | the imaginary unit; `I^n` reduces for integer `n` |
//! | `re(z)`, `im(z)`, `conj(z)` | real part, imaginary part, conjugate |
//! | `abs(z)`, `arg(z)` | modulus, principal argument in `(-pi, pi]` |
//!
//! The kernels split a term into real and imaginary parts structurally:
//! sums and products, integer powers, positive real bases with complex
//! exponents, and `exp`, `ln`, `sin`, `cos`, `sinh`, `cosh` of complex
//! arguments. `exp(I*t)` reduces to `cos(t) + I*sin(t)` when both sides
//! have closed values (`exp(I*pi) = -1`); otherwise Euler's formula is an
//! exploring rule in both directions.
//!
//! The real evaluator cannot see `I`, so complex results carry no numeric
//! witness. [`eval_complex`] evaluates terms over the complex numbers; the
//! complex-analysis kernels and the tests use it to check symbolic results.

pub mod analysis;
pub mod asymptotics;
pub mod branches;
pub mod conformal;
pub mod contour;
pub mod dynamics;
pub mod riemann;
pub mod transforms;

use std::collections::HashMap;

use num_complex::Complex64;

use crate::graph::Always;
use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::Facts;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpDescriptor;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::RuleError;
use crate::graph::RuleSet;
use crate::graph::SymbolId;
use crate::graph::Tier;
use crate::graph::op::core;
use crate::graph::rule::Installer;
use crate::rules::poly::best;

use super::calculus::calculus;
use super::elementary::elementary;

/// The complex-number rule set.
#[must_use]
pub fn complex() -> RuleSet {
    RuleSet::new("complex", install)
        .needs(elementary())
        .needs(calculus())
}

pub use crate::graph::ComplexEval;

/// Largest integer power expanded by the splitter.
const MAX_POWER: u64 = 32;

fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let unit = i.op(OpDescriptor::new("I", Arity::Fixed(0))
        .cost(1)
        .doc("`I`: the imaginary unit satisfying `I^2 = -1`."))?;
    let re = i.op(OpDescriptor::new("re", Arity::Fixed(1))
        .cost(8)
        .doc("`re(z)`: real part of complex number `z`."))?;
    let im = i.op(OpDescriptor::new("im", Arity::Fixed(1))
        .cost(8)
        .doc("`im(z)`: imaginary part of complex number `z`."))?;
    let conj = i.op(OpDescriptor::new("conj", Arity::Fixed(1))
        .cost(8)
        .doc("`conj(z)`: complex conjugate `re(z) - I * im(z)` of `z`."))?;
    let arg = i.op(OpDescriptor::new("arg", Arity::Fixed(1))
        .cost(8)
        .doc("`arg(z)`: principal argument (phase angle) of `z` in `(-pi, pi]`."))?;
    let complex_expand = i.op(OpDescriptor::new("complex_expand", Arity::Fixed(1))
        .flags(crate::graph::OpFlags::HEAVY)
        .cost(30)
        .doc("`complex_expand(z)`: expands `z` into explicit Cartesian form `re(z) + I*im(z)`."))?;
    let abs_sq = i.op(OpDescriptor::new("abs_sq", Arity::Fixed(1))
        .cost(8)
        .doc("`abs_sq(z)`: squared modulus `re(z)^2 + im(z)^2`."))?;
    let cis = i.op(OpDescriptor::new("cis", Arity::Fixed(1))
        .cost(8)
        .doc("`cis(t)`: unit complex number `cos(t) + I*sin(t)`."))?;
    let sign_complex = i.op(OpDescriptor::new("sign_complex", Arity::Fixed(1))
        .cost(8)
        .doc("`sign_complex(z)`: complex signum `z / abs(z)` on the unit circle."))?;
    let _to_polar = i.op(OpDescriptor::new("to_polar", Arity::Fixed(1))
        .cost(10)
        .doc("`to_polar(z)`: polar coordinates `list(abs(z), arg(z))`."))?;
    let from_polar = i.op(OpDescriptor::new("from_polar", Arity::Fixed(2))
        .cost(10)
        .doc("`from_polar(r, theta)`: reconstructs `r * (cos(theta) + I*sin(theta))`."))?;
    let _to_riemann_sphere = i.op(OpDescriptor::new("to_riemann_sphere", Arity::Fixed(1))
        .cost(12)
        .doc("`to_riemann_sphere(z)`: stereographic projection `(X, Y, Z)` onto the unit Riemann sphere `S^2`."))?;
    let from_riemann_sphere = i.op(OpDescriptor::new("from_riemann_sphere", Arity::Fixed(3))
        .cost(12)
        .doc("`from_riemann_sphere(X, Y, Z)`: inverse stereographic projection `(X + I*Y) / (1 - Z)` from `S^2` to `C`."))?;
    let chordal_distance = i.op(OpDescriptor::new("chordal_distance", Arity::Fixed(2))
        .cost(12)
        .doc("`chordal_distance(z1, z2)`: spherical chordal distance `2*abs(z1 - z2) / sqrt((1 + abs(z1)^2)*(1 + abs(z2)^2))` on `S^2`."))?;
    {
        let ops = i.graph().ops_mut();
        ops.set_attr(unit, Always(Facts::NONZERO));
        ops.set_attr(re, Always(Facts::REAL));
        ops.set_attr(im, Always(Facts::REAL));
        ops.set_attr(arg, Always(Facts::REAL));
        ops.set_attr(abs_sq, Always(Facts::NONNEGATIVE));
        ops.set_attr(chordal_distance, Always(Facts::NONNEGATIVE));
        ops.set_attr(
            cis,
            ComplexEval(|a| {
                a.first().map(|t| {
                    let i_unit = Complex64::new(0.0, 1.0);
                    (i_unit * *t).exp()
                })
            }),
        );
        ops.set_attr(
            abs_sq,
            ComplexEval(|a| a.first().map(|z| Complex64::new(z.norm_sqr(), 0.0))),
        );
        ops.set_attr(
            sign_complex,
            ComplexEval(|a| {
                a.first().map(|z| {
                    let n = z.norm();
                    if n == 0.0 {
                        Complex64::new(0.0, 0.0)
                    } else {
                        z / n
                    }
                })
            }),
        );
        ops.set_attr(
            from_polar,
            ComplexEval(|a| {
                let r = a.first()?;
                let t = a.get(1)?;
                let i_unit = Complex64::new(0.0, 1.0);
                Some(*r * (i_unit * *t).exp())
            }),
        );
        ops.set_attr(
            from_riemann_sphere,
            ComplexEval(|a| {
                let x = *a.first()?;
                let y = *a.get(1)?;
                let z = *a.get(2)?;
                let denom = Complex64::new(1.0, 0.0) - z;
                if denom.norm() < 1e-15 || !denom.is_finite() {
                    return None;
                }
                let i_unit = Complex64::new(0.0, 1.0);
                Some((x + i_unit * y) / denom)
            }),
        );
        ops.set_attr(
            chordal_distance,
            ComplexEval(|a| {
                let z1 = *a.first()?;
                let z2 = *a.get(1)?;
                Some(Complex64::new(
                    crate::kernels::complex_coords::chordal_distance(z1, z2),
                    0.0,
                ))
            }),
        );
    }
    if let Some(abs) = i.graph().ops().lookup("abs") {
        i.graph()
            .ops_mut()
            .set_attr(abs, Always(Facts::NONNEGATIVE));
    }
    register_complex_evals(i.graph());
    i.rewrites(
        Tier::Normalize,
        &[
            "complex/re-real: re(?x) => ?x if real(?x)",
            "complex/im-real: im(?x) => 0 if real(?x)",
            "complex/conj-real: conj(?x) => ?x if real(?x)",
            "complex/conj-conj: conj(conj(?x)) => ?x",
            // Conjugation is a field automorphism and commutes with every
            // function whose power series has real coefficients.
            "complex/conj-add: conj(?a + ?b) => conj(?a) + conj(?b)",
            "complex/conj-mul: conj(?a * ?b) => conj(?a) * conj(?b)",
            "complex/conj-pow: conj(?a ^ ?n) => conj(?a) ^ ?n if integer(?n)",
            "complex/conj-exp: conj(exp(?a)) => exp(conj(?a))",
            "complex/conj-sin: conj(sin(?a)) => sin(conj(?a))",
            "complex/conj-cos: conj(cos(?a)) => cos(conj(?a))",
            "complex/conj-sinh: conj(sinh(?a)) => sinh(conj(?a))",
            "complex/conj-cosh: conj(cosh(?a)) => cosh(conj(?a))",
            "complex/conj-tan: conj(tan(?a)) => tan(conj(?a))",
            "complex/conj-cot: conj(cot(?a)) => cot(conj(?a))",
            "complex/conj-sec: conj(sec(?a)) => sec(conj(?a))",
            "complex/conj-csc: conj(csc(?a)) => csc(conj(?a))",
            "complex/conj-tanh: conj(tanh(?a)) => tanh(conj(?a))",
            "complex/conj-coth: conj(coth(?a)) => coth(conj(?a))",
            "complex/conj-sech: conj(sech(?a)) => sech(conj(?a))",
            "complex/conj-csch: conj(csch(?a)) => csch(conj(?a))",
            "complex/conj-asin: conj(asin(?a)) => asin(conj(?a))",
            "complex/conj-acos: conj(acos(?a)) => acos(conj(?a))",
            "complex/conj-atan: conj(atan(?a)) => atan(conj(?a))",
            "complex/conj-asinh: conj(asinh(?a)) => asinh(conj(?a))",
            "complex/conj-acosh: conj(acosh(?a)) => acosh(conj(?a))",
            "complex/conj-atanh: conj(atanh(?a)) => atanh(conj(?a))",
            "complex/arg-positive: arg(?x) => 0 if positive(?x)",
            "complex/arg-negative: arg(?x) => pi if negative(?x)",
            "complex/abs-sq: abs_sq(?z) => re(?z)^2 + im(?z)^2",
            "complex/cis: cis(?t) => cos(?t) + I * sin(?t)",
            "complex/sign: sign_complex(?z) => ?z / abs(?z) if nonzero(?z)",
            "complex/to-polar: to_polar(?z) => list(abs(?z), arg(?z))",
            "complex/from-polar: from_polar(?r, ?t) => ?r * (cos(?t) + I * sin(?t))",
            "complex/to-sphere: to_riemann_sphere(?z) => list(2*re(?z)/(abs(?z)^2 + 1), 2*im(?z)/(abs(?z)^2 + 1), (abs(?z)^2 - 1)/(abs(?z)^2 + 1))",
            "complex/from-sphere: from_riemann_sphere(?x, ?y, ?z) => (?x + I * ?y) / (1 - ?z)",
            "complex/chordal: chordal_distance(?a, ?b) => 2 * abs(?a - ?b) / sqrt((1 + abs(?a)^2) * (1 + abs(?b)^2))",
        ],
    )?;
    i.rewrites(
        Tier::Explore,
        &["complex/euler: exp(I * ?x) <=> cos(?x) + I * sin(?x)"],
    )?;
    let ops = Ops::of(i.graph()).ok_or_else(|| {
        RuleError::Invalid {
            rule: "complex".into(),
            reason: "needs elementary",
        }
    })?;
    branches::install(i)?;
    riemann::install(i)?;
    contour::install(i)?;
    asymptotics::install(i)?;
    conformal::install(i)?;
    transforms::install(i)?;
    dynamics::install(i)?;
    analysis::install(i)?;
    i.kernel("complex/unit-power", Tier::Normalize, UnitPower { unit });
    i.kernel(
        "complex/parts",
        Tier::Normalize,
        Parts {
            ops,
            re,
            im,
            conj,
            arg,
        },
    );
    i.kernel(
        "complex/expand",
        Tier::Reduce,
        ComplexExpand {
            op: complex_expand,
            ops,
        },
    );
    i.kernel("complex/exp-i", Tier::Normalize, ExpI { ops });
    Ok(())
}

/// Operators the splitter needs.
#[derive(Copy, Clone, Debug)]
pub(crate) struct Ops {
    pub(crate) unit: OpId,
    pub(crate) re: OpId,
    pub(crate) im: OpId,
    pub(crate) exp: OpId,
    pub(crate) ln: OpId,
    pub(crate) sqrt: OpId,
    pub(crate) sin: OpId,
    pub(crate) cos: OpId,
    pub(crate) tan: OpId,
    pub(crate) cot: OpId,
    pub(crate) sec: OpId,
    pub(crate) csc: OpId,
    pub(crate) sinh: OpId,
    pub(crate) cosh: OpId,
    pub(crate) tanh: OpId,
    pub(crate) coth: OpId,
    pub(crate) sech: OpId,
    pub(crate) csch: OpId,
    pub(crate) asin: OpId,
    pub(crate) acos: OpId,
    pub(crate) atan: OpId,
    pub(crate) asinh: OpId,
    pub(crate) acosh: OpId,
    pub(crate) atanh: OpId,
    pub(crate) atan2: OpId,
    pub(crate) abs: OpId,
    pub(crate) conj: OpId,
    pub(crate) arg: OpId,
    pub(crate) pi: OpId,
}

impl Ops {
    pub(crate) fn of(graph: &Graph) -> Option<Self> {
        let get = |name: &str| graph.ops().lookup(name);
        Some(Self {
            unit: get("I")?,
            re: get("re")?,
            im: get("im")?,
            exp: get("exp")?,
            ln: get("ln")?,
            sqrt: get("sqrt")?,
            sin: get("sin")?,
            cos: get("cos")?,
            tan: get("tan")?,
            cot: get("cot")?,
            sec: get("sec")?,
            csc: get("csc")?,
            sinh: get("sinh")?,
            cosh: get("cosh")?,
            tanh: get("tanh")?,
            coth: get("coth")?,
            sech: get("sech")?,
            csch: get("csch")?,
            asin: get("asin")?,
            acos: get("acos")?,
            atan: get("atan")?,
            asinh: get("asinh")?,
            acosh: get("acosh")?,
            atanh: get("atanh")?,
            atan2: get("atan2")?,
            abs: get("abs")?,
            conj: get("conj")?,
            arg: get("arg")?,
            pi: get("pi")?,
        })
    }
}

/// Term builders shared by the complex kernels.
pub(crate) mod build {
    use super::Graph;
    use super::NodeId;
    use super::OpId;
    use super::core;

    pub fn add(
        graph: &mut Graph,
        terms: &[NodeId],
    ) -> NodeId {
        match terms {
            | [] => graph.int(0),
            | [only] => *only,
            | _ => graph.node(core::ADD, terms),
        }
    }

    pub fn mul(
        graph: &mut Graph,
        factors: &[NodeId],
    ) -> NodeId {
        match factors {
            | [] => graph.int(1),
            | [only] => *only,
            | _ => graph.node(core::MUL, factors),
        }
    }

    pub fn neg(
        graph: &mut Graph,
        x: NodeId,
    ) -> NodeId {
        let minus_one = graph.int(-1);
        mul(graph, &[minus_one, x])
    }

    pub fn sub(
        graph: &mut Graph,
        a: NodeId,
        b: NodeId,
    ) -> NodeId {
        let negated = neg(graph, b);
        add(graph, &[a, negated])
    }

    pub fn pow(
        graph: &mut Graph,
        base: NodeId,
        exponent: NodeId,
    ) -> NodeId {
        graph.node(core::POW, &[base, exponent])
    }

    pub fn powi(
        graph: &mut Graph,
        base: NodeId,
        exponent: i64,
    ) -> NodeId {
        let e = graph.int(exponent);
        pow(graph, base, e)
    }

    pub fn call(
        graph: &mut Graph,
        op: OpId,
        args: &[NodeId],
    ) -> NodeId {
        graph.node(op, args)
    }

    /// `a + b*I`.
    pub fn complex(
        graph: &mut Graph,
        unit: OpId,
        re: NodeId,
        im: NodeId,
    ) -> NodeId {
        let i = graph.node(unit, &[]);
        let imaginary = mul(graph, &[im, i]);
        add(graph, &[re, imaginary])
    }
}

use build::add;
use build::call;
use build::mul;
use build::neg;
use build::sub;

/// Real and imaginary parts of the best form of `z`, as (unsimplified)
/// real terms, when they can be read off its structure.
pub(crate) fn split(
    graph: &mut Graph,
    ops: Ops,
    z: NodeId,
) -> Option<(NodeId, NodeId)> {
    let term = best(graph, z)?;
    let mut memo = HashMap::new();
    split_term(graph, ops, term, &mut memo, 0)
}

#[allow(clippy::too_many_lines)]
fn split_term(
    graph: &mut Graph,
    ops: Ops,
    term: NodeId,
    memo: &mut HashMap<NodeId, Option<(NodeId, NodeId)>>,
    depth: usize,
) -> Option<(NodeId, NodeId)> {
    if let Some(known) = memo.get(&term) {
        return *known;
    }
    if depth > 64 {
        return None;
    }
    let zero = graph.int(0);
    let result = if graph.facts(term).has(Facts::REAL) {
        Some((term, zero))
    } else {
        let op = graph.op(term);
        let children = graph.children(term).to_vec();
        let mut part = |graph: &mut Graph, n: NodeId| split_term(graph, ops, n, memo, depth + 1);
        match (op, children.as_slice()) {
            | (op, []) if op == ops.unit => Some((zero, graph.int(1))),
            | (core::ADD, terms) => {
                let mut res = Vec::with_capacity(terms.len());
                let mut ims = Vec::with_capacity(terms.len());
                for &t in terms {
                    let (a, b) = part(graph, t)?;
                    res.push(a);
                    ims.push(b);
                }
                Some((add(graph, &res), add(graph, &ims)))
            },
            | (core::MUL, factors) => {
                let mut acc = (graph.int(1), zero);
                for &f in factors {
                    let next = part(graph, f)?;
                    acc = product(graph, acc, next);
                }
                Some(acc)
            },
            | (core::POW, &[base, exponent]) => {
                let n = graph.number_of(exponent).and_then(Number::to_i64);
                match n {
                    | Some(n) if n != 0 && n.unsigned_abs() <= MAX_POWER => {
                        let parts = part(graph, base)?;
                        let mut acc = (graph.int(1), zero);
                        for _ in 0..n.unsigned_abs() {
                            acc = product(graph, acc, parts);
                        }
                        if n < 0 {
                            // 1/(a + bi) = (a - bi)/(a^2 + b^2)
                            let (a, b) = acc;
                            let a2 = build::powi(graph, a, 2);
                            let b2 = build::powi(graph, b, 2);
                            let norm = add(graph, &[a2, b2]);
                            let inverse = build::powi(graph, norm, -1);
                            let re = mul(graph, &[a, inverse]);
                            let minus_b = neg(graph, b);
                            let im = mul(graph, &[minus_b, inverse]);
                            Some((re, im))
                        } else {
                            Some(acc)
                        }
                    },
                    | _ if graph.facts(base).has(Facts::POSITIVE) => {
                        // a^(x + iy) = a^x (cos(y ln a) + i sin(y ln a))
                        let (x, y) = part(graph, exponent)?;
                        let modulus = build::pow(graph, base, x);
                        let log = call(graph, ops.ln, &[base]);
                        let angle = mul(graph, &[y, log]);
                        Some(polar(graph, ops, modulus, angle))
                    },
                    | _ => {
                        let (a, b) = part(graph, base)?;
                        let (c, d) = part(graph, exponent)?;
                        let (lr, li) = split_ln(graph, ops, (a, b));
                        let (pr, pi_part) = product(graph, (c, d), (lr, li));
                        let modulus = call(graph, ops.exp, &[pr]);
                        Some(polar(graph, ops, modulus, pi_part))
                    },
                }
            },
            | (op, &[w]) if op == ops.exp => {
                let (a, b) = part(graph, w)?;
                let modulus = call(graph, ops.exp, &[a]);
                Some(polar(graph, ops, modulus, b))
            },
            | (op, &[w]) if op == ops.ln => {
                let (a, b) = part(graph, w)?;
                Some(split_ln(graph, ops, (a, b)))
            },
            | (op, &[w]) if op == ops.sqrt => {
                let (a, b) = part(graph, w)?;
                Some(split_sqrt(graph, ops, (a, b)))
            },
            | (op, &[w]) if op == ops.sin || op == ops.cos || op == ops.sinh || op == ops.cosh => {
                let (a, b) = part(graph, w)?;
                Some(split_trig_hyp(graph, ops, op, a, b))
            },
            | (op, &[w]) if op == ops.tan || op == ops.cot || op == ops.sec || op == ops.csc => {
                let (a, b) = part(graph, w)?;
                let (sa, ca) = split_trig_hyp(graph, ops, ops.sin, a, b);
                let (sb, cb) = split_trig_hyp(graph, ops, ops.cos, a, b);
                let one = (graph.int(1), zero);
                Some(if op == ops.tan {
                    cdiv(graph, (sa, ca), (sb, cb))
                } else if op == ops.cot {
                    cdiv(graph, (sb, cb), (sa, ca))
                } else if op == ops.sec {
                    cdiv(graph, one, (sb, cb))
                } else {
                    cdiv(graph, one, (sa, ca))
                })
            },
            | (op, &[w]) if op == ops.tanh || op == ops.coth || op == ops.sech || op == ops.csch =>
            {
                let (a, b) = part(graph, w)?;
                let (sa, ca) = split_trig_hyp(graph, ops, ops.sinh, a, b);
                let (sb, cb) = split_trig_hyp(graph, ops, ops.cosh, a, b);
                let one = (graph.int(1), zero);
                Some(if op == ops.tanh {
                    cdiv(graph, (sa, ca), (sb, cb))
                } else if op == ops.coth {
                    cdiv(graph, (sb, cb), (sa, ca))
                } else if op == ops.sech {
                    cdiv(graph, one, (sb, cb))
                } else {
                    cdiv(graph, one, (sa, ca))
                })
            },
            | (op, &[w]) if op == ops.atan => {
                // atan(a + ib) = 0.5*(atan2(a, 1-b) - atan2(-a, 1+b)) + (i/4)*ln(((1+b)^2+a^2)/((1-b)^2+a^2))
                let (a, b) = part(graph, w)?;
                let one = graph.int(1);
                let one_minus_b = sub(graph, one, b);
                let one_plus_b = add(graph, &[one, b]);
                let neg_a = neg(graph, a);
                let half = graph.num(Number::fraction(1, 2).unwrap_or_else(|| Number::from(0)));
                let quarter = graph.num(Number::fraction(1, 4).unwrap_or_else(|| Number::from(0)));
                let t1 = call(graph, ops.atan2, &[a, one_minus_b]);
                let t2 = call(graph, ops.atan2, &[neg_a, one_plus_b]);
                let diff_t = sub(graph, t1, t2);
                let re = mul(graph, &[half, diff_t]);
                let a2 = build::powi(graph, a, 2);
                let p_plus = build::powi(graph, one_plus_b, 2);
                let p_minus = build::powi(graph, one_minus_b, 2);
                let num = add(graph, &[p_plus, a2]);
                let den = add(graph, &[p_minus, a2]);
                let den_inv = build::powi(graph, den, -1);
                let ratio = mul(graph, &[num, den_inv]);
                let ln_ratio = call(graph, ops.ln, &[ratio]);
                let im = mul(graph, &[quarter, ln_ratio]);
                Some((re, im))
            },
            | (op, &[w]) if op == ops.atanh => {
                let (a, b) = part(graph, w)?;
                let one = graph.int(1);
                let one_plus_a = add(graph, &[one, a]);
                let one_minus_a = sub(graph, one, a);
                let neg_b = neg(graph, b);
                let half = graph.num(Number::fraction(1, 2).unwrap_or_else(|| Number::from(0)));
                let quarter = graph.num(Number::fraction(1, 4).unwrap_or_else(|| Number::from(0)));
                let b2 = build::powi(graph, b, 2);
                let p_plus = build::powi(graph, one_plus_a, 2);
                let p_minus = build::powi(graph, one_minus_a, 2);
                let num = add(graph, &[p_plus, b2]);
                let den = add(graph, &[p_minus, b2]);
                let den_inv = build::powi(graph, den, -1);
                let ratio = mul(graph, &[num, den_inv]);
                let ln_ratio = call(graph, ops.ln, &[ratio]);
                let re = mul(graph, &[quarter, ln_ratio]);
                let t1 = call(graph, ops.atan2, &[b, one_plus_a]);
                let t2 = call(graph, ops.atan2, &[neg_b, one_minus_a]);
                let diff_t = sub(graph, t1, t2);
                let im = mul(graph, &[half, diff_t]);
                Some((re, im))
            },
            | (op, &[w])
                if op == ops.asin || op == ops.acos || op == ops.asinh || op == ops.acosh =>
            {
                let (a, b) = part(graph, w)?;
                let (a2, b2) = product(graph, (a, b), (a, b));
                let one = graph.int(1);
                if op == ops.asinh {
                    // ln(z + sqrt(z^2 + 1))
                    let a2_p1 = add(graph, &[a2, one]);
                    let (sr, si) = split_sqrt(graph, ops, (a2_p1, b2));
                    let sum_r = add(graph, &[a, sr]);
                    let sum_i = add(graph, &[b, si]);
                    Some(split_ln(graph, ops, (sum_r, sum_i)))
                } else if op == ops.acosh {
                    // ln(z + sqrt(z-1)*sqrt(z+1))
                    let a_p1 = add(graph, &[a, one]);
                    let a_m1 = sub(graph, a, one);
                    let (p1r, p1i) = split_sqrt(graph, ops, (a_p1, b));
                    let (m1r, m1i) = split_sqrt(graph, ops, (a_m1, b));
                    let (pr, pi_p) = product(graph, (p1r, p1i), (m1r, m1i));
                    let sum_r = add(graph, &[a, pr]);
                    let sum_i = add(graph, &[b, pi_p]);
                    Some(split_ln(graph, ops, (sum_r, sum_i)))
                } else {
                    // asin(z) = -i * ln(i*z + sqrt(1 - z^2))
                    let one_m_a2 = sub(graph, one, a2);
                    let neg_b2 = neg(graph, b2);
                    let (sr, si) = split_sqrt(graph, ops, (one_m_a2, neg_b2));
                    let arg_r = sub(graph, sr, b);
                    let arg_i = add(graph, &[a, si]);
                    let (lr, li) = split_ln(graph, ops, (arg_r, arg_i));
                    if op == ops.asin {
                        Some((li, neg(graph, lr)))
                    } else {
                        let pi_node = call(graph, ops.pi, &[]);
                        let half =
                            graph.num(Number::fraction(1, 2).unwrap_or_else(|| Number::from(0)));
                        let hp = mul(graph, &[half, pi_node]);
                        Some((sub(graph, hp, li), lr))
                    }
                }
            },
            | (op, &[y, x]) if op == ops.atan2 => {
                // real for real arguments
                let (_, by) = part(graph, y)?;
                let (_, bx) = part(graph, x)?;
                let zero_part =
                    |graph: &Graph, n: NodeId| graph.number_of(n).is_some_and(Number::is_zero);
                (zero_part(graph, by) && zero_part(graph, bx)).then_some((term, zero))
            },
            | (op, &[w]) if op == ops.re => {
                let (a, _) = part(graph, w)?;
                Some((a, zero))
            },
            | (op, &[w]) if op == ops.im => {
                let (_, b) = part(graph, w)?;
                Some((b, zero))
            },
            | (op, &[w]) if op == ops.conj => {
                let (a, b) = part(graph, w)?;
                Some((a, neg(graph, b)))
            },
            | (op, &[w]) if op == ops.abs => {
                let (a, b) = part(graph, w)?;
                Some((modulus(graph, a, b), zero))
            },
            | (op, &[w]) if op == ops.arg => {
                let (a, b) = part(graph, w)?;
                Some((call(graph, ops.atan2, &[b, a]), zero))
            },
            | _ => None,
        }
    };
    memo.insert(term, result);
    result
}

/// `(a + bi) / (c + di) = ((ac + bd) + i(bc - ad)) / (c^2 + d^2)`.
fn cdiv(
    graph: &mut Graph,
    (a, b): (NodeId, NodeId),
    (c, d): (NodeId, NodeId),
) -> (NodeId, NodeId) {
    let c2 = build::powi(graph, c, 2);
    let d2 = build::powi(graph, d, 2);
    let denom = add(graph, &[c2, d2]);
    let inv = build::powi(graph, denom, -1);
    let ac = mul(graph, &[a, c]);
    let bd = mul(graph, &[b, d]);
    let bc = mul(graph, &[b, c]);
    let ad = mul(graph, &[a, d]);
    let sum_re = add(graph, &[ac, bd]);
    let diff_im = sub(graph, bc, ad);
    let re = mul(graph, &[sum_re, inv]);
    let im = mul(graph, &[diff_im, inv]);
    (re, im)
}

fn split_ln(
    graph: &mut Graph,
    ops: Ops,
    (a, b): (NodeId, NodeId),
) -> (NodeId, NodeId) {
    let m = modulus(graph, a, b);
    (call(graph, ops.ln, &[m]), call(graph, ops.atan2, &[b, a]))
}

fn split_sqrt(
    graph: &mut Graph,
    ops: Ops,
    (a, b): (NodeId, NodeId),
) -> (NodeId, NodeId) {
    let m = modulus(graph, a, b);
    let half = graph.num(Number::fraction(1, 2).unwrap_or_else(|| Number::from(0)));
    let r = build::pow(graph, m, half);
    let at = call(graph, ops.atan2, &[b, a]);
    let angle = mul(graph, &[half, at]);
    polar(graph, ops, r, angle)
}

fn split_trig_hyp(
    graph: &mut Graph,
    ops: Ops,
    op: OpId,
    a: NodeId,
    b: NodeId,
) -> (NodeId, NodeId) {
    if op == ops.sin {
        let sa = call(graph, ops.sin, &[a]);
        let cb = call(graph, ops.cosh, &[b]);
        let ca = call(graph, ops.cos, &[a]);
        let sb = call(graph, ops.sinh, &[b]);
        (mul(graph, &[sa, cb]), mul(graph, &[ca, sb]))
    } else if op == ops.cos {
        let sa = call(graph, ops.sin, &[a]);
        let sb = call(graph, ops.sinh, &[b]);
        let ca = call(graph, ops.cos, &[a]);
        let cb = call(graph, ops.cosh, &[b]);
        let im = mul(graph, &[sa, sb]);
        (mul(graph, &[ca, cb]), neg(graph, im))
    } else if op == ops.sinh {
        let sa = call(graph, ops.sinh, &[a]);
        let cb = call(graph, ops.cos, &[b]);
        let ca = call(graph, ops.cosh, &[a]);
        let sb = call(graph, ops.sin, &[b]);
        (mul(graph, &[sa, cb]), mul(graph, &[ca, sb]))
    } else {
        let ca = call(graph, ops.cosh, &[a]);
        let cb = call(graph, ops.cos, &[b]);
        let sa = call(graph, ops.sinh, &[a]);
        let sb = call(graph, ops.sin, &[b]);
        (mul(graph, &[ca, cb]), mul(graph, &[sa, sb]))
    }
}

/// `(a + bi)(c + di)`.
fn product(
    graph: &mut Graph,
    (a, b): (NodeId, NodeId),
    (c, d): (NodeId, NodeId),
) -> (NodeId, NodeId) {
    let ac = mul(graph, &[a, c]);
    let bd = mul(graph, &[b, d]);
    let ad = mul(graph, &[a, d]);
    let bc = mul(graph, &[b, c]);
    (sub(graph, ac, bd), add(graph, &[ad, bc]))
}

/// `r (cos t + i sin t)`.
fn polar(
    graph: &mut Graph,
    ops: Ops,
    r: NodeId,
    t: NodeId,
) -> (NodeId, NodeId) {
    let cos = call(graph, ops.cos, &[t]);
    let sin = call(graph, ops.sin, &[t]);
    (mul(graph, &[r, cos]), mul(graph, &[r, sin]))
}

/// `sqrt(a^2 + b^2)`.
fn modulus(
    graph: &mut Graph,
    a: NodeId,
    b: NodeId,
) -> NodeId {
    let a2 = build::powi(graph, a, 2);
    let b2 = build::powi(graph, b, 2);
    let sum = add(graph, &[a2, b2]);
    let half = graph.num(Number::fraction(1, 2).unwrap_or_else(|| Number::from(0)));
    build::pow(graph, sum, half)
}

/// Whether `node` contains one of `ops` anywhere.
fn mentions(
    graph: &Graph,
    node: NodeId,
    ops: &[OpId],
) -> bool {
    let mut stack = vec![node];
    let mut seen = std::collections::HashSet::new();
    while let Some(n) = stack.pop() {
        if !seen.insert(n) {
            continue;
        }
        if ops.contains(&graph.op(n)) {
            return true;
        }
        stack.extend_from_slice(graph.children(n));
    }
    false
}

/// `I^n` for integer `n`.
struct UnitPower {
    unit: OpId,
}

impl Kernel for UnitPower {
    fn ops(&self) -> Vec<OpId> {
        vec![core::POW]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let graph = &mut *cx.graph;
        let &[base, exponent] = graph.children(node) else {
            return Outcome::Pass;
        };
        if !graph
            .enodes(graph.find(base))
            .any(|n| graph.op(n) == self.unit)
        {
            return Outcome::Pass;
        }
        let Some(n) = graph.number_of(exponent).and_then(Number::to_i64) else {
            return Outcome::Pass;
        };
        let i = graph.node(self.unit, &[]);
        // Pinned: `I^3` and `-I` cost the same, and the power must not win.
        Outcome::Pinned(match n.rem_euclid(4) {
            | 0 => graph.int(1),
            | 1 => i,
            | 2 => graph.int(-1),
            | _ => neg(graph, i),
        })
    }

    fn revisit(&self) -> bool {
        true
    }
}

/// `re`, `im`, `conj`, `abs` and `arg` of terms whose parts are known.
struct Parts {
    ops: Ops,
    re: OpId,
    im: OpId,
    conj: OpId,
    arg: OpId,
}

impl Kernel for Parts {
    fn ops(&self) -> Vec<OpId> {
        vec![self.re, self.im, self.conj, self.arg, self.ops.abs]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let op = cx.graph.op(node);
        let &[z] = cx.graph.children(node) else {
            return Outcome::Pass;
        };
        let Some(term) = best(cx.graph, z) else {
            return Outcome::Pass;
        };
        // Symbols of unknown reality have no parts to read off.
        if cx.graph.as_symbol(term).is_some() {
            return Outcome::Pass;
        }
        if op == self.ops.abs && cx.graph.facts(term).has(Facts::REAL) {
            return Outcome::Pass;
        }
        let Some((a, b)) = split(cx.graph, self.ops, term) else {
            return Outcome::Pass;
        };
        let (a, b) = (cx.simplify(a), cx.simplify(b));
        let graph = &mut *cx.graph;
        let result = if op == self.re {
            a
        } else if op == self.im {
            b
        } else if op == self.conj {
            let minus_b = neg(graph, b);
            build::complex(graph, self.ops.unit, a, minus_b)
        } else if op == self.arg {
            // On the imaginary axis atan2 has no closed form here.
            let a_zero = graph.number_of(a).is_some_and(Number::is_zero);
            let b_sign = graph.facts(b);
            match (
                a_zero,
                b_sign.has(Facts::POSITIVE),
                b_sign.has(Facts::NEGATIVE),
            ) {
                | (true, true, _) | (true, _, true) => {
                    let pi = graph.ops().lookup("pi").map(|p| graph.node(p, &[]));
                    let Some(pi) = pi else {
                        return Outcome::Pass;
                    };
                    let half = graph.num(
                        Number::fraction(
                            if b_sign.has(Facts::POSITIVE) {
                                1
                            } else {
                                -1
                            },
                            2,
                        )
                        .unwrap_or_else(|| Number::from(0)),
                    );
                    mul(graph, &[half, pi])
                },
                | _ => call(graph, self.ops.atan2, &[b, a]),
            }
        } else {
            modulus(graph, a, b)
        };
        // Splitting a term that already is `re`/`im` of something gains
        // nothing and would only cycle.
        if mentions(graph, result, &[op]) && mentions(graph, term, &[self.ops.re, self.ops.im]) {
            return Outcome::Pass;
        }
        Outcome::Equal(result)
    }
}

struct ComplexExpand {
    op: OpId,
    ops: Ops,
}

impl Kernel for ComplexExpand {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let &[z] = cx.graph.children(node) else {
            return Outcome::Pass;
        };
        let Some((a, b)) = split(cx.graph, self.ops, z) else {
            return Outcome::Pass;
        };
        let (a, b) = (cx.simplify(a), cx.simplify(b));
        let expanded = build::complex(cx.graph, self.ops.unit, a, b);
        Outcome::Equal(expanded)
    }
}

/// `exp(z) = exp(a) (cos b + i sin b)` when the trigonometric values have
/// closed forms (`exp(I*pi) = -1`, `exp(I*pi/2) = I`).
struct ExpI {
    ops: Ops,
}

impl Kernel for ExpI {
    fn ops(&self) -> Vec<OpId> {
        vec![self.ops.exp]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let &[w] = cx.graph.children(node) else {
            return Outcome::Pass;
        };
        let Some(term) = best(cx.graph, w) else {
            return Outcome::Pass;
        };
        if !mentions(cx.graph, term, &[self.ops.unit]) {
            return Outcome::Pass;
        }
        let Some((a, b)) = split(cx.graph, self.ops, term) else {
            return Outcome::Pass;
        };
        let b = cx.simplify(b);
        let cos = call(cx.graph, self.ops.cos, &[b]);
        let sin = call(cx.graph, self.ops.sin, &[b]);
        let (cos, sin) = (cx.simplify(cos), cx.simplify(sin));
        if mentions(cx.graph, cos, &[self.ops.cos, self.ops.sin])
            || mentions(cx.graph, sin, &[self.ops.cos, self.ops.sin])
        {
            return Outcome::Pass;
        }
        let graph = &mut *cx.graph;
        let rotation = build::complex(graph, self.ops.unit, cos, sin);
        let a_is_zero = graph.number_of(a).is_some_and(Number::is_zero);
        if a_is_zero {
            return Outcome::Equal(rotation);
        }
        let scale = call(graph, self.ops.exp, &[a]);
        Outcome::Equal(mul(graph, &[scale, rotation]))
    }
}

/// Registers [`ComplexEval`] for the elementary functions and the parts
/// operators.
fn register_complex_evals(graph: &mut Graph) {
    macro_rules! unary {
        ($f:expr) => {
            ComplexEval(|a: &[Complex64]| a.first().map($f))
        };
    }
    let table: [(&str, ComplexEval); 29] = [
        ("I", ComplexEval(|_| Some(Complex64::new(0.0, 1.0)))),
        ("exp", unary!(|z: &Complex64| z.exp())),
        ("ln", unary!(|z: &Complex64| z.ln())),
        ("sqrt", unary!(|z: &Complex64| z.sqrt())),
        ("sin", unary!(|z: &Complex64| z.sin())),
        ("cos", unary!(|z: &Complex64| z.cos())),
        ("tan", unary!(|z: &Complex64| z.tan())),
        ("sinh", unary!(|z: &Complex64| z.sinh())),
        ("cosh", unary!(|z: &Complex64| z.cosh())),
        ("tanh", unary!(|z: &Complex64| z.tanh())),
        ("asin", unary!(|z: &Complex64| z.asin())),
        ("acos", unary!(|z: &Complex64| z.acos())),
        ("atan", unary!(|z: &Complex64| z.atan())),
        ("asinh", unary!(|z: &Complex64| z.asinh())),
        ("acosh", unary!(|z: &Complex64| z.acosh())),
        ("atanh", unary!(|z: &Complex64| z.atanh())),
        ("sec", unary!(|z: &Complex64| z.cos().inv())),
        ("csc", unary!(|z: &Complex64| z.sin().inv())),
        ("cot", unary!(|z: &Complex64| z.tan().inv())),
        ("sech", unary!(|z: &Complex64| z.cosh().inv())),
        ("csch", unary!(|z: &Complex64| z.sinh().inv())),
        ("coth", unary!(|z: &Complex64| z.tanh().inv())),
        ("abs", unary!(|z: &Complex64| Complex64::new(z.norm(), 0.0))),
        ("re", unary!(|z: &Complex64| Complex64::new(z.re, 0.0))),
        ("im", unary!(|z: &Complex64| Complex64::new(z.im, 0.0))),
        ("conj", unary!(|z: &Complex64| z.conj())),
        ("arg", unary!(|z: &Complex64| Complex64::new(z.arg(), 0.0))),
        (
            "pi",
            ComplexEval(|_| Some(Complex64::new(std::f64::consts::PI, 0.0))),
        ),
        (
            "E",
            ComplexEval(|_| Some(Complex64::new(std::f64::consts::E, 0.0))),
        ),
    ];
    for (name, eval) in table {
        if let Some(op) = graph.ops().lookup(name) {
            graph.ops_mut().set_attr(op, eval);
        }
    }
}

/// Evaluates the concrete term `node` over the complex numbers; see
/// [`Graph::eval_complex`].
#[must_use]
#[allow(clippy::implicit_hasher)] // delegates to `Graph::eval_complex`, which fixes the default hasher
pub fn eval_complex(
    graph: &Graph,
    node: NodeId,
    bindings: &HashMap<SymbolId, Complex64>,
) -> Option<Complex64> {
    graph.eval_complex(node, bindings)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Engine;
    use crate::rules::testing::simplify;

    fn run(src: &str) -> String {
        simplify(&[complex()], src)
    }

    #[test]
    fn unit_arithmetic() {
        assert_eq!(run("I^2"), "-1");
        assert_eq!(run("I*I*I"), "-I");
        assert_eq!(run("I^7"), "-I");
        assert_eq!(run("re((1 + I)*(1 - I))"), "2");
        assert_eq!(run("re((2 + 3*I)*(1 - I))"), "5");
        assert_eq!(run("im((2 + 3*I)*(1 - I))"), "1");
        assert_eq!(run("im(1/(1 + I))"), "-1/2");
        assert_eq!(run("conj(2 + 3*I)"), "2 - 3*I");
        assert_eq!(run("abs(3 + 4*I)"), "5");
        assert_eq!(run("arg(I)"), "1/2*pi");
        assert_eq!(run("arg(-2)"), "pi");
    }

    #[test]
    fn integer_power_outside_the_split_bound_stays_symbolic() {
        let got = run("re((1 + I)^(-9223372036854775808))");
        assert!(
            got.contains("re(") && got.contains("9223372036854775808"),
            "{got}"
        );
    }

    #[test]
    fn euler() {
        assert_eq!(run("exp(I*pi)"), "-1");
        assert_eq!(run("exp(I*pi/2)"), "I");
        assert_eq!(run("exp(2*pi*I)"), "1");
        assert_eq!(run("re(exp(I*x))"), "re(exp(x*I))", "x might not be real");
    }

    #[test]
    fn parts_of_real_symbols() {
        let mut g = Graph::new();
        let engine = Engine::install(&mut g, &[complex()]).unwrap_or_else(|e| panic!("{e}"));
        let x = g.interner_mut().symbol("x");
        let y = g.interner_mut().symbol("y");
        g.assume(x, Facts::REAL);
        g.assume(y, Facts::REAL);
        let cases = [
            ("re(exp(I*x))", "cos(x)"),
            ("im(exp(x + I*y))", "exp(x)*sin(y)"),
            ("re(x + I*y)", "x"),
            ("conj(x + I*y)", "x - y*I"),
            ("im(sin(x + I*y))", "cos(x)*sinh(y)"),
            ("abs(exp(I*x))", "1"),
        ];
        for (src, want) in cases {
            let root = g.parse(src).unwrap_or_else(|e| panic!("{e}"));
            engine.run(
                &mut g,
                &[root],
                &crate::graph::Env::symbolic(),
                &crate::graph::Saturate,
                &crate::graph::Budget::default(),
            );
            let got = crate::graph::Extractor::new(&g, &[root], &crate::graph::SizeCost)
                .build(&mut g, root)
                .map(|n| g.display(n));
            assert_eq!(got.as_deref(), Some(want), "{src}");
        }
    }

    #[test]
    fn complex_evaluation() {
        let mut g = Graph::new();
        assert!(Engine::install(&mut g, &[complex()]).is_ok());
        let root = g
            .parse("exp(I*z) * (1 + I)^2 + ln(z)")
            .unwrap_or_else(|e| panic!("{e}"));
        let z = g.interner_mut().symbol("z");
        let at = Complex64::new(0.3, -0.4);
        let bindings: HashMap<SymbolId, Complex64> = std::iter::once((z, at)).collect();
        let got = eval_complex(&g, root, &bindings).unwrap_or_default();
        let i = Complex64::new(0.0, 1.0);
        let want = (i * at).exp() * (1.0 + i) * (1.0 + i) + at.ln();
        assert!((got - want).norm() < 1e-12, "{got} vs {want}");
    }

    /// Splitting is an identity: both sides agree at random complex points
    /// for terms whose parts are computed.
    #[test]
    fn split_agrees_with_evaluation() {
        let mut g = Graph::new();
        assert!(Engine::install(&mut g, &[complex()]).is_ok());
        let ops = Ops::of(&g).unwrap_or_else(|| panic!("ops"));
        for src in [
            "exp(2 + 3*I)",
            "sin(1 - 2*I)",
            "cosh(0.5 + I)",
            "ln(-1 + I)",
            "(2 + I)^-3",
            "2^(1 + I)",
            "cos(I)*sinh(1 + I)",
        ] {
            let root = g.parse(src).unwrap_or_else(|e| panic!("{e}"));
            let (a, b) = split(&mut g, ops, root).unwrap_or_else(|| panic!("no split for {src}"));
            let none = HashMap::new();
            let whole = eval_complex(&g, root, &none).unwrap_or_default();
            let a = eval_complex(&g, a, &none).unwrap_or_default();
            let b = eval_complex(&g, b, &none).unwrap_or_default();
            assert!(a.im.abs() < 1e-12 && b.im.abs() < 1e-12, "{src}");
            assert!(
                (Complex64::new(a.re, b.re) - whole).norm() < 1e-12,
                "{src}: {whole} vs {a} + i{b}"
            );
        }
    }
}
