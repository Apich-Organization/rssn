//! Geometric (Clifford) algebra with symbolic coefficients.
//!
//! A multivector over the metric signature `(p, q, r)` — `p` basis vectors
//! squaring to `+1`, then `q` squaring to `-1`, then `r` squaring to `0` —
//! is the inert term
//!
//! ```text
//! mv(list(p, q, r), list(list(blade, coefficient), ...))
//! ```
//!
//! where `blade` is the bitmask of the basis vectors in the blade: bit `i`
//! stands for `e_{i+1}`, so `1 = e1`, `2 = e2`, `3 = e12`, `4 = e3`,
//! `5 = e13`, `6 = e23`, `7 = e123`, and `0` is the scalar. The blade
//! `e_{i1} e_{i2} ...` has the basis vectors in increasing order. The
//! canonical form has sorted blades, simplified coefficients and no zero
//! coefficients; the `mv` operator normalises itself into it. Coefficients
//! are arbitrary terms.
//!
//! A bare (non-`mv`) operand of a binary operation is a scalar of the
//! other operand's algebra.
//!
//! | operator | value |
//! |---|---|
//! | `mv_scalar(sig, s)`, `mv_vector(sig, list(x1, ..))`, `mv_blade(sig, mask, c)`, `mv_basis(sig, i)` (`e_i`, from 1), `mv_pseudoscalar(sig)` | constructors |
//! | `mv_add(A, B, ..)`, `mv_sub(A, B)`, `mv_scale(c, A)` | linear structure |
//! | `mv_coeff(A, mask)`, `mv_scalar_part(A)` | a coefficient (zero when absent) |
//! | `ga_gp(A, B)` | geometric product |
//! | `ga_wedge(A, B)` | outer product, `<A_r B_s>_{r+s}` |
//! | `ga_lcontract(A, B)`, `ga_rcontract(A, B)` | left contraction `<A_r B_s>_{s-r}`, right contraction `<A_r B_s>_{r-s}` |
//! | `ga_inner(A, B)` | Hestenes inner product `<A_r B_s>_{|r-s|}` (zero for scalar operands) |
//! | `ga_scalar_product(A, B)`, `ga_commutator(A, B)` | `<AB>_0` and `(AB - BA)/2` |
//! | `ga_grade(A, k)` | grade projection |
//! | `ga_reverse`, `ga_involute`, `ga_conjugate` | reverse `(-1)^{k(k-1)/2}`, grade involution `(-1)^k`, Clifford conjugate `(-1)^{k(k+1)/2}` |
//! | `ga_dual(A)`, `ga_undual(A)` | `A I⁻¹` and `A I`, `I` the pseudoscalar (for degenerate metrics, where `I` has no inverse, `ga_dual` is `A I`) |
//! | `ga_norm_squared(A)`, `ga_norm(A)` | `<A Ã>_0` and its square root |
//! | `ga_normalize(A)` | `A / ga_norm(A)` |
//! | `ga_inverse(A)` | `Ã / (A Ã)` for a versor, `A / A²` when `A²` is a scalar |
//! | `ga_exp(A)` | `exp(A)` when `A²` is a scalar |
//! | `ga_rotor(B, angle)` | the rotor `exp(-angle B̂ / 2)` of the plane of the bivector `B` (`B̂ = B/√(-B²)`); `R x R̃` rotates by `angle` from the first to the second basis vector of `B` |
//! | `ga_apply_versor(V, X)` | the sandwich `V X̂ V⁻¹` (`X̂` the grade involution for an odd versor, none for an even one) |

use std::collections::BTreeMap;

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
use crate::rules::poly::best;

use super::elementary::elementary;
use super::lie::add;
use super::lie::call;
use super::lie::exact_sqrt;
use super::lie::exp_even;
use super::lie::is_zero_literal;
use super::lie::mul;
use super::lie::pow_int;

/// The geometric algebra rule set.
#[must_use]
pub fn geometric_algebra() -> RuleSet {
    RuleSet::new("geometric_algebra", install).needs(elementary())
}

/// Largest number of basis vectors.
const MAX_DIM: u32 = 12;

type Sig = (u32, u32, u32);

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum Request {
    Scalar,
    Vector,
    Blade,
    Basis,
    Pseudoscalar,
    Add,
    Sub,
    Scale,
    Coeff,
    ScalarPart,
    Product(Product),
    Grade,
    Reverse,
    Involute,
    Conjugate,
    Dual,
    Undual,
    NormSquared,
    Norm,
    Normalize,
    Inverse,
    Exp,
    Rotor,
    Apply,
}

/// Which grades of the product of two blades are kept.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum Product {
    Geometric,
    Outer,
    Left,
    Right,
    Inner,
    Scalar,
    Commutator,
}

fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let mv = i.op(OpDescriptor::new("mv", Arity::Fixed(2)))?;
    i.kernel("geometric_algebra/mv", Tier::Normalize, Normalise { mv });
    for (name, arity, request) in [
        ("mv_scalar", Arity::Fixed(2), Request::Scalar),
        ("mv_vector", Arity::Fixed(2), Request::Vector),
        ("mv_blade", Arity::Fixed(3), Request::Blade),
        ("mv_basis", Arity::Fixed(2), Request::Basis),
        ("mv_pseudoscalar", Arity::Fixed(1), Request::Pseudoscalar),
        ("mv_add", Arity::Variadic, Request::Add),
        ("mv_sub", Arity::Fixed(2), Request::Sub),
        ("mv_scale", Arity::Fixed(2), Request::Scale),
        ("mv_coeff", Arity::Fixed(2), Request::Coeff),
        ("mv_scalar_part", Arity::Fixed(1), Request::ScalarPart),
        ("ga_gp", Arity::Fixed(2), Request::Product(Product::Geometric)),
        ("ga_wedge", Arity::Fixed(2), Request::Product(Product::Outer)),
        ("ga_lcontract", Arity::Fixed(2), Request::Product(Product::Left)),
        ("ga_rcontract", Arity::Fixed(2), Request::Product(Product::Right)),
        ("ga_inner", Arity::Fixed(2), Request::Product(Product::Inner)),
        ("ga_scalar_product", Arity::Fixed(2), Request::Product(Product::Scalar)),
        ("ga_commutator", Arity::Fixed(2), Request::Product(Product::Commutator)),
        ("ga_grade", Arity::Fixed(2), Request::Grade),
        ("ga_reverse", Arity::Fixed(1), Request::Reverse),
        ("ga_involute", Arity::Fixed(1), Request::Involute),
        ("ga_conjugate", Arity::Fixed(1), Request::Conjugate),
        ("ga_dual", Arity::Fixed(1), Request::Dual),
        ("ga_undual", Arity::Fixed(1), Request::Undual),
        ("ga_norm_squared", Arity::Fixed(1), Request::NormSquared),
        ("ga_norm", Arity::Fixed(1), Request::Norm),
        ("ga_normalize", Arity::Fixed(1), Request::Normalize),
        ("ga_inverse", Arity::Fixed(1), Request::Inverse),
        ("ga_exp", Arity::Fixed(1), Request::Exp),
        ("ga_rotor", Arity::Fixed(2), Request::Rotor),
        ("ga_apply_versor", Arity::Fixed(2), Request::Apply),
    ] {
        let op = i.op(OpDescriptor::new(name, arity).flags(OpFlags::HEAVY).cost(100))?;
        i.kernel(&format!("geometric_algebra/{name}"), Tier::Reduce, Ga { op, mv, request });
    }
    Ok(())
}

// ----------------------------------------------------------------------
// Multivectors
// ----------------------------------------------------------------------

/// A multivector: blade bitmask to coefficient.
#[derive(Clone, Debug)]
struct Mv {
    sig: Sig,
    terms: BTreeMap<u32, NodeId>,
}

/// Sign and metric factor (`0` when a degenerate vector is contracted) of
/// the product of two basis blades, and the resulting blade.
fn blade_product(
    sig: Sig,
    a: u32,
    b: u32,
) -> (i64, u32) {
    let mut sign = 1;
    for i in 0..32 {
        if (b >> i) & 1 == 1 && (a >> (i + 1)).count_ones() % 2 == 1 {
            sign = -sign;
        }
    }
    for i in 0..32 {
        if (a & b) >> i & 1 == 1 {
            let (p, q, _) = sig;
            if i >= p + q {
                return (0, a ^ b);
            }
            if i >= p {
                sign = -sign;
            }
        }
    }
    (sign, a ^ b)
}

fn read_sig(
    graph: &mut Graph,
    node: NodeId,
) -> Option<Sig> {
    let node = best(graph, node)?;
    if graph.op(node) != core::LIST {
        return None;
    }
    let kids = graph.children(node).to_vec();
    let &[p, q, r] = kids.as_slice() else {
        return None;
    };
    let get = |g: &Graph, n: NodeId| g.number_of(n).and_then(Number::to_i64).and_then(|v| u32::try_from(v).ok());
    let sig = (get(graph, p)?, get(graph, q)?, get(graph, r)?);
    (sig.0 + sig.1 + sig.2 <= MAX_DIM).then_some(sig)
}

fn sig_term(
    graph: &mut Graph,
    sig: Sig,
) -> NodeId {
    let parts: Vec<NodeId> = [sig.0, sig.1, sig.2].iter().map(|&v| graph.int(i64::from(v))).collect();
    graph.node(core::LIST, &parts)
}

/// Reads an `mv` term (merging repeated blades).
fn read_mv(
    cx: &mut Cx<'_>,
    mv_op: OpId,
    node: NodeId,
) -> Option<Mv> {
    let node = best(cx.graph, node)?;
    if cx.graph.op(node) != mv_op {
        return None;
    }
    let &[sig, terms] = cx.graph.children(node) else {
        return None;
    };
    let sig = read_sig(cx.graph, sig)?;
    let terms = best(cx.graph, terms)?;
    if cx.graph.op(terms) != core::LIST {
        return None;
    }
    let mut acc: BTreeMap<u32, Vec<NodeId>> = BTreeMap::new();
    for pair in cx.graph.children(terms).to_vec() {
        let pair = best(cx.graph, pair)?;
        let &[mask, coefficient] = cx.graph.children(pair) else {
            return None;
        };
        let mask = u32::try_from(cx.graph.number_of(mask)?.to_i64()?).ok()?;
        if mask >> (sig.0 + sig.1 + sig.2) != 0 {
            return None;
        }
        acc.entry(mask).or_default().push(coefficient);
    }
    Some(finish(cx, sig, acc))
}

/// Sums and simplifies the coefficients and drops the zero ones.
fn finish(
    cx: &mut Cx<'_>,
    sig: Sig,
    acc: BTreeMap<u32, Vec<NodeId>>,
) -> Mv {
    let mut terms = BTreeMap::new();
    for (blade, parts) in acc {
        let sum = add(cx.graph, &parts);
        let sum = cx.simplify(sum);
        if !is_zero_literal(cx.graph, sum) {
            terms.insert(blade, sum);
        }
    }
    Mv { sig, terms }
}

fn mv_term(
    cx: &mut Cx<'_>,
    mv_op: OpId,
    m: &Mv,
) -> NodeId {
    let pairs: Vec<NodeId> = m
        .terms
        .iter()
        .map(|(&blade, &c)| {
            let b = cx.graph.int(i64::from(blade));
            cx.graph.node(core::LIST, &[b, c])
        })
        .collect();
    let list = cx.graph.node(core::LIST, &pairs);
    let sig = sig_term(cx.graph, m.sig);
    cx.graph.node(mv_op, &[sig, list])
}

fn scalar_mv(
    sig: Sig,
    c: NodeId,
) -> Mv {
    Mv { sig, terms: BTreeMap::from([(0, c)]) }
}

/// Scales every coefficient of `m` by `factor`.
fn scaled(
    cx: &mut Cx<'_>,
    m: &Mv,
    factor: NodeId,
) -> Mv {
    let acc = m.terms.iter().map(|(&b, &c)| (b, vec![mul(cx.graph, &[factor, c])])).collect();
    finish(cx, m.sig, acc)
}

fn sum(
    cx: &mut Cx<'_>,
    parts: &[&Mv],
) -> Option<Mv> {
    let sig = parts.first()?.sig;
    if parts.iter().any(|m| m.sig != sig) {
        return None;
    }
    let mut acc: BTreeMap<u32, Vec<NodeId>> = BTreeMap::new();
    for m in parts {
        for (&b, &c) in &m.terms {
            acc.entry(b).or_default().push(c);
        }
    }
    Some(finish(cx, sig, acc))
}

fn product(
    cx: &mut Cx<'_>,
    a: &Mv,
    b: &Mv,
    keep: impl Fn(u32, u32, u32) -> bool,
) -> Option<Mv> {
    if a.sig != b.sig {
        return None;
    }
    let mut acc: BTreeMap<u32, Vec<NodeId>> = BTreeMap::new();
    for (&ba, &ca) in &a.terms {
        for (&bb, &cb) in &b.terms {
            let (sign, blade) = blade_product(a.sig, ba, bb);
            if sign == 0 || !keep(ba.count_ones(), bb.count_ones(), blade.count_ones()) {
                continue;
            }
            let s = cx.graph.int(sign);
            acc.entry(blade).or_default().push(mul(cx.graph, &[s, ca, cb]));
        }
    }
    Some(finish(cx, a.sig, acc))
}

fn gp(
    cx: &mut Cx<'_>,
    a: &Mv,
    b: &Mv,
) -> Option<Mv> {
    product(cx, a, b, |_, _, _| true)
}

/// Applies `sign(grade)` to every blade.
fn map_grades(
    cx: &mut Cx<'_>,
    m: &Mv,
    sign: impl Fn(u32) -> i64,
) -> Mv {
    let acc = m
        .terms
        .iter()
        .map(|(&b, &c)| {
            let s = cx.graph.int(sign(b.count_ones()));
            (b, vec![mul(cx.graph, &[s, c])])
        })
        .collect();
    finish(cx, m.sig, acc)
}

fn reverse(
    cx: &mut Cx<'_>,
    m: &Mv,
) -> Mv {
    map_grades(cx, m, |k| if (k * k.saturating_sub(1) / 2) % 2 == 0 { 1 } else { -1 })
}

fn involute(
    cx: &mut Cx<'_>,
    m: &Mv,
) -> Mv {
    map_grades(cx, m, |k| if k % 2 == 0 { 1 } else { -1 })
}

/// The coefficient of the scalar part.
fn scalar_part(
    cx: &mut Cx<'_>,
    m: &Mv,
) -> NodeId {
    m.terms.get(&0).copied().unwrap_or_else(|| cx.graph.int(0))
}

/// The scalar value of `m` when it has no other part.
fn as_scalar(
    cx: &mut Cx<'_>,
    m: &Mv,
) -> Option<NodeId> {
    m.terms.keys().all(|&b| b == 0).then(|| scalar_part(cx, m))
}

fn half_power(
    cx: &mut Cx<'_>,
    x: NodeId,
) -> Option<NodeId> {
    let half = cx.graph.num(Number::fraction(1, 2)?);
    let r = cx.graph.node(core::POW, &[x, half]);
    Some(cx.simplify(r))
}

fn inverse(
    cx: &mut Cx<'_>,
    a: &Mv,
) -> Option<Mv> {
    let reversed = reverse(cx, a);
    let aa = gp(cx, a, &reversed)?;
    let (numerator, norm) = if let Some(s) = as_scalar(cx, &aa) {
        (reversed, s)
    } else {
        let sq = gp(cx, a, a)?;
        (a.clone(), as_scalar(cx, &sq)?)
    };
    if is_zero_literal(cx.graph, norm) {
        return None;
    }
    let inv = pow_int(cx.graph, norm, -1);
    Some(scaled(cx, &numerator, inv))
}

fn parity(m: &Mv) -> Option<bool> {
    let mut odd = m.terms.keys().map(|b| b.count_ones() % 2 == 1);
    let first = odd.next()?;
    odd.all(|p| p == first).then_some(first)
}

fn rotor(
    cx: &mut Cx<'_>,
    b: &Mv,
    angle: NodeId,
) -> Option<Mv> {
    if b.terms.is_empty() || b.terms.keys().any(|k| k.count_ones() != 2) {
        return None;
    }
    let bb = gp(cx, b, b)?;
    let s = as_scalar(cx, &bb)?;
    let minus_one = cx.graph.int(-1);
    let neg_s = mul(cx.graph, &[minus_one, s]);
    let neg_s = cx.simplify(neg_s);
    let half = cx.graph.num(Number::fraction(1, 2)?);
    let t = mul(cx.graph, &[half, angle]);
    let t = cx.simplify(t);
    let euclidean = b.sig.1 == 0 && b.sig.2 == 0;
    let sign_known = cx.graph.number_of(s).map(Number::is_negative);
    let (trig, w) = if let Some(w) = exact_sqrt(cx.graph, neg_s).filter(|_| sign_known != Some(false)) {
        (true, w)
    } else if let Some(w) = exact_sqrt(cx.graph, s) {
        (false, w)
    } else if euclidean || sign_known == Some(true) {
        (true, half_power(cx, neg_s)?)
    } else if sign_known == Some(false) && !is_zero_literal(cx.graph, s) {
        (false, half_power(cx, s)?)
    } else if is_zero_literal(cx.graph, s) {
        let c = cx.graph.int(1);
        let k = mul(cx.graph, &[minus_one, t]);
        let k = cx.simplify(k);
        let scalar = scalar_mv(b.sig, c);
        let rest = scaled(cx, b, k);
        return sum(cx, &[&scalar, &rest]);
    } else {
        return None;
    };
    let (cos, sin) = if trig { ("cos", "sin") } else { ("cosh", "sinh") };
    let c = call(cx.graph, cos, &[t])?;
    let sn = call(cx.graph, sin, &[t])?;
    let w_inv = pow_int(cx.graph, w, -1);
    let k = mul(cx.graph, &[minus_one, sn, w_inv]);
    let k = cx.simplify(k);
    let scalar = scalar_mv(b.sig, cx.simplify(c));
    let rest = scaled(cx, b, k);
    sum(cx, &[&scalar, &rest])
}

const fn pseudoscalar(sig: Sig) -> u32 {
    (1 << (sig.0 + sig.1 + sig.2)) - 1
}

fn dual(
    cx: &mut Cx<'_>,
    a: &Mv,
    inverse: bool,
) -> Option<Mv> {
    let blade = pseudoscalar(a.sig);
    let one = cx.graph.int(1);
    let i = Mv { sig: a.sig, terms: BTreeMap::from([(blade, one)]) };
    let r = gp(cx, a, &i)?;
    if !inverse || a.sig.2 > 0 {
        return Some(r);
    }
    // I² = ±1 for a non-degenerate metric, so I⁻¹ = I / I².
    let (sign, _) = blade_product(a.sig, blade, blade);
    let s = cx.graph.int(sign);
    Some(scaled(cx, &r, s))
}

// ----------------------------------------------------------------------
// Kernels
// ----------------------------------------------------------------------

/// Brings a user-built `mv` into canonical form.
struct Normalise {
    mv: OpId,
}

impl Kernel for Normalise {
    fn ops(&self) -> Vec<OpId> {
        vec![self.mv]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let Some(m) = read_mv(cx, self.mv, node) else {
            return Outcome::Pass;
        };
        let canonical = mv_term(cx, self.mv, &m);
        if cx.graph.same(node, canonical) {
            Outcome::Pass
        } else {
            Outcome::Equal(canonical)
        }
    }
}

struct Ga {
    op: OpId,
    mv: OpId,
    request: Request,
}

impl Kernel for Ga {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let args = cx.graph.children(node).to_vec();
        match self.compute(cx, &args) {
            | Some(m) => Outcome::Equal(match m {
                | Out::Mv(m) => mv_term(cx, self.mv, &m),
                | Out::Term(t) => t,
            }),
            | None => Outcome::Pass,
        }
    }

    fn revisit(&self) -> bool {
        true
    }
}

enum Out {
    Mv(Mv),
    Term(NodeId),
}

fn int_arg(
    graph: &Graph,
    node: NodeId,
) -> Option<u32> {
    u32::try_from(graph.number_of(node)?.to_i64()?).ok()
}

impl Ga {
    /// Reads an operand; a non-`mv` term is a scalar of algebra `sig`.
    fn operand(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
        sig: Option<Sig>,
    ) -> Option<Mv> {
        if let Some(m) = read_mv(cx, self.mv, node) {
            return Some(m);
        }
        let c = Self::scalar(cx, node)?;
        let sig = sig?;
        Some(finish(cx, sig, BTreeMap::from([(0, vec![c])])))
    }

    /// A scalar coefficient: a term that is neither a list nor (still
    /// to be) a multivector.
    fn scalar(
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Option<NodeId> {
        let c = best(cx.graph, node)?;
        let mut stack = vec![c];
        while let Some(t) = stack.pop() {
            let op = cx.graph.op(t);
            if t == c && op == core::LIST {
                return None;
            }
            if cx.graph.ops().get(op).name.starts_with("mv") || cx.graph.ops().get(op).name.starts_with("ga_") {
                return None;
            }
            stack.extend_from_slice(cx.graph.children(t));
        }
        Some(c)
    }

    /// Both operands of a binary operation.
    fn pair(
        &self,
        cx: &mut Cx<'_>,
        a: NodeId,
        b: NodeId,
    ) -> Option<(Mv, Mv)> {
        let (x, y) = (read_mv(cx, self.mv, a), read_mv(cx, self.mv, b));
        let sig = x.as_ref().or(y.as_ref())?.sig;
        let x = match x {
            | Some(x) => x,
            | None => self.operand(cx, a, Some(sig))?,
        };
        let y = match y {
            | Some(y) => y,
            | None => self.operand(cx, b, Some(sig))?,
        };
        Some((x, y))
    }

    #[allow(clippy::too_many_lines)]
    fn compute(
        &self,
        cx: &mut Cx<'_>,
        args: &[NodeId],
    ) -> Option<Out> {
        let one = |cx: &mut Cx<'_>| -> Option<Mv> { self.operand(cx, *args.first()?, None) };
        match self.request {
            | Request::Scalar => {
                let sig = read_sig(cx.graph, args[0])?;
                let c = Self::scalar(cx, args[1])?;
                Some(Out::Mv(finish(cx, sig, BTreeMap::from([(0, vec![c])]))))
            },
            | Request::Vector => {
                let sig = read_sig(cx.graph, args[0])?;
                let v = best(cx.graph, args[1])?;
                if cx.graph.op(v) != core::LIST || cx.graph.children(v).len() > (sig.0 + sig.1 + sig.2) as usize {
                    return None;
                }
                let acc = cx.graph.children(v).to_vec().into_iter().enumerate().map(|(k, c)| (1_u32 << k, vec![c])).collect();
                Some(Out::Mv(finish(cx, sig, acc)))
            },
            | Request::Blade => {
                let sig = read_sig(cx.graph, args[0])?;
                let mask = int_arg(cx.graph, args[1])?;
                if mask >> (sig.0 + sig.1 + sig.2) != 0 {
                    return None;
                }
                Some(Out::Mv(finish(cx, sig, BTreeMap::from([(mask, vec![args[2]])]))))
            },
            | Request::Basis => {
                let sig = read_sig(cx.graph, args[0])?;
                let i = int_arg(cx.graph, args[1])?;
                if i == 0 || i > sig.0 + sig.1 + sig.2 {
                    return None;
                }
                let c = cx.graph.int(1);
                Some(Out::Mv(finish(cx, sig, BTreeMap::from([(1 << (i - 1), vec![c])]))))
            },
            | Request::Pseudoscalar => {
                let sig = read_sig(cx.graph, args[0])?;
                let c = cx.graph.int(1);
                Some(Out::Mv(finish(cx, sig, BTreeMap::from([(pseudoscalar(sig), vec![c])]))))
            },
            | Request::Add | Request::Sub => {
                let first = args.iter().find_map(|&a| read_mv(cx, self.mv, a))?;
                let mut parts = Vec::new();
                for &a in args {
                    parts.push(self.operand(cx, a, Some(first.sig))?);
                }
                if self.request == Request::Sub {
                    let minus_one = cx.graph.int(-1);
                    let last = parts.pop()?;
                    let negated = scaled(cx, &last, minus_one);
                    parts.push(negated);
                }
                let refs: Vec<&Mv> = parts.iter().collect();
                sum(cx, &refs).map(Out::Mv)
            },
            | Request::Scale => {
                let m = read_mv(cx, self.mv, args[1])?;
                let c = Self::scalar(cx, args[0])?;
                Some(Out::Mv(scaled(cx, &m, c)))
            },
            | Request::Coeff | Request::ScalarPart => {
                let m = one(cx)?;
                let mask = if self.request == Request::Coeff { int_arg(cx.graph, args[1])? } else { 0 };
                Some(Out::Term(m.terms.get(&mask).copied().unwrap_or_else(|| cx.graph.int(0))))
            },
            | Request::Product(kind) => {
                let (a, b) = self.pair(cx, args[0], args[1])?;
                let r = match kind {
                    | Product::Geometric => gp(cx, &a, &b)?,
                    | Product::Outer => product(cx, &a, &b, |ga, gb, gr| gr == ga + gb)?,
                    | Product::Left => product(cx, &a, &b, |ga, gb, gr| gb >= ga && gr == gb - ga)?,
                    | Product::Right => product(cx, &a, &b, |ga, gb, gr| ga >= gb && gr == ga - gb)?,
                    | Product::Inner => product(cx, &a, &b, |ga, gb, gr| ga > 0 && gb > 0 && gr == ga.abs_diff(gb))?,
                    | Product::Scalar => product(cx, &a, &b, |_, _, gr| gr == 0)?,
                    | Product::Commutator => {
                        let ab = gp(cx, &a, &b)?;
                        let ba = gp(cx, &b, &a)?;
                        let half = cx.graph.num(Number::fraction(1, 2)?);
                        let minus_half = cx.graph.num(Number::fraction(-1, 2)?);
                        let x = scaled(cx, &ab, half);
                        let y = scaled(cx, &ba, minus_half);
                        sum(cx, &[&x, &y])?
                    },
                };
                Some(Out::Mv(r))
            },
            | Request::Grade => {
                let m = read_mv(cx, self.mv, args[0])?;
                let k = int_arg(cx.graph, args[1])?;
                let terms = m.terms.iter().filter(|(b, _)| b.count_ones() == k).map(|(&b, &c)| (b, c)).collect();
                Some(Out::Mv(Mv { sig: m.sig, terms }))
            },
            | Request::Reverse => {
                let m = one(cx)?;
                Some(Out::Mv(reverse(cx, &m)))
            },
            | Request::Involute => {
                let m = one(cx)?;
                Some(Out::Mv(involute(cx, &m)))
            },
            | Request::Conjugate => {
                let m = one(cx)?;
                Some(Out::Mv(map_grades(cx, &m, |k| if (k * (k + 1) / 2) % 2 == 0 { 1 } else { -1 })))
            },
            | Request::Dual | Request::Undual => {
                let m = one(cx)?;
                dual(cx, &m, self.request == Request::Dual).map(Out::Mv)
            },
            | Request::NormSquared | Request::Norm | Request::Normalize => {
                let m = one(cx)?;
                let r = reverse(cx, &m);
                let p = product(cx, &m, &r, |_, _, gr| gr == 0)?;
                let ns = scalar_part(cx, &p);
                if self.request == Request::NormSquared {
                    return Some(Out::Term(ns));
                }
                let norm = half_power(cx, ns)?;
                if self.request == Request::Norm {
                    return Some(Out::Term(norm));
                }
                if is_zero_literal(cx.graph, norm) {
                    return None;
                }
                let inv = pow_int(cx.graph, norm, -1);
                Some(Out::Mv(scaled(cx, &m, inv)))
            },
            | Request::Inverse => {
                let m = one(cx)?;
                inverse(cx, &m).map(Out::Mv)
            },
            | Request::Exp => {
                let m = one(cx)?;
                let sq = gp(cx, &m, &m)?;
                let s = as_scalar(cx, &sq)?;
                let (c, k) = exp_even(cx, s)?;
                let scalar = scalar_mv(m.sig, c);
                let rest = scaled(cx, &m, k);
                sum(cx, &[&scalar, &rest]).map(Out::Mv)
            },
            | Request::Rotor => {
                let b = read_mv(cx, self.mv, args[0])?;
                rotor(cx, &b, args[1]).map(Out::Mv)
            },
            | Request::Apply => {
                let v = read_mv(cx, self.mv, args[0])?;
                let x = read_mv(cx, self.mv, args[1])?;
                let x = if parity(&v)? { involute(cx, &x) } else { x };
                let inv = inverse(cx, &v)?;
                let vx = gp(cx, &v, &x)?;
                gp(cx, &vx, &inv).map(Out::Mv)
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rules::testing::eval;
    use crate::rules::testing::reduce_with;
    use crate::rules::testing::simplify;

    fn run(src: &str) -> String {
        simplify(&[geometric_algebra()], src)
    }

    fn e(i: u32) -> String {
        format!("mv_basis(list(3, 0, 0), {i})")
    }

    fn v(x: &str, y: &str, z: &str) -> String {
        format!("mv_vector(list(3, 0, 0), list({x}, {y}, {z}))")
    }

    fn mv3(terms: &str) -> String {
        format!("mv(list(3, 0, 0), list({terms}))")
    }

    #[test]
    fn constructors_and_normal_form() {
        assert_eq!(run(&e(2)), mv3("list(2, 1)"));
        assert_eq!(run(&v("x", "0", "z")), mv3("list(1, x), list(4, z)"));
        assert_eq!(run("mv_scalar(list(3, 0, 0), 5)"), mv3("list(0, 5)"));
        assert_eq!(run("mv_blade(list(3, 0, 0), 5, a)"), mv3("list(5, a)"));
        assert_eq!(run("mv_pseudoscalar(list(3, 0, 0))"), mv3("list(7, 1)"));
        // user-built terms are merged, sorted and stripped of zeros
        assert_eq!(run("mv(list(3, 0, 0), list(list(2, y), list(1, x), list(2, -y), list(0, 0)))"), mv3("list(1, x)"));
        assert_eq!(run("mv(list(3, 0, 0), list(list(1, x), list(1, y)))"), mv3("list(1, x + y)"));
        // out of range blades are left alone
        assert!(!reduce_with(&[geometric_algebra()], "mv_blade(list(3, 0, 0), 8, 1)", &[]).1);
    }

    #[test]
    fn basis_products_follow_the_metric() {
        assert_eq!(run(&format!("ga_gp({}, {})", e(1), e(2))), mv3("list(3, 1)"));
        assert_eq!(run(&format!("ga_gp({}, {})", e(2), e(1))), mv3("list(3, -1)"));
        assert_eq!(run(&format!("ga_gp({}, {})", e(3), e(3))), mv3("list(0, 1)"));
        // e1 e2 e3 squares to -1 in Euclidean 3-space
        assert_eq!(run("ga_gp(mv_pseudoscalar(list(3, 0, 0)), mv_pseudoscalar(list(3, 0, 0)))"), mv3("list(0, -1)"));
        // Minkowski plane: e2^2 = -1; degenerate: e2^2 = 0
        assert_eq!(run("ga_gp(mv_basis(list(1, 1, 0), 2), mv_basis(list(1, 1, 0), 2))"), "mv(list(1, 1, 0), list(list(0, -1)))");
        assert_eq!(run("ga_gp(mv_basis(list(1, 0, 1), 2), mv_basis(list(1, 0, 1), 2))"), "mv(list(1, 0, 1), list())");
        // scalars promote
        assert_eq!(run(&format!("ga_gp(3, {})", e(2))), mv3("list(2, 3)"));
    }

    #[test]
    fn geometric_product_of_vectors_splits_into_dot_and_wedge() {
        let a = v("a1", "a2", "a3");
        let b = v("b1", "b2", "b3");
        let gp = run(&format!("ga_gp({a}, {b})"));
        let dot = run(&format!("ga_inner({a}, {b})"));
        let wedge = run(&format!("ga_wedge({a}, {b})"));
        assert_eq!(dot, mv3("list(0, a1*b1 + a2*b2 + a3*b3)"));
        assert_eq!(run(&format!("mv_add({dot}, {wedge})")), gp);
        assert_eq!(run(&format!("ga_wedge({}, {})", v("a", "b", "0"), v("c", "d", "0"))), mv3("list(3, a*d - b*c)"));
        // a ^ a = 0
        assert_eq!(run(&format!("ga_wedge({a}, {a})")), mv3(""));
    }

    #[test]
    fn product_is_associative() {
        let a = "mv(list(3, 0, 0), list(list(0, 2), list(1, 3), list(6, -1)))";
        let b = "mv(list(3, 0, 0), list(list(2, 5), list(3, -2), list(7, 7)))";
        let c = "mv(list(3, 0, 0), list(list(1, 4), list(5, 3), list(0, -6)))";
        let left = run(&format!("ga_gp(ga_gp({a}, {b}), {c})"));
        let right = run(&format!("ga_gp({a}, ga_gp({b}, {c}))"));
        assert_eq!(run(&format!("mv_sub({left}, {right})")), mv3(""));
    }

    #[test]
    fn contractions() {
        let e12 = "mv_blade(list(3, 0, 0), 3, 1)";
        // e1 | e12 = e2, e12 |_ e2 = e1, e2 | e1 = 0 (grade would be negative)
        assert_eq!(run(&format!("ga_lcontract({}, {e12})", e(1))), mv3("list(2, 1)"));
        assert_eq!(run(&format!("ga_rcontract({e12}, {})", e(2))), mv3("list(1, 1)"));
        assert_eq!(run(&format!("ga_lcontract({e12}, {})", e(1))), mv3(""));
        assert_eq!(run(&format!("ga_lcontract({}, {})", e(1), e(1))), mv3("list(0, 1)"));
        // Hestenes inner product of a scalar vanishes, of two bivectors is a scalar
        assert_eq!(run(&format!("ga_inner(2, {})", e(1))), mv3(""));
        assert_eq!(run(&format!("ga_inner({e12}, {e12})")), mv3("list(0, -1)"));
        assert_eq!(run(&format!("ga_scalar_product({e12}, {e12})")), mv3("list(0, -1)"));
        // commutator of e1, e2 is e12
        assert_eq!(run(&format!("ga_commutator({}, {})", e(1), e(2))), mv3("list(3, 1)"));
    }

    #[test]
    fn grades_and_involutions() {
        let m = "mv(list(3, 0, 0), list(list(0, s), list(1, a), list(3, b), list(7, c)))";
        assert_eq!(run(&format!("ga_grade({m}, 0)")), mv3("list(0, s)"));
        assert_eq!(run(&format!("ga_grade({m}, 2)")), mv3("list(3, b)"));
        assert_eq!(run(&format!("ga_grade({m}, 3)")), mv3("list(7, c)"));
        assert_eq!(run(&format!("ga_grade({m}, 1)")), run("mv_blade(list(3, 0, 0), 1, a)"));
        assert_eq!(run(&format!("ga_reverse({m})")), mv3("list(0, s), list(1, a), list(3, -b), list(7, -c)"));
        assert_eq!(run(&format!("ga_involute({m})")), mv3("list(0, s), list(1, -a), list(3, b), list(7, -c)"));
        assert_eq!(run(&format!("ga_conjugate({m})")), mv3("list(0, s), list(1, -a), list(3, -b), list(7, c)"));
        assert_eq!(run(&format!("mv_scalar_part({m})")), "s");
        assert_eq!(run(&format!("mv_coeff({m}, 3)")), "b");
        assert_eq!(run(&format!("mv_coeff({m}, 2)")), "0");
        // (AB)~ = B~ A~
        let (a, b) = (v("a1", "a2", "a3"), "mv_blade(list(3, 0, 0), 6, k)".to_owned());
        let lhs = run(&format!("ga_reverse(ga_gp({a}, {b}))"));
        let rhs = run(&format!("ga_gp(ga_reverse({b}), ga_reverse({a}))"));
        assert_eq!(lhs, rhs);
    }

    #[test]
    fn linear_structure() {
        assert_eq!(run(&format!("mv_add({}, {}, 2)", e(1), e(2))), mv3("list(0, 2), list(1, 1), list(2, 1)"));
        assert_eq!(run(&format!("mv_sub({}, {})", e(1), e(1))), mv3(""));
        assert_eq!(run(&format!("mv_scale(k, {})", v("x", "y", "z"))), mv3("list(1, k*x), list(2, k*y), list(4, k*z)"));
        // different algebras do not mix
        assert!(!reduce_with(&[geometric_algebra()], "mv_add(mv_basis(list(3, 0, 0), 1), mv_basis(list(2, 0, 0), 1))", &[]).1);
    }

    #[test]
    fn dual_in_three_dimensions() {
        // A* = A I^-1 with I^-1 = -e123
        assert_eq!(run(&format!("ga_dual({})", e(1))), mv3("list(6, -1)"));
        assert_eq!(run(&format!("ga_dual({})", e(3))), mv3("list(3, -1)"));
        assert_eq!(run("ga_dual(mv_scalar(list(3, 0, 0), 1))"), mv3("list(7, -1)"));
        assert_eq!(run("ga_dual(mv_pseudoscalar(list(3, 0, 0)))"), mv3("list(0, 1)"));
        let m = "mv(list(3, 0, 0), list(list(0, s), list(1, a), list(3, b), list(7, c)))";
        assert_eq!(run(&format!("ga_undual(ga_dual({m}))")), run(m));
        // x ^ x* is the volume element up to sign: e1 ^ (-e23) = -e123
        assert_eq!(run(&format!("ga_wedge({}, ga_dual({}))", e(1), e(1))), mv3("list(7, -1)"));
        // degenerate metric: M I
        assert_eq!(run("ga_dual(mv_basis(list(1, 0, 1), 1))"), "mv(list(1, 0, 1), list(list(2, 1)))");
    }

    #[test]
    fn norms_and_inverses() {
        let x = v("x", "y", "z");
        assert_eq!(run(&format!("ga_norm_squared({x})")), "x^2 + y^2 + z^2");
        assert_eq!(run(&format!("ga_norm({x})")), "(x^2 + y^2 + z^2)^(1/2)");
        assert_eq!(run(&format!("ga_norm({})", v("3", "4", "0"))), "5");
        assert_eq!(run(&format!("ga_normalize({})", v("3", "0", "4"))), mv3("list(1, 3/5), list(4, 4/5)"));
        assert_eq!(run(&format!("ga_inverse({})", v("2", "0", "0"))), mv3("list(1, 1/2)"));
        // v v^-1 = 1 for a symbolic vector
        let one = run(&format!("mv_scalar_part(ga_gp({x}, ga_inverse({x})))"));
        assert!((eval(&[geometric_algebra()], &one, &[("x", 1.0), ("y", 2.0), ("z", 3.0)]) - 1.0).abs() < 1e-12, "{one}");
        // a bivector squares to -|B|^2
        assert_eq!(run("ga_norm_squared(mv_blade(list(3, 0, 0), 6, 3))"), "9");
        // a null vector has no inverse
        assert!(!reduce_with(&[geometric_algebra()], "ga_inverse(mv_basis(list(1, 0, 1), 2))", &[]).1);
    }

    #[test]
    fn exponential_and_rotor() {
        assert_eq!(run("ga_exp(mv_blade(list(3, 0, 0), 3, t))"), mv3("list(0, cos(t)), list(3, sin(t))"));
        assert_eq!(run("ga_exp(mv_blade(list(1, 1, 0), 3, t))"), "mv(list(1, 1, 0), list(list(0, cosh(t)), list(list(3, sinh(t)))))".replace("list(list(3, sinh(t)))", "list(3, sinh(t))"));
        // rotor of the e12 plane
        assert_eq!(run("ga_rotor(mv_blade(list(3, 0, 0), 3, 1), t)"), mv3("list(0, cos(1/2*t)), list(3, -sin(1/2*t))"));
        // the rotor of a non-unit bivector is normalised
        assert_eq!(run("ga_rotor(mv_blade(list(3, 0, 0), 3, 2), t)"), run("ga_rotor(mv_blade(list(3, 0, 0), 3, 1), t)"));
        // R rotates e1 to cos t e1 + sin t e2
        assert_eq!(run(&format!("ga_apply_versor(ga_rotor(mv_blade(list(3, 0, 0), 3, 1), t), {})", e(1))), mv3("list(1, cos(t)), list(2, sin(t))"));
        // by a quarter turn e1 goes to e2
        let r = "ga_rotor(mv_blade(list(3, 0, 0), 3, 1), pi/2)";
        let x = run(&format!("mv_coeff(ga_apply_versor({r}, {}), 1)", e(1)));
        let y = run(&format!("mv_coeff(ga_apply_versor({r}, {}), 2)", e(1)));
        assert!(eval(&[geometric_algebra()], &x, &[]).abs() < 1e-12, "{x}");
        assert!((eval(&[geometric_algebra()], &y, &[]) - 1.0).abs() < 1e-12, "{y}");
        // e3 is untouched, and a rotor is a unit versor
        assert_eq!(run(&format!("ga_apply_versor(ga_rotor(mv_blade(list(3, 0, 0), 3, 1), t), {})", e(3))), mv3("list(4, 1)"));
        assert_eq!(run("ga_norm_squared(ga_rotor(mv_blade(list(3, 0, 0), 5, 1), t))"), "1");
        // not a bivector: refused
        assert!(!reduce_with(&[geometric_algebra()], &format!("ga_rotor({}, t)", e(1)), &[]).1);
    }

    #[test]
    fn versor_sandwich() {
        // reflection in the plane orthogonal to e1 (odd versor)
        assert_eq!(run(&format!("ga_apply_versor({}, {})", e(1), v("1", "1", "0"))), mv3("list(1, -1), list(2, 1)"));
        // a rotor times its inverse is one
        let r = "ga_rotor(mv_blade(list(3, 0, 0), 6, 1), a)";
        assert_eq!(run(&format!("ga_gp({r}, ga_inverse({r}))")), mv3("list(0, 1)"));
    }
}
