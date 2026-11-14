//! Real formulas over polynomial atoms: reading them from terms, deciding
//! them from a sign vector, and writing quantifier-free formulas back.
//!
//! An atom is always normalised to `p = 0`, `p != 0`, `p < 0`, `p <= 0`,
//! `p > 0` or `p >= 0` for a polynomial `p` over `Q` in a fixed list of
//! variables (the indices are the [`super::qe_poly::MPoly`] variable indices).
//! Boolean structure is `not`, `and`, `or`; `implies`, `iff` and `xor` are
//! expanded on input. Quantifiers are blocks of variables
//! (`exists`/`forall` with a symbol or a `list` prefix). Anything that is
//! not a polynomial comparison is rejected, which is what lets the Boolean
//! decision-diagram path keep its own semantics.

#![allow(
    clippy::missing_const_for_fn,
    clippy::must_use_candidate,
    clippy::option_if_let_else,
    clippy::redundant_pub_crate,
    clippy::single_call_fn
)]

use std::collections::BTreeSet;

use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::Zero;

use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::SymbolId;
use crate::graph::op::core;
use crate::rules::poly::repr::Gens;
use crate::rules::poly::repr::Limits as ReprLimits;
use crate::rules::poly::repr::Poly as ReprPoly;
use crate::rules::poly::repr::from_term;

use super::qe_poly::MPoly;

/// A comparison operator.
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum Cmp {
    /// `=`.
    Eq,
    /// `!=`.
    Ne,
    /// `<`.
    Lt,
    /// `<=`.
    Le,
    /// `>`.
    Gt,
    /// `>=`.
    Ge,
}

impl Cmp {
    /// Whether a polynomial with the given sign (at the sample point)
    /// satisfies the comparison against zero.
    #[must_use]
    pub(crate) fn holds_sign(
        self,
        sign: i32,
    ) -> bool {
        match self {
            | Self::Eq => sign == 0,
            | Self::Ne => sign != 0,
            | Self::Lt => sign < 0,
            | Self::Le => sign <= 0,
            | Self::Gt => sign > 0,
            | Self::Ge => sign >= 0,
        }
    }

    /// The comparison for a sign condition: `-1`, `0` or `1`.
    #[must_use]
    pub(crate) fn from_sign(sign: i32) -> Self {
        match sign {
            | -1 => Self::Lt,
            | 0 => Self::Eq,
            | _ => Self::Gt,
        }
    }
}

/// A polynomial atom.
#[derive(Clone, Debug)]
pub(crate) struct Atom {
    /// The comparison.
    pub cmp: Cmp,
    /// The polynomial compared against zero.
    pub poly: MPoly,
}

/// A formula over polynomial atoms.
#[derive(Clone, Debug)]
pub(crate) enum Formula {
    /// A truth constant.
    Truth(bool),
    /// A polynomial comparison.
    Atom(Atom),
    /// Negation.
    Not(Box<Self>),
    /// Conjunction.
    And(Vec<Self>),
    /// Disjunction.
    Or(Vec<Self>),
    /// Existential block.
    Exists(Vec<u32>, Box<Self>),
    /// Universal block.
    Forall(Vec<u32>, Box<Self>),
}

impl Formula {
    /// Whether the formula has no quantifier.
    #[must_use]
    #[allow(dead_code)] // used by the module's own tests
    pub(crate) fn is_quantifier_free(&self) -> bool {
        match self {
            | Self::Truth(_) | Self::Atom(_) => true,
            | Self::Not(x) => x.is_quantifier_free(),
            | Self::And(xs) | Self::Or(xs) => xs.iter().all(Self::is_quantifier_free),
            | Self::Exists(_, _) | Self::Forall(_, _) => false,
        }
    }

    /// Collects the variable indices of the atoms.
    pub(crate) fn atom_variables(
        &self,
        out: &mut BTreeSet<u32>,
    ) {
        match self {
            | Self::Truth(_) => {},
            | Self::Atom(atom) => out.extend(atom.poly.support()),
            | Self::Not(x) => x.atom_variables(out),
            | Self::And(xs) | Self::Or(xs) => {
                for x in xs {
                    x.atom_variables(out);
                }
            },
            | Self::Exists(_, body) | Self::Forall(_, body) => body.atom_variables(out),
        }
    }

    /// The free variable indices.
    ///
    /// Binding is tracked along the recursion: a variable bound by an inner
    /// quantifier is free again once that quantifier is left, so a global set
    /// difference would hide a free occurrence that shares a name with a bound
    /// one somewhere else in the formula.
    #[must_use]
    pub(crate) fn free_variables(&self) -> BTreeSet<u32> {
        let mut out = BTreeSet::new();
        self.free_variables_from(&BTreeSet::new(), &mut out);
        out
    }

    fn free_variables_from(
        &self,
        bound: &BTreeSet<u32>,
        out: &mut BTreeSet<u32>,
    ) {
        match self {
            | Self::Truth(_) => {},
            | Self::Atom(atom) => {
                out.extend(
                    atom.poly
                        .support()
                        .into_iter()
                        .filter(|v| !bound.contains(v)),
                );
            },
            | Self::Not(x) => x.free_variables_from(bound, out),
            | Self::And(xs) | Self::Or(xs) => {
                for x in xs {
                    x.free_variables_from(bound, out);
                }
            },
            | Self::Exists(vars, body) | Self::Forall(vars, body) => {
                let mut inner = bound.clone();
                inner.extend(vars.iter().copied());
                body.free_variables_from(&inner, out);
            },
        }
    }

    /// Truth of a quantifier-free formula, with each atom's polynomial
    /// looked up in `sign`.
    #[must_use]
    pub(crate) fn truth_with(
        &self,
        sign: &dyn Fn(&MPoly) -> i32,
    ) -> bool {
        match self {
            | Self::Truth(value) => *value,
            | Self::Atom(atom) => atom.cmp.holds_sign(sign(&atom.poly)),
            | Self::Not(x) => !x.truth_with(sign),
            | Self::And(xs) => xs.iter().all(|x| x.truth_with(sign)),
            | Self::Or(xs) => xs.iter().any(|x| x.truth_with(sign)),
            | Self::Exists(_, _) | Self::Forall(_, _) => false,
        }
    }

    /// Evaluates a quantifier-free formula at a rational assignment of all
    /// variables.
    #[must_use]
    pub(crate) fn eval_rational(
        &self,
        values: &[BigRational],
    ) -> Option<bool> {
        match self {
            | Self::Truth(value) => Some(*value),
            | Self::Atom(atom) => {
                let value = atom.poly.eval_rational(values)?;
                let sign = if value.is_zero() {
                    0
                } else if value.is_positive() {
                    1
                } else {
                    -1
                };
                Some(atom.cmp.holds_sign(sign))
            },
            | Self::Not(x) => Some(!x.eval_rational(values)?),
            | Self::And(xs) => {
                let mut result = true;
                for x in xs {
                    result &= x.eval_rational(values)?;
                }
                Some(result)
            },
            | Self::Or(xs) => {
                let mut result = false;
                for x in xs {
                    result |= x.eval_rational(values)?;
                }
                Some(result)
            },
            | Self::Exists(_, _) | Self::Forall(_, _) => None,
        }
    }

    /// Conjunction with constants and duplicates folded in.
    #[must_use]
    pub(crate) fn and(parts: Vec<Self>) -> Self {
        Self::junction(parts, true)
    }

    /// Disjunction with constants and duplicates folded in.
    #[must_use]
    pub(crate) fn or(parts: Vec<Self>) -> Self {
        Self::junction(parts, false)
    }

    fn junction(
        parts: Vec<Self>,
        conjunction: bool,
    ) -> Self {
        let mut kept: Vec<Self> = Vec::with_capacity(parts.len());
        let mut literals: BTreeSet<(u8, MPoly)> = BTreeSet::new();
        for part in parts {
            match part {
                | Self::Truth(true) if conjunction => {},
                | Self::Truth(false) if conjunction => return Self::Truth(false),
                | Self::Truth(false) => {},
                | Self::Truth(true) => return Self::Truth(true),
                | Self::Atom(atom) => {
                    let key = (cmp_code(atom.cmp), atom.poly.clone());
                    if literals.insert(key) {
                        kept.push(Self::Atom(atom));
                    }
                },
                | other => {
                    if !kept.iter().any(|k| same(k, &other)) {
                        kept.push(other);
                    }
                },
            }
        }
        match kept.len() {
            | 0 => Self::Truth(conjunction),
            | 1 => kept.pop().unwrap_or(Self::Truth(conjunction)),
            | _ if conjunction => Self::And(kept),
            | _ => Self::Or(kept),
        }
    }

    /// Negation with constants folded in.
    #[must_use]
    pub(crate) fn not(self) -> Self {
        match self {
            | Self::Truth(value) => Self::Truth(!value),
            | Self::Not(inner) => *inner,
            | Self::Atom(mut atom) => {
                atom.cmp = match atom.cmp {
                    | Cmp::Eq => Cmp::Ne,
                    | Cmp::Ne => Cmp::Eq,
                    | Cmp::Lt => Cmp::Ge,
                    | Cmp::Le => Cmp::Gt,
                    | Cmp::Gt => Cmp::Le,
                    | Cmp::Ge => Cmp::Lt,
                };
                Self::Atom(atom)
            },
            | other => Self::Not(Box::new(other)),
        }
    }
}

fn cmp_code(cmp: Cmp) -> u8 {
    match cmp {
        | Cmp::Eq => 0,
        | Cmp::Ne => 1,
        | Cmp::Lt => 2,
        | Cmp::Le => 3,
        | Cmp::Gt => 4,
        | Cmp::Ge => 5,
    }
}

/// Structural equality (formulas are small and this is used only for
/// deduplication).
fn same(
    a: &Formula,
    b: &Formula,
) -> bool {
    match (a, b) {
        | (Formula::Truth(x), Formula::Truth(y)) => x == y,
        | (Formula::Atom(x), Formula::Atom(y)) => x.cmp == y.cmp && x.poly == y.poly,
        | (Formula::Not(x), Formula::Not(y)) => same(x, y),
        | (Formula::And(xs), Formula::And(ys)) | (Formula::Or(xs), Formula::Or(ys)) => {
            xs.len() == ys.len() && xs.iter().zip(ys).all(|(x, y)| same(x, y))
        },
        | (Formula::Exists(xs, x), Formula::Exists(ys, y))
        | (Formula::Forall(xs, x), Formula::Forall(ys, y)) => xs == ys && same(x, y),
        | _ => false,
    }
}

/// The operators used to read and write formulas.
#[derive(Copy, Clone, Debug)]
pub(crate) struct FormulaOps {
    /// `and`.
    pub and: crate::graph::OpId,
    /// `or`.
    pub or: crate::graph::OpId,
    /// `not`.
    pub not: crate::graph::OpId,
    /// `implies`.
    pub implies: crate::graph::OpId,
    /// `iff`.
    pub iff: crate::graph::OpId,
    /// `xor`.
    pub xor: crate::graph::OpId,
    /// `exists`.
    pub exists: crate::graph::OpId,
    /// `forall`.
    pub forall: crate::graph::OpId,
    /// `lt`.
    pub lt: crate::graph::OpId,
    /// `le`.
    pub le: crate::graph::OpId,
    /// `gt`.
    pub gt: crate::graph::OpId,
    /// `ge`.
    pub ge: crate::graph::OpId,
    /// `ne`.
    pub ne: crate::graph::OpId,
    /// The truth literal.
    pub truth: crate::graph::OpId,
    /// The false literal.
    pub falsity: crate::graph::OpId,
}

impl FormulaOps {
    /// Looks the operators up in the graph. `None` when the formula
    /// operators are not installed (then the request is not a real formula
    /// request either).
    #[must_use]
    pub(crate) fn lookup(graph: &Graph) -> Option<Self> {
        let of = |name: &str| graph.ops().lookup(name);
        Some(Self {
            and: of("and")?,
            or: of("or")?,
            not: of("not")?,
            implies: of("implies")?,
            iff: of("iff")?,
            xor: of("xor")?,
            exists: of("exists")?,
            forall: of("forall")?,
            lt: of("lt")?,
            le: of("le")?,
            gt: of("gt")?,
            ge: of("ge")?,
            ne: of("ne")?,
            truth: of("true")?,
            falsity: of("false")?,
        })
    }

    /// The comparison for an operator, if any.
    #[must_use]
    pub(crate) fn comparison(
        self,
        op: crate::graph::OpId,
    ) -> Option<Cmp> {
        if op == self.lt {
            Some(Cmp::Lt)
        } else if op == self.le {
            Some(Cmp::Le)
        } else if op == self.gt {
            Some(Cmp::Gt)
        } else if op == self.ge {
            Some(Cmp::Ge)
        } else if op == self.ne {
            Some(Cmp::Ne)
        } else if op == core::EQ {
            Some(Cmp::Eq)
        } else {
            None
        }
    }

    /// The operator for a comparison.
    #[must_use]
    pub(crate) fn op_of(
        self,
        cmp: Cmp,
    ) -> crate::graph::OpId {
        match cmp {
            | Cmp::Eq => core::EQ,
            | Cmp::Ne => self.ne,
            | Cmp::Lt => self.lt,
            | Cmp::Le => self.le,
            | Cmp::Gt => self.gt,
            | Cmp::Ge => self.ge,
        }
    }
}

/// A parsed formula together with the variable table.
pub(crate) struct Parsed {
    /// The formula.
    pub formula: Formula,
    /// `variables[i]` is the symbol of variable index `i`.
    pub variables: Vec<SymbolId>,
    /// `nodes[i]` is the symbol leaf of variable index `i`.
    pub nodes: Vec<NodeId>,
    /// The operators.
    pub ops: FormulaOps,
}

/// Reads `root` as a real formula, adding the variables of `extra` to the
/// table even when they do not occur.
pub(crate) fn parse(
    graph: &mut Graph,
    root: NodeId,
    extra: &[NodeId],
) -> Option<Parsed> {
    let ops = FormulaOps::lookup(graph)?;
    let mut gens = Gens::default();
    let mut variables: Vec<SymbolId> = Vec::new();
    let mut nodes: Vec<NodeId> = Vec::new();
    let formula = read(graph, root, &ops, &mut gens, &mut variables, &mut nodes, 0)?;
    for &node in extra {
        let index = gens.index(graph, node);
        sync(graph, &gens, &mut variables, &mut nodes)?;
        let _ = index;
    }
    Some(Parsed {
        formula,
        variables,
        nodes,
        ops,
    })
}

/// Adds newly seen generators to the variable table, rejecting anything
/// that is not a symbol.
fn sync(
    graph: &Graph,
    gens: &Gens,
    variables: &mut Vec<SymbolId>,
    nodes: &mut Vec<NodeId>,
) -> Option<()> {
    while variables.len() < gens.len() {
        let index = u32::try_from(variables.len()).ok()?;
        let node = gens.node(index)?;
        let symbol = graph.as_symbol(node)?;
        variables.push(symbol);
        nodes.push(node);
    }
    Some(())
}

/// Converts both sides of a comparison to one polynomial.
fn comparison_poly(
    graph: &mut Graph,
    ops: FormulaOps,
    node: NodeId,
    gens: &mut Gens,
    variables: &mut Vec<SymbolId>,
    nodes: &mut Vec<NodeId>,
) -> Option<(Cmp, MPoly)> {
    let cmp = ops.comparison(graph.op(node))?;
    let &[lhs, rhs] = graph.children(node) else {
        return None;
    };
    let left = polynomial(graph, lhs, gens, variables, nodes)?;
    let right = polynomial(graph, rhs, gens, variables, nodes)?;
    let poly = left.sub(&right);
    Some((cmp, poly))
}

/// Converts a term to a polynomial over the shared generator table.
fn polynomial(
    graph: &mut Graph,
    node: NodeId,
    gens: &mut Gens,
    variables: &mut Vec<SymbolId>,
    nodes: &mut Vec<NodeId>,
) -> Option<MPoly> {
    let repr = from_term(graph, gens, node, ReprLimits::default())?;
    sync(graph, gens, variables, nodes)?;
    repr_to_mpoly(&repr, gens.len())
}

/// Converts a [`ReprPoly`] to an [`MPoly`].
#[must_use]
pub(crate) fn repr_to_mpoly(
    repr: &ReprPoly,
    variables: usize,
) -> Option<MPoly> {
    let mut out = MPoly::zero();
    for (mono, coeff) in repr.terms() {
        let value = coeff.to_rational()?;
        let mut mapped = Vec::with_capacity(mono.len());
        for &(generator, exponent) in mono {
            if usize::try_from(generator).ok()? >= variables {
                return None;
            }
            mapped.push((generator, exponent));
        }
        out = out.add(&MPoly::monomial(mapped, value));
    }
    Some(out)
}

const MAX_DEPTH: usize = 256;

fn read(
    graph: &mut Graph,
    node: NodeId,
    ops: &FormulaOps,
    gens: &mut Gens,
    variables: &mut Vec<SymbolId>,
    nodes: &mut Vec<NodeId>,
    depth: usize,
) -> Option<Formula> {
    if depth > MAX_DEPTH {
        return None;
    }
    let op = graph.op(node);
    let children = graph.children(node).to_vec();
    if let Some(value) = boolean_constant(graph, ops, node) {
        return Some(Formula::Truth(value));
    }
    if op == ops.and || op == ops.or {
        let mut parts = Vec::with_capacity(children.len());
        for child in children {
            parts.push(read(
                graph,
                child,
                ops,
                gens,
                variables,
                nodes,
                depth.saturating_add(1),
            )?);
        }
        return Some(if op == ops.and {
            Formula::And(parts)
        } else {
            Formula::Or(parts)
        });
    }
    if op == ops.not {
        let &[child] = children.as_slice() else {
            return None;
        };
        return Some(
            read(
                graph,
                child,
                ops,
                gens,
                variables,
                nodes,
                depth.saturating_add(1),
            )?
            .not(),
        );
    }
    if op == ops.implies || op == ops.iff {
        let &[a, b] = children.as_slice() else {
            return None;
        };
        let left = read(
            graph,
            a,
            ops,
            gens,
            variables,
            nodes,
            depth.saturating_add(1),
        )?;
        let right = read(
            graph,
            b,
            ops,
            gens,
            variables,
            nodes,
            depth.saturating_add(1),
        )?;
        return Some(if op == ops.implies {
            Formula::or(vec![left.not(), right])
        } else {
            Formula::or(vec![
                Formula::and(vec![left.clone(), right.clone()]),
                Formula::and(vec![left.not(), right.not()]),
            ])
        });
    }
    if op == ops.xor {
        return read_xor(graph, &children, ops, gens, variables, nodes, depth);
    }
    if op == ops.exists || op == ops.forall {
        let &[prefix, body] = children.as_slice() else {
            return None;
        };
        let bound_nodes: Vec<NodeId> = if graph.op(prefix) == core::LIST {
            graph.children(prefix).to_vec()
        } else {
            vec![prefix]
        };
        let mut indices = Vec::with_capacity(bound_nodes.len());
        for symbol in &bound_nodes {
            let index = gens.index(graph, *symbol);
            sync(graph, gens, variables, nodes)?;
            indices.push(index);
        }
        let body = read(
            graph,
            body,
            ops,
            gens,
            variables,
            nodes,
            depth.saturating_add(1),
        )?;
        // Drop variables that do not occur.
        let mut occurring = BTreeSet::new();
        body.atom_variables(&mut occurring);
        indices.retain(|i| occurring.contains(i));
        return Some(if indices.is_empty() {
            body
        } else if op == ops.exists {
            Formula::Exists(indices, Box::new(body))
        } else {
            Formula::Forall(indices, Box::new(body))
        });
    }
    if let Some((cmp, poly)) = comparison_poly(graph, *ops, node, gens, variables, nodes) {
        return Some(Formula::Atom(Atom { cmp, poly }));
    }
    None
}

fn read_xor(
    graph: &mut Graph,
    children: &[NodeId],
    ops: &FormulaOps,
    gens: &mut Gens,
    variables: &mut Vec<SymbolId>,
    nodes: &mut Vec<NodeId>,
    depth: usize,
) -> Option<Formula> {
    if children.len() > 8 {
        return None;
    }
    let mut parts = Vec::with_capacity(children.len());
    for &child in children {
        parts.push(read(
            graph,
            child,
            ops,
            gens,
            variables,
            nodes,
            depth.saturating_add(1),
        )?);
    }
    // Parity by the truth table over up to eight operands.
    let count = 1_usize.checked_shl(u32::try_from(parts.len()).ok()?)?;
    let mut alternatives = Vec::new();
    for mask in 0..count {
        if mask.count_ones() % 2 == 0 {
            continue;
        }
        let mut literals = Vec::with_capacity(parts.len());
        for (index, part) in parts.iter().enumerate() {
            literals.push(if mask >> index & 1 == 1 {
                part.clone()
            } else {
                part.clone().not()
            });
        }
        alternatives.push(Formula::and(literals));
    }
    Some(Formula::or(alternatives))
}

fn boolean_constant(
    graph: &Graph,
    ops: &FormulaOps,
    node: NodeId,
) -> Option<bool> {
    let op = graph.op(node);
    if op == ops.truth {
        return Some(true);
    }
    if op == ops.falsity {
        return Some(false);
    }
    match graph.payload(node) {
        | Some(crate::graph::Payload::Bool(value)) if op == core::LIT => Some(*value),
        | _ => None,
    }
}

/// Builds the term for a polynomial over the given variable nodes.
#[must_use]
pub(crate) fn polynomial_term(
    graph: &mut Graph,
    poly: &MPoly,
    nodes: &[NodeId],
) -> Option<NodeId> {
    let mut terms: Vec<NodeId> = Vec::with_capacity(poly.len());
    for (mono, coeff) in poly.terms() {
        let mut factors: Vec<NodeId> = Vec::new();
        if !coeff.is_one() || mono.is_empty() {
            factors.push(graph.num(Number::rat(coeff.clone())));
        }
        for &(variable, exponent) in mono {
            let base = *nodes.get(usize::try_from(variable).ok()?)?;
            if exponent == 1 {
                factors.push(base);
            } else {
                let power = graph.int(i64::from(exponent));
                factors.push(graph.node(core::POW, &[base, power]));
            }
        }
        terms.push(match factors.as_slice() {
            | [only] => *only,
            | _ => graph.node(core::MUL, &factors),
        });
    }
    Some(match terms.as_slice() {
        | [] => graph.int(0),
        | [only] => *only,
        | _ => graph.node(core::ADD, &terms),
    })
}

/// Builds an atom `poly cmp 0`.
#[must_use]
pub(crate) fn atom_term(
    graph: &mut Graph,
    ops: &FormulaOps,
    atom: &Atom,
    nodes: &[NodeId],
) -> Option<NodeId> {
    let poly = polynomial_term(graph, &atom.poly, nodes)?;
    let zero = graph.int(0);
    Some(graph.node(ops.op_of(atom.cmp), &[poly, zero]))
}

/// Builds the term for a formula.
#[must_use]
pub(crate) fn term(
    graph: &mut Graph,
    ops: &FormulaOps,
    formula: &Formula,
    nodes: &[NodeId],
) -> Option<NodeId> {
    match formula {
        | Formula::Truth(value) => Some(graph.lit(crate::graph::Payload::Bool(*value))),
        | Formula::Atom(atom) => atom_term(graph, ops, atom, nodes),
        | Formula::Not(inner) => {
            let inner = term(graph, ops, inner, nodes)?;
            Some(graph.node(ops.not, &[inner]))
        },
        | Formula::And(parts) => {
            let mut items = Vec::with_capacity(parts.len());
            for part in parts {
                items.push(term(graph, ops, part, nodes)?);
            }
            Some(match items.as_slice() {
                | [] => graph.lit(crate::graph::Payload::Bool(true)),
                | [only] => *only,
                | _ => graph.node(ops.and, &items),
            })
        },
        | Formula::Or(parts) => {
            let mut items = Vec::with_capacity(parts.len());
            for part in parts {
                items.push(term(graph, ops, part, nodes)?);
            }
            Some(match items.as_slice() {
                | [] => graph.lit(crate::graph::Payload::Bool(false)),
                | [only] => *only,
                | _ => graph.node(ops.or, &items),
            })
        },
        | Formula::Exists(_, _) | Formula::Forall(_, _) => None,
    }
}

/// Builds the term for an element of a number field, as a polynomial in the
/// field's generator (whose term is a `root_of`, a rational or a radical).
#[must_use]
pub(crate) fn elem_term(
    graph: &mut Graph,
    field: &super::qe_field::Field,
    value: &super::qe_field::Elem,
) -> Option<NodeId> {
    if let Some(constant) = value.constant() {
        return Some(graph.num(Number::rat(constant)));
    }
    let generator = field.root().to_term(graph)?;
    let mut terms = Vec::with_capacity(value.coeffs().len());
    for (power, coefficient) in value.coeffs().iter().enumerate() {
        if coefficient.is_zero() {
            continue;
        }
        if power == 0 {
            terms.push(graph.num(Number::rat(coefficient.clone())));
            continue;
        }
        let factor = graph.num(Number::rat(coefficient.clone()));
        let piece = if power == 1 {
            generator
        } else {
            let exponent = graph.int(i64::try_from(power).ok()?);
            graph.node(core::POW, &[generator, exponent])
        };
        terms.push(if coefficient.is_one() {
            piece
        } else {
            graph.node(core::MUL, &[factor, piece])
        });
    }
    Some(match terms.as_slice() {
        | [] => graph.int(0),
        | [only] => *only,
        | _ => graph.node(core::ADD, &terms),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Engine;
    use crate::graph::Graph;

    fn installed_graph() -> Graph {
        let mut graph = Graph::new();
        assert!(Engine::install(&mut graph, &crate::rules::standard()).is_ok());
        graph
    }

    fn read_source(source: &str) -> Parsed {
        let mut graph = installed_graph();
        let root = graph.parse(source).unwrap_or_else(|e| panic!("{e}"));
        parse(&mut graph, root, &[]).unwrap_or_else(|| panic!("parse"))
    }

    #[test]
    fn truth_constants() {
        assert!(read_source("true").formula.truth_with(&|_| 0));
        assert!(!read_source("false").formula.truth_with(&|_| 0));
    }

    #[test]
    fn atom_normal_form() {
        let parsed = read_source("lt(x^2 - 1, 0)");
        assert!(parsed.formula.is_quantifier_free());
        let Formula::Atom(atom) = &parsed.formula else {
            panic!()
        };
        assert_eq!(atom.cmp, Cmp::Lt);
        assert_eq!(atom.poly.degree_in(0), 2);
        assert_eq!(parsed.variables.len(), 1);
    }

    #[test]
    fn boolean_expansion() {
        let parsed = read_source("implies(lt(x, 0), gt(y, 1))");
        assert!(parsed.formula.is_quantifier_free());
        assert_eq!(parsed.variables.len(), 2);
        let signs = |poly: &MPoly| -> i32 {
            if poly == &MPoly::variable(0) {
                -1
            } else {
                1
            }
        };
        assert!(parsed.formula.truth_with(&signs));
    }

    #[test]
    fn quantifiers_are_read() {
        let parsed = read_source("exists(list(x, y), lt(x^2 + y^2, 1))");
        let Formula::Exists(vars, _) = &parsed.formula else {
            panic!()
        };
        assert_eq!(vars.len(), 2);
        assert!(parsed.formula.free_variables().is_empty());
    }
}
