//! The signature of every operator, and where it came from.
//!
//! # Why a derivation, and why it is trustworthy
//!
//! A thousand operators do not get a thousand hand-written signatures, and a
//! host that wants to generate typed bindings cannot do it from names. The
//! earlier answer was a blanket `(Real, …, Real) -> Real`, which was wrong for
//! almost every operator in a computer-algebra system: `sin(x)` with a symbol
//! `x` returns an expression, not a number.
//!
//! This module derives the signature from what the operator **declares about
//! itself**, and every rule below is a consequence of that declaration rather
//! than a guess about the name:
//!
//! | The operator declares | Therefore |
//! |---|---|
//! | a scalar evaluator (`eval`) | it accepts `n` reals and returns a real |
//! | `OpFlags::PREDICATE` | the result is a truth value |
//! | a `Binder` at index `k` | argument `k` is the bound *symbol*, not a value |
//! | `Arity::Variadic` | the last parameter repeats |
//! | `OpFlags::LEAF` | it carries a payload; it has no call signature |
//!
//! and, on top of that, **every operator that is not a leaf also has the
//! symbolic contract**: `n` terms in, one term out. An operator with a scalar
//! evaluator therefore has *two* overloads — `(Real, …) -> Real` for numbers
//! and `(Term, …) -> Term` for expressions — which is exactly what the engine
//! does: `sin(0.5)` folds to a number, `sin(x)` stays a term. The
//! [overload matcher](crate::kernel::overload) picks between them by argument
//! type, so a typed caller and the runtime cannot disagree.
//!
//! # Precedence
//!
//! 1. a signature the rule set **declared** on the descriptor;
//! 2. a correction from `signatures.overlay.json`, written by hand;
//! 3. the derivation below, reported as [`SignatureSource::Derived`], extended
//!    by the overloads the engine was *seen* to satisfy
//!    ([`SignatureSource::Witnessed`], see [`crate::kernel::sort`]);
//! 4. nothing, with a [`DynamicReason`] saying why.
//!
//! A derived signature is marked as such rather than as declared, because a
//! host is entitled to know which of the three it is looking at. It is
//! nonetheless *sound*: it states what the operator's own registration implies,
//! never what its name suggests.

use crate::graph::Arity;
use crate::graph::OpDescriptor;
use crate::graph::OpFlags;
use crate::kernel::manifest::SignatureSource;
use crate::kernel::manifest::overlay;
use crate::kernel::signature::DimExpr;
use crate::kernel::signature::DynamicReason;
use crate::kernel::signature::RepExpr;
use crate::kernel::signature::Signature;
use crate::kernel::signature::TypeExpr;
use crate::kernel::signature::Variadic;

/// Every signature of `desc`, and where they came from.
///
/// This is the single answer to "what does this operator accept", used by the
/// manifest, by `Session::overloads`, by `Session::resolve_overload` and by the
/// audit binary.
#[must_use]
pub fn signatures(desc: &OpDescriptor) -> (Vec<Signature>, SignatureSource) {
    if !desc.overloads.is_empty() {
        return (desc.overloads.clone(), SignatureSource::Declared);
    }
    if let Some(declared) = overlay().get(&*desc.name)
        && !declared.is_empty()
    {
        return (declared.clone(), SignatureSource::Overlay);
    }
    match derive(desc) {
        | Some(mut signatures) => {
            let witnessed = witnessed(desc, &signatures);
            if witnessed.is_empty() {
                (signatures, SignatureSource::Derived)
            } else {
                // Derived overloads keep their positions (index 0 stays the
                // primary overload); the witnessed ones follow.
                signatures.extend(witnessed);
                (signatures, SignatureSource::Witnessed)
            }
        },
        | None => (Vec::new(), SignatureSource::None),
    }
}

/// Where each overload of `desc` came from, in the order [`signatures`]
/// returns them.
///
/// An operator's [`SignatureSource`] names the *best* evidence it has; this
/// is the per-overload answer a binding generator needs to tell the derived
/// symbolic contract from an overload that was seen to hold.
#[must_use]
pub fn overload_sources(desc: &OpDescriptor) -> Vec<SignatureSource> {
    let (all, source) = signatures(desc);
    if source != SignatureSource::Witnessed {
        return vec![source; all.len()];
    }
    let derived = derive(desc).map_or(0, |d| d.len());
    (0..all.len())
        .map(|index| {
            if index < derived {
                SignatureSource::Derived
            } else {
                SignatureSource::Witnessed
            }
        })
        .collect()
}

/// The shipped witnessed signatures of `desc` that are coherent with
/// `derived` (see [`crate::kernel::sort`]).
///
/// The table is mined offline, so an entry can be stale with respect to the
/// registration it is read against (an operator changed arity, or gained a
/// declared signature). Coherence is therefore re-checked here: a witnessed
/// overload never makes a call ambiguous that the registration resolves.
fn witnessed(
    desc: &OpDescriptor,
    derived: &[Signature],
) -> Vec<Signature> {
    let mut kept: Vec<Signature> = derived.to_vec();
    let mut out = Vec::new();
    for entry in crate::kernel::sort::witnessed_for(desc) {
        // A refinement sharpens the result of an overload that exists; as an
        // overload of its own it would be ambiguous with that one.
        if !entry.refines && crate::kernel::sort::coherent(&kept, &entry.signature) {
            kept.push(entry.signature.clone());
            out.push(entry.signature.clone());
        }
    }
    out
}

/// The signatures the rules above imply, or `None` when nothing can be said.
#[must_use]
pub fn derive(desc: &OpDescriptor) -> Option<Vec<Signature>> {
    // A leaf carries a payload rather than arguments; it is not called.
    if desc.flags.has(OpFlags::LEAF) {
        return None;
    }
    let arity = parameter_count(desc);
    let bound = desc.binder.map(|binder| usize::from(binder.var));
    let mut out = Vec::with_capacity(2);
    // The numeric contract, when the operator has one.
    if desc.eval.is_some() {
        let ret = if desc.flags.has(OpFlags::PREDICATE) {
            TypeExpr::Bool
        } else {
            TypeExpr::Real
        };
        let mut signature = Signature::new(vec![TypeExpr::Real; arity], ret);
        if matches!(desc.arity, Arity::Variadic) {
            signature = signature.variadic(Variadic {
                min: 0,
                max: None,
                elem: TypeExpr::Real,
            });
        }
        out.push(signature);
    }
    // The symbolic contract: every non-leaf operator is reachable with
    // expressions, whether or not it also has a numeric evaluator.
    let mut symbolic = Signature::new(
        (0..arity)
            .map(|index| {
                if Some(index) == bound {
                    TypeExpr::Symbol
                } else {
                    term()
                }
            })
            .collect(),
        term(),
    );
    if matches!(desc.arity, Arity::Variadic) {
        symbolic = symbolic.variadic(Variadic {
            min: 0,
            max: None,
            elem: term(),
        });
    }
    out.push(symbolic);
    Some(out)
}

/// The type of a symbolic expression.
#[must_use]
pub fn term() -> TypeExpr {
    TypeExpr::Dyn(RepExpr {
        repr: "term".to_owned(),
        dim: Some(DimExpr::Any),
    })
}

/// Why an operator has no signature.
///
/// `NotYetSpecified` is the last resort, not the default: a leaf carries a
/// payload rather than arguments, and a rule-only operator has nothing to call.
#[must_use]
pub fn reason(desc: &OpDescriptor) -> DynamicReason {
    if desc.flags.has(OpFlags::LEAF) {
        return DynamicReason::OpaquePayload;
    }
    if matches!(desc.arity, Arity::Variadic) {
        return DynamicReason::Variadic;
    }
    if desc.eval.is_none() && desc.eval_dyn.is_none() {
        return DynamicReason::Request;
    }
    DynamicReason::NotYetSpecified
}

/// The number of parameters to declare for `desc`.
const fn parameter_count(desc: &OpDescriptor) -> usize {
    match desc.arity {
        | Arity::Fixed(n) => n as usize,
        | Arity::Variadic => 0,
    }
}

/// The first signature of `desc`, for callers that want just one.
#[must_use]
pub fn primary(desc: &OpDescriptor) -> Option<Signature> {
    signatures(desc).0.into_iter().next()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn numeric(
        name: &str,
        arity: Arity,
    ) -> OpDescriptor {
        OpDescriptor::new(name, arity).eval(|a| a.first().copied().unwrap_or(f64::NAN))
    }

    #[test]
    fn a_scalar_evaluator_gets_a_numeric_and_a_symbolic_overload() {
        let (signatures, source) = signatures(&numeric("demo_sin", Arity::Fixed(1)));
        assert_eq!(source, SignatureSource::Derived);
        assert_eq!(signatures.len(), 2);
        assert_eq!(signatures[0].params, vec![TypeExpr::Real]);
        assert_eq!(signatures[0].ret, TypeExpr::Real);
        assert_eq!(signatures[1].params, vec![term()]);
        assert_eq!(signatures[1].ret, term());
    }

    #[test]
    fn a_predicate_returns_a_truth_value() {
        let desc = OpDescriptor::new("demo_lt", Arity::Fixed(2))
            .flags(OpFlags::PREDICATE)
            .eval(|_| 1.0);
        let (signatures, _) = signatures(&desc);
        assert_eq!(signatures[0].ret, TypeExpr::Bool);
        assert_eq!(signatures[0].params.len(), 2);
    }

    #[test]
    fn a_request_takes_and_returns_terms() {
        let desc = OpDescriptor::new("demo_diff", Arity::Fixed(2)).flags(OpFlags::HEAVY);
        let (signatures, _) = signatures(&desc);
        assert_eq!(signatures.len(), 1, "a request has no numeric overload");
        assert_eq!(signatures[0].ret, term());
    }

    #[test]
    fn a_bound_variable_is_a_symbol_not_a_value() {
        let desc = OpDescriptor::new("demo_sum", Arity::Fixed(4))
            .flags(OpFlags::HEAVY)
            .binder(1, 1 << 2);
        let (signatures, _) = signatures(&desc);
        assert_eq!(
            signatures[0].params,
            vec![term(), TypeExpr::Symbol, term(), term()]
        );
    }

    #[test]
    fn a_variadic_operator_says_so() {
        let desc = OpDescriptor::new("demo_sum", Arity::Variadic).eval(|a| a.iter().sum());
        let (signatures, _) = signatures(&desc);
        for signature in &signatures {
            assert!(signature.variadic.is_some());
            assert!(signature.accepts_arity(7));
        }
    }

    #[test]
    fn a_leaf_has_no_signature_and_says_why() {
        let desc = OpDescriptor::new("demo_lit", Arity::Fixed(0)).flags(OpFlags::LEAF);
        assert_eq!(signatures(&desc).0.len(), 0);
        assert_eq!(reason(&desc), DynamicReason::OpaquePayload);
    }

    #[test]
    fn a_declared_overload_wins_over_the_derivation() {
        let declared = Signature::new(vec![TypeExpr::Int], TypeExpr::Bool);
        let desc = numeric("demo", Arity::Fixed(1)).signature(declared.clone());
        let (signatures, source) = signatures(&desc);
        assert_eq!(source, SignatureSource::Declared);
        assert_eq!(signatures, vec![declared]);
    }
}
