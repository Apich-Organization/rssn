//! Signature inference for operators that do not declare one.
//!
//! A thousand operators do not get a thousand hand-written signatures on day
//! one, and leaving them untyped would make the manifest useless for a host
//! that wants to generate bindings. So an operator without a declared
//! signature gets one inferred from what the descriptor already says:
//!
//! * a predicate is `(Scalar, ...) -> Bool`;
//! * a heavy operator — an unevaluated request — returns a term, so it is
//!   `(Scalar, ...) -> Dyn<Term>`;
//! * an operator with a scalar evaluator is `(Scalar, ...) -> Scalar`;
//! * a variadic one is variadic over `Scalar`.
//!
//! Inference is deliberately *conservative*: it never guesses that an
//! argument is a tensor, because a wrong tensor shape is worse than a
//! conservative `Scalar`. Operators whose arguments really are tensors are
//! corrected in [`signatures.overlay.json`], which is applied first and
//! always wins.
//!
//! Every inferred signature is marked
//! [`SignatureSource::Inferred`](super::manifest::SignatureSource), so a host
//! can tell a declaration from a guess.

use crate::graph::OpDescriptor;
use crate::graph::OpFlags;
use crate::kernel::manifest::SignatureSource;
use crate::kernel::manifest::overlay;
use crate::kernel::signature::DimExpr;
use crate::kernel::signature::RepExpr;
use crate::kernel::signature::Signature;
use crate::kernel::signature::TypeExpr;
use crate::kernel::signature::Variadic;

/// Every signature of `desc`, and where they came from.
///
/// This is the single answer to "what does this operator accept", used by the
/// manifest, by `Session::overloads` and by `Session::resolve_overload`. A
/// declared overload wins; otherwise the overlay and the inference rules
/// supply one; otherwise the operator is untyped and the manifest records a
/// [`DynamicReason`](crate::kernel::signature::DynamicReason) instead.
#[must_use]
pub fn signatures(desc: &OpDescriptor) -> (Vec<Signature>, SignatureSource) {
    if !desc.overloads.is_empty() {
        return (desc.overloads.clone(), SignatureSource::Declared);
    }
    match infer_signature_with_source(desc) {
        | Some((signature, source)) => (vec![signature], source),
        | None => (Vec::new(), SignatureSource::None),
    }
}

/// The signature of `desc`, from the overlay or inferred.
///
/// Returns `None` when nothing can be said: an operator with no evaluator,
/// no heavy flag and no declared overload is reduced by rules alone, and the
/// type of its result depends on which rule fires.
#[must_use]
pub fn infer_signature(desc: &OpDescriptor) -> Option<Signature> {
    infer_signature_with_source(desc).map(|(signature, _)| signature)
}

/// The signature of `desc`, and where it came from.
///
/// The distinction matters to a host: an overlay entry was written by a
/// person and can be trusted the way a declaration can, while a signature
/// derived from an `OpFlags::PREDICATE` bit is a guess that happens to be
/// right.
#[must_use]
pub fn infer_signature_with_source(
    desc: &OpDescriptor
) -> Option<(Signature, SignatureSource)> {
    if let Some(declared) = overlay().get(&*desc.name) {
        return declared
            .first()
            .cloned()
            .map(|signature| (signature, SignatureSource::Overlay));
    }
    if desc.flags.has(OpFlags::PREDICATE) {
        return Some((
            Signature::scalars(scalar_count(desc), TypeExpr::Bool),
            SignatureSource::Inferred,
        ));
    }
    if desc.flags.has(OpFlags::HEAVY) {
        return Some((
            Signature::scalars(
                scalar_count(desc),
                TypeExpr::Dyn(RepExpr {
                    repr: "term".to_owned(),
                    dim: Some(DimExpr::Any),
                }),
            ),
            SignatureSource::Inferred,
        ));
    }
    if desc.eval_dyn.is_some() {
        return Some((
            Signature::scalars(scalar_count(desc), TypeExpr::Real),
            SignatureSource::Inferred,
        ));
    }
    if desc.eval.is_some() {
        let signature = match desc.arity {
            | crate::graph::Arity::Fixed(n) => {
                Signature::scalars(usize::from(n), TypeExpr::Real)
            },
            | crate::graph::Arity::Variadic => Signature::scalars(0, TypeExpr::Real).variadic(
                Variadic {
                    min: 0,
                    max: None,
                    elem: TypeExpr::Real,
                },
            ),
        };
        return Some((signature, SignatureSource::Inferred));
    }
    None
}

/// The number of parameters to assume for a fixed-arity operator.
const fn scalar_count(desc: &OpDescriptor) -> usize {
    match desc.arity {
        | crate::graph::Arity::Fixed(n) => n as usize,
        | crate::graph::Arity::Variadic => 0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Arity;

    #[test]
    fn evaluators_infer_scalar_signatures() {
        let desc = OpDescriptor::new("demo_sin", Arity::Fixed(1)).eval(|a| {
            a.first().copied().unwrap_or(f64::NAN)
        });
        let signature = infer_signature(&desc).unwrap_or_else(|| panic!("no signature"));
        assert_eq!(signature.params, vec![TypeExpr::Real]);
        assert_eq!(signature.ret, TypeExpr::Real);
    }

    #[test]
    fn predicates_and_requests_infer_their_result_type() {
        let predicate = OpDescriptor::new("demo_lt", Arity::Fixed(2))
            .flags(OpFlags::PREDICATE)
            .eval(|_| 1.0);
        assert_eq!(
            infer_signature(&predicate).map(|s| s.ret),
            Some(TypeExpr::Bool)
        );
        let request = OpDescriptor::new("demo_diff", Arity::Fixed(2)).flags(OpFlags::HEAVY);
        assert!(matches!(
            infer_signature(&request).map(|s| s.ret),
            Some(TypeExpr::Dyn(_))
        ));
    }

    #[test]
    fn a_rule_only_operator_has_no_signature() {
        let desc = OpDescriptor::new("demo_rewrite_only", Arity::Fixed(1));
        assert!(infer_signature(&desc).is_none());
    }

    #[test]
    fn variadic_evaluators_are_variadic() {
        let desc = OpDescriptor::new("demo_sum", Arity::Variadic).eval(|a| a.iter().sum());
        let signature = infer_signature(&desc).unwrap_or_else(|| panic!("no signature"));
        assert!(signature.variadic.is_some());
        assert!(signature.accepts_arity(0));
        assert!(signature.accepts_arity(7));
    }
}
