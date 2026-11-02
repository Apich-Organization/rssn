//! Choosing between the overloads of one operator.
//!
//! # Why this exists
//!
//! `OpTable` lets several registrations share a name, and each carries its own
//! [`Signature`]. The dynamic call path used to pick the first one whose
//! *arity* matched, which is fine when the overloads differ in arity and wrong
//! when they differ in type: a generated typed front end can select
//! `Real -> Real` statically while the runtime executes `Tensor -> Tensor`.
//!
//! This module answers the question the runtime actually has — *which overload
//! do these arguments belong to* — and reports an ambiguity instead of
//! guessing.
//!
//! # The match lattice
//!
//! Every parameter is matched against an argument and classified:
//!
//! | Class | Meaning |
//! |---|---|
//! | `Exact` | the parameter type *is* the argument's type |
//! | `Widened` | a numeric widening: `Int` into `Real`, `Real` into `Complex` |
//! | `Shaped` | the same constructor with a looser argument (a dimension or an element type the signature does not pin down) |
//! | `Generic` | a type variable or `Dyn`: it accepts anything |
//! | `No` | the argument cannot be passed |
//!
//! A call is resolved by maximising the *least* match over the parameters, so
//! an overload that is exact everywhere beats one that is exact in one place
//! and generic in another. If two overloads tie, the call is ambiguous and is
//! refused — silently choosing one is how a typed caller and the runtime drift
//! apart.

use crate::kernel::signature::Signature;
use crate::kernel::signature::TypeExpr;
use crate::kernel::value::Value;

/// How one parameter accepts one argument.
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Match {
    /// The argument cannot be passed to this parameter.
    No,
    /// The parameter is a type variable or `Dyn`.
    Generic,
    /// The same constructor, with an argument the signature does not pin down.
    Shaped,
    /// A numeric widening: an `Int` passed where a `Real` is wanted.
    Widened,
    /// The parameter type is exactly the argument's type.
    Exact,
}

impl Match {
    /// Whether the argument can be passed at all.
    #[must_use]
    pub const fn accepted(self) -> bool {
        !matches!(self, Self::No)
    }

    /// The rank used to compare two overloads.
    #[must_use]
    const fn rank(self) -> u32 {
        match self {
            | Self::No => 0,
            | Self::Generic => 1,
            | Self::Shaped => 2,
            | Self::Widened => 3,
            | Self::Exact => 4,
        }
    }
}

/// The weaker of two matches, used to combine a constructor with its argument.
const fn weaker(
    a: Match,
    b: Match,
) -> Match {
    if a.rank() < b.rank() { a } else { b }
}

/// How the type `param` accepts the type `got`.
#[must_use]
pub fn matches_type(
    param: &TypeExpr,
    got: &TypeExpr,
) -> Match {
    if param == got {
        return Match::Exact;
    }
    match (param, got) {
        // A type variable or a `Dyn` accepts anything.
        | (TypeExpr::Var(_) | TypeExpr::Dyn(_), _) | (_, TypeExpr::Dyn(_)) => Match::Generic,
        // Numeric widening, one direction only: an `Int` is a `Real`, a
        // `Real` is a `Complex`, and none of the three narrows.
        | (TypeExpr::Real, TypeExpr::Int | TypeExpr::Float)
        | (TypeExpr::Complex, TypeExpr::Int | TypeExpr::Real | TypeExpr::Float) => Match::Widened,
        // The same constructor, with an argument the signature does not pin
        // down: a dimension, a rank, an element type.
        | (TypeExpr::Dim(a, _), TypeExpr::Dim(b, _))
        | (TypeExpr::Tensor { elem: a, .. }, TypeExpr::Tensor { elem: b, .. })
        | (TypeExpr::Sparse { elem: a, .. }, TypeExpr::Sparse { elem: b, .. }) => {
            weaker(Match::Shaped, matches_type(a, b))
        },
        | (TypeExpr::Function(pa, ra), TypeExpr::Function(pb, rb)) => {
            let params = pa
                .iter()
                .zip(pb)
                .map(|(x, y)| matches_type(x, y))
                .fold(Match::Exact, weaker);
            weaker(params, matches_type(ra, rb))
        },
        | _ => Match::No,
    }
}

/// How `param` accepts `value`.
#[must_use]
pub fn matches(
    param: &TypeExpr,
    value: &Value,
) -> Match {
    matches_type(param, &value.type_expr())
}

/// How well `signature` fits `args`.
///
/// The score is the *minimum* match over the parameters, so an overload is
/// only as good as its worst parameter. `None` means the arguments cannot be
/// passed at all.
#[must_use]
pub fn score(
    signature: &Signature,
    args: &[Value],
) -> Option<u32> {
    if !signature.accepts_arity(args.len()) {
        return None;
    }
    let mut worst = Match::Exact;
    for (index, value) in args.iter().enumerate() {
        // A variadic tail repeats its element type.
        let param = signature
            .params
            .get(index)
            .or_else(|| signature.variadic.as_ref().map(|variadic| &variadic.elem))?;
        let matched = matches(param, value);
        if !matched.accepted() {
            return None;
        }
        worst = weaker(worst, matched);
    }
    Some(worst.rank())
}

/// The outcome of choosing among the overloads of one operator.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Best {
    /// No overload accepts these arguments.
    None,
    /// Exactly one overload is the most specific.
    Unique(usize),
    /// Several overloads tie for most specific; their indices.
    Ambiguous(Vec<usize>),
}

/// Chooses among `signatures` for `args`.
#[must_use]
pub fn best(
    signatures: &[Signature],
    args: &[Value],
) -> Best {
    let scored: Vec<(usize, u32)> = signatures
        .iter()
        .enumerate()
        .filter_map(|(index, signature)| score(signature, args).map(|score| (index, score)))
        .collect();
    if scored.is_empty() {
        return Best::None;
    }
    let top = scored.iter().map(|(_, score)| *score).max().unwrap_or(0);
    let mut tied: Vec<usize> = scored
        .iter()
        .filter(|(_, score)| *score == top)
        .map(|(index, _)| *index)
        .collect();
    match tied.len() {
        | 1 => tied.pop().map_or(Best::None, Best::Unique),
        | _ => Best::Ambiguous(tied),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dim::UnitDef;
    use crate::kernel::signature::DimExpr;
    use crate::kernel::value::QuantityValue;
    use num_bigint::BigInt;
    use num_rational::BigRational;

    /// A `Real`: an exact rational, not a `Float`, so the parameter matches
    /// exactly rather than by widening.
    fn real() -> Value {
        Value::rational(BigRational::new(BigInt::from(1), BigInt::from(2)))
    }

    #[test]
    fn an_exact_parameter_beats_a_generic_one() {
        let exact = Signature::scalars(1, TypeExpr::Real);
        let generic = Signature::new(vec![TypeExpr::Var("T".to_owned())], TypeExpr::Real);
        let args = [real()];
        assert_eq!(score(&exact, &args), Some(4));
        assert_eq!(score(&generic, &args), Some(1));
        assert_eq!(best(&[generic, exact], &args), Best::Unique(1));
    }

    #[test]
    fn an_integer_widens_into_a_real_but_not_the_other_way() {
        // `Signature::scalars` always uses `Real` parameters, so the `Int`
        // parameter has to be spelled out.
        let real_param = Signature::new(vec![TypeExpr::Real], TypeExpr::Real);
        let int_param = Signature::new(vec![TypeExpr::Int], TypeExpr::Real);
        let int = [Value::int(1_i32)];
        assert_eq!(score(&real_param, &int), Some(3), "an int widens to a real");
        assert_eq!(score(&int_param, &int), Some(4), "an int is an int");
        let float = [Value::Float(1.0)];
        assert_eq!(score(&int_param, &float), None, "a float is not an int");
        assert_eq!(score(&real_param, &float), Some(3), "a float is a real");
    }

    #[test]
    fn a_dimension_pinned_by_the_signature_matches_exactly() {
        let length = Signature::new(
            vec![TypeExpr::Real.with_dim(DimExpr::base("Length"))],
            TypeExpr::Real,
        );
        let unit =
            UnitDef::linear("m", DimExpr::base("Length"), crate::dim::Scale::one()).to_value();
        let quantity = Value::Quantity(QuantityValue {
            value: Box::new(real()),
            unit,
        });
        assert_eq!(score(&length, &[quantity]), Some(4));
        // A plain real is not a length.
        assert_eq!(score(&length, &[real()]), None);
    }

    #[test]
    fn a_signature_that_leaves_the_element_type_open_is_generic() {
        // `T[Length]` pins the dimension but not the element type, so the
        // overload is only as specific as its weakest parameter.
        let open = Signature::new(
            vec![TypeExpr::Var("T".to_owned()).with_dim(DimExpr::base("Length"))],
            TypeExpr::Real,
        );
        let unit =
            UnitDef::linear("m", DimExpr::base("Length"), crate::dim::Scale::one()).to_value();
        let quantity = Value::Quantity(QuantityValue {
            value: Box::new(real()),
            unit,
        });
        assert_eq!(score(&open, &[quantity]), Some(1));
    }

    #[test]
    fn the_same_constructor_with_a_different_shape_is_shaped() {
        // A rank-1 tensor signature against a rank-2 value: the element type
        // matches exactly but the shape does not, which is `Shaped`.
        let vector = Signature::new(vec![TypeExpr::Real.vector(3_u64)], TypeExpr::Real);
        let matrix = Value::Tensor(crate::kernel::value::TensorValue {
            shape: vec![2, 2],
            data: vec![real(), real(), real(), real()],
        });
        assert_eq!(score(&vector, &[matrix]), Some(2));
    }

    #[test]
    fn a_tie_is_ambiguous_rather_than_silently_broken() {
        let a = Signature::scalars(1, TypeExpr::Real);
        let b = Signature::scalars(1, TypeExpr::Complex);
        assert!(matches!(
            best(&[a, b], &[real()]),
            Best::Ambiguous(ref indices) if indices.len() == 2
        ));
    }

    #[test]
    fn arity_still_filters_first() {
        let one = Signature::scalars(1, TypeExpr::Real);
        let two = Signature::scalars(2, TypeExpr::Real);
        assert_eq!(best(&[one, two], &[real()]), Best::Unique(0));
    }

    #[test]
    fn a_variadic_tail_repeats_its_element_type() {
        let variadic =
            Signature::scalars(0, TypeExpr::Real).variadic(crate::kernel::signature::Variadic {
                min: 0,
                max: None,
                elem: TypeExpr::Real,
            });
        assert_eq!(score(&variadic, &[real(), real(), real()]), Some(4));
        assert_eq!(
            score(&variadic, &[real(), Value::String("x".to_owned())]),
            None
        );
    }
}
