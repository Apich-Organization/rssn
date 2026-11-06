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

use crate::kernel::signature::DimExpr;
use crate::kernel::signature::NatExpr;
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

/// How a parameter's rank accepts the argument's.
///
/// Two literals are compared; anything else is unpinned and scores [`Match::Shaped`].
const fn rank_match(
    param: &NatExpr,
    got: &NatExpr,
) -> Match {
    match (param, got) {
        | (NatExpr::Lit(a), NatExpr::Lit(b)) => {
            if *a == *b {
                Match::Exact
            } else {
                Match::No
            }
        },
        | _ => Match::Shaped,
    }
}

/// How a parameter's extents accept the argument's.
fn shape_match(
    param: Option<&[NatExpr]>,
    got: Option<&[NatExpr]>,
) -> Match {
    let (Some(a), Some(b)) = (param, got) else {
        // One side does not pin the extents; the rank already compared.
        return Match::Shaped;
    };
    if a.len() != b.len() {
        return Match::No;
    }
    a.iter()
        .zip(b)
        .map(|(x, y)| rank_match(x, y))
        .fold(Match::Exact, weaker)
}

fn polymorphic_dimension(dim: &DimExpr) -> bool {
    match dim {
        | DimExpr::Mul(parts) => parts.iter().any(polymorphic_dimension),
        | DimExpr::Pow(inner, crate::kernel::signature::Exp::Rat(_)) | DimExpr::Inverse(inner) => {
            polymorphic_dimension(inner)
        },
        | DimExpr::Any | DimExpr::Var(_) | DimExpr::Pow(..) => true,
        | _ => false,
    }
}

/// Canonical factors for a fully concrete dimension. Named bases are kept as
/// exact identities; unknowns, symbolic exponents and invalid rationals never
/// acquire an invented value.
pub(crate) fn canonical_concrete_dimension(
    dim: &DimExpr
) -> Option<std::collections::BTreeMap<String, num_rational::BigRational>> {
    use num_rational::BigRational;
    use num_traits::{One, Zero};
    fn collect(
        dim: &DimExpr,
        power: &BigRational,
        factors: &mut std::collections::BTreeMap<String, BigRational>,
    ) -> Option<()> {
        match dim {
            | DimExpr::Dimless => {},
            | DimExpr::Base(name) => *factors.entry(name.clone()).or_default() += power.clone(),
            | DimExpr::Mul(parts) => {
                for part in parts {
                    collect(part, power, factors)?;
                }
            },
            | DimExpr::Inverse(inner) => collect(inner, &(-power.clone()), factors)?,
            | DimExpr::Pow(inner, crate::kernel::signature::Exp::Rat(r)) if r.denom != 0 => {
                let exponent = BigRational::new(r.numer.into(), r.denom.into());
                collect(inner, &(power * exponent), factors)?;
            },
            | _ => return None,
        }
        factors.retain(|_, exponent| !exponent.is_zero());
        Some(())
    }
    let mut factors = std::collections::BTreeMap::new();
    collect(dim, &BigRational::one(), &mut factors)?;
    Some(factors)
}

pub(crate) fn dimensions_equivalent(
    a: &DimExpr,
    b: &DimExpr,
) -> bool {
    if a == b && matches!(a, DimExpr::Dimless | DimExpr::Base(_)) {
        return true;
    }
    match (
        canonical_concrete_dimension(a),
        canonical_concrete_dimension(b),
    ) {
        | (Some(left), Some(right)) => left == right,
        | _ => false,
    }
}

/// How a parameter's dimension accepts the argument's.
///
/// `Any` and a variable accept anything; two concrete dimensions must be
/// equal. A length parameter does not accept a time argument.
fn dim_match(
    param: &DimExpr,
    got: &DimExpr,
) -> Match {
    match (param, got) {
        | (a, b) if polymorphic_dimension(a) || polymorphic_dimension(b) => Match::Generic,
        | (a, b) if a == b || dimensions_equivalent(a, b) => Match::Exact,
        | _ => Match::No,
    }
}

// A function's advertised schema is a promise, so unknown actual dimensions
// and extents cannot establish a concrete promise. `operator_on_left` tracks
// which side may contain the operator's existential variables as variance flips.
fn callback_schema_match(
    param: &TypeExpr,
    got: &TypeExpr,
    operator_on_left: bool,
) -> Match {
    fn nat(
        param: &NatExpr,
        got: &NatExpr,
        operator_on_left: bool,
    ) -> Match {
        let operator = if operator_on_left {
            param
        } else {
            got
        };
        let actual = if operator_on_left {
            got
        } else {
            param
        };
        if !matches!(actual, NatExpr::Lit(_)) {
            return Match::No;
        }
        if matches!(operator, NatExpr::Var(_)) {
            return Match::Generic;
        }
        match (param, got) {
            | (NatExpr::Lit(a), NatExpr::Lit(b)) if a == b => Match::Exact,
            | _ => Match::No,
        }
    }
    fn dimension(
        param: &DimExpr,
        got: &DimExpr,
        operator_on_left: bool,
    ) -> Match {
        let operator = if operator_on_left {
            param
        } else {
            got
        };
        let actual = if operator_on_left {
            got
        } else {
            param
        };
        if polymorphic_dimension(actual) {
            return Match::No;
        }
        if matches!(operator, DimExpr::Any) {
            return Match::No;
        }
        if polymorphic_dimension(operator) {
            return Match::Generic;
        }
        if dimensions_equivalent(param, got) {
            Match::Exact
        } else {
            Match::No
        }
    }
    match (param, got) {
        | (TypeExpr::Function(pa, ra), TypeExpr::Function(pb, rb)) => {
            if pa.len() != pb.len() {
                return Match::No;
            }
            let arguments = pa
                .iter()
                .zip(pb)
                .map(|(a, b)| callback_schema_match(b, a, !operator_on_left))
                .fold(Match::Exact, weaker);
            weaker(arguments, callback_schema_match(ra, rb, operator_on_left))
        },
        | (TypeExpr::Record(a), TypeExpr::Record(b)) => {
            if a.len() != b.len() {
                return Match::No;
            }
            let mut names = std::collections::BTreeSet::new();
            if a.iter().any(|(name, _)| !names.insert(name)) {
                return Match::No;
            }
            names.clear();
            if b.iter().any(|(name, _)| !names.insert(name)) {
                return Match::No;
            }
            a.iter()
                .map(|(name, field)| {
                    b.iter()
                        .find(|(other, _)| other == name)
                        .map_or(Match::No, |(_, other)| {
                            callback_schema_match(field, other, operator_on_left)
                        })
                })
                .fold(Match::Exact, weaker)
        },
        | (
            TypeExpr::Tensor {
                rank: ra,
                shape: sa,
                elem: ea,
            },
            TypeExpr::Tensor {
                rank: rb,
                shape: sb,
                elem: eb,
            },
        ) => {
            let shape = match (sa, sb) {
                | (None, _) => Match::Shaped,
                | (Some(a), Some(b)) if a.len() == b.len() => a
                    .iter()
                    .zip(b)
                    .map(|(a, b)| nat(a, b, operator_on_left))
                    .fold(Match::Exact, weaker),
                | _ => Match::No,
            };
            weaker(
                nat(ra, rb, operator_on_left),
                weaker(shape, callback_schema_match(ea, eb, operator_on_left)),
            )
        },
        | (TypeExpr::Dim(a, da), TypeExpr::Dim(b, db)) => weaker(
            dimension(da, db, operator_on_left),
            callback_schema_match(a, b, operator_on_left),
        ),
        | (TypeExpr::Unit(a), TypeExpr::Unit(b)) | (TypeExpr::Point(a), TypeExpr::Point(b)) => {
            dimension(a, b, operator_on_left)
        },
        | (TypeExpr::Dim(a, d), b) if dimensions_equivalent(d, &DimExpr::Dimless) => {
            callback_schema_match(a, b, operator_on_left)
        },
        | (a, TypeExpr::Dim(b, d)) if dimensions_equivalent(d, &DimExpr::Dimless) => {
            callback_schema_match(a, b, operator_on_left)
        },
        | (TypeExpr::Dim(a, d), b) if operator_on_left && matches!(d, DimExpr::Var(_)) => weaker(
            Match::Generic,
            callback_schema_match(a, b, operator_on_left),
        ),
        | (a, TypeExpr::Dim(b, d)) if !operator_on_left && matches!(d, DimExpr::Var(_)) => weaker(
            Match::Generic,
            callback_schema_match(a, b, operator_on_left),
        ),
        | (TypeExpr::Dyn(_), _) | (_, TypeExpr::Dyn(_)) => Match::No,
        | (TypeExpr::Var(_), _) if operator_on_left => Match::Generic,
        | (_, TypeExpr::Var(_)) if !operator_on_left => Match::Generic,
        | (TypeExpr::Var(_), _) | (_, TypeExpr::Var(_)) => Match::No,
        | _ => matches_type(param, got),
    }
}

/// How the type `param` accepts the type `got`.
#[must_use]
pub fn matches_type(
    param: &TypeExpr,
    got: &TypeExpr,
) -> Match {
    if matches!(
        (param, got),
        (TypeExpr::Function(..), TypeExpr::Function(..))
    ) {
        return callback_schema_match(param, got, true);
    }
    if let (TypeExpr::Record(a), TypeExpr::Record(b)) = (param, got) {
        if a.len() != b.len() {
            return Match::No;
        }
        let mut names = std::collections::BTreeSet::new();
        if a.iter().any(|(name, _)| !names.insert(name)) {
            return Match::No;
        }
        names.clear();
        if b.iter().any(|(name, _)| !names.insert(name)) {
            return Match::No;
        }
        return a
            .iter()
            .map(|(name, field)| {
                b.iter()
                    .find(|(other, _)| other == name)
                    .map_or(Match::No, |(_, other)| matches_type(field, other))
            })
            .fold(Match::Exact, weaker);
    }
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
        // down: a dimension, a rank, an element type. A pinned argument is
        // compared, not ignored: a rank-2 value must not match a rank-1
        // parameter, or the host's static lattice and this one disagree about
        // which overload the program chose.
        | (TypeExpr::Dim(a, da), TypeExpr::Dim(b, db)) => {
            weaker(dim_match(da, db), matches_type(a, b))
        },
        | (TypeExpr::Dim(inner, dim), other) if dimensions_equivalent(dim, &DimExpr::Dimless) => {
            matches_type(inner, other)
        },
        | (other, TypeExpr::Dim(inner, dim)) if dimensions_equivalent(dim, &DimExpr::Dimless) => {
            matches_type(other, inner)
        },
        // A unit is a dimension carrier: `quantity(value, unit)` declares its
        // unit parameter as `Unit(D)` and passes a `Value::Unit`, and without
        // this arm the typed overload scored `No` and lost to the all-`Dyn`
        // one.
        | (TypeExpr::Dim(_, da), TypeExpr::Unit(db) | TypeExpr::Point(db))
        | (TypeExpr::Unit(da) | TypeExpr::Point(da), TypeExpr::Dim(_, db)) => dim_match(da, db),
        | (TypeExpr::Unit(a), TypeExpr::Unit(b)) | (TypeExpr::Point(a), TypeExpr::Point(b)) => {
            dim_match(a, b)
        },
        // A dimension *variable* is polymorphic in the dimension: `Real[T]`
        // accepts a plain real (binding `T` to dimensionless) as well as any
        // dimensioned scalar. Without this arm a signature such as
        // `sde_wiener(Int, Int, Int, Real[T]) -> Tensor<2, Real[T^(1/2)]>`
        // could only accept a value whose type already *was* `Dim`, so the
        // dynamic call path refused every real argument. Concrete dimensions
        // are still compared by the arms above, and the element type is still
        // matched, so a `Real[T]` parameter does not accept a tensor.
        | (TypeExpr::Dim(inner, dim), _) if polymorphic_dimension(dim) => {
            weaker(Match::Generic, matches_type(inner, got))
        },
        // Contravariant function-input matching reverses the two types. An
        // operator's still-unbound call dimension is existential here, not a
        // demand that a concrete callback accept every physical dimension.
        | (_, TypeExpr::Dim(inner, dim)) if polymorphic_dimension(dim) => {
            weaker(Match::Generic, matches_type(param, inner))
        },
        | (
            TypeExpr::Tensor {
                rank: ra,
                shape: sa,
                elem: ea,
            },
            TypeExpr::Tensor {
                rank: rb,
                shape: sb,
                elem: eb,
            },
        ) => weaker(
            rank_match(ra, rb),
            weaker(
                shape_match(sa.as_deref(), sb.as_deref()),
                matches_type(ea, eb),
            ),
        ),
        | (
            TypeExpr::Sparse {
                format: fa,
                rank: ra,
                elem: ea,
            },
            TypeExpr::Sparse {
                format: fb,
                rank: rb,
                elem: eb,
            },
        ) => {
            let format = if fa == fb {
                Match::Exact
            } else {
                Match::No
            };
            weaker(format, weaker(rank_match(ra, rb), matches_type(ea, eb)))
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

/// The shared scoring fold: the *minimum* match over the parameters.
///
/// `match_at` classifies parameter `index` against the argument at that
/// position. Everything except that classification is identical for values
/// and for types, so both [`score`] and [`score_types`] go through here and
/// cannot disagree about arity, the variadic tail or the ranking.
fn score_fold<F>(
    signature: &Signature,
    arity: usize,
    mut match_at: F,
) -> Option<u32>
where
    F: FnMut(usize, &TypeExpr) -> Match,
{
    if !signature.accepts_arity(arity) {
        return None;
    }
    let mut worst = Match::Exact;
    for index in 0..arity {
        // A variadic tail repeats its element type.
        let param = signature
            .params
            .get(index)
            .or_else(|| signature.variadic.as_ref().map(|variadic| &variadic.elem))?;
        let matched = match_at(index, param);
        if !matched.accepted() {
            return None;
        }
        worst = weaker(worst, matched);
    }
    Some(worst.rank())
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
    score_fold(signature, args.len(), |index, param| {
        matches(param, &args[index])
    })
}

/// How well `signature` fits the static types `args`.
///
/// The type-level twin of [`score`]: a host that resolves an overload
/// statically gets the same number the dynamic call path would.
#[must_use]
pub fn score_types(
    signature: &Signature,
    args: &[TypeExpr],
) -> Option<u32> {
    score_fold(signature, args.len(), |index, param| {
        matches_type(param, &args[index])
    })
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
    choose(signatures.iter().map(|signature| score(signature, args)))
}

/// Chooses among `signatures` for the static types `args`.
///
/// The type-level twin of [`best`]: a host that resolves statically runs the
/// same lattice the dynamic call path runs, instead of rebuilding the fold
/// and letting the two drift apart.
#[must_use]
pub fn best_types(
    signatures: &[Signature],
    args: &[TypeExpr],
) -> Best {
    choose(
        signatures
            .iter()
            .map(|signature| score_types(signature, args)),
    )
}

/// The shared choice fold: keep the accepted scores, pick the highest, and
/// report a tie rather than guessing.
fn choose(scores: impl Iterator<Item = Option<u32>>) -> Best {
    let scored: Vec<(usize, u32)> = scores
        .enumerate()
        .filter_map(|(index, score)| score.map(|score| (index, score)))
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
    fn a_rank_mismatch_is_not_a_match() {
        // A rank-1 signature against a rank-2 value. The element type matches
        // exactly, but rank is not decoration: passing a matrix where a vector
        // is declared is a type error, not a looser match. Scoring this
        // `Shaped` let a host's static lattice and this one disagree about
        // which overload the program had chosen.
        let vector = Signature::new(vec![TypeExpr::Real.vector(3_u64)], TypeExpr::Real);
        let matrix = Value::Tensor(crate::kernel::value::TensorValue {
            shape: vec![2, 2],
            elem: TypeExpr::Real,
            data: vec![real(), real(), real(), real()],
        });
        assert_eq!(score(&vector, &[matrix]), None);
    }

    #[test]
    fn an_unpinned_rank_still_matches_any_rank() {
        // A signature that quantifies its rank accepts either, so the fix
        // above must not make every tensor signature rank-specific.
        let any = Signature::new(
            vec![TypeExpr::Tensor {
                rank: NatExpr::Var("r".to_owned()),
                shape: None,
                elem: Box::new(TypeExpr::Real),
            }],
            TypeExpr::Real,
        );
        let matrix = Value::Tensor(crate::kernel::value::TensorValue {
            shape: vec![2, 2],
            elem: TypeExpr::Real,
            data: vec![real(), real(), real(), real()],
        });
        assert_eq!(score(&any, &[matrix]), Some(2));
    }

    #[test]
    fn a_unit_matches_a_point_and_a_dimension_by_dimension() {
        // The unit-carrying arms. Before them, a unit-typed parameter that was
        // not *equal* to the argument scored `No`, so a `Real[Length]`
        // overload could never beat the all-`Dyn` fallback: `quantity(3, m)`
        // then ran the overload that loses the dimension. The association is
        // by measured dimension, and a different dimension is still refused.
        let length = DimExpr::base("Length");
        let time = DimExpr::base("Time");
        assert_eq!(
            matches_type(
                &TypeExpr::Real.with_dim(length.clone()),
                &TypeExpr::Unit(length.clone())
            ),
            Match::Exact
        );
        assert_eq!(
            matches_type(
                &TypeExpr::Real.with_dim(length.clone()),
                &TypeExpr::Point(length.clone())
            ),
            Match::Exact
        );
        assert_eq!(
            matches_type(
                &TypeExpr::Real.with_dim(length.clone()),
                &TypeExpr::Unit(time.clone())
            ),
            Match::No,
            "a length parameter must not accept a time unit"
        );
        assert_eq!(
            matches_type(
                &TypeExpr::Unit(length.clone()),
                &TypeExpr::Unit(length.clone())
            ),
            Match::Exact
        );
        assert_eq!(
            matches_type(
                &TypeExpr::Point(length.clone()),
                &TypeExpr::Point(length.clone())
            ),
            Match::Exact
        );
        assert_eq!(
            matches_type(&TypeExpr::Unit(length), &TypeExpr::Unit(time)),
            Match::No
        );
        // An unpinned unit parameter is generic over the dimension.
        assert_eq!(
            matches_type(
                &TypeExpr::Unit(DimExpr::Any),
                &TypeExpr::Unit(DimExpr::base("Length"))
            ),
            Match::Generic
        );
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
    fn a_dimension_variable_parameter_accepts_a_plain_scalar() {
        // `Real[T]` with a dimension variable binds `T`. It must accept a
        // plain real (`T = 1`), an exact integer and a dimensioned scalar —
        // otherwise a dim-generic overlay such as `sde_wiener`'s could not be
        // called through the dynamic ABI at all — but the element type is
        // still matched, so a tensor is not accepted.
        let generic = TypeExpr::Real.with_dim(DimExpr::Var("T".to_owned()));
        assert_eq!(matches_type(&generic, &TypeExpr::Real), Match::Generic);
        assert_eq!(matches_type(&generic, &TypeExpr::Int), Match::Generic);
        assert_eq!(
            matches_type(&generic, &TypeExpr::Real.with_dim(DimExpr::base("Time"))),
            Match::Generic
        );
        assert_eq!(
            matches_type(&generic, &TypeExpr::Real.vector(3_u64)),
            Match::No,
            "a dimensioned scalar must not accept a tensor"
        );
        // A concrete dimension is not polymorphic.
        let concrete = TypeExpr::Real.with_dim(DimExpr::base("Time"));
        assert_eq!(matches_type(&concrete, &TypeExpr::Real), Match::No);
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

    #[test]
    fn the_type_lattice_agrees_with_the_value_lattice() {
        // Every case the value-level `best` tests cover, run once over
        // `Value`s and once over the corresponding `TypeExpr`s.
        let exact = Signature::scalars(1, TypeExpr::Real);
        let generic = Signature::new(vec![TypeExpr::Var("T".to_owned())], TypeExpr::Real);
        let values = [real()];
        let types = [TypeExpr::Real];
        assert_eq!(
            best(&[generic.clone(), exact.clone()], &values),
            best_types(&[generic, exact], &types)
        );
        assert_eq!(
            best_types(&[Signature::scalars(1, TypeExpr::Real)], &types),
            Best::Unique(0)
        );

        // Widening: an `Int` prefers the `Int` parameter in both lattices.
        let real_param = Signature::new(vec![TypeExpr::Real], TypeExpr::Real);
        let int_param = Signature::new(vec![TypeExpr::Int], TypeExpr::Real);
        assert_eq!(
            best(
                &[real_param.clone(), int_param.clone()],
                &[Value::int(1_i32)]
            ),
            best_types(&[real_param, int_param], &[TypeExpr::Int])
        );

        // The tie case must be ambiguous in both, not silently broken.
        let a = Signature::scalars(1, TypeExpr::Real);
        let b = Signature::scalars(1, TypeExpr::Complex);
        let by_value = best(&[a.clone(), b.clone()], &values);
        let by_type = best_types(&[a, b], &types);
        assert_eq!(by_value, by_type);
        assert!(matches!(by_value, Best::Ambiguous(ref indices) if indices.len() == 2));

        // Arity filters first, and an impossible call has no candidate.
        let one = Signature::scalars(1, TypeExpr::Real);
        let two = Signature::scalars(2, TypeExpr::Real);
        assert_eq!(
            best(&[one.clone(), two.clone()], &values),
            best_types(&[one, two], &types)
        );
        let real_only = Signature::new(vec![TypeExpr::Real], TypeExpr::Real);
        assert_eq!(
            best(
                std::slice::from_ref(&real_only),
                &[Value::String("x".to_owned())]
            ),
            Best::None
        );
        assert_eq!(best_types(&[real_only], &[TypeExpr::String]), Best::None);
    }
}

#[cfg(test)]
mod concrete_dimension_equivalence_tests {
    use super::*;
    use crate::kernel::signature::{Exp, Rational};
    #[test]
    fn products_powers_cancellation_and_named_bases_match_canonically() {
        let length = DimExpr::Base("Length".into());
        let time = DimExpr::Base("Time".into());
        let square = DimExpr::Pow(Box::new(length.clone()), Exp::Rat(Rational::new(2, 1)));
        let product = DimExpr::Mul(vec![length.clone(), length.clone()]);
        let f = |dim| TypeExpr::Dim(Box::new(TypeExpr::Float), dim);
        assert_eq!(matches_type(&f(square), &f(product)), Match::Exact);
        let cancelled = DimExpr::Mul(vec![time.clone(), DimExpr::Inverse(Box::new(time))]);
        assert_eq!(matches_type(&f(cancelled), &TypeExpr::Float), Match::Exact);
        assert_eq!(
            matches_type(&f(length), &f(DimExpr::Base("length".into()))),
            Match::No
        );
        assert!(!dimensions_equivalent(&DimExpr::Any, &DimExpr::Dimless));
        assert!(!dimensions_equivalent(
            &DimExpr::Var("D".into()),
            &DimExpr::Dimless
        ));
    }
}
