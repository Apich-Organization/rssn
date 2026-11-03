//! # The dimensional layer
//!
//! This module is rssn's canonical view of dimensions and units, shared by
//! the run time and by anything that mirrors it statically.
//!
//! * [`DimTerm`] — a dimension as the engine computes with it: products of
//!   interned base dimensions and dimension variables with exact rational
//!   (or symbolic) exponents.
//! * [`DimArena`] — the interner that gives those bases and variables their
//!   ids, with snapshots so that a host can save and restore a numbering.
//! * [`Dimension`] — the owned, printable snapshot a caller sees.
//! * [`units`] — the one unit table: exact scales, affine offsets, aliases,
//!   prefixes, and a registry a host can extend.
//! * [`DimExpr`] — the *type-level* mirror of a dimension, written with
//!   names instead of ids, used in signatures and manifests.
//!
//! The implementation of the first three lives in
//! [`crate::graph::dim`] because dimension inference walks the term graph;
//! this module re-exports them so that a host depends on `rssn::dim` and
//! never on the graph kernel's internals.
//!
//! # Why names, not ids, at the boundary
//!
//! A [`BaseId`] is an index into one arena. Writing it to a manifest would
//! produce a document that only means something inside the process that
//! wrote it. Everything that leaves the process — a [`UnitDef`], a
//! signature, a manifest — therefore carries [`DimExpr`], which names the
//! bases, and converts to a [`DimTerm`] against the receiving arena.

pub mod constants;
pub mod units;

use num_rational::BigRational;

pub use crate::graph::dim::BaseId;
pub use crate::graph::dim::DimArena;
pub use crate::graph::dim::DimArenaSnapshot;
pub use crate::graph::dim::DimError;
pub use crate::graph::dim::DimOutcome;
pub use crate::graph::dim::DimRule;
pub use crate::graph::dim::DimTerm;
pub use crate::graph::dim::DimVarId;
pub use crate::graph::dim::Dimension;
pub use crate::graph::dim::DimensionMode;
pub use crate::graph::dim::SI_BASES;
pub use crate::graph::dim::UnitInfo;

pub use constants::PhysicalConstant;
pub use units::ResolvedUnit;
pub use units::Scale;
pub use units::StandardRegistry;
pub use units::UnitDef;
pub use units::UnitRegistry;
pub use units::UnitRegistrySnapshot;

pub use crate::kernel::signature::DimExpr;

use crate::kernel::signature::Exp as TypeExp;

impl DimExpr {
    /// This dimension as a term of `arena`.
    ///
    /// A base or variable name that the arena does not know is interned: the
    /// conversion is a *definition*, not a lookup, because the names in a
    /// manifest are the authority.
    ///
    /// # Errors
    /// Fails for [`DimExpr::Any`] and [`DimExpr::NonMetric`], which are not
    /// products of factors, and for a symbolic exponent on a product.
    pub fn to_term(
        &self,
        arena: &mut DimArena,
    ) -> Result<DimTerm, DimError> {
        let mut term = DimTerm::dimensionless();
        for (factor, exponent) in self.factors()? {
            // A rational exponent is the common case; a symbolic one (`x^n`)
            // is interned as an atom, which is what `DimTerm` stores.
            let rational = match &exponent {
                | TypeExp::Rat(r) => {
                    Some(BigRational::new(
                        num_bigint::BigInt::from(r.numer),
                        num_bigint::BigInt::from(r.denom),
                    ))
                },
                | TypeExp::Symbol(_) => None,
            };
            let symbol = match &exponent {
                | TypeExp::Symbol(name) => Some(name.clone()),
                | TypeExp::Rat(_) => None,
            };
            let factor_term = match factor {
                | Factor::Base(name) => {
                    let base = arena.define_base(&name);
                    match (&rational, &symbol) {
                        | (Some(power), _) => DimTerm::base(base, power.clone()),
                        | (None, Some(symbol)) => DimTerm::base_symbolic(base, symbol, arena),
                        | (None, None) => DimTerm::dimensionless(),
                    }
                },
                | Factor::Var(name) => {
                    let var = arena.intern_var(&name);
                    match (&rational, &symbol) {
                        | (Some(power), _) => DimTerm::var_pow(var, power.clone()),
                        | (None, Some(symbol)) => DimTerm::var_symbolic(var, symbol, arena),
                        | (None, None) => DimTerm::dimensionless(),
                    }
                },
            };
            term = term.mul(&factor_term);
        }
        Ok(term)
    }

    /// The dimension of `term` as a type-level expression.
    ///
    /// This is the exact inverse of [`DimExpr::to_term`] for every term that
    /// [`Exp`](crate::kernel::signature::Exp) can write: a rational exponent
    /// becomes `Exp::Rat`, an exponent of `-1` becomes [`DimExpr::Inverse`],
    /// and an exponent that is a single interned symbol becomes `Exp::Symbol`
    /// (`Length^n` stays `Length^n`).
    ///
    /// # Errors
    /// Fails for an exponent the type-level language cannot write — a
    /// compound symbolic exponent such as `Length^(2n)`, or a rational
    /// outside `i64` — instead of dropping the factor. The answer is checked
    /// by converting it back with [`DimExpr::to_term`] against a copy of
    /// `arena` and comparing terms, so this never returns a value that would
    /// rebuild a *different* term.
    ///
    /// The result is canonical rather than literally the input: factors are
    /// sorted by name with bases before variables, `Mul` is flat, and an
    /// exponent of `1` is the bare factor. A product written in that order
    /// round trips unchanged.
    pub fn try_from_term(
        term: &DimTerm,
        arena: &DimArena,
    ) -> Result<Self, DimError> {
        if term.is_any() {
            return Ok(Self::Any);
        }
        if term.is_container() {
            return Ok(Self::Container);
        }
        if term.is_non_metric() {
            return Ok(Self::NonMetric);
        }
        let view = term.to_dimension(arena);
        let mut factors: Vec<Self> =
            Vec::with_capacity(view.factors.len().saturating_add(view.unknown.len()));
        for (name, exponent) in view.factors {
            factors.push(factor_with_exponent(
                Self::base(&name),
                exponent_from_text(&exponent),
            ));
        }
        for (name, exponent) in view.unknown {
            factors.push(factor_with_exponent(
                Self::Var(name),
                exponent_from_text(&exponent),
            ));
        }
        let candidate = match factors.len() {
            | 0 => Self::Dimless,
            | 1 => factors.into_iter().next().unwrap_or(Self::Dimless),
            | _ => Self::Mul(factors),
        };
        // The printed exponent is the only window `DimTerm` offers onto an
        // exponent, and a compound one (`2n`) prints the same shape as a
        // symbol named `2*n`. Rebuilding the term is what tells them apart.
        let rebuilt = candidate.to_term(&mut arena.clone())?;
        if rebuilt == *term {
            return Ok(candidate);
        }
        let reason = match term.symbolic_factors(arena).into_iter().next() {
            | Some((name, exponent)) => {
                format!(
                    "the exponent {exponent} of {name} is not a single symbol, so Exp cannot write it"
                )
            },
            | None => String::from("the exponents of this term do not fit Exp's rational range"),
        };
        Err(DimError::Parse {
            text: term.describe(arena),
            reason,
        })
    }

    /// The dimension of `term` as a type-level expression, or [`DimExpr::Any`].
    ///
    /// The total form of [`DimExpr::try_from_term`]: the same value whenever
    /// that function succeeds. A term whose exponent `Exp` cannot write
    /// (`Length^(2n)`) becomes [`DimExpr::Any`] — "no type-level dimension is
    /// known" — never a wrong dimension such as `Length^n`. `Any` is not the
    /// input, and [`DimExpr::to_term`] refuses it, so converting back is an
    /// error rather than a quiet lie. Use `try_from_term` to see that failure
    /// where it happens.
    #[must_use]
    pub fn from_term(
        term: &DimTerm,
        arena: &DimArena,
    ) -> Self {
        Self::try_from_term(term, arena).unwrap_or(Self::Any)
    }

    /// The multiplicative factors of the expression, with their exponents.
    ///
    /// # Errors
    /// Fails for the three non-scalar variants.
    pub fn factors(&self) -> Result<Vec<(Factor, TypeExp)>, DimError> {
        let mut out = Vec::new();
        self.collect(&mut out, TypeExp::from(1_i64))?;
        Ok(out)
    }

    /// Flattens `self` into `out`, multiplying `exponent` into each factor's.
    ///
    /// The sign of a factor is carried by its exponent, so no separate
    /// polarity is threaded through the recursion.
    fn collect(
        &self,
        out: &mut Vec<(Factor, TypeExp)>,
        exponent: TypeExp,
    ) -> Result<(), DimError> {
        match self {
            | Self::Dimless => Ok(()),
            | Self::Base(name) => {
                out.push((Factor::Base(name.clone()), exponent));
                Ok(())
            },
            | Self::Var(name) => {
                out.push((Factor::Var(name.clone()), exponent));
                Ok(())
            },
            | Self::Mul(factors) => {
                for factor in factors {
                    factor.collect(out, exponent.clone())?;
                }
                Ok(())
            },
            | Self::Pow(base, inner) => {
                let combined = match (&exponent, inner) {
                    | (TypeExp::Rat(a), TypeExp::Rat(b)) => {
                        let n = i128::from(a.numer) * i128::from(b.numer);
                        let d = i128::from(a.denom) * i128::from(b.denom);
                        TypeExp::Rat(crate::kernel::signature::Rational::new(
                            i64::try_from(n).unwrap_or(i64::MAX),
                            i64::try_from(d).unwrap_or(1),
                        ))
                    },
                    // `(Length^n)^1` is `Length^n`. A symbolic exponent
                    // cannot be *scaled* by a rational — `Exp` has no product
                    // of a symbol and a rational — so anything else is refused
                    // with a message that says which case it was.
                    | (TypeExp::Rat(a), TypeExp::Symbol(symbol)) if a.numer == 1 && a.denom == 1 => {
                        TypeExp::Symbol(symbol.clone())
                    },
                    | (TypeExp::Rat(_), TypeExp::Symbol(_)) => {
                        return Err(DimError::Parse {
                            text: crate::kernel::signature::dim_name(self),
                            reason: String::from(
                                "a symbolic exponent cannot be scaled by a rational",
                            ),
                        });
                    },
                    | _ => {
                        return Err(DimError::Parse {
                            text: crate::kernel::signature::dim_name(self),
                            reason: String::from("nested symbolic exponents are not representable"),
                        });
                    },
                };
                base.collect(out, combined)
            },
            | Self::Inverse(base) => {
                let negated = match &exponent {
                    | TypeExp::Rat(r) => {
                        TypeExp::Rat(crate::kernel::signature::Rational::new(-r.numer, r.denom))
                    },
                    | TypeExp::Symbol(s) => TypeExp::Symbol(format!("-{s}")),
                };
                base.collect(out, negated)
            },
            | Self::Any | Self::NonMetric | Self::Container => {
                Err(DimError::Parse {
                    text: crate::kernel::signature::dim_name(self),
                    reason: String::from("not a product of factors"),
                })
            },
        }
    }
}

/// A factor of a type-level dimension: a base dimension or a variable.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Factor {
    /// A named base dimension.
    Base(String),
    /// A dimension variable.
    Var(String),
}

/// A factor raised to a type-level exponent.
///
/// `-1` is the reciprocal, which is the shape [`DimExpr::Inverse`] exists
/// for; any other exponent is a [`DimExpr::Pow`].
fn factor_with_exponent(
    factor: DimExpr,
    exponent: TypeExp,
) -> DimExpr {
    match exponent {
        | TypeExp::Rat(r) if r.numer == -1 && r.denom == 1 => DimExpr::Inverse(Box::new(factor)),
        | exp => factor.pow(exp),
    }
}

/// A printed exponent as a type-level exponent.
///
/// The printed form of a rational is `p`, `-p`, `p/q` or `-p/q`; anything
/// else is a symbol name, or a compound such as `2*n` that is kept as a
/// symbol here and refused by the caller's round-trip check.
fn exponent_from_text(text: &str) -> TypeExp {
    let (numer, denom) = match text.split_once('/') {
        | Some((numer, denom)) => (numer, denom),
        | None => (text, "1"),
    };
    match (numer.parse::<i64>(), denom.parse::<i64>()) {
        | (Ok(numer), Ok(denom)) if denom != 0 => {
            TypeExp::Rat(crate::kernel::signature::Rational::new(numer, denom))
        },
        | _ => TypeExp::Symbol(text.to_owned()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use num_bigint::BigInt;

    /// `expr` through `to_term` and back, on a fresh arena.
    fn round_trip(expr: &DimExpr) -> DimExpr {
        let mut arena = DimArena::new();
        let term = expr.to_term(&mut arena).unwrap_or_else(|e| panic!("{e}"));
        DimExpr::from_term(&term, &arena)
    }

    #[test]
    fn from_term_is_the_inverse_of_to_term() {
        let cases = [
            DimExpr::Dimless,
            DimExpr::base("Length"),
            DimExpr::Mul(vec![
                DimExpr::base("Length"),
                DimExpr::Inverse(Box::new(DimExpr::base("Time"))),
            ]),
            DimExpr::base("Length")
                .pow(TypeExp::Rat(crate::kernel::signature::Rational::new(1, 2))),
            DimExpr::Inverse(Box::new(DimExpr::base("Time"))),
            DimExpr::base("Length").pow(TypeExp::Symbol(String::from("n"))),
            DimExpr::Var(String::from("D")).pow(TypeExp::Symbol(String::from("n"))),
        ];
        for case in cases {
            assert_eq!(round_trip(&case), case, "{case:?}");
        }
    }

    #[test]
    fn inverse_comes_back_as_inverse() {
        let reciprocal = DimExpr::Inverse(Box::new(DimExpr::base("Time")));
        assert_eq!(round_trip(&reciprocal), reciprocal);
        // Not the old rewrite of `Inverse(x)` into `Pow(x, -1)`.
        assert_ne!(
            round_trip(&reciprocal),
            DimExpr::base("Time").pow(TypeExp::from(-1_i64))
        );
    }

    #[test]
    fn a_single_symbol_exponent_survives_and_a_compound_one_is_refused() {
        let mut arena = DimArena::new();
        let symbolic = DimExpr::base("Length").pow(TypeExp::Symbol(String::from("n")));
        let term = symbolic
            .to_term(&mut arena)
            .unwrap_or_else(|e| panic!("{e}"));
        // `Length^n` is one interned atom, and `Exp::Symbol` writes it exactly.
        assert_eq!(
            DimExpr::try_from_term(&term, &arena).unwrap_or_else(|e| panic!("{e}")),
            symbolic
        );

        // `Length^n * Length^n` is `Length^(2n)`. `Exp` has no product of a
        // symbol and a rational, so the honest outcome is a refusal — not a
        // silent `Length^n`, which is what the old code returned.
        let squared = term.mul(&term);
        assert!(DimExpr::try_from_term(&squared, &arena).is_err());
        assert_ne!(DimExpr::from_term(&squared, &arena), symbolic);
        assert_eq!(DimExpr::from_term(&squared, &arena), DimExpr::Any);
    }

    #[test]
    fn to_term_output_is_unchanged() {
        let mut arena = DimArena::new();
        let length = arena.define_base("Length");
        let time = arena.define_base("Time");
        let one = BigRational::new(BigInt::from(1), BigInt::from(1));
        let minus_one = BigRational::new(BigInt::from(-1), BigInt::from(1));

        let single = DimExpr::base("Length")
            .to_term(&mut arena)
            .unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(single, DimTerm::base(length, one.clone()));

        let product = DimExpr::Mul(vec![
            DimExpr::base("Length"),
            DimExpr::base("Time").pow(TypeExp::from(-1_i64)),
        ])
        .to_term(&mut arena)
        .unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(
            product.factors(&arena),
            vec![
                (String::from("Length"), one.clone()),
                (String::from("Time"), minus_one.clone()),
            ]
        );
        assert_eq!(
            product,
            DimTerm::base(length, one.clone()).mul(&DimTerm::base(time, minus_one))
        );

        let half = DimExpr::base("Length")
            .pow(TypeExp::Rat(crate::kernel::signature::Rational::new(1, 2)))
            .to_term(&mut arena)
            .unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(
            half,
            DimTerm::base(length, BigRational::new(BigInt::from(1), BigInt::from(2)))
        );

        let symbolic = DimExpr::base("Length")
            .pow(TypeExp::Symbol(String::from("n")))
            .to_term(&mut arena)
            .unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(
            symbolic.symbolic_factors(&arena),
            vec![(String::from("Length"), String::from("n"))]
        );
        assert_eq!(symbolic.describe(&arena), "Length^n");

        assert!(DimExpr::Any.to_term(&mut arena).is_err());
        assert!(DimExpr::NonMetric.to_term(&mut arena).is_err());
        assert!(DimExpr::Container.to_term(&mut arena).is_err());
    }

    #[test]
    fn symbolic_exponents_survive_the_bridge() {
        // `Length^n` is representable in `DimTerm` as an interned atom; the
        // bridge used to refuse it, which meant a static dimension with a
        // symbolic exponent could not be declared at all.
        let mut arena = DimArena::new();
        let expr = DimExpr::base("Length").pow(TypeExp::Symbol(String::from("n")));
        let term = expr.to_term(&mut arena).unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(term.describe(&arena), "Length^n");
        // Two occurrences of the same symbol are the same atom.
        let again = expr.to_term(&mut arena).unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(term, again);
    }

    #[test]
    fn non_metric_dimensions_do_not_become_terms() {
        assert!(DimExpr::Any.to_term(&mut DimArena::new()).is_err());
        assert!(DimExpr::NonMetric.to_term(&mut DimArena::new()).is_err());
    }

    #[test]
    fn user_bases_are_defined_on_demand() {
        let mut arena = DimArena::new();
        let term = DimExpr::base("Currency")
            .to_term(&mut arena)
            .unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(term.describe(&arena), "Currency");
    }
}
