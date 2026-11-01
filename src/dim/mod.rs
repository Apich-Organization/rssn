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
use num_traits::Zero;

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
                | TypeExp::Rat(r) => Some(BigRational::new(
                    num_bigint::BigInt::from(r.numer),
                    num_bigint::BigInt::from(r.denom),
                )),
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
    #[must_use]
    pub fn from_term(
        term: &DimTerm,
        arena: &DimArena,
    ) -> Self {
        if term.is_any() {
            return Self::Any;
        }
        if term.is_container() {
            return Self::Container;
        }
        if term.is_non_metric() {
            return Self::NonMetric;
        }
        let mut factors = Vec::new();
        for (name, exponent) in term.factors(arena) {
            factors.push(Self::base(&name).pow(rat_exp(&exponent)));
        }
        for (name, exponent) in term.vars(arena) {
            factors.push(Self::Var(name).pow(rat_exp(&exponent)));
        }
        match factors.len() {
            | 0 => Self::Dimless,
            | 1 => factors.into_iter().next().unwrap_or(Self::Dimless),
            | _ => Self::Mul(factors),
        }
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
                    | (TypeExp::Rat(a), TypeExp::Symbol(symbol))
                        if a.numer == 1 && a.denom == 1 =>
                    {
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
                    | TypeExp::Rat(r) => TypeExp::Rat(crate::kernel::signature::Rational::new(
                        -r.numer, r.denom,
                    )),
                    | TypeExp::Symbol(s) => TypeExp::Symbol(format!("-{s}")),
                };
                base.collect(out, negated)
            },
            | Self::Any | Self::NonMetric | Self::Container => Err(DimError::Parse {
                text: crate::kernel::signature::dim_name(self),
                reason: String::from("not a product of factors"),
            }),
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

/// A rational exponent as a type-level exponent.
fn rat_exp(value: &BigRational) -> TypeExp {
    if value.is_zero() {
        return TypeExp::from(0_i64);
    }
    let numer = i64::try_from(value.numer().clone()).unwrap_or(0);
    let denom = i64::try_from(value.denom().clone()).unwrap_or(1);
    TypeExp::Rat(crate::kernel::signature::Rational::new(numer, denom))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn type_level_dimensions_round_trip_through_terms() {
        let mut arena = DimArena::new();
        let cases = [
            DimExpr::Dimless,
            DimExpr::base("Length"),
            DimExpr::Mul(vec![
                DimExpr::base("Length"),
                DimExpr::base("Time").pow(TypeExp::from(-1_i64)),
            ]),
            DimExpr::base("Length").pow(TypeExp::Rat(
                crate::kernel::signature::Rational::new(1, 2),
            )),
        ];
        for case in cases {
            let term = case.to_term(&mut arena).unwrap_or_else(|e| panic!("{e}"));
            assert_eq!(DimExpr::from_term(&term, &arena), case, "{case:?}");
        }
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
