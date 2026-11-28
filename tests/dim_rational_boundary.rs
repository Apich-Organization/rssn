//! Regression tests for the runtime to type-level dimension conversion.

use num_bigint::BigInt;
use num_rational::BigRational;
use rssn::dim::BaseId;
use rssn::dim::DimArena;
use rssn::dim::DimExpr;
use rssn::dim::DimTerm;
use rssn::dim::DimVarId;
use rssn::dim::Dimension;
use rssn::kernel::signature::Exp;
use rssn::kernel::signature::Rational;

fn ratio(
    numer: i64,
    denom: i64,
) -> BigRational {
    BigRational::new(BigInt::from(numer), BigInt::from(denom))
}

#[test]
fn dimensional_rational_boundary_canonicalizes_each_kind_by_name() {
    let mut arena = DimArena::new();
    let input = DimExpr::Mul(vec![
        DimExpr::base("Zeta").pow(Exp::Rat(Rational::new(-2, 3))),
        DimExpr::Var("omega".to_owned()).pow(Exp::from(2)),
        DimExpr::Inverse(Box::new(DimExpr::base("Alpha"))),
        DimExpr::Var("Alpha".to_owned()),
    ]);
    let term = input.to_term(&mut arena).unwrap();
    let expected = DimExpr::Mul(vec![
        DimExpr::Inverse(Box::new(DimExpr::base("Alpha"))),
        DimExpr::base("Zeta").pow(Exp::Rat(Rational::new(-2, 3))),
        DimExpr::Var("Alpha".to_owned()),
        DimExpr::Var("omega".to_owned()).pow(Exp::from(2)),
    ]);
    let saved = arena.snapshot();
    let got = DimExpr::try_from_term(&term, &arena).unwrap();
    assert_eq!(got, expected);
    assert_eq!(saved, arena.snapshot());
    assert_eq!(got.to_term(&mut arena).unwrap(), term);
    assert_eq!(
        DimExpr::try_from_term(&DimTerm::dimensionless(), &arena).unwrap(),
        DimExpr::Dimless
    );
}

#[test]
fn dimensional_rational_boundary_handles_signed_limits_and_wide_fractions() {
    let mut arena = DimArena::new();
    let length = arena.base_id("Length").unwrap();
    for exp in [
        Rational::int(i64::MIN),
        Rational::int(i64::MAX),
        Rational::new(1, i64::MAX),
    ] {
        let term = DimTerm::base(length, ratio(exp.numer, exp.denom));
        let got = DimExpr::try_from_term(&term, &arena).unwrap();
        assert_eq!(got, DimExpr::base("Length").pow(Exp::Rat(exp)));
        assert_eq!(got.to_term(&mut arena).unwrap(), term);
    }
    let outside = BigInt::from(i64::MAX) + 1_i32;
    for exponent in [
        BigRational::from_integer(outside.clone()),
        BigRational::from_integer(-outside.clone() - 1),
        BigRational::new(BigInt::from(1), outside),
    ] {
        let term = DimTerm::base(length, exponent);
        let err = DimExpr::try_from_term(&term, &arena).unwrap_err();
        assert_eq!(err.code(), "parse");
        assert!(err.to_string().contains("rational range"));
        assert_eq!(DimExpr::from_term(&term, &arena), DimExpr::Any);
    }
}

#[test]
fn dimensional_rational_boundary_preserves_alias_ambiguity_and_arena_errors() {
    let mut arena = DimArena::new();
    let length = arena.base_id("Length").unwrap();
    let atom = DimTerm::base_symbolic(length, "n", &mut arena);
    let compound = atom.mul(&atom);
    let display = compound.symbolic_factors(&arena)[0].1.clone();
    let homonym = DimTerm::base_symbolic(length, &display, &mut arena);
    assert_eq!(compound.describe(&arena), homonym.describe(&arena));
    assert!(
        DimExpr::try_from_term(&compound, &arena)
            .unwrap_err()
            .to_string()
            .contains("not a single symbol")
    );
    let literal = DimExpr::base("Length").pow(Exp::Symbol(display));
    assert_eq!(DimExpr::try_from_term(&homonym, &arena).unwrap(), literal);
    assert!(
        DimExpr::try_from_term(&DimTerm::base_symbolic(length, "1/2", &mut arena), &arena).is_err()
    );

    let missing = DimTerm::base(BaseId(u32::MAX), ratio(1, 1));
    let missing_var = DimTerm::var_pow(DimVarId(u32::MAX), ratio(1, 1));
    for term in [missing, missing_var] {
        assert_eq!(
            DimExpr::try_from_term(&term, &arena).unwrap_err().code(),
            "parse"
        );
    }
    // The existing compatibility path can resolve the printed placeholder
    // when re-interning it happens to recreate this exact ID on an arena copy.
    let next = BaseId(u32::try_from(arena.bases().len()).unwrap());
    let next_dim = DimTerm::base(next, ratio(1, 1));
    assert_eq!(
        DimExpr::try_from_term(&next_dim, &arena).unwrap(),
        DimExpr::base("?")
    );
    assert!(arena.base_id("?").is_none());

    let duplicate = Dimension {
        factors: vec![("Length".to_owned(), "1".to_owned()); 2],
        ..Dimension::default()
    };
    let duplicate = DimTerm::from_dimension(&duplicate, &arena).unwrap();
    assert!(DimExpr::try_from_term(&duplicate, &arena).is_err());
}

mod session_value {
    use super::*;
    use rssn::api::Session;
    use rssn::graph::RuleError;
    use rssn::graph::RuleSet;
    use rssn::graph::rule::Installer;
    use rssn::kernel::OpError;
    use rssn::kernel::value::DynValue;
    use rssn::kernel::value::Value;

    fn populate(installer: &mut Installer<'_>) -> Result<(), RuleError> {
        installer.graph().with_dim_arena_mut(|arena| {
            let length = arena.base_id("Length").unwrap();
            for i in 0..2048 {
                arena.define_base(&format!("host_base_{i}"));
                DimExpr::Var(format!("host_var_{i}"))
                    .to_term(arena)
                    .unwrap();
                _ = DimTerm::base_symbolic(length, &format!("host_atom_{i}"), arena);
            }
        });
        Ok(())
    }

    #[test]
    fn dimensional_rational_boundary_keeps_the_session_value_term_and_conflict() {
        let session = Session::with_rules(&[
            rssn::rules::units(),
            RuleSet::new("dim_boundary_test_arena", populate),
        ])
        .unwrap();
        let magnitude = session.parse("quantity(3,m/s^2)").unwrap();
        let actual = magnitude.dimension().exact_factors().unwrap();
        assert_eq!(
            actual,
            vec![
                ("Length".to_owned(), ratio(1, 1)),
                ("Time".to_owned(), ratio(-2, 1))
            ]
        );
        let wrapped = |dim| {
            Value::Dyn(DynValue {
                value: Box::new(Value::Term(magnitude.handle())),
                dim,
            })
        };
        let expected = DimExpr::Mul(vec![
            DimExpr::base("Length"),
            DimExpr::base("Time").pow(Exp::from(-2)),
        ]);
        let reused = session.term_of(&wrapped(expected)).unwrap();
        assert_eq!(reused.id(), magnitude.id());
        let Value::Quantity(quantity) = reused.to_value().unwrap() else {
            panic!("expected a quantity")
        };
        assert_eq!(*quantity.value, Value::Int(BigInt::from(3)));
        let error = session
            .term_of(&wrapped(DimExpr::base("Time")))
            .expect_err("conflicting dimension");
        assert!(
            matches!(error, OpError::Domain(reason) if reason == "dynamic dimension conflicts with term dimension")
        );
    }
}
