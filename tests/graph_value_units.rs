//! Structural unit evidence across the host value / symbolic graph boundary.
use num_bigint::BigInt;
use num_rational::BigRational;
use rssn::api::Session;
use rssn::dim::Scale;
use rssn::dim::UnitDef;
use rssn::kernel::CallCtx;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::Exp;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::QuantityValue;
use rssn::kernel::value::TensorValue;

fn quantity(
    value: Value,
    name: &str,
    dim: DimExpr,
    scale: i64,
) -> Value {
    Value::Quantity(QuantityValue {
        value: Box::new(value),
        unit: UnitDef::linear(
            name,
            dim,
            Scale {
                rat: BigRational::from_integer(scale.into()),
                pi: 0,
            },
        )
        .to_value(),
    })
}
fn dynamic(
    value: Value,
    dim: DimExpr,
) -> Value {
    Value::Dyn(DynValue {
        value: Box::new(value),
        dim,
    })
}
fn call(
    session: &Session,
    name: &str,
    args: &[Value],
) -> Value {
    session
        .call_value(name, args, &CallCtx::new())
        .unwrap_or_else(|error| panic!("{name}: {error}"))
}
fn measured(value: &Value) -> (&Value, &rssn::kernel::value::UnitValue) {
    let Value::Quantity(q) = value else {
        panic!("expected quantity, got {value:?}")
    };
    (&q.value, &q.unit)
}

#[test]
fn conflicting_named_definitions_keep_physical_metadata_and_hash_cons() {
    let session = Session::new();
    let value = quantity(Value::int(3), "m", DimExpr::base("Time"), 7);
    let term = session.term_of(&value).unwrap();
    assert_eq!(term.id(), session.term_of(&value.clone()).unwrap().id());
    let report = term.dimension_report();
    assert!(report.errors.is_empty(), "{:?}", report.errors);
    assert_eq!(
        report.dimension.factors,
        vec![("Time".to_owned(), "1".to_owned())]
    );
    assert_eq!(session.value_of(term).unwrap(), value);
    let metre = call(
        &session,
        "quantity",
        &[
            Value::int(3),
            Value::Unit(rssn::dim::units::lookup("m").unwrap().to_value()),
        ],
    );
    assert_ne!(session.term_of(&metre).unwrap().id(), term.id());
}

#[test]
fn unit_arguments_preserve_the_supplied_definition_instead_of_looking_up_its_name() {
    let session = Session::new();
    let unit = UnitDef::linear("m", DimExpr::base("Time"), Scale::one())
        .display("clock")
        .to_value();
    let value = call(
        &session,
        "quantity",
        &[Value::int(3), Value::Unit(unit.clone())],
    );
    assert_eq!(measured(&value).1, &unit);
}

#[test]
fn anonymous_scaled_unit_converts_exactly_without_registration() {
    let session = Session::new();
    let value = quantity(Value::int(3), "private-length", DimExpr::base("Length"), 2);
    let unit = rssn::dim::units::lookup("m").unwrap().to_value();
    let result = call(&session, "convert", &[value, Value::Unit(unit)]);
    assert_eq!(
        measured(&result).0.as_rational(),
        Some(BigRational::from_integer(6.into()))
    );
    assert_eq!(measured(&result).1.name, "m");
}

#[test]
fn custom_base_units_survive_coherent_simplification_and_dimension_query() {
    let session = Session::new();
    let value = quantity(Value::int(4), "private-credit", DimExpr::base("Credit"), 3);
    let coherent = call(&session, "simplify_units", &[value]);
    assert_eq!(
        measured(&coherent).0.as_rational(),
        Some(BigRational::from_integer(12.into()))
    );
    assert_eq!(measured(&coherent).1.dimension, DimExpr::base("Credit"));
    assert_eq!(
        measured(&coherent).1.scale.rational(),
        BigRational::from_integer(1.into())
    );
    let dimension = call(&session, "dimension", &[coherent]);
    assert_eq!(dimension, Value::Symbol("Credit".to_owned()));
}

#[test]
fn identical_unit_names_with_distinct_scales_remain_distinct_product_factors() {
    let session = Session::new();
    let left = quantity(Value::int(2), "local", DimExpr::base("Credit"), 3);
    let right = quantity(Value::int(5), "local", DimExpr::base("Time"), 7);
    let product = call(&session, "mul", &[left, right]);
    let result = call(&session, "simplify_units", &[product]);
    assert_eq!(
        measured(&result).0.as_rational(),
        Some(BigRational::from_integer(210.into()))
    );
    let term = session.term_of(&result).unwrap();
    assert_eq!(
        term.dimension().factors,
        vec![
            ("Credit".to_owned(), "1".to_owned()),
            ("Time".to_owned(), "1".to_owned())
        ]
    );
}

#[test]
fn affine_by_value_units_are_coherent_in_their_own_dimension() {
    let session = Session::new();
    let unit = UnitDef::linear(
        "private-origin",
        DimExpr::base("Length"),
        Scale {
            rat: BigRational::from_integer(2.into()),
            pi: 0,
        },
    )
    .affine(BigRational::from_integer(10.into()))
    .to_value();
    let value = Value::Quantity(QuantityValue {
        value: Box::new(Value::int(3)),
        unit,
    });
    let metre = rssn::dim::units::lookup("m").unwrap().to_value();
    let result = call(&session, "convert", &[value, Value::Unit(metre)]);
    assert_eq!(
        measured(&result).0.as_rational(),
        Some(BigRational::from_integer(16.into()))
    );
}

#[test]
fn dynamic_dimension_is_materialized_and_nested_metadata_cannot_conflict() {
    let session = Session::new();
    let dim = DimExpr::base("Credit");
    let term = session
        .term_of(&dynamic(Value::int(5), dim.clone()))
        .unwrap();
    assert_eq!(
        term.dimension().factors,
        vec![("Credit".to_owned(), "1".to_owned())]
    );
    let back = session.value_of(term).unwrap();
    assert_eq!(
        measured(&back).0.as_rational(),
        Some(BigRational::from_integer(5.into()))
    );
    assert_eq!(measured(&back).1.dimension, dim);
    let q = quantity(Value::int(2), "private-credit", DimExpr::base("Credit"), 3);
    assert!(
        session
            .term_of(&dynamic(q.clone(), DimExpr::base("Time")))
            .is_err()
    );
    let same = session
        .term_of(&dynamic(q.clone(), DimExpr::base("Credit")))
        .unwrap();
    assert_eq!(session.value_of(same).unwrap(), q);
    assert!(
        session
            .term_of(&dynamic(
                dynamic(Value::int(1), DimExpr::base("Time")),
                DimExpr::base("Credit")
            ))
            .is_err()
    );
}

#[test]
fn graph_terms_inside_quantities_keep_the_real_session_handle() {
    let session = Session::new();
    let expression = session.parse("x + 1").unwrap();
    let value = quantity(
        Value::Term(session.handle(expression)),
        "private-credit",
        DimExpr::base("Credit"),
        1,
    );
    let term = session.term_of(&value).unwrap();
    let back = session.value_of(term).unwrap();
    assert_eq!(session.term_of(&back).unwrap().id(), term.id());
}

#[test]
fn unit_dimension_expressions_support_exact_powers_of_products() {
    let session = Session::new();
    let dimension = DimExpr::Mul(vec![DimExpr::base("Length"), DimExpr::base("Time")])
        .pow(Exp::Rat(rssn::kernel::signature::Rational::new(1, 2)));
    let value = quantity(Value::int(1), "root-product", dimension, 1);
    let term = session.term_of(&value).unwrap();
    assert_eq!(
        term.dimension().factors,
        vec![
            ("Length".to_owned(), "1/2".to_owned()),
            ("Time".to_owned(), "1/2".to_owned())
        ]
    );
}

#[test]
fn declared_tensor_element_dimensions_survive_graph_roundtrip() {
    let session = Session::new();
    let value = Value::Tensor(TensorValue {
        shape: vec![2, 1],
        data: vec![Value::int(3), Value::int(4)],
        elem: TypeExpr::Dim(Box::new(TypeExpr::Int), DimExpr::base("Credit")),
    });
    let term = session.term_of(&value).unwrap();
    let Value::Tensor(back) = session.value_of(term).unwrap() else {
        panic!("numeric quantity list must remain a tensor")
    };
    assert_eq!(back.shape, vec![2, 1]);
    assert_eq!(
        back.elem,
        TypeExpr::Dim(Box::new(TypeExpr::Int), DimExpr::base("Credit"))
    );
    assert_eq!(
        measured(&back.data[1]).0.as_rational(),
        Some(BigRational::from_integer(4.into()))
    );
}

#[test]
fn malformed_metadata_fails_without_panicking_or_erasing_the_schema() {
    let session = Session::new();
    let mut value = quantity(Value::int(1), "malformed", DimExpr::base("Length"), 1);
    let Value::Quantity(q) = &mut value else {
        unreachable!()
    };
    q.unit.scale.denom = BigInt::from(0);
    assert!(session.term_of(&value).is_err());
    assert!(
        session
            .term_of(&dynamic(Value::int(2), DimExpr::Var("D".to_owned())))
            .is_err()
    );
    let empty = Value::Tensor(TensorValue {
        shape: vec![0, 2],
        data: vec![],
        elem: TypeExpr::Real,
    });
    assert!(session.term_of(&empty).is_err());
}

#[test]
fn symbolic_core_overload_is_pinned_for_dimensioned_values() {
    let session = Session::new();
    let args = [
        quantity(Value::int(2), "private-credit", DimExpr::base("Credit"), 1),
        quantity(Value::int(3), "private-credit", DimExpr::base("Credit"), 1),
    ];
    let key = session.resolve_overload("add", &args).unwrap();
    let by_name = call(&session, "add", &args);
    let by_key = session.call_key(&key, &args, &CallCtx::new()).unwrap();
    assert_eq!(by_name, by_key);
    let sum = call(&session, "unify_expression", &[by_key]);
    assert_eq!(
        measured(&sum).0.as_rational(),
        Some(BigRational::from_integer(5.into()))
    );
    assert_eq!(measured(&sum).1.dimension, DimExpr::base("Credit"));
}

#[test]
fn symbolic_fallback_keeps_strict_dimension_checks_authoritative() {
    let session = Session::new();
    let args = [
        quantity(Value::int(2), "private-credit", DimExpr::base("Credit"), 1),
        quantity(Value::int(3), "private-time", DimExpr::base("Time"), 1),
    ];
    let context = CallCtx::new().dimensions(rssn::graph::DimensionMode::Strict);
    let result = session.call_value("add", &args, &context);
    assert!(
        matches!(result, Err(rssn::kernel::OpError::Dimension(_))),
        "{result:?}"
    );
}

#[test]
fn heterogeneous_unit_lists_do_not_invent_a_uniform_tensor_schema() {
    let session = Session::new();
    let args = [
        quantity(Value::int(2), "private-credit", DimExpr::base("Credit"), 1),
        quantity(Value::int(3), "private-time", DimExpr::base("Time"), 1),
    ];
    assert!(matches!(call(&session, "list", &args), Value::Term(_)));
}
