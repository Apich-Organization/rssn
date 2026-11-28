//! Native grid transfers preserve the input reference and session boundaries.

use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::TensorValue;

const LEVEL_SIZE: &str = "num_linear_geometric_level_size";
const RESTRICT: &str = "num_linear_geometric_restrict";
const PROLONG: &str = "num_linear_geometric_prolong";

fn call(
    session: &Session,
    name: &str,
    args: &[Value],
) -> Result<Value, OpError> {
    session.call_value(name, args, &CallCtx::new())
}

fn must_call(
    session: &Session,
    name: &str,
    args: &[Value],
) -> Value {
    call(session, name, args).unwrap_or_else(|error| panic!("{name}: {error:?}"))
}

fn ints(values: &[usize]) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![values.len()],
        elem: TypeExpr::Int,
        data: values.iter().map(|&value| Value::int(value)).collect(),
    })
}

fn floats(values: &[f64]) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![values.len()],
        elem: TypeExpr::Float,
        data: values.iter().copied().map(Value::Float).collect(),
    })
}

fn measured(
    value: f64,
    unit: &str,
) -> Value {
    Value::Dyn(DynValue {
        value: Box::new(Value::Float(value)),
        dim: DimExpr::base(unit),
    })
}

fn array(
    session: &Session,
    shape: &[usize],
    values: &[f64],
) -> Value {
    must_call(
        session,
        "num_linear_array_new",
        &[ints(shape), floats(values), measured(3.25, "Temperature")],
    )
}

fn hierarchy(
    session: &Session,
    dim: usize,
    side: usize,
    operator_scale: f64,
) -> Value {
    let response = must_call(
        session,
        "num_linear_geometric_new",
        &[
            Value::int(dim),
            Value::int(side),
            Value::Float(0.0),
            measured(operator_scale, "Length"),
        ],
    );
    let Value::Variant(variant) = response else {
        panic!("expected a hierarchy variant");
    };
    assert_eq!(variant.tag, "Some");
    let [handle] = variant.payload.as_slice() else {
        panic!("expected one hierarchy handle");
    };
    handle.clone()
}

fn field<'a>(
    v: &'a Value,
    key: &str,
) -> &'a Value {
    let Value::Record(record) = v else {
        panic!("expected packed array record");
    };
    &record.fields.iter().find(|(k, _)| k == key).unwrap().1
}

fn assert_array(
    session: &Session,
    value: &Value,
    expected: &[f64],
) {
    let packed = must_call(
        session,
        "num_linear_array_data",
        std::slice::from_ref(value),
    );
    assert_eq!(field(&packed, "shape"), &ints(&[expected.len()]));
    assert_eq!(field(&packed, "data"), &floats(expected));
    let reference = field(&packed, "reference");
    let measured_type = TypeExpr::Float.with_dim(DimExpr::base("Temperature"));
    assert!(
        rssn::kernel::overload::matches_type(&measured_type, &reference.type_expr()).accepted()
    );
    let Value::Dyn(reference) = reference else {
        panic!("expected a dimensioned, canonical Float reference");
    };
    assert_eq!(*reference.value, Value::Float(3.25));
}

#[test]
fn adjacent_level_transfers_preserve_values_and_the_input_reference() {
    let session = Session::new();
    let fine = array(&session, &[7], &[1.0, 4.0, 7.0, 10.0, 13.0, 16.0, 19.0]);

    for operator_scale in [2.0, 2.0_f64.powi(-700), 2.0_f64.powi(700)] {
        let mg = hierarchy(&session, 1, 7, operator_scale);
        for (level, size) in [(0usize, 7usize), (1, 3), (2, 1)] {
            assert_eq!(
                must_call(&session, LEVEL_SIZE, &[mg.clone(), Value::int(level)]),
                Value::int(size)
            );
        }
        let coarse = must_call(
            &session,
            RESTRICT,
            &[mg.clone(), Value::int(0usize), fine.clone()],
        );
        assert_array(&session, &coarse, &[4.0, 10.0, 16.0]);
        let prolong = must_call(
            &session,
            PROLONG,
            &[mg.clone(), Value::int(0usize), coarse.clone()],
        );
        assert_array(&session, &prolong, &[2.0, 4.0, 7.0, 10.0, 13.0, 16.0, 8.0]);
        let one = must_call(
            &session,
            RESTRICT,
            &[mg.clone(), Value::int(1usize), coarse],
        );
        assert_array(&session, &one, &[10.0]);
        let three = must_call(&session, PROLONG, &[mg, Value::int(1usize), one]);
        assert_array(&session, &three, &[5.0, 10.0, 5.0]);
    }

    let min_subnormal = f64::from_bits(1);
    let mg = hierarchy(&session, 1, 3, 2.0);
    let tiny = array(&session, &[3], &[min_subnormal; 3]);
    let coarse = must_call(&session, RESTRICT, &[mg, Value::int(0usize), tiny]);
    assert_array(&session, &coarse, &[min_subnormal]);
}

#[test]
fn native_transfer_rejects_missing_levels_wrong_shapes_and_foreign_handles() {
    let session = Session::new();
    let foreign = Session::new();
    let mg = hierarchy(&session, 1, 7, 2.0);
    let wrong_rank = array(&session, &[1, 7], &[0.0; 7]);
    let short = array(&session, &[6], &[0.0; 6]);
    let one = array(&session, &[1], &[1.0]);
    let fine = array(&session, &[7], &[0.0; 7]);
    let coarse = array(&session, &[3], &[0.0; 3]);
    let foreign_mg = hierarchy(&foreign, 1, 7, 2.0);
    let foreign_fine = array(&foreign, &[7], &[0.0; 7]);
    let foreign_coarse = array(&foreign, &[3], &[0.0; 3]);

    assert!(call(&session, LEVEL_SIZE, &[mg.clone(), Value::int(usize::MAX)]).is_err());
    assert!(
        call(
            &session,
            LEVEL_SIZE,
            &[foreign_mg.clone(), Value::int(0usize)]
        )
        .is_err()
    );
    assert!(call(&foreign, LEVEL_SIZE, &[mg.clone(), Value::int(0usize)]).is_err());
    for (op, level, value) in [
        (RESTRICT, 2, one.clone()),
        (PROLONG, 2, one.clone()),
        (RESTRICT, 0, wrong_rank),
        (RESTRICT, 0, short),
        (PROLONG, 0, one),
        (RESTRICT, 0, foreign_fine),
        (PROLONG, 0, foreign_coarse),
    ] {
        assert!(call(&session, op, &[mg.clone(), Value::int(level), value]).is_err());
    }
    assert!(
        call(
            &session,
            RESTRICT,
            &[foreign_mg, Value::int(0usize), fine.clone()]
        )
        .is_err()
    );
    assert!(call(&foreign, PROLONG, &[mg.clone(), Value::int(0usize), coarse]).is_err());

    // A valid call can reuse the same original hierarchy after rejected calls.
    let coarsened = must_call(&session, RESTRICT, &[mg, Value::int(0usize), fine]);
    assert_array(&session, &coarsened, &[0.0; 3]);
}

#[test]
fn largest_three_dimensional_transfer_reuses_the_reserved_hierarchy() {
    let session = Session::new();
    let mg = hierarchy(&session, 3, 31, 11.0);
    let size = 31usize.pow(3);
    let source = array(&session, &[size], &vec![1.0; size]);
    let expected_coarse = vec![1.0; 15usize.pow(3)];
    let expected_fine: Vec<_> = (0..size)
        .map(|i| {
            [i % 31, (i / 31) % 31, i / (31 * 31)]
                .iter()
                .map(|&axis| {
                    if axis == 0 || axis == 30 {
                        0.5
                    } else {
                        1.0
                    }
                })
                .product()
        })
        .collect();
    for _ in 0..2 {
        let coarse = must_call(
            &session,
            RESTRICT,
            &[mg.clone(), Value::int(0usize), source.clone()],
        );
        assert_array(&session, &coarse, &expected_coarse);
        let fine = must_call(&session, PROLONG, &[mg.clone(), Value::int(0usize), coarse]);
        assert_array(&session, &fine, &expected_fine);
    }
}
