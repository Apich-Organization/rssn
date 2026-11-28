//! Session and manifest coverage for prepared measured rectilinear grids.

use num_traits::ToPrimitive;
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::RepExpr;
use rssn::kernel::signature::ThreadSafety;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::QuantityValue;
use rssn::kernel::value::ScaleValue;
use rssn::kernel::value::TensorValue;
use rssn::kernel::value::UnitValue;

const NEW: &str = "num_native_multilinear_grid_new";
const EVAL: &str = "num_native_multilinear_grid_eval";
const KIND: &str = "rssn/native/interpolation/multilinear-grid/v1";

fn call(
    session: &Session,
    name: &str,
    args: &[Value],
) -> Result<Value, OpError> {
    session.call_value(name, args, &CallCtx::new())
}

fn vector(
    data: Vec<Value>,
    elem: TypeExpr,
) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![data.len()],
        data,
        elem,
    })
}

fn dynamic_vector(data: Vec<Value>) -> Value {
    vector(
        data,
        TypeExpr::Dyn(RepExpr {
            repr: "term".into(),
            dim: None,
        }),
    )
}

fn measured_vector(
    values: &[f64],
    dim: DimExpr,
) -> Value {
    vector(
        values.iter().copied().map(Value::Float).collect(),
        TypeExpr::Float.with_dim(dim),
    )
}

fn measured(
    value: f64,
    dimension: DimExpr,
) -> Value {
    Value::Dyn(DynValue {
        value: Box::new(Value::Float(value)),
        dim: dimension,
    })
}

fn centimeters(value: f64) -> Value {
    Value::Quantity(QuantityValue {
        value: Box::new(Value::Float(value)),
        unit: UnitValue {
            name: "cm".into(),
            display: "cm".into(),
            dimension: DimExpr::base("Length"),
            scale: ScaleValue {
                numer: 1.into(),
                denom: 100.into(),
                pi: 0,
            },
            offset: None,
            affine: false,
            aliases: vec![],
            prefixable: false,
            category: "length".into(),
            since: String::new(),
        },
    })
}

fn usage(session: &Session) -> (usize, usize) {
    let Value::Record(usage) = call(session, "num_native_usage", &[]).unwrap() else {
        panic!("expected a native usage record");
    };
    let get = |field: &str| {
        let (_, Value::Int(value)) = usage.fields.iter().find(|(name, _)| name == field).unwrap()
        else {
            panic!("expected an integer usage field");
        };
        value.to_usize().unwrap()
    };
    (get("objects"), get("bytes"))
}

#[test]
fn manifest_exposes_preparation_and_repeated_lookup() {
    let session = Session::new();
    let manifest = session.manifest();
    let create = manifest.op(NEW).unwrap();
    let evaluate = manifest.op(EVAL).unwrap();
    assert_eq!(create.category, "numerical_vector");
    for op in [create, evaluate] {
        assert_eq!(op.overloads.len(), 1);
        assert_eq!(op.overloads[0].params.len(), 2);
        assert_eq!(op.thread_safety, ThreadSafety::SessionBound);
        assert!(op.backend_support.interpreter);
        assert!(!op.doc.is_empty());
    }
    assert_eq!(create.overloads[0].ret, TypeExpr::Opaque(KIND.into()));
    assert!(matches!(&create.overloads[0].params[0], TypeExpr::Dyn(_)));
    assert!(matches!(
        &create.overloads[0].params[1],
        TypeExpr::Tensor { .. }
    ));
    assert_eq!(
        evaluate.overloads[0].params[0],
        TypeExpr::Opaque(KIND.into())
    );
    assert!(matches!(&evaluate.overloads[0].ret, TypeExpr::Dyn(_)));
}

#[test]
fn irregular_layout_and_separately_measured_axes_match_a_multiaffine_reference() {
    let session = Session::new();
    let input_x = vector(
        vec![centimeters(0.0), centimeters(200.0)],
        TypeExpr::Float.with_dim(DimExpr::base("Length")),
    );
    let input_t = measured_vector(&[0.0, 1.0, 4.0], DimExpr::base("Time"));
    let axes = dynamic_vector(vec![input_x, input_t]);
    let reference = |x, t| 5.0 + 2.0 * x - t + x * t;
    let values: Vec<_> = [0.0, 1.0, 4.0]
        .into_iter()
        .flat_map(|t| [0.0, 2.0].into_iter().map(move |x| reference(x, t)))
        .collect();
    let ordinate = DimExpr::base("Mass");
    let grid = call(
        &session,
        NEW,
        &[axes, measured_vector(&values, ordinate.clone())],
    )
    .unwrap();
    let before = usage(&session);
    assert_eq!(before.0, 1);
    let foreign = Session::new();
    for (x, t) in [(1.0, 0.5), (1.5, 3.0), (2.0, 4.0), (0.0, 0.0)] {
        let query = dynamic_vector(vec![
            centimeters(x * 100.0),
            measured(t, DimExpr::base("Time")),
        ]);
        let result = call(&session, EVAL, &[grid.clone(), query.clone()]).unwrap();
        let Value::Dyn(result) = result else {
            panic!("expected measured interpolation result");
        };
        assert_eq!(result.dim, ordinate);
        assert!((result.value.as_f64().unwrap() - reference(x, t)).abs() < 1e-12);
        assert!(call(&foreign, EVAL, &[grid.clone(), query]).is_err());
    }
    assert_eq!(usage(&session), before);
    for invalid in [
        dynamic_vector(vec![
            centimeters(100.0),
            measured(0.5, DimExpr::base("Length")),
        ]),
        dynamic_vector(vec![centimeters(100.0)]),
        dynamic_vector(vec![
            centimeters(300.0),
            measured(0.5, DimExpr::base("Time")),
        ]),
    ] {
        assert!(call(&session, EVAL, &[grid.clone(), invalid]).is_err());
    }
    assert_eq!(usage(&session), before);
    call(&session, "num_native_release", std::slice::from_ref(&grid)).unwrap();
    assert_eq!(usage(&session), (0, 0));
    let query = dynamic_vector(vec![
        centimeters(100.0),
        measured(0.5, DimExpr::base("Time")),
    ]);
    assert!(call(&session, EVAL, &[grid, query]).is_err());
}

#[test]
fn constructor_errors_leave_no_object_and_grid_charge_tracks_allocated_data() {
    let session = Session::new();
    let axes = |x: &[f64], y: &[f64]| {
        dynamic_vector(vec![
            measured_vector(x, DimExpr::base("Length")),
            measured_vector(y, DimExpr::base("Length")),
        ])
    };
    for (source, values) in [
        (axes(&[0.0, 2.0], &[0.0, 1.0, 4.0]), vec![1.0; 5]),
        (axes(&[0.0, 0.0], &[0.0, 1.0, 4.0]), vec![1.0; 6]),
        (axes(&[0.0, f64::INFINITY], &[0.0, 1.0, 4.0]), vec![1.0; 6]),
        (axes(&[0.0, 2.0], &[0.0, 1.0, 4.0]), vec![f64::NAN; 6]),
    ] {
        assert!(
            call(
                &session,
                NEW,
                &[source, measured_vector(&values, DimExpr::Dimless)]
            )
            .is_err()
        );
        assert_eq!(usage(&session), (0, 0));
    }
    let measure_object = |xs, ys, data: &[f64]| {
        let handle = call(
            &session,
            NEW,
            &[axes(xs, ys), measured_vector(data, DimExpr::Dimless)],
        )
        .unwrap();
        let (objects, bytes) = usage(&session);
        assert_eq!(objects, 1);
        assert_eq!(
            call(
                &session,
                EVAL,
                &[
                    handle.clone(),
                    measured_vector(&[1.0, 0.5], DimExpr::base("Length"))
                ],
            )
            .unwrap(),
            Value::Float(1.0)
        );
        call(&session, "num_native_release", &[handle]).unwrap();
        assert_eq!(usage(&session), (0, 0));
        bytes
    };
    let small = measure_object(&[0.0, 2.0], &[0.0, 1.0, 4.0], &[1.0; 6]);
    let large = measure_object(&[0.0, 1.0, 2.0], &[0.0, 1.0, 3.0, 4.0], &[1.0; 12]);
    assert!(large > small);
    assert!(large < small + 1024);
}
