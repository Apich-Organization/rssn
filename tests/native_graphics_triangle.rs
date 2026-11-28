//! Physical, bounded prepared triangle operators through the native session ABI.

#![allow(clippy::unwrap_used)]

use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::manifest::SignatureSource;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::NatExpr;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::TensorValue;
use rssn::kernel::value::VariantValue;

const PREPARE: &str = "num_native_graphics_prepared_triangle_new";
const INTERSECT: &str = "num_native_graphics_prepared_triangle_intersect";
const BATCH: &str = "num_native_graphics_prepared_triangle_intersect_parallel";
const POINT: &str = "num_native_graphics_point3_new";
const VECTOR: &str = "num_native_graphics_vector3_new";
const RAY: &str = "num_native_graphics_ray_new";
const RELEASE: &str = "num_native_release";
const POINT_KIND: &str = "rssn/native/discrete/graphics-point3d/v1";
const TRIANGLE_KIND: &str = "rssn/native/discrete/graphics-prepared-triangle/v1";

fn invoke(
    s: &Session,
    name: &str,
    args: &[Value],
) -> Result<Value, OpError> {
    s.call_value(name, args, &CallCtx::new())
}

fn call(
    s: &Session,
    name: &str,
    args: &[Value],
) -> Value {
    invoke(s, name, args).unwrap_or_else(|err| panic!("{name}: {err:?}"))
}

fn length() -> DimExpr {
    DimExpr::Base("Length".into())
}

fn measured(
    value: f64,
    dim: &DimExpr,
) -> Value {
    if *dim == DimExpr::Dimless {
        return Value::Float(value);
    }
    Value::Dyn(DynValue {
        dim: dim.clone(),
        value: Box::new(Value::Float(value)),
    })
}

fn meters(value: f64) -> Value {
    measured(value, &length())
}

fn some(value: Value) -> Value {
    Value::Variant(VariantValue {
        tag: "Some".into(),
        payload: vec![value],
    })
}

fn unbounded() -> Value {
    Value::Variant(VariantValue {
        tag: "None".into(),
        payload: vec![],
    })
}

fn optional(value: Value) -> Option<Value> {
    let Value::Variant(mut value) = value else {
        panic!("expected a native option");
    };
    match value.tag.as_str() {
        | "Some" => {
            assert_eq!(value.payload.len(), 1);
            value.payload.pop()
        },
        | "None" => {
            assert!(value.payload.is_empty());
            None
        },
        | other => panic!("unexpected option tag {other}"),
    }
}

fn field<'a>(
    record: &'a Value,
    name: &str,
) -> &'a Value {
    let Value::Record(record) = record else {
        panic!("expected a record");
    };
    record
        .fields
        .iter()
        .find(|(k, _)| k == name)
        .map(|(_, v)| v)
        .unwrap()
}

fn tensor<'a>(
    value: &'a Value,
    expected_shape: &[usize],
    elem: TypeExpr,
) -> &'a TensorValue {
    let Value::Tensor(tensor) = value else {
        panic!("expected a tensor");
    };
    assert_eq!(tensor.shape, expected_shape);
    assert_eq!(tensor.elem, elem);
    tensor
}

fn numeric(
    value: &Value,
    expected: &DimExpr,
) -> f64 {
    if let Value::Dyn(value) = value {
        assert_eq!(&value.dim, expected);
        return value.value.as_f64().unwrap();
    }
    assert_eq!(expected, &DimExpr::Dimless);
    value.as_f64().unwrap()
}

fn point(
    s: &Session,
    p: [f64; 3],
    dim: &DimExpr,
) -> Value {
    call(s, POINT, &p.map(|v| measured(v, dim)))
}

fn direction(
    s: &Session,
    components: [f64; 3],
    dim: &DimExpr,
) -> Value {
    call(s, VECTOR, &components.map(|v| measured(v, dim)))
}

fn origins(values: Vec<Value>) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![values.len()],
        elem: TypeExpr::Opaque(POINT_KIND.into()),
        data: values,
    })
}

fn usage(s: &Session) -> (f64, f64) {
    let usage = call(s, "num_native_usage", &[]);
    (
        numeric(field(&usage, "objects"), &DimExpr::Dimless),
        numeric(field(&usage, "bytes"), &DimExpr::Dimless),
    )
}

struct Scene {
    triangle: Value,
    triangle_bytes: f64,
    vertices: [Value; 3],
    origin: Value,
    direction: Value,
    ray: Value,
}

fn scene(
    s: &Session,
    dim: &DimExpr,
) -> Scene {
    let vertices = [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0], [0.0, 4.0, 0.0]].map(|p| point(s, p, dim));
    let [a, b, c] = vertices.clone();
    let before = usage(s);
    let triangle = optional(call(s, PREPARE, &[a, b, c, measured(4.0, dim)])).unwrap();
    let after = usage(s);
    let triangle_bytes = after.1 - before.1;
    assert_eq!(after.0, before.0 + 1.0);
    assert!(triangle_bytes >= 512.0); // The bank also charges its entry metadata.
    let origin = point(s, [1.0, 1.0, 4.0], dim);
    let direction = direction(s, [0.0, 0.0, -1.0], &DimExpr::Dimless);
    let ray = call(
        s,
        RAY,
        &[origin.clone(), direction.clone(), measured(16.0, dim)],
    );
    Scene {
        triangle,
        triangle_bytes,
        vertices,
        origin,
        direction,
        ray,
    }
}

#[test]
fn scalar_bounds_restore_physical_values_and_own_the_triangle() {
    let s = Session::new();
    let Scene {
        triangle,
        triangle_bytes,
        vertices,
        ray,
        ..
    } = scene(&s, &length());
    for vertex in vertices {
        assert_eq!(call(&s, RELEASE, &[vertex]), Value::Bool(true));
    }
    let before_hit = usage(&s);
    let hit = optional(call(
        &s,
        INTERSECT,
        &[
            triangle.clone(),
            ray.clone(),
            meters(4.0),
            some(meters(4.0)),
        ],
    ))
    .unwrap();
    assert_eq!(numeric(field(&hit, "t"), &length()), 4.0);
    let p = call(
        &s,
        "num_native_graphics_point3_components",
        &[field(&hit, "point").clone()],
    );
    let n = call(
        &s,
        "num_native_graphics_vector3_components",
        &[field(&hit, "normal").clone()],
    );
    let point = tensor(&p, &[3], TypeExpr::Float.with_dim(length()));
    assert_eq!(
        point
            .data
            .iter()
            .map(|v| numeric(v, &length()))
            .collect::<Vec<_>>(),
        [1.0, 1.0, 0.0]
    );
    let normal = tensor(&n, &[3], TypeExpr::Float.with_dim(DimExpr::Dimless));
    assert_eq!(
        normal
            .data
            .iter()
            .map(|v| numeric(v, &DimExpr::Dimless))
            .collect::<Vec<_>>(),
        [0.0, 0.0, 1.0]
    );
    for name in ["point", "normal"] {
        call(&s, RELEASE, &[field(&hit, name).clone()]);
    }
    assert_eq!(usage(&s), before_hit);
    assert!(
        optional(call(
            &s,
            INTERSECT,
            &[
                triangle.clone(),
                ray.clone(),
                meters(0.0),
                some(meters(3.0))
            ]
        ))
        .is_none()
    );
    let unbounded_hit = optional(call(
        &s,
        INTERSECT,
        &[triangle.clone(), ray.clone(), meters(0.0), unbounded()],
    ))
    .unwrap();
    for name in ["point", "normal"] {
        call(&s, RELEASE, &[field(&unbounded_hit, name).clone()]);
    }
    let before_release = usage(&s);
    assert_eq!(
        call(&s, RELEASE, std::slice::from_ref(&triangle)),
        Value::Bool(true)
    );
    assert_eq!(
        usage(&s),
        (before_release.0 - 1.0, before_release.1 - triangle_bytes)
    );
    assert!(invoke(&s, INTERSECT, &[triangle, ray, meters(0.0), unbounded()]).is_err());
}

#[test]
fn batch_preserves_shapes_units_and_budget_without_retaining_objects() {
    let s = Session::new();
    let scene = scene(&s, &length());
    let miss = point(&s, [5.0, 5.0, 4.0], &length());
    let edge = point(&s, [2.0, 2.0, 8.0], &length());
    let run = |values| {
        call(
            &s,
            BATCH,
            &[
                scene.triangle.clone(),
                origins(values),
                scene.direction.clone(),
                meters(0.0),
                unbounded(),
            ],
        )
    };
    let before = usage(&s);
    let measured_type = TypeExpr::Float.with_dim(length());
    let out = run(vec![scene.origin.clone(), miss, edge.clone()]);
    assert_eq!(
        tensor(field(&out, "hit_mask"), &[3], TypeExpr::Bool).data,
        [Value::Bool(true), Value::Bool(false), Value::Bool(true)]
    );
    assert_eq!(
        tensor(field(&out, "t"), &[3], measured_type.clone())
            .data
            .iter()
            .map(|v| numeric(v, &length()))
            .collect::<Vec<_>>(),
        [4.0, 0.0, 8.0]
    );
    let points = tensor(field(&out, "points"), &[3, 3], measured_type.clone());
    assert_eq!(
        points
            .data
            .iter()
            .map(|v| numeric(v, &length()))
            .collect::<Vec<_>>(),
        [1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 2.0, 2.0, 0.0]
    );
    let normals = tensor(field(&out, "normals"), &[3, 3], TypeExpr::Float);
    assert_eq!(
        normals
            .data
            .iter()
            .map(|v| numeric(v, &DimExpr::Dimless))
            .collect::<Vec<_>>(),
        [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
    );
    let clipped = call(
        &s,
        BATCH,
        &[
            scene.triangle.clone(),
            origins(vec![scene.origin.clone(), edge]),
            scene.direction.clone(),
            meters(0.0),
            some(meters(4.0)),
        ],
    );
    assert_eq!(
        tensor(field(&clipped, "hit_mask"), &[2], TypeExpr::Bool).data,
        [Value::Bool(true), Value::Bool(false)]
    );

    let empty = run(vec![]);
    assert!(
        tensor(field(&empty, "hit_mask"), &[0], TypeExpr::Bool)
            .data
            .is_empty()
    );
    assert!(
        tensor(field(&empty, "t"), &[0], measured_type.clone())
            .data
            .is_empty()
    );
    assert!(
        tensor(field(&empty, "points"), &[0, 3], measured_type.clone())
            .data
            .is_empty()
    );
    assert!(
        tensor(field(&empty, "normals"), &[0, 3], TypeExpr::Float)
            .data
            .is_empty()
    );

    let full = run(vec![scene.origin.clone(); 4096]);
    let mask = tensor(field(&full, "hit_mask"), &[4096], TypeExpr::Bool);
    assert!(mask.data.iter().all(|v| *v == Value::Bool(true)));
    let t = tensor(field(&full, "t"), &[4096], measured_type.clone());
    assert!(t.data.iter().all(|v| numeric(v, &length()) == 4.0));
    let p = tensor(field(&full, "points"), &[4096, 3], measured_type);
    let n = tensor(field(&full, "normals"), &[4096, 3], TypeExpr::Float);
    assert_eq!(
        mask.data.len() + t.data.len() + p.data.len() + n.data.len(),
        8 * 4096
    );
    let too_many = [
        scene.triangle,
        origins(vec![scene.origin; 4097]),
        scene.direction,
        meters(0.0),
        unbounded(),
    ];
    let error = invoke(&s, BATCH, &too_many).unwrap_err();
    assert!(matches!(error, OpError::Domain(message) if message.contains("4096")));
    assert_eq!(usage(&s), before);
}

#[test]
fn validates_common_dimensions_session_scope_and_optional_bound_tags() {
    let s = Session::new();
    let other = Session::new();
    let other_scene = scene(&other, &length());
    let scene = scene(&s, &length());
    let time = DimExpr::Base("Time".into());
    let [a, b, c] = scene.vertices;
    assert!(
        optional(call(
            &s,
            PREPARE,
            &[a.clone(), a.clone(), b.clone(), meters(1.0)]
        ))
        .is_none()
    );
    for invalid in [meters(0.0), meters(-4.0), measured(4.0, &time)] {
        assert!(invoke(&s, PREPARE, &[a.clone(), b.clone(), c.clone(), invalid]).is_err());
    }
    let foreign_point = point(&other, [1.0, 1.0, 0.0], &length());
    assert!(
        invoke(
            &s,
            PREPARE,
            &[a.clone(), b.clone(), foreign_point.clone(), meters(4.0)]
        )
        .is_err()
    );
    let time_point = point(&s, [1.0, 1.0, 0.0], &time);
    assert!(invoke(&s, PREPARE, &[a, b, time_point.clone(), meters(4.0)]).is_err());
    assert!(
        invoke(
            &s,
            INTERSECT,
            &[
                other_scene.triangle,
                scene.ray.clone(),
                meters(0.0),
                unbounded()
            ]
        )
        .is_err()
    );
    assert!(
        invoke(
            &s,
            INTERSECT,
            &[
                scene.triangle.clone(),
                other_scene.ray,
                meters(0.0),
                unbounded()
            ]
        )
        .is_err()
    );
    for (min, max) in [
        (measured(0.0, &time), unbounded()),
        (meters(0.0), some(measured(4.0, &time))),
        (meters(-1.0), unbounded()),
        (meters(5.0), some(meters(4.0))),
        (meters(0.0), some(meters(f64::INFINITY))),
        (
            meters(0.0),
            Value::Variant(VariantValue {
                tag: "Some".into(),
                payload: vec![],
            }),
        ),
    ] {
        assert!(
            invoke(
                &s,
                INTERSECT,
                &[scene.triangle.clone(), scene.ray.clone(), min, max]
            )
            .is_err()
        );
    }
    for invalid in [foreign_point, time_point] {
        assert!(
            invoke(
                &s,
                BATCH,
                &[
                    scene.triangle.clone(),
                    origins(vec![scene.origin.clone(), invalid]),
                    scene.direction.clone(),
                    meters(0.0),
                    unbounded()
                ]
            )
            .is_err()
        );
    }
    for invalid in [
        direction(&s, [0.0, 0.0, -1.0], &length()),
        direction(&s, [0.0, 0.0, -2.0], &DimExpr::Dimless),
    ] {
        assert!(
            invoke(
                &s,
                BATCH,
                &[
                    scene.triangle.clone(),
                    origins(vec![]),
                    invalid,
                    meters(0.0),
                    unbounded()
                ]
            )
            .is_err()
        );
    }
}

#[test]
fn a_shared_non_length_coordinate_dimension_is_preserved() {
    let s = Session::new();
    let dimension = DimExpr::Base("Time".into());
    let scene = scene(&s, &dimension);
    let hit = optional(call(
        &s,
        INTERSECT,
        &[
            scene.triangle,
            scene.ray,
            measured(0.0, &dimension),
            unbounded(),
        ],
    ))
    .unwrap();
    assert_eq!(numeric(field(&hit, "t"), &dimension), 4.0);
}

#[test]
fn manifest_declares_existing_handle_types_and_explicit_optional_bound() {
    let s = Session::new();
    let manifest = s.manifest();
    let opaque = TypeExpr::Opaque(TRIANGLE_KIND.into());
    for (name, source, arity) in [
        (PREPARE, "new", 4),
        (INTERSECT, "intersect", 4),
        (BATCH, "intersect_parallel", 5),
    ] {
        let method = manifest.op(name).unwrap();
        assert_eq!(method.signature_source, SignatureSource::Declared);
        assert_eq!(method.overloads.len(), 1);
        assert_eq!(method.overloads[0].params.len(), arity);
        let expected = format!("crate::kernels::computer_graphics::PreparedTriangle::{source}");
        assert_eq!(method.rust_path.as_deref(), Some(expected.as_str()));
        assert!(!method.doc.is_empty());
    }
    let option_of = |element| {
        TypeExpr::Variant(vec![
            ("None".into(), vec![]),
            ("Some".into(), vec![element]),
        ])
    };
    let prepare = &manifest.op(PREPARE).unwrap().overloads[0];
    assert_eq!(prepare.ret, option_of(opaque.clone()));
    let scalar = &manifest.op(INTERSECT).unwrap().overloads[0];
    let batch = &manifest.op(BATCH).unwrap().overloads[0];
    assert_eq!(scalar.params[0], opaque);
    assert_eq!(scalar.params[3], option_of(scalar.params[2].clone()));
    assert_eq!(batch.params[4], scalar.params[3]);
    assert_eq!(
        batch.params[1],
        TypeExpr::Tensor {
            rank: NatExpr::Lit(1),
            shape: None,
            elem: Box::new(TypeExpr::Opaque(POINT_KIND.into()))
        }
    );
}
