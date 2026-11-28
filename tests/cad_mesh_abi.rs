//! Checked measured triangle meshes over the native value ABI.
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::DimRuleInfo;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::QuantityValue;
use rssn::kernel::value::ScaleValue;
use rssn::kernel::value::TensorValue;
use rssn::kernel::value::UnitValue;

const FACES: [[usize; 3]; 4] = [[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]];

fn length(value: f64) -> Value {
    Value::Dyn(DynValue {
        value: Box::new(Value::Float(value)),
        dim: DimExpr::base("Length"),
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

fn vertices(
    points: &[[f64; 3]],
    measure: impl Fn(f64) -> Value,
) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![points.len(), 3],
        elem: TypeExpr::Float.with_dim(DimExpr::base("Length")),
        data: points.iter().flatten().copied().map(measure).collect(),
    })
}

fn indices(faces: &[[usize; 3]]) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![faces.len(), 3],
        elem: TypeExpr::Int,
        data: faces.iter().flatten().copied().map(Value::int).collect(),
    })
}

fn points(
    origin: f64,
    edge: f64,
) -> [[f64; 3]; 4] {
    let end = origin + edge;
    [
        [origin; 3],
        [end, origin, origin],
        [origin, end, origin],
        [origin, origin, end],
    ]
}

fn call(
    session: &Session,
    name: &str,
    args: &[Value],
) -> Result<Value, OpError> {
    session.call_value(&format!("num_cad_{name}"), args, &CallCtx::new())
}

fn one(
    session: &Session,
    name: &str,
    mesh: &Value,
) -> Value {
    call(session, name, std::slice::from_ref(mesh)).unwrap()
}

fn field(
    record: &Value,
    name: &str,
) -> Value {
    let Value::Record(record) = record else {
        panic!("expected record")
    };
    record
        .fields
        .iter()
        .find(|(key, _)| key == name)
        .unwrap()
        .1
        .clone()
}

fn volume(
    session: &Session,
    mesh: &Value,
    expected: f64,
) {
    let Value::Dyn(result) = one(session, "mesh_signed_volume", mesh) else {
        panic!("expected dimensional volume")
    };
    assert_eq!(result.dim, DimExpr::base("Length").pow(3));
    let measured = result.value.as_f64().unwrap();
    assert!(
        (measured / expected - 1.0).abs() < 1e-12,
        "{measured} vs {expected}"
    );
}

#[test]
fn measured_triangle_mesh_reuses_topology_and_returns_volume_units() {
    let session = Session::new();
    let mesh = call(
        &session,
        "mesh_new",
        &[vertices(&points(0.0, 100.0), centimeters), indices(&FACES)],
    )
    .unwrap();
    for _ in 0..2 {
        assert_eq!(one(&session, "mesh_vertex_count", &mesh), Value::int(4));
        assert_eq!(one(&session, "mesh_triangle_count", &mesh), Value::int(4));
        assert_eq!(one(&session, "mesh_is_closed", &mesh), Value::Bool(true));
        assert_eq!(one(&session, "mesh_is_manifold", &mesh), Value::Bool(true));
        assert_eq!(
            one(&session, "mesh_euler_characteristic", &mesh),
            Value::int(2)
        );
    }
    volume(&session, &mesh, 1.0 / 6.0);
    let arrays = one(&session, "mesh_arrays", &mesh);
    assert_eq!(
        field(&arrays, "vertices"),
        vertices(&points(0.0, 1.0), length)
    );
    assert_eq!(field(&arrays, "triangles"), indices(&FACES));
    let translated = call(
        &session,
        "mesh_new",
        &[vertices(&points(1e12, 1.0), length), indices(&FACES)],
    )
    .unwrap();
    volume(&session, &translated, 1.0 / 6.0);
    let inward = FACES.map(|[a, b, c]| [a, c, b]);
    let inward = call(
        &session,
        "mesh_new",
        &[vertices(&points(0.0, 1.0), length), indices(&inward)],
    )
    .unwrap();
    volume(&session, &inward, -1.0 / 6.0);
}

#[test]
fn repair_fills_a_hole_and_weld_accepts_small_geometry_and_tolerances() {
    let session = Session::new();
    let open = call(
        &session,
        "mesh_new",
        &[vertices(&points(0.0, 1.0), length), indices(&FACES[..3])],
    )
    .unwrap();
    assert_eq!(one(&session, "mesh_is_closed", &open), Value::Bool(false));
    assert!(call(&session, "mesh_signed_volume", std::slice::from_ref(&open)).is_err());
    let repaired = call(&session, "mesh_repair", &[open, length(0.0)]).unwrap();
    assert_eq!(field(&repaired, "filled_holes"), Value::int(1));
    volume(&session, &field(&repaired, "mesh"), 1.0 / 6.0);
    let mut duplicated = points(0.0, 1e-10).to_vec();
    duplicated.push(duplicated[0]);
    let mut faces = FACES;
    faces[0][0] = 4;
    let nanomesh = call(
        &session,
        "mesh_new",
        &[vertices(&duplicated, length), indices(&faces)],
    )
    .unwrap();
    let welded = call(&session, "mesh_weld", &[nanomesh, length(1e-200)]).unwrap();
    assert_eq!(field(&welded, "merged_vertices"), Value::int(1));
    let mesh = field(&welded, "mesh");
    assert_eq!(one(&session, "mesh_triangle_count", &mesh), Value::int(4));
    volume(&session, &mesh, 1e-30 / 6.0);
}

#[test]
fn stl_round_trips_byte_and_text_payloads_with_explicit_units() {
    let session = Session::new();
    let mesh = call(
        &session,
        "mesh_new",
        &[vertices(&points(0.0, 1.0), length), indices(&FACES)],
    )
    .unwrap();
    let unit = centimeters(1.0);
    let binary = call(&session, "stl_write_binary", &[mesh.clone(), unit.clone()]).unwrap();
    let Value::Tensor(bytes) = &binary else {
        panic!("expected byte tensor")
    };
    assert_eq!(bytes.shape, [284]);
    let ascii = call(
        &session,
        "stl_write_ascii",
        &[
            mesh.clone(),
            unit.clone(),
            Value::String("triangle model".into()),
        ],
    )
    .unwrap();
    let Value::String(text) = &ascii else {
        panic!("expected text")
    };
    assert!(text.contains("100.000000"));
    for (op, data) in [("stl_read", binary), ("stl_read_text", ascii)] {
        let read = call(&session, op, &[data, unit.clone()]).unwrap();
        assert_eq!(one(&session, "mesh_vertex_count", &read), Value::int(12));
        let repaired = call(&session, "mesh_weld", &[read, length(0.0)]).unwrap();
        assert_eq!(field(&repaired, "merged_vertices"), Value::int(8));
        volume(&session, &field(&repaired, "mesh"), 1.0 / 6.0);
    }
    assert!(
        call(
            &session,
            "stl_write_binary",
            &[mesh.clone(), length(1e-100)]
        )
        .is_err()
    );
    assert!(
        call(
            &session,
            "stl_write_ascii",
            &[mesh, unit, Value::String("bad\nlabel".into())]
        )
        .is_err()
    );
}

#[test]
fn native_mesh_rejects_wrong_units_bad_indices_unwound_volume_and_stale_scope() {
    let session = Session::new();
    let coordinates = vertices(&points(0.0, 1.0), length);
    let mut bad_indices = FACES;
    bad_indices[0][0] = 4;
    assert!(
        call(
            &session,
            "mesh_new",
            &[coordinates.clone(), indices(&bad_indices)]
        )
        .is_err()
    );
    let wrong = vertices(&points(0.0, 1.0), |value| {
        Value::Dyn(DynValue {
            value: Box::new(Value::Float(value)),
            dim: DimExpr::base("Time"),
        })
    });
    assert!(call(&session, "mesh_new", &[wrong, indices(&FACES)]).is_err());
    let mut winding = FACES;
    winding[0].swap(1, 2);
    let mismatched = call(
        &session,
        "mesh_new",
        &[coordinates.clone(), indices(&winding)],
    )
    .unwrap();
    assert_eq!(
        one(&session, "mesh_is_closed", &mismatched),
        Value::Bool(true)
    );
    assert!(call(&session, "mesh_signed_volume", &[mismatched]).is_err());
    let mesh = call(&session, "mesh_new", &[coordinates, indices(&FACES)]).unwrap();
    assert!(
        call(
            &Session::new(),
            "mesh_vertex_count",
            std::slice::from_ref(&mesh)
        )
        .is_err()
    );
    assert!(call(&session, "mesh_weld", &[mesh.clone(), Value::Float(0.0)]).is_err());
    assert_eq!(
        session
            .call_value(
                "num_native_release",
                std::slice::from_ref(&mesh),
                &CallCtx::new()
            )
            .unwrap(),
        Value::Bool(true)
    );
    assert!(call(&session, "mesh_vertex_count", &[mesh]).is_err());
}

#[test]
fn cad_manifest_declares_source_result_dimensions_and_owned_handles() {
    let manifest = Session::new().manifest();
    let ops = manifest
        .ops
        .iter()
        .filter(|entry| entry.name.starts_with("num_cad_"))
        .collect::<Vec<_>>();
    assert_eq!(
        ops.iter().filter(|entry| entry.rust_path.is_some()).count(),
        14
    );
    let volume = ops
        .iter()
        .find(|op| op.name == "num_cad_mesh_signed_volume")
        .unwrap();
    assert_eq!(
        volume.overloads[0].ret,
        TypeExpr::Float.with_dim(DimExpr::base("Length").pow(3))
    );
    assert!(matches!(volume.dim_rule, Some(DimRuleInfo::Custom(_))));
    for op in ops {
        assert_eq!(
            op.signature_source,
            rssn::kernel::manifest::SignatureSource::Declared
        );
        assert_eq!(op.overloads.len(), 1);
    }
}

#[test]
fn measured_planar_step_enforces_supported_geometry_and_budget() {
    const SOURCE: &str = include_str!("fixtures/cad_planar_tetra.step");
    let session = Session::new();
    let models = call(&session, "step_read_text", &[Value::String(SOURCE.into())]).unwrap();
    let Value::Tensor(models) = models else {
        panic!("expected STEP solid handles")
    };
    assert_eq!(models.shape, [1]);
    assert_eq!(
        one(&session, "mesh_triangle_count", &models.data[0]),
        Value::int(4)
    );
    volume(&session, &models.data[0], 1.0 / 6.0);
    let curved = SOURCE.replace("PLANE('',#121)", "CYLINDRICAL_SURFACE('',#121,1.)");
    assert_ne!(curved, SOURCE);
    assert!(call(&session, "step_read_text", &[Value::String(curved)]).is_err());
    let oversized = SOURCE.repeat(262_144 / SOURCE.len() + 1);
    assert!(call(&session, "step_read_text", &[Value::String(oversized)]).is_err());
}

#[test]
fn convex_mesh_boolean_returns_owned_measured_results() {
    let session = Session::new();
    let coords = vertices(&points(0.0, 1.0), length);
    let left = call(&session, "mesh_new", &[coords.clone(), indices(&FACES)]).unwrap();
    let right = call(&session, "mesh_new", &[coords.clone(), indices(&FACES)]).unwrap();
    for operation in ["union", "intersection"] {
        let output = call(
            &session,
            "mesh_boolean",
            &[left.clone(), right.clone(), Value::String(operation.into())],
        )
        .unwrap();
        assert_eq!(one(&session, "mesh_is_closed", &output), Value::Bool(true));
        volume(&session, &output, 1.0 / 6.0);
    }
    let difference = call(
        &session,
        "mesh_boolean",
        &[
            left.clone(),
            right.clone(),
            Value::String("difference".into()),
        ],
    )
    .unwrap();
    assert_eq!(
        one(&session, "mesh_triangle_count", &difference),
        Value::int(0)
    );
    assert!(
        call(
            &session,
            "mesh_boolean",
            &[left.clone(), right, Value::String("invalid".into())],
        )
        .is_err()
    );
    let open = call(&session, "mesh_new", &[coords, indices(&FACES[..3])]).unwrap();
    assert!(
        call(
            &session,
            "mesh_boolean",
            &[left, open, Value::String("union".into())],
        )
        .is_err()
    );
}
