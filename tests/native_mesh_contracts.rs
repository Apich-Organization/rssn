//! Typed mesh adapters and ordered triangle measurements.
use rssn::api::Session;
use rssn::execution::with_scoped_pool;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::overload::matches_type;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::TensorValue;
use rssn::kernels::meshgen::triangle_areas;

fn points<const N: usize>(rows: &[[f64; N]]) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![rows.len(), N],
        data: rows.iter().flatten().copied().map(Value::Float).collect(),
        elem: TypeExpr::Float.with_dim(DimExpr::base("Length")),
    })
}

fn indices<const N: usize>(rows: &[[usize; N]]) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![rows.len(), N],
        data: rows.iter().flatten().copied().map(Value::int).collect(),
        elem: TypeExpr::Int,
    })
}

fn check_measurements(
    value: Value,
    expected: &[f64],
    power: i64,
) {
    let Value::Tensor(result) = value else {
        panic!("expected a tensor of mesh measurements");
    };
    let target = TypeExpr::Float.with_dim(DimExpr::base("Length").pow(power));
    assert_eq!(result.shape, [expected.len()]);
    assert!(matches_type(&target, &result.elem).accepted());
    for (actual, expected) in result.data.iter().zip(expected) {
        let Value::Dyn(cell) = actual else {
            panic!("expected a measurement with a dimension witness");
        };
        assert_eq!(cell.value.as_f64(), Some(*expected));
        assert!(matches_type(&target, &TypeExpr::Float.with_dim(cell.dim.clone())).accepted());
    }
}

#[test]
fn mesh_measurements_keep_units_winding_and_cell_order() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let planar = points(&[[0.0, 0.0], [2.0, 0.0], [0.0, 3.0], [1.0, 0.0]]);
    let measurements = session
        .call_value(
            "num_meshgen_triangle_areas",
            &[planar.clone(), indices(&[[0, 3, 2], [0, 2, 1], [0, 1, 3]])],
            &ctx,
        )
        .unwrap();
    check_measurements(measurements, &[1.5, 3.0, 0.0], 2);
    let tessellation = session
        .call_value(
            "num_meshgen_delaunay_2d_incremental",
            &[points(&[[0.0, 0.0], [2.0, 0.0], [0.0, 3.0]])],
            &ctx,
        )
        .unwrap();
    let area = session
        .call_value("num_meshgen_triangle_areas", &[planar, tessellation], &ctx)
        .unwrap();
    check_measurements(area, &[3.0], 2);
    let solid = points(&[
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [0.0, 3.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -2.0],
    ]);
    let tetrahedra = indices(&[[0, 1, 2, 3], [0, 2, 1, 4]]);
    let volumes = session
        .call_value(
            "num_meshgen_tetra_volumes",
            &[solid, tetrahedra.clone()],
            &ctx,
        )
        .unwrap();
    check_measurements(volumes, &[1.0, 2.0], 3);
    let faces = session
        .call_value("num_meshgen_boundary_faces", &[tetrahedra], &ctx)
        .unwrap();
    assert_eq!(
        faces,
        indices(&[
            [1, 2, 3],
            [0, 3, 2],
            [0, 1, 3],
            [2, 1, 4],
            [0, 4, 1],
            [0, 2, 4]
        ])
    );
}

#[test]
fn mesh_contract_rejects_malformed_cells_and_preserves_empty_meshes() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let xy = points(&[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]);
    let area = "num_meshgen_triangle_areas";
    let angle = "num_meshgen_min_angle_deg";
    let faces = "num_meshgen_boundary_faces";
    assert!(matches!(
        session.call_value(area, &[xy.clone(), indices(&[[0, 1, 3]])], &ctx),
        Err(OpError::Domain(_))
    ));
    assert!(matches!(
        session.call_value(area, &[xy.clone(), indices(&[[0, 1]])], &ctx),
        Err(OpError::Domain(_))
    ));
    check_measurements(
        session
            .call_value(area, &[xy.clone(), indices::<3>(&[])], &ctx)
            .unwrap(),
        &[],
        2,
    );
    let measured_angle = session
        .call_value(angle, &[xy, indices(&[[0, 1, 2]])], &ctx)
        .unwrap();
    let Value::Float(degrees) = measured_angle else {
        panic!("expected a dimensionless angle in degrees");
    };
    assert!((degrees - 45.0).abs() < 1e-12);
    assert!(matches!(
        session.call_value(faces, &[indices(&[[0, 1, 1, 2]])], &ctx),
        Err(OpError::Domain(_))
    ));
    assert_eq!(
        session
            .call_value(faces, &[indices::<4>(&[])], &ctx)
            .unwrap(),
        indices::<3>(&[])
    );
}

#[test]
fn parallel_triangle_areas_preserve_order_and_bits() {
    let planar = [[0.0, 0.0], [2.0, 0.0], [0.0, 3.0], [1.0, 0.0]];
    let pattern = [[0, 1, 2], [0, 2, 3], [0, 1, 3], [0, 1, 7]];
    let cells: Vec<_> = pattern.into_iter().cycle().take(8_193).collect();
    let bit_pattern = |threads| {
        with_scoped_pool(Some(threads), || {
            triangle_areas(&planar, &cells)
                .into_iter()
                .map(f64::to_bits)
                .collect::<Vec<_>>()
        })
    };
    let serial = bit_pattern(1);
    assert_eq!(serial, bit_pattern(4));
    assert_eq!(&serial[..4], &[3.0_f64.to_bits(), 1.5_f64.to_bits(), 0, 0]);
}
