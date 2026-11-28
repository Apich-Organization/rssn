//! Measured planar adapters guard geometry, work, and the constraint chains.
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::TensorValue;
use rssn::kernels::meshgen::advancing_front_2d;
use rssn::kernels::meshgen::delaunay::constrained_delaunay_2d_mesh;

fn square() -> [[f64; 2]; 4] {
    [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]
}

#[test]
fn native_front_covers_a_unit_square_and_preserves_orientation() {
    let mesh = advancing_front_2d(&square(), 0.5).unwrap();
    let area: f64 = mesh
        .triangles
        .iter()
        .map(|&[i, j, k]| {
            let a = mesh.points.get(i).copied().unwrap();
            let b = mesh.points.get(j).copied().unwrap();
            let c = mesh.points.get(k).copied().unwrap();
            let cross = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
            assert!(cross > 0.0);
            0.5 * cross
        })
        .sum();
    assert!((area - 1.0).abs() < 1e-12);
}

#[test]
fn native_delaunay_keeps_four_split_constraint_edges() {
    let mesh = constrained_delaunay_2d_mesh(&square(), &[(0, 2), (1, 3)]).unwrap();
    assert_eq!(mesh.points.len(), 5);
    assert_eq!(mesh.triangles.len(), 4);
    for vertex in 0..4 {
        assert!(
            mesh.triangles
                .iter()
                .any(|triangle| triangle.contains(&vertex) && triangle.contains(&4))
        );
    }
}

#[test]
fn release_abi_preserves_length_and_limits_boundary_work() {
    let points = Value::Tensor(TensorValue {
        shape: vec![4, 2],
        elem: TypeExpr::Float.with_dim(DimExpr::base("Length")),
        data: square().into_iter().flatten().map(Value::Float).collect(),
    });
    let length = |v| {
        Value::Dyn(DynValue {
            value: Box::new(Value::Float(v)),
            dim: DimExpr::base("Length"),
        })
    };
    let ctx = CallCtx::new();
    let session = Session::new();
    assert_eq!(
        session.call_value(
            "num_meshgen_advancing_front_2d",
            &[points.clone(), length(0.001), length(1.0), Value::int(16)],
            &ctx
        ),
        Err(OpError::BudgetExceeded),
    );
    let result = session
        .call_value(
            "num_meshgen_advancing_front_2d",
            &[points, length(0.5), length(1.0), Value::int(256)],
            &ctx,
        )
        .unwrap();
    let Value::Record(result) = result else {
        panic!("expected mesh record")
    };
    let points = result
        .fields
        .iter()
        .find(|(name, _)| name == "points")
        .map(|(_, value)| value);
    assert!(
        matches!(points,Some(Value::Tensor(t)) if t.elem==TypeExpr::Float.with_dim(DimExpr::base("Length")))
    );
}
